#include "cuda_linalg.hpp"
#include <la.hpp>

using namespace ngla;



namespace ngla
{
  extern bool synckernels;
  using namespace ngs_cuda;

  using MicroTask = SparseCholeskyTM<double>::MicroTask;

  class DevSparseCholesky : public DevMatrix
  {
    size_t h, w, nused;
    Array<Dev<MicroTask>> microtasks;
    DevTable<int> micro_dependency, micro_dependency_trans;
    Array<int> host_incomingdep;

    // block i has dofs  [blocks[i], blocks[i+1])
    Array<Dev<int>> blocks; 

    // row-indices of non-zero entries
    // all row-indices within one block are identic, and stored just once
    Array<Dev<int>> rowindex2;

    // index-array to rowindex
    Array<Dev<size_t>> firstinrow_ri;

    // index-array to lfact
    Array<Dev<size_t>> firstinrow;

    // L-factor in compressed storage
    Array<Dev<double>> lfact;

    // diagonal 
    Array<Dev<double>> diag;

    // the reordering (original dofnr i -> order[i])
    Array<Dev<int>> order;

  public:
    DevSparseCholesky(const SparseCholeskyTM<double> & mat);
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    int VHeight() const override { return h; }
    int VWidth() const override { return w; }
  };



  static const auto a = []()
  {
    BaseMatrix::RegisterDeviceMatrixCreator
    (typeid(SparseCholesky<double>),
     [] (const BaseMatrix & bmat) -> shared_ptr<BaseMatrix>
    {
      auto & mat = dynamic_cast<const SparseCholeskyTM<double>&>(bmat);
      return make_shared<DevSparseCholesky>(mat);
    });
    return 0;
  } ();


  
  
  /* *************** kernels for SparseCholesky *********************** */

  __global__ void DeviceSparseCholeskyReorderKernel (FlatVector<Dev<double>> src, FlatVector<Dev<double>> dst, FlatArray<Dev<int>> order)
  {
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
    for (int i = tid; i < order.Size(); i += blockDim.x*gridDim.x)
      if (order[i] != -1)
        dst[order[i]] = src[i];
  }

  __global__ void DeviceSparseCholeskyReorderAddKernel (FlatVector<Dev<double>> src, FlatVector<Dev<double>> dst, FlatArray<Dev<int>> order, double s)
  {
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
    for (int i = tid; i < order.Size(); i += blockDim.x*gridDim.x)
      if (order[i] != -1)
        dst[i] += s*src[order[i]];
  }

  __global__ void DeviceSparseCholeskySolveLKernel (
        FlatTable<int> dependency, 
        FlatArray<Dev<int>> incomingdep, 
        FlatVector<Dev<double>> hy,
        int & cnt,
        FlatArray<Dev<MicroTask>> microtasks,
        FlatArray<Dev<int>> blocks,
        FlatArray<Dev<int>> rowindex2,
        FlatArray<Dev<size_t>> firstinrow_ri,
        FlatArray<Dev<size_t>> firstinrow,
        FlatArray<Dev<double>> lfact
      )
  {
    __shared__ int myjobs[16]; // max blockDim.y;   
    
    while (true)
       {
          if (threadIdx.x == 0)
             myjobs[threadIdx.y] = atomicAdd(&cnt, 1);
          __syncwarp();

          int myjob = myjobs[threadIdx.y];
          if (myjob >= dependency.Size())
              break;

          volatile int * n_deps = (int*)&incomingdep[myjob];
          if (threadIdx.x == 0)
            while(*n_deps);
          __syncwarp();

          // do the work for myjob....
          
          MicroTask task = microtasks[myjob];
          size_t blocknr = task.blocknr;
          auto range = T_Range<int>(blocks[blocknr], blocks[blocknr+1]);
          auto base = firstinrow_ri[range.First()] + range.Size()-1;
          auto ext_size =  firstinrow[range.First()+1]-firstinrow[range.First()] - range.Size()+1;
          auto extdofs = rowindex2.Range(base, base+ext_size);

          
          if ((task.type == MicroTask::L_BLOCK) || (task.type == MicroTask::LB_BLOCK))
            {
              if (threadIdx.x == 0)
                for (auto i : range)
                  {
                    size_t size = range.end()-i-1;
                    if (size == 0) continue;
                    
                    auto vlfact = lfact.Range(firstinrow[i], firstinrow[i]+size);
                    
                    double hyi = hy(i);
                    auto hyr = hy.Range(i+1, range.end());
                    for (size_t j = 0; j < hyr.Size(); j++)
                      hyr(j) -= vlfact[j] * hyi;
                  }
            }
          
          if ((task.type == MicroTask::B_BLOCK) || (task.type == MicroTask::LB_BLOCK))
            {
              if (extdofs.Size() != 0)
                if (threadIdx.x == 0)
                  {
                    auto myr = Range(extdofs).Split (task.bblock, task.nbblocks);
                    auto my_extdofs = extdofs.Range(myr);
                    
                    for (size_t j = 0; j < my_extdofs.Size(); j++)
                      {
                        double temp = 0.0;
                        for (auto i : range)
                          {
                            size_t first = firstinrow[i] + range.end()-i-1;
                            auto ext_lfact = lfact.Range(first, all_extdofs.Size());
                            
                            temp += Trans(ext_lfact[myr.begin()+j]) * hy(i);
                          }
                        atomicAdd ((double*)(&hy(my_extdofs[j])), -temp);
                      }
                  }
            }
          
          // myjob is done
          
          if (threadIdx.x == 0)
              for (int d : dependency[myjob])
                 atomicAdd((int*)&incomingdep[d], -1);
       }
  }
  


  DevSparseCholesky :: DevSparseCholesky(const SparseCholeskyTM<double> & mat)
    : h(mat.Height()), w(mat.Width()),
      microtasks(mat.GetMicroTasks()),
      micro_dependency(mat.GetMicroDependency()),
      micro_dependency_trans(mat.GetMicroDependencyTranspose()),
      host_incomingdep(mat.GetMicroDependency().Size()),
      blocks(mat.GetBlocks()),
      rowindex2(mat.GetRowIndex2()),
      firstinrow_ri(mat.GetFirstInRowRI()),
      firstinrow(mat.GetFirstInRow()),
      lfact(mat.GetLFact()),
      diag(mat.GetDiag()),
      order(mat.GetOrder())
  {
    auto hostdep = mat.GetMicroDependency();
  
    host_incomingdep = 0;
  
    bool directional = true;
    for (int i = 0; i < hostdep.Size(); i++)
      for (int d : hostdep[i])
        {
          if (d <= i) directional = false;
          host_incomingdep[d]++;
        }
    cout << "directional = " << (directional? "yes" : "no") << endl;
  }



  
  void DevSparseCholesky ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
  
    static Timer t("DevSparseCholesky::MultAdd");
    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);
    ux.UpdateDevice();
    uy.UpdateDevice();
    if (synckernels) cudaDeviceSynchronize();
    t.Start();
      
    cout << "MultAdd in DevSpasreCholesky" << endl;

    DevStackArray<double> mem_hx(x.Size());
    DevStackArray<double> mem_hy(y.Size());
    FlatVector<Dev<double>> hx(x.Size(), mem_hx.Data());
    FlatVector<Dev<double>> hy(y.Size(), mem_hy.Data());

    DeviceSparseCholeskyReorderKernel<<<512,256>>> (ux.FVDev(), hx, order);

    cout << "hx[:10] = " << endl << D2H(hx.Range(10)) << endl;

    Array<Dev<int>> incomingdep(host_incomingdep);

    Dev<int> * pcnt = Dev<int>::Malloc(1);
    pcnt->H2D(0);
    DeviceSparseCholeskySolveLKernel<<<512,dim3(32,8)>>>
      (
       micro_dependency, incomingdep, hx, *(int*)pcnt,
       microtasks, blocks, rowindex2, firstinrow_ri, firstinrow, lfact
       );
    Dev<int>::Free (pcnt);


    
    
    DeviceSparseCholeskyReorderAddKernel<<<512,256>>> (hy, uy.FVDev(), order, s);

    if (synckernels) cudaDeviceSynchronize();
    t.Stop();
      
    uy.InvalidateHost();
  }



}







