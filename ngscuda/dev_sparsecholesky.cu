#include "cuda_linalg.hpp"
#include "cuda_profiler.hpp"
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
    Array<int> host_incomingdep_trans;

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

  __global__ void DeviceSparseCholeskyMultDiagKernel 
          (FlatArray<Dev<double>> diag, FlatVector<Dev<double>> vec)
  {
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
    for (int i = tid; i < diag.Size(); i += blockDim.x*gridDim.x)
      vec(i) *= diag[i];
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
    DeviceBlockRegionTracer brt(gridDim.x*blockDim.y, gridDim.x*threadIdx.y + blockIdx.x, threadIdx.x);
    
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

     __threadfence();
         
          if ((task.type == MicroTask::L_BLOCK) || (task.type == MicroTask::LB_BLOCK))
            {
              DeviceRegionTracer rt(brt, 0, task.blocknr);
               for (int i = blocks[blocknr]; i < blocks[blocknr+1]-1; i++)
                  {
                    size_t size = range.end()-i-1;
                    
                    auto vlfact = lfact.Range(firstinrow[i]-i-1, firstinrow[i]+size);
                    
                    __threadfence_block();   

                    double hyi = hy(i);
                    for (int j = threadIdx.x+blocks[blocknr]; j < range.end(); j += blockDim.x)
                       if (j > i)
                          hy(j) -= vlfact[j] * hyi;
                    __threadfence_block();   
                  }
            }
          
          if ((task.type == MicroTask::B_BLOCK) || (task.type == MicroTask::LB_BLOCK))
            {
              if (extdofs.Size() != 0)
                  {
                    DeviceRegionTracer rt(brt, 1, task.blocknr);
                    auto myr = Range(extdofs).Split (task.bblock, task.nbblocks);
                    auto my_extdofs = extdofs.Range(myr);
                    
                    for (int j = threadIdx.x; j < my_extdofs.Size(); j+=blockDim.x)
                      {
                        double temp = 0.0;
                        for (auto i : range)
                          {
                            size_t first = firstinrow[i] + range.end()-i-1;
                            auto ext_lfact = lfact.Range(first, extdofs.Size());
                            
                            temp += ext_lfact[myr.begin()+j] * hy(i);
                          }
                        atomicAdd ((double*)(&hy(my_extdofs[j])), -temp);
                      }
                  }
            }
                      
          // myjob is done
          __syncwarp();
         __threadfence();
         
          if (threadIdx.x == 0)
              for (int d : dependency[myjob])
                 atomicAdd((int*)&incomingdep[d], -1);
       }
  }
  


  __global__ void DeviceSparseCholeskySolveLTransKernel (
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

          int myjob = dependency.Size()-1 - myjobs[threadIdx.y];
          if (myjob < 0)
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
    __threadfence();
          // TODO: needs transpose:
          if ((task.type == MicroTask::B_BLOCK) || (task.type == MicroTask::LB_BLOCK))
            {
              if (extdofs.Size() != 0)
                  {
                    auto myr = Range(extdofs).Split (task.bblock, task.nbblocks);
                    auto my_extdofs = extdofs.Range(myr);
                                          
                    // for (auto i : range)
                    for (int i = range.First()+threadIdx.x; i < range.Next(); i += blockDim.x)
                        {
                            size_t first = firstinrow[i] + range.end()-i-1;
                            auto ext_lfact = lfact.Range(first, first+extdofs.Size());
                            
                            double val = 0.0;
                            for (auto j : Range(my_extdofs))
                                val += ext_lfact[myr.begin()+j] * hy(my_extdofs[j]);
                                
                            atomicAdd ((double*)(&hy(i)), -val);
                        }
                  }
                
            }

          //       
          if ((task.type == MicroTask::L_BLOCK) || (task.type == MicroTask::LB_BLOCK))
            {
            /*
              if (threadIdx.x == 0)
                {
                if (range.Size() > 0) // for case [0,0)
                  for (size_t i = range.end()-1; i-- > range.begin(); )
                    {
                        size_t size = range.end()-i-1;
                        if (size == 0) continue;
                        FlatVector<Dev<double>> vlfact(size, &lfact[firstinrow[i]]);
                        auto hyr = hy.Range(i+1, range.end());

                        double hyi = hy(i);
                        for (size_t j = 0; j < vlfact.Size(); j++)
                           hyi -= vlfact(j) * hyr(j);
                        hy(i) = hyi;
                    }
                }
                */
                /*
              if (threadIdx.x == 0)
                {
                if (range.Size() > 0) // for case [0,0)
                  for (int i = range.end()-2; i >= range.begin(); i--)
                    {
                        size_t size = range.end()-i-1;
                        auto vlfact = lfact.Data()+firstinrow[i]; 
                        
                        auto hyr = hy.Range(i+1, range.end());

                        double hyi = hy(i);
                        for (size_t j = 0; j < size; j++)
                           hyi -= vlfact[j] * hyr(j);
                        hy(i) = hyi;
                    }
                }
                */
                
                
                /*
                 if (threadIdx.x == 0)
                {
                  auto hhy = hy.Range(range);
                  
                   for (int j = range.Size()-1; j >= 0; j--)
                      for (int i = 0; i < j; i++)
                      {
                        auto vlfact = lfact.Data()+firstinrow[range.First()+i]-i-1; 
                        hhy(i) -= vlfact[j] * hhy(j);                           
                      }   
                }
                */
                
                
                 auto hhy = hy.Range(range);
                 for (int j = range.Size()-1; j >= 0; j--)
                    {
                      for (int i = threadIdx.x; i < j; i += blockDim.x)
                        {
                          auto vlfact = lfact.Data()+firstinrow[range.First()+i]-i-1; 
                          hhy(i) -= vlfact[j] * hhy(j);                           
                        }   
                     __threadfence_block();   
                    }
            }
          

           __syncwarp();
    __threadfence();
         
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
      host_incomingdep_trans(mat.GetMicroDependency().Size()),      
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

    if (!directional) throw Exception("dependency-graph must be directional");
    
    for (int i = 0; i < hostdep.Size(); i++)
      host_incomingdep_trans[i] = hostdep[i].Size();
  }



  
  void DevSparseCholesky ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
  
    static Timer t("DevSparseCholesky::MultAdd");
    static Timer tL("DevSparseChol, L-fact");
    static Timer tLT("DevSparseChol, LT-fact");
    static Timer tR("DevSparseChol, Reorder");
    static Timer tRA("DevSparseChol, Reorder Add");
    static Timer tD("DevSparseChol, Diag");
    
    RegionTimer rt(t);
    
    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);
    ux.UpdateDevice();
    uy.UpdateDevice();
    if (synckernels) cudaDeviceSynchronize();
    
    DevStackArray<double> mem_hx(x.Size());
    FlatVector<Dev<double>> hx(x.Size(), mem_hx.Data());

    CudaRegionTimer rtR(tR);
    DeviceSparseCholeskyReorderKernel<<<512,256>>> (ux.FVDev(), hx, order);
    rtR.Stop();

    // cout << "reordered[:10] = " << endl << D2H(hx.Range(10)) << endl;

    // Array<Dev<int>> incomingdep(host_incomingdep);
    // Array<Dev<int>> incomingdep_trans(host_incomingdep_trans);    
    
    DevStackArray<int> incomingdep(host_incomingdep.Size());
    H2D(incomingdep, host_incomingdep);
    DevStackArray<int> incomingdep_trans(host_incomingdep_trans.Size());
    H2D(incomingdep_trans, host_incomingdep_trans);    


    CudaRegionTimer rtL(tL);
    Dev<int> * pcnt = Dev<int>::Malloc(2);
    pcnt->H2D(0);
    (pcnt+1)->H2D(0);

    DeviceSparseCholeskySolveLKernel<<<512,dim3(32,8)>>>
      (
       micro_dependency, incomingdep, hx, *(int*)pcnt,
       microtasks, blocks, rowindex2, firstinrow_ri, firstinrow, lfact
       );
    rtL.Stop();

    // cout << "SolveL[:10] = " << endl << D2H(hx.Range(10)) << endl;
    
    // TODO : Diag
    {
    CudaRegionTimer rtD(tD);
    DeviceSparseCholeskyMultDiagKernel<<<512,256>>> (diag, hx);
    }

    // cout << "SolveDL[:10] = " << endl << D2H(hx.Range(10)) << endl;
    // cout << "Norm DL = " << L2Norm(D2H(hx)) << endl;
    //  pcnt->H2D(0);

    {
    CudaRegionTimer rtLT(tLT);
    DeviceSparseCholeskySolveLTransKernel<<<512,dim3(32,8)>>>
      (
       micro_dependency_trans, incomingdep_trans, hx, *(int*)(pcnt+1),
       microtasks, blocks, rowindex2, firstinrow_ri, firstinrow, lfact
       );
    }

    // cout << "SolveLT[:10] = " << endl << D2H(hx.Range(10)) << endl;

    
    CudaRegionTimer rtRA(tRA);
    DeviceSparseCholeskyReorderAddKernel<<<512,256>>> (hx, uy.FVDev(), order, s);
    rtRA.Stop();

    Dev<int>::Free (pcnt);
    if (synckernels) cudaDeviceSynchronize();
      
    uy.InvalidateHost();
  }



}







