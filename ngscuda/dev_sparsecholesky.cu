#include <la.hpp>
#include "cuda_linalg.hpp"

using namespace ngla;



namespace ngla
{
  using namespace ngs_cuda;


  class DevSparseCholesky : public DevMatrix
  {
    double h, w;
    Array<Dev<SparseCholeskyTM<double>::MicroTask>> microtasks;
    DevTable<int> dependency;
    Array<int> host_incomingdep;
  public:
    DevSparseCholesky(const SparseCholeskyTM<double> & mat);
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    int VHeight() const override { return h; }
    int VWidth() const override { return w; }
  };



  auto a = []()->int
  {
    BaseMatrix::RegisterDeviceMatrixCreator
    (typeid(SparseCholesky<double>),
     [] (const BaseMatrix & bmat) -> shared_ptr<BaseMatrix>
    {
      auto & mat = dynamic_cast<const SparseCholeskyTM<double>&>(bmat);
      return make_shared<DevSparseCholesky>(mat);
    });
  } ();


  
  
  /* *************** kernels for SpasreCholesky *********************** */

  __global__ void DeviceSparseCholeskySolveLKernel (FlatTable<int> dependency, 
                                                    FlatArray<Dev<int>> incomingdep, 
                                                    FlatVector<Dev<double>> v,
                                                    int & cnt)
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

          // do the work ....


          if (threadIdx.x == 0)
              while (atomicAdd((int*)&incomingdep[myjob], 0) > 0);
          __syncwarp();
          
          if (threadIdx.x == 0)
              // for (int d : dependency[myjob])
              for (int j = 0; j < dependency[myjob].Size(); j++)
                 atomicAdd((int*)&incomingdep[dependency[myjob][j]], -1);
       }
  }
  
  void DeviceSparseCholeskySolveL (const DevTable<int> & dependency, 
                                   FlatArray<Dev<int>> incomingdep, 
                                   FlatVector<Dev<double>> v)
  {
    Dev<int> * pcnt = Dev<int>::Malloc(1);
    pcnt->H2D(0);
    DeviceSparseCholeskySolveLKernel<<<512,dim3(32,8)>>> (dependency, incomingdep, v, *(int*)pcnt);
    cout << "shared counter: " << pcnt->D2H();
    Dev<int>::Free (pcnt);
  }

  



  DevSparseCholesky :: DevSparseCholesky(const SparseCholeskyTM<double> & mat)
    : h(mat.Height()), w(mat.Width()),
      microtasks(mat.GetMicroTasks()),
      dependency(mat.GetMicroDependency()),
      host_incomingdep(mat.GetMicroDependency().Size())
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
    Array<Dev<int>> incomingdep(host_incomingdep);
    DeviceSparseCholeskySolveL (dependency, incomingdep, ux.FVDev());
    cout << "kernel is back" << endl;
    if (synckernels) cudaDeviceSynchronize();
    t.Stop();
      
    uy.InvalidateHost();
  }



}







