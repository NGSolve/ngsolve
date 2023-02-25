#include "cuda_linalg.hpp"
#include "cuda_profiler.hpp"
#include <la.hpp>

using namespace ngla;


namespace ngla
{
  class BlockJacobiCtr
  {
  public:
    SliceMatrix<Dev<double>> mat;
    Dev<int> * indices;
    BlockJacobiCtr() : mat(0,0,0,nullptr) { ; }
  };

  
  
  class DevBlockJacobiMatrix : public DevMatrix
  {
    double h, w;
    Array<Dev<int>> indices;
    Array<Dev<double>> matrices;
    Array<Dev<BlockJacobiCtr>> ctrstructs;
  public:
    DevBlockJacobiMatrix (const BlockJacobiPrecond<double> & mat);
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    // void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    int VHeight() const override { return h; }
    int VWidth() const override { return w; }
  };

  static const auto a = []()
  {
    BaseMatrix::RegisterDeviceMatrixCreator
    (typeid(BlockJacobiPrecond<double>),
     [] (const BaseMatrix & bmat) -> shared_ptr<BaseMatrix>
    {
      auto & mat = dynamic_cast<const BlockJacobiPrecond<double>&>(bmat);
      return make_shared<DevBlockJacobiMatrix>(mat);
    });
    return 0;
  } ();



  __global__ void BlockJacobiKernel (double s, FlatArray<Dev<BlockJacobiCtr>> ctrs, 
                                     BareVector<Dev<double>> x, BareVector<Dev<double>> y)
  {
    for (int i = blockIdx.x*blockDim.y+threadIdx.y; i < ctrs.Size(); i += gridDim.x*blockDim.y)
      {
        BlockJacobiCtr mv = ctrs[i];
        size_t h = mv.mat.Height();
        size_t w = mv.mat.Width();
     
        for (int r = threadIdx.x; r < h; r += blockDim.x)
          {
            double sum = 0;
            for (int c = 0; c < w; c++)
              sum += mv.mat(r,c) * x(mv.indices[c]);
            atomicAdd((double*)&y(mv.indices[r]), s*sum);
          }
      }
  }
    
  
  void DeviceBlockJacobi (double s, FlatArray<Dev<BlockJacobiCtr>> ctrs, 
                          BareVector<Dev<double>> x, BareVector<Dev<double>> y)
  {
    BlockJacobiKernel<<<512,dim3(16,16)>>> (s, ctrs, x, y);
  }
  


  DevBlockJacobiMatrix :: DevBlockJacobiMatrix (const BlockJacobiPrecond<double> & mat)
    : h(mat.Height()), w(mat.Width()),
      matrices(mat.MatrixData()), indices(mat.GetBlockTable()->AsArray())
  {
    const Array<FlatMatrix<double>> & inverses = mat.GetInverses();
        
    Array<BlockJacobiCtr> hostctrs(inverses.Size());
    Dev<double> * matptr = matrices.Data();
    Dev<int> * indexptr = indices.Data();
    for (size_t i = 0; i < inverses.Size(); i++)
    {
      size_t s = inverses[i].Height();
      new (&hostctrs[i].mat) SliceMatrix<Dev<double>> (s, s, s, matptr);
      hostctrs[i].indices = indexptr;
      matptr += s*s;
      indexptr += s;
    }
          
    ctrstructs = Array<Dev<BlockJacobiCtr>> (hostctrs);
  }


  void DevBlockJacobiMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DevBlockJacobi::MultAdd");
    UnifiedVectorWrapper ux(x);
    UnifiedVectorWrapper uy(y);
    ux.UpdateDevice();
    uy.UpdateDevice();
    if (synckernels) cudaDeviceSynchronize();
    t.Start();
      
    DeviceBlockJacobi (s, ctrstructs, ux.FVDev(), uy.FVDev());
      
    if (synckernels) cudaDeviceSynchronize();
    t.Stop();
      
    uy.InvalidateHost();
  }

  
}
