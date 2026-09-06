#ifndef CUDA_NGBLA
#define CUDA_NGBLA

#include <cuda_runtime.h>

#include <vector.hpp>
#include <matrix.hpp>

#include "cuda_ngstd.hpp"
#include "linalg_kernels.hpp"



namespace ngbla
{
  using namespace ngs_cuda;
    
  // template<> struct trivtrans<Dev<double>> { static constexpr bool value = true; };
  template<> struct is_scalar_type<Dev<double>> { static constexpr bool value = true; };
  

  template <typename T>  
  class Vector<Dev<T>> : public FlatVector<Dev<T>>
  { 
    using FlatVector<Dev<T>>::Size;
    using FlatVector<Dev<T>>::Data;

  public:
    Vector (Vector&) = delete;
    Vector (Vector&&v2)
      : FlatVector<Dev<T>>(v2.Size(), v2.Data())
    {
      v2.layout = { nullptr, 0 };      
    }
         
    Vector (size_t asize)
      : FlatVector<Dev<T>>(asize, Dev<T>::Malloc(asize)) { ; }
         
    Vector (FlatVector<T> vec)
      : FlatVector<Dev<T>>(vec.Size(), Dev<T>::Malloc(vec.Size()))
    {
      H2D(vec);
    }
    
    ~Vector()
    {
      Dev<T>::Free(Data());
    }


    template<typename TB>
    Vector & operator= (const Expr<TB> & v)
    {
      MatExpr<FlatVector<Dev<T>> >::operator= (v);
      return *this;
    }

    
    void D2H (FlatVector<T> vec) const
    {
      cudaMemcpy (vec.Data(), Data(), sizeof(T)*Size(), cudaMemcpyDeviceToHost);
    }

    void H2D (FlatVector<T> vec)
    {
      cudaMemcpy (Data(), vec.Data(), sizeof(T)*Size(), cudaMemcpyHostToDevice);
    }

    Vector<T> D2H() const
    {
      Vector<T> vh(Size());
      D2H (vh);
      return vh;
    }
  };
  
  inline Vector<double> D2H (FlatVector<Dev<double>> dvec)
  {
    Vector<double> hvec(dvec.Size());
    cudaMemcpy (hvec.Data(), dvec.Data(), sizeof(double)*hvec.Size(), cudaMemcpyDeviceToHost);
    return hvec;
  }

#ifdef OLDOLD
#ifdef __CUDACC__  
  template <typename TS, typename TD>
  __global__ void kernel_Assign (size_t n,  TD pod_dst, TS pod_src)
  {
    auto dst = *pod_dst;
    auto src = *pod_src;
    
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
    for (int i = tid; i < n; i += blockDim.x*gridDim.x)
      if (i < 5)
        dst(i) = src.S();
      else
        dst(i) = src.A()(i);
  }
#endif

  template <typename T>
  class AsPOD
  {
    std::array<char, sizeof(T)> data;
  public:
    AsPOD(const AsPOD&) = default;
    INLINE AsPOD (const T & adata)
    {
      char * pdata = (char*)(void*)&adata;
      for (int i = 0; i < sizeof(T); i++)
        data[i] = pdata[i];
    }

    INLINE const T & operator* () const
    {
      T * val = (T*)(void*)&data[0];
      return *val;
    }
    INLINE int operator[] (int i) const { return data[i]; }
  };
#endif
  
  

  template <typename TOP, typename T, typename TS, typename TDIST, typename TB>
  class assign_trait<TOP, VectorView<T,TS,TDIST>, TB,
                     enable_if_t < std::is_same<std::invoke_result_t<VectorView<T,TS,TDIST>,size_t>, Dev<double>&>::value, int>>
  {
  public:
    static INLINE VectorView<T,TS,TDIST> & Assign (MatExpr<VectorView<T,TS,TDIST>> & self, const Expr<TB> & v)
    {

#ifdef __CUDACC__ 
      
      ngs_cuda::DeviceParallelFor
        (self.Height(),
         [devself=self.Spec(), devv=v.Spec()] DEVICE_LAMBDA (auto tid) -> void
         {
           // devself(tid) = devv(tid);
           TOP()(devself(tid),devv(tid));
         });

#endif

      return self.Spec();
    }
  };    



  
    
  template <typename T>  
  class Matrix<Dev<T>> : public FlatMatrix<Dev<T>>
  { 
  public:
    using FlatMatrix<Dev<T>>::Height;
    using FlatMatrix<Dev<T>>::Width;
    using FlatMatrix<Dev<T>>::Data;

    Matrix (Matrix&) = delete;
    Matrix (Matrix&&) = default;
         
    Matrix (size_t h_, size_t w_)
      : FlatMatrix<Dev<T>>(h_, w_, Dev<T>::Malloc(h_*w_)) { ; }
         
    Matrix (FlatMatrix<T> mat)
      : FlatMatrix<Dev<T>>(mat.Height(), mat.Width(),
                           Dev<T>::Malloc(mat.Height()*mat.Width()))
    {
      H2D(mat);
    }
    
    ~Matrix()
    {
      Dev<T>::Free(Data());
    }
         
    void D2H (FlatMatrix<T> mat) const
    {
      cudaMemcpy (mat.Data(), Data(), sizeof(T)*Height()*Width(), cudaMemcpyDeviceToHost);
    }

    void H2D (FlatMatrix<T> mat)
    {
      cudaMemcpy (Data(), mat.Data(), sizeof(T)*Height()*Width(), cudaMemcpyHostToDevice);
    }

    Matrix<T> D2H() const
    {
      Matrix<T> mh(Height(), Width());
      D2H (mh);
      return mh;
    }
  };
  
  inline Matrix<double> D2H (SliceMatrix<Dev<double>> dmat)
  {
    Matrix<double> hmat(dmat.Height(), dmat.Width());
    for (size_t i = 0; i < hmat.Height(); i++)
      cudaMemcpy (&hmat(i,0), &dmat(i,0), sizeof(double)*hmat.Width(), cudaMemcpyDeviceToHost);
    return hmat;
  }
  
  inline Matrix<double,ColMajor> D2H (SliceMatrix<Dev<double>,ColMajor> dmat)
  {
    return Trans(D2H(Trans(dmat)));
  }

}

#endif
