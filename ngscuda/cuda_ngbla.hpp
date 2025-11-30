#ifndef CUDA_NGBLA
#define CUDA_NGBLA

#include <cuda_runtime.h>
#include <cublas_v2.h>

#include "cuda_ngstd.hpp"
#include "linalg_kernels.hpp"


namespace ngla
{
  cublasHandle_t Get_CuBlas_Handle ();
}


namespace ngbla
{
  using namespace ngs_cuda;
    
  // template<> struct trivtrans<Dev<double>> { static constexpr bool value = true; };
  template<> struct is_scalar_type<Dev<double>> { static constexpr bool value = true; };
  

  template <typename T>  
  class Vector<Dev<T>> : public FlatVector<Dev<T>>
  { 
    using FlatVector<Dev<T>>::size;
    using FlatVector<Dev<T>>::data;

  public:
    Vector (Vector&) = delete;
    Vector (Vector&&v2)
      : FlatVector<Dev<T>>(v2.Size(), v2.Data())
    {
      v2.data = nullptr;
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
      Dev<T>::Free(data);
    }

    void D2H (FlatVector<T> vec) const
    {
      cudaMemcpy (vec.Data(), data, sizeof(T)*size, cudaMemcpyDeviceToHost);
    }

    void H2D (FlatVector<T> vec)
    {
      cudaMemcpy (data, vec.Data(), sizeof(T)*size, cudaMemcpyHostToDevice);
    }

    Vector<T> D2H() const
    {
      Vector<T> vh(size);
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
  

  template <typename TOP, typename T, typename TS, typename TDIST, typename TB>
  class assign_trait<TOP, VectorView<T,TS,TDIST>, TB,
                     enable_if_t < std::is_same<std::invoke_result_t<VectorView<T,TS,TDIST>,size_t>, Dev<double>&>::value, int>>
  {
  public:
    static INLINE VectorView<T,TS,TDIST> & Assign (MatExpr<VectorView<T,TS,TDIST>> & self, const Expr<TB> & v)
    {

#if not defined(__CUDA_ARCH__)
      cout << "gpu expr template" << endl;
      cout << "res.size = " << self.Height() << "x" << self.Width() << endl;
      cout << "expr.size = " << v.Spec().Height() << "x" << v.Spec().Width() << endl;
      cout << "type = " << typeid(v.Spec()).name() << endl;
      cout << "expr = "; v.Spec().Dump(cout); cout << endl;
      cout << "sizeof(v) = " << sizeof(v.Spec()) << endl;
      cout << "host, S = " << v.Spec().S () << endl;
      cout << "triv copy = " << std::is_trivially_copyable<TB>::value << endl;
      cout << "triv copy dst = " << std::is_trivially_copyable<VectorView<T,TS,TDIST>>::value << endl;      
      auto vspec = v.Spec();
      int * ip = (int*)(void*)(&vspec);
      for (int i = 0; i < 6; i++)
        cout << "int[" << i << "] = " << ip[i] << endl;
      double * dp = (double*)(void*)(&vspec);
      for (int i = 0; i < 3; i++)
        cout << "double[" << i << "] = " << dp[i] << endl;
#endif      

      AsPOD src_pod{v.Spec()};
      AsPOD dst_pod{self.Spec()};
      
#if not defined(__CUDA_ARCH__)      
      for (int i = 0; i < sizeof(v.Spec()); i++)
        cout << "data[" << i << "] = " << src_pod[i] << endl;
#endif
      
      /*
      DeviceParallelFor
        (diag.Size(),
         [ddiag=diag.DevData(), dx=ux.DevData(), dy=uy.DevData()] DEVICE_LAMBDA (auto tid)
         {
           dy[tid] = ddiag[tid]*dx[tid];
         });
      */

#ifdef __CUDACC__
      ngs_cuda::DeviceParallelFor
        (self.Height(),
         [dst_pod, src_pod] DEVICE_LAMBDA (auto tid) -> void
         {
           TB v = *src_pod;
           auto devself = *dst_pod;
           // devself(tid) = v.S();
           // devself(tid) = src_pod[tid];
           // devself(tid) = 42; // devv(tid);
           // devself(tid) = v(tid);
           // devself(tid) = v.S();
           // devself(tid) = sizeof(devv);
           double * dp = (double*)(void*)(&v);
           /*
           if (tid < 3)
             devself(tid) = dp[tid];
           else
             // devself(tid) = sizeof(v.A());
             */
           if (tid == 0)
             devself(tid) = (char*)(&v.s) - (char*)(&v);
           // devself(tid) = offsetof(TB, a);
           else if (tid == 1)
             devself(tid) = offsetof(TB, s);
           else if (tid == 2)
             devself(tid) = sizeof(v);
           else if (tid == 3)
             devself(tid) = sizeof(v.a);
           else if (tid == 4)
             devself(tid) = sizeof(v.s);
             // devself(tid) = dp[2]*v.A()(tid);
    });
        
      /*
      ngs_cuda::DeviceParallelFor
        (self.Height(),
         [devself=self.Spec(), devv=v.Spec()] DEVICE_LAMBDA (auto tid) -> void
         {
           // devself(tid) = 42; // devv(tid);
           // devself(tid) = devv(tid);
           // devself(tid) = devv.S();
           // devself(tid) = sizeof(devv);
           double * dp = (double*)(void*)(&devv);           
           int * ip = (int*)(void*)(&devv);
           if (tid < 6)
             devself(tid) = ip[tid];
           else if (tid < 9)
             devself(tid) = dp[tid-6];             
           else
             devself(tid) = -1;
         });
      */


      // kernel_Assign<<<self.Height()/256+1,256>>> (self.Height(), AsPOD(self.Spec()), AsPOD(v.Spec()));

      
#endif
      
      return self.Spec();
    }
  };    



  
    
  template <typename T>  
  class Matrix<Dev<T>> : public FlatMatrix<Dev<T>>
  { 
    using FlatMatrix<Dev<T>>::h;
    using FlatMatrix<Dev<T>>::w;
    using FlatMatrix<Dev<T>>::data;

  public:
    Matrix (Matrix&) = delete;
    Matrix (Matrix&&m2)
      : FlatMatrix<Dev<T>>(m2.Height(), m2.Width(), m2.Data())
    {
      m2.data = nullptr;
    }
         
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
      Dev<T>::Free(data);
    }
         
    void D2H (FlatMatrix<T> mat) const
    {
      cudaMemcpy (mat.Data(), data, sizeof(T)*h*w, cudaMemcpyDeviceToHost);
    }

    void H2D (FlatMatrix<T> mat)
    {
      cudaMemcpy (data, mat.Data(), sizeof(T)*h*w, cudaMemcpyHostToDevice);
    }

    Matrix<T> D2H() const
    {
      Matrix<T> mh(h, w);
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
   
  
  
    
  template <ORDERING ORDA, ORDERING ORDB>  
  void CudaMultMatMat2 (SliceMatrix<Dev<double>, ORDA> a, SliceMatrix<Dev<double>,ORDB> b, 
                        SliceMatrix<Dev<double>, ORDERING::ColMajor> c,
                       double alpha, double beta)
  {
    cublasStatus_t stat =
      cublasDgemm(ngla::Get_CuBlas_Handle(), 
                  ORDA==ORDERING::RowMajor ? CUBLAS_OP_T : CUBLAS_OP_N, 
                  ORDB==ORDERING::RowMajor ? CUBLAS_OP_T : CUBLAS_OP_N, 
                  c.Height(), c.Width(), a.Width(),
                  &alpha, (double*)a.Data(), a.Dist(), (double*)b.Data(), b.Dist(),
                  &beta, (double*)c.Data(), c.Dist());
  }
  
  template <ORDERING ORDA, ORDERING ORDB>  
  void CudaMultMatMat2 (SliceMatrix<Dev<double>, ORDA> a, SliceMatrix<Dev<double>,ORDB> b, 
                        SliceMatrix<Dev<double>, ORDERING::RowMajor> c,
                       double alpha, double beta)
  {
    CudaMultMatMat2 (Trans(b), Trans(a), Trans(c), alpha, beta);
  }
   
    
  template <typename TA, typename TB, typename TC,
            enable_if_t<IsConvertibleToSliceMatrix<TA,Dev<double>>(),int> = 0,
            enable_if_t<IsConvertibleToSliceMatrix<TB,Dev<double>>(),int> = 0,
            enable_if_t<IsConvertibleToSliceMatrix<TC,Dev<double>>(),int> = 0>
  void MultMatMat (const TA & a, const TB & b, const TC & c, double alpha=1, double beta=0)
  {
    CudaMultMatMat2(make_SliceMatrix(a), make_SliceMatrix(b), make_SliceMatrix(c), alpha, beta);
  }
 

  template <typename TOP, typename T, typename TB1, typename TB2>
  class assign_trait<TOP, T, MultExpr<TB1,TB2>,
                     enable_if_t<IsConvertibleToSliceMatrix<T,Dev<double>>(),int>>
  {
  public:
    static INLINE T & Assign (MatExpr<T> & self, const Expr<MultExpr<TB1,TB2>> & v)
    {
      auto res = self.View();

      double alpha = std::is_same_v<TOP,typename MatExpr<T>::AsSub> ? -1 : 1;
      double beta = std::is_same_v<TOP,typename MatExpr<T>::As> ? 0 : 1;

      MultMatMat (v.Spec().A(), v.Spec().B(), self.Spec(), alpha, beta);
      return self.Spec();
    }
  };    
    
  template <typename TOP, typename T, typename TB1, typename TB2>
  class assign_trait<TOP, T, ScaleExpr<MultExpr<TB1,TB2>,double>,
                     enable_if_t<IsConvertibleToSliceMatrix<T,Dev<double>>(),int>>
  {
  public:
    static inline T & Assign (MatExpr<T> & self, const Expr<ScaleExpr<MultExpr<TB1,TB2>,double>> & v)
    {
      auto res = self.View();

      double alpha = is_same_v<TOP,typename MatExpr<T>::AsSub> ? -1 : 1;
      double beta = is_same_v<TOP,typename MatExpr<T>::As> ? 0 : 1;
      
      alpha *= v.View().S();

      MultMatMat (v.View().A().A(), v.View().A().B(), self.ViewRW(), alpha, beta);
      return self.Spec();
    }
  };    
    
  template <typename TOP, typename T, typename TB1, typename TB2>
  class assign_trait<TOP, T, MultExpr<ScaleExpr<TB1,double>,TB2>,
                     enable_if_t<IsConvertibleToSliceMatrix<T,Dev<double>>(),int>>
  {
  public:
    static inline T & Assign (MatExpr<T> & self, const Expr<MultExpr<ScaleExpr<TB1,double>,TB2>> & v)
    {
      auto res = self.View();

      double alpha = is_same_v<TOP,typename MatExpr<T>::AsSub> ? -1 : 1;
      double beta = is_same_v<TOP,typename MatExpr<T>::As> ? 0 : 1;
      
      alpha *= v.View().A().S();

      MultMatMat (v.View().A().A(), v.View().B(), self.ViewRW(), alpha, beta);
      return self.Spec();
    }
  };    
    
      
}

#endif
