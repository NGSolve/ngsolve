#ifndef CUDA_NGBLA
#define CUDA_NGBLA

#include "cuda_ngstd.hpp"
#include "linalg_kernels.hpp"

namespace ngbla
{
  using namespace ngs_cuda;
    
  // template<> struct trivtrans<Dev<double>> { static constexpr bool value = true; };
  template<> struct is_scalar_type<Dev<double>> { static constexpr bool value = true; };
  
  
    
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
      : FlatMatrix<Dev<T>>(h_, w_, nullptr)
    {
      this->data = Dev<T>::Malloc(h*w);
    }
         
    Matrix (FlatMatrix<T> mat)
    {
      h = mat.Height();
      w = mat.Width();
      if (auto err = cudaMalloc((void**)&data, h*w*sizeof(T)))
        throw Exception("UnifiedVector allocation error, ec="+ToString(err));

      cudaMemcpy (data, mat.Data(), sizeof(T)*h*w, cudaMemcpyHostToDevice);
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
    // double alpha = 1;
    // double beta = 0;
    cublasStatus_t stat = cublasDgemm(ngla::Get_CuBlas_Handle(), 
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
            typename enable_if<IsConvertibleToSliceMatrix<TA,Dev<double>>(),int>::type = 0,
            typename enable_if<IsConvertibleToSliceMatrix<TB,Dev<double>>(),int>::type = 0,
            typename enable_if<IsConvertibleToSliceMatrix<TC,Dev<double>>(),int>::type = 0>
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
    static INLINE T & Assign (MatExpr<T> & self, const Expr<ScaleExpr<MultExpr<TB1,TB2>,double>> & v)
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
    static INLINE T & Assign (MatExpr<T> & self, const Expr<MultExpr<ScaleExpr<TB1,double>,TB2>> & v)
    {
      auto res = self.View();

      double alpha = is_same_v<TOP,typename MatExpr<T>::AsSub> ? -1 : 1;
      double beta = is_same_v<TOP,typename MatExpr<T>::As> ? 0 : 1;
      
      alpha *= v.View().A().S();

      MultMatMat (v.View().A().A(), v.View().B(), self.ViewRW(), alpha, beta);
      return self.Spec();
    }
  };    
    
    /*
    // not working yet (is scalar val on host or device ?)
  template <>
  inline const FlatVector<Dev<double>> & FlatVector<Dev<double>> :: operator= (const Dev<double> & val) const
  {
      SetScalar (val, *this);
      return *this;
  }
  */
      
}

#endif
