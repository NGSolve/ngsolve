#ifndef CUDA_NGBLA
#define CUDA_NGBLA

namespace ngbla
{
    using namespace ngs_cuda;
    
  template<> struct trivtrans<Dev<double>> { static constexpr bool value = true; };

    
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
      if (auto err = cudaMalloc((void**)&data, h*w*sizeof(T)))
        throw Exception("UnifiedVector allocation error, ec="+ToString(err));
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
      cudaFree(data);
    }
         
    void D2H (FlatMatrix<T> mat)
    {
      cout << "D2H, myshape = " << h << "x" << w << " other shape = " << mat.Height() << "x" << mat.Width() << endl;
      cudaMemcpy (mat.Data(), data, sizeof(T)*h*w, cudaMemcpyDeviceToHost);
    }
  };
    
    
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
    
}

#endif
