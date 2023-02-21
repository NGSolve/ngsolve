#ifndef LINALG_KERNELS_HPP
#define LINALG_KERNELS_HPP


namespace ngs_cuda
{
  using namespace ngbla;



// own ngsolve cuda-kernels:
extern void SetScalar (double val, int n, double * dev_ptr);
extern void SetScalar (double val, ngbla::FlatVector<Dev<double>> vec);

extern void SetVector (double val, int n, double * x, double * y);
extern void MyDaxpy (double val, int n, double * x, double * y);


extern void MultDiagonal (int n, double * D, double * x, double * y);
extern void MultAddDiagonal (int n, double alpha, double * D, double * x, double * y);

// y = A * x
class MatVecData
{
public:
  SliceMatrix<Dev<double>> mat;
  // BareVector<Dev<double>> x, y;
  size_t offsetx, offsety;
  MatVecData() : mat(0,0,0,nullptr) /* , x(nullptr), y(nullptr) */ { ; }
};


  /* **************** BlockJacobi kernels ********************* */

  
extern void ManyMatVec (FlatArray<Dev<MatVecData>> matvecs, 
                        BareVector<Dev<double>> x, BareVector<Dev<double>> y); 


  class BlockJacobiCtr
  {
  public:
    SliceMatrix<Dev<double>> mat;
    Dev<int> * indices;
    BlockJacobiCtr() : mat(0,0,0,nullptr) { ; }
  };
  
  extern void DeviceBlockJacobi (double s, FlatArray<Dev<BlockJacobiCtr>> ctrs, 
                                 BareVector<Dev<double>> x, BareVector<Dev<double>> y); 


    
extern void ConstEBEKernelCopyIn (int numblocks, int bs, int * row_dnums, double * dev_ux, double * dev_hx);
extern void ConstEBEKernelCopyOut (int numblocks, int bs, int * col_dnums, double * dev_hy, double * dev_uy);
extern void ConstEBEKernelCopyInIdx (int numblocks, int * idx, int bs, int * row_dnums, double * dev_ux, double * dev_hx);
extern void ConstEBEKernelCopyOutIdx (int numblocks, int * idx, int bs, int * col_dnums, double * dev_hy, double * dev_uy);

  extern void DevBlockDiagonalMatrixSoAMultAddVecs (double s, int size, double * a, double * b, double * res);

  // for (i,j,k) in indices:
  //    res.Row(k) += s * a.Row(i) * b.Row(j)
  extern void DevBlockDiagonalMatrixSoAMultAddVecs (double s, FlatArray<Dev<int>>, 
                                                    SliceMatrix<Dev<double>> a, 
                                                    SliceMatrix<Dev<double>> b,
                                                    SliceMatrix<Dev<double>> res);


extern void DevProjectorMultAdd (double s, size_t size, const double * a, double * b, const unsigned char * bits, bool keep_values);
extern void DevProjectorProject (size_t size, double * a, const unsigned char * bits, bool keep_values);
}


#endif
