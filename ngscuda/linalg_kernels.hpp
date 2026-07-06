#ifndef LINALG_KERNELS_HPP
#define LINALG_KERNELS_HPP


namespace ngs_cuda
{
  using namespace ngbla;


 


// own ngsolve cuda-kernels:
extern void SetScalar (double val, int n, double * dev_ptr);
extern void SetScalar (double val, ngbla::FlatVector<Dev<double>> vec);

// y = val*x
extern void SetVector (double val, int n, Dev<double> * x, Dev<double> * y);
extern void MyDaxpy (double val, int n, Dev<double> * x, Dev<double> * y);

/*
extern void MultDiagonal (int n, double * D, double * x, double * y);
extern void MultAddDiagonal (int n, double alpha, double * D, double * x, double * y);
*/
    
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



    /*
extern void ConstEBEKernelCopyIn (int numblocks, int bs, int * row_dnums, double * dev_ux, double * dev_hx);
extern void ConstEBEKernelCopyOut (int numblocks, int bs, int * col_dnums, double * dev_hy, double * dev_uy);
*/

extern void ConstEBEKernelCopyInIdx (int numblocks, int * idx, int bs, int * row_dnums, double * dev_ux, double * dev_hx);
extern void ConstEBEKernelCopyOutIdx (int numblocks, int * idx, int bs, int * col_dnums, double * dev_hy, double * dev_uy);

// TFQMR scalar batch kernels — one kernel launch per scalar-update group between vector operarions
// Even step, after InnerProduct(r,v): compute alpha, neg_alpha, coeff
extern void TFQMREvenBatch1(double* rho, double* vtrstar, double* theta, double* eta,
                             double* alpha, double* neg_alpha, double* coeff);
// Even step, after InnerProduct(w,w): compute theta, c, tau, tau_sq, eta; copy rho->rho_last
extern void TFQMREvenTauBatch(double* wnorm_sq, double* tau_in, double* alpha_in, double* rho,
                               double* theta, double* c, double* tau_out, double* tau_sq,
                               double* eta, double* rho_last);
// Odd step, before d update: compute coeff from current theta, eta, alpha
extern void TFQMROddCoeff(double* theta, double* eta, double* alpha, double* coeff);
// Odd step, after InnerProduct(w,w): compute theta, c, tau, tau_sq, eta
extern void TFQMROddTauBatch(double* wnorm_sq, double* tau_in, double* alpha_in,
                              double* theta, double* c, double* tau_out, double* tau_sq, double* eta);
// Odd step, after InnerProduct(w,r): compute beta, beta_sq; update rho_last
extern void TFQMROddBeta(double* rho, double* rho_last, double* beta, double* beta_sq);
    
    /*
  extern void DevBlockDiagonalMatrixSoAMultAddVecs (double s, int size, double * a, double * b, double * res);

  // for (i,j,k) in indices:
  //    res.Row(k) += s * a.Row(i) * b.Row(j)
  extern void DevBlockDiagonalMatrixSoAMultAddVecs (double s, FlatArray<Dev<int>>, 
                                                    SliceMatrix<Dev<double>> a, 
                                                    SliceMatrix<Dev<double>> b,
                                                    SliceMatrix<Dev<double>> res);


extern void DevProjectorMultAdd (double s, size_t size, const double * a, double * b, const unsigned char * bits, bool keep_values);
extern void DevProjectorProject (size_t size, double * a, const unsigned char * bits, bool keep_values);
*/

}

// Convergence check for cudaGraphCondTypeWhile
extern void ConvergenceCheck(double* rz, double tol,
    cudaGraphConditionalHandle handle,
    int* iter_count, int maxsteps);

#endif
