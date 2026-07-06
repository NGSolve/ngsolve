#include <climits>
// #include <ngstd.hpp>
// #include <templates.hpp>
#include <ngs_stdcpp_include.hpp>
#include <vector.hpp>
#include <matrix.hpp>
// #include <bla.hpp>


#include "cuda_ngstd.hpp"

#include "linalg_kernels.hpp"


namespace ngs_cuda
{

// x = val
__global__ void SetScalarKernel (double val, int n, double * x)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < n; i += blockDim.x*gridDim.x)
    x[i] = val;
}

void SetScalar1 (double val, int n, double * x)
{
  static Timer t("CUDA::SetScalar");
  CudaRegionTimer rt(t);
  SetScalarKernel<<<512,256,0,ngs_cuda_stream>>> (val, n, x);
} 



void SetScalar (double val, int n, double * x)
{
  static Timer t("CUDA::SetScalar");
  CudaRegionTimer rt(t);

  auto lam = [val,n,x] DEVICE_LAMBDA (int tid) {
    x[tid] = val;
  };
  // CUDA_forall<<<512,256>>> (n, lam);
  DeviceParallelFor (n, lam);
}



__global__ void SetScalarKernelNew (double val, FlatVector<Dev<double>> vec)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < vec.Size(); i += blockDim.x*gridDim.x)
    vec(i) = val;
}

void SetScalar (double val, FlatVector<Dev<double>> vec)
{
  static Timer t("CUDA::SetScalar");
  CudaRegionTimer rt(t);
  SetScalar(val, vec.Size(), reinterpret_cast<double*>(vec.Data()));
  // SetScalarKernelNew<<<512,256>>> (val, vec);
} 





// y[i] = val * x[i]
__global__ void SetVectorKernel (double val, int n, Dev<double> * x, Dev<double> * y)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < n; i += blockDim.x*gridDim.x)
    y[i] = val * x[i];
}

void SetVector (double val, int n, Dev<double> * x, Dev<double> * y)
{
  static Timer t("CUDA::SetVector");
  CudaRegionTimer rt(t);
  SetVectorKernel<<<512,256,0,ngs_cuda_stream>>> (val, n, x, y);
} 


// y[i] += val * x[i]
__global__ void MyDaxpyKernel (double val, int n, double * x, double * y)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < n; i += blockDim.x*gridDim.x)
    y[i] += val * x[i];
}

void MyDaxpy1 (double val, int n, double * x, double * y)
{
// cout << "MyDaxpy" << endl;
  MyDaxpyKernel<<<512,256,0,ngs_cuda_stream>>> (val, n, x, y);
} 


void MyDaxpy (double val, int n, double * x, double * y)
{
// cout << "MyDaxpy2" << endl;
  DeviceParallelFor (n, 
    [val,x,y] DEVICE_LAMBDA (int tid) 
    {
       y[tid] += val*x[tid];
    });
}


// y = D * x
__global__ void MultDiagonalKernel (int n, double * D, double * x, double * y)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < n; i += blockDim.x*gridDim.x)
    y[i] = D[i] * x[i];
}

void MultDiagonal (int n, double * D, double * x, double * y)
{
  MultDiagonalKernel<<<512,256,0,ngs_cuda_stream>>> (n, D, x, y);
} 


// y += alpha D * x
__global__ void MultAddDiagonalKernel (int n, double alpha, double * D, double * x, double * y)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < n; i += blockDim.x*gridDim.x)
    y[i] += alpha * D[i] * x[i];
}

void MultAddDiagonal (int n, double alpha, double * D, double * x, double * y)
{
  MultAddDiagonalKernel<<<512,256,0,ngs_cuda_stream>>> (n, alpha, D, x, y);
} 


/* ***************** Many Mat-Vec kernels ******************** */


/*
// y = A * x
class MatVecData
{
    public:
  SliceMatrix<Dev<double>> mat;
  BareVector<Dev<double>> x, y;
    MatVecData() : mat(0,0,0,nullptr), x(nullptr), y(nullptr) { ; }
};
 */   
__global__ void ManyMatVecKernel (FlatArray<Dev<MatVecData>> matvecs, 
                        BareVector<Dev<double>> x, BareVector<Dev<double>> y)
{
  DeviceBlockRegionTracer brt(gridDim.x*blockDim.y, gridDim.x*threadIdx.y + blockIdx.x, threadIdx.x);
  for (int i = blockIdx.x*blockDim.y+threadIdx.y; i < matvecs.Size(); i += gridDim.x*blockDim.y)
  {
     DeviceRegionTracer rt(brt, 0, i);
     MatVecData mv = matvecs[i];
     size_t h = mv.mat.Height();
     size_t w = mv.mat.Width();
     
     auto myx = x.Range(mv.offsetx, mv.offsetx+w);
     auto myy = y.Range(mv.offsety, mv.offsety+h);
     
     for (int r = threadIdx.x; r < h; r += blockDim.x)
       {
          double sum = 0;
          for (int c = 0; c < w; c++)
            sum += mv.mat(r,c) * myx(c);
          myy(r) = sum;
          // atomicAdd((double*)&myy(r), sum);
       }
  }
}
    
void ManyMatVec (FlatArray<Dev<MatVecData>> matvecs, 
                        BareVector<Dev<double>> x, BareVector<Dev<double>> y)
{
  ManyMatVecKernel<<<512,dim3(16,16),0,ngs_cuda_stream>>> (matvecs, x, y);
}





/* ************** kernels for ConstantEBE Matrix ********************** */


__global__ void ConstEBEKernelCopyInKernel (int numblocks, int bs, int * row_dnums, double * dev_ux, double * dev_hx)
{
  int tid = threadIdx.x;
  for (int r = blockIdx.x*blockDim.y+threadIdx.y; r < numblocks; r+=gridDim.x*blockDim.y)
    for (int i = tid; i < bs; i += blockDim.x)
      dev_hx[r*bs+i] = dev_ux[row_dnums[r*bs+i]];
}

void ConstEBEKernelCopyIn (int numblocks, int bs, int * row_dnums, double * dev_ux, double * dev_hx)
{
   // ConstEBEKernelCopyInKernel<<<512,256>>> (numblocks, bs, row_dnums, dev_ux, dev_hx);
  ConstEBEKernelCopyInKernel<<<512,dim3(16,16),0,ngs_cuda_stream>>> (numblocks, bs, row_dnums, dev_ux, dev_hx);
}

__global__ void ConstEBEKernelCopyOutKernel (int numblocks, int bs, int *  col_dnums, double * dev_hy, double * dev_uy)
{
  int tid = threadIdx.x;

  for (int r = blockIdx.x*blockDim.y+threadIdx.y; r < numblocks; r+=gridDim.x*blockDim.y)  
    for (int i = tid; i < bs; i += blockDim.x)
      // dev_uy[col_dnums[r*bs+i]] += dev_hy[r*bs+i];
      atomicAdd((double*)dev_uy+col_dnums[r*bs+i], dev_hy[r*bs+i]);      
}

void ConstEBEKernelCopyOut (int numblocks, int bs, int * col_dnums, double * dev_hy, double * dev_uy)
{
  ConstEBEKernelCopyOutKernel<<<512,dim3(16,16),0,ngs_cuda_stream>>> (numblocks, bs, col_dnums, dev_hy, dev_uy);
}




__global__ void ConstEBEKernelCopyInIdxKernel (int numblocks, int * idx, int bs, int * row_dnums, double * dev_ux, double * dev_hx)
{
  int tid = threadIdx.x;
  // for (int r = blockIdx.x; r < numblocks; r+=gridDim.x)
  for (int r = blockIdx.x*blockDim.y+threadIdx.y; r < numblocks; r+=gridDim.x*blockDim.y)    
    for (int i = tid; i < bs; i += blockDim.x)
      dev_hx[r*bs+i] = dev_ux[row_dnums[idx[r]*bs+i]];
}

void ConstEBEKernelCopyInIdx (int numblocks, int * idx, int bs, int * row_dnums, double * dev_ux, double * dev_hx)
{
  ConstEBEKernelCopyInIdxKernel<<<512,dim3(16,16),0,ngs_cuda_stream>>> (numblocks, idx, bs, row_dnums, dev_ux, dev_hx);
}

__global__ void ConstEBEKernelCopyOutIdxKernel (int numblocks, int * idx, int bs, int *  col_dnums, double * dev_hy, double * dev_uy)
{
  int tid = threadIdx.x;

  // for (int r = blockIdx.x; r < numblocks; r+=gridDim.x)
  for (int r = blockIdx.x*blockDim.y+threadIdx.y; r < numblocks; r+=gridDim.x*blockDim.y)    
    for (int i = tid; i < bs; i += blockDim.x)
      dev_uy[col_dnums[idx[r]*bs+i]] += dev_hy[r*bs+i];
}

void ConstEBEKernelCopyOutIdx (int numblocks, int * idx, int bs, int * col_dnums, double * dev_hy, double * dev_uy)
{
  ConstEBEKernelCopyOutIdxKernel<<<512,dim3(16,16),0,ngs_cuda_stream>>> (numblocks, idx, bs, col_dnums, dev_hy, dev_uy);
}






/* ************** kernels for DevBlockDiagonalMatrixSoA Matrix ********************** */

__global__ void DevBlockDiagonalMatrixSoAMultAddVecsKernel (double s, int size, double * a, double * b, double * res)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < size; i += blockDim.x*gridDim.x)
    res[i] += s * a[i]*b[i];
  
}

void DevBlockDiagonalMatrixSoAMultAddVecs (double s, int size, double * a, double * b, double * res)
{    
  DevBlockDiagonalMatrixSoAMultAddVecsKernel<<<512,256,0,ngs_cuda_stream>>> (s, size, a, b, res);
}




__global__ void DevBlockDiagonalMatrixSoAMultAddVecsKernel (double s, FlatArray<Dev<int>> inds,
                                                            SliceMatrix<Dev<double>> a,
                                                            SliceMatrix<Dev<double>> b,
                                                            SliceMatrix<Dev<double>> res)
{
  // TODO: copy inds to shared memory ? 
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < res.Width(); i += blockDim.x*gridDim.x)
    for (int j = 0; j < inds.Size(); j+=3)
      {
        int rowa = inds[j];
        int rowb = inds[j+1];
        int rowres = inds[j+2];
        // res[i] += s * a[i]*b[i];
        res(rowres,i) += s * a(rowa,i) * b(rowb,i);
      }
}

// for (i,j,k) in indices:
//    res.Row(k) += s * a.Row(i) * b.Row(j)
void DevBlockDiagonalMatrixSoAMultAddVecs (double s, FlatArray<Dev<int>> inds, 
                                           SliceMatrix<Dev<double>> a, 
                                           SliceMatrix<Dev<double>> b,
                                           SliceMatrix<Dev<double>> res)
{
  DevBlockDiagonalMatrixSoAMultAddVecsKernel<<<512,256,0,ngs_cuda_stream>>> (s, inds, a, b, res);
}




/* ************** kernels for DevProjector Matrix ********************** */

__global__ void DevProjectorMultAddKernel1 (double s, size_t size, const double * a, 
                                            double * b, const unsigned char * bits)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < size; i += blockDim.x * gridDim.x)
  {
    unsigned char mask = (char(1) << (i % CHAR_BIT));
    unsigned int addr = i / CHAR_BIT;

    if (bits[addr] & mask)
      b[i] += s * a[i];
  }
}

__global__ void DevProjectorMultAddKernel2 (double s, size_t size, const double * a, 
                                            double * b, const unsigned char * bits)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < size; i += blockDim.x * gridDim.x)
  {
    unsigned char mask = (char(1) << (i % CHAR_BIT));
    unsigned int addr = i / CHAR_BIT;

    if (! (bits[addr] & mask))
      b[i] += s * a[i];
  }
}

void DevProjectorMultAdd (double s, size_t size, const double * a, double * b, 
                          const unsigned char * bits, bool keep_values)
{
  if (keep_values)
    DevProjectorMultAddKernel1<<<512,256,0,ngs_cuda_stream>>>(s, size, a, b, bits);
  else
    DevProjectorMultAddKernel2<<<512,256,0,ngs_cuda_stream>>>(s, size, a, b, bits);
}

__global__ void DevProjectorProjectKernel1 (size_t size, double * a, const unsigned char * bits)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < size; i += blockDim.x * gridDim.x)
  {
    unsigned char mask = (char(1) << (i % CHAR_BIT));
    unsigned int addr = i / CHAR_BIT;

    if (! (bits[addr] & mask))
      a[i] = 0.0;
  }
}

__global__ void DevProjectorProjectKernel2 (size_t size, double * a, const unsigned char * bits)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < size; i += blockDim.x * gridDim.x)
  {
    unsigned char mask = (char(1) << (i % CHAR_BIT));
    unsigned int addr = i / CHAR_BIT;

    if (bits[addr] & mask)
      a[i] = 0.0;
  }
}

void DevProjectorProject (size_t size, double * a, const unsigned char * bits, 
                          bool keep_values)
{
  if (keep_values)
    DevProjectorProjectKernel1<<<512,256,0,ngs_cuda_stream>>>(size, a, bits);
  else 
    DevProjectorProjectKernel2<<<512,256,0,ngs_cuda_stream>>>(size, a, bits);
}



//TFQMR scalar batch kernels
__global__ void TFQMREvenBatch1Kernel(double* rho, double* vtrstar,
                                       double* theta, double* eta,
                                       double* alpha, double* neg_alpha, double* coeff)
{
    double a   = *rho / *vtrstar;
    *alpha     = a;
    *neg_alpha = -a;
    *coeff     = (*theta) * (*theta) * (*eta) / a;
}

void TFQMREvenBatch1(double* rho, double* vtrstar, double* theta, double* eta,
                     double* alpha, double* neg_alpha, double* coeff)
{
    TFQMREvenBatch1Kernel<<<1,1,0,ngs_cuda_stream>>>(rho, vtrstar, theta, eta,
                                                      alpha, neg_alpha, coeff);
}

__global__ void TFQMREvenTauBatchKernel(double* wnorm_sq, double* tau_in, double* alpha_in,
                                         double* rho, double* theta, double* c,
                                         double* tau_out, double* tau_sq, double* eta, double* rho_last)
{
    double th  = sqrt(abs(*wnorm_sq)) / *tau_in;
    double cs  = 1.0 / sqrt(1.0 + th * th);
    double t   = *tau_in * th * cs;
    *theta     = th;
    *c         = cs;
    *tau_out   = t;
    *tau_sq    = t * t;
    *eta       = cs * cs * (*alpha_in);
    *rho_last  = *rho;
}

void TFQMREvenTauBatch(double* wnorm_sq, double* tau_in, double* alpha_in, double* rho,
                        double* theta, double* c, double* tau_out, double* tau_sq,
                        double* eta, double* rho_last)
{
    TFQMREvenTauBatchKernel<<<1,1,0,ngs_cuda_stream>>>(wnorm_sq, tau_in, alpha_in, rho,
                                                         theta, c, tau_out, tau_sq, eta, rho_last);
}

__global__ void TFQMROddCoeffKernel(double* theta, double* eta, double* alpha, double* coeff)
{
    *coeff = (*theta) * (*theta) * (*eta) / (*alpha);
}

void TFQMROddCoeff(double* theta, double* eta, double* alpha, double* coeff)
{
    TFQMROddCoeffKernel<<<1,1,0,ngs_cuda_stream>>>(theta, eta, alpha, coeff);
}

__global__ void TFQMROddTauBatchKernel(double* wnorm_sq, double* tau_in, double* alpha_in,
                                        double* theta, double* c,
                                        double* tau_out, double* tau_sq, double* eta)
{
    double th  = sqrt(abs(*wnorm_sq)) / *tau_in;
    double cs  = 1.0 / sqrt(1.0 + th * th);
    double t   = *tau_in * th * cs;
    *theta     = th;
    *c         = cs;
    *tau_out   = t;
    *tau_sq    = t * t;
    *eta       = cs * cs * (*alpha_in);
}

void TFQMROddTauBatch(double* wnorm_sq, double* tau_in, double* alpha_in,
                       double* theta, double* c, double* tau_out, double* tau_sq, double* eta)
{
    TFQMROddTauBatchKernel<<<1,1,0,ngs_cuda_stream>>>(wnorm_sq, tau_in, alpha_in,
                                                        theta, c, tau_out, tau_sq, eta);
}

__global__ void TFQMROddBetaKernel(double* rho, double* rho_last, double* beta, double* beta_sq)
{
    double b   = *rho / *rho_last;
    *beta      = b;
    *beta_sq   = b * b;
    *rho_last  = *rho;
}

void TFQMROddBeta(double* rho, double* rho_last, double* beta, double* beta_sq)
{
    TFQMROddBetaKernel<<<1,1,0,ngs_cuda_stream>>>(rho, rho_last, beta, beta_sq);
}


// Sets cudaGraphCondTypeWhile condition: 1 = continue, 0 = stop
// Also increments iter_count and stops when iter_count >= maxsteps
__global__ void ConvergenceCheckKernel(double* rz, double tol,
    cudaGraphConditionalHandle handle, int* iter_count, int maxsteps)
{
    int iter = ::atomicAdd(iter_count, 1) + 1;
    int converged  = (sqrt(abs(*rz)) <= tol);
    int max_reached = (iter >= maxsteps);
    cudaGraphSetConditional(handle,
        (converged || max_reached) ? 0 : 1);
}

void ConvergenceCheck(double* rz, double tol,
    cudaGraphConditionalHandle handle, int* iter_count, int maxsteps)
{
    ConvergenceCheckKernel<<<1, 1, 0, ngs_cuda_stream>>>(
        rz, tol, handle, iter_count, maxsteps);
}

}
