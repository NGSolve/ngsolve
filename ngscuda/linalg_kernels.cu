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
