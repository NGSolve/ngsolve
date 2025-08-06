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
  SetScalarKernel<<<512,256>>> (val, n, x);
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
  SetVectorKernel<<<512,256>>> (val, n, x, y);
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
  MyDaxpyKernel<<<512,256>>> (val, n, x, y);
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
  MultDiagonalKernel<<<512,256>>> (n, D, x, y);
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
  MultAddDiagonalKernel<<<512,256>>> (n, alpha, D, x, y);
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
  ManyMatVecKernel<<<512,dim3(16,16)>>> (matvecs, x, y);
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
   ConstEBEKernelCopyInKernel<<<512,dim3(16,16)>>> (numblocks, bs, row_dnums, dev_ux, dev_hx);
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
  ConstEBEKernelCopyOutKernel<<<512,dim3(16,16)>>> (numblocks, bs, col_dnums, dev_hy, dev_uy);
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
  ConstEBEKernelCopyInIdxKernel<<<512,dim3(16,16)>>> (numblocks, idx, bs, row_dnums, dev_ux, dev_hx);
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
  ConstEBEKernelCopyOutIdxKernel<<<512,dim3(16,16)>>> (numblocks, idx, bs, col_dnums, dev_hy, dev_uy);
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
  DevBlockDiagonalMatrixSoAMultAddVecsKernel<<<512,256>>> (s, size, a, b, res);
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
  DevBlockDiagonalMatrixSoAMultAddVecsKernel<<<512,256>>> (s, inds, a, b, res);
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
    DevProjectorMultAddKernel1<<<512,256>>>(s, size, a, b, bits);
  else
    DevProjectorMultAddKernel2<<<512,256>>>(s, size, a, b, bits);
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
    DevProjectorProjectKernel1<<<512,256>>>(size, a, bits);
  else 
    DevProjectorProjectKernel2<<<512,256>>>(size, a, bits);
}


}

