// x = val
__global__ void SetScalarKernel (double val, int n, double * x)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < n; i += blockDim.x*gridDim.x)
    x[i] = val;
}

void SetScalar (double val, int n, double * x)
{
  SetScalarKernel<<<512,256>>> (val, n, x);
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




/* ************** kernels for ConstantEBE Matrix ********************** */


__global__ void ConstEBEKernelCopyInKernel (int numblocks, int bs, int * row_dnums, double * dev_ux, double * dev_hx)
{
  int tid = threadIdx.x;
  for (int r = blockIdx.x*blockDim.y+blockIdx.y; r < numblocks; r+=gridDim.x*blockDim.y)
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

  // for (int r = blockIdx.x; r < numblocks; r+=gridDim.x)
  for (int r = blockIdx.x*blockDim.y+blockIdx.y; r < numblocks; r+=gridDim.x*blockDim.y)  
    for (int i = tid; i < bs; i += blockDim.x)
      dev_uy[col_dnums[r*bs+i]] += dev_hy[r*bs+i];
}

void ConstEBEKernelCopyOut (int numblocks, int bs, int * col_dnums, double * dev_hy, double * dev_uy)
{
  ConstEBEKernelCopyOutKernel<<<512,dim3(16,16)>>> (numblocks, bs, col_dnums, dev_hy, dev_uy);
}




__global__ void ConstEBEKernelCopyInIdxKernel (int numblocks, int * idx, int bs, int * row_dnums, double * dev_ux, double * dev_hx)
{
  int tid = threadIdx.x;
  // for (int r = blockIdx.x; r < numblocks; r+=gridDim.x)
  for (int r = blockIdx.x*blockDim.y+blockIdx.y; r < numblocks; r+=gridDim.x*blockDim.y)    
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
  for (int r = blockIdx.x*blockDim.y+blockIdx.y; r < numblocks; r+=gridDim.x*blockDim.y)    
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

