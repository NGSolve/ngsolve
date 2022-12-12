// x = val
__global__ void SetScalarKernel (double val, int n, double * x)
{
  int tid = threadIdx.x;
  for (int i = tid; i < n; i += blockDim.x)
    x[i] = val;
}

void SetScalar (double val, int n, double * x)
{
  SetScalarKernel<<<1,128>>> (val, n, x);
} 


// y = D * x
__global__ void MultDiagonalKernel (int n, double * D, double * x, double * y)
{
  int tid = threadIdx.x;
  for (int i = tid; i < n; i += blockDim.x)
    y[i] = D[i] * x[i];
}

void MultDiagonal (int n, double * D, double * x, double * y)
{
  MultDiagonalKernel<<<1,128>>> (n, D, x, y);
} 


// y += alpha D * x
__global__ void MultAddDiagonalKernel (int n, double alpha, double * D, double * x, double * y)
{
  int tid = threadIdx.x;
  for (int i = tid; i < n; i += blockDim.x)
    y[i] += alpha * D[i] * x[i];
}

void MultAddDiagonal (int n, double alpha, double * D, double * x, double * y)
{
  MultAddDiagonalKernel<<<1,128>>> (n, alpha, D, x, y);
} 




/* ************** kernels for ConstantEBE Matrix ********************** */

void ConstEBEKernelCopyInKernel (int numblocks, FlatTable<int> row_dnums, double * dev_ux, double * dev_hx)
{
  int tid = threadIdx.x;
  int wm = row_dnums[0].Size();
  
  for (int r = 0; r < numblocks; r++)
    {
      for (int i = tid; i < wm; i += blockDim.x)
        dev_hx[r*wm+i] = dev_ux[row_dnums[r][i]];
    }
}

void ConstEBEKernelCopyIn (int numblocks, const DevTable<int> & row_dnums, double * dev_ux, double * dev_hx)
{
  ConstEBEKernelCopyInKernel<<<1,128>>> (numblocks, row_dnums, dev_ux, dev_hx)
}

void ConstEBEKernelCopyOutKernel (int numblocks, FlatTable<int> col_dnums, double * dev_hy, double * dev_uy)
{
  int tid = threadIdx.x;
  int hm = col_dnums[0].Size();
  
  for (int r = 0; r < numblocks; r++)
    {
      for (int i = tid; i < hm; i += blockDim.x)
        dev_uy[col_dnums[r][i]] += dev_hx[r*hm+i];
    }
}

void ConstEBEKernelCopyOut (int numblocks, const DevTable<int> & col_dnums, double * dev_hy, double * dev_uy)
{
  ConstEBEKernelCopyOutKernel<<<1,128>>> (numblocks, col_dnums, dev_hy, dev_uy)
}
