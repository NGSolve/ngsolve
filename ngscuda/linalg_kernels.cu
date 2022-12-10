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



