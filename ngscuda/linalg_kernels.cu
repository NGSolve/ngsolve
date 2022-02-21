__global__ void SetScalarKernel (double val, int n, double * dev_ptr)
{
  int tid = threadIdx.x;
  for (int i = tid; i < n; i += blockDim.x)
    dev_ptr[i] = val;
}


void SetScalar (double val, int n, double * dev_ptr)
{
  SetScalarKernel<<<1,128>>> (val, n, dev_ptr);
} 


