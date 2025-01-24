/*********************************************************************/
/* File:   demo1.cpp                                                 */
/* Author: Joachim Schoeberl                                         */
/* Date:   11. Aug. 2014                                             */
/*********************************************************************/

#include <core/python_ngcore.hpp>
#include <la.hpp>
#include "cuda_linalg.hpp"

using namespace std;



// c[i] = a[i]+b[i]
__global__ void MyAddKernel1 (int n, double * a, double * b, double * c)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < n; i += blockDim.x*gridDim.x)
    c[i] += a[i]+b[i];
}

void Test1 (Vector<> & a, Vector<> & b, Vector<> & c)
{
  // c = a+b;
    
  double *dev_a, *dev_b, *dev_c;
  size_t n = a.Size();
      
  cudaMalloc (&dev_a, n*sizeof(double));
  cudaMalloc (&dev_b, n*sizeof(double));
  cudaMalloc (&dev_c, n*sizeof(double));
 
  cudaMemcpy (dev_a, a.Data(), n*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy (dev_b, b.Data(), n*sizeof(double), cudaMemcpyHostToDevice);
    
  MyAddKernel1<<<128,512>>> (n, dev_a, dev_b, dev_c);
    
  cudaMemcpy (c.Data(), dev_c, n*sizeof(double), cudaMemcpyDeviceToHost);
    
  cudaFree (dev_a);
  cudaFree (dev_b);
  cudaFree (dev_c);
}



__global__ void MyAddKernel2 (FlatArray<Dev<double>> a, FlatArray<Dev<double>> b, FlatArray<Dev<double>> c)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  for (int i = tid; i < a.Size(); i += blockDim.x*gridDim.x)
    c[i] += a[i]+b[i]; 
}

void Test2 (Vector<> & a, Vector<> & b, Vector<> & c)
{
  FlatArray<double> ha(a.Size(), a.Data());
  FlatArray<double> hb(b.Size(), b.Data());
  FlatArray<double> hc(c.Size(), c.Data());
    
  Array<Dev<double>> deva (ha);
  Array<Dev<double>> devb (hb);
  Array<Dev<double>> devc (hc.Size());
    
  MyAddKernel2<<<128,512>>> (deva, devb, devc);
    
  hc = D2H(devc);
}




void Test3 (Vector<> & a, Vector<> & b, Vector<> & c)
{
  FlatArray<double> ha(a.Size(), a.Data());
  FlatArray<double> hb(b.Size(), b.Data());
  FlatArray<double> hc(c.Size(), c.Data());
    
  Array<Dev<double>> deva (ha);
  Array<Dev<double>> devb (hb);
  Array<Dev<double>> devc (hc.Size());

  DeviceParallelFor
    (a.Size(), [a=FlatArray(deva), b=FlatArray(devb), c=FlatArray(devc)] DEVICE_LAMBDA  (auto tid)
       {
         c[tid] = a[tid]+b[tid];
       });
    
  hc = D2H(devc);
}




void ExportDemo(py::module & m)
{
   m.def("Test1", &Test1);
   m.def("Test2", &Test2);
   m.def("Test3", &Test3);
}



