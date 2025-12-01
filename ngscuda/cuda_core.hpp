#ifndef CUDA_CORE_HPP
#define CUDA_CORE_HPP

#include <cuda_runtime.h>

#include <core/array.hpp>
#include <core/exception.hpp>

namespace ngs_cuda
{


  // from CUDA C++ Programming Guide:
  // https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
#ifdef __CUDA_ARCH__
#if __CUDA_ARCH__ < 600
  inline __device__ double atomicAdd(double* address, double val)
  {
    unsigned long long int* address_as_ull =
      (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed,
                      __double_as_longlong(val +
                                           __longlong_as_double(assumed)));

      // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
  }
#endif
#endif



  // Kernel wrapper only available if we are compiling the current file with the cuda compiler
#ifdef __CUDACC__
   
  template<class F> __global__
  void CUDA_forall(int n, F f)
  {
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
    for (int i = tid; i < n; i += blockDim.x*gridDim.x)
      f(i);
  }

  template<class F> __global__
  void CUDA_forall2(int n, F f)
  {
    int tid = (blockIdx.x*blockDim.y+threadIdx.y)*blockDim.x+threadIdx.x;
    for (int i = tid; i < n; i += blockDim.x*blockDim.y*gridDim.x)
      f(i);
  }
    
#define DEVICE_LAMBDA __device__

  template <class F>
  inline void DeviceParallelFor (int n, F f)
  {
    // CUDA_forall<<<512,256>>> (n, f);
    CUDA_forall<<<n/256+1,256>>> (n, f);
    // CUDA_forall<<<4096,32>>> (n, f);           // slower
    // CUDA_forall2<<<512,dim3(16,16)>>> (n, f);  // same performance
  }   

#endif // __CUDACC__







  template <typename T>
  class Dev 
  {
    T data;
    
  public:
    __host__ __device__ Dev() = delete;
    
    static Dev<T> * Malloc(size_t size)
    {
      Dev<T> * ptr;
      if (auto err = cudaMalloc (&ptr, size*sizeof(T)))
        throw ngstd::Exception("cudaMalloc error, ec="+ngcore::ToString(err));
      return ptr;        
    }
    
    static void Free(Dev<T> * data)
    {
      cudaFree (data);   
    }


    
    
    T D2H() const
    {
      T res;
      cudaMemcpy (&res, &data, sizeof(T), cudaMemcpyDeviceToHost);
      return res;
    }

    void H2D (T val)

    {
      cudaMemcpy (this, &val, sizeof(T), cudaMemcpyHostToDevice);
    }

    void D2H (ngcore::FlatArray<T> hosta)
    {
      cudaMemcpy (hosta.Data(), &data, hosta.Size()*sizeof(T), cudaMemcpyDeviceToHost);
    }

    void H2D (ngcore::FlatArray<T> hosta)
    {
      cudaMemcpy (&data, hosta.Data(), hosta.Size()*sizeof(T), cudaMemcpyHostToDevice);
    }
    
    __host__ __device__ Dev<T> & operator= (T d2)
    {
#ifdef __CUDA_ARCH__      
      data = d2;
#else
      H2D(d2);
#endif
      return *this;
    }

    __host__ __device__ operator T() const
    {
#ifdef __CUDA_ARCH__
      return data;
#else
      return D2H();
#endif
    }
    
    __device__ const T& operator*() const { return data; }
    
    template <typename T2>
    __device__ auto & operator+= (T2 other) { data += other; return *this; }
    template <typename T2>
    __device__ auto & operator-= (T2 other) { data -= other; return *this; }
    template <typename T2>
    __device__ auto & operator*= (T2 other) { data *= other; return *this; }
  };


  template <typename T>
  Dev<T> * Host2Device (const T& val)
  {
    auto ptr = Dev<T>::Malloc(1);
    ptr->H2D(val);
    return ptr;
  }

  template <typename T>
  void Free(Dev<T> * dp)
  {
    Dev<T>::Free (dp);
  }  

  

}

#endif

