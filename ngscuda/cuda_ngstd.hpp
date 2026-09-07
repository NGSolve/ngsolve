#ifndef CUDA_NGSTD_HPP
#define CUDA_NGSTD_HPP

#include <cuda_runtime.h>
#include <ngstd.hpp>

#include "cuda_core.hpp"
#include "cuda_profiler.hpp"

namespace ngs_cuda
{
  using namespace ngstd;

  
  extern int gpu_clock;
  void InitCUDA (int verbose = 2);
  void WarmupCudaModule ();   // first launch of a kernel of this library

}

namespace ngcore {
  template <typename T>
  struct IsSafe<ngs_cuda::Dev<T>> {
    constexpr operator bool() const { return true; }
  };
}

namespace ngstd
{
  template <typename T>
  struct my_is_integral<ngs_cuda::Dev<T>> : my_is_integral<T>{};
}

namespace ngcore 
{
  using ngs_cuda::Dev;
  template <typename T>  
  class Array<Dev<T>> : public FlatArray<Dev<T>>
  {
  public:
    Array() = default;
    Array (size_t s)
      : FlatArray<Dev<T>>(s, Dev<T>::Malloc(s)) { } ;     
    Array (FlatArray<T> a2)
      : Array(a2.Size())
    {
      this->data->H2D(a2);
    }
    
    Array& operator= (Array<Dev<T>> && a2)
    {
      Swap (this->data, a2.data);
      Swap (this->size, a2.size);
      return *this;
    }
    
    Array& operator= (FlatArray<T> a2)
    {
      SetSize(a2.Size());
      this->data->H2D(a2);
      return *this;
    }
    
    void SetSize(size_t s)
    {
      if (this->Size() != s)
        {
          Dev<T>::Free(this->data);
          this->data = Dev<T>::Malloc(s);
          this->size = s;
        }
    }
    
    ~Array()
    {
      Dev<T>::Free(this->data);
    }
  };
}

#endif
