#ifndef CUDA_NGSTD_HPP
#define CUDA_NGSTD_HPP

#include <cuda_runtime.h>
#include <ngstd.hpp>

#include "cuda_profiler.hpp"

namespace ngs_cuda
{
  using namespace ngstd;

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

  
  extern int gpu_clock;
  void InitCUDA (int verbose = 2);


  template <typename T>
  class Dev 
  {
  public:
    T data;
    static Dev<T> * Malloc(size_t size)
    {
      Dev<T> * ptr;
      if (auto err = cudaMalloc (&ptr, size*sizeof(T)))
        throw Exception("cudaMalloc error, ec="+ToString(err));
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
      cudaMemcpy (&data, &val, sizeof(T), cudaMemcpyHostToDevice);
    }


    void D2H (FlatArray<T> hosta)
    {
      cudaMemcpy (hosta.Data(), &data, hosta.Size()*sizeof(T), cudaMemcpyDeviceToHost);
    }

    void H2D (FlatArray<T> hosta)
    {
      cudaMemcpy (&data, hosta.Data(), hosta.Size()*sizeof(T), cudaMemcpyHostToDevice);
    }

    
    __device__ Dev<T> & operator= (T d2) { data = d2; return *this; }
    __device__ operator T() const { return data; } 
    
    template <typename T2>
    __device__ auto & operator+= (T2 other) { data += other; return *this; }
    template <typename T2>
    __device__ auto & operator-= (T2 other) { data -= other; return *this; }
    template <typename T2>
    __device__ auto & operator*= (T2 other) { data *= other; return *this; }
  };
    


  /*
  template <typename T>
  class DevVar
  {
    T * ptr;
  public:
    DevVar()
    {
      cudaMalloc (&ptr, sizeof(T));
    }

    DevVar(T val)
    {
      cudaMalloc (&ptr, sizeof(T));
      cudaMemcpy (ptr, &val, sizeof(T), cudaMemcpyHostToDevice);
    }

    operator T () const
    {
      T tmp;
      cudaMemcpy (&tmp, ptr, sizeof(T), cudaMemcpyDeviceToHost);    
      return tmp;
    }

    T * DevPtr() const { return ptr; }
    T & DevRef() const { return *ptr; }

  };

  template <typename T>
  inline ostream & operator<< (ostream & ost, DevVar<T> & var)
  {
    ost << T(var);
    return ost;
  }
  */

    // TODO: Resize + error checking
  class DevStackMemory
  {
    char * data;
    char * stackptr;
  public:
    DevStackMemory (size_t s = 512*1024*1025)
      {
        cudaMalloc (&data, s);
        stackptr = data;
      }
    
    ~DevStackMemory ()
      {
        cudaFree (data);        
      }
    
    template <typename T>
      T * Alloc (size_t s)
    {
      char * tmp = stackptr;
      s *= sizeof(T);
      s = (s+255) & size_t(-256);
      stackptr += s;
      return reinterpret_cast<T*>(tmp);
    }
    
    void Free (void * ptr)
    {
      stackptr = reinterpret_cast<char*> (ptr);
    }
  };

  extern DevStackMemory stackmemory;

  template <typename T>
  class DevStackArray : public FlatArray<Dev<T>>
  {
  public:
    DevStackArray (size_t s)
      : FlatArray<Dev<T>> (s, (Dev<T>*)stackmemory.Alloc<T>(s))
      { ; } 
    ~DevStackArray ()
      {
        stackmemory.Free(this->data);
      }
    T * DevData () const { return (T*)this->data; }
  };
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


namespace ngs_cuda
{
  // use Array<Dev<T>> instead 
  template <typename T>
  class [[deprecated]] DevArray 
  {
    int size;
    T * dev_data;
  
  public: 
    DevArray (int asize)
    {
      size = asize;
      cudaMalloc((T**)&dev_data, size*sizeof(T));
    }

    DevArray (FlatArray<T> a2)
    {
      size = a2.Size();
      cudaMalloc((T**)&dev_data, size*sizeof(T));
      cudaMemcpy (dev_data, &a2[0], sizeof(T)*size, cudaMemcpyHostToDevice);
    }

    ~DevArray ()
    {
      cudaFree (dev_data);
    }

    T * DevPtr() { return dev_data; }

    DevArray & operator= (FlatArray<T> a2)
    {
      cudaMemcpy (dev_data, &a2[0], sizeof(T)*size, cudaMemcpyHostToDevice);
      return *this;
    }

    void D2H (FlatArray<T> a2) const
    {
      cudaMemcpy (&a2[0], dev_data, sizeof(T)*size, cudaMemcpyDeviceToHost);    
    }

    INLINE int Size() const { return size; }

    /*
    INLINE operator FlatArray<T> ()
    {
      return FlatArray<T> (size, dev_data);
    }
    */
    INLINE FlatArray<T> Dev() const
    {
      return FlatArray<T> (size, dev_data);
    }

    explicit INLINE operator Array<T> () const
    {
      Array<T> temp(size);
#ifdef __CUDA_ARCH__
      temp = FlatArray<T> (*this);
#else
      D2H (temp);
#endif
      return temp;
    }

    INLINE Array<T> Host() const
    {
      return Array<T> (*this);
    }
    
    T * DevData() const { return dev_data; }
  }; 



  template <typename T>
  inline Array<T> D2H (FlatArray<Dev<T>> deva)
  {
    Array<T> hosta(deva.Size());
    cudaMemcpy (hosta.Data(), deva.Data(), sizeof(T)*hosta.Size(), cudaMemcpyDeviceToHost);    
    return hosta;
  }

  template <typename T>
  inline void H2D (FlatArray<Dev<T>> deva, FlatArray<T> hosta)
  {
    cudaMemcpy (deva.Data(), hosta.Data(), sizeof(T)*hosta.Size(), cudaMemcpyHostToDevice);    
  }

  /*
    template <class T>
    class TableWrapper : public Table<T>
    {
    using Table<T>::size;
    using Table<T>::data;
    using Table<T>::index;
    public:
    INLINE TableWrapper (int asize, int * aindex, T * adata)
    // : Table<T> (0,0)
    { 
    size = asize;
    index = aindex;
    data = adata;
    }

    INLINE TableWrapper (const Table<T> & tab)
    // : Table<T> (0,0) 
    {
    const TableWrapper<T> & htab = static_cast<const TableWrapper<T>&> (tab);
    size = htab.size;
    data = htab.data;
    index = htab.index;
    }
    INLINE ~TableWrapper ()
    {
    data = NULL;
    index = NULL;
    }

    INLINE int SizeData() { return index[size]; }
    INLINE int* & Index() { return index; }
    INLINE T* & Data() { return data; }

    // HD const int * & Index() const { return index; }
    // HD const T * & Data() const { return data; }
    };
  */



  // only data at device, but index at host
  template <typename T>
  class DevDataTable
  {
    int size;
    size_t * index = nullptr;
    Dev<T> * dev_data = nullptr;
  
  public: 

    DevDataTable (FlatTable<T> t2)
    {
      size = t2.Size();
      if (size == 0) return;

      index = new size_t[size+1];
      for (int i = 0; i <= size; i++)
        index[i] = t2.IndexArray()[i];

      int sizedata = t2.AsArray().Size();
      dev_data = Dev<T>::Malloc(sizedata);
      cudaMemcpy (dev_data, t2.Data(), sizeof(T)*sizedata, cudaMemcpyHostToDevice);
    }

    ~DevDataTable ()
    {
      Dev<T>::Free (dev_data);
      delete [] index;
    }

    void D2H (FlatTable<T> & t2) const
    {
      int sizedata = t2.AsArray().Size();
      cudaMemcpy (&t2[0][0], dev_data, sizeof(T)*sizedata, cudaMemcpyDeviceToHost);    
    }

    operator FlatTable<Dev<T>> () const
    {
      return FlatTable<Dev<T>> (size, index, dev_data);
    }

    auto Index() const { return index; }
    auto DevData() const { return dev_data; }

    FlatArray<Dev<T>> Row(int i) const { return { index[i+1]-index[i], dev_data+index[i] }; }
    
    class Iterator
    {
      const DevDataTable & tab;
      size_t row;
    public:
      Iterator (const DevDataTable & _tab, size_t _row) : tab(_tab), row(_row) { ; }
      Iterator & operator++ () { ++row; return *this; }
      auto operator* () const { return tab.Row(row); }
      bool operator!= (const Iterator & it2) { return row != it2.row; }
    };

    Iterator begin() const { return Iterator(*this, 0); }
    Iterator end() const { return Iterator(*this, size); }
  }; 


  template <typename T>
  class DevTable
  {
    int size;
    Dev<size_t> * dev_index = nullptr;
    Dev<T> * dev_data = nullptr;
  
  public: 

    DevTable (FlatTable<T> t2)
    {
      size = t2.Size();
      if (size == 0) return;
      
      cudaMalloc((size_t**)&dev_index, (size+1)*sizeof(size_t));
      cudaMemcpy (dev_index, &t2.IndexArray()[0], sizeof(size_t)*(size+1), cudaMemcpyHostToDevice); 
      // cout << "res = " << cudaMemcpy (dev_index, t2.Index(), sizeof(int)*(size+1), cudaMemcpyHostToDevice) << endl;
    
      int sizedata = t2.AsArray().Size();
      cudaMalloc((int**)&dev_data, sizedata*sizeof(T));
      cudaMemcpy (dev_data, t2.Data(), sizeof(T)*sizedata, cudaMemcpyHostToDevice);
    }

    ~DevTable ()
    {
      cudaFree (dev_data);
      cudaFree (dev_index);
    }

    void D2H (FlatTable<T> & t2) const
    {
      int sizedata = t2.AsArray().Size();
      cudaMemcpy (&t2[0][0], dev_data, sizeof(T)*sizedata, cudaMemcpyDeviceToHost);    
    }

    operator FlatTable<T> () const
    {
      return FlatTable<T> (size, (size_t*)dev_index, (T*)dev_data);
    }

    size_t * DevIndex() const { return (size_t*)dev_index; }
    T * DevData() const { return (T*)dev_data; }

    FlatArray<Dev<T>> AsArray() const
    {
      return FlatArray<Dev<T>> ( dev_index[size].D2H(), dev_data );
    }    
  }; 






  class DevBitArray
  {
  protected:
    size_t size = 0;
    unsigned char * dev_data = nullptr;

  public:
    DevBitArray (size_t asize);
    DevBitArray (const ngcore::BitArray & ba);

    ~DevBitArray ();

    DevBitArray & operator= (const ngcore::BitArray &ba);

    size_t Size () const { return size; }
    auto Data () const { return dev_data; }

    void SetSize (size_t asize);

  private:
    size_t Addr (size_t i) const
    {
      return (i / CHAR_BIT);
    }
  };

}


#endif
