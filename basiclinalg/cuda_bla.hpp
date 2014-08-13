#ifdef CUDA

namespace ngs_cuda
{
  using namespace ngbla;


  template <typename T>
  class DevVector 
  {
    int size;
    T * dev_data;
  
  public: 
    DevVector (int asize)
    {
      size = asize;
      cudaMalloc((T**)&dev_data, size*sizeof(T));
    }

    DevVector (FlatVector<T> a2)
    {
      size = a2.Size();
      cudaMalloc((T**)&dev_data, size*sizeof(T));
      cudaMemcpy (dev_data, &a2[0], sizeof(T)*size, cudaMemcpyHostToDevice);
    }

    ~DevVector ()
    {
      cudaFree (dev_data);
    }

    T * Data() { return dev_data; }

    DevVector & operator= (FlatVector<T> a2)
    {
      cudaMemcpy (dev_data, &a2[0], sizeof(T)*size, cudaMemcpyHostToDevice);
      return *this;
    }

    void D2H (FlatVector<T> & a2)
    {
      cudaMemcpy (&a2[0], dev_data, sizeof(T)*size, cudaMemcpyDeviceToHost);    
    }

    INLINE int Size() const { return size; }

    INLINE operator FlatVector<T> ()
    {
      return FlatVector<T> (size, dev_data);
    }

    explicit INLINE operator Vector<T> ()
    {
      Vector<T> temp(size);
#ifdef __CUDA_ARCH__
      temp = FlatVector<T> (*this);
#else
      D2H (temp);
#endif
      return temp;
    }
  }; 

}

#endif
