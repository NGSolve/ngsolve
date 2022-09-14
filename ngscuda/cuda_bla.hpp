#ifdef CUDA

namespace ngs_cuda
{
  using namespace ngbla;


  template <typename T = double>
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

    void D2H (FlatVector<T> a2) const
    {
      cudaMemcpy (&a2[0], dev_data, sizeof(T)*size, cudaMemcpyDeviceToHost);    
    }

    INLINE int Size() const { return size; }

    /*
    INLINE operator FlatVector<T> ()
    {
      return FlatVector<T> (size, dev_data);
    }
    */

    INLINE FlatVector<T> Dev() const
    {
      return FlatVector<T> (size, dev_data);
    }

    explicit INLINE operator Vector<T> () const
    {
      Vector<T> temp(size);
#ifdef __CUDA_ARCH__
      temp = FlatVector<T> (*this);
#else
      D2H (temp);
#endif
      return temp;
    }

    INLINE Vector<T> Host() const
    {
      return Vector<T> (*this);
    }

  }; 
















  template <typename T = double>
  class DevMatrix 
  {
    int h, w;
    T * dev_data;
  
  public: 
    DevMatrix (int ah, int aw)
    {
      h = ah;
      w = aw;
      cudaMalloc((T**)&dev_data, h*w*sizeof(T));
    }

    DevMatrix (FlatMatrix<T> a2)
    {
      h = a2.Height();
      w = a2.Width();
      cudaMalloc((T**)&dev_data, h*w*sizeof(T));
      cudaMemcpy (dev_data, &a2(0,0), sizeof(T)*h*w, cudaMemcpyHostToDevice);
    }

    ~DevMatrix ()
    {
      cudaFree (dev_data);
    }

    T * Data() { return dev_data; }

    DevMatrix & operator= (FlatMatrix<T> a2)
    {
      cudaMemcpy (dev_data, &a2(0,0), sizeof(T)*h*w, cudaMemcpyHostToDevice);
      return *this;
    }

    void D2H (FlatMatrix<T> a2) const
    {
      cudaMemcpy (&a2(0,0), dev_data, sizeof(T)*h*w, cudaMemcpyDeviceToHost);    
    }

    INLINE int Height() const { return h; }
    INLINE int Width() const { return w; }

    /*
    INLINE operator FlatVector<T> ()
    {
      return FlatVector<T> (size, dev_data);
    }
    */

    INLINE FlatMatrix<T> Dev() const
    {
      return FlatMatrix<T> (h, w, dev_data);
    }

    explicit INLINE operator Matrix<T> () const
    {
      Matrix<T> temp(h,w);
#ifdef __CUDA_ARCH__
      temp = FlatMatrix<T> (*this);
#else
      D2H (temp);
#endif
      return std::move(temp);
    }

    INLINE Matrix<T> Host() const
    {
      return Matrix<T> (*this);
    }

  }; 




}

#endif
