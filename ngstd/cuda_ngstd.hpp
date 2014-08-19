#ifdef CUDA

#include <cuda_runtime.h>

namespace ngs_cuda
{
  using namespace ngstd;

  
  void InitCUDA (int verbose = 2);







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









  template <typename T>
  class DevArray 
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

  }; 






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


  template <typename T>
  class DevTable
  {
    int size;
    int * dev_index;
    T * dev_data;
  
  public: 

    DevTable (const Table<T> & t2)
    {
      size = t2.Size();
      cudaMalloc((int**)&dev_index, (size+1)*sizeof(int));
      cudaMemcpy (dev_index, &t2.IndexArray()[0], sizeof(int)*(size+1), cudaMemcpyHostToDevice); 
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
      return FlatTable<T> (size, dev_index, dev_data);
    }
  }; 










  /*
    template <typename T = double>
    class DevMatrix
    {
    int h;
    int w; 
    T * dev_data;
  
    public: 
    DevMatrix (int ah, int aw)
    {
    h = ah; w = aw;
    cudaMalloc((T**)&dev_data, h*w*sizeof(T));
    }

    DevMatrix (FlatMatrix<T> m2)
    {
    h = m2.Height();
    w = m2.Width();
    cudaMalloc((T**)&dev_data, h*w*sizeof(T));
    cudaMemcpy (dev_data, &m2(0,0), sizeof(T)*h*w, cudaMemcpyHostToDevice);
    }

    ~DevMatrix ()
    {
    cudaFree (dev_data);
    }

    T * Data() { return dev_data; }

    DevMatrix & operator= (FlatMatrix<T> m2)
    {
    cudaMemcpy (dev_data, &m2(0,0), sizeof(T)*h*w, cudaMemcpyHostToDevice);
    return *this;
    }

    void D2H (FlatMatrix<T> & m2)
    {
    cudaMemcpy (&m2(0,0), dev_data, sizeof(T)*h*w, cudaMemcpyDeviceToHost);    
    }

    operator FlatMatrix<T> ()
    {
    return FlatMatrix<T> (h, w, dev_data);
    }
    }; 
  */
}


#endif
