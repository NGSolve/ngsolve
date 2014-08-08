#ifdef CUDA

#include <cuda_runtime.h>

namespace ngs_cuda
{
  using namespace ngstd;

  
  void InitCUDA (int verbose = 2);


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

    T * Data() { return dev_data; }

    DevArray & operator= (FlatArray<T> a2)
    {
      cudaMemcpy (dev_data, &a2[0], sizeof(T)*size, cudaMemcpyHostToDevice);
      return *this;
    }

    void D2H (FlatArray<T> & a2)
    {
      cudaMemcpy (&a2[0], dev_data, sizeof(T)*size, cudaMemcpyDeviceToHost);    
    }

    operator FlatArray<T> ()
    {
      return FlatArray<T> (size, dev_data);
    }
  }; 







  template <class T>
  class TableWrapper : public Table<T>
  {
    using Table<T>::size;
    using Table<T>::data;
    using Table<T>::index;
  public:
    HD TableWrapper (int asize, int * aindex, T * adata)
    // : Table<T> (0,0)
    { 
      size = asize;
      index = aindex;
      data = adata;
    }

    HD TableWrapper (const Table<T> & tab)
    // : Table<T> (0,0) 
    {
      const TableWrapper<T> & htab = static_cast<const TableWrapper<T>&> (tab);
      size = htab.size;
      data = htab.data;
      index = htab.index;
    }
    HD ~TableWrapper ()
    {
      data = NULL;
      index = NULL;
    }

    HD int SizeData() { return index[size]; }
    HD int* & Index() { return index; }
    HD T* & Data() { return data; }
    // HD const int * & Index() const { return index; }
    // HD const T * & Data() const { return data; }
  };


  template <typename T>
  class DevTable
  {
    int size;
    int * dev_index;
    T * dev_data;
  
  public: 

    DevTable (TableWrapper<T> t2)
    {
      size = t2.Size();
      cudaMalloc((int**)&dev_index, (size+1)*sizeof(int));
      cudaMemcpy (dev_index, t2.Index(), sizeof(int)*(size+1), cudaMemcpyHostToDevice); 
      // cout << "res = " << cudaMemcpy (dev_index, t2.Index(), sizeof(int)*(size+1), cudaMemcpyHostToDevice) << endl;
    
      int sizedata = t2.SizeData();
      cudaMalloc((int**)&dev_data, sizedata*sizeof(T));
      cudaMemcpy (dev_data, t2.Data(), sizeof(T)*sizedata, cudaMemcpyHostToDevice);
    }

    ~DevTable ()
    {
      cudaFree (dev_data);
      cudaFree (dev_index);
    }

    operator TableWrapper<T> () const
    {
      return TableWrapper<T> (size, dev_index, dev_data);
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
