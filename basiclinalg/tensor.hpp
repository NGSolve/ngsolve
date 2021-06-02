#include <bla.hpp>





/*
  // example use of tensors:

int main ()
{
  Tensor<3> tensor(3,3,3);
  Tensor<3> tensor2(3,3,3);

  Matrix<> mat(3,3);
  mat = 0;
  for (int i = 0; i < mat.Height(); i++)
    mat(i,i) = i+1;

  cout << "mat = " << endl << mat << endl;

  tensor = 1.0;

  for (int j = 0; j < 3; j++)
    tensor2(j,STAR,STAR) = j * mat * tensor(STAR,STAR,j);

  cout << "tensor2 = " << endl << tensor2 << endl;

  cout << "sub = " << endl << tensor2(STAR,STAR,0) << endl;

  return 0;
}

*/


namespace ngbla
{


enum typestar { STAR };

template <int DIM, typename T = double, int DIMLIN = DIM> class FlatTensor;




template <int DIM, typename T, int LINDIM>
INLINE auto LargerTensor (FlatTensor<DIM,T,LINDIM> tensor, size_t as, size_t ad)
  -> FlatTensor<DIM+1,T,LINDIM> 
{
  return FlatTensor<DIM+1,T,LINDIM> (as, ad, tensor);
}

/*
template <typename T, int LINDIM>
INLINE auto LargerTensor (FlatTensor<0,T,LINDIM> tensor, int as, int ad)
  -> SliceVector<T> 
{
  return SliceVector<T> (as, ad, tensor.Data());
}
*/
template <typename T>
INLINE auto LargerTensor (FlatTensor<0,T,0> tensor, size_t as, size_t ad)
  -> SliceVector<T> 
{
  return SliceVector<T> (as, ad, tensor.Data());
}
template <typename T, int LINDIM>
INLINE auto LargerTensor (FlatTensor<0,T,LINDIM> tensor, size_t as, size_t ad)
  -> FlatVector<T> 
{
  return FlatVector<T> (as, tensor.Data());
}


template <typename T>
INLINE auto LargerTensor (SliceVector<T> vec, size_t as, size_t ad)
  -> DoubleSliceMatrix<T> 
{
  return DoubleSliceMatrix<T> (as, vec.Size(),
                               ad, vec.Dist(), vec.Data());
}

template <typename T>
INLINE auto LargerTensor (FlatVector<T> vec, int as, int ad)
  -> SliceMatrix<T> 
{
  return SliceMatrix<T> (as, vec.Size(),
                         ad, vec.Data());
}




template <typename T>
INLINE auto LargerTensor (DoubleSliceMatrix<T> mat, int as, int ad)
  -> FlatTensor<3,T,0> 
{
  FlatTensor<2,T,0> tens (mat.Height(), mat.DistRow(),
                          FlatTensor<1,T> (mat.Width(),mat.DistCol(),
                                           FlatTensor<0,T> (mat.Data())));
  
  return LargerTensor (tens , as, ad);
}

template <typename T>
INLINE auto LargerTensor (SliceMatrix<T> mat, int as, int ad)
  -> FlatTensor<3,T,0> 
{
  FlatTensor<2,T,1> tens (mat.Height(), mat.Dist(),
                          FlatTensor<1,T,1> (mat.Width(),1,
                                             FlatTensor<0,T,1> (mat.Data())));
  
  return LargerTensor (tens , as, ad);
}







template <int DIM, typename T, int LINDIM>
INLINE auto OffsetTensor (FlatTensor<DIM,T,LINDIM> tensor, int offset)
  -> FlatTensor<DIM,T,((LINDIM<DIM)?LINDIM:DIM)>
{
  tensor.Data() += offset;
  return tensor;
}

template <typename T>
INLINE auto OffsetTensor (SliceVector<T> vec, int offset)
  -> SliceVector<T> 
{
  return SliceVector<T> (vec.Size(), vec.Dist(), vec.Data()+offset);
}

template <typename T>
INLINE auto OffsetTensor (FlatVector<T> vec, int offset)
  -> FlatVector<T> 
{
  return FlatVector<T> (vec.Size(), vec.Data()+offset);
}

template <typename T>
INLINE auto OffsetTensor (DoubleSliceMatrix<T> mat, int offset)
  -> DoubleSliceMatrix<T> 
{
  return DoubleSliceMatrix<T> (mat.Height(),mat.Width(),
                               mat.DistRow(),mat.DistCol(),
                               mat.Data()+offset);
}

template <typename T>
INLINE auto OffsetTensor (SliceMatrix<T> mat, int offset)
  -> SliceMatrix<T> 
{
  return SliceMatrix<T> (mat.Height(),mat.Width(),
                         mat.Dist(),
                         &mat(0,0)+offset);
}


// a 0-tensor is treated as number:
template <typename T> auto ReduceTensor0 (T tensor) { return tensor; }
template <typename T, int DIMLIN> T & ReduceTensor0 (FlatTensor<0,T,DIMLIN> tensor) { return *tensor.Data(); } 


template <int DIM, typename T, int DIMLIN>
class FlatTensor
{
  size_t size;
  size_t dist;
  FlatTensor<DIM-1,T, DIMLIN> sub;


public: 
  FlatTensor () : size(0), dist(0) { ; }
  
  template <typename ... ARG>
  FlatTensor (size_t s, ARG ... args) : size(s), sub(args...) 
  {
    dist = sub.GetSize()*sub.GetDist(); 
  }

  template <typename ... ARG>
  FlatTensor (LocalHeap & lh, ARG ... args)
    : FlatTensor (args...)
  {
    size_t totsize = this->GetTotalSize();
    // TODO: why this call instead of lh.Alloc<T> as in FlatMatrix?
    this->Data() = new(lh) T[totsize];
  }
  
  template <typename ... ARG>
  FlatTensor (T * data, ARG ... args)
    : FlatTensor (args...)
  {
    this->Data() = data;
  }

  
  FlatTensor (size_t as, size_t ad, FlatTensor<DIM-1,T> asub) 
    : size(as), dist(ad), sub(asub) { ; }

  template <int DIMLIN2> 
  FlatTensor (FlatTensor<DIM,T,DIMLIN2> t2)
    : size(t2.GetSize()), dist(t2.GetDist()), sub(t2.GetSubTensor())
  {
    static_assert(DIMLIN <= DIMLIN2, "illegal tensor copy");
    ;
  }

  
  template <typename ... ARG>
  auto braces (size_t i, ARG ... args)
  {
    return OffsetTensor (sub.braces(args...), i*dist);
  }

  template <typename ... ARG>
  auto braces (typestar star, ARG ... args) 
  {
    return LargerTensor(sub.braces(args...), size, dist);
  }


  template <typename ... ARG>
  auto operator() (size_t i, ARG ... args)
    -> decltype (ReduceTensor0 (this->braces(i,args...)))
  {
    return ReduceTensor0 (braces(i,args...));
  }
 
  template <typename ... ARG>
  auto operator() (typestar star, ARG ... args) 
  {
    return braces(star,args...);
  }

  
  FlatTensor operator= (double d)
  {
    for (size_t i = 0; i < size; i++)
      GetSubTensor(i) = d;
    return *this;
  }

  size_t GetSize () const { return size; }
  size_t GetDist () const { return dist; }
  size_t GetTotalSize () const { return size*sub.GetTotalSize(); }

  auto GetSubTensor() const -> decltype(sub)
  { return sub; } 
  FlatTensor<DIM-1,T> GetSubTensor (size_t i) const 
  {
    return OffsetTensor (sub, i*dist);
  }

  T *& Data () { return sub.Data(); }
  T * Data () const { return sub.Data(); }
  
  
  
  template<typename ... ARG>
  INLINE void SetSize (size_t s, ARG ... args) throw ()
  {
    size = s;
    sub.SetSize(args...);
  }

  /// TODO: Still problems!
  /// copy data and sub pointers
  INLINE FlatTensor & Assign (const FlatTensor & m) throw()
  {
    sub.Assign(m.sub);
    dist = m.dist;
    size = m.size;
    return *this;
  }

  //TODO: is this feasible? Still problems!
  /// set size, and assign mem
  template<typename ... ARG>
  INLINE void AssignMemory (LocalHeap & lh, size_t s, ARG ... args) throw ()
  {
    FlatTensor tmp{lh, s, args...};
    Assign(tmp);
//    FlatTensor tmp{args...};
//    Assign(tmp);
//    size_t totsize = tmp->GetTotalSize();
//    this->Data() = lh.Alloc<T>(totsize);
  }
  
  /// set size, and assign mem
  template<typename ... ARG>
  INLINE void AssignMemory (T * mem, size_t s, ARG ... args) throw()
  {
    FlatTensor tmp{mem, s, args...};
    Assign(tmp);
//    FlatTensor tmp{args...};
//    this->Data() = mem;
  }
};



template <typename T, int LINDIM>
class FlatTensor<0,T,LINDIM>
{
  T * data;
public: 
  FlatTensor () { ; }
  FlatTensor (T * adata) : data(adata) { ; }
  template <int DIMLIN2>
  FlatTensor (FlatTensor<0,T,DIMLIN2> t2)
    : data(t2.Data()) { ; } 
              
  size_t GetSize () const { return 1; }
  size_t GetDist () const { return 1; }
  size_t GetTotalSize () const { return 1; }
  T *& Data () { return data; }
  T * Data () const { return data; }  
  FlatTensor<0,T,LINDIM> operator() () { return FlatTensor<0,T,LINDIM> (data); }
  FlatTensor<0,T,LINDIM> braces () { return FlatTensor<0,T,LINDIM> (data); }
  // const T & operator() () const { return *data; }
  // T & operator() () { return *data; }
  operator T  () const { return *data; }
  operator T& () { return *data; }
  T & operator= (double d) { *data = d; return *data; }
  T & operator-= (double d) { *data -= d; return *data; }
  T & operator+= (double d) { *data += d; return *data; }

  template<typename ... ARG>
  INLINE void SetSize (ARG ... args) throw ()
  {
  }

  INLINE void Assign (const FlatTensor& m) throw ()
  {
    this->data = m.data;
  }

//  //TODO: is this required?
//  /// set size, and assign mem
//  template<typename ... ARG>
//  INLINE void AssignMemory (LocalHeap & lh, size_t s, ARG ... args) throw ()
//  {
//    this->Data() = lh.Alloc<T>(this->GetSize());
//  }
//
//  /// set size, and assign mem
//  template<typename ... ARG>
//  INLINE void AssignMemory (T * mem, size_t s, ARG ... args) throw()
//  {
//    this->Data() = mem;
//  }
};













template <int DIM, typename T = double>
class Tensor : public FlatTensor<DIM,T>
{
public:
  
  template <typename ... ARG>
  Tensor (ARG ... args) : FlatTensor<DIM,T> (args...)
  {
    size_t totsize = this->GetTotalSize();
    this->Data() = new T[totsize];
  }
  
  ~Tensor() 
  {
    delete this->Data();
  }

  FlatTensor<DIM,T> operator= (double d)
  {
    return FlatTensor<DIM,T>::operator=(d);
  }
};







template <int DIM, typename T>
INLINE ostream & operator<< (ostream & ost, const FlatTensor<DIM,T> & tensor)
{
  ost << "tensor, dim = " << tensor.GetSize() 
      << ", dist = " << tensor.GetDist()
      << ", size = " << tensor.GetTotalSize() << endl;

  for (size_t i = 0; i < tensor.GetSize(); i++)
    {
      ost << "subtensor " << i;
      for (size_t j = 0; j < DIM-1; j++) ost << ",*";
      ost << ":" << endl;
      ost << tensor.GetSubTensor(i);
    }

  return ost;
}

template <typename T>
INLINE ostream & operator<< (ostream & ost, const FlatTensor<2,T> & tensor)
{
  for (size_t i = 0; i < tensor.GetSize(); i++)
    ost << tensor.GetSubTensor(i);
  return ost;
}

template <typename T>
INLINE ostream & operator<< (ostream & ost, const FlatTensor<1,T> & tensor)
{
  for (size_t i = 0; i < tensor.GetSize(); i++)
    ost << * (tensor.Data()+i*tensor.GetDist()) << " ";
  return ost << endl;
}


template <typename T>
INLINE ostream & operator<< (ostream & ost, const FlatTensor<0,T> & tensor)
{
  return ost << *tensor.Data();
}




}
