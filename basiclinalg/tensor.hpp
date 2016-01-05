#include <bla.hpp>
using namespace ngbla;





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





enum typestar { STAR };

template <int DIM, typename T = double, int DIMLIN = DIM> class FlatTensor;




template <int DIM, typename T, int LINDIM>
INLINE auto LargerTensor (FlatTensor<DIM,T,LINDIM> tensor, int as, int ad)
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
INLINE auto LargerTensor (FlatTensor<0,T,0> tensor, int as, int ad)
  -> SliceVector<T> 
{
  return SliceVector<T> (as, ad, tensor.Data());
}
template <typename T, int LINDIM>
INLINE auto LargerTensor (FlatTensor<0,T,LINDIM> tensor, int as, int ad)
  -> FlatVector<T> 
{
  return FlatVector<T> (as, tensor.Data());
}


template <typename T>
INLINE auto LargerTensor (SliceVector<T> vec, int as, int ad)
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






template <int DIM, typename T, int DIMLIN>
class FlatTensor
{
  int size;
  int dist;
  FlatTensor<DIM-1,T, DIMLIN> sub;


public: 
  FlatTensor () : size(0), dist(0) { ; }
  
  template <typename ... ARG>
  FlatTensor (int s, ARG ... args) : size(s), sub(args...) 
  {
    dist = sub.GetSize()*sub.GetDist(); 
  }

  template <typename ... ARG>
  FlatTensor (LocalHeap & lh, ARG ... args)
    : FlatTensor (args...)
  {
    int totsize = this->GetTotalSize();    
    this->Data() = new(lh) T[totsize];
  }  
  
  FlatTensor (int as, int ad, FlatTensor<DIM-1,T> asub) 
    : size(as), dist(ad), sub(asub) { ; }

  template <int DIMLIN2> 
  FlatTensor (FlatTensor<DIM,T,DIMLIN2> t2)
    : size(t2.GetSize()), dist(t2.GetDist()), sub(t2.GetSubTensor())
  {
    static_assert(DIMLIN <= DIMLIN2, "illegal tensor copy");
    ;
  }

  
  template <typename ... ARG>
  auto operator() (int i, ARG ... args) -> decltype ( OffsetTensor (sub(args...), i*dist) )
  {
    return OffsetTensor (sub(args...), i*dist);
  }

  template <typename ... ARG>
  auto operator() (typestar star, ARG ... args) 
    -> decltype (LargerTensor (sub(args...), size,dist) )
  {
    return LargerTensor(sub(args...), size, dist);
  }

  FlatTensor operator= (double d)
  {
    for (int i = 0; i < size; i++)
      GetSubTensor(i) = d;
    return *this;
  }

  int GetSize () const { return size; }
  int GetDist () const { return dist; }
  int GetTotalSize () const { return size*sub.GetTotalSize(); }

  auto GetSubTensor() const -> decltype(sub)
  { return sub; } 
  FlatTensor<DIM-1,T> GetSubTensor (int i) const 
  {
    return OffsetTensor (sub, i*dist);
  }

  T *& Data () { return sub.Data(); }
  T * Data () const { return sub.Data(); }
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
              
  int GetSize () const { return 1; }
  int GetDist () const { return 1; }
  int GetTotalSize () const { return 1; }
  T *& Data () { return data; }
  T * Data () const { return data; }  
  FlatTensor<0,T,LINDIM> operator() () { return FlatTensor<0,T,LINDIM> (data); }
  operator T  () const { return *data; }
  operator T& () { return *data; }
  void operator= (double d) { *data = d; }
  void operator-= (double d) { *data -= d; }
  void operator+= (double d) { *data += d; }
};













template <int DIM, typename T = double>
class Tensor : public FlatTensor<DIM,T>
{
public:
  
  template <typename ... ARG>
  Tensor (ARG ... args) : FlatTensor<DIM,T> (args...)
  {
    int totsize = this->GetTotalSize();
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

  for (int i = 0; i < tensor.GetSize(); i++)
    {
      ost << "subtensor " << i;
      for (int j = 0; j < DIM-1; j++) ost << ",*";
      ost << ":" << endl;
      ost << tensor.GetSubTensor(i);
    }

  return ost;
}

template <typename T>
INLINE ostream & operator<< (ostream & ost, const FlatTensor<2,T> & tensor)
{
  for (int i = 0; i < tensor.GetSize(); i++)
    ost << tensor.GetSubTensor(i);
  return ost;
}

template <typename T>
INLINE ostream & operator<< (ostream & ost, const FlatTensor<1,T> & tensor)
{
  for (int i = 0; i < tensor.GetSize(); i++)
    ost << * (tensor.Data()+i*tensor.GetDist()) << " ";
  return ost << endl;
}


template <typename T>
INLINE ostream & operator<< (ostream & ost, const FlatTensor<0,T> & tensor)
{
  return ost << *tensor.Data();
}




