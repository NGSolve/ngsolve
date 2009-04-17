#ifndef FILE_MATRIX_EXPR
#define FILE_MATRIX_EXPR

/**************************************************************************/
/* File:   matrix.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jan. 02                                                    */
/**************************************************************************/


template <int H, int W, typename T> class Mat;
template <typename T> class SliceMatrix;

///
extern void CheckMatRange(int h, int w, int i);
///
extern void CheckMatRange(int h, int w, int i, int j);



/**
   A simple matrix.
   Has height, width and data-pointer. 
   No memory allocation/deallocation. User must provide memory.
 */
template <typename T = double>
class FlatMatrix : public CMCPMatExpr<FlatMatrix<T> >
{
protected:
  /// the height
  int h;
  /// the width
  int w;
  /// the data
  T * data;
public:

  /// element type
  typedef T TELEM;
  /// scalar type of elements (double or Complex)
  typedef typename mat_traits<T>::TSCAL TSCAL;


  /// nothing done in default constructor
  FlatMatrix () throw () { ; }
  
  /// set height, width, and mem
  FlatMatrix (int ah, int aw, T * adata) throw ()
    : h(ah), w(aw), data(adata) { ; }
  
  /// set height = width, and mem
  FlatMatrix (int ah, T * adata) throw()
    : h(ah), w(ah), data(adata) { ; }

  /// allocates at local heap
  FlatMatrix (int ah, int aw, LocalHeap & lh) throw (LocalHeapOverflow)
    : h(ah), w(aw), data((T*)lh.Alloc(ah*aw*sizeof(T))) { ; }
  
  /// allocates at local heap
  FlatMatrix (int ah, LocalHeap & lh) throw (LocalHeapOverflow)
    : h(ah), w(ah), data((T*)lh.Alloc(ah*ah*sizeof(T))) { ; }
  
  /// copy constructor. copies pointers, not contents
  FlatMatrix (const FlatMatrix & m) throw () 
    : CMCPMatExpr<FlatMatrix> (), h(m.h), w(m.w) , data(m.data)
  { ; }
  

  /// useful to put FlatMatrix over other matrix
  template <typename T2>
  explicit FlatMatrix (const MatExpr<T2> & m)
    : h(m.Height()), w(m.Width()),
      data(const_cast<T*>(&m.Spec()(0,0)))  
  { ; }
  
  /// useful to put FlatMatrix over other Mat
  template <int H, int W>
  FlatMatrix (const Mat<H,W,TSCAL> & m) throw()
    : h(H), w(W), data(const_cast<T*>(&m(0,0)))
  { ; }

  /// do nothing 
  ~FlatMatrix () throw() { ; }

  /// set size, and assign mem
  void AssignMemory (int ah, int aw, LocalHeap & lh)  throw (LocalHeapOverflow)
  {
    h = ah;
    w = aw;
    data = (T*)lh.Alloc(h*w*sizeof(T));
  }
  
  /// set size, and assign mem
  void AssignMemory (int ah, int aw, T * mem) throw()
  {
    h = ah;
    w = aw;
    data = mem;
  }
  

  /// assign contents
  template<typename TBxx>
  const FlatMatrix & operator= (const Expr<TBxx> & m) const
  {
    return CMCPMatExpr<FlatMatrix>::operator= (m);
  }

  /// copy contents
  const FlatMatrix & operator= (const FlatMatrix & m) throw()
  {
    int hw = h*w;
    for (int i = 0; i < hw; i++)
      data[i] = m(i);
    return *this;
  }

  /// assign constant
  FlatMatrix & operator= (TSCAL s) throw()
  {
    int hw = h*w;
    for (int i = 0; i < hw; i++)
      data[i] = s; 
    return *this;
  }

  /// copy size and pointers
  FlatMatrix & Assign (const FlatMatrix & m) throw()
  {
    h = m.h;
    w = m.w;
    data = m.data;
    return *this;
  }


  /// access operator, linear access
  TELEM & operator() (int i) const 
  { 
#ifdef CHECK_RANGE
    CheckMatRange(h,w,i);
#endif
    return data[i]; 
  }

  /// access operator
  TELEM & operator() (int i, int j) const
  {
#ifdef CHECK_RANGE
    CheckMatRange(h,w,i,j);
#endif
    return data[i*w+j]; 
  }

  /*
  /// access operator, linear access
  const TELEM & operator() (int i) const
  {
#ifdef CHECK_RANGE
    CheckMatRange(h,w,i);
#endif
    return data[i]; 
  }

  /// access operator
  const TELEM & operator() (int i, int j) const
  {
#ifdef CHECK_RANGE
    CheckMatRange(h,w,i,j);
#endif
    return data[i*w+j]; 
  }
  */

  /// the height
  int Height () const throw() { return h; }

  /// the width
  int Width () const throw() { return w; }

  /*
  FlatVector<T> Row (int i)
  {
    return FlatVector<T> (w, &data[i*w]);
  }
  */

  const FlatVector<T> Row (int i) const
  {
    return FlatVector<T> (w, &data[i*w]);
  }

  const SliceVector<T> Col (int i)
  {
    return SliceVector<T> (h, w, &data[i]);
  }

  const FlatMatrix VRange (int first, int next) const
  {
    return FlatMatrix (next-first, w, data+first*w);
  }

  const SliceMatrix<T> HRange (int first, int next) const
  {
    return SliceMatrix<T> (h, next-first, w, data+first);
  }
  
};



/**
   A Matrix class with memory allocation/deallocation
 */
template <typename T = double>
class Matrix : public FlatMatrix<T>
{
public:

  /// element type
  typedef T TELEM;
  /// scalar type of elements (double or Complex)
  typedef typename mat_traits<T>::TSCAL TSCAL;

  /// default constructor
  Matrix () throw () : FlatMatrix<T> (0, 0) { ; }

  /// allocate matrix of size ah * ah
  Matrix (int ah) : FlatMatrix<T> (ah, new T[ah*ah]) { ; }

  /// allocate matrix of size ah * aw
  Matrix (int ah, int aw) : FlatMatrix<T> (ah, aw, new T[ah*aw]) { ; }

  /// delete memory
  ~Matrix() { delete [] this->data; }

  /// sets new size of matrix
  void SetSize(int ah, int aw)
  {
    if (this->h == ah && this->w == aw) return;
    delete [] this->data;
    this->h = ah;
    this->w = aw;
    this->data = new T[this->h * this->w];
  }

  /// sets new size of matrix
  void SetSize(int ah)
  {
    if (this->h == ah && this->w == ah) return;
    delete [] this->data;
    this->h = ah;
    this->w = ah;
    this->data = new T[this->h * this->w];
  }

  /// assign matrix, sizes must match
  template<typename TB>
  Matrix & operator= (const Expr<TB> & m)
  {
    CMCPMatExpr<FlatMatrix<T> >::operator= (m);
    return *this;
  }

  /// fill matrix with scalar
  Matrix & operator= (TSCAL s)
  {
    FlatMatrix<T>::operator= (s);
    return *this;
  }
};




/**
   A matrix of fixed size.
   Useful as entry type in system matrices,...
 */
template <int H, int W = H, typename T = double>
class Mat : public MatExpr<Mat<H,W,T> >
{
  T data[H*W];
public:
  typedef T TELEM;
  typedef typename mat_traits<T>::TSCAL TSCAL;
  typedef Vec<H, typename mat_traits<T>::TV_COL> TV_COL;
  typedef Vec<W, typename mat_traits<T>::TV_ROW> TV_ROW;
  enum { HEIGHT = H };
  enum { WIDTH  = W };

  /// do not initialize 
  Mat () throw () { ; }

  /// copy matrix
  Mat (const Mat & m) throw()
    : MatExpr<Mat> ()
  {
    (*this) = m;
  }

  /// assign values
  template<typename TB>
  Mat (const Expr<TB> & m)
  {
    MatExpr<Mat>::operator= (m);
  }

  /// fill with scalar
  Mat (TSCAL s) throw()     // explicit removed JS, Aug '07
  {
    for (int i = 0; i < H*W; i++)
      data[i] = s;
  }

  /// assign values
  template<typename TB>
  Mat & operator= (const Expr<TB> & m)
  {
    MatExpr<Mat>::operator= (m);
    return *this;
  }

  /// copy matrix
  Mat & operator= (const Mat & m) throw()
  {
    for (int i = 0; i < H*W; i++)
      data[i] = m.data[i];
    return *this;
  }
 
  /// fill values
  Mat & operator= (TSCAL s) throw()
  {
    for (int i = 0; i < H*W; i++)
      data[i] = s;
    return *this;
  }

  /// linear access
  TELEM & operator() (int i) { return data[i]; }
  /// access element
  TELEM & operator() (int i, int j) { return data[i*W+j]; }
  /// linear access
  const TELEM & operator() (int i) const { return data[i]; }
  /// access element
  const TELEM & operator() (int i, int j) const { return data[i*W+j]; }

  /// the height
  int Height () const throw() { return H; }
  /// the width
  int Width () const throw() { return W; }

  ///
  const FlatVec<W,T> Row (int i) 
  {
    return FlatVec<W,T> (&(*this)(i,0));
  }

  ///
  const FlatVec<W,const T> Row (int i) const
  {
    return FlatVec<W,const T> (&(*this)(i,0));
  }

  ///
  const FixSliceVector<W, T> Col (int i) 
  {
    return FixSliceVector<W,T> (H, &(*this)(0,i));
  }

};


template <class TB> ALWAYS_INLINE
Mat<3,3,double> operator+= (class Mat<3,3,double> & m, const Expr<TB> & exp)
{
  m(0,0) += exp.Spec()(0,0);
  m(0,1) += exp.Spec()(0,1);
  m(0,2) += exp.Spec()(0,2);
  m(1,0) += exp.Spec()(1,0);
  m(1,1) += exp.Spec()(1,1);
  m(1,2) += exp.Spec()(1,2);
  m(2,0) += exp.Spec()(2,0);
  m(2,1) += exp.Spec()(2,1);
  m(2,2) += exp.Spec()(2,2);
  return m;
}

template <class TB> ALWAYS_INLINE
Mat<2,2,double> operator+= (class Mat<2,2,double> & m, const Expr<TB> & exp)
{
  m(0,0) += exp.Spec()(0,0);
  m(0,1) += exp.Spec()(0,1);
  m(1,0) += exp.Spec()(1,0);
  m(1,1) += exp.Spec()(1,1);
  return m;
}












/**
   A Matrix with width known at compile time
   No memory allocation/deallocation. User must provide memory.
 */
template <int W, typename T = double>
class FlatMatrixFixWidth : public CMCPMatExpr<FlatMatrixFixWidth<W,T> >
{
protected:
  /// the data
  T * data;
  /// the height
  int h;
public:
  /// entry type
  typedef T TELEM;
  /// scalar type of entry
  typedef typename mat_traits<T>::TSCAL TSCAL;

  /// nothing done in default constructor
  FlatMatrixFixWidth () throw()
    : data(0), h(0) { ; }
  
  /// set height and mem
  FlatMatrixFixWidth (int ah, T * adata) throw()
    : data(adata), h(ah) { ; }

  /// allocates at local heap
  FlatMatrixFixWidth (int ah, LocalHeap & lh) throw (LocalHeapOverflow)
    : data((T*)lh.Alloc(ah*W*sizeof(T))), h(ah) { ; }
  
  /// copy constructor. copies pointers, not contents
  FlatMatrixFixWidth (const FlatMatrixFixWidth & m) throw()
    : CMCPMatExpr<FlatMatrixFixWidth> (), data(m.data), h(m.h)
  { ; }

  /// copy constructor. copies pointers, not contents
  FlatMatrixFixWidth (const FlatMatrix<TSCAL> & m) throw()
    : data(const_cast<T*> (&m(0))), h(m.Height())
  { ; }

  template <int H>
  FlatMatrixFixWidth (const Mat<H,W,TSCAL> & m) throw()
    : data(const_cast<T*>(&m(0,0))), h(H)
  { ; }  

  /// do nothing 
  ~FlatMatrixFixWidth () throw() { ; }

  /// set size, and assign mem
  void AssignMemory (int ah, LocalHeap & lh) throw (LocalHeapOverflow)
  {
    h = ah;
    data = (T*)lh.Alloc(h*W*sizeof(T));
  }
  
  /// set size, and assign mem
  void AssignMemory (int ah, T * mem) throw()
  {
    h = ah;
    data = mem;
  }
  

  /// assign contents
  template<typename TB>
  const FlatMatrixFixWidth & operator= (const Expr<TB> & m) const
  {
    return CMCPMatExpr<FlatMatrixFixWidth>::operator= (m);
  }

  /// copy contents
  FlatMatrixFixWidth & operator= (const FlatMatrixFixWidth & m) throw()
  {
    for (int i = 0; i < h*W; i++)
      data[i] = m(i);
    return *this;
  }

  /// assign constant
  FlatMatrixFixWidth & operator= (TSCAL s) throw()
  {
    for (int i = 0; i < h*W; i++)
      data[i] = s; 
    return *this;
  }

  /// copy size and pointers
  FlatMatrixFixWidth & Assign (const FlatMatrixFixWidth & m) throw()
  {
    h = m.h;
    data = m.data;
    return *this;
  }

  /*
  /// access operator, linear access
  TELEM & operator() (int i)
  { 
#ifdef CHECK_RANGE
    CheckMatRange(h,W,i);
#endif
    return data[i]; 
  }

  /// access operator
  TELEM & operator() (int i, int j) 
  {
#ifdef CHECK_RANGE
    CheckMatRange(h,W,i,j);
#endif
    return data[i*W+j]; 
  }
  */

  /// access operator, linear access
  TELEM & operator() (int i) const
  {
#ifdef CHECK_RANGE
    CheckMatRange(h,W,i);
#endif
    return data[i]; 
  }

  /// access operator
  TELEM & operator() (int i, int j) const
  {
#ifdef CHECK_RANGE
    CheckMatRange(h,W,i,j);
#endif
    return data[i*W+j]; 
  }

  /// the height
  int Height () const throw() { return h; }

  /// the width
  int Width () const throw() { return W; }

  /*
  ///
  FlatVec<W,T> Row (int i) 
  {
    return FlatVec<W,T> (&(*this)(i,0));
  }
  */

  ///
  const FlatVec<W,T> Row (int i) const
  {
    return FlatVec<W,T> (&(*this)(i,0));
  }

  const FixSliceVector<W,T> Col (int i)
  {
    return FixSliceVector<W,T> (h, &data[i]);
  }


  const FlatMatrixFixWidth VRange (int first, int next) const
  {
    return FlatMatrixFixWidth (next-first, data+first*W);
  }

  const SliceMatrix<T> HRange (int first, int next) const
  {
    return SliceMatrix<T> (h, next-first, W, data+first);
  }
};


/**
   A Matrix class with memory allocation/deallocation
 */
template <int W, typename T = double>
class MatrixFixWidth : public FlatMatrixFixWidth<W,T>
{
public:
  typedef typename mat_traits<T>::TSCAL TSCAL;

  MatrixFixWidth () { ; }

  /// allocate matrix of size ah * ah
  MatrixFixWidth (int ah) : FlatMatrixFixWidth<W,T> (ah, new T[ah*W]) { ; }

  /// delete memory
  ~MatrixFixWidth() { delete [] this->data; }

  /// sets new size of matrix
  void SetSize(int ah)
  {
    if (this->h == ah) return;
    delete [] this->data;
    this->h = ah;
    this->data = new T[this->h * W];
  }

  /// assign matrix, sizes must match
  template<typename TB>
  MatrixFixWidth & operator= (const Expr<TB> & m)
  {
    MatExpr<FlatMatrixFixWidth<W,T> >::operator= (m);
    return *this;
  }

  /// fill matrix with scalar
  MatrixFixWidth & operator= (TSCAL s)
  {
    FlatMatrixFixWidth<W,T>::operator= (s);
    return *this;
  }
};






















/**
   A Matrix which height is known at compile time
   No memory allocation/deallocation. User must provide memory.
   Matrix is stored colum-wise
 */
template <int H, typename T = double>
class FlatMatrixFixHeight : public MatExpr<FlatMatrixFixHeight<H,T> >
{
protected:
  /// the data
  T * data;
  /// the width
  int w;
public:
  ///
  typedef T TELEM;
  ///
  typedef typename mat_traits<T>::TSCAL TSCAL;

  /// 
  FlatMatrixFixHeight () throw()
    : data(0), w(0) { ; }
  
  /// set height and mem
  FlatMatrixFixHeight (int aw, T * adata) throw()
    : data(adata), w(aw) { ; }

  /// allocates at local heap
  FlatMatrixFixHeight (int aw, LocalHeap & lh) throw (LocalHeapOverflow)
    : data((T*)lh.Alloc(aw*H*sizeof(T))), w(aw) { ; }
  
  /// copy constructor. copies pointers, not contents
  FlatMatrixFixHeight (const FlatMatrixFixHeight & m)
    : data(m.data), w(m.w)
  { ; }
  

  /// do nothing 
  ~FlatMatrixFixHeight () throw() { ; }

  /// set size, and assign mem
  void AssignMemory (int aw, LocalHeap & lh) throw (LocalHeapOverflow)
  {
    w = aw;
    data = (T*)lh.Alloc(w*H*sizeof(T));
  }
  
  /// set size, and assign mem
  void AssignMemory (int aw, T * mem) 
  {
    w = aw;
    data = mem;
  }
  

  /// assign contents
  template<typename TB>
  FlatMatrixFixHeight & operator= (const Expr<TB> & m)
  {
    for (int j = 0; j < w; j++)
      for (int i = 0; i < H; i++)
	(*this)(i,j) = m.Spec()(i,j);
    return *this;
  }

  /// copy contents
  FlatMatrixFixHeight & operator= (const FlatMatrixFixHeight & m)
  {
    for (int i = 0; i < w*H; i++)
      data[i] = m(i);
    return *this;
  }

  /// assign constant
  FlatMatrixFixHeight & operator= (TSCAL s)
  {
    for (int i = 0; i < w*H; i++)
      data[i] = s; 
    return *this;
  }

  /// copy size and pointers
  FlatMatrixFixHeight & Assign (const FlatMatrixFixHeight & m)
  {
    w = m.w;
    data = m.data;
    return *this;
  }


  /// access operator, linear access
  TELEM & operator() (int i)
  { 
#ifdef CHECK_RANGE
    CheckMatRange(H,w,i);
#endif
    return data[i]; 
  }

  /// access operator
  TELEM & operator() (int i, int j) 
  {
#ifdef CHECK_RANGE
    CheckMatRange(H,w,i,j);
#endif
    return data[i+j*H]; 
  }


  /// access operator, linear access
  const TELEM & operator() (int i) const
  {
#ifdef CHECK_RANGE
    CheckMatRange(H,w,i);
#endif
    return data[i]; 
  }

  /// access operator
  const TELEM & operator() (int i, int j) const
  {
#ifdef CHECK_RANGE
    CheckMatRange(H,w,i,j);
#endif
    return data[i+j*H]; 
  }


  /// the height
  int Height () const { return H; }

  /// the width
  int Width () const { return w; }


  const FlatVec<H,T> Col (int i) 
  {
    return FlatVec<H,T> ( data+i*H ); 
  }


  const SliceVector<T> Row (int i)
  {
    return SliceVector<T> (w, H, &data[i]);
  }


  enum { IS_LINEAR = 0 };
  //  static bool IsLinear() { return false; }
};


/**
   A Matrix class with memory allocation/deallocation
 */
template <int H, typename T = double>
class MatrixFixHeight : public FlatMatrixFixHeight<H,T>
{
public:
  typedef typename mat_traits<T>::TSCAL TSCAL;

  /// allocate matrix of size ah * ah
  MatrixFixHeight (int aw) : FlatMatrixFixHeight<H,T> (aw, new T[aw*H]) { ; }

  /// delete memory
  ~MatrixFixHeight() { delete [] this->data; }

  /// sets new size of matrix
  void SetSize(int aw)
  {
    if (this->w == aw) return;
    delete [] this->data;
    this->w = aw;
    this->data = new T[H*this->w];
  }

  /// assign matrix, sizes must match
  template<typename TB>
  MatrixFixHeight & operator= (const Expr<TB> & m)
  {
    int k = 0;
    for (int j = 0; j < this->w; j++)
      for (int i = 0; i < H; i++)
	{
	  this->data[k] = m.Spec()(i,j);
	  k++;
	}
    return *this;
  }

  /// fill matrix with scalar
  MatrixFixHeight & operator= (TSCAL s)
  {
    FlatMatrixFixHeight<H,T>::operator= (s);
    return *this;
  }
};












template <typename T = double>
class SliceMatrix : public CMCPMatExpr<SliceMatrix<T> >
{
protected:
  /// the height
  int h;
  /// the width
  int w;
  /// the distance
  int dist;
  /// the data
  T * data;
public:

  /// element type
  typedef T TELEM;
  /// scalar type of elements (double or Complex)
  typedef typename mat_traits<T>::TSCAL TSCAL;
  enum { IS_LINEAR = 0 };

  /// set height, width, and mem
  SliceMatrix (int ah, int aw, int adist, T * adata) throw ()
    : h(ah), w(aw), dist(adist), data(adata) { ; }
  
  /// assign contents
  template<typename TB>
  const SliceMatrix & operator= (const Expr<TB> & m) const
  {
    return CMCPMatExpr<SliceMatrix>::operator= (m);
  }

  /// assign constant
  SliceMatrix & operator= (TSCAL s) throw()
  {
    for (int i = 0; i < h; i++)
      for (int j = 0; j < w; j++)
        data[i*dist+j] = s;
    return *this;
  }

  /// access operator
  TELEM & operator() (int i, int j) const
  {
#ifdef CHECK_RANGE
    CheckMatRange(h,w,i,j);
#endif
    return data[i*dist+j]; 
  }

  /// access operator, linear access
  TELEM & operator() (int i) const
  {
#ifdef CHECK_RANGE
    CheckMatRange(h,w,i);
#endif
    return data[i]; 
  }

  /// the height
  int Height () const throw() { return h; }

  /// the width
  int Width () const throw() { return w; }


  const SliceMatrix VRange (int first, int next) const
  {
    return SliceMatrix (next-first, w, dist, data+first*dist);
  }

  const SliceMatrix<T> HRange (int first, int next) const
  {
    return SliceMatrix<T> (h, next-first, dist, data+first);
  }

};


















/**
   Identity Matrix of fixed size
 */
template <int H>
class Id : public MatExpr<Id<H> >
{
public:
  typedef double TELEM;
  typedef double TSCAL;
  typedef Vec<H, double> TV_COL;
  typedef Vec<H, double> TV_ROW;
  enum { HEIGHT = H };
  enum { WIDTH  = H };
  // static bool IsLinear() { return 0; }
  enum { IS_LINEAR = 0 };

  /// do not initialize 
  Id () { ; }

  ///
  double operator() (int i) const
  { cerr << "id, linear access" << endl; return 0; }
  ///
  double operator() (int i, int j) const { return (i == j) ? 1 : 0; }

  /// the height
  int Height () const { return H; }
  /// the width
  int Width () const { return H; }
};






template<int H, int W, typename T>
inline std::ostream & operator<< (std::ostream & s, const Mat<H,W,T> & m)
{
  for (int i = 0; i < H*W; i++)
    s << " " << setw(7) << m(i);
  return s;
}



#endif
