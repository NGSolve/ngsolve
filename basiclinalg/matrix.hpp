#ifndef FILE_MATRIX_EXPR
#define FILE_MATRIX_EXPR

/**************************************************************************/
/* File:   matrix.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jan. 02                                                    */
/**************************************************************************/

namespace ngbla
{

  template <int H, int W, typename T> class Mat;
  template <typename T> class SliceMatrix;
  template <typename T> class DoubleSliceMatrix;


  ///
  extern void CheckMatRange(int h, int w, int i);
  ///
  extern void CheckMatRange(int h, int w, int i, int j);


  // compatibility (for a while)
  // #define VRange Rows
  // #define HRange Cols

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
    T * __restrict data;
  public:

    /// element type
    typedef T TELEM;
    typedef T& TREF;
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
    FlatMatrix (int ah, int aw, LocalHeap & lh) // throw (LocalHeapOverflow)
    // : h(ah), w(aw), data((T*)lh.Alloc(ah*aw*sizeof(T))) { ; }
      : h(ah), w(aw), data (lh.Alloc<T>(ah*aw)) { ; }
  
    /// allocates at local heap
    FlatMatrix (int ah, LocalHeap & lh) // throw (LocalHeapOverflow)
    // : h(ah), w(ah), data((T*)lh.Alloc(ah*ah*sizeof(T))) { ; }
      : h(ah), w(ah), data(lh.Alloc<T>(ah*ah)) { ; }
  
    /// copy constructor. copies pointers, not contents
    FlatMatrix (const FlatMatrix & m) throw () 
      : CMCPMatExpr<FlatMatrix> (), h(m.h), w(m.w) , data(m.data)
    { ; }
  
    /// allocate and compute 
    template<typename TB>
    FlatMatrix (const LocalHeapExpr<TB> & m2) 
    {
      h = m2.A().Height();
      w = m2.A().Width();
      LocalHeap & lh = m2.GetLocalHeap();
      data = lh.Alloc<T> (h*w);
      CMCPMatExpr<FlatMatrix<T> >::operator= (m2.A());
    }

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
    void AssignMemory (int ah, int aw, LocalHeap & lh)  // throw (LocalHeapOverflow)
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

    /*
    template <typename TA, typename TB>
    const FlatMatrix & operator= (const LapackProduct<TA, TB> & prod) const
    {
      LapackMult (prod.A(), prod.B(), *this);
      return *this;
    }
    */

    /// copy contents
    const FlatMatrix & operator= (const FlatMatrix & m) const 
    {
      int hw = h*w;
      for (int i = 0; i < hw; i++)
        data[i] = m(i);
      return *this;
    }

    /// assign constant
    const FlatMatrix & operator= (TSCAL s) const 
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
      return data[i*size_t(w)+j]; 
    }

    /// the height
    int Height () const throw() { return h; }

    /// the width
    int Width () const throw() { return w; }


    const FlatVector<T> Row (int i) const
    {
      return FlatVector<T> (w, &data[i*size_t(w)]);
    }

#ifdef FLATVECTOR_WITH_DIST
    const FlatVector<T> Col (int i) const
    {
      return FlatVector<T> (h, w, &data[i]);
    }
#else
    const SliceVector<T> Col (int i) const
    {
      return SliceVector<T> (h, w, &data[i]);
    }
#endif

    const SliceVector<T> Diag () const
    {
      return SliceVector<T> (h, w+1, &data[0]);
    }

    using CMCPMatExpr<FlatMatrix<T> >::Rows;
    using CMCPMatExpr<FlatMatrix<T> >::Cols;

    const FlatMatrix Rows (int first, int next) const
    {
      return FlatMatrix (next-first, w, data+first*w);
    }

    const SliceMatrix<T> Cols (int first, int next) const
    {
      return SliceMatrix<T> (h, next-first, w, data+first);
    }

    const FlatMatrix Rows (IntRange range) const
    {
      return FlatMatrix (range.Next()-range.First(), w, data+range.First()*w);
    }

    const SliceMatrix<T> Cols (IntRange range) const
    {
      return SliceMatrix<T> (h, range.Next()-range.First(), w, data+range.First());
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
    Matrix (int ah) : FlatMatrix<T> (ah, new T[size_t(ah)*size_t(ah)]) { ; }
    
    /// allocate matrix of size ah * aw
    Matrix (int ah, int aw) : FlatMatrix<T> (ah, aw, new T[size_t(ah)*size_t(aw)]) { ; }

    /// allocate and copy matrix  
    Matrix (const Matrix & m2) 
      : FlatMatrix<T> (m2.Height(), m2.Width(), new T[m2.Height()*m2.Width()]) 
    {
      FlatMatrix<T>::operator= (m2);
    }

    /// allocate and compute 
    template<typename TB>
    Matrix (const Expr<TB> & m2) 
      : FlatMatrix<T> (m2.Height(), m2.Width(), new T[m2.Height()*m2.Width()]) 
    {
      CMCPMatExpr<FlatMatrix<T> >::operator= (m2);
    }


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
      SetSize (m.Height(), m.Width());
      CMCPMatExpr<FlatMatrix<T> >::operator= (m);
      return *this;
    }

    /*
    template <typename TA, typename TB>
    Matrix & operator= (const LapackProduct<TA, TB> & prod) 
    {
      // SetSize (m.Height(), m.Width());
      LapackMult (prod.A(), prod.B(), *this);
      return *this;
    }
    */

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
    T data[(H*W>0) ? H*W : 1];
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
    const FixSliceVector<W,T> Col (int i) 
    {
      return FixSliceVector<W,T> (H, &(*this)(0,i));
    }

    ///
    const FixSliceVector<W,const T> Col (int i) const
    {
      return FixSliceVector<W,const T> (H, &(*this)(0,i));
    }

  };







  /**
     A diagonal matrix of fixed size.
  */
  template <int H, typename T = double>
  class DiagMat : public MatExpr<DiagMat<H,T> >
  {
    T data[(H>0) ? H : 1];
  public:
    typedef T TELEM;
    typedef typename mat_traits<T>::TSCAL TSCAL;
    typedef Vec<H, typename mat_traits<T>::TV_COL> TV_COL;
    typedef Vec<H, typename mat_traits<T>::TV_ROW> TV_ROW;
    enum { HEIGHT = H };
    enum { WIDTH  = H };

    /// do not initialize 
    DiagMat () { ; }

    /// copy matrix
    DiagMat (const DiagMat & m) 
      : MatExpr<DiagMat> ()
    {
      (*this) = m;
    }

    /// assign values
    template<typename TB>
    DiagMat (const Expr<TB> & m)
    {
      (*this) = m;
    }

    /// fill with scalar
    DiagMat (TSCAL s) throw()    
    {
      for (int i = 0; i < H; i++)
        data[i] = s;
    }

    /// assign values
    template<typename TB>
    DiagMat & operator= (const Expr<TB> & m)
    {
      for (int i = 0; i < H; i++)
        data[i] = m.Spec()(i,i);
      return *this;
    }

    /// copy matrix
    DiagMat & operator= (const DiagMat & m) throw()
    {
      for (int i = 0; i < H; i++)
        data[i] = m.data[i];
      return *this;
    }
 
    /// fill values
    DiagMat & operator= (TSCAL s) throw()
    {
      for (int i = 0; i < H; i++)
        data[i] = s;
      return *this;
    }

    /// linear access
    TELEM & operator() (int i) { return data[i]; }
    /// access element
    // TELEM & operator() (int i, int j) { return (i==j) ? data[i] : 0; }
    /// linear access
    TELEM operator() (int i) const { return data[i]; }
    TELEM operator[] (int i) const { return data[i]; }
    /// access element
    TELEM operator() (int i, int j) const { return (i==j) ? data[i] : 0; }

    /// the height
    int Height () const throw() { return H; }
    /// the width
    int Width () const throw() { return H; }
    
    enum { IS_LINEAR = 0 };
  };


















  /**
     A Matrix with width known at compile time
     No memory allocation/deallocation. User must provide memory.
  */
  template <int W, typename T = double>
  class FlatMatrixFixWidth : public CMCPMatExpr<FlatMatrixFixWidth<W,T> >
  {
  protected:
    /// the data
    T *  __restrict data;
    /// the height
    int h;
  public:
    /// entry type
    typedef T TELEM;
    typedef T& TREF;
    /// scalar type of entry
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /// nothing done in default constructor
    FlatMatrixFixWidth () throw() { ; }
  
    /// set height and mem
    FlatMatrixFixWidth (int ah, T * adata) throw()
      : data(adata), h(ah) { ; }

    /// allocates at local heap
    FlatMatrixFixWidth (int ah, LocalHeap & lh) // throw (LocalHeapOverflow)
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
    void AssignMemory (int ah, LocalHeap & lh) // throw (LocalHeapOverflow)
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
    const FlatMatrixFixWidth & operator= (const FlatMatrixFixWidth & m) const throw()
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

    ///
    operator const FlatMatrix<T>() const { return FlatMatrix<T> (h, W, data); }

    ///
    const FlatVec<W,T> Row (int i) const
    {
      return FlatVec<W,T> (&(*this)(i,0));
    }

    const FixSliceVector<W,T> Col (int i)
    {
      return FixSliceVector<W,T> (h, &data[i]);
    }

    using CMCPMatExpr<FlatMatrixFixWidth<W,T> >::Rows;
    using CMCPMatExpr<FlatMatrixFixWidth<W,T> >::Cols;

    const FlatMatrixFixWidth Rows (int first, int next) const
    {
      return FlatMatrixFixWidth (next-first, data+first*W);
    }

    const SliceMatrix<T> Cols (int first, int next) const
    {
      return SliceMatrix<T> (h, next-first, W, data+first);
    }

    const FlatMatrixFixWidth Rows (IntRange range) const
    {
      return Rows (range.First(), range.Next());
    }

    const SliceMatrix<T> Cols (IntRange range) const
    {
      return Cols (range.First(), range.Next());
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
    T * __restrict data;
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
    FlatMatrixFixHeight (int aw, LocalHeap & lh) // throw (LocalHeapOverflow)
      : data((T*)lh.Alloc(aw*H*sizeof(T))), w(aw) { ; }
  
    /// copy constructor. copies pointers, not contents
    FlatMatrixFixHeight (const FlatMatrixFixHeight & m)
      : data(m.data), w(m.w)
    { ; }
  

    /// do nothing 
    ~FlatMatrixFixHeight () throw() { ; }

    /// set size, and assign mem
    void AssignMemory (int aw, LocalHeap & lh) // throw (LocalHeapOverflow)
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


    const FlatVec<H,T> Col (int i) const
    {
      return FlatVec<H,T> ( data+i*H ); 
    }

    const SliceVector<T> Row (int i) const
    {
      return SliceVector<T> (w, H, &data[i]);
    }


    SubMatrixExpr<FlatMatrixFixHeight,true>
    Rows (int first, int next)
    { 
      return SubMatrixExpr<FlatMatrixFixHeight,true> (*this, first, 0, next-first, Width()); 
    }
    SubMatrixExpr<FlatMatrixFixHeight,true>
    Rows (IntRange range) 
    { 
      return Rows (range.First(), range.Next());
    }
    
    /*
    const DoubleSliceMatrix<T> Rows (int first, int next) const
    {
      return DoubleSliceMatrix<T> (next-first, w, 1, H, data+first);
    }

    const FlatMatrixFixHeight<H,T> Cols (int first, int next) const
    {
      return FlatMatrixFixHeight<H,T> (next-first, data+first*H);
    }

    const DoubleSliceMatrix<T> Rows (IntRange range) const
    {
      return Rows (range.First(), range.Next());
    }

    const FlatMatrixFixHeight<H,T> Cols (IntRange range) const
    {
      return Cols (range.First(), range.Next());
    }
    */


    const FlatMatrixFixWidth<H,T> Trans () const
    {
      return FlatMatrixFixWidth<H,T> (w, data);
    }



    enum { IS_LINEAR = 0 };
    enum { COL_MAJOR = 1 }; 
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
    T * __restrict data;
  public:

    /// element type
    typedef T TELEM;
    /// scalar type of elements (double or Complex)
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { IS_LINEAR = 0 };

    /// set height, width, and mem
    SliceMatrix (int ah, int aw, int adist, T * adata) throw ()
      : h(ah), w(aw), dist(adist), data(adata) { ; }
  
    SliceMatrix (const FlatMatrix<T> & mat)
      : h(mat.Height()), w(mat.Width()), dist(mat.Width()), data(&mat(0,0))
    { ; }

    /// assign contents
    template<typename TB>
    const SliceMatrix & operator= (const Expr<TB> & m) const
    {
      return CMCPMatExpr<SliceMatrix>::operator= (m);
    }

    /*
    template <typename TA, typename TB>
    const SliceMatrix & operator= (const LapackProduct<TA, TB> & prod) const
    {
      LapackMult (prod.A(), prod.B(), *this);
      return *this;
    }
    */

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

    /// 
    int Dist () const throw() { return dist; }

    const SliceMatrix Rows (int first, int next) const
    {
      return SliceMatrix (next-first, w, dist, data+first*dist);
    }

    const FlatVector<T> Row (int i) const
    {
      return FlatVector<T> (w, &data[i*size_t(dist)]);
    }

    const SliceVector<T> Col (int i) const
    {
      return SliceVector<T> (h, dist, &data[i]);
    }

    const SliceMatrix<T> Cols (int first, int next) const
    {
      return SliceMatrix<T> (h, next-first, dist, data+first);
    }

    const SliceMatrix Rows (IntRange range) const
    {
      return Rows (range.First(), range.Next());
    }

    const SliceMatrix<T> Cols (IntRange range) const
    {
      return Cols (range.First(), range.Next());
    }

  };





  template <typename T = double>
  class DoubleSliceMatrix : public CMCPMatExpr<DoubleSliceMatrix<T> >
  {
  protected:
    /// the height
    int h;
    /// the width
    int w;
    /// the distance
    int distr, distc;
    /// the data
    T * data;
  public:

    /// element type
    typedef T TELEM;
    /// scalar type of elements (double or Complex)
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { IS_LINEAR = 0 };

    /// set height, width, and mem
    DoubleSliceMatrix (int ah, int aw, int adistr, int adistc, T * adata) throw ()
      : h(ah), w(aw), distr(adistr), distc(adistc), data(adata) { ; }
  
    /// assign contents
    template<typename TB>
    const DoubleSliceMatrix & operator= (const Expr<TB> & m) const
    {
      return CMCPMatExpr<DoubleSliceMatrix>::operator= (m);
    }

    /// assign constant
    DoubleSliceMatrix & operator= (TSCAL s) throw()
    {
      for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++)
          data[i*distr+j*distc] = s;
      return *this;
    }

    /// access operator
    TELEM & operator() (int i, int j) const
    {
#ifdef CHECK_RANGE
      CheckMatRange(h,w,i,j);
#endif
      return data[i*distr+j*distc]; 
    }

    /// access operator, linear access
    TELEM & operator() (int i) const
    {
#ifdef CHECK_RANGE
      CheckMatRange(h,w,i);
#endif
      return data[i]; 
    }

    int Height () const throw() { return h; }
    int Width () const throw() { return w; }
    int DistRow () const { return distr; }
    int DistCol () const { return distc; }
    T * Data () const { return data; }

    /*
    const DoubleSliceMatrix Rows (int first, int next) const
    {
      return DoubleSliceMatrix (next-first, w, distr, distc, data+first*distr);
    }

    const DoubleSliceMatrix<T> Cols (int first, int next) const
    {
      return DoubleSliceMatrix<T> (h, next-first, distr, distc, data+first*distc);
    }

    const DoubleSliceMatrix Rows (IntRange range) const
    {
      return Rows(range.First(), range.Next());
    }

    const DoubleSliceMatrix<T> Cols (IntRange range) const
    {
      return Cols(range.First(), range.Next());
    }
    */

  };
















  template <class TM, class TSCAL>
  class Scalar2ElemMatrix 
  {
  public:
    const FlatMatrix<TSCAL> & mat;
    Scalar2ElemMatrix (const FlatMatrix<TSCAL> & amat) : mat(amat) { ; }

    enum { H = mat_traits<TM>::HEIGHT };
    enum { W = mat_traits<TM>::WIDTH };

    TM operator() (int i, int j) const
    {
      TM ret;
      for (int k = 0; k < H; k++)
	for (int l = 0; l < W; l++)
	  Access(ret, k,l) = mat(i*H+k, j*W+l);
      return ret;
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




  /// Variable size identity matrix
  class Identity : public MatExpr<Identity >
  {
    int size;
  public:
    typedef double TELEM;
    typedef double TSCAL;
    enum { IS_LINEAR = 0 };

    Identity (int s) : size(s) { ; }

    double operator() (int i) const
    { cerr << "Identity, linear access" << endl; return 0; }

    double operator() (int i, int j) const { return (i == j) ? 1 : 0; }

    int Height () const { return size; }
    int Width () const { return size; }
  };





  template<int H, int W, typename T>
  inline std::ostream & operator<< (std::ostream & s, const Mat<H,W,T> & m)
  {
    for (int i = 0; i < H*W; i++)
      s << " " << setw(7) << m(i);
    return s;
  }



  template <int H, int W, typename T>
  inline auto Trans (const Mat<H,W,T> & mat)
    -> Mat<W,H,decltype(Trans(mat(0,0)))>
  {
    Mat<W,H,decltype(Trans(mat(0,0)))> res;
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        res(j,i) = mat(i,j);
    return res;
  }

  template <int H, int W, int W2, typename T1, typename T2>
  inline auto operator* (const Mat<H,W,T1> & mat1, const Mat<W,W2,T2> & mat2) 
    -> Mat<H,W2,decltype(mat1(0,0)*mat2(0,0))>
  {
    typedef decltype(mat1(0,0)*mat2(0)) TRES;
    Mat<H,W2,TRES> res;
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W2; j++)
        {
          TRES sum(0);
          for (int k = 0; k < W; k++)
            sum += mat1(i,k) * mat2(k,j);
          res(i,j) = sum;
        }
    return res;
  }


  template <int H, int W, int W2, typename T1, typename T2>
  inline auto operator* (const Mat<H,W,T1> & mat, const Vec<W2,T2> & vec) 
    -> Vec<H, decltype(mat(0,0)*vec(0))>
  {
    typedef decltype(mat(0,0)*vec(0)) TRES;
    Vec<H, TRES> res = TRES(0);
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        res(i) += mat(i,j) * vec(j);
    return res;
  }

  template <int H, int W, typename T1, typename T2>
  inline auto operator* (const Mat<H,W,T1> & mat, const FlatVec<W,T2> & vec) 
    -> Vec<H, decltype(mat(0,0)*vec(0))>
  {
    typedef decltype(mat(0,0)*vec(0)) TRES;
    Vec<H, TRES> res = TRES(0);
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        res(i) += mat(i,j) * vec(j);
    return res;
  }


}


#ifdef PARALLEL
namespace ngstd
{
  template<int N, int M, typename T>
  class MPI_Traits<ngbla::Mat<N, M, T> >
  {
  public:
    static MPI_Datatype MPIType () 
    { 
      static MPI_Datatype MPI_T = 0;
      if (!MPI_T)
	{
	  int size = N * M;
	  MPI_Type_contiguous ( size, MPI_Traits<T>::MPIType(), &MPI_T);
	  MPI_Type_commit ( &MPI_T );
	}
      return MPI_T;
    }
  };
  

}
#endif

#endif
