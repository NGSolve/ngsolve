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
  template <typename T = double, ORDERING ORD = RowMajor> class SliceMatrix;
  // template <typename T> class SliceMatrixColMajor;
  template <typename T> class DoubleSliceMatrix;

  ///
  extern void CheckMatRange(int h, int w, int i);
  ///
  extern void CheckMatRange(int h, int w, int i, int j);



  /**
     A simple matrix.
     Has height, width and data-pointer. 
     No memory allocation/deallocation. User must provide memory.
  */
  template <typename T, ORDERING ORD>
  class FlatMatrix : public CMCPMatExpr<FlatMatrix<T,ORD> >
                     // class FlatMatrix : public CMCPMatExpr<FlatMatrix<T> >
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
    // typedef T& TREF;
    /// scalar type of elements (double or Complex)
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /// nothing done in default constructor
    INLINE FlatMatrix () { ; }
  
    /// set height, width, and mem
    INLINE FlatMatrix (int ah, int aw, T * adata) 
      : h(ah), w(aw), data(adata) { ; }
  
    /// set height = width, and mem
    INLINE FlatMatrix (int ah, T * adata) 
      : h(ah), w(ah), data(adata) { ; }

    /// allocates at local heap
    INLINE FlatMatrix (int ah, int aw, LocalHeap & lh) 
      : h(ah), w(aw), data (lh.Alloc<T>(ah*aw)) { ; }
  
    /// allocates at local heap
    INLINE FlatMatrix (int ah, LocalHeap & lh) 
      : h(ah), w(ah), data(lh.Alloc<T>(ah*ah)) { ; }
  
    /// copy constructor. copies pointers, not contents
    INLINE FlatMatrix (const FlatMatrix & m) throw () 
      : h(m.h), w(m.w) , data(m.data)
    { ; }
  
    /// allocate and compute 
    template<typename TB>
    INLINE FlatMatrix (const LocalHeapExpr<TB> & m2) 
    {
      h = m2.A().Height();
      w = m2.A().Width();
      LocalHeap & lh = m2.GetLocalHeap();
      data = lh.Alloc<T> (h*w);
      CMCPMatExpr<FlatMatrix<T> >::operator= (m2.A());
    }

    /// useful to put FlatMatrix over other matrix
    template <typename T2>
    INLINE explicit FlatMatrix (const MatExpr<T2> & m)
      : h(m.Height()), w(m.Width()),
        data(const_cast<T*>(&m.Spec()(0,0)))  
    { ; }
  
    /// useful to put FlatMatrix over other Mat
    template <int H, int W>
    INLINE FlatMatrix (const Mat<H,W,TSCAL> & m) throw()
      : h(H), w(W), data(const_cast<T*>(&m(0,0)))
    { ; }

    // do nothing 
    // ~FlatMatrix () throw() { ; }

    /// set size, and assign mem
    INLINE void AssignMemory (int ah, int aw, LocalHeap & lh)  
    {
      h = ah;
      w = aw;
      data = lh.Alloc<T>(h*w);
    }
  
    /// set size, and assign mem
    INLINE void AssignMemory (int ah, int aw, T * mem) throw()
    {
      h = ah;
      w = aw;
      data = mem;
    }
  

    /// assign contents
    template<typename TBxx>
    INLINE const FlatMatrix & operator= (const Expr<TBxx> & m) const
    {
      return CMCPMatExpr<FlatMatrix>::operator= (m);
    }

    /// copy contents
    INLINE const FlatMatrix & operator= (const FlatMatrix & m) const 
    {
      for (int i = 0; i < h*w; i++) data[i] = m(i);
      return *this;
    }

    /// assign constant
    INLINE const FlatMatrix & operator= (TSCAL s) const 
    {
      for (int i = 0; i < h*w; i++) data[i] = s; 
      return *this;
    }

    /// copy size and pointers
    INLINE FlatMatrix & Assign (const FlatMatrix & m) throw()
    {
      h = m.h;
      w = m.w;
      data = m.data;
      return *this;
    }


    /// access operator, linear access
    INLINE TELEM & operator() (int i) const 
    { 
#ifdef CHECK_RANGE
      CheckMatRange(h,w,i);
#endif
      return data[i]; 
    }

    /// access operator
    INLINE TELEM & operator() (int i, int j) const
    {
#ifdef CHECK_RANGE
      CheckMatRange(h,w,i,j);
#endif
      return data[i*size_t(w)+j]; 
    }

    /// the height
    INLINE int Height () const { return h; }

    /// the width
    INLINE int Width () const { return w; }

    INLINE const FlatVector<T> Row (int i) const
    {
      return FlatVector<T> (w, &data[i*size_t(w)]);
    }

#ifdef FLATVECTOR_WITH_DIST
    INLINE const FlatVector<T> Col (int i) const
    {
      return FlatVector<T> (h, w, &data[i]);
    }
#else
    INLINE const SliceVector<T> Col (int i) const
    {
      return SliceVector<T> (h, w, &data[i]);
    }
#endif

    INLINE const SliceVector<T> Diag () const
    {
      return SliceVector<T> (h, w+1, &data[0]);
    }

    using CMCPMatExpr<FlatMatrix<T> >::Rows;
    using CMCPMatExpr<FlatMatrix<T> >::Cols;

    INLINE const FlatMatrix Rows (int first, int next) const
    {
      return FlatMatrix (next-first, w, data+first*w);
    }

    INLINE const SliceMatrix<T> Cols (int first, int next) const
    {
      return SliceMatrix<T> (h, next-first, w, data+first);
    }

    INLINE const FlatMatrix Rows (IntRange range) const
    {
      return FlatMatrix (range.Next()-range.First(), w, data+range.First()*w);
    }

    INLINE const SliceMatrix<T> Cols (IntRange range) const
    {
      return SliceMatrix<T> (h, range.Next()-range.First(), w, data+range.First());
    }

    INLINE operator SliceMatrix<T> () const
    {
      return SliceMatrix<T> (h, w, w, data);
    }
  };





  template <typename T>
  class FlatMatrix<T,ColMajor> : public CMCPMatExpr<FlatMatrix<T,ColMajor> >
  {
  protected:
    int h;
    int w;
    T * __restrict data;
  public:
    enum { IS_LINEAR = 0 };

    /// element type
    typedef T TELEM;
    /// scalar type of elements (double or Complex)
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /// nothing done in default constructor
    INLINE FlatMatrix () { ; }
  
    /// set height, width, and mem
    INLINE FlatMatrix (int ah, int aw, T * adata) 
      : h(ah), w(aw), data(adata) { ; }
  
    /// set height = width, and mem
    INLINE FlatMatrix (int ah, T * adata) 
      : h(ah), w(ah), data(adata) { ; }

    /// allocates at local heap
    INLINE FlatMatrix (int ah, int aw, LocalHeap & lh) 
      : h(ah), w(aw), data (lh.Alloc<T>(ah*aw)) { ; }
  
    /// allocates at local heap
    INLINE FlatMatrix (int ah, LocalHeap & lh) 
      : h(ah), w(ah), data(lh.Alloc<T>(ah*ah)) { ; }
  
    /// copy constructor. copies pointers, not contents
    INLINE FlatMatrix (const FlatMatrix & m) throw () 
      : h(m.h), w(m.w) , data(m.data)
    { ; }
  
    /// allocate and compute 
    template<typename TB>
    INLINE FlatMatrix (const LocalHeapExpr<TB> & m2) 
    {
      h = m2.A().Height();
      w = m2.A().Width();
      LocalHeap & lh = m2.GetLocalHeap();
      data = lh.Alloc<T> (h*w);
      CMCPMatExpr<FlatMatrix<T> >::operator= (m2.A());
    }

    /// set size, and assign mem
    INLINE void AssignMemory (int ah, int aw, LocalHeap & lh)  
    {
      h = ah;
      w = aw;
      data = lh.Alloc<T>(h*w);
    }
  
    /// set size, and assign mem
    INLINE void AssignMemory (int ah, int aw, T * mem) throw()
    {
      h = ah;
      w = aw;
      data = mem;
    }

    /// assign contents
    template<typename TBxx>
    INLINE const FlatMatrix & operator= (const Expr<TBxx> & m) const
    {
      return CMCPMatExpr<FlatMatrix>::operator= (m);
    }

    /// copy contents
    INLINE const FlatMatrix & operator= (const FlatMatrix & m) const 
    {
      for (int i = 0; i < h*w; i++) data[i] = m(i);
      return *this;
    }

    /// assign constant
    INLINE const FlatMatrix & operator= (TSCAL s) const 
    {
      for (int i = 0; i < h*w; i++) data[i] = s; 
      return *this;
    }

    /// copy size and pointers
    INLINE FlatMatrix & Assign (const FlatMatrix & m) throw()
    {
      h = m.h;
      w = m.w;
      data = m.data;
      return *this;
    }

    /// access operator, linear access
    INLINE TELEM & operator() (int i) const 
    { 
#ifdef CHECK_RANGE
      CheckMatRange(h,w,i);
#endif
      return data[i]; 
    }

    /// access operator
    INLINE TELEM & operator() (int i, int j) const
    {
#ifdef CHECK_RANGE
      CheckMatRange(h,w,i,j);
#endif
      return data[j*size_t(h)+i]; 
    }

    /// the height
    INLINE int Height () const { return h; }

    /// the width
    INLINE int Width () const { return w; }


    INLINE const FlatVector<T> Col (int i) const
    {
      return FlatVector<T> (h, &data[i*size_t(h)]);
    }

    INLINE const SliceVector<T> Row (int i) const
    {
      return SliceVector<T> (w, h, &data[i]);
    }

    /*
    const SliceVector<T> Diag () const
    {
      return SliceVector<T> (h, w+1, &data[0]);
    }
    */

    INLINE const FlatMatrix Cols (int first, int next) const
    {
      return FlatMatrix (h, next-first, data+first*h);
    }

    INLINE const FlatMatrix Cols (IntRange range) const
    {
      return FlatMatrix (h, range.Size(), data+range.First()*h);
    }
    
    /*
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
    */

    INLINE operator SliceMatrix<T,ColMajor> () const
    {
      return SliceMatrix<T,ColMajor> (h, w, h, data);
    }
  };












  /**
     A Matrix class with memory allocation/deallocation
  */
  template <typename T, ORDERING ORD>
  class Matrix : public FlatMatrix<T,ORD>
  {
  public:

    /// element type
    typedef T TELEM;
    /// scalar type of elements (double or Complex)
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /// default constructor
    Matrix () : FlatMatrix<T> (0, 0) { ; }

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

    Matrix (initializer_list<initializer_list<T>> llist) 
      : FlatMatrix<T> (0,0,nullptr)
    {
      int h = llist.size();
      int w = 0;
      for (auto row : llist)
        w = max(w, int(row.size()));

      SetSize (h, w);
      (*this) = 0.0;

      int r = 0;
      for (auto row : llist)
        {
          int c = 0;
          for (auto col : row)
            (*this)(r,c++) = col;
          r++;
        }
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
    // enum { HEIGHT = H };
    // enum { WIDTH  = W };

    /// do not initialize 
    INLINE Mat () { ; }

    /// copy matrix
    INLINE Mat (const Mat & m) 
      : MatExpr<Mat> ()
    {
      (*this) = m;
    }

    /// assign values
    template<typename TB>
    INLINE Mat (const Expr<TB> & m)
    {
      MatExpr<Mat>::operator= (m);
    }

    /// fill with scalar
    INLINE Mat (TSCAL s) 
    {
      for (int i = 0; i < H*W; i++) 
        data[i] = s;
    }

    /// assign values
    template<typename TB>
    INLINE Mat & operator= (const Expr<TB> & m)
    {
      MatExpr<Mat>::operator= (m);
      return *this;
    }

    /// copy matrix
    INLINE Mat & operator= (const Mat & m) 
    {
      for (int i = 0; i < H*W; i++)
        data[i] = m.data[i];
      return *this;
    }
 
    /// fill values
    INLINE Mat & operator= (TSCAL s) 
    {
      for (int i = 0; i < H*W; i++)
        data[i] = s;
      return *this;
    }

    /// linear access
    INLINE TELEM & operator() (int i) { return data[i]; }
    /// access element
    INLINE TELEM & operator() (int i, int j) { return data[i*W+j]; }
    /// linear access
    INLINE const TELEM & operator() (int i) const { return data[i]; }
    /// access element
    INLINE const TELEM & operator() (int i, int j) const { return data[i*W+j]; }

    /// the height
    INLINE int Height () const { return H; }
    /// the width
    INLINE int Width () const { return W; }

    ///
    INLINE const FlatVec<W,T> Row (int i) 
    {
      return FlatVec<W,T> (&(*this)(i,0));
    }

    ///
    INLINE const FlatVec<W,const T> Row (int i) const
    {
      return FlatVec<W,const T> (&(*this)(i,0));
    }

    ///
    INLINE const FixSliceVector<W,T> Col (int i) 
    {
      return FixSliceVector<W,T> (H, &(*this)(0,i));
    }

    ///
    INLINE const FixSliceVector<W,const T> Col (int i) const
    {
      return FixSliceVector<W,const T> (H, &(*this)(0,i));
    }

  };

  template <int H, int W, typename T>
  class mat_traits<Mat<H,W,T>>
  {
  public:
    /// matrix element
    typedef T TELEM;
    /// field of matrix element
    typedef typename mat_traits<T>::TSCAL TSCAL;
    /// type of column vector
    typedef typename Mat<H,W,T>::TV_COL TV_COL;
    /// type of row vector
    typedef typename Mat<H,W,T>::TV_ROW TV_ROW;
    /// matrix height
    enum { HEIGHT = H };
    /// matrix with
    enum { WIDTH  = W  };
    ///
    enum { IS_COMPLEX = mat_traits<TSCAL>::IS_COMPLEX };
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
  template <int W, typename T = double, int DIST = W>
  class FlatMatrixFixWidth : 
    public CMCPMatExpr<FlatMatrixFixWidth<W,T,DIST>>
  {
  protected:
    /// the data
    T *  __restrict data;
    /// the height
    int h;
  public:
    /// entry type
    typedef T TELEM;
    /// scalar type of entry
    typedef typename mat_traits<T>::TSCAL TSCAL;

    enum { IS_LINEAR = (W == DIST) }; 

    /// nothing done in default constructor
    INLINE FlatMatrixFixWidth () { ; }
  
    /// set height and mem
    INLINE FlatMatrixFixWidth (int ah, T * adata) throw()
      : data(adata), h(ah) { ; }

    /// allocates at local heap
    INLINE FlatMatrixFixWidth (int ah, LocalHeap & lh) 
      : data(lh.Alloc<T>(ah*DIST)), h(ah) { ; }
  
    /// copy constructor. copies pointers, not contents
    INLINE FlatMatrixFixWidth (const FlatMatrixFixWidth & m)
      : data(m.data), h(m.h)
    { ; }

    /// copy constructor. copies pointers, not contents
    INLINE FlatMatrixFixWidth (const FlatMatrix<TSCAL> & m) 
      : data(const_cast<T*> (&m(0))), h(m.Height())
    { ; }

    template <int H>
    INLINE FlatMatrixFixWidth (const Mat<H,W,TSCAL> & m) 
      : data(const_cast<T*>(&m(0,0))), h(H)
    { ; }  

    /// set size, and assign mem
    INLINE void AssignMemory (int ah, LocalHeap & lh) 
    {
      h = ah;
      data = lh.Alloc<T> (h*DIST);
    }
  
    /// set size, and assign mem
    INLINE void AssignMemory (int ah, T * mem) throw()
    {
      h = ah;
      data = mem;
    }
  

    /// assign contents
    template<typename TB>
    INLINE const FlatMatrixFixWidth & operator= (const Expr<TB> & m) const
    {
      return CMCPMatExpr<FlatMatrixFixWidth>::operator= (m);
    }

    /// copy contents
    INLINE const FlatMatrixFixWidth & operator= (const FlatMatrixFixWidth & m) const throw()
    {
      if (W == DIST)
        for (int i = 0; i < h*W; i++)
          data[i] = m(i);
      else
        for (int i = 0; i < h; i++)
          for (int j = 0; j < W; j++)
            (*this)(i,j) = m(i,j);
        
      return *this;
    }

    /// assign constant
    INLINE const FlatMatrixFixWidth & operator= (TSCAL s) const throw()
    {
      if (W == DIST)
        for (int i = 0; i < h*W; i++)
          data[i] = s; 
      else
        for (int i = 0; i < h; i++)
          for (int j = 0; j < W; j++)
            (*this)(i,j) = s;
        
      return *this;
    }

    /// copy size and pointers
    INLINE const FlatMatrixFixWidth & Assign (const FlatMatrixFixWidth & m) throw()
    {
      h = m.h;
      data = m.data;
      return *this;
    }

    /// access operator, linear access
    INLINE TELEM & operator() (int i) const
    {
#ifdef CHECK_RANGE
      CheckMatRange(h,W,i);
#endif
      return data[i]; 
    }

    /// access operator
    INLINE TELEM & operator() (int i, int j) const
    {
#ifdef CHECK_RANGE
      CheckMatRange(h,W,i,j);
#endif
      return data[i*DIST+j]; 
    }

    /// the height
    INLINE int Height () const throw() { return h; }

    /// the width
    INLINE int Width () const throw() { return W; }

    ///
    INLINE operator const FlatMatrix<T>() const { return FlatMatrix<T> (h, W, data); }

    ///
    INLINE const FlatVec<W,T> Row (int i) const
    {
      return FlatVec<W,T> (&(*this)(i,0));
    }

    INLINE const FixSliceVector<DIST,T> Col (int i) const
    {
      return FixSliceVector<DIST,T> (h, &data[i]);
    }

    using CMCPMatExpr<FlatMatrixFixWidth<W,T,DIST> >::Rows;
    using CMCPMatExpr<FlatMatrixFixWidth<W,T,DIST> >::Cols;

    INLINE const FlatMatrixFixWidth Rows (int first, int next) const
    {
      return FlatMatrixFixWidth (next-first, data+first*DIST);
    }

    INLINE const SliceMatrix<T> Cols (int first, int next) const
    {
      return SliceMatrix<T> (h, next-first, DIST, data+first);
    }

    INLINE const FlatMatrixFixWidth Rows (IntRange range) const
    {
      return Rows (range.First(), range.Next());
    }

    INLINE const SliceMatrix<T> Cols (IntRange range) const
    {
      return Cols (range.First(), range.Next());
    }

    INLINE operator SliceMatrix<T> () const
    {
      return SliceMatrix<T> (h, W, DIST, data);
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

    INLINE MatrixFixWidth () { ; }

    /// allocate matrix of size ah * ah
    INLINE MatrixFixWidth (int ah) : FlatMatrixFixWidth<W,T> (ah, new T[ah*W]) { ; }

    /// delete memory
    INLINE ~MatrixFixWidth() { delete [] this->data; }

    /// sets new size of matrix
    INLINE void SetSize(int ah)
    {
      if (this->h == ah) return;
      delete [] this->data;
      this->h = ah;
      this->data = new T[this->h * W];
    }

    /// assign matrix, sizes must match
    template<typename TB>
    INLINE MatrixFixWidth & operator= (const Expr<TB> & m)
    {
      MatExpr<FlatMatrixFixWidth<W,T> >::operator= (m);
      return *this;
    }

    /// fill matrix with scalar
    INLINE MatrixFixWidth & operator= (TSCAL s)
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
  template <int H, typename T = double, int SLICE = H>
  class FlatMatrixFixHeight : public MatExpr<FlatMatrixFixHeight<H,T,SLICE> >
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
    FlatMatrixFixHeight () : data(0), w(0) { ; }
  
    /// set height and mem
    FlatMatrixFixHeight (int aw, T * adata) 
      : data(adata), w(aw) { ; }

    /// allocates at local heap
    FlatMatrixFixHeight (int aw, LocalHeap & lh) 
      : data(lh.Alloc<T>(aw*SLICE)), w(aw) { ; }
  
    /// copy constructor. copies pointers, not contents
    FlatMatrixFixHeight (const FlatMatrixFixHeight & m)
      : data(m.data), w(m.w)
    { ; }
  

    // do nothing 
    // ~FlatMatrixFixHeight () throw() { ; }
    
    /// set size, and assign mem
    void AssignMemory (int aw, LocalHeap & lh) 
    {
      w = aw;
      data = (T*)lh.Alloc(w*SLICE*sizeof(T));
    }
  
    /// set size, and assign mem
    void AssignMemory (int aw, T * mem) 
    {
      w = aw;
      data = mem;
    }
  

    /// assign contents
    template<typename TB>
    const FlatMatrixFixHeight & operator= (const Expr<TB> & m) const
    {
      for (int j = 0; j < w; j++)
        for (int i = 0; i < H; i++)
          (*this)(i,j) = m.Spec()(i,j);
      return *this;
    }

    /// copy contents
    const FlatMatrixFixHeight & operator= (const FlatMatrixFixHeight & m) const
    {
      if (H == SLICE)
        for (int i = 0; i < w*H; i++)
          data[i] = m(i);
      else
        for (int j = 0; j < w; j++)
          for (int i = 0; i < H; i++)
            (*this)(i,j) = m(i,j);
      
      return *this;
    }

    /// assign constant
    FlatMatrixFixHeight & operator= (TSCAL s)
    {
      if (H == SLICE)
        for (int i = 0; i < w*H; i++)
          data[i] = s; 
      else
        for (int j = 0; j < w; j++)
          for (int i = 0; i < H; i++)
            (*this)(i,j) = s;
        
      return *this;
    }

    /// copy size and pointers
    FlatMatrixFixHeight & Assign (const FlatMatrixFixHeight & m)
    {
      w = m.w;
      data = m.data;
      return *this;
    }

    /*
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
      return data[i+j*SLICE]; 
    }
    */

    /// access operator, linear access
    TELEM & operator() (int i) const
    {
#ifdef CHECK_RANGE
      CheckMatRange(H,w,i);
#endif
      return data[i]; 
    }

    /// access operator
    TELEM & operator() (int i, int j) const
    {
#ifdef CHECK_RANGE
      CheckMatRange(H,w,i,j);
#endif
      return data[i+j*SLICE]; 
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
      return SliceVector<T> (w, SLICE, &data[i]);
    }


    const FlatMatrixFixHeight Cols (int first, int next) const
    {
      return FlatMatrixFixHeight (next-first, data+first*SLICE);
    }

    const FlatMatrixFixHeight Cols (IntRange range) const
    {
      return FlatMatrixFixHeight (range.Size(), data+range.First()*SLICE);
    }

    const SliceMatrix<T,ColMajor>
    Rows (int first, int next) const
    { 
      return SliceMatrix<T,ColMajor> (next-first, w, SLICE, data+first);
    }

    const SliceMatrix<T,ColMajor>    
    Rows (IntRange range) const // -> decltype (this->Rows (range.First(), range.Next()))
    { 
      return Rows (range.First(), range.Next());
    }

    operator FlatMatrix<T,ColMajor> () const
    {
      static_assert(H==SLICE, "MatrixFixHeight to FlatMatrix, but slice != height");
      return FlatMatrix<T,ColMajor> (H, w, data);
    }

    operator SliceMatrix<T,ColMajor> () const
    {
      return SliceMatrix<T,ColMajor> (H, w, SLICE, data);
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












  template <typename T, ORDERING ORD>
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

    // 
    SliceMatrix() = delete;
    INLINE SliceMatrix(const SliceMatrix &) = default;
    INLINE SliceMatrix & operator= (const SliceMatrix &) = delete;

    /// set height, width, and mem
    INLINE SliceMatrix (int ah, int aw, int adist, T * adata) throw ()
      : h(ah), w(aw), dist(adist), data(adata) { ; }

    /*
    SliceMatrix (const FlatMatrix<T> & mat)
      : h(mat.Height()), w(mat.Width()), dist(mat.Width()), data(&mat(0,0))
    { ; }

    template<int W>
    SliceMatrix (const FlatMatrixFixWidth<W,T> & mat)
      : h(mat.Height()), w(mat.Width()), dist(mat.Width()), data(&mat(0,0))
    { ; }
    */

    template<int H, int W>
    INLINE SliceMatrix (Mat<H,W,T> & mat)
      : h(mat.Height()), w(mat.Width()), dist(mat.Width()), data(&mat(0,0))
    { ; }

    /// assign contents
    template<typename TB>
    INLINE const SliceMatrix & operator= (const Expr<TB> & m) const
    {
      return CMCPMatExpr<SliceMatrix>::operator= (m);
    }

  
    INLINE const SliceMatrix & operator= (const SliceMatrix & m) const
    {
      return CMCPMatExpr<SliceMatrix>::operator= (m);
    }

    /*
    template<int W>
    const SliceMatrix & operator= (const FlatMatrixFixWidth<W,T> & m) const
    {
      return CMCPMatExpr<SliceMatrix>::operator= (m);
    }
    template<int W>
    const SliceMatrix & operator= (const MatrixFixWidth<W,T> & m) const
    {
      return CMCPMatExpr<SliceMatrix>::operator= (m);
    }
    template <typename TA, typename TB>
    const SliceMatrix & operator= (const LapackProduct<TA, TB> & prod) const
    {
      LapackMult (prod.A(), prod.B(), *this);
      return *this;
    }
    */

    /// assign constant
    INLINE const SliceMatrix & operator= (TSCAL s) const
    {
      for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++)
          data[i*dist+j] = s;
      return *this;
    }

    /// access operator
    INLINE TELEM & operator() (int i, int j) const
    {
#ifdef CHECK_RANGE
      CheckMatRange(h,w,i,j);
#endif
      return data[i*dist+j]; 
    }

    /// access operator, linear access
    INLINE TELEM & operator() (int i) const
    {
#ifdef CHECK_RANGE
      CheckMatRange(h,w,i);
#endif
      return data[i]; 
    }

    /// the height
    INLINE int Height () const throw() { return h; }

    /// the width
    INLINE int Width () const throw() { return w; }

    /// 
    INLINE int Dist () const throw() { return dist; }

    INLINE const SliceMatrix Rows (int first, int next) const
    {
      return SliceMatrix (next-first, w, dist, data+first*dist);
    }

    INLINE const FlatVector<T> Row (int i) const
    {
      return FlatVector<T> (w, &data[i*size_t(dist)]);
    }

    INLINE const SliceVector<T> Col (int i) const
    {
      return SliceVector<T> (h, dist, &data[i]);
    }

    INLINE const SliceMatrix<T> Cols (int first, int next) const
    {
      return SliceMatrix<T> (h, next-first, dist, data+first);
    }

    INLINE const SliceMatrix Rows (IntRange range) const
    {
      return Rows (range.First(), range.Next());
    }

    INLINE const SliceMatrix<T> Cols (IntRange range) const
    {
      return Cols (range.First(), range.Next());
    }

  };




  template <typename T>
  class SliceMatrix<T,ColMajor> : public CMCPMatExpr<SliceMatrix<T,ColMajor> >
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
    enum { COL_MAJOR = 1 }; 


    /// set height, width, and mem
    SliceMatrix (int ah, int aw, int adist, T * adata) throw ()
      : h(ah), w(aw), dist(adist), data(adata) { ; }

    /// assign contents
    template<typename TB>
    INLINE const SliceMatrix & operator= (const Expr<TB> & m) const
    {
      return CMCPMatExpr<SliceMatrix>::operator= (m);
    }

    /// assign constant
    const SliceMatrix & operator= (TSCAL s) const throw()
    {
      for (int i = 0; i < w; i++)
        for (int j = 0; j < h; j++)
          data[i*dist+j] = s;
      return *this;
    }

    /// access operator
    TELEM & operator() (int i, int j) const
    {
      return data[j*dist+i]; 
    }

    /// access operator, linear access
    TELEM & operator() (int i) const
    {
      return data[i]; 
    }

    /// the height
    int Height () const throw() { return h; }

    /// the width
    int Width () const throw() { return w; }

    /// 
    int Dist () const throw() { return dist; }


    const FlatVector<T> Col (int i) const
    {
      return FlatVector<T> (h, &data[i*size_t(dist)]);
    }

    const SliceVector<T> Row (int i) const
    {
      return SliceVector<T> (w, dist, &data[i]);
    }


    /*
    const SliceMatrix Rows (int first, int next) const
    {
      return SliceMatrix (next-first, w, dist, data+first*dist);
    }
    const SliceMatrix Rows (IntRange range) const
    {
      return Rows (range.First(), range.Next());
    }
    */

    const SliceMatrix Cols (int first, int next) const
    {
      return SliceMatrix (h, next-first, dist, data+first*dist);
    }


    const SliceMatrix Cols (IntRange range) const
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




  INLINE
  FlatMatrix<double,ColMajor> Trans (FlatMatrix<double,RowMajor> mat)
  {
    return FlatMatrix<double,ColMajor> (mat.Width(), mat.Height(), &mat(0,0));
  }

  INLINE
  FlatMatrix<Complex,ColMajor> Trans (FlatMatrix<Complex,RowMajor> mat)
  {
    return FlatMatrix<Complex,ColMajor> (mat.Width(), mat.Height(), &mat(0,0));
  }

  INLINE
  FlatMatrix<double,RowMajor> Trans (FlatMatrix<double,ColMajor> mat)
  {
    return FlatMatrix<double,RowMajor> (mat.Width(), mat.Height(), &mat(0,0));
  }

  INLINE
  FlatMatrix<Complex,RowMajor> Trans (FlatMatrix<Complex,ColMajor> mat)
  {
    return FlatMatrix<Complex,RowMajor> (mat.Width(), mat.Height(), &mat(0,0));
  }




  INLINE 
  const SliceMatrix<double> Trans (SliceMatrix<double,ColMajor> mat)
  {
    return SliceMatrix<double> (mat.Width(), mat.Height(), mat.Dist(), &mat(0,0));
  }

  INLINE 
  const SliceMatrix<Complex> Trans (SliceMatrix<Complex,ColMajor> mat)
  {
    return SliceMatrix<Complex> (mat.Width(), mat.Height(), mat.Dist(), &mat(0,0));
  }

  template <int H, int DIST>
  INLINE const FlatMatrixFixWidth<H,double,DIST> Trans (FlatMatrixFixHeight<H,double,DIST> mat)
  {
    return FlatMatrixFixWidth<H,double,DIST> (mat.Width(), &mat(0,0));
  }

  template <int H, int DIST>
  INLINE const FlatMatrixFixWidth<H,Complex,DIST> Trans (FlatMatrixFixHeight<H,Complex,DIST> mat)
  {
    return FlatMatrixFixWidth<H,Complex,DIST> (mat.Width(), &mat(0,0));
  }






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
    enum { IS_LINEAR = 0 };

    /// nothing to do 
    Id () { ; }

    ///
    double operator() (int i) const
    {
      static_assert (true, "linear access of id");
      cerr << "id, linear access" << endl; return 0;
    }
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
  INLINE auto Trans (const Mat<H,W,T> & mat)
    -> Mat<W,H,decltype(Trans(mat(0,0)))>
  {
    Mat<W,H,decltype(Trans(mat(0,0)))> res;
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        res(j,i) = mat(i,j);
    return res;
  }
  



  template <int H, int W, int W2, typename T1, typename T2>
  INLINE auto operator* (const Mat<H,W,T1> & mat1, const Mat<W,W2,T2> & mat2) 
    -> Mat<H,W2,decltype( RemoveConst(mat1(0,0)*mat2(0,0)))>
  {
    typedef decltype( RemoveConst(mat1(0,0)*mat2(0))) TRES;
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
  INLINE auto operator* (const Mat<H,W,T1> & mat, const Vec<W2,T2> & vec) 
    -> Vec<H, decltype(RemoveConst(mat(0,0)*vec(0)))>
  {
    typedef decltype(RemoveConst(mat(0,0)*vec(0))) TRES;
    Vec<H, TRES> res = TRES(0);
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        res(i) += mat(i,j) * vec(j);
    return res;
  }

  template <int H, int W, typename T1, typename T2>
  INLINE auto operator* (const Mat<H,W,T1> & mat, const FlatVec<W,T2> & vec) 
    -> Vec<H, decltype(RemoveConst(mat(0,0)*vec(0)))>
  {
    typedef decltype(RemoveConst(mat(0,0)*vec(0))) TRES;
    Vec<H, TRES> res = TRES(0);
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        res(i) += mat(i,j) * vec(j);
    return res;
  }



  //
  //  Can we put a SliceMatrix over a matrix ? 
  //

  template <typename TMAT, typename T = double>
  class Is_Sliceable { public: enum { VAL = false };  };

  template <typename T>
  class Is_Sliceable<SliceMatrix<T>,T> { public: enum { VAL = true };  };
  
  template <typename T>
  class Is_Sliceable<FlatMatrix<T>,T> { public: enum { VAL = true };  };
  
  template <typename T>
  class Is_Sliceable<Matrix<T>,T> { public: enum { VAL = true };  };
  
  template <int W, typename T, int DIST>
  class Is_Sliceable<FlatMatrixFixWidth<W,T,DIST>,T> { public: enum { VAL = true };  };


  template <typename TMAT, typename T>
  class Is_Sliceable<const TMAT, T> { public: enum { VAL = Is_Sliceable<TMAT,T>::VAL };  };
  template <typename TMAT, typename T>
  class Is_Sliceable<TMAT&,T> { public: enum { VAL = Is_Sliceable<TMAT,T>::VAL };  };
  template <typename TMAT, typename T>
  class Is_Sliceable<TMAT&&,T> { public: enum { VAL = Is_Sliceable<TMAT,T>::VAL };  };



  template <int IsIt, typename TMAT>
  class slicetype
  {
  public:
    typedef const TMAT & TYPE;
  };
  
  template <typename TMAT>
  class slicetype<1,TMAT>
  {
  public:
    typedef SliceMatrix<typename TMAT::TELEM> TYPE;
  };


  //
  // if possible, a slicematrix is returned
  // if not, we return a reference to the orginal matrix
  // 
  
  template <typename T, typename TMAT>
  auto SliceIfPossible (const TMAT & mat) -> typename slicetype<Is_Sliceable<TMAT,T>::VAL,TMAT>::TYPE
  {
    // cout << "type(tmat) = " << typeid(TMAT).name() << endl;
    // cout << "sliceable = " << Is_Sliceable<TMAT,T>::VAL << endl;
    // cout << "return type = " << typeid (typename slicetype<Is_Sliceable<TMAT,T>::VAL,TMAT>::TYPE).name() << endl;
    return mat;
  }

  



}


#ifdef PARALLEL
namespace ngstd
{
  template<int N, int M, typename T>
  class MPI_Traits<ngbla::Mat<N, M, T> >
  {
  public:
    /// gets the MPI datatype
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
