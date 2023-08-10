#ifndef FILE_MATRIX_EXPR
#define FILE_MATRIX_EXPR

/**************************************************************************/
/* File:   matrix.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jan. 02                                                    */
/**************************************************************************/


#include "expr.hpp"
#include "vector.hpp"

namespace ngbla
{

  

  template <typename TS, typename TD, typename T>
  auto make_VectorView (TS size, TD dist, T * data)
  {
    if constexpr (is_same_v<TS, undefined_size> == true)
      {
        if constexpr (std::is_same_v<TD, IC<1>> == true)
          return BareVector(data, size);
        else
          return BareSliceVector(data, dist, size);
      }
    else
      {
        if constexpr (std::is_same_v<TD, IC<1>> == true)
          return FlatVector(size, data);
        else
          return SliceVector(size, dist, data);
      }
  }

  
  template <int H, int W, typename T> class Mat;
  // template <typename T = double, ORDERING ORD = RowMajor> class SliceMatrix;
  template <typename T = double, ORDERING ORD = RowMajor> class BareSliceMatrix;
  // template <typename T> class SliceMatrixColMajor;
  template <typename T> class DoubleSliceMatrix;

  /**
     A simple matrix.
     Has height, width and data-pointer. 
     No memory allocation/deallocation. User must provide memory.
  */
  template <typename T, ORDERING ORD>
  class FlatMatrix : public CMCPMatExpr<FlatMatrix<T,ORD> >
  {
  protected:
    /// the height
    size_t h;
    /// the width
    size_t w;
    /// the data
    T * __restrict data;
  public:

    /// element type
    typedef T TELEM;
    /// scalar type of elements (double or Complex)
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /// nothing done in default constructor
    INLINE FlatMatrix () = default; 
  
    /// set height, width, and mem
    INLINE FlatMatrix (size_t ah, size_t aw, T * adata) 
      : h(ah), w(aw), data(adata) { ; }
  
    /// set height = width, and mem
    INLINE FlatMatrix (size_t ah, T * adata) 
      : h(ah), w(ah), data(adata) { ; }

    /// allocates at local heap
    INLINE FlatMatrix (size_t ah, size_t aw, LocalHeap & lh) 
      : h(ah), w(aw), data (lh.Alloc<T>(ah*aw)) { ; }
  
    /// allocates at local heap
    INLINE FlatMatrix (size_t ah, LocalHeap & lh) 
      : h(ah), w(ah), data(lh.Alloc<T>(ah*ah)) { ; }
  
    /// copy constructor. copies pointers, not contents
    FlatMatrix (const FlatMatrix & m) = default;
    /// allocate and compute 
    template<typename TB>
    INLINE FlatMatrix (const LocalHeapExpr<TB> & m2) 
    {
      h = m2.A().Height();
      w = m2.A().Width();
      LocalHeap & lh = m2.GetLocalHeap();
      data = lh.Alloc<T> (h*w);
      CMCPMatExpr<FlatMatrix<T,ORD> >::operator= (m2.A());
    }

    /// useful to put FlatMatrix over other matrix
    template <typename T2>
    INLINE explicit FlatMatrix (const MatExpr<T2> & m)
      : h(m.Height()), w(m.Width()),
        data(const_cast<T*>(m.Spec().Data()))  
    {
      static_assert(ORD == m.Ordering());
    }
  
    /// useful to put FlatMatrix over other Mat
    template <int H, int W>
    INLINE FlatMatrix (const Mat<H,W,TSCAL> & m)
      : h(H), w(W), data(const_cast<T*>(m.Data()))
    { ; }

    // do nothing 
    // ~FlatMatrix () throw() { ; }

    /// set size, and assign mem
    INLINE void AssignMemory (size_t ah, size_t aw, LocalHeap & lh)  
    {
      h = ah;
      w = aw;
      data = lh.Alloc<T>(h*w);
    }
  
    /// set size, and assign mem
    INLINE void AssignMemory (size_t ah, size_t aw, T * mem) 
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
      // for (size_t i = 0; i < h*w; i++) data[i] = m(i);
      AsVector() = m.AsVector();
      return *this;
    }

    /// assign constant
    INLINE const FlatMatrix & operator= (TSCAL s) const 
    {
      // for (auto i : Range(h*w)) data[i] = s;
      AsVector() = s;
      return *this;
    }


    INLINE auto View() const { return FlatMatrix(*this); }     
    INLINE tuple<size_t, size_t> Shape() const { return { h, w }; }
    static constexpr auto Ordering() { return ORD; } 
    
    /// copy size and pointers
    INLINE FlatMatrix & Assign (const FlatMatrix & m) 
    {
      h = m.h;
      w = m.w;
      data = m.data;
      return *this;
    }


    /// access operator, linear access
    INLINE TELEM & operator() (size_t i) const 
    { 
      NETGEN_CHECK_RANGE(i, 0, Height()*Width());
      return data[i]; 
    }

    /// access operator
    INLINE TELEM & operator() (size_t i, size_t j) const
    {
      NETGEN_CHECK_RANGE(i, 0, Height());
      NETGEN_CHECK_RANGE(j, 0, Width());
      return data[i*w+j];
    }

    /// the height
    INLINE auto Height () const { return h; }
    /// the width
    INLINE auto Width () const { return w; }
    ///
    INLINE auto Dist () const { return w; }

    INLINE T * Data () const { return data; }
    
    INLINE const FlatVector<T> Row (size_t i) const
    {
      return FlatVector<T> (w, data+i*w);
    }

    INLINE const SliceVector<T> Col (size_t i) const
    {
      return SliceVector<T> (h, w, data+i);
    }

    INLINE const SliceVector<T> Diag () const
    {
      return SliceVector<T> (h, w+1, &data[0]);
    }

    const SliceVector<T> Diag (int offset) const
    {
      int dp = std::max(offset, 0);
      int dm = std::min(offset, 0);
      return SliceVector<T> (std::min(w-dp, h+dm), Dist()+1, data+dp-dm*Dist());
    }
    

    using CMCPMatExpr<FlatMatrix>::Rows;
    using CMCPMatExpr<FlatMatrix>::Cols;

    INLINE FlatMatrix Rows (size_t first, size_t next) const
    {
      return FlatMatrix (next-first, w, data+first*w);
    }

    INLINE SliceMatrix<T> Cols (size_t first, size_t next) const
    {
      return SliceMatrix<T> (h, next-first, w, data+first);
    }

    INLINE FlatMatrix Rows (IntRange range) const
    {
      return FlatMatrix (range.Next()-range.First(), w, data+range.First()*w);
    }

    INLINE SliceMatrix<T> Cols (IntRange range) const
    {
      return SliceMatrix<T> (h, range.Next()-range.First(), w, data+range.First());
    }

    BareSliceMatrix<T> RowSlice(size_t first, size_t adist) const
    {
      return BareSliceMatrix<T> (w*adist, data+first*w, DummySize( Height()/adist, w));
    }

    INLINE auto SplitRows (size_t split) const
    {
      return tuple(Rows(0,split), Rows(split, Height()));
    }

    INLINE auto SplitCols (size_t split) const
    {
      return tuple(Cols(0,split), Cols(split, Width()));
    }

    
    /*
    INLINE operator SliceMatrix<T> () const
    {
      return SliceMatrix<T> (h, w, w, data);
    }
    */
    auto AsVector() const
    {
      return FlatVector<T> (h*w, data);
    }

    auto Reshape(size_t h2, size_t w2)
    {
      return FlatMatrix{h2,w2,data};
    }
  };





  template <typename T>
  class FlatMatrix<T,ColMajor> : public CMCPMatExpr<FlatMatrix<T,ColMajor> >
  {
  protected:
    size_t h;
    size_t w;
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
    INLINE FlatMatrix (const FlatMatrix & m) 
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
      CMCPMatExpr<FlatMatrix<T,ColMajor> >::operator= (m2.A());
    }

    /// set size, and assign mem
    INLINE void AssignMemory (int ah, int aw, LocalHeap & lh)  
    {
      h = ah;
      w = aw;
      data = lh.Alloc<T>(h*w);
    }
  
    /// set size, and assign mem
    INLINE void AssignMemory (int ah, int aw, T * mem) 
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
      for (size_t i = 0; i < size_t(h)*size_t(w); i++) data[i] = m(i);
      return *this;
    }

    /// assign constant
    INLINE const FlatMatrix & operator= (TSCAL s) const 
    {
      for (size_t i = 0; i < size_t(h)*size_t(w); i++) data[i] = s; 
      return *this;
    }

    INLINE auto View() const { return FlatMatrix(*this); }     
    INLINE tuple<size_t, size_t> Shape() const { return { h, w }; }
    
    /// copy size and pointers
    INLINE FlatMatrix & Assign (const FlatMatrix & m) 
    {
      h = m.h;
      w = m.w;
      data = m.data;
      return *this;
    }

    /// access operator, linear access
    INLINE TELEM & operator() (size_t i) const 
    { 
      NETGEN_CHECK_RANGE(i, 0, Height()*Width());
      return data[i]; 
    }

    /// access operator
    INLINE TELEM & operator() (size_t i, size_t j) const
    {
      NETGEN_CHECK_RANGE(i, 0, Height());
      NETGEN_CHECK_RANGE(j, 0, Width());
      return data[j*size_t(h)+i]; 
    }

    /// the height
    INLINE size_t Height () const { return h; }

    /// the width
    INLINE size_t Width () const { return w; }

    INLINE size_t Dist () const { return h; }

    INLINE T * Data () const { return data; }
    
    INLINE const FlatVector<T> Col (size_t i) const
    {
      return FlatVector<T> (h, &data[i*size_t(h)]);
    }

    INLINE const SliceVector<T> Row (size_t i) const
    {
      return SliceVector<T> (w, h, &data[i]);
    }

    const SliceVector<T> Diag () const
    {
      return SliceVector<T> (h, w+1, &data[0]);
    }

    INLINE const FlatMatrix Cols (size_t first, size_t next) const
    {
      return FlatMatrix (h, next-first, data+first*h);
    }

    INLINE const FlatMatrix Cols (IntRange range) const
    {
      return FlatMatrix (h, range.Size(), data+range.First()*h);
    }


    INLINE const SliceMatrix<T,ColMajor> Rows (size_t first, size_t next) const
    {
      return SliceMatrix<T,ColMajor> (next-first, w, h, data+first);
    }

    INLINE const SliceMatrix<T,ColMajor> Rows (IntRange range) const
    {
      return SliceMatrix<T,ColMajor> (range.Size(), w, h, data+range.First());
    }

    INLINE auto SplitRows (size_t split) const
    {
      return tuple(Rows(0,split), Rows(split, Height()));
    }

    INLINE auto SplitCols (size_t split) const
    {
      return tuple(Cols(0,split), Cols(split, Width()));
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
    Matrix () : FlatMatrix<T,ORD> (0, 0) { ; }

    /// allocate matrix of size ah * ah
    Matrix (size_t ah) : FlatMatrix<T,ORD> (ah, new T[ah*ah]) { ; }
    
    /// allocate matrix of size ah * aw
    Matrix (size_t ah, size_t aw) : FlatMatrix<T,ORD> (ah, aw, new T[ah*aw]) { ; }

    /// allocate and copy matrix  
    INLINE Matrix (const Matrix & m2) 
      : FlatMatrix<T,ORD> (m2.Height(), m2.Width(), new T[m2.Height()*m2.Width()]) 
    {
      FlatMatrix<T,ORD>::operator= (m2);
    }
    
    /// move matrix
    INLINE Matrix (Matrix && m2)
      : FlatMatrix<T,ORD> (m2.h, m2.w, m2.data)
    { m2.data = nullptr; m2.w = 0; m2.h = 0; } 

    /// allocate and compute 
    template<typename TB>
    INLINE Matrix (const Expr<TB> & m2) 
      : FlatMatrix<T,ORD> (m2.Height(), m2.Width(), new T[m2.Height()*m2.Width()]) 
    {
      CMCPMatExpr<FlatMatrix<T,ORD> >::operator= (m2); // .View());   // does it work with view ? 
    }
    
    Matrix (initializer_list<initializer_list<T>> llist) 
      : FlatMatrix<T,ORD> (0,0,nullptr)
    {
      int h = llist.size();
      int w = 0;
      for (auto row : llist)
        w = std::max(w, int(row.size()));

      SetSize (h, w);
      (*this) = T(0.0);

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
    void SetSize(size_t ah, size_t aw)
    {
      if (this->h == ah && this->w == aw) return;
      delete [] this->data;
      this->h = ah;
      this->w = aw;
      this->data = new T[this->h * this->w];
    }

    /// sets new size of matrix
    void SetSize(size_t ah)
    {
      if (this->h == ah && this->w == ah) return;
      delete [] this->data;
      this->h = ah;
      this->w = ah;
      this->data = new T[this->h * this->w];
    }

    /// assign matrix, sizes must match
    template<typename TB>
    INLINE Matrix & operator= (const Expr<TB> & m) 
    { 
      SetSize (m.Height(), m.Width());
      CMCPMatExpr<FlatMatrix<T,ORD> >::operator= (m);
      return *this;
    }


    /// fill matrix with scalar
    Matrix & operator= (TSCAL s) 
    {
      FlatMatrix<T,ORD>::operator= (s);
      return *this;
    }


    /// fill matrix with scalar
    Matrix & operator= (const Matrix & m2) 
    {
      SetSize (m2.Height(), m2.Width());
      FlatMatrix<T,ORD>::operator= (m2);
      return *this;
    }

    Matrix & operator= (Matrix && m2) 
    {
      this->h = m2.h;      
      this->w = m2.w;
      Swap (this->data, m2.data);
      return *this;
    }

    template<typename M2>
    Matrix & operator+= (const M2 & m) 
    { 
      CMCPMatExpr<FlatMatrix<T,ORD> >::operator+= (m);
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
    // T data[(H*W>0) ? H*W : 1];
    HTArray<H*W,T> data;
  public:
    typedef T TELEM;
    typedef typename mat_traits<T>::TSCAL TSCAL;
    typedef Vec<H, typename mat_traits<T>::TV_COL> TV_COL;
    typedef Vec<W, typename mat_traits<T>::TV_ROW> TV_ROW;

    /// do not initialize 
    constexpr Mat () = default;
    Mat (const Mat &) = default;

    
    auto & HTData() const { return data; }                                    
    template <typename T2>
    Mat (const Mat<H,W,T2> & m2) : data(m2.HTData()) { ; }

    /// assign values
    template<typename TB>
    INLINE Mat (const Expr<TB> & m)
    {
      MatExpr<Mat>::operator= (m);
    }

    /// fill with scalar
    INLINE Mat (TSCAL s)     
    {
      for (size_t i = 0; i < H*W; i++) 
        data[i] = s;
    }

    Mat (initializer_list<initializer_list<T>> llist)
    {
      int r = 0;
      for (auto row : llist)
        {
          int c = 0;
          for (auto col : row)
            (*this)(r,c++) = col;
          r++;
        }
    }

    template <typename ...TTUP>
    Mat (tuple<TTUP...> tup)
    {
      constexpr int s = tuple_size<decltype(tup)>();
      Iterate<s> ([this, tup] (auto i) { this->data[i] = std::get<i>(tup); });
    }

    template <class... T2>
    Mat (T v0, T v1, T2... rest)
      : Mat(tuple(v0, v1, rest...)) { } 

    /// assign values
    template<typename TB>
    INLINE Mat & operator= (const Expr<TB> & m)
    {
      MatExpr<Mat>::operator= (m);
      return *this;
    }

    /// copy matrix
    Mat & operator= (const Mat &) = default;
    /*
    INLINE Mat & operator= (const Mat & m) 
    {
      for (int i = 0; i < H*W; i++)
        data[i] = m.data[i];
      return *this;
    }
    */
 
    /// fill values
    INLINE Mat & operator= (TSCAL s) 
    {
      for (size_t i = 0; i < H*W; i++)
        data[i] = s;
      return *this;
    }

    INLINE auto View() const { return Mat<H,W,const T>{*this}; }
    INLINE auto & ViewRW() { return *this; }
    INLINE auto Shape() const { return tuple(IC<H>(), IC<W>()); }
    
    INLINE T* Data() noexcept { return data.Ptr(); }
    INLINE const T* Data() const noexcept { return data.Ptr(); }
    /// linear access
    INLINE TELEM & operator() (size_t i) { return data[i]; }
    /// access element
    INLINE TELEM & operator() (size_t i, size_t j) { return data[i*W+j]; }
    /// linear access
    INLINE const TELEM & operator() (size_t i) const { return data[i]; }
    /// access element
    INLINE const TELEM & operator() (size_t i, size_t j) const { return data[i*W+j]; }

    /// the height
    INLINE static constexpr auto Height ()  { return IC<H>(); }
    /// the width
    INLINE static constexpr auto Width ()  { return IC<W>(); }

    ///
    INLINE const FlatVec<W,T> Row (size_t i) 
    {
      return FlatVec<W,T> (&(*this)(i,0));
    }

    ///
    INLINE const FlatVec<W,const T> Row (size_t i) const
    {
      return FlatVec<W,const T> (&(*this)(i,0));
    }

    ///
    INLINE auto Col (size_t i) 
    {
      // return FixSliceVector<W,T> (H, &(*this)(0,i));
      return FlatSliceVec<H,W,T> (&(*this)(0,i));
    }

    ///
    /*
    INLINE const FixSliceVector<W,const T> Col (size_t i) const
    {
      return FixSliceVector<W,const T> (H, &(*this)(0,i));
    }
    */
    INLINE const FlatSliceVec<H,W,const T> Col (size_t i) const
    {
      return FlatSliceVec<H,W,const T> (&(*this)(0,i));
    }

    auto AsVector() 
    {
      return FlatVec<H*W,T> (data.Ptr());
    }

    auto AsVector() const
    {
      return FlatVec<H*W,const T> (data.Ptr());
    }

    template <typename ARCHIVE>
    void DoArchive(ARCHIVE& ar)
    {
      ar.Do(&data[0], H*W);
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
    // enum { HEIGHT = H };
    static constexpr int HEIGHT = H;
    /// matrix with
    // enum { WIDTH  = W  };
    static constexpr int WIDTH = W;
    ///
    // enum { IS_COMPLEX = mat_traits<TSCAL>::IS_COMPLEX };
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
    DiagMat (TSCAL s) 
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
    DiagMat & operator= (const DiagMat & m) 
    {
      for (int i = 0; i < H; i++)
        data[i] = m.data[i];
      return *this;
    }
 
    /// fill values
    DiagMat & operator= (TSCAL s) 
    {
      for (int i = 0; i < H; i++)
        data[i] = s;
      return *this;
    }

    // auto View() const { return DiagMat<H,const T>(*this); }
    INLINE auto View() const { return *this; } 
    INLINE tuple<size_t, size_t> Shape() const { return { H,H }; }
    
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
    constexpr size_t Height () const { return H; }
    /// the width
    constexpr size_t  Width () const { return H; }
    
    enum { IS_LINEAR = 0 };
  };






  template <int DIM, typename T=double, int DIST=DIM, ORDERING ORD=RowMajor>
  class FlatMatrixFixInner :
    public CMCPMatExpr<FlatMatrixFixInner<DIM,T,DIST,ORD>>
  {
  protected:
    /// the data
    T *  __restrict data;
    /// the height
    size_t size;
    
  public:
    /// entry type
    typedef T TELEM;
    /// scalar type of entry
    typedef typename mat_traits<T>::TSCAL TSCAL;
    
    enum { IS_LINEAR = (DIM == DIST) && (ORD==RowMajor) }; 
    enum { COL_MAJOR = (ORD==ColMajor) };

    
    INLINE FlatMatrixFixInner () { ; }    
    INLINE FlatMatrixFixInner (size_t asize, T * adata) 
      : data(adata), size(asize) { ; }

    INLINE FlatMatrixFixInner (size_t asize, LocalHeap & lh) 
      : data(lh.Alloc<T>(asize*DIST)), size(asize) { ; }



    /// copy constructor. copies pointers, not contents
    INLINE FlatMatrixFixInner (FlatMatrix<T,ORD> m) 
      : data(m.Data()), size(ORD==RowMajor?m.Height():m.Width())
    { ; }

    template <int H>
    INLINE FlatMatrixFixInner (const Mat<H,DIM,TSCAL> & m) 
      : data(const_cast<T*>(m.Data())), size(H)
    {
      static_assert(ORD==RowMajor);
    }  


    
    /// set size, and assign mem
    INLINE void AssignMemory (size_t asize, LocalHeap & lh) 
    {
      size = asize;
      data = lh.Alloc<T> (size*DIST);
    }

    /// set size, and assign mem
    INLINE void AssignMemory (size_t asize, T * mem) 
    {
      size = asize;
      data = mem;
    }


    /// copy size and pointers
    INLINE auto & Assign (const FlatMatrixFixInner & m)
    {
      size = m.size;
      data = m.data;
      return *this;
    }
    
    
    /// assign contents
    template<typename TB>
    INLINE const auto & operator= (const Expr<TB> & m) const
    {
      return CMCPMatExpr<FlatMatrixFixInner>::operator= (m);
    }

    auto Size() const { return size; }

    auto Height() const
    {
      if constexpr (ORD==RowMajor)
        return size;
      else
        return IC<DIM>();
    }

    auto Width() const
    {
      if constexpr (ORD==RowMajor)
        return IC<DIM>();
      else
        return size;
    }

    INLINE constexpr auto Dist() { return DIST; }    
    
    INLINE auto View() const { return *this; } 
    INLINE auto Shape() const { return tuple(Height(), Width()); }

    INLINE T* Data() const noexcept { return data; }
    
    /// access operator, linear access
    INLINE TELEM & operator() (size_t i) const
    {
      static_assert (DIM==DIST);
      NETGEN_CHECK_RANGE(i, 0, Height()*Width());
      return data[i]; 
    }

    INLINE TELEM * Addr(size_t i, size_t j) const
    {
      NETGEN_CHECK_RANGE(i, 0, Height());
      NETGEN_CHECK_RANGE(j, 0, Width());
      if constexpr (ORD==RowMajor)      
        return data+i*DIST+j;
      else
        return data+j*DIST+i;
    }

    
    /// access operator
    INLINE TELEM & operator() (size_t i, size_t j) const
    {
      return *Addr(i,j);
      /*
      NETGEN_CHECK_RANGE(i, 0, Height());
      NETGEN_CHECK_RANGE(j, 0, Width());
      if constexpr (ORD==RowMajor)
        return data[i*DIST+j];
      else
        return data[j*DIST+i];
      */
    }




    /// copy contents
    INLINE const auto & operator= (const FlatMatrixFixInner & m) const 
    {
      if (DIM == DIST)
        for (size_t i = 0; i < size*DIM; i++)
          data[i] = m(i);
      else
        for (size_t i = 0; i < size; i++)
          for (size_t j = 0; j < DIM; j++)
            data[i*DIST+j] = m.data[i*DIST+j];
            // (*this)(i,j) = m(i,j);
        
      return *this;
    }

    /// assign constant
    INLINE const auto & operator= (TSCAL s) const 
    {
      if (DIM == DIST)
        for (size_t i = 0; i < size*DIM; i++)
          data[i] = s; 
      else
        for (size_t i = 0; i < size; i++)
          for (size_t j = 0; j < DIM; j++)
            data[i*DIST+j] = s;
            // (*this)(i,j) = s;
        
      return *this;
    }



    ///
    INLINE const auto Row (size_t i) const
    {
      if constexpr (ORD==RowMajor)
        return FlatVec<DIM,T> (Addr(i,0));
      else
        return SliceVector<T> (size, DIST, Addr(i,0));
    }

    INLINE const auto Col (size_t i) const
    {
      if constexpr (ORD==RowMajor)
        return FixSliceVector<DIST,T> (size, Addr(0,i));
      else
        return FlatVec<DIM,T> ( Addr(0,i) );       
    }

    using CMCPMatExpr<FlatMatrixFixInner<DIM,T,DIST,ORD>>::Rows;
    using CMCPMatExpr<FlatMatrixFixInner<DIM,T,DIST,ORD>>::Cols;

    INLINE const auto Rows (size_t first, size_t next) const
    {
      if constexpr (ORD==RowMajor)      
        return FlatMatrixFixInner (next-first, Addr(first, 0));
      else
        return SliceMatrix<T,ORD> (next-first, size, DIST, Addr(first,0));
    }

    INLINE const auto Cols (size_t first, size_t next) const
    {
      if constexpr (ORD==RowMajor)            
        return SliceMatrix<T,ORD> (size, next-first, DIST, Addr(0,first));
      else
        return FlatMatrixFixInner (next-first, Addr(0,first));        
    }

    INLINE const auto Rows (IntRange range) const
    {
      return Rows (range.First(), range.Next());
    }

    INLINE const auto Cols (IntRange range) const
    {
      return Cols (range.First(), range.Next());
    }

    
    
    ///
    INLINE operator const FlatMatrix<T,ORD>() const
    {
      static_assert(DIM==DIST);
      return FlatMatrix<T,ORD> (Height(), Width(), data);
    }
    
    operator SliceMatrix<T,ORD> () const
    {
      return SliceMatrix<T,ORD> (Height(), Width(), DIST, data);
    }
    
    operator BareSliceMatrix<T,ORD> () const
    {
      return BareSliceMatrix<T,ORD> (DIST, data, DummySize(Shape()));
    }
    
    
    auto AsVector() const
    {
      static_assert(DIM==DIST);
      return FlatVector<T> (Height()*Width(), data);
    }
    
    
  };




  template <int W, typename T = double, int DIST=W>
  using FlatMatrixFixWidth = FlatMatrixFixInner<W,T,DIST,RowMajor>;

  

  /**
     A Matrix class with memory allocation/deallocation
  */
  template <int W, typename T = double>
  class MatrixFixWidth : public FlatMatrixFixWidth<W,T>
  {
  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    INLINE MatrixFixWidth () : FlatMatrixFixWidth<W,T> (0, nullptr) { ; }

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
      FlatMatrixFixWidth<W,T>::operator= (m);
      return *this;
    }

    INLINE MatrixFixWidth & operator= (const MatrixFixWidth & m) 
    {
      FlatMatrixFixWidth<W,T>::operator= (m);
      return *this;
    }


    INLINE MatrixFixWidth & operator= (MatrixFixWidth && m)
    {
      this->h = m.h;
      Swap (this->data, m.data);
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

  template <int H, typename T = double, int DIST=H>
  using FlatMatrixFixHeight = FlatMatrixFixInner<H,T,DIST,ColMajor>;
  


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
  class SliceMatrix : public CMCPMatExpr<SliceMatrix<T,ORD>>
  {
  protected:
    /// the height
    size_t h;
    /// the width
    size_t w;
    /// the distance
    size_t dist;
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
    // INLINE SliceMatrix & operator= (const SliceMatrix &) = delete;

    /// set height, width, and mem
    INLINE SliceMatrix (size_t ah, size_t aw, size_t adist, T * adata) 
      : h(ah), w(aw), dist(adist), data(adata) { ; }

    SliceMatrix (FlatMatrix<T,ORD> mat)
      : h(mat.Height()), w(mat.Width()), dist(mat.Dist()), data(mat.Data())
    { ; }


    template<int W, int DIST>
    SliceMatrix (FlatMatrixFixInner<W,T,DIST,ORD> mat)
      : h(mat.Height()), w(mat.Width()), dist(mat.Dist()), data(mat.Data())
    { ; }


    template<int H, int W>
    INLINE SliceMatrix (Mat<H,W,T> & mat)
      : h(mat.Height()), w(mat.Width()), dist(mat.Width()), data(mat.Data())
    {
      static_assert(ORD==RowMajor);
    }

    template<int H, int W>    
    INLINE SliceMatrix(const Mat<H,W,const T> & mat) 
      : h(mat.Height()), w(mat.Width()), dist(mat.Width()), data(const_cast<T*>(mat.Data()))
    {
      static_assert(ORD==RowMajor);
    } 
    

    /// assign contents
    template<typename TB>
    INLINE const SliceMatrix & operator= (const Expr<TB> & m) const
    {
      return CMCPMatExpr<SliceMatrix<T,ORD>>::operator= (m);
    }

    INLINE const SliceMatrix & operator= (const SliceMatrix & m) const
    {
      return CMCPMatExpr<SliceMatrix<T,ORD>>::operator= (m);
    }

    
    /// assign constant
    INLINE const SliceMatrix & operator= (TSCAL s) const
    {
      /*
      if (w == 0) return *this;
      for (size_t i = 0; i < h; i++)
        {
          __assume (w > 0);
          for (size_t j = 0; j < w; j++)
            data[i*dist+j] = s;
        }
      */
      size_t hw = (ORD == RowMajor) ? w : h;
      size_t hh = (ORD == RowMajor) ? h : w;
      
      if (hw == 0) return *this;
      size_t i = 0, base = 0;
      for ( ; i+1 < hh; i+=2, base += 2*dist)
        {
          __assume (hw > 0);
          for (auto j : Range(hw))
            {
              data[base+j] = s;
              data[base+dist+j] = s;
            }
        }
      if (i < hh)
        {
          __assume (hw > 0);
          for (auto j : Range(hw))            
            data[base+j] = s;
        }
      return *this;
    }

    INLINE auto View() const { return *this; }     
    INLINE auto Shape() const { return tuple(h,w); }

    INLINE TELEM * Addr(size_t i, size_t j) const
    {
      if constexpr (ORD==RowMajor)      
        return data+i*dist+j;
      else
        return data+j*dist+i;
    }

    
    /// access operator
    INLINE TELEM & operator() (size_t i, size_t j) const
    {
      NETGEN_CHECK_RANGE(i, 0, Height());
      NETGEN_CHECK_RANGE(j, 0, Width());
      return *Addr(i,j);
    }

    /// access operator, linear access
    INLINE TELEM & operator() (size_t i) const
    {
      NETGEN_CHECK_RANGE(i, 0, Height()*Dist());
      return data[i]; 
    }

    INLINE auto RowDist() const
    {
      if constexpr (ORD==RowMajor)
        return IC<1>();
      else
        return dist;
    }

    INLINE auto ColDist() const
    {
      if constexpr (ORD==RowMajor)
        return dist;
      else
        return IC<1>();
    }

    
    INLINE auto Height () const { return h; }
    INLINE auto Width () const { return w; }
    INLINE auto Dist () const { return dist; }
    INLINE T* Data() const  { return data; }

    using CMCPMatExpr<SliceMatrix>::Rows;
    using CMCPMatExpr<SliceMatrix>::Cols;


    INLINE auto SplitRows (size_t split) const
    {
      return tuple(Rows(0,split), Rows(split, Height()));
    }

    INLINE auto SplitCols (size_t split) const
    {
      return tuple(Cols(0,split), Cols(split, Width()));
    }
    
    INLINE const auto Row (size_t i) const
    {
      return make_VectorView(Width(), RowDist(), Addr(i,0));
    }

    INLINE const auto  Col (size_t i) const
    {
      return make_VectorView(Height(), ColDist(), Addr(0,i));
    }

    INLINE const auto Cols (size_t first, size_t next) const
    {
      return SliceMatrix (h, next-first, dist, Addr(0,first));
    }

    INLINE const auto Rows (size_t first, size_t next) const
    {
      return SliceMatrix (next-first, w, dist, Addr(first, 0));
    }

    INLINE const auto Rows (IntRange range) const
    {
      return Rows (range.First(), range.Next());
    }

    INLINE const auto Cols (IntRange range) const
    {
      return Cols (range.First(), range.Next());
    }

    INLINE const SliceVector<T> Diag () const
    {
      if (ORD == RowMajor)
        return SliceVector<T> (h, dist+1, data);
      else
        return SliceVector<T> (w, dist+1, data);
    }

    INLINE const SliceVector<T> Diag (int offset) const
    {
      // return SliceVector<T> (h, dist+1, data);
      int dp = std::max(offset, 0);
      int dm = std::min(offset, 0);
      if (ORD == RowMajor)      
        return SliceVector<T> (min(w-dp, h+dm), dist+1, data+dp-dm*dist);      
      else
        return SliceVector<T> (min(w-dp, h+dm), dist+1, data-dm+dp*dist);
    }
  };






  
  
  template <typename T, ORDERING ORDER>
  SliceMatrix<T,ORDER> make_SliceMatrix (FlatMatrix<T,ORDER> mat) { return mat; }

  template <int DIM, typename T, int DIST, ORDERING ORDER>
  SliceMatrix<T,ORDER> make_SliceMatrix (FlatMatrixFixInner<DIM,T,DIST,ORDER> mat) { return mat; }

  template <int H, int W, typename T>
  SliceMatrix<T,RowMajor> make_SliceMatrix (const Mat<H,W,T> &mat) { return const_cast<Mat<H,W,T>&>(mat); }
  
  template <typename T, ORDERING ORDER>
  SliceMatrix<T,ORDER> make_SliceMatrix (SliceMatrix<T,ORDER> mat) { return mat; }





  template <typename T, ORDERING ORD> 
  class BareSliceMatrix : public CMCPMatExpr<BareSliceMatrix<T,ORD>>, DummySize
  {
  protected:
    /// the distance
    size_t dist;
    /// the data
    T * __restrict data;
  public:
    
    /// element type
    typedef T TELEM;
    /// scalar type of elements (double or Complex)
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { IS_LINEAR = 0 };

    // 
    BareSliceMatrix() : DummySize(0,0) { ; } // initialize with placement new later
    INLINE BareSliceMatrix(const BareSliceMatrix &) = default;

    BareSliceMatrix (FlatMatrix<T,ORD> mat)
      : DummySize(mat.Shape()), dist(mat.Dist()), data(mat.Data())
    { ; }

    BareSliceMatrix (SliceMatrix<T,ORD> mat)
      : DummySize(mat.Shape()), dist(mat.Dist()), data(mat.Data())
    { ; }

    template<int H, int W>
    BareSliceMatrix (Mat<H,W,T> & mat)
      : DummySize(mat.Shape()), dist(mat.Width()), data(mat.Data())
    {
      static_assert(ORD==RowMajor);
    }

    template<int H, int W>    
    INLINE BareSliceMatrix(const Mat<H,W,const T> & mat) 
      : DummySize(mat.Shape()), dist(mat.Width()), data(const_cast<T*>(mat.Data()))
    {
      static_assert(ORD==RowMajor);
    } 
    
    template<int W, int DIST>
    BareSliceMatrix (FlatMatrixFixInner<W,T,DIST,ORD> mat)
      : DummySize(mat.Shape()), dist(mat.Dist()), data(mat.Data()) { } 


    
    BareSliceMatrix (size_t adist, T * adata, DummySize ds) : DummySize(ds), dist(adist), data(adata) { ; } 
    
    BareSliceMatrix & operator= (const BareSliceMatrix & m) = delete;

    INLINE auto View() const { return *this; } 


    INLINE TELEM * Addr(size_t i, size_t j) const
    {
      if constexpr (ORD==RowMajor)      
        return data+i*dist+j;
      else
        return data+j*dist+i;
    }

    /// access operator
    INLINE TELEM & operator() (size_t i, size_t j) const
    {
      NETGEN_CHECK_RANGE(i, 0, Height());
      NETGEN_CHECK_RANGE(j, 0, Width());
      return *Addr(i,j);
    }

    /// access operator
    INLINE TELEM & operator() (size_t i) const
    {
      NETGEN_CHECK_RANGE(i, 0, Width()*Height());
      return data[i]; 
    }
    using DummySize::Height;
    using DummySize::Width;


    INLINE TELEM* Data() const { return data; }
    INLINE auto Shape() const { return tuple(Height(), Width()); }

    INLINE auto RowDist() const
    {
      if constexpr (ORD==RowMajor)
        return IC<1>();
      else
        return dist;
    }

    INLINE auto ColDist() const
    {
      if constexpr (ORD==RowMajor)
        return dist;
      else
        return IC<1>();
    }

    
    /// 
    INLINE size_t Dist () const  { return dist; }
    void IncPtr (size_t inc) { data += inc; } 
    const SliceMatrix<T,ORD> AddSize (size_t h, size_t w) const
    {
      NETGEN_CHECK_RANGE(h, Height(), Height()+1);
      NETGEN_CHECK_RANGE(w, Width(), Width()+1);
      return SliceMatrix<T,ORD> (h, w, dist, data);
    }
    
    INLINE const auto Rows (size_t first, size_t next) const
    {
      NETGEN_CHECK_RANGE(first, 0, first==next ? Height()+1 : Height()); // always allow Rows(0,0)
      NETGEN_CHECK_RANGE(next, 0, Height()+1);
      return BareSliceMatrix ( dist, Addr(first, 0), DummySize(next-first, Width()));
    }

    INLINE const auto Col (size_t col) const
    {
      NETGEN_CHECK_RANGE(col, 0, Width());
      return make_VectorView(Height(), ColDist(), Addr(0,col));
      /*
      if constexpr (ORD==RowMajor)
        return BareSliceVector<T> (Addr(0,col), dist, DummySize(Height()));
      else
        return BareVector<T> (Addr(0,col), DummySize(Height()));
      */
    }
    
    INLINE const auto Row (size_t i) const
    {
      NETGEN_CHECK_RANGE(i, 0, Height());
      return make_VectorView(Width(), RowDist(), Addr(i,0));
      /*
      if constexpr (ORD==RowMajor)
        return BareVector<T> (Addr(i,0), DummySize(Width()));
      else
        return BareSliceVector<T> (Addr(i,0), dist, DummySize(Width()));
      */
    }

    /*
    INLINE const FlatVector<T> Diag (size_t i) const
    {
      return SliceVector<T> (h, dist+1, data);
    }
    */
    
    INLINE const BareSliceMatrix Cols (size_t first, size_t next) const
    {
      NETGEN_CHECK_RANGE(first, 0, first==next ? Width()+1 : Width()); // always allow Cols(0,0)
      NETGEN_CHECK_RANGE(next, 0, Width()+1);
      return BareSliceMatrix (dist, Addr(0,first), DummySize(Height(), next-first));
    }

    INLINE const auto Rows (IntRange range) const
    {
      return Rows (range.First(), range.Next());
    }

    INLINE const auto Cols (IntRange range) const
    {
      return Cols (range.First(), range.Next());
    }
    /*
    INLINE const SliceVector<T> Diag () const
    {
      return SliceVector<T> (h, dist+1, &data[0]);
    }
    */
    const auto RowSlice(size_t first, size_t adist) const
    {
      static_assert (ORD==RowMajor);
      // NETGEN_CHECK_RANGE(first, 0, Height());  // too restrictive
      NETGEN_CHECK_RANGE(first, 0, adist);  
      return BareSliceMatrix<T,ORD> (dist*adist, data+first*dist, DummySize( Height()/adist, Width()));
    }
    
  };



  


  template <typename T = double>
  class DoubleSliceMatrix : public CMCPMatExpr<DoubleSliceMatrix<T> >
  {
  protected:
    /// the height
    size_t h;
    /// the width
    size_t w;
    /// the distance
    size_t distr, distc;
    /// the data
    T * data;
  public:

    /// element type
    typedef T TELEM;
    /// scalar type of elements (double or Complex)
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { IS_LINEAR = 0 };

    /// set height, width, and mem
    DoubleSliceMatrix (int ah, int aw, int adistr, int adistc, T * adata) 
      : h(ah), w(aw), distr(adistr), distc(adistc), data(adata) { ; }
  
    /// assign contents
    template<typename TB>
    const DoubleSliceMatrix & operator= (const Expr<TB> & m) const
    {
      return CMCPMatExpr<DoubleSliceMatrix>::operator= (m);
    }

    /// assign constant
    DoubleSliceMatrix & operator= (TSCAL s) 
    {
      for (size_t i = 0; i < h; i++)
        for (size_t j = 0; j < w; j++)
          data[i*distr+j*distc] = s;
      return *this;
    }

    INLINE auto View() const { return *this; }
    INLINE tuple<size_t,size_t> Shape() const { return { h, w }; }
    
    /// access operator
    TELEM & operator() (size_t i, size_t j) const
    {
      NETGEN_CHECK_RANGE(i, 0, Height());
      NETGEN_CHECK_RANGE(j, 0, Width());
      return data[i*distr+j*distc]; 
    }

    /// access operator, linear access
    TELEM & operator() (size_t i) const
    {
      NETGEN_CHECK_RANGE(i, 0, Height()*Width());
      return data[i]; 
    }

    size_t Height () const { return h; }
    size_t Width () const  { return w; }
    size_t DistRow () const { return distr; }
    size_t DistCol () const { return distc; }
    T * Data () const { return data; }

    const DoubleSliceMatrix Rows (size_t first, size_t next) const
    {
      return DoubleSliceMatrix (next-first, w, distr, distc, data+first*distr);
    }

    const DoubleSliceMatrix<T> Cols (size_t first, size_t next) const
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
  };

  /*
  template <typename T>
  struct trivtrans { static constexpr bool value = false; };
  template<> struct trivtrans<double> { static constexpr bool value = true; };
  template<> struct trivtrans<Complex> { static constexpr bool value = true; };
  template <int D, typename T>
  struct trivtrans<AutoDiff<D,T>> { static constexpr bool value = true; };
  template <int D, typename T>
  struct trivtrans<AutoDiffDiff<D,T>> { static constexpr bool value = true; };
  template <int D>
  struct trivtrans<SIMD<double,D>> { static constexpr bool value = true; };

  
  template <typename T>
  constexpr bool IsTrivialTranspose ()
  {
    return trivtrans<T>::value;
  }
  */

  template <typename T>
  constexpr bool IsTrivialTranspose () { return IsScalar<T>(); } 
  
  
  constexpr ORDERING operator! (ORDERING ordering)
  {
    return ordering == ColMajor ?  RowMajor : ColMajor;
  }



  template <typename T, ORDERING ord,
            typename enable_if<IsTrivialTranspose<T>(),int>::type = 0>
  INLINE auto Trans (FlatMatrix<T,ord> mat)
  {
    return FlatMatrix<T,!ord> (mat.Width(), mat.Height(), mat.Data());
  }

  template <typename T, ORDERING ord,
            typename enable_if<IsTrivialTranspose<T>(),int>::type = 0>
  INLINE auto Trans (const Matrix<T,ord> & mat)
  {
    return FlatMatrix<T,!ord> (mat.Width(), mat.Height(), mat.Data());
  }


  template <typename T, ORDERING ord,
            typename enable_if<IsTrivialTranspose<T>(),int>::type = 0>
  INLINE const SliceMatrix<T,!ord> Trans (SliceMatrix<T,ord> mat)
  {
    return SliceMatrix<T,!ord> (mat.Width(), mat.Height(), mat.Dist(), mat.Data());
  }

  

  template <typename T, ORDERING ord,
            typename enable_if<IsTrivialTranspose<T>(),int>::type = 0>
  INLINE const BareSliceMatrix<T,!ord> Trans (BareSliceMatrix<T,ord> mat)
  {
    return BareSliceMatrix<T,!ord> (mat.Dist(), mat.Data(), DummySize(mat.Width(), mat.Height()));
  }



  template <int DIM, typename T, int DIST, ORDERING ord,
            typename enable_if<IsTrivialTranspose<T>(),int>::type = 0>
  INLINE const auto Trans (FlatMatrixFixInner<DIM,T,DIST,ord> mat)
  {
    return FlatMatrixFixInner<DIM,T,DIST,!ord> (mat.Size(), mat.Data());
  }

  
  // should get rid of that ...
  template <int H, int DIST>
  INLINE const FlatMatrixFixWidth<H,double,DIST> Trans (FlatMatrixFixHeight<H,double,DIST> mat)
  {
    return FlatMatrixFixWidth<H,double,DIST> (mat.Width(), mat.Data());
  }

  template <int H, int DIST>
  INLINE const FlatMatrixFixWidth<H,Complex,DIST> Trans (FlatMatrixFixHeight<H,Complex,DIST> mat)
  {
    return FlatMatrixFixWidth<H,Complex,DIST> (mat.Width(), mat.Data());
  }






  template <class TM, class TSCAL>
  class Scalar2ElemMatrix 
  {
  public:
    const BareSliceMatrix<TSCAL> mat;
    Scalar2ElemMatrix (const BareSliceMatrix<TSCAL> amat) : mat(amat) { ; }

    /*
    enum { H = mat_traits<TM>::HEIGHT };
    enum { W = mat_traits<TM>::WIDTH };
    */
    enum { H = Height<TM>() };
    enum { W = Width<TM>() };
    
    TM operator() (size_t i, size_t j) const
    {
      /*
      TM ret;
      for (size_t k = 0; k < H; k++)
	for (size_t l = 0; l < W; l++)
	  Access(ret, k,l) = mat(i*H+k, j*W+l);
      return ret;
      */
      if constexpr (IsScalar<TM>())
        {
          return mat(i,j);
        }
      else
        {
          TM ret;
          for (size_t k = 0; k < H; k++)
            for (size_t l = 0; l < W; l++)
              ret(k,l) = mat(i*H+k, j*W+l);
          return ret;
        }
    }

    Scalar2ElemMatrix Rows(size_t first, size_t next) const
    {
      return Scalar2ElemMatrix(mat.Rows(H*first, H*next));
    }
  };






  /**
     Identity Matrix of fixed size
  */
  template <size_t H>
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

    INLINE auto View() const { return Id(); }
    INLINE tuple<size_t, size_t> Shape() const { return { H,H }; }    
    ///
    double operator() (int i) const
    {
      static_assert (true, "linear access of id");
      cerr << "id, linear access" << endl; return 0;
    }
    ///
    double operator() (int i, int j) const { return (i == j) ? 1 : 0; }

    /// the height
    constexpr size_t Height () const { return H; }
    /// the width
    constexpr size_t Width () const { return H; }
  };




  /// Variable size identity matrix
  class Identity : public MatExpr<Identity >
  {
    int size;
  public:
    typedef double TELEM;
    typedef double TSCAL;
    enum { IS_LINEAR = 0 };

    INLINE Identity (int s) : size(s) { ; }

    INLINE double operator() (int i) const
    { cerr << "Identity, linear access" << endl; return 0; }

    INLINE double operator() (int i, int j) const { return (i == j) ? 1 : 0; }
    INLINE auto View() const { return Identity(size); }
    INLINE tuple<size_t, size_t> Shape() const { return { size, size }; }    
    INLINE int Height () const { return size; }
    INLINE int Width () const { return size; }
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
    /*
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        res(j,i) = mat(i,j);
    */
    Iterate<H> ([&] (auto i) {
        Iterate<W> ([&] (auto j) {
            res(j.value,i.value) = mat(i.value, j.value);
          });
      });
    return res;
  }
  
  template <int H, int W, typename T>
  INLINE Mat<H,W,T> operator- (const Mat<H,W,T> & mat)
  {
    Mat<H,W,typename remove_const<T>::type> res;
    Iterate<H*W> ([&] (auto i) {
        res(i.value) = -mat(i.value);
      });
    return res;
  }


  
  template <int H, int W, typename T>
  INLINE Mat<H,W,T> operator+ (const Mat<H,W,T> & ma, const Mat<H,W,T> & mb)
  {
    Mat<H,W,typename remove_const<T>::type> res;
    Iterate<H*W> ([&] (auto i) {
        res(i.value) = ma(i.value) + mb(i.value);
      });
    return res;
  }
  
  template <int H, int W, typename T>
  INLINE Mat<H,W,T> operator- (const Mat<H,W,T> & ma, const Mat<H,W,T> & mb)
  {
    Mat<H,W,typename remove_const<T>::type> res;
    Iterate<H*W> ([&] (auto i) {
        res(i.value) = ma(i.value) - mb(i.value);
      });
    return res;
  }

  template <int H, int W, typename T>
  INLINE Mat<H,W,T> & operator+= (Mat<H,W,T> & ma, const Mat<H,W,T> & mb)
  {
    Iterate<H*W> ([&] (auto i) {
        ma(i.value) += mb(i.value);
      });
    return ma;
  }

  template <int H, int W, typename T>
  INLINE Mat<H,W,T> & operator-= (Mat<H,W,T> & ma, const Mat<H,W,T> & mb)
  {
    Iterate<H*W> ([&] (auto i) {
        ma(i.value) -= mb(i.value);
      });
    return ma;
  }

  

  template <int H, int W, typename T>
  INLINE Mat<H,W,T> operator* (T scal, const Mat<H,W,T> & mat)
  {
    Mat<H,W,typename remove_const<T>::type> res;
    Iterate<H*W> ([&] (auto i) {
        res(i.value) = scal * mat(i.value);
      });
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
    static_assert(W == W2, "Mat * Vec dimension mismatch!");
    typedef decltype(RemoveConst(mat(0,0)*vec(0))) TRES;
    Vec<H, TRES> res = TRES(0);
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        res(i) += mat(i,j) * vec(j);
    return res;
  }

  template <int H, int W, int W2, typename T1, typename T2>
  INLINE auto operator* (const Mat<H,W,T1> & mat, const FlatVec<W2,T2> & vec)
    -> Vec<H, decltype(RemoveConst(mat(0,0)*vec(0)))>
  {
    static_assert(W == W2, "Mat * FlatVec dimension mismatch!");
    typedef decltype(RemoveConst(mat(0,0)*vec(0))) TRES;
    Vec<H, TRES> res = TRES(0);
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        res(i) += mat(i,j) * vec(j);
    return res;
  }

  template <int H, int W, typename T1, typename T2>
  INLINE auto operator* (const Mat<H,W,T1> & mat, FlatVector<T2> vec) 
    -> Vec<H, decltype(RemoveConst(mat(0,0)*vec(0)))>
  {
    NETGEN_CHECK_RANGE(vec.Size(), W, W+1);
    typedef decltype(RemoveConst(mat(0,0)*vec(0))) TRES;
    Vec<H, TRES> res = TRES(0);
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        res(i) += mat(i,j) * vec(j);
    return res;
  }

  template <typename T>
  auto Det (Mat<0,0,T> m)
  {
    return T(1.0);
  }


  template <typename T>
  auto Det (Mat<1,1,T> m)
  {
    return m(0,0);
  }

  template <typename T>
  auto Det (Mat<2,2,T> m)
  {
    return m(0,0)*m(1,1) - m(0,1)*m(1,0); 
  }

  template <typename T>
  auto Det (Mat<3,3,T> m)
  {
    return 
      m(0) * (m(4) * m(8) - m(5) * m(7)) +
      m(1) * (m(5) * m(6) - m(3) * m(8)) +
      m(2) * (m(3) * m(7) - m(4) * m(6));
  }


  


  template <int H, int W, typename T>
  INLINE Mat<H,W,T> Adj (Mat<H,W,T> m)
  {
    cerr << "Adj<" << H << "," << W << "> not implemented" << endl;
    return m;
  }


  template <typename T>
  INLINE Mat<0,0,T> Adj (Mat<0,0,T> m)
  {
    return Mat<0,0,T>();
  }

  template <typename T>
  INLINE Mat<1,1,T> Adj (Mat<1,1,T> m)
  {
    Mat<1,1,T> adj;
    adj(0,0) = T(1); // m(0,0);
    return adj;
  }

  template <typename T>
  INLINE Mat<2,2,T> Adj (Mat<2,2,T> m)
  {
    Mat<2,2,T> adj;
    adj(0,0) = m(1,1);
    adj(0,1) = -m(0,1);
    adj(1,0) = -m(1,0);
    adj(1,1) = m(0,0);
    return adj;
  }


  template <typename T>
  INLINE Mat<3,3,T> Adj (Mat<3,3,T> m)
  {
    Mat<3,3,T> adj;
    adj(0,0) =  m(4)*m(8)-m(5)*m(7);
    adj(0,1) = -(m(1) * m(8) - m(2) * m(7));
    adj(0,2) =  m(1) * m(5) - m(2) * m(4);
    
    adj(1,0) =  m(5)*m(6)-m(3)*m(8);
    adj(1,1) =  m(0) * m(8) - m(2) * m(6);
    adj(1,2) = -(m(0) * m(5) - m(2) * m(3));
    
    adj(2,0) =  (m(3)*m(7)-m(4)*m(6));
    adj(2,1) = -(m(0) * m(7) - m(1) * m(6));
    adj(2,2) =  (m(0) * m(4) - m(1) * m(3));
    return adj;
  }



  template <int H, int W, typename T>
  INLINE Mat<H,W,T> Cof (Mat<H,W,T> m)
  {
    cerr << "Cof<" << H << "," << W << "> not implemented" << endl;
    return m;
  }

  template <typename T>
  INLINE Mat<0,0,T> Cof (Mat<0,0,T> m)
  {
    return Mat<0,0,T>();
  }

  template <typename T>
  INLINE Mat<1,1,T> Cof (Mat<1,1,T> m)
  {
    Mat<1,1,T> cof;
    cof(0,0) = T(1); // m(0,0);
    return cof;
  }

  template <typename T>
  INLINE Mat<2,2,T> Cof (Mat<2,2,T> m)
  {
    Mat<2,2,T> cof;
    cof(0,0) = m(1,1);
    cof(0,1) = -m(1,0);
    cof(1,0) = -m(0,1);
    cof(1,1) = m(0,0);
    return cof;
  }


  template <typename T>
  INLINE Mat<3,3,T> Cof (Mat<3,3,T> m)
  {
    Mat<3,3,T> cof;
    cof(0,0) =  m(1,1)*m(2,2)-m(2,1)*m(1,2);
    cof(0,1) = -m(1,0)*m(2,2)+m(2,0)*m(1,2);
    cof(0,2) =  m(1,0)*m(2,1)-m(2,0)*m(1,1);
    
    cof(1,0) = -m(0,1)*m(2,2)+m(2,1)*m(0,2); 
    cof(1,1) =  m(0,0)*m(2,2)-m(2,0)*m(0,2); 
    cof(1,2) = -m(0,0)*m(2,1)+m(2,0)*m(0,1);
    
    cof(2,0) =  m(0,1)*m(1,2)-m(1,1)*m(0,2); 
    cof(2,1) = -m(0,0)*m(1,2)+m(1,0)*m(0,2); 
    cof(2,2) =  m(0,0)*m(1,1)-m(1,0)*m(0,1);
    return cof;
  }


  template <typename T>
  INLINE Mat<4,4,T> Cof (Mat<4,4,T> m)
  {
    Mat<4,4,T> cof;
    cof(0,0) =   (m(1,1)*m(2,2)*m(3,3)+m(1,2)*m(2,3)*m(3,1)+m(1,3)*m(2,1)*m(3,2) - m(1,1)*m(3,2)*m(2,3) - m(2,1)*m(1,2)*m(3,3) - m(3,1)*m(2,2)*m(1,3));
    cof(0,1) =  -(m(1,0)*m(2,2)*m(3,3)+m(1,2)*m(2,3)*m(3,0)+m(1,3)*m(2,0)*m(3,2) - m(1,0)*m(3,2)*m(2,3) - m(2,0)*m(1,2)*m(3,3) - m(3,0)*m(2,2)*m(1,3));
    cof(0,2) =   (m(1,0)*m(2,1)*m(3,3)+m(1,1)*m(2,3)*m(3,0)+m(1,3)*m(2,0)*m(3,1) - m(1,0)*m(3,1)*m(2,3) - m(2,0)*m(1,1)*m(3,3) - m(3,0)*m(2,1)*m(1,3));
    cof(0,3) =  -(m(1,0)*m(2,1)*m(3,2)+m(1,1)*m(2,2)*m(3,0)+m(1,2)*m(2,0)*m(3,1) - m(1,0)*m(3,1)*m(2,2) - m(2,0)*m(1,1)*m(3,2) - m(3,0)*m(2,1)*m(1,2));

    cof(1,0) =  -(m(0,1)*m(2,2)*m(3,3)+m(0,2)*m(2,3)*m(3,1)+m(0,3)*m(2,1)*m(3,2) - m(0,1)*m(3,2)*m(2,3) - m(2,1)*m(0,2)*m(3,3) - m(3,1)*m(2,2)*m(0,3));
    cof(1,1) =   (m(0,0)*m(2,2)*m(3,3)+m(0,2)*m(2,3)*m(3,0)+m(0,3)*m(2,0)*m(3,2) - m(0,0)*m(3,2)*m(2,3) - m(2,0)*m(0,2)*m(3,3) - m(3,0)*m(2,2)*m(0,3));
    cof(1,2) =  -(m(0,0)*m(2,1)*m(3,3)+m(0,1)*m(2,3)*m(3,0)+m(0,3)*m(2,0)*m(3,1) - m(0,0)*m(3,1)*m(2,3) - m(2,0)*m(0,1)*m(3,3) - m(3,0)*m(2,1)*m(0,3));
    cof(1,3) =   (m(0,0)*m(2,1)*m(3,2)+m(0,1)*m(2,2)*m(3,0)+m(0,2)*m(2,0)*m(3,1) - m(0,0)*m(3,1)*m(2,2) - m(2,0)*m(0,1)*m(3,2) - m(3,0)*m(2,1)*m(0,2));

    cof(2,0) =   (m(0,1)*m(1,2)*m(3,3)+m(0,2)*m(1,3)*m(3,1)+m(0,3)*m(1,1)*m(3,2) - m(0,1)*m(3,2)*m(1,3) - m(1,1)*m(0,2)*m(3,3) - m(3,1)*m(1,2)*m(0,3));
    cof(2,1) =  -(m(0,0)*m(1,2)*m(3,3)+m(0,2)*m(1,3)*m(3,0)+m(0,3)*m(1,0)*m(3,2) - m(0,0)*m(3,2)*m(1,3) - m(1,0)*m(0,2)*m(3,3) - m(3,0)*m(1,2)*m(0,3));
    cof(2,2) =   (m(0,0)*m(1,1)*m(3,3)+m(0,1)*m(1,3)*m(3,0)+m(0,3)*m(1,0)*m(3,1) - m(0,0)*m(3,1)*m(1,3) - m(1,0)*m(0,1)*m(3,3) - m(3,0)*m(1,1)*m(0,3));
    cof(2,3) =  -(m(0,0)*m(1,1)*m(3,2)+m(0,1)*m(1,2)*m(3,0)+m(0,2)*m(1,0)*m(3,1) - m(0,0)*m(3,1)*m(1,2) - m(1,0)*m(0,1)*m(3,2) - m(3,0)*m(1,1)*m(0,2));

    cof(3,0) =  -(m(0,1)*m(1,2)*m(2,3)+m(0,2)*m(1,3)*m(2,1)+m(0,3)*m(1,1)*m(2,2) - m(0,1)*m(2,2)*m(1,3) - m(1,1)*m(0,2)*m(2,3) - m(2,1)*m(1,2)*m(0,3));
    cof(3,1) =   (m(0,0)*m(1,2)*m(2,3)+m(0,2)*m(1,3)*m(2,0)+m(0,3)*m(1,0)*m(2,2) - m(0,0)*m(2,2)*m(1,3) - m(1,0)*m(0,2)*m(2,3) - m(2,0)*m(1,2)*m(0,3));
    cof(3,2) =  -(m(0,0)*m(1,1)*m(2,3)+m(0,1)*m(1,3)*m(2,0)+m(0,3)*m(1,0)*m(2,1) - m(0,0)*m(2,1)*m(1,3) - m(1,0)*m(0,1)*m(2,3) - m(2,0)*m(1,1)*m(0,3));
    cof(3,3) =   (m(0,0)*m(1,1)*m(2,2)+m(0,1)*m(1,2)*m(2,0)+m(0,2)*m(1,0)*m(2,1) - m(0,0)*m(2,1)*m(1,2) - m(1,0)*m(0,1)*m(2,2) - m(2,0)*m(1,1)*m(0,2));
    return cof;
  }


  template <int H, int W, typename T>
  INLINE Mat<H,W,T> Inv (Mat<H,W,T> m)
  {
    return 1.0/Det(m) * Adj(m);
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
  // if not, we return a reference to the original matrix
  // 
  
  template <typename T, typename TMAT>
  auto SliceIfPossible (const TMAT & mat) -> typename slicetype<Is_Sliceable<TMAT,T>::VAL,TMAT>::TYPE
  {
    // cout << "type(tmat) = " << typeid(TMAT).name() << endl;
    // cout << "sliceable = " << Is_Sliceable<TMAT,T>::VAL << endl;
    // cout << "return type = " << typeid (typename slicetype<Is_Sliceable<TMAT,T>::VAL,TMAT>::TYPE).name() << endl;
    return mat;
  }

  


  template <int H, int W, typename SCAL, typename TANY>
  inline void AtomicAdd (Mat<H,W,SCAL> & x, TANY y)
  {
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        AtomicAdd (x(i,j), y(i,j));
  }

}


namespace ngstd
{
  template <typename ARCHIVE, int N, int M, typename T>
  inline auto & operator& (ARCHIVE & ar, ngbla::Mat<N,M,T> & m)
  {
    for (int i = 0; i < N*M; i++)
      ar & m(i);
    return ar;
  }
}


#ifdef PARALLEL
namespace ngcore
{
  template<int N, int M, typename T>
  class MPI_typetrait<ngbla::Mat<N, M, T> >
  {
  public:
    /// gets the MPI datatype
    static MPI_Datatype MPIType () 
    { 
      static MPI_Datatype MPI_T = 0;
      if (!MPI_T)
	{
	  int size = N * M;
	  MPI_Type_contiguous ( size, MPI_typetrait<T>::MPIType(), &MPI_T);
	  MPI_Type_commit ( &MPI_T );
	}
      return MPI_T;
    }
  };
  

}
#endif

#endif
