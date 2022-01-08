#ifndef FILE_VECTOR_EXPR
#define FILE_VECTOR_EXPR

/**************************************************************************/
/* File:   vector.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jan. 02                                                    */
/**************************************************************************/

namespace ngbla
{


  template <int S, class T> class Vec;
  template <int S, typename T> class FlatVec;
  template <class T> class SysVector;
  template <class T = double> class Vector;
  template <class T = double> class SliceVector;
  template <int DIST, typename T> class FixSliceVector;
  template <int S, int DIST, typename T> class FlatSliceVec;


  /**
     A simple vector.
     Has size and generic data-pointer. 
     No memory allocation/deallocation. User must provide memory.
  */
  template <typename T>
  class FlatVector : public CMCPMatExpr<FlatVector<T> > 
  {
  protected:
    /// vector size
    size_t size;
    /// the data
    T * __restrict data;
  public:
    /// element type
    typedef T TELEM;
    // typedef T& TREF;
    /// scalar of element type
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { IS_LINEAR = 1 };

    /// default constructor does nothing
    FlatVector () = default;
    ~FlatVector () = default;
    /// set size and mem
    INLINE FlatVector (size_t as, T * adata) : size(as), data(adata)
    { ; }

    ///
    FlatVector (const FlatVector &) = default;
    
    /// set size and mem
    INLINE FlatVector (size_t as, void * adata) : size(as), data(static_cast<TELEM*> (adata)) 
    {  ; }

    /// put FlatVector over fixed size vector
    template <int S>
    INLINE FlatVector (const Vec<S,TSCAL> & v)
      : size(v.Size()), data(const_cast<T*>(&v(0)))
    { ; }

    template <int S>
    INLINE FlatVector (const FlatVec<S,TSCAL> & v)
      : size(v.Size()), data(const_cast<T*>(&v(0)))
    { ; }

    /// allocate FlatVector on local heap
    INLINE FlatVector (size_t as, LocalHeap & lh) 
      : size(as), data(lh.Alloc<T> (as)) 
    { ; }

    /// put FlatVector over systemvector
    INLINE FlatVector (const SysVector<TSCAL> & sv)
      : size(sv.Size()*sv.BlockSize() / mat_traits<T>::VDIM), 
	data (sv(0))
    { ; }



    /// allocate and compute 
    template<typename TB>
    INLINE FlatVector (const LocalHeapExpr<TB> & m2) 
    {
      size = m2.A().Height();

      LocalHeap & lh = m2.GetLocalHeap();
      data = lh.Alloc<T> (size);

      // does not copy FlatVectors and don't know why
      // CMCPMatExpr<FlatVector<T> >::operator= (m2.A()); 

      for (size_t j = 0; j < size; j++) data[j] = m2.A()(j);
    }

    /// assign memory for vector on local heap
    INLINE void AssignMemory (size_t as, LocalHeap & lh) 
    {
      size = as;
      data = lh.Alloc<T>(size);
    }

    /// assign memory for vector
    void AssignMemory (size_t as, T * mem) 
    {
      size = as;
      data = mem;
    }

    /// copy vector. sizes must match
    INLINE const FlatVector & operator= (const FlatVector & v) const
    {
      NETGEN_CHECK_RANGE(v.Size(),0,Size()+1);
      for (auto i : ngstd::Range(size))
	data[i] = v(i);
      return *this;
    }

    template <int D, typename TSCAL2>
    INLINE const FlatVector & operator= (const Vec<D,TSCAL2> & v) const
    {
      NETGEN_CHECK_RANGE(D,0,Size()+1);
      for (int i = 0; i < D; i++)
	data[i] = v(i);
      return *this;
    }
    
    /// evaluate matrix expression
    template<typename TB>
    INLINE const FlatVector & operator= (const Expr<TB> & v) const
    {
      return CMCPMatExpr<FlatVector>::operator= (v);
    }

    /// assign constant value
    INLINE const FlatVector & operator= (TSCAL scal) const
    {
      for (auto i : Range())
	data[i] = scal; 
      return *this;
    }

    INLINE const FlatVector & operator= (initializer_list<T> list) const
    {
      NETGEN_CHECK_RANGE(list.size(),0,Size()+1);
      size_t cnt = 0;
      for (auto val : list)
        data[cnt++] = val;
      return *this;
    }

    /*  // prevents bla pattern matching
    template<typename TB>
    INLINE const FlatVector & operator+= (const Expr<TB> & v) const
    {
      if (TB::IS_LINEAR)
        for (auto i : Range())
          data[i] += v.Spec()(i);
      else
        for (auto i : Range())        
          data[i] += v.Spec()(i,0);
      return *this;
    }
    */

    /// constant element access
    INLINE TELEM & operator() (size_t i) const
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }
 
    INLINE RowsArrayExpr<FlatVector> operator() (FlatArray<int> rows) const
    { 
      return RowsArrayExpr<FlatVector> (*this, rows);
    }
    
    /// element access. index j is ignored
    INLINE TELEM & operator() (size_t i, size_t j) const
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }

    /// constant element access
    INLINE TELEM & operator[] (size_t i) const
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }

    // shape functions had a problem with icc v9.1
    // const CArray<T> Addr(int i) const { return CArray<T> (data+i*dist); }
    INLINE T * Addr(size_t i) const { return data+i; }
    /*
    const CArray<T> operator+(int i) const
    { return CArray<T> (data+i*dist); }
    */

    T * operator+(size_t i) const { return data+i; }

    /// sub-vector of size next-first, starting at first
    INLINE FlatVector Range (size_t first, size_t next) const
    { return FlatVector (next-first, data+first); }

    /// sub-vector given by range
    INLINE FlatVector Range (T_Range<size_t> range) const
    { return Range (range.First(), range.Next()); }
    

    /// vector size
    INLINE size_t Size () const { return size; }

    /// vector is matrix of height size
    INLINE size_t Height () const { return size; }

    /// vector is matrix of with 1
    INLINE constexpr size_t Width () const { return 1; }
    
    INLINE T_Range<size_t> Range () const
    { return T_Range<size_t> (0, size); }

    /// take a slice of the vector. Take elements first+i * dist. 
    INLINE const SliceVector<T> Slice (size_t first, size_t dist2) const
    {
      return SliceVector<T> (size/dist2, dist2, data+first);
    }

    INLINE FlatMatrix<T> AsMatrix (size_t h, size_t w)
    {
      return FlatMatrix<T> (h,w, data);
    }

    INLINE T * Data () const { return data; }

    // new for SysVectors:
    typedef FlatVector<T> TV_COL;
    typedef double TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };

    class Iterator
    {
      FlatVector vec;
      size_t ind;
    public:
      INLINE Iterator (FlatVector avec, size_t ai) : vec(avec), ind(ai) { ; }
      INLINE Iterator operator++ (int) { return Iterator(vec, ind++); }
      INLINE Iterator operator++ () { return Iterator(vec, ++ind); }
      INLINE TELEM operator*() const { return vec[ind]; }
      INLINE TELEM & operator*() { return vec[ind]; }
      INLINE bool operator != (Iterator d2) { return ind != d2.ind; }
      INLINE bool operator == (Iterator d2) { return ind == d2.ind; }
    };
    
    Iterator begin() const { return Iterator (*this, 0); }
    Iterator end() const { return Iterator (*this, size); }
  };







  template <int S, typename T>
  class FlatVector<Vec<S, T>> : public CMCPMatExpr<FlatVector<Vec<S, T>> > 
  {
  protected:
    /// vector size
    size_t size;
    /// the data
    T *  __restrict data;
  public:
    /// element type
    typedef Vec<S,T> TELEM;
    typedef FlatVec<S,T> TREF;
    /// scalar of element type
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /// default constructor does nothing
    FlatVector () { ; }  
    /// copy pointer
    // FlatVector (const FlatVector & v2) : size(v2.size), data(v2.data) { ; }
    INLINE FlatVector (const FlatVector<Vec<S,T>> & fv2) : size(fv2.Size()), data((T*)fv2.Data()) { ; }
    /// set size and mem
    FlatVector (size_t as, T * adata) : size(as), data(adata) { ; }
    
    /// set size and mem
    FlatVector (size_t as, void * adata) : size(as), data(static_cast<T*> (adata)) { ; }
    
    /*
    /// put FlatVector over fixed size vector
    template <int S>
    FlatVector (const Vec<S,TSCAL> & v)
      : s(v.Size()), data(const_cast<T*>(&v(0)))
    { ; }
    */

    /// allocate FlatVector on local heap
    FlatVector (size_t as, LocalHeap & lh) 
      : size(as), data((T*)lh.Alloc(size*S*sizeof(T))) { ; }

    /*
    /// put FlatVector over systemvector
    FlatVector (const SysVector<TSCAL> & sv)
      : s(sv.Size()*sv.BlockSize() / mat_traits<T>::VDIM), 
	data (sv(0))
    {
      ;
    }
    */

    /// assign memory for vector on local heap
    void AssignMemory (size_t as, LocalHeap & lh) 
    {
      size = as;
      data = (T*)lh.Alloc(size*S*sizeof(T));
    }

    /// assign memory for vector
    void AssignMemory (size_t as, T * mem) 
    {
      size = as;
      data = mem;
    }

    /// copy vector. sizes must match
    const FlatVector & operator= (const FlatVector & v) const
    {
      for (size_t i = 0; i < size; i++)
	(*this)(i) = v(i);
      return *this;
    }

    /// evaluate matrix expression
    template<typename TB>
    const FlatVector & operator= (const Expr<TB> & v) const
    {
      return CMCPMatExpr<FlatVector>::operator= (v);
    }

    /// assign constant value
    const FlatVector & operator= (TSCAL scal) const
    {
      for (size_t i = 0; i < size; i++)
	(*this)(i) = scal; 
      return *this;
    }

    template<typename TB>
    ALWAYS_INLINE const FlatVector & operator+= (const Expr<TB> & v) const
    {
      if (TB::IS_LINEAR)
        for (auto i : ngstd::Range(size))
          (*this)(i) += v.Spec()(i);
      else
	for (size_t i = 0; i < size; i++)
	  (*this)(i) += v.Spec()(i,0);
      return *this;
    }

    /// constant element access
    const FlatVec<S,T> operator() (size_t i) const
    {
      return FlatVec<S,T> (data+i*S); 
    }

    /// element access. index j is ignored
    const FlatVec<S,T> operator() (size_t i, size_t j) const
    {
      return FlatVec<S,T> (data+i*S); 
    }

    /// constant element access
    const FlatVec<S,T> operator[] (size_t i) const
    {
      return FlatVec<S,T> (data+i*S); 
    }

    RowsArrayExpr<FlatVector> operator() (FlatArray<int> rows) const
    { 
      return RowsArrayExpr<FlatVector> (*this, rows);
    }

    // shape functions had a problem with icc v9.1
    Vec<S,T> * Addr(size_t i) const
    {
      return static_cast<Vec<S,T>*> ((void*) (data+i*S)); 
    }

    /*
    const CArray<Vec<S,T> > Addr(int i) const
    {
      return CArray<Vec<S,T> > (static_cast<Vec<S,T>*> ((void*) (data+i*S))); 
    }
    */
    /*
    const CArray<T> Addr(int i) const
    {
      return CArray<T> (data+i*S); 
    }
    */

    /// sub-vector of size next-first, starting at first
    /* const */ FlatVector<Vec<S,T> > Range (size_t first, size_t next) const
    { return FlatVector<Vec<S,T> > (next-first, data+S*first); }

    /// sub-vector given by range
    /* const */ FlatVector<Vec<S,T> > Range (IntRange range) const
    { return FlatVector<Vec<S,T> > (range.Next()-range.First(), data+S*range.First()); }

    /// vector size
    INLINE size_t Size () const { return size; }

    /// vector is matrix of height size
    INLINE size_t Height () const { return size; }

    /// vector is matrix of with 1
    INLINE constexpr size_t Width () const { return 1; }


    INLINE SliceVector<T> Comp (size_t comp) const
    {
      return SliceVector<T> (size, S, data+comp);
    }

    /*
    /// take a slice of the vector. Take elements first+i * dist. 
    SliceVector<T> Slice (int first, int dist)
    {
      return SliceVector<T> (s/dist, dist, data+first);
    }

    /// take a slice of the vector. Take elements first+i * dist. 
    const SliceVector<T> Slice (int first, int dist) const
    {
      return SliceVector<T> (s/dist, dist, data+first);
    }
    */

    /// access to data
    const void * Data () const { return static_cast<const void*>(data); }
    /// access to data
    void * Data () { return static_cast<void*>(data); }

    // new for SysVectors:
    // typedef FlatVector<T> TV_COL;
    // typedef double TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
  };





  

  /**
     A Vector class with memory allocation/deallocation
  */
  template <typename T>
  class Vector : public FlatVector<T>
  {
    // using FlatVector<T>::data;
    // using FlatVector<T>::size;
  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /// default constructor
    Vector () : FlatVector<T> (0, (T*)nullptr) { ; }

    /// allocate vector
    explicit Vector (size_t as) : FlatVector<T> (as, new T[as]) { ; }


    /// allocate and copy matrix  
    Vector (const Vector & v2) 
      : FlatVector<T> (v2.Size(), new T[v2.Size()]) 
    {
      FlatVector<T>::operator= (v2);
    }

    Vector (Vector && v2)
      : FlatVector<T> (v2.size, v2.data)
    { v2.data = nullptr; v2.size = 0; } 
    
    /// allocate and compute 
    template<typename TB>
    Vector (const Expr<TB> & v2) 
      : FlatVector<T> (v2.Height(), new T[v2.Height()]) 
    {
      FlatVector<T>::operator= (v2);
    }


    Vector (initializer_list<T> list) 
      : FlatVector<T> (list.size(), new T[list.size()])
    {
      /*
      int cnt = 0;
      for (auto i = list.begin(); i < list.end(); i++, cnt++)
        data[cnt] = *i;
      */
      size_t cnt = 0;
      for (auto val : list)
        (*this)[cnt++] = val;
    }


    /// deallocate vector
    ~Vector() { delete [] this->data; }

    /// set vector to constant values
    Vector & operator= (TSCAL scal)
    {
      FlatVector<T>::operator= (scal);
      return *this;
    }

    /// set vector size
    void SetSize(size_t as)
    {
      if (this->size == as) return;
      delete [] this->data;
      this->size = as;
      if (as != 0)
        this->data = new T[this->size];
      else
        this->data = nullptr;
    }

    /// evaluate matrix expression
    template<typename TB>
    Vector & operator= (const Expr<TB> & v)
    {
      MatExpr<FlatVector<T> >::operator= (v);
      return *this;
    }

    Vector & operator= (const Vector & v2)
    {
      FlatVector<T>::operator= (static_cast<FlatVector<T> >(v2));
      // MatExpr<FlatVector<T> >::operator= (v2);  // does not work, we don't know why
      return *this;
    }

    Vector & operator= (Vector && v2)
    {
      this->size = v2.size;
      Swap (this->data, v2.data);
      return *this;
    }

    Vector & operator= (initializer_list<T> list) 
    {
      SetSize (list.size());
      size_t cnt = 0;
      for (auto val : list)
        (*this)[cnt++] = val;
      return *this;
    }

  };





  template <int S, typename T>
  class Vector<Vec<S,T> > : public FlatVector<Vec<S,T> >
  {
  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /// default constructor
    Vector () : FlatVector<Vec<S,T> > (0, (T*)0) { ; }

    /// allocate vector
    explicit Vector (size_t as) : FlatVector<Vec<S,T> > (as, new T[as*S]) { ; }


    /// allocate and copy matrix  
    Vector (const Vector & v2) 
      : FlatVector<Vec<S,T> > (v2.Size(), new T[S*v2.Size()]) 
    {
      FlatVector<Vec<S,T> >::operator= (v2);
    }
    
    /// allocate and compute 
    template<typename TB>
    Vector (const Expr<TB> & v2) 
      : FlatVector<Vec<S,T> > (v2.Height(), new T[S*v2.Height()]) 
    {
      FlatVector<Vec<S,T> >::operator= (v2);
    }

    /// deallocate vector
    ~Vector() { delete [] this->data; }

    /// set vector to constant values
    Vector & operator= (TSCAL scal)
    {
      FlatVector<Vec<S,T> >::operator= (scal);
      return *this;
    }

    /// set vector size
    void SetSize(size_t as)
    {
      if (this->size == as) return;
      delete [] this->data;
      this->size = as;
      this->data = new T[S*this->size];
    }

    /// evaluate matrix expression
    template<typename TB>
    Vector & operator= (const Expr<TB> & v)
    {
      MatExpr<FlatVector<Vec<S,T> > >::operator= (v);
      return *this;
    }

    Vector & operator= (const Vector & v2)
    {
      FlatVector<Vec<S,T> >::operator= (static_cast<FlatVector<Vec<S,T> > >(v2));
      // MatExpr<FlatVector<Vec<S,T> > >::operator= (v2);  // does not work, we don't know why
      return *this;
    }

  };







  /**
     A Vector class with memory allocation/deallocation.  At
     compile-time, a certain amount of vector entries is defined. If the
     dynamic size fits into this memory, the vector is allocated in this memory.
     Otherwise, dynamic memory is allocated.
  */
  template <int S, typename T = double>
  class VectorMem : public FlatVector<T>
  {
    T mem[S];
  public:
    /// the scalar type
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /** allocate vector. 
	If the dynamic size fits into the static size, use static memory. Otherwise use dynamic alloation
    */
    explicit INLINE VectorMem (size_t as) : FlatVector<T> (as, (as <= S) ? 
                                                           static_cast<void*> (&mem[0]) : 
                                                           static_cast<void*> (new T[as])) { ; }

    /// deallocates dynamic memory
    INLINE ~VectorMem() { if (this->Size() > S) delete [] this->data; }

    /// assigns constant value
    INLINE VectorMem & operator= (TSCAL scal)
    {
      FlatVector<T>::operator= (scal);
      return *this;
    }

    /// evaluates matrix expression
    template<typename TB>
    INLINE VectorMem & operator= (const Expr<TB> & v)
    {
      MatExpr<FlatVector<T> >::operator= (v);
      return *this;
    }
  };















  // A system vector (not completely functional)
  template <typename T = double>
  class FlatSysVector : public MatExpr<FlatSysVector<T> >
  {
  protected:
    size_t s;
    size_t blocksize;
    T *  __restrict data;
  public:
    typedef FlatVector<T> TELEM;
    typedef typename mat_traits<T>::TSCAL TSCAL;
    
    INLINE FlatSysVector (size_t as, size_t bs, T * adata) 
      : s(as), blocksize(bs), data(adata) { ; }
  
    INLINE FlatSysVector (size_t as, size_t bs, LocalHeap & lh) 
      : s(as), blocksize(bs), data (new (lh) T[as*bs])
    { ; }
  
    /*
      FlatSysVector (const FlatSysVector<TSCAL> & sv)
      : s(sv.Size()),
      blocksize(sv.BlockSize()),
      data (sv(0))
      {
      ;
      }
    */
  
    INLINE FlatSysVector & operator= (const FlatSysVector & v)
    {
      for (size_t i = 0; i < this->s * this->blocksize; i++)
	data[i] = v.data[i];
      return *this;
    }
  
    template<typename TB>
    INLINE FlatSysVector & operator= (const Expr<TB> & v)
    {
      return MatExpr<FlatSysVector>::operator= (v);
    }
  
    INLINE FlatSysVector & operator= (TSCAL s)
    {
      for (size_t i = 0; i < this->s*this->blocksize; i++)
	data[i] = s;
      return *this;
    }

  

    INLINE TELEM operator() (size_t i) 
    {
      return FlatVector<T> (blocksize, &data[i*blocksize]); 
    }

    INLINE const TELEM operator() (size_t i) const
    {
      return FlatVector<T> (blocksize, &data[i*blocksize]); 
    }

    INLINE const TELEM operator() (size_t i, size_t j) const
    {
      return FlatVector<T> (blocksize, &data[i*blocksize]); 
    }

    INLINE SliceVector<T> Comp (size_t comp) const
    {
      return SliceVector<T> (s, blocksize, data+comp);
    }

    INLINE FlatSysVector<T> Range(size_t first, size_t last)
    { return FlatSysVector<T> (last-first, blocksize, data+(first*blocksize)); }
    INLINE /* const */ FlatSysVector<T> Range(size_t first, size_t last) const
    { return FlatSysVector<T> (last-first, blocksize, data+(first*blocksize)); }

    INLINE size_t Size () const { return s; }
    INLINE size_t Height () const { return s; }
    INLINE constexpr size_t Width () const { return 1; }
  };



  /**
     A vector of fixed size.
     Useful as entry type in system vectors.
  */
  template <int S, typename T = double>
  class Vec : public MatExpr<Vec<S,T> > // , protected BaseVec
  {
    /// the values
    // T data[S];
    HTArray<S,T> data;
  public:
    /// type of the elements
    typedef T TELEM;
    /// is the element double or complex ?
    typedef typename mat_traits<T>::TSCAL TSCAL;
    /// a vec is a S times 1 matrix, the according colume vector
    typedef Vec<S, typename mat_traits<T>::TV_COL> TV_COL;
    /// a vec is a S times 1 matrix, the according row vector
    typedef Vec<1, typename mat_traits<T>::TV_ROW> TV_ROW;

    enum { SIZE = S };
    /// height of matrix
    enum { HEIGHT = S };
    /// with of matrix
    enum { WIDTH  = 1 };

    /// constructor, no initialization
    INLINE Vec () { ; }
    /// copy vector
    INLINE Vec (const Vec & v) : MatExpr<Vec> ()
    {
      for (size_t i = 0; i < S; i++)
	data[i] = v.data[i];
    }

    /// initialize with values
    INLINE Vec (const TELEM & scal)
    {
      for (size_t i = 0; i < S; i++)
	data[i] = scal;
    }

    /// initialize with expression
    template<typename TB>
    INLINE Vec (const Expr<TB> & v)
    {
      (*this) = v;
    }

    template <int D>
    INLINE Vec(FlatSliceVec<S,D,T> fsv)
    {
      for (int i = 0; i < S; i++)
        data[i] = fsv(i);
    }

    // Helper function for variadic constructor
    template <int I, class... T2>
    void Set(const TELEM &v, T2... rest)
    {
      data[I] = v;
      Set<I+1>(rest...);
    }

    template <int I>
    void Set(const TELEM &v)
    {
      data[I] = v;
    }

    /*
    template <class... T2>
    Vec(const TELEM &v, T2... rest) {
      static_assert(S==1+sizeof...(rest),"Vec<S> ctor with wrong number of arguments called");
      Set<0>(v, rest...);
    }
    */

    template <class... T2,
              typename enable_if<S==1+sizeof...(T2),int>::type=0>
    Vec(const TELEM &v, T2... rest) {
      // static_assert(S==1+sizeof...(rest),"Vec<S> ctor with wrong number of arguments called");
      Set<0>(v, rest...);
    }
  
    /// copy vector
    INLINE Vec & operator= (const Vec & v)
    {
      for (size_t i = 0; i < S; i++)
	data[i] = v.data[i];
      return *this;
    }
  

    /// assign scalar value
    INLINE Vec & operator= (const TELEM & scal)
    {
      for (size_t i = 0; i < S; i++)
	data[i] = scal;
      return *this;
    }

    /// assign expression
    template<typename TB>
    INLINE Vec & operator= (const Expr<TB> & v)
    {
      for (size_t i = 0; i < S; i++)
	data[i] = v.Spec()(i,0);
      return *this;
    }

    /// access vector
    INLINE TELEM & operator() (size_t i) 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }

    /// access vector
    INLINE const TELEM & operator() (size_t i) const 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }


    /// access vector
    INLINE TELEM & operator[] (size_t i) 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }

    /// access vector
    INLINE const TELEM & operator[] (size_t i) const 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }

    /// access vector
    INLINE TELEM & operator() (size_t i, size_t j) 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }

    /// access vector
    INLINE const TELEM & operator() (size_t i, size_t j) const 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }

    /// vector size
    INLINE constexpr size_t Size () const { return S; }
    /// corresponding matrix height
    INLINE constexpr size_t Height () const { return S; }
    /// corresponding matrix with
    INLINE constexpr size_t Width () const { return 1; }
    INLINE T* Data() { return data.Ptr(); }

    INLINE /* const */ FlatVector<const T> Range(size_t first, size_t next) const
    { return FlatVector<const T> (next-first, data+first); }

    INLINE /* const */ FlatVector<T> Range(size_t first, size_t next) 
    { return FlatVector<T> (next-first, data+first); }
  };



  template <typename T>
  class Vec<0,T>  : public MatExpr<Vec<0,T> > 
  {
  public:
    INLINE Vec () { ; }
    INLINE Vec (const Vec &d) { ; }
    INLINE Vec (T d) { ; }
    template<typename TB>
    INLINE Vec (const Expr<TB> & v) {;}
    INLINE constexpr size_t Size() const { return 0; }
    INLINE constexpr size_t Height() const { return 0; }
    INLINE constexpr size_t Width() const { return 1; }
    template<typename TB>
    INLINE Vec & operator= (const Expr<TB> & v) { return *this;}
    INLINE Vec & operator= (const T & /* scal */) { return *this; } 
    
    INLINE T & operator[] (int i) const  { return *(T*)(void*)(this); }
    INLINE T & operator() (int i) const  { return *(T*)(void*)(this); }
    INLINE T & operator() (int i, int j) const  { return *(T*)(void*)(this); }
  };
  /*
  template <typename T>
  class Vec<0,T>  : public MatExpr<Vec<0,T> > 
  {
    T dummy;
  public:
    INLINE Vec () { ; }
    INLINE Vec (T d) { ; }
    template<typename TB>
    INLINE Vec (const Expr<TB> & v) {;}

    template<typename TB>
    INLINE Vec & operator= (const Expr<TB> & v) { return *this;}

    INLINE T operator[] (int i) const  { return 0.0; }
    INLINE T operator() (int i) const  { return 0.0; }
    INLINE T operator() (int i, int j) const  { return 0.0; }

    INLINE T & operator[] (int i)  { return dummy; }
    INLINE T & operator() (int i)  { return dummy; }
    INLINE T & operator() (int i, int j)  { return dummy; }
  };
  */


  template <int S, typename T>
  class mat_traits<Vec<S,T>>
  {
  public:
    /// matrix element
    typedef T TELEM;
    /// field of matrix element
    typedef typename mat_traits<T>::TSCAL TSCAL;
    /// type of column vector
    typedef Vec<S, typename mat_traits<T>::TV_COL> TV_COL;
    /// a vec is a S times 1 matrix, the according row vector
    typedef Vec<1, typename mat_traits<T>::TV_ROW> TV_ROW;
    /// matrix height
    enum { HEIGHT = S };
    /// matrix with
    enum { WIDTH  = 1  };
    ///
    enum { IS_COMPLEX = mat_traits<TSCAL>::IS_COMPLEX };
  };


  
  template <typename T> struct ConstVecSize { static constexpr int VSIZE = -1; }; 
  template <int S, typename T> struct ConstVecSize<Vec<S,T>> { static constexpr int VSIZE = S; };
  template <int S, typename T> struct ConstVecSize<FlatVec<S,T>> { static constexpr int VSIZE = S; };
  template <int S, int D, typename T> struct ConstVecSize<FlatSliceVec<S,D,T>> { static constexpr int VSIZE = S; };
  template <typename T>
  constexpr auto ConstVectorSize() { return ConstVecSize<T>::VSIZE; }
      
  /// cross product of 3-vectors
  template <typename TA, typename TB,
            // std::enable_if_t<ConstVectorSize<TA>() == 3, bool> = true,
            // std::enable_if_t<ConstVectorSize<TB>() == 3, bool> = true>
            std::enable_if_t<ConstVecSize<TA>::VSIZE == 3, bool> = true,
            std::enable_if_t<ConstVecSize<TB>::VSIZE == 3, bool> = true>
  INLINE auto Cross (const TA & a, const TB & b)
  {
    typedef decltype (a(0)*b(0)) T;
    return Vec<3,T>({ a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0) });
  }

  /*
  /// cross product
  template <typename S>
  INLINE Vec<3,S> Cross (const Vec<3,S> & a, const Vec<3,S> & b)
  {
    return Vec<3,S>({ a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0) });
  }
  */


  template <typename S>
  INLINE Vec<1,S> Cross (const Vec<2,S> & a, const Vec<2,S> & b)
  {
    return Vec<1,S> ( { a(0) * b(1) - a(1) * b(0) } );
    /*
    Vec<1,S> prod;
    prod(0) = a(0) * b(1) - a(1) * b(0);
    return prod;
    */
  }

  template <typename S>
  INLINE Vec<0,S> Cross (const Vec<1,S> & a, const Vec<1,S> & b)
  {
    return Vec<0,S>();
  }
  
  /// output vector
  template<int S, typename T>
  inline ostream & operator<< (ostream & ost, const Vec<S,T> & v)
  {
    for (int i = 0; i < S; i++)
      ost << " " << setw(7) << v(i);
    return ost;
  }


  template<int S, typename TB>
  INLINE Vec<S> & operator+= (Vec<S> & v, const Expr<TB> & v2)
  {
    for (int i = 0; i < S; i++)
      v(i) += v2.Spec()(i,0);
    return v;
  }








  /**
     A pointer to a vector of fixed size.
  */
  template <int S, typename T = double>
  class FlatVec : public CMCPMatExpr<FlatVec<S,T> > 
  {
    /// the values
    T *  __restrict data;
  public:
    /// type of the elements
    typedef T TELEM;
    /// is the element double or complex ?
    typedef typename mat_traits<T>::TSCAL TSCAL;
    /// a vec is a S times 1 matrix, the according colume vector
    typedef Vec<S, typename mat_traits<T>::TV_COL> TV_COL;
    /// a vec is a S times 1 matrix, the according row vector
    typedef Vec<1, typename mat_traits<T>::TV_ROW> TV_ROW;

    enum { SIZE = S };
    /// height of matrix
    enum { HEIGHT = S };
    /// with of matrix
    enum { WIDTH  = 1 };

    /// constructor, no initialization
    INLINE FlatVec (T * adata) : data(adata) { ; }

    /// constructor, no initialization
    INLINE FlatVec (Vec<S,T> & v2) : data(&v2(0)) { ; }


    /// copy vector
    INLINE const FlatVec & operator= (const FlatVec & v) const
    {
      for (int i = 0; i < S; i++)
	data[i] = v.data[i];
      return *this;
    }

    /// assign scalar value
    INLINE const FlatVec & operator= (TSCAL scal) const
    {
      for (int i = 0; i < S; i++)
	data[i] = scal;
      return *this;
    }

    /// assign expression
    template<typename TB>
    INLINE const FlatVec & operator= (const Expr<TB> & v) const
    {
      for (int i = 0; i < S; i++)
	data[i] = v.Spec()(i,0);
      return *this;
    }

    template<typename TB>
    INLINE const FlatVec & operator+= (const Expr<TB> & v) const
    {
      for (int i = 0; i < S; i++)
	data[i] += v.Spec()(i,0);
      return *this;
    }

    /// access vector
    INLINE TELEM & operator() (int i) const 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }

    /// access vector
    INLINE TELEM & operator[] (int i) const 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }

    /// access vector
    INLINE TELEM & operator() (int i, int j) const 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i]; 
    }

    INLINE /* const */ FlatVector<T> Range(int first, int next) const
    { return FlatVector<T> (next-first, data+first); }

    /// vector size
    INLINE constexpr int Size () const { return S; }
    /// corresponding matrix height
    INLINE constexpr int Height () const { return S; }
    /// corresponding matrix with
    INLINE constexpr int Width () const { return 1; }
  };

  /// output vector.
  template<int S, typename T>
  inline ostream & operator<< (ostream & ost, const FlatVec<S,T> & v)
  {
    for (int i = 0; i < S; i++)
      ost << " " << setw(7) << v(i);
    return ost;
  }





  /**
     A pointer to a vector of fixed size.
  */
  template <int S, int D, typename T = double>
  class FlatSliceVec : public CMCPMatExpr<FlatSliceVec<S,D,T> > 
  {
    /// the values
    T *  __restrict data;
  public:
    /// type of the elements
    typedef T TELEM;
    /// is the element double or complex ?
    typedef typename mat_traits<T>::TSCAL TSCAL;
    /// a vec is a S times 1 matrix, the according colume vector
    typedef Vec<S, typename mat_traits<T>::TV_COL> TV_COL;
    /// a vec is a S times 1 matrix, the according row vector
    typedef Vec<1, typename mat_traits<T>::TV_ROW> TV_ROW;

    enum { SIZE = S };
    /// height of matrix
    enum { HEIGHT = S };
    /// with of matrix
    enum { WIDTH  = 1 };

    /// constructor
    INLINE FlatSliceVec (T * adata) : data(adata) { ; }

    /// copy vector
    INLINE auto operator= (const FlatSliceVec & v) const
    {
      for (int i = 0; i < S; i++)
	data[i*D] = v.data[i*D];
      return *this;
    }
    
    /// assign scalar value
    INLINE auto operator= (TSCAL scal) const
    {
      for (int i = 0; i < S; i++)
	data[i*D] = scal;
      return *this;
    }

    /// assign expression
    template<typename TB>
    INLINE auto operator= (const Expr<TB> & v) const
    {
      for (int i = 0; i < S; i++)
	data[i*D] = v.Spec()(i,0);
      return *this;
    }

    template<typename TB>
    INLINE auto operator+= (const Expr<TB> & v) const
    {
      for (int i = 0; i < S; i++)
	data[i*D] += v.Spec()(i,0);
      return *this;
    }

    operator Vec<S,T>() const
    {
      Vec<S,T> ret;
      for (int i = 0; i < S; i++)
        ret(i) = data[i*D];
      return ret;
    }
    
    /// access vector
    INLINE TELEM & operator() (int i) const 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*D]; 
    }

    /// access vector
    INLINE TELEM & operator[] (int i) const 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*D]; 
    }

    /// access vector
    INLINE TELEM & operator() (int i, int j) const 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*D]; 
    }

    // INLINE /* const */ FlatVector<T> Range(int first, int next) const
    // { return FlatVector<T> (next-first, data+first); }

    /// vector size
    INLINE constexpr int Size () const { return S; }
    /// corresponding matrix height
    INLINE constexpr int Height () const { return S; }
    /// corresponding matrix with
    INLINE constexpr int Width () const { return 1; }
  };

  /// output vector.
  template<int S, int D, typename T>
  inline ostream & operator<< (ostream & ost, const FlatSliceVec<S,D,T> & v)
  {
    for (int i = 0; i < S; i++)
      ost << " " << setw(7) << v(i*D);
    return ost;
  }







  
  // Template helper for the SliceVector from Vec constructor
  // to prevent range checks from triggering on call to v(0) for empty v
  template <typename T, int D>
  struct SliceVecFromVecHelper
  {
    static T *ptr(Vec<D,T> & v) { return &v(0); }
  };

  template <typename T>
  struct SliceVecFromVecHelper<T, 0>
  {
    static T *ptr(Vec<0,T> &) { return nullptr; }
  };

  /**
     A vector with non-linear data access.
     Has size and generic data-pointer. 
     No memory allocation/deallocation. User must provide memory.
  */
  template <typename T>
  class SliceVector : public CMCPMatExpr<SliceVector<T> > 
  {
  protected:
    /// vector size
    size_t s;
    /// distance between entries
    size_t dist;
    /// the data
    T *  __restrict data;
  public:
    /// the entry type
    typedef T TELEM;
    /// the scalar type of the vector
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /// element access is not linear
    enum { IS_LINEAR = 0 };

    /// set size, distance and memory
    INLINE SliceVector (size_t as, size_t ad, T * adata) 
      : s(as), dist(ad), data(adata) { ; }
    
    /// SV from FlatVector
    INLINE SliceVector (FlatVector<T> fv)
      : s(fv.Size()), dist(1), data((T*)fv.Data()) { ; }

    /// SV from Vec
    template <int D>
    INLINE SliceVector (Vec<D,T> & v)
      : s(D), dist(1), data(SliceVecFromVecHelper<T, D>::ptr(v)) { ; }

    ///
    template <int D>
    INLINE SliceVector (FixSliceVector<D,T> fsv)
      : s(fsv.Size()), dist(D), data(&fsv(0)) { ; }


    /// evaluates matrix expression
    template<typename TB>
    INLINE const SliceVector & operator= (const Expr<TB> & v) const
    {
      return CMCPMatExpr<SliceVector>::operator= (v);
    }

    /// assigns constant value
    INLINE const SliceVector & operator= (TSCAL scal) const
    {
      for (size_t i = 0; i < s; i++)
        data[i*size_t(dist)] = scal; 
      return *this;
    }

    /// copies contents of vector
    INLINE const SliceVector & operator= (const SliceVector & v2) const
    {
      for (size_t i = 0; i < s; i++)
	data[i*size_t(dist)] = v2(i);
      return *this;
    }


    template<typename TB>
    INLINE const SliceVector & operator+= (const Expr<TB> & v) const
    {
      if (TB::IS_LINEAR)
	for (size_t i = 0; i < s; i++)
	  data[i*size_t(dist)] += v.Spec()(i);
      else
	for (size_t i = 0; i < s; i++)
	  data[i*size_t(dist)] += v.Spec()(i,0);
      return *this;
    }



    /// access element
    INLINE TELEM & operator() (size_t i) 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*size_t(dist)]; 
    }

    /// access element
    INLINE TELEM & operator() (size_t i) const
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*size_t(dist)]; 
    }

    /// access element, index j is unused
    INLINE TELEM & operator() (size_t i, size_t j) const
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*size_t(dist)]; 
    }

    /// access element, index j is unused
    INLINE TELEM & operator() (size_t i, size_t j) 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*size_t(dist)]; 
    }

    /// access element
    INLINE TELEM & operator[] (size_t i) 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*size_t(dist)]; 
    }

    /// access element
    INLINE TELEM & operator[] (size_t i) const
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*dist]; 
    }

    
    INLINE SliceVector<T> operator+(int i) const { return SliceVector<T> (s-i, dist, data+i*size_t(dist)); }

    INLINE TELEM * Addr (size_t i) const
    {
      return data+i*size_t(dist);
    }

    /// vector size
    INLINE size_t Size () const { return s; }
    INLINE size_t Dist () const { return dist; }
    /// vector is a matrix of height size
    INLINE size_t Height () const { return s; }
    /// vector is a matrix of width 1
    INLINE constexpr size_t Width () const { return 1; }

    INLINE T * Data () const { return data; }

    INLINE const SliceVector<T> Range (size_t first, size_t next) const
    {
      return SliceVector<T> (next-first, dist, data+first*dist);
    }

    INLINE const SliceVector<T> Range (IntRange range) const
    {
      return Range (range.First(), range.Next());
    }

    INLINE const SliceVector<T> Slice (size_t first, size_t adist) const
    {
      return SliceVector<T> (s/adist, dist*adist, data+first*dist);
    }

    class Iterator
    {
      SliceVector vec;
      size_t ind;
    public:
      INLINE Iterator (SliceVector avec, size_t ai) : vec(avec), ind(ai) { ; }
      INLINE Iterator operator++ (int) { return Iterator(vec, ind++); }
      INLINE Iterator operator++ () { return Iterator(vec, ++ind); }
      INLINE TELEM operator*() const { return vec[ind]; }
      INLINE TELEM & operator*() { return vec[ind]; }
      INLINE bool operator != (Iterator d2) { return ind != d2.ind; }
      INLINE bool operator == (Iterator d2) { return ind == d2.ind; }
    };

    Iterator begin() const { return Iterator (*this, 0); }
    Iterator end() const { return Iterator (*this, this->Size()); }
  };

  
#ifdef NETGEN_ENABLE_CHECK_RANGE
  // Record with and height for classes that usually have no such information
  class DummySize {
    size_t height;
    size_t width;
  public:
    size_t Height() const { return height; }
    size_t Width() const { return width; }
    DummySize( size_t aheight, size_t awidth=1 ) :
      height(aheight), width(awidth) {;}
  };
#else 
  class DummySize {
  public:
    DummySize( size_t aheight, size_t awidth=1 ) {}
  protected:
    static INLINE size_t Height() { return 0; }
    static INLINE size_t Width() { return 0; }
  };
#endif



  
  template <class T = double>
  class BareVector : public CMCPMatExpr<BareVector<T> >, DummySize
  {
    T * __restrict data;
  public:
#ifdef NETGEN_ENABLE_CHECK_RANGE
    using DummySize::Width;
    using DummySize::Height;
#endif
    BareVector(T * _data) : DummySize(0,0), data(_data) { ; }
    BareVector(FlatVector<T> vec) : DummySize( vec.Size() ), data(&vec(0)) { ; }

    template <int D, typename TSCAL2>
    INLINE const BareVector & operator= (const Vec<D,TSCAL2> & v) const
    {
      for (int i = 0; i < D; i++)
	data[i] = v(i);
      return *this;
    }
    
    FlatVector<T> AddSize(size_t size) const
    {
      NETGEN_CHECK_RANGE(size, Height(), Height()+1);
      return FlatVector<T> (size, data);
    }
    
    T & operator() (size_t i) const
    {
      NETGEN_CHECK_RANGE(i, 0, Height());
      return data[i];
    }
    T & operator() (size_t i, size_t j) const{ return (*this)(i); }
    T & operator[] (size_t i) const { return (*this)(i); }

    /// sub-vector of size next-first, starting at first
    INLINE auto Range (size_t first, size_t next) const
    {
      NETGEN_CHECK_RANGE(first, 0, Height());
      NETGEN_CHECK_RANGE(next, 0, Height()+1);
      return FlatVector (next-first, data+first);
    }

    /// sub-vector given by range
    INLINE auto Range (T_Range<size_t> range) const
    {
      NETGEN_CHECK_RANGE(range.First(), 0, Height());
      NETGEN_CHECK_RANGE(range.Next(), 0, Height()+1);
      return Range (range.First(), range.Next());
    }
  };




  
  template <class T = double>
  class BareSliceVector : public CMCPMatExpr<BareSliceVector<T> >, DummySize
  {
    T * __restrict data;
    size_t dist;
#ifdef NETGEN_ENABLE_CHECK_RANGE
    BareSliceVector(T * _data, size_t _dist, DummySize dsize) : DummySize(dsize), data(_data), dist(_dist) { ; }
#else
    BareSliceVector(T * _data, size_t _dist) : DummySize(0,0), data(_data), dist(_dist) { ; }
#endif
  public:
    typedef typename mat_traits<T>::TSCAL TSCAL;
#ifdef NETGEN_ENABLE_CHECK_RANGE
    using DummySize::Width;
    using DummySize::Height;
#endif
    BareSliceVector(SliceVector<T> vec) : DummySize( vec.Size() ), data(vec.Data()), dist(vec.Dist()) { ; }
    template <int D>
    BareSliceVector(FixSliceVector<D,T> vec) : DummySize( vec.Size() ), data(vec.Data()), dist(D)  { ; }
    BareSliceVector(FlatVector<T> vec) : DummySize( vec.Size() ), data(vec.Data()), dist(1)  { ; }
    template <int D>
    BareSliceVector(Vec<D,T> & vec) :  DummySize( vec.Size() ), data(vec.Data()), dist(1) { ; }
    BareSliceVector(const BareSliceVector &) = default;
    BareSliceVector & operator= (const BareSliceVector&) = delete;
    size_t Dist () const { return dist; }
    T* Data() const { return data; }

    [[deprecated("Use Range(0,size) instead!")]]                
    SliceVector<T> AddSize(size_t size) const
    {
      NETGEN_CHECK_RANGE(size, Height(), Height()+1);
      return SliceVector<T> (size, dist, data);
    }
    
    T & operator() (size_t i) const
    {
      NETGEN_CHECK_RANGE(i, 0, Height());
      return data[i*dist];
    }
    T & operator() (size_t i, size_t j) const { return (*this)(i); }
    T & operator[] (size_t i) const { return (*this)(i); }
    BareSliceVector<T> operator+(size_t i) const
    {
#ifdef NETGEN_ENABLE_CHECK_RANGE
      return BareSliceVector<T> (data+i*dist, dist, Height()-i);
#else
      return BareSliceVector<T> (data+i*dist, dist);
#endif
    }
    T * Addr (size_t i) const { return data+i*dist; }
    SliceVector<T> Range (size_t first, size_t next) const
    {
      return SliceVector<T> (next-first, dist, data+first*dist);
    }
    SliceVector<T> Range (T_Range<size_t> range) const
    {
      return Range(range.First(), range.Next());
    }    
    BareSliceVector Slice (size_t first, size_t adist) const
    {
#ifdef NETGEN_ENABLE_CHECK_RANGE
      return BareSliceVector (data+first*dist, dist*adist, Height()/adist );
#else
      return BareSliceVector (data+first*dist, dist*adist);
#endif
    }

  };







  /**
     A vector with non-linear data access.
     Has size and generic data-pointer. 
     No memory allocation/deallocation. User must provide memory.
  */
  template <int DIST, typename T = double>
  class FixSliceVector : public CMCPMatExpr<FixSliceVector<DIST, T> > 
  {
  protected:
    /// vector size
    size_t s;
    /// the data
    T *  __restrict data;
  public:
    /// the entry type
    typedef T TELEM;
    /// the scalar type of the vector
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /// element access is not linear
    enum { IS_LINEAR = 0 };

    /// set size, distance and memory
    FixSliceVector (size_t as, T * adata) 
      : s(as), data(adata) { ; }

    /// evaluates matrix expression
    template<typename TB>
    FixSliceVector & operator= (const Expr<TB> & v)
    {
      return MatExpr<FixSliceVector>::operator= (v);
    }

    /// assigns constant value
    const FixSliceVector & operator= (TSCAL scal) const
    {
      for (size_t i = 0; i < s; i++)
	data[i*DIST] = scal; 
      return *this;
    }

    /// copies contents of vector
    INLINE const FixSliceVector & operator= (const FixSliceVector & v2) const
    {
      for (size_t i = 0; i < s; i++)
	data[i*DIST] = v2(i);
      return *this;
    }

    template<typename TB>
    INLINE const FixSliceVector & operator= (const Expr<TB> & v) const
    {
      if (TB::IS_LINEAR)
	for (size_t i = 0; i < s; i++)
	  data[i*DIST] = v.Spec()(i);
      else
	for (size_t i = 0; i < s; i++)
	  data[i*DIST] = v.Spec()(i,0);
      return *this;
    }


    template<typename TB>
    INLINE const FixSliceVector & operator+= (const Expr<TB> & v) const
    {
      if (TB::IS_LINEAR)
	for (size_t i = 0; i < s; i++)
	  data[i*DIST] += v.Spec()(i);
      else
	for (size_t i = 0; i < s; i++)
	  data[i*DIST] += v.Spec()(i,0);
      return *this;
    }



    /// access element
    TELEM & operator() (size_t i) 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*DIST]; 
    }

    /// access element
    TELEM & operator() (size_t i) const
    {
      NETGEN_CHECK_RANGE(i,0,s);
      return data[i*DIST]; 
    }

    /// access element, index j is unused
    TELEM & operator() (size_t i, size_t j) const
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*DIST]; 
    }

    /// access element, index j is unused
    TELEM & operator() (size_t i, size_t j) 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*DIST]; 
    }

    /// access element
    TELEM & operator[] (size_t i) 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*DIST]; 
    }

    /// access element
    TELEM & operator[] (size_t i) const
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*DIST]; 
    }

    /*
    TELEM * Addr (int i) const
    {
      return data+i*DIST;
    }
    */

    INLINE T * Data () const { return data; }

    /// vector size
    size_t Size () const { return s; }

    /// vector is a matrix of height size
    size_t Height () const { return s; }
    /// vector is a matrix of width 1
    size_t constexpr Width () const { return 1; }

    const FixSliceVector Range (size_t first, size_t next) const
    {
      return FixSliceVector (next-first, data+first*DIST);
    }

    const FixSliceVector Range (IntRange range) const
    {
      return Range (range.First(), range.Next());
    }


    const SliceVector<T> Slice (size_t first, size_t adist) const
    {
      return SliceVector<T> (s/adist, DIST*adist, data+first*DIST);
    }
  };










  template <class TV, class TSCAL> class Scalar2ElemVector
  {
  public:
    const FlatVector<TSCAL> & vec;
    Scalar2ElemVector (const FlatVector<TSCAL> & avec) : vec(avec) { ; }

    enum { H = mat_traits<TV>::HEIGHT };

    FlatVec<H,TSCAL> operator() (int i) const
    {
      return FlatVec<H,TSCAL> (&vec(i*H));
    }

  };
  

  template <class TSCAL> class Scalar2ElemVector<TSCAL,TSCAL>
  {
  public:
    const FlatVector<TSCAL> & vec;
    Scalar2ElemVector (const FlatVector<TSCAL> & avec) : vec(avec) { ; }


    const TSCAL & operator() (int i) const
    {
      return vec(i);
    }

    TSCAL & operator() (int i)
    {
      return vec(i);
    }

  };



  class SparseVector
  {
    size_t size;
    ClosedHashTable<size_t, double> data;
  public:
    SparseVector (size_t asize) : size(asize), data(10) { }
    size_t Size() const { return size; }
    double operator[] (size_t i) const { return data[i]; }
    double & operator[] (size_t i) { return data[i]; }
    double InnerProduct (FlatVector<> v2) const
    {
      double sum = 0;
      for (auto [i,v] : data)
        sum += v * v2(i);
      return sum;
    }
    friend ostream & operator<< (ostream & ost, const SparseVector & sv);
  };
  
  inline ostream & operator<< (ostream & ost, const SparseVector & sv) {
    for (auto [i,v] : sv.data)
      ost << i << ": " << v << ", ";
    return ost;
  }


  template <class T>
  class mat_traits<FlatVector<T> >
  {
  public:
    typedef T TELEM;
    typedef T TSCAL;
  };




  /*
  template <int S, typename T>
  INLINE auto operator* (double a, const Vec<S,T> & vec) 
    -> Vec<S, decltype(RemoveConst(a*vec(0)))>
  {
    typedef decltype(RemoveConst(a*vec(0))) TRES;
    Vec<S, TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a * vec(i);
    return res;
  }
  */
  
  template <int S, typename T>
  INLINE auto operator* (double a, const Vec<S,T> & vec)
  {
    typedef decltype(RemoveConst(a*vec(0))) TRES;
    Vec<S, TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a * vec(i);
    return res;
  }
  
  template <int S, typename T>
  INLINE auto operator* (Complex a, const Vec<S,T> & vec) 
  {
    typedef decltype(RemoveConst(a*vec(0))) TRES;
    Vec<S, TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a * vec(i);
    return res;
  }

  // all other cases ...
  template <int S, typename T,
            typename enable_if<!is_convertible_v<T,Complex>,int>::type=0>
  INLINE auto operator* (T a, const Vec<S,T> & vec) 
  {
    typedef decltype(RemoveConst(a*vec(0))) TRES;
    Vec<S, TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a * vec(i);
    return res;
  }
  
  template <int S, typename T>
  INLINE auto operator+ (const Vec<S,T> & a, const Vec<S,T> & b) 
  {
    typedef decltype(RemoveConst(a(0))) TRES;    
    Vec<S,TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a(i)+b(i);
    return res;
  }

  template <int S, typename T>
  INLINE auto operator- (const Vec<S,T> & a, const Vec<S,T> & b) 
  {
    typedef decltype(RemoveConst(a(0))) TRES;        
    Vec<S,TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a(i)-b(i);
    return res;
  }

  
  template <int S, typename T>
  INLINE auto operator* (double a, FlatVec<S,T> vec) 
    -> Vec<S, decltype(RemoveConst(a*vec(0)))>
  {
    typedef decltype(RemoveConst(a*vec(0))) TRES;
    Vec<S, TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a * vec(i);
    return res;
  }


  template <int S, typename T>
  INLINE auto operator* (Complex a, FlatVec<S,T> vec) 
    -> Vec<S, decltype(RemoveConst(a*vec(0)))>
  {
    typedef decltype(RemoveConst(a*vec(0))) TRES;
    Vec<S, TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a * vec(i);
    return res;
  }


  template <int S, typename T>
  INLINE auto operator+ (FlatVec<S,T> x, FlatVec<S,T> y) -> Vec<S,T>
  {
    Vec<S,T> tmp = x;
    tmp += y;
    return tmp;
  }

  template <int S, typename T>
  INLINE auto operator- (FlatVec<S,T> x, FlatVec<S,T> y) -> Vec<S,T>
  {
    Vec<S,T> tmp = x;
    tmp -= y;
    return tmp;
  }

  
  
  template <int DIM, typename SCAL, typename TANY>
  inline void AtomicAdd (Vec<DIM,SCAL> & x, TANY y)
  {
    for (int i = 0; i < DIM; i++)
      AtomicAdd (x(i), y(i));
  }
  
  template <int DIM, typename SCAL, typename TANY>
  inline void AtomicAdd (FlatVec<DIM,SCAL> x, TANY y)
  {
    for (int i = 0; i < DIM; i++)
      AtomicAdd (x(i), y(i));
  }
  

  
}

namespace ngstd
{
  template <int S, typename T>
  inline Archive & operator& (Archive & ar, ngbla::Vec<S,T> & v)
  {
    for (int i = 0; i < S; i++)
      ar & v(i);
    return ar;
  }
}


#ifdef PARALLEL
namespace ngcore
{
  template<int S, typename T>
  class MPI_typetrait<ngbla::Vec<S, T> >
  {
  public:
    /// gets the MPI datatype
    static MPI_Datatype MPIType () 
    { 
      static MPI_Datatype MPI_T = 0;
      if (!MPI_T)
	{
	  MPI_Type_contiguous ( S, MPI_typetrait<T>::MPIType(), &MPI_T);
	  MPI_Type_commit ( &MPI_T );
	}
      return MPI_T;
    }
  };
}
#endif

#endif
