#ifndef FILE_VECTOR_EXPR
#define FILE_VECTOR_EXPR

/**************************************************************************/
/* File:   vector.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jan. 02                                                    */
/**************************************************************************/

#include "expr.hpp"
#include <core/hashtable.hpp>   // for SparseVector

namespace ngbla
{

  template <typename T=double, typename TS=size_t, typename TDIST=IC<1>> class VectorView;


  template <typename T, typename ...Args>
  INLINE void SetVector (typename mat_traits<T>::TSCAL val, VectorView<T,Args...> vec) NETGEN_NOEXCEPT
  {
    for (size_t i : Range(vec))
      vec[i] = val;
  }
  

  
  template <typename T = double>
  using FlatVector = VectorView<T,size_t,IC<1>>;

  template <typename T, typename TELEM=typename T::TELEM>
  constexpr bool IsConvertibleToFlatVector () {
    return is_constructible_v<FlatVector<TELEM>,T>;
  }

  template <typename T, typename TELEM=typename T::TELEM>
  auto make_FlatVector (const T & v) {
    return FlatVector<TELEM> (v);
  }
  

  
  template <typename T = double>
  using SliceVector = VectorView<T,size_t,size_t>;

  template <typename T, typename TELEM=typename T::TELEM>
  constexpr bool IsConvertibleToSliceVector () {
    return is_constructible_v<SliceVector<TELEM>,T>;    
  }

  template <typename T, typename TELEM=typename T::TELEM>
  auto make_SliceVector (const T & v) {
    return SliceVector<TELEM> (v);
  }
  
  


  template <typename T = double>
  using BareSliceVector = VectorView<T,undefined_size,size_t>;
  
  template <typename T, typename TELEM=typename T::TELEM>
  constexpr bool IsConvertibleToBareSliceVector () {
    return is_constructible_v<BareSliceVector<TELEM>,T>;        
  }

  template <typename T, typename TELEM=typename T::TELEM>
  auto make_BareSliceVector (const T & v) {
    return BareSliceVector<TELEM> (v);
  }


  
  template <typename T = double>
  using BareVector = VectorView<T,undefined_size,IC<1>>;

  template <typename T, typename TELEM=typename T::TELEM>
  constexpr bool IsConvertibleToBareVector () {
    return is_constructible_v<BareVector<TELEM>,T>;            
  }

  template <typename T, typename TELEM=typename T::TELEM>
  auto make_BareVector (const T & v) {
    return BareVector<TELEM> (v);
  }
  
  
  template <typename T, typename TS>
  using LinearVector = VectorView<T,TS,IC<1>>;


  


  
  template <int S, class T> class Vec;

  template <int S, typename T = double>
  using FlatVec = VectorView<T,IC<S>,IC<1>>;
  
  template <int S, int D, typename T = double>
  using FlatSliceVec = VectorView<T,IC<S>,IC<D>>;
  
  template <int D, typename T = double>
  using FixSliceVec = VectorView<T,size_t,IC<D>>;  


  
  template <class T> class SysVector;
  template <class T = double> class Vector;
  template <int DIST, typename T> class FixSliceVector;

#ifdef WIN32
  // #pragma warning( disable : 4848) // support for standard attribute 'no_unique_address' in C++17 and earlier is a vendor extension
  // #define NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]

  // Disable this optimization with MSVC because it creates inconsistent results with different versions.
  // Following code returns 8 for compilers up to MSVC version 19.31, and returns 16 from version 19.32, see https://godbolt.org/z/v6P5vsq1M
  #define NO_UNIQUE_ADDRESS
  /*
  struct Base{};
  struct Empty {};

  template <typename T, typename TS>
  class Class : public Base{
      T* data;
      [[msvc::no_unique_address]] TS size;
  };

  int f() {
      return sizeof(Class<void, Empty>);
  }
  */

#else
  #define NO_UNIQUE_ADDRESS [[no_unique_address]]
#endif

  template <typename T, typename TS, typename TDIST>
  class VectorView : public MatExpr<VectorView<T,TS,TDIST>>
  {
  protected:
    typedef MatExpr<VectorView> BASE;
    T * __restrict data;
    NO_UNIQUE_ADDRESS TS size;
    NO_UNIQUE_ADDRESS TDIST dist;
  public:
    typedef T TELEM;
    typedef typename mat_traits<T>::TSCAL TSCAL;
    typedef TDIST type_dist;
    
    /// linear element access ? 
    static constexpr bool IsLinear() { return std::is_same<type_dist,IC<1>>(); }
    
    INLINE VectorView () = default;
    INLINE VectorView (const VectorView&) = default;
    INLINE VectorView (VectorView&&) = default;
    
    template <typename T2, typename TS2, typename TDIST2,
              enable_if_t<is_convertible<T2*,T*>::value, int> =0,
              enable_if_t<is_constructible<TS,TS2>::value, int> =0,
              enable_if_t<is_constructible<TDIST,TDIST2>::value, int> =0>
    INLINE VectorView (const VectorView<T2,TS2,TDIST2> & v2)
      : data{v2.Data()}, size{TS(v2.Size())}, dist{TDIST(v2.Dist())} { }
    
    INLINE explicit VectorView (T * adata)
      : data(adata)
    {
      ; // static_assert(std::is_same<type_dist,IC<1>>());
    } 
    INLINE VectorView (TS asize, T * adata)
      : data(adata), size(asize)
    {
      static_assert(std::is_same<type_dist,IC<1>>());
    } 
    INLINE VectorView (TS asize, TDIST adist, T * adata)
      : data(adata), size(asize), dist(adist) { } 

    /// allocate FlatVector on local heap
    INLINE VectorView (size_t as, LocalHeap & lh) 
      :  data(lh.Alloc<T> (as)), size(as), dist(IC<1>()) { }

    template <typename EXPR>
    INLINE VectorView(LocalHeapExpr<EXPR>&& lhe)
      : VectorView(lhe.Height(), lhe.GetLocalHeap())
    {
      *this = lhe.A();
    }

    /// put FlatVector over fixed size vector

    template <int S>
    INLINE VectorView (Vec<S,T> & v)
      : data(v.Data()), size(v.Size()), dist(IC<1>()) { }

    
    /*
    template <int S>
    INLINE VectorView (Vec<S,T> & v)
      : data(const_cast<T*>(v.Data())), size(v.Size()), dist(IC<1>()) { }
    */
    template <int S>
    INLINE VectorView (const Vec<S,T> & v)
      : data(const_cast<T*>(v.Data())), size(v.Size()), dist(IC<1>()) { }
    /*
    template <int S>
    INLINE VectorView (Vec<S,const T> & v)
      : data(const_cast<T*>(v.Data())), size(v.Size()), dist(IC<1>()) { }
    template <int S>
    INLINE VectorView (const Vec<S,const T> & v)
      : data(const_cast<T*>(v.Data())), size(v.Size()), dist(IC<1>()) { } 
    */
    
    INLINE auto Size() const { return size; }
    INLINE auto Dist() const { return dist; }
    INLINE auto Shape() const { return tuple(size); }
    

    INLINE auto Height () const { return size; }
    INLINE auto Width () const { return IC<1>(); }

    INLINE auto Range () const { return IntRange (0, size); }

    INLINE T * Addr(size_t i) const { return data+i*dist; }
    INLINE T * Data() const { return data; }

    INLINE auto View() const { return VectorView(*this); }         
    INLINE auto ViewRW() { return this->View(); }    



    /// assign memory for vector on local heap
    INLINE void AssignMemory (size_t as, LocalHeap & lh) 
    {
      size = as;
      dist = IC<1>();
      data = lh.Alloc<T>(size);
    }

    /// assign memory for vector
    void AssignMemory (size_t as, T * mem) 
    {
      size = as;
      dist = IC<1>(); 
      data = mem;
    }
    
    INLINE auto & operator= (const VectorView & v)
    {
      /*
      auto cs = CombinedSize(this->Size(), v.Size());      
      for (size_t i = 0; i < cs; i++)
      data[i*dist] = v(i);
      */
      this->template Assign<typename BASE::As> (v);      
      return *this;
    }
    INLINE auto & operator= (VectorView && v)
    {
      /*
      auto cs = CombinedSize(this->Size(), v.Size());      
      for (size_t i = 0; i < cs; i++)
        data[i*dist] = v(i);
      */
      this->template Assign<typename BASE::As> (v);
      return *this;
    }

    /// copy vector. sizes must match
    template <typename ...Args>
    INLINE auto & operator= (const VectorView<Args...> & v) 
    {
      // BASE::operator= (v);
      this->template Assign<typename BASE::As> (v);      
      return *this;
    }
    
    template<typename TB>
    INLINE auto & operator= (const Expr<TB> & v) 
    {
      // BASE::operator= (v);
      this->template Assign<typename BASE::As> (v);            
      return *this;
    }

    /// assign constant value
    INLINE auto & operator= (TSCAL scal) const
    {
      SetVector (scal, *this);
      return *this;
    }
    
    template <int D, typename TSCAL2>
    INLINE auto & operator= (const Vec<D,TSCAL2> & v) 
    {
      NETGEN_CHECK_RANGE(D,0,Size()+1);
      for (int i = 0; i < D; i++)
	data[i*dist] = v(i);
      return *this;
    }
    
    /// access element
    INLINE TELEM & operator[] (size_t i) 
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*dist]; 
    }

    INLINE const TELEM & operator[] (size_t i) const
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*dist]; 
    }

    
    /// constant element access
    INLINE TELEM & operator() (size_t i) const
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*dist]; 
    }

    /// element access. index j is ignored
    INLINE TELEM & operator() (size_t i, size_t j) const
    {
      NETGEN_CHECK_RANGE(i,0,Size());
      return data[i*dist]; 
    }
    
    RowsArrayExpr<VectorView> operator() (FlatArray<int> rows) const
    { 
      return RowsArrayExpr<VectorView> (*this, rows);
    }

    INLINE auto Range (size_t first, size_t next) const
    {
      NETGEN_CHECK_RANGE(next,first,Size()+1);
      return VectorView<T,size_t,TDIST> (next-first, dist, data+first*dist);
    }    
    
    INLINE auto Range (IntRange range) const
    {
      return Range (range.First(), range.Next());
    }

    /*
    template <int S> 
    INLINE auto Range (IC<S> r) const
    {
      return VectorView<T,IC<S>,TDIST> (r, Dist(), Data());
    }
    */

    INLINE auto Slice(size_t first, size_t dist2) const
    {
      // return VectorView<T,decltype(declval<TS>()/size_t()), decltype(declval<TDIST>()*size_t())> (size/dist2, dist2*dist, Addr(first));
      return VectorView<T,decltype(declval<TS>()/size_t()), decltype(declval<TDIST>()*size_t())> ( (size-first+dist2-1)/dist2, dist2*dist, Addr(first));
    }
    
    INLINE auto operator+(int i) const { return VectorView(size-i, dist, data+i*dist); }

    INLINE auto RemoveConst() const { return VectorView<typename remove_const<T>::type,TS,TDIST>(size, dist, const_cast<typename remove_const<T>::type*> (data)); }
    
    INLINE auto AsMatrix (size_t h, size_t w) const
    {
      // todo: checking
      static_assert(std::is_same<TDIST,IC<1>>());
      return FlatMatrix<T> (h,w, data);
    }

    class Iterator
    {
      VectorView vec;
      size_t ind;
    public:
      INLINE Iterator (VectorView avec, size_t ai) : vec(avec), ind(ai) { ; }
      INLINE Iterator operator++ (int) { return Iterator(vec, ind++); }
      INLINE Iterator operator++ () { return Iterator(vec, ++ind); }
      INLINE TELEM operator*() const { return vec[ind]; }
      INLINE TELEM & operator*() { return vec[ind]; }
      INLINE bool operator != (Iterator d2) { return ind != d2.ind; }
      INLINE bool operator == (Iterator d2) { return ind == d2.ind; }
    };
    
    INLINE Iterator begin() const { return Iterator (*this, 0); }
    INLINE Iterator end() const { return Iterator (*this, size); }
    
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


    /// allocate and copy vector 
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
                                                           &mem[0] : new T[as]) { ; }

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

  

    INLINE auto operator() (size_t i) 
    {
      return FlatVector<T> (blocksize, &data[i*blocksize]); 
    }

    INLINE auto operator() (size_t i) const
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
    static constexpr bool IsLinear() { return true; } 

    
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
    /*
    INLINE Vec (const Vec & v) : MatExpr<Vec> ()
    {
      for (size_t i = 0; i < S; i++)
	data[i] = v.data[i];
    }
    */
    
    Vec (const Vec &) = default;
    auto & HTData() const { return data; }                                    
    template <typename T2>
    Vec (const Vec<S,T2> & v2) : data(v2.HTData()) { ; }

    
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

    /*
    template <int D>
    INLINE Vec(FlatSliceVec<S,D,T> fsv)
    {
      for (int i = 0; i < S; i++)
        data[i] = fsv(i);
    }
    */
    template <typename T2, typename S2, typename D2>
    INLINE Vec(VectorView<T2,S2,D2> fsv)
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

    template <class... T2,
              enable_if_t<S==1+sizeof...(T2),bool> = true>
    Vec(const TELEM &v, T2... rest) {
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
	data[i] = v.Spec()(i);
      return *this;
    }

    // auto View() const { return FlatVec(const_cast<Vec&>(*this)); }
    // auto View() const { return Vec(*this); }
    INLINE auto View() const { return Vec<S,const T>{*this}; }    
    INLINE auto Shape() const { return tuple { IC<S>() }; }
    
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
    INLINE auto Size () const { return IC<S>(); }
    /// corresponding matrix height
    INLINE constexpr size_t Height () const { return S; }
    /// corresponding matrix with
    INLINE constexpr size_t Width () const { return 1; }
    INLINE T* Data() { return data.Ptr(); }
    INLINE const T* Data() const { return data.Ptr(); }

    INLINE /* const */ FlatVector<const T> Range(size_t first, size_t next) const
    { return FlatVector<const T> (next-first, data+first); }

    INLINE /* const */ FlatVector<T> Range(size_t first, size_t next) 
    { return FlatVector<T> (next-first, data+first); }

    const T * begin() const { return data.Ptr(); }
    const T * end() const { return data.Ptr()+S; }
    T * begin() { return data.Ptr(); }
    T * end() { return data.Ptr()+S; }
  };



  template <typename T>
  class Vec<0,T>  : public MatExpr<Vec<0,T> > 
  {
  public:
    static constexpr bool IsLinear() { return true; }     
    INLINE Vec () { }
    INLINE Vec (const Vec &d) { }
    INLINE Vec (T d) { }
    template<typename TB>
    INLINE Vec (const Expr<TB> & v) {;}
    INLINE auto View() const { return Vec(*this); }
    INLINE auto Shape() const { return tuple(IC<0>()); }    
    INLINE constexpr size_t Size() const { return 0; }
    INLINE constexpr size_t Height() const { return 0; }
    INLINE constexpr size_t Width() const { return 1; }
    template<typename TB>
    INLINE Vec & operator= (const Expr<TB> & v) { return *this;}
    INLINE Vec & operator= (const T & /* scal */) { return *this; } 
    
    INLINE T & operator[] (int i) const  { return *(T*)(void*)(this); }
    INLINE T & operator() (int i) const  { return *(T*)(void*)(this); }
    INLINE T & operator() (int i, int j) const  { return *(T*)(void*)(this); }
    INLINE T* Data() { return nullptr; }
    INLINE const T* Data() const { return nullptr; }
    INLINE FlatVector<const T> Range(size_t first, size_t next) const
    { return FlatVector<const T> (next-first, (T*)nullptr); }
    INLINE FlatVector<T> Range(size_t first, size_t next) 
    { return FlatVector<T> (next-first, (T*)nullptr); }
  };


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
    // enum { IS_COMPLEX = mat_traits<TSCAL>::IS_COMPLEX };
  };


  
  template <typename T> struct ConstVecSize { static constexpr int VSIZE = -1; }; 
  template <int S, typename T> struct ConstVecSize<Vec<S,T>> { static constexpr int VSIZE = S; };
  template <int S, typename T> struct ConstVecSize<FlatVec<S,T>> { static constexpr int VSIZE = S; };
  template <int S, int D, typename T> struct ConstVecSize<FlatSliceVec<S,D,T>> { static constexpr int VSIZE = S; };
  template <typename T>
  constexpr auto ConstVectorSize() { return ConstVecSize<T>::VSIZE; }
      
  /// cross product of 3-vectors
  template <typename TA, typename TB,
            std::enable_if_t<ConstVecSize<TA>::VSIZE == 3, bool> = true,
            std::enable_if_t<ConstVecSize<TB>::VSIZE == 3, bool> = true>
  INLINE auto Cross (const TA & a, const TB & b)
  {
    typedef decltype (a(0)*b(0)) T;
    return Vec<3,T>({ a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0) });
  }

  template <typename S>
  INLINE Vec<1,S> Cross (const Vec<2,S> & a, const Vec<2,S> & b)
  {
    return Vec<1,S> ( { a(0) * b(1) - a(1) * b(0) } );
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


  template<int S, typename T, typename TB>
  INLINE Vec<S,T> & operator+= (Vec<S,T> & v, const Expr<TB> & v2)
  {
    for (int i = 0; i < S; i++)
      v(i) += v2.Spec()(i,0);
    return v;
  }

  template<int S, typename T, typename TB>
  INLINE Vec<S,T> & operator-= (Vec<S,T> & v, const Expr<TB> & v2)
  {
    for (int i = 0; i < S; i++)
      v(i) -= v2.Spec()(i,0);
    return v;
  }






  
  /// output vector.
  template<int S, typename T>
  inline ostream & operator<< (ostream & ost, const FlatVec<S,T> & v)
  {
    for (int i = 0; i < S; i++)
      ost << " " << setw(7) << v(i);
    return ost;
  }





  template <class TV, class TSCAL> class Scalar2ElemVector
  {
  public:
    FlatVector<TSCAL> vec;
    Scalar2ElemVector (FlatVector<TSCAL> avec) : vec(avec) { }

    enum { H = mat_traits<TV>::HEIGHT };

    FlatVec<H,TSCAL> operator() (int i) const
    {
      return FlatVec<H,TSCAL> (&vec(i*H));
    }

  };
  

  template <class TSCAL> class Scalar2ElemVector<TSCAL,TSCAL>
  {
  public:
    FlatVector<TSCAL> vec;
    Scalar2ElemVector (FlatVector<TSCAL> avec) : vec(avec) { }


    const TSCAL & operator() (int i) const
    {
      return vec(i);
    }

    TSCAL & operator() (int i)
    {
      return vec(i);
    }

  };


  template <typename T=double>
  class SparseVector
  {
    size_t size;
    ClosedHashTable<size_t, T> data;
  public:
    SparseVector (size_t asize) : size(asize), data(10) { }
    size_t Size() const { return size; }
    T operator[] (size_t i) const { return data[i]; }
    T & operator[] (size_t i) { return data[i]; }
    T InnerProduct (FlatVector<T> v2) const
    {
      T sum = 0;
      for (auto [i,v] : data)
        sum += v * v2(i);
      return sum;
    }
    auto & Data() const { return data; }
  };

  template <typename T>  
  inline ostream & operator<< (ostream & ost, const SparseVector<T> & sv) {
    for (auto [i,v] : sv.Data())
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


  template <int S, typename T>
  INLINE auto operator* (double a, const Vec<S,T> & vec)
  {
    // typedef decltype(RemoveConst(a*vec(0))) TRES;
    typedef typename std::remove_const<decltype(a*std::declval<T>())>::type TRES;
    Vec<S, TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a * vec(i);
    return res;
  }

  template <int S, typename T>
  INLINE auto operator* (int a, const Vec<S,T> & vec) { return double(a)*vec; }
  
  template <int S, typename T>
  INLINE auto operator* (Complex a, const Vec<S,T> & vec) 
  {
    // typedef decltype(RemoveConst(a*vec(0))) TRES;
    typedef typename std::remove_const<decltype(a*std::declval<T>())>::type TRES;    
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
    // typedef decltype(RemoveConst(a*vec(0))) TRES;
    typedef typename std::remove_const<decltype(a*std::declval<T>())>::type TRES;    
    Vec<S, TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a * vec(i);
    return res;
  }
  
  template <int S, typename T>
  INLINE auto operator+ (const Vec<S,T> & a, const Vec<S,T> & b) 
  {
    // typedef decltype(RemoveConst(a(0))) TRES;
    typedef typename std::remove_const<T>::type TRES;        
    Vec<S,TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a(i)+b(i);
    return res;
  }

  template <int S, typename T>
  INLINE auto operator- (const Vec<S,T> & a, const Vec<S,T> & b) 
  {
    // typedef decltype(RemoveConst(a(0))) TRES;
    typedef typename std::remove_const<T>::type TRES;    
    Vec<S,TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a(i)-b(i);
    return res;
  }

  
  template <int S, typename T>
  INLINE auto operator* (double a, FlatVec<S,T> vec) 
  // -> Vec<S, decltype(RemoveConst(a*vec(0)))>
  {
    // typedef decltype(RemoveConst(a*vec(0))) TRES;
    typedef typename std::remove_const<decltype(a*std::declval<T>())>::type TRES;
    Vec<S, TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a * vec(i);
    return res;
  }


  template <int S, typename T>
  INLINE auto operator* (Complex a, FlatVec<S,T> vec) 
  // -> Vec<S, decltype(RemoveConst(a*vec(0)))>
  {
    // typedef decltype(RemoveConst(a*vec(0))) TRES;
    typedef typename std::remove_const<decltype(a*std::declval<T>())>::type TRES;    
    Vec<S, TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a * vec(i);
    return res;
  }

  template <int S, int D, typename T>
  INLINE auto operator* (double a, FlatSliceVec<S,D,T> vec) 
  // -> Vec<S, decltype(RemoveConst(a*vec(0)))>
  {
    // typedef decltype(RemoveConst(a*vec(0))) TRES;
    typedef typename std::remove_const<decltype(a*std::declval<T>())>::type TRES;        
    Vec<S, TRES> res;
    for (int i = 0; i < S; i++)
      res(i) = a * vec(i);
    return res;
  }


  template <int S, int D, typename T>
  INLINE auto operator* (Complex a, FlatSliceVec<S,D,T> vec) 
  // -> Vec<S, decltype(RemoveConst(a*vec(0)))>
  {
    // typedef decltype(RemoveConst(a*vec(0))) TRES;
    typedef typename std::remove_const<decltype(a*std::declval<T>())>::type TRES;            
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
  template <typename ARCHIVE, int S, typename T>
  inline auto & operator& (ARCHIVE & ar, ngbla::Vec<S,T> & v)
  {
    for (int i = 0; i < S; i++)
      ar & v(i);
    return ar;
  }
}


namespace ngcore
{
  template<typename T> struct MPI_typetrait;
  
  template<int S, typename T>
  struct MPI_typetrait<ngbla::Vec<S, T> > {
    static auto MPIType () {
      return MPI_typetrait<std::array<T,S>>::MPIType();
    }
  };
}

#endif
