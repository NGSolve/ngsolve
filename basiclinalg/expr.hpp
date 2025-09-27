#ifndef FILE_NGBLA_EXPR
#define FILE_NGBLA_EXPR

/**************************************************************************/
/* File:   expr.hpp                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jan. 02                                                    */
/**************************************************************************/

#include <core/array.hpp>
#include <core/localheap.hpp>
#include <core/exception.hpp>

#include <cstddef>
#include <ngs_stdcpp_include.hpp>  // for INLINE
#include "complex_wrapper.hpp"


#if defined(NETGEN_ENABLE_CHECK_RANGE) && !defined(__CUDA_ARCH__)
#define NETGEN_CHECK_SHAPE(a,b) \
  { if(a.Shape() != b.Shape()) \
      ngcore::ThrowException(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t: shapes don't match"); }
#else // defined(NETGEN_ENABLE_CHECK_RANGE) && !defined(__CUDA_ARCH__)
#define NETGEN_CHECK_SHAPE(a,b)
#endif // defined(NETGEN_ENABLE_CHECK_RANGE) && !defined(__CUDA_ARCH__)


template <typename T>
struct SafeIndex
{
  T i;
  INLINE SafeIndex(T ai) : i(ai) { };
  INLINE operator T() const { return i; }
  INLINE auto operator++() { return ++i; }
  INLINE auto operator++(int) { return i++; }
};

namespace ngcore
{
  template <typename T>
  struct IsSafe<SafeIndex<T>> {
    constexpr operator bool() const { return true; } };
} // namespace ngcore

/*
namespace std {
  template <typename T>  
  struct is_integral<SafeIndex<T>> {
    static constexpr bool value = true;
  };
}
*/

namespace ngbla
{
  using namespace std;
  using namespace ngcore;
  using namespace ngstd;

  
  enum ORDERING { ColMajor, RowMajor };


  struct unused_dist
  {
    unused_dist () = default;
    unused_dist (size_t d) { };
    template <int S>
    unused_dist (IC<S> d) { };
  };
  
  template <typename T = double, ORDERING ORD = RowMajor, typename TH=size_t, typename TW=size_t, typename TDIST=size_t>
  class MatrixView;

  template <typename T = double, ORDERING ORD = RowMajor>
  using FlatMatrix = MatrixView<T,ORD,size_t, size_t, unused_dist>;

  
  // template <typename T = double, ORDERING ORD = RowMajor> class FlatMatrix;
  template <typename T = double, ORDERING ORD = RowMajor> class Matrix;

  
  template <int H, int W, typename T> class Mat;
  template <int H, typename T> class DiagMat;
  template <int S, typename T> class Vec;


  template <typename T>
  struct is_scalar_type { static constexpr bool value = false; };
  
  template <typename T>
  constexpr bool IsScalar ()
  {
    return is_scalar_type<T>::value;
  }




  
  /*
    Matrix expression templates
  */

  
  template <typename TM, enable_if_t<!IsScalar<TM>(),bool> = true> 
  inline auto Access (const TM & mat, int i, int j)
  {
    return mat(i,j);
  }

  template <typename TM, enable_if_t<IsScalar<TM>(),bool> = true>  
  inline auto Access (const TM & mat, int i, int j)
  {
    return mat;
  }




  template <typename T> struct is_scalar_type;
    
  template<> struct is_scalar_type<int> { static constexpr bool value = true; };  
  template<> struct is_scalar_type<double> { static constexpr bool value = true; };
  template<> struct is_scalar_type<Complex> { static constexpr bool value = true; };
  


  /**
     Trait to obtain vector and scalar types for given matrix types.
     Is specified for double, Complex, AutoDiff<doube>, AutoDiff<Complex>
  */




  template <class T>  class mat_traits;

  template <class T>
  class mat_traits
  {
  public:
    /// matrix element
    typedef T TELEM;
    /// field of matrix element
    typedef T TSCAL;
    /// type of column vector
    typedef T TV_COL;
    /// type of row vector
    typedef T TV_ROW;
    /// matrix height
    // enum { HEIGHT = 1 };
    /// matrix with
    // enum { WIDTH  = 1  };
    // static constexpr int HEIGHT = 1;
    // static constexpr int WIDTH = 1;
    ///
    // enum { IS_COMPLEX = 0 };
  };

  template <class T>
  class mat_traits<const T> : public mat_traits<T> { };
  template <class T>
  class mat_traits<T&> : public mat_traits<T> { };


  /*
  template <int D>
  class mat_traits<ngcore::INT<D> >
  {
  public:
    typedef int TELEM;
    typedef int TSCAL;
    typedef int TV_COL;
    typedef int TV_ROW;
    // enum { HEIGHT = D };
    // enum { WIDTH = 1 };
    // enum { IS_COMPLEX = 0 };
  };
  */

  template <>
  class mat_traits<Complex>
  {
  public:
    typedef Complex TELEM;
    typedef Complex TSCAL;
    typedef Complex TV_COL;
    typedef Complex TV_ROW;
    // static constexpr int HEIGHT = 1;
    // static constexpr int WIDTH = 1;
    // enum { IS_COMPLEX = 1 };
  };



  /// Height of matrix
  template <class TM> 
  inline auto Height (const TM & m)
  {
    return m.Height();
  }

  /// Width of matrix
  template <class TM> 
  inline auto Width (const TM & m)
  {
    return m.Width();
  }

  template <> inline constexpr auto Height<double> (const double&) { return 1; }
  template <> inline constexpr auto Height<Complex> (const Complex&) { return 1; }
  template <> inline constexpr auto Width<double> (const double&) { return 1; }
  template <> inline constexpr auto Width<Complex> (const Complex&) { return 1; }

  /*
  template <class TM> 
  inline constexpr size_t Height () { return Height(TM()); }
  template <class TM> 
  inline constexpr size_t Width () { return Width(TM()); }
  */

  template <class TM> 
  inline constexpr auto Height () { return TM::Height(); }
  template <class TM> 
  inline constexpr auto Width () { return TM::Width(); }

  template <> inline constexpr auto Height<double> () { return 1; }
  template <> inline constexpr auto Height<Complex> () { return 1; }
  template <> inline constexpr auto Width<double> () { return 1; }
  template <> inline constexpr auto Width<Complex> () { return 1; }

  
  template <class TM> 
  inline constexpr bool IsComplex () { return IsComplex<typename mat_traits<TM>::TSCAL>(); }
  template <> inline constexpr bool IsComplex<double> () { return false; }
  template <> inline constexpr bool IsComplex<Complex> () { return true; }  

  

  template <class TA> class RowsArrayExpr;
  template <class TA> class ColsArrayExpr;
  template <class TA> class SubMatrixExpr;
  template <class TA> class RowExpr;
  template <class TA> class ColExpr;




  /**
     Expr is the base class for all matrix template expressions.
     Barton and Nackman Trick for template polymorphism, function Spec.

     provides Height and Width of matrix.
     IsLinear allows linear matrix element access.
  */

#ifdef NETGEN_ENABLE_CHECK_RANGE
  struct undefined_size
  {
    size_t size;
    
    undefined_size() = default;
    constexpr undefined_size(size_t s) : size(s) { }
    template <int S>
    explicit constexpr undefined_size(IC<S> s) : size(s) { }
    explicit constexpr operator size_t() const { return size; }
    explicit constexpr operator int() const { return size; }    
    explicit constexpr operator ptrdiff_t() const { return size; }
  };

    
  inline ostream & operator<< (ostream & ost, undefined_size s) { ost << "undefined("<<size_t(s)<<")"; return ost; }
  inline constexpr auto operator/ (undefined_size ud, size_t i) { return undefined_size(size_t(ud)/i); }
  inline constexpr auto operator- (undefined_size ud, size_t i) { return undefined_size(size_t(ud)-i); }
  inline constexpr auto operator+ (undefined_size ud, size_t i) { return undefined_size(size_t(ud)+i); }
  inline constexpr auto operator* (undefined_size ud, size_t i) { return undefined_size(size_t(ud)*i); }
  inline constexpr auto operator* (size_t i, undefined_size ud) { return undefined_size(size_t(ud)*i); }    
  inline constexpr bool operator< (size_t i, undefined_size ud) { return i < size_t(ud); }
  inline constexpr bool operator< (undefined_size ud, size_t i) { return size_t(ud) < i; }
  inline constexpr bool operator>= (size_t i, undefined_size ud) { return i >= size_t(ud); }
  inline constexpr bool operator>= (undefined_size ud, size_t i) { return size_t(ud) >= i; }
  inline constexpr bool operator== (size_t i, undefined_size ud) { return i == size_t(ud); }
  inline constexpr bool operator== (undefined_size ud, size_t i) { return size_t(ud) == i; }
  inline constexpr bool operator== (undefined_size ud, undefined_size ud2) { return size_t(ud) == size_t(ud2); }  
  inline constexpr bool operator!= (size_t i, undefined_size ud) { return i != size_t(ud); }
  inline constexpr bool operator!= (undefined_size ud, size_t i) { return size_t(ud) != i; }
  inline constexpr bool operator!= (undefined_size ud, undefined_size ud2) { return size_t(ud) != size_t(ud2); }  

#else
  struct undefined_size
    {
      undefined_size() = default;
      undefined_size(size_t s) { }
      template <int S>
      explicit constexpr undefined_size(IC<S> s) { }
  };
  
  inline ostream & operator<< (ostream & ost, undefined_size s) { ost << "undefined"; return ost; }
  inline auto operator/ (undefined_size ud, size_t i) { return ud; }
  inline auto operator- (undefined_size ud, size_t i) { return ud; }
  inline auto operator+ (undefined_size ud, size_t i) { return ud; }
#endif

  
  INLINE constexpr auto CombinedSize(undefined_size s1, undefined_size s2) {
    NETGEN_CHECK_SAME(size_t(s1), size_t(s2)); return s1; }
  INLINE constexpr auto CombinedSize(undefined_size s1, size_t s2) {
    NETGEN_CHECK_SAME(size_t(s1), size_t(s2)); return s2; }  
  INLINE constexpr auto CombinedSize(size_t s1, undefined_size s2) {
    NETGEN_CHECK_SAME(size_t(s1), size_t(s2)); return s1; }  
  INLINE constexpr auto CombinedSize(size_t s1, size_t s2) {
    NETGEN_CHECK_SAME(size_t(s1), size_t(s2)); return s1; }
  template <int S1> INLINE constexpr auto CombinedSize(IC<S1> s1, undefined_size s2) {
    NETGEN_CHECK_SAME(size_t(s1), size_t(s2)); return s1; }  
  template <int S1> INLINE constexpr auto CombinedSize(IC<S1> s1, size_t s2) {
    NETGEN_CHECK_SAME(size_t(s1), size_t(s2)); return s1; }  
  template <int S1, int S2> INLINE constexpr auto CombinedSize(IC<S1> s1, IC<S2> s2) {
    NETGEN_CHECK_SAME(size_t(s1), size_t(s2)); return s1; }  
  template <int S2> INLINE constexpr auto CombinedSize(undefined_size s1, IC<S2> s2) {
    NETGEN_CHECK_SAME(size_t(s1), size_t(s2)); return s2; }  
  template <int S2> INLINE constexpr auto CombinedSize(size_t s1, IC<S2> s2) {
     NETGEN_CHECK_SAME(size_t(s1), size_t(s2)); return s2; }  

  template <typename T1, typename T2>
  INLINE constexpr auto CombinedSize(tuple<T1> tup1, tuple<T2> tup2)
  { return tuple(CombinedSize(get<0>(tup1), get<0>(tup2))); }

  template <typename T11, typename T12, typename T21, typename T22>
  INLINE constexpr auto CombinedSize(tuple<T11,T12> tup1, tuple<T21,T22> tup2)
  { return tuple(CombinedSize(get<0>(tup1), get<0>(tup2)),
                 CombinedSize(get<1>(tup1), get<1>(tup2))); }



    
  
  
  template <typename T>
  class Expr 
  {
  public:
    constexpr Expr () = default;

    /// cast to specific type
    INLINE T & Spec() { return static_cast<T&> (*this); }

    /// cast to specific type
    INLINE const T & Spec() const { return static_cast<const T&> (*this); }

    INLINE auto View() const { return static_cast<const T&> (*this).View(); }
    INLINE decltype(auto) ViewRW() { return static_cast<T&>(*this).ViewRW(); }
    INLINE auto Shape() const { return Spec().T::Shape(); }


    INLINE auto Height() const { return Spec().T::Height(); }
    INLINE auto Width() const { return Spec().T::Width(); }
    

    void Dump (ostream & ost) const { Spec().T::Dump(ost); }


    INLINE auto Row (size_t r) const
    {
      // return RowExpr<const T> (static_cast<const T&> (*this), r);
      return RowExpr<const T> (this->View(), r);
    }

    INLINE auto Col (size_t r) const
    {
      // return ColExpr<const T> (static_cast<const T&> (*this), r);
      return ColExpr<const T> (this->View(), r);
    }


    INLINE SubMatrixExpr<T>
    Rows (size_t first, size_t next) const
    { 
      // return SubMatrixExpr<T> (static_cast<T&> (*this), first, 0, next-first, Width());
      return SubMatrixExpr<T> (this->View(), first, 0, next-first, Width()); 
    }

    INLINE SubMatrixExpr<T>
    Cols (size_t first, size_t next) const
    { 
      // return SubMatrixExpr<T> (static_cast<T&> (*this), 0, first, Height(), next-first);
      return SubMatrixExpr<T> (this->View(), 0, first, Height(), next-first);
    }

    INLINE SubMatrixExpr<T>
    Rows (IntRange range) const
    { 
      return Rows (range.First(), range.Next());
    }

    INLINE SubMatrixExpr<T>
    Cols (IntRange range) const
    { 
      return Cols (range.First(), range.Next());
    }


    INLINE RowsArrayExpr<T>
    Rows (FlatArray<int> rows) const
    { 
      return RowsArrayExpr<T> (static_cast<const T&> (*this), rows); 
    }

    INLINE ColsArrayExpr<T>
    Cols (FlatArray<int> cols) const
    { 
      return ColsArrayExpr<T> (static_cast<const T&> (*this), cols); 
    }

  };









  /**
     Caller knows that matrix expression is a symmetric matrix.
     Thus, only one half of the matrix needs to be computed.
  */
  template <typename T>
  class SymExpr : public Expr<SymExpr<T> >
  {
    T a;
  public:

    SymExpr (T aa) : a(aa) { ; }

    INLINE auto operator() (size_t i) const { return a(i); }
    INLINE auto operator() (size_t i, size_t j) const { return a(i,j); }
    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return a.Width(); }
    
    auto View() const { return *this; }
    auto Shape() const { return a.Shape(); }
    
    static constexpr bool IsLinear() { return T::IsLinear(); }
    void Dump (ostream & ost) const
    { ost << "Sym ("; a.Dump(ost); ost << ")"; }
  };


  /**
     Declare that matrix expression is symmetric
  */
  template <typename T>
  inline SymExpr<T> Symmetric (const Expr<T> & a)
  {
    return SymExpr<T> (a.View());
  }




  template <typename TA>
  class LocalHeapExpr : public Expr<LocalHeapExpr<TA> >
  {
    const TA & a;
    LocalHeap * lh;
  public:
    INLINE LocalHeapExpr (const TA & aa, LocalHeap & alh) : a(aa), lh(&alh) { ; }
    INLINE const TA & A() const { return a; }
    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return a.Width(); }
    INLINE LocalHeap & GetLocalHeap() const { return *lh; }
  };
  
  template <typename TA>
  INLINE LocalHeapExpr<TA> operator| (const Expr<TA> & a, LocalHeap & lh)
  {
    return LocalHeapExpr<TA> (a.Spec(), lh);
  }


  template <class TA> class MatExpr;
  template <class TA, class TB> class MultExpr;
  template <class TA> class MinusExpr;
  template <class TA> class TransExpr;
  template <class TA, class TS> class ScaleExpr;


  
  template <typename TOP, typename T, typename TB, typename Enable=int>
  class assign_trait
  {
  public:
    static INLINE T & Assign (MatExpr<T> & self, const Expr<TB> & v)
    {
      // NETGEN_CHECK_RANGE(self.Height(), v.Height(), v.Height()+1);
      // NETGEN_CHECK_RANGE(self.Width(), v.Width(), v.Width()+1);
      // NETGEN_CHECK_SHAPE(self.Spec(), v);


      auto src = v.View();
      decltype(auto) dest = self.ViewRW();

      auto h = CombinedSize (src.Height(), dest.Height()); // checks if same
      auto w = CombinedSize (src.Width(), dest.Width());   // checks if same
      
      if (T::COL_MAJOR)
        {
          if (h > 0)
            for (size_t j = 0; j < w; j++)
              for (size_t i = 0; i < h; i++)
                TOP()(dest(i,j), src(i,j));
          return self.Spec();
        }


      if constexpr (TB::IsLinear())
	{
	  if (T::IsLinear())
	    {
	      auto hw = h*w;
              for (SafeIndex<size_t> i : Range(hw))  
                TOP()(dest(i), src(i));
	    }
	  else
	    {
              if (w > 0)
                for (SafeIndex<size_t> i = 0, k = 0; i < h; i++)
                  for (SafeIndex<size_t> j = 0; j < w; j++, k++)
                    TOP() (dest(i,j), src(k));
	    }
	}
      else
	{
          if (w > 0)
            {
              if (T::IsLinear())
                for (SafeIndex<size_t> i = 0, k = 0; i < h; i++)
                  for (SafeIndex<size_t> j = 0; j < w; j++, k++)
                    TOP() (dest(k), src(i,j));                    
              else
                {
                  for (SafeIndex<size_t> i = 0; i < h; i++)
                    for (SafeIndex<size_t> j = 0; j < w; j++)
                      TOP() (dest(i,j), src(i,j));                        
                }
            }
        }
      return self.Spec();
     
    }
  };
  


  

  
  /**
     The base class for matrices.
  */
  template <class T>
  class MatExpr : public Expr<T>
  {
  public:
 
    constexpr MatExpr () = default;

    using Expr<T>::Spec;
    using Expr<T>::Height;
    using Expr<T>::Width;

    enum { COL_MAJOR = 0 };  // matrix is stored col-major

    void Dump (ostream & ost) const { ost << "Matrix"; }



    template<typename TOP, typename TB>
    INLINE auto & Assign (const Expr<TB> & v)
    {
      return assign_trait<TOP, T, TB>::Assign (*this, v);
    }


    class As 
    {
    public:
      template <typename T1, typename T2> 
      INLINE void operator() (T1 && v1, const T2 & v2) { v1 = v2; }
      static constexpr bool IsPos() { return true; } 
      static constexpr bool IsAdd() { return false; } 
    };
    class AsAdd 
    {
    public:
      template <typename T1, typename T2> 
      INLINE void operator() (T1 && v1, const T2 & v2) { v1 += v2; }
      static constexpr bool IsPos() { return true; } 
      static constexpr bool IsAdd() { return true; } 
    };
    class AsSub 
    {
    public:
      template <typename T1, typename T2> 
      INLINE void operator() (T1 && v1, const T2 & v2) { v1 -= v2; }
      static constexpr bool IsPos() { return false; } 
      static constexpr bool IsAdd() { return true; } 
    };
	


    

    /*
      TODO: move to traits ... 
    
    // x += s*y
    template <typename OP, typename TA, 
              enable_if_t<std::is_same<OP,AsAdd>::value,bool> = true,
              enable_if_t<is_constructible_v<SliceVector<double>,typename pair<T,TA>::first_type>,bool> = true,
              enable_if_t<is_constructible_v<SliceVector<double>,TA>,bool> = true>
    INLINE T & Assign (const Expr<ScaleExpr<TA,double>> & scaled)
    {
      AddVector (scaled.View().S(),
                 make_SliceVector(scaled.View().A()),
                 make_SliceVector(this->Spec()));
      return Spec();
    }

    
    // x += s*(m*y)
    template <typename OP, typename TA, typename TB,
              enable_if_t<std::is_same_v<OP,AsAdd>,bool> = true,
              enable_if_t<IsConvertibleToSliceMatrix<TA,double>(),bool> = true,
              enable_if_t<is_convertible_v<TB,FlatVector<double>>,bool> = true,
              enable_if_t<is_convertible<typename pair<T,TB>::first_type,FlatVector<double>>::value,bool> = true>
    INLINE T & Assign (const Expr<ScaleExpr<MultExpr<TA, TB>,double>> & prod)
    {
      MultAddMatVec (prod.Spec().S(),
                     make_SliceMatrix(prod.Spec().A().A()),
                     prod.Spec().A().B(),
                     Spec());
      return Spec();
    }

    // x += (s*m)*y
    template <typename OP, typename TA, typename TB,
              typename enable_if<std::is_same<OP,AsAdd>::value,int>::type = 0,
              typename enable_if<IsConvertibleToSliceMatrix<TA,double>(),int>::type = 0,
              typename enable_if<is_convertible<TB,FlatVector<double>>::value,int>::type = 0,
              typename enable_if<is_convertible<typename pair<T,TB>::first_type,FlatVector<double>>::value,int>::type = 0>
    INLINE T & Assign (const Expr<MultExpr<ScaleExpr<TA,double>, TB>> & prod)
    {
      MultAddMatVec (prod.Spec().A().S(),
                     make_SliceMatrix(prod.Spec().A().A()),
                     prod.Spec().B(),
                     Spec());
      return Spec();
    }
    */

    
    template<typename TB>
    INLINE T & operator= (const Expr<TB> & v)
    {
      Assign<As> (v);
      return Spec();
    }

    INLINE T & operator= (const T & v)
    {
      Assign<As> (v);
      return Spec();
    }

    template<typename TB>
    INLINE T & operator+= (const Expr<TB> & v)
    {
      Assign<AsAdd> (v);
      return Spec();
    }

    template<typename TB>
    INLINE MatExpr<T> & operator-= (const Expr<TB> & v)
    {
      Assign<AsSub> (v);
      return Spec();
    }




    template<typename TB>
    INLINE T & operator+= (const Expr<SymExpr<TB> > & v)
    {
      NETGEN_CHECK_RANGE(Height(), v.Height(), v.Height()+1);
      NETGEN_CHECK_RANGE(Width(), v.Width(), v.Width()+1);
      size_t h = Height();
      for (size_t i = 0; i < h; i++)
	{
	  for (size_t j = 0; j < i; j++)
	    {
	      double val = v.Spec()(i,j);
	      Spec()(i,j) += val;
	      Spec()(j,i) += val;
	    }
	  Spec()(i,i) += v.Spec()(i,i);
	}
      return Spec();
    }



    template <class SCAL2>
    INLINE T & operator*= (SCAL2 s)
    {
      if (T::IsLinear())
	{
	  size_t hw = Height() * Width();
	  for (size_t i = 0; i < hw; i++)
	    Spec()(i) *= s;
	}
      else
	for (size_t i = 0; i < Height(); i++)
	  for (size_t j = 0; j < Width(); j++)
	    Spec()(i,j) *= s;
	
      return Spec();
    }

    template <class SCAL2>
    INLINE T & operator/= (SCAL2 s)
    {
      return (*this) *= (1./s);
    }
  };







  /* *************************** SumExpr **************************** */

  /**
     Sum of 2 matrix expressions
  */

  template <class TA, class TB> 
  class SumExpr : public Expr<SumExpr<TA,TB> >
  {
    TA a;
    TB b;
  public:

    static constexpr bool IsLinear() { return TA::IsLinear() && TB::IsLinear(); }     
    
    INLINE SumExpr (TA aa, TB ab) : a(aa), b(ab) { ; }

    template <typename ...I>
    INLINE auto operator() (I... i) const { return a(i...)+b(i...); }    
    
    INLINE auto Height() const { return CombinedSize(a.Height(), b.Height()); }
    INLINE auto Width() const { return CombinedSize(a.Width(), b.Width()); }

    INLINE auto View() const { return SumExpr(a,b); }
    INLINE auto Shape() const { return CombinedSize(a.Shape(), b.Shape()); }     
    
    void Dump (ostream & ost) const
    { ost << "("; a.Dump(ost); ost << ") + ("; b.Dump(ost); ost << ")"; }
  };

  template <typename TA, typename TB>
  INLINE auto
  operator+ (const Expr<TA> & a, const Expr<TB> & b)
  {
    NETGEN_CHECK_SHAPE(a, b);
    return SumExpr(a.View(), b.View());
  }




  /* *************************** SubExpr **************************** */


  /**
     Matrix-expr minus Matrix-expr
  */

  template <class TA, class TB> 
  class SubExpr : public Expr<SubExpr<TA,TB> >
  {
    TA a;
    TB b;
  public:

    static constexpr bool IsLinear() { return TA::IsLinear() && TB::IsLinear(); }     
    
    INLINE SubExpr (TA aa, TB ab) : a(aa), b(ab) { ; }

    template <typename ...I>
    INLINE auto operator() (I... i) const { return a(i...)-b(i...); }
    
    INLINE auto View() const { return SubExpr(a,b); }
    INLINE auto Shape() const { return CombinedSize(a.Shape(), b.Shape()); }
    
    INLINE auto Height() const { return CombinedSize(a.Height(), b.Height()); }
    INLINE auto Width() const { return CombinedSize(a.Width(), b.Width()); }
  };


  template <typename TA, typename TB>
  inline auto operator- (const Expr<TA> & a, const Expr<TB> & b)
  {
    NETGEN_CHECK_SHAPE(a, b);    
    return SubExpr(a.View(), b.View());
  }







  /* *************************** MinusExpr **************************** */


  /**
     minus Matrix-expr
  */

  template <class TA>
  class MinusExpr : public Expr<MinusExpr<TA> >
  {
    TA a;
  public:
    MinusExpr (TA aa) : a(aa) { ; }

    template <typename ...I>
    INLINE auto operator() (I... i) const { return -a(i...); }
    
    INLINE auto View() const { return MinusExpr(a); }
    INLINE auto Shape() const { return a.Shape(); }
    
    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return a.Width(); }
    INLINE TA A() const { return a; }

    static constexpr bool IsLinear() { return TA::IsLinear(); } 
  };

  template <typename TA>
  INLINE auto operator- (const Expr<TA> & a)
  {
    return MinusExpr (a.View());
  }


  /* *************************** PW_Mult_Expr **************************** */

  template <class TA, class TB>
  class PW_Mult_Expr : public Expr<PW_Mult_Expr<TA,TB> >
  {
    TA a;
    TB b;
  public:
    static constexpr bool IsLinear() { return TA::IsLinear() && TB::IsLinear(); }     

    INLINE PW_Mult_Expr (TA aa, TB ab) : a(aa), b(ab) { ; }

    INLINE auto operator() (size_t i) const { return a(i)*b(i); }
    INLINE auto operator() (size_t i, size_t j) const { return a(i,j)*b(i,j); }

    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return a.Width(); }
    
    INLINE auto View() const { return *this; }
    INLINE auto Shape() const { return CombinedSize(a.Shape(), b.Shape()); }
    void Dump (ostream & ost) const
    { ost << "("; a.Dump(ost); ost << ") + ("; b.Dump(ost); ost << ")"; }
  };

  template <typename TA, typename TB>
  INLINE auto pw_mult (const Expr<TA> & a, const Expr<TB> & b)
  {
    NETGEN_CHECK_SHAPE(a, b);    
    return PW_Mult_Expr (a.View(), b.View());
  }


  /* *************************** PW_Inv_Expr **************************** */

  template <class TA>
  class PW_Inv_Expr : public Expr<PW_Inv_Expr<TA> >
  {
    TA a;
  public:
    static constexpr bool IsLinear() { return TA::IsLinear(); }     
    INLINE PW_Inv_Expr (TA aa) : a(aa) { ; }

    INLINE auto operator() (size_t i) const { return 1.0/a(i); }
    INLINE auto operator() (size_t i, size_t j) const { return 1.0/a(i,j); }

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }
    INLINE auto View () const { return *this; }
    void Dump (ostream & ost) const
    { ost << "1/("; a.Dump(ost); ost << ")"; }
  };

  template <typename TA>
  INLINE auto pw_inv (const Expr<TA> & a)
  {
    return PW_Inv_Expr (a.View());
  }



  /* *************************** ScaleExpr **************************** */


  /**
     Scalar times Matrix-expr
  */
  template <class TA, class TS> 
  class ScaleExpr : public Expr<ScaleExpr<TA,TS> >
  {
    TA a;
    TS s;
  public:
    static constexpr bool IsLinear() { return TA::IsLinear(); }

    INLINE ScaleExpr (TA aa, TS as) : a(aa), s(as) { ; }
    
    // INLINE auto operator() (size_t i) const { return s * a(i); }
    // INLINE auto operator() (size_t i, size_t j) const { return s * a(i,j); }

    template <typename ...I>
    INLINE auto operator() (I... i) const { return s*a(i...); }

    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return a.Width(); }

    INLINE auto View() const { return *this; }
    INLINE auto Shape() const { return a.Shape(); }         
    
    INLINE TA A() const { return a; }
    INLINE TS S() const { return s; }
    
    void Dump (ostream & ost) const
    { ost << "Scale, s=" << s << " * "; a.Dump(ost);  }
  };


  template <typename TS, typename TA,
            typename enable_if<IsScalar<TS>(),int>::type = 0>
  INLINE auto operator* (TS s, const Expr<TA> & a)
  {
    return ScaleExpr (a.View(), s);
  }



  /* ************************* MultExpr ************************* */


  /**
     Matrix-expr timex Matrix-expr
  */
  template <class TA, class TB> class MultExpr : public Expr<MultExpr<TA,TB> >
  {
    TA a;
    TB b;
  public:

    INLINE MultExpr (TA aa, TB ab) : a(aa), b(ab) { ; }

    /*
    INLINE auto operator() (size_t i) const
    { return operator()(i,0); }  

    INLINE auto operator() (size_t i, size_t j) const // -> decltype (a(0,0)*b(0,0))
    { 
      size_t wa = a.Width();

      if (wa >= 1)
	{
	  auto sum = a(i,0) * b(0,j);
	  for (size_t k = 1; k < wa; k++)
	    sum += a(i,k) * b(k,j);
          return sum;
	}

      decltype (a(0,0)*b(0,0)) sum (0);
      return sum;
    }
    */
    template <typename ...J>
    INLINE auto operator() (size_t i, J... j) const
    { 
      size_t wa = a.Width();

      if (wa >= 1)
	{
	  auto sum = a(i,SafeIndex(0)) * b(0,j...);
	  for (size_t k = 1; k < wa; k++)
	    sum += a(i,SafeIndex(k)) * b(k,j...);
          return sum;
	}

      decltype (a(0,0)*b(0,j...)) sum (0);
      return sum;
    }

    INLINE auto View() const { return MultExpr(a,b); }
    INLINE auto Shape() const
    {
      if constexpr (tuple_size<decltype(b.Shape())>() == 1)
        return tuple(get<0>(a.Shape()));
      else
        return tuple(get<0>(a.Shape()), get<1>(b.Shape()));
      /*
        // too complicated ? 
      return tuple_cat(tuple(get<0>(a.Shape())),
                       std::apply([](auto&&, const auto&... args) {return std::tie(args...);}, b.Shape()) );
      */
    }
    
    INLINE TA A() const { return a; }
    INLINE TB B() const { return b; }
    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return b.Width(); }
    static constexpr bool IsLinear() { return false; }     
  };


  template <int H, typename SCALA, class TB> class MultExpr<DiagMat<H,SCALA>,TB> 
    : public Expr<MultExpr<DiagMat<H,SCALA>,TB> >
  {
    DiagMat<H,SCALA> a;
    TB b;
  public:

    MultExpr (DiagMat<H,SCALA> aa, TB ab) : a(aa), b(ab) { ; }

    INLINE auto operator() (size_t i) const { return a[i] * b(i); }  
    INLINE auto operator() (size_t i, size_t j) const { return a[i] * b(i,j); }

    INLINE auto View() const { return *this; }
    INLINE auto Shape() const
    {
      typedef decltype(b.Shape()) TBSHAPE;
      if constexpr (tuple_size<TBSHAPE>() == 1)
                     return tuple<size_t> (H);
      else
        return tuple<size_t,size_t> (H, b.Width());
    }
    
    INLINE const auto A() const { return a; }
    INLINE const auto B() const { return b; }
    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return b.Width(); }

    static constexpr bool IsLinear() { return false; }         
  };


  template <typename TA, typename TB>
  INLINE auto operator* (const Expr<TA> & a, const Expr<TB> & b)
  {
    return MultExpr (a.View(), b.View());
  }


  /* ************************** Trans *************************** */

  template <typename TA,
            typename enable_if<IsScalar<TA>(),int>::type = 0>
  INLINE auto Trans (TA a) { return a; } 
  

  /**
     Transpose of Matrix-expr
  */
  template <class TA> class TransExpr : public MatExpr<TransExpr<TA> >
  {
    TA a;
  public:
    INLINE TransExpr (TA aa) : a(aa) { ; }

    INLINE auto Height() const { return a.Width(); }
    INLINE auto Width() const { return a.Height(); }

    INLINE auto operator() (size_t i, size_t j) const { return Trans (a(j,i)); }
    INLINE auto operator() (size_t i) const { return Trans(a(0,0)); }
    // auto Row (int i) const -> decltype (a.Col(i)) { return a.Col(i); }
    // auto Col (int i) const -> decltype (a.Row(i)) { return a.Row(i); }

    INLINE auto View() const { return *this; }
    INLINE auto Shape() const { return tuple (a.Width(), a.Height()); }

    static constexpr bool IsLinear() { return false; }     
    INLINE const TA & A() const { return a; }
  };


  /// Transpose 
  template <typename TA>
  INLINE auto Trans (const Expr<TA> & a)
  {
    return TransExpr (a.View());
  }

  template <typename TA, typename TB>
  INLINE auto Trans (const Expr<MultExpr<TA,TB>> & expr)
  {
    return Trans(expr.Spec().B()) * Trans(expr.Spec().A());
  }

  template <typename TA>
  INLINE auto Trans (const Expr<TransExpr<TA>> & expr)
  {
    return expr.Spec().A();
  }
  
  /* ************************* Real/Imag ************************ */
  
  INLINE double Real(double a) { return a; }
  INLINE double Imag(double a) { return 0; }

  INLINE double Real(Complex a) { return a.real(); }
  INLINE double Imag(Complex a) { return a.imag(); }
  
  template <class TA>
  class RealExpr : public Expr<RealExpr<TA> >
  {
    TA a;
  public:
    RealExpr (TA aa) : a(aa) { ; }

    INLINE auto operator() (size_t i) const { return Real(a(i)); }
    INLINE auto operator() (size_t i, size_t j) const { return Real(a(i,j)); }
    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return a.Width(); }
    INLINE auto View() const { return *this; }
    INLINE auto Shape() const { return a.Shape(); }
    static constexpr bool IsLinear() { return TA::IsLinear(); } 
  };

  template <typename TA>
  INLINE auto Real(const Expr<TA> & a)
  {
    return RealExpr(a.View());
  }


  template <class TA>
  class ImagExpr : public Expr<ImagExpr<TA> >
  {
    TA a;
  public:
    ImagExpr (TA aa) : a(aa) { ; }

    INLINE auto operator() (size_t i) const { return Imag(a(i)); }
    INLINE auto operator() (size_t i, size_t j) const { return Imag(a(i,j)); }
    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return a.Width(); }
    INLINE auto View() const { return *this; }
    INLINE auto Shape() const { return a.Shape(); }    
    static constexpr bool IsLinear() { return TA::IsLinear(); }     
  };

  template <typename TA>
  INLINE auto Imag(const Expr<TA> & a)
  {
    return ImagExpr (a.View());
  }


  
  /* ************************* SubMatrix ************************ */

  template <class TA> 
  class SubMatrixExpr : public MatExpr<SubMatrixExpr<TA> >
  {
    TA a;
    size_t first_row, first_col;
    size_t height, width;
  public:
    SubMatrixExpr (TA aa, size_t fr, size_t fc, size_t ah, size_t aw) 
      : a(aa), first_row(fr), first_col(fc), height(ah), width(aw) { ; }

    INLINE size_t Height() const { return height; }
    INLINE size_t Width() const { return width; }

    // auto operator() (size_t i, size_t j) { return a(i+first_row, j+first_col); }
    // auto operator() (size_t i) { return a(i+first_row); }
    INLINE decltype(auto) operator() (size_t i, int j) const  { return a(i+first_row, j+first_col); }
    INLINE decltype(auto) operator() (size_t i) const { return a(i+first_row); }

    typedef typename TA::TELEM TELEM;
    typedef typename TA::TSCAL TSCAL;

    static constexpr bool IsLinear() { return false; }
    
    enum { COL_MAJOR = TA::COL_MAJOR };

    template<typename TB>
    INLINE const SubMatrixExpr & operator= (const Expr<TB> & m) 
    {
      MatExpr<SubMatrixExpr<TA> >::operator= (m);
      return *this;
    }

    INLINE auto View() const { return SubMatrixExpr(a, first_row, first_col, height, width); }
    INLINE auto ViewRW() { return SubMatrixExpr(a, first_row, first_col, height, width); }
    INLINE tuple<size_t,size_t> Shape() const { return { height, width }; }
  };


  template <class TA> 
  class RowExpr : public MatExpr<RowExpr<TA> >
  {
    TA a;
    size_t row;
  public:
    RowExpr (TA aa, size_t r)
      : a(aa), row(r) { ; }

    INLINE size_t Height() const { return 1; }
    INLINE size_t Width() const { return a.Width(); }

    INLINE auto operator() (size_t i, size_t j) -> decltype(a(0,0)) { return a(row,i); }
    INLINE auto operator() (size_t i) -> decltype(a(0,0)) { return a(row,i); }
    INLINE auto operator() (size_t i, size_t j) const { return a(row,i); }
    INLINE auto operator() (size_t i) const { return a(row,i); }

    static constexpr bool IsLinear() { return false; }
    
    template<typename TB>
    INLINE const RowExpr & operator= (const Expr<TB> & m) 
    {
      MatExpr<RowExpr<TA> >::operator= (m);
      return *this;
    }

    INLINE auto View() const { return *this; }
    INLINE tuple<size_t> Shape() const { return a.Width(); }
  };



  template <class TA> 
  class ColExpr : public MatExpr<ColExpr<TA> >
  {
    TA a;
    size_t col;
  public:
    ColExpr (TA aa, size_t c)
      : a(aa), col(c) { ; }

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return 1; }

    INLINE auto operator() (size_t i, size_t j) -> decltype(a(0,0)) { return a(i,col); }
    INLINE auto operator() (size_t i) -> decltype(a(0,0)) { return a(i,col); }
    INLINE auto operator() (size_t i, size_t j) const { return a(i,col); }
    INLINE auto operator() (size_t i) const { return a(i,col); }

    static constexpr bool IsLinear() { return false; }
    
    template<typename TB>
    INLINE const ColExpr & operator= (const Expr<TB> & m) 
    {
      MatExpr<ColExpr<TA> >::operator= (m);
      return *this;
    }

    INLINE auto View() const { return *this; }
    INLINE tuple<size_t> Shape() const { return a.Height(); }
  };





  /* ************************* RowsArray ************************ */

  /**
     RowsArray
  */
  template <class TA> class RowsArrayExpr : public MatExpr<RowsArrayExpr<TA> >
  {
    const TA & a;
    FlatArray<int> rows;
  public:
    typedef typename TA::TELEM TELEM;
    typedef typename TA::TSCAL TSCAL;
    static constexpr bool IsLinear() { return false; }
    
    INLINE RowsArrayExpr (const TA & aa, FlatArray<int> arows) : a(aa), rows(arows) { ; }

    INLINE auto Height() const { return rows.Size(); }
    INLINE auto Width() const { return a.Width(); }

    INLINE auto operator() (size_t i, size_t j) const-> decltype(a(rows[i],j)) { return a(rows[i], j); }
    INLINE auto operator() (size_t i) const-> decltype(a(rows[i])) { return a(rows[i]); }

    INLINE auto Row (size_t i) const { return a.Row(rows[i]); }
    INLINE auto View() const { return RowsArrayExpr(a, rows); }
    INLINE auto ViewRW() { return RowsArrayExpr(a, rows); }
    INLINE auto Shape() const
    {
      typedef decltype(a.Shape()) TASHAPE;
      if constexpr (tuple_size<TASHAPE>() == 1)
                     return tuple<size_t> (rows.Size());
      else
        return tuple<size_t,size_t> (rows.Size(), a.Width());
    }
    
    template<typename TB>
    INLINE const RowsArrayExpr & operator= (const Expr<TB> & m) 
    {
      MatExpr<RowsArrayExpr<TA> >::operator= (m);
      return *this;
    }

    INLINE const RowsArrayExpr & operator= (const RowsArrayExpr & m) 
    {
      MatExpr<RowsArrayExpr<TA> >::operator= (m);
      return *this;
    }
  };
  

  /**
     ColsArray
  */
  template <class TA> class ColsArrayExpr : public MatExpr<ColsArrayExpr<TA> >
  {
    const TA & a;
    FlatArray<int> cols;
  public:
    typedef typename TA::TELEM TELEM;
    static constexpr bool IsLinear() { return false; } 
    typedef typename TA::TSCAL TSCAL;

    INLINE ColsArrayExpr (const TA & aa, FlatArray<int> acols) : a(aa), cols(acols) { ; }

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return cols.Size(); }

    INLINE auto operator() (size_t i, size_t j) const -> decltype(a(i, cols[j]))  { return a(i, cols[j]); }
    INLINE auto operator() (size_t i) const -> decltype(a(i, cols[0]))  { return a(i, cols[0]); }

    INLINE auto View() const { return *this; }
    INLINE auto ViewRW() { return *this; }
    INLINE tuple<size_t,size_t> Shape() const { return { a.Height(), cols.Size() }; }
    
    template<typename TB>
    INLINE const ColsArrayExpr & operator= (const Expr<TB> & m) 
    {
      MatExpr<ColsArrayExpr<TA> >::operator= (m);
      return *this;
    }

  };
  



  /* ************************* Conjugate *********************** */ 


  INLINE double Conj (double a)
  {
    return a;
  }

  INLINE Complex Conj (Complex a)
  {
    return conj(a);
  }

  /**
     Conjugate of Matrix-expr
  */
  // Attention NOT transpose !!! elemwise conjugate !!! 
  template <class TA> class ConjExpr : public Expr<ConjExpr<TA> >
  {
    TA a;
  public:
    // typedef typename TA::TELEM TELEM;
    // typedef typename TA::TSCAL TSCAL;

    INLINE ConjExpr (TA aa) : a(aa) { ; }

    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return a.Width(); }
    
    INLINE auto View() const { return *this; }
    INLINE auto Shape() const { return a.Shape(); }
    
    INLINE auto operator() (size_t i, size_t j) const { return Conj(a(i,j)); }
    INLINE auto operator() (size_t i) const { return Conj(a(i)); }

    static constexpr bool IsLinear() { return TA::IsLinear(); }         
  };


  /// Conjugate
  template <typename TA>
  INLINE auto Conj (const Expr<TA> & a)
  {
    if constexpr (IsComplex<decltype(a.Spec()(0))>())
      return ConjExpr (a.View());
    else
      return a.View();
  }



  /* ************************* Truncate ************************* */

  INLINE double Truncate (double v, double eps = 1e-12)
  {
    return (fabs(v) < eps) ? 0 : v;
  }

  INLINE Complex Truncate (Complex v, double eps = 1e-12)
  {
    return Complex(Truncate(Real(v), eps), Truncate(Imag(v),eps));
  }
  
  template <class TA> class TruncateExpr : public Expr<TruncateExpr<TA> >
  {
    TA a;
    double eps;
  public:
    INLINE TruncateExpr (TA aa, double aeps) : a(aa), eps(aeps) { ; }

    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return a.Width(); }
 
    INLINE auto operator() (size_t i, size_t j) const { return Truncate(a(i,j), eps); }
    INLINE auto operator() (size_t i) const { return Truncate(a(i), eps); }
    INLINE auto View() const { return *this; }

    static constexpr bool IsLinear() { return TA::IsLinear(); }     
  };

  /// Conjugate
  template <typename TA>
  INLINE TruncateExpr<TA>
  Truncate (const Expr<TA> & a, double eps = 1e-12)
  {
    return TruncateExpr<TA> (a.View(), eps);
  }
  

  /* ************************* InnerProduct ********************** */

  template <typename TA, typename TB,
            typename enable_if<IsScalar<TA>(),int>::type = 0,
            typename enable_if<IsScalar<TB>(),int>::type = 0>
  INLINE auto InnerProduct (TA a, TB b) { return a*b; } 



  
  /**
     Inner product
  */

  template <class TA, class TB>
  INLINE auto InnerProduct (const Expr<TA> & a, const Expr<TB> & b)
    -> decltype (InnerProduct(a.Spec()(0), b.Spec()(0)))
  {
    auto h = CombinedSize(a.Height(), b.Height());
    auto w = CombinedSize(a.Width(), b.Width());
      
    if (h*w == 0) return 0; 

    auto sum = InnerProduct (a.Spec()(0), b.Spec()(0));
    for (size_t i = 1; i < h*w; i++)
      sum += InnerProduct (a.Spec()(i), b.Spec()(i));
    return sum;
  }


  /* ************************* OuterProduct ********************** */
  
  template <class TA, class TB>
  class OuterProductExpr : public Expr<OuterProductExpr<TA,TB>>
  {
    TA a;
    TB b;
  public:
    OuterProductExpr (TA aa, TB ab) : a(aa), b(ab) { ; }

    // INLINE auto operator() (size_t i) const { return a[i] * b(i); }  
    INLINE auto operator() (size_t i, size_t j) const { return a[i] * b[j]; }

    INLINE auto View() const { return *this; }
    INLINE auto Shape() const
    {
      return tuple<size_t,size_t> (a.Size(), b.Width());
    }
    
    INLINE const auto A() const { return a; }
    INLINE const auto B() const { return b; }
    INLINE auto Height() const { return a.Size(); }
    INLINE auto Width() const { return b.Size(); }

    static constexpr bool IsLinear() { return false; }         
  };


  template <typename TA, typename TB>
  INLINE auto OuterProduct (const Expr<TA> & a, const Expr<TB> & b)
  {
    return OuterProductExpr (a.View(), b.View());
  }



  /* **************************** Trace **************************** */


  /**
     Calculates the trace of a matrix expression.
  */
  template <class TA>
  INLINE auto Trace (const Expr<TA> & a) // -> decltype (a.Spec()(0,0))
  {
    typedef decltype( RemoveConst(a.Spec()(0,0)) ) TRES;
    TRES sum = 0;    
    for (size_t i = 0; i < Height(a); i++)
      sum += a.Spec()(i,i);
    return sum;
  }

  /* **************************** L2Norm **************************** */

  /// Euclidean norm squared
  INLINE double L2Norm2 (double v)
  {
    return v*v;
  }

  INLINE double L2Norm2 (Complex v)
  {
    return v.real()*v.real()+v.imag()*v.imag();
  }

  template <class TA>
  INLINE auto L2Norm2 (const Expr<TA> & v) 
  {
    decltype(L2Norm2(v.Spec()(0))) sum = 0.0;
    if (TA::IsLinear())
      for (size_t i = 0; i < v.Height()*v.Width(); i++)
	sum += L2Norm2 (v.Spec()(i));  
    else
      for (size_t i = 0; i < v.Height(); i++)
	for (size_t j = 0; j < v.Width(); j++)
	  sum += L2Norm2 (v.Spec()(i,j));  
    
    return sum;
  }

  /**
     Calculates the Euclidean norm
  */
  template <class TA>
  INLINE auto L2Norm (const Expr<TA> & v) // -> decltype(L2Norm2(v))
  {
    return sqrt (L2Norm2(v));
  }



  /* **************************** MaxNorm **************************** */

  /// Euclidean norm squared
  INLINE double MaxNorm (double v)
  {
    return fabs(v);
  }

  INLINE double MaxNorm (Complex v)
  {
    return fabs(v);
  }


  template <class TA>
  INLINE double MaxNorm (const Expr<TA> & v)
  {
    double sum = 0;

    if constexpr (TA::IsLinear())
      for (size_t i = 0; i < v.Height()*v.Width(); i++)
	sum = max(sum, MaxNorm( v.Spec()(i)) );  
    else
      for (size_t i = 0; i < v.Height(); i++)
	for (size_t j = 0; j < v.Width(); j++)
	  sum = max(sum, MaxNorm ( v.Spec()(i,j)) );  
    
    return sum;
  }


  /* *************************** Output ****************************** */



  /// Print matrix-expr
  template<typename T>
  ostream & operator<< (ostream & s, const Expr<T> & v)
  { 
    int width = s.width();
    if (width == 0) width = 8;
    s.width(0);
    for (size_t i = 0; i < v.Height(); i++)
      {
	for (size_t j = 0 ; j < v.Width(); j++)
	  s << " " << setw(width-1) << v.Spec()(i,j);
	s << endl;
      }
    return s;
  }




  
  template <typename TA,
            enable_if_t<IsScalar<TA>(),bool> = true>
  INLINE auto Inv (TA val) { return 1.0/val; } 
  


  /* ********************** Determinant *********************** */


  /**
     Calculates the determinant of a Matrix.
  */

  template <class T>
  inline typename T::TELEM Det (const MatExpr<T> & m)
  {
    const T & sm = m.Spec();
    switch (sm.Height())
      {
      case 1: 
	{
	  return sm(0,0); 
	}
      case 2:
	{
	  return ( sm(0,0)*sm(1,1) - sm(0,1)*sm(1,0) ); 
	}
      case 3:
	{
	  return 
	    sm(0) * (sm(4) * sm(8) - sm(5) * sm(7)) +
	    sm(1) * (sm(5) * sm(6) - sm(3) * sm(8)) +
	    sm(2) * (sm(3) * sm(7) - sm(4) * sm(6));
	}
      default:
	{
	  cerr << "general det not implemented" << endl;
	}
      }

    return typename T::TELEM (0);  
  }

}

#endif
