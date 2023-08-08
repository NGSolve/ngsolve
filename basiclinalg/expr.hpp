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
#include <core/hashtable.hpp>   // for INT

#include "complex_wrapper.hpp"

namespace ngbla
{
  using namespace std;
  using namespace ngcore;
  using namespace ngstd;

  
  enum ORDERING { ColMajor, RowMajor };

  template <typename T = double, ORDERING ORD = RowMajor> class FlatMatrix;
  template <typename T = double, ORDERING ORD = RowMajor> class Matrix;

  template <int H, int W, typename T> class Mat;
  template <int H, typename T> class DiagMat;
  template <int S, typename T> class Vec;

  template <typename T = double, ORDERING ORD = RowMajor> class SliceMatrix;
  template <typename T, ORDERING ORD> class BareSliceMatrix;
  template <class T = double> class FlatVector;
  template <class T = double> class BareVector;  
  template <class T = double> class SliceVector;
  template <class T = double> class BareSliceVector;

  
  template <typename T, typename TELEM=typename T::TELEM>
  constexpr bool IsConvertibleToSliceMatrix ()
  {
    return is_convertible_v<T,SliceMatrix<TELEM, RowMajor>> ||
      is_convertible_v<T,SliceMatrix<TELEM, ColMajor>>;
  }

  template <typename T, typename TELEM=typename T::TELEM>
  constexpr bool IsConvertibleToBareSliceMatrix ()
  {
    return is_convertible_v<T,BareSliceMatrix<TELEM, RowMajor>> ||
      is_convertible_v<T,BareSliceMatrix<TELEM, ColMajor>>;
  }



  namespace detail {
    template <typename T>
    struct test_conv_flatvector {
      template<typename T2>
      static constexpr auto check(T2*) -> decltype(FlatVector(T2()));
      template<typename>
      static constexpr std::false_type check(...);
      
      using type = decltype(check<T>(nullptr));
      static constexpr bool value = !std::is_same<type, std::false_type>::value;
    };
  }
  
  template <typename T>
  constexpr bool IsConvertibleToFlatVector()
  {
    return detail::test_conv_flatvector<T>::value;
  }


  namespace detail {
    template <typename T>
    struct test_conv_slicevector {
      template<typename T2>
      static constexpr auto check(T2*) -> decltype(SliceVector(T2()));
      template<typename>
      static constexpr std::false_type check(...);
      
      using type = decltype(check<T>(nullptr));
      static constexpr bool value = !std::is_same<type, std::false_type>::value;
    };
  }
  
  template <typename T>
  constexpr bool IsConvertibleToSliceVector()
  {
    return detail::test_conv_slicevector<T>::value;
  }

  

  namespace detail {
    template <typename T>
    struct test_conv_bareslicevector {
      template<typename T2>
      static constexpr auto check(T2*) -> decltype(BareSliceVector(T2()));
      template<typename>
      static constexpr std::false_type check(...);
      
      using type = decltype(check<T>(nullptr));
      static constexpr bool value = !std::is_same<type, std::false_type>::value;
    };
  }
  
  template <typename T>
  constexpr bool IsConvertibleToBareSliceVector()
  {
    return detail::test_conv_bareslicevector<T>::value;
  }

  



  
  template <bool ADD, bool POS, ORDERING orda, ORDERING ordb>
  void NgGEMM (SliceMatrix<double,orda> a, SliceMatrix<double, ordb> b, SliceMatrix<double> c);

  template <bool ADD, bool POS, ORDERING orda, ORDERING ordb>
  void NgGEMM (SliceMatrix<double,orda> a, SliceMatrix<double, ordb> b, SliceMatrix<double,ColMajor> c);
  
  
  

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

  /*
  template <typename TVEC> 
  inline typename TVEC::TELEM & Access (TVEC & vec, int nr)
  {
    return vec(nr);
  }

  template <typename TVEC> 
  inline typename TVEC::TELEM Access (const TVEC & vec, int nr)
  {
    return vec(nr);
  }

  inline double & Access (double & vec, int nr)
  {
    return vec;
  }

  inline double Access (const double & vec, int nr)
  {
    return vec;
  }

  inline Complex & Access (Complex & vec, int nr)
  {
    return vec;
  }

  inline Complex Access (const Complex & vec, int nr)
  {
    return vec;
  }
  */
  
  /*
  template <typename TM> 
  inline typename TM::TELEM & Access (TM & mat, int i, int j)
  {
    return mat(i,j);
  }
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

  
  /*
  inline double & Access (double & mat, int i, int j)
  {
    return mat;
  }
  */

  /*
  inline double Access (const double & mat, int i, int j)
  {
    return mat;
  }
  */
  /*
  inline Complex & Access (Complex & mat, int i, int j)
  {
    return mat;
  }
  */
  /*
  inline Complex Access (const Complex & mat, int i, int j)
  {
    return mat;
  }
  */
  



  template <typename T> struct is_scalar_type;
    
  template<> struct is_scalar_type<int> { static constexpr bool value = true; };  
  template<> struct is_scalar_type<double> { static constexpr bool value = true; };
  template<> struct is_scalar_type<Complex> { static constexpr bool value = true; };
  


  /**
     Trait to obtain vector and scalar types for given matrix types.
     Is specified for double, Complex, AutoDiff<doube>, AutoDiff<Complex>
  */




  /*
  template <class T>
  class mat_traits
  {
  public:
    /// matrix element
    typedef typename T::TELEM TELEM;
    /// field of matrix element
    typedef typename T::TSCAL TSCAL;
    /// type of column vector
    typedef typename T::TV_COL TV_COL;
    /// type of row vector
    typedef typename T::TV_ROW TV_ROW;
    /// matrix height
    enum { HEIGHT = T::HEIGHT };
    /// matrix with
    enum { WIDTH  = T::WIDTH  };
    ///
    enum { IS_COMPLEX = mat_traits<TSCAL>::IS_COMPLEX };
  };
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


  /*
  template <int D, typename SCAL>
  class mat_traits<AutoDiff<D,SCAL> >
  {
  public:
    typedef AutoDiff<D,SCAL> TELEM;
    typedef AutoDiff<D,SCAL> TSCAL;
    typedef AutoDiff<D,SCAL> TV_COL;
    typedef AutoDiff<D,SCAL> TV_ROW;
    // static constexpr int HEIGHT = 1;
    // static constexpr int WIDTH = 1;
    // enum { IS_COMPLEX = mat_traits<SCAL>::IS_COMPLEX };
  };
  */



  /// Height of matrix
  template <class TM> 
  inline size_t Height (const TM & m)
  {
    return m.Height();
  }

  /// Width of matrix
  template <class TM> 
  inline size_t Width (const TM & m)
  {
    return m.Width();
  }

  template <> inline constexpr size_t Height<double> (const double&) { return 1; }
  template <> inline constexpr size_t Height<Complex> (const Complex&) { return 1; }
  template <> inline constexpr size_t Width<double> (const double&) { return 1; }
  template <> inline constexpr size_t Width<Complex> (const Complex&) { return 1; }

  /*
  template <class TM> 
  inline constexpr size_t Height () { return Height(TM()); }
  template <class TM> 
  inline constexpr size_t Width () { return Width(TM()); }
  */

  template <class TM> 
  inline constexpr size_t Height () { return TM::Height(); }
  template <class TM> 
  inline constexpr size_t Width () { return TM::Width(); }

  template <> inline constexpr size_t Height<double> () { return 1; }
  template <> inline constexpr size_t Height<Complex> () { return 1; }
  template <> inline constexpr size_t Width<double> () { return 1; }
  template <> inline constexpr size_t Width<Complex> () { return 1; }

  
  template <class TM> 
  inline constexpr bool IsComplex () { return IsComplex<typename mat_traits<TM>::TSCAL>(); }
  template <> inline constexpr bool IsComplex<double> () { return false; }
  template <> inline constexpr bool IsComplex<Complex> () { return true; }  

  
  

  /// Complex to double assignment called
  class Complex2RealException : public Exception
  {
  public:
    Complex2RealException()
      : Exception("Assignment of Complex 2 Real") { ; }
  };



  template <typename TO>
  inline TO ConvertTo (double f)
  {
    return TO(f);
  }

  template <typename TO>
  inline TO ConvertTo (Complex f)
  {
    return TO(f);
  }

  /*
  template <typename TO>
  inline TO ConvertTo (const AutoDiff<1, Complex> & f)
  {
    return TO(f);
  }
  */

  template <>
  inline double ConvertTo (Complex f)
  {
    throw Complex2RealException();
  }


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

  struct undefined_size { };
  inline auto operator/ (undefined_size ud, size_t i) { return undefined_size(); }
  inline auto operator- (undefined_size ud, size_t i) { return undefined_size(); }
  
  INLINE auto CombinedSize(undefined_size s1, undefined_size s2) { return undefined_size(); }
  INLINE auto CombinedSize(undefined_size s1, size_t s2) { return s2; }  
  INLINE auto CombinedSize(size_t s1, undefined_size s2) { return s1; }  
  INLINE auto CombinedSize(size_t s1, size_t s2) { return s1; }
  template <int S1> INLINE auto CombinedSize(IC<S1> s1, undefined_size s2) { return s1; }  
  template <int S1> INLINE auto CombinedSize(IC<S1> s1, size_t s2) { return s1; }  
  template <int S1, int S2> INLINE auto CombinedSize(IC<S1> s1, IC<S2> s2) { return s1; }  
  template <int S2> INLINE auto CombinedSize(undefined_size s1, IC<S2> s2) { return s2; }  
  template <int S2> INLINE auto CombinedSize(size_t s1, IC<S2> s2) { return s2; }  

  template <typename T1, typename T2>
  INLINE auto CombinedSize(tuple<T1> tup1, tuple<T2> tup2)
  { return tuple(CombinedSize(get<0>(tup1), get<0>(tup2))); }

  template <typename T11, typename T12, typename T21, typename T22>
  INLINE auto CombinedSize(tuple<T11,T12> tup1, tuple<T21,T22> tup2)
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
    // INLINE auto Shape() const { return static_cast<const T&> (*this).Shape(); }
    INLINE auto Shape() const { return Spec().T::Shape(); }


    /// height
    INLINE auto Height() const { return Spec().T::Height(); }
    INLINE auto Width() const { return Spec().T::Width(); }

    // auto Shape() const { return std::tuple(Height(), Width()); }

    // INLINE auto operator() (int i) const -> decltype (this->Spec()(i)) { return this->Spec()(i); }
    // INLINE auto operator() (int i, int j) const -> decltype (this->Spec()(i,j)) { return this->Spec()(i,j); }

    void Dump (ostream & ost) const { Spec().T::Dump(ost); }


    INLINE auto Row (size_t r) const
    {
      return RowExpr<const T> (static_cast<const T&> (*this), r);
    }

    INLINE auto Col (size_t r) const
    {
      return ColExpr<const T> (static_cast<const T&> (*this), r);
    }


    INLINE SubMatrixExpr<T>
    Rows (size_t first, size_t next)
    { 
      return SubMatrixExpr<T> (static_cast<T&> (*this), first, 0, next-first, Width()); 
    }

    INLINE SubMatrixExpr<T>
    Cols (size_t first, size_t next) 
    { 
      return SubMatrixExpr<T> (static_cast<T&> (*this), 0, first, Height(), next-first);
    }

    INLINE SubMatrixExpr<T>
    Rows (IntRange range) 
    { 
      return Rows (range.First(), range.Next());
    }

    INLINE SubMatrixExpr<T>
    Cols (IntRange range) 
    { 
      return Cols (range.First(), range.Next());
    }


    INLINE RowsArrayExpr<T>
    Rows (FlatArray<int> rows) 
    { 
      return RowsArrayExpr<T> (static_cast<const T&> (*this), rows); 
    }

    INLINE ColsArrayExpr<T>
    Cols (FlatArray<int> cols) 
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
    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }
    
    auto View() const { return *this; }
    auto Shape() const { return a.Shape(); }
    
    enum { IS_LINEAR = T::IS_LINEAR };
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
    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }
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
      NETGEN_CHECK_RANGE(self.Height(), v.Height(), v.Height()+1);
      NETGEN_CHECK_RANGE(self.Width(), v.Width(), v.Width()+1);
      NETGEN_CHECK_SHAPE(self.Spec(), v);

      /*
      if constexpr (std::is_same_v<TOP,typename MatExpr<T>::As> && 
                    is_convertible_v<TB,FlatVector<typename T::TELEM>> && 
                    is_convertible_v<T,FlatVector<typename T::TELEM>>)
        {
          CopyVector(FlatVector<typename T::TELEM>(v.Spec()), FlatVector<typename T::TELEM>(self.Spec()));
          return self.Spec();
        }
      */
      if constexpr (std::is_same_v<TOP,typename MatExpr<T>::As> && 
                    IsConvertibleToFlatVector<TB>() && 
                    IsConvertibleToFlatVector<T>())
        {
          CopyVector(BareVector(FlatVector(v.Spec())), FlatVector(self.Spec()));
          return self.Spec();
        }

      /*
      else if constexpr (std::is_same_v<TOP,typename MatExpr<T>::As> && 
                         is_convertible_v<TB,SliceVector<typename T::TELEM>> && 
                         is_convertible_v<T,SliceVector<typename T::TELEM>>)
        {
          CopyVector(SliceVector<typename T::TELEM>(v.Spec()), SliceVector<typename T::TELEM>(self.Spec()));
          return self.Spec();
        }
      */
      else if constexpr (std::is_same_v<TOP,typename MatExpr<T>::As> &&
                         IsConvertibleToSliceVector<TB>() && 
                         IsConvertibleToSliceVector<T>())
        {
          CopyVector(BareSliceVector(SliceVector(v.Spec())), SliceVector(self.Spec()));
          return self.Spec();
        }

      auto src = v.View();
      decltype(auto) dest = self.ViewRW();
      
      if (T::COL_MAJOR)
        {
	  size_t h = self.Height();
	  size_t w = self.Width();

          if (h > 0)
            for (size_t j = 0; j < w; j++)
              for (size_t i = 0; i < h; i++)
                // TOP()(Spec()(i,j), v.Spec()(i,j));
                // TOP()(Spec()(i,j), v.View()(i,j));
                TOP()(dest(i,j), src(i,j));
          return self.Spec();
        }


      if (TB::IS_LINEAR)
	{
	  if (T::IS_LINEAR)
	    {
	      auto hw = self.Height() * self.Width();
              for (auto i : Range(hw))  // int i = 0; i < hw; i++)
                // TOP()(Spec()(i),v.Spec()(i));
                TOP()(dest(i), src(i));
	    }
	  else
	    {
	      size_t h = self.Height();
	      size_t w = self.Width();
              if (w > 0)
                for (size_t i = 0, k = 0; i < h; i++)
                  for (size_t j = 0; j < w; j++, k++)
                    // TOP() (Spec()(i,j), v.Spec()(k));
                    // TOP() (Spec()(i,j), v.View()(k));
                    TOP() (dest(i,j), src(k));
	    }
	}
      else
	{
	  size_t h = self.Height();
	  size_t w = self.Width();
          if (w > 0)
            {
              if (T::IS_LINEAR)
                for (size_t i = 0, k = 0; i < h; i++)
                  for (size_t j = 0; j < w; j++, k++)
                    // TOP() (Spec()(k), v.Spec()(i,j));
                    // TOP() (Spec()(k), v.View()(i,j));
                    TOP() (dest(k), src(i,j));                    
              else
                {
                  for (size_t i = 0; i < h; i++)
                    for (size_t j = 0; j < w; j++)
                      {
                        // TOP() (Spec()(i,j), v.View()(i,j));
                        TOP() (dest(i,j), src(i,j));                        
                      }
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

    enum { IS_LINEAR = 1 };  // row-major continuous storage (dist=width)
    enum { COL_MAJOR = 0 };  // matrix is stored col-major

    void Dump (ostream & ost) const { ost << "Matrix"; }



    template<typename TOP, typename TB>
    INLINE T & Assign (const Expr<TB> & v)
    {
      return assign_trait<TOP, T, TB>::Assign (*this, v);
    }


    class As 
    {
    public:
      template <typename T1, typename T2> 
      INLINE void operator() (T1 && v1, const T2 & v2) { v1 = v2; }
    };
    class AsAdd 
    {
    public:
      template <typename T1, typename T2> 
      INLINE void operator() (T1 && v1, const T2 & v2) { v1 += v2; }
    };
    class AsSub 
    {
    public:
      template <typename T1, typename T2> 
      INLINE void operator() (T1 && v1, const T2 & v2) { v1 -= v2; }
    };
	


    
    /*
    template <typename OP, typename TA,
              typename enable_if<std::is_same<OP,As>::value,int>::type = 0,
              typename enable_if<is_convertible<TA,SliceVector<typename TA::TELEM>>::value,int>::type = 0>
    INLINE T & Assign (const Expr<TA> & src) 
    {
      CopyVector(src.Spec(), this->Spec());
    }
    */


    
    template <typename OP, typename TA, typename TB,
              enable_if_t<IsConvertibleToSliceMatrix<TA,double>(),bool> = true,
              enable_if_t<IsConvertibleToSliceMatrix<TB,double>(),bool> = true,
              enable_if_t<IsConvertibleToSliceMatrix<typename pair<T,TB>::first_type,double>(),bool> = true>
    INLINE T & Assign (const Expr<MultExpr<TA, TB>> & prod) 
    {
      constexpr bool ADD = std::is_same<OP,AsAdd>::value || std::is_same<OP,AsSub>::value;
      constexpr bool POS = std::is_same<OP,As>::value || std::is_same<OP,AsAdd>::value;
      
      NgGEMM<ADD,POS> (make_SliceMatrix(prod.View().A()),
                       make_SliceMatrix(prod.View().B()),
                       make_SliceMatrix(Spec()));
      return Spec();
    }


    template <typename OP, typename TA, typename TB,
              typename enable_if<IsConvertibleToSliceMatrix<TA,double>(),int>::type = 0,
              typename enable_if<IsConvertibleToSliceMatrix<TB,double>(),int>::type = 0,
              typename enable_if<IsConvertibleToSliceMatrix<typename pair<T,TB>::first_type,double>(),int>::type = 0>    
    INLINE T & Assign (const Expr<MultExpr<MinusExpr<TA>, TB>> & prod) 
    {
      constexpr bool ADD = std::is_same<OP,AsAdd>::value || std::is_same<OP,AsSub>::value;
      constexpr bool POS = std::is_same<OP,As>::value || std::is_same<OP,AsAdd>::value;
      
      NgGEMM<ADD,!POS> (make_SliceMatrix(prod.View().A().A()),
                        make_SliceMatrix(prod.View().B()),
                        make_SliceMatrix(Spec()));
      return Spec();
    }

 

    
    // x += s*y
    template <typename OP, typename TA, 
              enable_if_t<std::is_same<OP,AsAdd>::value,bool> = true,
              enable_if_t<is_convertible_v<TA,SliceVector<double>>,bool> = true,
              enable_if_t<is_convertible<typename pair<T,TA>::first_type,SliceVector<double>>::value,bool> = true>
    INLINE T & Assign (const Expr<ScaleExpr<TA,double>> & scaled)
    {
      AddVector (scaled.View().S(),
                 SliceVector<typename T::TELEM>(scaled.View().A()),
                 SliceVector<typename T::TELEM>(this->Spec()));
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
    
    // rank 1 update
    template <typename OP, typename TA, typename TB,
              typename enable_if<is_convertible_v<TA,FlatVector<double>>,int>::type = 0,
              typename enable_if<is_convertible_v<TB,FlatVector<double>>,int>::type = 0,
              typename enable_if<IsConvertibleToSliceMatrix<typename pair<T,TB>::first_type,double>(),int>::type = 0>
    INLINE T & Assign (const Expr<MultExpr<TA, TransExpr<TB>>> & prod)
    {
      constexpr bool ADD = std::is_same<OP,AsAdd>::value || std::is_same<OP,AsSub>::value;
      constexpr bool POS = std::is_same<OP,As>::value || std::is_same<OP,AsAdd>::value;

      auto veca = prod.Spec().A();
      auto mata = FlatMatrix<typename TA::TELEM>(veca.Height(), 1, veca.Data());
      auto vecb = prod.Spec().B().A();
      auto matb = FlatMatrix<typename TB::TELEM>(1, vecb.Height(), vecb.Data());
      
      NgGEMM<ADD,POS> (SliceMatrix<typename TA::TELEM>(mata),
                       SliceMatrix<typename TB::TELEM>(matb), Spec());
      return Spec();
    }


    
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


    /*
    template <typename TA, typename TB,
              typename enable_if<IsConvertibleToSliceMatrix<TA,double>(),int>::type = 0,
              typename enable_if<IsConvertibleToSliceMatrix<TB,double>(),int>::type = 0>
    INLINE T & operator= (const Expr<MultExpr<TA, TB>> & prod) 
    {
      cout << "using fast" << endl;
      NgGEMM<false,true> (make_SliceMatrix(prod.Spec().A()),
                          make_SliceMatrix(prod.Spec().B()),
                          make_SliceMatrix(Spec()));
      return Spec();
    }
    */





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
      if (T::IS_LINEAR)
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











  /**
     The base class for matrices.
     Constant-Means-Constat-Pointer
     matrix-values may be changed by const methods
  */
  template <class T>
  class CMCPMatExpr : public MatExpr<T>
  {
  public:
    // int Height() const { return Spec().T::Height(); }
    // int Width() const { return Spec().T::Width(); }

    INLINE CMCPMatExpr () { ; }

    using MatExpr<T>::Spec;
    using MatExpr<T>::Height;
    using MatExpr<T>::Width;

    // T & Spec() { return static_cast<T&> (*this); }
    // const T & Spec() const { return static_cast<const T&> (*this); }
    // enum { IS_LINEAR = 1 };

    template<typename TB>
    INLINE const T & operator= (const Expr<TB> & v) const
    {
      const_cast<CMCPMatExpr*> (this) -> MatExpr<T>::operator= (v);
      return Spec();
    }

    INLINE auto ViewRW() { return this->View(); }    

    
    template<typename TB>
    INLINE const T & operator+= (const Expr<TB> & v) const
    {
      const_cast<CMCPMatExpr*> (this) -> MatExpr<T>::operator+= (v);
      return Spec();
    }

    template<typename TB>
    INLINE const T & operator+= (const Expr<SymExpr<TB> > & v) const
    {
      const_cast<CMCPMatExpr*> (this) -> MatExpr<T>::operator+= (v);
      return Spec();
    }

    template<typename TB>
    INLINE const T & operator-= (const Expr<TB> & v) const
    {
      const_cast<CMCPMatExpr*> (this) -> MatExpr<T>::operator-= (v);
      return Spec();
    }

    template <class SCAL2>
    INLINE const T & operator*= (SCAL2 s) const
    {
      const_cast<CMCPMatExpr*> (this) -> MatExpr<T>::operator*= (s);
      return Spec();
    }

    template <class SCAL2>
    INLINE const T & operator/= (SCAL2 s) const 
    {
      return (*this) *= (1.0/s);
    }

    SubMatrixExpr<const T>
    INLINE Rows (size_t first, size_t next) const
    { 
      return SubMatrixExpr<const T> (static_cast<const T&> (*this), first, 0, next-first, Width()); 
    }

    SubMatrixExpr<const T>
    INLINE Cols (size_t first, size_t next) const
    { 
      return SubMatrixExpr<const T> (static_cast<const T&> (*this), 0, first, Height(), next-first);
    }

    SubMatrixExpr<const T>
    INLINE Rows (IntRange range) const
    { 
      return Rows (range.First(), range.Next());
    }

    SubMatrixExpr<const T>
    INLINE Cols (IntRange range) const
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

    enum { IS_LINEAR = TA::IS_LINEAR && TB::IS_LINEAR };
    
    INLINE SumExpr (TA aa, TB ab) : a(aa), b(ab) { ; }

    INLINE auto operator() (size_t i) const { return a(i)+b(i); }
    INLINE auto operator() (size_t i, size_t j) const { return a(i,j)+b(i,j); }

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

    enum { IS_LINEAR = TA::IS_LINEAR && TB::IS_LINEAR };
    
    INLINE SubExpr (TA aa, TB ab) : a(aa), b(ab) { ; }

    INLINE auto operator() (size_t i) const { return a(i)-b(i); }
    INLINE auto operator() (size_t i, size_t j) const { return a(i,j)-b(i,j); }

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

    INLINE auto operator() (size_t i) const { return -a(i); }
    INLINE auto operator() (size_t i, size_t j) const { return -a(i,j); }
    INLINE auto View() const { return MinusExpr(a); }
    INLINE auto Shape() const { return a.Shape(); }
    
    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return a.Width(); }
    INLINE TA A() const { return a; }
    enum { IS_LINEAR = TA::IS_LINEAR };
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
    enum { IS_LINEAR = TA::IS_LINEAR && TB::IS_LINEAR };

    INLINE PW_Mult_Expr (TA aa, TB ab) : a(aa), b(ab) { ; }

    INLINE auto operator() (size_t i) const { return a(i)*b(i); }
    INLINE auto operator() (size_t i, size_t j) const { return a(i,j)*b(i,j); }

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }
    
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
    enum { IS_LINEAR = TA::IS_LINEAR };

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
    enum { IS_LINEAR = TA::IS_LINEAR };

    INLINE ScaleExpr (TA aa, TS as) : a(aa), s(as) { ; }

    INLINE auto operator() (size_t i) const { return s * a(i); }
    INLINE auto operator() (size_t i, size_t j) const { return s * a(i,j); }

    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return a.Width(); }

    INLINE auto View() const { return *this; }
    INLINE auto Shape() const { return a.Shape(); }         
    
    INLINE TA A() const { return a; }
    INLINE TS S() const { return s; }
    
    void Dump (ostream & ost) const
    { ost << "Scale, s=" << s << " * "; a.Dump(ost);  }
  };


  /*
  template <typename TA>
  INLINE auto operator* (double b, const Expr<TA> & a)
  {
    return ScaleExpr (a.View(), b);
  }

  template <typename TA>
  INLINE auto operator* (Complex b, const Expr<TA> & a)
  {
    return ScaleExpr (a.View(), b);
  }
  
  template <int D, typename TAD, typename TA>
  INLINE auto operator* (const AutoDiff<D,TAD> & b, const Expr<TA> & a)
  {
    return ScaleExpr (a.View(), b );
  }

  template <int N, typename TA>
  INLINE auto operator* (SIMD<double,N> b, const Expr<TA> & a)
  {
    return ScaleExpr (a.View(), b);
  }
  
  template <int N, typename TA>
  INLINE auto operator* (SIMD<Complex,N> b, const Expr<TA> & a)
  {
    return ScaleExpr (a.View(), b);
  }
  */

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

    INLINE auto View() const { return MultExpr(a,b); }
    INLINE auto Shape() const
    {
      typedef decltype(b.Shape()) TBSHAPE;
      if constexpr (tuple_size<TBSHAPE>() == 1)
        return tuple(get<0>(a.Shape()));
        // return tuple<size_t> (a.Height());
      else
        return tuple(get<0>(a.Shape()), get<1>(b.Shape()));
        // return tuple<size_t,size_t> (a.Height(), b.Width());
    }
    
    INLINE TA A() const { return a; }
    INLINE TB B() const { return b; }
    INLINE auto Height() const { return a.Height(); }
    INLINE auto Width() const { return b.Width(); }
    enum { IS_LINEAR = 0 };
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
    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return b.Width(); }
    enum { IS_LINEAR = 0 };
  };


  template <typename TA, typename TB>
  INLINE auto operator* (const Expr<TA> & a, const Expr<TB> & b)
  {
    return MultExpr (a.View(), b.View());
  }


  /* ************************** Trans *************************** */

  /*
  INLINE double Trans (double a) { return a; }
  INLINE Complex Trans (Complex a) { return a; }
  template<int D, typename TAD>
  INLINE AutoDiff<D,TAD> Trans (const AutoDiff<D,TAD> & a) { return a; }
  */
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

    INLINE size_t Height() const { return a.Width(); }
    INLINE size_t Width() const { return a.Height(); }

    INLINE auto operator() (size_t i, size_t j) const { return Trans (a(j,i)); }
    INLINE auto operator() (size_t i) const { return Trans(a(0,0)); }
    // auto Row (int i) const -> decltype (a.Col(i)) { return a.Col(i); }
    // auto Col (int i) const -> decltype (a.Row(i)) { return a.Row(i); }

    INLINE auto View() const { return *this; }
    INLINE tuple<size_t,size_t> Shape() const { return { a.Width(), a.Height() }; }
    enum { IS_LINEAR = 0 };

    INLINE const TA & A() const { return a; }
  };


  /// Transpose 
  template <typename TA>
  INLINE auto Trans (const Expr<TA> & a)
  {
    return TransExpr (a.View());
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
    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }
    INLINE auto View() const { return *this; }
    INLINE auto Shape() const { return a.Shape(); }
    enum { IS_LINEAR = TA::IS_LINEAR };
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
    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }
    INLINE auto View() const { return *this; }
    INLINE auto Shape() const { return a.Shape(); }    
    enum { IS_LINEAR = TA::IS_LINEAR };
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
    TA & a;
    size_t first_row, first_col;
    size_t height, width;
  public:
    SubMatrixExpr (TA & aa, size_t fr, size_t fc, size_t ah, size_t aw) 
      : a(aa), first_row(fr), first_col(fc), height(ah), width(aw) { ; }

    INLINE size_t Height() const { return height; }
    INLINE size_t Width() const { return width; }

    // auto operator() (size_t i, size_t j) { return a(i+first_row, j+first_col); }
    // auto operator() (size_t i) { return a(i+first_row); }
    INLINE decltype(auto) operator() (size_t i, int j) const  { return a(i+first_row, j+first_col); }
    INLINE decltype(auto) operator() (size_t i) const { return a(i+first_row); }

    typedef typename TA::TELEM TELEM;
    typedef typename TA::TSCAL TSCAL;
    
    enum { IS_LINEAR = 0 };
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
    TA & a;
    size_t row;
  public:
    RowExpr (TA & aa, size_t r)
      : a(aa), row(r) { ; }

    INLINE size_t Height() const { return 1; }
    INLINE size_t Width() const { return a.Width(); }

    INLINE auto operator() (size_t i, size_t j) -> decltype(a(0,0)) { return a(row,i); }
    INLINE auto operator() (size_t i) -> decltype(a(0,0)) { return a(row,i); }
    INLINE auto operator() (size_t i, size_t j) const { return a(row,i); }
    INLINE auto operator() (size_t i) const { return a(row,i); }

    enum { IS_LINEAR = 0 };
    
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
    TA & a;
    size_t col;
  public:
    ColExpr (TA & aa, size_t c)
      : a(aa), col(c) { ; }

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return 1; }

    INLINE auto operator() (size_t i, size_t j) -> decltype(a(0,0)) { return a(i,col); }
    INLINE auto operator() (size_t i) -> decltype(a(0,0)) { return a(i,col); }
    INLINE auto operator() (size_t i, size_t j) const { return a(i,col); }
    INLINE auto operator() (size_t i) const { return a(i,col); }

    enum { IS_LINEAR = 0 };
    
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
    // typedef typename TA::TSCAL TSCAL;

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
    enum { IS_LINEAR = 0 };

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
    // typedef typename TA::TSCAL TSCAL;

    INLINE ColsArrayExpr (const TA & aa, FlatArray<int> acols) : a(aa), cols(acols) { ; }

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return cols.Size(); }

    INLINE auto operator() (size_t i, size_t j) const -> decltype(a(i, cols[j]))  { return a(i, cols[j]); }
    INLINE auto operator() (size_t i) const -> decltype(a(i, cols[0]))  { return a(i, cols[0]); }

    enum { IS_LINEAR = 0 };

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

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }
    
    INLINE auto View() const { return *this; }
    INLINE auto Shape() const { return a.Shape(); }
    
    INLINE auto operator() (size_t i, size_t j) const { return Conj(a(i,j)); }
    INLINE auto operator() (size_t i) const { return Conj(a(i)); }

    enum { IS_LINEAR = 0 };
  };


  /// Conjugate
  template <typename TA>
  INLINE auto Conj (const Expr<TA> & a)
  {
    return ConjExpr (a.View()); 
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

    INLINE size_t Height() const { return a.Height(); }
    INLINE size_t Width() const { return a.Width(); }
 
    INLINE auto operator() (size_t i, size_t j) const { return Truncate(a(i,j), eps); }
    INLINE auto operator() (size_t i) const { return Truncate(a(i), eps); }
    INLINE auto View() const { return *this; }
    enum { IS_LINEAR = TA::IS_LINEAR };
  };

  /// Conjugate
  template <typename TA>
  INLINE TruncateExpr<TA>
  Truncate (const Expr<TA> & a, double eps = 1e-12)
  {
    return TruncateExpr<TA> (a.View(), eps);
  }
  

  /* ************************* InnerProduct ********************** */

  /*
  INLINE double InnerProduct (double a, double b) {return a * b;}
  INLINE Complex InnerProduct (Complex a, Complex b) {return a * b;}
  INLINE Complex InnerProduct (double a, Complex b) {return a * b;}
  INLINE Complex InnerProduct (Complex a, double b) {return a * b;}
  INLINE int InnerProduct (int a, int b) {return a * b;}

  template <int DIM>
  AutoDiff<DIM> InnerProduct (AutoDiff<DIM> a, AutoDiff<DIM> b) {return a * b;}

  template <int N>
  auto InnerProduct (SIMD<double,N> a, SIMD<double,N> b) { return a*b; }
  template <int N>
  auto InnerProduct (SIMD<Complex,N> a, SIMD<double,N> b) { return a*b; }
  template <int N>
  auto InnerProduct (SIMD<double,N> a, SIMD<Complex,N> b) { return a*b; }
  template <int N>
  auto InnerProduct (SIMD<Complex,N> a, SIMD<Complex,N> b) { return a*b; }
  */
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
    if (a.Height()*a.Width() == 0) return 0; 

    auto sum = InnerProduct (a.Spec()(0), b.Spec()(0));
    for (size_t i = 1; i < a.Height()*a.Width(); i++)
      sum += InnerProduct (a.Spec()(i), b.Spec()(i));
    return sum;
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
    if (TA::IS_LINEAR)
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

    if (TA::IS_LINEAR)
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



  
  extern NGS_DLL_HEADER void CalcLU (SliceMatrix<double> A, FlatArray<int> p);
  extern NGS_DLL_HEADER void InverseFromLU (SliceMatrix<double> A, FlatArray<int> p);
  extern NGS_DLL_HEADER void SolveFromLU (SliceMatrix<double> A, FlatArray<int> p, SliceMatrix<double,ColMajor> X);
  extern NGS_DLL_HEADER void SolveTransFromLU (SliceMatrix<double> A, FlatArray<int> p, SliceMatrix<double,ColMajor> X);
  



  
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
