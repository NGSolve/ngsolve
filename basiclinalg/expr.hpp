#ifndef FILE_NGBLA_EXPR
#define FILE_NGBLA_EXPR

/**************************************************************************/
/* File:   expr.hpp                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jan. 02                                                    */
/**************************************************************************/



namespace ngbla
{

  template <int H, int W, typename T> class Mat;
  template <int S, typename T> class Vec;


  /*
    Matrix expression templates
  */


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






  template <typename TM> 
  inline typename TM::TELEM & Access (TM & mat, int i, int j)
  {
    return mat(i,j);
  }

  template <typename TM> 
  inline typename TM::TELEM Access (const TM & mat, int i, int j)
  {
    return mat(i,j);
  }

  inline double & Access (double & mat, int i, int j)
  {
    return mat;
  }

  inline double Access (const double & mat, int i, int j)
  {
    return mat;
  }

  inline Complex & Access (Complex & mat, int i, int j)
  {
    return mat;
  }

  inline Complex Access (const Complex & mat, int i, int j)
  {
    return mat;
  }










  /**
     Trait to obtain vector and scalar types for given matrix types.
     Is specified for double, Complex, AutoDiff<doube>, AutoDiff<Complex>
  */
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


  // some compiler thinks it needs Mat<int>...
  template <>
  class mat_traits<Mat<1,1,int> >
  {
  public:
    typedef int TELEM;
    typedef int TSCAL;
    typedef int TV_COL;
    typedef int TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
    enum { IS_COMPLEX = 0 };
  };

  template <>
  class mat_traits<int>
  {
  public:
    typedef int TELEM;
    typedef int TSCAL;
    typedef int TV_COL;
    typedef int TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
    enum { IS_COMPLEX = 0 };
  };


  template <>
  class mat_traits<double>
  {
  public:
    typedef double TELEM;
    typedef double TSCAL;
    typedef double TV_COL;
    typedef double TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
    enum { IS_COMPLEX = 0 };
  };

  template <>
  class mat_traits<Complex>
  {
  public:
    typedef Complex TELEM;
    typedef Complex TSCAL;
    typedef Complex TV_COL;
    typedef Complex TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
    enum { IS_COMPLEX = 1 };
  };

  template <int D, typename SCAL>
  class mat_traits<AutoDiff<D,SCAL> >
  {
  public:
    typedef AutoDiff<D,SCAL> TELEM;
    typedef AutoDiff<D,SCAL> TSCAL;
    typedef AutoDiff<D,SCAL> TV_COL;
    typedef AutoDiff<D,SCAL> TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
    enum { IS_COMPLEX = mat_traits<SCAL>::IS_COMPLEX };
  };

  template <int D>
  class mat_traits<AutoDiffDiff<D> >
  {
  public:
    typedef AutoDiffDiff<D> TELEM;
    typedef AutoDiffDiff<D> TSCAL;
    typedef AutoDiffDiff<D> TV_COL;
    typedef AutoDiffDiff<D> TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
    enum { IS_COMPLEX = false };
  };


  template <>
  class mat_traits<const int>
  {
  public:
    typedef int TELEM;
    typedef int TSCAL;
    typedef int TV_COL;
    typedef int TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
    enum { IS_COMPLEX = 0 };
  };

  template <>
  class mat_traits<const double>
  {
  public:
    typedef double TELEM;
    typedef double TSCAL;
    typedef double TV_COL;
    typedef double TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
    enum { IS_COMPLEX = 0 };
  };

  template <>
  class mat_traits<const Complex>
  {
  public:
    typedef Complex TELEM;
    typedef Complex TSCAL;
    typedef Complex TV_COL;
    typedef Complex TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
    enum { IS_COMPLEX = 1 };
  };


  /// matrix type from column and row vectors
  template <typename TV_COL, typename TV_ROW>
  class mat_from_vecs
  {
    enum { HEIGHT = mat_traits<TV_COL>::HEIGHT };
    enum { WIDTH  = mat_traits<TV_ROW>::HEIGHT };
    typedef typename mat_from_vecs<typename TV_COL::TELEM, typename TV_ROW::TELEM>::TMAT TELEM;
    typedef Mat<HEIGHT,WIDTH,TELEM> TMAT;
  };

  template <> class mat_from_vecs<double,double> { typedef double TMAT; };
  template <> class mat_from_vecs<double,Complex> { typedef Complex TMAT; };
  template <> class mat_from_vecs<Complex,double> { typedef Complex TMAT; };
  template <> class mat_from_vecs<Complex,Complex> { typedef Complex TMAT; };



  /// matrix type of product
  template <typename TA, typename TB>
  class mat_prod_type
  {
  public:
    enum { HEIGHT = mat_traits<TA>::HEIGHT };
    enum { WIDTH  = mat_traits<TB>::WIDTH };
    typedef typename mat_prod_type<typename TA::TELEM, 
                                   typename TB::TELEM>::TMAT TELEM;
    typedef Mat<HEIGHT,WIDTH,TELEM> TMAT;
  };

  template <int S, typename T> class mat_prod_type<double, Vec<S,T> > { public: typedef Vec<S,T> TMAT; };
  template <int S, typename T> class mat_prod_type<Complex, Vec<S,T> > { public: typedef Vec<S,T> TMAT; };
  template <> class mat_prod_type<double,double> { public: typedef double TMAT; };
  template <> class mat_prod_type<double,Complex> { public: typedef Complex TMAT; };
  template <> class mat_prod_type<Complex,double> { public: typedef Complex TMAT; };
  template <> class mat_prod_type<Complex,Complex> { public: typedef Complex TMAT; };

  template <int D, typename TA, typename TB>
  class mat_prod_type<AutoDiff<D, TA>, TB> 
  { public: typedef AutoDiff<D, typename mat_prod_type<TA,TB>::TMAT> TMAT; };

  template <int D, typename TA, typename TB>
  class mat_prod_type<TB, AutoDiff<D, TA> > 
  { public: typedef AutoDiff<D, typename mat_prod_type<TA,TB>::TMAT> TMAT; };

  template <int D, typename TA, int E, typename TB>
  class mat_prod_type<AutoDiff<D, TA>, AutoDiff<E,TB> > 
  { public: typedef AutoDiff<D, typename mat_prod_type<TA,TB>::TMAT> TMAT; };

  /// matrix type of sum (important for double+Complex)
  template <typename TA, typename TB>
  class mat_sum_type
  {
  public:
    enum { HEIGHT = mat_traits<TA>::HEIGHT };
    enum { WIDTH  = mat_traits<TA>::WIDTH };
    typedef typename mat_sum_type<typename TA::TELEM, 
                                  typename TB::TELEM>::TMAT TELEM;
    typedef Mat<HEIGHT,WIDTH,TELEM> TMAT;
  };

  template <> class mat_sum_type<double,double> { public: typedef double TMAT; };
  template <> class mat_sum_type<double,Complex> { public: typedef Complex TMAT; };
  template <> class mat_sum_type<Complex,double> { public: typedef Complex TMAT; };
  template <> class mat_sum_type<Complex,Complex> { public: typedef Complex TMAT; };

  template <int D, typename TA, typename TB>
  class mat_sum_type<AutoDiff<D, TA>, TB> 
  { public: typedef AutoDiff<D, typename mat_sum_type<TA,TB>::TMAT> TMAT; };


  template <int D, typename TA, typename TB>
  class mat_sum_type<TB, AutoDiff<D, TA> > 
  { public: typedef AutoDiff<D, typename mat_sum_type<TA,TB>::TMAT> TMAT; };


  template <int D, typename TA, int E, typename TB>
  class mat_sum_type<AutoDiff<D, TA>, AutoDiff<E,TB> > 
  { public: typedef AutoDiff<D, typename mat_sum_type<TA,TB>::TMAT> TMAT; };







  /// matrix type of scale
  template <typename TM, typename TS>
  class mat_scale_type
  {
  public:
    enum { HEIGHT = mat_traits<TM>::HEIGHT };
    enum { WIDTH  = mat_traits<TM>::WIDTH };
    typedef typename mat_scale_type<typename TM::TELEM, TS>::TMAT TELEM;
    typedef Mat<HEIGHT,WIDTH,TELEM> TMAT;
  };

  template <> class mat_scale_type<double,double> { public: typedef double TMAT; };
  template <> class mat_scale_type<double,Complex> { public: typedef Complex TMAT; };
  template <> class mat_scale_type<Complex,double> { public: typedef Complex TMAT; };
  template <> class mat_scale_type<Complex,Complex> { public: typedef Complex TMAT; };

  template <int D, typename TA, typename TB>
  class mat_scale_type<TB, AutoDiff<D, TA> > 
  { public: typedef AutoDiff<D, TA> TMAT; };

  template <int D, typename TA, typename TB>
  class mat_scale_type<AutoDiff<D, TA>, TB> 
  { public: typedef AutoDiff<D, TA> TMAT; };

  template <int D, typename TA, int E, typename TB>
  class mat_scale_type<AutoDiff<D, TA>, AutoDiff<E,TB> > 
  { public: typedef AutoDiff<D, typename mat_scale_type<TA,TB>::TMAT> TMAT; };




  /// Height of matrix
  template <class TM> 
  inline int Height (const TM & m)
  {
    return m.Height();
  }

  /// Width of matrix
  template <class TM> 
  inline int Width (const TM & m)
  {
    return m.Width();
  }

  template <> inline int Height<double> (const double&) { return 1; }
  template <> inline int Height<Complex> (const Complex&) { return 1; }
  template <> inline int Width<double> (const double&) { return 1; }
  template <> inline int Width<Complex> (const Complex&) { return 1; }


  /// Complex to double assignment called
  class Complex2RealException : public Exception
  {
  public:
    Complex2RealException()
      : Exception("Assignment of Complex 2 Real") { ; }
  };

  /*
  /// converts Complex to type T. Throw exception if illegal
  template<class T>
  inline T ReduceComplex (Complex x)
  {
    return x;
  }

  template<> 
  inline double ReduceComplex<double> (Complex x)
  {
    throw Complex2RealException();
    }
  */

  /*
  /// converts Complex to type T.
  template<class T>
  inline T ReduceComplex2 (Complex x)
  {
    return x;
  }

  template<> 
  inline double ReduceComplex2<double> (Complex x)
  {
    return x.real()+x.imag();//abs(x);
  }
  */

  /*
  /// converts type T to double. Throw exception if illegal
  template<class T>
  inline double ReduceToReal (T x)
  {
    return x;
  }

  template<>
  inline double ReduceToReal (Complex x)
  {
    throw Complex2RealException();
  }
  */




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


  template <typename TO>
  inline TO ConvertTo (const AutoDiff<1, Complex> & f)
  {
    return TO(f);
  }


  template <>
  inline double ConvertTo (Complex f)
  {
    throw Complex2RealException();
  }











  



#ifdef CHECK_RANGE
  /// matrices do not fit for matrix-matrix operation
  class MatrixNotFittingException : public Exception
  {
  public:
    MatrixNotFittingException (const string & where, int h1, int w1,
                               int h2, int w2)
      : Exception ("") 
    {
      stringstream str;
      str << where << " mat1 = " << h1 << " x " << w1 << ", mat2 = " << h2 << " x " << w2 << endl;
      Append (str.str());
    }
  };
#endif


  /**
     Expr is the base class for all matrix template expressions.
     Barton and Nackman Trick for template polymorphism, function Spec.

     provides Height and Width of matrix.
     IsLinear allows linear matrix element access.
  */

  template <typename T>
  class Expr 
  {
  public:
    Expr () { ; }

    typedef T TConv;  
    /// cast to specific type
    T & Spec() { return static_cast<T&> (*this); }
    /// cast to specific type
    const T & Spec() const { return static_cast<const T&> (*this); }

    /// height 
    int Height() const { return Spec().T::Height(); }
    int Width() const { return Spec().T::Width(); }

	void Dump (ostream & ost) const { Spec().T::Dump(ost); }
  };



  /**
     Caller knows that matrix expression is a symmetric matrix.
     Thus, only one half of the matrix needs to be computed.
  */
  template <typename T>
  class SymExpr : public Expr<SymExpr<T> >
  {
    const T a;
  public:
    typedef typename T::TELEM TELEM;
    typedef typename mat_traits<TELEM>::TSCAL TSCAL;

    SymExpr (const T & aa) : a(aa) { ; }

    TELEM operator() (int i) const { return a(i); }
    TELEM operator() (int i, int j) const { return a(i,j); }
    int Height() const { return a.Height(); }
    int Width() const { return a.Width(); }
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
    return SymExpr<T> (static_cast<const T&> (a));
  }



  /**
     Reference to matrix expression.
     T is matrix class.
  */
  template <typename T>
  class RefMatExpr : public Expr<RefMatExpr<T> >
  {
  protected:
    const T & data;
  public:
    typedef typename T::TELEM TELEM;
    typedef typename T::TSCAL TSCAL;

    RefMatExpr (const T & d) : data(d) { ; }

    TELEM operator() (int i) const  { return data(i); }
    TELEM operator() (int i, int j) const  { return data(i, j); }

    int Height() const { return data.T::Height(); }
    int Width() const { return data.T::Width(); }

    enum { IS_LINEAR = T::IS_LINEAR };

    void Dump (ostream & ost) const { ost << "Matrix-Ref(" << data << ")"; }
  };



  /**
     The base class for matrices.
  */
  template <class T>
  class MatExpr : public Expr<T>
  {
  public:
    typedef RefMatExpr<T> TConv;


//	typedef typename mat_traits<T>::TELEM TELEM;
//    typedef typename mat_traits<TELEM>::TSCAL TSCAL;

//	typedef typename T::TSCAL TSCAL;

    int Height() const { return Spec().T::Height(); }
    int Width() const { return Spec().T::Width(); }

    MatExpr () : Expr<T>() { ; }

    T & Spec() { return static_cast<T&> (*this); }
    const T & Spec() const { return static_cast<const T&> (*this); }

    enum { IS_LINEAR = 1 };

	void Dump (ostream & ost) const { ost << "Matrix"; }
	
    template<typename TB>
    T & operator= (const Expr<TB> & v)
    {
#ifdef CHECK_RANGE
      if (Height() != v.Height() || Width() != v.Width())
        {
          throw MatrixNotFittingException ("operator=", 
                                           Height(), Width(),
                                           v.Height(), v.Width());
        }
#endif

      if (TB::IS_LINEAR)
        {
          int hw = Expr<T>::Height() * Expr<T>::Width();
          for (int i = 0; i < hw; i++)
            Spec()(i) = v.Spec()(i);
        }
      else
        {
          int h = Expr<T>::Height();
          int w = Expr<T>::Width();

          if (T::IS_LINEAR)
            for (int i = 0, k = 0; i < h; i++)
              for (int j = 0; j < w; j++, k++)
                Spec()(k) = v.Spec()(i,j);
          else
            for (int i = 0; i < h; i++)
              for (int j = 0; j < w; j++)
                Spec()(i,j) = v.Spec()(i,j);
        }
      return Spec();
    }


    template<typename TB>
    T & operator+= (const Expr<TB> & v)
    {
#ifdef CHECK_RANGE
      if (Height() != v.Height() || Width() != v.Width())
        throw MatrixNotFittingException ("operator+=", 
                                         Height(), Width(),
                                         v.Height(), v.Width());
#endif

      if (TB::IS_LINEAR)
        {
          if (T::IS_LINEAR)
            {
              int hw = Height() * Width(); 
              for (int i = 0; i < hw; i++)
                Spec()(i) += v.Spec()(i);
            }
          else
            {
              int h = Height();
              int w = Width();
              for (int i = 0, k = 0; i < h; i++)
                for (int j = 0; j < w; j++, k++)
                  Spec()(i,j) += v.Spec()(k);
            }
        }
      else
        {
          int h = Height();
          int w = Width();

          if (T::IS_LINEAR)
            for (int i = 0, k = 0; i < h; i++)
              for (int j = 0; j < w; j++, k++)
                Spec()(k) += v.Spec()(i,j);
          else
            for (int i = 0; i < h; i++)
              for (int j = 0; j < w; j++)
                Spec()(i,j) += v.Spec()(i,j);
        }
      return Spec();
    }




    template<typename TB>
    T & operator+= (const Expr<SymExpr<TB> > & v)
    {
#ifdef CHECK_RANGE
      if (Height() != v.Height() || Width() != v.Width())
        throw MatrixNotFittingException ("operator+=", 
                                         Height(), Width(),
                                         v.Height(), v.Width());
#endif
      int h = Height();
      for (int i = 0; i < h; i++)
        {
          for (int j = 0; j < i; j++)
            {
              double val = v.Spec()(i,j);
              Spec()(i,j) += val;
              Spec()(j,i) += val;
            }
          Spec()(i,i) += v.Spec()(i,i);
        }
      return Spec();
    }



    MatExpr<T> & operator+= (double scal)
    {
      int hw = Height() * Width();
      for (int i = 0; i < hw; i++)
        Spec()(i) += scal;
	  return *this;
    }
  

    template<typename TB>
    MatExpr<T> & operator-= (const Expr<TB> & v)
    {
#ifdef CHECK_RANGE
      if (Height() != v.Height() || Width() != v.Width())
        throw MatrixNotFittingException ("operator-=", 
                                         Height(), Width(),
                                         v.Height(), v.Width());
#endif
      if (TB::IS_LINEAR)
        {
          if (T::IS_LINEAR)
            {
              int hw = Height() * Width(); 
              for (int i = 0; i < hw; i++)
                Spec()(i) -= v.Spec()(i);
            }
          else
            {
              int h = Height();
              int w = Width();
              for (int i = 0, k = 0; i < h; i++)
                for (int j = 0; j < w; j++, k++)
                  Spec()(i,j) -= v.Spec()(k);
            }
        }
      else
        {
          int h = Height();
          int w = Width();

          if (T::IS_LINEAR)
            for (int i = 0, k = 0; i < h; i++)
              for (int j = 0; j < w; j++, k++)
                Spec()(k) -= v.Spec()(i,j);
          else
            for (int i = 0; i < h; i++)
              for (int j = 0; j < w; j++)
                Spec()(i,j) -= v.Spec()(i,j);
        }
      return Spec();
    }


    template <class SCAL2>
    T & operator*= (const SCAL2 & s)
    {
      int hw = Height() * Width();
      for (int i = 0; i < hw; i++)
        Spec()(i) *= s;
      return Spec();
    }

    template <class SCAL2>
    T & operator/= (const SCAL2 & s)
    {
      return (*this) *= (1./s);
    }
  };









  /**
     The base class for matrices.
     Constant-Means-Constat-Pointer
     matrix-values may be chaneged by const methods
  */
  template <class T>
  class CMCPMatExpr : public MatExpr<T>
  {
  public:
    typedef RefMatExpr<T> TConv;

    int Height() const { return Spec().T::Height(); }
    int Width() const { return Spec().T::Width(); }

    CMCPMatExpr () { ; }

    T & Spec() { return static_cast<T&> (*this); }
    const T & Spec() const { return static_cast<const T&> (*this); }

    // enum { IS_LINEAR = 1 };

    template<typename TB>
    const T & operator= (const Expr<TB> & v) const
    {
#ifdef CHECK_RANGE
      if (Height() != v.Height() || Width() != v.Width())
        {
          throw MatrixNotFittingException ("operator=", 
                                           Height(), Width(),
                                           v.Height(), v.Width());
        }
#endif

      if (TB::IS_LINEAR)
        {
          if (T::IS_LINEAR)
            {
              int hw = Expr<T>::Height() * Expr<T>::Width();
              for (int i = 0; i < hw; i++)
                Spec()(i) = v.Spec()(i);
            }
          else
            {
              int h = Expr<T>::Height();
              int w = Expr<T>::Width();
              for (int i = 0, k = 0; i < h; i++)
                for (int j = 0; j < w; j++, k++)
                  Spec()(i,j) = v.Spec()(k);
            }
        }
      else
        {
          int h = Expr<T>::Height();
          int w = Expr<T>::Width();

          if (T::IS_LINEAR)
            for (int i = 0, k = 0; i < h; i++)
              for (int j = 0; j < w; j++, k++)
                Spec()(k) = v.Spec()(i,j);
          else
            for (int i = 0; i < h; i++)
              for (int j = 0; j < w; j++)
                Spec()(i,j) = v.Spec()(i,j);
        }
      return Spec();
    }


    template<typename TB>
    const T & operator+= (const Expr<TB> & v) const
    {
#ifdef CHECK_RANGE
      if (Height() != v.Height() || Width() != v.Width())
        throw MatrixNotFittingException ("operator+=", 
                                         Height(), Width(),
                                         v.Height(), v.Width());
#endif

      if (TB::IS_LINEAR)
        {
          if (T::IS_LINEAR)
            {
              int hw = Height() * Width(); 
              for (int i = 0; i < hw; i++)
                Spec()(i) += v.Spec()(i);
            }
          else
            {
              int h = Height();
              int w = Width();
              for (int i = 0, k = 0; i < h; i++)
                for (int j = 0; j < w; j++, k++)
                  Spec()(i,j) += v.Spec()(k);
            }
        }
      else
        {
          int h = Height();
          int w = Width();

          if (T::IS_LINEAR)
            for (int i = 0, k = 0; i < h; i++)
              for (int j = 0; j < w; j++, k++)
                Spec()(k) += v.Spec()(i,j);
          else
            {
              for (int i = 0; i < h; i++)
                for (int j = 0; j < w; j++)
                  Spec()(i,j) += v.Spec()(i,j);
            }
        }
      return Spec();
    }


    template<typename TB>
    const T & operator+= (const Expr<SymExpr<TB> > & v) const
    {
#ifdef CHECK_RANGE
      if (Height() != v.Height() || Width() != v.Width())
        throw MatrixNotFittingException ("operator+=", 
                                         Height(), Width(),
                                         v.Height(), v.Width());
#endif
      typedef typename mat_traits<typename TB::TELEM>::TSCAL TSCAL;

      int h = Height();
      for (int i = 0; i < h; i++)
        {
          for (int j = 0; j < i; j++)
            {
              TSCAL val = v.Spec()(i,j);
              Spec()(i,j) += val;
              Spec()(j,i) += val;
            }
          Spec()(i,i) += v.Spec()(i,i);
        }
      return Spec();
    }



    const T & operator+= (double scal) const
    {
      int hw = Height() * Width();
      for (int i = 0; i < hw; i++)
        Spec()(i) += scal;
      return Spec();
    }
  

    template<typename TB>
    const T & operator-= (const Expr<TB> & v) const
    {
#ifdef CHECK_RANGE
      if (Height() != v.Height() || Width() != v.Width())
        throw MatrixNotFittingException ("operator-=", 
                                         Height(), Width(),
                                         v.Height(), v.Width());
#endif

      if (TB::IS_LINEAR)
        {
          if (T::IS_LINEAR)
            {
              int hw = Height() * Width(); 
              for (int i = 0; i < hw; i++)
                Spec()(i) -= v.Spec()(i);
            }
          else
            {
              int h = Height();
              int w = Width();
              for (int i = 0, k = 0; i < h; i++)
                for (int j = 0; j < w; j++, k++)
                  Spec()(i,j) -= v.Spec()(k);
            }
        }
      else
        {
          int h = Height();
          int w = Width();

          if (T::IS_LINEAR)
            for (int i = 0, k = 0; i < h; i++)
              for (int j = 0; j < w; j++, k++)
                Spec()(k) -= v.Spec()(i,j);
          else
            {
              for (int i = 0; i < h; i++)
                for (int j = 0; j < w; j++)
                  Spec()(i,j) -= v.Spec()(i,j);
            }
        }
      return Spec();
    }


    template <class SCAL2>
    const T & operator*= (const SCAL2 & s) const
    {
      int hw = Height() * Width();
      for (int i = 0; i < hw; i++)
        Spec()(i) *= s;
      return Spec();
    }

    template <class SCAL2>
    const T & operator/= (const SCAL2 & s) const 
    {
      return (*this) *= (1.0/s);
    }
  };

















  /* *************************** SumExpr **************************** */

  /**
     Sum of 2 matrix expressions
  */


  template <class TA, class TB> 
  class SumExpr : public Expr<SumExpr<TA,TB> >
  {
    const TA a;
    const TB b;
  public:
    typedef typename mat_sum_type<typename TA::TELEM,
                                  typename TB::TELEM>::TMAT TELEM;
    typedef typename mat_traits<TELEM>::TSCAL TSCAL;

    SumExpr (const TA & aa, const TB & ab) : a(aa), b(ab) 
    { ; }

    TELEM operator() (int i) const { return a(i)+b(i); }
    TELEM operator() (int i, int j) const { return a(i,j)+b(i,j); }
    int Height() const { return a.Height(); }
    int Width() const { return a.Width(); }

    enum { IS_LINEAR = TA::IS_LINEAR && TB::IS_LINEAR };
    void Dump (ostream & ost) const
    { ost << "("; a.Dump(ost); ost << ") + ("; b.Dump(ost); ost << ")"; }
  };

  template <typename TA, typename TB>
  inline SumExpr<typename TA::TConv, typename TB::TConv>
  operator+ (const Expr<TA> & a, const Expr<TB> & b)
  {
    return SumExpr<typename TA::TConv, typename TB::TConv> (static_cast <const TA&> (a), static_cast <const TB&> (b));
  }



  /* *************************** SubExpr **************************** */


  /**
     Matrix-expr minus Matrix-expr
  */

  template <class TA, class TB> 
  class SubExpr : public Expr<SubExpr<TA,TB> >
  {
    const TA a;
    const TB b;
  public:
    // typedef typename TA::TELEM TELEM;
    // typedef typename TA::TSCAL TSCAL;

    typedef typename mat_sum_type<typename TA::TELEM,
                                  typename TB::TELEM>::TMAT TELEM;
    typedef typename mat_traits<TELEM>::TSCAL TSCAL;



    SubExpr (const TA & aa, const TB & ab) : a(aa), b(ab) { ; }

    TELEM operator() (int i) const { return a(i)-b(i); }
    TELEM operator() (int i, int j) const { return a(i,j)-b(i,j); }
    int Height() const { return a.Height(); }
    int Width() const { return a.Width(); }
    // static bool IsLinear() { return TA::IsLinear() && TB::IsLinear(); }
    enum { IS_LINEAR = TA::IS_LINEAR && TB::IS_LINEAR };
  };


  template <typename TA, typename TB>
  inline SubExpr<typename TA::TConv, typename TB::TConv>
  operator- (const Expr<TA> & a, const Expr<TB> & b)
  {
    return SubExpr<typename TA::TConv, typename TB::TConv> 
      (static_cast <const TA&> (a), static_cast <const TB&> (b));
  }







  /* *************************** MinusExpr **************************** */


  /**
     minus Matrix-expr
  */

  template <class TA>
  class MinusExpr : public Expr<MinusExpr<TA> >
  {
    const TA a;
  public:
    typedef typename TA::TELEM TELEM;
    typedef typename TA::TSCAL TSCAL;

    MinusExpr (const TA & aa) : a(aa) { ; }

    TELEM operator() (int i) const { return -a(i); }
    TELEM operator() (int i, int j) const { return -a(i,j); }
    int Height() const { return a.Height(); }
    int Width() const { return a.Width(); }

    enum { IS_LINEAR = TA::IS_LINEAR };
  };



  template <typename TA>
  inline MinusExpr<typename TA::TConv>
  operator- (const Expr<TA> & a)
  {
    return MinusExpr<typename TA::TConv> (static_cast<const TA&> (a) );
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
    typedef typename mat_scale_type<typename TA::TELEM, TS>::TMAT TELEM;
    typedef typename mat_traits<TELEM>::TSCAL TSCAL;
    enum { IS_LINEAR = TA::IS_LINEAR };

    ScaleExpr (const TA & aa, TS as) : a(aa), s(as) { ; }

    TELEM operator() (int i) const { return s * a(i); }
    TELEM operator() (int i, int j) const { return s * a(i,j); }

    int Height() const { return a.Height(); }
    int Width() const { return a.Width(); }
    //  static bool IsLinear() { return TA::IsLinear(); }
    void Dump (ostream & ost) const
    { ost << "Scale, s=" << s << " * "; a.Dump(ost);  }
  };

  template <typename TA>
  inline ScaleExpr<typename TA::TConv, double> 
  operator* (double b, const Expr<TA> & a)
  {
    return ScaleExpr<typename TA::TConv, double> (static_cast<const TA&> (a), b);
  }

  template <typename TA>
  inline ScaleExpr<TA, Complex> 
  operator* (Complex b, const Expr<TA> & a)
  {
    return ScaleExpr<TA, Complex> (static_cast<const TA&> (a), b);
  }

  template <int D, typename TAD, typename TA>
  inline ScaleExpr<typename TA::TConv, AutoDiff<D,TAD> > 
  operator* (const AutoDiff<D,TAD> & b, const Expr<TA> & a)
  {
    return ScaleExpr<typename TA::TConv, AutoDiff<D,TAD> > (static_cast<const TA&> (a), b );
  }




  /* ************************* MultExpr ************************* */


  /**
     Matrix-expr timex Matrix-expr
  */
  template <class TA, class TB> class MultExpr : public Expr<MultExpr<TA,TB> >
  {
    const TA a;
    const TB b;
  public:
    typedef typename mat_prod_type<typename TA::TELEM, 
                                   typename TB::TELEM>::TMAT TELEM;
    typedef typename mat_traits<TELEM>::TSCAL TSCAL;

    MultExpr (const TA & aa, const TB & ab) : a(aa), b(ab) { ; }

    TELEM operator() (int i) const { return operator()(i,0); }  //   TELEM(TSCAL(0)); }  JS, 310508
    TELEM operator() (int i, int j) const
    { 
      int wa = a.Width();
      TELEM sum;
      if (wa >= 1)
        {
          sum = a(i,0) * b(0,j);
          for (int k = 1; k < wa; k++)
            sum += a(i,k) * b(k, j);
        }
      else
        sum = TSCAL(0);
    
      return sum;
    }

    int Height() const { return a.Height(); }
    int Width() const { return b.Width(); }
    enum { IS_LINEAR = 0 };
  };


  template <typename TA, typename TB>
  inline MultExpr<typename TA::TConv, typename TB::TConv>
  operator* (const Expr<TA> & a, const Expr<TB> & b)
  {
    return MultExpr<typename TA::TConv, typename TB::TConv> 
      (static_cast <const TA&> (a), static_cast <const TB&> (b));
  }


  /* ************************** Trans *************************** */


  /**
     Transpose of Matrix-expr
  */
  template <class TA> class TransExpr : public Expr<TransExpr<TA> >
  {
    const TA a;
  public:
    typedef typename TA::TELEM TELEM;

    TransExpr (const TA & aa) : a(aa) { ; }

    int Height() const { return a.Width(); }
    int Width() const { return a.Height(); }

    TELEM operator() (int i, int j) const { return Trans (a(j,i)); }
    TELEM operator() (int i) const { return 0; }

    enum { IS_LINEAR = 0 };
  };


  /// Transpose 
  template <typename TA>
  inline TransExpr<typename TA::TConv>
  Trans (const Expr<TA> & a)
  {
    return TransExpr<typename TA::TConv> (static_cast <const TA&> (a));
  }

  inline double Trans (double a) { return a; }
  inline Complex Trans (Complex a) { return a; }
  template<int D, typename TAD>
  inline AutoDiff<D,TAD> Trans (const AutoDiff<D,TAD> & a) { return a; }


  /* ************************* Conjugate *********************** */ 


  inline double Conj (double a)
  {
    return a;
  }

  inline Complex Conj (Complex a)
  {
    return conj(a);
  }

  /**
     Conjugate of Matrix-expr
  */
  // Attention NOT transpose !!! elemwise conjugate !!! 
  template <class TA> class ConjExpr : public Expr<ConjExpr<TA> >
  {
    const TA a;
  public:
    typedef typename TA::TELEM TELEM;

    ConjExpr (const TA & aa) : a(aa) { ; }

    int Height() const { return a.Height(); }
    int Width() const { return a.Width(); }
 
    TELEM operator() (int i, int j) const { return Conj(a(i,j)); }
    TELEM operator() (int i) const { return Conj(a(i)); }

    enum { IS_LINEAR = 0 };
  };


  /// Conjugate
  template <typename TA>
  inline ConjExpr<typename TA::TConv>
  Conj (const Expr<TA> & a)
  {
    return ConjExpr<typename TA::TConv> (static_cast <const TA&> (a));
  }

  template<int D, typename TAD>
  inline AutoDiff<D,TAD> Conj (const AutoDiff<D,TAD> & a) 
  { 
    AutoDiff<D,TAD> b; 
    b.Value() = conj(a.Value()); 
  
    for(int i=0;i<D;i++) 
      b.DValue(i) = conj(a.DValue(i)); 

    return b;
  }



  /* ************************* ChangeSize ********************** */


  template <class TA>
  class ChangeSizeExpr : public Expr<ChangeSizeExpr<TA> >
  {
    const TA a;
    const int height;
    const int width;
  public:
    typedef typename TA::TELEM TELEM;
    typedef typename TA::TSCAL TSCAL;

    ChangeSizeExpr (const TA & aa, const int h, const int w) : a(aa),height(h),width(w) { ; }

    TELEM operator() (int i) const { return a(i); }
    TELEM operator() (int i, int j) const 
    { 
      //if(i >= a.Height() || j>= a.Width())
      //  return 0;
      return a(i,j); 
    }
    int Height() const { return height; }
    int Width() const { return width; }

    enum { IS_LINEAR = TA::IS_LINEAR };
  };

  template <typename TA>
  inline ChangeSizeExpr<typename TA::TConv>
  ChangeSize (const Expr<TA> & a, const int h, const int w)
  {
    return ChangeSizeExpr<typename TA::TConv> (static_cast<const TA&> (a),h,w );
  }


  /* ************************* InnerProduct ********************** */


  inline double InnerProduct (double a, double b)
  {
    return a * b;
  }

  inline Complex InnerProduct (Complex a, Complex b)
  {
    return a * b;
  }


  /**
     Inner product
  */
  template <class TA, class TB>
  // inline typename TA::TSCAL InnerProduct (const MatExpr<TA> & a, const MatExpr<TB> & b)
  inline typename TA::TSCAL InnerProduct (const Expr<TA> & a, const Expr<TB> & b)
  {
    typedef typename TA::TSCAL TSCAL;

    if (a.Height() == 0) return TSCAL(0);

    TSCAL sum = InnerProduct (a.Spec()(0), b.Spec()(0));
    for (int i = 1; i < a.Height(); i++)
      sum += InnerProduct (a.Spec()(i), b.Spec()(i));

    return sum;
  }

  /*
  // shouldn't be necessary with InnerProduct (Expr<>,Expr<>)
  template <class TA, class TB>
  inline typename TA::TSCAL InnerProduct (const MatExpr<TA> & a, const ConjExpr<TB> & b)
  {
  typedef typename TA::TSCAL TSCAL;
  TSCAL sum = 0;
  for (int i = 0; i < a.Height(); i++)
  sum += InnerProduct (a.Spec()(i), b.Spec()(i));
   
  return sum;
  }

  template <class TA, class TB>
  inline typename TA::TSCAL InnerProduct (const ConjExpr<TA> & a, const MatExpr<TB> & b)
  {
  typedef typename TA::TSCAL TSCAL;
  TSCAL sum = 0;

  for (int i = 0; i < a.Height(); i++)
  sum += InnerProduct (a.Spec()(i), b.Spec()(i));
   
  return sum;
  }


  template <class TA, class TB>
  inline typename TA::TSCAL InnerProduct (const ConjExpr<TA> & a, const ConjExpr<TB> & b)
  {
  typedef typename TA::TSCAL TSCAL;
  TSCAL sum = 0;
  for (int i = 0; i < a.Height(); i++)
  sum += InnerProduct (a.Spec()(i), b.Spec()(i));
     
  return sum;
  }
  */





  /* **************************** Trace **************************** */


  /**
     Calculates the trace of a matrix expression.
  */
  template <class TA>
  inline typename TA::TELEM Trace (const Expr<TA> & a)
  {
    typedef typename TA::TELEM TELEM;
    TELEM sum = 0;
    for (int i = 0; i < Height(a); i++)
      sum += a.Spec()(i,i);
    return sum;
  }

  /* **************************** L2Norm **************************** */

  /// Euklidean norm squared
  inline double L2Norm2 (double v)
  {
    return v*v;
  }

  inline double L2Norm2 (Complex v)
  {
    return v.real()*v.real()+v.imag()*v.imag();
  }

  template<int D, typename SCAL>
  inline double L2Norm2 (const AutoDiff<D,SCAL> & x) throw()
  {
    return L2Norm2(x.Value());
  }

  template <class TA>
  inline double L2Norm2 (const Expr<TA> & v)
  {
    typedef typename TA::TSCAL TSCAL;
    double sum = 0;
    if (TA::IS_LINEAR)
      for (int i = 0; i < v.Height()*v.Width(); i++)
        sum += L2Norm2 (v.Spec()(i));  
    else
      for (int i = 0; i < v.Height(); i++)
        for (int j = 0; j < v.Width(); j++)
          sum += L2Norm2 (v.Spec()(i,j));  
    
    return sum;
  }

  /**
     Calculates the Euklidean norm
  */
  template <class TA>
  inline double L2Norm (const Expr<TA> & v)
  {
    return sqrt (L2Norm2(v));
  }

  /* *************************** Output ****************************** */



  /// Print matrix-expr
  template<typename T>
  std::ostream & operator<< (std::ostream & s, const Expr<T> & v)
  {
    for (int i = 0; i < v.Height(); i++)
      {
        for (int j = 0 ; j < v.Width(); j++)
          s << " " << setw(7) << v.Spec()(i,j);
        s << endl;
      }
    return s;
  }



  /* **************************** Inverse *************************** */

//  template <typename T> class FlatMatrix;
  template <typename T> class Matrix;
  template <int H, int W, typename T> class Mat;



  inline void CalcInverse (const double & m, double & inv)
  {
    inv = 1 / m;
  }

  inline void CalcInverse (const Complex & m, Complex & inv)
  {
    inv = 1.0 / m;
  }

  template <int H, int W, typename T, typename TINV>
  inline void CalcInverse (const Mat<H,W,T> & m, TINV & inv)
  {
    switch (H)
      {
      case 1:
        {
          inv(0,0) = 1.0 / m(0,0);
          return;
        }
      case 2:
        {
          T idet = 1.0 / (m(0,0) * m(1,1) - m(0,1) * m(1,0));
          inv(0,0) = idet * m(1,1);
          inv(0,1) = -idet * m(0,1);
          inv(1,0) = -idet * m(1,0);
          inv(1,1) = idet * m(0,0);
          return;
        }
      case 3:
        {
          T h0 = m(4)*m(8)-m(5)*m(7);
          T h1 = m(5)*m(6)-m(3)*m(8);
          T h2 = m(3)*m(7)-m(4)*m(6);
          T det = m(0) * h0 + m(1) * h1 + m(2) * h2;
          T idet = 1.0 / det;

          inv(0,0) =  idet * h0; 
          inv(0,1) = -idet * (m(1) * m(8) - m(2) * m(7));
          inv(0,2) =  idet * (m(1) * m(5) - m(2) * m(4));
	
          inv(1,0) =  idet * h1; 
          inv(1,1) =  idet * (m(0) * m(8) - m(2) * m(6));
          inv(1,2) = -idet * (m(0) * m(5) - m(2) * m(3));
	
          inv(2,0) =  idet * h2; 
          inv(2,1) = -idet * (m(0) * m(7) - m(1) * m(6));
          inv(2,2) =  idet * (m(0) * m(4) - m(1) * m(3));
          return;
        }
      default:
        {
          FlatMatrix<T> fm(m);
          FlatMatrix<T> finv(inv);
          CalcInverse (fm, finv);
        }
      }
  }


  template <typename T, typename TINV>
  inline void CalcInverse (const Mat<2,2,T> & m, Mat<2,2,TINV> & inv)
  {
    T idet = 1.0 / (m(0,0) * m(1,1) - m(0,1) * m(1,0));
    inv(0,0) = idet * m(1,1);
    inv(0,1) = -idet * m(0,1);
    inv(1,0) = -idet * m(1,0);
    inv(1,1) = idet * m(0,0);
  }




  template <int H, int W, typename T>
  inline Mat<H,W,T> Inv (const Mat<H,W,T> & m)
  {
    Mat<H,W,T> inv;
    CalcInverse (m, inv);
    return inv;
  }



  /* ********************** Determinant *********************** */


  /**
     Calculates the determinant of a Matrix.
  */
  template <class T>
  inline typename T::TSCAL Det (const MatExpr<T> & m)
  {
    const T & sm = m.Spec();
    switch (sm.Height())
      {
      case 1: 
        {
          return sm(0,0); 
          break;
        }
      case 2:
        {
          return ( sm(0,0)*sm(1,1) - sm(0,1)*sm(1,0) ); 
          break;
        }
      case 3:
        {
          return 
            sm(0) * (sm(4) * sm(8) - sm(5) * sm(7)) +
            sm(1) * (sm(5) * sm(6) - sm(3) * sm(8)) +
            sm(2) * (sm(3) * sm(7) - sm(4) * sm(6));
          break;
        }
      default:
        {
          cerr << "general det not implemented" << endl;
        }
      }
    return typename T::TSCAL(0);
  }
}

#endif
