#ifndef FILE_RECURSIVE_POL
#define FILE_RECURSIVE_POL

/*********************************************************************/
/* File:   recursive_pol.hpp                                         */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


namespace ngfem
{

  /*
    Recursive Polynomials
  */



  // Use Lambda function with square-bracket assignment

  template <typename FUNC>
  class SBLambdaElement
  {
    FUNC f;
    int i;
  public:
    INLINE SBLambdaElement (FUNC af, int hi) : f(af), i(hi) { ; }
    template <typename VAL>
    INLINE VAL operator= (VAL v) { f(i, v); return v; }
  };

  template <typename FUNC>
  class Class_SBLambda
  {
    FUNC func;
    int offset;
  public:
    INLINE Class_SBLambda (FUNC f, int ao = 0) : func(f), offset(ao) { ; }
    INLINE SBLambdaElement<FUNC> operator[] (int i) const { return SBLambdaElement<FUNC> (func, offset+i); }
    INLINE Class_SBLambda<FUNC> operator+ (int i) const { return Class_SBLambda<FUNC> (func, offset+i); }
  };

  template <typename FUNC> 
  INLINE const Class_SBLambda<FUNC> SBLambda (FUNC f)
  {
    return Class_SBLambda<FUNC> (f);
  }
    






  /// a helper class for fixed order evaluation
  template <class REC, int N>
  class CEvalFO
  {
  public:
    template <class S, class T>
    ALWAYS_INLINE static void Eval (S x, T && values, S & p1, S & p2) 
    {
      S p3;
      CEvalFO<REC,N-1>::Eval (x, values, p2, p3);
      values[N] = p1 = ( REC::A(N) * x + REC::B(N)) * p2 + REC::C(N) * p3;
    }


    template <class S, class Sc, class T>
    ALWAYS_INLINE static void EvalMult (S x, Sc c, T && values, S & p1, S & p2) 
    {
      S p3;
      CEvalFO<REC,N-1>::EvalMult (x, c, values, p2, p3);
      values[N] = p1 = ( REC::A(N) * x + REC::B(N)) * p2 + REC::C(N) * p3;
    }


    template <class S, class Sy, class T>
    ALWAYS_INLINE static void EvalScaled (S x, Sy y, T && values, S & p1, S & p2) 
    {
      S p3;
      CEvalFO<REC,N-1>::EvalScaled (x, y, values, p2, p3);
      values[N] = p1 = ( REC::A(N) * x + REC::B(N) * y) * p2 + REC::C(N)*(y*y) * p3;
    }


    template <class S, class Sy, class Sc, class T>
    ALWAYS_INLINE static void EvalScaledMult (S x, Sy y, Sc c, T && values, S & p1, S & p2) 
    {
      S p3;
      CEvalFO<REC,N-1>::EvalScaledMult (x, y, c, values, p2, p3);
      values[N] = p1 = ( REC::A(N) * x + REC::B(N) * y) * p2 + REC::C(N)*(y*y) * p3;
    }
  };


  template <class REC>
  class CEvalFO<REC, -1>
    {
    public:
      template <class S, class T>
      ALWAYS_INLINE static void Eval (S x, T && values, S & /* p1 */, S & /* p2 */) 
      { ; }

      template <class S, class Sc, class T>
      ALWAYS_INLINE static void EvalMult (S x, Sc c, T && values, S & /* p1 */, S & /* p2 */) 
      { ; }


      template <class S, class Sy, class T>
      ALWAYS_INLINE static void EvalScaled (S x, Sy y, T && values, S & /* p1 */, S & /* p2 */) 
      { ; }

      template <class S, class Sy, class Sc, class T>
      ALWAYS_INLINE static void EvalScaledMult (S x, Sy y, Sc c, T && values, S & /* p1 */, S & /* p2 */) 
      { ; }

    };


  template <class REC>
  class CEvalFO<REC, 0>
    {
    public:
      template <class S, class T>
      ALWAYS_INLINE static void Eval (S x, T && values, S & p1, S & /* p2 */) 
      {
	values[0] = p1 = REC::P0(x);
      }

      template <class S, class Sc, class T>
      ALWAYS_INLINE static void EvalMult (S x, Sc c, T && values, S & p1, S & /* p2 */) 
      {
	values[0] = p1 = c * REC::P0(x);
      }


      template <class S, class Sy, class T>
      ALWAYS_INLINE static void EvalScaled (S x, Sy y, T && values, S & p1, S & /* p2 */) 
      {
	values[0] = p1 = REC::P0(x);
      }

      template <class S, class Sy, class Sc, class T>
      ALWAYS_INLINE static void EvalScaledMult (S x, Sy y, Sc c, T && values, S & p1, S & /* p2 */) 
      {
	values[0] = p1 = c * REC::P0(x);
      }

    };

  template <class REC>
  class CEvalFO<REC, 1>
  {
  public:
    template <class S, class T>
    ALWAYS_INLINE static void Eval (S x, T && values, S & p1, S & p2) 
    {
      values[0] = p2 = REC::P0(x);
      values[1] = p1 = REC::P1(x);
    }

    template <class S, class Sc, class T>
    ALWAYS_INLINE static void EvalMult (S x, Sc c, T && values, S & p1, S & p2) 
    {
      values[0] = p2 = c * REC::P0(x);
      values[1] = p1 = c * REC::P1(x);
    }

    template <class S, class Sy, class T>
    ALWAYS_INLINE static void EvalScaled (S x, Sy y, T && values, S & p1, S & p2) 
    {
      values[0] = p2 = REC::P0(x);
      values[1] = p1 = REC::P1(x);
    }

    template <class S, class Sy, class Sc, class T>
    ALWAYS_INLINE static void EvalScaledMult (S x, Sy y, Sc c, T && values, S & p1, S & p2) 
    {
      values[0] = p2 = c * REC::P0(x);
      values[1] = p1 = c * REC::P1(x);
    }
  };
  




  // P_i = (a_i x + b_i) P_{i-1} + c_i P_{i-2}
  template<class REC>
  class RecursivePolynomial
  {
  public:

    template <class S>
    ALWAYS_INLINE static S EvalNext (int i, S x, S & p1, S & p2)
    {
      switch (i)
        {
        case 0: return p1 = REC::P0(x);
        case 1: return p1 = REC::P1(x);
        default:
          {
            if (REC::ZERO_B)
              {
                p1 *= REC::C(i);
                p1 += REC::A(i) * x * p2;
                return p1;
              }
            else
              {
                p1 *= REC::C(i);
                p1 += (REC::A(i) * x + REC::B(i)) * p2;
                return p1;
              }
          }
        }
    }


    template <class S, class Sc>
    ALWAYS_INLINE static S EvalNextMult (int i, S x, Sc c, S & p1, S & p2)
    {
      switch (i)
        {
        case 0: return p1 = c * REC::P0(x);
        case 1: return p1 = c * REC::P1(x);
        default: return EvalNext (i, x, p1, p2);
        }
    }


    template <class S, class Sy>
    ALWAYS_INLINE static void EvalScaledNext (int i, S x, Sy y, S & p1, S & p2)
    {
      if (REC::ZERO_B)
        {
          S pnew = REC::A(i) * x * p1 + REC::C(i) * (y*y)*p2;
          p2 = p1;
          p1 = pnew;
        }
      else
        {
          S pnew = (REC::A(i) * x + REC::B(i) * y) * p1 + REC::C(i) * (y*y)*p2;
          p2 = p1;
          p1 = pnew;
        }
    }

    template <class S, class Sy>
    ALWAYS_INLINE static S EvalScaledNext2 (int i, S x, Sy y, S & p1, S & p2)
    {
      switch (i)
        {
        case 0: return p1 = REC::P0(x);
        case 1: return p1 = REC::P1(x);
        default:
          {
            if (REC::ZERO_B)
              {
                p1 *= REC::C(i) *y*y;
                p1 += REC::A(i) * x * p2;
                return p1;
              }
            else
              {
                p1 *= REC::C(i);
                p1 += (REC::A(i) * x + REC::B(i) * y) * p2;
                return p1;
              }
          }
        }
    }




  public:


    template <class S, class T>
    static void Eval (int n, S x, T && values) 
    {
      EvalMult (n, x, 1.0, values);
    }

    template <class S, class Sc, class T>
    static void EvalMult (int n, S x, Sc c, T && values) 
    {
      S p1, p2;

      if (n < 0) return;

      values[0] = EvalNextMult(0, x, c, p1, p2);
      if (n < 1) return;

      values[1] = EvalNextMult(1, x, c, p2, p1);
      if (n < 2) return;
      /*
      values[2] = EvalNext(2, x, p1, p2);
      if (n < 3) return;

      values[3] = EvalNext(3, x, p2, p1);
      if (n < 4) return;

      values[4] = EvalNext(4, x, p1, p2);
      if (n < 5) return;

      values[5] = EvalNext(5, x, p2, p1);
      if (n < 6) return;

      values[6] = EvalNext(6, x, p1, p2);
      if (n < 7) return;
     
      values[7] = EvalNext(7, x, p2, p1);
      if (n < 8) return;

      values[8] = EvalNext(8, x, p1, p2);
      if (n < 9) return;

      values[9] = EvalNext(9, x, p2, p1);
      if (n < 10) return;

      values[10] = EvalNext(10, x, p1, p2);
      if (n < 11) return;

      values[11] = EvalNext(11, x, p2, p1);
      if (n < 12) return;
      */
      int i = 2;
      for ( ; i < n; i+=2)
	{	
	  values[i] = EvalNext (i, x, p1, p2);
	  values[i+1] = EvalNext (i+1, x, p2, p1);
	}
      if (i <= n)
        values[i] = EvalNext (i, x, p1, p2);
    }





    template <class S, class Sy, class T>
    static void EvalScaled (int n, S x, Sy y, T && values)
    {
      EvalScaledMult (n, x, y, 1.0, values);
    }

    template <class S, class Sy, class Sc, class T>
    static void EvalScaledMult (int n, S x, Sy y, Sc c, T && values)
    {
      S p1, p2;

      if (n < 0) return;

      values[0] = p2 = c * REC::P0(x);
      if (n < 1) return;

      values[1] = p1 = c * REC::P1(x);
      if (n < 2) return;

      EvalScaledNext(2, x, y, p1, p2);
      values[2] = p1;
      if (n < 3) return;

      EvalScaledNext(3, x, y, p1, p2);
      values[3] = p1;
      if (n < 4) return;

      EvalScaledNext(4, x, y, p1, p2);
      values[4] = p1;
      if (n < 5) return;

      EvalScaledNext(5, x, y, p1, p2);
      values[5] = p1;
      if (n < 6) return;

      EvalScaledNext(6, x, y, p1, p2);
      values[6] = p1;
      if (n < 7) return;

      EvalScaledNext(7, x, y, p1, p2);
      values[7] = p1;
      if (n < 8) return;

      EvalScaledNext(8, x, y, p1, p2);
      values[8] = p1;

      for (int i = 9; i <= n; i++)
	{	
	  EvalScaledNext (i, x, y, p1, p2);
	  values[i] = p1;
	}
    }






    template <int N, class S, class T>
    ALWAYS_INLINE static void EvalFO (S x, T && values) 
    {
      S p1, p2;
      CEvalFO<REC, N>::Eval (x, values, p1, p2);
    }

    template <int N, class S, class Sc, class T>
    ALWAYS_INLINE static void EvalMultFO (S x, Sc c, T && values) 
    {
      S p1, p2;
      CEvalFO<REC, N>::EvalMult (x, c, values, p1, p2);
    }

    template <int N, class S, class Sy, class T>
    ALWAYS_INLINE static void EvalScaledFO (S x, Sy y, T && values) 
    {
      S p1, p2;
      CEvalFO<REC, N>::EvalScaled (x, y, values, p1, p2);
    }

    template <int N, class S, class Sy, class Sc, class T>
    ALWAYS_INLINE static void EvalScaledMultFO (S x, Sy y, Sc c,T && values) 
    {
      S p1, p2;
      CEvalFO<REC, N>::EvalScaledMult (x, y, c, values, p1, p2);
    }

    enum { ZERO_B = 0 };
  };







  // P_i = ( (a_i x + b_i) P_{i-1} + c_i P_{i-2} ) / d_i
  template<class REC>
  class RecursivePolynomialNonStatic
  {
  public:

    ALWAYS_INLINE void ABC (int i, double & a, double & b, double & c) const   
    {
      a = static_cast<const REC&>(*this).A(i);
      b = static_cast<const REC&>(*this).B(i);
      c = static_cast<const REC&>(*this).C(i); 
      double d = 1.0 / static_cast<const REC&>(*this).D(i);
      a *= d; b *= d; c *= d;
    }

    template <class S>
    ALWAYS_INLINE S EvalNext (int i, S x, S & p1, S & p2)
    {
      switch (i)
        {
        case 0: return p1 = static_cast<const REC&>(*this).P0(x);
        case 1: return p1 = static_cast<const REC&>(*this).P1(x);
        default:
          {
            if (REC::ZERO_B)
              {
                p1 *= static_cast<const REC&>(*this).C(i);
                p1 += static_cast<const REC&>(*this).A(i) * x * p2;
                p1 *= 1.0 / static_cast<const REC&>(*this).D(i);
                return p1;
              }
            else
              {
                double a, b, c;
                static_cast<const REC&> (*this).ABC (i, a, b, c);
                
                p1 *= c;
                p1 += (a * x + b) * p2;
                return p1;
              }
          }
        }
    }


    template <class S, class Sc>
    ALWAYS_INLINE S EvalNextMult (int i, S x, Sc c, S & p1, S & p2)
    {
      switch (i)
        {
        case 0: return p1 = c * static_cast<const REC&>(*this).P0(x);
        case 1: return p1 = c * static_cast<const REC&>(*this).P1(x);
        default: return EvalNext (i, x, p1, p2);
        }
    }


    template <class S, class Sy>
    ALWAYS_INLINE void EvalScaledNext (int i, S x, Sy y, S & p1, S & p2)
    {
      if (REC::ZERO_B)
        {
          S pnew = static_cast<const REC&>(*this).A(i) * x * p1 + static_cast<const REC&>(*this).C(i) * (y*y)*p2;
          pnew *= 1.0 / static_cast<const REC&>(*this).D(i);
          p2 = p1;
          p1 = pnew;
        }
      else
        {
          S pnew = (static_cast<const REC&>(*this).A(i) * x + static_cast<const REC&>(*this).B(i) * y) * p1 + static_cast<const REC&>(*this).C(i) * (y*y)*p2;
          pnew *= 1.0 / static_cast<const REC&>(*this).D(i);
          p2 = p1;
          p1 = pnew;
        }
    }


  public:


    template <class S, class T>
    INLINE void Eval (int n, S x, T && values) 
    {
      EvalMult (n, x, 1.0, values);
    }

    template <class S, class Sc, class T>
    INLINE void EvalMult (int n, S x, Sc c, T && values) 
    {
      S p1, p2;

      if (n < 0) return;

      values[0] = EvalNextMult(0, x, c, p1, p2);
      if (n < 1) return;

      values[1] = EvalNextMult(1, x, c, p2, p1);
      if (n < 2) return;

      for (int i = 2; i <= n; i+=2)
	{	
	  values[i] = EvalNext (i, x, p1, p2);
          if (i == n) break;
	  values[i+1] = EvalNext (i+1, x, p2, p1);
	}
    }



    template <class S, class Sy, class T>
    void EvalScaled (int n, S x, Sy y, T && values)
    {
      EvalScaledMult (n, x, y, 1.0, values);
    }

    template <class S, class Sy, class Sc, class T>
    void EvalScaledMult (int n, S x, Sy y, Sc c, T && values)
    {
      S p1, p2;

      if (n < 0) return;

      values[0] = p2 = c * static_cast<const REC&>(*this).P0(x);
      if (n < 1) return;

      values[1] = p1 = c * static_cast<const REC&>(*this).P1(x);
      if (n < 2) return;

      EvalScaledNext(2, x, y, p1, p2);
      values[2] = p1;
      if (n < 3) return;

      EvalScaledNext(3, x, y, p1, p2);
      values[3] = p1;
      if (n < 4) return;

      EvalScaledNext(4, x, y, p1, p2);
      values[4] = p1;
      if (n < 5) return;

      EvalScaledNext(5, x, y, p1, p2);
      values[5] = p1;
      if (n < 6) return;

      EvalScaledNext(6, x, y, p1, p2);
      values[6] = p1;
      if (n < 7) return;

      EvalScaledNext(7, x, y, p1, p2);
      values[7] = p1;
      if (n < 8) return;

      EvalScaledNext(8, x, y, p1, p2);
      values[8] = p1;

      for (int i = 9; i <= n; i++)
	{	
	  EvalScaledNext (i, x, y, p1, p2);
	  values[i] = p1;
	}
    }

    enum { ZERO_B = 0 };
  };









  /* ******************** Legendre Polynomial ****************** */

  class LegendrePolynomial_Old : public RecursivePolynomial<LegendrePolynomial_Old>
  {
  public:
    LegendrePolynomial_Old () { ; }

    template <class S, class T>
    inline LegendrePolynomial_Old (int n, S x, T && values)
    {
      Eval (n, x, values);
    }

    template <class S>
    static ALWAYS_INLINE S P0(S x)  { return S(1.0); }
    template <class S>
    static ALWAYS_INLINE S P1(S x)  { return x; }
    
    static ALWAYS_INLINE double A (int i) { return 2.0-1.0/i; }
    static ALWAYS_INLINE double B (int i) { return 0; }
    static ALWAYS_INLINE double C (int i) { return 1.0/i-1.0; }
    enum { ZERO_B = 1 };
  };
    

  class LegendrePolynomial : public RecursivePolynomial<LegendrePolynomial>
  {
    static Array< double[2] > coefs;
    
  public:
    LegendrePolynomial () { ; }

    template <class S, class T>
    inline LegendrePolynomial (int n, S x, T && values)
    { 
      Eval (n, x, values);
    }

    static void Calc (int n);

    template <class S>
    static ALWAYS_INLINE S P0(S x)  { return S(1.0); }
    template <class S>
    static ALWAYS_INLINE S P1(S x)  { return x; }
    
    static ALWAYS_INLINE double A (int i) { return coefs[i][0]; } // 2.0-1.0/i; 
    static ALWAYS_INLINE double B (int i) { return 0; }
    static ALWAYS_INLINE double C (int i) { return coefs[i][1]; } // 1.0/i-1.0; 
    enum { ZERO_B = 1 };
  };


    

  /* ******************** Jacobi Polynomial  (with fixed alpha, beta)  ************ */

  template <int al, int be>
  class JacobiPolynomialFix : public RecursivePolynomial<JacobiPolynomialFix<al,be> >
  {
  public:
    JacobiPolynomialFix () { ; }

    template <class S, class T>
    inline JacobiPolynomialFix (int n, S x, T && values)
    {
      Eval (n, x, values);
    }

    // typedef RecursivePolynomial<JacobiPolynomialFix<al,be> > BASE;
    // using BASE::EvalFO;
    
    template <class S>
    static ALWAYS_INLINE S P0(S x) { return S(1.0); }
    template <class S>
    static ALWAYS_INLINE S P1(S x) { return 0.5 * (2*(al+1)+(al+be+2)*(x-1)); }
      
    static ALWAYS_INLINE double A (int i) 
    { i--; return (2.0*i+al+be)*(2*i+al+be+1)*(2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    static ALWAYS_INLINE double B (int i)
    { i--; return (2.0*i+al+be+1)*(al*al-be*be) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    static ALWAYS_INLINE double C (int i) 
    { i--; return -2.0*(i+al)*(i+be) * (2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
  };



  class JacobiPolynomial2 : public RecursivePolynomialNonStatic<JacobiPolynomial2>
  {
    double al, be;
  public:
    JacobiPolynomial2 (double a, double b) : al(a), be(b) { ; }

    template <class S, class T>
    inline JacobiPolynomial2 (int n, S x, double alpha, double beta, T && values)
      : al(alpha), be(beta)
    {
      Eval (n, x, values);
    }

    template <class S>
    ALWAYS_INLINE S P0(S x) const { return S(1.0); }
    template <class S>
    ALWAYS_INLINE S P1(S x) const { return 0.5 * (2*(al+1)+(al+be+2)*(x-1)); }
      
    ALWAYS_INLINE double A (int i) const
    { i--; return (2.0*i+al+be)*(2*i+al+be+1)*(2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    ALWAYS_INLINE double B (int i) const
    { i--; return (2.0*i+al+be+1)*(al*al-be*be) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    ALWAYS_INLINE double C (int i) const
    { i--; return -2.0*(i+al)*(i+be) * (2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    ALWAYS_INLINE double D (int i) const { return 1; }
  };




  class JacobiPolynomial3 : public RecursivePolynomialNonStatic<JacobiPolynomial3>
  {
    double al, be;
    Vec<4> pola;
    Vec<4> polc;
    Vec<4> pold;
    Vec<2> polb;
  public:
    JacobiPolynomial3 (double a, double b) : al(a), be(b) 
    { 
      ;
    }

    template <class S, class T>
    inline JacobiPolynomial3 (int n, S x, double alpha, double beta, T && values)
      : al(alpha), be(beta)
    {
      double ab = al+be;
      // A(i) = ( pola(3) i^3 + pola(2) i^2 + pola(1) i + pola(0) )  / pold
      // ...

      pola(0) = ab*(ab+1)*(ab+2);
      pola(1) = 6*ab*ab+12*ab+4;
      pola(2) = 12*ab+12;
      pola(3) = 8;

      polb(0) = (ab+1)*(al*al-be*be);
      polb(1) = 2*(al*al-be*be);

      polc(0) = -2*al*be*(ab+2);
      polc(1) = -2*ab*(ab+2)-4*al*be;
      polc(2) = -6*ab-4;
      polc(3) = -4;

      pold(0) = 2*(ab+1)*ab;
      pold(1) = 2*ab*ab+8*ab+4;
      pold(2) = 6*ab+8;
      pold(3) = 4;

      Eval (n, x, values);
    }

    template <class S>
    ALWAYS_INLINE S P0(S x) const { return S(1.0); }
    template <class S>
    ALWAYS_INLINE S P1(S x) const { return 0.5 * (2*(al+1)+(al+be+2)*(x-1)); }

    ALWAYS_INLINE double A (int i) const
    { i--; return ( pola(0) + i * (pola(1) + i * (pola(2) + i * pola(3))) ); }
    ALWAYS_INLINE double B (int i) const
    { i--; return ( polb(0) + i * polb(1) ); }
    ALWAYS_INLINE double C (int i) const
    { i--; return ( polc(0) + i * (polc(1) + i * (polc(2) + i * polc(3))) ); }

    ALWAYS_INLINE double D (int i) const
    { i--; return ( pold(0) + i * (pold(1) + i * (pold(2) + i * pold(3))) ); }
  };











  /*
  class JacobiPolynomial2 : public RecursivePolynomial<JacobiPolynomial2>
  {
    double al, be;
  public:
    JacobiPolynomial2 (double aal, double abe) : al(aal), be(abe) { ; }

    template <class S>
    S P0(S x) const { return S(1.0); }
    template <class S>
    S P1(S x) const { return 0.5 * (2*(al+1)+(al+be+2)*(x-1)); }
      
    double A (int i) const 
    { i--; return (2.0*i+al+be)*(2*i+al+be+1)*(2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    double B (int i) const
    { i--; return (2.0*i+al+be+1)*(al*al-be*be) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    double C (int i) const 
    { i--; return -2.0*(i+al)*(i+be) * (2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
  };
  */





  template <class S, class Sc, class T>
  inline void LegendrePolynomialMult (int n, S x, Sc c , T && values)
  {
    LegendrePolynomial leg;
    leg.EvalMult (n, x, c, values);
  }




  template <class S, class T>
  inline void JacobiPolynomial (int n, S x, double alpha, double beta, T && values)
  {
    S p1 = 1.0, p2 = 0.0, p3;

    if (n >= 0) 
      values[0] = p2 = 1.0;
    if (n >= 1) 
      values[1] = p1 = 0.5 * (2*(alpha+1)+(alpha+beta+2)*(x-1));

    for (int i  = 1; i < n; i++)
      {
	p3 = p2; p2=p1;
	p1 =
	  1.0 / ( 2 * (i+1) * (i+alpha+beta+1) * (2*i+alpha+beta) ) *
	  ( 
	   ( (2*i+alpha+beta+1)*(alpha*alpha-beta*beta) + 
	     (2*i+alpha+beta)*(2*i+alpha+beta+1)*(2*i+alpha+beta+2) * x) 
	   * p2
	   - 2*(i+alpha)*(i+beta) * (2*i+alpha+beta+2) * p3
	   );
	values[i+1] = p1;
      }
  }





  template <class S, class Sc, class T>
  inline void JacobiPolynomialMult (int n, S x, double alpha, double beta, Sc c, T && values)
  {
    S p1 = c, p2 = 0.0, p3;

    if (n >= 0) 
      p2 = values[0] = c;
    if (n >= 1) 
      p1 = values[1] = 0.5 * c * (2*(alpha+1)+(alpha+beta+2)*(x-1));

    for (int i  = 1; i < n; i++)
      {
	p3 = p2; p2=p1;
	p1 =
	  1.0 / ( 2 * (i+1) * (i+alpha+beta+1) * (2*i+alpha+beta) ) *
	  ( 
	   ( (2*i+alpha+beta+1)*(alpha*alpha-beta*beta) + 
	     (2*i+alpha+beta)*(2*i+alpha+beta+1)*(2*i+alpha+beta+2) * x) 
	   * p2
	   - 2*(i+alpha)*(i+beta) * (2*i+alpha+beta+2) * p3
	   );
	values[i+1] = p1;
      }
  }







  template <class S, class St, class T>
  inline void ScaledJacobiPolynomial (int n, S x, St t, double alpha, double beta, T && values)
  {
    /*
      S p1 = 1.0, p2 = 0.0, p3;

      if (n >= 0) values[0] = 1.0;
    */

    S p1 = 1.0, p2 = 0.0, p3;

    if (n >= 0) 
      p2 = values[0] = 1.0;
    if (n >= 1) 
      p1 = values[1] = 0.5 * (2*(alpha+1)*t+(alpha+beta+2)*(x-t));

    for (int i=1; i < n; i++)
      {
	p3 = p2; p2=p1;
	p1 =
	  1.0 / ( 2 * (i+1) * (i+alpha+beta+1) * (2*i+alpha+beta) ) *
	  ( 
	   ( (2*i+alpha+beta+1)*(alpha*alpha-beta*beta) * t + 
	     (2*i+alpha+beta)*(2*i+alpha+beta+1)*(2*i+alpha+beta+2) * x) 
	   * p2
	   - 2*(i+alpha)*(i+beta) * (2*i+alpha+beta+2) * t * t * p3
	   );
	values[i+1] = p1;
      }
  }






  template <class S, class St, class Sc, class T>
  inline void ScaledJacobiPolynomialMult (int n, S x, St t, double alpha, double beta, Sc c, T && values)
  {
    /*
      S p1 = 1.0, p2 = 0.0, p3;
      if (n >= 0) values[0] = 1.0;
    */

    S p1 = c, p2 = 0.0, p3;

    if (n >= 0) 
      p2 = values[0] = c;
    if (n >= 1) 
      p1 = values[1] = 0.5 * c * (2*(alpha+1)*t+(alpha+beta+2)*(x-t));

    for (int i=1; i < n; i++)
      {
	p3 = p2; p2=p1;
	p1 =
	  1.0 / ( 2 * (i+1) * (i+alpha+beta+1) * (2*i+alpha+beta) ) *
	  ( 
	   ( (2*i+alpha+beta+1)*(alpha*alpha-beta*beta) * t + 
	     (2*i+alpha+beta)*(2*i+alpha+beta+1)*(2*i+alpha+beta+2) * x) 
	   * p2
	   - 2*(i+alpha)*(i+beta) * (2*i+alpha+beta+2) * t * t * p3
	   );
	values[i+1] = p1;
      }
  }















  template <int n, int i, int alpha0, int beta>
  class DubinerJacobiPolynomialsFO
  {
  public:
    template <class S, class T>
    static void Eval (S x, T && values)
    {
      DubinerJacobiPolynomialsFO<n, i-1, alpha0, beta>::Eval (x, values);

      S p1, p2;
      CEvalFO<JacobiPolynomialFix<alpha0+2*i, beta>, n-i>::Eval (x, values.Row(i), p1, p2);

      // why not ???
      // JacobiPolynomialFix<alpha0+2*i, beta>::EvalFO<n-i> (x, values.Row(i));


      // typedef JacobiPolynomialFix<alpha0+2*i, beta> JAC;
      // JAC::EvalFO<n-i> (x, values.Row(i));
      // JacobiPolynomialFix<1,2>::EvalFO<n-i> (hx, hmat.Row(i));
    }


    template <class S, class St, class T>
    static void EvalScaled (S x, St t, T && values)
    {
      DubinerJacobiPolynomialsFO<n, i-1, alpha0, beta>::EvalScaled (x, t, values);
      // JacobiPolynomialFO<n-i, alpha0+2*i, beta>::Eval (x, values.Row(i));

      JacobiPolynomialFix<alpha0+2*i, beta> jac;
      S p1, p2;
      CEvalFO<JacobiPolynomialFix<alpha0+2*i, beta>, n-i>::EvalScaled (x, t, values.Row(i), p1, p2);
    }    
  };
  
  template <int n, int alpha0, int beta>
  class DubinerJacobiPolynomialsFO<n, -1, alpha0, beta>
  {
  public:
    template <class S, class T>
    static void Eval (S x, T && values)
    { ; }
    template <class S, class St, class T>
    static void EvalScaled (S x, St t, T && values)
    { ; }
  };
  
  


  template <int n, int i, int alpha0, int beta>
  class DubinerJacobiPolynomialsPowFO
  {
  public:
    template <class S, class T>
    static S Eval (S x, T && values)
    {
      S power = DubinerJacobiPolynomialsPowFO<n, i-1, alpha0, beta>::Eval (x, values);
      S p1, p2;
      CEvalFO<JacobiPolynomialFix<alpha0+2*i, beta>, n-i>::EvalMult (x, power, values.Row(i), p1, p2);
      return power * (1-x)/2;
    }


    template <class S, class St, class T>
    static S EvalScaled (S x, St t, T && values)
    {
      S power = DubinerJacobiPolynomialsPowFO<n, i-1, alpha0, beta>::EvalScaled (x, t, values);
      S p1, p2;
      CEvalFO<JacobiPolynomialFix<alpha0+2*i, beta>, n-i>::EvalScaledMult (x, t, power, values.Row(i), p1, p2);
      return power * (1-x)/2;
    }
  };
  
  template <int n, int alpha0, int beta>
  class DubinerJacobiPolynomialsPowFO<n, -1, alpha0, beta>
  {
  public:
    template <class S, class T>
    static S Eval (S x, T && values)
    { return 1.0; }
    template <class S, class St, class T>
    static S EvalScaled (S x, St t, T && values)
    { return 1.0; }
  };
















#ifdef OLD
  template <int ALPHA0, int BETA, class S, class T>
  void DubinerJacobiPolynomials2 (int n, S x, T && values)
  {
    switch (n)
      {
      case 0: DubinerJacobiPolynomialsFO<0, 0, ALPHA0, BETA>::Eval (x, values); break;
      case 1: DubinerJacobiPolynomialsFO<1, 1, ALPHA0, BETA>::Eval (x, values); break;
      case 2: DubinerJacobiPolynomialsFO<2, 2, ALPHA0, BETA>::Eval (x, values); break;
      case 3: DubinerJacobiPolynomialsFO<3, 3, ALPHA0, BETA>::Eval (x, values); break;
      case 4: DubinerJacobiPolynomialsFO<4, 4, ALPHA0, BETA>::Eval (x, values); break;
      case 5: DubinerJacobiPolynomialsFO<5, 5, ALPHA0, BETA>::Eval (x, values); break;
      case 6: DubinerJacobiPolynomialsFO<6, 6, ALPHA0, BETA>::Eval (x, values); break;
      case 7: DubinerJacobiPolynomialsFO<7, 7, ALPHA0, BETA>::Eval (x, values); break;
      case 8: DubinerJacobiPolynomialsFO<8, 8, ALPHA0, BETA>::Eval (x, values); break;
      case 9: DubinerJacobiPolynomialsFO<9, 9, ALPHA0, BETA>::Eval (x, values); break;
      case 10: DubinerJacobiPolynomialsFO<10, 10, ALPHA0, BETA>::Eval (x, values); break;
      default: DubinerJacobiPolynomials1 (n, x, ALPHA0, BETA, values);
      }
  }
#endif





  // compute J_j^{2i+alpha0, beta} (x),  for i+j <= n

  template <class S, class T>
  void DubinerJacobiPolynomials1 (int n, S x, int alpha0, int beta, T && values)
  {
    for (int i = 0; i <= n; i++)
      JacobiPolynomial (n-i, x, alpha0+2*i, beta, values.Row(i));
  }


  template <class S, class St, class T>
  void DubinerJacobiPolynomialsScaled1 (int n, S x, St t, int alpha0, int beta, T && values)
  {
    for (int i = 0; i <= n; i++)
      ScaledJacobiPolynomial (n-i, x, t, alpha0+2*i, beta, values.Row(i));
  }







  template <int ALPHA0, int BETA, int DIAG, int ORDER = DIAG>
  class DubinerJacobiPolynomialsDiag
  {
  public:
    template<class S, class Thelp, class T>
    ALWAYS_INLINE static void Step (S x, Thelp & help, T & values)
    {
      DubinerJacobiPolynomialsDiag<ALPHA0, BETA, DIAG, ORDER-1>::Step (x, help, values);
      typedef JacobiPolynomialFix<ALPHA0+2*(DIAG-ORDER), BETA> REC;
      
      if (ORDER == 0)
	help[DIAG-ORDER][0] = REC::P0(x);
      else if (ORDER == 1)
	{
	  help[DIAG-ORDER][0] = REC::P1(x);
	  help[DIAG-ORDER][1] = REC::P0(x);
	}
      else
	{
          Swap (help[DIAG-ORDER][0], help[DIAG-ORDER][1]);
	  REC::EvalNext (ORDER, x, help[DIAG-ORDER][0], help[DIAG-ORDER][1]);
	}
      values(DIAG-ORDER, ORDER) = help[DIAG-ORDER][0];
    }
  };

  template <int ALPHA0, int BETA, int DIAG>
  class DubinerJacobiPolynomialsDiag<ALPHA0, BETA, DIAG, -1>
  {
  public:
    template<class S, class Thelp, class T>
    ALWAYS_INLINE static void Step (S x, Thelp & help, T & values) {;}
  };

  
  template <int ALPHA0, int BETA, class S, class T>
  void DubinerJacobiPolynomials (int n, S x, T && values)
  {
    S help[20][2];
    if (n < 0) return;
    DubinerJacobiPolynomialsDiag<ALPHA0, BETA, 0>::Step (x, help, values);

    if (n < 1) return;
    DubinerJacobiPolynomialsDiag<ALPHA0, BETA, 1>::Step (x, help, values);

    if (n < 2) return;
    DubinerJacobiPolynomialsDiag<ALPHA0, BETA, 2>::Step (x, help, values);

    if (n < 3) return;
    DubinerJacobiPolynomialsDiag<ALPHA0, BETA, 3>::Step (x, help, values);

    if (n < 4) return;
    DubinerJacobiPolynomialsDiag<ALPHA0, BETA, 4>::Step (x, help, values);

    if (n < 5) return;
    DubinerJacobiPolynomialsDiag<ALPHA0, BETA, 5>::Step (x, help, values);

    if (n < 6) return;
    DubinerJacobiPolynomialsDiag<ALPHA0, BETA, 6>::Step (x, help, values);

    if (n < 7) return;
    DubinerJacobiPolynomialsDiag<ALPHA0, BETA, 7>::Step (x, help, values);

    if (n < 8) return;
    DubinerJacobiPolynomialsDiag<ALPHA0, BETA, 8>::Step (x, help, values);

    if (n < 9) return;
    DubinerJacobiPolynomialsDiag<ALPHA0, BETA, 9>::Step (x, help, values);

    if (n < 10) return;
    DubinerJacobiPolynomialsDiag<ALPHA0, BETA, 10>::Step (x, help, values);

    if (n < 11) return;

    DubinerJacobiPolynomials1 (n, x, ALPHA0, BETA, values);
  }













  template <int ALPHA0, int BETA, int DIAG, int ORDER = DIAG>
  class DubinerJacobiPolynomialsScaledDiag
  {
  public:
    template<class S, class St, class Thelp, class T>
    ALWAYS_INLINE static void Step (S x, St t, Thelp & help, T & values)
    {
      DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, DIAG, ORDER-1>::Step (x, t, help, values);
      typedef JacobiPolynomialFix<ALPHA0+2*(DIAG-ORDER), BETA> REC;
      
      if (ORDER == 0)
	help[DIAG-ORDER][0] = REC::P0(x);
      else if (ORDER == 1)
	{
	  help[DIAG-ORDER][0] = REC::P1(x);
	  help[DIAG-ORDER][1] = REC::P0(x);
	}
      else
	{
	  REC::EvalScaledNext (ORDER, x, t, help[DIAG-ORDER][0], help[DIAG-ORDER][1]);
	}
      values(DIAG-ORDER, ORDER) = help[DIAG-ORDER][0];
    }
  };

  template <int ALPHA0, int BETA, int DIAG>
  class DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, DIAG, -1>
  {
  public:
    template<class S, class St, class Thelp, class T>
    ALWAYS_INLINE static void Step (S x, St t, Thelp & help, T & values) {;}
  };

  
  template <int ALPHA0, int BETA, class S, class St, class T>
  void DubinerJacobiPolynomialsScaled (int n, S x, St t, T && values)
  {
    S help[20][2];
    if (n < 0) return;
    DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, 0>::Step (x, t, help, values);

    if (n < 1) return;
    DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, 1>::Step (x, t, help, values);

    if (n < 2) return;
    DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, 2>::Step (x, t, help, values);

    if (n < 3) return;
    DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, 3>::Step (x, t, help, values);

    if (n < 4) return;
    DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, 4>::Step (x, t, help, values);

    if (n < 5) return;
    DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, 5>::Step (x, t, help, values);

    if (n < 6) return;
    DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, 6>::Step (x, t, help, values);

    if (n < 7) return;
    DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, 7>::Step (x, t, help, values);

    if (n < 8) return;
    DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, 8>::Step (x, t, help, values);

    if (n < 9) return;
    DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, 9>::Step (x, t, help, values);

    if (n < 10) return;
    DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, 10>::Step (x, t, help, values);

    if (n < 11) return;

    DubinerJacobiPolynomialsScaled1 (n, x, t, ALPHA0, BETA, values);
  }










  inline int TrigIndex(int i, int j, int n)
  {
    return j + i*(2*n+3-i)/2;
  }



  template <int ALPHA0, int BETA, int DIAG, int ORDER = DIAG>
  class DubinerJacobiPolynomialsDiag_Linear
  {
  public:
    template<class S, class Thelp, class T>
    ALWAYS_INLINE static void Step (S x, Thelp & help, T & values)
    {
      DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, DIAG, ORDER-1>::Step (x, help, values);
      typedef JacobiPolynomialFix<ALPHA0+2*(DIAG-ORDER), BETA> REC;
      
      if (ORDER == 0)
	help(DIAG-ORDER,0) *= REC::P0(x);
      else if (ORDER == 1)
	{
	  help(DIAG-ORDER,1) = help(DIAG-ORDER,0);
	  help(DIAG-ORDER,0) *= REC::P1(x) / REC::P0(x);
	}
      else
	{
          Swap (help(DIAG-ORDER,0), help(DIAG-ORDER,1));
	  REC::EvalNext (ORDER, x, help(DIAG-ORDER,0), help(DIAG-ORDER,1));
	}
      values[TrigIndex(ORDER,DIAG-ORDER,help.Height()-1)] = help(DIAG-ORDER,0);
    }
  };

  template <int ALPHA0, int BETA, int DIAG>
  class DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, DIAG, -1>
  {
  public:
    template<class S, class Thelp, class T>
    ALWAYS_INLINE static void Step (S x, Thelp & help, T & values) {;}
  };

  
  template<class S, class Thelp, class T>
  void DubinerJacobiPolynomialsDiag_Linear1 (double alpha0, double beta0, 
					     int diag, S x, Thelp & help, T & values)
  {
    for (int order = 0; order <= diag; order++)
      {
	double al = alpha0 + 2*(diag-order);
	double be = beta0;

	if (order == 0)
	  help(diag-order,0) *= 1; // REC::P0(x);
	else if (order == 1)
	  {
	    help(diag-order,1) = help(diag-order,0);
	    help(diag-order,0) *= 0.5 * (2*(al+1)+(al+be+2)*(x-1));
	  }
	else
	  {
	    // REC::EvalNext (ORDER, x, help(DIAG-ORDER,0), help(DIAG-ORDER,1));
	    // S pnew = (REC::A(i) * x + REC::B(i)) * p1 + REC::C(i) * p2;

	    S p1 = help(diag-order,0);
	    S p2 = help(diag-order,1);
	    
	    int i = order-1;
	    S pnew =
	      1.0 / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be) ) *
	      ( 
	       ( (2*i+al+be+1)*(al*al-be*be) + 
		 (2*i+al+be)*(2*i+al+be+1)*(2*i+al+be+2) * x) 
	       * p1
	       - 2*(i+al)*(i+be) * (2*i+al+be+2) * p2
	       );

	    help(diag-order,0) = pnew;
	    help(diag-order,1) = p1;
	  }
	values[TrigIndex(order,diag-order,help.Height()-1)] = help(diag-order,0);
      
      }
  }

  template <int ALPHA0, int BETA, class S, class Sy, class Sc, class T>
  void DubinerJacobiPolynomials_Linear (int n, S x, Sy y, Sc c, T & values)
  {
    ArrayMem<S,40> hmem(2*(n+1));
    FlatMatrixFixWidth<2,S> help(n+1, &hmem[0]);

    LegendrePolynomial leg;
    leg.EvalScaledMult (n, 2*y+x-1, 1-x, c, help.Col(0));

    x = 2*x-1;
    if (n >= 0) DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, 0>::Step (x, help, values);
    if (n >= 1) DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, 1>::Step (x, help, values);
    if (n >= 2) DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, 2>::Step (x, help, values);
    if (n >= 3) DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, 3>::Step (x, help, values);
    if (n >= 4) DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, 4>::Step (x, help, values);
    if (n >= 5) DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, 5>::Step (x, help, values);
    if (n >= 6) DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, 6>::Step (x, help, values);
    if (n >= 7) DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, 7>::Step (x, help, values);
    if (n >= 8) DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, 8>::Step (x, help, values);
    
    for (int i = 9; i <= n; i++)
      DubinerJacobiPolynomialsDiag_Linear1 (ALPHA0, BETA, i, x, help, values);
  }









  template <int ALPHA0, int BETA, int DIAG, int ORDER = DIAG>
  class DubinerJacobiPolynomialsScaledDiag_Linear
  {
  public:
    template<class S, class St, class Thelp, class T>
    ALWAYS_INLINE static void Step (S x, St t, Thelp & help, T & values)
    {
      DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, DIAG, ORDER-1>::Step (x, t, help, values);
      typedef JacobiPolynomialFix<ALPHA0+2*(DIAG-ORDER), BETA> REC;
      
      if (ORDER == 0)
	help(DIAG-ORDER,0) *= REC::P0(x);
      else if (ORDER == 1)
	{
	  help(DIAG-ORDER,1) = help(DIAG-ORDER,0);
	  help(DIAG-ORDER,0) *= REC::P1(x) / REC::P0(x);
	}
      else
	{
	  REC::EvalScaledNext (ORDER, x, t, help(DIAG-ORDER,0), help(DIAG-ORDER,1));
	}
      values[TrigIndex(ORDER,DIAG-ORDER,help.Height()-1)] = help(DIAG-ORDER,0);
    }
  };

  template <int ALPHA0, int BETA, int DIAG>
  class DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, DIAG, -1>
  {
  public:
    template<class S, class St, class Thelp, class T>
    ALWAYS_INLINE static void Step (S x, St t, Thelp & help, T & values) {;}
  };

  
  template<class S, class St, class Thelp, class T>
  void DubinerJacobiPolynomialsScaledDiag_Linear1 (double alpha0, double beta0, 
						   int diag, S x, St t, Thelp & help, T & values)
  {
    for (int order = 0; order <= diag; order++)
      {
	double al = alpha0 + 2*(diag-order);
	double be = beta0;

	if (order == 0)
	  help(diag-order,0) *= 1; // REC::P0(x);
	else if (order == 1)
	  {
	    help(diag-order,1) = help(diag-order,0);
	    help(diag-order,0) *= 0.5 * (2*(al+1)+(al+be+2)*(x-1));
	  }
	else
	  {
	    S p1 = help(diag-order,0);
	    S p2 = help(diag-order,1);
	    
	    int i = order-1;
	    S pnew =
	      1.0 / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be) ) *
	      ( 
	       ( (2*i+al+be+1)*(al*al-be*be) * t + 
		 (2*i+al+be)*(2*i+al+be+1)*(2*i+al+be+2) * x) 
	       * p1
	       - 2*(i+al)*(i+be) * (2*i+al+be+2) * t*t*p2
	       );

	    help(diag-order,0) = pnew;
	    help(diag-order,1) = p1;
	  }
	values[TrigIndex(order,diag-order,help.Height()-1)] = help(diag-order,0);
      
      }
  }

  template <int ALPHA0, int BETA, class S, class Sy, class Sc, class St, class T>
  void DubinerJacobiPolynomialsScaled_Linear (int n, S x, Sy y, St t, Sc c, T & values)
  {
    ArrayMem<S,40> hmem(2*(n+1));
    FlatMatrixFixWidth<2,S> help(n+1, &hmem[0]);

    LegendrePolynomial leg;
    leg.EvalScaledMult (n, 2*y+x-1, t-x, c, help.Col(0));

    x = 2*x-1;
    if (n >= 0) DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, 0>::Step (x, t, help, values);
    if (n >= 1) DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, 1>::Step (x, t, help, values);
    if (n >= 2) DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, 2>::Step (x, t, help, values);
    if (n >= 3) DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, 3>::Step (x, t, help, values);
    if (n >= 4) DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, 4>::Step (x, t, help, values);
    if (n >= 5) DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, 5>::Step (x, t, help, values);
    if (n >= 6) DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, 6>::Step (x, t, help, values);
    if (n >= 7) DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, 7>::Step (x, t, help, values);
    if (n >= 8) DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, 8>::Step (x, t, help, values);

    for (int i = 9; i <= n; i++)
      DubinerJacobiPolynomialsScaledDiag_Linear1 (ALPHA0, BETA, i, x, t, help, values);
  }




  class DubinerBasis
  {
  public:
    template <class S, class T>
    static void Eval (int n, S x, S y, T && values)
    {
      EvalMult (n, x, y, 1.0, values);
    }

    template <class S, class Sc, class T>
    static void EvalMult (int n, S x, S y, Sc c, T && values)
    {
      DubinerJacobiPolynomials_Linear<1,0> (n, x, y, c, values);
    }

    template <class S, class St, class Sc, class T>
    static void EvalScaledMult (int n, S x, S y, St t, Sc c, T && values)
    {
      DubinerJacobiPolynomialsScaled_Linear<1,0> (n, x, y, t, c, values);
    }
  };










  template <class S, class T>
  inline void GegenbauerPolynomial (int n, S x, double lam, T && values)
  {
    S p1 = 1.0, p2 = 0.0, p3;

    if (n >= 0)
      values[0] = 1.0;

    for (int j=1; j<=n; j++)
      {
	p3=p2; p2=p1;
	p1=( 2.0*(j+lam-1.0) * x * p2 - (j+2*lam-2.0) * p3) / j;
	values[j] = p1;
      }  
  }


  /**
     Integrated Legendre polynomials on (-1,1)

     value[0] = -1
     value[1] = x
     value[i] (x) = \int_{-1}^x P_{i-1} (s) ds  for i >= 2

     WARNING: is not \int P_i
  */
  template <class S, class T>
  inline void IntegratedLegendrePolynomial (int n, S x, T && values)
  {
    S p1 = -1.0;
    S p2 = 0.0; 
    S p3;

    if (n >= 0)
      values[0] = (S)-1.0;

    for (int j=1; j<=n; j++)
      {
	p3=p2; p2=p1;
	p1=( (2*j-3) * x * p2 - (j-3) * p3) / j;
	values[j] = p1;
    
      }
  }







  /**
     Hermite polynomials H_i, orthogonal w.r.t. \int_R exp(-x*x)
     up to order n --> n+1 values
   
     H_l =  2 x H_{l-1} - 2 (l-1)  H_{l-2}

     P_0 = 1
     P_1 = 2*x
     P_2 = 4*x*x - 2
     P_2 = 8*x*x*x - 12 x
  */

  template <class S, class T>
  inline void HermitePolynomial (int n, S x, T && values)
  {
    S p1, p2, p3;
  
    p2 = 0;
    if (n >= 0)
      p1 = values[0] = 1.0;
    for (int j=1; j<=n; j++)
      {
	p3 = p2; p2 = p1;
	p1 = 2*x*p2 - 2*(j-1)*p3;
	values[j] = p1;
      }
  }












  /**
     Compute triangle edge-shape functions

     functions vanish on upper two edges

     x,y: coordinates in triangle (-1, 0), (1, 0), (0, 1)

     f_i (x, 0) = IntegratedLegendrePol_i (x)

     f_i ... pol of order i


     Monomial extension:
  */
  template <class Sx, class Sy, class T>
  inline void TriangleExtensionMonomial (int n, Sx x, Sy y, T & values)
  {
    Sx p1 = -1.0, p2 = 0.0, p3;
    values[0] = -1.0;
    Sy fy = (1-y)*(1-y);
    for (int j=1; j<=n; j++)
      {
	p3=p2; p2=p1;
	p1=( (2*j-3) * x * p2 - (j-3) * fy * p3) / j;
	values[j] = p1;
      }    
  }

  template <class Sx, class Sy, class T>
  inline void DiffTriangleExtensionMonomial (int n, Sx x, Sy y, T & values)
  {
    Array<AutoDiff<2> > ad_values(n+1);
    AutoDiff<2> ad_x(x, 0);
    AutoDiff<2> ad_y(y, 1);

    TriangleExtensionMonomial (n, ad_x, ad_y, ad_values);

    for (int i = 0; i <= n; i++)
      for (int j = 0; j < 2; j++)
	values(i,j) = ad_values[i].DValue(j);
  }



  /**
     Extension is the optimal averaging extension:
  */
  template <class Sx, class Sy, class T>
  inline void TriangleExtensionJacobi (int n, Sx x, Sy y, T & values)
  {
    if ( (1-y) != 0.0)
      {
	int j;

	JacobiPolynomial (n-2, x / (1-y), 2, 2, values);
	Sy fac = (1.-x-y) * (1.+x-y);
	for (j = 0; j <= n-2; j++)
	  {
	    values[j] *= fac;
	    fac *= 1-y;
	  }
	for (j = n; j >= 2; j--)
	  values[j] = values[j-2];
	if (n >= 0) values[0] = 0;
	if (n >= 1) values[1] = 0;
      }
    else
      {
	for (int j = 0; j <= n; j++)
	  values[j] = 0;
      }
  }

  template <class Sx, class Sy, class T>
  inline void DiffTriangleExtensionJacobi (int n, Sx x, Sy y, T & values)
  {
    Array<AutoDiff<2> > ad_values(n+1);
    AutoDiff<2> ad_x(x, 0);
    AutoDiff<2> ad_y(y, 1);

    TriangleExtensionJacobi (n, ad_x, ad_y, ad_values);
    for (int i = 0; i <= n; i++)
      for (int j = 0; j < 2; j++)
	values(i,j) = ad_values[i].DValue(j);
  }







  /**
     Extension is the optimal averaging extension:
  */
  template <class Sx, class Sy, class T>
  inline void TriangleExtensionOpt (int n, Sx x, Sy y, T & values)
  {
    if (y < 1e-10)
      {
	IntegratedLegendrePolynomial (n, x, values);
      }
    else
      {
	Array<Sx> ge1(n+2);
	Array<Sx> ge2(n+2);
	Array<Sx> ge3(n+2);
	Array<Sx> ge4(n+2);

	GegenbauerPolynomial (n+1, Sx(-1.0), -1.5, ge1);
	GegenbauerPolynomial (n+1, x-y, -1.5, ge2);
	GegenbauerPolynomial (n+1, x+y, -1.5, ge3);
	GegenbauerPolynomial (n+1, Sx(1.0), -1.5, ge4);
 
	for (int i = 0; i <= n; i++)
	  values[i] = 1.0/3.0 *
	    (  (2*y/(1+x+y)/(1+x+y) - y/2) * ge1[i+1]  +
	       (-1/(2*y) + 2*y/(1-x+y)/(1-x+y)) * ge2[i+1] +
	       (1/(2*y) - 2*y/(1+x+y)/(1+x+y) ) * ge3[i+1] +
	       (-2*y/(1-x+y)/(1-x+y) + y/2 ) * ge4[i+1] );
      }
  }

  template <class Sx, class Sy, class T>
  inline void DiffTriangleExtensionOpt (int n, Sx x, Sy y, T & values)
  {
    Array<AutoDiff<2> > ad_values(n+1);
    AutoDiff<2> ad_x(x, 0);
    AutoDiff<2> ad_y(y, 1);

    if (y < 1e-10)
      {
	values = 0.;
      }
    else
      {
	TriangleExtensionOpt (n, ad_x, ad_y, ad_values);

	for (int i = 0; i <= n; i++)
	  for (int j = 0; j < 2; j++)
	    values(i,j) = ad_values[i].DValue(j);
      }
  }



  template <class S1, class S2, class S3>
  inline void StdOp (S1 & v1, const S2 & tt, const S3 & v2, double fac)
  {
    v1 = fac * (v1*tt - v2) + v2;
    // v1 = fac * (v1*tt) + (1-fac) * v2;
  }

  template <int D>
  inline void StdOp (AutoDiff<D> & v1, const AutoDiff<D> & tt, const AutoDiff<D> & v2, double fac)
  {
    for (int j = 0; j < D; j++)
      v1.DValue(j) = fac * (v1.DValue(j) * tt.Value() + v1.Value() * tt.DValue(j) - v2.DValue(j)) + v2.DValue(j);
    v1.Value() = fac * (v1.Value()*tt.Value()-v2.Value()) + v2.Value();
  }


  /* 
     E_i(x,y) = P_i(x/t) * t^i 
  */ 
  template <class Sx, class St, class T>
  inline void ScaledLegendrePolynomial (int n, Sx x, St t, T && values)
  {
    // St tt = t*t;
    St tt = sqr(t);

    Sx p1, p2;

    if (n < 0) return;
    values[0] = p2 = 1.0;
    if (n < 1) return;
    values[1] = p1 = x;
    if (n < 2) return;

    for (int j=2; j < n; j+=2)
      {
	/*
	  double invj = 1.0/j;
	  p2 *= (invj-1) * tt;
	  p2 += (2-invj) * x * p1;
	  values[j]   = p2; 

	  double invj2 = 1.0/(j+1);
	  p1 *= (invj2-1) * tt;
	  p1 += (2-invj2) * x * p2;
	  values[j+1] = p1; 
	*/

	StdOp (p2, tt, x*p1, 1.0/j-1);
	values[j]   = p2; 
	StdOp (p1, tt, x*p2, 1.0/(j+1)-1);
	values[j+1] = p1; 
      }

    if (n % 2 == 0)
      {
	double invn = 1.0/n;
	values[n] = (2-invn)*x*p1 - (1-invn) * tt*p2;
      }


    /*
      if (n < 0) return;
      values[0] = 1.0;
      if (n < 1) return;
      values[1] = x;
      if (n < 2) return;
      values[2] = p2 = 1.5 * x * x - 0.5 * tt;
      if (n < 3) return;
      values[3] = p1 =  (5.0/3.0) * p2 * x - (2.0/3.0) * tt * x;
      if (n < 4) return;

      for (int j=4; j < n; j+=2)
      {
      double invj = 1.0/j;
      p2 *= (invj-1) * tt;
      p2 += (2-invj) * x * p1;
      values[j]   = p2; 

      invj = 1.0/(j+1);
      p1 *= (invj-1) * tt;
      p1 += (2-invj) * x * p2;
      values[j+1] = p1; 
      }

      if (n % 2 == 0)
      {
      double invn = 1.0/n;
      values[n] = (2-invn)*x*p1 - (invn-1) * tt*p2;
      }
    */





    /*
      Sx p1 = 1.0, p2 = 0.0, p3;
  
      if (n>=0) values[0] = 1.0;

      for (int j=1; j<=n; j++)
      {
      p3=p2; p2=p1;
      p1=((2.0*j-1.0) * x*p2 - tt*(j-1.0)*p3)/j;
      values[j] = p1;
      }
    */
  }



  /* 
     E_i(x,y) = c * P_i(x/t) * t^i 
  */ 
  template <class Sx, class St, class Sc, class T>
  inline void ScaledLegendrePolynomialMult (int n, Sx x, St t, Sc c, T && values)
  {
    St tt = sqr(t);
    Sx p1, p2;

    if (n < 0) return;
    values[0] = p2 = c;
    if (n < 1) return;
    values[1] = p1 = c * x;
    if (n < 2) return;

    for (int j=2; j < n; j+=2)
      {
	StdOp (p2, tt, x*p1, 1.0/j-1);
	values[j] = p2;
	StdOp (p1, tt, x*p2, 1.0/(j+1)-1);
	values[j+1] = p1;
	/*
	  double invj = 1.0/j;
	  p2 *= (invj-1) * tt;
	  p2 += (2-invj) * x * p1;
	  values[j]   = p2; 

	  invj = 1.0/(j+1);
	  p1 *= (invj-1) * tt;
	  p1 += (2-invj) * x * p2;
	  values[j+1] = p1; 
	*/
      }

    if (n % 2 == 0)
      {
	StdOp (p2, tt, x*p1, 1.0/n-1);
	values[n] = p2;

	//      double invn = 1.0/n;
	//      values[n] = (2-invn)*x*p1 - (1-invn) * tt*p2;
      }
  }



















  template <class T> 
  inline void DiffScaledLegendrePolynomial (int n, double x, double t, T && values)
  {
    ArrayMem<AutoDiff<2>,10> ad_values(n+1);
    AutoDiff<2> ad_x(x, 0);
    AutoDiff<2> ad_t(t, 1);

    ScaledLegendrePolynomial(n, ad_x, ad_t, ad_values);

    for (int i = 0; i <= n; i++)
      for (int j = 0; j < 2; j++)
	values(i,j) = ad_values[i].DValue(j);
  }


  template <class Sx, class St, class T>
  inline void ScaledIntegratedLegendrePolynomial (int n, Sx x, St t, T & values)
  {
    Sx p1 = -1.0;
    Sx p2 = 0.0; 
    Sx p3;
    St tt = t*t;
    if (n >= 0)
      values[0] = -1.0;

    for (int j=1; j<=n; j++)
      {
	p3=p2; p2=p1;
	p1=((2.0*j-3.0) * x*p2 - tt*(j-3.0)*p3)/j;
	values[j] = p1;
      }
  }






  template <class T> 
  inline void ScaledLegendrePolynomialandDiff(int n, double x, double t, T & P, T & Px, T & Pt)  
  {
    /*
      ArrayMem<AutoDiff<2>,10> ad_values(n+1);
      AutoDiff<2> ad_x(x, 0);
      AutoDiff<2> ad_t(t, 1);

      ScaledLegendrePolynomial(n, ad_x, ad_t, ad_values);

      for (int i = 0; i <= n; i++)
      {
      P[i] = ad_values[i].Value();
      Px[i] = ad_values[i].DValue(0);
      Pt[i] = ad_values[i].DValue(1);
      }
    */
    if(n>=0) 
      {
	P[0] = 1.; 
	Px[0] = 0.; 
	Pt[0] = 0.; 
	if(n>=1) 
	  {
	    P[1] = x;
	    Px[1] = 1.; 
	    Pt[1] = 0.; 
	  } 

	double px0 = 0., px1 = 0., px2 =1.;
	double sqt = t*t; 
	for(int l = 2; l<=n; l++) 
	  { 
	    px0=px1; 
	    px1=px2;  
	    px2=  ( (2*l-1)*x*px1 - l*sqt*px0)/(l-1); 
	    Px[l] = px2; 
	    Pt[l] = -t*px1;
	    P[l] = (x*px2-sqt*px1)/l;
	  }
      }

  }
	  
  template <class T> 
  inline void LegendrePolynomialandDiff(int n, double x,  T & P, T & Px) 
  {
    /*
      ArrayMem<AutoDiff<1>,10> ad_values(n+1);
      AutoDiff<1> ad_x(x, 0);
      LegendrePolynomial(n, ad_x, ad_values);

      for (int i = 0; i <= n; i++)
      {
      P[i] = ad_values[i].Value();
      Px[i] = ad_values[i].DValue(0);
      }
 
      (*testout) << "P = " << endl << P << endl
      << "Px = " << endl << Px << endl;
    */
    if(n>=0) 
      {
	P[0] = 1.; 
	Px[0] = 0.;  
	if(n>=1) 
	  {
	    P[1] = x;
	    Px[1] = 1.;  
	  }       
	double px0 = 0., px1 = 0., px2 =1.;
	for(int l = 2; l<=n; l++) 
	  { 
	    px0=px1; 
	    px1=px2;  

	    px2=  ( (2*l-1)*x*px1 - l*px0)/(l-1); 
	    Px[l] = px2; 
	    P[l] = (x*px2 - px1)/l;
	  }
      }
  }




#ifdef OLD


  /*
    u(0) = 0, u(1) = 1,
    min \int_0^1 (1-x)^{DIM-1} (u')^2 dx

    representation as
    \sum c_i P_i(2x-1)
  */

  template <int DIM>
  class VertexExtensionOptimal
  {
    enum { SIZE = 50 };
    static double coefs[SIZE][SIZE];
    static bool initialized;
  public:

    VertexExtensionOptimal ();

    template <typename Tx>
    inline static Tx Calc (int p, Tx x)
    {
      Tx p1 = 1.0, p2 = 0.0, p3;
      Tx sum = 0;

      x = 2*x-1;

      if (p >= 0)
	sum += coefs[0][p];

      for (int j=1; j<=p; j++)
	{
	  p3 = p2; p2 = p1;
	  p1 = ((2.0*j-1.0)*x*p2 - (j-1.0)*p3) / j;
	  sum += coefs[j][p] * p1;
	}

      return sum;
    }

    inline static double CalcDeriv (int p, double x)
    {
      AutoDiff<1> p1 = 1.0, p2 = 0.0, p3;
      AutoDiff<1> sum = 0;
      AutoDiff<1> adx (x, 0);  // \nabla adx = e_0

      adx = 2.0*adx-1;

      if (p >= 0)
	sum += coefs[0][p];

      for (int j=1; j<=p; j++)
	{
	  p3 = p2; p2 = p1;
	  p1 = ((2.0*j-1.0)*adx*p2 - (j-1.0)*p3) / j;
	  sum += coefs[j][p] * p1;
	}

      return sum.DValue(0);
    }


    /*
    // Based on Jacobi pols, 3D only
    template <typename Tx>
    inline static Tx Calc (int p, Tx x)
    {
    ArrayMem<Tx,20> jacpol(p+1);

    JacobiPolynomial (p, 2*x-1, 1, -1, jacpol);
    
    Tx sum = 0;
    for (int j = 0; j <= p; j++)
    sum += coefs[j][p] * jacpol[j];

    return sum;
    }

    inline static double CalcDeriv (int p, double x)
    {
    ArrayMem<double,20> jacpol(p+1);

    JacobiPolynomial (p, 2*x-1, 2, 0, jacpol);
    
    double sum = 0;
    for (int j = 1; j <= p; j++)
    sum += coefs[j][p] * 0.5 * (j+1) * jacpol[j-1];

    return sum;
    }
    */
  };




  /*
    u(-1) = 0, u(1) = 1,
    min \int_-1^1 (1-x)^{DIM-1} (u')^2 dx
  */

  template <class Tx, class T>
  inline void LowEnergyVertexPolynomials2D  (int n, Tx x, T & values)
  {
    JacobiPolynomial (n, x, 0, -1, values);
    Tx sum1 = 0.0, sum2 = 0.0;
    for (int i = 1; i <= n; i++)
      {
	sum1 += 1.0/i;
	sum2 += values[i] / i;
	values[i] = sum2/sum1;
      }
    values[0] = 1;
  }

  template <class Tx, class T>
  inline void LowEnergyVertexPolynomials3D  (int n, Tx x, T & values)
  {
    JacobiPolynomial (n, x, 1, -1, values);
    Tx sum = 0.0;
    for (int i = 1; i <= n; i++)
      {
	sum += (2.0*i+1)/(i+1) * values[i];
	values[i] = 1.0/(i*(i+2)) * sum;
      }
    values[0] = 1;
  }





  class VertexStandard
  {
  public:
    template <typename Tx>
    inline static Tx Calc (int p, Tx x)
    {
      return x;
    }

    inline static double CalcDeriv (int p, double x)
    {
      return 1;
    }
  };

#endif


  /**
     Compute triangle edge-shape functions

     functions vanish on upper two edges

     x,y: coordinates in triangle (-1, 0), (1, 0), (0, 1)

     f_i-2 (x, 0) = IntegratedLegendrePol_i (x)

     f_i-2 ... pol of order i, max order = n

     Monomial extension:
  */

  class IntegratedLegendreMonomialExt
  {
    enum { SIZE = 1000 };
    static double coefs[SIZE][2];
  
  public:

    static void CalcCoeffs ()
    {
      for (int j = 1; j < SIZE; j++)
	{
	  coefs[j][0] = double(2*j-3)/j;
	  coefs[j][1] = double(j-3)/j;
	}
    }

    template <class Sx, class Sy, class T>
    inline static int CalcScaled (int n, Sx x, Sy y, T & values)
    {
      Sy fy = y*y;
      Sx p3 = 0;
      Sx p2 = -1;
      Sx p1 = x;

      for (int j=2; j<=n; j++)
	{
	  p3=p2; p2=p1;
	  p1 = double(2*j-3)/j  * x * p2 - double(j-3)/j * fy * p3;
	  values[j-2] = p1;
	}     

      return n-1;
    }


    template <int n, class Sx, class Sy, class T>
    inline static int CalcScaled (Sx x, Sy y, T & values)
    {
      Sy fy = y*y;
      Sx p3 = 0;
      Sx p2 = -1;
      Sx p1 = x;

      for (int j=2; j<=n; j++)
	{
	  p3=p2; p2=p1;
	  p1 = double(2*j-3)/j  * x * p2 - double(j-3)/j * fy * p3;
	  values[j-2] = p1;
	}     

      return n-1;
    }





    template <class Sx, class Sy, class T>
    inline static int CalcTrigExt (int n, Sx x, Sy y, T & values)
    {
      Sy fy = (1-y)*(1-y);
      Sx p3 = 0;
      Sx p2 = -1;
      Sx p1 = x;

      for (int j=2; j<=n; j++)
	{
	  p3=p2; p2=p1;
	  p1 = double(2*j-3)/j  * x * p2 - double(j-3)/j * fy * p3;
	  values[j-2] = p1;
	}     

      return n-1;
    }

    template <class Sx, class Sy, class Sf, class T>
    inline static int CalcTrigExtMult (int n, Sx x, Sy y, Sf fac, T & values)
    {
      Sy fy = (1-y)*(1-y);
      Sx p3 = 0;
      Sx p2 = -fac;
      Sx p1 = x * fac;

      for (int j=2; j<=n; j++)
	{
	  p3=p2; p2=p1;
	  // p1=( (2*j-3) * x * p2 - (j-3) * fy * p3) / j;
	  p1 = double(2*j-3)/j  * x * p2 - double(j-3)/j * fy * p3;
	  // p1= coefs[j][0] * x * p2 - coefs[j][1] * fy * p3;
	  values[j-2] = p1;
	}     

      return n-1;
    }




    template <class T>
    inline static int CalcTrigExtDeriv (int n, double x, double y, T & values)
    {
      double fy = (1-y)*(1-y);
      double p3 = 0, p3x = 0, p3y = 0;
      double p2 = -1, p2x = 0, p2y = 0;
      double p1 = x, p1x = 1, p1y = 0;

      for (int j=2; j<=n; j++)
	{
	  p3=p2; p3x = p2x; p3y = p2y;
	  p2=p1; p2x = p1x; p2y = p1y;
	  double c1 = (2.0*j-3) / j;
	  double c2 = (j-3.0) / j;
	
	  p1  = c1 * x * p2 - c2 * fy * p3;
	  p1x = c1 * p2 + c1 * x * p2x - c2 * fy * p3x;
	  p1y = c1 * x * p2y - (c2 * 2 * (y-1) * p3 + c2 * fy * p3y);
	  values (j-2, 0) = p1x;
	  values (j-2, 1) = p1y;
	}    
      return n-1;
    }


    template <class Sx, class T>
    inline static int Calc (int n, Sx x, T & values)
    {
      Sx p3 = 0;
      Sx p2 = -1;
      Sx p1 = x;

      for (int j=2; j<=n; j++)
	{
	  p3=p2; p2=p1;
	  p1=( (2*j-3) * x * p2 - (j-3) * p3) / j;
	  values[j-2] = p1;
	}
      return n-1;
    }

    template <class Sx, class Sf, class T>
    inline static int CalcMult (int n, Sx x, Sf fac, T & values)
    {
      Sx p3 = 0;
      Sx p2 = -fac;
      Sx p1 = x*fac;

      for (int j=2; j<=n; j++)
	{
	  p3=p2; p2=p1;
	  p1=( (2*j-3) * x * p2 - (j-3) * p3) / j;
	  values[j-2] = p1;
	}
      return n-1;
    }



    template <class T>
    inline static int CalcDeriv (int n, double x, T & values)
    {
      double p1 = 1.0, p2 = 0.0, p3;

      for (int j=1; j<=n-1; j++)
	{
	  p3 = p2; p2 = p1;
	  p1 = ((2.0*j-1.0)*x*p2 - (j-1.0)*p3) / j;
	  values[j-1] = p1;
	}
      return n-1;
    }


  };






  /*   Conversion of orthogonal polynomials */

  // given: \sum u_i P^{al,0}
  // find:  \sum v_i P^{al-1, 0}
  template <class T>
  void ConvertJacobiReduceAlpha (int n, int alpha, T & inout)
  {
    for (int i = n; i >= 1; i--)
      {
	double val = inout(i) / (i+alpha);
	inout(i) = (2*i+alpha) * val;
	inout(i-1) += i * val;
      }
  }

  // given: \sum u_i (1-x) P^{al+1,0}   0 <= i <  n
  // find:  \sum v_i P^{al, 0}          0 <= i <= n
  template <class T>
  void ConvertJacobiReduceAlphaFactor (int n, double alpha, T & inout) 
  {
    inout(n) = 0;
    for (int i = n; i > 0; i--)
      {
	double val = inout(i-1) / (i+alpha/2);
	inout(i-1) = (i+alpha) * val;
	inout(i) -= i * val;
      }
  }





  // Differentiate Jacobi
  // (P_i^alpha)' 

  template <class T>
  void DifferentiateJacobi (int n, double alpha, T & inout)
  {
    for (int i = n; i >= 1; i--)
      {
	double val = inout(i);
	inout(i-1) -= alpha*(2*i+alpha-1) / ( (i+alpha)*(2*i+alpha-2)) * val;
	if (i > 1)
	  inout(i-2) += (i-1)*(2*i+alpha) / ( (i+alpha)*(2*i+alpha-2))  * val;

	inout(i) *= (2*i+alpha)*(2*i+alpha-1) / ( 2 * (i+alpha) );
      }
    for (int i = 0; i < n; i++)
      inout(i) = inout(i+1);
    inout(n) = 0;
  }


  template <class T>
  void DifferentiateJacobiTrans (int n, double alpha, T & inout)
  {
    for (int i = n-1; i >= 0; i--)
      inout(i+1) = inout(i);
    inout(0) = 0;

    for (int i = 1; i <= n; i++)
      {
	inout(i) *= (2*i+alpha)*(2*i+alpha-1) / ( 2 * (i+alpha) );

	inout(i) -= alpha*(2*i+alpha-1) / ( (i+alpha)*(2*i+alpha-2)) * inout(i-1);
	if (i > 1)
	  inout(i) += (i-1)*(2*i+alpha) / ( (i+alpha)*(2*i+alpha-2))  * inout(i-2);
      }
  }


  template <class T>
  void DifferentiateLegendre (int n, T & inout)
  {
    for (int i = n; i >= 1; i--)
      {
	if (i > 1) inout(i-2) += inout(i);
	inout(i) *= 2*i-1;
      }
    for (int i = 0; i < n; i++)
      inout(i) = inout(i+1);
    inout(n) = 0;
  }

  template <class T>
  void DifferentiateLegendreTrans (int n, T & inout)
  {
    for (int i = n-1; i >= 0; i--)
      inout(i+1) = inout(i);
    inout(0) = 0;

    for (int i = 1; i <= n; i++)
      {
	inout(i) *= (2*i-1);
	if (i > 1)
	  inout(i) += inout(i-2);
      }
  }






  class ConvertJacobi
  {
    typedef double d2[2];
    NGS_DLL_HEADER static Array<d2*> coefs_increasealpha;
    NGS_DLL_HEADER static Array<d2*> coefs_reducealpha;
    NGS_DLL_HEADER static Array<d2*> coefs_reducealphafac;

  public:
    ConvertJacobi ();
    ~ConvertJacobi ();

    template <class T>
    static void ReduceAlpha (int n, int alpha, T & inout)  // alpha of input
    {
      d2 * c = coefs_reducealpha[alpha];

      double val = inout[n];
      for (int i = n; i >= 1; i--)
	{
	  inout[i] = c[i][1] * val;
	  val = c[i][0] * val + inout[i-1];
	}
      inout[0] = val;

      /*
	for (int i = n; i >= 1; i--)
	{
	inout[i-1] += c[i][0] * inout[i];
	inout[i] *= c[i][1];
	}    
      */
    }

    template <class T>
    static void ReduceAlphaFactor (int n, int alpha, T & inout) // alpha of output
    {
      d2 * c = coefs_reducealphafac[alpha];

      inout(n) = c[n][0] * inout[n-1];
      for (int i = n-1; i >= 1; i--)
	inout[i] = c[i+1][1] * inout[i] + c[i][0] * inout[i-1];
      inout[0] = c[1][1] * inout[0];

      /*
	for (int i = n; i >= 1; i--)
	{
	inout[i] += c[i][0] * inout[i-1];
	inout[i-1] *= c[i][1];
	}    
      */
    }


    template <class T>
    static void ReduceAlphaTrans (int n, int alpha, T & inout)  // alpha of input
    {
      d2 * c = coefs_reducealpha[alpha];

      /*
	for (int i = n; i >= 1; i--)
	{
	inout[i-1] += c[i][0] * inout[i];
	inout[i] *= c[i][1];
	}    
      */
      for (int i = 1; i <= n; i++)
	{
	  inout[i] *= c[i][1];
	  inout[i] += c[i][0] * inout[i-1];
	}    

    }

    template <class T>
    static void ReduceAlphaFactorTrans (int n, int alpha, T & inout) // alpha of output
    {
      d2 * c = coefs_reducealphafac[alpha];
      /*
	for (int i = n; i >= 1; i--)
	{
	inout[i] += c[i][0] * inout[i-1];
	inout[i-1] *= c[i][1];
	}    
      */
      for (int i = 1; i <= n; i++)
	{
	  inout[i-1] *= c[i][1];
	  inout[i-1] += c[i][0] * inout[i];
	}    
      inout[n] = 0;
    }

















    // reduce, fac
    // P_i^alpha(x) (1-x)/2 = c_i^{alpha} P_i^{\alpha-1} + hatc_i^{alpha} P_{i+1}^{\alpha-1}
    static double c (int i, int alpha) { return double(i+alpha)/double(2*i+alpha+1); }
    static double hatc (int i, int alpha) { return -double(i+1)/double(2*i+alpha+1); }

    // increase alpha
    // P_i^alpha(x)  = d_i^{alpha} P_i^{\alpha+1} + hatd_i^{alpha} P_{i-1}^{\alpha+1}
    static double d (int i, int alpha) { return double(i+alpha+1)/double(2*i+alpha+1); }
    static double hatd (int i, int alpha) { return -double(i)/double(2*i+alpha+1); }

    // decrease alpha
    // P_i^alpha(x)  = e_i^{alpha} P_i^{\alpha-1} + hate_i^{alpha} P_{i-1}^{\alpha}
    static double e (int i, int alpha) { return double(2*i+alpha)/double(i+alpha); }
    static double hate (int i, int alpha) { return double(i)/double(i+alpha); }

  
    static Array<d2*> coefs_c, coefs_d, coefs_e;

    // alpha,beta of input,
    // order of input
    // reduce alpha, reduce beta, reduce factor (1-x)/2 (1-y)/2
    template <class T>
    static void TriangularReduceFactor (int order, int alpha, int beta, T & inout)
    {
      for (int i = 0; i <= order+1; i++)
	inout(i, order+1-i) = 0;

      for (int j = order; j >= 0; j--)
	{
	  d2 * calpha = coefs_c[alpha+2*j];
	  d2 * cbeta = coefs_c[beta];
	  d2 * dalpha = coefs_d[alpha+2*j];

	  for (int i = order-j; i >= 0; i--)
	    {
	      double val = inout(i,j);
	      inout(i,j)    = val * calpha[i][0] * cbeta[j][0];
	      inout(i+1,j) += val * calpha[i][1] * cbeta[j][0];
	      inout(i,j+1) += val * dalpha[i][0] * cbeta[j][1];
	      if (i > 0)
		inout(i-1,j+1) += val * dalpha[i][1] * cbeta[j][1];
	    }

	  /*
	    for (int i = order-j; i >= 0; i--)
	    {
	    double val = inout(i,j);
	    inout(i,j)   = val * c(i,alpha+2*j) * c(j,beta);
	    inout(i+1,j) += val * hatc(i, alpha+2*j) * c(j, beta);
	    inout(i,j+1) += val * d(i, alpha+2*j) * hatc (j,beta);
	    if (i > 0)
	    inout(i-1,j+1) += val * hatd(i, alpha+2*j) * hatc(j,beta);
	    }
	  */
	}
    }

    // alpha,beta of input,
    // const alpha, dec beta
    // order is constant
    template <class T, class S>
    static void TriangularReduceBeta (int order, int alpha, int beta, T & inout, S & hv)
    {
      d2 * ebeta  = coefs_e[beta];

      for (int j = order; j > 0; j--)
	{
	  d2 * calpha = coefs_c[alpha+2*j];

	  hv = 0.0;
	  for (int i = 0; i <= order-j; i++)
	    {
	      double val = inout(i,j);
	      inout(i,j)  = val * ebeta[j][0]; 

	      hv(i)   += val * calpha[i][0] * ebeta[j][1];
	      hv(i+1) += val * calpha[i][1] * ebeta[j][1];
	    }

	  ReduceAlpha (order-j+1, alpha+2*j-1, hv);
	  inout.Col(j-1) += hv;
	}
    }  
  };












  // template meta-programming

  template <int N, int AL>
  class TReduceAlpha
  { 
  public:
    enum { FLOP = TReduceAlpha<N-1,AL>::FLOP + 2 };

    template <class T>
    static ALWAYS_INLINE void Do (T & inout)
    {
      inout[N-1] += double(N)/double(N+AL) * inout[N];
      inout[N] *= double(2*N+AL)/double(N+AL);
      TReduceAlpha<N-1,AL>::Do(inout);
    }

    template <class T>
    static ALWAYS_INLINE void Trans (T & inout)
    {
      TReduceAlpha<N-1,AL>::Trans(inout);
      inout[N] *= double(2*N+AL)/double(N+AL);
      inout[N] += double(N)/double(N+AL) * inout[N-1];
    }  
  };

  template <int AL>
  class TReduceAlpha<0,AL> 
  { 
  public:
    enum { FLOP = 0 };
    template <class T>
    static void Do (T & inout) { ; }

    template <class T>
    static void Trans (T & inout) { ; }
  };





  template <int N, int AL, int HIGHEST = 1>
  class TReduceAlphaFactor
  { 
  public:
    enum { FLOP = TReduceAlphaFactor<N-1,AL>::FLOP + 2 };

    template <class T>
    static  ALWAYS_INLINE void Do (T & inout)
    {
      if (HIGHEST)
	inout[N] = double(-N)/double(2*N+AL) * inout[N-1];
      else
	inout[N] += double(-N)/double(2*N+AL) * inout[N-1];
      inout[N-1] *= double(N+AL)/double(2*N+AL);
      TReduceAlphaFactor<N-1,AL,0>::Do(inout);
    }

    template <class T>
    static  ALWAYS_INLINE void Trans (T & inout)
    {
      TReduceAlphaFactor<N-1,AL,0>::Trans(inout);
      inout[N-1] *= double(N+AL)/double(2*N+AL);
      inout[N-1] += double(-N)/double(2*N+AL) * inout[N];
    }
  };

  template <int AL, int HIGHEST>
  class TReduceAlphaFactor<0,AL,HIGHEST> 
  { 
  public:
    enum { FLOP = 0 };

    template <class T>
    static inline void Do (T & inout)
    { 
      if (HIGHEST) inout[0] = 0.0;
    }

    template <class T>
    static inline void Trans (T & inout) { ; }
  };






  /*
    template <class T>
    void DifferentiateJacobi (int n, double alpha, T & inout)
    {
    for (int i = n; i >= 1; i--)
    {
    double val = inout(i);
    inout(i-1) -= alpha*(2*i+alpha-1) / ( (i+alpha)*(2*i+alpha-2)) * val;
    if (i > 1)
    inout(i-2) += (i-1)*(2*i+alpha) / ( (i+alpha)*(2*i+alpha-2))  * val;

    inout(i) *= (2*i+alpha)*(2*i+alpha-1) / ( 2 * (i+alpha) );
    }
    for (int i = 0; i < n; i++)
    inout(i) = inout(i+1);
    inout(n) = 0;
    }
  */



  template <int N, int AL>
  class TDifferentiateJacobi
  { 
  public:
    enum { FLOP = TDifferentiateJacobi<N-1,AL>::FLOP + 3 };

    // (P_i^Al)' = c1 P_i^AL + c2 (P_{i-1}^AL)' + c3 (P_{i-2}^AL)' 

    static double c1() { return double ( (2*N+AL)*(2*N+AL-1) ) / double( 2 * (N+AL) ); }
    static double c2() { return -double (AL*(2*N+AL-1)) / double( (N+AL)*(2*N+AL-2)); }
    static double c3() { return double((N-1)*(2*N+AL)) / double( (N+AL)*(2*N+AL-2)); }
    /*
      enum { c2 = -double (AL*(2*N+AL-1)) / double( (N+AL)*(2*N+AL-2)) };
      enum { c3 = double((N-1)*(2*N+AL)) / double( (N+AL)*(2*N+AL-2)) };
    */
    template <class T>
    static  ALWAYS_INLINE void Do (T & inout)
    {
      double val = inout(N);
      inout(N-1) += c2() * val;
      if (N > 1)
	inout(N-2) += c3() * val;

      inout(N) *= c1();
      /*
	double val = inout(N);
	inout(N-1) -= double (AL*(2*N+AL-1)) / double( (N+AL)*(2*N+AL-2)) * val;
	if (N > 1)
	inout(N-2) += double((N-1)*(2*N+AL)) / double( (N+AL)*(2*N+AL-2))  * val;

	inout(N) *= double ( (2*N+AL)*(2*N+AL-1) ) / double( 2 * (N+AL) );
      */
      TDifferentiateJacobi<N-1,AL>::Do(inout);

      inout(N-1) = inout(N);
      inout(N) = 0;
    }

    template <class T>
    static  ALWAYS_INLINE void Trans (T & inout)
    {
      inout(N) = inout(N-1);
      inout(N-1) = 0;

      TDifferentiateJacobi<N-1,AL>::Trans(inout);

      inout(N) *= double ((2*N+AL)*(2*N+AL-1)) / ( 2 * (N+AL) );

      inout(N) -= double (AL*(2*N+AL-1)) / ( (N+AL)*(2*N+AL-2)) * inout(N-1);
      if (N > 1)
	inout(N) += double ((N-1)*(2*N+AL)) / ( (N+AL)*(2*N+AL-2))  * inout(N-2);

      /*      
	      double val = inout(N);
	      inout(N-1) -= double (AL*(2*N+AL-1)) / double( (N+AL)*(2*N+AL-2)) * val;
	      if (N > 1)
	      inout(N-2) += double((N-1)*(2*N+AL)) / double( (N+AL)*(2*N+AL-2))  * val;

	      inout(N) *= double ( (2*N+AL)*(2*N+AL-1) ) / double( 2 * (N+AL) );
      */

    }


  };

  template <int AL>
  class TDifferentiateJacobi<0,AL> 
  { 
  public:
    enum { FLOP = 0 };

    template <class T>
    static inline void Do (T & inout) { ; }

    template <class T>
    static inline void Trans (T & inout) { ; }
  };









  template <int N>
  class TDifferentiateLegendre
  { 
  public:
    enum { FLOP = TDifferentiateLegendre<N-1>::FLOP + 3 };

    // (P_i^Al)' = c1 P_i^AL + c2 (P_{i-1}^AL)' + c3 (P_{i-2}^AL)' 

    template <class T>
    static  ALWAYS_INLINE void Do (T & inout)
    {
      double val = inout(N);
      if (N > 1) inout(N-2) += val;
    
      inout(N) *= double (2*N-1);
      TDifferentiateLegendre<N-1>::Do(inout);

      inout(N-1) = inout(N);
      inout(N) = 0;
    }

    template <class T>
    static  ALWAYS_INLINE void Trans (T & inout)
    {
      inout(N) = inout(N-1);
      inout(N-1) = 0;

      TDifferentiateLegendre<N-1>::Trans(inout);

      inout(N) *= double (2*N-1);
      if (N > 1) inout(N) += inout(N-2);
    }

  };

  template <>
  class TDifferentiateLegendre<0> 
  { 
  public:
    enum { FLOP = 0 };

    template <class T>
    static inline void Do (T & inout) { ; }

    template <class T>
    static inline void Trans (T & inout) { ; }
  };




















  template <int N, int J, int I, int AL, int BE>
  class TTriangleReduceFactorCol
  {
    // reduce, fac
    // P_i^alpha(x) (1-x)/2 = c_i^{alpha} P_i^{\alpha-1} + hatc_i^{alpha} P_{i+1}^{\alpha-1}
    static double c (int i, int alpha) { return double(i+alpha)/double(2*i+alpha+1); }
    static double hatc (int i, int alpha) { return -double(i+1)/double(2*i+alpha+1); }

    // increase alpha
    // P_i^alpha(x)  = d_i^{alpha} P_i^{\alpha+1} + hatd_i^{alpha} P_{i-1}^{\alpha+1}
    static double d (int i, int alpha) { return double(i+alpha+1)/double(2*i+alpha+1); }
    static double hatd (int i, int alpha) { return -double(i)/double(2*i+alpha+1); }

  public:
    enum { FLOP = TTriangleReduceFactorCol<N,J,I-1,AL,BE>::FLOP + 4 };

    template <class T>
    static  ALWAYS_INLINE void Do (T & inout)
    {
      double val    = inout(I,J);

      inout(I,J)    = val * c(I,AL+2*J) * c(J,BE);
      inout(I+1,J) += val * hatc(I, AL+2*J) * c(J, BE);
      inout(I,J+1) += val * d(I, AL+2*J) * hatc (J, BE);
      if (I > 0)
	inout(I-1,J+1) += val * hatd(I, AL+2*J) * hatc(J,BE);

      TTriangleReduceFactorCol<N,J,I-1,AL,BE> ::Do(inout);
    }


    template <class T>
    static  ALWAYS_INLINE void Trans (T & inout)
    {
      TTriangleReduceFactorCol<N,J,I-1,AL,BE> ::Trans(inout);

      double val = 
	c(I,AL+2*J) * c(J,BE) * inout(I,J)
	+ hatc(I, AL+2*J) * c(J, BE) * inout(I+1,J)
	+ d(I, AL+2*J) * hatc (J, BE) * inout(I,J+1);

      if (I > 0)
	val += hatd(I, AL+2*J) * hatc(J,BE) * inout(I-1,J+1);

      inout(I,J) = val;
    }

  };


  template <int N, int J, int AL, int BE>
  class TTriangleReduceFactorCol<N,J,-1,AL,BE>
  {
  public:
    enum { FLOP = 0 };
    template <class T>
    static void Do (T & inout) { ; }
    template <class T>
    static void Trans (T & inout) { ; }
  };


  template <int N, int J, int AL, int BE>
  class TTriangleReduceFactor
  {
  public:
    enum { FLOP = TTriangleReduceFactor<N,J-1,AL,BE>::FLOP + 
	   TTriangleReduceFactorCol<N,J,N-J,AL,BE>::FLOP };

    template <class T>
    static  ALWAYS_INLINE void Do (T & inout)
    {
      if (J == N)
	for (int i = 0; i <= N+1; i++)
	  inout(i, N+1-i) = 0;
      
      TTriangleReduceFactorCol<N,J,N-J,AL,BE> ::Do(inout);
      TTriangleReduceFactor<N,J-1,AL,BE>  ::Do(inout);
    }

    template <class T>
    static  ALWAYS_INLINE void Trans (T & inout)
    {
      /*
	if (J == N)
	for (int i = 0; i <= N+1; i++)
	inout(i, N+1-i) = 0;
      */
      TTriangleReduceFactor<N,J-1,AL,BE>  ::Trans(inout);
      TTriangleReduceFactorCol<N,J,N-J,AL,BE> ::Trans(inout);
    }
  };


  template <int N, int AL, int BE>
  class TTriangleReduceFactor<N,-1,AL,BE>
  {
  public:
    enum { FLOP = 0 };
    template <class T>
    static void Do (T & inout) { ; }
    template <class T>
    static void Trans (T & inout) { ; }
  };

 






  /*


  template <int N, int J, int AL, int BE>
  class TTriangleReduce
  {
  // reduce, fac
  // P_i^alpha(x) (1-x)/2 = c_i^{alpha} P_i^{\alpha-1} + hatc_i^{alpha} P_{i+1}^{\alpha-1}
  static double c (int i, int alpha) { return double(i+alpha)/double(2*i+alpha+1); }
  static double hatc (int i, int alpha) { return -double(i+1)/double(2*i+alpha+1); }

  // decrease alpha
  // P_i^alpha(x)  = e_i^{alpha} P_i^{\alpha-1} + hate_i^{alpha} P_{i-1}^{\alpha}
  static double e (int i, int alpha) { return double(2*i+alpha)/double(i+alpha); }
  static double hate (int i, int alpha) { return double(i)/double(i+alpha); }


  public:
  enum { FLOP = TTriangleReduce<N,J-1,AL,BE>::FLOP + 3*(N-J+1) +
  TReduceAlpha<N-J+1, AL+2*J-1>::FLOP
  };


  template <class T>
  static  ALWAYS_INLINE void Do (T & inout)
  {
  Vec<N-J+2> hv;
  hv(0) = 0.0;

  for (int i = 0; i <= N-J; i++)
  {
  double val = inout(i,J);
  inout(i,J)  = val * e(J,BE);
  hv(i)   += val * c(i, AL+2*J) * hate(J, BE);
  hv(i+1)  = val * hatc(i, AL+2*J) * hate(J, BE);
  }
      
  TReduceAlpha<N-J+1, AL+2*J-1>::Do(hv);
      
  for (int i = 0; i <= N-J+1; i++)
  inout(i,J-1) += hv(i);

  TTriangleReduce<N,J-1,AL,BE>  ::Do(inout);

  TReduceAlpha<N-J, AL+2*J>::Do(inout.Col(J));
  }




  template <class T>
  static  ALWAYS_INLINE void Trans (T & inout)
  {
  Vec<N-J+2> hv;
    
  TReduceAlpha<N-J, AL+2*J>::Trans(inout.Col(J));

  TTriangleReduce<N,J-1,AL,BE>::Trans(inout);

  for (int i = 0; i <= N-J+1; i++)
  hv(i) = inout(i,J-1);

  TReduceAlpha<N-J+1, AL+2*J-1>::Trans(hv);    

  for (int i = 0; i <= N-J; i++)
  {
  inout(i,J) = 
  e(J,BE) * inout(i,J) 
  + c(i, AL+2*J) * hate(J, BE) * hv(i)
  + hatc(i, AL+2*J) * hate(J, BE) * hv(i+1);
  } 
  }
  };


  template <int N, int AL, int BE>
  class TTriangleReduce<N,0,AL,BE>
  {
  public:
  enum { FLOP = 0 };
  template <class T>
  static void Do (T & inout) 
  { 
  TReduceAlpha<N, AL>::Do(inout.Col(0));
  }
  
  template <class T>
  static void Trans (T & inout) 
  {
  TReduceAlpha<N, AL>::Trans(inout.Col(0));
  }
  };
  */








  template <int N, int I, int J, int AL, int BE>
  class TTriangleReduceLoop2New
  {
  public:
    template <class T>
    static  ALWAYS_INLINE void Do (T & inout)
    {
      if (BE+I == 0) cout << "is 0" << endl;
      double fac = 1.0 / ( (BE + I)*(4 + AL-1 + 2*I-2 + J-1) );
      double val = inout(J,I);
      inout(J,I) = fac * (BE + 2*I)*(5 + AL-1 + 2*I-2 + 2*J-2) * val;
      if (I >= 1)
	{
	  inout(J,I-1) += fac * I*(3 + AL-1 + 2*I-2 + J-1) * val;
	  inout(1+J,I-1) -= fac * I*(2 + J-1) * val;
	}
      if (J >= 1)
	inout(J-1, I) += fac * (BE + I)*(1 + J-1) * val;

      TTriangleReduceLoop2New<N,I,J-1,AL,BE>::Do(inout);
    }


    template <class T>
    static  ALWAYS_INLINE void Trans (T & inout)
    {
      TTriangleReduceLoop2New<N,I,J-1,AL,BE>::Trans(inout);

      double fac = 1.0 / ( (BE + I)*(4 + AL-1 + 2*I-2 + J-1) );
      double val = inout(J,I) * ( fac * (BE + 2*I)*(5 + AL-1 + 2*I-2 + 2*J-2) );
      if (I >= 1)
	{
	  val += inout(J,I-1) * ( fac * I*(3 + AL-1 + 2*I-2 + J-1) );
	  val -= inout(1+J,I-1) * ( fac * I*(2 + J-1));
	}
      if (J >= 1)
	val += inout(J-1, I) * (fac * (BE + I)*(1 + J-1));
      
      inout(J,I) = val;
    }

  };

  template <int N, int AL, int I, int BE>
  class TTriangleReduceLoop2New<N,I,-1,AL,BE>
  {
  public:
    enum { FLOP = 0 };
    template <class T>
    static void Do (T & inout) { ; }
  
    template <class T>
    static void Trans (T & inout) { ; }
  };



  template <int N, int I, int AL, int BE>
  class TTriangleReduceNew
  {
  public:
    enum { FLOP = TTriangleReduceNew<N,I-1,AL,BE>::FLOP + 4*(N-I+1) };


    template <class T>
    static  ALWAYS_INLINE void Do (T & inout)
    {
      TTriangleReduceLoop2New<N,I,N-I,AL,BE>::Do(inout);
      TTriangleReduceNew<N,I-1,AL,BE>::Do(inout);
    }

    template <class T>
    static  ALWAYS_INLINE void Trans (T & inout)
    {
      TTriangleReduceNew<N,I-1,AL,BE>::Trans(inout);
      TTriangleReduceLoop2New<N,I,N-I,AL,BE>::Trans(inout);
    }
  };

  template <int N, int AL, int BE>
  class TTriangleReduceNew<N,-1,AL,BE>
  {
  public:
    enum { FLOP = 0 };

    template <class T>
    static void Do (T & inout) { ; }
  
    template <class T>
    static void Trans (T & inout) { ; }
  };





  /*

  template <int N, int I, int AL, int BE>
  class TTriangleReduceNew
  {
  public:

  template <class T>
  static  ALWAYS_INLINE void Do (T & inout)
  {
  for (int J = N-I; J >= 0; J--)
  {
  double fac = 1.0 / ( (2 + BE-1 + I-1)*(4 + AL-1 + 2*I-2 + J-1) );
  double val = inout(J,I);
  inout(J,I) = fac * (3 + BE-1 + 2*I-2)*(5 + AL-1 + 2*I-2 + 2*J-2) * val;
  if (I >= 1)
  {
  inout(1+J,I-1) += -fac * I*(2 + J-1) * val;
  inout(J,I-1) += fac * I*(3 + AL-1 + 2*I-2 + J-1) * val;
  }
  if (J >= 1)
  inout(J-1, I) += fac * (2 + BE-1 + I-1)*(1 + J-1) * val;
  }

  TTriangleReduceNew<N,I-1,AL,BE>::Do(inout);
  }


  template <class T>
  static  ALWAYS_INLINE void Trans (T & inout)
  {
  TTriangleReduceNew<N,I-1,AL,BE>::Trans(inout);

  for (int J = 0; J <= N-I; J++)
  {
  double fac = 1.0 / ( (2 + BE-1 + I-1)*(4 + AL-1 + 2*I-2 + J-1) );
  double val = inout(J,I) * ( fac * (3 + BE-1 + 2*I-2)*(5 + AL-1 + 2*I-2 + 2*J-2) );
  if (I >= 1)
  {
  val += inout(J,I-1) * ( fac * I*(3 + AL-1 + 2*I-2 + J-1) );
  val -= inout(1+J,I-1) * ( fac * I*(2 + J-1));
  }
          
  if (J >= 1)
  val += inout(J-1, I) * (fac * (2 + BE-1 + I-1)*(1 + J-1));
          
  inout(J,I) = val;
  }
  }

  };

  template <int N, int AL, int BE>
  class TTriangleReduceNew<N,-1,AL,BE>
  {
  public:
  enum { FLOP = 0 };
  template <class T>
  static void Do (T & inout) 
  { ; }
  
  template <class T>
  static void Trans (T & inout) 
  { ; }
  };
  */












}


#endif
