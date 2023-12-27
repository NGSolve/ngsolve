#ifndef FILE_RECURSIVE_POL
#define FILE_RECURSIVE_POL

/*********************************************************************/
/* File:   recursive_pol.hpp                                         */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#include <bla.hpp>


namespace ngfem
{
  using namespace ngbla;

  
  /*
    Recursive Polynomials
  */


  // Use Lambda function with square-bracket assignment
  template <typename TI, typename FUNC>
  class SBLambdaElement
  {
    FUNC f;
    TI i;
  public:
    INLINE SBLambdaElement (const SBLambdaElement & e2) = default; 
    INLINE SBLambdaElement (FUNC af, TI hi) : f(af), i(hi) { ; }
    template <typename VAL>
    INLINE VAL operator= (VAL v) { f(i, v); return v; }
  };

  template <typename TI, typename FUNC>
  class Class_SBLambda
  {
    FUNC func;
    TI offset;
  public:
    INLINE Class_SBLambda (const Class_SBLambda & l2) = default; 
    INLINE Class_SBLambda (FUNC f, TI ao) : func(f), offset(ao) { ; }

    template <typename TI2>    
    INLINE auto operator[] (TI2 i) const -> SBLambdaElement<decltype(offset+i),FUNC> 
    { 
      return SBLambdaElement<decltype(offset+i),FUNC> (func, offset+i); 
    }

    template <typename TI2>
    INLINE auto operator+ (TI2 i) const -> Class_SBLambda<decltype(offset+i),FUNC>
    { return Class_SBLambda<decltype(offset+i),FUNC> (func, offset+i); }

    template <typename TI2>
    INLINE auto Addr (TI2 i) const -> Class_SBLambda<decltype(offset+i),FUNC>
    {return Class_SBLambda<decltype(offset+i),FUNC> (func, offset+i); }
  };

  template <typename FUNC> 
  INLINE const Class_SBLambda<IC<0>,FUNC> SBLambda (FUNC f)
  {
    return Class_SBLambda<IC<0>,FUNC> (f, IC<0>());
  }




  template <typename FUNC, typename FUNC2>
  class Class_SBLambdaDuo
  {
    FUNC func;
    FUNC2 func2;
    int offset;
  public:
    INLINE Class_SBLambdaDuo (const Class_SBLambdaDuo & l2) : func(l2.func), func2(l2.func2), offset(l2.offset) { ; }
    INLINE Class_SBLambdaDuo (FUNC f, FUNC2 f2, int ao = 0) : func(f), func2(f2), offset(ao) { ; }
    INLINE SBLambdaElement<int,FUNC> operator[] (int i) const { return SBLambdaElement<int,FUNC> (func, offset+i); }
    INLINE Class_SBLambdaDuo<FUNC, FUNC2> operator+ (int i) const { return Class_SBLambdaDuo<FUNC,FUNC2> (func, func2, offset+i); }
    INLINE Class_SBLambdaDuo<FUNC, FUNC2> Addr (int i) const { return Class_SBLambdaDuo<FUNC,FUNC2> (func, func2, offset+i); }
    template <typename VAL1, typename VAL2>
    void operator() (int i1, VAL1 v1, int i2, VAL2 v2) const { func2(i1, v1, i2, v2); }
  };

  template <typename FUNC, typename FUNC2> 
  INLINE const Class_SBLambdaDuo<FUNC, FUNC2> SBLambdaDuo (FUNC f, FUNC2 f2)
  {
    return Class_SBLambdaDuo<FUNC,FUNC2> (f, f2);
  }





  /// a helper class for fixed order evaluation
  template <class REC, int N>
  class CEvalFO
  {
  public:
    template <class S, class T>
    INLINE static void Eval (S x, T && values, S & p1, S & p2) 
    {
      S p3;
      CEvalFO<REC,N-1>::Eval (x, values, p2, p3);
      if (REC::ZERO_B)
	values[N] = p1 = ( REC::CalcA(N) * x ) * p2 + REC::CalcC(N) * p3;
      else
	values[N] = p1 = ( REC::CalcA(N) * x + REC::CalcB(N)) * p2 + REC::CalcC(N) * p3;
    }


    template <class S, class Sc, class T>
    INLINE static void EvalMult (S x, Sc c, T && values, S & p1, S & p2) 
    {
      S p3;
      CEvalFO<REC,N-1>::EvalMult (x, c, values, p2, p3);
      if (REC::ZERO_B)
	values[N] = p1 = ( REC::CalcA(N) * x ) * p2 + REC::CalcC(N) * p3;
      else
	values[N] = p1 = ( REC::CalcA(N) * x + REC::CalcB(N)) * p2 + REC::CalcC(N) * p3;
    }
    

    template <class S, class Sy, class T>
    INLINE static void EvalScaled (S x, Sy y, T && values, S & p1, S & p2) 
    {
      S p3;
      CEvalFO<REC,N-1>::EvalScaled (x, y, values, p2, p3);
      values[N] = p1 = ( REC::CalcA(N) * x + REC::CalcB(N) * y) * p2 + REC::CalcC(N)*(y*y) * p3;
    }


    template <class S, class Sy, class Sc, class T>
    INLINE static void EvalScaledMult (S x, Sy y, Sc c, T && values, S & p1, S & p2) 
    {
      S p3;
      CEvalFO<REC,N-1>::EvalScaledMult (x, y, c, values, p2, p3);
      values[N] = p1 = ( REC::CalcA(N) * x + REC::CalcB(N) * y) * p2 + REC::CalcC(N)*(y*y) * p3;
    }
  };


  template <class REC>
  class CEvalFO<REC, -1>
    {
    public:
      template <class S, class T>
      INLINE static void Eval (S x, T && values, S & /* p1 */, S & /* p2 */) 
      { ; }

      template <class S, class Sc, class T>
      INLINE static void EvalMult (S x, Sc c, T && values, S & /* p1 */, S & /* p2 */) 
      { ; }


      template <class S, class Sy, class T>
      INLINE static void EvalScaled (S x, Sy y, T && values, S & /* p1 */, S & /* p2 */) 
      { ; }

      template <class S, class Sy, class Sc, class T>
      INLINE static void EvalScaledMult (S x, Sy y, Sc c, T && values, S & /* p1 */, S & /* p2 */) 
      { ; }

    };


  template <class REC>
  class CEvalFO<REC, 0>
    {
    public:
      template <class S, class T>
      INLINE static void Eval (S x, T && values, S & p1, S & /* p2 */) 
      {
	values[0] = p1 = REC::P0(x);
      }

      template <class S, class Sc, class T>
      INLINE static void EvalMult (S x, Sc c, T && values, S & p1, S & /* p2 */) 
      {
	values[0] = p1 = c * REC::P0(x);
      }


      template <class S, class Sy, class T>
      INLINE static void EvalScaled (S x, Sy y, T && values, S & p1, S & /* p2 */) 
      {
	values[0] = p1 = REC::P0(x);
      }

      template <class S, class Sy, class Sc, class T>
      INLINE static void EvalScaledMult (S x, Sy y, Sc c, T && values, S & p1, S & /* p2 */) 
      {
	values[0] = p1 = c * REC::P0(x);
      }

    };

  template <class REC>
  class CEvalFO<REC, 1>
  {
  public:
    template <class S, class T>
    INLINE static void Eval (S x, T && values, S & p1, S & p2) 
    {
      values[0] = p2 = REC::P0(x);
      values[1] = p1 = REC::P1(x);
    }

    template <class S, class Sc, class T>
    INLINE static void EvalMult (S x, Sc c, T && values, S & p1, S & p2) 
    {
      values[0] = p2 = c * REC::P0(x);
      values[1] = p1 = c * REC::P1(x);
    }

    template <class S, class Sy, class T>
    INLINE static void EvalScaled (S x, Sy y, T && values, S & p1, S & p2) 
    {
      values[0] = p2 = REC::P0(x);
      values[1] = p1 = REC::P1(x,y);
    }

    template <class S, class Sy, class Sc, class T>
    INLINE static void EvalScaledMult (S x, Sy y, Sc c, T && values, S & p1, S & p2) 
    {
      values[0] = p2 = c * REC::P0(x);
      values[1] = p1 = c * REC::P1(x,y);
    }
  };
  




  // P_i = (a_i x + b_i) P_{i-1} + c_i P_{i-2}
  template<class REC>
  class RecursivePolynomial
  {
  public:
    template <typename TINT, class S>
    INLINE static S EvalNext2 (TINT i, S x, S & p1, S & p2)
    {
      if (REC::ZERO_B)
        {
          // S pnew = REC::A(i) * x * p1 + REC::C(i) * p2; // maybe ?
          // S pnew = FMA(REC::C(i), p2, REC::A(i) * x * p1);  // bad
          S pnew = FMA(REC::A(i) * x, p1, REC::C(i) * p2);  // good
          p2 = p1;
          p1 = pnew;
        }
      else
        {
          S pnew = (REC::A(i) * x + REC::B(i)) * p1 + REC::C(i) * p2;
          p2 = p1;
          p1 = pnew;
        }
      return p1;
    }

#ifdef OLD
    template <class S>
    INLINE static S EvalNext (int i, S x, S & p1, S & p2)
    {
      if (i == 0) 
        {
          p1 = REC::P0(x);
          return p1;
        }
      if (i == 1) 
        {
          p2 = p1;
          p1 = REC::P1(x);
          return p1;
        }
      return EvalNext2 (i, x, p1, p2);
      /*
      if (REC::ZERO_B)
        {
          S pnew = REC::A(i) * x * p1 + REC::C(i) * p2;
          p2 = p1;
          p1 = pnew;
        }
      else
        {
          S pnew = (REC::A(i) * x + REC::B(i)) * p1 + REC::C(i) * p2;
          p2 = p1;
          p1 = pnew;
        }
      return p1;
      */
    }
#endif



    template <class S>
    INLINE static S EvalNextTicTac (int i, S x, S & p1, S & p2)
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
    
    template <typename TI, class S>
    INLINE static S EvalNextTicTac2 (TI i, S x, S & p1, S & p2)
    {
      if (REC::ZERO_B)
	{
          /*
	  p1 *= REC::C(i);
	  p1 += REC::A(i) * x * p2;
	  return p1;
          */
          p1 = FMA (REC::A(i)*x, p2, REC::C(i)*p1);
          return p1;
	}
      else
	{
	  p1 *= REC::C(i);
	  p1 += (REC::A(i) * x + REC::B(i)) * p2;
	  return p1;
	}
    }

    template <typename TI, class S, class Sy>
    INLINE static S EvalScaledNextTicTac2 (TI i, S x, Sy y, S & p1, S & p2)
    {
      // p1 = (y*y) * REC::C(i) * p1;
      if (REC::ZERO_B)
        {
          // p1 = (y*y) * REC::C(i) * p1;
          // p1 += REC::A(i) * x * p2;
          double a = REC::A(i);
          double c = REC::C(i);
          p1 = FMA(a*x, p2, c*(y*y)*p1);
        }
      else
        {
          p1 = (y*y) * REC::C(i) * p1;
          p1 += (REC::A(i) * x + REC::B(i) * y) * p2;
        }
      return p1;
    }

    /*
    template <class S, class Sc>
    INLINE static S EvalNextMultTicTac (int i, S x, Sc c, S & p1, S & p2)
    {
      switch (i)
        {
        case 0: return p1 = c * REC::P0(x);
        case 1: return p1 = c * REC::P1(x);
        default: return EvalNextTicTac2 (i, x, p1, p2);
        }
    }
    */

    template <class S, class Sy>
    INLINE static void EvalScaledNext (int i, S x, Sy y, S & p1, S & p2)
    {
      if (i == 0) 
        {
          p1 = REC::P0(x);
          return;
        }
      if (i == 1) 
        {
          p2 = p1;
          p1 = REC::P1(x,y);
          return;
        }

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

    template <typename TINT, class S, class Sy>
    INLINE static S EvalScaledNext2 (TINT i, S x, Sy y, S & p1, S & p2)
    {
      if (REC::ZERO_B)
        {
          // S pnew = REC::A(i) * x * p1 + REC::C(i) * (y*y)*p2;
          S pnew = FMA( REC::A(i)*x, p1, REC::C(i)*(y*y)*p2);
          p2 = p1;
          p1 = pnew;
        }
      else
        {
          S pnew = (REC::A(i) * x + REC::B(i) * y) * p1 + REC::C(i) * (y*y)*p2;
          p2 = p1;
          p1 = pnew;
        }
      return p1;
      /*
      // switch (i)
        {
          // case 0: return p1 = REC::P0(x);
          // case 1: return p1 = REC::P1(x);
          // default:
          {
            if (REC::ZERO_B)
              {
                p1 *= REC::C(i) *(y*y);
                p1 += REC::A(i) * x * p2;
                return p1;
              }
            else
              {
                p1 *= REC::C(i) *(y*y);
                p1 += (REC::A(i) * x + REC::B(i) * y) * p2;
                return p1;
              }
          }
        }
      */
    }

    template <class S, class Sy>
    static INLINE S P1(S x, Sy y)  { return REC::P1(x); }

  public:

    template <typename TI, class S, class T>
    INLINE static void Eval (TI n, S x, T && values) 
    {
      // if (n < 0) return;
      EvalMult (n, x, 1.0, values);
    }

    template <typename TI, class S, class Sc, class T>
    INLINE static void EvalMult (TI n, S x, Sc c, T && values) 
    {
      /*
      S p1, p2;
      if (n < 0) return;

      values[0] = EvalNextMultTicTac(0, x, c, p1, p2);
      if (n < 1) return;

      values[1] = EvalNextMultTicTac(1, x, c, p2, p1);
      if (n < 2) return;

      for (int i = 2; i < n; i+=2)
	{	
	  values[i] = EvalNextTicTac2 (i, x, p1, p2);
	  values[i+1] = EvalNextTicTac2 (i+1, x, p2, p1);
	}
      if ( (n&1) == 0 )
        values[n] = EvalNextTicTac2 (n, x, p1, p2);
      */



      S p1(c * REC::P0(x));
      S p2 = c * REC::P1(x);

      TI i = 0;
      for ( ; i < n; i+=2)
	{	
          values[i] = p1;
          values[i+1] = p2;
          
	  EvalNextTicTac2 (i+2, x, p1, p2);
	  EvalNextTicTac2 (i+3, x, p2, p1);
	}
      if (i == n)
        values[n] = p1;
    }


    template <typename TI, class S, class Sc, class T1, class T2>
    INLINE static void EvalMult (TI n, S x, Sc c, const Class_SBLambdaDuo<T1,T2> & values) 
    {
      /*
      S p1, p2;
      if (n < 0) return;
      values[0] = EvalNextMultTicTac(0, x, c, p1, p2);
      if (n < 1) return;
      values[1] = EvalNextMultTicTac(1, x, c, p2, p1);
      if (n < 2) return;
      */
      
      S p1(c * REC::P0(x));
      S p2 = c * REC::P1(x);

      if (n < 0) return;
      if (n == 0) { values [0] = p1; return; }
      values (0, p1, 1, p2);
      TI i = 2;
      for ( ; i < n; i+=2)
	{	
	  S v1  = EvalNextTicTac2 (i, x, p1, p2);
	  S v2 = EvalNextTicTac2 (i+1, x, p2, p1);
	  values (i, v1, i+1, v2);
	}
      if (i <= n)
        values[i] = EvalNextTicTac2 (i, x, p1, p2);
    }


    template <int N, class S, class Sc, class T>
    INLINE static void EvalMult (IC<N> n, S x, Sc c, T && values) 
    {
      S p1(c * REC::P0(x));
      S p2 = c * REC::P1(x);

      int i = 0;
      for ( ; i < n; i+=2)
	{	
          values[i] = p1;
          values[i+1] = p2;
          
	  EvalNextTicTac2 (i+2, x, p1, p2);
	  EvalNextTicTac2 (i+3, x, p2, p1);
	}
      if (i == n)
        values[n] = p1;

      /*
        // needs checking, code looks the same
      S p2(c * REC::P0(x));
      S p1(c * REC::P1(x));
      
      Iterate<N+1> ([&] (auto i)
                  {
                    values[i] = p2;
                    EvalNext2 (i+2, x, p1, p2);
                  });
      */
    }


    template <typename TI, class S, class T>
    INLINE static void Eval1Assign (TI n, S x, T && values)
    {
      EvalMult1Assign (n, x, 1.0, values);
    }

    template <typename TI, class S, class Sc, class T>
    INLINE static void EvalMult1Assign (TI n, S x, Sc c, T && values)
    {
      S p1 = c*REC::P1(x), p2(c * REC::P0(x));
      for (TI i = 0; i <= n; i++)
        {
	  values[i] = p2;
          EvalNext2 (i+2, x, p1, p2);
        }
    }



    template <typename TI, class S, class Sy, class T>
    INLINE static void EvalScaled (TI n, S x, Sy y, T && values)
    {
      if (n < 0) return;
      EvalScaledMult (n, x, y, 1.0, values);
    }

    template <typename TI, class S, class Sy, class Sc, class T>
    INLINE static void EvalScaledMult (TI n, S x, Sy y, Sc c, T && values)
    {
      /*
      S p1, p2;
      
      // if (n < 0) return;

      values[0] = p2 = c * REC::P0(x);

      if (n < 1) return;
      values[1] = p1 = c * REC::P1(x);

      for (int i = 2; i <= n; i++)
	{	
	  EvalScaledNext2 (i, x, y, p1, p2);
	  values[i] = p1;
	}
      */


      S p1(c * REC::P0(x));
      S p2(c * REC::P1(x,y));

      TI i = 0;
      for ( ; i < n; i+=2)
	{	
          values[i] = p1;
          values[i+1] = p2;
          
	  EvalScaledNextTicTac2 (i+2, x, y, p1, p2);
	  EvalScaledNextTicTac2 (i+3, x, y, p2, p1);
	}
      if (i == n)
        values[n] = p1;


      // EvalScaledMult1Assign (n, x, y, c, values);
    }


    template <int N, class S, class Sy, class Sc, class T>
    HD INLINE static void EvalScaledMult (IC<N> n, S x, Sy y, Sc c, T && values)
    {
      S p1(c*REC::P1(x,y)), p2(c * REC::P0(x));
      Iterate<N+1> ([&] (auto i) LAMBDA_INLINE
                  {
                    values[i] = p2;
                    EvalScaledNext2 (i+2, x, y, p1, p2);
                  });
    }  

    
    template <typename TI, class S, class Sy, class T>
    INLINE static void EvalScaled1Assign (TI n, S x, Sy y, T && values)
    {
      // if (n  < 0) return;
      EvalScaledMult1Assign (n, x, y, 1.0, values);
    }

    template <typename TI, class S, class Sy, class Sc, class T>
    INLINE static void EvalScaledMult1Assign (TI n, S x, Sy y, Sc c, T && values)
    {
      S p1(c*REC::P1(x,y)), p2(c * REC::P0(x));
      if (n < 0) return;

      TI i = 0;
      while (true)
        {
	  values[i] = p2;
          if (i == n) break;
          EvalScaledNext2 (i+2, x, y, p1, p2);
          i++;
        }

      /*
      TI i = 0;
      goto lab1;
      while (i <= n)
        {
          EvalScaledNext2 (i+1, x, y, p1, p2);
        lab1:
	  values[i] = p2;
          i++;
        }
      */
      
    }


// C++14 only
#if __cplusplus > 201103L
    template <int N, class S, class T>
    INLINE static void Eval (IC<N> n, S x, T && values) 
    {
      // S p1, p2;
      // CEvalFO<REC, N>::Eval (x, values, p1, p2);
      
      S p1 = REC::P1(x), p2 = REC::P0(x);
      Iterate<n+1> ( [&] (auto i) LAMBDA_INLINE
        {
	  values[i] = p2;
          RecursivePolynomial<REC>::EvalNext2 (IC<i+2>(), x, p1, p2);
        } );
    }
    
    template <int N, class S, class Sy, class T>
    INLINE static void EvalScaled (IC<N> n,
                                   S x, Sy y, T && values) 
    {
      S p1(REC::P1(x,y)), p2(REC::P0(x));
      Iterate<n+1> ( [&] (auto i) LAMBDA_INLINE
        {
          // cout << "eval scaled, type(i) = " << typeid(i).name() << endl;
	  values[i] = p2;
          RecursivePolynomial<REC>::EvalScaledNext2 (i+IC<2>(), x, y, p1, p2);
        });
    }
#endif

    template <int N, class S, class T>
    INLINE static void EvalFO (S x, T && values) 
    {
      S p1, p2;
      CEvalFO<REC, N>::Eval (x, values, p1, p2);
    }

    template <int N, class S, class Sc, class T>
    INLINE static void EvalMultFO (S x, Sc c, T && values) 
    {
      S p1, p2;
      CEvalFO<REC, N>::EvalMult (x, c, values, p1, p2);
    }

    template <int N, class S, class Sy, class T>
    INLINE static void EvalScaledFO (S x, Sy y, T && values) 
    {
      S p1, p2;
      CEvalFO<REC, N>::EvalScaled (x, y, values, p1, p2);
    }

    template <int N, class S, class Sy, class Sc, class T>
    INLINE static void EvalScaledMultFO (S x, Sy y, Sc c,T && values) 
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

    template <typename TI>
    INLINE void ABC (TI i, double & a, double & b, double & c) const   
    {
      a = static_cast<const REC&>(*this).A(i);
      b = static_cast<const REC&>(*this).B(i);
      c = static_cast<const REC&>(*this).C(i); 
      double d = 1.0 / static_cast<const REC&>(*this).D(i);
      a *= d; b *= d; c *= d;
    }

    INLINE const REC & Cast() const { return static_cast<const REC &> (*this); }

    template <typename TI, class S>
    INLINE S EvalNext2 (TI i, S x, S & p1, S & p2) const
    {
      if (REC::ZERO_B)
	{
	  S pnew = Cast().A(i) * x * p1 + Cast().C(i)*p2;
	  p2 = p1;
	  p1 = pnew;
	  return p1;
	}
      else
	{
          // S pnew = (Cast().A(i) * x + Cast().B(i)) * p1 + Cast().C(i)*p2;
          // S pnew = FMA(Cast().C(i), p2, (Cast().A(i) * x + Cast().B(i)) * p1); // bad
          S pnew = FMA( FMA(Cast().A(i), x, S(Cast().B(i))), p1, Cast().C(i) * p2);  // good          
          /*
          S axpb = (Cast().A(i) * x + Cast().B(i));
          S hv = Cast().C(i)*p2;
          S pnew = hv + axpb * p1;
          */
	  p2 = p1;
	  p1 = pnew;
	  return p1;
	}
    }


    /*
    template <class S>
    INLINE S EvalNext (int i, S x, S & p1, S & p2) const
    {
      switch (i)
        {
        case 0: return p1 = Cast().P0(x);
        case 1: p2 = p1; return p1 = Cast().P1(x);
        default:
          {
            if (REC::ZERO_B)
              {
                S pnew = Cast().A(i) * x * p1 + Cast().C(i)*p2;
                p2 = p1;
                p1 = pnew;
                return p1;
              }
            else
              {
                S pnew = (Cast().A(i) * x + Cast().B(i)) * p1 + Cast().C(i)*p2;
                p2 = p1;
                p1 = pnew;
                return p1;
              }
          }
        }
    }
    */


    template <class S, class Sc>
    INLINE S EvalNextMult (int i, S x, Sc c, S & p1, S & p2)
    {
      switch (i)
        {
        case 0: return p1 = c * Cast().P0(x);
        case 1: p2 = p1; return p1 = c * Cast().P1(x);
        default:
          {
            S pnew = REC::ZERO_B ?
              Cast().A(i) * x * p1 + Cast().C(i)*p2 :
              (Cast().A(i) * x + Cast().B(i)) * p1 + Cast().C(i)*p2;
            p2 = p1; p1 = pnew;
            return p1;
          }
        }
    }


    template <typename TI, class S>
    INLINE S EvalNextTicTac2 (TI i, S x, S & p1, S & p2) const
    {
      double a,b,c;
      static_cast<const REC&> (*this).ABC (i, a, b, c);
      if (REC::ZERO_B)
	{
	  p1 *= c;
	  p1 += a*x*p2;
	  return p1;
	}
      else
	{
          /*
	  p1 *= c;
	  p1 += (a * x + b) * p2;
          */
          p1 = FMA(FMA(a,x,S(b)), p2, c*p1);
	  return p1;
	}
    }

    /*
    template <class S>
    INLINE S EvalNextTicTac (int i, S x, S & p1, S & p2) const
    {
      switch (i)
        {
        case 0: return p1 = static_cast<const REC&>(*this).P0(x);
        case 1: return p1 = static_cast<const REC&>(*this).P1(x);
        default:
          {
            double a,b,c;
            static_cast<const REC&> (*this).ABC (i, a, b, c);
            
            if (REC::ZERO_B)
              {
                p1 *= c;
                p1 += a*x*p2;
                return p1;
              }
            else
              {
                p1 *= c;
                p1 += (a * x + b) * p2;
                return p1;
              }
          }
        }
    }
    */

    /*
    template <class S, class Sc>
    INLINE S EvalNextMultTicTac (int i, S x, Sc c, S & p1, S & p2) const
    {
      switch (i)
        {
        case 0: return p1 = c * static_cast<const REC&>(*this).P0(x);
        case 1: return p1 = c * static_cast<const REC&>(*this).P1(x);
        default: return EvalNextTicTac2 (i, x, p1, p2);
        }
    }
    */

    template <typename TI, class S, class Sy>
    INLINE S EvalScaledNext2 (TI i, S x, Sy y, S & p1, S & p2) const
    {
      if (REC::ZERO_B)
        {
          S pnew = Cast().A(i) * x * p1 + Cast().C(i) * (y*y)*p2;
          // pnew *= 1.0 / Cast().D(i);
          p2 = p1;
          p1 = pnew;
        }
      else
        {
          S pnew = (Cast().A(i) * x + Cast().B(i) * y) * p1 + Cast().C(i) * (y*y)*p2;
          // pnew *= 1.0 / static_cast<const REC&>(*this).D(i);
          p2 = p1;
          p1 = pnew;
        }
      return p1;
    }

    /*
    // new
    template <class S, class Sy, class Tc>
    INLINE S EvalScaledMultNext (int i, S x, Sy y, Tc c, S & p1, S & p2) const
    {
      switch (i)
        {
        case 0: 
          return p1 = c * static_cast<const REC&>(*this).P0(x);
        case 1: 
          p2 = p1; return p1 = c * static_cast<const REC&>(*this).P1(x,y);
        default:
          {
            if (REC::ZERO_B)
              {
                S pnew = Cast().A(i) * x * p1 + Cast().C(i) * (y*y)*p2;
                pnew *= 1.0 / Cast().D(i);
                p2 = p1;
                p1 = pnew;
              }
            else
              {
                S pnew = (Cast().A(i) * x + Cast().B(i) * y) * p1 + Cast().C(i) * (y*y)*p2;
                pnew *= 1.0 / static_cast<const REC&>(*this).D(i);
                p2 = p1;
                p1 = pnew;
              }
            return p1;
          }
        }
    }
    */


  public:

    template <typename TI, class S, class T>
    INLINE void Eval (TI n, S x, T && values) const
    {
      if (n < 0) return;
      EvalMult (n, x, 1.0, values);
    }

    template <typename TI, class S, class T>
    INLINE void Eval1Assign (TI n, S x, T && values) const
    {
      // if (n < 0) return;
      EvalMult1Assign (n, x, 1.0, values);
    }

    template <typename TI, class S, class Sc, class T>
    INLINE void EvalMult (TI n, S x, Sc c, T && values) const
    {
      /*
      S p1, p2;

      if (n < 0) return;

      values[0] = EvalNextMultTicTac(0, x, c, p1, p2);
      if (n < 1) return;

      values[1] = EvalNextMultTicTac(1, x, c, p2, p1);
      if (n < 2) return;

      // for (int i = 2; i <= n; i++)
      // values[i] = EvalNext2 (i, x, p1, p2);

      for (int i = 2; i <= n; i+=2)
	{	
	  values[i] = EvalNextTicTac2 (i, x, p1, p2);
          if (i == n) break;
	  values[i+1] = EvalNextTicTac2 (i+1, x, p2, p1);
	}
      */

      /*
      // leading !!!!
      // if (n < 0) return;
      S p1(c * static_cast<const REC&>(*this).P0(x));
      values[0] = p1;
      if (n < 1) return;

      S p2(c * static_cast<const REC&>(*this).P1(x));
      values[1] = p2;

      for (int i = 2; i <= n; i+=2)
	{	
	  values[i] = EvalNextTicTac2 (i, x, p1, p2);
          if (i == n) break;
	  values[i+1] = EvalNextTicTac2 (i+1, x, p2, p1);
	}
      */


      // if (n < 0) return;
      S p1(c * static_cast<const REC&>(*this).P0(x));
      values[0] = p1;
      if (n < 1) return;

      S p2(c * static_cast<const REC&>(*this).P1(x));
      values[1] = p2;

      TI i = 2;
      if ( (n&1) == 0)
        {
          values[i] = EvalNext2 (i, x, p2, p1);
          i++;
        }
      for ( ; i < n; i+=2)
	{	
	  values[i] = EvalNextTicTac2 (i, x, p1, p2);
	  values[i+1] = EvalNextTicTac2 (i+1, x, p2, p1);
	}


      
      /*
      // if (n < 0) return;
      S p1(c * static_cast<const REC&>(*this).P0(x));
      values[0] = p1;
      if (n < 1) return;

      S p2(c * static_cast<const REC&>(*this).P1(x));
      values[1] = p2;
      // if (n < 2) return;

      for (int i = 2; i <= n; i++)
        values[i] = EvalNext2 (i, x, p2, p1);
      */
      
      /*
      S p1(c * static_cast<const REC&>(*this).P0(x));
      S p2(c * static_cast<const REC&>(*this).P1(x));
      TI i = 0;
      for ( ; i < n; i+=2)
	{	
          values[i] = p1;
	  EvalNextTicTac2 (i+2, x, p1, p2);
          values[i+1] = p2;
          EvalNextTicTac2 (i+3, x, p2, p1);
	}
      if (i == n)
        values[n] = p1;
      */

      
      /*
      S p1(c * static_cast<const REC&>(*this).P0(x));
      S p2(c * static_cast<const REC&>(*this).P1(x));
      TI i = 0;
      for ( ; i <= n; i+=2)
	{	
          values[i] = p1;
          if (i == n) break;
	  EvalNextTicTac2 (i+2, x, p1, p2);
          values[i+1] = p2;
          EvalNextTicTac2 (i+3, x, p2, p1);
	}
      */
    }

    template <int N, class S, class Sc, class T>
    INLINE void EvalMult (IC<N> n, S x, Sc c, T && values) const
    {
      // S p1(c*REC::P1(x)), p2(c * REC::P0(x));

      S p2(c * static_cast<const REC&>(*this).P0(x));
      S p1(c * static_cast<const REC&>(*this).P1(x));
      
      Iterate<N+1> ([&] (auto i)
                  {
                    values[i] = p2;
                    this->EvalNext2 (i+2, x, p1, p2);
                  });
    }

    template <int N, class S, class Sc, class T1, class T2>
    INLINE void EvalMult (IC<N> n, S x, Sc c, Class_SBLambdaDuo<T1,T2> && values) const
    {
      S p1(c * static_cast<const REC&>(*this).P0(x));
      S p2(c * static_cast<const REC&>(*this).P1(x));

      size_t i = 0;
      for ( ; i < n; i+=2)
	{	
          values (i, p1, i+1, p2);
	  EvalNextTicTac2 (i+2, x, p1, p2);
	  EvalNextTicTac2 (i+3, x, p2, p1);
	}
      if (i == n)
        values[n] = p1;
    }
    


    
    template <typename TI, class S, class Sc, class T>
    INLINE void EvalMult1Assign (TI n, S x, Sc c, T && values) const
    {
      /*
      S p1(0.0), p2(0.0); // initialize for suppressing warning
      for (int i = 0; i <= n; i++)
        values[i] = EvalNextMult (i, x, c, p1, p2);
      */

      S p1(c * static_cast<const REC&>(*this).P1(x)), 
        p2(c * static_cast<const REC&>(*this).P0(x));
      for (TI i = 0; i <= n; i++)
        {
	  values[i] = p2;
          EvalNext2 (i+2, x, p1, p2);
        }

    }
    
    template <typename TI, class S, class Sc, class T1, class T2>
    INLINE void EvalMult (TI n, S x, Sc c, Class_SBLambdaDuo<T1,T2> && values) const
    {
      S p1(c * static_cast<const REC&>(*this).P0(x));
      S p2(c * static_cast<const REC&>(*this).P1(x));

      TI i = 0;
      for ( ; i < n; i+=2)
	{	
          values (i, p1, i+1, p2);
	  EvalNextTicTac2 (i+2, x, p1, p2);
	  EvalNextTicTac2 (i+3, x, p2, p1);
	}
      if (i == n)
        values[n] = p1;
    }

      
    template <typename TI, class S, class Sy, class T>
    INLINE void EvalScaled (TI n, S x, Sy y, T && values) const
    {
      EvalScaledMult (n, x, y, 1.0, values);
    }
    
    template <class S, class Sy>
    INLINE S P1(S x, Sy y) const 
    { 
      return static_cast<const REC&>(*this).P1(x);
    }

    template <typename TI, class S, class Sy, class Sc, class T>
    INLINE void EvalScaledMult (TI n, S x, Sy y, Sc c, T && values) const
    {
      S p1, p2;

      if (n < 0) return;

      values[0] = p2 = c * static_cast<const REC&>(*this).P0(x);
      if (n < 1) return;

      values[1] = p1 = c * static_cast<const REC&>(*this).P1(x,y);

      for (TI i = 2; i <= n; i++)
        values[i] = EvalScaledNext2 (i, x, y, p1, p2);
    }

    template <int N, class S, class Sy, class Sc, class T>
    INLINE void EvalScaledMult (IC<N> n, S x, Sy y, Sc c, T && values) const
    {
      S p1(c*Cast().P1(x,y)), p2(c * Cast().P0(x));
      Iterate<N+1> ([&] (auto i) LAMBDA_INLINE
                  {
                    values[i] = p2;
                    this->EvalScaledNext2 (i+2, x, y, p1, p2);
                  });
    }  
    
    template <typename TI, class S, class Sy, class Sc, class T>
    INLINE void EvalScaledMult1Assign (TI n, S x, Sy y, Sc c, T && values) const
    {
      /*
      S p1(0.0), p2(0.0);
      for (int i = 0; i <= n; i++)
        values[i] = EvalScaledMultNext (i, x, y, c, p1, p2);
      */


      S p1 = c*static_cast<const REC&>(*this).P1(x,y), 
        p2 = c * static_cast<const REC&>(*this).P0(x);
      for (TI i = 0; i <= n; i++)
        {
	  values[i] = p2;
          EvalScaledNext2 (i+2, x, y, p1, p2);
        }
      
    }

    template <typename TI, class S>
    INLINE S CalcHighestOrder (TI n, S x) const
    {
      return CalcHighestOrderMult (n, x, 1.0);
    }
    
    template <typename TI, class S, class Sc>
    INLINE S CalcHighestOrderMult (TI n, S x, Sc c) const
    {
      S p1(c * static_cast<const REC&>(*this).P1(x)), 
        p2(c * static_cast<const REC&>(*this).P0(x));
      for (TI i = 0; i < n; i++)
        EvalNext2 (i+2, x, p1, p2);
      return p2;
    }    
    
    enum { ZERO_B = 0 };
  };









  /* ******************** Legendre Polynomial ****************** */

  class LegendrePolynomial_CalcCoefficient : public RecursivePolynomial<LegendrePolynomial_CalcCoefficient>
  {
  public:
    INLINE LegendrePolynomial_CalcCoefficient () { ; }

    template <class S, class T>
    INLINE LegendrePolynomial_CalcCoefficient (int n, S x, T && values)
    {
      Eval (n, x, values);
    }

    template <class S>
    static INLINE S P0(S x)  { return S(1.0); }
    template <class S>
    static INLINE S P1(S x)  { return x; }
    template <class S, class Sy>
    static INLINE S P1(S x, Sy y)  { return P1(x); }

    static INLINE double A (int i) { return 2.0-1.0/i; }
    static INLINE double B (int i) { return 0; }
    static INLINE double C (int i) { return 1.0/i-1.0; }
    enum { ZERO_B = 1 };
  };



  class NGS_DLL_HEADER LegendrePolynomial : public RecursivePolynomial<LegendrePolynomial>
  {
    static Array< Vec<2> > coefs;
    
  public:
    INLINE LegendrePolynomial () { ; }

    template <class S, class T>
    INLINE LegendrePolynomial (int n, S x, T && values)
    { 
      Eval (n, x, values);
    }

    static void Calc (int n);
    static FlatArray<Vec<2>> GetCoefs() { return coefs; }

    template <class S>
    static INLINE double P0(S x)  { return 1.0; }
    template <class S>
    static INLINE S P1(S x)  { return x; }
    template <class S, class Sy>
    static INLINE S P1(S x, Sy y)  { return P1(x); }


    static INLINE double CalcA (int i) { return 2.0-1.0/i; }
    static INLINE double CalcB (int i) { return 0; }
    static INLINE double CalcC (int i) { return 1.0/i-1.0; }

    template <typename TI>
      static INLINE double A (TI i) { return coefs[i][0]; } 
    template <typename TI>
      static INLINE double B (TI i) { return 0; }
    template <typename TI>    
      static INLINE double C (TI i) { return coefs[i][1]; } 
    
    template <int I>
    static INLINE double A (IC<I> i) { return 2.0-1.0/i; }
    template <int I>
    static INLINE double B (IC<I> i) { return 0; }
    template <int I>
    static INLINE double C (IC<I> i) { return 1.0/i-1.0; }

    
    enum { ZERO_B = 1 };
  };


  class LegendrePolynomialNonStatic : public RecursivePolynomialNonStatic<LegendrePolynomialNonStatic>
  {
    FlatArray<Vec<2>> coefs;
  public:
    LegendrePolynomialNonStatic ()
      : coefs(LegendrePolynomial::GetCoefs()) { ; }

    template <class S>
    INLINE double P0(S x) const { return 1.0; }
    template <class S>
    INLINE S P1(S x) const { return x; }
    template <class S, class Sy>
    INLINE S P1(S x, Sy y) const { return P1(x); }
    
    template <typename TI>
    INLINE double A (TI i) const { return coefs[i][0]; } 
    template <typename TI>
    INLINE double B (TI i) const { return 0; }
    template <typename TI>    
    INLINE double C (TI i) const { return coefs[i][1]; } 
  };
  

#ifdef __CUDACC__    
  extern __device__ Vec<2> * intlegnobubble_coefs;
#endif

  class IntLegNoBubble : public RecursivePolynomial<IntLegNoBubble>
  {
    static Array< Vec<2> > coefs;

  public:
    IntLegNoBubble () { ; }
    
    template <class S, class T>
    inline IntLegNoBubble (int n, S x, T && values)
    {
      Eval (n, x, values);
    }
    
    static void Calc (int n);

    template <class S>
    static INLINE double P0(S x)  { return -0.5; }
    template <class S>
    static INLINE S P1(S x)  { return -0.5*x; }
    template <class S, class Sy>
    static INLINE S P1(S x, Sy y)  { return P1(x); }

    static INLINE double A (int i) { return coefs[i][0]; } 
    static INLINE double B (int i) { return 0; }
    static INLINE double C (int i) { return coefs[i][1]; } 
    
    static INLINE double CalcA (int i) { i+=2; return (2*i-3)/double(i); }
    static INLINE double CalcB (int i) { return 0; }
    static INLINE double CalcC (int i) { i+=2; return -(i-3.0)/i; }
    enum { ZERO_B = 1 };
  };
    


  
  class IntegratedLegendrePolynomial : public RecursivePolynomial<IntegratedLegendrePolynomial>
  {
    // static Array< double[2] > coefs;
    
  public:
    IntegratedLegendrePolynomial () { ; }

    template <class S, class T>
    inline IntegratedLegendrePolynomial (int n, S x, T && values)
    { 
      Eval (n, x, values);
    }

    static void Calc (int n);

    template <class S>
    static INLINE S P0(S x)  { return S(-1.0); }
    template <class S>
    static INLINE S P1(S x)  { return x; }
    
    static INLINE double A (int i) { return (2*i-3)/double(i); }
    static INLINE double B (int i) { return 0; }
    static INLINE double C (int i) { return -(i-3.0)/i; }
    enum { ZERO_B = 1 };
  };



  class ChebyPolynomial : public RecursivePolynomial<ChebyPolynomial>
  {
  public:
    ChebyPolynomial () { ; }

    template <class S, class T>
    inline ChebyPolynomial (int n, S x, T && values)
    { 
      Eval (n, x, values);
    }

    template <class S>
    static INLINE double P0(S x)  { return 1.0; }
    template <class S>
    static INLINE S P1(S x)  { return x; }
    template <class S, class Sy>
    static INLINE S P1(S x, Sy y)  { return P1(x); }

    static INLINE double A (int i) { return 2; } 
    static INLINE double B (int i) { return 0; }
    static INLINE double C (int i) { return -1; }

    static INLINE double CalcA (int i) { return 2; } 
    static INLINE double CalcB (int i) { return 0; }
    static INLINE double CalcC (int i) { return -1; }

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
    static INLINE S P0(S x) { return S(1.0); }
    template <class S>
    // static INLINE S P1(S x) { return 0.5 * (2*(al+1)+(al+be+2)*(x-1)); }
    static INLINE S P1(S x) { return 0.5*(2*(al+1)-(al+be+2)) + 0.5*(al+be+2)*x; }
    
    template <class S>
    // static INLINE S P1(S x, S y) { return 0.5 * (2*(al+1)*y+(al+be+2)*(x-y)); }
    static INLINE S P1(S x, S y) { return 0.5*(al+be+2)*x+0.5*(al-be)*y; }
      
    static INLINE double CalcA (int i) 
    { i--; return (2.0*i+al+be)*(2*i+al+be+1)*(2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    static INLINE double CalcB (int i)
    { i--; return (2.0*i+al+be+1)*(al*al-be*be) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    static INLINE double CalcC (int i) 
    { i--; return -2.0*(i+al)*(i+be) * (2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }

    static INLINE double A (int i) { return CalcA (i); }
    static INLINE double B (int i) { return CalcB (i); }
    static INLINE double C (int i) { return CalcC (i); }
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
    INLINE S P0(S x) const { return S(1.0); }
    template <class S>
    INLINE S P1(S x) const { return 0.5 * (2*(al+1)+(al+be+2)*(x-1)); }
    template <class S>
    INLINE S P1(S x, S y) const { return 0.5 * (2*(al+1)*y+(al+be+2)*(x-y)); }
      
    INLINE double A (int i) const
    { i--; return (2.0*i+al+be)*(2*i+al+be+1)*(2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    INLINE double B (int i) const
    { i--; return (2.0*i+al+be+1)*(al*al-be*be) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    INLINE double C (int i) const
    { i--; return -2.0*(i+al)*(i+be) * (2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    INLINE double D (int i) const { return 1; }
  };




#ifdef __CUDACC__    
  // __device__ Vec<3> * jacobialpha_coefs;
  extern __device__ double jacobialpha_coefs[100][100][4];
  extern __device__ int jacobialpha_maxn;
#endif

  
  class JacobiPolynomialAlpha : public RecursivePolynomialNonStatic<JacobiPolynomialAlpha>
  {
  public:
    /*
    static Array< Vec<4> > coefs;
    static int maxnp;
    */
    static constexpr size_t maxnp = 128;
    static constexpr size_t maxalpha = 128;
    static Vec<4> coefs[maxnp*maxalpha];
    size_t n2;
    Vec<4> * coefsal;
  public:
    INLINE JacobiPolynomialAlpha (size_t a) 
    { 
      // offset = alpha*maxnp;
      coefsal = &coefs[a*maxnp];
      n2 = 2*maxnp;
      // coefsal = &jacobialpha_coefs[alpha*(jacobialpha_maxn+1)];      
    }

    template <class S, class T>
    INLINE JacobiPolynomialAlpha (int n, S x, int a, T && values)
    { 
      // offset = alpha*(maxn+1);
      coefsal = &coefs[a*maxnp];
      n2 = 2*maxnp;
      Eval (n, x, values);
    }

    void IncAlpha2 () { coefsal += n2; }
    static void Calc (int n, int maxalpha);

    template <class S>
    INLINE double P0(S x) const { return 1.0; }
    template <class S>
    INLINE S P1(S x) const 
    { 
      return FMA (coefsal[1][0],x,S(coefsal[1][1])); 
    }
    template <class S, class Sy>
    INLINE S P1(S x, Sy y) const 
    { 
      return FMA(coefsal[1][0],x,coefsal[1][1]*y);
    }


    template <typename TI>
    INLINE const double & A (TI i) const { return coefsal[i][0]; } 
    template <typename TI>
    INLINE const double & B (TI i) const { return coefsal[i][1]; } 
    template <typename TI>
    INLINE const double & C (TI i) const { return coefsal[i][2]; } 

    template <typename TI>
    INLINE void ABC (TI i, double & a, double & b, double & c) const   
    {
      a = A(i);
      b = B(i);
      c = C(i);
    }

    enum { ZERO_B = 0 };


    static INLINE double CalcA (int i, double al, double be) 
    { i--; return (2.0*i+al+be)*(2*i+al+be+1)*(2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    static INLINE double CalcB (int i, double al, double be)
    { i--; return (2.0*i+al+be+1)*(al*al-be*be) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    static INLINE double CalcC (int i, double al, double be)
    { i--; return -2.0*(i+al)*(i+be) * (2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }

    static INLINE double D (int i) { return 1; }
  };



  /*
// see Beuchler+Schoeberl : New shape functions for trigs, formula (2.11)
class IntJacobiPolynomialAlpha : public RecursivePolynomialNonStatic<IntJacobiPolynomialAlpha>
  {
  public:
    int alpha;

  public:
    IntJacobiPolynomialAlpha (int a) : alpha(a) 
    { 
      ;
    }

    template <class S, class T>
    inline IntJacobiPolynomialAlpha (int n, S x, int a, T && values)
      : alpha(a)
    { 
      Eval (n, x, values);
    }

    template <class S>
    INLINE double P0(S x) const { return 1.0; }

    // d/dx  ( (x+1) P1) =  0.5*(alpha+2) x + 0.5 * alpha;
    // (x+1) P1 = 0.25*(alpha+2) x^2 + 0.5 * alpha * x  - 0.25*(alpha+2) (-1)^2 - 0.5*alpha*(-1)
    // = 0.25*(alpha+2) (x-1)*(x+1)  + 0.5 * alpha (x+1)
    // P1 = 0.25*(alpha+2)(x-1) + 0.5 * alpha
    
    template <class S>
    INLINE S P1(S x) const 
    { 
      return 0.25 * (alpha+2)*x + 0.25*(alpha-2);
    }

    template <class S, class Sy>
    INLINE S P1(S x, Sy y) const 
    { 
      return 0.25 * (alpha+2)*x + 0.25*(alpha-2) * y;
    }

    INLINE double A (int i) const { return CalcA (i, alpha, 0); }
    INLINE double B (int i) const { return CalcB (i, alpha, 0); }
    INLINE double C (int i) const { return CalcC (i, alpha, 0); }

    enum { ZERO_B = 0 };

    static double CalcA (int n, double al, double be) 
    { return (2.0*n+al-1)*(2*n+al-2)*(2*n+al) / ( (2*n+2) * (n+al) * (2*n+al-2) ); }
    static double CalcB (int n, double al, double be)
    { return (2.0*n+al-1)*al * (al-2) / ( (2*n+2) * (n+al) * (2*n+al-2) ); }
    static double CalcC (int n, double al, double be)
    { return - 2*(n-1)*(n+al-2)*(2*n+al) / ( (2*n+2) * (n+al) * (2*n+al-2) ); }

    INLINE double D (int i) const { return 1; }
  };
  */




// see Beuchler+Schoeberl : New shape functions for trigs, formula (2.11)
// \int_(-1)^x P^(alpha,0)(s) ds   /  (1-s)
class IntegratedJacobiPolynomialAlpha : public RecursivePolynomialNonStatic<IntegratedJacobiPolynomialAlpha>
  {
  public:
    static Array< Vec<4> > coefs;
    int alpha;
    static int maxn;
    Vec<4> * coefsal;

  public:
    IntegratedJacobiPolynomialAlpha (int a) : alpha(a) 
    { 
      coefsal = &coefs[alpha*(maxn+1)];
    }

    template <class S, class T>
    inline IntegratedJacobiPolynomialAlpha (int n, S x, int a, T && values)
      : alpha(a)
    { 
      coefsal = &coefs[alpha*(maxn+1)];
      Eval (n, x, values);
    }

    static void Calc (int n, int maxalpha);

    template <class S>
    INLINE double P0(S x) const { return 1.0; }

    // d/dx  ( (x+1) P1) =  0.5*(alpha+2) x + 0.5 * alpha;
    // (x+1) P1 = 0.25*(alpha+2) x^2 + 0.5 * alpha * x  - 0.25*(alpha+2) (-1)^2 - 0.5*alpha*(-1)
    // = 0.25*(alpha+2) (x-1)*(x+1)  + 0.5 * alpha (x+1)
    // P1 = 0.25*(alpha+2)(x-1) + 0.5 * alpha
    
    template <class S>
    INLINE S P1(S x) const 
    { 
      return coefsal[1][0]*x+coefsal[1][1]; 
      // return 0.25 * (alpha+2)*x + 0.25*(alpha-2);
    }

    template <class S, class Sy>
    INLINE S P1(S x, Sy y) const 
    { 
      return coefsal[1][0]*x+coefsal[1][1]*y; 
      // return 0.25 * (alpha+2)*x + 0.25*(alpha-2) * y;
    }


    INLINE double A (int i) const { return coefsal[i][0]; } 
    INLINE double B (int i) const { return coefsal[i][1]; } 
    INLINE double C (int i) const { return coefsal[i][2]; } 
    /*
    INLINE double A (int i) const { return CalcA (i, alpha, 0); }
    INLINE double B (int i) const { return CalcB (i, alpha, 0); }
    INLINE double C (int i) const { return CalcC (i, alpha, 0); }
    */
    enum { ZERO_B = 0 };

    static double CalcA (int n, double al, double be) 
    { return (2.0*n+al-1)*(2*n+al-2)*(2*n+al) / ( (2*n+2) * (n+al) * (2*n+al-2) ); }
    static double CalcB (int n, double al, double be)
    { return (2.0*n+al-1)*al * (al-2) / ( (2*n+2) * (n+al) * (2*n+al-2) ); }
    static double CalcC (int n, double al, double be)
    { return - 2*(n-1)*(n+al-2)*(2*n+al) / ( (2*n+2) * (n+al) * (2*n+al-2) ); }

    INLINE double D (int i) const { return 1; }
  };










  class JacobiPolynomialNew : public RecursivePolynomialNonStatic<JacobiPolynomialNew>
  {
  protected:
    double al, be;

  public:
    INLINE JacobiPolynomialNew (double alpha, double beta) : al(alpha), be(beta) { ; }

    template <class S>
    INLINE double P0(S x) const { return 1.0; }
    template <class S>
    INLINE S P1(S x) const 
    { 
      return 0.5 * (2*(al+1)+(al+be+2)*(x-1));
    }
    template <class S>
    INLINE S P1(S x, S y) const 
    { 
      return 0.5 * (2*(al+1)*y+(al+be+2)*(x-y));
    }

    enum { ZERO_B = 0 };

    INLINE double A (int i) const
    { i--; return (2.0*i+al+be)*(2*i+al+be+1)*(2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    INLINE double B (int i) const
    { i--; return (2.0*i+al+be+1)*(al*al-be*be) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }
    INLINE double C (int i) const
    { i--; return -2.0*(i+al)*(i+be) * (2*i+al+be+2) / ( 2 * (i+1) * (i+al+be+1) * (2*i+al+be)); }

    INLINE double D (int i) const { return 1; }
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
  INLINE void LegendrePolynomialMult (int n, S x, Sc c , T && values)
  {
    LegendrePolynomial leg;
    leg.EvalMult (n, x, c, values);
  }




  template <class S, class T>
  inline void JacobiPolynomial (int n, S x, double alpha, double beta, T && values)
  {
    S p1(1.0), p2(0.0), p3;

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
    S p1 = c, p2(0.0), p3;

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

    S p1(1.0), p2(0.0), p3;

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

    S p1 = c, p2(0.0), p3;

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
      return power * (1.0-x)*0.5;
    }


    template <class S, class St, class T>
    static S EvalScaled (S x, St t, T && values)
    {
      S power = DubinerJacobiPolynomialsPowFO<n, i-1, alpha0, beta>::EvalScaled (x, t, values);
      S p1, p2;
      CEvalFO<JacobiPolynomialFix<alpha0+2*i, beta>, n-i>::EvalScaledMult (x, t, power, values.Row(i), p1, p2);
      return power * (1.0-x)*0.5;
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
    INLINE static void Step (S x, Thelp & help, T & values)
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
	  REC::EvalNextTicTac (ORDER, x, help[DIAG-ORDER][0], help[DIAG-ORDER][1]);
	}
      values(DIAG-ORDER, ORDER) = help[DIAG-ORDER][0];
    }
  };

  template <int ALPHA0, int BETA, int DIAG>
  class DubinerJacobiPolynomialsDiag<ALPHA0, BETA, DIAG, -1>
  {
  public:
    template<class S, class Thelp, class T>
    INLINE static void Step (S x, Thelp & help, T & values) {;}
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
    INLINE static void Step (S x, St t, Thelp & help, T & values)
    {
      DubinerJacobiPolynomialsScaledDiag<ALPHA0, BETA, DIAG, ORDER-1>::Step (x, t, help, values);
      typedef JacobiPolynomialFix<ALPHA0+2*(DIAG-ORDER), BETA> REC;
      
      if (ORDER == 0)
	help[DIAG-ORDER][0] = REC::P0(x);
      else if (ORDER == 1)
	{
	  help[DIAG-ORDER][0] = REC::P1(x,t);
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
    INLINE static void Step (S x, St t, Thelp & help, T & values) {;}
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
    INLINE static void Step (S x, Thelp & help, T & values)
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
	  REC::EvalNextTicTac (ORDER, x, help(DIAG-ORDER,0), help(DIAG-ORDER,1));
	}
      values[TrigIndex(ORDER,DIAG-ORDER,help.Height()-1)] = help(DIAG-ORDER,0);


      /*
	 .// hat aber nichts gebracht ...
      if (ORDER == 0)
	help(DIAG-ORDER,0) *= REC::P0(x);
      else if (ORDER == 1)
	{
	  // help(DIAG-ORDER,1) = help(DIAG-ORDER,0);
	  help(DIAG-ORDER,1) = REC::P1(x) / REC::P0(x) * help(DIAG-ORDER,0);
	}
      else
	{
          // Swap (help(DIAG-ORDER,0), help(DIAG-ORDER,1));
	  if ( ((DIAG-ORDER) % 2) == 0 )
	    REC::EvalNext (ORDER, x, help(DIAG-ORDER,0), help(DIAG-ORDER,1));
	  else
	    REC::EvalNext (ORDER, x, help(DIAG-ORDER,1), help(DIAG-ORDER,0));
	}
      if ( (DIAG-ORDER) % 2 )
	values[TrigIndex(ORDER,DIAG-ORDER,help.Height()-1)] = help(DIAG-ORDER,1);
      else
	values[TrigIndex(ORDER,DIAG-ORDER,help.Height()-1)] = help(DIAG-ORDER,0);
      */
    }
  };

  template <int ALPHA0, int BETA, int DIAG>
  class DubinerJacobiPolynomialsDiag_Linear<ALPHA0, BETA, DIAG, -1>
  {
  public:
    template<class S, class Thelp, class T>
    INLINE static void Step (S x, Thelp & help, T & values) {;}
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
    INLINE static void Step (S x, St t, Thelp & help, T & values)
    {
      DubinerJacobiPolynomialsScaledDiag_Linear<ALPHA0, BETA, DIAG, ORDER-1>::Step (x, t, help, values);
      typedef JacobiPolynomialFix<ALPHA0+2*(DIAG-ORDER), BETA> REC;
      
      if (ORDER == 0)
	help(DIAG-ORDER,0) *= REC::P0(x);
      else if (ORDER == 1)
	{
	  help(DIAG-ORDER,1) = help(DIAG-ORDER,0);
	  help(DIAG-ORDER,0) *= REC::P1(x,t) / REC::P0(x);
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
    INLINE static void Step (S x, St t, Thelp & help, T & values) {;}
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



  /*
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
  */




  class DubinerBasis
  {
  public:
    template <typename TI, class S, class T>
    INLINE static void Eval (TI n, S x, S y, T && values)
    {
      EvalMult (n, x, y, 1, values);
    }

    template <typename TI, class S, class Sc, class T>
    INLINE static void EvalMult (TI n, S x, S y, Sc c, T && values)
    {
      LegendrePolynomial leg;
      TI ii = 0;
      JacobiPolynomialAlpha jac(1);      
      leg.EvalScaledMult1Assign (n, y-(1-x-y), 1-x, c,
            SBLambda ([&] (TI i, S val) LAMBDA_INLINE 
                   {
                     jac.EvalMult (n-i, 2*x-1, val, values+ii);
                     ii += n-i+1;
                     jac.IncAlpha2();
                   }));
    }

    template <class S, class Sc, class T>
    INLINE static void EvalScaled (int n, S x, S y, Sc t, T && values)
    {
      EvalScaledMult (n, x, y, t, 1, values);
    }

    template <class S, class St, class Sc, class T>
    INLINE static void EvalScaledMult (int n, S x, S y, St t, Sc c, T && values)
    {
      LegendrePolynomial leg;
      int ii = 0;
      leg.EvalScaledMult1Assign (n, y-(t-x-y), t-x, c,
          SBLambda ([&] (int i, S val) LAMBDA_INLINE  // clang
                   {
                     JacobiPolynomialAlpha jac(1+2*i);
                     jac.EvalScaledMult1Assign (n-i, 2*x-t, t, val, values+ii);
                     ii += n-i+1;
                   }));
    }


    // evaluate basis functions of hightest order only
    template <typename TI, class S, class T>
    INLINE static void EvalHighestOrder (TI n, S x, S y, T && values)
    {
      EvalHighestOrderMult (n, x, y, 1, values);
    }

    template <typename TI, class S, class Sc, class T>
    INLINE static void EvalHighestOrderMult (TI n, S x, S y, Sc c, T && values)
    {
      LegendrePolynomial leg;
      TI ii = 0;
      JacobiPolynomialAlpha jac(1);      
      leg.EvalScaledMult1Assign (n, y-(1-x-y), 1-x, c,
            SBLambda ([&] (TI i, S val) LAMBDA_INLINE 
                   {
                     values[ii++] = jac.CalcHighestOrderMult(n-i, 2*x-1, val);
                     jac.IncAlpha2();
                   }));
    }

  };


  // orthogonal w.r.t. cubic bubble
  class DubinerBasisOrthoBub
  {
  public:
    template <typename TI, class S, class T>
    INLINE static void Eval (TI n, S x, S y, T && values)
    {
      EvalMult (n, x, y, 1, values);
    }

    template <typename TI, class S, class Sc, class T>
    INLINE static void EvalMult (TI n, S x, S y, Sc c, T && values)
    {
      JacobiPolynomialFix<1,1> leg;
      TI ii = 0;
      leg.EvalScaledMult1Assign (n, y-(1-x-y), 1-x, c,
            SBLambda ([&] (TI i, S val) LAMBDA_INLINE 
                   {
                     IntegratedJacobiPolynomialAlpha jac(4+2*i);                     
                     jac.EvalMult (n-i, 2*x-1, val, values+ii);
                     ii += n-i+1;
                   }));
    }

    template <class S, class Sc, class T>
    INLINE static void EvalScaled (int n, S x, S y, Sc t, T && values)
    {
      EvalScaledMult (n, x, y, t, 1, values);
    }

    template <class S, class St, class Sc, class T>
    INLINE static void EvalScaledMult (int n, S x, S y, St t, Sc c, T && values)
    {
      JacobiPolynomialFix<1,1> leg;
      int ii = 0;
      leg.EvalScaledMult1Assign (n, y-(t-x-y), t-x, c,
          SBLambda ([&] (int i, S val) LAMBDA_INLINE  // clang
                   {
                     IntegratedJacobiPolynomialAlpha jac(4+2*i);
                     jac.EvalScaledMult1Assign (n-i, 2*x-t, t, val, values+ii);
                     ii += n-i+1;
                   }));
    }

    
    template <class T>
    static size_t CalcNormInv (size_t n, T && norminv)
    {
      size_t ii = 0;
      for (size_t i = 0; i <= n; i++)
	for (size_t j = 0; j <= n-i; j++)
	  norminv[ii++] =
	    0.5 * (5+2*i+2*j)*(4+2*i+j)*(j+1) * (2*i+3)*(2*i+4) / (i+1);
      return ii;
    }


    /*
    // evaluate basis functions of hightest order only
    template <typename TI, class S, class T>
    INLINE static void EvalHighestOrder (TI n, S x, S y, T && values)
    {
      EvalHighestOrderMult (n, x, y, 1, values);
    }

    template <typename TI, class S, class Sc, class T>
    INLINE static void EvalHighestOrderMult (TI n, S x, S y, Sc c, T && values)
    {
      LegendrePolynomial leg;
      TI ii = 0;
      JacobiPolynomialAlpha jac(1);      
      leg.EvalScaledMult1Assign (n, y-(1-x-y), 1-x, c,
            SBLambda ([&] (TI i, S val) LAMBDA_INLINE 
                   {
                     values[ii++] = jac.CalcHighestOrderMult(n-i, 2*x-1, val);
                     jac.IncAlpha2();
                   }));
    }
    */
  };



  class DubinerBasis3D
  {
  public:
    template <typename TI, class S, class T>
    INLINE static void Eval (TI n, S x, S y, S z, T && values)
    {
      EvalMult (n, x, y, z, 1, values);
    }

    template <typename TI, class S, class Sc, class T>
    INLINE static void EvalMult (TI n, S x, S y, S z, Sc c, T && values)
    {
    size_t ii = 0;
    S lam4 = 1.0 - x-y-z;
    LegendrePolynomial leg;
    JacobiPolynomialAlpha jac1(1);    
    leg.EvalScaledMult1Assign 
      (n, z-lam4, z+lam4, c,
       SBLambda ([&](size_t k, S polz) LAMBDA_INLINE
                 {
                   // JacobiPolynomialAlpha jac(2*k+1);
                   JacobiPolynomialAlpha jac2(2*k+2);
 
                   jac1.EvalScaledMult1Assign
                     (n-k, y-z-lam4, 1-x, polz, 
                      SBLambda ([&] (size_t j, S polsy) LAMBDA_INLINE
                                {
                                  jac2.EvalMult1Assign(n-k-j, 2*x - 1, polsy, values+ii);
                                  ii += n-k-j+1;
                                  jac2.IncAlpha2();
                                }));
                   jac1.IncAlpha2();
                 }));
    }


    // evaluate basis functions of hightest order only
    template <typename TI, class S, class T>
    INLINE static void EvalHighestOrder (TI n, S x, S y, S z, T && values)
    {
      EvalHighestOrderMult (n, x, y, z, 1, values);
    }

    template <typename TI, class S, class Sc, class T>
    INLINE static void EvalHighestOrderMult (TI n, S x, S y, S z, Sc c, T && values)
    {
    size_t ii = 0;
    S lam4 = 1.0 - x-y-z;
    LegendrePolynomial leg;
    JacobiPolynomialAlpha jac1(1);    
    leg.EvalScaledMult1Assign 
      (n, z-lam4, z+lam4, c,
       SBLambda ([&](size_t k, S polz) LAMBDA_INLINE
                 {
                   JacobiPolynomialAlpha jac2(2*k+2);
 
                   jac1.EvalScaledMult1Assign
                     (n-k, y-z-lam4, 1-x, polz, 
                      SBLambda ([&] (size_t j, S polsy) LAMBDA_INLINE
                                {
                                  values[ii++] =
                                    jac2.CalcHighestOrderMult(n-k-j, 2*x - 1, polsy);
                                  jac2.IncAlpha2();
                                }));
                   jac1.IncAlpha2();
                 }));
    }
    
  };



  class DubinerBasis3DOrthoBub
  {
  public:
    template <typename TI, class S, class T>
    INLINE static void Eval (TI n, S x, S y, S z, T && values)
    {
      EvalMult (n, x, y, z, 1, values);
    }

    template <typename TI, class S, class Sc, class T>
    INLINE static void EvalMult (TI n, S x, S y, S z, Sc c, T && values)
    {
    size_t ii = 0;
    S lam4 = 1.0 - x-y-z;
    // JacobiPolynomialAlpha jac1(1);
    JacobiPolynomialFix<1,1> leg;
    leg.EvalScaledMult1Assign 
      (n, z-lam4, z+lam4, c,
       SBLambda ([&](size_t k, S polz) LAMBDA_INLINE
                 {
                   // JacobiPolynomialAlpha jac(2*k+1);
                   // JacobiPolynomialAlpha jac2(2*k+2);
		   
		   IntegratedJacobiPolynomialAlpha jac1(4+2*k);                     
                   jac1.EvalScaledMult1Assign
                     (n-k, y-z-lam4, 1-x, polz, 
                      SBLambda ([&] (size_t j, S polsy) LAMBDA_INLINE
                                {
				  IntegratedJacobiPolynomialAlpha jac2(6+2*k+2*j);                     
                                  jac2.EvalMult1Assign(n-k-j, 2*x - 1, polsy, values+ii);
                                  ii += n-k-j+1;
                                }));
                 }));
    }

    template <class T>
    static size_t CalcNormInv (size_t n, T && norminv)
    {
      size_t ii = 0;
      for (size_t i = 0; i <= n; i++)
	for (size_t j = 0; j <= n-i; j++)
	  for (size_t k = 0; k <= n-i-j; k++)
	    norminv[ii++] =
	      0.5 * (7+2*i+2*j+2*k)*(6+2*i+2*j+k)*(k+1) *
	      (5+2*i+2*j)*(4+2*i+j)*(j+1) *
	      (2*i+3)*(2*i+4) / (i+1);
      return ii;
    }
    /*
    // evaluate basis functions of hightest order only
    template <typename TI, class S, class T>
    INLINE static void EvalHighestOrder (TI n, S x, S y, S z, T && values)
    {
      EvalHighestOrderMult (n, x, y, z, 1, values);
    }

    template <typename TI, class S, class Sc, class T>
    INLINE static void EvalHighestOrderMult (TI n, S x, S y, S z, Sc c, T && values)
    {
    size_t ii = 0;
    S lam4 = 1.0 - x-y-z;
    LegendrePolynomial leg;
    JacobiPolynomialAlpha jac1(1);    
    leg.EvalScaledMult1Assign 
      (n, z-lam4, z+lam4, c,
       SBLambda ([&](size_t k, S polz) LAMBDA_INLINE
                 {
                   JacobiPolynomialAlpha jac2(2*k+2);
 
                   jac1.EvalScaledMult1Assign
                     (n-k, y-z-lam4, 1-x, polz, 
                      SBLambda ([&] (size_t j, S polsy) LAMBDA_INLINE
                                {
                                  values[ii++] =
                                    jac2.CalcHighestOrderMult(n-k-j, 2*x - 1, polsy);
                                  jac2.IncAlpha2();
                                }));
                   jac1.IncAlpha2();
                 }));
    }
    */
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
  INLINE void IntegratedLegendrePolynomial_Old (int n, S x, T && values)
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
  INLINE void LegendrePolynomialandDiff(int n, double x,  T & P, T & Px) 
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
    inline static int CalcScaled (int n, Sx x, Sy y, T && values)
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
    inline static int CalcScaled (Sx x, Sy y, T && values)
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
    inline static int CalcTrigExt (int n, Sx x, Sy y, T && values)
    {
      Sy fy = (1-y)*(1-y);
      Sx p3(0);
      Sx p2(-1);
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
    inline static int CalcTrigExtMult (int n, Sx x, Sy y, Sf fac, T && values)
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
    inline static int CalcTrigExtDeriv (int n, double x, double y, T && values)
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
      Sx p3(0.0);
      Sx p2(-1.0);
      Sx p1 = x;

      for (int j=2; j<=n; j++)
	{
	  p3=p2; p2=p1;
	  p1=( (2*j-3)/double(j) * x * p2 - (j-3)/double(j) * p3); //  / double(j);
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
    inline static int CalcDeriv (int n, double x, T && values)
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
    static INLINE void Do (T & inout)
    {
      inout[N-1] += double(N)/double(N+AL) * inout[N];
      inout[N] *= double(2*N+AL)/double(N+AL);
      TReduceAlpha<N-1,AL>::Do(inout);
    }

    template <class T>
    static INLINE void Trans (T & inout)
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
    static  INLINE void Do (T & inout)
    {
      if (HIGHEST)
	inout[N] = double(-N)/double(2*N+AL) * inout[N-1];
      else
	inout[N] += double(-N)/double(2*N+AL) * inout[N-1];
      inout[N-1] *= double(N+AL)/double(2*N+AL);
      TReduceAlphaFactor<N-1,AL,0>::Do(inout);
    }

    template <class T>
    static  INLINE void Trans (T & inout)
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
    static  INLINE void Do (T & inout)
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
    static  INLINE void Trans (T & inout)
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
    static  INLINE void Do (T & inout)
    {
      double val = inout(N);
      if (N > 1) inout(N-2) += val;
    
      inout(N) *= double (2*N-1);
      TDifferentiateLegendre<N-1>::Do(inout);

      inout(N-1) = inout(N);
      inout(N) = 0;
    }

    template <class T>
    static  INLINE void Trans (T & inout)
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
    static  INLINE void Do (T & inout)
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
    static  INLINE void Trans (T & inout)
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
    static  INLINE void Do (T & inout)
    {
      if (J == N)
	for (int i = 0; i <= N+1; i++)
	  inout(i, N+1-i) = 0;
      
      TTriangleReduceFactorCol<N,J,N-J,AL,BE> ::Do(inout);
      TTriangleReduceFactor<N,J-1,AL,BE>  ::Do(inout);
    }

    template <class T>
    static  INLINE void Trans (T & inout)
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
  static  INLINE void Do (T & inout)
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
  static  INLINE void Trans (T & inout)
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
    static  INLINE void Do (T & inout)
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
    static  INLINE void Trans (T & inout)
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
    static  INLINE void Do (T & inout)
    {
      TTriangleReduceLoop2New<N,I,N-I,AL,BE>::Do(inout);
      TTriangleReduceNew<N,I-1,AL,BE>::Do(inout);
    }

    template <class T>
    static  INLINE void Trans (T & inout)
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
  static  INLINE void Do (T & inout)
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
  static  INLINE void Trans (T & inout)
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


#include "recursive_pol_trig.hpp"
#include "recursive_pol_tet.hpp"

#endif
