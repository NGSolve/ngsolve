#ifndef FILE_RECURSIVE_POL_TRIG
#define FILE_RECURSIVE_POL_TRIG


#include "recursive_pol.hpp"

namespace ngfem
{

  /**
     Computes face shape functions for triangles.
   
  */
  class TrigShapesInnerLegendre
  {
  public:
    /// compute shape functions in factored form $\varphi_{ij} = u_i v_j$
    template <typename Sx, typename Sy, typename T>
    INLINE static void CalcSplitted (int n, Sx x, Sy y, T & val1, T & val2)
    {
      LegendrePolynomial leg;

      Sx bub1 = (1-x-y)*(1+x-y);
      leg.EvalScaledMult (n-3, x, 1-y, bub1, val1); 
      leg.EvalMult (n-3, 2*y-1, y, val2);
    }


    template <int n, typename Sx, typename Sy, typename T>
    INLINE static void CalcSplitted (Sx x, Sy y, T & val1, T & val2)
    {
      Sx bub1 = (1-x-y)*(1+x-y);
      // ScaledLegendrePolynomialMult (n-3, x, 1-y, bub1, val1); 
      // LegendrePolynomialFO<n-3>::EvalScaledMult (x, 1-y, bub1, val1); 
      LegendrePolynomial leg;
      leg.EvalScaledMultFO<n-3> (x, 1-y, bub1, val1);

      // LegendrePolynomialMult (n-3, 2*y-1, y, val2);
      // LegendrePolynomialFO<n-3>::EvalMult (2*y-1, y, val2);
      leg.EvalMultFO<n-3> (2*y-1, y, val2);
    }


    /// computes all shape functions
    template < typename Sx, typename Sy, typename T>
    static int Calc (int n, Sx x, Sy y, T & values)
    {
      ArrayMem<Sx, 20> polx(n-2), poly(n-2);

      ScaledLegendrePolynomial (n-3, x, 1-y, polx);
      LegendrePolynomial (n-3, 2*y-1, poly);
      Sx bub = y * (1-x-y) * (1+x-y);

      int ii = 0;
      for (int i = 0; i <= n-3; i++)
	for (int j = 0; j <= n-3-i; j++)
	  values[ii++] = bub * polx[i] * poly[j];

      return ii;
    }

    template <int n, typename Sx, typename Sy, typename T>
    static int Calc (Sx x, Sy y, T & values)
    {
      // ArrayMem<Sx, 20> polx(n-2), poly(n-2);
      Sx polx[n], poly[n];

      /*
	ScaledLegendrePolynomial (n-3, x, 1-y, polx);
	LegendrePolynomial (n-3, 2*y-1, poly);
	Sx bub = y * (1-x-y) * (1+x-y);
      */
      CalcSplitted<n> (x, y, polx, poly);
      int ii = 0;
      for (int i = 0; i <= n-3; i++)
	for (int j = 0; j <= n-3-i; j++)
	  values[ii++] = /* bub * */ polx[i] * poly[j];

      return ii;
    }





    /// computes all shape functions
    template <typename Sx, typename Sy, typename Sf, typename T>
    static int CalcMult (int n, Sx x, Sy y, Sf & fac, T & values)
    {
      ArrayMem<Sx, 20> polx(n-2), poly(n-2);

      ScaledLegendrePolynomial (n-3, x, 1-y, polx);
      LegendrePolynomial (n-3, 2*y-1, poly);

      Sx bub = fac * y * (1-x-y) * (1+x-y);

      int ii = 0;
      for (int i = 0; i <= n-3; i++)
	for (int j = 0; j <= n-3-i; j++)
	  values[ii++] = bub * polx[i] * poly[j];

      return ii;
    }


  };


  /**
     Compute triangular face shape functions.
     Shape functions are $L_2$-orthogonal Jacobi polynomials
  */
  class TrigShapesInnerJacobi
  {
  public:
    /// computes all base functions
    template <typename Sx, typename Sy, typename T>
    static int Calc (int n, Sx x, Sy y, T & values)
    {
      int ii = 0;
      ArrayMem<Sx, 20> polx(n+1), poly(n+1);

      Sx bub = y * (1-x-y) * (1+x-y);
      ScaledJacobiPolynomial (n-3, x, 1-y, 2, 2, polx);

      for (int ix = 0; ix <= n-3; ix++)
	{
	  JacobiPolynomial (n-3, 2*y-1, 2*ix+5, 2, poly);
	  for (int j = 0; j <= n-3-ix; j++)
	    values[ii++] = bub * polx[ix] * poly[j];
	}
      return ii;
    }


  };


  //*************************************MONOMIAL EXTENSION***********************************************************

    /**
       Compute triangle edge-shape functions.

       functions vanish on upper two edges

       x,y: coordinates in triangle (-1, 0), (1, 0), (0, 1)

       f_i (x, 0) = IntegratedLegendrePol_i (x)

       f_i ... pol of order i

       Monomial extension:
    */
  class TrigExtensionMonomial
  {
  public:
    /// computes function on triangle
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
	  p1=( (2*j-3) * x * p2 - (j-3) * fy * p3) / j;
	  values[j-2] = p1;
	}
      return n-1;
    }


    /// computes derivates on triangle, values must be $N \times 2$ matrix
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

    /// computes values on edge
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

    /// computes derivatives on edge
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


  // ***********************  OPTIMAL EXTENSION  *****************************************
  /**
     Evaluete optimal edge-shape function.
     on triangle (with coordinates (-1,0),(1,0),(0,1))
     $$
     \varphi_i:=\frac{1}{y}\int_{x-y}^{x+y}C_{i-1}^{-1/2}(s)ds
     $$
     vanishes on other both edges
     Computed by recurrence
     $$
     {\varphi}_i:=
     p_n(x,y) = a_n p_{n-4}} + b_n x p_{n-3} + (c_n + d_n (x*x-y*y)) p_{n-2}} + e_n x p_{n-1}}
     $$
  */
  class TrigExtensionOptimal
  {
    enum { SIZE = 1000 };
    static double coefs[SIZE][5];
    static bool initialized;
  public:

    TrigExtensionOptimal ();

    template <class Sx, class Sy, class T>
    inline static int CalcTrigExt (int p, Sx x, Sy y, T & values)
    {
      Sx p1(0.0), p2(0.0), p3(0.0), p4(0.0), p5(0.0);
      Sx bub = (1.+x-y)*(x+y-1);
      switch(p)
	{
	default:
	case 5:

	  {
	    p5 = -1./24. * bub*x * (-9.+21.*x*x-14.*y+35.*y*y);
	    values[3] = p5;
	  }

	case 4:

	  {
	    p4 = -0.125 * bub * (-1.+5.*x*x-2.*y+3.*y*y);
	    values[2] = p4;
	  }

	case 3:

	  {
	    p3 = -0.5 * bub * x;
	    values[1] = p3;
	  }

	case 2:

	  {
	    p2 = -0.5 * bub;
	    values[0] = p2;
	  }
	case 1:
	case 0:
	  ;
	}

      for (int j = 6; j <= p; j++)
	{
	  p1 = 
	    -coefs[j-6][0] * p2 
	    + coefs[j-6][1] * x * p3 
	    - (coefs[j-6][2]+coefs[j-6][3]*(x*x-y*y)) * p4 
	    + coefs[j-6][4] * x * p5;
	  p2 = p3; p3 = p4; p4 = p5; p5 = p1;
	  values[j-2] = p1;
	}

      return p-1;
    }

    template <class T>
    inline static int CalcTrigExtDeriv (int n, double x, double y, T & values)
    {
      ArrayMem<AutoDiff<2>,10> ad_values(n-1);
      AutoDiff<2> ad_x(x, 0);
      AutoDiff<2> ad_y(y, 1);

      CalcTrigExt (n, ad_x, ad_y, ad_values);

      for (int i = 0; i < n-1; i++)
	for (int j = 0; j < 2; j++)
	  values(i,j) = ad_values[i].DValue(j);
      return n-1;
    }


    template <class Sx, class T>
    inline static int Calc (int n, Sx x, T & values)
    {
      Sx y = 0.0;
      return CalcTrigExt (n, x, y, values);
    }

    template <class T>
    inline static int CalcDeriv (int n, double x, T & values)
    {
      return CalcTrigExtDeriv (n, x, 0.0, values);
    }
  };



  // ************************ MINIMAL EXTENSION *************************


  /**
     .......
  */
  class TrigExtensionMin
  {
  public:
    template <class Sx, class Sy, class T>
    inline static int CalcTrigExt (int n, Sx x, Sy y, T & values)
    {
#ifdef ABC
      //    TrigExtensionOptimal::CalcTrigExt (n, x, y, values);

      ArrayMem<Sx, 20> polx(n+1);

      const IntegrationRule & rule =
	GetIntegrationRules().SelectIntegrationRule (ET_SEGM, n+2);

      Sx bub = (1-x-y) * (1+x-y);
      for (int i = 0; i <= n-2; i++)
	values[i] = 0;

      for (int i = 0; i < rule.GetNIP(); i++)
	{
	  const IntegrationPoint & ip = rule.GetIP(i);

	  Sx hx = x + y * (2.0*ip(0)-1.0);
			 
	  JacobiPolynomial (n-2, hx, 1, 1, polx);
	  // LegendrePolynomial (n-2, hx, polx);

	  Sx fac = bub * ip.Weight(); //  / (z*z);

	  for (int j = 0; j <= n-2; j++)
	    values[j] += fac * polx[j];
	}



      // IntegratedLegendreMonomialExt::CalcTrigExt (n, x, y, values);

      /*
	Sy fy = (1-y)*(1-y);
	Sx p3 = 0;
	Sx p2 = -1;
	Sx p1 = x;

	for (int j=2; j<=n; j++)
	{
        p3=p2; p2=p1;
        p1=( (2*j-3) * x * p2 - (j-3) * fy * p3) / j;
	values[j-2] = p1;
	}
      */

      /*
	Sx fac = 1;
	for (int j = n-2; j >= 0; j--)
	{
        values[j] = fac * values[j];
        fac = fac * (1-y);
	}
      */

      for (int j = 0; j < n-3; j++)
	values[j] *= VertexExtensionOptimal<2>::Calc (n-2-j, 1-y);

      return n-1;
#endif

      if (y == 0.0) y += 1e-8;
      ArrayMem<Sx, 20> pol1(n+2), pol2(n+2);
      JacobiPolynomial (n+1, x-y, 1, 1, pol1);
      JacobiPolynomial (n+1, x+y, 1, 1, pol2);
      for (int i = 0; i <= n-2; i++)
	values[i] = (pol2[i+1]-pol1[i+1]) / (2*y);

      /*
	ScaledJacobiPolynomial (n-2, x, 1-y, 2, 2, values);
	for (int i = 0; i < n-1; i++)
	values[i] *= 0.5 * (i+4) * (1-x*x);
      */



      LowEnergyVertexPolynomials2D (n, 1-2*y, pol1);
      for (int i = 0; i <= n-2; i++)
	values[i] *= pol1[n-2-i] * (1-x-y) * (1+x-y);

      return n-1;
    }

    template <class T>
    inline static int CalcTrigExtDeriv (int n, double x, double y, T & values)
    {
      ArrayMem<AutoDiff<2>,10> ad_values(n-1);
      AutoDiff<2> ad_x(x, 0);
      AutoDiff<2> ad_y(y, 1);

      CalcTrigExt (n, ad_x, ad_y, ad_values);

      for (int i = 0; i < n-1; i++)
	for (int j = 0; j < 2; j++)
	  values(i,j) = ad_values[i].DValue(j);
      return n-1;
    }



    template <class Sx, class T>
    inline static int Calc (int n, Sx x, T & values)
    {
      /*
	Sx y = 0.0;
	CalcTrigExt (n, x, y, values);
      */

      JacobiPolynomial (n-2, x, 2, 2, values);
      for (int i = 0; i < n-1; i++)
	values[i] *= 0.5 * (i+4) * (1-x*x);

      return n-1;

    }

    template <class T>
    inline static int CalcDeriv (int n, double x, T & values)
    {
      ArrayMem<AutoDiff<1>,10> ad_values(n-1);
      AutoDiff<1> ad_x(x, 0);

      Calc (n, ad_x, ad_values);

      for (int i = 0; i < n-1; i++)
	values[i] = ad_values[i].DValue(0);
      return n-1;

      // return CalcTrigExtDeriv (n, x, 0.0, values);
    }
  };


}

#endif
