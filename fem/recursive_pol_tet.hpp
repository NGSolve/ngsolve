#ifndef FILE_RECURSIVE_POL_TET
#define FILE_RECURSIVE_POL_TET

/*********************************************************************/
/* File:   recursive_pol.hpp                                         */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

#include "recursive_pol.hpp"

namespace ngfem
{

  class TetShapesInnerLegendre
  {
  public:
    template <typename Sx, typename Sy, typename Sz, typename T>
    static INLINE int Calc (int n, Sx x, Sy y, Sz z, T && values)
    {
      STACK_ARRAY(Sx, tmp, 3*n+2);
      Sx * polx = &tmp[0];
      Sx * poly = &tmp[n+1];
      // Sx * polz = &tmp[2*n+2];

      Sx bub = (1-x-y-z) * (1+x-y-z) * y * z ; 

      LegendrePolynomial::EvalScaledMult (n-4, x, 1-y-z, bub, polx);
      LegendrePolynomial::EvalScaled (n-4, 2*y - (1-z) , (1-z),poly); 
      
      // LegendrePolynomial (n-4, 2*z-1, polz);

      /*
      int ii = 0;
      for (int i = 0; i <= n-4; i++)
	for (int j = 0; j <= n-4-i; j++)
          {
            Sx hp = polx[i]*poly[j];
            for (int k = 0; k <= n-4-i-j; k++)
              values[ii++] = hp * polz[k];
          }
      */
      int ii = 0;
      for (int i = 0; i <= n-4; i++)
	for (int j = 0; j <= n-4-i; j++)
	  {
	    LegendrePolynomial::EvalMult (n-4-i-j, 2*z-1, polx[i]*poly[j], values+ii);
	    ii += n-3-i-j;
	  }
      return ii;
    }

    template <typename Sx, typename Sy, typename Sz, typename T1, typename T2, typename T3>
    static INLINE void CalcSplitted (int n, Sx x, Sy y, Sz z, T1 & val1, T2 & val2, T3 & val3)
    {
      Sx bub1 = sqr (1-y-z) - sqr(x);
      // ScaledLegendrePolynomialMult (n-4, x, 1-y-z, bub1, val1);
      // ScaledLegendrePolynomialMult (n-4, 2*y - (1-z) , (1-z), y, val2); 
      LegendrePolynomial::EvalScaledMult (n-4, x, 1-y-z, bub1, val1);
      LegendrePolynomial::EvalScaledMult (n-4, 2*y - (1-z) , (1-z), y, val2); 

      LegendrePolynomial leg;
      leg.EvalMult (n-4, 2*z-1, z, val3);
    }
  };



  class TetShapesInnerJacobi
  {
  public:
    template <typename Sx, typename Sy, typename Sz, typename T>
    static int Calc (int n, Sx x, Sy y, Sz z, T & values)
    {
      int ii = 0, j, k;
      ArrayMem<Sx, 20> polx(n+1), poly(n+1), polz(n+1);
      
      Sx bub = y * z * (1-x-y-z) * (1+x-y-z);
      ScaledJacobiPolynomial (n-4, x, (1-y-z), 2, 2, polx);
    
      for (int ix = 0; ix <= n-4; ix++)
	{
	  ScaledJacobiPolynomial (n-4, (2*y-1+z),(1-z), 2*ix+5, 2, poly);
	  JacobiPolynomial (n-4, 2*z-1, 2*ix+5, 2, polz);
	
	  for (j = 0; j <= n-4-ix; j++)
	    for (k = 0; k <= n-4-ix-j; k++)
	      values[ii++] = bub * polx[ix] * poly[j] * polz[k];
	}
      return ii;
    }
  };

  class TetShapesFaceLegendre
  {
  public:
    template <typename Sx, typename Sy, typename Sz, typename T>
    static int Calc (int n, Sx x, Sy y, Sz z, T && values)
    {
      ArrayMem<Sx, 20> polx(n+1), poly(n+1);
    
      ScaledLegendrePolynomial (n-3, x, 1-y-z, polx);
      ScaledLegendrePolynomial (n-3, 2*y-(1-z),(1-z), poly);
      Sx bub = (1-x-y-z) * (1+x-y-z)*y; 

      int ii = 0;
      for (int i = 0; i <= n-3; i++)
	for (int j = 0; j <= n-3-i; j++)
	  values[ii++] = bub * polx[i] * poly[j];
 
      return ii;
    }
    template <typename Sx, typename Sy, typename Sz, typename T>
    static INLINE void CalcSplitted (int n, Sx x, Sy y, Sz z, T && val1, T && val2)
    {
      Sx bub1 = sqr (1-y-z) - sqr (x); 
      /*
      ScaledLegendrePolynomialMult (n-3, x, 1-y-z, bub1, val1);
      ScaledLegendrePolynomialMult (n-3, 2*y-(1-z),(1-z), y, val2);
      */
      LegendrePolynomial::EvalScaledMult (n-3, x, 1-y-z, bub1, val1);
      LegendrePolynomial::EvalScaledMult (n-3, 2*y-(1-z),(1-z), y, val2);
    }
  };




  class TetShapesFaceJacobi
  {
  public:
    template <typename Sx, typename Sy, typename Sz, typename T>
    static int Calc (int n, Sx x, Sy y, Sz z, T && values)
    {
      int ii = 0;
      ArrayMem<Sx, 20> polx(n+1), poly(n+1);

      Sx bub = y * (1-x-y-z) * (1+x-y-z);
      ScaledJacobiPolynomial (n-3, x, 1-y-z, 2, 2, polx);
    
      for (int ix = 0; ix <= n-3; ix++)
	{
	  ScaledJacobiPolynomial (n-3, 2*y-1+z, 1-z, 2*ix+5, 2, poly);
	  for (int j = 0; j <= n-3-ix; j++)
	    values[ii++] = bub * polx[ix] * poly[j];
	}
      return ii;
    }
  };


#ifdef NOTAVAILABLE
  class TetShapesFaceOpt1
  {
  public:
    template <typename Sx, typename Sy, typename Sz, typename T>
    static int Calc (int n, Sx x, Sy y, Sz z, T && values)
    {
      int nd = Calc1 (n, x, y, z, values);
      ArrayMem<Sx, 100> hvalues(nd);

      Sx lam1 = 0.5 * (1 + x - y - z);
      Sx lam2 = 0.5 * (1 - x - y - z);
      Sx lam3 = y;
      Sx lam4 = z;

      Sx hlam1 = 0;
      Sx hlam2 = lam2;
      Sx hlam3 = lam3;
      Sx hlam4 = lam4 + lam1;

      Sx hx = hlam1 - hlam2;
      Sx hy = hlam3;
      Sx hz = hlam4;

      Sx frac;
      if (hlam4 < 1e-12)
	frac = 0.0;
      else
	frac = lam4 / hlam4;

      Calc1 (n, hx, hy, hz, hvalues);
      for (int i = 0; i < nd; i++)
	values[i] -= frac * hvalues[i];

      return nd;
    }
  
    template <typename Sx, typename Sy, typename Sz, typename T>
    static int Calc1 (int n, Sx x, Sy y, Sz z, T && values)
    {
      int nd = Calc2 (n, x, y, z, values);
      ArrayMem<Sx, 100> hvalues(nd);

      Sx lam1 = 0.5 * (1 + x - y - z);
      Sx lam2 = 0.5 * (1 - x - y - z);
      Sx lam3 = y;
      Sx lam4 = z;

      Sx hlam1 = lam1;
      Sx hlam2 = 0;
      Sx hlam3 = lam3;
      Sx hlam4 = lam4 + lam2;

      Sx hx = hlam1 - hlam2;
      Sx hy = hlam3;
      Sx hz = hlam4;

      Sx frac;
      if (hlam4 < 1e-12)
	frac = 0.0;
      else
	frac = lam4 / hlam4;

      Calc2 (n, hx, hy, hz, hvalues);
      for (int i = 0; i < nd; i++)
	values[i] -= frac * hvalues[i];
    
      return nd;
    }

    template <typename Sx, typename Sy, typename Sz, typename T>
    static int Calc2 (int n, Sx x, Sy y, Sz z, T && values)
    {
      Sx * hp = &values[0];
      int nd = Calc3 (n, x, y, z, hp);
      ArrayMem<Sx, 100> hvalues(nd);

      Sx lam1 = 0.5 * (1 + x - y - z);
      Sx lam2 = 0.5 * (1 - x - y - z);
      Sx lam3 = y;
      Sx lam4 = z;

      Sx hlam1 = lam1;
      Sx hlam2 = lam2;
      Sx hlam3 = 0;
      Sx hlam4 = lam4 + lam3;

      Sx hx = hlam1 - hlam2;
      Sx hy = hlam3;
      Sx hz = hlam4;

      Sx frac;
      if (hlam4 < 1e-12)
	frac = 0.0;
      else
	frac = lam4 / hlam4;

      hp = &hvalues[0];
      Calc3 (n, hx, hy, hz, hp);
      for (int i = 0; i < nd; i++)
	values[i] -= frac * hvalues[i];
    
      return nd;
    }






    template <typename Sx, typename Sy, typename Sz, typename T>
    static int Calc3 (int n, Sx x, Sy y, Sz z, T && values)
    {
      int ii = 0;
      ArrayMem<Sx, 20> polx(n+1), poly(n+1);

      const IntegrationRule & rule = SelectIntegrationRule (ET_TRIG, n+2);

      for (int ix = 0; ix <= n-3; ix++)
	for (int j = 0; j <= n-3-ix; j++)
	  values[ii++] = 0;
      for (int i = 0; i < rule.GetNIP(); i++)
	{
	  ii = 0;
	  const IntegrationPoint & ip = rule[i];

	  Sx hx = x + z * (-1+2*ip(0)+ip(1));
	  Sy hy = y + z * ip(1);

	  Sx bub = hy * (1-hx-hy) * (1+hx-hy);
	  ScaledJacobiPolynomial (n-3, hx, 1-hy, 2, 2, polx);

	  Sx fac = 2 * bub * ip.Weight(); //  / (z*z);

	  for (int ix = 0; ix <= n-3; ix++)
	    {
	      ScaledJacobiPolynomial (n-3, 2*hy-1, 1, 2*ix+5, 2, poly);
	      for (int j = 0; j <= n-3-ix; j++)
		values[ii++] += fac * polx[ix] * poly[j];
	    }
	}
      return ii;
    }
  };



  class TetShapesFaceOpt2
  {
  public:
    template <typename Sx, typename Sy, typename Sz, typename T>
    static int Calc (int n, Sx x, Sy y, Sz z, T & values)
    {
      int ii = 0, i, j;
      ArrayMem<Sx, 20> polx(n+1), poly(n+1);

      const IntegrationRule & rule = SelectIntegrationRule (ET_TRIG, n+2);

      for (int ix = 0; ix <= n-3; ix++)
	for (j = 0; j <= n-3-ix; j++)
	  values[ii++] = 0;
      for (i = 0; i < rule.GetNIP(); i++)
	{
	  ii = 0;
	  const IntegrationPoint & ip = rule[i];

	  Sx hx = x + z * (-1+2*ip(0)+ip(1));
	  Sy hy = y + z * ip(1);
          
	  //Sx bub = hy * (1-hx-hy) * (1+hx-hy);
	  ScaledJacobiPolynomial (n-3, hx, 1-hy, 2, 2, polx);

	  Sx fac = 2  * ip.Weight(); //  / (z*z);

	  for (int ix = 0; ix <= n-3; ix++)
	    {
	      ScaledJacobiPolynomial (n-3, 2*hy-1, 1, 2*ix+5, 2, poly);
	      for (j = 0; j <= n-3-ix; j++)
		values[ii++] += fac * polx[ix] * poly[j];
	    }
	}
      int jj = ii;
      ArrayMem<Sx, 100> hvalues(jj);
      for (int i = 0; i < ii; i++)
	hvalues[i] = values[i];
      ii = 0; jj = 0;
      for (int k = 0; k <= n-3; k++)
	{
	  for (int m = 0; m <= n-3-k; m++)
	    {
	      values[ii++] = y*(1-x-y-z)*(1+x-y-z)*hvalues[jj++];
	    }
	}

      return ii;
    }
  };

#endif
  








  /*
 /// compute shape
 virtual void CalcShape (const IntegrationPoint & ip,
 FlatVector<> shape) const
 {
 double lam1 = ip(0);
 double lam2 = 1-ip(0)-ip(1);
 double lam3 = ip(1);
 double bub = lam1*lam2*lam3;

 ArrayMem<double, 50> polx(order-2);
 ArrayMem<double, 50> poly(order-2);

 double fac = 1;
 for (int i = 0; i <= order-3; i++)
 {
 polx[i] *= fac;
 fac *= (1-lam3);
 }

 int ii = 0;
 for (int ix = 0; ix <= order-3; ix++)
 {
 JacobiPolynomial (order-3, lam3-(lam1+lam2), 2*ix+5, 2, poly);
 for (int iy = 0; iy <= order-3-ix; iy++)
 shape(ii++) = bub*polx[ix] * poly[iy];
 }
 }

  */


}


#endif
