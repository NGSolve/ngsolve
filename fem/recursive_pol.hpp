#ifndef FILE_RECURSIVE_POL
#define FILE_RECURSIVE_POL

/*********************************************************************/
/* File:   recursive_pol.hpp                                         */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

/*
  Recursive Polynomials
*/



/**
   Legendre polynomials P_i on (-1,1), 
   up to order n --> n+1 values
   
   l P_l = (2 l - 1) x P_{l-1} - (l-1) P_{l-2}

   P_0 = 1
   P_1 = x
   P_2 = -1/2 + 3/2 x^2
*/
template <class S, class T>
inline void LegendrePolynomial (int n, S x, T & values)
{
  S p1, p2;

  if (n < 0) return;
  values[0] = 1.0;
  if (n < 1) return;
  values[1] = x;
  if (n < 2) return;
  values[2] = p2 = 1.5 * sqr(x) - 0.5;
  if (n < 3) return;
  values[3] = p1 = ( (5.0/3.0) * p2 - (2.0/3.0)) * x;
  if (n < 4) return;

  for (int j=4; j < n; j+=2)
    {
      double invj = 1.0 / j;
      p2 *= (invj-1);
      p2 += (2-invj) * x * p1;
      values[j] = p2; 

      invj = 1.0 / (j+1);
      p1 *= (invj-1);
      p1 += (2-invj) * x * p2;
      values[j+1] = p1; 
    }

  if (n % 2 == 0)
    {
      double invn = 1.0 / n;
      values[n] = (2-invn)*x*p1 - (1-invn)*p2;
    }

  /*
    S p1 = 1.0, p2 = 0.0, p3;

    if (n >= 0)
    values[0] = 1.0;

    for (int j=1; j<=n; j++)
    {
    p3 = p2; p2 = p1;
    p1 = ((2.0*j-1.0)*x*p2 - (j-1.0)*p3) / j;
    values[j] = p1;
    }
  */
}





/**
   c P_i 
*/
template <class S, class Sc, class T>
inline void LegendrePolynomialMult (int n, S x, Sc c , T & values)
{
  S p1, p2;

  if (n < 0) return;
  values[0] = c;
  if (n < 1) return;
  values[1] = c * x;
  if (n < 2) return;
  values[2] = p2 = c * (1.5 * sqr(x) - 0.5);
  if (n < 3) return;
  values[3] = p1 = ( (5.0/3.0) * p2 - (2.0/3.0) * c) * x;
  if (n < 4) return;

  for (int j=4; j < n; j+=2)
    {
      double invj = 1.0 / j;
      p2 *= (invj-1);
      p2 += (2-invj) * x * p1;
      values[j] = p2; 

      invj = 1.0 / (j+1);
      p1 *= (invj-1);
      p1 += (2-invj) * x * p2;
      values[j+1] = p1; 
    }

  if (n % 2 == 0)
    {
      double invn = 1.0 / n;
      values[n] = (2-invn)*x*p1 - (1-invn)*p2;
    }

  /*
    S p1 = 1.0, p2 = 0.0, p3;

    if (n >= 0)
    values[0] = 1.0;

    for (int j=1; j<=n; j++)
    {
    p3 = p2; p2 = p1;
    p1 = ((2.0*j-1.0)*x*p2 - (j-1.0)*p3) / j;
    values[j] = p1;
    }
  */
}











  template <int n>
  class LegendrePolynomialFO
  {
  public:
    template <class S, class T>
    static void Eval (S x, T & values)
    {
      LegendrePolynomialFO<n-1>::Eval (x, values);
      values[n] = (2.0*n-1)/n * x * values[n-1] - (n-1.0)/n * values[n-2];    
    }

    template <class S, class Sc, class T>
    static void EvalMult (S x, Sc c, T & values)
    {
      LegendrePolynomialFO<n-1>::EvalMult (x, c, values);
      values[n] = (2.0*n-1)/n * x * values[n-1] - (n-1.0)/n * values[n-2];    
    }
  };


  template <> class LegendrePolynomialFO<-1>
  {
  public:
    template <class S, class T>
    static void Eval (S x, T & values)
    { ; }

    template <class S, class Sc, class T>
    static void EvalMult (S x, Sc c, T & values)
    { ; }

  };

  template <> class LegendrePolynomialFO<0>
  {
  public:
    template <class S, class T>
    static void Eval (S x, T & values)
    {
      values[0] = 1;
    }

    template <class S, class Sc, class T>
    static void EvalMult (S x, Sc c, T & values)
    {
      values[0] = c;
    }

  };

  template <> class LegendrePolynomialFO<1>
  {
  public:
    template <class S, class T>
    static void Eval (S x, T & values)
    {
      values[0] = 1;
      values[1] = x;
    }

    template <class S, class Sc, class T>
    static void EvalMult (S x, Sc c, T & values)
    {
      values[0] = c;
      values[1] = c*x;
    }
  };
















template <class S, class T>
inline void GegenbauerPolynomial (int n, S x, double lam, T & values)
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
inline void IntegratedLegendrePolynomial (int n, S x, T & values)
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
   Derived Legendre polynomials on (-1,1)
   value[i] (x) = d/dx P_{i+1}
*/
template <class S, class T>
inline void DerivedLegendrePolynomial (int n, S x, T & values)
{
  GegenbauerPolynomial<S,T> (n, x, 1.5, values);
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
inline void HermitePolynomial (int n, S x, T & values)
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
inline void ScaledLegendrePolynomial (int n, Sx x, St t, T & values)
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
inline void ScaledLegendrePolynomialMult (int n, Sx x, St t, Sc c, T & values)
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
inline void DiffScaledLegendrePolynomial (int n, double x, double t, T & values)
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
      p1=((2.0*j-3.0) * x*p2 - t*t*(j-3.0)*p3)/j;
      values[j] = p1;
    }
}



template <class S, class T>
inline void JacobiPolynomial (int n, S x, double alpha, double beta, T & values)
{
  S p1 = 1.0, p2 = 0.0, p3;

  if (n >= 0) 
    p2 = values[0] = 1.0;
  if (n >= 1) 
    p1 = values[1] = 0.5 * (2*(alpha+1)+(alpha+beta+2)*(x-1));

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
inline void JacobiPolynomialMult (int n, S x, double alpha, double beta, Sc c, T & values)
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
inline void ScaledJacobiPolynomial (int n, S x, St t, double alpha, double beta, T & values)
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
inline void ScaledJacobiPolynomialMult (int n, S x, St t, double alpha, double beta, Sc c, T & values)
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







class ConvertJacobi
{
  typedef double d2[2];
  static Array<d2*> coefs_increasealpha;
  static Array<d2*> coefs_reducealpha;
  static Array<d2*> coefs_reducealphafac;

public:
  ConvertJacobi ();

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














#endif
