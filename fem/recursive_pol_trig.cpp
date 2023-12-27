/*********************************************************************/
/* File:   recurisve_pol_trig.cpp                                    */
/* Author: Almedin Becirovic                                         */
/* Date:   14. Apr. 2003                                             */
/*********************************************************************/






// #include <fem.hpp>
#include "recursive_pol_trig.hpp"

namespace ngfem
{
  using namespace ngfem;


  double TrigExtensionOptimal :: coefs[SIZE][5];
  bool TrigExtensionOptimal :: initialized = 0;


  TrigExtensionOptimal :: TrigExtensionOptimal ()
  {
    //    cout<<"TrigEXOpt"<<endl;
    if (!initialized)
      {
	// cout << "calculate TrigExtOpt coefficeints " << endl;
	int i = 0;
	for (double n=6; n< SIZE; n++)
         {
           double help1 = n*(n+1)*(-7+2*n);
	   double help2 = -5+2*n;

           coefs[i][0] = (-6+n)*(-5+n)*(-3+2*n)/help1;
           coefs[i][1] = 2*(-5+n)*(-7+2*n)*(-3+2*n)/help1;

	   coefs[i][2] = help2*(3-10*n+2*n*n)/help1;
           coefs[i][3] = help2*(21-20*n+4*n*n)/help1;

	   coefs[i][4] = 2*n*(-7+2*n)*(-3+2*n)/help1;

	   i++;
        }
	initialized = 1;
      }
  }


  TrigExtensionOptimal dummy_textopt;

 


  /*
  template <typename Sx, typename Sy, typename Sz, typename T>
  int TetShapesFaceOpt1 :: Calc3 (int n, Sx x, Sy y, Sz z, T & values)
  {
    int ii = 0, i, j, k;
    ArrayMem<Sx, 20> polx(n+1), poly(n+1);

    const IntegrationRule & rule = 
      GetIntegrationRules().SelectIntegrationRule (ET_TRIG, n);
    
    for (int ix = 0; ix <= n-3; ix++)
      for (j = 0; j <= n-3-ix; j++)
	values[ii++] = 0;
    for (i = 0; i < rule.GetNIP(); i++)
      {
	ii = 0;
	const IntegrationPoint & ip = rule.GetIP(i);
	
	Sx hx = x + z * (-1+2*ip(0)+ip(1));
	Sy hy = y + z * ip(1);
	
	Sx bub = hy * (1-hx-hy) * (1+hx-hy);
	ScaledJacobiPolynomial (n-3, x, 1-y, 2, 2, polx);

	Sx fac = 2 * bub * ip.Weight(); //  / (z*z);
	
	for (int ix = 0; ix <= n-3; ix++)
	  {
	    ScaledJacobiPolynomial (n-3, 2*y-1, 1, 2*ix+5, 2, poly);
	    for (j = 0; j <= n-3-ix; j++)
	      values[ii++] += fac * polx[ix] * poly[j];
	  }
      }
    return ii;
  }


  template <>
  int TetShapesFaceOpt1 :: Calc3<double, double, double, double*> (int n, double x, double y, double z, double *& values);

  template <>
  int TetShapesFaceOpt1 :: Calc3<AutoDiff<3>, AutoDiff<3>, AutoDiff<3>, AutoDiff<3>*> (int n, AutoDiff<3> x, AutoDiff<3> y, AutoDiff<3> z,  AutoDiff<3> *& values);
  */

}
