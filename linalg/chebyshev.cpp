/**************************************************************************/
/* File:   chebyshev.cpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   May 2001                                                       */
/**************************************************************************/

/* 

Chebyshev iteration
  
*/ 

// #include <la.hpp>
#include "chebyshev.hpp"

namespace ngla
{
  using namespace ngla;

  ChebyshevIteration :: ChebyshevIteration
  (const BaseMatrix & aa, const BaseMatrix & ac, int asteps)
  {
    a = &aa;
    c = &ac;
    steps = asteps;
  }
  
  void ChebyshevIteration :: 
  SetBounds (double almin, double almax)
  {
    lmin = almin;
    lmax = almax;
  }
    
  void  ChebyshevIteration :: 
  Mult (const BaseVector & b, BaseVector & x) const
  {
    int m;
    double sigma, kappa;

    auto xold = b.CreateVector();
    auto xoold = b.CreateVector();
    auto w = b.CreateVector();
    auto r = b.CreateVector();

    xold = 0;

    x = (*c) * b;
    
    if (fabs (1-lmax) > 1e-7)
      {
	x *= 2 / (2-lmin-lmax);
	
	kappa = (1-lmin) / (1-lmax);
	sigma = 2;
	for (m = 1; m <= steps; m++)
	  {
	    //	cout << "it " << m << " ";
	    sigma = 4 / (4 - sqr( (1-1/kappa) / (1+1/kappa)) * sigma);
	    
	    // a->Residuum (x, b, r);
	    // c->Mult (r, w);
	    r = b - (*a) * x;
	    w = (*c) * r;

	    //	cout << "err = " << (r * w) << endl;

	    xoold = xold;
	    xold = x;
	    x += (2 / (2-lmin-lmax)) * w;
	    x *= sigma;
	    x += (1-sigma) * xoold;
	    /*	    
	    xoold.Set (1, xold);
	    xold.Set (1, x);
	    x.Add (2 / (2-lmin-lmax), w);
	    x *= sigma;
	    x.Add (1-sigma, xoold);
	    */
	  }
      }
  }
}
