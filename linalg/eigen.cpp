/**************************************************************************/
/* File:   cg.cc                                                          */
/* Author: Joachim Schoeberl                                              */
/* Date:   5. Jul. 96                                                     */
/**************************************************************************/

/* 

Conjugate Gradient Soler
  
*/ 

#include "eigen.hpp"

namespace ngla
{
  using namespace ngla;

  EigenSystem :: EigenSystem (const BaseMatrix & aa)
  {
    SetMatrix (aa);
    c = NULL;
    SetPrecision (1e-10);
    SetMaxSteps (200);
  }



  EigenSystem :: EigenSystem (const BaseMatrix & aa,
			      const BaseMatrix & ac)
  {
    SetMatrix (aa);
    SetPrecond (ac);
    SetPrecision (1e-10);
    SetMaxSteps (200);
  }

  void EigenSystem :: SetMatrix (const BaseMatrix & aa) 
  {
    a = &aa;
  }

  void EigenSystem :: SetPrecond (const BaseMatrix & ac)
  {
    c = &ac;
  }

  void EigenSystem :: SetMaxSteps (int amaxsteps)
  {
    maxsteps = amaxsteps;
  }

  void EigenSystem :: SetPrecision (double aprec)
  {
    prec = aprec;
  }






 
  int EigenSystem :: Calc ()
  {
    int retval = 0;

    auto v = a->CreateColVector();
    auto u = a->CreateColVector();
    auto w = a->CreateColVector();
    auto cv = a->CreateColVector();

    int it;
    double gamma, err;

    (*testout).precision(10);

    ai.SetSize(0);
    bi.SetSize(0);

    bool is_real = (dynamic_cast< S_BaseVector<double> *>(&*cv));// ||
		    //dynamic_cast< S_BaseVector<float> *>(&cv));

    /*    
    for(int i=0; i<v.Size(); i++)
      {
	if(is_real)
	  v.FVDouble()(i) = double (rand()) / RAND_MAX;
	else
	  v.FVComplex()(i) = double (rand()) / RAND_MAX;
      }
    */
    if (is_real)
      {
        for(int i = 0; i < v.FVDouble().Size(); i++)
          v.FVDouble()(i) = double (rand()) / RAND_MAX;
      }
    else
      for(int i = 0; i < v.FVComplex().Size(); i++)
        v.FVComplex()(i) = double (rand()) / RAND_MAX;

    v.Distribute();
    v.Cumulate();
      
    cv = 0;

    if (c)
      *cv = (*c) * *v;
    else
      *cv = *v;

    if(is_real)
      gamma = InnerProduct (*cv, *v);
    else
      gamma = (S_InnerProduct<Complex>(*cv,*v)).real();

    if (gamma < -1e-20)
      {
	cerr << "Eigensystem: Preconditioner negative" << endl;
	(*testout) << "Eigensystem: Preconditioner negative" << endl;
	(*testout) << "|v| = " << L2Norm(*v) << endl;
	(*testout) << "|cv| = " << L2Norm(*cv) << endl;
	(*testout) << "(v,cv) = " << gamma << endl;
	return 1;
      }

    gamma = sqrt (gamma);
    *v /= gamma;
    *cv /= gamma;
  
    *u = (*a) * *cv;
    it = 0;
    bi.Append (0);
    err = 1;

    while (err > prec && it < maxsteps)
      {
	//      cout << "." << flush;
	cout << IM(3) << "eigen-it " << it << "/" << maxsteps << endl;
	it++;
	double cvu;
	if(is_real)
	  cvu = InnerProduct (*cv, *u);
	else
	  cvu = (S_InnerProduct<Complex>(*cv,*u)).real();

	ai.Append (cvu);

	*w = *u - ai[it-1] * *v;

	if (c)
	  *cv = (*c) * *w;
	else
	  *cv = *w;

	if(is_real)
	  gamma = InnerProduct (*cv, *w);
	else
	  gamma = (S_InnerProduct<Complex>(*cv,*w)).real();

	if (gamma < -1e-20)
	  {
	    (*testout) << "Eigensystem: Preconditioner negative" << endl;
	    (*testout) << "|w| = " << L2Norm(*w) << endl;
	    (*testout) << "|cw| = " << L2Norm(*cv) << endl;
	    (*testout) << "(w,cw) = " << gamma << endl;
	    cerr << "Eigensystem: Preconditioner negative" << endl;
	    retval = 1;
	    break;
	  }
       
	if (gamma < 1e-40) break;
	gamma = sqrt (gamma);

	bi.Append (gamma);
	*cv /= gamma;
	*w /= gamma;
      
	*u = (*a) * *cv;
	*u -= gamma * *v;
	*v = *w;
	err *= fabs (gamma / ai[it-1]);

	//(*testout) << "|w| = " << L2Norm(w) << endl;
	//(*testout) << "|cw| = " << L2Norm(cv) << endl;
	//(*testout) << "(w,cw) = " << gamma << endl;
	//(*testout) << "err " << err << " prec " << prec << endl;
      }

    
      (*testout) << "it = " << it << " maxstep = " << maxsteps  
	   << " prec = " << prec << "  err = " << err  << endl; 
      /*
	//  << " gamma = " << gamma << endl;

      cout << endl;
      cout << "its = " << it << endl;
    */
    if (it >= maxsteps)
      {
	cout << IM(1) << "maxsteps " << maxsteps << " exceeded !!" << endl;
	retval = 2;
      }

    return retval;
  }


  int EigenSystem :: NumEigenValues () const
  {
    return ai.Size();
  }

  double EigenSystem :: MaxEigenValue () const
  {
    return EigenValue (NumEigenValues());
  }

  /*

  Eigenwertberechnung einer symmetrischen Tridiagonalmatrix:

  a1  b2
  b2  a2 b3
  b3 a3 ..
  .. .. bn
  bn an


  n ..... Matrixgroesse
  lami .. welcher Eigenwert (1 .. kleinster, n .. groesster
  alpha . Diagonale
  beta .. Nebendiagonale (beta(1) muss 0 sein !)


  Algorithmus nach:

  W. Barth, R.S. Martin, J.H. Wilkinson
  Numerische Mathematik 9, pp 386 - 393 (1967)

  */


  double EigenSystem::EigenValue (int nr) const
  {
    int n = ai.Size();
    int i, zbelow;
    double xl, xu, x, q;

    xu = 0;
    for (i = 0; i < n; i++)
      {
	q = fabs (ai[i]) + fabs (bi[i]) +
	  ((i < n-1) ? fabs (bi[i+1]) : 0);
	if (q > xu) xu = q;
      }

    xl = -xu;

    while (xu - xl > 1e-15 * fabs(xu) &&
	   xu - xl > 1e-100)
      {
	x = 0.5 * (xu + xl);

	zbelow = 0;
	q = 1;

	for (i = 0; i < n; i++)
	  {
	    if ( fabs (q) > 1e-100 )
	      q = ai[i] - x - bi[i] * bi[i] / q;
	    else
	      q = ai[i] - x - fabs(bi[i]) * 1e100;

	    if (q < 0) zbelow++;
	  }

	if (zbelow < nr)  
	  xl = x;
	else
	  xu = x;
      }

    return 0.5 * (xl + xu);
  }



  void EigenSystem :: PrintEigenValues (ostream & ost) const
  {
    for (int i = 1; i <= NumEigenValues(); i++)
      ost << "lam(" << i << ") = " << EigenValue(i) << endl;
  }

}
