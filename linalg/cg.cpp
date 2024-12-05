/**************************************************************************/
/* File:   cg.cpp                                                         */
/* Author: Joachim Schoeberl                                              */
/* Date:   5. Jul. 96                                                     */
/**************************************************************************/

/* 

  Conjugate Gradient Soler
  
*/ 

#include "cg.hpp"
#include "sparsematrix.hpp"

namespace ngla
{
  inline double Abs (const double & v)
  {
    return fabs (v);
  }

  inline double Abs (const Complex & v)
  {
    return std::abs (v);
  }


  KrylovSpaceSolver :: KrylovSpaceSolver ()
  {
    //      SetSymmetric();
    
    a = 0;  
    c = 0;
    SetPrecision (1e-10);
    SetMaxSteps (200); 
    SetInitialize (1);
    printrates = 0;
    useseed = false;
  }
  

  KrylovSpaceSolver :: KrylovSpaceSolver (shared_ptr<BaseMatrix> aa)
  {
    //  SetSymmetric();
    
    SetMatrix (aa);
    c = NULL;
    SetPrecision (1e-10);
    SetMaxSteps (200);
    SetInitialize (1);
    printrates = 0;
    useseed = false;
  }



  KrylovSpaceSolver :: KrylovSpaceSolver (shared_ptr<BaseMatrix> aa, shared_ptr<BaseMatrix> ac)
  {
    //  SetSymmetric();
    
    SetMatrix (aa);
    SetPrecond (ac);
    SetPrecision (1e-8);
    SetMaxSteps (200);
    SetInitialize (1);
    printrates = 0;
    useseed = false;
  }

    template <class SCAL>
  void BruteInnerProduct(const BaseVector & a, const BaseVector & b, Vector<SCAL> & result, const int start = 0)
  {
    const SCAL * pa;
    const SCAL * pb;
    int i;

    for(int i=start; i<result.Size(); i++)
      result[i] = 0;

    
    if(start == 0)
      for(i=0, pa = (SCAL*)(a.Memory()), pb = (SCAL*)(b.Memory()); i<a.Size()*result.Size(); i++,pa++,pb++)
	result[i%result.Size()] += (*pa)*(*pb);
    else
      {
	pa = (SCAL*)(a.Memory());
	pb = (SCAL*)(b.Memory());
	for(i=0; i<a.Size();i++)
	  {
	    pa += start;
	    pb += start;
	
	    for(int j=start; j<result.Size(); j++)
	      {
		result[j] += (*pa)*(*pb);
		pa++;
		pb++;
	      }
	  }
      }

  }


  template <class SCAL>
  void BruteInnerProduct2(const BaseVector & a, const BaseVector & b, Vector<SCAL> & result, const int start)
  {
    const SCAL * pa;
    const SCAL * pb;
    int i;

    for(int i=start; i<result.Size(); i++)
      result[i] = 0;

    pa = (SCAL*)(a.Memory());
    pb = (SCAL*)(b.Memory());
    for(i=0; i<a.Size();i++)
      {
	pb += start;

	for(int j=start; j<result.Size(); j++)
	  {
	    result[j] += (*pa)*(*pb);
	    pb++;
	  }
	pa++;
      }
      
  }

  template <class IPTYPE>
  void CGSolver<IPTYPE> :: MultiMult (const BaseVector & f, BaseVector & u, const int dim) const
  {
    try
      {
	// Solve A u = f
        BaseStatusHandler::SetThreadPercentage(0);

	auto d = f.CreateVector();
	auto w = f.CreateVector();
	auto s = f.CreateVector();

	int n = 0;
	Vector<SCAL> al(dim), be(dim), wd(dim), wdn(dim), kss(dim);
	double err;

	if (initialize)
	  {
	    u = 0.0;
	    d = f;
	  }
	else
	  {
	    d = f - (*a) * u;
	  }
	if (c)
	  w = (*c) * d;
	else
	  w = d;

	s = w;
	
	BruteInnerProduct(w,d,wdn);	 

	if (printrates) cout << IM(1) << "0 " << sqrt(L2Norm(wdn)) << endl;
	if (L2Norm(wdn) == 0.0) wdn = 1;	

	if(stop_absolute)
	  err = prec * prec;
	else
	  err = prec * prec * L2Norm (wdn);
	
	double lwstart = log(L2Norm(wdn));
	double lerr = log(err);
	

	while (n++ < maxsteps && L2Norm(wdn) > err && !(BaseStatusHandler::ShouldTerminate()))
	  {
	    w = (*a) * s;

	    wd = wdn;

	    BruteInnerProduct(s,w,kss);
	   
	    //(*testout) << "INNERPROD kss " <<kss << endl;
	    if (L2Norm(kss) == 0.0) break;
	    
	    for(int i = 0; i<dim; i++)
	      al[i] = wd[i] / kss[i];
	    
	    SCAL * pl;
	    const SCAL * pr;

	    int i;

	    for(pl = (SCAL*)(u.Memory()), pr = (SCAL*)(s.Memory()), i=0; i<dim*u.Size(); i++,pl++,pr++)
	      *pl += al[i%dim]*(*pr);
	      
	    for(pl = (SCAL*)(d.Memory()), pr = (SCAL*)(w.Memory()), i=0; i<dim*u.Size(); i++,pl++,pr++)
	      *pl -= al[i%dim]*(*pr);
	      

	    //u += al * s;
	    //d -= al * w;

	    if (c)
	      w = (*c) * d;
	    else
	      w = d;

	    BruteInnerProduct(w,d,wdn);

	    //(*testout) << "wdn " << wdn << endl;
	    
	    for(int i = 0; i<dim; i++)
	      be[i] = wdn[i] / wd[i];
	    
	    for(pl = (SCAL*)(s.Memory()), pr = (SCAL*)(w.Memory()), i=0; i<dim*s.Size(); i++,pl++,pr++)
	      *pl = (*pl)*be[i%dim] + *pr;

	    //s *= be;
	    //s += w;

	    if (printrates ) cout << IM(1) << n << " " << sqrt(L2Norm (wdn)) << endl;
            BaseStatusHandler::SetThreadPercentage(100.*max2(double(n)/double(maxsteps),
						(lwstart-log(L2Norm(wdn)))/(lwstart-lerr)));
	  } 
	
	const_cast<int&> (steps) = n;
	
        /*
	delete &d;
	delete &w;
	delete &s;
        */
      }

    catch (Exception & e)
      {
	e.Append ("in caught in CGSolver::Mult\n");
	throw;
      }
    catch (exception & e)
      {
	throw Exception(e.what() +
			string ("\ncaught in CGSolver::Mult\n"));
      }
  }


  template <class IPTYPE>
  void CGSolver<IPTYPE> :: MultiMultSeed (const BaseVector & f, BaseVector & u, const int dim) const
  {
    try
      {
	// Solve A u = f
        BaseStatusHandler::SetThreadPercentage(0);
 
	SCAL * pl;
	const SCAL * pr;
	int i;

	auto d = f.CreateVector();

	BaseMatrix * smalla;
        /*
	if(dynamic_cast< const SparseMatrixSymmetricTM<SCAL> *>(a))
	  smalla = new SparseMatrixSymmetric<SCAL,SCAL>(*dynamic_cast< const SparseMatrixSymmetricTM<SCAL> *>(a));
	else
        */
        if (dynamic_cast< const SparseMatrixTM<SCAL> *>(a.get()))
	  smalla = new SparseMatrix<SCAL,SCAL>(*dynamic_cast< const SparseMatrixTM<SCAL> *>(a.get()));
	else
	  throw Exception("Assumption about bilinearform wrong.");


	//BaseVector & aux1 = (smalla) ? d : *f.CreateVector();
	//BaseVector & aux2 = (smalla) ? d : *f.CreateVector();
	

	VVector<SCAL> w(f.Size());
	VVector<SCAL> d_reduced(f.Size());
	VVector<SCAL> s(f.Size());

	int n = 0;

	SCAL be,wd,wdn,kss;
	Vector<SCAL> al(dim);
	Array<double> err(dim);

	if (initialize)
	  {
	    u = 0.0;
	    d = f;
	  }
	else
	  {
	    d = f - (*a) * u;
	  }

		
	double lwstart;
	double lerr;
	


	for(int seed = dim-1; seed >= 0; seed--)
	  {
	    
	    pr = (SCAL*)(d.Memory());
	    pr += seed;

	    for(i=0, pl = (SCAL*)(d_reduced.Memory()); i<d.Size(); i++, pl++)
	      {
		(*pl) = (*pr);
		pr += dim;
	      }
	    
	    
	   
	    if (c)
	      w = (*c) * d_reduced;
	    else
	      w = d_reduced;

	    if(stop_absolute)
	      err[seed] = prec * prec;
	    else
	      err[seed] = prec * prec * Abs (S_InnerProduct<SCAL>(w,d_reduced));
	  }


	for(int seed = 0; seed < dim; seed++)
	  {
	    (*testout) << "seed " << seed << endl;

	    if(seed > 0)
	      {
		pr = (SCAL*)(d.Memory());
		pr += seed;

		for(i=0, pl = (SCAL*)(d_reduced.Memory()); i<d.Size(); i++, pl++)
		  {
		    (*pl) = (*pr);
		    pr += dim;
		  }
		
		
		
		if (c)
		  w = (*c) * d_reduced;
		else
		  w = d_reduced;
	      }
	    
	    s = w;	    
	    
	    wdn = S_InnerProduct<SCAL>(w,d_reduced);
	    
	    
	    if (printrates ) cout << IM(1) << n << " (block " << seed+1 << ") " << sqrt (Abs (wdn)) << endl;
	    if(Abs(wdn) == 0.0) wdn = 1;

	    lwstart = log(Abs(wdn));
	    lerr = log(err[seed]);
	    


	    while (n++ < maxsteps && Abs(wdn) > err[seed] && !(BaseStatusHandler::ShouldTerminate()))
	      {
		//if(smalla)
		w = (*smalla)  * s;
		/*
		else
		  {
		    pl = (SCAL*)(aux1.Memory());
		    pr = (SCAL*)(s.Memory());
		    for(i=0; i<s.Size(); i++)
		      {
			for(int j=0; j<dim; j++)
			  {
			    *pl = *pr;
			    pl++;
			  }
			pr++;
		      }
		    aux2 = (*a) * aux1;
		    pl = (SCAL*)(w.Memory());
		    pr = (SCAL*)(aux2.Memory());
		    for(i=0; i<s.Size(); i++)
		      {
			*pl = *pr;
			pl++;
			pr += dim;
		      }
		  }
		*/

		//w = (*a) * s;
		
		wd = wdn;
		
		kss = S_InnerProduct<IPTYPE> (s, w);
		if (kss == 0.0) break;
		

		BruteInnerProduct2(s,d,al,seed+1);
		al[seed] = wd;
		
		for(i=seed; i<dim; i++)
		  al[i] /= kss;

		
		
		//(*testout) << "al " << al << endl;
		
		pl = (SCAL*)(u.Memory());
		pr = (SCAL*)(s.Memory());
		for(i=0; i<u.Size(); i++)
		  {
		    pl += seed;

		    for(int j=seed; j<dim; j++)
		      {
			*pl += al[j]*(*pr);
			pl++;
		      }
		    pr++;
		  }
		
		pl = (SCAL*)(d.Memory());
		pr = (SCAL*)(w.Memory());
		for(i=0; i<d.Size(); i++)
		  {
		    pl += seed;

		    for(int j=seed; j<dim; j++)
		      {
			*pl -= al[j]*(*pr);
			pl++;
		      }
		    pr++;
		  }
				
		//u += al * s;
		//d -= al * w;


		
		pr = (SCAL*)(d.Memory());
		pr += seed;

		for(i=0, pl = (SCAL*)(d_reduced.Memory()); i<d.Size(); i++, pl++)
		  {
		    *pl = *pr;
		    pr += dim;
		  }

		
		if (c)
		  w = (*c) * d_reduced;
		else
		  w = d_reduced;

		wdn = S_InnerProduct<IPTYPE> (d_reduced, w);

		be = wdn/wd;
		
		s *= be;
		s += w;

		if (printrates ) cout << IM(1) << n << " (block " << seed+1 << ") " << sqrt (Abs (wdn)) << endl;
                BaseStatusHandler::SetThreadPercentage(100.*max2(double(n)/double(maxsteps),
						    (lwstart-log(Abs(wdn)))/(lwstart-lerr)));
	      } 
	  }
	const_cast<int&> (steps) = n;
	
	/*
	if(!smalla)
	  {
	    delete &aux1;
	    delete &aux2;
	  }
	*/
	delete smalla;
      }

    catch (Exception & e)
      {
	e.Append ("in caught in CGSolver::Mult\n");
	throw;
      }
    catch (exception & e)
      {
	throw Exception(e.what() +
			string ("\ncaught in CGSolver::Mult\n"));
      }
  }


  template <class IPTYPE>
  void CGSolver<IPTYPE> :: Mult (const BaseVector & f, BaseVector & u) const
  {
    static Timer timer ("CG solver");
    RegionTimer reg (timer);

    int dim = 1;

    if(dynamic_cast<VVector< Vec<2, SCAL> >* >(&u))
      dim = 2;
    else if(dynamic_cast<VVector< Vec<3, SCAL> >* >(&u))
      dim = 3;
    else if(dynamic_cast<VVector< Vec<4, SCAL> >* >(&u))
      dim = 4;
    else if(dynamic_cast<VVector< Vec<5, SCAL> >* >(&u))
      dim = 5;
    else if(dynamic_cast<VVector< Vec<6, SCAL> >* >(&u))
      dim = 6;
    else if(dynamic_cast<VVector< Vec<7, SCAL> >* >(&u))
      dim = 7;
    else if(dynamic_cast<VVector< Vec<8, SCAL> >* >(&u))
      dim = 8;
    /*
    else if(dynamic_cast<VVector< Vec<9, SCAL> >* >(&u))
      dim = 9;
    else if(dynamic_cast<VVector< Vec<10, SCAL> >* >(&u))
      dim = 10;
    else if(dynamic_cast<VVector< Vec<11, SCAL> >* >(&u))
      dim = 11;
    else if(dynamic_cast<VVector< Vec<12, SCAL> >* >(&u))
      dim = 12;
    else if(dynamic_cast<VVector< Vec<13, SCAL> >* >(&u))
      dim = 13;
    else if(dynamic_cast<VVector< Vec<14, SCAL> >* >(&u))
      dim = 14;
    else if(dynamic_cast<VVector< Vec<15, SCAL> >* >(&u))
      dim = 15;
    */
    //cout << "useseed: " << useseed << " dim: " << dim << endl;

    if(useseed && dim != 1)
      {
	MultiMultSeed(f,u,dim);
	//MultiMult(f,u,dim);
	return;
      }
 
    
    try
      {
	// Solve A u = f
        BaseStatusHandler::SetThreadPercentage(0);
 
        auto w = u.CreateVector();
        auto s = u.CreateVector();
        auto d = f.CreateVector();
        auto as = f.CreateVector();
        
	int n = 0;
	SCAL al, be, wd, wdn, kss;
	double err;
	if (initialize)
	  {
	    u = 0.0;
	    d = f;
	  }
	else
	  {
	    d = f - (*a) * u;
	  }

	if (c)
	  w = (*c) * d;
	else
	  w = d;

	s = w;
	wdn = S_InnerProduct<IPTYPE> (w,d);

	if (printrates) cout << IM(1) << "0 " << sqrt(Abs(wdn)) << endl;
	if (wdn == 0.0) wdn = 1;	

	if(stop_absolute)
	  err = prec * prec;
	else
	  err = prec * prec * Abs (wdn);
	
	double lwstart = log(Abs(wdn));
	double lerr = log(err);
	
	while (n++ < maxsteps && Abs(wdn) > err && !(BaseStatusHandler::ShouldTerminate()))
	  {
	    as = (*a) * s;
	    wd = wdn;
	    kss = S_InnerProduct<IPTYPE> (s, as);
	    if (kss == 0.0) break;
	    
	    al = wd / kss;
	    u += al * s;
	    d -= al * as;
            
	    if (c)
	      w = (*c) * d;
	    else
	      w = d;
	    wdn = S_InnerProduct<IPTYPE> (d, w);

	    be = wdn / wd;
	    
	    s *= be;
	    s += w;

	    if (printrates ) cout << IM(1) << n << " " << sqrt (Abs (wdn)) << endl;
            BaseStatusHandler::SetThreadPercentage(100.*max2(double(n)/double(maxsteps),
						(lwstart-log(Abs(wdn)))/(lwstart-lerr)));
	  } 
	
	const_cast<int&> (steps) = n;
      }

    catch (Exception & e)
      {
	e.Append ("in caught in CGSolver::Mult\n");
	throw;
      }
    catch (exception & e)
      {
	throw Exception(e.what() +
			string ("\ncaught in CGSolver::Mult\n"));
      }
  }





  template <class IPTYPE>
  void BiCGStabSolver<IPTYPE> :: Mult (const BaseVector & f, BaseVector & u) const
  {
    
    try
      {
	// Solve A u = f
        BaseStatusHandler::SetThreadPercentage(0);
 
	auto r = f.CreateVector();
	auto r_tilde = f.CreateVector();
	auto p = f.CreateVector();
	auto p_tilde = f.CreateVector();
	auto s = f.CreateVector();
	auto s_tilde = f.CreateVector();
	auto t = f.CreateVector();
	auto v = f.CreateVector();

	int n = 0;
	SCAL rho_old, rho_new, beta, alpha, omega;
	double err, err_i;

	if (initialize)
	  {
	    u = 0.0;
	    r = f;
	  }
	else
	  {
	    r = f - (*a) * u;
	  }
	r_tilde = r;

	rho_new = S_InnerProduct<IPTYPE>(r_tilde, r);
	p = r;
	if (c)
	  p_tilde = (*c) * p;
	else
	  p_tilde = p;

	v = (*a) * p_tilde;
	alpha = rho_new / S_InnerProduct<IPTYPE> (r_tilde, v);
	s = r;
	s -= alpha * v;

	err_i = L2Norm(s);
	if (c)
	  s_tilde = (*c) * s;
	else
	  s_tilde = s;

	t = (*a) * s_tilde;

	omega = S_InnerProduct<IPTYPE> (t, s) / S_InnerProduct<IPTYPE> (t, t);
	u += alpha * p_tilde + omega * s_tilde;
	r = s;
	r -= omega * t;

	err_i = L2Norm(r);
	if (printrates) cout << IM(1) << "0 " << err_i << endl;


	if(stop_absolute)
	  err = prec * prec;
	else
	  err = prec * prec * err_i;
	
	double lwstart = log(err_i);
	double lerr = log(err);
	

	while (n++ < maxsteps && err_i > err && !(BaseStatusHandler::ShouldTerminate()))
	  {
	    rho_old = rho_new;
	    rho_new = S_InnerProduct<IPTYPE>(r_tilde, r);
	    beta = (rho_new / rho_old ) * ( alpha / omega );
	    p = r;
	    p += beta * p;
	    p -= beta*omega * v;

	    if (c)
	      p_tilde = (*c) * p;
	    else
	      p_tilde = p;
	    
	    v = (*a) * p_tilde;
	    alpha = rho_new / S_InnerProduct<IPTYPE> (r_tilde, v);
	    s = r;
	    s -= alpha * v;

	    err_i = L2Norm(s);
	    u += alpha * p_tilde;
	    
	    if ( err_i < err )
	      {
		break;
	      }

	    if (c)
	      s_tilde = (*c) * s;
	    else
	      s_tilde = s;

	    t = (*a) * s_tilde;
	    
	    omega = S_InnerProduct<IPTYPE> (t, s) / S_InnerProduct<IPTYPE> (t, t);
	    u +=  omega * s_tilde;
	    r = s;
	    r -= omega * t;

	    err_i = L2Norm(r);

	    if (printrates ) cout << IM(1) << n << " " << err_i << endl;
            BaseStatusHandler::SetThreadPercentage(100.*max2(double(n)/double(maxsteps),
						(lwstart-log(err_i))/(lwstart-lerr)));
	  } 
	
	const_cast<int&> (steps) = n;
      }

    catch (Exception & e)
      {
	e.Append ("in caught in BiCGStabSolver::Mult\n");
	throw;
      }
    catch (exception & e)
      {
	throw Exception(e.what() +
			string ("\ncaught in BiCGStabSolver::Mult\n"));
      }
  }




  template <class IPTYPE>
  void SimpleIterationSolver<IPTYPE> :: Mult (const BaseVector & f, BaseVector & u) const
  {

  try
      {
	// Solve A u = f
        BaseStatusHandler::SetThreadPercentage(0);
 
	auto d = f.CreateVector();
	auto w = f.CreateVector();

	int n = 0;
	double err, err0;

	if (initialize)
	  {
	    u = 0.0;
	    d = f;
	  }
	else
	  {
	    d = f - (*a) * u;
	  }


        err = err0 = 1;

	while (n++ < maxsteps && err > prec * err0)
          {
            d = f - (*a) * u;

            if (c)
              w = (*c) * d;
            else
              w = d;

            u += tau * w;

            err = Abs (S_InnerProduct<IPTYPE> (w, d));
            if (n == 1) err0 = err;

	    if (printrates ) cout << IM(1) << n << " " << sqrt (err) << endl;
          }

	const_cast<int&> (steps) = n;
      }

    catch (Exception & e)
      {
	e.Append ("in caught in SimpleIterationSolver::Mult\n");
	throw;
      }
    catch (exception & e)
      {
	throw Exception(e.what() +
			string ("\ncaught in SimpleIterationSolver::Mult\n"));
      }
  }





















  template <class IPTYPE>
  void GMRESSolver<IPTYPE> :: Mult (const BaseVector & f, BaseVector & x) const
  {
    // from Wikipedia

    try
      {
	// Solve A u = f

	auto v = f.CreateVector();
	auto av = f.CreateVector();
	auto r = f.CreateVector();
	auto w = f.CreateVector();
	auto hv = f.CreateVector();

        Array<AutoVector> vi(maxsteps);
        Matrix<SCAL> h(maxsteps+1, maxsteps);
        Matrix<SCAL> h2(maxsteps+1, maxsteps);
        Vector<SCAL> gammai(maxsteps), ci(maxsteps), si(maxsteps);


        h = SCAL(0.0);
        h2 = SCAL(0.0);

	if (initialize)
	  {
	    x = 0.0;
	    r = f;
	  }
	else
	  {
	    r = f - (*a) * x;
	  }

	if (c)
          {
            hv = (*c) * r;
            r = hv;
          }


        double norm = r.L2Norm();
        v = (1.0/sqrt(S_InnerProduct<IPTYPE>(r,r))) * r;

        gammai(0) = norm;

	if (printrates) cout << IM(1) << "0 " << norm << endl;
	
	double err;
	if(stop_absolute)
	  err = prec;
	else
	  err = prec * Abs (norm);
	
	int j = -1;
	while (j++ < maxsteps-2 && norm > err)
	  {
            vi[j].AssignPointer (f.CreateVector());
            vi[j] = v;

            av = (*a) * v;
            if (c)
              {
                hv = (*c) * av;
                av = hv;
              }

            for (int i = 0; i <= j; i++)
              h2(i,j) = h(i,j) = S_InnerProduct<IPTYPE> (*vi[i], av);

            w = av;
            for (int i = 0; i <= j; i++)
              w -= h(i,j) * (*vi[i]);

            v = (1.0 / sqrt (S_InnerProduct<IPTYPE> (w, w))) * w;
            h2(j+1,j) = h(j+1,j) = S_InnerProduct<IPTYPE> (v, av);

            for (int i = 0; i < j; i++)
              {
                SCAL hi = h(i,j), hip = h(i+1, j);
                h(i,j)   = ci(i+1) * hi + si(i+1) * hip;
                h(i+1,j) = si(i+1) * hi - ci(i+1) * hip;
              }
            SCAL beta = sqrt ( sqr(h(j,j)) + sqr(h(j+1,j)));
            si(j+1) = h(j+1,j) / beta;
            ci(j+1) = h(j,j) / beta;
            h(j,j) = beta;
            gammai(j+1) = si(j+1) * gammai(j);
            gammai(j) = ci(j+1) * gammai(j);
            
	    if (printrates ) cout << IM(1) << j 
                                  << " ci = " << ci(j+1) 
                                  << " si = " << si(j+1) 
                                  << " gammi = " << gammai(j) << endl;


            norm = fabs (gammai(j));
          }
        
        j--;
        cout << IM(5) << "gmres - Triangular matrix" << endl << h.Rows(0,j+2).Cols(0,j+2) << endl;
        Vector<SCAL> y(maxsteps);
        for (int i = j; i >= 0; i--)
          {
            SCAL sum = gammai(i);
            for (int k = i+1; k <= j; k++)
              sum -= h(i,k) * y(k);
            y(i) = sum / h(i,i);
          }

        for (int i = 0; i <= j; i++)
          x += y(i) * *vi[i];

	const_cast<int&> (steps) = j;
	
        /*
        *testout << "h2 = " << endl << h2 << endl;

        for (int k = 0; k < 10; k++)
          for (int l = 0; l < 10; l++)
            *testout << "< v(" << k << ") , v(" << l << ") > = " 
                     << S_InnerProduct<IPTYPE> (*vi[k], *vi[l]) << endl;
        
        for (int k = 0; k < 10; k++)
          {
            hv = (*a) * (*vi[k]);
            av = (*c) * hv;
            for (int l = 0; l < 10; l++)
              *testout << "< Av(" << k << ") , v(" << l << ") > = " 
                       << S_InnerProduct<IPTYPE> (av, *vi[l]) << endl;
          }


        Matrix<SCAL> hs(j+1,j+1), hsinv(j+1,j+1);
        Vector<SCAL> rs(j+1), us(j+1);
        for (int i = 0; i <= j; i++)
          for (int k = 0; k <= j; k++)
            hs(i,k) = h2(i,k);

        CalcInverse (hs, hsinv);
        rs = SCAL(0.0);
        rs(0) = 1.0;
        us = hsinv * rs;
        
        x = 0.0;
        for (int i = 0; i <= j; i++)
          x += us(i) * *vi[i];
        */
      }

    catch (Exception & e)
      {
	e.Append ("in caught in GMRESSolver::Mult\n");
	throw;
      }
    catch (exception & e)
      {
	throw Exception(e.what() +
			string ("\ncaught in GMRESSolver::Mult\n"));
      }
  }









//*****************************************************************
// Iterative template routine -- QMR
//
// QMR.h solves the unsymmetric linear system Ax = b using the
// Quasi-Minimal Residual method following the algorithm as described
// on p. 24 in the SIAM Templates book.
//
//   -------------------------------------------------------------
//   return value     indicates
//   ------------     ---------------------
//        0           convergence within max_iter iterations
//        1           no convergence after max_iter iterations
//                    breakdown in:
//        2             rho
//        3             beta
//        4             gamma
//        5             delta
//        6             ep
//        7             xi
//   -------------------------------------------------------------
//   
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax=b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************



template <class SCAL>
void QMRSolver<SCAL> :: Mult (const BaseVector & b, BaseVector & x) const
{
  try
    {
      cout << IM(1) << "QMR called" << endl;
      double resid;
      SCAL rho, rho_1, xi, gamma, gamma_1, theta, theta_1, eta, delta, ep=1.0, beta;
      

      auto r = b.CreateVector();
      auto v_tld = b.CreateVector();
      auto y = b.CreateVector();
      auto w_tld = b.CreateVector();
      auto z = b.CreateVector();
      auto v = b.CreateVector();
      auto w = b.CreateVector();
      auto y_tld = b.CreateVector();
      auto z_tld = b.CreateVector();
      auto p = b.CreateVector();
      auto q = b.CreateVector();
      auto p_tld = b.CreateVector();
      auto d = b.CreateVector();
      auto s = b.CreateVector();

      double normb = b.L2Norm();


      if (initialize)
	x = 0;


      r = b - (*a) * x;

      if (normb == 0.0)
	normb = 1;
      
      cout.precision(12);
      
      // 
      double tol = prec;
      int max_iter = maxsteps;
      
      if ((resid = r.L2Norm() / normb) <= tol) {
	tol = resid;
	max_iter = 0;
	((int&)status) = 0;
	return;
      }
  
      v_tld = r;

      // use preconditioner c1
      if (c)
	y = (*c) * v_tld;
      else
	y = v_tld;

      rho = y.L2Norm();
      
      w_tld = r;

      if (c2) 
	z = Transpose (*c2) * w_tld; 
      // z = (*c2) * w_tld; 
      else
	z = w_tld;
      
      xi = z.L2Norm();

      gamma = 1.0;
      eta = -1.0;
      theta = 0.0;
      ((int&)steps) = 0;


      for (int i = 1; i <= max_iter; i++) 
	{

	  ((int&)steps) = i;  
	  
	  if (rho == 0.0)
	    {
	      (*testout) << "QMR: breakdown in rho" << endl;
	      ((int&)status) = 2;
	      return;                        // return on breakdown
	    }
	  
	  if (xi == 0.0)
	    {
	      (*testout) << "QMR: breakdown in xi" << endl;
	      ((int&)status) = 7;
	      return;                        // return on breakdown
	    }

	  v = (1.0/rho) * v_tld;
	  y /= rho;

	  w = (1.0/xi) * w_tld;
	  z /= xi;


	  delta = S_InnerProduct<SCAL> (z, y);
	  if (delta == 0.0)
	    {
	      (*testout) << "QMR: breakdown in delta" << endl;
	      ((int&)status) = 5;
	      return;                        // return on breakdown
	    }

	  
	  if (c2) 
	    y_tld = (*c2) * y;
	  else
	    y_tld = y;

	  
	  if (c)
	    z_tld = Transpose (*c) * z;
	  // z_tld = (*c) * z;
	  else
	    z_tld = z;

	  if (i > 1) 
	    {
	      //  p = y_tld - (xi(0) * delta(0) / ep(0)) * p;
	      //  q = z_tld - (rho(0) * delta(0) / ep(0)) * q;
	      p *= (-xi * delta / ep);
	      p += y_tld;
	      q *= (-rho * delta / ep);
	      q += z_tld;
	    } 
	  else 
	    {
	      p = y_tld;
	      q = z_tld;
	    }
	  
	  p_tld = (*a) * p;
	  ep = S_InnerProduct<SCAL> (q, p_tld);

	  if (ep == 0.0)
	    {
	      (*testout) << "QMR: breakdown in ep" << endl;
	      ((int&)status) = 6;
	      return;                        // return on breakdown
	    }

	  beta = ep / delta;
	  if (beta == 0.0)
	    {
	      (*testout) << "QMR: breakdown in beta" << endl;
	      ((int&)status) = 3;
	      return;                        // return on breakdown
	    }

	  v_tld = p_tld;
	  v_tld -= beta * v;

	  if (c)
	    y = (*c) * v_tld;
	  else
	    y = v_tld;


	  rho_1 = rho;
	  rho = y.L2Norm();

	  w_tld = Transpose(*a) * q;
	  w_tld -= beta * w;
	  
	  if (c2) 
	    z = Transpose (*c2) * w_tld;
	  // z = (*c2) * w_tld;
	  else
	    z = w_tld;
	  
	  xi = z.L2Norm();
	  
	  gamma_1 = gamma;
	  theta_1 = theta;
	  
	  theta = rho / (gamma_1 * Abs(beta));    // abs (beta) ???
	  gamma = 1.0 / sqrt(1.0 + theta * theta);
	  
	  if (gamma == 0.0)
	    {
	      (*testout) << "QMR: breakdown in gamma" << endl;
	      ((int&)status) = 4;
	      return;                        // return on breakdown
	    }
	  
	  eta = -eta * rho_1 * gamma * gamma / 
	    (beta * gamma_1 * gamma_1);

	  if (i > 1) 
	    {
	      // d = eta(0) * p + (theta_1(0) * theta_1(0) * gamma(0) * gamma(0)) * d;
	      // s = eta(0) * p_tld + (theta_1(0) * theta_1(0) * gamma(0) * gamma(0)) * s;
	      d *= (theta_1 * theta_1 * gamma * gamma);
	      d += eta * p;
	      s *= (theta_1 * theta_1 * gamma * gamma);
	      s += eta * p_tld;
	    } 
	  else 
	    {
	      d = eta * p;
	      s = eta * p_tld;
	    }
	  
	  x += d;
	  r -= s;

	  if ( printrates ) cout << IM(1) << i << " " << r.L2Norm() << endl;
	  
	  if ((resid = r.L2Norm() / normb) <= tol) {
	    tol = resid;
	    max_iter = i;
	    ((int&)status) = 0;
	    return;
	  }
	}
      
      /*
      (*testout) << "no convergence" << endl;

      (*testout) << "res = " << endl << r << endl;
      (*testout) << "x = " << endl << x << endl;
      (*testout) << "b = " << endl << b << endl;
      */
      tol = resid;
      ((int&)status) = 1;
      return;                            // no convergence
    }

  

  catch (Exception & e)
    {
      e.Append ("in caught in QMRSolver::Mult\n"); 
      throw;
    }
  catch (exception & e)
    {
      throw Exception(e.what() +
		      string ("\ncaught in QMRSolver::Mult\n"));
    }
}
  
 
  
  template class CGSolver<double>;
  template class CGSolver<Complex>;
  template class CGSolver<ComplexConjugate>;
  template class CGSolver<ComplexConjugate2>;
  template class BiCGStabSolver<double>;
  template class BiCGStabSolver<Complex>;
  template class BiCGStabSolver<ComplexConjugate>;
  template class BiCGStabSolver<ComplexConjugate2>;
  template class SimpleIterationSolver<double>;
  template class SimpleIterationSolver<Complex>;
  template class SimpleIterationSolver<ComplexConjugate>;
  template class SimpleIterationSolver<ComplexConjugate2>;
  template class QMRSolver<double>;
  template class QMRSolver<Complex>;
  template class QMRSolver<ComplexConjugate>;
  template class QMRSolver<ComplexConjugate2>;
  template class GMRESSolver<double>;
  template class GMRESSolver<Complex>;
  template class GMRESSolver<ComplexConjugate>;
  template class GMRESSolver<ComplexConjugate2>;


}
