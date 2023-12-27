/*********************************************************************/
/* File:   recurisve_pol.cpp                                         */
/* Author: Almedin Becirovic, JS                                     */
/* Date:   14. Apr. 2003                                             */
/*********************************************************************/



// #include <fem.hpp>
#include <recursive_pol.hpp>

namespace ngfem
{

  Array< Vec<2> > LegendrePolynomial :: coefs;
  
  void LegendrePolynomial :: Calc (int n)
  {
    static mutex calclegendre_mutex;
    if (coefs.Size() > n) return;

    {
      lock_guard<mutex> guard(calclegendre_mutex);
      if (coefs.Size() <= n)
        {
          coefs.SetSize (n+1);
          
          coefs[0][0] = 1;
          coefs[1][1] = 1;
          for (int i = 1; i <= n; i++)
            {
              coefs[i][0] = CalcA(i);  // (2.0*i-1)/i;
              coefs[i][1] = CalcC(i);  // (1.0-i)/i;
            }
        }
    }
  }


  Array< Vec<2> > IntLegNoBubble :: coefs;

  void IntLegNoBubble :: Calc (int n)
  {
    static mutex calcintlegnobub_mutex;
    if (coefs.Size() > n) return;

    {
      lock_guard<mutex> guard(calcintlegnobub_mutex);
      if (coefs.Size() <= n)
        {
          coefs.SetSize (n+1);
          coefs[0][0] = coefs[0][1] = 1e10;
          // coefs[0][0] = -0.5;
          // coefs[1][1] = -0.5;
          for (int i = 1; i <= n; i++)
            {
              coefs[i][0] = CalcA(i);
              coefs[i][1] = CalcC(i);
            }
        }
    }
  }






  // int JacobiPolynomialAlpha :: maxnp;
  // int JacobiPolynomialAlpha :: maxalpha;
  // Array< Vec<4> > JacobiPolynomialAlpha :: coefs;
  Vec<4> JacobiPolynomialAlpha :: coefs[maxnp*maxalpha];


  void JacobiPolynomialAlpha :: Calc (int n, int alpha)
  {

    /*
    if (coefs.Size() < (n+1)*(alpha+1))
      {
        coefs.SetSize ((n+1)*(alpha+1));

        for (int a = 0; a <= alpha; a++)
          {
            for (int i = 1; i <= n; i++)
              {
                coefs[a*(n+1)+i][0] = CalcA (i, a, 0);
                coefs[a*(n+1)+i][1] = CalcB (i, a, 0);
                coefs[a*(n+1)+i][2] = CalcC (i, a, 0);
              }
            // ALWAYS_INLINE S P1(S x) const { return 0.5 * (2*(al+1)+(al+be+2)*(x-1)); }
            double al = a, be = 0;
            coefs[a*(n+1)+1][0] = 0.5 * (al+be+2);
            coefs[a*(n+1)+1][1] = 0.5 * (2*(al+1)-(al+be+2));
            coefs[a*(n+1)+1][2] = 0.0;
          }

        maxnp = n+1;
      }
    */


    // compile-time fixed alpha and maxnp
    for (int a = 0; a < alpha; a++)
      {
        for (int i = 1; i < maxnp; i++)
          {
            coefs[a*maxnp+i][0] = CalcA (i, a, 0);
            coefs[a*maxnp+i][1] = CalcB (i, a, 0);
            coefs[a*maxnp+i][2] = CalcC (i, a, 0);
          }
        // ALWAYS_INLINE S P1(S x) const { return 0.5 * (2*(al+1)+(al+be+2)*(x-1)); }
        double al = a, be = 0;
        coefs[a*maxnp+1][0] = 0.5 * (al+be+2);
        coefs[a*maxnp+1][1] = 0.5 * (2*(al+1)-(al+be+2));
        coefs[a*maxnp+1][2] = 0.0;
      }

  }



  int IntegratedJacobiPolynomialAlpha :: maxn;
  // int IntegratedJacobiPolynomialAlpha :: maxalpha;
  Array< Vec<4> > IntegratedJacobiPolynomialAlpha :: coefs;


  void IntegratedJacobiPolynomialAlpha :: Calc (int n, int alpha)
  {
    if (coefs.Size() < (n+1)*(alpha+1))
      {
        coefs.SetSize ((n+1)*(alpha+1));

        for (int a = 0; a <= alpha; a++)
          {
            for (int i = 1; i <= n; i++)
              {
                coefs[a*(n+1)+i][0] = CalcA (i, a, 0);
                coefs[a*(n+1)+i][1] = CalcB (i, a, 0);
                coefs[a*(n+1)+i][2] = CalcC (i, a, 0);
              }
            // ALWAYS_INLINE S P1(S x) const { return 0.5 * (2*(al+1)+(al+be+2)*(x-1)); }
            double alpha = a;
            coefs[a*(n+1)+1][0] = 0.25 * (alpha+2);
            coefs[a*(n+1)+1][1] = 0.25 * (alpha-2);
            coefs[a*(n+1)+1][2] = 0.0;
          }

        maxn = n;
      }
  }




  class InitRecPol
  {
  public:
    InitRecPol()
    {
      LegendrePolynomial::Calc(1000);
      IntLegNoBubble::Calc(1000);
      JacobiPolynomialAlpha::Calc(100,100); 
      IntegratedJacobiPolynomialAlpha::Calc(100,100); 
    }
  };

  InitRecPol init;













#ifdef OLD
  
  template <int DIM>
  double VertexExtensionOptimal<DIM> :: coefs[SIZE][SIZE];

  template <int DIM>
  bool VertexExtensionOptimal<DIM> :: initialized = 0;


  template <int DIM>
  VertexExtensionOptimal<DIM> :: VertexExtensionOptimal ()
  {
    if (!initialized)
      {
	// ofstream tout ("test1.out");
	int i, j, k;
	Matrix<double> a(SIZE, SIZE);
	Matrix<double> b(2, SIZE);

	try
	  {
	    const IntegrationRule & ir =
	      GetIntegrationRules().SelectIntegrationRule (ET_SEGM, 2*SIZE+2);

	    Vector<double> diff_leg(SIZE); 

	    a = 0;
	    b = 0;

	    for (i = 0; i < ir.GetNIP(); i++)
	      {
		double x = ir[i](0);
		DerivedLegendrePolynomial (SIZE-2, 2*x-1, diff_leg);
		for (j = SIZE-1; j >= 1; j--)
		  diff_leg(j) = 2 * diff_leg(j-1);
		diff_leg(0) = 0;

		double fac = ir[i].Weight();
		for (j = 0; j < DIM-1; j++)
		  fac *= (1-x);
	    
		for (j = 0; j < SIZE; j++)
		  for (k = 0; k < SIZE; k++)
		    a(j,k) += fac * diff_leg(j) * diff_leg(k);
	      }

	    double fac = 1;
	    for (j = 0; j < SIZE; j++)
	      {
		b(0,j) = fac;
		b(1,j) = 1;
		fac *= -1;
	      }


	    // compute coefficients
	    for (i = 1; i < SIZE; i++)
	      {
		Matrix<double> mat(i+3, i+3);
		Vector<double> rhs(i+3), sol(i+3);
		mat = 0;
		for (j = 0; j <= i; j++)
		  for (k = 0; k <= i; k++)
		    mat(j,k) = a(j,k);

		for (j = 0; j <= i; j++)
		  {
		    mat(i+1, j) = mat(j,i+1) = b(0,j);
		    mat(i+2, j) = mat(j,i+2) = b(1,j);
		  }
	    
		rhs = 0;
		rhs(i+2) = 1;

		// stabilize A
		for (j = 0; j <= i; j++)
		  {
		    for (k = 0; k <= i; k++)
		      mat(j,k) += b(0,j) * b(0,j) + b(1,j) * b(1,k);
		    rhs(j) += b(1,j);
		  }
	    
		CholeskyFactors<double> inv(mat);
		inv.Mult (rhs, sol);

		// tout << "sol = " << sol << endl;

		for (j = 0; j <= i; j++)
		  coefs[j][i] = sol(j);
	      }
	

	    /*
	    // 	    // testing
	    ofstream out ("vext.dat");
	    ofstream dout ("dvext.dat");
	    for (double x = 0; x <= 1; x += 0.01)
	      {
		out << x << " ";
		for (i = 1; i <= 10; i++)
		  out << Calc (i, x) << " ";
		out << endl;
		
		dout << x << " ";
		for (i = 1; i <= 10; i++)
		  dout << CalcDeriv (i, x) << " ";
		dout << endl;
	      }
	    */
	    initialized = 1;
	  }
	catch (const Exception& e)
	  {
	    cout << "Caught in VertexExtensionOptimal:" << e.What() << endl;
	  }
      }
}

#endif

  /*

  template <int DIM>
  VertexExtensionOptimal<DIM> :: VertexExtensionOptimal ()
  {
    if (!initialized)
      {

	ofstream tout;
	if (DIM == 3)
	  tout.open("test3d.out");
	else
	  tout.open("test2d.out");

	int i, j, k;
	Matrix<double> a(SIZE, SIZE);
	Matrix<double> b(2, SIZE);


	ofstream jacout ("jacobi.out");
	Vector<> jac(11);
	for (double x = -1; x <= 1; x += 0.01)
	  {
	    JacobiPolynomial (10, x, 1, -1, jac);
	    for (int j = 0; j <= 10; j++)
	      jacout << " " << jac(j) << " ";
	    jacout << endl;
	  }

	try
	  {
	    const IntegrationRule & ir =
	      GetIntegrationRules().SelectIntegrationRule (ET_SEGM, 2*SIZE+2);

	    Vector<double> diff_leg(SIZE); 
	    Vector<double> help(SIZE); 

	    a = 0;
	    b = 0;

	    for (i = 0; i < ir.GetNIP(); i++)
	      {
		double x = ir.GetIP(i)(0);

		JacobiPolynomial (SIZE-2, 2*x-1, 2, 0, help);
		diff_leg(0) = 0;
		for (j = 1; j < SIZE; j++)
		  diff_leg(j) = 0.5 * (j+1) * 2.0 * help(j-1);

		double fac = ir.GetIP(i).Weight();
		for (j = 0; j < DIM-1; j++)
		  fac *= (1-x);
	    
		for (j = 0; j < SIZE; j++)
		  for (k = 0; k < SIZE; k++)
		    a(j,k) += fac * diff_leg(j) * diff_leg(k);
	      }

	    JacobiPolynomial (SIZE-1, -1.0, 1, -1, help);
	    for (j = 0; j < SIZE; j++)
	      b(0,j) = help(j);
	    JacobiPolynomial (SIZE-1, 1.0, 1, -1, help);
	    for (j = 0; j < SIZE; j++)
	      b(1,j) = help(j);


	    tout << "mata = " << a << endl;
	    tout << "matb = " << b << endl;

	    tout << "diag (A), 1/diag(A),   bi*bi/aii " << endl;
	    for (j = 1; j < SIZE; j++)
	      tout << a(j,j) << ",    " << 1.0/a(j,j) << "    " << b(1,j)*b(1,j)/a(j,j) << endl;

	    // compute coefficients
	    for (i = 1; i < SIZE; i++)
	      {
		Matrix<double> mat(i+3, i+3);
		Vector<double> rhs(i+3), sol(i+3);
		mat = 0;
		for (j = 0; j <= i; j++)
		  for (k = 0; k <= i; k++)
		    mat(j,k) = a(j,k);

		for (j = 0; j <= i; j++)
		  {
		    mat(i+1, j) = mat(j,i+1) = b(0,j);
		    mat(i+2, j) = mat(j,i+2) = b(1,j);
		  }
	    
		// stabilize A
		rhs = 0;
		rhs(i+2) = 1;

		tout << "mat,1 = " << mat << endl;
		tout << "rhs,1 = " << rhs << endl;

		for (j = 0; j <= i; j++)
		  {
		    for (k = 0; k <= i; k++)
		      mat(j,k) += b(0,j) * b(0,k) + b(1,j) * b(1,k);
		    rhs(j) += b(1,j);
		  }
	    
		tout << "mat,2 = " << mat << endl;
		tout << "rhs,2 = " << rhs << endl;

		CholeskyFactors<double> inv(mat);
		inv.Mult (rhs, sol);

		tout << "sol = " << sol << endl;

		for (j = 0; j <= i; j++)
		  coefs[j][i] = sol(j);

		double schur = 0;
		for (j = 1; j <= i; j++)
		  schur += b(1,j) * b(1,j) / a(j,j);
		double lam = 1.0 / schur;

		tout << "schur = " << schur << ", lam = " << lam << endl;

		sol(0) = 0;
		for (j = 1; j <= i; j++)
		  sol(j) = (2*j+1) / ( (j+1.0)*(i*i+2*i));
		tout << "sol,2 = " << sol << endl;
	      }
	

	    // testing
	    ofstream out ("vext.dat");
	    ofstream dout ("dvext.dat");
	    for (double x = 0; x <= 1.0001; x += 0.01)
	      {
		out << x << " ";
		for (i = 1; i <= 10; i++)
		  out << Calc (i, x) << " ";
		out << endl;
		
		dout << x << " ";
		for (i = 1; i <= 10; i++)
		  dout << CalcDeriv (i, x) << " ";
		dout << endl;
	      }
	    
	    initialized = 1;
	  }
	catch (const Exception& e)
	  {
	    cout << "Caught in VertexExtensionOptimal:" << e.What() << endl;
	  }
      }
  }
  */



  /*
  template <int DIM>
  VertexExtensionOptimal<DIM> :: VertexExtensionOptimal ()
  {
    if (!initialized)
      {

	ofstream tout;
	if (DIM == 3)
	  tout.open("test3d.out");
	else
	  tout.open("test2d.out");

	int i, j, k;
	Matrix<double> a(SIZE, SIZE);
	Matrix<double> b(2, SIZE);


	ofstream jacout ("jacobi.out");
	Vector<> jac(11);
	for (double x = -1; x <= 1; x += 0.01)
	  {
	    JacobiPolynomial (10, x, 1, -1, jac);
	    for (int j = 0; j <= 10; j++)
	      jacout << " " << jac(j) << " ";
	    jacout << endl;
	  }

	try
	  {
	    const IntegrationRule & ir =
	      GetIntegrationRules().SelectIntegrationRule (ET_SEGM, 2*SIZE+4);

	    Vector<double> diff_leg(SIZE); 
	    Vector<double> help(SIZE); 
	    Vector<double> help2(SIZE); 

	    a = 0;
	    b = 0;

	    for (i = 0; i < ir.GetNIP(); i++)
	      {
		double x = ir.GetIP(i)(0);
		double hx = 2*x-1;

		JacobiPolynomial (SIZE-1, hx, 1, 1, help);
		JacobiPolynomial (SIZE-1, hx, 2, 2, help2);

		diff_leg(0) = 1;
		for (j = 1; j < SIZE; j++)
		  diff_leg(j) = help(j) + (1+hx)*0.5*(j+3)*help2(j-1);

		
		double fac = ir.GetIP(i).Weight();
		for (j = 0; j < DIM-1; j++)
		  fac *= (1-x);
	    
		for (j = 0; j < SIZE; j++)
		  for (k = 0; k < SIZE; k++)
		    a(j,k) += fac * diff_leg(j) * diff_leg(k);
	      }

	    JacobiPolynomial (SIZE-1, -1.0, 1, -1, help);
	    for (j = 0; j < SIZE; j++)
	      b(0,j) = help(j);
	    JacobiPolynomial (SIZE-1, 1.0, 1, -1, help);
	    for (j = 0; j < SIZE; j++)
	      b(1,j) = help(j);


	    tout << "mata = " << a << endl;
	    tout << "matb = " << b << endl;

	    tout << "diag (A), 1/diag(A),   bi*bi/aii " << endl;
	    for (j = 1; j < SIZE; j++)
	      tout << a(j,j) << ",    " << 1.0/a(j,j) << "    " << b(1,j)*b(1,j)/a(j,j) << endl;

	    // compute coefficients
	    for (i = 1; i < SIZE; i++)
	      {
		Matrix<double> mat(i+3, i+3);
		Vector<double> rhs(i+3), sol(i+3);
		mat = 0;
		for (j = 0; j <= i; j++)
		  for (k = 0; k <= i; k++)
		    mat(j,k) = a(j,k);

		for (j = 0; j <= i; j++)
		  {
		    mat(i+1, j) = mat(j,i+1) = b(0,j);
		    mat(i+2, j) = mat(j,i+2) = b(1,j);
		  }
	    
		// stabilize A
		rhs = 0;
		rhs(i+2) = 1;

		tout << "mat,1 = " << mat << endl;
		tout << "rhs,1 = " << rhs << endl;

		for (j = 0; j <= i; j++)
		  {
		    for (k = 0; k <= i; k++)
		      mat(j,k) += b(0,j) * b(0,k) + b(1,j) * b(1,k);
		    rhs(j) += b(1,j);
		  }
	    
		tout << "mat,2 = " << mat << endl;
		tout << "rhs,2 = " << rhs << endl;

		CholeskyFactors<double> inv(mat);
		inv.Mult (rhs, sol);

		tout << "sol = " << sol << endl;

		for (j = 0; j <= i; j++)
		  coefs[j][i] = sol(j);

		double schur = 0;
		for (j = 1; j <= i; j++)
		  schur += b(1,j) * b(1,j) / a(j,j);
		double lam = 1.0 / schur;

		tout << "schur = " << schur << ", lam = " << lam << endl;

		sol(0) = 0;
		for (j = 1; j <= i; j++)
		  sol(j) = (2*j+1) / ( (j+1.0)*(i*i+2*i));
		tout << "sol,2 = " << sol << endl;
	      }
	

	    // testing
	    ofstream out ("vext.dat");
	    ofstream dout ("dvext.dat");
	    for (double x = 0; x <= 1.0001; x += 0.01)
	      {
		out << x << " ";
		for (i = 1; i <= 10; i++)
		  out << Calc (i, x) << " ";
		out << endl;
		
		dout << x << " ";
		for (i = 1; i <= 10; i++)
		  dout << CalcDeriv (i, x) << " ";
		dout << endl;
	      }
	    
	    initialized = 1;
	  }
	catch (const Exception& e)
	  {
	    cout << "Caught in VertexExtensionOptimal:" << e.What() << endl;
	  }
      }
  }
  */



  /*
  VertexExtensionOptimal<2> dummy_vextopt2;
  VertexExtensionOptimal<3> dummy_vextopt3;
  */



  Array<ConvertJacobi::d2*> ConvertJacobi::coefs_reducealpha;
  Array<ConvertJacobi::d2*> ConvertJacobi::coefs_reducealphafac;
  Array<ConvertJacobi::d2*> ConvertJacobi::coefs_c;
  Array<ConvertJacobi::d2*> ConvertJacobi::coefs_d;
  Array<ConvertJacobi::d2*> ConvertJacobi::coefs_e;


  ConvertJacobi :: ConvertJacobi()
  {
    int n = 200;

    coefs_reducealpha.SetSize(n);
    for (int al = 0; al < n; al++)
      {
        coefs_reducealpha[al] = new d2[n];
        for (int i = 0; i < n; i++)
          {
            coefs_reducealpha[al][i][0] = i / double(i+al);
            coefs_reducealpha[al][i][1] = (2*i+al)/double(i+al);
          }
      }


    coefs_reducealphafac.SetSize(n);
    for (int al = 0; al < n; al++)
      {
        coefs_reducealphafac[al] = new d2[n];
        for (int i = 1; i < n; i++)
          {
            coefs_reducealphafac[al][i][0] = -i / (i+0.5*al) / 2;
            coefs_reducealphafac[al][i][1] = (i+al) / (i+0.5*al) / 2;
          }
      }


    coefs_c.SetSize(n);
    for (int al = 0; al < n; al++)
      {
        coefs_c[al] = new d2[n];
        for (int i = 0; i < n; i++)
          {
            coefs_c[al][i][0] = c(i, al);
            coefs_c[al][i][1] = hatc(i,al);
          }
      }

    coefs_d.SetSize(n);
    for (int al = 0; al < n; al++)
      {
        coefs_d[al] = new d2[n];
        for (int i = 0; i < n; i++)
          {
            coefs_d[al][i][0] = d(i, al);
            coefs_d[al][i][1] = hatd(i,al);
          }
      }

    coefs_e.SetSize(n);
    for (int al = 0; al < n; al++)
      {
        coefs_e[al] = new d2[n];
        for (int i = 0; i < n; i++)
          {
            coefs_e[al][i][0] = e(i, al);
            coefs_e[al][i][1] = hate(i,al);
          }
      }


  }

  ConvertJacobi :: ~ConvertJacobi ()
  {
    for (int i = 0; i < coefs_reducealpha.Size(); i++)
      {
        delete [] coefs_reducealpha[i];
        delete [] coefs_reducealphafac[i];
        delete [] coefs_c[i];
        delete [] coefs_d[i];
        delete [] coefs_e[i];
      }
  }

  ConvertJacobi init_convjac;





  double IntegratedLegendreMonomialExt :: coefs[SIZE][2];



  namespace recursive_pol_cpp
  {
    class Init
    {
      public:
      Init ()
      {
	IntegratedLegendreMonomialExt :: CalcCoeffs();


        /*
        int n = 10;
        Vector<double> ui(n+2), jac(n+1), jacm(n+2);
        double x = 0.5;
        JacobiPolynomial (n, x, 5, 0, jac);
        JacobiPolynomial (n+1, x, 4, 0, jacm);
        ui = 0; ui(5) = 1;
        
        cout << "jac_i = " << jac << endl;
        cout << "jac(x) = "<< InnerProduct (ui.Range(0,n+1), jac) * (1-x) << endl;

        ConvertJacobiReduceAlphaFactor (n, 4, ui);
        cout << " =?= jac(x) = "<< InnerProduct (ui, jacm) << endl;
        */
      }
    };

    Init init;
  }

  void testit(FlatVector<> v)
  {
    ConvertJacobi::ReduceAlphaFactor (3, 10, v);
  }

  void testit210(FlatVector<> v)
  {
    // Vector<> v(100);
    ConvertJacobi::ReduceAlphaFactor (10, 10, v);
  }

  void testit2a(FlatVector<> v)
  {
    // Vector<> v(100);
    for (int i = 0; i < 3; i++)
      ConvertJacobi::ReduceAlphaFactor (5-i, 2*i, v);
  }

  template<int N>
  void testit2brec (FlatVector<> v)
  {
    ConvertJacobi::ReduceAlphaFactor (5-N, 2*N, v);
    testit2brec<N-1> (v);
  }
  template<>
  void testit2brec<0> (FlatVector<> v)
  { ; }


  void testit2bOrig(FlatVector<> v)
  {
    testit2brec<5> (v);
  }
  /*
  void testit2b1(FlatVector<> v)
  {
    testit2brec<1> (v);
  }
  void testit2b2(FlatVector<> v)
  {
    testit2brec<2> (v);
  }
  */
  void testit2b3(FlatVector<> v)
  {
    testit2brec<3> (v);
  }

  
  void testit3(int n, FlatMatrix<> m)
  {
    ConvertJacobi::TriangularReduceFactor (1, 10, 8, m);
  }
  void testit3b(int n, FlatMatrix<> m)
  {
    ConvertJacobi::TriangularReduceFactor (2, 10, 8, m);
  }



  void testit4(int n, int al, int be)
  {
    Matrix<> m(n, n);
    Vector<> v(n);
    ConvertJacobi::TriangularReduceBeta (3, al, be, m, v);
  }




  void testit5(FlatVector<> v)
  {
    TReduceAlpha<5,8>::Do(v);
  }


}
