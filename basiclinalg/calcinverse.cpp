#include <bla.hpp>

namespace ngbla
{


  void test()
  {
    Matrix<Complex> a, b, c;
 
    LapackMult (a, b, c);
    LapackMult (Trans(a), b, c);
    LapackMult (a, Trans(b), c);
  }


  
  inline double abs (double a)
  {
    return std::fabs(a);
  }

  inline double abs (Complex a)
  {
    return std::abs(a);
  }

#ifdef USE_GMP

  inline double abs (mpq_class a)
  {
    return std::fabs(a.get_d());
  }

  void CalcInverse (mpq_class a, mpq_class & ia)
  {
    ia = 1 / a;
  }
#endif


  template <int N, int N2, typename SCAL>
  inline double abs (Mat<N,N2,SCAL> & m)
  {
    double sum = 0;
    for (int i = 0; i < N; i++)
      sum += abs(m(i,i));
    return sum;
  }


  template <class T2>
  // void T_CalcInverse (const FlatMatrix<T> m, FlatMatrix<T2> inv)
  void T_CalcInverse (FlatMatrix<T2> inv)
  {
    // Gauss - Jordan - algorithm
    // Algorithm of Stoer, Einf. i. d. Num. Math, S 145
    // int n = m.Height();
    int n = inv.Height();

    ngstd::Array<int> p(n);   // pivot-permutation

    // inv = m;

    for (int j = 0; j < n; j++) p[j] = j;
    
    for (int j = 0; j < n; j++)
      {
	// pivot search
	double maxval = abs(inv(j,j));
	int r = j;

	for (int i = j+1; i < n; i++)
	  if (abs (inv(j, i)) > maxval)
	    {
	      r = i;
	      maxval = abs (inv(j, i));
	    }
      
	if (maxval < 1e-20)
	  {
	    throw Exception ("Inverse matrix: Matrix singular");
	  }

	// exchange rows
	if (r > j)
	  {
	    for (int k = 0; k < n; k++)
	      swap (inv(k, j), inv(k, r));
	    swap (p[j], p[r]);
	  }
      

	// transformation
	

	T2 hr;
	CalcInverse (inv(j,j), hr);
	for (int i = 0; i < n; i++)
	  {
	    T2 h = hr * inv(j, i);
	    inv(j, i) = h;
	  }
	inv(j,j) = hr;


	for (int k = 0; k < n; k++)
	  if (k != j)
	    {
	      T2 help = inv(n*k+j);
	      for (int i = 0; i < j; i++)
		{
		  T2 h = help * inv(n*j+i); 
		  inv(n*k+i) -= h;
		}
	      for (int i = j+1; i < n; i++)
		{
		  T2 h = help * inv(n*j+i); 
		  inv(n*k+i) -= h;
		}

	      T2 h = inv(k,j) * hr;   
	      inv(k,j) = -h;
	    }
      }

    // row exchange
  
    Vector<T2> hv(n);
    for (int i = 0; i < n; i++)
      {
	for (int k = 0; k < n; k++) hv(p[k]) = inv(k, i);
	for (int k = 0; k < n; k++) inv(k, i) = hv(k);
      }
  }











  void CalcSchurComplement (const FlatMatrix<double> a, 
			    FlatMatrix<double> s,
			    const BitArray & used,
			    LocalHeap & lh)
  {
    // cout << "Calc SchurComplement" << endl;

    HeapReset hr(lh);

    int n = a.Height(), n_used = 0, n_unused = 0;
    for (int i = 0; i < n; i++)
      if (used[i]) n_used++; else n_unused++;

    FlatMatrix<> b1(n_unused, n_used, lh);
    FlatMatrix<> hb1(n_unused, n_used, lh);
    FlatMatrix<> b2(n_used, n_unused, lh);
    FlatMatrix<> c(n_unused, n_unused, lh);
    FlatMatrix<> hc(n_unused, n_unused, lh);

    //(*testout) << "used = " << endl << used << endl;
    int cnt_usedi = 0;
    for (int i = 0; i < n; i++)
      {
	int cnt_usedj = 0;
	for (int j = 0; j < n; j++)
	  {
	    if (used[i] && used[j])
	      s(cnt_usedi, cnt_usedj) = a(i,j);
	    else if (!used[i] && used[j])
	      b1(i-cnt_usedi, cnt_usedj) = a(i,j);
	    else if (used[i] && !used[j])
	      b2(cnt_usedi, j-cnt_usedj) = a(i,j);
	    else 
	      c(i-cnt_usedi, j-cnt_usedj) = a(i,j);
		
	    if (used[j]) cnt_usedj++;
	  }
	if (used[i]) cnt_usedi++;
      }
    
    /*
    (*testout) << "Schur, a = " << endl << a << endl;
    (*testout) << "Schur, s = " << endl << s << endl;
    (*testout) << "Schur, b1 = " << endl << b1 << endl;
    (*testout) << "Schur, b2 = " << endl << b2 << endl;
    (*testout) << "Schur, c = " << endl << c << endl;
    */

    LapackInverse (c);

    // LapackMult (c, b1, hb1);
    hb1 = c * b1 | Lapack;

    // LapackMultAddAB (b2, hb1, -1, s);
    s -= b2 * hb1 | Lapack;

    // (*testout) << "Schur, cinv = " << endl << c << endl;
    // (*testout) << "Schur = " << endl << s << endl;
  }



#ifdef LAPACK
  template <>
  void CalcInverse (FlatMatrix<Complex> inv)
  {
    LapackInverse (inv);
  }
#endif



  template <class T2>
  extern void CalcInverse (FlatMatrix<T2> inv)
  {
    T_CalcInverse (inv);
  }

  extern NGS_DLL_HEADER void CalcInverse (FlatMatrix<double> inv)
  {
    T_CalcInverse (inv);
  }


#ifdef USE_GMP
  template void CalcInverse (FlatMatrix<mpq_class> inv);
#endif


  template void CalcInverse (FlatMatrix<Mat<1,1,double> > inv);
  template void CalcInverse (FlatMatrix<Mat<2,2,double> > inv);
  template void CalcInverse (FlatMatrix<Mat<3,3,double> > inv);
  template void CalcInverse (FlatMatrix<Mat<4,4,double> > inv);
  template void CalcInverse (FlatMatrix<Mat<5,5,double> > inv);
  template void CalcInverse (FlatMatrix<Mat<6,6,double> > inv);
  template void CalcInverse (FlatMatrix<Mat<7,7,double> > inv);
  template void CalcInverse (FlatMatrix<Mat<8,8,double> > inv);


#ifndef LAPACK
  template void CalcInverse (FlatMatrix<Complex> inv);
#endif

  template void CalcInverse (FlatMatrix<Mat<1,1,Complex> > inv);
  template void CalcInverse (FlatMatrix<Mat<2,2,Complex> > inv);
  template void CalcInverse (FlatMatrix<Mat<3,3,Complex> > inv);
  template void CalcInverse (FlatMatrix<Mat<4,4,Complex> > inv);
  template void CalcInverse (FlatMatrix<Mat<5,5,Complex> > inv);
  template void CalcInverse (FlatMatrix<Mat<6,6,Complex> > inv);
  template void CalcInverse (FlatMatrix<Mat<7,7,Complex> > inv);
  template void CalcInverse (FlatMatrix<Mat<8,8,Complex> > inv);
}
