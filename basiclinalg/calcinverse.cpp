#include <bla.hpp>

namespace ngbla
{
  using namespace ngbla;
  
  inline double abs (double a)
  {
    return std::fabs(a);
  }

  inline double abs (Complex a)
  {
    return std::abs(a);
  }

  template <int N, int N2, typename SCAL>
  inline double abs (Mat<N,N2,SCAL> & m)
  {
    double sum = 0;
    for (int i = 0; i < N; i++)
      sum += abs(m(i,i));
    return sum;
  }


  // #define COLWISE
#ifdef COLWISE 

  // slow memory access

  template <class T, class T2>
  void CalcInverse (const FlatMatrix<T> m, FlatMatrix<T2> inv)
  {
    //	  Gauss - Jordan - algorithm
  
    int n = m.Height();

    T hr;

    ngstd::ARRAY<int> p(n);   // pivot-permutation
    Vector<T> hv(n);
   
    inv = m;

    // Algorithm of Stoer, Einf. i. d. Num. Math, S 145

    for (int j = 0; j < n; j++)
      p[j] = j;
    
    for (int j = 0; j < n; j++)
      {
	// pivot search

	double maxval = abs(inv(j,j));
	int r = j;

	for (int i = j+1; i < n ;i++)
	  if (abs (inv(i, j)) > maxval)
	    {
	      r = i;
	      maxval = abs (inv(i, j));
	    }
      
	if (maxval < 1e-20)
	  {
	    throw Exception ("Inverse matrix: Matrix singular");
	  }

	// exchange rows
	if (r > j)
	  {
	    for (int k = 0; k < n; k++)
	      swap (inv(j,k), inv(r,k));
	    swap (p[j], p[r]);
	  }
      

	// transformation

	CalcInverse (inv(j,j), hr);
	for (int i = 0; i < n; i++)
	  {
	    T  h = inv(i,j) * hr;
	    inv(i,j) = h;
	  }
	inv(j,j) = hr;


	for (int k = 0; k < n; k++)
	  if (k != j)
	    {
	      
	      /*
	      for (i = 0; i < n; i++)
		if (i != j)
		  {
		    T h = inv(n*i+j) * inv(n*j+k);
		    inv(n*i+k) -= h;
		  }
	      */
	      T help = inv(n*j+k);
	      for (int i = 0; i < j; i++)
		{
		  T h = inv(n*i+j) * help; // inv(n*j+k);
		  inv(n*i+k) -= h;
		}
	      for (int i = j+1; i < n; i++)
		{
		  T h = inv(n*i+j) * help; // inv(n*j+k);
		  inv(n*i+k) -= h;
		}

	      T h = hr * inv(j,k);   
	      inv(j,k) = -h;
	    }
      }

    // col exchange
  
    for (int i = 0; i < n; i++)
      {
	for (int k = 0; k < n; k++)
	  hv(p[k]) = inv(i, k);
	for (int k = 0; k < n; k++)
	  inv(i, k) = hv(k);
      }
  }

#else

  template <class T, class T2>
  void CalcInverse (const FlatMatrix<T> m, FlatMatrix<T2> inv)
  {
    // Gauss - Jordan - algorithm
    // Algorithm of Stoer, Einf. i. d. Num. Math, S 145
    int n = m.Height();

    ngstd::ARRAY<int> p(n);   // pivot-permutation

    T hr;
    inv = m;

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

	CalcInverse (inv(j,j), hr);
	for (int i = 0; i < n; i++)
	  {
	    T  h = hr * inv(j, i);
	    inv(j, i) = h;
	  }
	inv(j,j) = hr;


	for (int k = 0; k < n; k++)
	  if (k != j)
	    {
	      /*
	      for (int i = 0; i < n; i++)
		if (i != j)
		  {
		    T h = inv(n*j+i) * inv(n*k+j);  // exchange both ?
		    inv(n*k+i) -= h;
		  }
	      */

	      T help = inv(n*k+j);
	      for (int i = 0; i < j; i++)
		{
		  T h = help * inv(n*j+i); 
		  inv(n*k+i) -= h;
		}
	      for (int i = j+1; i < n; i++)
		{
		  T h = help * inv(n*j+i); 
		  inv(n*k+i) -= h;
		}

	      T h = inv(k,j) * hr;   
	      inv(k,j) = -h;
	    }
      }

    // row exchange
  
    Vector<T> hv(n);
    for (int i = 0; i < n; i++)
      {
	for (int k = 0; k < n; k++) hv(p[k]) = inv(k, i);
	for (int k = 0; k < n; k++) inv(k, i) = hv(k);
      }
  }



#endif













  void CalcSchurComplement (const FlatMatrix<double> a, 
			    FlatMatrix<double> s,
			    const BitArray & used,
			    LocalHeap & lh)
  {
    // cout << "Calc SchurComplement" << endl;

    void * heapp = lh.GetPointer();

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

#ifdef LAPACK
    LapackInverse (c);
    LapackMultAB (c, b1, hb1);
    // hb1 = c * b1;
    // s -= b2 * hb1;
    LapackMultAddAB (b2, hb1, -1, s);
#else
    FlatMatrix<> c_inv(n_unused, n_unused, lh);
    CalcInverse (c, c_inv);
    c = c_inv;

    hb1 = c * b1;
    s -= b2 * hb1;
#endif

    lh.CleanUp (heapp);

    // (*testout) << "Schur, cinv = " << endl << c << endl;


    // (*testout) << "Schur = " << endl << s << endl;
  }



#ifdef LAPACK
  template <>
  void CalcInverse (const FlatMatrix<Complex> m, 
                    FlatMatrix<Complex> inv)
  {
    inv = m;
    LapackInverse (inv);
  }
#endif




  template void CalcInverse (const FlatMatrix<double> m, 
			     FlatMatrix<double> inv);
  template void CalcInverse (const FlatMatrix<Mat<1,1,double> > m, 
			     FlatMatrix<Mat<1,1,double> > inv);

  template void CalcInverse (const FlatMatrix<Mat<2,2,double> > m, 
			     FlatMatrix<Mat<2,2,double> > inv);
  template void CalcInverse (const FlatMatrix<Mat<3,3,double> > m, 
			     FlatMatrix<Mat<3,3,double> > inv);
  template void CalcInverse (const FlatMatrix<Mat<4,4,double> > m, 
			     FlatMatrix<Mat<4,4,double> > inv);
  template void CalcInverse (const FlatMatrix<Mat<5,5,double> > m, 
			     FlatMatrix<Mat<5,5,double> > inv);
  template void CalcInverse (const FlatMatrix<Mat<6,6,double> > m, 
			     FlatMatrix<Mat<6,6,double> > inv);
  template void CalcInverse (const FlatMatrix<Mat<7,7,double> > m, 
			     FlatMatrix<Mat<7,7,double> > inv);
  template void CalcInverse (const FlatMatrix<Mat<8,8,double> > m, 
			     FlatMatrix<Mat<8,8,double> > inv);


#ifndef LAPACK
  template void CalcInverse (const FlatMatrix<Complex> m, 
                             FlatMatrix<Complex> inv);
#endif

  template void CalcInverse (const FlatMatrix<Mat<1,1,Complex> > m, 
			     FlatMatrix<Mat<1,1,Complex> > inv);
  template void CalcInverse (const FlatMatrix<Mat<2,2,Complex> > m, 
			     FlatMatrix<Mat<2,2,Complex> > inv);
  template void CalcInverse (const FlatMatrix<Mat<3,3,Complex> > m, 
			     FlatMatrix<Mat<3,3,Complex> > inv);
  template void CalcInverse (const FlatMatrix<Mat<4,4,Complex> > m, 
			     FlatMatrix<Mat<4,4,Complex> > inv);
  template void CalcInverse (const FlatMatrix<Mat<5,5,Complex> > m, 
			     FlatMatrix<Mat<5,5,Complex> > inv);
  template void CalcInverse (const FlatMatrix<Mat<6,6,Complex> > m, 
			     FlatMatrix<Mat<6,6,Complex> > inv);
  template void CalcInverse (const FlatMatrix<Mat<7,7,Complex> > m, 
			     FlatMatrix<Mat<7,7,Complex> > inv);
  template void CalcInverse (const FlatMatrix<Mat<8,8,Complex> > m, 
			     FlatMatrix<Mat<8,8,Complex> > inv);
}
