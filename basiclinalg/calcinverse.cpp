#include <bla.hpp>

namespace ngbla
{
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
  void T_CalcInverse (FlatMatrix<T2> inv)
  {
    // static Timer t("CalcInverse");
    // RegionTimer reg(t);

    // Gauss - Jordan - algorithm
    // Algorithm of Stoer, Einf. i. d. Num. Math, S 145
    // int n = m.Height();

    int n = inv.Height();

    ngstd::ArrayMem<int,100> p(n);   // pivot-permutation
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
      
        double rest = 0.0;
        for (int i = j+1; i < n; i++)
          rest += abs(inv(r, i));
	if (maxval < 1e-20*rest)
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
	      T2 h = help * hr;   

	      for (int i = 0; i < n; i++)
		{
		  T2 h = help * inv(n*j+i); 
		  inv(n*k+i) -= h;
		}

	      inv(k,j) = -h;
	    }
      }

    // row exchange
  
    VectorMem<100,T2> hv(n);
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
    if (s.Height() == 0) return;
    if (s.Height() == a.Height())
      {
        s = a;
        return;
      }

    HeapReset hr(lh);

    int n = a.Height();
    Array<int> used_dofs(n, lh);
    Array<int> unused_dofs(n, lh);
    used_dofs.SetSize(0);
    unused_dofs.SetSize(0);
    for (int i = 0; i < n; i++)
      if (used[i])
        used_dofs.Append(i);
      else
        unused_dofs.Append(i);

    s = a.Rows(used_dofs).Cols(used_dofs);
    FlatMatrix<> b1 = a.Rows(unused_dofs).Cols(used_dofs) | lh;
    FlatMatrix<> b2 = a.Rows(used_dofs).Cols(unused_dofs) | lh;
    FlatMatrix<> c = a.Rows(unused_dofs).Cols(unused_dofs) | lh;
    FlatMatrix<> hb1 (b1.Height(), b1.Width(), lh);

    if (n > 10)
      {
        LapackInverse (c);
        hb1 = c * b1 | Lapack;
        s -= b2 * hb1 | Lapack;
      }
    else
      {
        CalcInverse (c);
        hb1 = c * b1;
        s -= b2 * hb1;
      }
  }



#ifdef LAPACK
  template <>
  void CalcInverse (FlatMatrix<Complex> inv, INVERSE_LIB il)
  {
    LapackInverse (inv);
  }
#endif



  template <class T2>
  extern void CalcInverse (FlatMatrix<T2> inv, INVERSE_LIB il)
  {
    T_CalcInverse (inv);
  }

  extern NGS_DLL_HEADER void CalcInverse (FlatMatrix<double> inv, INVERSE_LIB il)
  {
    if (il == INVERSE_LIB::INV_CHOOSE)
      il = inv.Height() >= 100 ? INVERSE_LIB::INV_LAPACK : INVERSE_LIB::INV_NGBLA;

    if (il == INVERSE_LIB::INV_NGBLA_QR)
      {
        QRFactorizationInPlace (inv);
        InverseFromQR (inv);
      }
    
    else if (il == INVERSE_LIB::INV_NGBLA_LU)
      {
        // Matrix save = inv;
        ArrayMem<int,100> p(inv.Height());
        CalcLU (inv, p);
        InverseFromLU (inv, p); // has issues with numerical stability

        /*
          // stable as lapack degetri (but not as stable as old CalcInverse)
        Matrix<double,ColMajor> X = Identity(inv.Height());
        SolveFromLU (inv, p, X);
        inv = X;
        */
        
        /*
        double err = L2Norm(inv*save-Identity(inv.Height()));
        if (err > 0.1)
          {
            cout << "inverse, err = " << err << endl;

            inv = save;
            T_CalcInverse (inv);
            double err2 = L2Norm(inv*save-Identity(inv.Height()));
            cout << "old inverse, err = " << err2 << endl;
          }
        */
      }
#ifdef LAPACK
    else if (il == INVERSE_LIB::INV_LAPACK)
      LapackInverse(inv);
#endif
    else
      T_CalcInverse (inv);
  }


#ifdef USE_GMP
  template void CalcInverse (FlatMatrix<mpq_class> inv);
#endif


#ifndef LAPACK
  template void CalcInverse (FlatMatrix<Complex> inv, INVERSE_LIB il);
#endif



#if MAX_SYS_DIM >= 1
  template void CalcInverse (FlatMatrix<Mat<1,1,double> > inv, INVERSE_LIB il);
  template void CalcInverse (FlatMatrix<Mat<1,1,Complex> > inv, INVERSE_LIB il);
#endif
#if MAX_SYS_DIM >= 2
  template void CalcInverse (FlatMatrix<Mat<2,2,double> > inv, INVERSE_LIB il);
  template void CalcInverse (FlatMatrix<Mat<2,2,Complex> > inv, INVERSE_LIB il);
#endif
#if MAX_SYS_DIM >= 3
  template void CalcInverse (FlatMatrix<Mat<3,3,double> > inv, INVERSE_LIB il);
  template void CalcInverse (FlatMatrix<Mat<3,3,Complex> > inv, INVERSE_LIB il);
#endif
#if MAX_SYS_DIM >= 4
  template void CalcInverse (FlatMatrix<Mat<4,4,double> > inv, INVERSE_LIB il);
  template void CalcInverse (FlatMatrix<Mat<4,4,Complex> > inv, INVERSE_LIB il);
#endif
#if MAX_SYS_DIM >= 5
  template void CalcInverse (FlatMatrix<Mat<5,5,Complex> > inv, INVERSE_LIB il);
  template void CalcInverse (FlatMatrix<Mat<5,5,double> > inv, INVERSE_LIB il);
#endif
#if MAX_SYS_DIM >= 6
  template void CalcInverse (FlatMatrix<Mat<6,6,double> > inv, INVERSE_LIB il);
  template void CalcInverse (FlatMatrix<Mat<6,6,Complex> > inv, INVERSE_LIB il);
#endif
#if MAX_SYS_DIM >= 7
  template void CalcInverse (FlatMatrix<Mat<7,7,double> > inv, INVERSE_LIB il);
  template void CalcInverse (FlatMatrix<Mat<7,7,Complex> > inv, INVERSE_LIB il);
#endif
#if MAX_SYS_DIM >= 8
  template void CalcInverse (FlatMatrix<Mat<8,8,double> > inv, INVERSE_LIB il);
  template void CalcInverse (FlatMatrix<Mat<8,8,Complex> > inv, INVERSE_LIB il);
#endif

}
