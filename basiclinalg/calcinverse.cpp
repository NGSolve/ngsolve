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


  template <int N, int N2, typename SCAL>
  inline double abs (Mat<N,N2,SCAL> & m)
  {
    double sum = 0;
    for (int i = 0; i < N; i++)
      sum += abs(m(i,i));
    return sum;
  }


  template <class T2>
  void T_CalcInverse (SliceMatrix<T2> inv)
  {
    // static Timer t("CalcInverse");
    // RegionTimer reg(t);

    // Gauss - Jordan - algorithm
    // Algorithm of Stoer, Einf. i. d. Num. Math, S 145
    // int n = m.Height();

    int n = inv.Height();

    ArrayMem<int,100> p(n);   // pivot-permutation
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
	      T2 help = inv(k,j);
	      T2 h = help * hr;   

	      for (int i = 0; i < n; i++)
		{
		  T2 h = help * inv(j,i); 
		  inv(k,i) -= h;
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



  template <>
  void CalcInverse (SliceMatrix<Complex> inv, INVERSE_LIB il)
  {
    // crossover vs lapack: ~450 Accelerate (M4), ~60 MKL (AVX512)
    if (il == INVERSE_LIB::INV_CHOOSE)
      il = inv.Height() >= 150 ? INVERSE_LIB::INV_LAPACK : INVERSE_LIB::INV_NGBLA_LU;

    if (il == INVERSE_LIB::INV_NGBLA_LU)
      {
        ArrayMem<int,100> p(inv.Height());
        CalcLU (inv, p);
        InverseFromLU (inv, p);
      }
#ifdef LAPACK
    else if (il == INVERSE_LIB::INV_LAPACK)
      LapackInverse (inv);
#endif
    else
      T_CalcInverse (inv);
  }

  /*
  extern NGS_DLL_HEADER
  void CalcInverse (FlatMatrix<float> inv, INVERSE_LIB il)
  {
    T_CalcInverse (inv);    
  }
  */
  template void CalcInverse (SliceMatrix<float> inv, INVERSE_LIB il);  


  template <class T2>
  extern void CalcInverse (SliceMatrix<T2> inv, INVERSE_LIB il)
  {
    T_CalcInverse (inv);
  }

  extern NGS_DLL_HEADER void CalcInverse (SliceMatrix<double> inv, INVERSE_LIB il)
  {
    if (il == INVERSE_LIB::INV_CHOOSE)
      il = inv.Height() >= 100 ? INVERSE_LIB::INV_LAPACK : INVERSE_LIB::INV_NGBLA_LU;

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
        InverseFromLU (inv, p);

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


  



#define INST(N) \
  template void CalcInverse (SliceMatrix<Mat<N,N,double> > inv, INVERSE_LIB il); \
  template void CalcInverse (SliceMatrix<Mat<N,N,Complex> > inv, INVERSE_LIB il);
  NGS_FOR_SYS_DIMS(INST)
#undef INST

}
