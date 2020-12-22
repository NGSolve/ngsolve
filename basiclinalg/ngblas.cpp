#define COMPILE_NGBLAS

#include <bla.hpp>


namespace ngbla
{

  

  int dgemm(char *transa, char *transb, integer *m, integer *
		  n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
		  doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
		  integer *ldc)
  {
    return dgemm_ (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
  }

  int zgemm(char *transa, char *transb, integer *m, integer *
		    n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, 
		    doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *
		    c__, integer *ldc)
  {
    return zgemm_ (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
  }

  int dger(integer *m, integer *n, doublereal *alpha,
                   doublereal *x, integer *incx, doublereal *y, integer *incy,
                   doublereal *a, integer *lda)
  {
    return dger_ (m, n, alpha, y, incx, x, incy, a, lda);
  }

  int dgetri(integer* n, double* a, integer* lda, integer* ipiv,
             double* hwork, integer* lwork, integer* info)
  {
    return dgetri_(n,a,lda,ipiv,hwork,lwork,info);
  }

  int dgetrf(integer* n, integer* m, double* a, integer* lda, integer* ipiv, integer* info)
  {
    return dgetrf_(n,m,a,lda,ipiv,info);
  }

#include "matkernel.hpp"
  // template <OPERATION OP>
  constexpr OPERATION AddOp(OPERATION OP) { return (OP == ADD || OP == SET) ? ADD : SUB; }


  /* ***************************** Copy Matrix *********************** */
  // copy matrix
  /*
    with AVX2 we get ~3 GF, independent of version 1 or version 2
    prefetch not yet tested 
  */


  /*
  template <int SW>
  inline void CopyMatrixIn (size_t h, size_t w,
                            double * ps, size_t dists,
                            SIMD<double,SW> * pd, size_t distd)
  {
    // constexpr int SW = SIMD<double>::Size();
    SIMD<mask64> mask(w % SW);

    for (size_t i = 0; i < h; i++, pd += distd, ps += dists)
      {
        size_t js = 0, jd=0;
        for ( ; js+SW <= w; js+=SW, jd++)
          pd[jd] = SIMD<double>(ps+js);
        SIMD<double>(ps+js, mask).Store((double*) (pd+jd));
        }
        }
  */
  template <int SW>
  inline void CopyMatrixIn (size_t h, size_t w,
                            double * ps, size_t dists,
                            SIMD<double,SW> * pd, size_t distd)
  {
    SIMD<mask64> mask(w % SW);
    
    for (size_t i = 0; i < h; i++, pd += distd, ps += dists)
      {
        auto ps2 = ps;
        auto pd2 = pd;
        
        size_t js = 0; 
        for ( ; js+SW <= w; js+=SW, pd2++, ps2+=SW)
          *pd2 = SIMD<double>(ps2);
        SIMD<double>(ps2, mask).Store((double*) (pd2), mask);
      }
  }
  


  
  /* ************************ Matrix * Vector ************************** */


  template <int SX>
  void MultMatVecShort (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    KernelMatVec<SX,SET> (y.Size(), &a(0), a.Dist(), &x(0), &y(0));
  }


  template <int SX>
  void MultAddMatVecShort (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    KernelAddMatVec<SX> (s, y.Size(), &a(0), a.Dist(), &x(0), &y(0));
  }


  NGS_DLL_HEADER void MultMatVec_intern (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    // constexpr int SW = SIMD<double>::Size();
    size_t h = y.Size();
    size_t w = x.Size();
    size_t i = 0;

    double * pa = &a(i,0);
    for ( ; i+8 <= h; i+=8, pa += 8*a.Dist())
      {
        // SIMD<double,4> sum1, sum2;
        // tie(sum1, sum2) = MatKernelScalAB<8,1> (w, pa, a.Dist(), &x(0), 0);
        auto [sum1, sum2] = MatKernelScalAB<8,1> (w, pa, a.Dist(), &x(0), 0);
        sum1.Store(&y(i));        
        sum2.Store(&y(i+4));        
      }
    
    if (i+4 <= h)
      {
        // SIMD<double,4> sum;
        // tie(sum) = MatKernelScalAB<4,1> (w, pa, a.Dist(), &x(0), 0);
        auto [sum] = MatKernelScalAB<4,1> (w, pa, a.Dist(), &x(0), 0);
        sum.Store(&y(i));
        i += 4;
        pa += 4*a.Dist();
      }

    if (i+2 <= h)
      {
        auto scal = MatKernelScalAB<2,1> (w, pa, a.Dist(), &x(0), 0);
        SIMD<double,2> sum(get<0>(scal), get<1>(scal));
        sum.Store(&y(i));
        i += 2;
        pa += 2*a.Dist();
      }

    if (i+1 <= h)
      {
        auto scal = MatKernelScalAB<1,1> (w, pa, a.Dist(), &x(0), 0);
        y(i) = get<0>(scal);
      }

  }


  NGS_DLL_HEADER void MultAddMatVec_intern (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    y += s * a.AddSize(y.Size(),x.Size()) * x;
  }

  pmult_matvec dispatch_matvec[];
  /*=
    {
      &MultMatVecShort<0>, &MultMatVecShort<1>, &MultMatVecShort<2>, &MultMatVecShort<3>,
      &MultMatVecShort<4>, &MultMatVecShort<5>, &MultMatVecShort<6>, &MultMatVecShort<7>,
      &MultMatVecShort<8>, &MultMatVecShort<9>, &MultMatVecShort<10>, &MultMatVecShort<11>,
      &MultMatVecShort<12>, &MultMatVecShort<13>, &MultMatVecShort<14>, &MultMatVecShort<15>,
      &MultMatVecShort<16>, &MultMatVecShort<17>, &MultMatVecShort<18>, &MultMatVecShort<19>,
      &MultMatVecShort<20>, &MultMatVecShort<21>, &MultMatVecShort<22>, &MultMatVecShort<23>,
      &MultMatVecShort<24>
    };
  */
  
  auto init_matvec = [] ()
  {
    Iterate<std::size(dispatch_matvec)-1> ([&] (auto i)
    { dispatch_matvec[i] = &MultMatVecShort<i>; });
    dispatch_matvec[std::size(dispatch_matvec)-1] = &MultMatVec_intern;
    return 1;
  }();
  
  
  /*
  template <template <int> typename FUNC, typename T>
  void InitDispatchArray (T * ap)
  {
    cout << "array size = " << std::size(*ap) << endl;
  }
  
  auto myinit = [] ()
  {
    cout << "init" << endl;
    InitDispatchArray<MultMatVecShort> (dispatch_matvec);
    return 0;
  };
  int dummy_myinit = myinit();
  */
  

 pmultadd_matvec dispatch_addmatvec[25] =
    {
      &MultAddMatVecShort<0>, &MultAddMatVecShort<1>, &MultAddMatVecShort<2>, &MultAddMatVecShort<3>,
      &MultAddMatVecShort<4>, &MultAddMatVecShort<5>, &MultAddMatVecShort<6>, &MultAddMatVecShort<7>,
      &MultAddMatVecShort<8>, &MultAddMatVecShort<9>, &MultAddMatVecShort<10>, &MultAddMatVecShort<11>,
      &MultAddMatVecShort<12>, &MultAddMatVecShort<13>, &MultAddMatVecShort<14>, &MultAddMatVecShort<15>,
      &MultAddMatVecShort<16>, &MultAddMatVecShort<17>, &MultAddMatVecShort<18>, &MultAddMatVecShort<19>,
      &MultAddMatVecShort<20>, &MultAddMatVecShort<21>, &MultAddMatVecShort<22>, &MultAddMatVecShort<23>,
      &MultAddMatVecShort<24>
    };

  

  // ************************** transpose Mat * vec ***************

  
  template <int SX>
  void MultMatTransVecShort (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    MatKernelDaxpy<1, SX, SET> (y.Size(), &x(0), 1, &a(0), a.Dist(), &y(0), 1);
  }
  


  NGS_DLL_HEADER void MultMatTransVec_intern (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    constexpr int SW = SIMD<double>::Size();
    size_t h = x.Size();
    size_t w = y.Size();
    size_t dist = a.Dist();

    size_t i = 0;
    for ( ; i+SW <= w; i+= SW)
      {
        SIMD<double> s0(0), s1(0), s2(0), s3(0);
        size_t j = 0;
        double * pa = &a(0,i);
        for ( ; j+4 <= h; j += 4, pa += 4*dist)
          {
            s0 += SIMD<double>(x(j)) * SIMD<double>(pa);
            s1 += SIMD<double>(x(j+1)) * SIMD<double>(pa+dist);
            s2 += SIMD<double>(x(j+2)) * SIMD<double>(pa+2*dist);
            s3 += SIMD<double>(x(j+3)) * SIMD<double>(pa+3*dist);
          }
        for ( ; j+2 <= h; j += 2, pa += 2*dist)
          {
            s0 += SIMD<double>(x(j)) * SIMD<double>(pa);
            s1 += SIMD<double>(x(j+1)) * SIMD<double>(pa+dist);
          }
        for ( ; j+1 <= h; j += 1, pa += dist)
          s2 += SIMD<double>(x(j)) * SIMD<double>(pa);
        SIMD<double> sum = (s0+s1)+(s2+s3);
        sum.Store(&y(i));
      }
    
    if (i < w)
      {
        SIMD<mask64> mask(w % SW);
        SIMD<double> s0(0), s1(0), s2(0), s3(0);
        size_t j = 0;
        double * pa = &a(0,i);
        for ( ; j+4 <= h; j += 4, pa += 4*dist)
          {
            s0 += SIMD<double>(x(j)) * SIMD<double>(pa, mask);
            s1 += SIMD<double>(x(j+1)) * SIMD<double>(pa+dist, mask);
            s2 += SIMD<double>(x(j+2)) * SIMD<double>(pa+2*dist, mask);
            s3 += SIMD<double>(x(j+3)) * SIMD<double>(pa+3*dist, mask);
          }
        for ( ; j+2 <= h; j += 2, pa += 2*dist)
          {
            s0 += SIMD<double>(x(j)) * SIMD<double>(pa, mask);
            s1 += SIMD<double>(x(j+1)) * SIMD<double>(pa+dist, mask);
          }
        for ( ; j+1 <= h; j += 1, pa += dist)
          s2 += SIMD<double>(x(j)) * SIMD<double>(pa, mask);
        SIMD<double> sum = (s0+s1)+(s2+s3);
        sum.Store(&y(i), mask);
      }
  }


  
  // typedef void REGCALL (*pmult_mattransvec)(BareSliceMatrix<>, FlatVector<>, FlatVector<>);  
  pmult_mattransvec dispatch_mattransvec[13] =
    {
      &MultMatTransVecShort<0>,
      &MultMatTransVecShort<1>,
      &MultMatTransVecShort<2>,
      &MultMatTransVecShort<3>,
      &MultMatTransVecShort<4>,
      &MultMatTransVecShort<5>,
      &MultMatTransVecShort<6>,
      &MultMatTransVecShort<7>,
      &MultMatTransVecShort<8>,
      &MultMatTransVecShort<9>,
      &MultMatTransVecShort<10>,
      &MultMatTransVecShort<11>,
      &MultMatTransVecShort<12>
    };
  








  template <int SX>
  void MultAddMatTransVecShort (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    double hx[max(1,SX)];
    for (size_t i = 0; i < SX; i++)
      hx[i] = s*x(i);
    MatKernelDaxpy<1, SX, ADD> (y.Size(), &hx[0], 1, &a(0), a.Dist(), &y(0), 1);
  }
  


  NGS_DLL_HEADER void MultAddMatTransVec_intern (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    constexpr int SW = SIMD<double>::Size();
    size_t h = x.Size();
    size_t w = y.Size();
    size_t dist = a.Dist();

    size_t i = 0;
    for ( ; i+SW <= w; i+= SW)
      {
        SIMD<double> s0(0), s1(0), s2(0), s3(0);
        size_t j = 0;
        double * pa = &a(0,i);
        for ( ; j+4 <= h; j += 4, pa += 4*dist)
          {
            s0 += SIMD<double>(x(j)) * SIMD<double>(pa);
            s1 += SIMD<double>(x(j+1)) * SIMD<double>(pa+dist);
            s2 += SIMD<double>(x(j+2)) * SIMD<double>(pa+2*dist);
            s3 += SIMD<double>(x(j+3)) * SIMD<double>(pa+3*dist);
          }
        for ( ; j+2 <= h; j += 2, pa += 2*dist)
          {
            s0 += SIMD<double>(x(j)) * SIMD<double>(pa);
            s1 += SIMD<double>(x(j+1)) * SIMD<double>(pa+dist);
          }
        for ( ; j+1 <= h; j += 1, pa += dist)
          s2 += SIMD<double>(x(j)) * SIMD<double>(pa);
        SIMD<double> sum = (s0+s1)+(s2+s3);
        sum *= s;
        sum += SIMD<double>(&y(i));
        sum.Store(&y(i));
      }
    
    if (i < w)
      {
        SIMD<mask64> mask(w % SW);
        SIMD<double> s0(0), s1(0), s2(0), s3(0);
        size_t j = 0;
        double * pa = &a(0,i);
        for ( ; j+4 <= h; j += 4, pa += 4*dist)
          {
            s0 += SIMD<double>(x(j)) * SIMD<double>(pa, mask);
            s1 += SIMD<double>(x(j+1)) * SIMD<double>(pa+dist, mask);
            s2 += SIMD<double>(x(j+2)) * SIMD<double>(pa+2*dist, mask);
            s3 += SIMD<double>(x(j+3)) * SIMD<double>(pa+3*dist, mask);
          }
        for ( ; j+2 <= h; j += 2, pa += 2*dist)
          {
            s0 += SIMD<double>(x(j)) * SIMD<double>(pa, mask);
            s1 += SIMD<double>(x(j+1)) * SIMD<double>(pa+dist, mask);
          }
        for ( ; j+1 <= h; j += 1, pa += dist)
          s2 += SIMD<double>(x(j)) * SIMD<double>(pa, mask);
        SIMD<double> sum = (s0+s1)+(s2+s3);
        sum *= s;
        sum += SIMD<double>(&y(i), mask);
        sum.Store(&y(i), mask);
      }
  }


  
  // typedef void REGCALL (*pmult_mattransvec)(BareSliceMatrix<>, FlatVector<>, FlatVector<>);  
  pmultadd_mattransvec dispatch_addmattransvec[13] =
    {
      &MultAddMatTransVecShort<0>,
      &MultAddMatTransVecShort<1>,
      &MultAddMatTransVecShort<2>,
      &MultAddMatTransVecShort<3>,
      &MultAddMatTransVecShort<4>,
      &MultAddMatTransVecShort<5>,
      &MultAddMatTransVecShort<6>,
      &MultAddMatTransVecShort<7>,
      &MultAddMatTransVecShort<8>,
      &MultAddMatTransVecShort<9>,
      &MultAddMatTransVecShort<10>,
      &MultAddMatTransVecShort<11>,
      &MultAddMatTransVecShort<12>
    };
  









  

  // ************************** Mult Add transpose Mat * vec, indirect ***************


  template <int SX>
  void MultAddMatTransVecShortI (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y,
                                 FlatArray<int> ind)
  {
    KernelAddMatTransVecI<SX> (s, ind.Size(), &a(0,0), a.Dist(), &x(0), &y(0), &ind[0]);
  }

  NGS_DLL_HEADER void MultAddMatTransVecIndirect_intern (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y,
                                                         FlatArray<int> ind)
  {
    constexpr int BS = 24;
    size_t i = 0;
    for (i = 0; i+BS <= y.Size(); i += BS)
      MultAddMatTransVecIndirect (s, a.Cols(i, i+BS), x, y.Range(i,i+BS), ind);
    MultAddMatTransVecIndirect (s, a.Cols(i, y.Size()), x, y.Range(i,y.Size()), ind);
    
    /*
    constexpr int SW = SIMD<double>::Size();
    size_t h = ind.Size();
    size_t w = y.Size();
    size_t dist = a.Dist();

    size_t i = 0;
    for ( ; i+SW <= w; i+= SW)
      {
        SIMD<double> s0(0), s1(0), s2(0), s3(0);
        size_t j = 0;
        double * pa = &a(0,i);
        for ( ; j+4 <= h; j += 4, pa += 4*dist)
          {
            s0 += SIMD<double>(x(ind[j])) * SIMD<double>(pa);
            s1 += SIMD<double>(x(ind[j+1])) * SIMD<double>(pa+dist);
            s2 += SIMD<double>(x(ind[j+2])) * SIMD<double>(pa+2*dist);
            s3 += SIMD<double>(x(ind[j+3])) * SIMD<double>(pa+3*dist);
          }
        for ( ; j+2 <= h; j += 2, pa += 2*dist)
          {
            s0 += SIMD<double>(x(ind[j])) * SIMD<double>(pa);
            s1 += SIMD<double>(x(ind[j+1])) * SIMD<double>(pa+dist);
          }
        for ( ; j+1 <= h; j += 1, pa += dist)
          s2 += SIMD<double>(x(ind[j])) * SIMD<double>(pa);
        SIMD<double> sum = (s0+s1)+(s2+s3);
        sum = SIMD<double>(&y(i)) + SIMD<double>(s) * sum;
        sum.Store(&y(i));
      }
    
    if (i < w)
      {
        SIMD<mask64> mask(w % SW);
        SIMD<double> s0(0), s1(0), s2(0), s3(0);
        size_t j = 0;
        double * pa = &a(0,i);
        for ( ; j+4 <= h; j += 4, pa += 4*dist)
          {
            s0 += SIMD<double>(x(ind[j])) * SIMD<double>(pa, mask);
            s1 += SIMD<double>(x(ind[j+1])) * SIMD<double>(pa+dist, mask);
            s2 += SIMD<double>(x(ind[j+2])) * SIMD<double>(pa+2*dist, mask);
            s3 += SIMD<double>(x(ind[j+3])) * SIMD<double>(pa+3*dist, mask);
          }
        for ( ; j+2 <= h; j += 2, pa += 2*dist)
          {
            s0 += SIMD<double>(x(ind[j])) * SIMD<double>(pa, mask);
            s1 += SIMD<double>(x(ind[j+1])) * SIMD<double>(pa+dist, mask);
          }
        for ( ; j+1 <= h; j += 1, pa += dist)
          s2 += SIMD<double>(x(ind[j])) * SIMD<double>(pa, mask);
        SIMD<double> sum = (s0+s1)+(s2+s3);
        sum = SIMD<double>(&y(i), mask) + SIMD<double>(s) * sum;        
        sum.Store(&y(i), mask);
      }
    */
  }

  pmultadd_mattransvecind dispatch_addmattransvecI[25] =
    {
      &MultAddMatTransVecShortI<0>,
      &MultAddMatTransVecShortI<1>,
      &MultAddMatTransVecShortI<2>,
      &MultAddMatTransVecShortI<3>,
      &MultAddMatTransVecShortI<4>,
      &MultAddMatTransVecShortI<5>,
      &MultAddMatTransVecShortI<6>,
      &MultAddMatTransVecShortI<7>,
      &MultAddMatTransVecShortI<8>,
      &MultAddMatTransVecShortI<9>,
      &MultAddMatTransVecShortI<10>,
      &MultAddMatTransVecShortI<11>,
      &MultAddMatTransVecShortI<12>,
      &MultAddMatTransVecShortI<13>,
      &MultAddMatTransVecShortI<14>,
      &MultAddMatTransVecShortI<15>,
      &MultAddMatTransVecShortI<16>,
      &MultAddMatTransVecShortI<17>,
      &MultAddMatTransVecShortI<18>,
      &MultAddMatTransVecShortI<19>,
      &MultAddMatTransVecShortI<20>,
      &MultAddMatTransVecShortI<21>,
      &MultAddMatTransVecShortI<22>,
      &MultAddMatTransVecShortI<23>,
      &MultAddMatTransVecShortI<24>
    };

  
  /* *********************** C = A * B ********************************* */
  
  // b.Width() = W * SIMD
  template <int W>
  INLINE void MatKernel2MultAB(size_t ha, size_t wa,
                               BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    size_t r = 0;
    size_t da = a.Dist();
    size_t dc = c.Dist();
    double * pa = &a(0,0);
    double * pc = &c(0,0);
    for ( ; r+4 <= ha; r += 4, pa += 4*da, pc += 4*dc)
      MatKernelMultAB<4,W> (wa, pa, da, &b(0,0), b.Dist(), pc, dc);
    switch (ha-r)
      {
      case 0: break;
      case 1:
        MatKernelMultAB<1,W> (wa, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      case 2:
        MatKernelMultAB<2,W> (wa, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      case 3:
        MatKernelMultAB<3,W> (wa, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      default:
        ;
      }
    return;
  }

  INLINE void MatKernel2MultABMask(SIMD<mask64> mask, size_t ha, size_t wa, BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    size_t r = 0;
    size_t da = a.Dist();
    size_t dc = c.Dist();
    double * pa = &a(0,0);
    double * pc = &c(0,0);
    for ( ; r+4 <= ha; r += 4, pa += 4*da, pc += 4*dc)
      MatKernelMultABMask<4,SET> (wa, mask, pa, da, &b(0,0), b.Dist(), pc, dc);
    switch (ha-r)
      {
      case 0: break;
      case 1:
        MatKernelMultABMask<1,SET> (wa, mask, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      case 2:
        MatKernelMultABMask<2,SET> (wa, mask, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      case 3:
        MatKernelMultABMask<3,SET> (wa, mask, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      default:
        ;
      }
    
  }

  /*
  // c = a * b
  void MultMatMat (size_t ha, size_t wa, size_t wb,
                   BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    size_t k = 0;
    constexpr size_t SW = SIMD<double>::Size();
    for ( ; k+3*SW <= wb; k += 3*SW)
      MatKernel2MultAB<3>(ha, wa, a, b.Cols(k,k+3*SW), c.Cols(k,k+3*SW));
    for ( ; k+SW <= wb; k += SW)
      MatKernel2MultAB<1>(ha, wa, a, b.Cols(k,k+SW), c.Cols(k,k+SW));

    if (k < wb)
      MatKernel2MultABMask(SIMD<mask64>(wb-k), ha, wa, a, b.Cols(k,k+SW), c.Cols(k,k+SW));
  }
  */

  template <size_t H, OPERATION OP, typename TB>
  INLINE void MatKernel2AddAB (size_t hb, size_t wb, double * pa, size_t da, TB * pb, size_t db, double * pc, size_t dc)
  {
    constexpr size_t SW = SIMD<double>::Size();
    constexpr size_t SWdTB = sizeof(SIMD<double>)/sizeof(TB);
    size_t l = 0, lb = 0;
    for ( ; l+3*SW <= wb; l += 3*SW, lb += 3*SWdTB)
      MatKernelMultAB<H,3,OP> (hb, pa, da, pb+lb, db, pc+l, dc);
    for ( ; l+SW <= wb; l += SW, lb += SWdTB)
      MatKernelMultAB<H,1,OP> (hb, pa, da, pb+lb, db, pc+l, dc);
    if (l < wb)
      MatKernelMultABMask<H,OP>(hb, SIMD<mask64>(wb-l), pa, da, pb+lb, db, pc+l, dc);
  }

  /*
  void MultMatMat_intern (size_t ha, size_t wa, size_t wb,
                          BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
#ifdef __AVX512F__
    constexpr size_t HA = 6;
#else
    constexpr size_t HA = 4;
#endif
    
    // blockwise B, fits into L2 cache
    constexpr size_t BBH = 128;
    constexpr size_t BBW = 96;
    constexpr size_t SW = SIMD<double>::Size();
    alignas(64) SIMD<double> bb[BBH*BBW/SW];

    double *pb = &b(0);
    for (size_t i = 0; i < wa; i += BBH, pb += BBH*b.Dist())
      for (size_t j = 0; j < wb; j += BBW)
        {
          size_t hbi = min2(BBH, wa-i);
          size_t wbi = min2(BBW, wb-j);
          CopyMatrixIn (hbi, wbi, pb+j, b.Dist(), &bb[0], BBW/SW);
          double * pa = &a(0)+i;
          double * pc = &c(0)+j;

          if (i == 0)
            {
              size_t k = 0;
              for ( ; k+HA <= ha; k += HA, pa += HA*a.Dist(), pc += HA * c.Dist())
                // MatKernel2AddAB<4,SET> (hbi, wbi, pa, a.Dist(),  (double*)&bb[0], BBW, pc, c.Dist());
                MatKernel2AddAB<HA,SET> (hbi, wbi, pa, a.Dist(),  &bb[0], BBW/SW, pc, c.Dist());
              switch (ha-k)
                {
                case 0: break;
                case 1: MatKernel2AddAB<1,SET> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                case 2: MatKernel2AddAB<2,SET> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                case 3: MatKernel2AddAB<3,SET> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                case 4:
                  if (HA > 4)
                    MatKernel2AddAB<4,SET> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist());
                  break;
                case 5:
                  if (HA > 5)
                    MatKernel2AddAB<5,SET> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist());
                  break;
                default: ; 
                }
              
            }
          else
            {
              size_t k = 0;
              for ( ; k+HA <= ha; k += HA, pa += HA*a.Dist(), pc += HA * c.Dist())
                MatKernel2AddAB<HA,ADD> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist());
              switch (ha-k)
                {
                case 0: break;
                case 1: MatKernel2AddAB<1,ADD> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                case 2: MatKernel2AddAB<2,ADD> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                case 3: MatKernel2AddAB<3,ADD> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                case 4:
                  if (HA > 4)
                    MatKernel2AddAB<4,ADD> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist());
                  break;
                case 5:
                  if (HA > 5)
                    MatKernel2AddAB<5,ADD> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist());
                  break;
                default: ;
                }
            }
        }
  }
  */



  template <size_t BBH, OPERATION OP>
  void REGCALL MultMatMat_intern2_SlimB (size_t ha, size_t wa, size_t wb,
                                         BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    double * pa0 = &a(0);
    size_t dista = a.Dist();
    double * pb = &b(0);
    size_t distb = b.Dist();

    
    constexpr size_t SW = SIMD<double>::Size();
    alignas(64) SIMD<double> bb[BBH];
#ifdef __AVX512F__
    constexpr size_t HA = 6;
#else
    constexpr size_t HA = 4;
#endif

    double * pc = &c(0,0);
    for (size_t j = 0; j+SW <= wb; j+=SW, pb += SW, pc += SW)
      {
        for (size_t k = 0; k < wa; k++)
          bb[k] = SIMD<double> (pb+k*distb);

        double * pc1 = pc;
        double * pa1 = pa0;
        size_t k = 0;
        for ( ; k+2*HA <= ha; k += 2*HA, pc1 += 2*HA*c.Dist(), pa1 += 2*HA*dista)
          MatKernelMultAB<2*HA, 1, OP> (wa, pa1, dista, bb, 1, pc1, c.Dist());
        for ( ; k+HA <= ha; k += HA, pc1 += HA*c.Dist(), pa1 += HA*dista)
          MatKernelMultAB<HA, 1, OP> (wa, pa1, dista, bb, 1, pc1, c.Dist());
        for ( ; k+1 <= ha; k += 1, pc1 += c.Dist(), pa1 += dista)
          MatKernelMultAB<1, 1, OP> (wa, pa1, dista, bb, 1, pc1, c.Dist());
      }

    
    if (wb % SW != 0)
      {
        SIMD<mask64> mask(wb%SW);
        for (size_t k = 0; k < wa; k++)
          bb[k] = SIMD<double> (pb+k*distb, mask);
        
        size_t k = 0;
        double * pc1 = pc;        
        for ( ; k+HA <= ha; k += HA, pc1 += HA*c.Dist())
          MatKernelMultABMask<HA, OP> (wa, mask, pa0+k*dista, dista, bb, 1, pc1, c.Dist());
        for ( ; k+1 <= ha; k += 1, pc1 += c.Dist())
          MatKernelMultABMask<1, OP> (wa, mask, pa0+k*dista, dista, bb, 1, pc1, c.Dist());
      }
  } 


  template <size_t BBH, OPERATION OP>
  void  REGCALL MultMatMat_intern2 (size_t ha, size_t wa, size_t wb,
                                     BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    if (wb < 3*SIMD<double>::Size())
      {
        MultMatMat_intern2_SlimB<BBH,OP> (ha, wa, wb, a, b, c);
        return;
      }

    double * pa0 = &a(0);
    size_t dista = a.Dist();
    double * pb = &b(0);
    size_t distb = b.Dist();
    
#ifdef __AVX512F__
    constexpr size_t HA = 6;
#else
    constexpr size_t HA = 4;
#endif
    
    // blockwise B, fits into L2 cache
    // constexpr size_t BBH = 128;
    constexpr size_t BBW = 96;
    constexpr size_t SW = SIMD<double>::Size();
    alignas(64) SIMD<double> bb[BBH*BBW/SW];

    for (size_t j = 0; j < wb; j += BBW)
      {
        size_t hbi = wa;
        size_t wbi = min2(BBW, wb-j);
        CopyMatrixIn (hbi, wbi, pb+j, distb, &bb[0], BBW/SW);

        double * pa = pa0;
        double * pc = &c(0)+j;
        
        size_t k = 0;
        for ( ; k+HA <= ha; k += HA, pa += HA*dista, pc += HA * c.Dist())
          MatKernel2AddAB<HA,OP> (hbi, wbi, pa, dista,  &bb[0], BBW/SW, pc, c.Dist());
        switch (ha-k)
          {
          case 0: break;
          case 1: MatKernel2AddAB<1,OP> (hbi, wbi, pa, dista, &bb[0], BBW/SW, pc, c.Dist()); break;
          case 2: MatKernel2AddAB<2,OP> (hbi, wbi, pa, dista, &bb[0], BBW/SW, pc, c.Dist()); break;
          case 3: MatKernel2AddAB<3,OP> (hbi, wbi, pa, dista, &bb[0], BBW/SW, pc, c.Dist()); break;
          case 4:
            if (HA > 4)
              MatKernel2AddAB<4,OP> (hbi, wbi, pa, dista, &bb[0], BBW/SW, pc, c.Dist());
            break;
          case 5:
            if (HA > 5)
              MatKernel2AddAB<5,OP> (hbi, wbi, pa, dista, &bb[0], BBW/SW, pc, c.Dist());
            break;
          default: ; 
          }
      }
  }
 

  template <size_t WA, OPERATION OP=SET> 
  REGCALL void MultMatMat_intern2_ShortSum (size_t ha, size_t wb,
                                           BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    if (WA <= 6) 
      MatKernelShortSum2<WA,OP> (ha, wb, a.Data(), a.Dist(), b.Data(), b.Dist(), c.Data(), c.Dist());
    else
      MatKernelShortSum<WA,OP> (ha, wb, a.Data(), a.Dist(), b.Data(), b.Dist(), c.Data(), c.Dist());
  }

  template <size_t WA, OPERATION OP=SET> 
  REGCALL void MultMatMat_intern2_ShortSumW (size_t ha, size_t wa, size_t wb,
                                             BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    if (WA <= 6) 
      MatKernelShortSum2<WA,OP> (ha, wb, a.Data(), a.Dist(), b.Data(), b.Dist(), c.Data(), c.Dist());
    else
      MatKernelShortSum<WA,OP> (ha, wb, a.Data(), a.Dist(), b.Data(), b.Dist(), c.Data(), c.Dist());
  }




  
  void MultMatMat_intern (size_t ha, size_t wa, size_t wb,
                          BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    constexpr size_t BBH = 128;
    if (wa <= BBH)
      {
        if (wb < 3*SIMD<double>::Size())
          MultMatMat_intern2_SlimB<BBH,SET> (ha, wa, wb, a, b, c);
        else
          MultMatMat_intern2<BBH,SET> (ha, wa, wb, a, b, c);
      }
    else
      {
        MultMatMat_intern2<BBH,SET> (ha, BBH, wb, a, b, c);    

        for (size_t i = BBH; i < wa; i += BBH)
          {
            a.IncPtr(BBH);
            b.IncPtr(BBH*b.Dist());
            size_t hbi = min2(BBH, wa-i);        
            MultMatMat_intern2<BBH,ADD> (ha, hbi, wb, a, b, c);
          }
      }
  }

  void REGCALL MinusMultAB_intern (size_t ha, size_t wa, size_t wb,
                                   BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    constexpr size_t BBH = 128;
    if (wa <= BBH)
      {
        if (wb < 3*SIMD<double>::Size())
          MultMatMat_intern2_SlimB<BBH,SETNEG> (ha, wa, wb, a, b, c);
        else
          MultMatMat_intern2<BBH,SETNEG> (ha, wa, wb, a, b, c);
      }
    else
      {
        MultMatMat_intern2<BBH,SETNEG> (ha, BBH, wb, a, b, c);    

        for (size_t i = BBH; i < wa; i += BBH)
          {
            a.IncPtr(BBH);
            b.IncPtr(BBH*b.Dist());
            size_t hbi = min2(BBH, wa-i);        
            MultMatMat_intern2<BBH,SUB> (ha, hbi, wb, a, b, c);
          }
      }
  }

  
  /*
  void MinusMultAB_intern (size_t ha, size_t wa, size_t wb,
                           BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    constexpr size_t BBH = 128;
    if (wb < 3*SIMD<double>::Size())
      MultMatMat_intern2_SlimB<BBH,SETNEG> (ha, wa, wb, a, b, c);
    else
      {
        for (size_t i = 0; i < wa; i += BBH, a.IncPtr(BBH), b.IncPtr(BBH*b.Dist()))
          {
            size_t hbi = min2(BBH, wa-i);
            if (i == 0)
              MultMatMat_intern2<BBH,SETNEG> (ha, hbi, wb, a, b, c);
            else
              MultMatMat_intern2<BBH,SUB> (ha, hbi, wb, a, b, c);              
          }
      }
  }
  */
  
  void REGCALL AddAB_intern (size_t ha, size_t wa, size_t wb,
                             BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    constexpr size_t BBH = 128;
    if (wa <= BBH && wb < 3*SIMD<double>::Size())
      MultMatMat_intern2_SlimB<BBH,ADD> (ha, wa, wb, a, b, c);
    else
      for (size_t i = 0; i < wa; i += BBH, a.IncPtr(BBH), b.IncPtr(BBH*b.Dist()))
        {
          size_t hbi = min2(BBH, wa-i);        
          MultMatMat_intern2<BBH,ADD> (ha, hbi, wb, a, b, c);
        }
  }

  void REGCALL SubAB_intern (size_t ha, size_t wa, size_t wb,
                             BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    constexpr size_t BBH = 128;
    if (wa <= BBH && wb < 3*SIMD<double>::Size())
      MultMatMat_intern2_SlimB<BBH,SUB> (ha, wa, wb, a, b, c);
    else
      for (size_t i = 0; i < wa; i += BBH, a.IncPtr(BBH), b.IncPtr(BBH*b.Dist()))
        {
          size_t hbi = min2(BBH, wa-i);        
          MultMatMat_intern2<BBH,SUB> (ha, hbi, wb, a, b, c);
        }
  }


  pmultABW dispatch_multAB[];
  auto init_multAB = [] ()
  {
    Iterate<std::size(dispatch_multAB)-1> ([&] (auto i)
    { dispatch_multAB[i] = &MultMatMat_intern2_ShortSumW<i,SET>; });
    // Iterate<std::size(dispatch_multAB)-1> ([&] (auto i)
    // { dispatch_multAB[i] = &MultMatMat_intern; });
    dispatch_multAB[std::size(dispatch_multAB)-1] = &MultMatMat_intern;
    return 1;
  }();

  pmultABW dispatch_minusmultAB[];
  auto init_minusmultAB = [] ()
  {
    Iterate<std::size(dispatch_minusmultAB)-1> ([&] (auto i)
    { dispatch_minusmultAB[i] = &MultMatMat_intern2_ShortSumW<i,SETNEG>; });
    dispatch_minusmultAB[std::size(dispatch_minusmultAB)-1] = &MinusMultAB_intern;
    return 1;
  }();

  pmultABW dispatch_addAB[];
  auto init_addAB = [] ()
  {
    Iterate<std::size(dispatch_addAB)-1> ([&] (auto i)
    { dispatch_addAB[i] = &MultMatMat_intern2_ShortSumW<i,ADD>; });
    dispatch_addAB[std::size(dispatch_addAB)-1] = &AddAB_intern;
    return 1;
  }();

  pmultABW dispatch_subAB[];
  auto init_subAB = [] ()
  {
    Iterate<std::size(dispatch_subAB)-1> ([&] (auto i)
    { dispatch_subAB[i] = &MultMatMat_intern2_ShortSumW<i,SUB>; });
    dispatch_subAB[std::size(dispatch_subAB)-1] = &SubAB_intern;
    return 1;
  }();


  


  /*

    // was not fast ...

  template <size_t HA,size_t WAREST>
  INLINE void MultMatMat_SmallA_intern2 (size_t wa, size_t wb,
                                  double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
  {
    MatKernelDaxpy<HA,WAREST,SET> (wb, pa, da, pb, db, pc, dc);
    for (size_t j = WAREST ; j+4 <= wa; j += 4)
      MatKernelDaxpy<HA,4,ADD> (wb, pa+j, da, pb+j*db, db, pc, dc);
  }
  
  template <size_t WAREST>
  void REGCALL MultMatMat_SmallA_intern3 (size_t ha, size_t wa, size_t wb,
                                            BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    double * pa = &a(0);
    size_t da = a.Dist();
    double * pb = &b(0);
    size_t db = b.Dist();
    double * pc = &c(0);
    size_t dc = c.Dist();

    size_t i = 0;
    for ( ; i+3 <= ha; i += 3, pa += 3*da, pc += 3*dc)
      MultMatMat_SmallA_intern2<3,WAREST> (wa, wb, pa, da, pb, db, pc, dc);
    switch (ha-i)
      {
      case 1:
        MultMatMat_SmallA_intern2<1,WAREST> (wa, wb, pa, da, pb, db, pc, dc);
        break;
      case 2:
        MultMatMat_SmallA_intern2<2,WAREST> (wa, wb, pa, da, pb, db, pc, dc);
        break;
      }
  }

  void REGCALL MultMatMat_SmallA_intern (size_t ha, size_t wa, size_t wb,
                                           BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    switch (wa & 3)
      {
      case 1: MultMatMat_SmallA_intern3<1> (ha, wa, wb, a, b, c); break;
      case 2: MultMatMat_SmallA_intern3<2> (ha, wa, wb, a, b, c); break;
      case 3: MultMatMat_SmallA_intern3<3> (ha, wa, wb, a, b, c); break;
      case 0: MultMatMat_SmallA_intern3<4> (ha, wa, wb, a, b, c); break;
      default:
        __assume(false); 
      }
        
  }
  */


  
  /* ********************* C = A * B  with B is SIMD **************************** */

  // b.Width() = W * SIMD
  template <int W>
  INLINE void MatKernel2MultAB(size_t ha, size_t wa,
                               BareSliceMatrix<> a,
                               BareSliceMatrix<SIMD<double>> b, BareSliceMatrix<SIMD<double>> c)
  {
    size_t r = 0;
    size_t da = a.Dist();
    size_t dc = c.Dist();
    double * pa = &a(0,0);
    SIMD<double> * pc = &c(0,0);
    for ( ; r+4 <= ha; r += 4, pa += 4*da, pc += 4*dc)
      MatKernelAlignedMultAB<4,W> (wa, pa, da, &b(0,0), b.Dist(), pc, dc);
    switch (ha-r)
      {
      case 0: break;
      case 1:
        MatKernelAlignedMultAB<1,W> (wa, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      case 2:
        MatKernelAlignedMultAB<2,W> (wa, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      case 3:
        MatKernelAlignedMultAB<3,W> (wa, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      default:
        ;
      }
    return;
  }

  // c = a * b
  void MultMatMat_intern (size_t ha, size_t wa, size_t wb,
                          BareSliceMatrix<> a, BareSliceMatrix<SIMD<double>> b, BareSliceMatrix<SIMD<double>> c)
  {
    size_t k = 0;
    constexpr size_t SW = SIMD<double>::Size();
    for ( ; k+3 <= wb; k += 3)
      MatKernel2MultAB<3>(ha, wa, a, b.Cols(k,k+3), c.Cols(k,k+3));
    for ( ; k+SW <= wb; k += SW)
      MatKernel2MultAB<1>(ha, wa, a, b.Cols(k,k+SW), c.Cols(k,k+SW));
  }

  /* ******************************* A^T B *************************************** */

  template <size_t WA>
  INLINE void MultAtBSmallWA2 (size_t ha, size_t wb, BareSliceMatrix<double> a, BareSliceMatrix<double> b,
                               BareSliceMatrix<double> c)
  {
    constexpr size_t SW = SIMD<double>::Size();    

    size_t da = a.Dist();
    size_t db = b.Dist();
    size_t dc = c.Dist();

    size_t j = 0;
    double * pc0 = &c(0);
    for ( ; j+2*SW <= wb; j+=2*SW, pc0+=2*SW)
      {
        Vec<WA, SIMD<double>> sum0;
        Vec<WA, SIMD<double>> sum1;
        for (size_t i = 0; i < WA; i++)
          {
            sum0[i] = SIMD<double> (0);
            sum1[i] = SIMD<double> (0);
          }
        
        double * pa = &a(0);
        double * pb = &b(j);
        __assume(ha > 0);
        for (size_t k = 0; k < ha; k++, pa += da, pb += db)
          {
            SIMD<double> bjk0(pb);
            SIMD<double> bjk1(pb+SW);
            
            for (size_t i = 0; i < WA; i++)
              {
                SIMD<double> ai(pa[i]);
                FMAasm (bjk0, ai, sum0[i]);
                FMAasm (bjk1, ai, sum1[i]);
              }
          }

        double * pc = pc0;
        for (size_t i = 0; i < WA; i++, pc += dc)
          {
            sum0[i].Store (pc);
            sum1[i].Store (pc+SW);
          }
      }

    for ( ; j+SW <= wb; j+=SW, pc0+=SW)
      {
        Vec<WA, SIMD<double>> sum;
        for (size_t i = 0; i < WA; i++)
          sum[i] = SIMD<double> (0);
        
        double * pa = &a(0);
        double * pb = &b(j);
        __assume(ha > 0);
        for (size_t k = 0; k < ha; k++, pa += da, pb += db)
          {
            SIMD<double> bjk(pb);
            for (size_t i = 0; i < WA; i++)
              // sum[i] += bjk*pa[i];
              FMAasm (bjk, SIMD<double>(pa[i]), sum[i]);
          }

        double * pc = pc0;
        for (size_t i = 0; i < WA; i++, pc += dc)
          sum[i].Store (pc);
      }

    if(j==wb)
        return;
    
    SIMD<mask64> mask(wb-j);
    std::array<SIMD<double>,WA> sum;
    for (size_t i = 0; i < WA; i++)
      sum[i] = SIMD<double> (0);
    
    double * pa = &a(0);
    double * pb = &b(j);
    __assume(ha > 0);    
    for (size_t k = 0; k < ha; k++, pa += da, pb += db)
      {
        SIMD<double> bi(pb, mask);
        for (size_t i = 0; i < WA; i++)
          sum[i] += bi*pa[i];
      }

    double * pc = &c(j);    
    for (size_t i = 0; i < WA; i++, pc += dc)
      sum[i].Store (pc, mask);
  }



  template <size_t WA, OPERATION OP>
  void REGCALL MultAtBSmallWA (size_t ha, size_t /* wa */, size_t wb, BareSliceMatrix<double> a, BareSliceMatrix<double> b,
                               BareSliceMatrix<double> c)

  {
    if (WA <= 6 && OP == SET)
      {
        MultAtBSmallWA2<WA> (ha, wb, a, b, c);
        return;
      }
    MatKernelAtB_SmallWA<WA,OP> (ha, wb, a.Data(), a.Dist(), b.Data(), b.Dist(), c.Data(), c.Dist());
  }

  // template <> pmultABW dispatch_atb<false,true>[];

  template <OPERATION OP, size_t MAXHA>
  void REGCALL MultAtB_intern2 (size_t ha, size_t wa, size_t wb,
                                BareSliceMatrix<double> a, BareSliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    // c.AddSize(a.Width(), b.Width()) = 1.0 * Trans(a) * b;  // avoid recursion
    // return;

    constexpr size_t bs = SIMD<double>::Size();   

    alignas(64) SIMD<double> mem[MAXHA];
    
    size_t i = 0;
    
    for ( ; i+bs <= wa; i += bs, a.IncPtr(bs), c.IncPtr(bs*c.Dist()))
      {
        CopyMatrixIn (ha, bs, a.Data(), a.Dist(), &mem[0], 1);
        MatKernelAtB_SmallWA2<bs,OP> (ha, wb,
                                      (double*)&mem[0], bs,
                                      &b(0), b.Dist(),
                                      &c(0), c.Dist());
      }

    if (i == wa) return;
    dispatch_atb<OP==ADD || OP==SUB, OP==SET || OP==ADD>::ptrs[wa-i] (ha, wa-i, wb, a, b, c);    
  }

  
  template <OPERATION OP>
  void REGCALL MultAtB_intern (size_t ha, size_t wa, size_t wb, BareSliceMatrix<double> a, BareSliceMatrix<double> b,
                               BareSliceMatrix<double> c)
  {
    constexpr size_t BBW = 120;
    constexpr size_t ABH = 120;

    alignas(64) SIMD<double> mem[ABH*BBW/SIMD<double>::Size()];
    FlatMatrix<> bbm(ABH, BBW, (double*)&mem);
    
    for (size_t i = 0; i < ha; i += ABH)
      {
        IntRange ri(i, i+min(ha-i, ABH));
        auto bi = b.Rows(ri);
        auto ai = a.Rows(ri);

        for (size_t j = 0; j < wb; j += BBW)
          {
            IntRange rj(j, j+min(wb-j, BBW));

            auto bij = bi.Cols(rj);
            CopyMatrixIn (ri.Size(), rj.Size(), bij.Data(), bij.Dist(), &mem[0], BBW/SIMD<double>::Size());

            if (i > 0 || (OP == ADD || OP == SUB))
              MultAtB_intern2<AddOp(OP),ABH> (ri.Size(), wa, rj.Size(), ai, bbm, c.Cols(rj));
            else
              MultAtB_intern2<OP,ABH> (ri.Size(), wa, rj.Size(), ai, bbm, c.Cols(rj));
          }
      }
  }
  

  template <bool ADD, bool POS>
  pmultABW dispatch_atb<ADD,POS>::ptrs[];

  auto init_atb = [] ()
  {
    Iterate<std::size(dispatch_atb<false,true>::ptrs)-1> ([&] (auto i)
    {
      dispatch_atb<false,true>::ptrs[i] = &MultAtBSmallWA<i,SET>;
      dispatch_atb<true,true>::ptrs[i] = &MultAtBSmallWA<i,ADD>;
      dispatch_atb<false,false>::ptrs[i] = &MultAtBSmallWA<i,SETNEG>;
      dispatch_atb<true,false>::ptrs[i] = &MultAtBSmallWA<i,SUB>;
    });
    dispatch_atb<false,true>::ptrs[std::size(dispatch_atb<false,true>::ptrs)-1] = &MultAtB_intern<SET>;
    dispatch_atb<true,true>::ptrs[std::size(dispatch_atb<true,true>::ptrs)-1] = &MultAtB_intern<ADD>;
    dispatch_atb<false,false>::ptrs[std::size(dispatch_atb<false,false>::ptrs)-1] = &MultAtB_intern<SETNEG>;
    dispatch_atb<true,false>::ptrs[std::size(dispatch_atb<true,false>::ptrs)-1] = &MultAtB_intern<SUB>;
    return 1;
  }();
  
  /* ***************************** A * B^T *************************************** */

  template <int SX, OPERATION OP>
  void REGCALL MultABtSmallWA (size_t ah, size_t bh, BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    double * pa = a.Data();
    double * pc = c.Data();
    for (size_t i = 0; i < ah; i++, pa += a.Dist(), pc += c.Dist())
      KernelMatVec<SX,OP> (bh, b.Data(), b.Dist(), pa, pc);
  }

  /*
  pfunc_abt dispatch_abt[25] =
    { &MultABtSmallWA<0>, &MultABtSmallWA<1>, &MultABtSmallWA<2>, &MultABtSmallWA<3>,
      &MultABtSmallWA<4>, &MultABtSmallWA<5>, &MultABtSmallWA<6>, &MultABtSmallWA<7>,
      &MultABtSmallWA<8>, &MultABtSmallWA<9>, &MultABtSmallWA<10>, &MultABtSmallWA<11>,
      &MultABtSmallWA<12>, &MultABtSmallWA<13>, &MultABtSmallWA<14>, &MultABtSmallWA<15>,
      &MultABtSmallWA<16>, &MultABtSmallWA<17>, &MultABtSmallWA<18>, &MultABtSmallWA<19>,
      &MultABtSmallWA<20>, &MultABtSmallWA<21>, &MultABtSmallWA<22>, &MultABtSmallWA<23>,
      &MultABtSmallWA<24>
    };
  */

  pfunc_abt dispatch_abt[];
  auto init_abt = [] ()
  {
    Iterate<std::size(dispatch_abt)> ([&] (auto i)
    { dispatch_abt[i] = &MultABtSmallWA<i,SET>; });
    // dispatch_matvec[std::size(dispatch_matvec)-1] = &MultMatVec_intern;
    return 1;
  }();

  pfunc_abt dispatch_addabt[];
  auto init_addabt = [] ()
  {
    Iterate<std::size(dispatch_abt)> ([&] (auto i)
    { dispatch_addabt[i] = &MultABtSmallWA<i,ADD>; });
    // dispatch_matvec[std::size(dispatch_matvec)-1] = &MultMatVec_intern;
    return 1;
  }();
  
  
  template <typename TAB, typename FUNC>
  INLINE void TAddABt4 (size_t wa, size_t hc, size_t wc,
                        TAB * pa, size_t da, TAB * pb, size_t db, double * pc, size_t dc,
                        FUNC func)
  {
#ifdef __AVX512F__
    constexpr size_t HA = 6;
#else
    constexpr size_t HA = 3;
#endif
    
    TAB * pb0 = pb;
    size_t i = 0;
    for ( ; i+HA <= hc; i += HA, pa += HA*da, pc += HA*dc)
      {
        TAB * pb = pb0;
        size_t j = 0;
        for ( ; j+4 <= wc; j += 4, pb += 4*db)
          {
            auto scal = MatKernelScalAB<HA,4>(wa, pa, da, pb, db);
            Iterate<HA> ([&] (auto i) {
                double * pci = pc+i.value*dc+j;
                auto si = func (SIMD<double,4>(pci), get<i.value>(scal));
                si.Store(pci);
              });
          }
        for ( ; j+2 <= wc; j += 2, pb += 2*db)
          {
            auto scal = MatKernelScalAB<HA,2>(wa, pa, da, pb, db);
            Iterate<HA> ([&] (auto i) {
                double * pci = pc+i.value*dc+j;
                auto si = func (SIMD<double,2>(pci), get<i.value>(scal));
                si.Store(pci);
              });
          }
        for ( ; j < wc; j++, pb += db)
          {
            auto scal = MatKernelScalAB<HA,1>(wa, pa, da, pb, db);
            Iterate<HA> ([&] (auto i) {
                double * pci = pc+i.value*dc+j;
                auto si = func (*pci, get<i.value>(scal));
                *pci = si;
              });
          }
      }
    for ( ; i < hc; i ++, pa += da, pc += dc)
      {
        double * pc1 = pc;
        TAB * pb = pb0;
        size_t j = 0;
        for ( ; j+4 <= wc; j += 4, pb += 4*db)
          {
            auto scal = MatKernelScalAB<1,4>(wa, pa, da, pb, db);
            auto s1 = func (SIMD<double,4>(pc1+j), get<0>(scal));
            s1.Store(pc1+j);
          }
        for ( ; j < wc; j++, pb += db)
          {
            auto scal = MatKernelScalAB<1,1>(wa, pa, da, pb, db);
            auto s1 = func (pc1[j], get<0>(scal));
            pc1[j] = s1;
          }
      }
  }
  

  template <typename TAB, typename FUNC>
  void TAddABt2 (size_t wa, size_t ha, size_t hb,
                 TAB * pa, size_t da, TAB * pb, size_t db, double * pc, size_t dc,
                 FUNC func)
  {
    constexpr size_t bsa = 96; // height a
    constexpr size_t bsb = 32; // height b    
    for (size_t i = 0; i < ha; i += bsa, pa += bsa*da, pc += bsa*dc)
      {
        size_t hha = min2(bsa, ha-i);
        TAB * hpb = pb;
        for (size_t j = 0; j < hb; j += bsb, hpb += bsb*db)
          TAddABt4 (wa, hha, min2(bsb, hb-j),
                    pa, da, hpb, db, pc+j, dc, func);
        
      }
  }

  
  template <typename FUNC>
  void TAddABt1 (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c,
                FUNC func)
  {
    constexpr size_t bs = 256; // inner-product loop
    size_t wa = a.Width();
    double *pa = a.Data();
    double *pb = b.Data();
    double *pc = c.Data();
    for (size_t i = 0; i < wa; i += bs, pa+=bs, pb+=bs)
      TAddABt2 (min2(bs,wa-i), a.Height(), b.Height(),
                pa, a.Dist(), pb, b.Dist(), pc, c.Dist(), func);
  }

  void MultABt_intern (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    // c = a * Trans(b);

    constexpr size_t bs = 256;
    size_t wa = a.Width();

    TAddABt2 (min2(bs, wa), a.Height(), b.Height(),
              &a(0), a.Dist(), &b(0), b.Dist(), &c(0), c.Dist(),
              [] (auto c, auto ab) { return ab; });

    if (wa > bs)
      TAddABt1 (a.Cols(bs, wa), b.Cols(bs, wa), c, [] (auto c, auto ab) { return c+ab; });        
  }

  void MinusMultABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    // c = -a * Trans(b);
    
    constexpr size_t bs = 256;
    size_t wa = a.Width();

    TAddABt2 (min2(bs, wa), a.Height(), b.Height(),
              &a(0), a.Dist(), &b(0), b.Dist(), &c(0), c.Dist(),
              [] (auto c, auto ab) { return -ab; });

    if (wa > bs)
      TAddABt1 (a.Cols(bs, wa), b.Cols(bs, wa), c, [] (auto c, auto ab) { return c-ab; });        
  }

  
  void AddABt_intern (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    // c += a * Trans(b);
    TAddABt1 (a, b, c, [] (auto c, auto ab) { return c+ab; });
  }

  void SubABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    // c -= a * Trans(b);
    TAddABt1 (a, b, c, [] (auto c, auto ab) { return c-ab; });
  }






  /* ***************************** A * B^T, A,B SIMD *********************************** */

  
  template <typename FUNC>
  void TAddABt1 (SliceMatrix<SIMD<double>> a, SliceMatrix<SIMD<double>> b, BareSliceMatrix<double> c,
                 FUNC func)
  {
    constexpr size_t bs = 256; // inner-product loop
    size_t wa = a.Width();
    SIMD<double> *pa = a.Data();
    SIMD<double> *pb = b.Data();
    double *pc = c.Data();
    for (size_t i = 0; i < wa; i += bs, pa+=bs, pb+=bs)
      TAddABt2 (min2(bs,wa-i), a.Height(), b.Height(),
                pa, a.Dist(), pb, b.Dist(), pc, c.Dist(), func);
  }
  
  void AddABt (SliceMatrix<SIMD<double>> a, SliceMatrix<SIMD<double>> b, BareSliceMatrix<double> c)
  {
    // c += a * Trans(b);
    TAddABt1 (a, b, c, [] (auto c, auto ab) { return c+ab; });
  }

  void SubABt (SliceMatrix<SIMD<double>> a, SliceMatrix<SIMD<double>> b, BareSliceMatrix<double> c)
  {
    // c -= a * Trans(b);
    TAddABt1 (a, b, c, [] (auto c, auto ab) { return c-ab; });
  }



  /* *********************** AddABt-Sym ************************ */

  template <typename TAB, typename FUNC>
  INLINE void TAddABt4Sym (size_t wa, size_t hc, size_t wc,
                           TAB * pa, size_t da, TAB * pb, size_t db, double * pc, size_t dc,
                           FUNC func)
  {
#ifdef __AVX512F__
    constexpr size_t HA = 6;
#else
    constexpr size_t HA = 3;
#endif
    
    TAB * pb0 = pb;
    size_t i = 0;
    for ( ; i+HA <= hc; i += HA, pa += HA*da, pc += HA*dc)
      {
        TAB * pb = pb0;
        size_t j = 0;
        for ( ; j+4 <= i+HA; j += 4, pb += 4*db)
          {
            auto scal = MatKernelScalAB<HA,4>(wa, pa, da, pb, db);
            Iterate<HA> ([&] (auto i) {
                double * pci = pc+i.value*dc+j;
                auto si = func (SIMD<double,4>(pci), get<i.value>(scal));
                si.Store(pci);
              });
          }
        for ( ; j < i+HA; j++, pb += db)
          {
            auto scal = MatKernelScalAB<HA,1>(wa, pa, da, pb, db);
            Iterate<HA> ([&] (auto i) {
                double * pci = pc+i.value*dc+j;
                auto si = func (*pci, get<i.value>(scal));
                *pci = si;
              });
          }
      }
    for ( ; i < hc; i ++, pa += da, pc += dc)
      {
        double * pc1 = pc;
        TAB * pb = pb0;
        size_t j = 0;
        for ( ; j+3 <= i; j += 4, pb += 4*db)
          {
            auto scal = MatKernelScalAB<1,4>(wa, pa, da, pb, db);
            auto s1 = func (SIMD<double,4>(pc1+j), get<0>(scal));
            s1.Store(pc1+j);
          }
        for ( ; j <= i; j++, pb += db)
          {
            auto scal = MatKernelScalAB<1,1>(wa, pa, da, pb, db);
            auto s1 = func (pc1[j], get<0>(scal));
            pc1[j] = s1;
          }
      }
  }


  void AddABtSym (SliceMatrix<double> a,
                  SliceMatrix<double> b,
                  BareSliceMatrix<double> c)
  {
    TAddABt4Sym(a.Width(), a.Height(), b.Height(),
                a.Data(), a.Dist(), b.Data(), b.Dist(), c.Data(), c.Dist(),
                [] (auto c, auto ab) { return c+ab; });
  }


  void AddABtSym (SliceMatrix<SIMD<double>> a,
                  SliceMatrix<SIMD<double>> b,
                  BareSliceMatrix<double> c)
  {
    TAddABt4Sym(a.Width(), a.Height(), b.Height(),
                a.Data(), a.Width(), b.Data(), b.Width(), c.Data(), c.Dist(),
                [] (auto c, auto ab) { return c+ab; });
    /*
    AddABtSym (SliceMatrix<double> (AFlatMatrix<double>(a)),
               SliceMatrix<double> (AFlatMatrix<double>(b)), c);
    */
  }



  
  
  /* *************************** copied from symbolicintegrator, needs some rework ***** */



  

  void AddABt (FlatMatrix<SIMD<Complex>> a,
               FlatMatrix<SIMD<Complex>> b,
               SliceMatrix<Complex> c)
  {
    for (size_t i = 0; i < c.Height(); i++)
      for (size_t j = 0; j < c.Width(); j++)
        {
          SIMD<Complex> sum(0.0);
          for (size_t k = 0; k < a.Width(); k++)
            sum += a(i,k) * b(j,k);
          c(i,j) += HSum(sum);
        }
  }
  
  void AddABtSym (FlatMatrix<SIMD<Complex>> a,
                  FlatMatrix<SIMD<Complex>> b,
                  SliceMatrix<Complex> c)
  {
    AddABt (a, b, c);
  }
  /*
  void AddABt (FlatMatrix<SIMD<double>> a,
               FlatMatrix<SIMD<Complex>> b,
               SliceMatrix<Complex> c)
  {
    size_t i = 0;
    for ( ; i < c.Height()-1; i+=2)
      for (size_t j = 0; j < c.Width(); j++)
        {
          SIMD<Complex> sum1(0.0);
          SIMD<Complex> sum2(0.0);
          for (size_t k = 0; k < a.Width(); k++)
            {
              sum1 += a(i,k) * b(j,k);
              sum2 += a(i+1,k) * b(j,k);
            }
          c(i,j) += HSum(sum1);
          c(i+1,j) += HSum(sum2);
        }
    
    if (i < c.Height())
      for (size_t j = 0; j < c.Width(); j++)
        {
          SIMD<Complex> sum(0.0);
          for (size_t k = 0; k < a.Width(); k++)
            sum += a(i,k) * b(j,k);
          c(i,j) += HSum(sum);
        }
  }
  */


      
  void AddABt1 (SliceMatrix<SIMD<double>> a,
                SliceMatrix<SIMD<Complex>> b,
                SliceMatrix<Complex> c)
  {
    size_t i = 0;
    size_t wa = a.Width();
    size_t da = a.Dist();
    size_t db = b.Dist();
    if (wa == 0) return;
    
    for ( ; i+1 < c.Height(); i+=2)
      {
        auto pa1 = &a(i,0);
        auto pa2 = pa1 + da;
        auto pb1 = &b(0,0);
        size_t j = 0;
        for ( ; j+1 < c.Width(); j+=2, pb1 += 2*db)
          // for ( ; j+1 < c.Width(); j+=1, pb1 += db)
          {
            auto pb2 = pb1 + db;
            
            SIMD<Complex> sum11(0.0);
            SIMD<Complex> sum21(0.0);
            SIMD<Complex> sum12(0.0);
            SIMD<Complex> sum22(0.0);
            __assume (wa > 0);
            for (size_t k = 0; k < wa; k++)
              {
                sum11 += pa1[k] * pb1[k];
                sum21 += pa2[k] * pb1[k];
                sum12 += pa1[k] * pb2[k];
                sum22 += pa2[k] * pb2[k];
              }

            Complex s11, s21, s12, s22;
            std::tie(s11,s21) = HSum(sum11, sum21);
            std::tie(s12,s22) = HSum(sum12, sum22);
            c(i,j) += s11;
            c(i,j+1) += s12;
            c(i+1,j) += s21;
            c(i+1,j+1) += s22;
          }
        if (j < c.Width())
          {
            SIMD<Complex> sum1(0.0);
            SIMD<Complex> sum2(0.0);
            __assume (wa > 0);
            for (size_t k = 0; k < wa; k++)
              {
                sum1 += pa1[k] * pb1[k];
                sum2 += pa2[k] * pb1[k];
              }

            Complex s1, s2;
            std::tie(s1,s2) = HSum(sum1, sum2);
            c(i,j) += s1;
            c(i+1,j) += s2;
          }
      }
    
    if (i < c.Height())
      for (size_t j = 0; j < c.Width(); j++)
        {
          SIMD<Complex> sum(0.0);
          for (size_t k = 0; k < wa; k++)
            sum += a(i,k) * b(j,k);
          c(i,j) += HSum(sum);
        }
  }

  Timer timer_addabtdc ("AddABt-double-complex");
  Timer timer_addabtcd ("AddABt-complex-double");
  Timer timer_addabtdcsym ("AddABt-double-complex, sym");

  // block and pack B
  template <size_t K>
  void AddABt2 (SliceMatrix<SIMD<double>> a,
                SliceMatrix<SIMD<Complex>> b,
                SliceMatrix<Complex> c)
  {
    constexpr size_t bs = 32;
    SIMD<Complex> memb[bs*K];
    // M * K * sizeof(SIMD<Complex>) = 32 * 64 * 64 = 128 KB
    for (size_t k = 0; k < b.Height(); k+= bs)
      {
        size_t k2 = min2(k+bs, b.Height());
        FlatMatrix<SIMD<Complex>> tempb(k2-k, b.Width(), &memb[0]);
        tempb = b.Rows(k,k2);
        AddABt1 (a, tempb, c.Cols(k,k2));
      }
  }
  
  void AddABt (SliceMatrix<SIMD<double>> a,
               SliceMatrix<SIMD<Complex>> b,
               SliceMatrix<Complex> c)
  {
    ThreadRegionTimer reg(timer_addabtdc, TaskManager::GetThreadId());
    NgProfiler::AddThreadFlops(timer_addabtdc, TaskManager::GetThreadId(),
                               a.Height()*b.Height()*a.Width()*2*SIMD<double>::Size());
    constexpr size_t bs = 64;
    for (size_t k = 0; k < a.Width(); k+=bs)
      {
        size_t k2 = min2(k+bs, a.Width());
        AddABt2<bs> (a.Cols(k,k2), b.Cols(k,k2), c);
      }
  }



  void AddABt (SliceMatrix<SIMD<Complex>> a, SliceMatrix<SIMD<double>> b, SliceMatrix<Complex> c)
  {
    ThreadRegionTimer reg(timer_addabtcd, TaskManager::GetThreadId());
    NgProfiler::AddThreadFlops(timer_addabtcd, TaskManager::GetThreadId(),
                               a.Height()*b.Height()*a.Width()*2*SIMD<double>::Size());

    for (size_t i = 0; i < c.Height(); i++)
      for (size_t j = 0; j < c.Width(); j++)
        {
          SIMD<Complex> sum = 0.0;
          auto rowa = a.Row(i);
          auto rowb = b.Row(j);
          for (size_t k = 0; k < a.Width(); k++)
            sum += rowa(k)*rowb(k);
          c(i,j) += HSum(sum);
        }
  }

  
  
  void AddABtSym (FlatMatrix<SIMD<double>> a,
                  FlatMatrix<SIMD<Complex>> b,
                  SliceMatrix<Complex> c)
  {
    size_t ha = a.Height();
    size_t bs = 192;
    if (ha > bs)
      {
        AddABtSym(a.Rows(0,bs), b.Rows(0,bs), c.Rows(0,bs).Cols(0,bs));
        AddABt(a.Rows(bs,ha), b.Rows(0,bs), c.Rows(bs,ha).Cols(0,bs));
        AddABtSym(a.Rows(bs,ha), b.Rows(bs,ha), c.Rows(bs,ha).Cols(bs,ha));
        return;
      }
    
    bs = 96;
    if (ha > bs)
      {
        AddABtSym(a.Rows(0,bs), b.Rows(0,bs), c.Rows(0,bs).Cols(0,bs));
        AddABt(a.Rows(bs,ha), b.Rows(0,bs), c.Rows(bs,ha).Cols(0,bs));
        AddABtSym(a.Rows(bs,ha), b.Rows(bs,ha), c.Rows(bs,ha).Cols(bs,ha));
        return;
      }

    bs = 48;
    if (ha > bs)
      {
        AddABtSym(a.Rows(0,bs), b.Rows(0,bs), c.Rows(0,bs).Cols(0,bs));
        AddABt(a.Rows(bs,ha), b.Rows(0,bs), c.Rows(bs,ha).Cols(0,bs));
        AddABtSym(a.Rows(bs,ha), b.Rows(bs,ha), c.Rows(bs,ha).Cols(bs,ha));
        return;
      }
    bs = 24;
    if (ha > bs)
      {
        AddABtSym(a.Rows(0,bs), b.Rows(0,bs), c.Rows(0,bs).Cols(0,bs));
        AddABt(a.Rows(bs,ha), b.Rows(0,bs), c.Rows(bs,ha).Cols(0,bs));
        AddABtSym(a.Rows(bs,ha), b.Rows(bs,ha), c.Rows(bs,ha).Cols(bs,ha));
        return;
      }
    
    ThreadRegionTimer reg(timer_addabtdcsym, TaskManager::GetThreadId());
    NgProfiler::AddThreadFlops(timer_addabtdcsym, TaskManager::GetThreadId(),
                               a.Height()*b.Height()*a.Width()*8);
    
    // AddABt (a, b, c);
    size_t da = a.Width();
    size_t db = b.Width();
    size_t wa = a.Width();
    // size_t ha = a.Height();
    size_t hb = b.Height();
    size_t dc = c.Dist();
    if (wa == 0) return;
    
    size_t i = 0;
    for ( ; i+1 < ha; i+=2)
      {
        auto pa1 = &a(i,0);
        auto pa2 = pa1 + da;
        auto pb1 = &b(0,0);
        auto pc = &c(i,0);

        for (size_t j = 0; j <= i; j+=2, pb1 += 2*db)
          {
            auto pb2 = pb1 + db;
            
            SIMD<Complex> sum11(0.0);
            SIMD<Complex> sum21(0.0);
            SIMD<Complex> sum12(0.0);
            SIMD<Complex> sum22(0.0);

            __assume (wa > 0);
            for (size_t k = 0; k < wa; k++)
              {
                sum11 += pa1[k] * pb1[k];
                sum21 += pa2[k] * pb1[k];
                sum12 += pa1[k] * pb2[k];
                sum22 += pa2[k] * pb2[k];
              }

            Complex s11, s21, s12, s22;
            std::tie(s11,s12) = HSum(sum11, sum12);
            std::tie(s21,s22) = HSum(sum21, sum22);

            pc[j] += s11;
            pc[j+1] += s12;            
            pc[j+dc] += s21;
            pc[j+dc+1] += s22;            
          }
      }
    
    if (i < ha)
      for (size_t j = 0; j < hb; j++)
        {
          SIMD<Complex> sum(0.0);
          for (size_t k = 0; k < wa; k++)
            sum += a(i,k) * b(j,k);
          c(i,j) += HSum(sum);
        }
  }
  
  void AddABt (FlatMatrix<SIMD<double>> a,
               FlatMatrix<SIMD<double>> b,
               SliceMatrix<Complex> c)
  {
    constexpr size_t M = 92;
    constexpr size_t N = 64;
    double mem[M*N];
    for (size_t i = 0; i < a.Height(); i += M)
      {
        size_t i2 = min2(a.Height(), i+M);
        for (size_t j = 0; j < b.Height(); j += N)
          {
            size_t j2 = min2(b.Height(), j+N);
            FlatMatrix<double> tempc(i2-i, j2-j, &mem[0]);
            tempc = 0.0;
            AddABt (a.Rows(i,i2), b.Rows(j,j2), tempc);
            c.Rows(i,i2).Cols(j,j2) += tempc;
          }
      }
  }
  
  void AddABtSym (FlatMatrix<SIMD<double>> a,
                  FlatMatrix<SIMD<double>> b,
                  SliceMatrix<Complex> c)
  {
    constexpr size_t N = 92;
    double mem[N*N];
    for (size_t i = 0; i < a.Height(); i += N)
      {
        size_t i2 = min2(a.Height(), i+N);
        for (size_t j = 0; j < i; j += N)
          {
            size_t j2 = min2(b.Height(), j+N);
            FlatMatrix<double> tempc(i2-i, j2-j, &mem[0]);
            tempc = 0.0;
            AddABt (a.Rows(i,i2), b.Rows(j,j2), tempc);
            c.Rows(i,i2).Cols(j,j2) += tempc;
          }
        // j == i
        FlatMatrix<double> tempc(i2-i, i2-i, &mem[0]);
        tempc = 0.0;
        AddABtSym (a.Rows(i,i2), b.Rows(i,i2), tempc);
        c.Rows(i,i2).Cols(i,i2) += tempc;
      }
  }
  
  void AddABt (SliceMatrix<double> a,
               SliceMatrix<double> b,
               SliceMatrix<Complex> c)
  {
    constexpr size_t M = 92;
    constexpr size_t N = 64;
    double mem[M*N];
    for (size_t i = 0; i < a.Height(); i += M)
      {
        size_t i2 = min2(a.Height(), i+M);
        for (size_t j = 0; j < b.Height(); j += N)
          {
            size_t j2 = min2(b.Height(), j+N);
            FlatMatrix<double> tempc(i2-i, j2-j, &mem[0]);
            tempc = 0.0;
            AddABt (a.Rows(i,i2), b.Rows(j,j2), tempc);
            c.Rows(i,i2).Cols(j,j2) += tempc;
          }
      }
  }
  
  void AddABtSym (SliceMatrix<double> a,
                  SliceMatrix<double> b,
                  SliceMatrix<Complex> c)
  {
    constexpr size_t N = 92;
    double mem[N*N];
    for (size_t i = 0; i < a.Height(); i += N)
      {
        size_t i2 = min2(a.Height(), i+N);
        for (size_t j = 0; j < i; j += N)
          {
            size_t j2 = min2(b.Height(), j+N);
            FlatMatrix<double> tempc(i2-i, j2-j, &mem[0]);
            tempc = 0.0;
            AddABt (a.Rows(i,i2), b.Rows(j,j2), tempc);
            c.Rows(i,i2).Cols(j,j2) += tempc;
          }
        // j == i
        FlatMatrix<double> tempc(i2-i, i2-i, &mem[0]);
        tempc = 0.0;
        AddABtSym (a.Rows(i,i2), b.Rows(i,i2), tempc);
        c.Rows(i,i2).Cols(i,i2) += tempc;
      }
  }


  /* ************************** SubAtDB ***************************** */

  static constexpr size_t NA = 128;
  //static constexpr size_t NB = 96;
  static constexpr size_t NK = 128;


  void MyTranspose (SliceMatrix<> a, SliceMatrix<> b)
  {
    size_t j = 0;
    size_t ha = a.Height();
    size_t wa = a.Width();
    size_t da = a.Dist();
    size_t db = b.Dist();
    for ( ; j+4 <= wa; j+=4)
      {
        size_t i = 0;
        for ( ; i+4 <= ha; i+=4)
          {
            double * pa = &a(i,j);
            double * pb = &b(j,i);
            SIMD<double,4> a0(pa);
            SIMD<double,4> a1(pa+1*da);
            SIMD<double,4> a2(pa+2*da);
            SIMD<double,4> a3(pa+3*da);
            SIMD<double,4> b0, b1, b2, b3;
            SIMDTranspose(a0,a1,a2,a3, b0,b1,b2,b3);
            b0.Store(pb);
            b1.Store(pb+1*db);
            b2.Store(pb+2*db);
            b3.Store(pb+3*db);
          }
        for ( ; i < ha; i++)
          {
            double * pa = &a(i,j);
            double * pb = &b(j,i);
            pb[0] = pa[0];
            pb[db] = pa[1];
            pb[2*db] = pa[2];
            pb[3*db] = pa[3];
          }
      }
    for ( ; j < wa; j++)
      b.Row(j) = a.Col(j);
  }

  // scale every row from a .. 
  void MyTransposeScaleNeg (SliceMatrix<> a, SliceMatrix<> b,
                            SliceVector<> d)
  {
    size_t ha = a.Height();
    size_t wa = a.Width();
    size_t da = a.Dist();
    size_t db = b.Dist();
    size_t j = 0;
    for ( ; j+4 <= ha; j+=4)
      {
        SIMD<double,4> di(-d(j), -d(j+1), -d(j+2), -d(j+3));
        size_t i = 0;
        double * pa = &a(j,0);
        double * pb = &b(0,j);
        for ( ; i+4 <= wa; i+=4, pa += 4, pb += 4*db)
          {
            SIMD<double,4> a0(pa);
            SIMD<double,4> a1(pa+1*da);
            SIMD<double,4> a2(pa+2*da);
            SIMD<double,4> a3(pa+3*da);
            SIMD<double,4> b0, b1, b2, b3;
            SIMDTranspose(a0,a1,a2,a3, b0,b1,b2,b3);
            (b0*di).Store(pb);
            (b1*di).Store(pb+1*db);
            (b2*di).Store(pb+2*db);
            (b3*di).Store(pb+3*db);
          }
        for ( ; i < wa; i++, pa++, pb+=db)
          {
            SIMD<double,4> b0(pa[0], pa[1*da], pa[2*da], pa[3*da]);
            (b0*di).Store(pb);
          }
      }
    for ( ; j < ha; j++)
      b.Col(j) = (-d(j)) * a.Row(j);
  }


  
  void SubAtDB_BP (SliceMatrix<double> a,
                   SliceVector<double> diag,
                   SliceMatrix<double> b, SliceMatrix<double> c)
  {
    // constexpr size_t SW = SIMD<double>::Size();
    alignas (64) double mema[NA*NK];
    // SIMD<double> memb[3*NK];
    // size_t na = a.Width();
    // size_t nb = b.Width();
    // size_t ha = a.Height();
    
    // loca = Trans(a);
    // for (size_t i = 0; i < loca.Width(); i++)
    // loca.Col(i) *= -diag(i);
    // c += loca * b;
    // return;

#ifdef __AVX512F__
    // constexpr size_t HA = 6;
#else
    // constexpr size_t HA = 4;
#endif

    // size_t da = NA;
    // size_t db = b.Dist();
    // double * pc = c.Data();

    SliceMatrix<> loca(a.Width(), a.Height(), NA, &mema[0]);
    MyTransposeScaleNeg (a, loca, diag);
    c += loca * b;
    return;
    
    /*    
    size_t j = 0;
    for ( ; j+3*SW <= nb; j += 3*SW)
      {
        for (size_t i = 0; i < b.Height(); i++)
          {
            memb[3*i  ] = SIMD<double>(&b(i,j));
            memb[3*i+1] = SIMD<double>(&b(i,j+SW));
            memb[3*i+2] = SIMD<double>(&b(i,j+2*SW));
          }
        
        double * pc =&c(0,j);
        double * pa = &mema[0];
        size_t k = 0;
        for ( ; k+HA <= na; k += HA, pa += HA*da, pc += HA * c.Dist())
          MatKernelMultAB<HA,3,ADD> (ha, pa, da, &memb[0], 3, pc, c.Dist());
        switch (na-k) 
          {
          case 0: break;
          case 1: MatKernelMultAB<1,3,ADD> (ha, pa, da, &memb[0], 3, pc, c.Dist()); break;
          case 2: MatKernelMultAB<2,3,ADD> (ha, pa, da, &memb[0], 3, pc, c.Dist()); break;
          case 3: MatKernelMultAB<3,3,ADD> (ha, pa, da, &memb[0], 3, pc, c.Dist()); break;
          case 4:
            if (HA > 4)
              MatKernelMultAB<4,3,ADD> (ha, pa, da, &memb[0], 3, pc, c.Dist());
            break;
          case 5:
            if (HA > 5)
              MatKernelMultAB<5,3,ADD> (ha, pa, da, &memb[0], 3, pc, c.Dist());
            break;
          default: ;
          }
      }

    if (j == nb) return;
    SliceMatrix<> locb(b.Height(), nb-j, 3*SW, (double*)&memb[0]);
    locb = b.Cols(j, nb);
    pc =&c(0,j);
    double * pa = &mema[0];    
    size_t k = 0;
    for ( ; k+HA <= na; k += HA, pa += HA*da, pc += HA * c.Dist())
      MatKernel2AddAB<HA,ADD> (ha, nb-j, pa, da, &locb(0), 3*SW, pc, c.Dist());
    switch (na-k) 
      {
      case 0: break;
      case 1: MatKernel2AddAB<1,ADD> (ha, nb-j, pa, da, &locb(0), 3*SW, pc, c.Dist()); break;
      case 2: MatKernel2AddAB<2,ADD> (ha, nb-j, pa, da, &locb(0), 3*SW, pc, c.Dist()); break;
      case 3: MatKernel2AddAB<3,ADD> (ha, nb-j, pa, da, &locb(0), 3*SW, pc, c.Dist()); break;
      case 4:
        if (HA > 4)
          MatKernel2AddAB<4,ADD> (ha, nb-j, pa, da, &locb(0), 3*SW, pc, c.Dist());
        break;
      case 5:
        if (HA > 5)
          MatKernel2AddAB<5,ADD> (ha, nb-j, pa, da, &locb(0), 3*SW, pc, c.Dist());
        break;
      default: ;
      }
*/

    
    /*
    size_t k = 0;
    for ( ; k+HA <= na; k += HA, pa += HA*da, pc += HA * c.Dist())
      MatKernel2AddAB<HA,ADD> (ha, nb, pa, da, &b(0), db, pc, c.Dist());
    switch (na-k) 
      {
      case 0: break;
      case 1: MatKernel2AddAB<1,ADD> (ha, nb, pa, da, &b(0), db, pc, c.Dist()); break;
      case 2: MatKernel2AddAB<2,ADD> (ha, nb, pa, da, &b(0), db, pc, c.Dist()); break;
      case 3: MatKernel2AddAB<3,ADD> (ha, nb, pa, da, &b(0), db, pc, c.Dist()); break;
      case 4:
        if (HA > 4)
          MatKernel2AddAB<4,ADD> (ha, nb, pa, da, &b(0), db, pc, c.Dist()); break;
      case 5:
        if (HA > 5)
          MatKernel2AddAB<5,ADD> (ha, nb, pa, da, &b(0), db, pc, c.Dist()); break;
      default: ;
      }
    */
  }
  
  void SubAtDB_PM (SliceMatrix<double> a,
                   SliceVector<double> diag,
                   SliceMatrix<double> b, SliceMatrix<double> c)
  {
    for (size_t i = 0; i < a.Height(); i += NK)
      {
        size_t i2 = min2(i+NK, a.Height());
        SubAtDB_BP (a.Rows(i,i2), diag.Range(i,i2), b.Rows(i,i2), c);
      }
  }
  
  void SubAtDB (SliceMatrix<double> a,
                SliceVector<double> diag,
                SliceMatrix<double> b, SliceMatrix<double> c)
  {
    for (size_t i = 0; i < a.Width(); i += NA)
      {
        size_t i2 = min2(i+NA, a.Width());
        SubAtDB_PM (a.Cols(i,i2), diag, b, c.Rows(i,i2));
      }
  }

  // ************************ SubADBt ******************************
  
  void SubADBt (SliceMatrix<double> a,
                SliceVector<double> diag,
                SliceMatrix<double> b, SliceMatrix<double> c)
  {
    static Timer t("SubADBt"); RegionTimer r(t);
    t.AddFlops(diag.Size()*c.Height()*c.Width());
    constexpr size_t N = 128;
    constexpr size_t M = 128;
    double memb[N*M];

    for (size_t i = 0; i < b.Height(); i += N)
      {
        size_t i2 = min2(i+N, b.Height());
        for (size_t j = 0; j < b.Width(); j += M)
          {
            size_t j2 = min2(j+M, b.Width());
            SliceMatrix<> hb(j2-j, i2-i, N, &memb[0]);
            hb = Trans(b.Rows(i,i2).Cols(j,j2));
            for (size_t k : Range(j2-j))
              hb.Row(k) *= -diag(k+j);
            c.Cols(i,i2) += a.Cols(j,j2) * hb;
          }
      }

    /*
    for (size_t i = 0; i < a.Height(); i++)
      for (size_t j = 0; j < b.Height(); j++)
        for (size_t k = 0; k < a.Width(); k++)
          c(i,j) -= diag(k) * a(i,k) * b(j,k);
    */
  }


  void ScaleCols (SliceMatrix<double,RowMajor> a, BareSliceVector<double> diag)
  {
    static Timer t("ScaleCols, RowMajor"); RegionTimer r(t);
    t.AddFlops(a.Height()*a.Width());

    for (size_t i = 0; i < a.Width(); i++)
      a.Col(i) *= diag(i);
  }

  void ScaleCols (SliceMatrix<double,ColMajor> a, BareSliceVector<double> diag)
  {
    static Timer t("ScaleCols, ColMajor"); RegionTimer r(t);
    t.AddFlops(a.Height()*a.Width());
    
    for (size_t i = 0; i < a.Width(); i++)
      a.Col(i) *= diag(i);
  }

  
  
  // ************************************** Complex ADB^t *********************
  
  /*
  
  void CopyMatrixInScaleRows (size_t h, size_t w,
                              Complex * ps, size_t dists,
                              Complex * pd, size_t distd,
                              Complex * pscale, size_t distscale)
  {
    for (size_t i = 0; i < h; i++, ps += dists, pd += distd, pscale += distscale)
      {
        Complex scale = *pscale;
        for (size_t j = 0; j < w; j++)
          pd[j] = scale * ps[j];
      }
  }
  
  */
  
  void CopyMatrixInScaleRows (size_t h, size_t w,
                              Complex * ps, size_t dists,
                              Complex * pd, size_t distd,
                              Complex * pscale, size_t distscale)
  {
    for (size_t i = 0; i < h; i++, ps += dists, pd += distd, pscale += distscale)
      {
        SIMD<Complex> scale (*pscale);
        size_t j = 0;
        size_t WS = SIMD<double>::Size();
        for ( ; j+4*WS <= w; j+=4*WS)
          {
            SIMD<Complex> val1, val2, val3, val4;
            val1.LoadFast(ps+j);
            val2.LoadFast(ps+j+WS);
            val3.LoadFast(ps+j+2*WS);
            val4.LoadFast(ps+j+3*WS);
            val1 = val1 * scale;
            val2 = val2 * scale;
            val3 = val3 * scale;
            val4 = val4 * scale;
            val1.StoreFast(pd+j);
            val2.StoreFast(pd+j+WS);
            val3.StoreFast(pd+j+2*WS);
            val4.StoreFast(pd+j+3*WS);
          }
        for ( ; j+WS <= w; j+=WS)
          {
            SIMD<Complex> val;
            val.LoadFast(ps+j);
            val = val * scale;
            val.StoreFast(pd+j);
          }
        SIMD<Complex> val;
        val.LoadFast(ps+j, w-j);
        val = val * scale;
        val.StoreFast(pd+j, w-j);
      }
  }  



  

  
  void KernelScal4x4Trans (Complex * pa, size_t da,
                           Complex * pb, size_t db,
                           Complex * pc, size_t dc,
                           size_t ninner)
  {
    SIMD<Complex> sum1, sum2, sum3, sum4;
    sum1.LoadFast (pc);
    sum2.LoadFast (pc+dc);
    sum3.LoadFast (pc+2*dc);
    sum4.LoadFast (pc+3*dc);
    for (size_t i = 0; i < ninner; i++, pa += da, pb += db)
      {
        SIMD<Complex> b1;
        b1.LoadFast(pb);
        sum1 = sum1 - SIMD<Complex> (pa[0]) * b1;
        sum2 = sum2 - SIMD<Complex> (pa[1]) * b1;
        sum3 = sum3 - SIMD<Complex> (pa[2]) * b1;
        sum4 = sum4 - SIMD<Complex> (pa[3]) * b1;
      }
    sum1.StoreFast(pc);
    sum2.StoreFast(pc+dc);
    sum3.StoreFast(pc+2*dc);
    sum4.StoreFast(pc+3*dc);
  }

  void KernelScal4x4Trans (Complex * pa, size_t da,
                           Complex * pb, size_t db,
                           Complex * pc, size_t dc,
                           size_t ninner, int mask)
  {
    SIMD<Complex> sum1, sum2, sum3, sum4;
    sum1.LoadFast (pc, mask);
    sum2.LoadFast (pc+dc, mask);
    sum3.LoadFast (pc+2*dc, mask);
    sum4.LoadFast (pc+3*dc, mask);
    for (size_t i = 0; i < ninner; i++, pa += da, pb += db)
      {
        SIMD<Complex> b1;
        b1.LoadFast(pb, mask);
        sum1 = sum1 - SIMD<Complex> (pa[0]) * b1;
        sum2 = sum2 - SIMD<Complex> (pa[1]) * b1;
        sum3 = sum3 - SIMD<Complex> (pa[2]) * b1;
        sum4 = sum4 - SIMD<Complex> (pa[3]) * b1;
      }
    sum1.StoreFast(pc, mask);
    sum2.StoreFast(pc+dc, mask);
    sum3.StoreFast(pc+2*dc, mask);
    sum4.StoreFast(pc+3*dc, mask);
  }
  

   
  void KernelScal1x4Trans (Complex * pa, size_t da,
                           Complex * pb, size_t db,
                           Complex * pc, size_t dc,
                           size_t ninner)
  {
    SIMD<Complex> sum1;
    sum1.LoadFast (pc);
    for (size_t i = 0; i < ninner; i++, pa += da, pb += db)
      {
        SIMD<Complex> b1;
        b1.LoadFast(pb);
        sum1 = sum1 - SIMD<Complex> (*pa) * b1;
      }
    sum1.StoreFast(pc);
  }
  
  void KernelScal1x4Trans (Complex * pa, size_t da,
                           Complex * pb, size_t db,
                           Complex * pc, size_t dc,
                           size_t ninner, int mask)
  {
    SIMD<Complex> sum1;
    sum1.LoadFast (pc, mask);
    for (size_t i = 0; i < ninner; i++, pa += da, pb += db)
      {
        SIMD<Complex> b1;
        b1.LoadFast(pb, mask);
        sum1 = sum1 - SIMD<Complex> (*pa) * b1;
      }
    sum1.StoreFast(pc, mask);
  }
  
  void MySubAtDB_BB (
                      Complex * pa, size_t da,
                      Complex * pb, size_t db,
                      Complex * pc, size_t dc,
                      size_t na, size_t nb, size_t ninner
                      )
  {
    size_t WS = SIMD<double>::Size();
    size_t i = 0;
    for ( ; i+4 <= na; i+=4, pa += 4, pc += 4*dc)
      {
        size_t j = 0;
        for ( ; j+WS <= nb; j+=WS)
          KernelScal4x4Trans (pa, da, pb+j, db, pc+j, dc, ninner);
        if (j < nb)
          KernelScal4x4Trans (pa, da, pb+j, db, pc+j, dc, ninner, nb-j);          
          /* maybe better for AVX, but not portable to 512
          {
          for ( ; j < nb; j++)
            Complex tmpc[4] = { pc[j], pc[dc+j], pc[2*dc+j], pc[3*dc+j] };
            KernelScal1x4Trans (pb+j, db, pa, da, tmpc, 1, ninner);
            pc[j] = tmpc[0];
            pc[dc+j] = tmpc[1];
            pc[2*dc+j] = tmpc[2];
            pc[3*dc+j] = tmpc[3];
          }
          */
      }
    for ( ; i < na; i++, pa ++, pc += dc)
      {
        size_t j = 0;
        for ( ; j+WS <= nb; j+=WS)
          KernelScal1x4Trans (pa, da, pb+j, db, pc+j, dc, ninner);
        if (j < nb)
          KernelScal1x4Trans (pa, da, pb+j, db, pc+j, dc, ninner, nb-j);
      }
  }

  /*

  void MySubAtDB_BB (
                     Complex * pa, size_t da,
                     Complex * pb, size_t db,
                     Complex * pc, size_t dc,
                     size_t na, size_t nb, size_t ninner
                     )
  {
  // SliceMatrix<Complex> a(ninner, na, da, pa);
  // SliceMatrix<Complex> b(ninner, nb, db, pb);
  // SliceMatrix<Complex> c(na, nb, dc, pc);
  // c -= Trans(a) * b; //  | Lapack;

    for (size_t i = 0; i < na; i++)
      for (size_t j = 0; j < nb; j++)
        {
          Complex sum = pc[i*dc+j];
          for (size_t k = 0; k < ninner; k++)
            sum -= pa[k*da+i] * pb[k*db+j];
          pc[i*dc+j] = sum;
        }
  }
  */


  constexpr size_t CNA = 32;
  constexpr size_t CNB = 32;
  constexpr size_t CNK = 32;
  
  void MySubAtDB_BP (SliceMatrix<Complex> a,
                     SliceVector<Complex> diag,
                     SliceMatrix<Complex> b, SliceMatrix<Complex> c)
  {
    // alignas (64) Complex mema[CNA*CNK];   // slow !!!
    alignas (64) double mema[2*CNA*CNK];
    size_t na = a.Width();
    size_t nb = b.Width();
    size_t k = a.Height();
    
    CopyMatrixInScaleRows (k, na,
                           a.Data(), a.Dist(), (Complex*)&mema[0], CNA,
                           diag.Data(), diag.Dist());

    size_t i = 0;
    constexpr size_t bs = CNB;
    for ( ; i+bs <= nb; i += bs)
      MySubAtDB_BB ((Complex*)mema, CNA, &b(size_t(0),i), b.Dist(), &c(size_t(0),i), c.Dist(), na, bs,k);
    if (i < nb)
      MySubAtDB_BB ((Complex*)mema, CNA, &b(size_t(0),i), b.Dist(), &c(size_t(0),i), c.Dist(), na, nb-i, k);    
  }
  
  void MySubAtDB_PM (SliceMatrix<Complex> a,
                     SliceVector<Complex> diag,
                     SliceMatrix<Complex> b, SliceMatrix<Complex> c)
  {
    size_t k = a.Height();
    size_t i = 0;
    constexpr size_t bs = CNK;
    for ( ; i+bs <= k; i += bs) 
      MySubAtDB_BP (a.Rows(i,i+bs), diag.Range(i,i+bs), b.Rows(i,i+bs), c);
    if (i < k)
      MySubAtDB_BP (a.Rows(i,k), diag.Range(i,k), b.Rows(i,k), c);      
  }
  
  void SubAtDB (SliceMatrix<Complex> a,
                SliceVector<Complex> diag,
                SliceMatrix<Complex> b, SliceMatrix<Complex> c)
  {
    size_t na = a.Width();
    size_t i = 0;
    constexpr size_t bs = CNA;
    for ( ; i+bs <= na; i += bs)
      MySubAtDB_PM (a.Cols(i,i+bs), diag, b, c.Rows(i,i+bs));
    if (i < na)
      MySubAtDB_PM (a.Cols(i,na), diag, b, c.Rows(i,na));
  }


  void PairwiseInnerProduct (size_t n, FlatArray<double*> x, FlatArray<double*> y, BareSliceMatrix<double> ip) {

    #ifdef __AVX512F__
        constexpr size_t HA = 6;
    #else
      constexpr size_t HA = 3;
    #endif

    int x_size = x.Size();
    int y_size = y.Size();

    size_t j = 0;
    for (; j+HA <= x_size; j+=HA) {

      size_t k = 0;
      for (; k+4 <= y_size; k+=4) {

        auto res = MultiVecScalAB<HA, 4>(n, &x[j], &y[k]);

        Iterate<HA> ([&] (auto i) {
          auto sum = get<i.value>(res);
          sum.Store(&ip(j + i.value, k));
        });

      }

      for (; k+2 <= y_size; k+=2) {

        auto res = MultiVecScalAB<HA, 2>(n, &x[j], &y[k]);

        Iterate<HA> ([&] (auto i) {
          auto sum = get<i.value>(res);
          sum.Store(&ip(j + i.value, k));
        });

      }

      for (; k < y_size; k++) {

        auto res = MultiVecScalAB<HA, 1>(n, &x[j], &y[k]);

        Iterate<HA> ([&] (auto i) {
          auto sum = get<i.value>(res);
          ip(j + i.value, k) = sum;
        });

      }
    }

    for ( ; j < x_size; j++)
    {

      size_t k = 0;
      for (; k+4 <= y_size; k+=4) {

        auto res = MultiVecScalAB<1, 4>(n, &x[j], &y[k]);

        auto sum0 = get<0>(res);
        sum0.Store(&ip(j, k));
      }

      for (; k < y_size; k++) {

	       auto res = MultiVecScalAB<1, 1>(n, &x[j], &y[k]);

         auto sum0 = get<0>(res);
         ip(j, k) = sum0;
      }
    }

  }


  void PairwiseInnerProduct (size_t n, FlatArray<Complex*> x, FlatArray<Complex*> y, BareSliceMatrix<Complex> ip, bool conj) {

    #ifdef __AVX512F__
        constexpr size_t HA = 6;
    #else
      constexpr size_t HA = 3;
    #endif

    int x_size = x.Size();
    int y_size = y.Size();

    if (conj) {
      size_t j = 0;
      for (; j+HA <= x_size; j+=HA) {

        size_t k = 0;
        for (; k+2 <= y_size; k+=2) {
          MultiVecScalC<HA, 2, 1>(n, &x[j], &y[k], &ip(j,k), ip.Dist());
        }
        for (; k < y_size; k++) {
          MultiVecScalC<HA, 1, 1>(n, &x[j], &y[k], &ip(j,k), ip.Dist());
        }

      }

      for ( ; j < x_size; j++)
      {

        size_t k = 0;
        for (; k+8 <= y_size; k+=8) {
          MultiVecScalC<1, 8, 1>(n, &x[j], &y[k], &ip(j,k), ip.Dist());
        }
        for (; k+4 <= y_size; k+=4) {
          MultiVecScalC<1, 4, 1>(n, &x[j], &y[k], &ip(j,k), ip.Dist());
        }
        for (; k < y_size; k++) {
          MultiVecScalC<1, 1, 1>(n, &x[j], &y[k], &ip(j,k), ip.Dist());
        }
      }
    }
    else {

      size_t j = 0;
      for (; j+HA <= x_size; j+=HA) {

        size_t k = 0;
        for (; k+2 <= y_size; k+=2) {
          MultiVecScalC<HA, 2, 0>(n, &x[j], &y[k], &ip(j,k), ip.Dist());
        }
        for (; k < y_size; k++) {
          MultiVecScalC<HA, 1, 0>(n, &x[j], &y[k], &ip(j,k), ip.Dist());
        }
      }

      for ( ; j < x_size; j++)
      {

        size_t k = 0;
        for (; k+8 <= y_size; k+=8) {
          MultiVecScalC<1, 8, 0>(n, &x[j], &y[k], &ip(j,k), ip.Dist());
        }
        for (; k+4 <= y_size; k+=4) {
          MultiVecScalC<1, 4, 0>(n, &x[j], &y[k], &ip(j,k), ip.Dist());
        }
        for (; k < y_size; k++) {
          MultiVecScalC<1, 1, 0>(n, &x[j], &y[k], &ip(j,k), ip.Dist());
        }
      }
    }

  }


    // x_i += sum_j a(j,i) y_j
  void MultiVectorAdd (size_t n, FlatArray<double*> x, FlatArray<double*> y, BareSliceMatrix<double> a) {

    constexpr int Hx = 6;
    constexpr int Hy = 6;

    size_t i = 0;
    for (; i + Hy <= y.Size(); i += Hy) {

      size_t j = 0;
      for (; j + Hx <= x.Size(); j+=Hx) {
        MultiScaleAdd<Hx,Hy>(n, x+j, y+i, &a(i,j), a.Dist());
      }
      for (; j + 2 <= x.Size(); j+=2) {
        MultiScaleAdd<2,Hy>(n, x+j, y+i, &a(i,j), a.Dist());
      }
      for (; j < x.Size(); j++) {
        MultiScaleAdd<1,Hy>(n, x+j, y+i, &a(i,j), a.Dist());
      }

    }
    for (; i + 2 <= y.Size(); i += 2) {
      size_t j = 0;
      for (; j + Hx <= x.Size(); j+=Hx) {
        MultiScaleAdd<Hx,2>(n, x+j, y+i, &a(i,j), a.Dist());
      }
      for (; j + 2 <= x.Size(); j+=2) {
        MultiScaleAdd<2,2>(n, x+j, y+i, &a(i,j), a.Dist());
      }
      for (; j < x.Size(); j++) {
        MultiScaleAdd<1,2>(n, x+j, y+i, &a(i,j), a.Dist());
      }
    }
    for (; i < y.Size(); i++) {

      size_t j = 0;
      for (; j + Hx <= x.Size(); j+=Hx) {
        MultiScaleAdd<Hx,1>(n, x+j, y+i, &a(i,j), a.Dist());
      }
      for (; j + 2 <= x.Size(); j+=2) {
        MultiScaleAdd<2,1>(n, x+j, y+i, &a(i,j), a.Dist());
      }
      for (; j < x.Size(); j++) {
        MultiScaleAdd<1,1>(n, x+j, y+i, &a(i,j), a.Dist());
      }

    }
  }


    void MultiVectorAdd (size_t n, FlatArray<Complex*> x, FlatArray<Complex*> y, BareSliceMatrix<Complex> a) {

    constexpr int Hx = 3;
    constexpr int Hy = 4;

    size_t i = 0;
    for (; i + Hy <= y.Size(); i += Hy) {

      size_t j = 0;
      for (; j + Hx <= x.Size(); j+=Hx) {
        MultiScaleAddC<Hx,Hy>(n, x+j, y+i, &a(i,j), a.Dist());
      }
      for (; j + 2 <= x.Size(); j+=2) {
        MultiScaleAddC<2,Hy>(n, x+j, y+i, &a(i,j), a.Dist());
      }
      for (; j < x.Size(); j++) {
        MultiScaleAddC<1,Hy>(n, x+j, y+i, &a(i,j), a.Dist());
      }

    }
    for (; i < y.Size(); i++) {

      size_t j = 0;
      for (; j + Hx <= x.Size(); j+=Hx) {
        MultiScaleAddC<Hx,1>(n, x+j, y+i, &a(i,j), a.Dist());
      }
      for (; j + 2 <= x.Size(); j+=2) {
        MultiScaleAddC<2,1>(n, x+j, y+i, &a(i,j), a.Dist());
      }
      for (; j < x.Size(); j++) {
        MultiScaleAddC<1,1>(n, x+j, y+i, &a(i,j), a.Dist());
      }

    }

  }




  

  /**************** timings *********************** */

  extern void MultUL (SliceMatrix<> A);
  extern void LapackSVD (SliceMatrix<> A,
                         SliceMatrix<double, ColMajor> U,
                         SliceMatrix<double, ColMajor> V);

  
  list<tuple<string,double>> Timing (int what, size_t n, size_t m, size_t k, bool lapack, size_t maxits)
  {
    if (what < 0)
      {
        cout << "Available options timings are:\n"
          "-1 .. this help\n"
          "0 ... run all timings\n"
          "1 ... A = B,   A,B = n*m,   A = aligned, fixed dist\n"
          "2 ... A = 0,   A = n*m,     but sliced\n"
          "3 ... A = B^t, A = n*m, \n"
          "5 ... y = A*x,   A = n*m\n"
          "6 ... y = A^t*x,   A = n*m\n"
          "7 ... y += A^t*x(ind),   A = n*m\n"
          "10 .. C = A * B,   A=n*m, B=m*k, C=n*k\n"
          "11 .. C += A * B,   A=n*m, B=m*k, C=n*k\n"
          // "20 .. C = A * B    A=n*m, B=n*k', C=n*k', k'=round(k), B aligned\n"
          "20 .. X = T * X       T=n*n triangular, X=n*m "
          "21 .. X = T^-1 * X     T=n*n triangular, X=n*m "
          "22 .. T^-1             T=n*n triangular"
          "50 .. C += A * B^t,   A=n*k, B=m*k, C=n*m\n"
          "51 .. C += A * B^t,   A=n*k, B=m*k, C=n*m,  A,B aligned\n"
          "52 .. C = A * B^t,   A=n*k, B=m*k, C=n*m\n"
          "60 .. C -= A^t * D B,  A=n*k, B=n*m, C = k*m, D=diag\n"
          "61 .. C = A^t B,  A=n*k, B=n*m, C = k*m\n"
          "70 .. C += A B^t,  A=n*k, B=m*k, C = n*m, A,B SIMD\n"
          "100.. MultAddKernel  C += A * B,  A=4*n, B=n*3SW\n"
          "101.. MultAddKernel  C += A * B,  A=4*n, B=n*3SW, B aligned\n"
          "110.. MultAddKernel2  C += A * B,  A=4*n, B=n*m, m multiple of 3*SW\n"
          "111.. MultAddKernel2  C += A * B,  A=4*n, B=n*m, m multiple of 3*SW, B aligned\n"
          "150.. ScalKernel     C = A * B^t,  A=4*n, B = 3*n\n"
          "151.. ScalKernel     C = A * B^t,  A=4*n, B = 3*n\n, A,B aligned\n"
          "200.. CalcInverse        A = nxn\n"
          "201.. CalcInverse by LU  A = nxn\n"          
          "205.. LDL                A = nxn\n"
          "210.. CalcInverseLapack  A = nxn\n"
          "300.. CalcSVD            A = nxn\n"
             << endl;
        return list<tuple<string,double>>();
      }

    list<tuple<string,double>> timings;
    constexpr int SW = SIMD<double>::Size();
    if (what == 0 || what == 1)
      {
        // A = B
        constexpr size_t WA = 128;
        if (m > WA)
          {
            m = WA;
            cout << "max width = " << WA << endl;
          }
        Matrix<> b(n,m);
        STACK_ARRAY(SIMD<double>, mema, n*WA/SIMD<double>::Size());
        FlatMatrix<SIMD<double>> a(n,WA/SIMD<double>::Size(),&mema[0]);
        b = 1;
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("Copy matrix, packed dest");
          t.Start();
          for (size_t j = 0; j < its; j++)
            CopyMatrixIn(n,m, &b(0,0), m, &a(0,0), a.Width());
          t.Stop();
          cout << "Lapack GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Copy matrix, packed dest", 1e-9 * n*m*its / t.GetTime()));
        }
      }


    if (what == 0 || what == 2)
      {
        // A = 0
        Matrix<> a(n,m);
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("Zero matrix, packed dest");
          t.Start();
          for (size_t j = 0; j < its; j++)
            a.Rows(0,n).Cols(0,m) = j;
          t.Stop();
          cout << "Zero matrix GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Zero matrix", 1e-9 * n*m*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 3)
      {
        // A = B^t
        Matrix<> a(n,m), b(m,n);
        b = 1;
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("Matrix Transpose");
          t.Start();
          for (size_t j = 0; j < its; j++)
            TransposeMatrix(b, a);
          t.Stop();
          cout << "Lapack GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Transpose matrix", 1e-9 * tot*its / t.GetTime()));
        }
      }


    
    if (what == 0 || what == 5)
      {
        // y = A*x
        Matrix<> a(n,m);
        Vector<> x(m), y(n);
        a = 1; x = 2;
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          MultMatVec(a,x,y);
          if (L2Norm(a*x-y) > 1e-8)
            throw Exception("MultMatVec is faulty");
          Timer t("y = A*x");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MultMatVec(a,x,y);
          t.Stop();
          cout << "MultMatVec GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatVec", 1e-9 * n*m*its / t.GetTime()));
        }
        {
          Timer t("y = A*x, Lapack");
          t.Start();
          for (size_t j = 0; j < its; j++)
            LapackMultAx (a, x, y);
          t.Stop();
          cout << "MultMatVec Lapack GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatVecLapack", 1e-9 * n*m*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 6)
      {
        // y = A*x
        Matrix<> a(n,m);
        Vector<> x(n), y(m);
        a = 1; x = 2;
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("y = A*x");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MultMatTransVec(a,x,y);
          t.Stop();
          cout << "MultMatTransVec GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatVec", 1e-9 * n*m*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 7)
      {
        // y = A*x
        Matrix<> a(n,m);
        Vector<> x(1000), y(m);
        Array<int> index(n);
        for (size_t i = 0; i < n; i++)
          index[i] = (17*i)%1000;
        a = 1; x = 2; y = 0;
        double tot = n*m;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("y = A*x");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MultAddMatTransVecIndirect(1, a,x,y, index);
          t.Stop();
          cout << "MultAddMatTransVecIndirect GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultAddMatVecIndirect", 1e-9 * n*m*its / t.GetTime()));
        }
      }

    
    if (what == 0 || what == 10)
      {
        // C=A*B
        Matrix<> a(n,m), b(m,k), c(n,k);
        a = 1; b = 2;
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < m; j++)
            a(i,j) = sin(i+1) * cos(j);
        for (size_t i = 0; i < m; i++)
          for (size_t j = 0; j < k; j++)
            b(i,j) = cos(i+3) * cos(j);
        
        double tot = n*m*k;
        size_t its = 1e9 / tot + 1;
        // MultMatMat(a,b,c);
        c = a * b;
        double err = L2Norm(a*b-c);
        if (err > 1e-8)
          throw Exception("MultMatMat is faulty");
        
        {
          Timer t("C = A*B");
          t.Start();
          if (!lapack)
            for (size_t j = 0; j < its; j++)
              // MultMatMat(a,b,c);
              c = a*b;
          else
            for (size_t j = 0; j < its; j++)
              c = a*b | Lapack;
          t.Stop();
          cout << "MultMatMat GFlops = " << 1e-9 * n*m*k*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatMat", 1e-9 * n*m*k*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 11)
      {
        // C=A*B
        Matrix<> a(n,m), b(m,k), c(n,k);
        a = 1; b = 2;
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < m; j++)
            a(i,j) = sin(i+1) * cos(j);
        for (size_t i = 0; i < m; i++)
          for (size_t j = 0; j < k; j++)
            b(i,j) = cos(i+3) * cos(j);
        c = 0.0;
        double tot = n*m*k;
        size_t its = 1e9 / tot + 1;
        // MultMatMat(a,b,c);
        {
          Timer t("C += A*B");
          t.Start();
          if (!lapack)
            for (size_t j = 0; j < its; j++)
              c += a*b;
          else
            for (size_t j = 0; j < its; j++)
              c += a*b | Lapack;
          t.Stop();
          cout << "MultMatMat GFlops = " << 1e-9 * n*m*k*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatMat", 1e-9 * n*m*k*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 20)
      {
        Matrix<> a(n,n), b(n,m);
        a = 1; b = 2;
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < n; j++)
            a(i,j) = sin(i+1) * cos(j);
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < m; j++)
            b(i,j) = cos(i+3) * cos(j);
        Matrix<> saveb = b;
        
        double tot = n*n*m/2;
        size_t its = 1e9 / tot + 1;

        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularMult<LowerLeft> (a, b);
            }
          t.Stop();
          cout << "TriangularMult<L> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularMult<L>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }
        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularMult<UpperRight> (a, b);
            }
          t.Stop();
          cout << "TriangularMult<R> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularMult<R>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }

        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularMult<LowerLeft,Normalized> (a, b);
            }
          t.Stop();
          cout << "TriangularMult<L,N> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularMult<L,N>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }
        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularMult<UpperRight,Normalized> (a, b);
            }
          t.Stop();
          cout << "TriangularMult<R,N> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularMult<R,N>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }

      }
    if (what == 0 || what == 21)
      {
        Matrix<> a(n,n), b(n,m);
        a = 1; b = 2;
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < n; j++)
            a(i,j) = sin(i+1) * cos(j);
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < m; j++)
            b(i,j) = cos(i+3) * cos(j);
        Matrix<> saveb = b;
        
        double tot = n*n*m/2;
        size_t its = 1e9 / tot + 1;
        // MultMatMat(a,b,c);
        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularSolve<LowerLeft> (a, b);
            }
          t.Stop();
          cout << "TriangularSolve<L> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularSolve<L>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }
        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularSolve<UpperRight> (a, b);
            }
          t.Stop();
          cout << "TriangularSolve<R> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularSolve<R>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }

        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularSolve<LowerLeft, Normalized> (a, b);
            }
          t.Stop();
          cout << "TriangularSolve<L,N> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularSolve<L,N>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }
        {
          Timer t("X = L * X");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              b = saveb;
              TriangularSolve<UpperRight, Normalized> (a, b);
            }
          t.Stop();
          cout << "TriangularSolve<R,N> GFlops = " << 1e-9 * n*n*m/2*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularSolve<R,N>", 1e-9 * n*n*m/2*its / t.GetTime()));
        }


      }




    if (what == 0 || what == 22)
      {
        // T^{-1}
        Matrix<> a(n,n);
        a = 1; 
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j < n; j++)
            a(i,j) = sin(i+1) * cos(j);
        Matrix<> savea = a;
        
        double tot = n*n*n/6;
        size_t its = 1e9 / tot + 1;
        if (its > maxits) its = maxits;
        {
          Timer t("L^-1");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              // a = savea;
              TriangularInvert<LowerLeft> (a);
              TriangularInvert<LowerLeft> (a);
            }
          t.Stop();
          cout << "TriangularInvert<L> GFlops = " << 1e-9 * 2*n*n*n/6*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularInvert<L>", 1e-9 * 2*n*n*n/6*its / t.GetTime()));
        }
        {
          Timer t("R^-1");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              // a = savea;
              TriangularInvert<UpperRight> (a);
              TriangularInvert<UpperRight> (a);
            }
          t.Stop();
          cout << "TriangularInvert<R> GFlops = " << 1e-9 * 2*n*n*n/6*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularInvert<R>", 1e-9 * 2*n*n*n/6*its / t.GetTime()));
        }
        {
          Timer t("L^-1");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              // a = savea;
              TriangularInvert<LowerLeft, Normalized> (a);
              TriangularInvert<LowerLeft, Normalized> (a);
            }
          t.Stop();
          cout << "TriangularInvert<LN> GFlops = " << 1e-9 * 2*n*n*n/6*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularInvert<LN>", 1e-9 * 2*n*n*n/6*its / t.GetTime()));
        }
        {
          Timer t("R^-1");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              // a = savea;
              TriangularInvert<UpperRight, Normalized> (a);
              TriangularInvert<UpperRight, Normalized> (a);
            }
          t.Stop();
          cout << "TriangularInvert<RN> GFlops = " << 1e-9 * 2*n*n*n/6*its / t.GetTime() << endl;
          timings.push_back(make_tuple("TriangularInvert<RN>", 1e-9 * 2*n*n*n/6*its / t.GetTime()));
        }

        {
          Timer t("UL");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              a = savea;
              MultUL (a);
            }
          t.Stop();
          cout << "MultUL GFlops = " << 1e-9 * n*n*n/3*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultUL", 1e-9 * n*n*n/3*its / t.GetTime()));
        }
      }





    
    
    if (what == 0 || what == 50)
      {
        // C=A*B^t
        Matrix<> a(n,k), b(m,k), c(n,m);
        a = 1; b = 2;
        c = 0.0;        
        double tot = n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            // AddABt(a,b,c);
            c += a * Trans(b);
          t.Stop();
          cout << "AddABt GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddABt", 1e-9 * tot *its / t.GetTime()));
        }
      }

    if (what == 0 || what == 51)
      {
        // C=A*B^t
        if (k % SW != 0)
          cout << "k should be a multiple of " << SW << endl;
        size_t ks = k/SW;
        Matrix<SIMD<double>> a(n,ks), b(m,ks);
        Matrix<> c(n,m);
        a = SIMD<double>(1); b = SIMD<double>(2);
        c = 0.0;
        double tot = n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            AddABt(SliceMatrix<double> (a.Height(), SW*a.Width(), SW*a.Width(), (double*)&a(0)),
                   SliceMatrix<double> (b.Height(), SW*b.Width(), SW*b.Width(), (double*)&b(0)),
                   // SliceMatrix<double> (AFlatMatrix<double>(b)),
                   c);
          t.Stop();
          cout << "AddABt GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddABt", 1e-9 * tot *its / t.GetTime()));
        }
      }

    if (what == 0 || what == 52)
      {
        // C=A*B^t
        Matrix<> a(n,k), b(m,k), c(n,m);
        a = 1; b = 2;
        c = 0.0;        
        double tot = n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          if (!lapack)
            for (size_t j = 0; j < its; j++)
              c = a * Trans(b);
          else
            for (size_t j = 0; j < its; j++)
              c = a * Trans(b) | Lapack;
          t.Stop();
          cout << "AddABt GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddABt", 1e-9 * tot *its / t.GetTime()));
        }
      }

    
    if (what == 0 || what == 60)
      {
        // C=A*B^t
        Matrix<> a(n,k), b(n,m), c(k,m);
        Vector<> d(n);
        a = 1, b = 1, d = 2;
        c = 0.0;
        double tot = n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C -= A^t*D*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            SubAtDB(a, d, b, c);
          t.Stop();
          cout << "AddAtDB GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddAtDB", 1e-9 * tot *its / t.GetTime()));
        }
      }

    if (what == 0 || what == 61)
      {
        // C=A^t*B
        Matrix<> a(n,k), b(n,m), c(k,m);
        for (size_t i = 0; i < a.Height(); i++)
          for (size_t j = 0; j < a.Width(); j++)
            a(i,j) = sin(i+1) * cos(j);
        for (size_t i = 0; i < b.Height(); i++)
          for (size_t j = 0; j < b.Width(); j++)
            b(i,j) = cos(i+3) * cos(j);
        
        c = 0.0;
        MultAtB (a,b,c);
        if (n <= 1000)
          {
            double err = L2Norm(Trans(a)*b-c);
            if (err > 1e-8)
              throw Exception("MultAtB is faulty");
          }
        double tot = n*m*k;
        size_t its = 5e9 / tot + 1;
        {
          Timer t("C = A^t*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MultAtB(a, b, c);
          t.Stop();
          cout << "MultAtB GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultAtB", 1e-9 * tot *its / t.GetTime()));
        }
        {
          Timer t("C = A^t*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            c = Trans(a)*b;
          t.Stop();
          cout << "C=A^t*B GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("A^T*B", 1e-9 * tot *its / t.GetTime()));
        }
        {
          Timer t("C = A^t*B, block block");
          constexpr size_t BS = 96;
          if (n >= BS && m >= BS && k >= BS)
            {
              tot = BS*BS*BS;
              its = 5e9 / tot+1;
              auto ba = a.Rows(BS).Cols(BS);
              // auto bb = b.Rows(BS).Cols(BS);
              auto bc = c.Rows(BS).Cols(BS);
              // Matrix ba(BS,BS);
              Matrix bb(BS,BS);
              t.Start();
              for (size_t j = 0; j < its; j++)
                MultAtB_intern2<SET,BS> (ba.Height(), ba.Width(), bb.Width(), ba, bb, bc);
              t.Stop();
              cout << "MultAtB - block GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
              timings.push_back(make_tuple("MultAtB - block", 1e-9 * tot *its / t.GetTime()));
            }
        }
      }
    
    if (what == 0 || what == 70)
      {
        // C=A*B^t
        if (k % SW != 0)
          cout << "k should be a multiple of " << SW << endl;
        size_t ks = k/SW;
        Matrix<SIMD<double>> a(n,ks), b(m,ks);
        Matrix<> c(n,m);
        a = SIMD<double>(1); b = SIMD<double>(2);
        c = 0.0;
        double tot = n*m*k;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C += A*Bt, sym");
          t.Start();
          for (size_t j = 0; j < its; j++)
            AddABtSym(a, b, c);
          t.Stop();
          cout << "AddABt, sym GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddABt", 1e-9 * tot *its / t.GetTime()));
        }
      }

    
    if (what == 0 || what == 100)
      {
        // C=A*B
        Matrix<> a(4,n), b(n,3*SW), c(4,3*SW);
        a = 1; b = 2; c = 0;
        double tot = n*4*3*SW;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MatKernelMultAB<4,3,ADD>(n,&a(0), a.Width(), &b(0), b.Width(), &c(0), c.Width());
          t.Stop();
          cout << "MatKernelAddAB 3x4 = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelAddAB", 1e-9 * tot*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 101)
      {
        // C=A*B
        Matrix<> a(4,n), c(4,3*SW);
        Matrix<SIMD<double>> b(n, 3);
        a = 1; b = SIMD<double>(2); c = 0;
        double tot = n*4*3*SW;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            MatKernelMultAB<4,3,ADD>(n,&a(0), a.Width(), &b(0), b.Width(), &c(0), c.Width());
          t.Stop();
          cout << "MatKernelAddAB 3x4, algined GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelAddAB aligned", 1e-9 * tot*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 110)
      {
        // C=A*B
        if (m % (3*SW) != 0)
          cout << "m should be a multiple of 3*SIMD::Size" << endl;
        Matrix<> a(4,n), b(n,m), c(4,m);
        a = 1; b = 2; c = 0;
        double tot = n*4*m;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            for (size_t i = 0; i+3*SW <= m; i += 3*SW)
              MatKernelMultAB<4,3,ADD>(n,&a(0), a.Width(), &b(i), b.Width(), &c(i), c.Width());
          t.Stop();
          cout << "MatKernel2AddAB 3x4 = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelAddAB", 1e-9 * tot*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 111)
      {
        // C=A*B
        if (m % (3*SW) != 0)
          cout << "m should be a multiple of 3*SIMD::Size" << endl;
        Matrix<> a(4,n), c(4,m);
        Matrix<SIMD<double>> b(n, m/SW);
        a = 1; b = SIMD<double>(2); c = 0;
        double tot = n*4*m;
        size_t its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            for (size_t i = 0; i+3*SW <= m; i += 3*SW)            
              MatKernelMultAB<4,3,ADD>(n,&a(0), a.Width(), &b(i/SW), b.Width(), &c(i), c.Width());
          t.Stop();
          cout << "MatKernel2AddAB 3x4, algined GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelAddAB aligned", 1e-9 * tot*its / t.GetTime()));
        }
      }




    if (what == 0 || what == 150)
      {
        // C=A*B
        Matrix<> a(4,n), b(4,n), c(3,4);
        a = 1; b = 2; c = 0;
        double tot = n*4*3;
        size_t its = 1e10 / tot + 1;
        SIMD<double,4> sum(0);
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              auto res = MatKernelScalAB<3,4>(n,&a(0), a.Width(), &b(0), b.Width());
              sum += get<0>(res) + get<1>(res) + get<2>(res);
            }
          t.Stop();
          cout << sum;
          cout << "MatKernelScalAB 4x3 = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelScalAB 4x3", 1e-9 * tot*its / t.GetTime()));
        }
      }
    

    if (what == 0 || what == 151)
      {
        // C=A*B
        Matrix<SIMD<double>> a(4,n), b(4,n);
        Matrix<> c(3,4);
        a = SIMD<double>(1); b = SIMD<double>(2); c = 0;
        double tot = n*4*3*SW;
        size_t its = 1e10 / tot + 1;
        SIMD<double,4> sum(0);
        {
          Timer t("C = A*B");
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              auto res = MatKernelScalAB<3,4>(n,&a(0), a.Width(), &b(0), b.Width());
              sum += get<0>(res) + get<1>(res) + get<2>(res);
            }
          t.Stop();
          cout << sum;
          cout << "MatKernelScalAB, simd 4x3 = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelScalAB, simd 4x3", 1e-9 * tot*its / t.GetTime()));
        }
      }
    


    if (what == 0 || what == 200)
      {
        // CalcInverse
        Matrix<> a(n,n);
        a = 1;
        a.Diag() = 10000;
        double tot = n*n*n;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("Inv(A)");
          t.Start();
          for (size_t j = 0; j < its; j++)
            CalcInverse(a, INVERSE_LIB::INV_NGBLA);
          t.Stop();
          cout << "Inv(A) GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Inv(A)", 1e-9 * tot *its / t.GetTime()));
        }
      }

    if (what == 0 || what == 201)
      {
        // CalcInverse
        Matrix<> a(n,n);
        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++)
            a(i,j) = cos(i+j);
        // a = 1;
        // a.Diag() = 1.1;
        double tot = n*n*n;
        size_t its = 1e9 / tot + 1;
        if (its > maxits) its = maxits;
        {
          Timer t("Inv(A)");
          t.Start();
          for (size_t j = 0; j < its; j++)
            CalcInverse(a, INVERSE_LIB::INV_NGBLA_LU);
          t.Stop();
          cout << "Inv(A) GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Inv(A)", 1e-9 * tot *its / t.GetTime()));
        }

        
        {
          Timer t("CalcLU");
          Array<int> p(n);
          Matrix<> ha = a;
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              ha = a;
              CalcLU(ha, p);
            }
          t.Stop();
          cout << "CalcLU GFlops = " << 1e-9 * tot/3*its / t.GetTime() << endl;
          timings.push_back(make_tuple("CalcLU", 1e-9 * tot/3 *its / t.GetTime()));
        }

        {
          Timer t("InvFromLU");
          Array<int> p(n);
          CalcLU(a, p);          
          Matrix<> ha = a;
          t.Start();
          for (size_t j = 0; j < its; j++)
            {
              ha = a;
              InverseFromLU(ha, p);
            }
          t.Stop();
          cout << "InvFromLU GFlops = " << 1e-9 * tot*2/3*its / t.GetTime() << endl;
          timings.push_back(make_tuple("InvFromLU", 1e-9 * tot*2/3 *its / t.GetTime()));
        }



      }

    

    if (what == 0 || what == 205)
      {
        // CalcInverse
        Matrix<double,ColMajor> a(n,n);
        a = 1;
        Trans(a).Diag() = 10000;
        double tot = n*n*n;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("Inv(A)");
          t.Start();
          for (size_t j = 0; j < its; j++)
            CalcLDL (SliceMatrix<double,ColMajor> (a));
          t.Stop();
          cout << "Inv(A) GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Inv(A)", 1e-9 * tot *its / t.GetTime()));
        }
      }

    
     if (what == 0 || what == 210)
      {
        // CalcInverse
        Matrix<> a(n,n);
        a = 1;
        a.Diag() = 10000;
        double tot = n*n*n;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("Inv(A)");
          t.Start();
          for (size_t j = 0; j < its; j++)
            LapackInverse(a);
          t.Stop();
          cout << "LapackInv(A) GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("LapackInv(A)", 1e-9 * tot *its / t.GetTime()));
        }
      }

     if (what == 0 || what == 211)
      {
        // CalcInverse
        Matrix<> a(n,n);
        a = 1;
        a.Diag() = 10000;
        double tot = n*n*n;
        size_t its = 1e9 / tot + 1;
        {
          Timer t("Inv(A)");
          t.Start();
          for (size_t j = 0; j < its; j++)
            LapackInverseSPD(a);
          t.Stop();
          cout << "LapackInv(A) GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("LapackInv(A)", 1e-9 * tot *its / t.GetTime()));
        }
      }
     
     if (what == 0 || what == 300)
      {
        // CalcSVD
        Matrix<> a(n,n);
        Matrix<double,ColMajor> U(n,n), V(n,n);
        for (int i = 0; i < a.Height(); i++)
          for (int j = 0; j < a.Width(); j++)
            a(i,j) = rand() / double(RAND_MAX);
        
        Matrix aorig = a;
        double tot = 5 * n*n*n;
        size_t its = 1e9 / tot + 1;
        if (!lapack)
          {
            Timer t("CalcSVD");
            t.Start();
            for (size_t j = 0; j < its; j++)
              {
                a = aorig;
                CalcSVD(a, U, V);
            }
            t.Stop();
            cout << "CalcSVD GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
            timings.push_back(make_tuple("CalcSVD", 1e-9 * tot *its / t.GetTime()));
          }
        else
          {
            Timer t("LapackSVD");
            t.Start();
            for (size_t j = 0; j < its; j++)
              {
                a = aorig;
                LapackSVD(a, U, V);
              }
            t.Stop();
            cout << "LapackSVD GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
            timings.push_back(make_tuple("LapackSVD", 1e-9 * tot *its / t.GetTime()));
          }
      }

    
    return timings;
  }


#if defined __AVX512F__

  double MatKernelMaskedScalAB (size_t n,
				double * pa, size_t da,
				double * pb, size_t db,
				const BitArray & ba)
  {
    SIMD<double,8> sum0 = 0.0;
    SIMD<double,8> sum1 = 0.0;
    int i(0), i0(0);
    auto bad = ba.Data();
    for ( ; i+16 <= n; i += 16, i0 += 2)
      {
	unsigned char mask0 = bad[i0];
	SIMD<mask64> m0 = GetMaskFromBits (unsigned(mask0));
	unsigned char mask1 = bad[i0+1];
	SIMD<mask64> m1 = GetMaskFromBits (unsigned(mask1));
	sum0 = If (m0, sum0+SIMD<double,8>(pa+i)*SIMD<double,8> (pb + i), sum0);
	sum1 = If (m1, sum1+SIMD<double,8>(pa+i+8)*SIMD<double,8>(pb + i + 8), sum1);
      } // n < i + 16
    if (i + 8 <= n)
      {
	unsigned char mask = bad[i0];
	SIMD<mask64> m0 = GetMaskFromBits (unsigned(mask));
	sum0 = If (m0, sum0+SIMD<double,8>(pa + i)*SIMD<double,8> (pb + i), sum0);
	i += 8;
      } // n < i + 8
    double hsum = HSum(sum0+sum1);
    for ( ; i < n; i++ )
      {
	if (ba.Test(i)) {
	  hsum += pa[i] * pb[i];
	}
      }
    return hsum;
  }

#elif defined __AVX__

  double MatKernelMaskedScalAB (size_t n,
				double * pa, size_t da,
				double * pb, size_t db,
				const BitArray & ba)
  {
    SIMD<double,4> sum0 = 0.0;
    SIMD<double,4> sum1 = 0.0;
    int i(0), i0(0);
    auto bad = ba.Data();
    for ( ; i+8 <= n; i += 8, i0++)
      {
	unsigned char mask = bad[i0];
	SIMD<mask64> m0 = GetMaskFromBits (unsigned(mask));
	SIMD<mask64> m1 = GetMaskFromBits (unsigned(mask) / 16);
	sum0 = If (m0, sum0+SIMD<double,4>(pa+i)*SIMD<double,4> (pb + i), sum0);
	sum1 = If (m1, sum1+SIMD<double,4>(pa+i+4)*SIMD<double,4>(pb + i + 4), sum1);
      } // n < i + 8
    if (i + 4 <= n)
      {
	unsigned char mask = bad[i0];
	SIMD<mask64> m0 = GetMaskFromBits (unsigned(mask));
	sum0 = If (m0, sum0+SIMD<double,4>(pa+i)*SIMD<double,4> (pb+i), sum0);
	i += 4;
      } // n < i + 4

    double ssum = HSum(sum0 + sum1);
    for ( ; i < n; i++ )
      if (ba.Test(i)) 
        ssum += pa[i] * pb[i];
    return ssum;
  }

#elif defined NETGEN_ARCH_AMD64

  double MatKernelMaskedScalAB (size_t n,
				double * pa, size_t da,
				double * pb, size_t db,
				const BitArray & ba)
  {
    SIMD<double,2> sum0 = 0.0;
    SIMD<double,2> sum1 = 0.0;
    SIMD<double,2> sum2 = 0.0;
    SIMD<double,2> sum3 = 0.0;
    int i(0), i0(0);
    auto bad = ba.Data();
    for ( ; i + 8 <= n; i += 8, i0++)
      {
	unsigned char mask = bad[i0];
	SIMD<mask64,2> m0 = GetMaskFromBits (unsigned(mask));
	SIMD<mask64,2> m1 = GetMaskFromBits (unsigned(mask) / 4);   // shift by 2
	SIMD<mask64,2> m2 = GetMaskFromBits (unsigned(mask) / 16);  // shift by 4
	SIMD<mask64,2> m3 = GetMaskFromBits (unsigned(mask) / 64);  // shift by 6
	sum0 = If (m0, sum0+SIMD<double,2>(pa+i)*SIMD<double,2> (pb+i), sum0);
	sum1 = If (m1, sum1+SIMD<double,2>(pa+i+2)*SIMD<double,2>(pb+i+2), sum1);
	sum2 = If (m2, sum2+SIMD<double,2>(pa+i+4)*SIMD<double,2>(pb+i+4), sum2);
	sum3 = If (m3, sum3+SIMD<double,2>(pa+i+6)*SIMD<double,2>(pb+i+6), sum3);
      } // n < i+8
    unsigned char mask = bad[i0];
    int shift = 1;
    if (i+4 <= n)
      {
	SIMD<mask64,2> m0 = GetMaskFromBits (unsigned(mask));
	SIMD<mask64,2> m1 = GetMaskFromBits (unsigned(mask) / 4); // shift by 2
	sum0 = If (m0, sum0+SIMD<double,2>(pa+i)*SIMD<double,2> (pb+i), sum0);
	sum1 = If (m1, sum1+SIMD<double,2>(pa+i+2)*SIMD<double,2>(pb+i+2), sum1);
	i += 4;
	shift = 16;
      } // n < i+4
    if (i+2 <= n)
      {
	SIMD<mask64,2> m0 = GetMaskFromBits (unsigned(mask) / shift);
	sum0 = If (m0, sum0+SIMD<double,2>(pa+i)*SIMD<double,2> (pb+i), sum0);
	i += 2;
      } // n < i+2
    sum0 += sum2;
    sum1 += sum3;
    double sum = HSum(sum0+sum1);
    for ( ; i < n; i++ )
      {
	if (ba.Test(i)) {
	  sum += pa[i]*pb[i];
	}
      }
    return sum;
  }

#else // ifdef AVX512/AVX/SSE

  double MatKernelMaskedScalAB (size_t n,
				double * pa, size_t da,
				double * pb, size_t db,
				const BitArray & ba)
  {
    double vhsum[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    int i(0);
    for ( ; i+8 <= n; i += 8)
      {
	for (int j = 0; j < 8; j++)
	  {
	    double hprod = pa[i+j]*pb[i+j];
	    if (ba.Test(i+j))
	      vhsum[j] += hprod;
	  }
      }
    for ( ; i < n; i++)
      if (ba.Test(i))
	vhsum[0] += pa[i]*pb[i];
    for (int j = 1; j < 8; j++)
      vhsum[0] += vhsum[j];
    return vhsum[0];
  }

#endif  // ifdef AVX512/AVX/SSE
  
}

