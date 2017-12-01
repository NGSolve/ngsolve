#include <bla.hpp>



namespace ngbla
{

#include "matkernel.hpp"


  /* ***************************** Copy Matrix *********************** */
  // copy matrix
  /*
    with AVX2 we get ~3 GF, independent of version 1 or version 2
    prefetch not yet tested 
   */
  inline void CopyMatrixIn (size_t h, size_t w,
                            double * ps, size_t dists,
                            SIMD<double> * pd, size_t distd)
  {
    constexpr int SW = SIMD<double>::Size();
    SIMD<mask64> mask(w % SW);
    // bool rest = (w % SW) != 0;
    auto pd0 = pd;
    for (size_t i = 0; i < h; i++, pd += distd, ps += dists)
      {
        // memcpy (pd, ps, sizeof(double)*w);

        size_t js = 0, jd=0;
        for ( ; js+SW <= w; js+=SW, jd++)
          {
            pd[jd] = SIMD<double>(ps+js);
            // Prefetch (psnext+j,  _MM_HINT_T1);
          }
        SIMD<double>(ps+js, mask).Store((double*) (pd+jd), mask);
         /*
        if (rest)
          pd[jd] = SIMD<double>(ps+js, mask);
         */
        /*
        size_t js = 0, jd=0;
        for ( ; js+4*SW <= w; js+=4*SW, jd+=4)
          {
            SIMD<double> v0(ps+js);
            SIMD<double> v1(ps+js+1*SW);
            SIMD<double> v2(ps+js+2*SW);
            SIMD<double> v3(ps+js+3*SW);
            pd[jd+0] = v0;
            pd[jd+1] = v1;
            pd[jd+2] = v2;
            pd[jd+3] = v3;
            // Prefetch (psnext+j,  _MM_HINT_T1);
          }
        for ( ; js+SW <= w; js+=SW, jd+=1)
          {
            SIMD<double> v0(ps+js);
            pd[jd+0] = v0;
            // Prefetch (psnext+j,  _MM_HINT_T1);
          }
        SIMD<double>(ps+js, mask).Store((double*) (pd+jd), mask);
        ps = psnext;
        */
      }
  }

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
    // constexpr size_t SWdTB = sizeof(TB)/sizeof(double);
    constexpr size_t SWdTB = sizeof(SIMD<double>)/sizeof(TB);
    size_t l = 0, lb = 0;
    for ( ; l+3*SW <= wb; l += 3*SW, lb += 3*SWdTB)
      MatKernelMultAB<H,3,OP> (hb, pa, da, pb+lb, db, pc+l, dc);
    for ( ; l+SW <= wb; l += SW, lb += SWdTB)
      MatKernelMultAB<H,1,OP> (hb, pa, da, pb+lb, db, pc+l, dc);
    if (l < wb)
      MatKernelMultABMask<H,OP>(hb, SIMD<mask64>(wb-l), pa, da, pb+lb, db, pc+l, dc);
  }
  
  void MultMatMat (size_t ha, size_t wa, size_t wb,
                   BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    // blockwise B, fits into L2 cache
    constexpr size_t BBH = 128;
    constexpr size_t BBW = 96;
    constexpr size_t SW = SIMD<double>::Size();
    alignas(64) SIMD<double> bb[BBH*BBW/SW];

    // c.AddSize(ha, wb) = 0.0;

    double *pb = &b(0);
    for (size_t i = 0; i < wa; i += BBH, pb += BBH*b.Dist())
      for (size_t j = 0; j < wb; j += BBW)
        {
          size_t hbi = min2(BBH, wa-i);
          size_t wbi = min2(BBW, wb-j);
          // cout << "copy B source = " << endl << SliceMatrix<> (hbi, wbi, b.Dist(), pb+j) << endl;
          CopyMatrixIn (hbi, wbi, pb+j, b.Dist(), &bb[0], BBW/SW);
          // cout << "copy B dest = " << endl << SliceMatrix<> (hbi, wbi, BBW, (double*)&bb[0]) << endl;          
          double * pa = &a(0)+i;
          double * pc = &c(0)+j;

          if (i == 0)
            {
              size_t k = 0;
              for ( ; k+4 <= ha; k += 4, pa += 4*a.Dist(), pc += 4 * c.Dist())
                // MatKernel2AddAB<4,SET> (hbi, wbi, pa, a.Dist(),  (double*)&bb[0], BBW, pc, c.Dist());
                MatKernel2AddAB<4,SET> (hbi, wbi, pa, a.Dist(),  &bb[0], BBW/SW, pc, c.Dist());
              switch (ha-k)
                {
                case 0: break;
                case 1: MatKernel2AddAB<1,SET> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                case 2: MatKernel2AddAB<2,SET> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                case 3: MatKernel2AddAB<3,SET> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                default: ; 
                }
              
            }
          else
            {
              size_t k = 0;
              for ( ; k+4 <= ha; k += 4, pa += 4*a.Dist(), pc += 4 * c.Dist())
                MatKernel2AddAB<4,ADD> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist());
              switch (ha-k)
                {
                case 0: break;
                case 1: MatKernel2AddAB<1,ADD> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                case 2: MatKernel2AddAB<2,ADD> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                case 3: MatKernel2AddAB<3,ADD> (hbi, wbi, pa, a.Dist(), (double*)&bb[0], BBW, pc, c.Dist()); break;
                default: ; 
                }
            }
        }
  }

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
  void MultMatMat (size_t ha, size_t wa, size_t wb,
                   BareSliceMatrix<> a, BareSliceMatrix<SIMD<double>> b, BareSliceMatrix<SIMD<double>> c)
  {
    size_t k = 0;
    constexpr size_t SW = SIMD<double>::Size();
    for ( ; k+3 <= wb; k += 3)
      MatKernel2MultAB<3>(ha, wa, a, b.Cols(k,k+3), c.Cols(k,k+3));
    for ( ; k+SW <= wb; k += SW)
      MatKernel2MultAB<1>(ha, wa, a, b.Cols(k,k+SW), c.Cols(k,k+SW));
  }




  /* ***************************** A * B^T *************************************** */

  template <typename FUNC>
  INLINE void TAddABt4 (size_t wa, size_t hc, size_t wc,
                        double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc,
                        FUNC func)
  {
    double * pb0 = pb;
    size_t i = 0;
    for ( ; i+3 <= hc; i += 3, pa += 3*da, pc += 3*dc)
      {
        double * pc1 = pc;
        double * pc2 = pc1 + dc;
        double * pc3 = pc2 + dc;
        double * pb = pb0;
        size_t j = 0;
        for ( ; j+4 <= wc; j += 4, pb += 4*db)
          {
            auto scal = MatKernelScalAB<3,4>(wa, pa, da, pb, db);
            auto s1 = func (SIMD<double,4>(pc1+j), get<0>(scal));
            auto s2 = func (SIMD<double,4>(pc2+j), get<1>(scal));
            auto s3 = func (SIMD<double,4>(pc3+j), get<2>(scal));
            s1.Store(pc1+j);
            s2.Store(pc2+j);
            s3.Store(pc3+j);
          }
        for ( ; j < wc; j++, pb += db)
          {
            auto scal = MatKernelScalAB<3,1>(wa, pa, da, pb, db);
            auto s1 = func (pc1[j], get<0>(scal));
            auto s2 = func (pc2[j], get<1>(scal));
            auto s3 = func (pc3[j], get<2>(scal));
            pc1[j] = s1;
            pc2[j] = s2;
            pc3[j] = s3;
          }
      }
    for ( ; i < hc; i ++, pa += da, pc += dc)
      {
        double * pc1 = pc;
        double * pb = pb0;
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
  

  template <typename FUNC>
  void TAddABt3 (size_t wa, size_t ha, size_t hb,
                 double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc,
                 FUNC func)
  {
    constexpr size_t bs = 32;  // height b
    for (size_t i = 0; i < hb; i += bs, pb += bs*db, pc += bs)
      TAddABt4 (wa, ha, min2(bs, hb-i),
                pa, da, pb, db, pc, dc, func);
  }

  template <typename FUNC>
  void TAddABt2 (size_t wa, size_t ha, size_t hb,
                 double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc,
                 // SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c,
                 FUNC func)
  {
    constexpr size_t bs = 96; // height a
    for (size_t i = 0; i < ha; i += bs, pa += bs*da, pc += bs*dc)
      TAddABt3(wa, min2(bs, ha-i), hb,
               pa, da, pb, db, pc, dc, func);
  }

  
  template <typename FUNC>
  void TAddABt1 (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c,
                FUNC func)
  {
    constexpr size_t bs = 256; // inner-product loop
    size_t wa = a.Width();
    double *pa = &a(0);
    double *pb = &b(0);
    double *pc = &c(0);
    for (size_t i = 0; i < wa; i += bs, pa+=bs, pb+=bs)
      TAddABt2 (min2(bs,wa-i), a.Height(), b.Height(),
                pa, a.Dist(), pb, b.Dist(), pc, c.Dist(), func);
  }
  
  void AddABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    // c += a * Trans(b);
    TAddABt1 (a, b, c, [] (auto c, auto ab) { return c+ab; });
  }

  void SubABt (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c)
  {
    // c -= a * Trans(b);
    TAddABt1 (a, b, c, [] (auto c, auto ab) { return c-ab; });
  }



  list<tuple<string,double>> Timing (int what, size_t n, size_t m, size_t k)
  {
    if (what < 0)
      {
        cout << "Available options timings are:\n"
          "-1 .. this help\n"
          "0 ... run all timings\n"
          "1 ... A = B,   A,B = n*m,   A = aligned, fixed dist\n"
          "2 ... A = 0,   A = n*m,     but sliced\n"
          "10 .. C = A * B,   A=n*m, B=m*k, C=n*k\n"
          "20 .. C = A * B    A=n*m, B=n*k', C=n*k', k'=round(k), B aligned\n"
          "100.. MultAddKernel  C += A * B,  C=4x12\n"
          "101.. MultAddKernel  C += A * B,  C=4x12\n,  B aligned"
             << endl;
        return list<tuple<string,double>>();
      }

    list<tuple<string,double>> timings;

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
        int its = 1e10 / tot + 1;
        {
          Timer t("Copy matrix, packed dest");
          t.Start();
          for (int j = 0; j < its; j++)
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
        int its = 1e10 / tot + 1;
        {
          Timer t("Zero matrix, packed dest");
          t.Start();
          for (int j = 0; j < its; j++)
            a.Rows(0,n).Cols(0,m) = j;
          t.Stop();
          cout << "Zero matrix GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("Zero matrix", 1e-9 * n*m*its / t.GetTime()));
        }
      }


    
    if (what == 0 || what == 10)
      {
        // C=A*B
        Matrix<> a(n,m), b(m,k), c(n,k);
        a = 1; b = 2;
        double tot = n*m*k;
        int its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (int j = 0; j < its; j++)
            MultMatMat(a,b,c);
          t.Stop();
          cout << "MultMatMat GFlops = " << 1e-9 * n*m*k*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatMat", 1e-9 * n*m*k*its / t.GetTime()));
        }
      }
    
    if (what == 0 || what == 100)
      {
        // C=A*B
        Matrix<> a(4,n), b(n,12), c(4,12);
        a = 1; b = 2; c = 0;
        double tot = n*4*12;
        int its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (int j = 0; j < its; j++)
            MatKernelMultAB<4,3,ADD>(n,&a(0), a.Width(), &b(0), b.Width(), &c(0), c.Width());
          t.Stop();
          cout << "Lapack GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelAddAB", 1e-9 * tot*its / t.GetTime()));
        }
      }

    if (what == 0 || what == 101)
      {
        // C=A*B
        Matrix<> a(4,n), c(4,12);
        Matrix<SIMD<double>> b(n, 3);
        a = 1; b = SIMD<double>(2); c = 0;
        double tot = n*4*12;
        int its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (int j = 0; j < its; j++)
            MatKernelMultAB<4,3,ADD>(n,&a(0), a.Width(), &b(0), b.Width(), &c(0), c.Width());
          t.Stop();
          cout << "Lapack GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MatKernelAddAB", 1e-9 * tot*its / t.GetTime()));
        }
      }

    return timings;
  }

  
}

