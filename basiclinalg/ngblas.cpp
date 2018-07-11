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
    // auto pd0 = pd;
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

  /* ************************ Matrix * Vector ************************** */

  void MultMatVec (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    constexpr int SW = SIMD<double>::Size();
    size_t h = y.Size();
    size_t w = x.Size();
    size_t i = 0;
    SIMD<mask64> mask(w % SW);

    double * pa = &a(i,0);
    for ( ; i+8 <= h; i+=8, pa += 8*a.Dist())
      {
        auto scal = MatKernelScalAB<8,1> (w, pa, a.Dist(), &x(0), 0);
        SIMD<double,4> sum1(get<0>(scal), get<1>(scal), get<2>(scal), get<3>(scal));
        SIMD<double,4> sum2(get<4>(scal), get<5>(scal), get<6>(scal), get<7>(scal));
        sum1.Store(&y(i));        
        sum2.Store(&y(i+4));        
      }
    
    for ( ; i+4 <= h; i+=4, pa += 4*a.Dist())
      {
        /*
        SIMD<double> s0(0), s1(0), s2(0), s3(0);
        double * pa0 = &a(i,0);
        double * pa1 = pa0 + a.Dist();
        double * pa2 = pa0 + 2*a.Dist();
        double * pa3 = pa0 + 3*a.Dist();
        size_t j = 0;
        for ( ; j+SW <= w; j+=SW)
          {
            SIMD<double> xi(&x(j));
            s0 += xi * SIMD<double>(pa0+j);
            s1 += xi * SIMD<double>(pa1+j);
            s2 += xi * SIMD<double>(pa2+j);
            s3 += xi * SIMD<double>(pa3+j);
          }
        if (j < w)
          {
            SIMD<double> xi(&x(j), mask);
            s0 += xi * SIMD<double>(pa0+j, mask);
            s1 += xi * SIMD<double>(pa1+j, mask);
            s2 += xi * SIMD<double>(pa2+j, mask);
            s3 += xi * SIMD<double>(pa3+j, mask);
          }
        SIMD<double,4> sum = HSum(s0,s1,s2,s3);
        sum.Store(&y(i));
        */
        auto scal = MatKernelScalAB<4,1> (w, pa, a.Dist(), &x(0), 0);
        SIMD<double,4> sum(get<0>(scal), get<1>(scal), get<2>(scal), get<3>(scal));
        sum.Store(&y(i));        
      }

    for ( ; i+2 <= h; i+=2, pa += 2*a.Dist())
      {
        SIMD<double> s0(0), s1(0);
        double * pa0 = pa;
        double * pa1 = pa0 + a.Dist();
        size_t j = 0;
        for ( ; j+SW <= w; j+=SW)
          {
            SIMD<double> xi(&x(j));
            s0 += xi * SIMD<double>(pa0+j);
            s1 += xi * SIMD<double>(pa1+j);
          }
        if (j < w)
          {
            SIMD<double> xi(&x(j), mask);
            s0 += xi * SIMD<double>(pa0+j, mask);
            s1 += xi * SIMD<double>(pa1+j, mask);
          }
        // SIMD<double,2> sum = HSum(s0,s1);
        auto scal = HSum(s0,s1);
        SIMD<double,2> sum(get<0>(scal), get<1>(scal));
        sum.Store(&y(i));
      }

    for ( ; i+1 <= h; i++)
      {
        SIMD<double> s0(0), s1(0);
        double * pa0 = pa;
        size_t j = 0;
        for ( ; j+SW <= w; j+=SW)
          {
            SIMD<double> xi(&x(j));
            s0 += xi * SIMD<double>(pa0+j);
          }
        if (j < w)
          {
            SIMD<double> xi(&x(j), mask);
            s0 += xi * SIMD<double>(pa0+j, mask);
          }
        y(i) = HSum(s0);
      }

  }


  void MultMatTransVec (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
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




  /* ***************************** A * B^T *************************************** */

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
    SIMD<double> *pa = &a(0);
    SIMD<double> *pb = &b(0);
    double *pc = &c(0);
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
                &a(0), a.Dist(), &b(0), b.Dist(), &c(0), c.Dist(),
                [] (auto c, auto ab) { return c+ab; });
  }


  void AddABtSym (SliceMatrix<SIMD<double>> a,
                  SliceMatrix<SIMD<double>> b,
                  BareSliceMatrix<double> c)
  {
    TAddABt4Sym(a.Width(), a.Height(), b.Height(),
                &a(0), a.Width(), &b(0), b.Width(), &c(0), c.Dist(),
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
                               a.Height()*b.Height()*a.Width()*8);
    constexpr size_t bs = 64;
    for (size_t k = 0; k < a.Width(); k+=bs)
      {
        size_t k2 = min2(k+bs, a.Width());
        AddABt2<bs> (a.Cols(k,k2), b.Cols(k,k2), c);
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
  static constexpr size_t NB = 96;
  static constexpr size_t NK = 128;


  /*
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
            Transpose(a0,a1,a2,a3, b0,b1,b2,b3);
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
  */

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
    constexpr size_t SW = SIMD<double>::Size();
    alignas (64) double mema[NA*NK];
    SIMD<double> memb[3*NK];
    size_t na = a.Width();
    size_t nb = b.Width();
    size_t ha = a.Height();
    
    // loca = Trans(a);
    // for (size_t i = 0; i < loca.Width(); i++)
    // loca.Col(i) *= -diag(i);
    // c += loca * b;
    // return;

#ifdef __AVX512F__
    constexpr size_t HA = 6;
#else
    constexpr size_t HA = 4;
#endif

    size_t da = NA;
    // size_t db = b.Dist();
    double * pc = &c(0);

    SliceMatrix<> loca(a.Width(), a.Height(), NA, &mema[0]);
    MyTransposeScaleNeg (a, loca, diag);

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
                           &a(0,0), a.Dist(), (Complex*)&mema[0], CNA,
                           &diag(0), diag.Dist());

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


  

  /**************** timings *********************** */

  
  list<tuple<string,double>> Timing (int what, size_t n, size_t m, size_t k)
  {
    if (what < 0)
      {
        cout << "Available options timings are:\n"
          "-1 .. this help\n"
          "0 ... run all timings\n"
          "1 ... A = B,   A,B = n*m,   A = aligned, fixed dist\n"
          "2 ... A = 0,   A = n*m,     but sliced\n"
          "5 ... y = A*x,   A = n*m\n"
          "10 .. C = A * B,   A=n*m, B=m*k, C=n*k\n"
          // "20 .. C = A * B    A=n*m, B=n*k', C=n*k', k'=round(k), B aligned\n"
          "50 .. C += A * B^t,   A=n*k, B=m*k, C=n*m\n"
          "51 .. C += A * B^t,   A=n*k, B=m*k, C=n*m,  A,B aligned\n"
          "60 .. C -= A^t * D B,  A=n*k, B=n*m, C = k*m, D=diag\n"
          "100.. MultAddKernel  C += A * B,  A=4*n, B=n*3SW\n"
          "101.. MultAddKernel  C += A * B,  A=4*n, B=n*3SW, B aligned\n"
          "110.. MultAddKernel2  C += A * B,  A=4*n, B=n*m, m multiple of 3*SW\n"
          "111.. MultAddKernel2  C += A * B,  A=4*n, B=n*m, m multiple of 3*SW, B aligned\n"
          "150.. ScalKernel     C = A * B^t,  A=4*n, B = 3*n\n"
          "151.. ScalKernel     C = A * B^t,  A=4*n, B = 3*n\n, A,B aligned"
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

    if (what == 0 || what == 5)
      {
        // y = A*x
        Matrix<> a(n,m);
        Vector<> x(m), y(n);
        a = 1; x = 2;
        double tot = n*m;
        int its = 1e10 / tot + 1;
        {
          Timer t("y = A*x");
          t.Start();
          for (int j = 0; j < its; j++)
            MultMatVec(a,x,y);
          t.Stop();
          cout << "MultMatVec GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatVec", 1e-9 * n*m*its / t.GetTime()));
        }
        {
          Timer t("y = A*x, Lapack");
          t.Start();
          for (int j = 0; j < its; j++)
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
        int its = 1e10 / tot + 1;
        {
          Timer t("y = A*x");
          t.Start();
          for (int j = 0; j < its; j++)
            MultMatTransVec(a,x,y);
          t.Stop();
          cout << "MultMatTransVec GFlops = " << 1e-9 * n*m*its / t.GetTime() << endl;
          timings.push_back(make_tuple("MultMatVec", 1e-9 * n*m*its / t.GetTime()));
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

    if (what == 0 || what == 50)
      {
        // C=A*B^t
        Matrix<> a(n,k), b(m,k), c(n,m);
        a = 1; b = 2;
        c = 0.0;        
        double tot = n*m*k;
        int its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (int j = 0; j < its; j++)
            AddABt(a,b,c);
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
        int its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (int j = 0; j < its; j++)
            AddABt(SliceMatrix<double> (AFlatMatrix<double>(a)),
                   SliceMatrix<double> (AFlatMatrix<double>(b)),c);
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
        int its = 1e10 / tot + 1;
        {
          Timer t("C -= A^t*D*B");
          t.Start();
          for (int j = 0; j < its; j++)
            SubAtDB(a, d, b, c);
          t.Stop();
          cout << "AddAtDB GFlops = " << 1e-9 * tot*its / t.GetTime() << endl;
          timings.push_back(make_tuple("AddAtDB", 1e-9 * tot *its / t.GetTime()));
        }
      }

    
    if (what == 0 || what == 100)
      {
        // C=A*B
        Matrix<> a(4,n), b(n,3*SW), c(4,3*SW);
        a = 1; b = 2; c = 0;
        double tot = n*4*3*SW;
        int its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (int j = 0; j < its; j++)
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
        int its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (int j = 0; j < its; j++)
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
        int its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (int j = 0; j < its; j++)
            for (int i = 0; i+3*SW <= m; i += 3*SW)
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
        int its = 1e10 / tot + 1;
        {
          Timer t("C = A*B");
          t.Start();
          for (int j = 0; j < its; j++)
            for (int i = 0; i+3*SW <= m; i += 3*SW)            
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
        int its = 1e10 / tot + 1;
        SIMD<double,4> sum(0);
        {
          Timer t("C = A*B");
          t.Start();
          for (int j = 0; j < its; j++)
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
        int its = 1e10 / tot + 1;
        SIMD<double,4> sum(0);
        {
          Timer t("C = A*B");
          t.Start();
          for (int j = 0; j < its; j++)
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
    



    
    return timings;
  }

  
}

