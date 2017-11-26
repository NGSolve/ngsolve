#include <bla.hpp>



namespace ngbla
{

#include "matkernel.hpp"

  
  INLINE void MatKernel2MultABMask(SIMD<mask64> mask, size_t ha, size_t wa, BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c)
  {
    size_t r = 0;
    size_t da = a.Dist();
    size_t dc = c.Dist();
    double * pa = &a(0,0);
    double * pc = &c(0,0);
    for ( ; r+4 <= ha; r += 4, pa += 4*da, pc += 4*dc)
      MatKernelMultABMask<4> (wa, mask, pa, da, &b(0,0), b.Dist(), pc, dc);
    switch (ha-r)
      {
      case 0: break;
      case 1:
        MatKernelMultABMask<1> (wa, mask, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      case 2:
        MatKernelMultABMask<2> (wa, mask, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      case 3:
        MatKernelMultABMask<3> (wa, mask, pa, da, &b(0,0), b.Dist(), pc, dc);
        break;
      default:
        ;
      }
    
  }
  
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
  void TAddABt3 (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c,
                 FUNC func)
  {
    constexpr size_t bs = 32;  // height b
    size_t hb = b.Height();
    size_t i = 0;
    double * pa = &a(0,0);
    double * pb = &b(0,0);
    for (size_t i = 0; i < hb; i += bs, pb += bs*b.Dist())
      {
        size_t wc = min2(bs, hb-i);
        TAddABt4 (a.Width(), c.Height(), wc, pa, a.Dist(), pb, b.Dist(), &c(i,0), c.Dist(), func);
      }
  }

  template <typename FUNC>
  void TAddABt2 (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c,
                 FUNC func)
  {
    constexpr size_t bs = 96;  // height a
    size_t ha = a.Height();
    size_t i = 0;
    for (size_t i = 0; i < ha; i += bs)
      {
        size_t i2 = i+bs;
        if (i2 > ha) i2 = ha;
        TAddABt3(a.Rows(i,i+i2), b, c.Rows(i,i2), func);
      }
  }

  
  template <typename FUNC>
  void TAddABt1 (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c,
                FUNC func)
  {
    constexpr size_t bs = 256; // inner-product loop
    size_t wa = a.Width();
    for (size_t i = 0; i < wa; i += bs)
      {
        size_t i2 = i + bs;
        if (i2 > wa) i2 = wa;
        TAddABt2 (a.Cols(i, i2), b.Cols(i,i2), c, func);
      }
  }
  
  void AddABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    // c += a * Trans(b);
    // return;
    TAddABt1 (a, b, c.AddSize(a.Height(),b.Height()), [] (auto c, auto ab) { return c+ab; });
  }

  void SubABt (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c)
  {
    // c -= a * Trans(b);
    // return;
    TAddABt1 (a, b, c, [] (auto c, auto ab) { return c-ab; });
  }

  
}

