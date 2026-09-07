#include <bla.hpp>

namespace ngbla
{

  namespace
  {
    template <bool ADD, bool POS, typename T>
    INLINE void Func (T & value, const T & sum)
    {
      if constexpr (ADD)
        {
          if constexpr (POS)
            value += sum;
          else
            value -= sum;
        }
      else
        {
          if constexpr (POS)
            value = sum;
          else
            value = -sum;
        }
    }

    INLINE void CopyMatrix32(size_t height, size_t width,
                             BareSliceMatrix<float,ColMajor> src,
                             BareSliceMatrix<float,RowMajor> dest)
    {
      constexpr size_t BS = 32;
      for (size_t ib = 0; ib < height; ib += BS)
        for (size_t jb = 0; jb < width; jb += BS)
          {
            size_t iend = min(height,ib+BS);
            size_t jend = min(width,jb+BS);
            size_t j = jb;
            for ( ; j+4 <= jend; j += 4)
              {
                size_t i = ib;
                for ( ; i+4 <= iend; i += 4)
                  {
                    SIMD<float,4> a0(src.Addr(i,j));
                    SIMD<float,4> a1(src.Addr(i,j+1));
                    SIMD<float,4> a2(src.Addr(i,j+2));
                    SIMD<float,4> a3(src.Addr(i,j+3));
                    auto [t0,t1] = Unpack(a0,a1);
                    auto [t2,t3] = Unpack(a2,a3);
                    SIMD<float,4>(t0.Lo(),t2.Lo()).Store(dest.Addr(i,j));
                    SIMD<float,4>(t1.Lo(),t3.Lo()).Store(dest.Addr(i+1,j));
                    SIMD<float,4>(t0.Hi(),t2.Hi()).Store(dest.Addr(i+2,j));
                    SIMD<float,4>(t1.Hi(),t3.Hi()).Store(dest.Addr(i+3,j));
                  }
                for ( ; i < iend; i++)
                  for (size_t jj = 0; jj < 4; jj++)
                    dest(i,j+jj) = src(i,j+jj);
              }
            for ( ; j < jend; j++)
              for (size_t i = ib; i < iend; i++)
                dest(i,j) = src(i,j);
          }
    }

  }

#include "matkernel32.hpp"

  /* *********************** GEMM float ******************************** */
  template <size_t H, size_t W, bool ADD, bool POS, ORDERING OA, ORDERING OB>
  INLINE void NgGEMMFloatWidths(size_t aw, size_t bw,
                                BareSliceMatrix<float,OA> a,
                                BareSliceMatrix<float,OB> b,
                                BareSliceMatrix<float,RowMajor> c,
                                size_t & j)
  {
    for ( ; j+W <= bw; j += W)
      MatKernelFloat<H,W,ADD,POS,OA,OB>(aw, a, b.Cols(j,j+W), c.Cols(j,j+W));

    if constexpr (W > 1)
      NgGEMMFloatWidths<H,W/2,ADD,POS,OA,OB>(aw, bw, a, b, c, j);
  }

  template <size_t H, bool ADD, bool POS, ORDERING OA, ORDERING OB>
  INLINE void NgGEMMFloatRows(size_t aw, size_t bw, size_t first_row,
                             BareSliceMatrix<float,OA> a,
                             BareSliceMatrix<float,OB> b,
                             BareSliceMatrix<float,RowMajor> c)
  {
    constexpr size_t MAX_NR = OA == RowMajor && OB == ColMajor ? NGBLAS32_FLOAT_INNER_MAX_NR : NGBLAS32_FLOAT_MAX_NR;
    size_t j = 0;
    NgGEMMFloatWidths<H,MAX_NR,ADD,POS,OA,OB> (aw, bw, a.Rows(first_row,first_row+H), b, c.Rows(first_row,first_row+H), j);
  }

  template <bool ADD, bool POS, ORDERING OA, ORDERING OB>
  void NgGEBP32 (size_t ah, size_t aw, size_t bw,
                 BareSliceMatrix<float,OA> a,
                 BareSliceMatrix<float,OB> b,
                 BareSliceMatrix<float,RowMajor> c)
  {
    constexpr size_t H = OA == RowMajor && OB == ColMajor ? NGBLAS32_FLOAT_INNER_MAX_MR : 4;
    size_t i = 0;
    for ( ; i+H <= ah; i += H)
      NgGEMMFloatRows<H,ADD,POS,OA,OB> (aw, bw, i, a, b, c);

    switch (ah-i)
      {
      case 3:
        NgGEMMFloatRows<3,ADD,POS,OA,OB> (aw, bw, i, a, b, c);
        break;
      case 2:
        NgGEMMFloatRows<2,ADD,POS,OA,OB> (aw, bw, i, a, b, c);
        break;
      case 1:
        NgGEMMFloatRows<1,ADD,POS,OA,OB> (aw, bw, i, a, b, c);
        break;
      }
  }

  template <bool ADD, bool POS, ORDERING OA, ORDERING OB>
  void NgGEMM32Intern (size_t ah, size_t aw, size_t bw,
                       BareSliceMatrix<float,OA> a,
                       BareSliceMatrix<float,OB> b,
                       BareSliceMatrix<float,RowMajor> c)
  {
    constexpr size_t BM = 72;
    constexpr size_t BK = 72;
    constexpr size_t BN = 512;

    if constexpr (OB == ColMajor)
      {
        Matrix<float> packed_b(BK,min(BN,bw));
        for (size_t j = 0; j < bw; j += BN)
          {
            size_t jmax = min(bw,j+BN);
            for (size_t k = 0; k < aw; k += BK)
              {
                size_t kmax = min(aw,k+BK);
                auto packed = packed_b.Rows(0,kmax-k).Cols(0,jmax-j);
                CopyMatrix32 (kmax-k,jmax-j,b.Rows(k,kmax).Cols(j,jmax),packed);
                if (k == 0)
                  NgGEMM32Intern<ADD,POS>(ah,kmax-k,jmax-j, AsBareSliceMatrix(a.Cols(k,kmax)), AsBareSliceMatrix(packed), AsBareSliceMatrix(c.Cols(j,jmax)));
                else
                  NgGEMM32Intern<true,POS>(ah,kmax-k,jmax-j, AsBareSliceMatrix(a.Cols(k,kmax)), AsBareSliceMatrix(packed), AsBareSliceMatrix(c.Cols(j,jmax)));
              }
          }
        return;
      }

    if constexpr (OA == ColMajor && OB == RowMajor)
      {
        Matrix<float> packed_a(BM,BK);
        for (size_t i = 0; i < ah; i += BM)
          {
            size_t imax = min(ah,i+BM);
            for (size_t k = 0; k < aw; k += BK)
              {
                size_t kmax = min(aw,k+BK);
                auto packed = packed_a.Rows(0,imax-i).Cols(0,kmax-k);
                CopyMatrix32 (imax-i,kmax-k,a.Rows(i,imax).Cols(k,kmax),packed);
                for (size_t j = 0; j < bw; j += BN)
                  {
                    size_t jmax = min(bw,j+BN);
                    if (k == 0)
                      NgGEBP32<ADD,POS,RowMajor,RowMajor>(imax-i,kmax-k,jmax-j,packed, b.Rows(k,kmax).Cols(j,jmax), c.Rows(i,imax).Cols(j,jmax));
                    else
                      NgGEBP32<true,POS,RowMajor,RowMajor>(imax-i,kmax-k,jmax-j,packed, b.Rows(k,kmax).Cols(j,jmax), c.Rows(i,imax).Cols(j,jmax));
                  }
              }
          }
        return;
      }

    for (size_t i = 0; i < ah; i += BM)
      {
        size_t imax = min(ah,i+BM);
        for (size_t j = 0; j < bw; j += BN)
          {
            size_t jmax = min(bw,j+BN);
            size_t kmax = min(aw,BK);
            NgGEBP32<ADD,POS,OA,OB>(imax-i,kmax,jmax-j,a.Rows(i,imax).Cols(0,kmax), b.Rows(0,kmax).Cols(j,jmax), c.Rows(i,imax).Cols(j,jmax));
            for (size_t k = BK; k < aw; k += BK)
              {
                size_t kmax = min(aw,k+BK);
                NgGEBP32<true,POS,OA,OB> (imax-i,kmax-k,jmax-j,a.Rows(i,imax).Cols(k,kmax), b.Rows(k,kmax).Cols(j,jmax), c.Rows(i,imax).Cols(j,jmax));
              }
          }
      }
  }

  /* *********************** GEMM float small k (A width) ******************************** */

  template <size_t AW, bool ADD, bool POS, ORDERING OA, ORDERING OB>
  void NgGEMMFuncSmall32 (size_t ah, size_t bw,
                          BareSliceMatrix<float,OA> a,
                          BareSliceMatrix<float,OB> b,
                          BareSliceMatrix<float,RowMajor> c)
  {
    MatKernelFloatSmallK<AW,ADD,POS,OA,OB> (ah,bw,a,b,c);
  }

  /* *********************** GEMM Complex32 small k (A width) ******************************** */

  template <size_t AW, bool ADD, bool POS, ORDERING OA, ORDERING OB>
  void NgGEMMFuncSmall32 (size_t ah, size_t bw,
                          BareSliceMatrix<Complex32,OA> a,
                          BareSliceMatrix<Complex32,OB> b,
                          BareSliceMatrix<Complex32,RowMajor> c)
  {
    if (ah == 4 && bw == 4)
      {
        MatKernelComplex32<4,4,ADD,POS,OA,OB> (AW,a,b,c);
        return;
      }
    MatKernelComplex32SmallK<AW,ADD,POS,OA,OB> (ah,bw,a,b,c);
  }


  /* *********************** GEMM Complex32 ******************************** */
  template <size_t H, size_t W, bool ADD, bool POS, ORDERING OA, ORDERING OB>
  INLINE void NgGEMMComplex32Widths(size_t aw, size_t bw,
                                    BareSliceMatrix<Complex32,OA> a,
                                    BareSliceMatrix<Complex32,OB> b,
                                    BareSliceMatrix<Complex32,RowMajor> c,
                                    size_t & j)
  {
    for ( ; j+W <= bw; j += W)
      MatKernelComplex32<H,W,ADD,POS,OA,OB> (aw, a, b.Cols(j,j+W), c.Cols(j,j+W));

    if constexpr (W > 1)
      NgGEMMComplex32Widths<H,W/2,ADD,POS,OA,OB> (aw, bw, a, b, c, j);
  }

  template <size_t H, bool ADD, bool POS, ORDERING OA, ORDERING OB>
  INLINE void NgGEMMComplex32Rows(size_t aw, size_t bw, size_t first_row,
                                  BareSliceMatrix<Complex32,OA> a,
                                  BareSliceMatrix<Complex32,OB> b,
                                  BareSliceMatrix<Complex32,RowMajor> c)
  {
    size_t j = 0;
    NgGEMMComplex32Widths<H,NGBLAS32_COMPLEX_MAX_NR,ADD,POS,OA,OB>
      (aw, bw, a.Rows(first_row,first_row+H), b, c.Rows(first_row,first_row+H), j);
  }

  template <bool ADD, bool POS, ORDERING OA, ORDERING OB>
  void NgGEBP32 (size_t ah, size_t aw, size_t bw,
                 BareSliceMatrix<Complex32,OA> a,
                 BareSliceMatrix<Complex32,OB> b,
                 BareSliceMatrix<Complex32,RowMajor> c)
  {
    size_t i = 0;
    for ( ; i+4 <= ah; i += 4)
      NgGEMMComplex32Rows<4,ADD,POS,OA,OB> (aw, bw, i, a, b, c);

    switch (ah-i)
      {
      case 3:
        NgGEMMComplex32Rows<3,ADD,POS,OA,OB> (aw, bw, i, a, b, c);
        break;
      case 2:
        NgGEMMComplex32Rows<2,ADD,POS,OA,OB> (aw, bw, i, a, b, c);
        break;
      case 1:
        NgGEMMComplex32Rows<1,ADD,POS,OA,OB> (aw, bw, i, a, b, c);
        break;
      }
  }

  template <bool ADD, bool POS, ORDERING OA, ORDERING OB>
  void NgGEMM32Intern (size_t ah, size_t aw, size_t bw,
                       BareSliceMatrix<Complex32,OA> a,
                       BareSliceMatrix<Complex32,OB> b,
                       BareSliceMatrix<Complex32,RowMajor> c)
  {
    constexpr size_t BM = 72;
    constexpr size_t BK = 72;
    constexpr size_t BN = 512;
    for (size_t i = 0; i < ah; i += BM)
      {
        size_t imax = min(ah,i+BM);
        for (size_t j = 0; j < bw; j += BN)
          {
            size_t jmax = min(bw,j+BN);
            size_t kmax = min(aw,BK);
            NgGEBP32<ADD,POS,OA,OB> (imax-i,kmax,jmax-j,a.Rows(i,imax).Cols(0,kmax), b.Rows(0,kmax).Cols(j,jmax),c.Rows(i,imax).Cols(j,jmax));
            for (size_t k = BK; k < aw; k += BK)
              {
                size_t kmax = min(aw,k+BK);
                NgGEBP32<true,POS,OA,OB> (imax-i,kmax-k,jmax-j,a.Rows(i,imax).Cols(k,kmax), b.Rows(k,kmax).Cols(j,jmax),c.Rows(i,imax).Cols(j,jmax));
              }
          }
      }
  }

  /* *********************** GEMM intern instantiations ******************************** */
  template <typename T, bool ADD, bool POS, ORDERING OA, ORDERING OB>
  void NgGEMM32Intern (size_t ah, size_t aw, size_t bw,
                      BareSliceMatrix<T,OA> a,
                      BareSliceMatrix<T,OB> b,
                      BareSliceMatrix<T,RowMajor> c)
  {
    NgGEMM32Intern<ADD,POS> (ah, aw, bw, a, b, c);
  }
  template NGS_DLL_HEADER void NgGEMM32Intern<float,false,false,RowMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,false,true,RowMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,true,false,RowMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,true,true,RowMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,false,false,RowMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,false,true,RowMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,true,false,RowMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,true,true,RowMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,false,false,ColMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,false,true,ColMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,true,false,ColMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,true,true,ColMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,false,false,ColMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,false,true,ColMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,true,false,ColMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<float,true,true,ColMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,ColMajor>, BareSliceMatrix<float,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,false,false,RowMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,false,true,RowMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,true,false,RowMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,true,true,RowMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,false,false,RowMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,false,true,RowMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,true,false,RowMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,true,true,RowMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,false,false,ColMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,false,true,ColMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,true,false,ColMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,true,true,ColMajor,RowMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,false,false,ColMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,false,true,ColMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,true,false,ColMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>);
  template NGS_DLL_HEADER void NgGEMM32Intern<Complex32,true,true,ColMajor,ColMajor>(size_t, size_t, size_t, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,ColMajor>, BareSliceMatrix<Complex32,RowMajor>);

  /* *********************** GEMM small k (A width) dispatch ******************************** */
  template <> NGS_DLL_HEADER pmatmat<float,RowMajor,RowMajor> dispatch_matmat<float,false,false,RowMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,RowMajor,RowMajor> dispatch_matmat<float,false,true,RowMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,RowMajor,RowMajor> dispatch_matmat<float,true,false,RowMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,RowMajor,RowMajor> dispatch_matmat<float,true,true,RowMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,RowMajor,ColMajor> dispatch_matmat<float,false,false,RowMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,RowMajor,ColMajor> dispatch_matmat<float,false,true,RowMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,RowMajor,ColMajor> dispatch_matmat<float,true,false,RowMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,RowMajor,ColMajor> dispatch_matmat<float,true,true,RowMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,ColMajor,RowMajor> dispatch_matmat<float,false,false,ColMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,ColMajor,RowMajor> dispatch_matmat<float,false,true,ColMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,ColMajor,RowMajor> dispatch_matmat<float,true,false,ColMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,ColMajor,RowMajor> dispatch_matmat<float,true,true,ColMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,ColMajor,ColMajor> dispatch_matmat<float,false,false,ColMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,ColMajor,ColMajor> dispatch_matmat<float,false,true,ColMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,ColMajor,ColMajor> dispatch_matmat<float,true,false,ColMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<float,ColMajor,ColMajor> dispatch_matmat<float,true,true,ColMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,RowMajor,RowMajor> dispatch_matmat<Complex32,false,false,RowMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,RowMajor,RowMajor> dispatch_matmat<Complex32,false,true,RowMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,RowMajor,RowMajor> dispatch_matmat<Complex32,true,false,RowMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,RowMajor,RowMajor> dispatch_matmat<Complex32,true,true,RowMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,RowMajor,ColMajor> dispatch_matmat<Complex32,false,false,RowMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,RowMajor,ColMajor> dispatch_matmat<Complex32,false,true,RowMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,RowMajor,ColMajor> dispatch_matmat<Complex32,true,false,RowMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,RowMajor,ColMajor> dispatch_matmat<Complex32,true,true,RowMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,ColMajor,RowMajor> dispatch_matmat<Complex32,false,false,ColMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,ColMajor,RowMajor> dispatch_matmat<Complex32,false,true,ColMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,ColMajor,RowMajor> dispatch_matmat<Complex32,true,false,ColMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,ColMajor,RowMajor> dispatch_matmat<Complex32,true,true,ColMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,ColMajor,ColMajor> dispatch_matmat<Complex32,false,false,ColMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,ColMajor,ColMajor> dispatch_matmat<Complex32,false,true,ColMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,ColMajor,ColMajor> dispatch_matmat<Complex32,true,false,ColMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmat<Complex32,ColMajor,ColMajor> dispatch_matmat<Complex32,true,true,ColMajor,ColMajor>[9];


  template <typename T, bool ADD, bool POS, ORDERING OA, ORDERING OB>
  void InitMatMat32 ()
  {
    Iterate<std::size(dispatch_matmat<T,ADD,POS,OA,OB>)> ([&] (auto i)
    {
      dispatch_matmat<T,ADD,POS,OA,OB>[i] = static_cast<pmatmat<T,OA,OB>> (&NgGEMMFuncSmall32<i,ADD,POS,OA,OB>);
    });
  }

  template <typename T, ORDERING OA, ORDERING OB>
  void InitMatMat32 ()
  {
    InitMatMat32<T,false,false,OA,OB>();
    InitMatMat32<T,false,true,OA,OB>();
    InitMatMat32<T,true,false,OA,OB>();
    InitMatMat32<T,true,true,OA,OB>();
  }

  template <typename T>
  void InitMatMat32 ()
  {
    InitMatMat32<T,RowMajor,RowMajor>();
    InitMatMat32<T,RowMajor,ColMajor>();
    InitMatMat32<T,ColMajor,RowMajor>();
    InitMatMat32<T,ColMajor,ColMajor>();
  }

  auto init_matmat32 = [] ()
  {
    InitMatMat32<float>();
    InitMatMat32<Complex32>();
    return 1;
  }();
}
