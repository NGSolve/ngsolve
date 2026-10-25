#ifndef FILE_MATKERNEL32
#define FILE_MATKERNEL32

#if NETGEN_DEFAULT_SIMD_SIZE==1
#include "matkernel32_1.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==2
#include "matkernel32_2.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==4
#include "matkernel32_4.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==8
#include "matkernel32_8.hpp"
#endif


/*
  Template kernels for float and Complex32, replacing generated ones.
  Rules as in matkernel.hpp: Vec/Mat accumulators with compile-time
  Get/Set, forceinline, no C arrays.
  Unrolled loops are comma folds over index sequences instead of KFor
  lambdas: MSVC compiles every instantiated lambda operator() on its own,
  which tripled the compile time of ngblas32.cpp.
  k loops bump pointers, rows of A and C are indexed (as generated).
*/

#ifndef NGK_INLINE
#ifdef _MSC_VER
#define NGK_INLINE __forceinline
#define NGK_LAMBDA_INLINE [[msvc::forceinline]]
#elif defined(__clang__)
// nodebug: the inline chains otherwise triple the debug info
#define NGK_INLINE NETGEN_INLINE __attribute__((nodebug))
#define NGK_LAMBDA_INLINE NETGEN_LAMBDA_INLINE __attribute__((nodebug))
#else
#define NGK_INLINE NETGEN_INLINE
#define NGK_LAMBDA_INLINE NETGEN_LAMBDA_INLINE
#endif
#endif

// MSVC also compiles every inline function standalone (no /Zc:inline):
// the SmallK panel loops go out of line once instead of into each caller
#ifdef _MSC_VER
#define K32_PANEL_INLINE __declspec(noinline)
#else
#define K32_PANEL_INLINE NGK_INLINE
#endif

template <int N> using KSeq32 = std::make_integer_sequence<int,N>;

// wide SIMD<float,N> are composed of native ones, FMA and FMAComplex
// act per native part, as in the generated kernels
constexpr int K32_NATIVE = 2*SIMD<double>::Size();

template <int N>
NGK_INLINE SIMD<float,N> KFMA32 (float a, SIMD<float,N> b, SIMD<float,N> c)
{
  if constexpr (N <= K32_NATIVE) return FMA(a, b, c);
  else return { KFMA32(a, b.Lo(), c.Lo()), KFMA32(a, b.Hi(), c.Hi()) };
}

template <int N>
NGK_INLINE void KFMAC32 (SIMD<float,N> a, SIMD<float,N> b, SIMD<float,N> & c)
{
  if constexpr (N <= K32_NATIVE) FMAComplex(a, b, c);
  else
    {
      auto lo = c.Lo(), hi = c.Hi();
      KFMAC32 (a.Lo(), b.Lo(), lo);
      KFMAC32 (a.Hi(), b.Hi(), hi);
      c = { lo, hi };
    }
}

// N floats at distance d
template <int N>
NGK_INLINE SIMD<float,N> KGather32 (float * p, size_t d)
{
  if constexpr (N == 1) return SIMD<float,1>(*p);
  else return SIMD<float,N>(KGather32<N/2>(p, d), KGather32<N/2>(p+N/2*d, d));
}

// N Complex32 at distance d
template <int N>
NGK_INLINE SIMD<float,2*N> KGatherC32 (Complex32 * p, size_t d)
{
  if constexpr (N == 1) return SIMD<float,2>((float*)(void*)p);
  else return SIMD<float,2*N>(KGatherC32<N/2>(p, d), KGatherC32<N/2>(p+N/2*d, d));
}

// one Complex32 repeated N times
template <int N>
NGK_INLINE SIMD<float,2*N> KRepeat32 (SIMD<float,2> a)
{
  if constexpr (N == 1) return a;
  else
    {
      auto h = KRepeat32<N/2>(a);
      return SIMD<float,2*N>(h, h);
    }
}

// W entries of B's row, lane distance dl
template <int W, ORDERING OB>
NGK_INLINE SIMD<float,W> KLoadB32 (float * p, size_t dl)
{
  if constexpr (OB == RowMajor) return SIMD<float,W>(p);
  else return KGather32<W>(p, dl);
}

template <int W, ORDERING OB>
NGK_INLINE SIMD<float,2*W> KLoadB32 (Complex32 * p, size_t dl)
{
  if constexpr (OB == RowMajor) return SIMD<float,2*W>((float*)(void*)p);
  else return KGatherC32<W>(p, dl);
}

// C (+/-)= v, Complex32 as pairs of floats
template <bool ADD, bool POS, typename TS, typename T>
NGK_INLINE void KStore32 (TS v, T * pc)
{
  float * p = (float*)(void*)pc;
  if constexpr (ADD)
    v = POS ? TS(p) + v : TS(p) - v;
  else if constexpr (!POS)
    v = -v;
  v.Store(p);
}

template <bool ADD, bool POS, typename T>
NGK_INLINE void KFunc32 (T & c, T v)
{
  if constexpr (ADD)
    {
      if constexpr (POS) c += v;
      else c -= v;
    }
  else if constexpr (POS)
    c = v;
  else
    c = -v;
}

// row and k distances of A, k and lane distances of B
template <ORDERING OA, ORDERING OB, typename T>
NGK_INLINE auto KDists32 (BareSliceMatrix<T,OA> a, BareSliceMatrix<T,OB> b)
{
  size_t da = OA == RowMajor ? a.Dist() : 1;
  size_t dka = OA == RowMajor ? 1 : a.Dist();
  size_t dkb = OB == RowMajor ? b.Dist() : 1;
  size_t dl = OB == RowMajor ? 1 : b.Dist();
  return std::tuple(da, dka, dkb, dl);
}

// a * b + c for one entry of A
template <typename TW>
NGK_INLINE TW KUpdate32 (float * pa, TW b, TW c)
{
  return KFMA32 (*pa, b, c);
}

template <typename TW>
NGK_INLINE TW KUpdate32 (Complex32 * pa, TW b, TW c)
{
  KFMAC32 (KRepeat32<TW::Size()/2> (SIMD<float,2>((float*)(void*)pa)), b, c);
  return c;
}

// row-block helpers with only the template parameters they need: fewer instantiations

template <int H, typename TW, int... I>
NGK_INLINE void KZero32 (Vec<H,TW> & sum, std::integer_sequence<int,I...>)
{
  (sum.template Set<I> (TW(0.0f)), ...);
}

// sum(i) (+)= a(i) * b for H rows of A at distance da
template <bool FIRST, int H, typename T, typename TW, int... I>
NGK_INLINE void KCol32 (Vec<H,TW> & sum, TW b, T * pa, size_t da, std::integer_sequence<int,I...>)
{
  if constexpr (FIRST)
    (sum.template Set<I> (pa[I*da] * b), ...);
  else
    (sum.template Set<I> (KUpdate32 (pa+I*da, b, sum.template Get<I>())), ...);
}

// C (+/-)= sum for H rows of C
template <bool ADD, bool POS, int H, typename T, typename TW, int... I>
NGK_INLINE void KStoreRows32 (Vec<H,TW> & sum, T * pc, size_t dc, std::integer_sequence<int,I...>)
{
  (KStore32<ADD,POS> (sum.template Get<I>(), pc+I*dc), ...);
}



/* ****************** MatKernelFloat, MatKernelComplex32 ****************** */

// C (+/-)= A * B, outer product, one wide SIMD per row of C
template <typename T, size_t H, size_t W, bool ADD, bool POS, ORDERING OA, ORDERING OB>
NGK_INLINE void K32AB (size_t n, BareSliceMatrix<T,OA> a, BareSliceMatrix<T,OB> b, BareSliceMatrix<T,RowMajor> c)
{
  constexpr int L = std::is_same_v<T,Complex32> ? 2 : 1;   // floats per entry
  typedef SIMD<float,L*W> TW;
  auto [da, dka, dkb, dl] = KDists32 (a, b);
  T * pa = a.Data(), * pb = b.Data();
  Vec<H,TW> sum;
  KZero32 (sum, KSeq32<H>());
  for (size_t k = 0; k < n; k++, pa += dka, pb += dkb)
    KCol32<false> (sum, KLoadB32<W,OB> (pb, dl), pa, da, KSeq32<H>());
  KStoreRows32<ADD,POS> (sum, c.Data(), c.Dist(), KSeq32<H>());
}

// C (+/-)= A * B, A RowMajor, B ColMajor: inner products along k, entry IJ = (IJ/W, IJ%W)
template <size_t H, size_t W, bool ADD, bool POS>
struct KFloatABt
{
  static constexpr int VS = 2*SIMD<double>::Size();
  typedef SIMD<float,VS> TS;
  typedef Mat<H,W,TS> TSUM;
  typedef Mat<H,W,float> TVAL;

  template <int... IJ>
  static NGK_INLINE void Zero (TSUM & sum, std::integer_sequence<int,IJ...>)
  {
    (sum.template Set<IJ/W,IJ%W> (TS(0.0f)), ...);
  }

  template <int... I, int... J, int... IJ>
  static NGK_INLINE void Step (TSUM & sum, float * pa, size_t da, float * pb, size_t db,
                               std::integer_sequence<int,I...>, std::integer_sequence<int,J...>,
                               std::integer_sequence<int,IJ...>)
  {
    Vec<H,TS> va;
    Vec<W,TS> vb;
    (va.template Set<I> (TS(pa+I*da)), ...);
    (vb.template Set<J> (TS(pb+J*db)), ...);
    (sum.template Set<IJ/W,IJ%W> (FMA(va.template Get<IJ/W>(), vb.template Get<IJ%W>(), sum.template Get<IJ/W,IJ%W>())), ...);
  }

  // horizontal sums, scalar k tail, store
  template <int... IJ>
  static NGK_INLINE void Finish (TSUM & sum, size_t k, size_t n, float * pa, size_t da, float * pb, size_t db,
                                 float * pc, size_t dc, std::integer_sequence<int,IJ...>)
  {
    TVAL val;
    (val.template Set<IJ/W,IJ%W> (HSum(sum.template Get<IJ/W,IJ%W>())), ...);
    for ( ; k < n; k++, pa++, pb++)
      (val.template Set<IJ/W,IJ%W> (val.template Get<IJ/W,IJ%W>() + pa[(IJ/W)*da]*pb[(IJ%W)*db]), ...);
    (KFunc32<ADD,POS> (pc[(IJ/W)*dc+IJ%W], val.template Get<IJ/W,IJ%W>()), ...);
  }

  static NGK_INLINE void Run (size_t n, BareSliceMatrix<float,RowMajor> a, BareSliceMatrix<float,ColMajor> b, BareSliceMatrix<float,RowMajor> c)
  {
    float * pa = a.Data(), * pb = b.Data();
    size_t da = a.Dist(), db = b.Dist();
    TSUM sum;
    Zero (sum, KSeq32<H*W>());
    size_t k = 0;
    for ( ; k+VS <= n; k += VS, pa += VS, pb += VS)
      Step (sum, pa, da, pb, db, KSeq32<H>(), KSeq32<W>(), KSeq32<H*W>());
    Finish (sum, k, n, pa, da, pb, db, c.Data(), c.Dist(), KSeq32<H*W>());
  }
};

template <size_t H, size_t W, bool ADD, bool POS,
          ORDERING OA, ORDERING OB>
NGK_INLINE void MatKernelFloat
(size_t nk,
 BareSliceMatrix<float,OA> a,
 BareSliceMatrix<float,OB> b,
 BareSliceMatrix<float,RowMajor> c)
{
  if constexpr (OA == RowMajor && OB == ColMajor)
    KFloatABt<H,W,ADD,POS>::Run (nk, a, b, c);
  else
    K32AB<float,H,W,ADD,POS,OA,OB> (nk, a, b, c);
}

template <size_t H, size_t W, bool ADD, bool POS,
          ORDERING OA, ORDERING OB>
NGK_INLINE void MatKernelComplex32
(size_t nk,
 BareSliceMatrix<Complex32,OA> a,
 BareSliceMatrix<Complex32,OB> b,
 BareSliceMatrix<Complex32,RowMajor> c)
{
  K32AB<Complex32,H,W,ADD,POS,OA,OB> (nk, a, b, c);
}



/* ************ MatKernelFloatSmallK, MatKernelComplex32SmallK ************ */

// sum = A * B for H rows of A, k-major as the generated kernels;
// float: the first product initializes (K0 = 1), complex: zero init (K0 = 0)
template <int K0, int H, int AW, typename T, typename TW, int... K>
NGK_INLINE void KSmallKRows32 (Vec<H,TW> & sum, Vec<AW,TW> & b, T * pa, size_t da, size_t dk,
                               std::integer_sequence<int,K...>)
{
  if constexpr (K0 == 1)
    KCol32<true> (sum, b.template Get<0>(), pa, da, KSeq32<H>());
  else
    KZero32 (sum, KSeq32<H>());
  (KCol32<false> (sum, b.template Get<K0+K>(), pa+(K0+K)*dk, da, KSeq32<H>()), ...);
}

// AW rows of B, W entries each
template <int W, ORDERING OB, int AW, typename T, typename TW, int... K>
NGK_INLINE void KLoadPanel32 (Vec<AW,TW> & bp, T * pb, size_t dkb, size_t dl, std::integer_sequence<int,K...>)
{
  (bp.template Set<K> (KLoadB32<W,OB> (pb+K*dkb, dl)), ...);
}

// C (+/-)= A * B with fixed small k = AW, B panel kept in registers
template <typename T, size_t AW, bool ADD, bool POS, ORDERING OA, ORDERING OB>
struct K32SmallK
{
  static constexpr int L = std::is_same_v<T,Complex32> ? 2 : 1;
  static constexpr int K0 = (L == 1 && AW > 0) ? 1 : 0;

  // panels of W columns, HR rows at once
  template <int W, int HR>
  static K32_PANEL_INLINE void Panels (size_t ah, size_t bw, size_t & j,
                                       const BareSliceMatrix<T,OA> & a, const BareSliceMatrix<T,OB> & b,
                                       const BareSliceMatrix<T,RowMajor> & c)
  {
    typedef SIMD<float,L*W> TW;
    auto [da, dk, dkb, dl] = KDists32 (a, b);
    size_t dc = c.Dist();
    for ( ; j+W <= bw; j += W)
      {
        Vec<AW,TW> bp;
        KLoadPanel32<W,OB> (bp, b.Data()+j*dl, dkb, dl, KSeq32<AW>());
        // rows by index: bumped pointers cost clang and MSVC register copies here
        T * pa = a.Data(), * pc = c.Data()+j;
        size_t i = 0;
        for ( ; i+HR <= ah; i += HR)
          {
            Vec<HR,TW> sum;
            KSmallKRows32<K0> (sum, bp, pa+i*da, da, dk, KSeq32<AW-K0>());
            KStoreRows32<ADD,POS> (sum, pc+i*dc, dc, KSeq32<HR>());
          }
        for ( ; i < ah; i++)
          {
            Vec<1,TW> sum;
            KSmallKRows32<K0> (sum, bp, pa+i*da, da, dk, KSeq32<AW-K0>());
            KStoreRows32<ADD,POS> (sum, pc+i*dc, dc, KSeq32<1>());
          }
      }
  }

  // panel widths W0, W0/2, ..., 1, called only if a panel fits (calls are not free on MSVC)
  template <int W0, int HR, int... S>
  static NGK_INLINE void Widths (size_t ah, size_t bw, size_t & j,
                                 BareSliceMatrix<T,OA> a, BareSliceMatrix<T,OB> b, BareSliceMatrix<T,RowMajor> c,
                                 std::integer_sequence<int,S...>)
  {
    ((j+(W0>>S) <= bw ? Panels<(W0>>S),HR> (ah, bw, j, a, b, c) : void()), ...);
  }

  static NGK_INLINE void Run (size_t ah, size_t bw, BareSliceMatrix<T,OA> a, BareSliceMatrix<T,OB> b, BareSliceMatrix<T,RowMajor> c)
  {
    constexpr int WMAX = K32_NATIVE / L;   // entries per native vector
    size_t j = 0;
#if defined(__AVX512F__) || defined(__arm64__) || defined(__aarch64__)
    if constexpr (AW <= 2)
      if (4*WMAX <= bw)
        Panels<4*WMAX,4> (ah, bw, j, a, b, c);
    Widths<2*WMAX,4> (ah, bw, j, a, b, c, KSeq32<std::bit_width(unsigned(2*WMAX))>());
#else
    Widths<WMAX,2> (ah, bw, j, a, b, c, KSeq32<std::bit_width(unsigned(WMAX))>());
#endif
  }
};

template <size_t AW, bool ADD, bool POS,
          ORDERING OA, ORDERING OB>
inline void MatKernelFloatSmallK
(size_t ah, size_t bw,
 BareSliceMatrix<float,OA> a,
 BareSliceMatrix<float,OB> b,
 BareSliceMatrix<float,RowMajor> c)
{
  K32SmallK<float,AW,ADD,POS,OA,OB>::Run (ah, bw, a, b, c);
}

template <size_t AW, bool ADD, bool POS,
          ORDERING OA, ORDERING OB>
inline void MatKernelComplex32SmallK
(size_t ah, size_t bw,
 BareSliceMatrix<Complex32,OA> a,
 BareSliceMatrix<Complex32,OB> b,
 BareSliceMatrix<Complex32,RowMajor> c)
{
  K32SmallK<Complex32,AW,ADD,POS,OA,OB>::Run (ah, bw, a, b, c);
}

#endif
