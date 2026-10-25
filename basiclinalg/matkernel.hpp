#if NETGEN_DEFAULT_SIMD_SIZE==1
#include "matkernel_1.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==2
#include "matkernel_2.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==4
#include "matkernel_4.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==8
#include "matkernel_8.hpp"
#endif


/*
  Template kernels, replacing generated ones.
  Accumulators and per-step vectors live in Mat<H,W,SIMD> accessed by
  compile-time Get/Set, so compilers keep them in registers.
  Rules found with MSVC 17.x/18.x, gcc and clang:
  - no local C arrays: MSVC then copies registers every iteration
  - no capturing lambdas: in big callers MSVC leaves their closure
    objects out of line; bodies are stateless, state is passed by KFor
  - everything __forceinline on MSVC, plain inline is not reliable
*/

#ifdef _MSC_VER
#define NGK_INLINE __forceinline
#define NGK_LAMBDA_INLINE [[msvc::forceinline]]
#else
#define NGK_INLINE NETGEN_INLINE
#define NGK_LAMBDA_INLINE NETGEN_LAMBDA_INLINE
#endif

// f(IC<I>(), args...) for I = I0 ... N-1
template <int I, int N, typename F, typename... ARGS>
NGK_INLINE void KFor (F f, ARGS &... args)
{
  if constexpr (I < N)
    {
      f(IC<I>(), args...);
      KFor<I+1,N> (f, args...);
    }
}

template <typename T>
NGK_INLINE SIMD<double> KLoad (T * p)
{
  if constexpr (std::is_same_v<T,double>) return SIMD<double>(p);
  else return *p;
}



/* ********************* MatKernelMultAB ********************** */

#ifndef NGCORE_HAVE_WIDE_FMA
#error "netgen too old: KMultAB needs the fused FMA for wide SIMD<double,N> (ngcore simd_generic.hpp)"
#endif

// C (op)= A * B, one wide SIMD<double,W*SW> per row of C
template <size_t H, size_t W, OPERATION OP, typename TB>
struct KMultAB
{
  static constexpr int SW = SIMD<double>::Size();
  typedef SIMD<double,W*SW> TW;
  typedef Vec<H,TW> TSUM;

  static NGK_INLINE void Step (TSUM & sum, double * pa, size_t da, TB * pb)
  {
    TW b((double*)pb);
    KFor<0,H> ([] (auto i, TSUM & sum, TW & b, double *& pa, size_t & da) NGK_LAMBDA_INLINE
      {
        TW a(pa[i*da]);
        if constexpr (OP == SET || OP == ADD)
          sum.template Set<i> (FMA(a, b, sum.template Get<i>()));
        else
          sum.template Set<i> (FNMA(a, b, sum.template Get<i>()));
      }, sum, b, pa, da);
  }

  static NGK_INLINE void Run (size_t n, double * pa, size_t da, TB * pb, size_t db, double * pc, size_t dc)
  {
    // C rows by a bumped pointer: index arithmetic costs MSVC registers
    TSUM sum;
    double * pci = pc;
    KFor<0,H> ([] (auto i, TSUM & sum, double *& pci, size_t & dc) NGK_LAMBDA_INLINE
      {
        if constexpr (OP == SET || OP == SETNEG)
          sum.template Set<i> (TW(0.0));
        else
          sum.template Set<i> (TW(pci));
        pci += dc;
      }, sum, pci, dc);

    for (size_t k = 0; k < n; k++, pa++, pb += db)
      Step (sum, pa, da, pb);

    KFor<0,H> ([] (auto i, TSUM & sum, double *& pc, size_t & dc) NGK_LAMBDA_INLINE
      { sum.template Get<i>().Store(pc); pc += dc; }, sum, pc, dc);
  }
};

template <size_t H, size_t W, OPERATION OP, typename TB>
NGK_INLINE void MatKernelMultAB_impl (size_t n, double * pa, size_t da, TB * pb, size_t db, double * pc, size_t dc)
{
  KMultAB<H,W,OP,TB>::Run (n, pa, da, pb, db, pc, dc);
}

template <size_t H, size_t W, OPERATION OP>
NGK_INLINE void MatKernelMultAB
(size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  MatKernelMultAB_impl<H,W,OP> (n, pa, da, pb, db, pc, dc);
}

template <size_t H, size_t W, OPERATION OP>
NGK_INLINE void MatKernelMultAB
(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, double * pc, size_t dc)
{
  MatKernelMultAB_impl<H,W,OP> (n, pa, da, pb, db, pc, dc);
}



/* ****************** MatKernelAlignedMultAB ******************* */

// C = A * B, B and C as SIMD rows
template <size_t H, size_t W>
struct KAlignedMultAB
{
  typedef Mat<H,W,SIMD<double>> TSUM;
  typedef Mat<1,W,SIMD<double>> TB_ROW;

  static NGK_INLINE void Step (TSUM & sum, double * pa, size_t da, SIMD<double> * pb)
  {
    TB_ROW b;
    KFor<0,W> ([] (auto j, TB_ROW & b, SIMD<double> *& pb) NGK_LAMBDA_INLINE
      { b.template Set<0,j> (pb[j]); }, b, pb);
    KFor<0,H> ([] (auto i, TSUM & sum, TB_ROW & b, double *& pa, size_t & da) NGK_LAMBDA_INLINE
      {
        SIMD<double> a(pa[i*da]);
        KFor<0,W> ([] (auto j, auto i, TSUM & sum, TB_ROW & b, SIMD<double> & a) NGK_LAMBDA_INLINE
          { sum.template Set<i,j> (FMA(a, b.template Get<0,j>(), sum.template Get<i,j>())); }, i, sum, b, a);
      }, sum, b, pa, da);
  }

  static NGK_INLINE void Run (size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, SIMD<double> * pc, size_t dc)
  {
    TSUM sum;
    KFor<0,H> ([] (auto i, TSUM & sum) NGK_LAMBDA_INLINE
      {
        KFor<0,W> ([] (auto j, auto i, TSUM & sum) NGK_LAMBDA_INLINE
          { sum.template Set<i,j> (SIMD<double>(0.0)); }, i, sum);
      }, sum);

    for (size_t k = 0; k < n; k++, pa++, pb += db)
      Step (sum, pa, da, pb);

    KFor<0,H> ([] (auto i, TSUM & sum, SIMD<double> *& pc, size_t & dc) NGK_LAMBDA_INLINE
      {
        KFor<0,W> ([] (auto j, auto i, TSUM & sum, SIMD<double> *& pc) NGK_LAMBDA_INLINE
          { pc[j] = sum.template Get<i,j>(); }, i, sum, pc);
        pc += dc;
      }, sum, pc, dc);
  }
};

template <size_t H, size_t W>
NGK_INLINE void MatKernelAlignedMultAB_impl
(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, SIMD<double> * pc, size_t dc)
{
  KAlignedMultAB<H,W>::Run (n, pa, da, pb, db, pc, dc);
}

template <size_t H, size_t W>
NGK_INLINE void MatKernelAlignedMultAB
(size_t n, double * pa, size_t da, SIMD<double> * pb, size_t db, SIMD<double> * pc, size_t dc)
{
  MatKernelAlignedMultAB_impl<H,W> (n, pa, da, pb, db, pc, dc);
}



/* ******************** MatKernelMultABMask ******************** */

// masked last column block: C(:,0:mask) (op)= A * B
template <size_t H, OPERATION OP, typename TB>
struct KMultABMask
{
  typedef Vec<H,SIMD<double>> TSUM;

  static NGK_INLINE void Step (TSUM & sum, SIMD<mask64> mask, double * pa, size_t da, TB * pb)
  {
    SIMD<double> b((double*)pb, mask);
    KFor<0,H> ([] (auto i, TSUM & sum, SIMD<double> & b, double *& pa, size_t & da) NGK_LAMBDA_INLINE
      {
        SIMD<double> a(pa[i*da]);
        if constexpr (OP == SET || OP == ADD)
          sum.template Set<i> (FMA(a, b, sum.template Get<i>()));
        else
          sum.template Set<i> (FNMA(a, b, sum.template Get<i>()));
      }, sum, b, pa, da);
  }

  static NGK_INLINE void Run (size_t n, SIMD<mask64> mask, double * pa, size_t da, TB * pb, size_t db, double * pc, size_t dc)
  {
    TSUM sum;
    double * pci = pc;
    KFor<0,H> ([] (auto i, TSUM & sum, SIMD<mask64> & mask, double *& pci, size_t & dc) NGK_LAMBDA_INLINE
      {
        if constexpr (OP == SET || OP == SETNEG)
          sum.template Set<i> (SIMD<double>(0.0));
        else
          sum.template Set<i> (SIMD<double>(pci, mask));
        pci += dc;
      }, sum, mask, pci, dc);

    for (size_t k = 0; k < n; k++, pa++, pb += db)
      Step (sum, mask, pa, da, pb);

    KFor<0,H> ([] (auto i, TSUM & sum, SIMD<mask64> & mask, double *& pc, size_t & dc) NGK_LAMBDA_INLINE
      { sum.template Get<i>().Store(pc, mask); pc += dc; }, sum, mask, pc, dc);
  }
};

template <size_t H, OPERATION OP, typename TB>
NGK_INLINE void MatKernelMultABMask_impl (size_t n, SIMD<mask64> mask, double * pa, size_t da, TB * pb, size_t db, double * pc, size_t dc)
{
  KMultABMask<H,OP,TB>::Run (n, mask, pa, da, pb, db, pc, dc);
}

template <size_t H, OPERATION OP>
NGK_INLINE void MatKernelMultABMask
(size_t n, SIMD<mask64> mask, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  MatKernelMultABMask_impl<H,OP> (n, mask, pa, da, pb, db, pc, dc);
}

template <size_t H, OPERATION OP>
NGK_INLINE void MatKernelMultABMask
(size_t n, SIMD<mask64> mask, double * pa, size_t da, SIMD<double> * pb, size_t db, double * pc, size_t dc)
{
  MatKernelMultABMask_impl<H,OP> (n, mask, pa, da, pb, db, pc, dc);
}



/* ********************** MatKernelScalAB ********************** */

// inner products of rows of A with rows of B, H x W results
// ALLFMA=false keeps mul+add for H*W < 12 (shorter dependency chain)
template <size_t H, size_t W, bool ALLFMA, typename T>
struct KScalAB
{
  typedef Mat<H,W,SIMD<double>> TSUM;
  static constexpr bool simded = !std::is_same_v<T,double>;
  static constexpr size_t SW = simded ? 1 : SIMD<double>::Size();

  // a is loaded where it is used: compilers fold it into the multiply
  static NGK_INLINE void Step (TSUM & sum, T * pa, size_t da, T * pb, size_t db)
  {
    KFor<0,W> ([] (auto j, TSUM & sum, T *& pa, size_t & da, T *& pb, size_t & db) NGK_LAMBDA_INLINE
      {
        SIMD<double> b = KLoad(pb+j*db);
        KFor<0,H> ([] (auto i, auto j, TSUM & sum, T *& pa, size_t & da, SIMD<double> & b) NGK_LAMBDA_INLINE
          {
            if constexpr (ALLFMA || H*W >= 12)
              sum.template Set<i,j> (FMA(KLoad(pa+i*da), b, sum.template Get<i,j>()));
            else
              sum.template Set<i,j> (sum.template Get<i,j>() + KLoad(pa+i*da)*b);
          }, j, sum, pa, da, b);
      }, sum, pa, da, pb, db);
  }

  static NGK_INLINE void Rest (TSUM & sum, size_t r, double * pa, size_t da, double * pb, size_t db)
  {
    SIMD<mask64> mask(r);
    Vec<H,SIMD<double>> a;   // masked loads cannot be folded into the fma
    KFor<0,H> ([] (auto i, Vec<H,SIMD<double>> & a, double *& pa, size_t & da, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      { a.template Set<i> (SIMD<double>(pa+i*da, mask)); }, a, pa, da, mask);
    KFor<0,W> ([] (auto j, TSUM & sum, Vec<H,SIMD<double>> & a, double *& pb, size_t & db, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> b(pb+j*db, mask);
        KFor<0,H> ([] (auto i, auto j, TSUM & sum, Vec<H,SIMD<double>> & a, SIMD<double> & b) NGK_LAMBDA_INLINE
          { sum.template Set<i,j> (FMA(a.template Get<i>(), b, sum.template Get<i,j>())); }, j, sum, a, b);
      }, sum, a, pb, db, mask);
  }

  template <size_t I>
  static NGK_INLINE auto RowSum (TSUM & sum)
  {
    return [] <size_t... J> (std::index_sequence<J...>, TSUM & sum) NGK_LAMBDA_INLINE
      { return HSum(sum.template Get<I,J>()...); } (std::make_index_sequence<W>(), sum);
  }

  template <size_t... I>
  static NGK_INLINE auto Result (TSUM & sum, std::index_sequence<I...>)
  {
    if constexpr (W == 1 && H % 4 == 0)
      return Concat(make_tuple(HSum(sum.template Get<4*I,0>(), sum.template Get<4*I+1,0>(),
                                    sum.template Get<4*I+2,0>(), sum.template Get<4*I+3,0>())...));
    else
      return make_tuple(RowSum<I>(sum)...);
  }

  static NGK_INLINE auto Run (size_t n, T * pa, size_t da, T * pb, size_t db)
  {
    TSUM sum;
    KFor<0,H> ([] (auto i, TSUM & sum) NGK_LAMBDA_INLINE
      { KFor<0,W> ([] (auto j, auto i, TSUM & sum) NGK_LAMBDA_INLINE
        { sum.template Set<i,j> (SIMD<double>(0.0)); }, i, sum); }, sum);

    for (size_t k = 0; k < n/SW; k++, pa += SW, pb += SW)
      Step (sum, pa, da, pb, db);
    if constexpr (!simded)
      if (size_t r = n % SW)
        Rest (sum, r, pa, da, pb, db);

    if constexpr (W == 1 && H % 4 == 0)
      return Result (sum, std::make_index_sequence<H/4>());
    else
      return Result (sum, std::make_index_sequence<H>());
  }
};

template <size_t H, size_t W, bool ALLFMA = false, typename T>
NGK_INLINE auto MatKernelScalAB_impl (size_t n, T * pa, size_t da, T * pb, size_t db)
{
  return KScalAB<H,W,ALLFMA,T>::Run (n, pa, da, pb, db);
}

template <size_t H, size_t W>
NGK_INLINE auto MatKernelScalAB (size_t n, double * pa, size_t da, double * pb, size_t db)
{
  return MatKernelScalAB_impl<H,W> (n, pa, da, pb, db);
}

template <size_t H, size_t W>
NGK_INLINE auto MatKernelScalAB (size_t n, SIMD<double> * pa, size_t da, SIMD<double> * pb, size_t db)
{
  return MatKernelScalAB_impl<H,W> (n, pa, da, pb, db);
}



/* ********************* MultiVector kernels ******************** */
/* rows given by pointer arrays, row pointers are loaded once      */

template <size_t N, typename T>
using KPtrs = Vec<N,T*>;

template <size_t N, typename T>
NGK_INLINE KPtrs<N,T> KLoadPtrs (T ** pp)
{
  KPtrs<N,T> p;
  KFor<0,N> ([] (auto i, KPtrs<N,T> & p, T **& pp) NGK_LAMBDA_INLINE
    { p.template Set<i> (pp[i]); }, p, pp);
  return p;
}

template <size_t H, size_t W, typename TM>
NGK_INLINE void KZero (TM & sum)
{
  KFor<0,H> ([] (auto i, TM & sum) NGK_LAMBDA_INLINE
    { KFor<0,W> ([] (auto j, auto i, TM & sum) NGK_LAMBDA_INLINE
      { sum.template Set<i,j> (SIMD<double>(0.0)); }, i, sum); }, sum);
}

// C = A * B^t, inner products of rows
template <size_t H, size_t W>
struct KMultiVecScalAB
{
  typedef Mat<H,W,SIMD<double>> TSUM;
  static constexpr size_t SW = SIMD<double>::Size();

  template <bool MASK>
  static NGK_INLINE void Step (TSUM & sum, KPtrs<H,double> & pa, KPtrs<W,double> & pb, size_t i, SIMD<mask64> mask)
  {
    KFor<0,W> ([] (auto j, TSUM & sum, KPtrs<H,double> & pa, KPtrs<W,double> & pb, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> b = MASK ? SIMD<double>(pb.template Get<j>()+i, mask) : SIMD<double>(pb.template Get<j>()+i);
        KFor<0,H> ([] (auto k, auto j, TSUM & sum, KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask, SIMD<double> & b) NGK_LAMBDA_INLINE
          {
            SIMD<double> a = MASK ? SIMD<double>(pa.template Get<k>()+i, mask) : SIMD<double>(pa.template Get<k>()+i);
            if constexpr (MASK || H*W >= 12)
              sum.template Set<k,j> (FMA(a, b, sum.template Get<k,j>()));
            else
              sum.template Set<k,j> (sum.template Get<k,j>() + a*b);
          }, j, sum, pa, i, mask, b);
      }, sum, pa, pb, i, mask);
  }

  template <size_t I>
  static NGK_INLINE auto RowSum (TSUM & sum)
  {
    return [] <size_t... J> (std::index_sequence<J...>, TSUM & sum) NGK_LAMBDA_INLINE
      { return HSum(sum.template Get<I,J>()...); } (std::make_index_sequence<W>(), sum);
  }

  template <size_t... I>
  static NGK_INLINE auto Result (TSUM & sum, std::index_sequence<I...>)
  {
    if constexpr (W == 1 && H % 4 == 0)
      return make_tuple(HSum(sum.template Get<4*I,0>(), sum.template Get<4*I+1,0>(),
                             sum.template Get<4*I+2,0>(), sum.template Get<4*I+3,0>())...);
    else
      return make_tuple(RowSum<I>(sum)...);
  }

  static NGK_INLINE auto Run (size_t n, double ** ppa, double ** ppb)
  {
    KPtrs<H,double> pa = KLoadPtrs<H>(ppa);
    KPtrs<W,double> pb = KLoadPtrs<W>(ppb);
    TSUM sum;
    KZero<H,W> (sum);
    SIMD<mask64> mask(SW);
    size_t i = 0;
    for ( ; i+SW <= n; i += SW)
      Step<false> (sum, pa, pb, i, mask);
    if (size_t r = n % SW)
      Step<true> (sum, pa, pb, i, SIMD<mask64>(r));

    if constexpr (W == 1 && H % 4 == 0)
      return Result (sum, std::make_index_sequence<H/4>());
    else
      return Result (sum, std::make_index_sequence<H>());
  }
};

template <size_t H, size_t W>
NGK_INLINE auto MultiVecScalAB_impl (size_t n, double ** ppa, double ** ppb)
{
  return KMultiVecScalAB<H,W>::Run (n, ppa, ppb);
}


// C = A * B^t (or conj(B)^t), complex rows
template <size_t H, size_t W, bool CONJ>
struct KMultiVecScalC
{
  typedef Mat<H,W,SIMD<double>> TSUM;
  static constexpr size_t SW = SIMD<double>::Size();

  static NGK_INLINE void Step (TSUM & sum, KPtrs<H,double> & pa, KPtrs<W,double> & pb, size_t i, SIMD<double> & conj)
  {
    Vec<H,SIMD<double>> a;
    KFor<0,H> ([] (auto k, Vec<H,SIMD<double>> & a, KPtrs<H,double> & pa, size_t & i) NGK_LAMBDA_INLINE
      { a.template Set<k> (SIMD<double>(pa.template Get<k>()+i)); }, a, pa, i);
    KFor<0,W> ([] (auto j, TSUM & sum, Vec<H,SIMD<double>> & a, KPtrs<W,double> & pb, size_t & i, SIMD<double> & conj) NGK_LAMBDA_INLINE
      {
        SIMD<double> b(pb.template Get<j>()+i);
        if constexpr (CONJ) b *= conj;
        KFor<0,H> ([] (auto k, auto j, TSUM & sum, Vec<H,SIMD<double>> & a, SIMD<double> & b) NGK_LAMBDA_INLINE
          {
            SIMD<double> s = sum.template Get<k,j>();
            FMAComplex (a.template Get<k>(), b, s);
            sum.template Set<k,j> (s);
          }, j, sum, a, b);
      }, sum, a, pb, i, conj);
  }

  static NGK_INLINE void Rest (TSUM & sum, KPtrs<H,double> & pa, KPtrs<W,double> & pb, size_t i, SIMD<double> & conj, SIMD<mask64> mask)
  {
    Vec<H,SIMD<double>> are, aim;
    KFor<0,H> ([] (auto k, Vec<H,SIMD<double>> & are, Vec<H,SIMD<double>> & aim,
                   KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> a(pa.template Get<k>()+i, mask);
        auto [re, im] = Unpack(a, a);
        are.template Set<k> (re);
        aim.template Set<k> (im);
      }, are, aim, pa, i, mask);
    KFor<0,W> ([] (auto j, TSUM & sum, Vec<H,SIMD<double>> & are, Vec<H,SIMD<double>> & aim,
                   KPtrs<W,double> & pb, size_t & i, SIMD<double> & conj, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> b(pb.template Get<j>()+i, mask);
        if constexpr (CONJ) b *= conj;
        SIMD<double> bswap = SwapPairs(b);
        KFor<0,H> ([] (auto k, auto j, TSUM & sum, Vec<H,SIMD<double>> & are, Vec<H,SIMD<double>> & aim,
                       SIMD<double> & b, SIMD<double> & bswap) NGK_LAMBDA_INLINE
          {
            SIMD<double> t = aim.template Get<k>() * bswap;
            sum.template Set<k,j> (sum.template Get<k,j>() + FMAddSub(are.template Get<k>(), b, t));
          }, j, sum, are, aim, b, bswap);
      }, sum, are, aim, pb, i, conj, mask);
  }

  static NGK_INLINE void Run (size_t n, Complex ** ppa, Complex ** ppb, Complex * pc, size_t dc)
  {
    KPtrs<H,double> pa = KLoadPtrs<H>((double**)ppa);
    KPtrs<W,double> pb = KLoadPtrs<W>((double**)ppb);
    SIMD<double> conj([](int i) { return (i%2==0) ? 1. : -1.; });
    TSUM sum;
    KZero<H,W> (sum);
    size_t i = 0;
    for ( ; i+SW <= 2*n; i += SW)
      Step (sum, pa, pb, i, conj);
    if (int r = (2*n) % SW)
      Rest (sum, pa, pb, i, conj, SIMD<mask64>(r));

    KFor<0,H> ([] (auto k, TSUM & sum, Complex *& pc, size_t & dc) NGK_LAMBDA_INLINE
      {
        KFor<0,W> ([] (auto j, auto k, TSUM & sum, Complex *& pc) NGK_LAMBDA_INLINE
          {
            SIMD<double,2> hsum = HSum128(sum.template Get<k,j>());
            pc[j].real(hsum[0]);
            pc[j].imag(hsum[1]);
          }, k, sum, pc);
        pc += dc;
      }, sum, pc, dc);
  }
};

template <size_t H, size_t W, bool CONJ>
NGK_INLINE void MultiVecScalC_impl (size_t n, Complex ** ppa, Complex ** ppb, Complex * pc, size_t dc)
{
  KMultiVecScalC<H,W,CONJ>::Run (n, ppa, ppb, pc, dc);
}


// A[i] += sum_j c(i,j) * B[j],  c(i,j) = pc[i+j*dc]
template <size_t H, size_t W>
struct KMultiScaleAdd
{
  typedef Mat<H,W,double> TC;
  typedef Vec<H,SIMD<double>> TA;
  static constexpr size_t SW = SIMD<double>::Size();

  template <bool MASK>
  static NGK_INLINE void Step (TC & c, KPtrs<H,double> & pa, KPtrs<W,double> & pb, size_t i, SIMD<mask64> mask)
  {
    TA a;
    KFor<0,H> ([] (auto k, TA & a, KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      { a.template Set<k> (MASK ? SIMD<double>(pa.template Get<k>()+i, mask) : SIMD<double>(pa.template Get<k>()+i)); }, a, pa, i, mask);
    KFor<0,W> ([] (auto j, TA & a, TC & c, KPtrs<W,double> & pb, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> b = MASK ? SIMD<double>(pb.template Get<j>()+i, mask) : SIMD<double>(pb.template Get<j>()+i);
        KFor<0,H> ([] (auto k, auto j, TA & a, TC & c, SIMD<double> & b) NGK_LAMBDA_INLINE
          { a.template Set<k> (a.template Get<k>() + c.template Get<k,j>() * b); }, j, a, c, b);
      }, a, c, pb, i, mask);
    KFor<0,H> ([] (auto k, TA & a, KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        if constexpr (MASK) a.template Get<k>().Store(pa.template Get<k>()+i, mask);
        else a.template Get<k>().Store(pa.template Get<k>()+i);
      }, a, pa, i, mask);
  }

  static NGK_INLINE void Run (size_t n, double ** ppa, double ** ppb, double * pc, size_t dc)
  {
    TC c;
    KFor<0,H> ([] (auto k, TC & c, double *& pc, size_t & dc) NGK_LAMBDA_INLINE
      { KFor<0,W> ([] (auto j, auto k, TC & c, double *& pc, size_t & dc) NGK_LAMBDA_INLINE
        { c.template Set<k,j> (pc[k+j*dc]); }, k, c, pc, dc); }, c, pc, dc);
    KPtrs<H,double> pa = KLoadPtrs<H>(ppa);
    KPtrs<W,double> pb = KLoadPtrs<W>(ppb);
    SIMD<mask64> mask(SW);
    size_t i = 0;
    for ( ; i+SW <= n; i += SW)
      Step<false> (c, pa, pb, i, mask);
    if (size_t r = n % SW)
      Step<true> (c, pa, pb, i, SIMD<mask64>(r));
  }
};

template <size_t H, size_t W>
NGK_INLINE void MultiScaleAdd_impl (size_t n, double ** ppa, double ** ppb, double * pc, size_t dc)
{
  KMultiScaleAdd<H,W>::Run (n, ppa, ppb, pc, dc);
}


// complex A[i] += sum_j c(i,j) * B[j]
template <size_t H, size_t W>
struct KMultiScaleAddC
{
  typedef Mat<H,W,SIMD<double>> TC;
  typedef Vec<H,SIMD<double>> TA;
  static constexpr size_t SW = SIMD<double>::Size();

  template <bool MASK>
  static NGK_INLINE void Step (TC & cre, TC & cim, KPtrs<H,double> & pa, KPtrs<W,double> & pb, size_t i, SIMD<mask64> mask)
  {
    TA a;
    KFor<0,H> ([] (auto k, TA & a, KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      { a.template Set<k> (MASK ? SIMD<double>(pa.template Get<k>()+i, mask) : SIMD<double>(pa.template Get<k>()+i)); }, a, pa, i, mask);
    KFor<0,W> ([] (auto j, TA & a, TC & cre, TC & cim, KPtrs<W,double> & pb, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> b = MASK ? SIMD<double>(pb.template Get<j>()+i, mask) : SIMD<double>(pb.template Get<j>()+i);
        SIMD<double> bswap = SwapPairs(b);
        KFor<0,H> ([] (auto k, auto j, TA & a, TC & cre, TC & cim, SIMD<double> & b, SIMD<double> & bswap) NGK_LAMBDA_INLINE
          {
            SIMD<double> t = cim.template Get<k,j>() * bswap;
            a.template Set<k> (a.template Get<k>() + FMAddSub(cre.template Get<k,j>(), b, t));
          }, j, a, cre, cim, b, bswap);
      }, a, cre, cim, pb, i, mask);
    KFor<0,H> ([] (auto k, TA & a, KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        if constexpr (MASK) a.template Get<k>().Store(pa.template Get<k>()+i, mask);
        else a.template Get<k>().Store(pa.template Get<k>()+i);
      }, a, pa, i, mask);
  }

  static NGK_INLINE void Run (size_t n, Complex ** ppa, Complex ** ppb, Complex * pc, size_t dc)
  {
    TC cre, cim;
    KFor<0,H> ([] (auto k, TC & cre, TC & cim, Complex *& pc, size_t & dc) NGK_LAMBDA_INLINE
      { KFor<0,W> ([] (auto j, auto k, TC & cre, TC & cim, Complex *& pc, size_t & dc) NGK_LAMBDA_INLINE
        {
          cre.template Set<k,j> (SIMD<double>(pc[k+j*dc].real()));
          cim.template Set<k,j> (SIMD<double>(pc[k+j*dc].imag()));
        }, k, cre, cim, pc, dc); }, cre, cim, pc, dc);
    KPtrs<H,double> pa = KLoadPtrs<H>((double**)ppa);
    KPtrs<W,double> pb = KLoadPtrs<W>((double**)ppb);
    SIMD<mask64> mask(SW);
    size_t i = 0;
    for ( ; i+SW <= 2*n; i += SW)
      Step<false> (cre, cim, pa, pb, i, mask);
    if (size_t r = (2*n) % SW)
      Step<true> (cre, cim, pa, pb, i, SIMD<mask64>(r));
  }
};

template <size_t H, size_t W>
NGK_INLINE void MultiScaleAddC_impl (size_t n, Complex ** ppa, Complex ** ppb, Complex * pc, size_t dc)
{
  KMultiScaleAddC<H,W>::Run (n, ppa, ppb, pc, dc);
}

template <size_t H, size_t W>
NGK_INLINE auto MultiVecScalAB (size_t n, double ** ppa, double ** ppb)
{
  return MultiVecScalAB_impl<H,W> (n, ppa, ppb);
}

template <size_t H, size_t W, bool conjugate>
NGK_INLINE void MultiVecScalC (size_t n, Complex ** ppa, Complex ** ppb, Complex * pc, size_t dc)
{
  MultiVecScalC_impl<H,W,conjugate> (n, ppa, ppb, pc, dc);
}

template <size_t H, size_t W>
NGK_INLINE void MultiScaleAdd (size_t n, double ** ppa, double ** ppb, double * pc, size_t dc)
{
  MultiScaleAdd_impl<H,W> (n, ppa, ppb, pc, dc);
}

template <size_t H, size_t W>
NGK_INLINE void MultiScaleAddC (size_t n, Complex ** ppa, Complex ** ppb, Complex * pc, size_t dc)
{
  MultiScaleAddC_impl<H,W> (n, ppa, ppb, pc, dc);
}



/* ************** ShortSum, AtB_SmallWA, Daxpy kernels ************** */

#if defined(__GNUC__) || defined(__clang__)
#define NGK_UNROLL4 _Pragma("GCC unroll 4")
#else
#define NGK_UNROLL4
#endif

// sum (op)= a*b, fused or as separate mul and add
template <OPERATION OP, bool FUSED, typename T>
NGK_INLINE T KUpdate (T a, T b, T sum)
{
  constexpr bool neg = OP == SUB || OP == SETNEG;
  if constexpr (FUSED)
    return neg ? FNMA(a, b, sum) : FMA(a, b, sum);
  else
    return neg ? sum - a*b : sum + a*b;
}

// C = A * B, A is ha x WA (short), C rows by panels of NV vectors of type TS
template <size_t WA, OPERATION OP>
struct KShortSum
{
  static constexpr int SW = SIMD<double>::Size();

  // one panel of C: NV vectors per row, the last one masked if MASK
  template <typename TS, int NV, bool MASK, bool FUSED, bool UNROLL>
  static NGK_INLINE void Panel (size_t ha, double * pa, size_t da, double * pb, size_t db,
                                double * pc, size_t dc, SIMD<mask64> mask)
  {
    constexpr int R = TS::Size();
    Mat<WA,NV,TS> b;
    KFor<0,int(WA)> ([] (auto k, Mat<WA,NV,TS> & b, double *& pb, size_t & db, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        KFor<0,NV> ([] (auto v, auto k, Mat<WA,NV,TS> & b, double *& pb, size_t & db, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
          {
            if constexpr (MASK && v == NV-1)
              b.template Set<k,v> (TS(pb+k*db+v*R, mask));
            else
              b.template Set<k,v> (TS(pb+k*db+v*R));
          }, k, b, pb, db, mask);
      }, b, pb, db, mask);

    __assume(ha>0);
    auto row = [] (Mat<WA,NV,TS> & b, double * pa2, double * pc2, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        Mat<1,NV,TS> sum;
        KFor<0,NV> ([] (auto v, Mat<1,NV,TS> & sum, double *& pc2, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
          {
            if constexpr (OP == SET || OP == SETNEG)
              sum.template Set<0,v> (TS(0.0));
            else if constexpr (MASK && v == NV-1)
              sum.template Set<0,v> (TS(pc2+v*R, mask));
            else
              sum.template Set<0,v> (TS(pc2+v*R));
          }, sum, pc2, mask);
        KFor<0,int(WA)> ([] (auto k, Mat<1,NV,TS> & sum, Mat<WA,NV,TS> & b, double *& pa2) NGK_LAMBDA_INLINE
          {
            KFor<0,NV> ([] (auto v, auto k, Mat<1,NV,TS> & sum, Mat<WA,NV,TS> & b, double *& pa2) NGK_LAMBDA_INLINE
              { sum.template Set<0,v> (KUpdate<OP,FUSED> (TS(pa2[k]), b.template Get<k,v>(), sum.template Get<0,v>())); },
              k, sum, b, pa2);
          }, sum, b, pa2);
        KFor<0,NV> ([] (auto v, Mat<1,NV,TS> & sum, double *& pc2, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
          {
            if constexpr (MASK && v == NV-1)
              sum.template Get<0,v>().Store(pc2+v*R, mask);
            else
              sum.template Get<0,v>().Store(pc2+v*R);
          }, sum, pc2, mask);
      };
    double * pa2 = pa, * pc2 = pc;
    if constexpr (UNROLL)
      {
        NGK_UNROLL4
        for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
          row (b, pa2, pc2, mask);
      }
    else
      for (size_t j = 0; j < ha; j++, pa2 += da, pc2 += dc)
        row (b, pa2, pc2, mask);
  }

  static NGK_INLINE void Run (size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
  {
    SIMD<mask64> mask(SW);
    for (size_t i = 0; i+SW <= wb; i += SW, pb += SW, pc += SW)
      Panel<SIMD<double>,1,false,true,true> (ha, pa, da, pb, db, pc, dc, mask);
    size_t rest = wb % SW;
    if (rest == 0) return;
    Panel<SIMD<double>,1,true,false,true> (ha, pa, da, pb, db, pc, dc, SIMD<mask64>(rest));
  }

  template <int R>
  static NGK_INLINE void Rest2 (size_t rest, size_t ha, double * pa, size_t da, double *& pb, size_t db, double *& pc, size_t dc)
  {
    if constexpr (R >= 1)
      {
        if constexpr (R <= SW)
          if (rest & R)
            {
              Panel<SIMD<double,R>,1,false,false,true> (ha, pa, da, pb, db, pc, dc, SIMD<mask64>(SW));
              pc += R;
              pb += R;
            }
        Rest2<R/2> (rest, ha, pa, da, pb, db, pc, dc);
      }
  }

  // B width unrolled by 2
  static NGK_INLINE void Run2 (size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
  {
    SIMD<mask64> mask(SW);
    for (size_t i = 0; i+2*SW <= wb; i += 2*SW, pb += 2*SW, pc += 2*SW)
      Panel<SIMD<double>,2,false,true,false> (ha, pa, da, pb, db, pc, dc, mask);
    size_t rest = wb % (2*SW);
    if (rest == 0) return;
    Rest2<8> (rest, ha, pa, da, pb, db, pc, dc);
  }
};

template <size_t WA, OPERATION OP>
NGK_INLINE void MatKernelShortSum_impl (size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  KShortSum<WA,OP>::Run (ha, wb, pa, da, pb, db, pc, dc);
}

template <size_t WA, OPERATION OP>
NGK_INLINE void MatKernelShortSum2_impl (size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  KShortSum<WA,OP>::Run2 (ha, wb, pa, da, pb, db, pc, dc);
}


// C = A^t * B, A is ha x WA (short), all updates fused
template <size_t WA, OPERATION OP>
struct KAtB_SmallWA
{
  static constexpr int SW = SIMD<double>::Size();

  template <int NV, bool MASK>
  static NGK_INLINE void Panel (size_t ha, double * pa, size_t da, double * pb, size_t db,
                                double * pc, size_t dc, SIMD<mask64> mask)
  {
    typedef Mat<WA,NV,SIMD<double>> TSUM;
    TSUM sum;
    double * pc2 = pc;
    KFor<0,int(WA)> ([] (auto k, TSUM & sum, double *& pc2, size_t & dc, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        KFor<0,NV> ([] (auto v, auto k, TSUM & sum, double *& pc2, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
          {
            if constexpr (OP == SET || OP == SETNEG)
              sum.template Set<k,v> (SIMD<double>(0.0));
            else if constexpr (MASK)
              sum.template Set<k,v> (SIMD<double>(pc2+v*SW, mask));
            else
              sum.template Set<k,v> (SIMD<double>(pc2+v*SW));
          }, k, sum, pc2, mask);
        pc2 += dc;
      }, sum, pc2, dc, mask);

    __assume(ha>0);
    auto row = [] (TSUM & sum, double * pa2, double * pb2, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        Mat<1,NV,SIMD<double>> b;
        KFor<0,NV> ([] (auto v, Mat<1,NV,SIMD<double>> & b, double *& pb2, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
          {
            if constexpr (MASK) b.template Set<0,v> (SIMD<double>(pb2+v*SW, mask));
            else b.template Set<0,v> (SIMD<double>(pb2+v*SW));
          }, b, pb2, mask);
        KFor<0,int(WA)> ([] (auto k, TSUM & sum, Mat<1,NV,SIMD<double>> & b, double *& pa2) NGK_LAMBDA_INLINE
          {
            KFor<0,NV> ([] (auto v, auto k, TSUM & sum, Mat<1,NV,SIMD<double>> & b, double *& pa2) NGK_LAMBDA_INLINE
              { sum.template Set<k,v> (KUpdate<OP,true> (b.template Get<0,v>(), SIMD<double>(pa2[k]), sum.template Get<k,v>())); },
              k, sum, b, pa2);
          }, sum, b, pa2);
      };
    double * pa2 = pa, * pb2 = pb;
    NGK_UNROLL4
    for (size_t j = 0; j < ha; j++, pa2 += da, pb2 += db)
      row (sum, pa2, pb2, mask);

    pc2 = pc;
    KFor<0,int(WA)> ([] (auto k, TSUM & sum, double *& pc2, size_t & dc, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        KFor<0,NV> ([] (auto v, auto k, TSUM & sum, double *& pc2, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
          {
            if constexpr (MASK) sum.template Get<k,v>().Store(pc2+v*SW, mask);
            else sum.template Get<k,v>().Store(pc2+v*SW);
          }, k, sum, pc2, mask);
        pc2 += dc;
      }, sum, pc2, dc, mask);
  }

  template <int NV0>
  static NGK_INLINE void Run (size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
  {
    SIMD<mask64> mask(SW);
    size_t i = 0;
    if constexpr (NV0 > 1)
      for ( ; i+NV0*SW <= wb; i += NV0*SW, pb += NV0*SW, pc += NV0*SW)
        Panel<NV0,false> (ha, pa, da, pb, db, pc, dc, mask);
    for ( ; i+SW <= wb; i += SW, pb += SW, pc += SW)
      Panel<1,false> (ha, pa, da, pb, db, pc, dc, mask);
    size_t rest = wb % SW;
    if (rest == 0) return;
    Panel<1,true> (ha, pa, da, pb, db, pc, dc, SIMD<mask64>(rest));
  }
};

template <size_t WA, OPERATION OP>
NGK_INLINE void MatKernelAtB_SmallWA_impl (size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  KAtB_SmallWA<WA,OP>::template Run<1> (ha, wb, pa, da, pb, db, pc, dc);
}

template <size_t WA, OPERATION OP>
NGK_INLINE void MatKernelAtB_SmallWA2_impl (size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  KAtB_SmallWA<WA,OP>::template Run<3> (ha, wb, pa, da, pb, db, pc, dc);
}


// C (op)= A * B, A is H x W (small, kept in registers), B is W x n
template <size_t H, size_t W, OPERATION OP>
struct KDaxpy
{
  static constexpr int SW = SIMD<double>::Size();
  typedef Mat<H,W,SIMD<double>> TA;
  typedef Vec<H,SIMD<double>> TC;

  template <bool MASK>
  static NGK_INLINE void Step (TA & a, KPtrs<W,double> & pb, KPtrs<H,double> & pc, size_t i, SIMD<mask64> mask)
  {
    TC c;
    KFor<0,int(H)> ([] (auto k, TC & c, KPtrs<H,double> & pc, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        if constexpr (OP == SET || OP == SETNEG)
          c.template Set<k> (SIMD<double>(0.0));
        else if constexpr (MASK)
          c.template Set<k> (SIMD<double>(pc.template Get<k>()+i, mask));
        else
          c.template Set<k> (SIMD<double>(pc.template Get<k>()+i));
      }, c, pc, i, mask);
    KFor<0,int(W)> ([] (auto j, TA & a, TC & c, KPtrs<W,double> & pb, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> b = MASK ? SIMD<double>(pb.template Get<j>()+i, mask) : SIMD<double>(pb.template Get<j>()+i);
        KFor<0,int(H)> ([] (auto k, auto j, TA & a, TC & c, SIMD<double> & b) NGK_LAMBDA_INLINE
          { c.template Set<k> (KUpdate<OP,!MASK> (a.template Get<k,j>(), b, c.template Get<k>())); }, j, a, c, b);
      }, a, c, pb, i, mask);
    KFor<0,int(H)> ([] (auto k, TC & c, KPtrs<H,double> & pc, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        if constexpr (MASK) c.template Get<k>().Store(pc.template Get<k>()+i, mask);
        else c.template Get<k>().Store(pc.template Get<k>()+i);
      }, c, pc, i, mask);
  }

  static NGK_INLINE void Run (size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
  {
    TA a;
    KFor<0,int(H)> ([] (auto k, TA & a, double *& pa, size_t & da) NGK_LAMBDA_INLINE
      { KFor<0,int(W)> ([] (auto j, auto k, TA & a, double *& pa, size_t & da) NGK_LAMBDA_INLINE
        { a.template Set<k,j> (SIMD<double>(pa[k*da+j])); }, k, a, pa, da); }, a, pa, da);
    KPtrs<H,double> pci;
    KFor<0,int(H)> ([] (auto k, KPtrs<H,double> & pci, double *& pc, size_t & dc) NGK_LAMBDA_INLINE
      { pci.template Set<k> (pc+k*dc); }, pci, pc, dc);
    KPtrs<W,double> pbj;
    KFor<0,int(W)> ([] (auto j, KPtrs<W,double> & pbj, double *& pb, size_t & db) NGK_LAMBDA_INLINE
      { pbj.template Set<j> (pb+j*db); }, pbj, pb, db);
    SIMD<mask64> mask(SW);
    size_t i = 0;
    for ( ; i+SW <= n; i += SW)
      Step<false> (a, pbj, pci, i, mask);
    Step<true> (a, pbj, pci, i, SIMD<mask64>(n%SW));
  }
};

template <size_t H, size_t W, OPERATION OP>
NGK_INLINE void MatKernelDaxpy_impl (size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  KDaxpy<H,W,OP>::Run (n, pa, da, pb, db, pc, dc);
}

template <size_t WA, OPERATION OP>
NGK_INLINE void MatKernelShortSum (size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  MatKernelShortSum_impl<WA,OP> (ha, wb, pa, da, pb, db, pc, dc);
}

template <size_t WA, OPERATION OP>
NGK_INLINE void MatKernelShortSum2 (size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  MatKernelShortSum2_impl<WA,OP> (ha, wb, pa, da, pb, db, pc, dc);
}

template <size_t WA, OPERATION OP>
NGK_INLINE void MatKernelAtB_SmallWA (size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  MatKernelAtB_SmallWA_impl<WA,OP> (ha, wb, pa, da, pb, db, pc, dc);
}

template <size_t WA, OPERATION OP>
NGK_INLINE void MatKernelAtB_SmallWA2 (size_t ha, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  MatKernelAtB_SmallWA2_impl<WA,OP> (ha, wb, pa, da, pb, db, pc, dc);
}

template <size_t H, size_t W, OPERATION OP>
NGK_INLINE void MatKernelDaxpy (size_t n, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
{
  MatKernelDaxpy_impl<H,W,OP> (n, pa, da, pb, db, pc, dc);
}



/* ****************** MatVec kernels, fixed width WA ****************** */

// y (op)= A * x,  MODE: SET, ADD, or KSCALE for y += s * A * x
constexpr int KSCALE = 16;

template <size_t WA, int MODE>
struct KMatVec
{
  static constexpr int SW = SIMD<double>::Size();
  static constexpr int NF = WA / SW;
  static constexpr int R = WA % SW;
  static constexpr bool R1 = SW == 4 && R == 1;   // scalar tail
  static constexpr bool R2 = SW == 4 && R == 2;   // half-vector tail
  static constexpr bool RM = R && !R1 && !R2;     // masked tail
  typedef Vec<NF,SIMD<double>> TX;
  typedef std::conditional_t<R1, double, std::conditional_t<R2, SIMD<double,2>, SIMD<double>>> TXT;

  // row sums of H rows, SW-vectors per row
  template <int H>
  static NGK_INLINE Vec<H,SIMD<double>> Rows (double * pa, size_t da, TX & xs, TXT & xt, SIMD<mask64> & mask)
  {
    typedef Vec<H,SIMD<double>> TS;
    TS sum;
    KFor<0,H> ([] (auto k, TS & sum) NGK_LAMBDA_INLINE { sum.template Set<k> (SIMD<double>(0.0)); }, sum);
    KFor<0,NF> ([] (auto b, TS & sum, double *& pa, size_t & da, TX & xs) NGK_LAMBDA_INLINE
      {
        KFor<0,H> ([] (auto k, auto b, TS & sum, double *& pa, size_t & da, TX & xs) NGK_LAMBDA_INLINE
          { sum.template Set<k> (sum.template Get<k>() + SIMD<double>(pa+k*da+b*SW) * xs.template Get<b>()); },
          b, sum, pa, da, xs);
      }, sum, pa, da, xs);
    KFor<0,H> ([] (auto k, TS & sum, double *& pa, size_t & da, TXT & xt, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        constexpr size_t o = NF*SW;
        if constexpr (R1 && H < 4)
          sum.template Set<k> (sum.template Get<k>() + SIMD<double,4> (pa[k*da+o] * xt, 0,0,0));
        else if constexpr (R2)
          sum.template Set<k> (sum.template Get<k>() + SIMD<double,4> (SIMD<double,2>(pa+k*da+o) * xt, SIMD<double,2>(0.0)));
        else if constexpr (RM)
          sum.template Set<k> (sum.template Get<k>() + SIMD<double>(pa+k*da+o, mask) * xt);
      }, sum, pa, da, xt, mask);
    return sum;
  }

  static NGK_INLINE void Run (double s, size_t ha, double * pa, size_t da, double * x, double * y)
  {
    TX xs;
    KFor<0,NF> ([] (auto b, TX & xs, double *& x) NGK_LAMBDA_INLINE
      { xs.template Set<b> (SIMD<double>(x+b*SW)); }, xs, x);
    SIMD<mask64> mask(R);
    TXT xt{};
    if constexpr (R1) xt = x[NF*SW];
    if constexpr (R2) xt = SIMD<double,2>(x+NF*SW);
    if constexpr (RM) xt = SIMD<double>(x+NF*SW, mask);

    size_t i = 0;
    for ( ; i+4 <= ha; i += 4, pa += 4*da)
      {
        auto sum = Rows<4> (pa, da, xs, xt, mask);
        SIMD<double,4> vsum = HSum(sum.template Get<0>(), sum.template Get<1>(), sum.template Get<2>(), sum.template Get<3>());
        if constexpr (R1)
          vsum += xt * SIMD<double,4> (pa[0*da+NF*SW], pa[1*da+NF*SW], pa[2*da+NF*SW], pa[3*da+NF*SW]);
        if constexpr (MODE == KSCALE)
          {
            vsum *= SIMD<double,4>(s);
            vsum += SIMD<double,4>(y+i);
          }
        else if constexpr (MODE == ADD)
          vsum += SIMD<double,4>(y+i);
        vsum.Store(y+i);
      }
    if (ha & 2)
      {
        auto sum = Rows<2> (pa, da, xs, xt, mask);
        SIMD<double,2> vsum = HSum(sum.template Get<0>(), sum.template Get<1>());
        if constexpr (MODE == KSCALE)
          {
            vsum *= SIMD<double,2>(s);
            vsum += SIMD<double,2>(y+i);
          }
        else if constexpr (MODE == ADD)
          vsum += SIMD<double,2>(y+i);
        vsum.Store(y+i);
        i += 2; pa += 2*da;
      }
    if (ha & 1)
      {
        auto sum = Rows<1> (pa, da, xs, xt, mask);
        if constexpr (MODE == KSCALE)
          y[i] += s*HSum(sum.template Get<0>());
        else if constexpr (MODE == ADD)
          y[i] += HSum(sum.template Get<0>());
        else
          y[i] = HSum(sum.template Get<0>());
      }
  }
};

template <size_t WA, OPERATION OP>
NGK_INLINE void KernelMatVec_impl (size_t ha, double * pa, size_t da, double * x, double * y)
{
  KMatVec<WA,OP>::Run (0.0, ha, pa, da, x, y);
}

template <size_t WA>
NGK_INLINE void KernelAddMatVec_impl (double s, size_t ha, double * pa, size_t da, double * x, double * y)
{
  KMatVec<WA,KSCALE>::Run (s, ha, pa, da, x, y);
}


// y(0:WA) += s * sum_i A(i,:) * x[ind[i]],  U independent accumulators
template <size_t WA>
struct KAddMatTransVecI
{
  static constexpr int SW = SIMD<double>::Size();
  static constexpr int NF = WA / SW;
  static constexpr int R = WA % SW;
  static constexpr int U = NF <= 1 ? 4 : (NF <= 4 ? 2 : 1);
  static constexpr bool R1 = R == 1;
  static constexpr bool R2 = R == 2 && SW == 4;
  typedef std::conditional_t<R1, SIMD<double,1>, std::conditional_t<R2, SIMD<double,2>, SIMD<double>>> TR;
  typedef Mat<NF,U,SIMD<double>> TY;
  typedef Mat<1,U,TR> TYR;

  template <int J>
  static NGK_INLINE void Row (TY & sy, TYR & syr, double * pa, double xv, SIMD<mask64> & mask)
  {
    SIMD<double> sx(xv);
    KFor<0,NF> ([] (auto b, TY & sy, double *& pa, SIMD<double> & sx) NGK_LAMBDA_INLINE
      { sy.template Set<b,J> (sy.template Get<b,J>() + SIMD<double>(pa+b*SW) * sx); }, sy, pa, sx);
    if constexpr (R1)
      syr.template Set<0,J> (syr.template Get<0,J>() + TR(pa+NF*SW) * sx[0]);
    else if constexpr (R2)
      syr.template Set<0,J> (syr.template Get<0,J>() + TR(pa+NF*SW) * sx.Lo());
    else if constexpr (R > 0)
      syr.template Set<0,J> (syr.template Get<0,J>() + SIMD<double>(pa+NF*SW, mask) * sx);
  }

  // remaining r < U rows go to accumulators r, r-1, ..., 1
  template <int J>
  static NGK_INLINE void Rest (size_t r, TY & sy, TYR & syr, double *& pa, size_t da, double * x, int * ind, size_t & i, SIMD<mask64> & mask)
  {
    if constexpr (J >= 1)
      {
        if (r >= J)
          {
            Row<J> (sy, syr, pa, x[ind[i]], mask);
            pa += da; i++;
          }
        Rest<J-1> (r, sy, syr, pa, da, x, ind, i, mask);
      }
  }

  static NGK_INLINE void Run (double s, size_t ha, double * pa, size_t da, double * x, double * y, int * ind)
  {
    TY sy;
    KZero<NF,U> (sy);
    TYR syr;
    KFor<0,U> ([] (auto j, TYR & syr) NGK_LAMBDA_INLINE { syr.template Set<0,j> (TR(0.0)); }, syr);
    SIMD<mask64> mask(R);

    size_t i;
    for (i = 0; i+U <= ha; i += U)
      KFor<0,U> ([] (auto j, TY & sy, TYR & syr, double *& pa, size_t & da, double *& x, int *& ind, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
        {
          Row<j> (sy, syr, pa, x[ind[i+j]], mask);
          pa += da;
        }, sy, syr, pa, da, x, ind, i, mask);
    Rest<U-1> (ha-i, sy, syr, pa, da, x, ind, i, mask);

    if constexpr (U > 1)
      {
        KFor<0,NF> ([] (auto b, TY & sy) NGK_LAMBDA_INLINE
          {
            if constexpr (U == 2)
              sy.template Set<b,0> (sy.template Get<b,0>() + sy.template Get<b,1>());
            else
              sy.template Set<b,0> (sy.template Get<b,0>() + (sy.template Get<b,1>() + sy.template Get<b,2>() + sy.template Get<b,3>()));
          }, sy);
        if constexpr (R > 0)
          {
            if constexpr (U == 2)
              syr.template Set<0,0> (syr.template Get<0,0>() + syr.template Get<0,1>());
            else
              syr.template Set<0,0> (syr.template Get<0,0>() + (syr.template Get<0,1>() + syr.template Get<0,2>() + syr.template Get<0,3>()));
          }
      }

    KFor<0,NF> ([] (auto b, TY & sy, double & s, double *& y) NGK_LAMBDA_INLINE
      { (s * sy.template Get<b,0>() + SIMD<double>(y+b*SW)).Store(y+b*SW); }, sy, s, y);
    if constexpr (R1 || R2)
      (s * syr.template Get<0,0>() + TR(y+NF*SW)).Store(y+NF*SW);
    else if constexpr (R > 0)
      (s * syr.template Get<0,0>() + SIMD<double>(y+NF*SW, mask)).Store(y+NF*SW, mask);
  }
};

template <size_t WA>
NGK_INLINE void KernelAddMatTransVecI_impl (double s, size_t ha, double * pa, size_t da, double * x, double * y, int * ind)
{
  KAddMatTransVecI<WA>::Run (s, ha, pa, da, x, y, ind);
}

template <size_t WA, OPERATION OP>
NGK_INLINE void KernelMatVec (size_t ha, double * pa, size_t da, double * x, double * y)
{
  KernelMatVec_impl<WA,OP> (ha, pa, da, x, y);
}

template <size_t WA>
NGK_INLINE void KernelAddMatVec (double s, size_t ha, double * pa, size_t da, double * x, double * y)
{
  KernelAddMatVec_impl<WA> (s, ha, pa, da, x, y);
}

template <size_t WA>
NGK_INLINE void KernelAddMatTransVecI (double s, size_t ha, double * pa, size_t da, double * x, double * y, int * ind)
{
  KernelAddMatTransVecI_impl<WA> (s, ha, pa, da, x, y, ind);
}



/* *************** Triangular kernels, DIM x DIM blocks *************** */

// X = T X or X = T^-1 X, T is lower-left or upper-right, X has wx columns
template <TRIG_SIDE SIDE, TRIG_NORMAL NORM, ORDERING ORD, int DIM>
struct KTriangular
{
  static constexpr int SW = SIMD<double>::Size();
  static constexpr bool LOWER = SIDE == LowerLeft;
  static constexpr bool NORMALIZED = NORM == Normalized;
  typedef Mat<DIM,DIM,SIMD<double>> TT;
  typedef Vec<DIM,SIMD<double>> TX;

  template <int I, int J> static constexpr bool Used ()
  { return LOWER ? (J < I || (J == I && !NORMALIZED)) : (J > I || (J == I && !NORMALIZED)); }

  static NGK_INLINE void LoadT (TT & t, double * pt, size_t dt)
  {
    KFor<0,DIM> ([] (auto i, TT & t, double *& pt, size_t & dt) NGK_LAMBDA_INLINE
      {
        KFor<0,DIM> ([] (auto j, auto i, TT & t, double *& pt, size_t & dt) NGK_LAMBDA_INLINE
          {
            if constexpr (Used<i,j>())
              t.template Set<i,j> (SIMD<double> (ORD == RowMajor ? pt[i*dt+j] : pt[i+dt*j]));
          }, i, t, pt, dt);
      }, t, pt, dt);
  }

  // row I of the operation, rows processed in the generator's order
  template <bool SOLVE, int I>
  static NGK_INLINE void Row (TT & t, TX & inv, TX & x)
  {
    SIMD<double> xi = x.template Get<I>();
    if constexpr (!SOLVE && !NORMALIZED)
      xi *= t.template Get<I,I>();
    KFor<0,DIM> ([] (auto j, TT & t, TX & x, SIMD<double> & xi) NGK_LAMBDA_INLINE
      {
        if constexpr (j != I && Used<I,j>())
          {
            if constexpr (SOLVE) xi -= t.template Get<I,j>() * x.template Get<j>();
            else xi += t.template Get<I,j>() * x.template Get<j>();
          }
      }, t, x, xi);
    if constexpr (SOLVE && !NORMALIZED)
      xi *= inv.template Get<I>();
    x.template Set<I> (xi);
  }

  template <bool SOLVE>
  static NGK_INLINE void Op (TT & t, TX & inv, TX & x)
  {
    // mult: lower bottom-up, upper top-down; solve the reverse
    constexpr bool down = SOLVE ? LOWER : !LOWER;
    KFor<0,DIM> ([] (auto k, TT & t, TX & inv, TX & x) NGK_LAMBDA_INLINE
      {
        constexpr int I = down ? int(k) : DIM-1-int(k);
        Row<SOLVE,I> (t, inv, x);
      }, t, inv, x);
  }

  template <bool SOLVE, bool MASK>
  static NGK_INLINE void Step (TT & t, TX & inv, double * px, size_t dx, SIMD<mask64> mask)
  {
    TX x;
    KFor<0,DIM> ([] (auto k, TX & x, double *& px, size_t & dx, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      { x.template Set<k> (MASK ? SIMD<double>(px+k*dx, mask) : SIMD<double>(px+k*dx)); }, x, px, dx, mask);
    Op<SOLVE> (t, inv, x);
    KFor<0,DIM> ([] (auto k, TX & x, double *& px, size_t & dx, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        constexpr bool skip = NORMALIZED && (LOWER ? k == 0 : k == DIM-1);
        if constexpr (!skip)
          {
            if constexpr (MASK) x.template Get<k>().Store(px+k*dx, mask);
            else x.template Get<k>().Store(px+k*dx);
          }
      }, x, px, dx, mask);
  }

  template <bool SOLVE>
  static NGK_INLINE void Run (size_t wx, double * pt, size_t dt, double * px, size_t dx)
  {
    TT t;
    LoadT (t, pt, dt);
    TX inv;
    if constexpr (SOLVE && !NORMALIZED)
      KFor<0,DIM> ([] (auto k, TT & t, TX & inv) NGK_LAMBDA_INLINE
        { inv.template Set<k> (1.0 / t.template Get<k,k>()); }, t, inv);
    SIMD<mask64> mask(SW);
    size_t i = 0;
    for ( ; i+SW <= wx; i += SW)
      Step<SOLVE,false> (t, inv, px+i, dx, mask);
    size_t rest = wx % SW;
    if (rest == 0) return;
    Step<SOLVE,true> (t, inv, px+i, dx, SIMD<mask64>(rest));
  }
};

template <TRIG_SIDE SIDE, TRIG_NORMAL NORM, ORDERING ORD, int DIM>
NGK_INLINE void KernelTriangularMult_impl (size_t wx, double * pt, size_t dt, double * px, size_t dx)
{
  KTriangular<SIDE,NORM,ORD,DIM>::template Run<false> (wx, pt, dt, px, dx);
}

template <TRIG_SIDE SIDE, TRIG_NORMAL NORM, ORDERING ORD, int DIM>
NGK_INLINE void KernelTriangularSolve_impl (size_t wx, double * pt, size_t dt, double * px, size_t dx)
{
  KTriangular<SIDE,NORM,ORD,DIM>::template Run<true> (wx, pt, dt, px, dx);
}


// Y (op)= T X, fused updates, a normalized diagonal counts as 1
template <TRIG_SIDE SIDE, TRIG_NORMAL NORM, ORDERING ORD, OPERATION OP, int DIM>
struct KTriangularMultXY
{
  static constexpr int SW = SIMD<double>::Size();
  static constexpr bool LOWER = SIDE == LowerLeft;
  static constexpr bool NORMALIZED = NORM == Normalized;
  typedef Mat<DIM,DIM,SIMD<double>> TT;
  typedef Vec<DIM,SIMD<double>> TX;

  static NGK_INLINE void LoadT (TT & t, double * pt, size_t dt)
  {
    KFor<0,DIM> ([] (auto i, TT & t, double *& pt, size_t & dt) NGK_LAMBDA_INLINE
      {
        KFor<0,DIM> ([] (auto j, auto i, TT & t, double *& pt, size_t & dt) NGK_LAMBDA_INLINE
          {
            if constexpr (j == i && NORMALIZED)
              t.template Set<i,j> (SIMD<double>(1.0));
            else if constexpr (LOWER ? j <= i : j >= i)
              t.template Set<i,j> (SIMD<double> (ORD == RowMajor ? pt[i*dt+j] : pt[i+dt*j]));
          }, i, t, pt, dt);
      }, t, pt, dt);
  }

  template <bool MASK>
  static NGK_INLINE void Step (TT & t, double * px, size_t dx, double * py, size_t dy, SIMD<mask64> mask)
  {
    TX x, y;
    KFor<0,DIM> ([] (auto k, TX & x, TX & y, double *& px, size_t & dx, double *& py, size_t & dy, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        x.template Set<k> (MASK ? SIMD<double>(px+k*dx, mask) : SIMD<double>(px+k*dx));
        if constexpr (OP == SET || OP == SETNEG)
          y.template Set<k> (SIMD<double>(0.0));
        else
          y.template Set<k> (MASK ? SIMD<double>(py+k*dy, mask) : SIMD<double>(py+k*dy));
      }, x, y, px, dx, py, dy, mask);
    // lower bottom-up, upper top-down, as generated
    KFor<0,DIM> ([] (auto k, TT & t, TX & x, TX & y) NGK_LAMBDA_INLINE
      {
        constexpr int I = LOWER ? DIM-1-int(k) : int(k);
        SIMD<double> yi = y.template Get<I>();
        IC<I> ic;
        KFor<0,DIM> ([] (auto j, auto i, TT & t, TX & x, SIMD<double> & yi) NGK_LAMBDA_INLINE
          {
            if constexpr (LOWER ? j <= i : j >= i)
              yi = KUpdate<OP,true> (t.template Get<i,j>(), x.template Get<j>(), yi);
          }, ic, t, x, yi);
        y.template Set<I> (yi);
      }, t, x, y);
    KFor<0,DIM> ([] (auto k, TX & y, double *& py, size_t & dy, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        if constexpr (MASK) y.template Get<k>().Store(py+k*dy, mask);
        else y.template Get<k>().Store(py+k*dy);
      }, y, py, dy, mask);
  }

  static NGK_INLINE void Run (size_t wx, double * pt, size_t dt, double * px, size_t dx, double * py, size_t dy)
  {
    TT t;
    LoadT (t, pt, dt);
    SIMD<mask64> mask(SW);
    size_t i = 0;
    for ( ; i+SW <= wx; i += SW)
      Step<false> (t, px+i, dx, py+i, dy, mask);
    size_t rest = wx % SW;
    if (rest == 0) return;
    Step<true> (t, px+i, dx, py+i, dy, SIMD<mask64>(rest));
  }
};

template <TRIG_SIDE SIDE, TRIG_NORMAL NORM, ORDERING ORD, OPERATION OP, int DIM>
NGK_INLINE void KernelTriangularMultXY_impl (size_t wx, double * pt, size_t dt, double * px, size_t dx, double * py, size_t dy)
{
  KTriangularMultXY<SIDE,NORM,ORD,OP,DIM>::Run (wx, pt, dt, px, dx, py, dy);
}

template <TRIG_SIDE SIDE, TRIG_NORMAL NORM, ORDERING ORD, int DIM>
NGK_INLINE void KernelTriangularMult (size_t wx, double * pt, size_t dt, double * px, size_t dx)
{
  KernelTriangularMult_impl<SIDE,NORM,ORD,DIM> (wx, pt, dt, px, dx);
}

template <TRIG_SIDE SIDE, TRIG_NORMAL NORM, ORDERING ORD, int DIM>
NGK_INLINE void KernelTriangularSolve (size_t wx, double * pt, size_t dt, double * px, size_t dx)
{
  KernelTriangularSolve_impl<SIDE,NORM,ORD,DIM> (wx, pt, dt, px, dx);
}

template <TRIG_SIDE SIDE, TRIG_NORMAL NORM, ORDERING ORD, OPERATION OP, int DIM>
NGK_INLINE void KernelTriangularMultXY (size_t wx, double * pt, size_t dt, double * px, size_t dx, double * py, size_t dy)
{
  KernelTriangularMultXY_impl<SIDE,NORM,ORD,OP,DIM> (wx, pt, dt, px, dx, py, dy);
}
