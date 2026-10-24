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

// C (op)= A * B, C is H x W SIMDs
template <size_t H, size_t W, OPERATION OP, typename TB>
struct KMultAB
{
  typedef Mat<H,W,SIMD<double>> TSUM;
  typedef Mat<1,W,SIMD<double>> TB_ROW;
  static constexpr int SW = SIMD<double>::Size();
  static constexpr int SWB = std::is_same_v<TB,double> ? SW : 1;

  static NGK_INLINE void Step (TSUM & sum, double * pa, size_t da, TB * pb)
  {
    TB_ROW b;
    KFor<0,W> ([] (auto j, TB_ROW & b, TB *& pb) NGK_LAMBDA_INLINE
      { b.template Set<0,j> (KLoad(pb+j*SWB)); }, b, pb);
    KFor<0,H> ([] (auto i, TSUM & sum, TB_ROW & b, double *& pa, size_t & da) NGK_LAMBDA_INLINE
      {
        SIMD<double> a(pa[i*da]);
        KFor<0,W> ([] (auto j, auto i, TSUM & sum, TB_ROW & b, SIMD<double> & a) NGK_LAMBDA_INLINE
          {
            if constexpr (OP == SET || OP == ADD)
              sum.template Set<i,j> (FMA(a, b.template Get<0,j>(), sum.template Get<i,j>()));
            else
              sum.template Set<i,j> (FNMA(a, b.template Get<0,j>(), sum.template Get<i,j>()));
          }, i, sum, b, a);
      }, sum, b, pa, da);
  }

  static NGK_INLINE void Run (size_t n, double * pa, size_t da, TB * pb, size_t db, double * pc, size_t dc)
  {
    // C rows by a bumped pointer: index arithmetic costs MSVC registers
    TSUM sum;
    double * pci = pc;
    KFor<0,H> ([] (auto i, TSUM & sum, double *& pci, size_t & dc) NGK_LAMBDA_INLINE
      {
        KFor<0,W> ([] (auto j, auto i, TSUM & sum, double *& pci) NGK_LAMBDA_INLINE
          {
            if constexpr (OP == SET || OP == SETNEG)
              sum.template Set<i,j> (SIMD<double>(0.0));
            else
              sum.template Set<i,j> (SIMD<double>(pci+j*SW));
          }, i, sum, pci);
        pci += dc;
      }, sum, pci, dc);

    for (size_t k = 0; k < n; k++, pa++, pb += db)
      Step (sum, pa, da, pb);

    KFor<0,H> ([] (auto i, TSUM & sum, double *& pc, size_t & dc) NGK_LAMBDA_INLINE
      {
        KFor<0,W> ([] (auto j, auto i, TSUM & sum, double *& pc) NGK_LAMBDA_INLINE
          { sum.template Get<i,j>().Store(pc+j*SW); }, i, sum, pc);
        pc += dc;
      }, sum, pc, dc);
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
  typedef Mat<H,1,SIMD<double>> TSUM;

  static NGK_INLINE void Step (TSUM & sum, SIMD<mask64> mask, double * pa, size_t da, TB * pb)
  {
    SIMD<double> b((double*)pb, mask);
    KFor<0,H> ([] (auto i, TSUM & sum, SIMD<double> & b, double *& pa, size_t & da) NGK_LAMBDA_INLINE
      {
        SIMD<double> a(pa[i*da]);
        if constexpr (OP == SET || OP == ADD)
          sum.template Set<i,0> (FMA(a, b, sum.template Get<i,0>()));
        else
          sum.template Set<i,0> (FNMA(a, b, sum.template Get<i,0>()));
      }, sum, b, pa, da);
  }

  static NGK_INLINE void Run (size_t n, SIMD<mask64> mask, double * pa, size_t da, TB * pb, size_t db, double * pc, size_t dc)
  {
    TSUM sum;
    double * pci = pc;
    KFor<0,H> ([] (auto i, TSUM & sum, SIMD<mask64> & mask, double *& pci, size_t & dc) NGK_LAMBDA_INLINE
      {
        if constexpr (OP == SET || OP == SETNEG)
          sum.template Set<i,0> (SIMD<double>(0.0));
        else
          sum.template Set<i,0> (SIMD<double>(pci, mask));
        pci += dc;
      }, sum, mask, pci, dc);

    for (size_t k = 0; k < n; k++, pa++, pb += db)
      Step (sum, mask, pa, da, pb);

    KFor<0,H> ([] (auto i, TSUM & sum, SIMD<mask64> & mask, double *& pc, size_t & dc) NGK_LAMBDA_INLINE
      { sum.template Get<i,0>().Store(pc, mask); pc += dc; }, sum, mask, pc, dc);
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
    Mat<H,1,SIMD<double>> a;   // masked loads cannot be folded into the fma
    KFor<0,H> ([] (auto i, Mat<H,1,SIMD<double>> & a, double *& pa, size_t & da, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      { a.template Set<i,0> (SIMD<double>(pa+i*da, mask)); }, a, pa, da, mask);
    KFor<0,W> ([] (auto j, TSUM & sum, Mat<H,1,SIMD<double>> & a, double *& pb, size_t & db, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> b(pb+j*db, mask);
        KFor<0,H> ([] (auto i, auto j, TSUM & sum, Mat<H,1,SIMD<double>> & a, SIMD<double> & b) NGK_LAMBDA_INLINE
          { sum.template Set<i,j> (FMA(a.template Get<i,0>(), b, sum.template Get<i,j>())); }, j, sum, a, b);
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
using KPtrs = Mat<N,1,T*>;

template <size_t N, typename T>
NGK_INLINE KPtrs<N,T> KLoadPtrs (T ** pp)
{
  KPtrs<N,T> p;
  KFor<0,N> ([] (auto i, KPtrs<N,T> & p, T **& pp) NGK_LAMBDA_INLINE
    { p.template Set<i,0> (pp[i]); }, p, pp);
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
        SIMD<double> b = MASK ? SIMD<double>(pb.template Get<j,0>()+i, mask) : SIMD<double>(pb.template Get<j,0>()+i);
        KFor<0,H> ([] (auto k, auto j, TSUM & sum, KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask, SIMD<double> & b) NGK_LAMBDA_INLINE
          {
            SIMD<double> a = MASK ? SIMD<double>(pa.template Get<k,0>()+i, mask) : SIMD<double>(pa.template Get<k,0>()+i);
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
    Mat<H,1,SIMD<double>> a;
    KFor<0,H> ([] (auto k, Mat<H,1,SIMD<double>> & a, KPtrs<H,double> & pa, size_t & i) NGK_LAMBDA_INLINE
      { a.template Set<k,0> (SIMD<double>(pa.template Get<k,0>()+i)); }, a, pa, i);
    KFor<0,W> ([] (auto j, TSUM & sum, Mat<H,1,SIMD<double>> & a, KPtrs<W,double> & pb, size_t & i, SIMD<double> & conj) NGK_LAMBDA_INLINE
      {
        SIMD<double> b(pb.template Get<j,0>()+i);
        if constexpr (CONJ) b *= conj;
        KFor<0,H> ([] (auto k, auto j, TSUM & sum, Mat<H,1,SIMD<double>> & a, SIMD<double> & b) NGK_LAMBDA_INLINE
          {
            SIMD<double> s = sum.template Get<k,j>();
            FMAComplex (a.template Get<k,0>(), b, s);
            sum.template Set<k,j> (s);
          }, j, sum, a, b);
      }, sum, a, pb, i, conj);
  }

  static NGK_INLINE void Rest (TSUM & sum, KPtrs<H,double> & pa, KPtrs<W,double> & pb, size_t i, SIMD<double> & conj, SIMD<mask64> mask)
  {
    Mat<H,1,SIMD<double>> are, aim;
    KFor<0,H> ([] (auto k, Mat<H,1,SIMD<double>> & are, Mat<H,1,SIMD<double>> & aim,
                   KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> a(pa.template Get<k,0>()+i, mask);
        auto [re, im] = Unpack(a, a);
        are.template Set<k,0> (re);
        aim.template Set<k,0> (im);
      }, are, aim, pa, i, mask);
    KFor<0,W> ([] (auto j, TSUM & sum, Mat<H,1,SIMD<double>> & are, Mat<H,1,SIMD<double>> & aim,
                   KPtrs<W,double> & pb, size_t & i, SIMD<double> & conj, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> b(pb.template Get<j,0>()+i, mask);
        if constexpr (CONJ) b *= conj;
        SIMD<double> bswap = SwapPairs(b);
        KFor<0,H> ([] (auto k, auto j, TSUM & sum, Mat<H,1,SIMD<double>> & are, Mat<H,1,SIMD<double>> & aim,
                       SIMD<double> & b, SIMD<double> & bswap) NGK_LAMBDA_INLINE
          {
            SIMD<double> t = aim.template Get<k,0>() * bswap;
            sum.template Set<k,j> (sum.template Get<k,j>() + FMAddSub(are.template Get<k,0>(), b, t));
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
  typedef Mat<H,1,SIMD<double>> TA;
  static constexpr size_t SW = SIMD<double>::Size();

  template <bool MASK>
  static NGK_INLINE void Step (TC & c, KPtrs<H,double> & pa, KPtrs<W,double> & pb, size_t i, SIMD<mask64> mask)
  {
    TA a;
    KFor<0,H> ([] (auto k, TA & a, KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      { a.template Set<k,0> (MASK ? SIMD<double>(pa.template Get<k,0>()+i, mask) : SIMD<double>(pa.template Get<k,0>()+i)); }, a, pa, i, mask);
    KFor<0,W> ([] (auto j, TA & a, TC & c, KPtrs<W,double> & pb, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> b = MASK ? SIMD<double>(pb.template Get<j,0>()+i, mask) : SIMD<double>(pb.template Get<j,0>()+i);
        KFor<0,H> ([] (auto k, auto j, TA & a, TC & c, SIMD<double> & b) NGK_LAMBDA_INLINE
          { a.template Set<k,0> (a.template Get<k,0>() + c.template Get<k,j>() * b); }, j, a, c, b);
      }, a, c, pb, i, mask);
    KFor<0,H> ([] (auto k, TA & a, KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        if constexpr (MASK) a.template Get<k,0>().Store(pa.template Get<k,0>()+i, mask);
        else a.template Get<k,0>().Store(pa.template Get<k,0>()+i);
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
  typedef Mat<H,1,SIMD<double>> TA;
  static constexpr size_t SW = SIMD<double>::Size();

  template <bool MASK>
  static NGK_INLINE void Step (TC & cre, TC & cim, KPtrs<H,double> & pa, KPtrs<W,double> & pb, size_t i, SIMD<mask64> mask)
  {
    TA a;
    KFor<0,H> ([] (auto k, TA & a, KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      { a.template Set<k,0> (MASK ? SIMD<double>(pa.template Get<k,0>()+i, mask) : SIMD<double>(pa.template Get<k,0>()+i)); }, a, pa, i, mask);
    KFor<0,W> ([] (auto j, TA & a, TC & cre, TC & cim, KPtrs<W,double> & pb, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        SIMD<double> b = MASK ? SIMD<double>(pb.template Get<j,0>()+i, mask) : SIMD<double>(pb.template Get<j,0>()+i);
        SIMD<double> bswap = SwapPairs(b);
        KFor<0,H> ([] (auto k, auto j, TA & a, TC & cre, TC & cim, SIMD<double> & b, SIMD<double> & bswap) NGK_LAMBDA_INLINE
          {
            SIMD<double> t = cim.template Get<k,j>() * bswap;
            a.template Set<k,0> (a.template Get<k,0>() + FMAddSub(cre.template Get<k,j>(), b, t));
          }, j, a, cre, cim, b, bswap);
      }, a, cre, cim, pb, i, mask);
    KFor<0,H> ([] (auto k, TA & a, KPtrs<H,double> & pa, size_t & i, SIMD<mask64> & mask) NGK_LAMBDA_INLINE
      {
        if constexpr (MASK) a.template Get<k,0>().Store(pa.template Get<k,0>()+i, mask);
        else a.template Get<k,0>().Store(pa.template Get<k,0>()+i);
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
