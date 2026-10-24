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
