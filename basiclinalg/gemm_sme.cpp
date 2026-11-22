// double GEMM on the ARM SME unit (Apple M4 and later).
// Compiled with the normal flags; only the streaming kernels carry SME target attributes.
// If the CPU has SME, InitGEMM_SME replaces the general entries of dispatch_gemm.

#include <bla.hpp>

#if defined(__aarch64__) && defined(__clang__) && (defined(__APPLE__) || defined(__linux__))
#if __has_include(<arm_sme.h>)
#define NGS_SME_GEMM
#endif
#endif

#ifdef NGS_SME_GEMM
#include <arm_sme.h>
#ifdef __APPLE__
#include <sys/sysctl.h>
#else
#include <sys/auxv.h>
#ifndef HWCAP2_SME
#define HWCAP2_SME (1ul << 23)
#endif
#ifndef HWCAP2_SME_F64F64
#define HWCAP2_SME_F64F64 (1ul << 25)
#endif
#endif
#endif

namespace ngbla
{
  void SubAtDB_intern (SliceMatrix<double> a, SliceVector<double> diag,
                     SliceMatrix<double> b, SliceMatrix<double> c);
  extern void (*subatdb_hook) (SliceMatrix<double>, SliceVector<double>,
                               SliceMatrix<double>, SliceMatrix<double>);
  void InitGEMM_SME ();
  // mode 0: NEON only, 1: SME where profitable, 2: SME whenever possible; returns old mode
  NGS_DLL_HEADER int SetGEMM_SME (int mode);
}


#ifndef NGS_SME_GEMM

namespace ngbla
{
  void InitGEMM_SME () { }
  int SetGEMM_SME (int) { return 0; }
}

#else

namespace ngbla
{
  namespace
  {
#pragma clang attribute push (__attribute__((target("sme,sme-f64f64"))), apply_to = function)

    // C (m x n) op= A^T B, A: k x m, B: k x n, all row-major
    // C block of 16 x 32 in 8 ZA tiles: tiles 0-3 rows 0-7, tiles 4-7 rows 8-15
    // DIAG: C op= A^T diag(d) B, with the product a*d rounded first
    template <bool ADD, bool POS, bool DIAG = false>
    __arm_new("za") __arm_locally_streaming
    void AtB_Kernel (size_t k, size_t m, size_t n,
                     const double * a, size_t da, const double * b, size_t db,
                     double * c, size_t dc, const double * d = nullptr, size_t dd = 0)
    {
      const size_t vl = svcntd();

      for (size_t i = 0; i < m; i += 2*vl)
        {
          svbool_t pa0 = svwhilelt_b64_u64(i, m);
          svbool_t pa1 = svwhilelt_b64_u64(i+vl, m);
          size_t rows = m-i < 2*vl ? m-i : 2*vl;
          size_t rows0 = rows < vl ? rows : vl;
          size_t rows1 = rows - rows0;
          const bool two_a = rows1 > 0;

          for (size_t j = 0; j < n; j += 4*vl)
            {
              svbool_t pb0 = svwhilelt_b64_u64(j, n);
              svbool_t pb1 = svwhilelt_b64_u64(j+vl, n);
              svbool_t pb2 = svwhilelt_b64_u64(j+2*vl, n);
              svbool_t pb3 = svwhilelt_b64_u64(j+3*vl, n);
              const size_t nbv = (n-j+vl-1) / vl;
              double * c0 = c + i*dc + j;
              double * c1 = c0 + vl*dc;

              if constexpr (!ADD)
                svzero_za();
              else
                {
                  for (size_t r = 0; r < rows0; r++)
                    {
                      svld1_hor_za64(0, r, pb0, c0+r*dc);
                      svld1_hor_za64(1, r, pb1, c0+r*dc+vl);
                      svld1_hor_za64(2, r, pb2, c0+r*dc+2*vl);
                      svld1_hor_za64(3, r, pb3, c0+r*dc+3*vl);
                    }
                  for (size_t r = 0; r < rows1; r++)
                    {
                      svld1_hor_za64(4, r, pb0, c1+r*dc);
                      svld1_hor_za64(5, r, pb1, c1+r*dc+vl);
                      svld1_hor_za64(6, r, pb2, c1+r*dc+2*vl);
                      svld1_hor_za64(7, r, pb3, c1+r*dc+3*vl);
                    }
                }

#define MOP(t, pa, pb, va, vb) \
  if constexpr (POS) svmopa_za64_f64_m(t, pa, pb, va, vb); else svmops_za64_f64_m(t, pa, pb, va, vb)

              const double * al = a+i;
              const double * bl = b+j;
              for (size_t l = 0; l < k; l++, al += da, bl += db)
                {
                  svfloat64_t a0 = svld1_f64(pa0, al);
                  svfloat64_t a1 = svld1_f64(pa1, al+vl);
                  if constexpr (DIAG)
                    {
                      svfloat64_t dl = svdup_f64(d[l*dd]);
                      a0 = svmul_f64_x(pa0, a0, dl);
                      a1 = svmul_f64_x(pa1, a1, dl);
                    }
                  svfloat64_t b0 = svld1_f64(pb0, bl);
                  MOP(0, pa0, pb0, a0, b0);
                  if (two_a) { MOP(4, pa1, pb0, a1, b0); }
                  if (nbv > 1)
                    {
                      svfloat64_t b1 = svld1_f64(pb1, bl+vl);
                      MOP(1, pa0, pb1, a0, b1);
                      if (two_a) { MOP(5, pa1, pb1, a1, b1); }
                    }
                  if (nbv > 2)
                    {
                      svfloat64_t b2 = svld1_f64(pb2, bl+2*vl);
                      svfloat64_t b3 = svld1_f64(pb3, bl+3*vl);
                      MOP(2, pa0, pb2, a0, b2);
                      MOP(3, pa0, pb3, a0, b3);
                      if (two_a)
                        {
                          MOP(6, pa1, pb2, a1, b2);
                          MOP(7, pa1, pb3, a1, b3);
                        }
                    }
                }
#undef MOP

              for (size_t r = 0; r < rows0; r++)
                {
                  svst1_hor_za64(0, r, pb0, c0+r*dc);
                  svst1_hor_za64(1, r, pb1, c0+r*dc+vl);
                  svst1_hor_za64(2, r, pb2, c0+r*dc+2*vl);
                  svst1_hor_za64(3, r, pb3, c0+r*dc+3*vl);
                }
              for (size_t r = 0; r < rows1; r++)
                {
                  svst1_hor_za64(4, r, pb0, c1+r*dc);
                  svst1_hor_za64(5, r, pb1, c1+r*dc+vl);
                  svst1_hor_za64(6, r, pb2, c1+r*dc+2*vl);
                  svst1_hor_za64(7, r, pb3, c1+r*dc+3*vl);
                }
            }
        }
    }

    // dst (cols x rows) = Trans(src (rows x cols)), through ZA tiles 0-3
    __arm_new("za") __arm_locally_streaming
    void Transpose_Kernel (size_t rows, size_t cols, const double * src, size_t ds,
                           double * dst, size_t dd)
    {
      const size_t vl = svcntd();
      for (size_t i = 0; i < rows; i += vl)
        {
          svbool_t pr = svwhilelt_b64_u64(i, rows);
          size_t nr = std::min(rows-i, vl);
          const double * si = src + i*ds;
          for (size_t j = 0; j < cols; j += 4*vl)
            {
              svbool_t pc0 = svwhilelt_b64_u64(j, cols);
              svbool_t pc1 = svwhilelt_b64_u64(j+vl, cols);
              svbool_t pc2 = svwhilelt_b64_u64(j+2*vl, cols);
              svbool_t pc3 = svwhilelt_b64_u64(j+3*vl, cols);
              for (size_t r = 0; r < nr; r++)
                {
                  svld1_hor_za64(0, r, pc0, si+r*ds+j);
                  svld1_hor_za64(1, r, pc1, si+r*ds+j+vl);
                  svld1_hor_za64(2, r, pc2, si+r*ds+j+2*vl);
                  svld1_hor_za64(3, r, pc3, si+r*ds+j+3*vl);
                }
              size_t nc = std::min(cols-j, 4*vl);
              double * dj = dst + j*dd + i;
              for (size_t c = 0; c < std::min(nc, vl); c++)
                svst1_ver_za64(0, c, pr, dj+c*dd);
              for (size_t c = vl; c < std::min(nc, 2*vl); c++)
                svst1_ver_za64(1, c-vl, pr, dj+c*dd);
              for (size_t c = 2*vl; c < std::min(nc, 3*vl); c++)
                svst1_ver_za64(2, c-2*vl, pr, dj+c*dd);
              for (size_t c = 3*vl; c < nc; c++)
                svst1_ver_za64(3, c-3*vl, pr, dj+c*dd);
            }
        }
    }

#pragma clang attribute pop

    // keep the A and B panels in cache; k-blocks after the first accumulate
    constexpr size_t KB = 512, MB = 256, NB = 512;

    template <bool ADD, bool POS, bool DIAG = false>
    void AtB_Blocked (size_t k, size_t m, size_t n,
                      const double * a, size_t da, const double * b, size_t db,
                      double * c, size_t dc, const double * d = nullptr, size_t dd = 0)
    {
      for (size_t l = 0; l < k; l += KB)
        {
          size_t kb = std::min(k-l, KB);
          for (size_t i = 0; i < m; i += MB)
            {
              size_t mb = std::min(m-i, MB);
              for (size_t j = 0; j < n; j += NB)
                {
                  size_t nb = std::min(n-j, NB);
                  if (l == 0)
                    AtB_Kernel<ADD,POS,DIAG> (kb, mb, nb, a+i, da, b+j, db, c+i*dc+j, dc, d, dd);
                  else
                    AtB_Kernel<true,POS,DIAG> (kb, mb, nb, a+l*da+i, da, b+l*db+j, db, c+i*dc+j, dc, d+l*dd, dd);
                }
            }
        }
    }


    // Outer products need rows of A^T (k x m) and of B (k x n). PACKA: a is A (m x k),
    // PACKB: b is B^T (n x k); such operands are transposed panel by panel into a buffer.
    template <bool ADD, bool POS, bool PACKA, bool PACKB>
    void Packed_Blocked (size_t k, size_t m, size_t n,
                         const double * a, size_t da, const double * b, size_t db,
                         double * c, size_t dc)
    {
      constexpr size_t KBT = 256;
      thread_local std::unique_ptr<double[]> buf;
      thread_local size_t bufsize = 0;
      size_t mp = (std::min(m, MB)+7) & ~size_t(7);
      size_t np = (std::min(n, NB)+7) & ~size_t(7);
      size_t need = std::min(k, KBT) * (mp+np);
      if (bufsize < need)
        {
          buf.reset(new double[need+8]);
          bufsize = need;
        }
      // 64-byte aligned panels: unaligned SME loads cost up to 13% at small sizes
      double * pa = reinterpret_cast<double*>((reinterpret_cast<uintptr_t>(buf.get()) + 63) & ~uintptr_t(63));
      double * pb = pa + std::min(k, KBT)*mp;

      for (size_t l = 0; l < k; l += KBT)
        {
          size_t kb = std::min(k-l, KBT);
          auto block = [&] (size_t i, size_t j, const double * ai, size_t dai, const double * bj, size_t dbj)
            {
              size_t mb = std::min(m-i, MB), nb = std::min(n-j, NB);
              if (l == 0)
                AtB_Kernel<ADD,POS> (kb, mb, nb, ai, dai, bj, dbj, c+i*dc+j, dc);
              else
                AtB_Kernel<true,POS> (kb, mb, nb, ai, dai, bj, dbj, c+i*dc+j, dc);
            };

          if constexpr (PACKB)
            for (size_t j = 0; j < n; j += NB)
              {
                Transpose_Kernel (std::min(n-j, NB), kb, b+j*db+l, db, pb, np);
                for (size_t i = 0; i < m; i += MB)
                  if constexpr (PACKA)
                    {
                      Transpose_Kernel (std::min(m-i, MB), kb, a+i*da+l, da, pa, mp);
                      block (i, j, pa, mp, pb, np);
                    }
                  else
                    block (i, j, a+l*da+i, da, pb, np);
              }
          else
            for (size_t i = 0; i < m; i += MB)
              {
                Transpose_Kernel (std::min(m-i, MB), kb, a+i*da+l, da, pa, mp);
                for (size_t j = 0; j < n; j += NB)
                  block (i, j, pa, mp, b+l*db+j, db);
              }
        }
    }


    template <bool ADD, bool POS, ORDERING OA, ORDERING OB>
    pgemm<double,OA,OB> neon_gemm = nullptr;

    int sme_mode = 0;

    inline uint64_t Ticks ()
    {
      uint64_t t;
      asm volatile ("mrs %0, cntvct_el0" : "=r"(t));
      return t;
    }

    // The SME unit powers down when idle: a call after t us without SME pays about
    // min(0.5 + t/3, 9) us extra. t is taken from the cadence of SME-worthy GEMMs on this
    // thread, so a stream switches to SME if it would stay warm. Times in ns, M4 Pro.
    double ns_per_tick = 41.67;
    thread_local uint64_t gemm_last = 0;   // end of this thread's last SME-worthy GEMM

    // 0: NEON, 1: NEON but SME-worthy, 2: SME
    // kernels that round differently from NEON (ABt, SubAtDB) decide by size only,
    // so runs stay reproducible
    template <bool BITIDENTICAL>
    int UseSME (size_t m, size_t n, size_t k)
    {
      if (m < 16 || n < 16) return 0;
      if (sme_mode == 2) return 2;
      double w = double(m)*n*k;
      if constexpr (!BITIDENTICAL)
        return w >= 4e5 ? 2 : 0;
      double tneon = 0.0335*w;
      double tsme = 100 + 0.0041*w;
      if (tneon < tsme + 500) return 0;
      double wake = std::min(500 + (Ticks()-gemm_last)*ns_per_tick/3, 9000.0);
      return tneon > tsme + wake ? 2 : 1;
    }

    inline void Stamp (int use)
    {
      if (use) gemm_last = Ticks();
    }

    template <bool ADD, bool POS>
    void REGCALL AtB_SME (size_t ah, size_t aw, size_t bw,
                          BareSliceMatrix<double,ColMajor> a, BareSliceMatrix<double,RowMajor> b,
                          BareSliceMatrix<double,RowMajor> c) NETGEN_NOEXCEPT
    {
      int use = UseSME<true> (ah, bw, aw);
      if (use < 2)
        (*neon_gemm<ADD,POS,ColMajor,RowMajor>) (ah, aw, bw, a, b, c);
      else
        AtB_Blocked<ADD,POS> (aw, ah, bw, a.Data(), a.Dist(), b.Data(), b.Dist(), c.Data(), c.Dist());
      Stamp (use);
    }

    template <bool ADD, bool POS>
    void REGCALL ABt_SME (size_t ah, size_t aw, size_t bw,
                          BareSliceMatrix<double,RowMajor> a, BareSliceMatrix<double,ColMajor> b,
                          BareSliceMatrix<double,RowMajor> c) NETGEN_NOEXCEPT
    {
      int use = UseSME<false> (ah, bw, aw);
      if (use < 2)
        (*neon_gemm<ADD,POS,RowMajor,ColMajor>) (ah, aw, bw, a, b, c);
      else
        Packed_Blocked<ADD,POS,true,true> (aw, ah, bw, a.Data(), a.Dist(), b.Data(), b.Dist(), c.Data(), c.Dist());
      Stamp (use);
    }

    template <bool ADD, bool POS>
    void REGCALL AB_SME (size_t ah, size_t aw, size_t bw,
                         BareSliceMatrix<double,RowMajor> a, BareSliceMatrix<double,RowMajor> b,
                         BareSliceMatrix<double,RowMajor> c) NETGEN_NOEXCEPT
    {
      int use = UseSME<false> (ah, bw, aw);
      if (use < 2)
        (*neon_gemm<ADD,POS,RowMajor,RowMajor>) (ah, aw, bw, a, b, c);
      else
        Packed_Blocked<ADD,POS,true,false> (aw, ah, bw, a.Data(), a.Dist(), b.Data(), b.Dist(), c.Data(), c.Dist());
      Stamp (use);
    }

    template <bool ADD, bool POS>
    void REGCALL AtBt_SME (size_t ah, size_t aw, size_t bw,
                           BareSliceMatrix<double,ColMajor> a, BareSliceMatrix<double,ColMajor> b,
                           BareSliceMatrix<double,RowMajor> c) NETGEN_NOEXCEPT
    {
      int use = UseSME<false> (ah, bw, aw);
      if (use < 2)
        (*neon_gemm<ADD,POS,ColMajor,ColMajor>) (ah, aw, bw, a, b, c);
      else
        Packed_Blocked<ADD,POS,false,true> (aw, ah, bw, a.Data(), a.Dist(), b.Data(), b.Dist(), c.Data(), c.Dist());
      Stamp (use);
    }

    // c -= Trans(a) diag(d) b, the SparseCholesky update
    void SubAtDB_SME (SliceMatrix<double> a, SliceVector<double> diag,
                      SliceMatrix<double> b, SliceMatrix<double> c)
    {
      size_t k = a.Height(), m = a.Width(), n = b.Width();
      int use = UseSME<false> (m, n, k);
      if (use < 2)
        SubAtDB_intern (a, diag, b, c);
      else
        AtB_Blocked<true,false,true> (k, m, n, a.Data(), a.Dist(), b.Data(), b.Dist(),
                                      c.Data(), c.Dist(), diag.Data(), diag.Dist());
    }

    template <bool ADD, bool POS, ORDERING OA, ORDERING OB>
    void SwitchEntry (bool on, pgemm<double,OA,OB> sme)
    {
      constexpr size_t N = gemm_dispatch_size<OA,OB>;
      auto & entry = dispatch_gemm<ADD,POS,OA,OB>[N-1];
      auto & neon = neon_gemm<ADD,POS,OA,OB>;
      if (!neon) neon = entry;
      entry = on ? sme : neon;
    }

    template <bool ADD, bool POS>
    void Switch (bool on)
    {
      SwitchEntry<ADD,POS,ColMajor,RowMajor> (on, &AtB_SME<ADD,POS>);
      SwitchEntry<ADD,POS,RowMajor,ColMajor> (on, &ABt_SME<ADD,POS>);
      SwitchEntry<ADD,POS,RowMajor,RowMajor> (on, &AB_SME<ADD,POS>);
      SwitchEntry<ADD,POS,ColMajor,ColMajor> (on, &AtBt_SME<ADD,POS>);
    }

    bool HaveSME ()
    {
#ifdef __APPLE__
      auto feature = [] (const char * name)
      {
        int val = 0;
        size_t len = sizeof(val);
        return sysctlbyname(name, &val, &len, nullptr, 0) == 0 && val;
      };
      return feature("hw.optional.arm.FEAT_SME") && feature("hw.optional.arm.FEAT_SME_F64F64");
#else
      unsigned long hw2 = getauxval(AT_HWCAP2);
      return (hw2 & HWCAP2_SME) && (hw2 & HWCAP2_SME_F64F64);
#endif
    }

  }


  int SetGEMM_SME (int mode)
  {
    int old = sme_mode;
    if (!HaveSME()) mode = 0;
    bool on = mode > 0;
    Switch<false,false> (on);
    Switch<false,true> (on);
    Switch<true,false> (on);
    Switch<true,true> (on);
    subatdb_hook = on ? &SubAtDB_SME : nullptr;
    sme_mode = mode;
    return old;
  }

  void InitGEMM_SME ()
  {
    uint64_t freq;
    asm volatile ("mrs %0, cntfrq_el0" : "=r"(freq));
    ns_per_tick = 1e9 / freq;
    const char * env = getenv("NGS_SME");
    SetGEMM_SME (env ? atoi(env) : 1);
  }
}

#endif
