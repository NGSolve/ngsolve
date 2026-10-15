/*********************************************************************/
/* File:   mptools_device.cpp                                        */
/* Purpose: the FMM passes (translations, P2M, L2P, M2P, P2P) on the */
/*          common GPU layer                                         */
/*********************************************************************/

#include "mptools.hpp"
#include "intop_parameters.hpp"
#include "fmmoperator.hpp"
#include <gpuwrapper.hpp>
#include <gpukernel.hpp>
#include <array>
#include <climits>
#include <map>
#include <mutex>

namespace ngsbem
{
  // embedded from the kernel headers by ngstd/embed_source.cmake
  extern const std::string code_rotatey_gpu;
  extern const std::string code_shiftz_gpu;
  extern const std::string code_transform_gpu;
  extern const std::string code_p2m_gpu;
  extern const std::string code_l2p_gpu;
  extern const std::string code_m2l_fused_gpu;

  namespace
  {
    // the real type of the device objects of one tree: float, or double (fp64) for Complex host expansions
    size_t RealBytes (bool fp64) { return fp64 ? sizeof(double) : sizeof(float); }
    size_t EntryBytes (bool is_complex, bool fp64) { return (is_complex ? 2 : 1)*RealBytes (fp64); }

    // a buffer of reals or complex numbers of that type (the kernels' REAL, Complex<REAL>); sizes count its elements
    struct RealBuffer
    {
      std::shared_ptr<ngs_gpu::Buffer> buf;
      bool fp64 = false, is_complex = false;

      explicit operator bool () const { return bool(buf); }
      size_t ElemBytes () const { return EntryBytes (is_complex, fp64); }
      size_t Size () const { return buf->Size() / ElemBytes(); }
      operator ngs_gpu::KernelArg () const
      {
        return ngs_gpu::KernelArg (buf).WithType ({ is_complex ? ngs_gpu::ArgType::Complex : ngs_gpu::ArgType::Float,
                                                    unsigned(ElemBytes()) });
      }
    };

    RealBuffer NewRealBuffer (const std::shared_ptr<ngs_gpu::Device> & device, size_t n, bool fp64, bool is_complex)
    {
      return { device->NewBuffer (n*EntryBytes (is_complex, fp64)), fp64, is_complex };
    }

    struct SqrtTables
    {
      RealBuffer buffer;
      int stride = 0;
    };

    std::string ReplaceAll (std::string src, const std::string & token, const std::string & value)
    {
      for (size_t p = src.find(token); p != std::string::npos; p = src.find(token, p+value.size()))
        src.replace(p, token.size(), value);
      return src;
    }

    // compiled once per device and source; published to the cache only once everything succeeded
    std::vector<std::shared_ptr<ngs_gpu::Kernel>>
    GetKernels (const std::shared_ptr<ngs_gpu::Device> & device, const std::string & source,
                std::initializer_list<const char*> names)
    {
      static std::map<std::pair<ngs_gpu::Device*, std::string>,
                      std::pair<std::shared_ptr<ngs_gpu::Library>, std::vector<std::shared_ptr<ngs_gpu::Kernel>>>> cache;
      auto key = std::make_pair (device.get(), source);
      if (auto it = cache.find(key); it != cache.end())
        return it->second.second;
      auto library = device->CompileSource (ngs_gpu::code_gpukernel + source);
      std::vector<std::shared_ptr<ngs_gpu::Kernel>> kernels;
      for (auto name : names)
        kernels.push_back (library->GetKernel (name));
      cache[key] = { library, kernels };
      return kernels;
    }

    // work-items per group: 1 on the CPU device (it runs a thread per work-item), else up to cap
    unsigned GroupThreads (const std::shared_ptr<ngs_gpu::Device> & device, size_t cap)
    {
      return device == ngs_gpu::GetCpuDevice() ? 1u : unsigned(std::min<size_t> (cap, device->MaxThreadsPerGroup()));
    }

    unsigned KernelThreads (const ngs_gpu::Kernel & kernel, unsigned t)
    {
      size_t kmax = kernel.Info(t).max_threads_per_group;
      return kmax > 0 && kmax < t ? unsigned(kmax) : t;
    }

    // work-items per group of a kernel whose group runs loops of up to width items between its barriers: width
    // rounded up to the SIMD width, at least one SIMD group, at most what the kernel allows (1 on the CPU device)
    struct WidthThreads
    {
      unsigned simd = 1, limit = 1;
      WidthThreads () = default;
      WidthThreads (const std::shared_ptr<ngs_gpu::Device> & device, const ngs_gpu::Kernel & kernel)
        : simd (unsigned(std::max<size_t> (device->SimdWidth(), 1))),
          limit (KernelThreads (kernel, GroupThreads (device, device->MaxThreadsPerGroup()))) { }
      unsigned operator() (size_t width) const
      { return unsigned(std::min<size_t> (std::max<size_t> ((width+simd-1)/simd, 1)*simd, limit)); }
    };

    // host on the device, at least one element (a valid kernel argument when empty)
    template <typename T>
    ngs_gpu::TypedBuffer<T> UploadBuffer (const std::shared_ptr<ngs_gpu::Device> & device, FlatArray<T> host)
    {
      ngs_gpu::TypedBuffer<T> buf = device->NewBuffer<T> (std::max<size_t> (host.Size(), 1));
      if (host.Size()) buf.H2D (host.Data(), host.Size());
      return buf;
    }

    // host values in double (T double or Complex) rounded once to the real type, uploaded; at least one element
    template <typename T>
    RealBuffer UploadReal (const std::shared_ptr<ngs_gpu::Device> & device, FlatArray<T> host, bool fp64)
    {
      constexpr bool is_complex = !std::is_same_v<T,double>;
      RealBuffer buf = NewRealBuffer (device, std::max<size_t> (host.Size(), 1), fp64, is_complex);
      if (!host.Size()) return buf;
      if (fp64)
        buf.buf->H2D (host.Data(), host.Size()*sizeof(T));
      else
        {
          using F = std::conditional_t<is_complex, std::complex<float>, float>;
          Array<F> rounded(host.Size());
          for (size_t i = 0; i < host.Size(); i++)
            rounded[i] = F(host[i]);
          buf.buf->H2D (rounded.Data(), rounded.Size()*sizeof(F));
        }
      return buf;
    }

    // the points p[3j..3j+2] as hi and lo parts xyz[6j..6j+5] of the real type (the device P2P); lo is 0 in double
    RealBuffer UploadHiLo (const std::shared_ptr<ngs_gpu::Device> & device, FlatArray<double> p, bool fp64)
    {
      Array<double> xyz(2*p.Size());
      for (size_t j = 0; j < p.Size()/3; j++)
        for (int k = 0; k < 3; k++)
          {
            double hi = fp64 ? p[3*j+k] : double(float(p[3*j+k]));
            xyz[6*j+k] = hi;
            xyz[6*j+3+k] = p[3*j+k] - hi;
          }
      return UploadReal (device, FlatArray<double>(xyz), fp64);
    }

    // the entries [first, first+n) of a buffer of complex C: f(dst, off, cnt) writes [first+off, first+off+cnt) into dst
    template <typename C, typename F>
    void FillBuffer (const RealBuffer & buf, size_t first, size_t n, F f)
    {
      buf.buf->H2D (n*sizeof(C), first*sizeof(C), [&] (void * dst, size_t off, size_t nb)
                    { f ((C*)dst, off/sizeof(C), nb/sizeof(C)); });
    }

    // f(data) with the first n entries of a buffer of complex C: in place if the host sees them (the CPU device, shared
    // memory), else a copy
    template <typename C, typename F>
    void ReadBuffer (const RealBuffer & buf, size_t n, F f)
    {
      if (auto host = buf.buf->HostData<C>())
        f (static_cast<const C*> (host));
      else
        {
          Array<C> tmp(n);
          buf.buf->D2H (tmp.Data(), n*sizeof(C));
          f (static_cast<const C*> (tmp.Data()));
        }
    }

    // a kernel source with REAL the real type; double only on a device with float64 (not Metal)
    std::string RealSource (const std::shared_ptr<ngs_gpu::Device> & device, const std::string & source, bool fp64)
    {
      if (fp64 && !device->HasFloat64())
        throw Exception("fmm device: " + device->Name() + " has no float64 for the fp64 FMM kernels, pass fp32=True");
      return ReplaceAll (source, "REAL", fp64 ? "double" : "float");
    }

    // rows and lanes per tile of the gather/scatter (bem/transform_gpu_kernel.h)
    constexpr size_t TRANSFORM_TILE = 32;

    auto GetTransformKernels (const std::shared_ptr<ngs_gpu::Device> & device, bool fp64)
    {
      return GetKernels (device, "#define TRANSFORM_TILE " + ToString(TRANSFORM_TILE) + "\n"
                         + RealSource (device, code_transform_gpu, fp64),
                         { "transform_gather", "transform_scatter", "fmm_zero" });
    }

    // S->R units with max(os,ot) up to FusedMaxOrder run fused (bem/m2l_fused_gpu_kernel.h), FUSED_LANES lanes per group
    constexpr int FUSED_M2L_MAX_ORDER = 27;
    constexpr int FUSED_LANES = 2;

    // the fused kernel's group memory for FUSED_ROWS = width^2: two buffers of FUSED_LANES expansions, within the 32 KiB
    // of Metal
    constexpr size_t FUSED_GROUP_BYTES = 32*1024;
    constexpr size_t FusedGroupBytes (int width, size_t real_bytes) { return 2*2*FUSED_LANES*size_t(width)*width*real_bytes; }
    static_assert (FusedGroupBytes (FUSED_M2L_MAX_ORDER+1, sizeof(float)) <= FUSED_GROUP_BYTES);

    // 27 in float, 21 in double
    int FusedMaxOrder (bool fp64)
    {
      int p = FUSED_M2L_MAX_ORDER;
      while (FusedGroupBytes (p+1, RealBytes (fp64)) > FUSED_GROUP_BYTES) p--;
      return p;
    }

    // one library per device, real type and row bound FUSED_ROWS >= (max(os,ot)+1)^2: in float 18 KiB of group memory up
    // to order 23, 25 KiB up to 27; in double 30 KiB up to 21
    auto GetFusedKernels (const std::shared_ptr<ngs_gpu::Device> & device, int max_order, bool fp64)
    {
      int width = std::min (max_order <= 23 ? 24 : FUSED_M2L_MAX_ORDER+1, FusedMaxOrder (fp64)+1);
      std::string defs = "#define FUSED_LANES " + ToString(FUSED_LANES) + "\n#define FUSED_ROWS " + ToString(sqr (width))
        + "\n#define TRANSFORM_TILE " + ToString(TRANSFORM_TILE) + "\n";
      return GetKernels (device, defs + RealSource (device, code_transform_gpu + code_m2l_fused_gpu, fp64),
                         { "rot_fold", "m2l_fused" });
    }

    // sqrt(k), 1/sqrt(k), 1/sqrt((2n+1)(2n+3)) for k,n < stride, in double and cast once as the CPU presqrt tables;
    // one buffer per device and real type, replaced for a larger order (rotations prepared earlier keep theirs alive)
    SqrtTables GetSqrtTables (const std::shared_ptr<ngs_gpu::Device> & device, int order, bool fp64)
    {
      static std::map<std::pair<ngs_gpu::Device*,bool>, SqrtTables> cache;
      int needed = 2*order+4;   // the kernels read indices up to 2*order+3
      auto key = std::make_pair (device.get(), fp64);
      if (auto it = cache.find(key); it != cache.end() && it->second.stride >= needed)
        return it->second;
      SqrtTables tables;
      tables.stride = needed;
      Array<double> host(3*size_t(needed));
      for (int k = 0; k < needed; k++)
        {
          host[k] = sqrt(double(k));
          host[needed+k] = k > 0 ? 1/sqrt(double(k)) : 0.0;
          host[2*needed+k] = 1/sqrt((2*k+1.0)*(2*k+3.0));
        }
      tables.buffer = UploadReal (device, FlatArray<double>(host), fp64);
      cache[key] = tables;   // only after the upload succeeded
      return tables;
    }

    // SphericalHarmonics::CalcAmn/CalcBmn in closed form, evaluated as the presqrt tables (sqrt of the exact integer
    // products, same operation order), so float(...) is the float CPU value; unlike CalcBmn without a table bound
    double ShiftZCoefAmn (int m, int n)
    {
      if (m < 0) m = -m;
      if (n < m) return 0;
      return std::sqrt(double(n+1+m))*std::sqrt(double(n+1-m)) * (1.0/std::sqrt(double(int64_t(2*n+1)*(2*n+3))));
    }
    double ShiftZCoefBmn (int m, int n)
    {
      if (m >= n || -m > n) return 0;
      double sgn = m >= 0 ? 1 : -1;
      int64_t k = int64_t(n)-m-1;   // sqrt(k(k+1)) = presqrt.sqrt_n_np1[k]
      return sgn * std::sqrt(double(k*(k+1))) * (1.0/std::sqrt(double(int64_t(2*n-1)*(2*n+1))));
    }

    // Amn, 1/Amn, Bmn(-m,n), Bmn(m-1,n), 1/Bmn(-m,n) for m <= max_m, n < stride, correctly rounded from the closed forms,
    // per device and real type; a larger request replaces them by the maximum of both sizes in each dimension
    struct Coef2DTables
    {
      RealBuffer buffer;
      int max_m = -1, stride = 0;
    };
    Coef2DTables GetCoef2DTables (const std::shared_ptr<ngs_gpu::Device> & device, int max_m, int stride, bool fp64)
    {
      static std::map<std::pair<ngs_gpu::Device*,bool>, Coef2DTables> cache;
      auto key = std::make_pair (device.get(), fp64);
      auto it = cache.find(key);
      if (it != cache.end() && it->second.max_m >= max_m && it->second.stride >= stride)
        return it->second;
      if (it != cache.end())
        {
          // grow, unless the merged table would exceed the kernels' int indexing (then the requested size, which was checked)
          int mm = std::max (max_m, it->second.max_m), ms = std::max (stride, it->second.stride);
          if (5*(size_t(mm)+1)*size_t(ms) <= INT_MAX)
            { max_m = mm; stride = ms; }
        }
      Coef2DTables tables;
      tables.max_m = max_m; tables.stride = stride;
      size_t block = size_t(max_m+1)*stride;
      Array<double> host(5*block);
      for (int m = 0; m <= max_m; m++)
        for (int n = 0; n < stride; n++)
          {
            size_t i = size_t(m)*stride+n;
            double amn = ShiftZCoefAmn(m,n), bneg = ShiftZCoefBmn(-m,n), bpos = ShiftZCoefBmn(m-1,n);
            host[i] = amn;
            host[block+i] = amn != 0 ? 1/amn : 0.0;
            host[2*block+i] = bneg;
            host[3*block+i] = bpos;
            host[4*block+i] = bneg != 0 ? 1/bneg : 0.0;
          }
      tables.buffer = UploadReal (device, FlatArray<double>(host), fp64);
      cache[key] = tables;
      return tables;
    }

    // the ShiftZ kernels index with int: a layer (H*W), the boundary buffer ((M+1)*(H+4W)) and the coefficient
    // tables (5*(M+1)*(H+1)); for os = ot the tables are the tightest, os = ot <= 14653
    void CheckShiftZOrders (int os, int ot)
    {
      size_t H = size_t(os)+ot+1, W = size_t(std::max(os,ot))+1, M = size_t(std::min(os,ot));
      if (H*W > INT_MAX || (M+1)*(H+4*W) > INT_MAX || 5*(M+1)*(H+1) > INT_MAX)
        throw Exception("FMMDevicePlan: os=" + ToString(os) + ", ot=" + ToString(ot)
                        + " exceed the int indexing of the device kernels (os = ot up to 14653)");
    }

    // NormalizedLegendreFunctions' recursion coefficients, CalcAmn and CalcBmn in double, in the layout of
    // p2m_gpu_kernel.h
    Array<double> P2MTables ()
    {
      constexpr int W = FMMDeviceP2M::max_order+1;
      size_t amn = 2*W*W+2*W, bmn = amn+(W+1)*(W+1);
      Array<double> tab(bmn+(W+1)*(2*W+1));
      tab = 0.0;
      for (int m = 0; m < W; m++)
        {
          for (int n = m+2; n < W; n++)
            {
              tab[n*W+m] = (2*n-1) / sqrt((n-m+0.0)*(n+m));
              tab[W*W+n*W+m] = sqrt((n+m-1.0)*(n-m-1.0)) / sqrt((n-m+0.0)*(n+m));
            }
          tab[2*W*W+m] = m > 0 ? sqrt((2*m-1.0)/(2*m)) : 0.0;
          tab[2*W*W+W+m] = sqrt(2*m+1.0);
        }
      for (int n = 0; n <= W; n++)
        for (int m = -W; m <= W; m++)
          {
            if (m >= 0) tab[amn+n*(W+1)+m] = ShiftZCoefAmn(m, n);
            tab[bmn+n*(2*W+1)+m+W] = ShiftZCoefBmn(m, n);
          }
      return tab;
    }

    // the group size bound of fmm_l2p and fmm_m2p (their group memory), defined for the kernels as L2P_THREADS
    constexpr unsigned L2P_THREADS = 64;

    struct LeafKernels
    {
      std::shared_ptr<ngs_gpu::Kernel> p2m, l2p, m2p, p2p, p2p_grad;
      unsigned p2m_threads = 1, eval_threads = 1;
      RealBuffer tables;
    };

    // per device, real type and kappa type: the P2M, L2P, M2P and P2P kernels for FMMDeviceP2M::max_order and their
    // tables. A GPU group has up to 256 P2M work-items, each owning (max_order+1)^2 / p2m_threads coefficients, and up to
    // L2P_THREADS others. Their group memory in double: fmm_p2m and fmm_l2p with a complex kappa 48 KiB, fmm_m2p 32 KiB.
    LeafKernels GetLeafKernels (const std::shared_ptr<ngs_gpu::Device> & device, bool complex_kappa, bool fp64)
    {
      constexpr int W = FMMDeviceP2M::max_order+1;
      LeafKernels k;
      k.p2m_threads = GroupThreads (device, 256);
      auto compile = [&] ()
      {
        std::string defs = "#define P2M_PMAX " + ToString(W-1) + "\n#define P2M_KMAX "
          + ToString((W*W+k.p2m_threads-1)/k.p2m_threads) + "\n#define P2M_STAGE 8\n#define L2P_THREADS "
          + ToString(L2P_THREADS) + "\n";
        std::string source = ReplaceAll (defs + code_p2m_gpu + code_l2p_gpu, "KAPPA_T",
                                         complex_kappa ? "Complex<REAL>" : "REAL");
        auto ks = GetKernels (device, RealSource (device, source, fp64),
                              { "fmm_p2m", "fmm_l2p", "fmm_m2p", "fmm_p2p", "fmm_p2p_grad" });
        k.p2m = ks[0]; k.l2p = ks[1]; k.m2p = ks[2]; k.p2p = ks[3]; k.p2p_grad = ks[4];
      };
      compile();
      if (unsigned t = KernelThreads (*k.p2m, k.p2m_threads); t < k.p2m_threads)
        {
          k.p2m_threads = t;      // fewer work-items, each with more coefficients
          compile();
          if (unsigned t2 = KernelThreads (*k.p2m, k.p2m_threads); t2 < k.p2m_threads)
            throw Exception("fmm device P2M: the device allows only " + ToString(t2) + " work-items per group");
        }
      k.eval_threads = GroupThreads (device, L2P_THREADS);
      for (auto & kern : { k.l2p, k.m2p, k.p2p, k.p2p_grad })
        k.eval_threads = KernelThreads (*kern, k.eval_threads);
      static std::map<std::pair<ngs_gpu::Device*,bool>, RealBuffer> tables;
      auto & t = tables[std::make_pair (device.get(), fp64)];
      if (!t)
        t = UploadReal (device, FlatArray<double>(P2MTables()), fp64);
      k.tables = t;
      return k;
    }

    ngs_gpu::KernelArg RealArg (double v, bool fp64)
    {
      return fp64 ? ngs_gpu::KernelArg (v) : ngs_gpu::KernelArg (float(v));
    }

    struct LeafKappa
    {
      std::complex<double> kappa, khat = 1.0;
      bool complex_kappa = false;
      LeafKappa () = default;
      LeafKappa (std::variant<double,std::complex<double>> k)
        : kappa (std::visit ([] (auto v) { return std::complex<double>(v); }, k)),
          complex_kappa (std::holds_alternative<std::complex<double>>(k))
      {
        if (complex_kappa && std::abs(kappa) > 0) khat = kappa / std::abs(kappa);
      }
      // a KAPPA_T value of the leaf kernels: REAL, or Complex<REAL> for a complex kappa
      ngs_gpu::KernelArg Arg (std::complex<double> v, bool fp64) const
      {
        if (!complex_kappa) return RealArg (v.real(), fp64);
        return fp64 ? ngs_gpu::KernelArg (v) : ngs_gpu::KernelArg (std::complex<float>(v));
      }
    };

    void CheckOrder (const char * where, int order)
    {
      // the kernels index the trafo blocks and the coefficient rows with ints
      if (order < 1 || (size_t(order)+2)*(2*size_t(order)+1) > INT_MAX)
        throw Exception(std::string(where) + ": order must be >= 1 and (order+2)*(2*order+1) <= INT_MAX, got "
                        + ToString(order));
    }

    // scratch reals per degree: the (order+1) x (2*order+1) trafo block plus 2*order+1 Dmn values
    int RotateYWorkStride (int order) { return (order+2)*(2*order+1); }

    // rows per item of rotatey_apply_b (bem/rotatey_gpu_kernel.h)
    constexpr int ROTATEY_ROWS = 4;

    // lanes per block of the multi-group applies: up to this many consecutive lanes of one translation group
    constexpr int DEVICE_LANE_BLOCK = 32;

    constexpr size_t FMM_DEVICE_LANE_BUFFER_BYTES = size_t(16) << 20;

    // bytes per workspace buffer: 16 MiB on the CPU device, else 1/64 of the device memory within [16 MiB, 1 GiB]
    // (1 GiB keeps the kernels' int offsets into the buffers valid)
    size_t FMMDeviceBudget (const std::shared_ptr<ngs_gpu::Device> & device)
    {
      if (device == ngs_gpu::GetCpuDevice()) return FMM_DEVICE_LANE_BUFFER_BYTES;
      return std::clamp (device->MemoryBytes() / 64, FMM_DEVICE_LANE_BUFFER_BYTES, size_t(1) << 30);
    }

    // the caps per unit (lanes, groups, RotateY degrees, ShiftZ layer slots): what fits 16 MiB, beyond that what fits
    // the budget up to this many parallel lane blocks / groups / degrees / layers per compute unit (more are no faster)
    constexpr size_t FMM_DEVICE_FILL = 4;

    // lanes per unit: each of the two lane buffers (rows x lanes Complex<REAL>) stays within budget bytes
    size_t FMMDeviceLaneCap (int max_order, bool fp64, size_t budget)
    {
      size_t rows = sqr (size_t(max_order)+1);
      return std::max (size_t(1), budget / (EntryBytes (true, fp64)*rows));
    }

    // groups per unit: their ShiftZ construction storage (layer 0, the boundaries and one layer slot per
    // group, which grow with G) stays within about budget bytes; at least one group, at most 65535 (the grid y dimension)
    size_t ShiftZGroupCap (int os, int ot, bool complex_trafo, bool fp64, size_t budget)
    {
      size_t H = size_t(os)+ot+1, W = size_t(std::max(os,ot))+1, M = size_t(std::min(os,ot));
      size_t per_group = (2*H*W + (M+1)*(H+4*W))*EntryBytes (complex_trafo, fp64);
      return std::max (size_t(1), std::min (size_t(65535), budget / per_group));
    }

    // reals of a rotation's folded blocks in a fused unit (FoldOffset(order+1) in bem/m2l_fused_gpu_kernel.h)
    size_t FoldReals (int order)
    {
      size_t n = size_t(order)+1;
      return n*(n+1)*(2*n+1)/3-2;
    }

    // groups per fused unit: all its layer slots within 16 MiB, the least ShiftZNeeds reserves, and each rotation's
    // folded blocks within a lane buffer (they share the workspace, see EnqueueFusedUnit)
    size_t FusedGroupCap (int os, int ot, bool fp64)
    {
      size_t H = size_t(os)+ot+1, W = size_t(std::max(os,ot))+1, M = size_t(std::min(os,ot));
      size_t layers = std::max (M, size_t(1))*H*W*EntryBytes (true, fp64);
      size_t fold = RealBytes (fp64)*std::max (FoldReals (os), FoldReals (ot));
      return std::max (size_t(1), std::min ({ ShiftZGroupCap (os, ot, true, fp64, FMM_DEVICE_LANE_BUFFER_BYTES),
                                              FMM_DEVICE_LANE_BUFFER_BYTES / layers, FMM_DEVICE_LANE_BUFFER_BYTES / fold }));
    }

    // SphericalExpansion::Scale of a host expansion of the real type (its basis_scale), singular or regular
    double ExpansionScale (bool singular, double abs_kappa, double rtyp, bool fp64)
    {
      if (fp64)
        return singular ? SphericalExpansion<Singular,Vec<1,Complex>>(-1, abs_kappa, rtyp).Scale()
                        : SphericalExpansion<Regular,Vec<1,Complex>>(-1, abs_kappa, rtyp).Scale();
      return singular ? SphericalExpansion<Singular,Vec<1,Complex32>>(-1, abs_kappa, rtyp).Scale()
                      : SphericalExpansion<Regular,Vec<1,Complex32>>(-1, abs_kappa, rtyp).Scale();
    }

    /*
      RotateY of a unit of G >= 1 rotations of one order, each with its own theta. Their lanes are the columns
      of in and out, group after group: (order+1)^2 rows, row-major, row n*(n+1)+m = Coef(n,m). Each lane gets
      the bits of its group run alone. Its workspace holds RotateYWorkStride(order) reals per group and degree of
      a chunk; launches sharing it must go to one in-order queue: the construct of a chunk of degrees overwrites
      the scratch that the previous apply reads.
    */
    struct DeviceRotation
    {
      int order;
      int groups, lanes;                     // G rotations; the columns of in/out, group after group
      RealBuffer legendre;                   // [g*(order+2)^2 + n*(order+2)+m] = normalized P_n^m(cos theta_g), odd m negated for theta_g < 0
      RealBuffer scal;                       // [3g..3g+2] = sin(theta_g), (1-cos)/2, (1+cos)/2, computed in double
      ngs_gpu::TypedBuffer<int> blocks;      // the unit's LaneBlocks
      int nblocks;
      RealBuffer tables;                     // sqrt(k), 1/sqrt(k), 1/sqrt((2n+1)(2n+3)); shared per device, see the kernel
      int table_stride;                      // length of each of the three tables
      unsigned threads;                      // work-items per group of the apply and the fold
      WidthThreads construct_threads;        // of a construct chunk: from the row width 2n+1 of its largest degree
      std::shared_ptr<ngs_gpu::Kernel> construct, apply;
    };

    // the lane blocks of a unit (group, first lane, lanes <= block_lanes), group after group, uploaded once
    // and shared by the unit's RotateY and ShiftZ applies (DEVICE_LANE_BLOCK) or its fused kernel (FUSED_LANES)
    struct LaneBlocks
    {
      ngs_gpu::TypedBuffer<int> buffer;
      int count = 0, lanes = 0;
    };

    LaneBlocks MakeLaneBlocks (const std::shared_ptr<ngs_gpu::Device> & device, FlatArray<int> lanes, int block_lanes)
    {
      size_t total = 0, nblocks = 0;
      for (size_t g = 0; g < lanes.Size(); g++)
        {
          if (lanes[g] < 1)
            throw Exception("FMMDevicePlan: group " + ToString(g) + " has no lanes");
          total += lanes[g];
          nblocks += (size_t(lanes[g])+block_lanes-1) / block_lanes;
        }
      if (lanes.Size() < 1 || total > INT_MAX || 3*nblocks > INT_MAX)
        throw Exception("FMMDevicePlan: a unit needs >= 1 groups, and its lanes must fit an int");
      Array<int> blocks(3*nblocks);
      size_t lane0 = 0, b = 0;
      for (size_t g = 0; g < lanes.Size(); g++)
        {
          for (int s0 = 0; s0 < lanes[g]; s0 += block_lanes, b++)
            {
              blocks[3*b] = int(g);
              blocks[3*b+1] = int(lane0)+s0;
              blocks[3*b+2] = std::min (block_lanes, lanes[g]-s0);
            }
          lane0 += lanes[g];
        }
      LaneBlocks lb;
      lb.buffer = UploadBuffer (device, FlatArray<int>(blocks));
      lb.count = int(nblocks);
      lb.lanes = int(total);
      return lb;
    }

    // the host side of one rotation, in double as RotateYImpl in mptools.cpp; rounded once to the real type at packing
    struct RotationHost
    {
      std::array<double,3> scal;    // sin(theta), (1-cos)/2, (1+cos)/2
      Array<double> legendre;       // [n*(order+2)+m] = normalized P_n^m(cos theta), odd m negated for theta < 0
    };

    RotationHost ComputeRotation (double theta, int order)
    {
      CheckOrder ("fmm device RotateY", order);
      if (!std::isfinite(theta))
        throw Exception("fmm device RotateY: theta is not finite");
      RotationHost rot;
      // leg(m,n) = normalized P_n^m(cos theta)
      double s = sin(theta);
      double c = cos(theta);
      Matrix<> leg(order+2, order+2);
      NormalizedLegendreFunctions (order+1, order+1, c, leg);
      if (theta < 0)
        for (int m = 1; m <= order+1; m += 2)
          leg.Row(m) *= -1;
      rot.scal = { s, (1-c)/2, (1+c)/2 };

      // the kernel reads legendre[n*(order+2)+m]: the transpose of leg
      int width = order+2;
      rot.legendre.SetSize (size_t(width)*width);
      for (int n = 0; n < width; n++)
        for (int m = 0; m < width; m++)
          rot.legendre[size_t(n)*width+m] = leg(m,n);
      return rot;
    }

    // a unit of G rotations (hosts[g]) on the lanes of blocks; share: a unit of the same rotations whose packed
    // Legendre and scalar tables this one reuses (the units of one group cut by the lane cap)
    std::shared_ptr<DeviceRotation>
    PackRotation (const std::shared_ptr<ngs_gpu::Device> & device, int order, FlatArray<const RotationHost*> hosts,
                  const LaneBlocks & blocks, const DeviceRotation * share, bool fp64)
    {
      size_t G = hosts.Size(), width = size_t(order)+2;
      // the kernels address group g's table at g*(order+2)^2 and its scratch at (g*nslots+k)*work_stride with ints;
      // the group index is the y grid dimension of rotatey_construct_g (65535 on CUDA)
      if (G < 1 || G > 65535 || G*width*width > INT_MAX || G*size_t(RotateYWorkStride (order)) > INT_MAX)
        throw Exception("fmm device RotateY: " + ToString(G) + " rotations of order " + ToString(order)
                        + " exceed the int indexing of the device kernels (1 <= G <= 65535, G*(order+2)^2 and "
                        "G*(order+2)*(2*order+1) must fit an int)");
      auto rot = std::make_shared<DeviceRotation>();
      rot->order = order;
      rot->groups = int(G);
      rot->lanes = blocks.lanes;
      rot->blocks = blocks.buffer;
      rot->nblocks = blocks.count;
      if (share)
        {
          rot->legendre = share->legendre;
          rot->scal = share->scal;
        }
      else
        {
          Array<double> legendre(G*width*width), scal(3*G);
          for (size_t g = 0; g < G; g++)
            {
              for (size_t i = 0; i < width*width; i++)
                legendre[g*width*width+i] = hosts[g]->legendre[i];
              for (size_t i = 0; i < 3; i++)
                scal[3*g+i] = hosts[g]->scal[i];
            }
          rot->legendre = UploadReal (device, FlatArray<double>(legendre), fp64);
          rot->scal = UploadReal (device, FlatArray<double>(scal), fp64);
        }

      auto tables = GetSqrtTables (device, order, fp64);
      rot->tables = tables.buffer;
      rot->table_stride = tables.stride;

      rot->threads = GroupThreads (device, 64);

      auto kernels = GetKernels (device, "#define ROTATEY_ROWS " + ToString(ROTATEY_ROWS) + "\n"
                                 + RealSource (device, code_rotatey_gpu, fp64),
                                 { "rotatey_construct_g", "rotatey_apply_b" });
      rot->construct = kernels[0];
      rot->apply = kernels[1];
      rot->construct_threads = WidthThreads (device, *rot->construct);
      return rot;
    }

    // parts per lane block of the multi-group applies: about 8 items per work-item for a full block, so that the few blocks
    // of a high-order unit still fill a GPU; within the int grid and loop index; 1 with one work-item per group
    int ApplyParts (size_t rows, unsigned threads, int nblocks)
    {
      if (threads <= 1 || nblocks < 1) return 1;
      size_t items = rows*DEVICE_LANE_BLOCK;
      if (items >= size_t(INT_MAX)) return 1;
      size_t parts = (items + 8*size_t(threads) - 1) / (8*size_t(threads));
      parts = std::min ({ parts, size_t(INT_MAX) / size_t(nblocks), (size_t(INT_MAX) - items) / threads });
      return int(std::max (parts, size_t(1)));
    }

    // RotateY scratch reals for a unit of G rotations of this order: chunks of 1 .. order degrees by the caps
    size_t RotateYNeeds (int order, size_t G, bool fp64, size_t budget, size_t fill)
    {
      size_t stride = size_t(RotateYWorkStride (order));
      auto degrees = [&] (size_t bytes) { return std::max (size_t(1), bytes / (RealBytes (fp64)*G*stride)); };
      size_t chunk = std::min (size_t(order), std::max (degrees (FMM_DEVICE_LANE_BUFFER_BYTES),
                                                        std::min (degrees (budget), std::max (size_t(1), fill / G))));
      return G*chunk*stride;
    }

    bool RotateYFits (const DeviceRotation & rot, const RealBuffer & ws)
    {
      return ws && ws.Size() >= size_t(rot.groups)*size_t(RotateYWorkStride (rot.order));
    }

    // after the checks, per chunk of degrees (nslots per group; the kernels address degree k of group g at
    // (g*nslots+k)*work_stride: ints) its construct into ws, then use(first, degrees, nslots, work_stride)
    template <typename F>
    void EnqueueRotateYChunks (ngs_gpu::Queue & queue, const DeviceRotation & rot, const RealBuffer & ws, F use)
    {
      int order = rot.order, work_stride = RotateYWorkStride (order);
      size_t G = rot.groups;
      if (!RotateYFits (rot, ws))
        throw Exception("fmm device RotateY: workspace too small for " + ToString(G) + " rotations of order "
                        + ToString(order));
      int nslots = int(std::min ({ size_t(order), ws.Size() / (G*work_stride), size_t(INT_MAX) / (G*work_stride) }));
      for (int first = 1; first <= order; first += nslots)
        {
          int degrees = std::min (nslots, order-first+1);
          queue.Launch (*rot.construct, ngs_gpu::Dim3(unsigned(degrees), unsigned(G)),
                        ngs_gpu::Dim3(rot.construct_threads (2*size_t(first+degrees-1)+1)),
                        { rot.legendre, rot.tables, rot.scal, ws, order, first, work_stride, rot.table_stride, nslots });
          use (first, degrees, nslots, work_stride);
        }
    }

    // in and out: distinct, (order+1)^2 * lanes entries each (the caller checks); launches only
    void EnqueueRotateY (ngs_gpu::Queue & queue, const DeviceRotation & rot,
                         const RealBuffer & ws, const RealBuffer & in, const RealBuffer & out)
    {
      EnqueueRotateYChunks (queue, rot, ws, [&] (int first, int degrees, int nslots, int work_stride)
      {
        int parts = ApplyParts ((2*size_t(first+degrees-1)+1+ROTATEY_ROWS-1)/ROTATEY_ROWS, rot.threads, rot.nblocks);
        queue.Launch (*rot.apply, ngs_gpu::Dim3(unsigned(rot.nblocks)*unsigned(parts), unsigned(degrees)), ngs_gpu::Dim3(rot.threads),
                      { ws, rot.blocks, in, out, first, work_stride, nslots, rot.lanes, parts });
      });
    }

    // the folded blocks of all degrees (rot_fold) into fold, groups*FoldReals(order) reals; launches only
    void EnqueueRotateYFold (ngs_gpu::Queue & queue, const DeviceRotation & rot, const RealBuffer & ws,
                             ngs_gpu::Kernel & fold_kernel, const ngs_gpu::KernelArg & fold)
    {
      EnqueueRotateYChunks (queue, rot, ws, [&] (int first, int degrees, int nslots, int work_stride)
      {
        queue.Launch (fold_kernel, ngs_gpu::Dim3(unsigned(degrees), unsigned(rot.groups)), ngs_gpu::Dim3(rot.threads),
                      { ws, fold, rot.order, first, work_stride, nslots });
      });
    }
  }


  /* ******************************* ShiftZ ******************************* */

  namespace
  {
    /*
      ShiftZ of a unit of G >= 1 translation groups, os -> ot, S->S, S->R or R->R, real or complex kappa; they
      share the variant, orders, kappa and radii, each has its own z. Layout, lanes and the in-order queue as
      RotateY; rows of out with |m| > os are zeroed.
      ComputeShiftZGroup mirrors SphericalExpansion::ShiftZ in double and rounds once: the O(p) seed column,
      powscale and the recursion scalars. The recursion coefficients are the correctly rounded 2D tables
      (GetCoef2DTables): 1D sqrt-table products or an in-kernel sqrt lose accuracy at high order.
      In float the recursion runs at the balanced scales of SphericalExpansion::ShiftZ (Metal flushes float
      subnormals), converted exactly by fin(l) and fout(n) in the apply; double keeps the expansion scales.
    */
    struct DeviceShiftZ
    {
      FMMTranslation variant;
      int os, ot;
      bool complex_trafo;                    // S->R or complex kappa: Complex<REAL> layers, else REAL
      bool fp64;
      int groups, lanes;                     // G translation groups; the columns of in/out, group after group
      // per group g, uniform strides (H = os+ot+1, W = max(os,ot)+1):
      RealBuffer seed;                             // TRAFO[g*H+l]: layer-0 first column
      RealBuffer pow_mant;                         // g*(W+os+ot+2) + [powscale mantissa, k < W | fin(l), l <= os | fout mantissa, n <= ot]
      ngs_gpu::TypedBuffer<int> pow_exp;           // g*(W+ot+1) + [powscale exponent, k < W | fout exponent, n <= ot]
      RealBuffer scal;                             // g*9 + the recursion scalars, in the order of ComputeShiftZGroup
      ngs_gpu::TypedBuffer<int> blocks;            // the unit's LaneBlocks
      int nblocks;
      RealBuffer ctab;                             // Amn, 1/Amn, Bmn(-m,n), Bmn(m-1,n), 1/Bmn(-m,n) tables, shared per device
      int cstride, cblock;                         // entry m*cstride+n, block size cblock (see the kernel source)
      unsigned threads_boundaries, threads_layers; // WidthThreads of the column length H, of the layer row width W
      unsigned threads_apply;                      // work-items per group of the apply
      std::shared_ptr<ngs_gpu::Kernel> boundaries, layers, apply;
    };

    struct ShiftZWorkspace
    {
      // per group of a unit (H = os+ot+1, W = max(os,ot)+1, M = min(os,ot)): layer 0 (H*W entries), the boundaries
      // ((M+1)*(H+4W)) and K >= 1 layer slots of H*W
      std::shared_ptr<ngs_gpu::Buffer> layer0, boundary, layers;
    };

    // group g of DeviceShiftZ on the host, all in double; rounded once to the real type at packing
    struct ShiftZGroupHost
    {
      bool complex_trafo;
      std::vector<Complex> cseed;                // complex layers: the H seed entries
      std::vector<double> rseed;                 // real layers: the H seed entries
      std::vector<double> mant;
      std::vector<int> expo;
      std::array<double,9> scal;
    };

    ShiftZGroupHost ComputeShiftZGroup (FMMTranslation variant, std::variant<double,std::complex<double>> kappa_v, double z,
                                        double source_rtyp, double target_rtyp, int os, int ot, bool fp64)
    {
      if (!std::isfinite(z) || !(source_rtyp > 0) || !(target_rtyp > 0))
        throw Exception("fmm device ShiftZ: z must be finite and the radii positive");
      bool is_ss = variant == FMMTranslation::SS, is_sr = variant == FMMTranslation::SR, is_rr = variant == FMMTranslation::RR;
      bool complex_kappa = std::holds_alternative<std::complex<double>>(kappa_v);
      std::complex<double> kappa = complex_kappa ? std::get<std::complex<double>>(kappa_v) : std::complex<double>(std::get<double>(kappa_v));
      if (!std::isfinite(std::abs(kappa)))
        throw Exception("fmm device ShiftZ: kappa is not finite");

      ShiftZGroupHost sh;
      sh.complex_trafo = is_sr || complex_kappa;

      // as SphericalExpansion::ShiftZ, in double; the source is Singular for S->S and S->R, the target for S->S only
      double abs_kappa = std::abs(kappa);
      double scale = ExpansionScale (!is_rr, abs_kappa, source_rtyp, fp64);
      double inv_scale = 1/scale;
      double tscale = ExpansionScale (is_ss, abs_kappa, target_rtyp, fp64);
      double inv_tscale = 1/tscale;
      double hankel_prefactor = 1;
      if (is_sr && abs_kappa > 0 && abs_kappa < 1)
        hankel_prefactor = abs_kappa;
      double inv_hankel_prefactor = 1/hankel_prefactor;
      if (is_sr && (!(abs_kappa*std::abs(z) > 0) || !std::isfinite(abs_kappa*std::abs(z))))
        throw Exception("fmm device ShiftZ: S->R needs a finite kappa*z != 0");
      // float at the balanced scales (see DeviceShiftZ): fin(l) = (1/hp) rin^l <= 1/hp, fout(n) = rout^n <= 1
      double rin = 1, rout = 1;
      if (!fp64 && is_sr)
        {
          double lambda = std::max ({ std::ldexp(1.0, BalancedScaleExponent (os+ot, kappa*std::abs(z))),
                                      scale, inv_tscale });
          rin = scale/lambda; rout = inv_tscale/lambda;
          scale = lambda; inv_scale = 1/lambda; tscale = 1/lambda; inv_tscale = lambda;
        }
      if (!fp64 && is_ss)
        {
          rin = scale/tscale;
          scale = tscale; inv_scale = inv_tscale;
        }
      if (!fp64 && is_rr)
        {
          rout = scale/tscale;
          tscale = scale; inv_tscale = inv_scale;
        }
      if (is_rr)
        {
          std::swap (scale, inv_tscale);
          std::swap (tscale, inv_scale);
          z = -z;
        }
      int H = os+ot+1, W = std::max(os,ot)+1;

      // the seed column as SphericalExpansion::ShiftZ, in the real or complex double arithmetic of the CPU trafo_type
      auto make_seed = [&] (auto kappa_t, auto & out)
      {
        using TT = typename std::decay_t<decltype(out)>::value_type;
        Vector<TT> seed(H);
        seed = TT(0.0);
        if constexpr (std::is_same_v<TT,Complex>)
          if (is_sr) SphericalHankel1 (os+ot, kappa_t*std::abs(z), inv_tscale, seed, hankel_prefactor);
        if (!is_sr) SphericalBessel (os+ot, kappa_t*std::abs(z), tscale, seed);
        for (int l = 0; l < H; l++)
          out.push_back (seed(l) * (z < 0 && l % 2 ? -1.0 : 1.0) * sqrt(2*l+1));
      };
      if (!sh.complex_trafo) make_seed (kappa.real(), sh.rseed);
      else if (complex_kappa) make_seed (kappa, sh.cseed);
      else make_seed (kappa.real(), sh.cseed);
      if (is_sr)
        for (Complex v : sh.cseed)
          {
            if (!fp64) v = Complex(std::complex<float>(v));
            if (std::isfinite(v.real()) && std::isfinite(v.imag())) continue;
            std::string real = fp64 ? "double" : "float";
            throw Exception("fmm device ShiftZ: the S->R seed overflows " + real + " (order os+ot=" + ToString(os+ot)
                            + " too high for " + real + " at kappa*|z|=" + ToString(abs_kappa*std::abs(z))
                            + "; a low-frequency high-order translation" + (fp64 ? ")" : ", where fp64 breaks down too)"));
          }

      sh.scal = { scale, tscale, inv_scale, inv_tscale, scale*tscale, scale*inv_tscale, scale*scale, inv_tscale*inv_tscale,
                  inv_scale*inv_tscale };

      // powscale(i) = (-scale*tscale)^i (balanced S->R: (-1)^i) and fout(n) can leave the range of REAL: a mantissa in
      // [1,2) and an int exponent, zero stays zero. fp64: fin = 1/hp, fout = 2^0, the unscaled apply's bits.
      sh.mant.resize (W + os+1 + ot+1);
      sh.expo.resize (W + ot+1);
      double scaleprod = -scale*tscale, prod = 1;
      for (int i = 0; i < W; i++, prod *= scaleprod)
        {
          int e = 0;
          double f = prod != 0 ? std::frexp (prod, &e) : 0.0;
          sh.mant[i] = 2*f;
          sh.expo[i] = prod != 0 ? e-1 : 0;
        }
      for (int l = 0; l <= os; l++)
        sh.mant[W+l] = rin == 1 ? inv_hankel_prefactor : std::exp2(std::log2(inv_hankel_prefactor) + l*std::log2(rin));
      for (int n = 0; n <= ot; n++)
        {
          double lf = n*std::log2(rout), fe = std::floor(lf);
          sh.mant[W+os+1+n] = std::exp2(lf-fe);
          sh.expo[W+n] = int(fe);
        }
      return sh;
    }

    // one unit of G translation groups (hosts[g]) on the lanes of blocks: the per-group tables packed with uniform
    // strides, one upload per table kind (bem/shiftz_gpu_kernel.h)
    std::shared_ptr<DeviceShiftZ>
    PackShiftZ (const std::shared_ptr<ngs_gpu::Device> & device, FMMTranslation variant, int os, int ot,
                FlatArray<const ShiftZGroupHost*> hosts, const LaneBlocks & blocks, bool fp64)
    {
      size_t G = hosts.Size();
      size_t H = size_t(os)+ot+1, W = size_t(std::max(os,ot))+1, M = size_t(std::min(os,ot));
      size_t nmant = W+os+ot+2, nexpo = W+ot+1;
      // the kernels address the groups' tables and workspace with ints (g*H*W, g*(M+1)*(H+4W), ...); the group
      // index is the y grid dimension of shiftz_layers_g (65535 on CUDA)
      if (G > 65535 || G*H*W > INT_MAX || G*(M+1)*(H+4*W) > INT_MAX || G*nmant > INT_MAX)
        throw Exception("fmm device ShiftZ: " + ToString(G) + " groups of os=" + ToString(os) + ", ot=" + ToString(ot)
                        + " exceed the int indexing of the device kernels (G <= 65535, G*(os+ot+1)*(max(os,ot)+1) and "
                        "G*(min(os,ot)+1)*(os+ot+1+4*(max(os,ot)+1)) must fit an int)");
      bool complex_trafo = hosts[0]->complex_trafo;

      auto sh = std::make_shared<DeviceShiftZ>();
      sh->variant = variant;
      sh->os = os; sh->ot = ot;
      sh->complex_trafo = complex_trafo;
      sh->fp64 = fp64;
      sh->groups = int(G);
      sh->lanes = blocks.lanes;
      sh->blocks = blocks.buffer;
      sh->nblocks = blocks.count;

      Array<double> mant(G*nmant), scal(G*9);
      Array<int> expo(G*nexpo);
      Array<Complex> cseed(complex_trafo ? G*H : 0);
      Array<double> rseed(complex_trafo ? 0 : G*H);
      for (size_t g = 0; g < G; g++)
        {
          const ShiftZGroupHost & h = *hosts[g];
          for (size_t i = 0; i < H; i++)
            {
              if (complex_trafo) cseed[g*H+i] = h.cseed[i];
              else rseed[g*H+i] = h.rseed[i];
            }
          for (size_t i = 0; i < nmant; i++) mant[g*nmant+i] = h.mant[i];
          for (size_t i = 0; i < nexpo; i++) expo[g*nexpo+i] = h.expo[i];
          for (size_t i = 0; i < 9; i++) scal[g*9+i] = h.scal[i];
        }

      sh->seed = complex_trafo ? UploadReal (device, FlatArray<Complex>(cseed), fp64)
                               : UploadReal (device, FlatArray<double>(rseed), fp64);
      sh->pow_mant = UploadReal (device, FlatArray<double>(mant), fp64);
      sh->pow_exp = UploadBuffer (device, FlatArray<int>(expo));
      sh->scal = UploadReal (device, FlatArray<double>(scal), fp64);

      // recursion coefficients: m <= min(os,ot), n <= os+ot+1 (Bmn(m-1,l+1) with l <= os+ot)
      auto tables = GetCoef2DTables (device, int(M), int(H)+1, fp64);
      sh->ctab = tables.buffer;
      sh->cstride = tables.stride;
      sh->cblock = (tables.max_m+1)*tables.stride;

      sh->threads_apply = GroupThreads (device, 128);

      // one library per device, real type and layer scalar type (REAL, Complex<REAL>)
      std::string source = ReplaceAll (code_shiftz_gpu, "TRAFO", complex_trafo ? "Complex<REAL>" : "REAL");
      auto kernels = GetKernels (device, RealSource (device, source, fp64),
                                 { "shiftz_boundaries_g", "shiftz_layers_g", "shiftz_apply_b" });
      sh->boundaries = kernels[0];
      sh->layers = kernels[1];
      sh->apply = kernels[2];
      sh->threads_boundaries = WidthThreads (device, *sh->boundaries) (H);
      sh->threads_layers = WidthThreads (device, *sh->layers) (W);
      return sh;
    }

    // bytes of a ShiftZWorkspace for units of G groups of (os, ot): K = 1 .. max(M,1) layer slots per group by the caps
    struct ShiftZNeeds
    {
      size_t layer0 = 0, boundary = 0, layers = 0;

      ShiftZNeeds () = default;
      ShiftZNeeds (int os, int ot, size_t G, bool complex_trafo, bool fp64, size_t budget, size_t fill)
      {
        size_t H = size_t(os)+ot+1, W = size_t(std::max(os,ot))+1, M = size_t(std::min(os,ot));
        size_t elem = EntryBytes (complex_trafo, fp64);
        size_t layer_bytes = H*W*elem;
        auto fit = [&] (size_t bytes) { return std::max (size_t(1), bytes / (G*layer_bytes)); };
        size_t slots = std::min (std::max (M, size_t(1)), std::max (fit (FMM_DEVICE_LANE_BUFFER_BYTES),
                                                                    std::min (fit (budget), std::max (size_t(1), fill / G))));
        layer0 = G*layer_bytes;
        boundary = G*(M+1)*(H+4*W)*elem;
        layers = G*slots*layer_bytes;
      }
      void Max (const ShiftZNeeds & b)
      {
        layer0 = std::max (layer0, b.layer0); boundary = std::max (boundary, b.boundary);
        layers = std::max (layers, b.layers);
      }
      bool operator== (const ShiftZNeeds & b) const
      { return layer0 == b.layer0 && boundary == b.boundary && layers == b.layers; }
    };

    bool ShiftZFits (const DeviceShiftZ & sh, const ShiftZWorkspace & ws)
    {
      size_t G = sh.groups;
      size_t H = size_t(sh.os)+sh.ot+1, W = size_t(std::max(sh.os,sh.ot))+1, M = size_t(std::min(sh.os,sh.ot));
      size_t elem = EntryBytes (sh.complex_trafo, sh.fp64);
      return ws.layer0 && ws.boundary && ws.layers
        && ws.layer0->Size() >= G*H*W*elem && ws.boundary->Size() >= G*(M+1)*(H+4*W)*elem
        && ws.layers->Size() >= G*H*W*elem;
    }

    // after the checks: layer 0 and the boundaries of every group; returns K, the layer slots per group of ws
    int EnqueueShiftZBoundaries (ngs_gpu::Queue & queue, const DeviceShiftZ & sh, const ShiftZWorkspace & ws)
    {
      using ngs_gpu::Dim3;
      int os = sh.os, ot = sh.ot, G = sh.groups;
      size_t H = size_t(os)+ot+1, W = size_t(std::max(os,ot))+1, M = size_t(std::min(os,ot));
      size_t layer_bytes = H*W*EntryBytes (sh.complex_trafo, sh.fp64);
      if (!ShiftZFits (sh, ws))
        throw Exception("fmm device ShiftZ: workspace too small for " + ToString(G) + " groups of os=" + ToString(os)
                        + ", ot=" + ToString(ot));

      queue.Launch (*sh.boundaries, Dim3(unsigned(G)), Dim3(sh.threads_boundaries),
                    { ws.layer0, ws.boundary, sh.seed, sh.ctab, sh.pow_mant, sh.pow_exp, sh.scal,
                      os, ot, int(sh.variant), sh.cstride, sh.cblock });

      // K layer slots per group; the kernels address slot k of group g at (g*K+k)*H*W: ints
      return int(std::min ({ std::max (M, size_t(1)), ws.layers->Size() / (G*layer_bytes), size_t(INT_MAX) / (G*H*W) }));
    }

    // layers m0 .. m0+nlayers-1 >= 1 of every group, rebuilt from their boundaries into slots 0 .. nlayers-1 of K
    void EnqueueShiftZLayers (ngs_gpu::Queue & queue, const DeviceShiftZ & sh, const ShiftZWorkspace & ws,
                              int m0, int nlayers, int K)
    {
      queue.Launch (*sh.layers, ngs_gpu::Dim3(unsigned(nlayers), unsigned(sh.groups)), ngs_gpu::Dim3(sh.threads_layers),
                    { ws.boundary, ws.layers, sh.ctab, sh.pow_mant, sh.pow_exp, sh.scal,
                      sh.os, sh.ot, int(sh.variant), sh.cstride, sh.cblock, m0, K });
    }

    // in and out: distinct, (os+1)^2 and (ot+1)^2 * lanes entries (the caller checks); launches only
    void EnqueueShiftZ (ngs_gpu::Queue & queue, const DeviceShiftZ & sh,
                        const ShiftZWorkspace & ws, const RealBuffer & in, const RealBuffer & out)
    {
      using ngs_gpu::Dim3;
      int os = sh.os, ot = sh.ot, dim = sh.lanes;
      size_t M = size_t(std::min(os,ot));
      int K = EnqueueShiftZBoundaries (queue, sh, ws);
      size_t slots = size_t(K);
      int variant = int(sh.variant);
      // layer 0 (in layer0), then chunks of K layers rebuilt from their boundaries, each applied over the lane blocks
      auto apply = [&] (int m0, int nlayers)
      {
        // the chunk's first layer has the most items per lane: ot-m0+1
        int parts = ApplyParts (size_t(ot-m0+1), sh.threads_apply, sh.nblocks);
        queue.Launch (*sh.apply, Dim3(unsigned(sh.nblocks)*unsigned(parts), unsigned(nlayers)), Dim3(sh.threads_apply),
                      { ws.layer0, ws.layers, sh.pow_mant, sh.pow_exp, sh.blocks, in, out,
                        os, ot, variant, m0, K, dim, parts });
      };
      apply (0, 1);
      for (size_t m0 = 1; m0 <= M; m0 += slots)
        {
          int nlayers = int(std::min (slots, M-m0+1));
          EnqueueShiftZLayers (queue, sh, ws, int(m0), nlayers, K);
          apply (int(m0), nlayers);
        }
    }
  }


  /* ******************************* Transform units ******************************* */

  namespace
  {
    /*
      The batched translation of the CPU FMM (ProcessDynamicBatchSS/SR/RR) for a unit of G >= 1 translation groups
      of one shape (variant, os, ot, kappa, both radii), each with its own len and theta, lanes group after group:
        work = RotateZFlip(phi_b, flip_b) of src_pool[src_idx[b]]       (transform_gather)
        work = RotateY(theta), ShiftZ(-len), RotateY(-theta)           (over the groups of the unit)
        tgt_pool[tgt_idx[b]] += RotateZFlip(-phi_b, flip_b) of work      (transform_scatter)
      A pool holds expansions one after another, expansion i at i*(order+1)^2. PrepareUnits and ReserveNeeds
      allocate and upload, so never inside a recording; EnqueueTransformUnit only launches. A fused S->R unit runs
      the five steps in one kernel (m2l_fused) on its tables, see EnqueueFusedUnit.
    */
    struct DeviceTransformUnit
    {
      bool fp64;
      int os, ot, lanes;                             // lanes = B
      int max_src_idx, max_tgt_idx;                  // for the pool size checks
      int phase_stride;                              // max(os,ot)+1
      RealBuffer phase;                              // [b*phase_stride+k] = e^{i k phi_b}
      ngs_gpu::TypedBuffer<int> flip, src_idx, tgt_idx;  // [b]
      int groups;                                    // G >= 1, lanes consecutive in group order
      std::shared_ptr<DeviceRotation> rotate, rotate_back;   // RotateY(theta_g), order os; RotateY(-theta_g), order ot; G groups each
      std::shared_ptr<DeviceShiftZ> shift;           // ShiftZ(-len_g) of every group, os -> ot, one unit
      unsigned threads;                              // work-items per group of gather/scatter
      std::shared_ptr<ngs_gpu::Kernel> gather, scatter;
      // a fused unit: its lane blocks of up to FUSED_LANES lanes, the kernels and their work-items per group
      std::shared_ptr<ngs_gpu::Kernel> fold, m2l;
      LaneBlocks fused_blocks;
      unsigned fused_threads = 1;
    };

    struct TransformWorkspace
    {
      // (max(os,ot)+1)^2 * B entries each, for every unit it serves; a fused unit's folded rotations instead
      RealBuffer work0, work1;
      RealBuffer rotate;                             // at least G*RotateYWorkStride(max(os,ot)) reals for every unit
      ShiftZWorkspace shift;
    };

    // a unit's piece of one group: the lanes [first, first+lanes) of the lane arrays
    struct UnitPiece { size_t group, first, lanes; };

    // prev: the unit prepared before this one; if both consist of one piece of the same group (a group cut by the
    // lane cap), this one reuses its packed rotation tables
    std::shared_ptr<DeviceTransformUnit>
    PrepareUnit (const std::shared_ptr<ngs_gpu::Device> & device, FMMTranslation variant, int os, int ot,
                 FlatArray<UnitPiece> pieces, const std::vector<ShiftZGroupHost> & hosts,
                 const std::vector<RotationHost> & rotate, const std::vector<RotationHost> & rotate_back,
                 FlatArray<double> phi, FlatArray<int> flip, FlatArray<int> src_idx, FlatArray<int> tgt_idx,
                 const DeviceTransformUnit * prev, size_t prev_group, bool fused, bool fp64)
    {
      size_t first = pieces[0].first, nb = 0;
      for (auto & p : pieces) nb += p.lanes;
      size_t rows_s = size_t(os+1)*(os+1), rows_t = size_t(ot+1)*(ot+1);
      size_t rows = std::max (rows_s, rows_t);
      // the kernels index the unit's columns (ii*B+b) and the (2n+1)*B items of a degree with ints
      if (rows*nb > INT_MAX || (2*size_t(std::max(os,ot))+1)*nb > INT_MAX)
        throw Exception("FMMDevicePlan: (max(os,ot)+1)^2 * B must fit an int");
      int B = int(nb);

      // indices; the kernels address pool entries with ints, the scatter as 2*index on Metal (complex as two reals)
      int max_src = -1, max_tgt = -1;
      for (size_t b = first; b < first+nb; b++)
        {
          if (src_idx[b] < 0 || tgt_idx[b] < 0)
            throw Exception("FMMDevicePlan: negative source or target index at element " + ToString(b));
          max_src = std::max (max_src, src_idx[b]);
          max_tgt = std::max (max_tgt, tgt_idx[b]);
        }
      if ((size_t(max_src)+1)*rows_s > INT_MAX || 2*(size_t(max_tgt)+1)*rows_t > INT_MAX)
        throw Exception("FMMDevicePlan: pool offsets exceed the int indexing of the kernels "
                        "((max_src_idx+1)*(os+1)^2 and 2*(max_tgt_idx+1)*(ot+1)^2 must fit an int)");

      auto unit = std::make_shared<DeviceTransformUnit>();
      unit->fp64 = fp64;
      unit->os = os; unit->ot = ot; unit->lanes = B;
      unit->max_src_idx = max_src; unit->max_tgt_idx = max_tgt;

      // RotateY(theta_g), ShiftZ(-len_g) and RotateY(-theta_g) of all groups, on one lane-block list
      size_t G = pieces.Size();
      Array<const ShiftZGroupHost*> shift_hosts(G);
      Array<const RotationHost*> rot_hosts(G), back_hosts(G);
      Array<int> lanes(G);
      for (size_t i = 0; i < G; i++)
        {
          auto & p = pieces[i];
          shift_hosts[i] = &hosts[p.group];
          rot_hosts[i] = &rotate[p.group];
          back_hosts[i] = &rotate_back[p.group];
          lanes[i] = int(p.lanes);
        }
      bool share = prev && G == 1 && prev->groups == 1 && prev_group == pieces[0].group;
      auto blocks = MakeLaneBlocks (device, lanes, fused ? FUSED_LANES : DEVICE_LANE_BLOCK);
      unit->groups = int(G);
      unit->rotate = PackRotation (device, os, rot_hosts, blocks, share ? prev->rotate.get() : nullptr, fp64);
      unit->rotate_back = PackRotation (device, ot, back_hosts, blocks, share ? prev->rotate_back.get() : nullptr, fp64);
      unit->shift = PackShiftZ (device, variant, os, ot, shift_hosts, blocks, fp64);

      // e^{i k phi}, k <= max(os,ot), by the recurrence of RotateZFlip in double, rounded once
      int K = std::max(os,ot)+1;
      unit->phase_stride = K;
      Array<Complex> phase(size_t(B)*K);
      for (int b = 0; b < B; b++)
        {
          Complex e(cos(phi[first+b]), sin(phi[first+b]));
          Complex prod = 1.0;
          for (int k = 0; k < K; k++)
            {
              phase[size_t(b)*K+k] = prod;
              prod *= e;
            }
        }
      unit->phase = UploadReal (device, FlatArray<Complex>(phase), fp64);
      unit->flip = UploadBuffer (device, flip.Range(first, first+B));
      unit->src_idx = UploadBuffer (device, src_idx.Range(first, first+B));
      unit->tgt_idx = UploadBuffer (device, tgt_idx.Range(first, first+B));

      // gather/scatter work on tiles in group memory (bem/transform_gpu_kernel.h): 256 work-items per tile group, or what
      // the compiled kernels allow if a backend reports less (their loops stride over the tile, any count works)
      auto kernels = GetTransformKernels (device, fp64);
      unit->gather = kernels[0];
      unit->scatter = kernels[1];
      unit->threads = GroupThreads (device, 256);
      for (auto * k : { unit->gather.get(), unit->scatter.get() })
        unit->threads = KernelThreads (*k, unit->threads);
      if (fused)
        {
          auto fk = GetFusedKernels (device, std::max (os, ot), fp64);
          unit->fold = fk[0];
          unit->m2l = fk[1];
          unit->fused_blocks = blocks;
          unit->fused_threads = KernelThreads (*unit->m2l, GroupThreads (device, 256));
        }
      return unit;
    }

    // units of whole lane ranges: at most FMMDeviceLaneCap lanes and ShiftZGroupCap groups each; fused S->R units at
    // most FusedGroupCap groups, their lanes bounded only by the int indexing (they need no lane buffers)
    std::vector<std::shared_ptr<DeviceTransformUnit>>
    PrepareUnits (std::shared_ptr<ngs_gpu::Device> device, FMMTranslation variant,
                  std::variant<double,std::complex<double>> kappa,
                  double source_rtyp, double target_rtyp, int os, int ot,
                  FlatArray<double> len, FlatArray<double> theta, FlatArray<int> group_lanes,
                  FlatArray<double> phi, FlatArray<int> flip, FlatArray<int> src_idx, FlatArray<int> tgt_idx, bool fused,
                  bool fp64)
    {
      if (os < 1 || ot < 1)
        throw Exception("FMMDevicePlan: orders must be >= 1 (the rotations need it), got os=" + ToString(os)
                        + ", ot=" + ToString(ot));
      CheckShiftZOrders (os, ot);   // the int bounds of the ShiftZ kernels and GetCoef2DTables
      size_t G = len.Size(), nb = phi.Size();
      for (size_t b = 0; b < nb; b++)
        if (!std::isfinite(phi[b]) || (flip[b] != 0 && flip[b] != 1))
          throw Exception("FMMDevicePlan: element " + ToString(b) + " needs a finite phi and flip 0 or 1");
      // the caps (FMM_DEVICE_FILL); a lane block or a group is one GPU group of the applies or the ShiftZ boundaries
      int order = std::max (os, ot);
      size_t budget = FMMDeviceBudget (device), fill = FMM_DEVICE_FILL*device->ComputeUnits();
      bool complex_trafo = variant == FMMTranslation::SR || std::holds_alternative<std::complex<double>>(kappa);
      auto lanes = [&] (size_t bytes) { return FMMDeviceLaneCap (order, fp64, bytes); };
      auto groups = [&] (size_t bytes) { return ShiftZGroupCap (os, ot, complex_trafo, fp64, bytes); };
      size_t lane_cap = fused ? size_t(INT_MAX) / sqr (size_t(order)+1)
        : std::max (lanes (FMM_DEVICE_LANE_BUFFER_BYTES), std::min (lanes (budget), fill*DEVICE_LANE_BLOCK));
      size_t group_cap = fused ? FusedGroupCap (os, ot, fp64)
        : std::max (groups (FMM_DEVICE_LANE_BUFFER_BYTES), std::min (groups (budget), fill));

      // per group once, on the host: the rotations as the CPU batch calls them, RotateY(theta) and RotateY(-theta),
      // and ShiftZ(-len); packed and uploaded per unit
      std::vector<ShiftZGroupHost> hosts;
      std::vector<RotationHost> rotate, rotate_back;
      for (size_t g = 0; g < G; g++)
        {
          rotate.push_back (ComputeRotation (theta[g], os));
          rotate_back.push_back (ComputeRotation (-theta[g], ot));
          hosts.push_back (ComputeShiftZGroup (variant, kappa, -len[g], source_rtyp, target_rtyp, os, ot, fp64));
        }

      std::vector<std::shared_ptr<DeviceTransformUnit>> units;
      Array<UnitPiece> pieces;
      size_t unit_lanes = 0, lane = 0, prev_group = 0;
      auto flush = [&] ()
      {
        if (pieces.Size() == 0) return;
        units.push_back (PrepareUnit (device, variant, os, ot, pieces, hosts, rotate, rotate_back,
                                      phi, flip, src_idx, tgt_idx, units.empty() ? nullptr : units.back().get(), prev_group,
                                      fused, fp64));
        prev_group = pieces[0].group;
        pieces.SetSize0();
        unit_lanes = 0;
      };
      for (size_t g = 0; g < G; g++)
        for (size_t done = 0; done < size_t(group_lanes[g]); )
          {
            if (unit_lanes == lane_cap || pieces.Size() == group_cap)
              flush();
            size_t take = std::min (size_t(group_lanes[g])-done, lane_cap-unit_lanes);
            pieces.Append (UnitPiece{ g, lane, take });
            unit_lanes += take; lane += take; done += take;
          }
      flush();
      return units;
    }

    // entries of each lane buffer for the unit: its lanes' rows, or a fused unit's folded rotations (as reals)
    size_t WorkNeeds (const DeviceTransformUnit & unit)
    {
      if (unit.m2l)
        return (size_t(unit.groups)*std::max (FoldReals (unit.os), FoldReals (unit.ot))+1) / 2;
      return sqr(size_t(std::max (unit.os, unit.ot))+1)*size_t(unit.lanes);
    }

    // what units need of a workspace; the maximum over units is what a shared workspace needs
    struct TransformNeeds
    {
      size_t rotate = 0;             // RotateY scratch reals; 0: no unit
      size_t work = 0;               // entries of each lane buffer
      ShiftZNeeds shift;

      void Add (const DeviceTransformUnit & unit, size_t budget, size_t fill)
      {
        rotate = std::max ({ rotate, RotateYNeeds (unit.os, size_t(unit.groups), unit.fp64, budget, fill),
                             RotateYNeeds (unit.ot, size_t(unit.groups), unit.fp64, budget, fill) });
        work = std::max (work, WorkNeeds (unit));
        shift.Max (ShiftZNeeds (unit.os, unit.ot, size_t(unit.groups), unit.shift->complex_trafo, unit.fp64, budget, fill));
      }
      void Max (const TransformNeeds & b)
      {
        rotate = std::max (rotate, b.rotate);
        work = std::max (work, b.work);
        shift.Max (b.shift);
      }
      bool operator== (const TransformNeeds & b) const
      { return rotate == b.rotate && work == b.work && shift == b.shift; }
      bool Empty () const { return rotate == 0; }
    };

    std::shared_ptr<TransformWorkspace> ReserveNeeds (const std::shared_ptr<ngs_gpu::Device> & device, const TransformNeeds & needs,
                                                      bool fp64)
    {
      auto ws = std::make_shared<TransformWorkspace>();
      ws->rotate = NewRealBuffer (device, needs.rotate, fp64, false);
      ws->shift = { device->NewBuffer (needs.shift.layer0), device->NewBuffer (needs.shift.boundary),
                    device->NewBuffer (needs.shift.layers) };
      ws->work0 = NewRealBuffer (device, needs.work, fp64, true);
      ws->work1 = NewRealBuffer (device, needs.work, fp64, true);
      return ws;
    }

    /*
      A fused unit: the folded rotation blocks of every group into the lane buffers, which it doesn't use (rot_fold
      after each construct chunk), all its ShiftZ layers into the layer slots (the group cap keeps them within the
      workspace), then m2l_fused over its lane blocks. Only the last launch writes the targets.
    */
    void EnqueueFusedUnit (ngs_gpu::Queue & queue, const DeviceTransformUnit & unit, const TransformWorkspace & ws,
                           const RealBuffer & src_pool, const RealBuffer & tgt_pool)
    {
      auto & sh = *unit.shift;
      RealBuffer fold { ws.work0.buf, unit.fp64, false }, fold_back { ws.work1.buf, unit.fp64, false };   // as reals
      EnqueueRotateYFold (queue, *unit.rotate, ws.rotate, *unit.fold, fold);
      EnqueueRotateYFold (queue, *unit.rotate_back, ws.rotate, *unit.fold, fold_back);
      int K = EnqueueShiftZBoundaries (queue, sh, ws.shift);
      int M = std::min (unit.os, unit.ot);
      if (K < M)
        throw Exception("fmm device fused M2L: the workspace holds " + ToString(K) + " of the " + ToString(M) + " layers");
      EnqueueShiftZLayers (queue, sh, ws.shift, 1, M, K);
      queue.Launch (*unit.m2l, ngs_gpu::Dim3(unsigned(unit.fused_blocks.count)), ngs_gpu::Dim3(unit.fused_threads),
                    { unit.fused_blocks.buffer, src_pool, unit.src_idx, unit.phase, unit.flip, fold, fold_back,
                      ws.shift.layer0, ws.shift.layers, sh.pow_mant, sh.pow_exp, unit.tgt_idx, tgt_pool,
                      unit.os, unit.ot, unit.phase_stride, K });
    }

    void EnqueueTransformUnit (ngs_gpu::Queue & queue, const DeviceTransformUnit & unit,
                               const TransformWorkspace & ws, const RealBuffer & src_pool, const RealBuffer & tgt_pool)
    {
      using ngs_gpu::Dim3;
      int os = unit.os, ot = unit.ot, B = unit.lanes;
      size_t rows_s = size_t(os+1)*(os+1), rows_t = size_t(ot+1)*(ot+1);
      // every check (EnqueueRotateY/ShiftZ included) precedes the scatter: a failing call never touches the targets
      if (!ws.rotate || !ws.work0 || !ws.work1 || ws.work0.Size() < WorkNeeds (unit) || ws.work1.Size() < WorkNeeds (unit))
        throw Exception("fmm device transform: workspace too small for os=" + ToString(os) + ", ot=" + ToString(ot)
                        + ", B=" + ToString(B) + ", G=" + ToString(unit.groups));
      if (!src_pool || !tgt_pool || src_pool.Size() < (size_t(unit.max_src_idx)+1)*rows_s
          || tgt_pool.Size() < (size_t(unit.max_tgt_idx)+1)*rows_t)
        throw Exception("fmm device transform: source or target pool smaller than the unit's largest index");
      if (unit.m2l)
        return EnqueueFusedUnit (queue, unit, ws, src_pool, tgt_pool);

      // gather/scatter: one group per tile of TRANSFORM_TILE rows (x) x TRANSFORM_TILE lanes (y); the lane tiles are the
      // y grid dimension (65535 on CUDA)
      auto tiles = [] (size_t n) { return unsigned((n+TRANSFORM_TILE-1) / TRANSFORM_TILE); };
      if (tiles(size_t(B)) > 65535)
        throw Exception("fmm device transform: B=" + ToString(B) + " lanes exceed the grid of the gather/scatter tiles (65535*"
                        + ToString(TRANSFORM_TILE) + ")");

      queue.Launch (*unit.gather, Dim3(tiles(rows_s), tiles(size_t(B))), Dim3(unit.threads),
                    { src_pool, unit.src_idx, unit.phase, unit.flip, ws.work0, os, B, unit.phase_stride });
      EnqueueRotateY (queue, *unit.rotate, ws.rotate, ws.work0, ws.work1);
      EnqueueShiftZ (queue, *unit.shift, ws.shift, ws.work1, ws.work0);
      EnqueueRotateY (queue, *unit.rotate_back, ws.rotate, ws.work0, ws.work1);
      queue.Launch (*unit.scatter, Dim3(tiles(rows_t), tiles(size_t(B))), Dim3(unit.threads),
                    { ws.work1, unit.tgt_idx, unit.phase, unit.flip, tgt_pool, ot, B, unit.phase_stride });
    }

    bool GpuRegistered ()
    {
      return ngs_gpu::GetDevice() && ngs_gpu::GetDevice() != ngs_gpu::GetCpuDevice();
    }
  }

  FMM_Device ResolveFMMDevice (const std::string & device, bool fp32, bool use_fmm)
  {
    if (device != "auto" && device != "gpu" && device != "cpu" && device != "none")
      throw Exception("fmm device='" + device + "': valid values are 'auto', 'gpu', 'cpu', 'none'");

    if (!use_fmm || device == "none")
      return FMM_Device::None;

    // "cpu" never touches GetDevice(): it needs no GPU at all, and the CPU device has float64
    if (device == "cpu")
      return FMM_Device::Cpu;

    if (device == "gpu")
      {
        // GetDevice() lazily runs the device creator (e.g. cuInit): "gpu" propagates its failure
        if (!GpuRegistered())
          throw Exception("fmm device='gpu': no GPU device registered (import ngsolve.gpu with a Metal or CUDA build)");
        if (!fp32 && !ngs_gpu::GetDevice()->HasFloat64())
          throw Exception("fmm device='gpu': " + ngs_gpu::GetDevice()->Name()
                          + " has no float64 for the fp64 FMM kernels, pass fp32=True");
        return FMM_Device::Gpu;
      }

    // "auto": a failing device creator (e.g. cuInit without a driver) means no GPU
    bool usable = false;
    try
      {
        usable = GpuRegistered() && (fp32 || ngs_gpu::GetDevice()->HasFloat64());
      }
    catch (...) { }
    return usable ? FMM_Device::Gpu : FMM_Device::None;
  }


  /* *********************** FMM passes on the device (mptools_device.hpp) *********************** */

  namespace
  {
    // one device call at a time: the transfers and the launches all use the device's default queue
    std::mutex fmm_device_mutex;

    Timer<> & AllocTimer () { static Timer t("fmm device alloc"); return t; }
    Timer<> & PrepareTimer () { static Timer t("fmm device prepare"); return t; }
    Timer<> & ExecuteTimer () { static Timer t("fmm device execute"); return t; }
    Timer<> & RecordTimer () { static Timer t("fmm device record"); return t; }
    Timer<> & UploadTimer () { static Timer t("fmm device upload"); return t; }
    Timer<> & DownloadTimer () { static Timer t("fmm device download"); return t; }
    Timer<> & MatvecTimer () { static Timer t("fmm device matvec"); return t; }

    std::shared_ptr<ngs_gpu::Device> GetFMMDevice (FMM_Device device)
    {
      if (device == FMM_Device::Cpu)
        return ngs_gpu::GetCpuDevice();
      if (device == FMM_Device::Gpu)
        {
          if (!GpuRegistered())
            throw Exception("fmm device='gpu': no GPU device registered (import ngsolve.gpu with a Metal or CUDA build)");
          return ngs_gpu::GetDevice();
        }
      throw Exception("fmm device execution without a device (FMM_Device::None)");
    }
  }

  struct FMMDeviceTree::Impl
  {
    struct Pool
    {
      int order;
      RealBuffer buffer;   // slots * (order+1)^2 complex entries
    };
    std::shared_ptr<ngs_gpu::Device> device;
    bool fp64 = false;     // double pools (Complex host expansions), else float (Complex32)
    std::vector<Pool> pools;

    const Pool & Get (int pool) const
    {
      if (pool < 0 || size_t(pool) >= pools.size())
        throw Exception("FMMDeviceTree: no pool " + ToString(pool));
      return pools[pool];
    }

    template <typename C>
    void Upload (int ipool, size_t first, size_t n, const std::function<void(C*, size_t, size_t)> & fill)
    {
      RegionTimer reg(UploadTimer());
      std::lock_guard<std::mutex> lock(fmm_device_mutex);
      auto & pool = Get (ipool);
      if (pool.buffer.fp64 != std::is_same_v<C,Complex>)
        throw Exception("FMMDeviceTree::Upload: the pool holds the other real type");
      if (first > pool.buffer.Size() || n > pool.buffer.Size()-first)
        throw Exception("FMMDeviceTree::Upload: range exceeds the pool");
      if (n == 0) return;
      FillBuffer<C> (pool.buffer, first, n, [&] (C * dst, size_t off, size_t cnt) { fill (dst, first+off, cnt); });
    }

    template <typename C>
    void Download (int ipool, size_t n, const std::function<void(const C*, size_t, size_t)> & consume)
    {
      RegionTimer reg(DownloadTimer());
      std::lock_guard<std::mutex> lock(fmm_device_mutex);
      auto & pool = Get (ipool);
      if (pool.buffer.fp64 != std::is_same_v<C,Complex>)
        throw Exception("FMMDeviceTree::Download: the pool holds the other real type");
      if (n > pool.buffer.Size())
        throw Exception("FMMDeviceTree::Download: range exceeds the pool");
      if (n == 0) return;
      ReadBuffer<C> (pool.buffer, n, [&] (const C * src) { consume (src, 0, n); });
    }
  };

  FMMDeviceTree :: FMMDeviceTree (FMM_Device device, bool fp64)
    : impl(std::make_unique<Impl>())
  {
    impl->device = GetFMMDevice (device);
    impl->fp64 = fp64;
  }

  FMMDeviceTree :: ~FMMDeviceTree () = default;

  int FMMDeviceTree :: AddPool (int order, size_t slots)
  {
    RegionTimer reg(AllocTimer());
    if (order < 0 || slots < 1 || slots > size_t(INT_MAX))
      throw Exception("FMMDeviceTree::AddPool: order must be >= 0 and 1 <= slots <= INT_MAX");
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    Impl::Pool pool { order, NewRealBuffer (impl->device, slots * sqr(size_t(order)+1), impl->fp64, true) };
    impl->pools.push_back (pool);
    return int(impl->pools.size())-1;
  }

  void FMMDeviceTree :: Upload (int ipool, size_t first, size_t n, const std::function<void(Complex32*, size_t, size_t)> & fill)
  {
    impl->Upload (ipool, first, n, fill);
  }

  void FMMDeviceTree :: Upload (int ipool, size_t first, size_t n, const std::function<void(Complex*, size_t, size_t)> & fill)
  {
    impl->Upload (ipool, first, n, fill);
  }

  void FMMDeviceTree :: Download (int ipool, size_t n, const std::function<void(const Complex32*, size_t, size_t)> & consume)
  {
    impl->Download (ipool, n, consume);
  }

  void FMMDeviceTree :: Download (int ipool, size_t n, const std::function<void(const Complex*, size_t, size_t)> & consume)
  {
    impl->Download (ipool, n, consume);
  }

  namespace
  {
    // the entries [first[e], first[e]+number[e]) within items, and the slots of p (comps per node) they use
    void CheckLevel (const FMMDeviceTree::Impl::Pool & p, int comps, size_t count, FlatArray<int> first,
                     FlatArray<int> number, FlatArray<int> slot, size_t items, const char * where)
    {
      size_t rows = sqr(size_t(p.order)+1);
      if (p.order > FMMDeviceP2M::max_order || count == 0 || count > size_t(INT_MAX) || first.Size() != number.Size()
          || slot.Size() < first.Size() || p.buffer.Size() > size_t(INT_MAX))
        throw Exception(std::string(where) + ": invalid level");
      for (size_t e = 0; e < first.Size(); e++)
        if (first[e] < 0 || number[e] < 1 || size_t(first[e])+number[e] > items)
          throw Exception(std::string(where) + ": an entry exceeds its range");
      for (int sl : slot)
        if (sl < 0 || (size_t(sl)+1)*comps*rows > p.buffer.Size())
          throw Exception(std::string(where) + ": a slot exceeds the pool");
    }
  }


  struct FMMDeviceP2M::Impl
  {
    struct Level
    {
      RealBuffer pool;
      int order, leaves;
      ngs_gpu::KernelArg scale, kappa_over_s, kappa_times_s;
      ngs_gpu::TypedBuffer<int> first, count, slot;
    };
    std::shared_ptr<ngs_gpu::Device> device;
    bool fp64 = false;            // the tree's real type
    LeafKernels k;
    int kind, comps;
    LeafKappa kappa;
    Array<int> order;             // source j is point order[j]
    RealBuffer geo, dir, xyz;
    RealBuffer charge, dipole, dipvec;
    std::vector<Level> levels;
    std::map<std::pair<int,int>, std::pair<int,int>> sources;   // (pool, slot) -> the leaf's sources

    // launches only (recordable): one per level, a group per (leaf, component)
    void Enqueue (ngs_gpu::Queue & queue) const
    {
      for (auto & lev : levels)
        queue.Launch (*k.p2m, ngs_gpu::Dim3(lev.leaves, comps), ngs_gpu::Dim3(k.p2m_threads),
                      { lev.first, lev.count, lev.slot, geo, dir, charge, dipole, dipvec, k.tables, lev.pool,
                        lev.order, lev.scale, kappa.Arg (kappa.kappa, fp64), lev.kappa_over_s,
                        lev.kappa_times_s, kappa.Arg (kappa.khat, fp64), kind, comps });
    }
  };

  FMMDeviceP2M :: FMMDeviceP2M (const FMMDeviceTree & tree, int kind, int comps,
                                std::variant<double,std::complex<double>> kappa,
                                FlatArray<int> order, FlatArray<double> geo, FlatArray<double> dir, FlatArray<double> xyz)
    : impl(std::make_unique<Impl>())
  {
    RegionTimer reg(AllocTimer());
    size_t nsrc = order.Size();
    if (kind < 1 || kind > 5 || comps < 1 || comps > 65535 || nsrc == 0 || geo.Size() != 5*nsrc || dir.Size() != 3*nsrc
        || xyz.Size() != 3*nsrc || 6*nsrc > size_t(INT_MAX) || 3*nsrc*comps > size_t(INT_MAX))
      throw Exception("FMMDeviceP2M: invalid sources");
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    auto & im = *impl;
    im.device = tree.impl->device;
    im.fp64 = tree.impl->fp64;
    im.kappa = kappa;
    im.k = GetLeafKernels (im.device, im.kappa.complex_kappa, im.fp64);
    im.kind = kind; im.comps = comps;
    im.order = order;
    im.geo = UploadReal (im.device, geo, im.fp64);
    im.dir = UploadReal (im.device, (kind & 2) ? dir : FlatArray<double>(), im.fp64);
    im.xyz = UploadHiLo (im.device, xyz, im.fp64);
    im.charge = NewRealBuffer (im.device, (kind & 1) ? nsrc*comps : 1, im.fp64, true);
    im.dipole = NewRealBuffer (im.device, (kind & 2) ? nsrc*comps : 1, im.fp64, true);
    im.dipvec = NewRealBuffer (im.device, (kind & 4) ? 3*nsrc*comps : 1, im.fp64, true);
  }

  FMMDeviceP2M :: ~FMMDeviceP2M () = default;

  void FMMDeviceP2M :: AddLevel (const FMMDeviceTree & tree, int pool, double scale,
                                 FlatArray<int> first, FlatArray<int> count, FlatArray<int> slot)
  {
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    auto & im = *impl;
    auto & p = tree.impl->Get (pool);
    size_t leaves = first.Size();
    if (tree.impl->device != im.device || tree.impl->fp64 != im.fp64 || slot.Size() != leaves)
      throw Exception("FMMDeviceP2M::AddLevel: invalid level");
    CheckLevel (p, im.comps, leaves, first, count, slot, im.order.Size(), "FMMDeviceP2M::AddLevel");
    im.levels.push_back ({ p.buffer, p.order, int(leaves), RealArg (scale, im.fp64), im.kappa.Arg (im.kappa.kappa/scale, im.fp64),
                           im.kappa.Arg (im.kappa.kappa*scale, im.fp64), UploadBuffer (im.device, first),
                           UploadBuffer (im.device, count), UploadBuffer (im.device, slot) });
    for (size_t k = 0; k < leaves; k++)
      im.sources[{ pool, slot[k] }] = { first[k], count[k] };
  }

  std::pair<int,int> FMMDeviceP2M :: Sources (int pool, int slot) const
  {
    auto it = impl->sources.find ({ pool, slot });
    if (it == impl->sources.end())
      throw Exception("FMMDeviceP2M::Sources: no device leaf in this slot");
    return it->second;
  }

  void FMMDeviceP2M :: Upload (const FMMDeviceSources & sources)
  {
    RegionTimer reg(UploadTimer());
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    auto & im = *impl;
    size_t comps = im.comps, n = im.order.Size()*comps;
    // w values per source and component, in the tree's real type
    auto fill = [&] (FlatArray<Complex> values, const RealBuffer & buf, size_t w)
    {
      if (values.Size() != sources.pts.Size()*comps*w)
        throw Exception("FMMDeviceP2M::Upload: the sources do not match the plan");
      auto copy = [&] (auto * dst, size_t first, size_t cnt)
      {
        for (size_t e = first; e < first+cnt; e++)
          dst[e-first] = values[size_t(im.order[e/(comps*w)])*comps*w + e%(comps*w)];
      };
      if (im.fp64)
        FillBuffer<Complex> (buf, 0, n*w, copy);
      else
        FillBuffer<Complex32> (buf, 0, n*w, copy);
    };
    if (sources.kind != im.kind || size_t(sources.comps) != comps)
      throw Exception("FMMDeviceP2M::Upload: the sources do not match the plan");
    if (im.kind & 1) fill (sources.charge, im.charge, 1);
    if (im.kind & 2) fill (sources.dipole, im.dipole, 1);
    if (im.kind & 4) fill (sources.dipvec, im.dipvec, 3);
  }


  struct FMMDeviceL2P::Impl
  {
    struct Level
    {
      RealBuffer pool;
      int order, count;           // count: leaves (L2P) or target entries (M2P)
      ngs_gpu::KernelArg scale;
      ngs_gpu::KernelArg ks, kis;   // kappa times and over the expansion's Scale(), for the gradient
      ngs_gpu::TypedBuffer<int> target, first, number, slot;
      RealBuffer geo;
    };
    struct NearField
    {
      int entries, kind, branch;
      ngs_gpu::KernelArg kappa;
      ngs_gpu::TypedBuffer<int> target, list, lcount, ranges;
      RealBuffer sxyz, sdir;
      RealBuffer charge, dipole, dipvec;
    };
    std::shared_ptr<ngs_gpu::Device> device;
    bool fp64 = false;            // the tree's real type
    LeafKernels k;
    int comps, mode, outs;        // outs: the entries per target and component (FMMDeviceTargets)
    LeafKappa kappa;
    Array<int> order;             // target j is point order[j]
    RealBuffer geo, xyz;
    RealBuffer out;
    std::vector<Level> l2p, m2p;
    std::vector<NearField> p2p;

    // launches only (recordable): L2P writes every target, then M2P adds, a pool after the other, then P2P (the
    // gradient near field from the charges, as the host's EvaluateDeriv)
    void Enqueue (ngs_gpu::Queue & queue) const
    {
      using ngs_gpu::Dim3;
      unsigned t = k.eval_threads;
      auto khat = kappa.Arg (kappa.khat, fp64);
      for (auto & lev : l2p)
        queue.Launch (*k.l2p, Dim3(lev.count, comps), Dim3(t),
                      { lev.first, lev.number, lev.slot, geo, k.tables, lev.pool, out, lev.order, lev.scale,
                        khat, lev.ks, lev.kis, mode, comps });
      for (auto & lev : m2p)
        queue.Launch (*k.m2p, Dim3((lev.count+t-1)/t, comps), Dim3(t),
                      { lev.target, lev.first, lev.number, lev.slot, lev.geo, k.tables, lev.pool, out,
                        lev.order, lev.scale, khat, lev.ks, lev.kis, mode, comps, lev.count });
      for (auto & nf : p2p)
        {
          if (mode & 1)
            queue.Launch (*k.p2p, Dim3((nf.entries+t-1)/t, comps), Dim3(t),
                          { nf.target, nf.list, nf.lcount, nf.ranges, xyz, nf.sxyz, nf.sdir, nf.charge, nf.dipole, nf.dipvec,
                            out, nf.kappa, nf.kind, nf.branch, comps, nf.entries, outs });
          if ((mode & 2) && (nf.kind & 1))
            queue.Launch (*k.p2p_grad, Dim3((nf.entries+t-1)/t, comps), Dim3(t),
                          { nf.target, nf.list, nf.lcount, nf.ranges, xyz, nf.sxyz, nf.charge, out,
                            nf.kappa, comps, nf.entries, outs, mode & 1 });
        }
    }
  };

  FMMDeviceL2P :: FMMDeviceL2P (const FMMDeviceTree & tree, int comps, int mode,
                                std::variant<double,std::complex<double>> kappa,
                                FlatArray<int> order, FlatArray<double> geo, FlatArray<double> xyz)
    : impl(std::make_unique<Impl>())
  {
    RegionTimer reg(AllocTimer());
    size_t ntg = order.Size(), outs = FMMDeviceTargetOuts (mode);
    if (comps < 1 || comps > 65535 || outs == 0 || ntg == 0 || geo.Size() != 5*ntg || xyz.Size() != 3*ntg
        || 6*ntg > size_t(INT_MAX) || ntg*comps*outs > size_t(INT_MAX))
      throw Exception("FMMDeviceL2P: invalid targets");
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    auto & im = *impl;
    im.device = tree.impl->device;
    im.fp64 = tree.impl->fp64;
    im.kappa = kappa;
    im.k = GetLeafKernels (im.device, im.kappa.complex_kappa, im.fp64);
    im.comps = comps; im.mode = mode; im.outs = int(outs);
    im.order = order;
    im.geo = UploadReal (im.device, geo, im.fp64);
    im.xyz = UploadHiLo (im.device, xyz, im.fp64);
    im.out = NewRealBuffer (im.device, ntg*comps*outs, im.fp64, true);
  }

  FMMDeviceL2P :: ~FMMDeviceL2P () = default;

  void FMMDeviceL2P :: AddLevel (const FMMDeviceTree & tree, int pool, double scale,
                                 FlatArray<int> first, FlatArray<int> count, FlatArray<int> slot)
  {
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    auto & im = *impl;
    auto & p = tree.impl->Get (pool);
    if (tree.impl->device != im.device || tree.impl->fp64 != im.fp64 || slot.Size() != first.Size())
      throw Exception("FMMDeviceL2P::AddLevel: invalid level");
    CheckLevel (p, im.comps, first.Size(), first, count, slot, im.order.Size(), "FMMDeviceL2P::AddLevel");
    // the local expansion's Scale() is 1/scale
    im.l2p.push_back ({ p.buffer, p.order, int(first.Size()), RealArg (scale, im.fp64),
                        im.kappa.Arg (im.kappa.kappa/scale, im.fp64), im.kappa.Arg (im.kappa.kappa*scale, im.fp64), {},
                        UploadBuffer (im.device, first), UploadBuffer (im.device, count), UploadBuffer (im.device, slot), {} });
  }

  void FMMDeviceL2P :: AddM2P (const FMMDeviceTree & tree, int pool, double scale, FlatArray<int> target,
                               FlatArray<int> first, FlatArray<int> count, FlatArray<int> slot, FlatArray<double> geo)
  {
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    auto & im = *impl;
    auto & p = tree.impl->Get (pool);
    if (tree.impl->device != im.device || tree.impl->fp64 != im.fp64 || target.Size() != first.Size()
        || geo.Size() != 5*slot.Size() || geo.Size() > size_t(INT_MAX))
      throw Exception("FMMDeviceL2P::AddM2P: invalid pairs");
    for (int j : target)
      if (j < 0 || size_t(j) >= im.order.Size())
        throw Exception("FMMDeviceL2P::AddM2P: a target is out of range");
    CheckLevel (p, im.comps, target.Size(), first, count, slot, slot.Size(), "FMMDeviceL2P::AddM2P");
    im.m2p.push_back ({ p.buffer, p.order, int(target.Size()), RealArg (scale, im.fp64),
                        im.kappa.Arg (im.kappa.kappa*scale, im.fp64), im.kappa.Arg (im.kappa.kappa/scale, im.fp64),
                        UploadBuffer (im.device, target), UploadBuffer (im.device, first), UploadBuffer (im.device, count),
                        UploadBuffer (im.device, slot), UploadReal (im.device, geo, im.fp64) });
  }

  void FMMDeviceL2P :: AddP2P (const FMMDeviceP2M & p2m, FlatArray<int> target, FlatArray<int> list, FlatArray<int> lcount,
                               FlatArray<int> ranges)
  {
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    auto & im = *impl;
    auto & src = *p2m.impl;
    size_t entries = target.Size(), nsrc = src.order.Size(), nranges = ranges.Size()/2;
    if (src.device != im.device || src.fp64 != im.fp64 || src.comps != im.comps || src.kappa.kappa != im.kappa.kappa
        || src.kappa.complex_kappa != im.kappa.complex_kappa || entries == 0 || entries > size_t(INT_MAX)
        || list.Size() != entries || lcount.Size() != entries || ranges.Size() % 2 || ranges.Size() > size_t(INT_MAX))
      throw Exception("FMMDeviceL2P::AddP2P: invalid entries");
    // the host's EvaluateDeriv adds only charges (dipoles add nothing)
    if (!(im.mode & 1) && !(src.kind & 1))
      return;
    for (size_t r = 0; r < nranges; r++)
      if (ranges[2*r] < 0 || ranges[2*r+1] < 1 || size_t(ranges[2*r])+ranges[2*r+1] > nsrc)
        throw Exception("FMMDeviceL2P::AddP2P: a range exceeds the sources");
    std::vector<char> owner (im.order.Size(), 0);
    for (size_t e = 0; e < entries; e++)
      {
        if (target[e] < 0 || size_t(target[e]) >= im.order.Size() || list[e] < 0 || lcount[e] < 1
            || size_t(list[e])+lcount[e] > nranges)
          throw Exception("FMMDeviceL2P::AddP2P: an entry exceeds the targets or the ranges");
        if (owner[target[e]]++)
          throw Exception("FMMDeviceL2P::AddP2P: two entries share a target");
      }
    double ak = std::abs(src.kappa.kappa);
    int branch = ak < 1e-12 ? 0 : ak < 1e-8 ? 1 : 2;      // SingularMLExpansion::Node::Evaluate for charges
    im.p2p.push_back ({ int(entries), src.kind, branch, src.kappa.Arg (src.kappa.kappa, im.fp64),
                        UploadBuffer (im.device, target), UploadBuffer (im.device, list), UploadBuffer (im.device, lcount),
                        UploadBuffer (im.device, ranges), src.xyz, src.dir, src.charge, src.dipole, src.dipvec });
  }

  void FMMDeviceL2P :: Download (const FMMDeviceTargets & targets)
  {
    RegionTimer reg(DownloadTimer());
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    auto & im = *impl;
    size_t n = size_t(im.comps)*im.outs;
    if (targets.comps != im.comps || targets.mode != im.mode || targets.values.Size() != targets.pts.Size()*n)
      throw Exception("FMMDeviceL2P::Download: the targets do not match the plan");
    auto copy = [&] (const auto * src)
    {
      for (size_t j = 0; j < im.order.Size(); j++)
        for (size_t c = 0; c < n; c++)
          targets.values[size_t(im.order[j])*n+c] = src[j*n+c];
    };
    if (im.fp64)
      ReadBuffer<Complex> (im.out, im.order.Size()*n, copy);
    else
      ReadBuffer<Complex32> (im.out, im.order.Size()*n, copy);
  }

  struct FMMDevicePlan::Impl
  {
    struct Group
    {
      FMMTranslation variant;
      std::variant<double,std::complex<double>> kappa;
      double source_rtyp, target_rtyp;
      int os, ot;
      double len, theta;
      RealBuffer src, tgt;
      std::vector<double> phi;
      std::vector<int> flip, src_idx, tgt_idx;
    };
    struct Zero
    {
      RealBuffer pool;
      size_t first, count;
    };
    struct Pass
    {
      std::vector<Zero> zeros;      // run before the pass's groups
      std::vector<Group> groups;
    };
    std::shared_ptr<ngs_gpu::Device> device;
    bool fp64 = false;                   // the real type of all its pools
    std::vector<Pass> passes = std::vector<Pass>(1);
    std::shared_ptr<FMMDeviceP2M> p2m;   // before the first pass
    std::shared_ptr<FMMDeviceL2P> l2p;   // after the last pass

    void UseDevice (const std::shared_ptr<ngs_gpu::Device> & dev, bool real64)
    {
      if (!device)
        {
          device = dev;
          fp64 = real64;
        }
      else if (dev != device || real64 != fp64)
        throw Exception("FMMDevicePlan: the pools live on different devices or in different real types");
    }
    void UseDevice (const FMMDeviceTree & tree) { UseDevice (tree.impl->device, tree.impl->fp64); }

    // everything prepared and uploaded, ready to enqueue: no host transfer from here on
    struct Ready
    {
      struct Unit
      {
        std::shared_ptr<DeviceTransformUnit> unit;
        RealBuffer src, tgt;
      };
      std::shared_ptr<ngs_gpu::Device> device;
      bool fp64 = false;
      std::shared_ptr<ngs_gpu::Kernel> zero;
      std::shared_ptr<FMMDeviceP2M> p2m;
      std::shared_ptr<FMMDeviceL2P> l2p;
      std::vector<std::vector<Zero>> zeros;       // per pass
      std::vector<std::vector<Unit>> units;       // per pass, in order
      TransformNeeds needs;                       // of these units
      std::shared_ptr<TransformWorkspace> workspace;   // one for all passes (Record: shared with another program)
      TransformNeeds workspace_needs;             // what the workspace was reserved for
    };

    static bool SameUnit (const Group & a, const Group & b)
    {
      return a.variant == b.variant && a.os == b.os && a.ot == b.ot
        && a.source_rtyp == b.source_rtyp && a.target_rtyp == b.target_rtyp && a.kappa == b.kappa
        && a.src.buf == b.src.buf && a.tgt.buf == b.tgt.buf;
    }

    // reserve: a workspace for these units alone (Run); Record reserves its own, possibly shared one
    Ready Prepare (bool reserve) const
    {
      RegionTimer reg(PrepareTimer());
      Ready ready;
      ready.device = device;
      ready.fp64 = fp64;
      if (!ready.device) return ready;
      ready.zero = GetTransformKernels (ready.device, fp64)[2];
      ready.p2m = p2m;
      ready.l2p = l2p;
      size_t budget = FMMDeviceBudget (ready.device), fill = FMM_DEVICE_FILL*ready.device->ComputeUnits();
      for (auto & pass : passes)
        {
          ready.zeros.push_back (pass.zeros);
          ready.units.emplace_back();

          // runs of consecutive groups of one key, so the lanes and the scatter order stay in the plan's order
          std::vector<std::vector<const Group*>> runs;
          for (auto & g : pass.groups)
            {
              if (runs.empty() || !SameUnit (*runs.back()[0], g)
                  || g.src.buf == g.tgt.buf)   // a unit gathers all lanes before it scatters
                runs.emplace_back();
              runs.back().push_back (&g);
            }

          // each run packed and uploaded once per table kind, split into units of whole lane ranges by the caps
          for (auto & run : runs)
            {
              const Group & k = *run[0];
              size_t nb = 0;
              for (auto * g : run) nb += g->phi.size();
              Array<double> len(run.size()), theta(run.size()), phi(nb);
              Array<int> lanes(run.size()), flip(nb), src_idx(nb), tgt_idx(nb);
              size_t lane = 0;
              for (size_t i = 0; i < run.size(); i++)
                {
                  auto & g = *run[i];
                  if (g.phi.size() > size_t(INT_MAX))
                    throw Exception("FMMDevicePlan: a group has more than INT_MAX lanes");
                  len[i] = g.len; theta[i] = g.theta; lanes[i] = int(g.phi.size());
                  for (size_t b = 0; b < g.phi.size(); b++, lane++)
                    {
                      phi[lane] = g.phi[b]; flip[lane] = g.flip[b];
                      src_idx[lane] = g.src_idx[b]; tgt_idx[lane] = g.tgt_idx[b];
                    }
                }
              // fused: S->R of low order between distinct pools (a group reads its sources while others scatter)
              bool fused = k.variant == FMMTranslation::SR && std::max (k.os, k.ot) <= FusedMaxOrder (fp64)
                && k.src.buf != k.tgt.buf;
              for (auto & u : PrepareUnits (ready.device, k.variant, k.kappa, k.source_rtyp, k.target_rtyp,
                                            k.os, k.ot, len, theta, lanes, phi, flip, src_idx, tgt_idx, fused, fp64))
                {
                  ready.needs.Add (*u, budget, fill);
                  ready.units.back().push_back ({ u, k.src, k.tgt });
                }
            }
        }
      if (reserve && !ready.needs.Empty())
        {
          ready.workspace = ReserveNeeds (ready.device, ready.needs, fp64);
          ready.workspace_needs = ready.needs;
        }
      return ready;
    }

    // launches only (recordable): per pass the zero fills, then the units in order, all on the one workspace
    static void Enqueue (ngs_gpu::Queue & queue, const Ready & ready)
    {
      using ngs_gpu::Dim3;
      bool cpu = ready.device == ngs_gpu::GetCpuDevice();
      unsigned threads = GroupThreads (ready.device, 256);
      if (ready.p2m)
        ready.p2m->impl->Enqueue (queue);
      for (size_t p = 0; p < ready.units.size(); p++)
        {
          for (auto & z : ready.zeros[p])
            {
              if (z.count == 0) continue;
              unsigned groups = cpu ? 1u : unsigned(std::min<size_t> ((z.count+threads-1)/threads, 4096));
              queue.Launch (*ready.zero, Dim3(groups), Dim3(threads), { z.pool, int(z.first), int(z.count) });
            }
          for (auto & u : ready.units[p])
            EnqueueTransformUnit (queue, *u.unit, *ready.workspace, u.src, u.tgt);
        }
      if (ready.l2p)
        ready.l2p->impl->Enqueue (queue);
    }

    static std::shared_ptr<ngs_gpu::Program> RecordReady (const Ready & ready)
    {
      auto queue = ready.device->DefaultQueue();
      // process-wide until EndRecording: any transfer throws, other threads' launches on this queue are captured
      queue->BeginRecording();
      try
        {
          Enqueue (*queue, ready);
        }
      catch (...)
        {
          queue->EndRecording();              // never leave the process in recording mode (every transfer would throw)
          throw;
        }
      return queue->EndRecording();
    }
  };

  struct FMMDeviceProgram::Impl
  {
    std::shared_ptr<ngs_gpu::Program> program;
    // what the launches reference: the prepared units and the workspace (the recording keeps the buffers alive,
    // this keeps the host structures too and allows a re-recording against a new shared workspace)
    FMMDevicePlan::Impl::Ready ready;
  };

  FMMDeviceProgram :: FMMDeviceProgram () : impl(std::make_unique<Impl>()) { }
  FMMDeviceProgram :: ~FMMDeviceProgram () = default;

  void FMMDeviceProgram :: Replay ()
  {
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    if (!impl->program || impl->program->Size() == 0) return;
    RegionTimer rexec(ExecuteTimer());
    impl->program->Run();                  // Queue::Replay on the recording queue
    impl->program->GetQueue()->Finish();
  }

  FMMDevicePlan :: FMMDevicePlan () : impl(std::make_unique<Impl>()) { }
  FMMDevicePlan :: ~FMMDevicePlan () = default;

  void FMMDevicePlan :: BeginPass ()
  {
    impl->passes.emplace_back();
  }

  void FMMDevicePlan :: AddZero (const FMMDeviceTree & tree, int pool, size_t first, size_t count)
  {
    auto & p = tree.impl->Get (pool);
    if (first > p.buffer.Size() || count > p.buffer.Size()-first || first+count > size_t(INT_MAX))
      throw Exception("FMMDevicePlan::AddZero: range exceeds the pool or the int indexing");
    impl->UseDevice (tree);
    impl->passes.back().zeros.push_back ({ p.buffer, first, count });
  }

  void FMMDevicePlan :: SetP2M (std::shared_ptr<FMMDeviceP2M> p2m)
  {
    impl->UseDevice (p2m->impl->device, p2m->impl->fp64);
    impl->p2m = std::move(p2m);
  }

  void FMMDevicePlan :: SetL2P (std::shared_ptr<FMMDeviceL2P> l2p)
  {
    impl->UseDevice (l2p->impl->device, l2p->impl->fp64);
    impl->l2p = std::move(l2p);
  }

  void FMMDevicePlan :: AddGroup (FMMTranslation variant, std::variant<double,std::complex<double>> kappa,
                                  double source_rtyp, double target_rtyp, int os, int ot, double len, double theta,
                                  const FMMDeviceTree & src, int src_pool, const FMMDeviceTree & tgt, int tgt_pool,
                                  FlatArray<double> phi, FlatArray<int> flip, FlatArray<int> src_idx, FlatArray<int> tgt_idx)
  {
    auto & sp = src.impl->Get (src_pool);
    auto & tp = tgt.impl->Get (tgt_pool);
    if (sp.order != os || tp.order != ot)
      throw Exception("FMMDevicePlan::AddGroup: the orders do not match the pools");
    size_t lanes = phi.Size();
    if (flip.Size() != lanes || src_idx.Size() != lanes || tgt_idx.Size() != lanes)
      throw Exception("FMMDevicePlan::AddGroup: phi, flip, src_idx and tgt_idx need the same size");
    if (lanes == 0) return;
    impl->UseDevice (src);
    impl->UseDevice (tgt);
    Impl::Group g { variant, kappa, source_rtyp, target_rtyp, os, ot, len, theta, sp.buffer, tp.buffer,
                    std::vector<double>(phi.begin(), phi.end()), std::vector<int>(flip.begin(), flip.end()),
                    std::vector<int>(src_idx.begin(), src_idx.end()), std::vector<int>(tgt_idx.begin(), tgt_idx.end()) };
    impl->passes.back().groups.push_back (std::move(g));
  }

  void FMMDevicePlan :: Run ()
  {
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    // all host work and uploads first: on a GPU every H2D waits for the queue
    auto ready = impl->Prepare(true);
    if (!ready.device) return;

    // one in-order queue, from this thread; the prepared units and the workspace live until Finish
    RegionTimer reg(ExecuteTimer());
    auto queue = ready.device->DefaultQueue();
    Impl::Enqueue (*queue, ready);
    queue->Finish();
  }

  std::shared_ptr<FMMDeviceProgram> FMMDevicePlan :: Record (FMMDeviceProgram * share)
  {
    std::lock_guard<std::mutex> lock(fmm_device_mutex);
    auto prog = std::make_shared<FMMDeviceProgram>();
    auto ready = impl->Prepare(false);      // every allocation and upload before the recording starts
    if (!ready.device) return prog;

    // one workspace for this program and share, sized for both; if share's is not exactly that, a new one is
    // reserved and share is re-recorded against it, which frees the old one
    FMMDeviceProgram::Impl * other = (share && share->impl->program) ? share->impl.get() : nullptr;
    if (other && (other->ready.device != ready.device || other->ready.fp64 != ready.fp64))
      throw Exception("FMMDevicePlan::Record: the shared program runs on another device or in another real type");
    TransformNeeds target = ready.needs;
    if (other) target.Max (other->ready.needs);
    std::shared_ptr<TransformWorkspace> ws;
    bool rerecord = false;
    if (other && other->ready.workspace && other->ready.workspace_needs == target)
      ws = other->ready.workspace;
    else if (!target.Empty())
      {
        RegionTimer reg(PrepareTimer());
        ws = ReserveNeeds (ready.device, target, ready.fp64);
        rerecord = other && !other->ready.needs.Empty();
      }
    ready.workspace = ws;
    ready.workspace_needs = target;

    RegionTimer reg(RecordTimer());
    if (rerecord)
      {
        auto rd = other->ready;             // the same prepared units, the new workspace
        rd.workspace = ws;
        rd.workspace_needs = target;
        auto program = Impl::RecordReady (rd);
        other->program = program;
        other->ready = std::move(rd);
      }
    prog->impl->program = Impl::RecordReady (ready);
    prog->impl->ready = std::move(ready);
    return prog;
  }


  shared_ptr<FMMDeviceCache> MakeFMMDeviceCache ()
  {
    return std::make_shared<FMMDeviceCache>();
  }

  bool FMMDeviceMatvec (const FMM_Parameters & params)
  {
    auto * cache = params.device_cache.get();
    if (!cache || !cache->fast || !params.device_sources || !params.device_targets) return false;
    RegionTimer reg(MatvecTimer());
    cache->singular.p2m->Upload (*params.device_sources);
    cache->singular.program->Replay();
    cache->regular.program->Replay();
    cache->regular.l2p->Download (*params.device_targets);
    return true;
  }
}
