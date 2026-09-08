/*********************************************************************/
/* File:   device_blockjacobi.cpp                                    */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   3. Sep. 2026                                              */
/*********************************************************************/

#define FILE_DEVICE_BLOCKJACOBI_CPP

#include <la.hpp>
#include <gpukernel.hpp>
#include <climits>

namespace ngla
{
  using namespace ngs_gpu;

  namespace
  {
    string Substitute (string src, const string & token, const string & value)
    {
      for (size_t p = src.find(token); p != string::npos; p = src.find(token, p+value.size()))
        src.replace (p, token.size(), value);
      return src;
    }

    /*
      Entry (r,c) of the effective block matrix is mat[c*bs+r] if strided,
      mat[r*bs+c] otherwise - the host packs the inverses column- or
      row-major and folds the transpose into the flag.

      Small blocks: lanes consecutive work-items share a block and take
      its rows in turns. The atomic variant is for block tables where a
      dof occurs more than once.
    */
    const char * lane_source = R"RAW(

      #define BJ_BODY(ACCUMULATE)                                            \
        int lane = int(LOCAL_ID_X) & (lanes-1);                              \
        int item = int(GLOBAL_ID_X) / lanes;                                 \
        if (item >= nitems) return;                                          \
        int block = blocklist[item];                                         \
        int first = blockfirst[block];                                       \
        int bs = blockfirst[block+1] - first;                                \
        GLOBAL_PTR(const int) idx = indices + first;                         \
        GLOBAL_PTR(const SCAL) mat = mats + matfirst[block];                 \
        for (int r = lane; r < bs; r += lanes)                               \
          {                                                                  \
            SCAL sum = 0;                                                    \
            if (strided)                                                     \
              for (int c = 0; c < bs; c++) sum += mat[c*bs+r] * x[idx[c]];   \
            else                                                             \
              for (int c = 0; c < bs; c++) sum += mat[r*bs+c] * x[idx[c]];   \
            ACCUMULATE;                                                      \
          }

      KERNEL(bj_mult, GLOBAL_IN(int,blocklist), GLOBAL_IN(int,blockfirst), GLOBAL_IN(int,indices),
                      GLOBAL_IN(int,matfirst), GLOBAL_IN(SCAL,mats), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                      VALUE(SCAL,s), VALUE(int,lanes), VALUE(int,strided), VALUE(int,nitems))
      {
        BJ_BODY(y[idx[r]] += s*sum)
      }

      KERNEL(bj_mult_atomic, GLOBAL_IN(int,blocklist), GLOBAL_IN(int,blockfirst), GLOBAL_IN(int,indices),
                             GLOBAL_IN(int,matfirst), GLOBAL_IN(SCAL,mats), GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y),
                             VALUE(SCAL,s), VALUE(int,lanes), VALUE(int,strided), VALUE(int,nitems))
      {
        BJ_BODY(ATOMIC_ADD(&y[idx[r]], s*sum))
      }

    )RAW";


    /*
      Large blocks: one warp per block, the block's x entries staged in
      group memory ($MAXBS per warp), rows taken in turns by the 32 lanes,
      up to four rows per lane at once for independent load chains.
      Compiled per maximal block size.
    */
    const char * warp_source = R"RAW(

      #define ACC_PLAIN(row,val)   y[idx[row]] += s*(val)
      #define ACC_ATOMIC(row,val)  ATOMIC_ADD(&y[idx[row]], s*(val))

      // entry (r,c) of the block
      #define BJ_M(r,c)  (strided ? mat[(c)*bs+(r)] : mat[(r)*bs+(c)])

      #define BJ_WARP_BODY(ACC)                                              \
        uint tid = LOCAL_ID_X;                                               \
        int lane = int(tid & 31);                                            \
        int warp = int(tid >> 5);                                            \
        $XS_DECL                                                             \
        int item = int(GROUP_ID_X)*$WARPS + warp;                            \
        if (item >= nitems) return;                                          \
        int block = blocklist[item];                                         \
        int first = blockfirst[block];                                       \
        int bs = blockfirst[block+1] - first;                                \
        GLOBAL_PTR(const int) idx = indices + first;                         \
        GLOBAL_PTR(const SCAL) mat = mats + matfirst[block];                 \
        $XS_LOAD                                                             \
        int r = lane;                                                        \
        for ( ; r + 96 < bs; r += 128)                                       \
          {                                                                  \
            SCAL s0 = 0, s1 = 0, s2 = 0, s3 = 0;                             \
            for (int c = 0; c < bs; c++)                                     \
              {                                                              \
                SCAL xc = $XC;                                               \
                s0 += BJ_M(r,c) * xc;    s1 += BJ_M(r+32,c) * xc;            \
                s2 += BJ_M(r+64,c) * xc; s3 += BJ_M(r+96,c) * xc;            \
              }                                                              \
            ACC(r, s0); ACC(r+32, s1); ACC(r+64, s2); ACC(r+96, s3);         \
          }                                                                  \
        for ( ; r + 32 < bs; r += 64)                                        \
          {                                                                  \
            SCAL s0 = 0, s1 = 0;                                             \
            for (int c = 0; c < bs; c++)                                     \
              {                                                              \
                SCAL xc = $XC;                                               \
                s0 += BJ_M(r,c) * xc; s1 += BJ_M(r+32,c) * xc;               \
              }                                                              \
            ACC(r, s0); ACC(r+32, s1);                                       \
          }                                                                  \
        for ( ; r < bs; r += 32)                                             \
          {                                                                  \
            SCAL s0 = 0;                                                     \
            for (int c = 0; c < bs; c++) s0 += BJ_M(r,c) * $XC;              \
            ACC(r, s0);                                                      \
          }

      KERNEL(bj_warp, GLOBAL_IN(int,blocklist), GLOBAL_IN(int,blockfirst), GLOBAL_IN(int,indices),
                      GLOBAL_IN(int,matfirst), GLOBAL_IN(SCAL,mats), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                      VALUE(SCAL,s), VALUE(int,strided), VALUE(int,nitems))
      {
        BJ_WARP_BODY(ACC_PLAIN)
      }

      KERNEL(bj_warp_atomic, GLOBAL_IN(int,blocklist), GLOBAL_IN(int,blockfirst), GLOBAL_IN(int,indices),
                             GLOBAL_IN(int,matfirst), GLOBAL_IN(SCAL,mats), GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y),
                             VALUE(SCAL,s), VALUE(int,strided), VALUE(int,nitems))
      {
        BJ_WARP_BODY(ACC_ATOMIC)
      }

    )RAW";


    template <typename T>
    class DeviceBlockJacobiKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<ngs_gpu::Queue> queue;
      shared_ptr<Kernel> mult, mult_atomic;
      unsigned groupsize;

      DeviceBlockJacobiKernels (shared_ptr<Device> adevice)
        : device(adevice)
      {
        if constexpr (is_same_v<T,double>)
          if (!device->HasFloat64())
            throw Exception("DeviceBlockJacobi<double> on "+device->Name()+
                            ", which has no fp64 - use DeviceBlockJacobi<float>");

        string scal = is_same_v<T,double> ? "double" : "float";
        library = device->CompileSource (string(code_gpukernel) +
                                         Substitute (lane_source, "SCAL", scal));
        mult        = library->GetKernel ("bj_mult");
        mult_atomic = library->GetKernel ("bj_mult_atomic");
        queue = device->DefaultQueue();
        groupsize = (device->SimdWidth() > 1) ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
      }

      // a power of two up to maxlanes, matching the average block size
      int ChooseLanes (double avg, int maxlanes) const
      {
        size_t simd = device->SimdWidth();
        if (simd <= 1) return 1;
        int lanes = 1;
        while (lanes < avg && size_t(2*lanes) <= min<size_t>(simd, maxlanes) && size_t(2*lanes) <= groupsize)
          lanes *= 2;
        return lanes;
      }

      static const DeviceBlockJacobiKernels & Get()
      {
        static mutex mtx;
        static shared_ptr<DeviceBlockJacobiKernels> cached;

        auto dev = GetGpuDevice();
        auto lock = lock_guard<mutex>(mtx);
        if (!cached || cached->device != dev)
          cached = make_shared<DeviceBlockJacobiKernels> (dev);
        return *cached;
      }
    };
  }


  template <typename T>
  class DeviceBlockJacobiWarpKernels
  {
    shared_ptr<Library> library;
  public:
    shared_ptr<Kernel> mult, mult_atomic;
    int maxbs, warps;
    bool staged;

    DeviceBlockJacobiWarpKernels (shared_ptr<Device> device, int amaxbs)
      : maxbs(amaxbs)
    {
      // group memory for the staged x, metal has 32 KB
      size_t per_warp = size_t(maxbs)*sizeof(T);
      staged = per_warp <= 24*1024;
      warps = staged ? int(min<size_t> (8, (24*1024)/per_warp)) : 8;
      if (auto e = getenv("NGS_BJ_WARPS")) warps = atoi(e);

      string code = Substitute (warp_source, "SCAL", is_same_v<T,double> ? "double" : "float");
      code = Substitute (code, "$WARPS", ToString(warps));
      if (staged)
        {
          code = Substitute (code, "$XS_DECL", "SHARED(SCAL, xs_all, " + ToString(warps*maxbs) + "); "
                             "LOCAL_PTR(SCAL) xs = xs_all + warp*" + ToString(maxbs) + ";");
          code = Substitute (code, "$XS_LOAD", "for (int c = lane; c < bs; c += 32) xs[c] = x[idx[c]]; SIMD_BARRIER();");
          code = Substitute (code, "$XC", "xs[c]");
        }
      else
        {
          code = Substitute (code, "$XS_DECL", "");
          code = Substitute (code, "$XS_LOAD", "");
          code = Substitute (code, "$XC", "x[idx[c]]");
        }
      code = Substitute (code, "SCAL", is_same_v<T,double> ? "double" : "float");

      library = device->CompileSource (string(code_gpukernel) + code);
      mult        = library->GetKernel ("bj_warp");
      mult_atomic = library->GetKernel ("bj_warp_atomic");

      if (getenv("NGS_BJ_INFO"))
        cout << "bj_warp maxbs=" << maxbs << " warps=" << warps << " staged=" << staged
             << " " << mult->Info(warps*32) << endl;
    }

    // maxbs rounded up to a power of two, so few variants get compiled
    static shared_ptr<const DeviceBlockJacobiWarpKernels> Get (shared_ptr<Device> device, int maxbs)
    {
      static mutex mtx;
      static map<pair<Device*,int>, shared_ptr<const DeviceBlockJacobiWarpKernels>> cache;

      if (device->SimdWidth() != 32) return nullptr;
      int cls = 32;
      while (cls < maxbs) cls *= 2;
      auto lock = lock_guard<mutex>(mtx);
      auto & entry = cache[{device.get(), cls}];
      if (!entry)
        entry = make_shared<DeviceBlockJacobiWarpKernels> (device, cls);
      return entry;
    }
  };



  template <typename T>
  template <typename TM>
  DeviceBlockJacobi<T> :: DeviceBlockJacobi (const BlockJacobiPrecond<TM> & pre)
  {
    height = pre.Height();
    width = pre.Width();
    const Table<int> & blocktable = *pre.GetBlockTable();
    nblocks = blocktable.Size();
    const auto & inverses = pre.GetInverses();

    const auto & kern = DeviceBlockJacobiKernels<T>::Get();
    device = kern.device;
    queue = kern.queue;
    memtype = PreferredMemType();

    // offsets into the concatenated index and matrix arrays
    Array<int> blockfirst(nblocks+1), matfirst(nblocks+1);
    size_t nind = 0, nmat = 0;
    for (size_t i = 0; i < nblocks; i++)
      {
        blockfirst[i] = int(nind);
        matfirst[i] = int(nmat);
        size_t bs = blocktable[i].Size();
        nind += bs;
        nmat += bs*bs;
        if (nmat > size_t(INT_MAX))
          throw Exception("DeviceBlockJacobi: block inverses too large for 32-bit offsets");
      }
    blockfirst[nblocks] = int(nind);
    matfirst[nblocks] = int(nmat);

    // the table is already sorted per block, its data is contiguous
    FlatArray<int> indices = blocktable.AsArray();
    if (indices.Size() != nind)
      throw Exception("DeviceBlockJacobi: block table is not contiguous");

    // symmetric inverses: the transpose is the matrix itself
    symmetric = true;
    for (size_t i = 0; i < nblocks && symmetric; i++)
      {
        const auto & inv = inverses[i];
        size_t bs = blocktable[i].Size();
        double scale = 0, asym = 0;
        for (size_t r = 0; r < bs; r++)
          for (size_t c = 0; c < bs; c++)
            {
              scale = max (scale, double(fabs(inv(r,c))));
              asym = max (asym, double(fabs(inv(r,c)-inv(c,r))));
            }
        if (asym > 1e-10*scale) symmetric = false;
      }

    // column-major gives coalesced reads across the lanes of a warp
    colmajor = true;
    if (auto e = getenv("NGS_BJ_COLMAJOR")) colmajor = atoi(e);

    Array<T> mats(nmat);
    for (size_t i = 0; i < nblocks; i++)
      {
        size_t bs = blocktable[i].Size();
        T * dst = mats.Data() + matfirst[i];
        if (bs == 0) continue;
        const auto & inv = inverses[i];
        for (size_t r = 0; r < bs; r++)
          for (size_t c = 0; c < bs; c++)
            dst[colmajor ? c*bs+r : r*bs+c] = T(inv(r,c));
      }

    // disjoint blocks can be accumulated without atomics
    Array<int> count(width);
    count = 0;
    overlapping = false;
    for (int d : indices)
      if (count[d]++ > 0) { overlapping = true; break; }

    // size classes: empty blocks are dropped, large ones take a warp each
    int large_from = 16;
    if (auto e = getenv("NGS_BJ_LARGE")) large_from = atoi(e);
    if (device->SimdWidth() != 32) large_from = INT_MAX;
    Array<int> small, large;
    size_t nind_small = 0;
    maxbs = 0;
    for (size_t i = 0; i < nblocks; i++)
      {
        int bs = blocktable[i].Size();
        if (bs == 0) continue;
        if (bs > large_from) { large.Append(i); maxbs = max(maxbs, bs); }
        else { small.Append(i); nind_small += bs; }
      }
    nsmall = small.Size();
    nlarge = large.Size();
    lanes = nsmall ? kern.ChooseLanes (double(nind_small)/nsmall, large_from) : 1;
    if (nlarge)
      warpkern = DeviceBlockJacobiWarpKernels<T>::Get (device, maxbs);

    cout << IM(7) << "DeviceBlockJacobi<" << (is_same_v<T,double> ? "double" : "float")
         << "> nblocks = " << nblocks << " (" << nsmall << " small, " << nlarge << " large, max " << maxbs << ")"
         << ", indices = " << nind << ", matrix entries = " << nmat
         << ", lanes = " << lanes << (overlapping ? ", overlapping" : ", disjoint")
         << (symmetric ? ", symmetric" : "") << endl;

    dev_blockfirst = device->NewBuffer<int> (nblocks+1, MemType::Device);
    dev_matfirst   = device->NewBuffer<int> (nblocks+1, MemType::Device);
    dev_indices    = device->NewBuffer<int> (max<size_t>(nind,1), MemType::Device);
    dev_mats       = device->NewBuffer<T> (max<size_t>(nmat,1), MemType::Device);
    dev_small      = device->NewBuffer<int> (max<size_t>(nsmall,1), MemType::Device);
    dev_large      = device->NewBuffer<int> (max<size_t>(nlarge,1), MemType::Device);

    dev_blockfirst.H2D (blockfirst.Data(), nblocks+1);
    dev_matfirst.H2D (matfirst.Data(), nblocks+1);
    dev_indices.H2D (indices.Data(), nind);
    dev_mats.H2D (mats.Data(), nmat);
    if (nsmall) dev_small.H2D (small.Data(), nsmall);
    if (nlarge) dev_large.H2D (large.Data(), nlarge);
  }


  template <typename T>
  void DeviceBlockJacobi<T> :: Launch (const BaseVector & x, BaseVector & y, T s, bool trans) const
  {
    if (x.Size() != width || y.Size() != height)
      throw Exception("DeviceBlockJacobi::MultAdd - size mismatch");
    if (nsmall + nlarge == 0) return;

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);

    const auto & kern = DeviceBlockJacobiKernels<T>::Get();
    int strided = (colmajor != (trans && !symmetric)) ? 1 : 0;

    if (nsmall)
      {
        size_t items = nsmall * lanes;
        unsigned groups = (items + kern.groupsize-1) / kern.groupsize;
        queue->Launch (overlapping ? *kern.mult_atomic : *kern.mult,
                       Dim3(groups), Dim3(kern.groupsize),
                       { KernelArg(dev_small), KernelArg(dev_blockfirst), KernelArg(dev_indices), KernelArg(dev_matfirst),
                         KernelArg(dev_mats), ux.DevArgRO(), uy.DevArgRW(),
                         KernelArg(s), KernelArg(int(lanes)), KernelArg(strided), KernelArg(int(nsmall)) });
      }
    if (nlarge)
      {
        unsigned groups = (nlarge + warpkern->warps-1) / warpkern->warps;
        queue->Launch (overlapping ? *warpkern->mult_atomic : *warpkern->mult,
                       Dim3(groups), Dim3(32*warpkern->warps),
                       { KernelArg(dev_large), KernelArg(dev_blockfirst), KernelArg(dev_indices), KernelArg(dev_matfirst),
                         KernelArg(dev_mats), ux.DevArgRO(), uy.DevArgRW(),
                         KernelArg(s), KernelArg(strided), KernelArg(int(nlarge)) });
      }
  }

  template <typename T>
  void DeviceBlockJacobi<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceBlockJacobi::MultAdd"); RegionTimer reg(t);
    Launch (x, y, T(s), false);
  }

  template <typename T>
  void DeviceBlockJacobi<T> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceBlockJacobi::MultTransAdd"); RegionTimer reg(t);
    Launch (x, y, T(s), true);
  }




  template <typename T>
  BaseMatrix::OperatorInfo DeviceBlockJacobi<T> :: GetOperatorInfo () const
  {
    return { string("DeviceBlockJacobi<") + (is_same_v<T,double> ? "double" : "float")
             + "> (blocks=" + ToString(nblocks) + ")", height, width };
  }

  template <typename T>
  ostream & DeviceBlockJacobi<T> :: Print (ostream & ost) const
  {
    ost << "DeviceBlockJacobi<" << (is_same_v<T,double> ? "double" : "float")
        << ">, height = " << height << ", nblocks = " << nblocks
        << " (" << nsmall << " small, " << nlarge << " large, max " << maxbs << ")"
        << ", on " << device->Name() << endl;
    return ost;
  }


  template class DeviceBlockJacobi<double>;
  template class DeviceBlockJacobi<float>;
  template DeviceBlockJacobi<double>::DeviceBlockJacobi (const BlockJacobiPrecond<double> &);
  template DeviceBlockJacobi<double>::DeviceBlockJacobi (const BlockJacobiPrecond<float> &);
  template DeviceBlockJacobi<float>::DeviceBlockJacobi (const BlockJacobiPrecond<double> &);
  template DeviceBlockJacobi<float>::DeviceBlockJacobi (const BlockJacobiPrecond<float> &);
}
