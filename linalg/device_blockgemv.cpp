/*********************************************************************/
/* File:   device_blockgemv.cpp                                      */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   9. Sep. 2026                                              */
/*********************************************************************/

#define FILE_DEVICE_BLOCKGEMV_CPP

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
      Entry (r,c) of block b is mat[c*nout+r]; the transposed view reads
      the stored matrix with in/out swapped, i.e. mat[r*nin+c] (strided).

      Small blocks: lanes consecutive work-items share a block and take
      its output rows in turns.
    */
    const char * lane_source = R"RAW(

      #define BG_M(r,c)  (strided ? mat[(r)*nin+(c)] : mat[(c)*nout+(r)])

      #define BG_BODY(ACC)                                                   \
        int lane = int(LOCAL_ID_X) & (lanes-1);                              \
        int item = int(GLOBAL_ID_X) / lanes;                                 \
        if (item >= nitems) return;                                          \
        int block = blocklist[item];                                         \
        int i0 = infirst[block], nin = infirst[block+1] - i0;                \
        int o0 = outfirst[block], nout = outfirst[block+1] - o0;             \
        GLOBAL_PTR(const int) in = inidx + i0;                               \
        GLOBAL_PTR(const int) out = outidx + o0;                             \
        GLOBAL_PTR(const SCAL) mat = mats + matfirst[block];                 \
        for (int r = lane; r < nout; r += lanes)                             \
          {                                                                  \
            SCAL sum = 0;                                                    \
            for (int c = 0; c < nin; c++) sum += BG_M(r,c) * x[in[c]];       \
            ACC(out[r], s*sum);                                              \
          }

      #define ACC_PLAIN(d,v)   y[d] += v
      #define ACC_ATOMIC(d,v)  ATOMIC_ADD(&y[d], v)

      KERNEL(bg_lane, GLOBAL_IN(int,blocklist), GLOBAL_IN(int,infirst), GLOBAL_IN(int,outfirst),
                      GLOBAL_IN(int,matfirst), GLOBAL_IN(int,inidx), GLOBAL_IN(int,outidx),
                      GLOBAL_IN(SCAL,mats), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                      VALUE(SCAL,s), VALUE(int,lanes), VALUE(int,strided), VALUE(int,nitems))
      { BG_BODY(ACC_PLAIN) }

      KERNEL(bg_lane_atomic, GLOBAL_IN(int,blocklist), GLOBAL_IN(int,infirst), GLOBAL_IN(int,outfirst),
                             GLOBAL_IN(int,matfirst), GLOBAL_IN(int,inidx), GLOBAL_IN(int,outidx),
                             GLOBAL_IN(SCAL,mats), GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y),
                             VALUE(SCAL,s), VALUE(int,lanes), VALUE(int,strided), VALUE(int,nitems))
      { BG_BODY(ACC_ATOMIC) }

    )RAW";

    /*
      Large blocks: one group per block ($WARPS warps), the block's x
      staged in group memory ($MAXIN entries), output rows over the
      group's lanes, up to four at once for independent load chains.
      Compiled per maximal input size.
    */
    const char * group_source = R"RAW(

      #define BG_M(r,c)  (strided ? mat[(r)*nin+(c)] : mat[(c)*nout+(r)])

      #define BG_GROUP_BODY(ACC)                                             \
        constexpr int gsize = 32*$WARPS;                                     \
        int tid = int(LOCAL_ID_X);                                           \
        SHARED(SCAL, xs, $MAXIN);                                            \
        int item = int(GROUP_ID_X);                                          \
        if (item >= nitems) return;                                          \
        int block = blocklist[item];                                         \
        int i0 = infirst[block], nin = infirst[block+1] - i0;                \
        int o0 = outfirst[block], nout = outfirst[block+1] - o0;             \
        GLOBAL_PTR(const int) in = inidx + i0;                               \
        GLOBAL_PTR(const int) out = outidx + o0;                             \
        GLOBAL_PTR(const SCAL) mat = mats + matfirst[block];                 \
        for (int c = tid; c < nin; c += gsize) xs[c] = x[in[c]];             \
        BARRIER();                                                           \
        int r = tid;                                                         \
        for ( ; r + 3*gsize < nout; r += 4*gsize)                            \
          {                                                                  \
            SCAL s0 = 0, s1 = 0, s2 = 0, s3 = 0;                             \
            for (int c = 0; c < nin; c++)                                    \
              {                                                              \
                SCAL xc = xs[c];                                             \
                s0 += BG_M(r,c) * xc;         s1 += BG_M(r+gsize,c) * xc;    \
                s2 += BG_M(r+2*gsize,c) * xc; s3 += BG_M(r+3*gsize,c) * xc;  \
              }                                                              \
            ACC(out[r], s*s0); ACC(out[r+gsize], s*s1);                      \
            ACC(out[r+2*gsize], s*s2); ACC(out[r+3*gsize], s*s3);            \
          }                                                                  \
        for ( ; r + gsize < nout; r += 2*gsize)                              \
          {                                                                  \
            SCAL s0 = 0, s1 = 0;                                             \
            for (int c = 0; c < nin; c++)                                    \
              {                                                              \
                SCAL xc = xs[c];                                             \
                s0 += BG_M(r,c) * xc; s1 += BG_M(r+gsize,c) * xc;            \
              }                                                              \
            ACC(out[r], s*s0); ACC(out[r+gsize], s*s1);                      \
          }                                                                  \
        for ( ; r < nout; r += gsize)                                        \
          {                                                                  \
            SCAL s0 = 0;                                                     \
            for (int c = 0; c < nin; c++) s0 += BG_M(r,c) * xs[c];           \
            ACC(out[r], s*s0);                                               \
          }

      #define ACC_PLAIN(d,v)   y[d] += v
      #define ACC_ATOMIC(d,v)  ATOMIC_ADD(&y[d], v)

      KERNEL(bg_group, GLOBAL_IN(int,blocklist), GLOBAL_IN(int,infirst), GLOBAL_IN(int,outfirst),
                       GLOBAL_IN(int,matfirst), GLOBAL_IN(int,inidx), GLOBAL_IN(int,outidx),
                       GLOBAL_IN(SCAL,mats), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                       VALUE(SCAL,s), VALUE(int,strided), VALUE(int,nitems))
      { BG_GROUP_BODY(ACC_PLAIN) }

      KERNEL(bg_group_atomic, GLOBAL_IN(int,blocklist), GLOBAL_IN(int,infirst), GLOBAL_IN(int,outfirst),
                              GLOBAL_IN(int,matfirst), GLOBAL_IN(int,inidx), GLOBAL_IN(int,outidx),
                              GLOBAL_IN(SCAL,mats), GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y),
                              VALUE(SCAL,s), VALUE(int,strided), VALUE(int,nitems))
      { BG_GROUP_BODY(ACC_ATOMIC) }

    )RAW";


    // the lane kernels, one library per device
    template <typename T>
    class LaneKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<Kernel> mult, mult_atomic;
      unsigned groupsize;

      LaneKernels (shared_ptr<Device> adevice)
        : device(adevice)
      {
        if constexpr (is_same_v<T,double>)
          if (!device->HasFloat64())
            throw Exception("DeviceBlockGemv<double> on "+device->Name()+", which has no fp64");
        library = device->CompileSource (string(code_gpukernel) +
                                         Substitute (lane_source, "SCAL", is_same_v<T,double> ? "double" : "float"));
        mult        = library->GetKernel ("bg_lane");
        mult_atomic = library->GetKernel ("bg_lane_atomic");
        groupsize = (device->SimdWidth() > 1) ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
      }

      static const LaneKernels & Get (const shared_ptr<Device> & dev)
      {
        static mutex mtx;
        static shared_ptr<LaneKernels> cached;
        auto lock = lock_guard<mutex>(mtx);
        if (!cached || cached->device != dev)
          cached = make_shared<LaneKernels> (dev);
        return *cached;
      }
    };
  }


  template <typename T>
  class DeviceBlockGemvKernels
  {
    shared_ptr<Library> library;
  public:
    shared_ptr<Kernel> mult, mult_atomic;
    int maxin, warps;

    DeviceBlockGemvKernels (shared_ptr<Device> device, int amaxin)
      : maxin(amaxin)
    {
      // group memory for the staged x, metal has 32 KB
      warps = min (8, (maxin+31)/32);
      if (auto e = getenv("NGS_BG_WARPS")) warps = atoi(e);
      string code = Substitute (group_source, "SCAL", is_same_v<T,double> ? "double" : "float");
      code = Substitute (code, "$WARPS", ToString(warps));
      code = Substitute (code, "$MAXIN", ToString(maxin));
      library = device->CompileSource (string(code_gpukernel) + code);
      mult        = library->GetKernel ("bg_group");
      mult_atomic = library->GetKernel ("bg_group_atomic");
      if (getenv("NGS_BG_INFO"))
        cout << "bg_group maxin=" << maxin << " warps=" << warps << " " << mult->Info(warps*32) << endl;
    }

    static bool Supported (const Device & device, int maxin)
    { return device.SimdWidth() == 32 && size_t(maxin)*sizeof(T) <= 24*1024; }

    // maxin rounded up to a power of two, so few variants get compiled
    static shared_ptr<const DeviceBlockGemvKernels> Get (shared_ptr<Device> device, int maxin)
    {
      static mutex mtx;
      static map<pair<Device*,int>, shared_ptr<const DeviceBlockGemvKernels>> cache;
      int cls = 32;
      while (cls < maxin) cls *= 2;
      auto lock = lock_guard<mutex>(mtx);
      auto & entry = cache[{device.get(), cls}];
      if (!entry)
        entry = make_shared<DeviceBlockGemvKernels> (device, cls);
      return entry;
    }
  };



  template <typename T>
  DeviceBlockGemv<T> :: DeviceBlockGemv (shared_ptr<Device> adevice, const BlockGemvBuilder<T> & b,
                                         size_t width, size_t height)
    : device(adevice), strided(false)
  {
    queue = device->DefaultQueue();
    LaneKernels<T>::Get (device);   // checks fp64 support
    nblocks = b.NBlocks();
    if (b.mats.Size() > size_t(INT_MAX))
      throw Exception("DeviceBlockGemv: block matrices too large for 32-bit offsets");

    nin.SetSize(nblocks); nout.SetSize(nblocks);
    for (size_t i = 0; i < nblocks; i++)
      {
        nin[i] = b.infirst[i+1]-b.infirst[i];
        nout[i] = b.outfirst[i+1]-b.outfirst[i];
      }

    static Timer tup("DeviceBlockGemv ctor upload");
    auto upload = [&] (auto & buf, const auto & arr)
    {
      using TB = typename std::remove_reference_t<decltype(arr)>::value_type;
      buf = device->NewBuffer<TB> (max<size_t>(arr.Size(),1), MemType::Device);
      if (arr.Size()) buf.H2D (arr.Data(), arr.Size());
    };
    tup.Start();
    upload (dev_infirst, b.infirst);
    upload (dev_outfirst, b.outfirst);
    upload (dev_matfirst, b.matfirst);
    upload (dev_inidx, b.inidx);
    upload (dev_outidx, b.outidx);
    upload (dev_mats, b.mats);
    tup.Stop();

    Classify ();
  }

  // size classes and atomic / plain accumulation for the current in/out roles
  template <typename T>
  void DeviceBlockGemv<T> :: Classify ()
  {
    static Timer t("DeviceBlockGemv classify"), tdown("DeviceBlockGemv classify outidx D2H"),
      tatomic("DeviceBlockGemv classify atomic scan");
    RegionTimer reg(t);

    int large_from = 16;
    if (auto e = getenv("NGS_BG_LARGE")) large_from = atoi(e);

    Array<int> small, large;
    size_t nout_small = 0;
    maxin = 0;
    for (size_t i = 0; i < nblocks; i++)
      {
        if (nout[i] > large_from && DeviceBlockGemvKernels<T>::Supported (*device, nin[i]))
          { large.Append(i); maxin = max(maxin, nin[i]); }
        else
          { small.Append(i); nout_small += nout[i]; }
      }
    nsmall = small.Size();
    nlarge = large.Size();

    // a power of two up to large_from, matching the average number of output rows
    lanes = 1;
    double avg = nsmall ? double(nout_small)/nsmall : 1;
    if (device->SimdWidth() > 1)
      while (lanes < avg && 2*lanes <= large_from) lanes *= 2;

    kern_large = nlarge ? DeviceBlockGemvKernels<T>::Get (device, maxin) : nullptr;

    // output dofs in more than one block need atomic accumulation
    Array<int> outidx (dev_outidx.Size());
    tdown.Start();
    dev_outidx.D2H (outidx.Data(), outidx.Size());
    tdown.Stop();
    tatomic.Start();
    size_t total = 0;
    for (size_t i = 0; i < nblocks; i++) total += nout[i];
    int maxdof = -1;
    for (size_t i = 0; i < total; i++) maxdof = max (maxdof, outidx[i]);
    BitArray seen (maxdof+1);
    seen.Clear();
    atomic = false;
    for (size_t i = 0; i < total && !atomic; i++)
      {
        if (seen.Test(outidx[i])) atomic = true;
        seen.SetBit (outidx[i]);
      }
    tatomic.Stop();

    dev_small = device->NewBuffer<int> (max<size_t>(nsmall,1), MemType::Device);
    dev_large = device->NewBuffer<int> (max<size_t>(nlarge,1), MemType::Device);
    if (nsmall) dev_small.H2D (small.Data(), nsmall);
    if (nlarge) dev_large.H2D (large.Data(), nlarge);
  }

  template <typename T>
  shared_ptr<DeviceBlockGemv<T>> DeviceBlockGemv<T> :: Transpose () const
  {
    auto t = shared_ptr<DeviceBlockGemv> (new DeviceBlockGemv);
    t->device = device;
    t->queue = queue;
    t->nblocks = nblocks;
    t->strided = !strided;
    t->nin = nout;
    t->nout = nin;
    t->dev_infirst = dev_outfirst;
    t->dev_outfirst = dev_infirst;
    t->dev_matfirst = dev_matfirst;
    t->dev_inidx = dev_outidx;
    t->dev_outidx = dev_inidx;
    t->dev_mats = dev_mats;
    t->Classify ();
    return t;
  }

  template <typename T>
  string DeviceBlockGemv<T> :: Info () const
  {
    return ToString(nblocks) + " blocks (" + ToString(nsmall) + " small, " + ToString(nlarge) + " large, max in "
      + ToString(maxin) + "), lanes = " + ToString(lanes) + (atomic ? ", atomic" : ", disjoint")
      + (strided ? ", transposed" : "");
  }

  template <typename T>
  void DeviceBlockGemv<T> :: MultAdd (T s, KernelArg x, KernelArg y) const
  {
    const auto & lane = LaneKernels<T>::Get (device);
    if (nsmall)
      {
        size_t items = nsmall * lanes;
        unsigned groups = (items + lane.groupsize-1) / lane.groupsize;
        queue->Launch (atomic ? *lane.mult_atomic : *lane.mult, Dim3(groups), Dim3(lane.groupsize),
                       { KernelArg(dev_small), KernelArg(dev_infirst), KernelArg(dev_outfirst), KernelArg(dev_matfirst),
                         KernelArg(dev_inidx), KernelArg(dev_outidx), KernelArg(dev_mats), x, y,
                         KernelArg(s), KernelArg(lanes), KernelArg(int(strided)), KernelArg(int(nsmall)) });
      }
    if (nlarge)
      queue->Launch (atomic ? *kern_large->mult_atomic : *kern_large->mult, Dim3(nlarge), Dim3(32*kern_large->warps),
                     { KernelArg(dev_large), KernelArg(dev_infirst), KernelArg(dev_outfirst), KernelArg(dev_matfirst),
                       KernelArg(dev_inidx), KernelArg(dev_outidx), KernelArg(dev_mats), x, y,
                       KernelArg(s), KernelArg(int(strided)), KernelArg(int(nlarge)) });
  }

  template class DeviceBlockGemv<double>;
  template class DeviceBlockGemv<float>;
}
