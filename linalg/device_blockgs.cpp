/*********************************************************************/
/* File:   device_blockgs.cpp                                        */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   8. Sep. 2026                                              */
/*********************************************************************/

#define FILE_DEVICE_BLOCKGS_CPP

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
      One group per block of the current color, one warp per 32 rows (at
      most 8). The residual of each block row is a reduction over its csr
      entries by lpr lanes (from the average row length), into group memory; then x_block += inv * res with the
      column-major inverse (coalesced across lanes), up to four rows per
      lane at once. Blocks of one color share no coupled dofs, so the x
      updates of a launch are independent.
    */
    const char * kernel_source = R"RAW(

      #define BGS_M(r,c)  mat[(c)*bs+(r)]

      KERNEL(bgs_sweep, GLOBAL_IN(int,blocklist), GLOBAL_IN(int,blockfirst), GLOBAL_IN(int,indices),
                        GLOBAL_IN(int,matfirst), GLOBAL_IN(SCAL,mats),
                        GLOBAL_IN(int,firsti), GLOBAL_IN(int,colnr), GLOBAL_IN(SCAL,values),
                        GLOBAL_IN(SCAL,b), GLOBAL(SCAL,x), VALUE(int,nitems))
      {
        constexpr int warps = $WARPS;
        constexpr int gsize = 32*warps;
        int tid = int(LOCAL_ID_X);
        int lane = tid & 31;
        int warp = tid >> 5;
        SHARED(SCAL, hx, $MAXBS);

        int item = int(GROUP_ID_X);
        if (item >= nitems) return;
        int block = blocklist[item];
        int first = blockfirst[block];
        int bs = blockfirst[block+1] - first;
        GLOBAL_PTR(const int) idx = indices + first;
        GLOBAL_PTR(const SCAL) mat = mats + matfirst[block];

        // residual of the block rows, lpr lanes per row; the shuffles need the whole warp
        constexpr int lpr = $LPR, rpw = 32/lpr;
        int sublane = lane % lpr;
        for (int j0 = warp*rpw; j0 < bs; j0 += warps*rpw)
          {
            int j = j0 + lane/lpr;
            bool active = j < bs;
            SCAL sum = 0;
            if (active)
              {
                int row = idx[j];
                for (int k = firsti[row]+sublane; k < firsti[row+1]; k += lpr)
                  sum += values[k] * x[colnr[k]];
              }
            for (int o = lpr/2; o > 0; o >>= 1)
              sum += SIMD_SHUFFLE_XOR(sum, o);
            if (active && sublane == 0) hx[j] = b[idx[j]] - sum;
          }
        BARRIER();

        // x_block += inv * hx, rows over the group's lanes
        int r = tid;
        for ( ; r + 3*gsize < bs; r += 4*gsize)
          {
            SCAL s0 = 0, s1 = 0, s2 = 0, s3 = 0;
            for (int c = 0; c < bs; c++)
              {
                SCAL xc = hx[c];
                s0 += BGS_M(r,c) * xc;         s1 += BGS_M(r+gsize,c) * xc;
                s2 += BGS_M(r+2*gsize,c) * xc; s3 += BGS_M(r+3*gsize,c) * xc;
              }
            x[idx[r]] += s0; x[idx[r+gsize]] += s1; x[idx[r+2*gsize]] += s2; x[idx[r+3*gsize]] += s3;
          }
        for ( ; r + gsize < bs; r += 2*gsize)
          {
            SCAL s0 = 0, s1 = 0;
            for (int c = 0; c < bs; c++)
              {
                SCAL xc = hx[c];
                s0 += BGS_M(r,c) * xc; s1 += BGS_M(r+gsize,c) * xc;
              }
            x[idx[r]] += s0; x[idx[r+gsize]] += s1;
          }
        for ( ; r < bs; r += gsize)
          {
            SCAL s0 = 0;
            for (int c = 0; c < bs; c++) s0 += BGS_M(r,c) * hx[c];
            x[idx[r]] += s0;
          }
      }

    )RAW";
  }


  template <typename T>
  class DeviceBlockGSKernels
  {
    shared_ptr<Library> library;
  public:
    shared_ptr<Kernel> sweep;
    int maxbs, warps, lpr;

    DeviceBlockGSKernels (shared_ptr<Device> device, int amaxbs, int alpr)
      : maxbs(amaxbs), lpr(alpr)
    {
      // group memory for the block residual, metal has 32 KB
      if (size_t(maxbs)*sizeof(T) > 24*1024)
        throw Exception ("DeviceBlockGaussSeidel: block size " + ToString(maxbs) + " too large");
      warps = min (8, (maxbs+31)/32);
      if (auto e = getenv("NGS_BGS_WARPS")) warps = atoi(e);

      string code = Substitute (kernel_source, "SCAL", is_same_v<T,double> ? "double" : "float");
      code = Substitute (code, "$WARPS", ToString(warps));
      code = Substitute (code, "$MAXBS", ToString(maxbs));
      code = Substitute (code, "$LPR", ToString(lpr));
      library = device->CompileSource (string(code_gpukernel) + code);
      sweep = library->GetKernel ("bgs_sweep");

      if (getenv("NGS_BGS_INFO"))
        cout << "bgs_sweep maxbs=" << maxbs << " warps=" << warps << " lanes/row=" << lpr
             << " " << sweep->Info(warps*32) << endl;
    }

    // lanes per csr row: a power of two, about a quarter of the average row length
    static int LanesPerRow (double avgnnz)
    {
      int lpr = 1;
      while (lpr < 32 && 4*lpr < avgnnz) lpr *= 2;
      if (auto e = getenv("NGS_BGS_LPR")) lpr = atoi(e);
      return lpr;
    }

    // maxbs rounded up to a power of two, so few variants get compiled
    static shared_ptr<const DeviceBlockGSKernels> Get (shared_ptr<Device> device, int maxbs, int lpr)
    {
      static mutex mtx;
      static map<tuple<Device*,int,int>, shared_ptr<const DeviceBlockGSKernels>> cache;

      if constexpr (is_same_v<T,double>)
        if (!device->HasFloat64())
          throw Exception("DeviceBlockGaussSeidel<double> on "+device->Name()+", which has no fp64");
      if (device->SimdWidth() != 32)
        throw Exception("DeviceBlockGaussSeidel needs a simd width of 32, "+device->Name()+" has "
                        + ToString(device->SimdWidth()));
      int cls = 32;
      while (cls < maxbs) cls *= 2;
      auto lock = lock_guard<mutex>(mtx);
      auto & entry = cache[{device.get(), cls, lpr}];
      if (!entry)
        entry = make_shared<DeviceBlockGSKernels> (device, cls, lpr);
      return entry;
    }
  };



  template <typename T>
  template <typename TM>
  DeviceBlockGaussSeidel<T> :: DeviceBlockGaussSeidel (const BlockJacobiPrecond<TM> & pre)
  {
    height = pre.Height();
    const Table<int> & blocktable = *pre.GetBlockTable();
    const Table<int> & coloring = pre.GetBlockColoring();
    nblocks = blocktable.Size();
    const auto & inverses = pre.GetInverses();

    device = GetGpuDevice();
    queue = device->DefaultQueue();
    memtype = PreferredMemType();
    devmat = make_shared<DeviceSparseMatrix<T>> (*pre.GetMatrix());

    // offsets into the concatenated index and matrix arrays
    Array<int> blockfirst(nblocks+1), matfirst(nblocks+1);
    size_t nind = 0, nmat = 0;
    maxbs = 0;
    for (size_t i = 0; i < nblocks; i++)
      {
        blockfirst[i] = int(nind);
        matfirst[i] = int(nmat);
        int bs = blocktable[i].Size();
        nind += bs;
        nmat += size_t(bs)*bs;
        maxbs = max (maxbs, bs);
        if (nmat > size_t(INT_MAX))
          throw Exception("DeviceBlockGaussSeidel: block inverses too large for 32-bit offsets");
      }
    blockfirst[nblocks] = int(nind);
    matfirst[nblocks] = int(nmat);

    FlatArray<int> indices = blocktable.AsArray();
    if (indices.Size() != nind)
      throw Exception("DeviceBlockGaussSeidel: block table is not contiguous");

    Array<T> mats(nmat);
    for (size_t i = 0; i < nblocks; i++)
      {
        size_t bs = blocktable[i].Size();
        T * dst = mats.Data() + matfirst[i];
        const auto & inv = inverses[i];
        for (size_t r = 0; r < bs; r++)
          for (size_t c = 0; c < bs; c++)
            dst[c*bs+r] = T(inv(r,c));
      }

    // non-empty blocks, color after color, small ones (one warp) before large ones
    Array<int> colorblocks;
    colorfirst.SetSize (coloring.Size()+1);
    colorsplit.SetSize (coloring.Size());
    for (size_t c = 0; c < coloring.Size(); c++)
      {
        colorfirst[c] = colorblocks.Size();
        for (int b : coloring[c])
          if (blocktable[b].Size() && blocktable[b].Size() <= 32) colorblocks.Append (b);
        colorsplit[c] = colorblocks.Size();
        for (int b : coloring[c])
          if (blocktable[b].Size() > 32) colorblocks.Append (b);
      }
    colorfirst[coloring.Size()] = colorblocks.Size();

    int lpr = DeviceBlockGSKernels<T>::LanesPerRow (double(pre.GetMatrix()->NZE()) / max<size_t>(height,1));
    kern_small = DeviceBlockGSKernels<T>::Get (device, 32, lpr);
    kern_large = (maxbs > 32) ? DeviceBlockGSKernels<T>::Get (device, maxbs, lpr) : nullptr;

    cout << IM(7) << "DeviceBlockGaussSeidel<" << (is_same_v<T,double> ? "double" : "float")
         << "> nblocks = " << nblocks << " (" << colorblocks.Size() << " non-empty), colors = " << coloring.Size()
         << ", max block " << maxbs << ", matrix entries = " << nmat << endl;

    dev_blockfirst  = device->NewBuffer<int> (nblocks+1, MemType::Device);
    dev_matfirst    = device->NewBuffer<int> (nblocks+1, MemType::Device);
    dev_indices     = device->NewBuffer<int> (max<size_t>(nind,1), MemType::Device);
    dev_mats        = device->NewBuffer<T> (max<size_t>(nmat,1), MemType::Device);
    dev_colorblocks = device->NewBuffer<int> (max<size_t>(colorblocks.Size(),1), MemType::Device);

    dev_blockfirst.H2D (blockfirst.Data(), nblocks+1);
    dev_matfirst.H2D (matfirst.Data(), nblocks+1);
    if (nind) dev_indices.H2D (indices.Data(), nind);
    if (nmat) dev_mats.H2D (mats.Data(), nmat);
    if (colorblocks.Size()) dev_colorblocks.H2D (colorblocks.Data(), colorblocks.Size());
  }


  template <typename T>
  void DeviceBlockGaussSeidel<T> :: Sweep (BaseVector & x, const BaseVector & b, bool backward) const
  {
    if (x.Size() != height || b.Size() != height)
      throw Exception("DeviceBlockGaussSeidel::Smooth - size mismatch");

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> ub(b, memtype);
    auto xarg = ux.DevArgRW(), barg = ub.DevArgRO();

    auto launch = [&] (const DeviceBlockGSKernels<T> & kern, int first, int nitems)
    {
      if (!nitems) return;
      queue->Launch (*kern.sweep, Dim3(nitems), Dim3(32*kern.warps),
                     { KernelArg(dev_colorblocks, first), KernelArg(dev_blockfirst), KernelArg(dev_indices),
                       KernelArg(dev_matfirst), KernelArg(dev_mats),
                       KernelArg(devmat->DevFirstI()), KernelArg(devmat->DevColNr()), KernelArg(devmat->DevValues()),
                       barg, xarg, KernelArg(nitems) });
    };

    int ncolors = colorfirst.Size()-1;
    for (int ci = 0; ci < ncolors; ci++)
      {
        int c = backward ? ncolors-1-ci : ci;
        launch (*kern_small, colorfirst[c], colorsplit[c]-colorfirst[c]);
        if (kern_large)
          launch (*kern_large, colorsplit[c], colorfirst[c+1]-colorsplit[c]);
      }
  }

  template <typename T>
  void DeviceBlockGaussSeidel<T> :: Smooth (BaseVector & x, const BaseVector & b, int steps) const
  {
    static Timer t("DeviceBlockGaussSeidel::Smooth"); RegionTimer reg(t);
    for (int k = 0; k < steps; k++)
      Sweep (x, b, false);
  }

  template <typename T>
  void DeviceBlockGaussSeidel<T> :: SmoothBack (BaseVector & x, const BaseVector & b, int steps) const
  {
    static Timer t("DeviceBlockGaussSeidel::SmoothBack"); RegionTimer reg(t);
    for (int k = 0; k < steps; k++)
      Sweep (x, b, true);
  }

  template <typename T>
  void DeviceBlockGaussSeidel<T> :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceBlockGaussSeidel::Mult"); RegionTimer reg(t);
    y = 0.0;
    Sweep (y, x, false);
    Sweep (y, x, true);
  }

  template <typename T>
  void DeviceBlockGaussSeidel<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    auto tmp = CreateColVector();
    Mult (x, *tmp);
    y += s * *tmp;
  }


  template <typename T>
  BaseMatrix::OperatorInfo DeviceBlockGaussSeidel<T> :: GetOperatorInfo () const
  {
    return { string("DeviceBlockGaussSeidel<") + (is_same_v<T,double> ? "double" : "float")
             + "> (blocks=" + ToString(nblocks) + ", colors=" + ToString(colorfirst.Size()-1) + ")",
             height, height };
  }

  template <typename T>
  ostream & DeviceBlockGaussSeidel<T> :: Print (ostream & ost) const
  {
    ost << "DeviceBlockGaussSeidel<" << (is_same_v<T,double> ? "double" : "float")
        << ">, height = " << height << ", nblocks = " << nblocks << ", colors = " << colorfirst.Size()-1
        << ", max block " << maxbs << ", on " << device->Name() << endl;
    return ost;
  }


  template <typename TM>
  shared_ptr<BaseMatrix> CreateDeviceBlockGaussSeidelT (const BlockJacobiPrecond<TM> & pre)
  {
    if (!ngs_gpu::HasDevice() || GetGpuDevice()->SimdWidth() != 32)
      return nullptr;
    if constexpr (is_same_v<TM,double>)
      if (GetGpuDevice()->HasFloat64())
        return make_shared<DeviceBlockGaussSeidel<double>> (pre);
    return make_shared<DeviceBlockGaussSeidel<float>> (pre);
  }

  shared_ptr<BaseMatrix> CreateDeviceBlockGaussSeidel (const BlockJacobiPrecond<double> & pre)
  { return CreateDeviceBlockGaussSeidelT (pre); }
  shared_ptr<BaseMatrix> CreateDeviceBlockGaussSeidel (const BlockJacobiPrecond<float> & pre)
  { return CreateDeviceBlockGaussSeidelT (pre); }

  template class DeviceBlockGaussSeidel<double>;
  template class DeviceBlockGaussSeidel<float>;
  template DeviceBlockGaussSeidel<double>::DeviceBlockGaussSeidel (const BlockJacobiPrecond<double> &);
  template DeviceBlockGaussSeidel<double>::DeviceBlockGaussSeidel (const BlockJacobiPrecond<float> &);
  template DeviceBlockGaussSeidel<float>::DeviceBlockGaussSeidel (const BlockJacobiPrecond<double> &);
  template DeviceBlockGaussSeidel<float>::DeviceBlockGaussSeidel (const BlockJacobiPrecond<float> &);
}
