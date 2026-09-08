/*********************************************************************/
/* File:   device_ebeconst.cpp                                       */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   7. Sep. 2026                                              */
/*********************************************************************/

#define FILE_DEVICE_EBECONST_CPP

#include <la.hpp>
#include <gpukernel.hpp>
#include <tinybla.hpp>

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

    template <int N> int RoundUp (int i) { return N * ((i+N-1) / N); }

    /*
      Both kernels see the element matrix as B (nin x nout, row-major with
      leading dimension ldb, zero padded to multiples of 8): B = M^T for
      Mult, B = M for MultTrans, so y[out[r]] += sum_c B(c,r) x[in[c]].

      Lane kernel, any backend: lanes consecutive work-items share a block
      and take its output entries in turns. beta 0 stores, otherwise
      accumulates; the atomic variant accumulates only (output dofs shared
      between blocks).
    */
    const char * lane_source = R"RAW(

      #define EBE_BODY(ACCUMULATE)                                                  \
        int lane = int(LOCAL_ID_X) & (lanes-1);                                     \
        int block = int(GLOBAL_ID_X) / lanes;                                       \
        if (block >= nblocks) return;                                               \
        GLOBAL_PTR(const int) bin  = dnin + block*nin;                              \
        GLOBAL_PTR(const int) bout = dnout + block*nout;                            \
        for (int r = lane; r < nout; r += lanes)                                    \
          {                                                                         \
            SCAL v = 0;                                                             \
            for (int c = 0; c < nin; c++) v += bmat[c*ldb+r] * x[bin[c]];           \
            int d = bout[r];                                                        \
            ACCUMULATE;                                                             \
          }

      KERNEL(ebe_mult, GLOBAL_IN(SCAL,bmat), GLOBAL_IN(int,dnin), GLOBAL_IN(int,dnout),
                       GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                       VALUE(SCAL,s), VALUE(SCAL,beta), VALUE(int,nin), VALUE(int,nout),
                       VALUE(int,ldb), VALUE(int,lanes), VALUE(int,nblocks))
      {
        EBE_BODY(y[d] = (beta == SCAL(0)) ? s*v : beta*y[d] + s*v)
      }

      KERNEL(ebe_mult_atomic, GLOBAL_IN(SCAL,bmat), GLOBAL_IN(int,dnin), GLOBAL_IN(int,dnout),
                              GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y),
                              VALUE(SCAL,s), VALUE(SCAL,beta), VALUE(int,nin), VALUE(int,nout),
                              VALUE(int,ldb), VALUE(int,lanes), VALUE(int,nblocks))
      {
        EBE_BODY(ATOMIC_ADD(&y[d], s*v))
      }

    )RAW";


    /*
      Gemm kernel, simd width 32: a group takes BS_ELS blocks. Their x
      entries are gathered into group memory (BS_ELS x NIN8), multiplied
      with B in 8x8 warp tiles (tinybla WarpMatrix, tensor units where
      available), and scattered from group memory into y. The shape is
      compiled in: one library per (scalar, nin, nout, index layout);
      consecutive dofs (block*n + j) skip the index table.
    */
    const char * gemm_source = R"RAW(
      using namespace tinybla;

      #define EBE_GEMM_BODY(ACCUMULATE)                                                 \
        constexpr uint nin = $NIN, nout = $NOUT, nin8 = $NIN8, nout8 = $NOUT8;          \
        constexpr uint bs_els = $BS_ELS, warps = $WARPS;                                \
        constexpr uint el_tiles = bs_els/8, out_tiles = nout8/8;                        \
        uint tid = LOCAL_ID_X;                                                          \
        uint bdim = GROUP_SIZE_X;                                                       \
        uint warp = tid/32;                                                             \
        uint baseel = GROUP_ID_X*bs_els;                                                \
        SHARED_2D(SCAL, elx, bs_els, nin8);                                             \
        SHARED_2D(SCAL, ely, bs_els, nout8);                                            \
        for (uint i = tid; i < bs_els*nin8; i += bdim)                                  \
          {                                                                             \
            uint r = i / nin8, c = i % nin8, el = baseel + r;                           \
            elx[r][c] = (c < nin && el < uint(nblocks)) ? x[$IN_INDEX] : SCAL(0);       \
          }                                                                             \
        BARRIER();                                                                      \
        auto mat_elx = MakeBareMatrix<RowMajor>(elx);                                   \
        auto mat_ely = MakeBareMatrix<RowMajor>(ely);                                   \
        auto mat_b = MakeBareMatrix<RowMajor>(bmat, nout8);                             \
        for (uint t = warp; t < el_tiles*out_tiles; t += warps)                         \
          {                                                                             \
            uint et = t % el_tiles, ot = t / el_tiles;                                  \
            WarpMatrix<8,8,SCAL> sum = 0;                                               \
            sum.AddMM<nin8> (mat_elx.SubMatrix(8*et,0), mat_b.SubMatrix(0,8*ot), tid);  \
            sum.Store (mat_ely.SubMatrix(8*et,8*ot), tid);                              \
          }                                                                             \
        BARRIER();                                                                      \
        for (uint i = tid; i < bs_els*nout; i += bdim)                                  \
          {                                                                             \
            uint r = i / nout, c = i % nout, el = baseel + r;                           \
            if (el < uint(nblocks))                                                     \
              {                                                                         \
                int d = $OUT_INDEX;                                                     \
                SCAL v = ely[r][c];                                                     \
                ACCUMULATE;                                                             \
              }                                                                         \
          }

      KERNEL(ebe_gemm, GLOBAL_IN(SCAL,bmat), GLOBAL_IN(int,dnin), GLOBAL_IN(int,dnout),
                       GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                       VALUE(SCAL,s), VALUE(SCAL,beta), VALUE(int,nblocks))
      {
        EBE_GEMM_BODY(y[d] = (beta == SCAL(0)) ? s*v : beta*y[d] + s*v)
      }

      KERNEL(ebe_gemm_atomic, GLOBAL_IN(SCAL,bmat), GLOBAL_IN(int,dnin), GLOBAL_IN(int,dnout),
                              GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y),
                              VALUE(SCAL,s), VALUE(SCAL,beta), VALUE(int,nblocks))
      {
        EBE_GEMM_BODY(ATOMIC_ADD(&y[d], s*v))
      }

    )RAW";


    template <typename T>
    class DeviceEBEKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<ngs_gpu::Queue> queue;
      shared_ptr<Kernel> mult, mult_atomic;
      unsigned groupsize;

      DeviceEBEKernels (shared_ptr<Device> adevice)
        : device(adevice)
      {
        if constexpr (is_same_v<T,double>)
          if (!device->HasFloat64())
            throw Exception("DeviceConstantEBEMatrix<double> on "+device->Name()+
                            ", which has no fp64 - use DeviceConstantEBEMatrix<float>");

        string scal = is_same_v<T,double> ? "double" : "float";
        library = device->CompileSource (string(code_gpukernel) +
                                         Substitute (lane_source, "SCAL", scal));
        mult        = library->GetKernel ("ebe_mult");
        mult_atomic = library->GetKernel ("ebe_mult_atomic");
        queue = device->DefaultQueue();
        groupsize = (device->SimdWidth() > 1) ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
      }

      // a power of two up to the simd width, matching the number of output rows
      int ChooseLanes (int nout) const
      {
        size_t simd = device->SimdWidth();
        if (simd <= 1) return 1;
        int lanes = 1;
        while (lanes < nout && size_t(2*lanes) <= simd && size_t(2*lanes) <= groupsize)
          lanes *= 2;
        return lanes;
      }

      static const DeviceEBEKernels & Get()
      {
        static mutex mtx;
        static shared_ptr<DeviceEBEKernels> cached;

        auto dev = GetGpuDevice();
        auto lock = lock_guard<mutex>(mtx);
        if (!cached || cached->device != dev)
          cached = make_shared<DeviceEBEKernels> (dev);
        return *cached;
      }
    };

    // dofs of block i are i*n .. i*n+n-1
    bool IsConsecutive (FlatTable<int> table, int n)
    {
      for (size_t i = 0; i < table.Size(); i++)
        for (int j = 0; j < n; j++)
          if (table[i][j] != int(i*n+j)) return false;
      return true;
    }

    // the dof tables are rectangular by construction, the kernel relies on it
    template <typename T>
    void UploadTable (FlatTable<int> table, int rowsize, shared_ptr<Device> device,
                      TypedBuffer<int> & buf, const char * name)
    {
      std::vector<int> flat;
      flat.reserve (table.Size()*rowsize);
      for (size_t i = 0; i < table.Size(); i++)
        {
          if (int(table[i].Size()) != rowsize)
            throw Exception (string("DeviceConstantEBEMatrix: ") + name + " table is not rectangular");
          for (int d : table[i]) flat.push_back (d);
        }
      buf = device->template NewBuffer<int> (max<size_t>(flat.size(),1), MemType::Device);
      buf.H2D (flat.data(), flat.size());
    }
  }


  template <typename T>
  class DeviceEBEGemmKernels
  {
    shared_ptr<Library> library;
  public:
    shared_ptr<Kernel> mult, mult_atomic;
    int nin, nout, bs_els, warps;
    bool in_consecutive, out_consecutive;

    // group memory for the batch, metal has 32 KB
    static int BatchSize (int nin8, int nout8)
    {
      int bs = 64;
      while (bs > 8 && size_t(bs)*(nin8+nout8)*sizeof(T) > 24*1024) bs /= 2;
      return bs;
    }

    static bool Supported (const Device & device, int nin, int nout)
    {
      if (getenv("NGS_EBE_LANES")) return false;   // force the lane kernel
      return device.SimdWidth() == 32 && nin > 0 && nout > 0
        && size_t(8)*(RoundUp<8>(nin)+RoundUp<8>(nout))*sizeof(T) <= 24*1024;
    }

    DeviceEBEGemmKernels (shared_ptr<Device> device, int anin, int anout, bool ain_consec, bool aout_consec)
      : nin(anin), nout(anout), in_consecutive(ain_consec), out_consecutive(aout_consec)
    {
      int nin8 = RoundUp<8>(nin), nout8 = RoundUp<8>(nout);
      bs_els = BatchSize (nin8, nout8);
      if (auto e = getenv("NGS_EBE_BS")) bs_els = atoi(e);
      int tiles = (bs_els/8) * (nout8/8);
      warps = min (tiles, 8);
      if (auto e = getenv("NGS_EBE_WARPS")) warps = atoi(e);

      string code = Substitute (gemm_source, "SCAL", is_same_v<T,double> ? "double" : "float");
      code = Substitute (code, "$NIN8", ToString(nin8));
      code = Substitute (code, "$NOUT8", ToString(nout8));
      code = Substitute (code, "$NIN", ToString(nin));
      code = Substitute (code, "$NOUT", ToString(nout));
      code = Substitute (code, "$BS_ELS", ToString(bs_els));
      code = Substitute (code, "$WARPS", ToString(warps));
      code = Substitute (code, "$IN_INDEX", in_consecutive ? "el*nin+c" : "dnin[el*nin+c]");
      code = Substitute (code, "$OUT_INDEX", out_consecutive ? "int(el*nout+c)" : "dnout[el*nout+c]");

      library = device->CompileSource (string(code_gpukernel) + code_tinybla + code);
      mult        = library->GetKernel ("ebe_gemm");
      mult_atomic = library->GetKernel ("ebe_gemm_atomic");

      if (getenv("NGS_EBE_INFO"))
        cout << "ebe_gemm " << nin << "x" << nout << ": bs_els=" << bs_els << " warps=" << warps
             << " consecutive in/out=" << in_consecutive << "/" << out_consecutive
             << " " << mult->Info(warps*32) << endl;
    }

    // null if the gemm kernel is not available for this shape on this device
    static shared_ptr<const DeviceEBEGemmKernels> Get (shared_ptr<Device> device, int nin, int nout,
                                                        bool in_consec, bool out_consec)
    {
      static mutex mtx;
      static map<tuple<Device*,int,int,bool,bool>, shared_ptr<const DeviceEBEGemmKernels>> cache;

      if (!Supported (*device, nin, nout)) return nullptr;
      auto lock = lock_guard<mutex>(mtx);
      auto & entry = cache[{device.get(), nin, nout, in_consec, out_consec}];
      if (!entry)
        entry = make_shared<DeviceEBEGemmKernels> (device, nin, nout, in_consec, out_consec);
      return entry;
    }
  };



  template <typename T>
  template <typename TM>
  DeviceConstantEBEMatrix<T> :: DeviceConstantEBEMatrix (const ConstantElementByElementMatrix<TM> & mat)
    : memtype (PreferredMemType())
  {
    const auto & kern = DeviceEBEKernels<T>::Get();
    device = kern.device;
    queue = kern.queue;

    height = mat.Height();
    width = mat.Width();
    FlatMatrix<TM> m = mat.GetMatrix();
    hm = m.Height();
    wm = m.Width();
    nblocks = mat.GetRowDNums().Size();

    int hm8 = RoundUp<8>(hm), wm8 = RoundUp<8>(wm);
    std::vector<T> bt (size_t(wm8)*hm8, T(0)), b (size_t(hm8)*wm8, T(0));
    for (int i = 0; i < hm; i++)
      for (int j = 0; j < wm; j++)
        {
          bt[size_t(j)*hm8+i] = T(m(i,j));
          b[size_t(i)*wm8+j] = T(m(i,j));
        }
    dev_mat = device->template NewBuffer<T> (max<size_t>(bt.size(),1), MemType::Device);
    dev_mat.H2D (bt.data(), bt.size());
    dev_mat_trans = device->template NewBuffer<T> (max<size_t>(b.size(),1), MemType::Device);
    dev_mat_trans.H2D (b.data(), b.size());

    UploadTable<T> (mat.GetRowDNums(), wm, device, dev_rowdnums, "row");
    UploadTable<T> (mat.GetColDNums(), hm, device, dev_coldnums, "col");

    disjoint_rows = (mat.GetRowColoring().Size() == 0);
    disjoint_cols = (mat.GetColColoring().Size() == 0);
    onto_cols = disjoint_cols && mat.GetColDNums().AsArray().Size() == height;
    onto_rows = disjoint_rows && mat.GetRowDNums().AsArray().Size() == width;

    lanes = kern.ChooseLanes (hm);
    lanes_trans = kern.ChooseLanes (wm);
    bool rows_consec = IsConsecutive (mat.GetRowDNums(), wm);
    bool cols_consec = IsConsecutive (mat.GetColDNums(), hm);
    gemm = DeviceEBEGemmKernels<T>::Get (device, wm, hm, rows_consec, cols_consec);
    gemm_trans = DeviceEBEGemmKernels<T>::Get (device, hm, wm, cols_consec, rows_consec);
  }


  template <typename T>
  void DeviceConstantEBEMatrix<T> :: Launch (const BaseVector & x, BaseVector & y, T s, T beta, bool trans) const
  {
    if (x.Size() != (trans ? height : width) || y.Size() != (trans ? width : height))
      throw Exception("DeviceConstantEBEMatrix::Mult - size mismatch");
    if (nblocks == 0) return;

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);

    bool atomic = trans ? !disjoint_rows : !disjoint_cols;
    const auto & bmat = trans ? dev_mat_trans : dev_mat;
    const auto & dnin = trans ? dev_coldnums : dev_rowdnums;
    const auto & dnout = trans ? dev_rowdnums : dev_coldnums;
    auto yarg = (beta == T(0) && !atomic) ? uy.DevArgW() : uy.DevArgRW();

    if (const auto & g = trans ? gemm_trans : gemm)
      {
        unsigned groups = (nblocks + g->bs_els-1) / g->bs_els;
        queue->Launch (atomic ? *g->mult_atomic : *g->mult,
                       Dim3(groups), Dim3(g->warps*32),
                       { KernelArg(bmat), KernelArg(dnin), KernelArg(dnout),
                         ux.DevArgRO(), yarg,
                         KernelArg(s), KernelArg(beta), KernelArg(int(nblocks)) });
        return;
      }

    const auto & kern = DeviceEBEKernels<T>::Get();
    int nin = trans ? hm : wm, nout = trans ? wm : hm;
    int ldb = RoundUp<8>(nout);
    int l = trans ? lanes_trans : lanes;
    size_t items = nblocks * l;
    unsigned groups = (items + kern.groupsize-1) / kern.groupsize;
    queue->Launch (atomic ? *kern.mult_atomic : *kern.mult,
                   Dim3(groups), Dim3(kern.groupsize),
                   { KernelArg(bmat), KernelArg(dnin), KernelArg(dnout),
                     ux.DevArgRO(), yarg,
                     KernelArg(s), KernelArg(beta), KernelArg(nin), KernelArg(nout),
                     KernelArg(ldb), KernelArg(l), KernelArg(int(nblocks)) });
  }

  template <typename T>
  void DeviceConstantEBEMatrix<T> :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceConstantEBEMatrix::Mult"); RegionTimer reg(t);
    if (onto_cols) Launch (x, y, T(1), T(0), false);
    else { y = 0.0; Launch (x, y, T(1), T(1), false); }
  }

  template <typename T>
  void DeviceConstantEBEMatrix<T> :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceConstantEBEMatrix::MultTrans"); RegionTimer reg(t);
    if (onto_rows) Launch (x, y, T(1), T(0), true);
    else { y = 0.0; Launch (x, y, T(1), T(1), true); }
  }

  template <typename T>
  void DeviceConstantEBEMatrix<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceConstantEBEMatrix::MultAdd"); RegionTimer reg(t);
    Launch (x, y, T(s), T(1), false);
  }

  template <typename T>
  void DeviceConstantEBEMatrix<T> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceConstantEBEMatrix::MultTransAdd"); RegionTimer reg(t);
    Launch (x, y, T(s), T(1), true);
  }




  template <typename T>
  BaseMatrix::OperatorInfo DeviceConstantEBEMatrix<T> :: GetOperatorInfo () const
  {
    return { string("DeviceConstantEBEMatrix<") + (is_same_v<T,double> ? "double" : "float")
             + "> (blocks=" + ToString(nblocks) + ", " + ToString(hm) + "x" + ToString(wm) + ")",
             height, width };
  }

  template <typename T>
  ostream & DeviceConstantEBEMatrix<T> :: Print (ostream & ost) const
  {
    ost << "DeviceConstantEBEMatrix<" << (is_same_v<T,double> ? "double" : "float")
        << ">, " << height << " x " << width << ", blocks = " << nblocks
        << ", element matrix " << hm << " x " << wm << ", on " << device->Name()
        << (gemm ? ", gemm kernel" : ", lane kernel") << endl;
    return ost;
  }


  template <typename SCAL>
  shared_ptr<BaseMatrix> ConstantElementByElementMatrix<SCAL> :: CreateDeviceMatrix () const
  {
    if constexpr (is_same_v<SCAL,double> || is_same_v<SCAL,float>)
      if (ngs_gpu::HasDevice())
        {
          if constexpr (is_same_v<SCAL,double>)
            if (GetGpuDevice()->HasFloat64())
              return make_shared<DeviceConstantEBEMatrix<double>> (*this);
          return make_shared<DeviceConstantEBEMatrix<float>> (*this);
        }
    return BaseMatrix::CreateDeviceMatrix();
  }

  template shared_ptr<BaseMatrix> ConstantElementByElementMatrix<double>::CreateDeviceMatrix () const;
  template shared_ptr<BaseMatrix> ConstantElementByElementMatrix<float>::CreateDeviceMatrix () const;
  template shared_ptr<BaseMatrix> ConstantElementByElementMatrix<Complex>::CreateDeviceMatrix () const;

  template class DeviceConstantEBEMatrix<double>;
  template class DeviceConstantEBEMatrix<float>;
  template DeviceConstantEBEMatrix<double>::DeviceConstantEBEMatrix (const ConstantElementByElementMatrix<double>&);
  template DeviceConstantEBEMatrix<double>::DeviceConstantEBEMatrix (const ConstantElementByElementMatrix<float>&);
  template DeviceConstantEBEMatrix<float>::DeviceConstantEBEMatrix (const ConstantElementByElementMatrix<double>&);
  template DeviceConstantEBEMatrix<float>::DeviceConstantEBEMatrix (const ConstantElementByElementMatrix<float>&);
}
