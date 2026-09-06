/*********************************************************************/
/* File:   device_ebeconst.cpp                                       */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   7. Sep. 2026                                              */
/*********************************************************************/

#define FILE_DEVICE_EBECONST_CPP

#include <la.hpp>
#include <gpukernel.hpp>

namespace ngla
{
  using namespace ngs_gpu;

  namespace
  {
    string SubstituteScal (string src, const string & value)
    {
      const string token = "SCAL";
      for (size_t p = src.find(token); p != string::npos; p = src.find(token, p))
        src.replace (p, token.size(), value);
      return src;
    }

    /*
      lanes consecutive work-items share a block and take its output rows
      in turns: row r of the block is sum_c mat(r,c) * x[in[c]], written to
      y[out[r]]. trans swaps the dof tables and reads mat column-wise.
      beta 0 stores, otherwise accumulates; the atomic variant accumulates
      only (output dofs shared between blocks).
    */
    const char * kernel_source = R"RAW(

      #define EBE_BODY(ACCUMULATE)                                                  \
        int lane = int(LOCAL_ID_X) & (lanes-1);                                     \
        int block = int(GLOBAL_ID_X) / lanes;                                       \
        if (block >= nblocks) return;                                               \
        int nin  = trans ? hm : wm;                                                 \
        int nout = trans ? wm : hm;                                                 \
        GLOBAL_PTR(const int) in  = trans ? coldnums + block*hm : rowdnums + block*wm; \
        GLOBAL_PTR(const int) out = trans ? rowdnums + block*wm : coldnums + block*hm; \
        for (int r = lane; r < nout; r += lanes)                                    \
          {                                                                         \
            SCAL sum = 0;                                                           \
            if (trans)                                                              \
              for (int c = 0; c < nin; c++) sum += mat[c*wm+r] * x[in[c]];          \
            else                                                                    \
              for (int c = 0; c < nin; c++) sum += mat[r*wm+c] * x[in[c]];          \
            ACCUMULATE;                                                             \
          }

      KERNEL(ebe_mult, GLOBAL_IN(SCAL,mat), GLOBAL_IN(int,rowdnums), GLOBAL_IN(int,coldnums),
                       GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                       VALUE(SCAL,s), VALUE(SCAL,beta), VALUE(int,hm), VALUE(int,wm),
                       VALUE(int,lanes), VALUE(int,trans), VALUE(int,nblocks))
      {
        EBE_BODY(y[out[r]] = (beta == SCAL(0)) ? s*sum : beta*y[out[r]] + s*sum)
      }

      KERNEL(ebe_mult_atomic, GLOBAL_IN(SCAL,mat), GLOBAL_IN(int,rowdnums), GLOBAL_IN(int,coldnums),
                              GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y),
                              VALUE(SCAL,s), VALUE(SCAL,beta), VALUE(int,hm), VALUE(int,wm),
                              VALUE(int,lanes), VALUE(int,trans), VALUE(int,nblocks))
      {
        EBE_BODY(ATOMIC_ADD(&y[out[r]], s*sum))
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
                                         SubstituteScal (kernel_source, scal));
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

    std::vector<T> hmat (size_t(hm)*wm);
    for (int i = 0; i < hm; i++)
      for (int j = 0; j < wm; j++)
        hmat[size_t(i)*wm+j] = T(m(i,j));
    dev_mat = device->template NewBuffer<T> (max<size_t>(hmat.size(),1), MemType::Device);
    dev_mat.H2D (hmat.data(), hmat.size());

    UploadTable<T> (mat.GetRowDNums(), wm, device, dev_rowdnums, "row");
    UploadTable<T> (mat.GetColDNums(), hm, device, dev_coldnums, "col");

    disjoint_rows = (mat.GetRowColoring().Size() == 0);
    disjoint_cols = (mat.GetColColoring().Size() == 0);
    onto_cols = disjoint_cols && mat.GetColDNums().AsArray().Size() == height;
    onto_rows = disjoint_rows && mat.GetRowDNums().AsArray().Size() == width;

    lanes = kern.ChooseLanes (hm);
    lanes_trans = kern.ChooseLanes (wm);
  }


  template <typename T>
  void DeviceConstantEBEMatrix<T> :: Launch (const BaseVector & x, BaseVector & y, T s, T beta, bool trans) const
  {
    if (x.Size() != (trans ? height : width) || y.Size() != (trans ? width : height))
      throw Exception("DeviceConstantEBEMatrix::Mult - size mismatch");
    if (nblocks == 0) return;

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);

    const auto & kern = DeviceEBEKernels<T>::Get();
    bool atomic = trans ? !disjoint_rows : !disjoint_cols;
    int l = trans ? lanes_trans : lanes;
    size_t items = nblocks * l;
    unsigned groups = (items + kern.groupsize-1) / kern.groupsize;
    queue->Launch (atomic ? *kern.mult_atomic : *kern.mult,
                   Dim3(groups), Dim3(kern.groupsize),
                   { KernelArg(dev_mat), KernelArg(dev_rowdnums), KernelArg(dev_coldnums),
                     ux.DevArgRO(), (beta == T(0) && !atomic) ? uy.DevArgW() : uy.DevArgRW(),
                     KernelArg(s), KernelArg(beta), KernelArg(hm), KernelArg(wm),
                     KernelArg(l), KernelArg(int(trans)), KernelArg(int(nblocks)) });
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
  AutoVector DeviceConstantEBEMatrix<T> :: CreateRowVector () const
  {
    return make_unique<DeviceVector<T>> (width, memtype);
  }

  template <typename T>
  AutoVector DeviceConstantEBEMatrix<T> :: CreateColVector () const
  {
    return make_unique<DeviceVector<T>> (height, memtype);
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
        << ", element matrix " << hm << " x " << wm << ", on " << device->Name() << endl;
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
