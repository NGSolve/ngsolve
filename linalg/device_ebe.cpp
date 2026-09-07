#define FILE_DEVICE_EBE_CPP
#include <la.hpp>
#include <gpukernel.hpp>
#include <climits>

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

    const char * kernel_source = R"RAW(

      KERNEL(ebe_multadd, GLOBAL_IN(int,rowfirst), GLOBAL_IN(int,colfirst), GLOBAL_IN(int,matfirst),
                          GLOBAL_IN(int,rowidx), GLOBAL_IN(int,colidx), GLOBAL_IN(SCAL,mats),
                          GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y),
                          VALUE(SCAL,s), VALUE(int,lanes), VALUE(int,trans), VALUE(int,nel))
      {
        int lane = int(LOCAL_ID_X) & (lanes-1);
        int el = int(GLOBAL_ID_X) / lanes;
        if (el >= nel) return;
        int r0 = rowfirst[el], nr = rowfirst[el+1] - r0;
        int c0 = colfirst[el], nc = colfirst[el+1] - c0;
        GLOBAL_PTR(const SCAL) mat = mats + matfirst[el];
        GLOBAL_PTR(const int) in  = trans ? rowidx + r0 : colidx + c0;
        GLOBAL_PTR(const int) out = trans ? colidx + c0 : rowidx + r0;
        int nin = trans ? nr : nc, nout = trans ? nc : nr;
        for (int r = lane; r < nout; r += lanes)
          {
            SCAL sum = 0;
            if (trans)
              for (int c = 0; c < nin; c++) sum += mat[c*nc+r] * x[in[c]];
            else
              for (int c = 0; c < nin; c++) sum += mat[r*nc+c] * x[in[c]];
            ATOMIC_ADD(&y[out[r]], s*sum);
          }
      }

    )RAW";

    template <typename T>
    class DeviceEBEKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<ngs_gpu::Queue> queue;
      shared_ptr<Kernel> multadd;
      unsigned groupsize;

      DeviceEBEKernels (shared_ptr<Device> adevice)
        : device(adevice)
      {
        if constexpr (is_same_v<T,double>)
          if (!device->HasFloat64())
            throw Exception("DeviceEBEMatrix<double> on "+device->Name()+
                            ", which has no fp64 - use DeviceEBEMatrix<float>");
        string scal = is_same_v<T,double> ? "double" : "float";
        library = device->CompileSource (string(code_gpukernel) +
                                         SubstituteScal (kernel_source, scal));
        multadd = library->GetKernel ("ebe_multadd");
        queue = device->DefaultQueue();
        groupsize = (device->SimdWidth() > 1) ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
      }

      // a power of two up to the simd width, matching the average number of output rows
      int ChooseLanes (double nout) const
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

    // reserved but unused element slots carry a -1 entry, as in the host MultAdd
    template <typename TM>
    bool ElementUsed (const ElementByElementMatrix<TM> & mat, size_t i)
    {
      auto rdi = mat.GetElementRowDNums(i);
      auto cdi = mat.GetElementColumnDNums(i);
      return rdi.Size() && cdi.Size() && rdi[0] != -1 && cdi[0] != -1;
    }
  }


  template <typename T>
  template <typename TM>
  DeviceEBEMatrix<T> :: DeviceEBEMatrix (const ElementByElementMatrix<TM> & mat)
    : memtype (PreferredMemType())
  {
    const auto & kern = DeviceEBEKernels<T>::Get();
    device = kern.device;
    queue = kern.queue;
    height = mat.Height();
    width = mat.Width();

    // offsets into the concatenated dof and matrix arrays, used elements only
    size_t nrow = 0, ncol = 0, nmat = 0;
    nel = 0;
    for (size_t i = 0; i < mat.GetNumElMats(); i++)
      if (ElementUsed (mat, i))
        {
          nel++;
          size_t nr = mat.GetElementRowDNums(i).Size(), nc = mat.GetElementColumnDNums(i).Size();
          nrow += nr;
          ncol += nc;
          nmat += nr*nc;
          if (nmat > size_t(INT_MAX))
            throw Exception("DeviceEBEMatrix: element matrices too large for 32-bit offsets");
        }

    Array<int> rowfirst(nel+1), colfirst(nel+1), matfirst(nel+1);
    Array<int> rowidx(nrow), colidx(ncol);
    Array<T> mats(nmat);
    size_t el = 0;
    rowfirst[0] = colfirst[0] = matfirst[0] = 0;
    for (size_t i = 0; i < mat.GetNumElMats(); i++)
      {
        if (!ElementUsed (mat, i)) continue;
        auto rdi = mat.GetElementRowDNums(i);
        auto cdi = mat.GetElementColumnDNums(i);
        auto m = mat.GetElementMatrix(i);
        if (m.Height() != rdi.Size() || m.Width() != cdi.Size())
          throw Exception("DeviceEBEMatrix: element matrix does not match its dof lists");
        rowidx.Range(rowfirst[el], rowfirst[el]+rdi.Size()) = rdi;
        colidx.Range(colfirst[el], colfirst[el]+cdi.Size()) = cdi;
        T * dst = mats.Data() + matfirst[el];
        for (size_t r = 0; r < rdi.Size(); r++)
          for (size_t c = 0; c < cdi.Size(); c++)
            dst[r*cdi.Size()+c] = T(m(r,c));
        rowfirst[el+1] = rowfirst[el] + int(rdi.Size());
        colfirst[el+1] = colfirst[el] + int(cdi.Size());
        matfirst[el+1] = matfirst[el] + int(rdi.Size()*cdi.Size());
        el++;
      }

    lanes = kern.ChooseLanes (nel ? double(nrow)/nel : 1);
    lanes_trans = kern.ChooseLanes (nel ? double(ncol)/nel : 1);
    cout << IM(7) << "DeviceEBEMatrix<" << (is_same_v<T,double> ? "double" : "float")
         << "> elements = " << nel << ", matrix entries = " << nmat
         << ", lanes = " << lanes << "/" << lanes_trans << endl;

    dev_rowfirst = device->NewBuffer<int> (nel+1, MemType::Device);
    dev_colfirst = device->NewBuffer<int> (nel+1, MemType::Device);
    dev_matfirst = device->NewBuffer<int> (nel+1, MemType::Device);
    dev_rowidx   = device->NewBuffer<int> (max<size_t>(nrow,1), MemType::Device);
    dev_colidx   = device->NewBuffer<int> (max<size_t>(ncol,1), MemType::Device);
    dev_mats     = device->NewBuffer<T> (max<size_t>(nmat,1), MemType::Device);

    dev_rowfirst.H2D (rowfirst.Data(), nel+1);
    dev_colfirst.H2D (colfirst.Data(), nel+1);
    dev_matfirst.H2D (matfirst.Data(), nel+1);
    dev_rowidx.H2D (rowidx.Data(), nrow);
    dev_colidx.H2D (colidx.Data(), ncol);
    dev_mats.H2D (mats.Data(), nmat);
  }

  template <typename T>
  void DeviceEBEMatrix<T> :: Launch (const BaseVector & x, BaseVector & y, T s, bool trans) const
  {
    if (x.Size() != (trans ? height : width) || y.Size() != (trans ? width : height))
      throw Exception("DeviceEBEMatrix::MultAdd - size mismatch");
    if (nel == 0) return;
    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    const auto & kern = DeviceEBEKernels<T>::Get();
    int l = trans ? lanes_trans : lanes;
    size_t items = nel * l;
    unsigned groups = (items + kern.groupsize-1) / kern.groupsize;
    queue->Launch (*kern.multadd, Dim3(groups), Dim3(kern.groupsize),
                   { KernelArg(dev_rowfirst), KernelArg(dev_colfirst), KernelArg(dev_matfirst),
                     KernelArg(dev_rowidx), KernelArg(dev_colidx), KernelArg(dev_mats),
                     ux.DevArgRO(), uy.DevArgRW(),
                     KernelArg(s), KernelArg(l), KernelArg(int(trans)), KernelArg(int(nel)) });
  }

  template <typename T>
  void DeviceEBEMatrix<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceEBEMatrix::MultAdd"); RegionTimer reg(t);
    Launch (x, y, T(s), false);
  }

  template <typename T>
  void DeviceEBEMatrix<T> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceEBEMatrix::MultTransAdd"); RegionTimer reg(t);
    Launch (x, y, T(s), true);
  }

  template <typename T>
  AutoVector DeviceEBEMatrix<T> :: CreateRowVector () const
  {
    return make_unique<DeviceVector<T>> (width, memtype);
  }

  template <typename T>
  AutoVector DeviceEBEMatrix<T> :: CreateColVector () const
  {
    return make_unique<DeviceVector<T>> (height, memtype);
  }

  template <typename T>
  BaseMatrix::OperatorInfo DeviceEBEMatrix<T> :: GetOperatorInfo () const
  {
    return { string("DeviceEBEMatrix<") + (is_same_v<T,double> ? "double" : "float")
             + "> (elements=" + ToString(nel) + ")", height, width };
  }

  template <typename T>
  ostream & DeviceEBEMatrix<T> :: Print (ostream & ost) const
  {
    ost << "DeviceEBEMatrix<" << (is_same_v<T,double> ? "double" : "float")
        << ">, " << height << " x " << width << ", elements = " << nel
        << ", on " << device->Name() << endl;
    return ost;
  }


  template <typename SCAL>
  shared_ptr<BaseMatrix> ElementByElementMatrix<SCAL> :: CreateDeviceMatrix () const
  {
    if constexpr (is_same_v<SCAL,double> || is_same_v<SCAL,float>)
      if (ngs_gpu::HasDevice())
        {
          if constexpr (is_same_v<SCAL,double>)
            if (GetGpuDevice()->HasFloat64())
              return make_shared<DeviceEBEMatrix<double>> (*this);
          return make_shared<DeviceEBEMatrix<float>> (*this);
        }
    return BaseMatrix::CreateDeviceMatrix();
  }

  template shared_ptr<BaseMatrix> ElementByElementMatrix<double>::CreateDeviceMatrix () const;
  template shared_ptr<BaseMatrix> ElementByElementMatrix<float>::CreateDeviceMatrix () const;
  template shared_ptr<BaseMatrix> ElementByElementMatrix<Complex>::CreateDeviceMatrix () const;

  template class DeviceEBEMatrix<double>;
  template class DeviceEBEMatrix<float>;
  template DeviceEBEMatrix<double>::DeviceEBEMatrix (const ElementByElementMatrix<double>&);
  template DeviceEBEMatrix<double>::DeviceEBEMatrix (const ElementByElementMatrix<float>&);
  template DeviceEBEMatrix<float>::DeviceEBEMatrix (const ElementByElementMatrix<double>&);
  template DeviceEBEMatrix<float>::DeviceEBEMatrix (const ElementByElementMatrix<float>&);
}
