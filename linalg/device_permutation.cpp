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

    const char * kernel_source = R"RAW(

      KERNEL(perm_mult, GLOBAL_IN(int,ind), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                        VALUE(SCAL,s), VALUE(SCAL,beta), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i >= n) return;
        SCAL val = (ind[i] >= 0) ? s*x[ind[i]] : SCAL(0);
        y[i] = (beta == SCAL(0)) ? val : beta*y[i] + val;
      }

      KERNEL(perm_multtrans, GLOBAL_IN(int,ind), GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y),
                             VALUE(SCAL,s), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n && ind[i] >= 0) ATOMIC_ADD (&y[ind[i]], s * x[i]);
      }

    )RAW";

    template <typename T>
    class PermKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<ngs_gpu::Queue> queue;
      shared_ptr<Kernel> mult, multtrans;
      unsigned groupsize;

      PermKernels (shared_ptr<Device> adevice) : device(adevice)
      {
        if constexpr (is_same_v<T,double>)
          if (!device->HasFloat64())
            throw Exception("DevicePermutationMatrix: double vectors on "+device->Name()+
                            ", which has no fp64");

        string scal = is_same_v<T,double> ? "double" : "float";
        library = device->CompileSource (string(code_gpukernel) +
                                         SubstituteScal (kernel_source, scal));
        mult = library->GetKernel ("perm_mult");
        multtrans = library->GetKernel ("perm_multtrans");
        queue = device->DefaultQueue();
        groupsize = (device->SimdWidth() > 1) ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
      }

      static const PermKernels & Get()
      {
        static mutex mtx;
        static shared_ptr<PermKernels> cached;
        auto dev = GetGpuDevice();
        auto lock = lock_guard<mutex>(mtx);
        if (!cached || cached->device != dev)
          cached = make_shared<PermKernels> (dev);
        return *cached;
      }
    };

    template <typename T>
    void LaunchPerm (const TypedBuffer<int> & dev_ind, size_t height, MemType memtype,
                     bool trans, T s, T beta, const BaseVector & x, BaseVector & y)
    {
      DeviceVectorWrapper<T> ux(x, memtype);
      DeviceVectorWrapper<T> uy(y, memtype);

      const auto & kern = PermKernels<T>::Get();
      unsigned groups = (height + kern.groupsize-1) / kern.groupsize;
      if (trans)
        kern.queue->Launch (*kern.multtrans, Dim3(groups), Dim3(kern.groupsize),
                            { KernelArg(dev_ind), ux.DevArgRO(), uy.DevArgRW(),
                              KernelArg(s), KernelArg(int(height)) });
      else
        kern.queue->Launch (*kern.mult, Dim3(groups), Dim3(kern.groupsize),
                            { KernelArg(dev_ind), ux.DevArgRO(),
                              beta == T(0) ? uy.DevArgW() : uy.DevArgRW(),
                              KernelArg(s), KernelArg(beta), KernelArg(int(height)) });
    }

    bool IsFloatVector (const BaseVector & v)
    { return dynamic_cast<const DeviceVector<float>*> (&v) != nullptr; }
  }


  class DevicePermutationMatrix : public BaseMatrix
  {
    size_t height, width;
    MemType memtype;
    shared_ptr<ngs_gpu::Device> device;
    shared_ptr<ngs_gpu::Queue> queue;
    ngs_gpu::TypedBuffer<int> dev_ind;

    void Launch (bool trans, double s, double beta,
                 const BaseVector & x, BaseVector & y) const;
    AutoVector CreateVec (size_t size) const;

  public:
    DevicePermutationMatrix (const PermutationMatrix & mat);

    int VHeight() const override { return height; }
    int VWidth() const override { return width; }
    bool IsComplex() const override { return false; }

    void Mult (const BaseVector & x, BaseVector & y) const override;
    void MultTrans (const BaseVector & x, BaseVector & y) const override;
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    AutoVector CreateRowVector () const override { return CreateVec (width); }
    AutoVector CreateColVector () const override { return CreateVec (height); }

    VecFormat RowFormat () const override { return VecFormat(width).OnDevice(memtype); }
    VecFormat ColFormat () const override { return VecFormat(height).OnDevice(memtype); }

    BaseMatrix::OperatorInfo GetOperatorInfo () const override
    { return { "DevicePermutationMatrix", height, width }; }
  };


  DevicePermutationMatrix :: DevicePermutationMatrix (const PermutationMatrix & mat)
    : height(mat.Height()), width(mat.Width()), memtype(PreferredMemType())
  {
    if (height > size_t(INT_MAX) || width > size_t(INT_MAX))
      throw Exception("DevicePermutationMatrix: matrix too large for 32-bit indices");

    device = GetGpuDevice();
    queue = device->DefaultQueue();

    static Timer tup("DevicePermutationMatrix ctor upload");
    RegionTimer reg(tup);

    FlatArray<size_t> ind = mat.GetIndices();
    dev_ind = device->NewBuffer<int> (max<size_t>(height,1), MemType::Device);
    dev_ind.Fill (height, [ind] (int * dst, size_t off, size_t n)
      {
        ParallelFor (n, [dst, off, ind] (size_t i)
          { dst[i] = (ind[off+i] == size_t(-1)) ? -1 : int(ind[off+i]); });
      });
  }

  void DevicePermutationMatrix :: Launch (bool trans, double s, double beta,
                                          const BaseVector & x, BaseVector & y) const
  {
    if (height == 0) return;
    if ((trans ? y.Size() : x.Size()) != width ||
        (trans ? x.Size() : y.Size()) != height)
      throw Exception("DevicePermutationMatrix - size mismatch");

    if (!device->HasFloat64() || (IsFloatVector(x) && IsFloatVector(y)))
      LaunchPerm<float> (dev_ind, height, memtype, trans, float(s), float(beta), x, y);
    else
      LaunchPerm<double> (dev_ind, height, memtype, trans, s, beta, x, y);
  }

  AutoVector DevicePermutationMatrix :: CreateVec (size_t size) const
  {
    if (device->HasFloat64())
      return make_unique<DeviceVector<double>> (size, memtype);
    return make_unique<DeviceVector<float>> (size, memtype);
  }

  void DevicePermutationMatrix :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DevicePermutationMatrix::Mult"); RegionTimer reg(t);
    Launch (false, 1.0, 0.0, x, y);        // rows without a source stay zero
  }

  void DevicePermutationMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DevicePermutationMatrix::MultAdd"); RegionTimer reg(t);
    Launch (false, s, 1.0, x, y);
  }

  void DevicePermutationMatrix :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DevicePermutationMatrix::MultTrans"); RegionTimer reg(t);
    y = 0.0;                               // the scatter only touches ind[i]
    Launch (true, 1.0, 1.0, x, y);
  }

  void DevicePermutationMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DevicePermutationMatrix::MultTransAdd"); RegionTimer reg(t);
    Launch (true, s, 1.0, x, y);
  }


  shared_ptr<BaseMatrix> PermutationMatrix :: CreateDeviceMatrix () const
  {
    if (ngs_gpu::HasDevice())
      return make_shared<DevicePermutationMatrix> (*this);
    return BaseMatrix::CreateDeviceMatrix();
  }
}
