/*********************************************************************/
/* File:   device_diagonalmatrix.cpp                                 */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   4. Sep. 2026                                              */
/*********************************************************************/

#define FILE_DEVICE_DIAGONALMATRIX_CPP

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

      // y[i] = beta*y[i] + s * d[i] * x[i]; beta 0 overwrites without reading y
      KERNEL(diag_mult, GLOBAL_IN(SCAL,d), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                        VALUE(SCAL,s), VALUE(SCAL,beta), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n) y[i] = (beta == SCAL(0)) ? s * d[i] * x[i] : beta*y[i] + s * d[i] * x[i];
      }


      // the projector: a 0/1 diagonal kept as a bit mask, packed into 32 bit words.
      #define NGS_PROJ_SET(i)  ((bits[(i) >> 5] >> ((i) & 31)) & 1)

      KERNEL(proj_mult, GLOBAL_IN(int,bits), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                        VALUE(int,keep), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n) y[i] = (NGS_PROJ_SET(i) == keep) ? x[i] : SCAL(0);
      }

      KERNEL(proj_multadd, GLOBAL_IN(int,bits), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                           VALUE(SCAL,s), VALUE(int,keep), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n && NGS_PROJ_SET(i) == keep) y[i] += s * x[i];
      }

      // in place: drop everything the projector does not keep
      KERNEL(proj_project, GLOBAL_IN(int,bits), GLOBAL(SCAL,x),
                           VALUE(int,keep), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n && NGS_PROJ_SET(i) != keep) x[i] = SCAL(0);
      }


      // block-diagonal SoA: y(j,i) = beta*y(j,i) + s * sum_k a(aind[k],i) * x(xind[k],i)
      // for the entries k of row j, one work-item per block index i
      KERNEL(blockdiag_soa, GLOBAL_IN(SCAL,a), GLOBAL_IN(int,first), GLOBAL_IN(int,aind), GLOBAL_IN(int,xind),
                            GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                            VALUE(SCAL,s), VALUE(SCAL,beta), VALUE(int,nrows), VALUE(int,blocks))
      {
        int i = int(GLOBAL_ID_X);
        if (i >= blocks) return;
        for (int j = 0; j < nrows; j++)
          {
            SCAL sum = SCAL(0);
            for (int k = first[j]; k < first[j+1]; k++)
              sum += a[aind[k]*blocks+i] * x[xind[k]*blocks+i];
            int yi = j*blocks+i;
            y[yi] = (beta == SCAL(0)) ? s*sum : beta*y[yi] + s*sum;
          }
      }

    )RAW";


    template <typename T>
    class DeviceDiagonalKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<ngs_gpu::Queue> queue;
      shared_ptr<Kernel> mult;
      shared_ptr<Kernel> proj_mult, proj_multadd, proj_project;
      shared_ptr<Kernel> blockdiag_soa;
      unsigned groupsize;

      DeviceDiagonalKernels (shared_ptr<Device> adevice)
        : device(adevice)
      {
        string scal = is_same_v<T,double> ? "double" : "float";
        library = device->CompileSource (string(code_gpukernel) +
                                         SubstituteScal (kernel_source, scal));
        mult = library->GetKernel ("diag_mult");
        proj_mult    = library->GetKernel ("proj_mult");
        proj_multadd = library->GetKernel ("proj_multadd");
        proj_project = library->GetKernel ("proj_project");
        blockdiag_soa = library->GetKernel ("blockdiag_soa");
        queue = device->DefaultQueue();
        groupsize = (device->SimdWidth() > 1) ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
      }

      static const DeviceDiagonalKernels & Get()
      {
        static mutex mtx;
        static shared_ptr<DeviceDiagonalKernels> cached;

        auto dev = GetGpuDevice();
        auto lock = lock_guard<mutex>(mtx);
        if (!cached || cached->device != dev)
          cached = make_shared<DeviceDiagonalKernels> (dev);
        return *cached;
      }
    };
  }





  class DeviceProjector : public Projector
  {
    MemType memtype;
    TypedBuffer<int> dev_bits;

    void Launch (Kernel & kernel, size_t n,
                 const std::vector<KernelArg> & args) const
    {
      const auto & kern = DeviceDiagonalKernels<double>::Get();
      unsigned groups = (n + kern.groupsize-1) / kern.groupsize;
      kern.queue->Launch (kernel, Dim3(groups), Dim3(kern.groupsize), args);
    }

  public:
    DeviceProjector (const Projector & proj)
      : Projector (proj.Mask(), proj.KeepValues()), memtype (PreferredMemType())
    {
      auto & ba = *Mask();
      size_t nwords = (ba.Size() + 31) / 32;
      std::vector<int> words (max<size_t>(nwords,1), 0);
      for (size_t i = 0; i < ba.Size(); i++)
        if (ba.Test(i)) words[i >> 5] |= (1 << (i & 31));
      dev_bits = DeviceDiagonalKernels<double>::Get().device->NewBuffer<int> (words.size(), MemType::Device);
      dev_bits.H2D (words.data(), words.size());
    }

    void Mult (const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("DeviceProjector::Mult"); RegionTimer reg(t);
      size_t n = Mask()->Size();
      if (n == 0) return;
      DeviceVectorWrapper<double> ux(x, memtype);
      DeviceVectorWrapper<double> uy(y, memtype);
      Launch (*DeviceDiagonalKernels<double>::Get().proj_mult, n,
              { dev_bits, ux.DevArgRO(), uy.DevArgW(),
                KernelArg(int(KeepValues())), KernelArg(int(n)) });
    }

    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("DeviceProjector::MultAdd"); RegionTimer reg(t);
      size_t n = Mask()->Size();
      if (n == 0) return;
      DeviceVectorWrapper<double> ux(x, memtype);
      DeviceVectorWrapper<double> uy(y, memtype);
      Launch (*DeviceDiagonalKernels<double>::Get().proj_multadd, n,
              { dev_bits, ux.DevArgRO(), uy.DevArgRW(),
                KernelArg(s), KernelArg(int(KeepValues())), KernelArg(int(n)) });
    }

    void Project (BaseVector & x) const override
    {
      static Timer t("DeviceProjector::Project"); RegionTimer reg(t);
      size_t n = Mask()->Size();
      if (n == 0) return;
      DeviceVectorWrapper<double> ux(x, memtype);
      Launch (*DeviceDiagonalKernels<double>::Get().proj_project, n,
              { dev_bits, ux.DevArgRW(),
                KernelArg(int(KeepValues())), KernelArg(int(n)) });
    }

    xbool IsSymmetric () const override { return true; }
    void MultTrans (const BaseVector & x, BaseVector & y) const override { Mult (x, y); }
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override { MultAdd (s, x, y); }

    AutoVector CreateRowVector () const override
    { return make_unique<DeviceVector<double>> (Mask()->Size(), memtype); }
    AutoVector CreateColVector () const override
    { return make_unique<DeviceVector<double>> (Mask()->Size(), memtype); }
    VecFormat RowFormat () const override { return DeviceVectorFormat<double> (Mask()->Size(), memtype); }
    VecFormat ColFormat () const override { return RowFormat(); }

    BaseMatrix::OperatorInfo GetOperatorInfo () const override
    { return { "DeviceProjector", Mask()->Size(), Mask()->Size() }; }
  };


  shared_ptr<BaseMatrix> Projector :: CreateDeviceMatrix () const
  {
    if (ngs_gpu::HasDevice())
      return make_shared<DeviceProjector> (*this);
    return BaseMatrix::CreateDeviceMatrix();
  }


  template <typename T>
  template <typename TS>
  DeviceDiagonalMatrix<T> :: DeviceDiagonalMatrix (FlatVector<TS> adiag)
    : diag (adiag.Size(), PreferredMemType())
  {
    // fills the host side, uploaded on first device use
    diag.FVScal() = adiag;
  }


  template <typename T>
  void DeviceDiagonalMatrix<T> :: Launch (const BaseVector & x, BaseVector & y, T s, T beta) const
  {
    size_t n = diag.Size();
    if (x.Size() != n || y.Size() != n)
      throw Exception("DeviceDiagonalMatrix::Mult - size mismatch");
    if (n == 0) return;

    DeviceVectorWrapper<T> ux(x, diag.GetMemType());
    DeviceVectorWrapper<T> uy(y, diag.GetMemType());

    const auto & kern = DeviceDiagonalKernels<T>::Get();
    unsigned groups = (n + kern.groupsize-1) / kern.groupsize;
    kern.queue->Launch (*kern.mult, Dim3(groups), Dim3(kern.groupsize),
                        { diag.DevArgRO(), ux.DevArgRO(), beta == T(0) ? uy.DevArgW() : uy.DevArgRW(),
                          KernelArg(s), KernelArg(beta), KernelArg(int(n)) });
  }

  template <typename T>
  void DeviceDiagonalMatrix<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceDiagonalMatrix::MultAdd"); RegionTimer reg(t);
    Launch (x, y, T(s), T(1));
  }

  template <typename T>
  void DeviceDiagonalMatrix<T> :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceDiagonalMatrix::Mult"); RegionTimer reg(t);
    Launch (x, y, T(1), T(0));
  }


  template <typename T>
  AutoVector DeviceDiagonalMatrix<T> :: CreateRowVector () const
  {
    return make_unique<DeviceVector<T>> (diag.Size(), diag.GetMemType());
  }

  template <typename T>
  AutoVector DeviceDiagonalMatrix<T> :: CreateColVector () const
  {
    return make_unique<DeviceVector<T>> (diag.Size(), diag.GetMemType());
  }

  template <typename T>
  VecFormat DeviceDiagonalMatrix<T> :: RowFormat () const
  { return DeviceVectorFormat<T> (diag.Size(), diag.GetMemType()); }
  template <typename T>
  VecFormat DeviceDiagonalMatrix<T> :: ColFormat () const
  { return DeviceVectorFormat<T> (diag.Size(), diag.GetMemType()); }

  template <typename T>
  BaseMatrix::OperatorInfo DeviceDiagonalMatrix<T> :: GetOperatorInfo () const
  {
    return { string("DeviceDiagonalMatrix<") + (is_same_v<T,double> ? "double" : "float") + ">",
             diag.Size(), diag.Size() };
  }

  template <typename T>
  ostream & DeviceDiagonalMatrix<T> :: Print (ostream & ost) const
  {
    ost << "DeviceDiagonalMatrix<" << (is_same_v<T,double> ? "double" : "float")
        << ">, size = " << diag.Size() << endl;
    return ost;
  }


  template <typename T>
  DeviceBlockDiagonalMatrixSoA<T> :: DeviceBlockDiagonalMatrixSoA (const BlockDiagonalMatrixSoA & mat)
    : memtype (PreferredMemType())
  {
    FlatTensor<3> blockdiag = mat.GetBlockDiag();
    dimy = blockdiag.GetSize();
    dimx = blockdiag.GetSubTensor().GetSize();
    blocks = blockdiag.GetSubTensor().GetSubTensor().GetSize();

    auto device = DeviceDiagonalKernels<T>::Get().device;
    size_t n = size_t(dimy)*dimx*blocks;
    dev_data = device->template NewBuffer<T> (max<size_t>(n,1), MemType::Device);
    if constexpr (is_same_v<T,double>)
      dev_data.H2D (blockdiag.Data(), n);
    else
      {
        std::vector<T> tmp (n);
        for (size_t i = 0; i < n; i++) tmp[i] = T(blockdiag.Data()[i]);
        dev_data.H2D (tmp.data(), n);
      }

    // row lists: a-row i*dimx+j holds block (i,j)
    auto upload = [&] (FlatTable<int> table, bool trans,
                       ngs_gpu::TypedBuffer<int> & bfirst, ngs_gpu::TypedBuffer<int> & baind,
                       ngs_gpu::TypedBuffer<int> & bxind)
    {
      std::vector<int> hfirst, haind, hxind;
      hfirst.push_back(0);
      for (size_t j = 0; j < table.Size(); j++)
        {
          for (int k : table[j])
            {
              haind.push_back (trans ? k*dimx+int(j) : int(j)*dimx+k);
              hxind.push_back (k);
            }
          hfirst.push_back (int(haind.size()));
        }
      bfirst = device->template NewBuffer<int> (hfirst.size(), MemType::Device);
      baind  = device->template NewBuffer<int> (max<size_t>(haind.size(),1), MemType::Device);
      bxind  = device->template NewBuffer<int> (max<size_t>(hxind.size(),1), MemType::Device);
      bfirst.H2D (hfirst.data(), hfirst.size());
      baind.H2D (haind.data(), haind.size());
      bxind.H2D (hxind.data(), hxind.size());
    };
    upload (mat.GetSparseMatrix(), false, first, aind, xind);
    upload (mat.GetSparseMatrixTrans(), true, firstT, aindT, xindT);
  }

  template <typename T>
  void DeviceBlockDiagonalMatrixSoA<T> :: Launch (const BaseVector & x, BaseVector & y, T s, T beta, bool trans) const
  {
    int nrows = trans ? dimx : dimy, ncols = trans ? dimy : dimx;
    if (x.Size() != size_t(ncols)*blocks || y.Size() != size_t(nrows)*blocks)
      throw Exception("DeviceBlockDiagonalMatrixSoA::Mult - size mismatch");
    if (blocks == 0) return;

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);

    const auto & kern = DeviceDiagonalKernels<T>::Get();
    unsigned groups = (blocks + kern.groupsize-1) / kern.groupsize;
    kern.queue->Launch (*kern.blockdiag_soa, Dim3(groups), Dim3(kern.groupsize),
                        { dev_data,
                          trans ? firstT : first, trans ? aindT : aind, trans ? xindT : xind,
                          ux.DevArgRO(), beta == T(0) ? uy.DevArgW() : uy.DevArgRW(),
                          KernelArg(s), KernelArg(beta), KernelArg(nrows), KernelArg(blocks) });
  }

  template <typename T>
  void DeviceBlockDiagonalMatrixSoA<T> :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceBlockDiagonalMatrixSoA::Mult"); RegionTimer reg(t);
    Launch (x, y, T(1), T(0), false);
  }

  template <typename T>
  void DeviceBlockDiagonalMatrixSoA<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceBlockDiagonalMatrixSoA::MultAdd"); RegionTimer reg(t);
    Launch (x, y, T(s), T(1), false);
  }

  template <typename T>
  void DeviceBlockDiagonalMatrixSoA<T> :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceBlockDiagonalMatrixSoA::MultTrans"); RegionTimer reg(t);
    Launch (x, y, T(1), T(0), true);
  }

  template <typename T>
  void DeviceBlockDiagonalMatrixSoA<T> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceBlockDiagonalMatrixSoA::MultTransAdd"); RegionTimer reg(t);
    Launch (x, y, T(s), T(1), true);
  }

  template <typename T>
  AutoVector DeviceBlockDiagonalMatrixSoA<T> :: CreateRowVector () const
  {
    return make_unique<DeviceVector<T>> (size_t(dimx)*blocks, memtype);
  }

  template <typename T>
  AutoVector DeviceBlockDiagonalMatrixSoA<T> :: CreateColVector () const
  {
    return make_unique<DeviceVector<T>> (size_t(dimy)*blocks, memtype);
  }

  template <typename T>
  VecFormat DeviceBlockDiagonalMatrixSoA<T> :: RowFormat () const
  { return DeviceVectorFormat<T> (size_t(dimx)*blocks, memtype); }
  template <typename T>
  VecFormat DeviceBlockDiagonalMatrixSoA<T> :: ColFormat () const
  { return DeviceVectorFormat<T> (size_t(dimy)*blocks, memtype); }

  template <typename T>
  BaseMatrix::OperatorInfo DeviceBlockDiagonalMatrixSoA<T> :: GetOperatorInfo () const
  {
    return { string("DeviceBlockDiagonalMatrixSoA<") + (is_same_v<T,double> ? "double" : "float") + ">",
             size_t(dimy)*blocks, size_t(dimx)*blocks };
  }

  template <typename T>
  ostream & DeviceBlockDiagonalMatrixSoA<T> :: Print (ostream & ost) const
  {
    ost << "DeviceBlockDiagonalMatrixSoA<" << (is_same_v<T,double> ? "double" : "float")
        << ">, blocks = " << blocks << ", dim = " << dimy << " x " << dimx << endl;
    return ost;
  }


  shared_ptr<BaseMatrix> BlockDiagonalMatrixSoA :: CreateDeviceMatrix () const
  {
    if (ngs_gpu::HasDevice())
      {
        if (GetGpuDevice()->HasFloat64())
          return make_shared<DeviceBlockDiagonalMatrixSoA<double>> (*this);
        return make_shared<DeviceBlockDiagonalMatrixSoA<float>> (*this);
      }
    return BaseMatrix::CreateDeviceMatrix();
  }


  template class DeviceDiagonalMatrix<double>;
  template class DeviceDiagonalMatrix<float>;
  template class DeviceBlockDiagonalMatrixSoA<double>;
  template class DeviceBlockDiagonalMatrixSoA<float>;
  template DeviceDiagonalMatrix<double>::DeviceDiagonalMatrix (FlatVector<double>);
  template DeviceDiagonalMatrix<double>::DeviceDiagonalMatrix (FlatVector<float>);
  template DeviceDiagonalMatrix<float>::DeviceDiagonalMatrix (FlatVector<double>);
  template DeviceDiagonalMatrix<float>::DeviceDiagonalMatrix (FlatVector<float>);
}
