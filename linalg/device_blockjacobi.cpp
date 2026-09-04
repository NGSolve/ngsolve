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
    string SubstituteScal (string src, const string & value)
    {
      const string token = "SCAL";
      for (size_t p = src.find(token); p != string::npos; p = src.find(token, p))
        src.replace (p, token.size(), value);
      return src;
    }

    /*
      lanes consecutive work-items share a block and take its rows in
      turns. trans applies the transposed inverse. The atomic variant is
      for block tables where a dof occurs more than once.
    */
    const char * kernel_source = R"RAW(

      #define BJ_BODY(ACCUMULATE)                                            \
        int lane = int(LOCAL_ID_X) & (lanes-1);                              \
        int block = int(GLOBAL_ID_X) / lanes;                                \
        if (block >= nblocks) return;                                        \
        int first = blockfirst[block];                                       \
        int bs = blockfirst[block+1] - first;                                \
        GLOBAL_PTR(const int) idx = indices + first;                         \
        GLOBAL_PTR(const SCAL) mat = mats + matfirst[block];                 \
        for (int r = lane; r < bs; r += lanes)                               \
          {                                                                  \
            SCAL sum = 0;                                                    \
            if (trans)                                                       \
              for (int c = 0; c < bs; c++) sum += mat[c*bs+r] * x[idx[c]];   \
            else                                                             \
              for (int c = 0; c < bs; c++) sum += mat[r*bs+c] * x[idx[c]];   \
            ACCUMULATE;                                                      \
          }

      KERNEL(bj_mult, GLOBAL_IN(int,blockfirst), GLOBAL_IN(int,indices), GLOBAL_IN(int,matfirst),
                      GLOBAL_IN(SCAL,mats), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                      VALUE(SCAL,s), VALUE(int,lanes), VALUE(int,trans), VALUE(int,nblocks))
      {
        BJ_BODY(y[idx[r]] += s*sum)
      }

      KERNEL(bj_mult_atomic, GLOBAL_IN(int,blockfirst), GLOBAL_IN(int,indices), GLOBAL_IN(int,matfirst),
                             GLOBAL_IN(SCAL,mats), GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y),
                             VALUE(SCAL,s), VALUE(int,lanes), VALUE(int,trans), VALUE(int,nblocks))
      {
        BJ_BODY(ATOMIC_ADD(&y[idx[r]], s*sum))
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
                                         SubstituteScal (kernel_source, scal));
        mult        = library->GetKernel ("bj_mult");
        mult_atomic = library->GetKernel ("bj_mult_atomic");
        queue = device->DefaultQueue();
        groupsize = (device->SimdWidth() > 1) ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
      }

      // a power of two up to the simd width, matching the average block size
      int ChooseLanes (size_t nindices, size_t nblocks) const
      {
        size_t simd = device->SimdWidth();
        if (simd <= 1 || nblocks == 0) return 1;
        double avg = double(nindices) / nblocks;
        int lanes = 1;
        while (lanes < avg && size_t(2*lanes) <= simd && size_t(2*lanes) <= groupsize)
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

    // inverses row-major, block after block, values converted to T
    Array<T> mats(nmat);
    for (size_t i = 0; i < nblocks; i++)
      {
        size_t bs = blocktable[i].Size();
        T * dst = mats.Data() + matfirst[i];
        if (bs == 0) continue;
        const auto & inv = inverses[i];
        for (size_t r = 0; r < bs; r++)
          for (size_t c = 0; c < bs; c++)
            dst[r*bs+c] = T(inv(r,c));
      }

    // disjoint blocks can be accumulated without atomics
    Array<int> count(width);
    count = 0;
    overlapping = false;
    for (int d : indices)
      if (count[d]++ > 0) { overlapping = true; break; }

    lanes = kern.ChooseLanes (nind, nblocks);

    cout << IM(7) << "DeviceBlockJacobi<" << (is_same_v<T,double> ? "double" : "float")
         << "> nblocks = " << nblocks << ", indices = " << nind << ", matrix entries = " << nmat
         << ", lanes = " << lanes << (overlapping ? ", overlapping" : ", disjoint") << endl;

    dev_blockfirst = device->NewBuffer<int> (nblocks+1, MemType::Device);
    dev_matfirst   = device->NewBuffer<int> (nblocks+1, MemType::Device);
    dev_indices    = device->NewBuffer<int> (max<size_t>(nind,1), MemType::Device);
    dev_mats       = device->NewBuffer<T> (max<size_t>(nmat,1), MemType::Device);

    dev_blockfirst.H2D (blockfirst.Data(), nblocks+1);
    dev_matfirst.H2D (matfirst.Data(), nblocks+1);
    dev_indices.H2D (indices.Data(), nind);
    dev_mats.H2D (mats.Data(), nmat);
  }


  template <typename T>
  void DeviceBlockJacobi<T> :: Launch (const BaseVector & x, BaseVector & y, T s, bool trans) const
  {
    if (x.Size() != width || y.Size() != height)
      throw Exception("DeviceBlockJacobi::MultAdd - size mismatch");
    if (nblocks == 0) return;

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);

    const auto & kern = DeviceBlockJacobiKernels<T>::Get();
    size_t items = nblocks * lanes;
    unsigned groups = (items + kern.groupsize-1) / kern.groupsize;
    queue->Launch (overlapping ? *kern.mult_atomic : *kern.mult,
                   Dim3(groups), Dim3(kern.groupsize),
                   { KernelArg(dev_blockfirst), KernelArg(dev_indices), KernelArg(dev_matfirst),
                     KernelArg(dev_mats), ux.DevArgRO(), uy.DevArgRW(),
                     KernelArg(s), KernelArg(int(lanes)), KernelArg(int(trans)), KernelArg(int(nblocks)) });
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
  AutoVector DeviceBlockJacobi<T> :: CreateRowVector () const
  {
    return make_unique<DeviceVector<T>> (width, memtype);
  }

  template <typename T>
  AutoVector DeviceBlockJacobi<T> :: CreateColVector () const
  {
    return make_unique<DeviceVector<T>> (height, memtype);
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
