#include <la.hpp>
#include "cuda_linalg.hpp"
#include "cuda_device.hpp"

namespace ngla
{

  UnifiedVector :: UnifiedVector (int asize)
    : DeviceVector<double> (asize, MemType::Device)
  { }

  UnifiedVector :: UnifiedVector (const BaseVector& vec)
    : DeviceVector<double> (vec, MemType::Device)
  { }

  UnifiedVector :: UnifiedVector (const UnifiedVector& vec)
    : DeviceVector<double> (vec, MemType::Device)
  { }

  UnifiedVector :: UnifiedVector (UnifiedVector && vec)
    : DeviceVector<double> (vec, MemType::Device)
  { }


  Dev<double> * UnifiedVector :: DevData() const
  {
    return (Dev<double>*)ngs_cuda::BufferDevPtr(*devbuffer) + devoffset;
  }

  FlatVector<Dev<double>> UnifiedVector :: FVDev() const
  {
    UpdateDevice();
    InvalidateHost();
    return { Size(), DevData() };
  }

  FlatVector<Dev<double>> UnifiedVector :: FVDevRO() const
  {
    UpdateDevice();
    return { Size(), DevData() };
  }

  FlatVector<Complex> UnifiedVector :: FVComplex () const
  {
    throw Exception ("unified complex not yet supported");
  }


  BaseVector & UnifiedVector :: operator= (double d)
  {
    return SetScalar (d);
  }

  BaseVector & UnifiedVector :: operator= (const BaseVector & v2)
  {
    if (auto dv2 = dynamic_cast<const DeviceVector<double>*> (&v2))
      if (dv2->IsDevUptodate())
        {
          auto src = (double*)ngs_cuda::BufferDevPtr(*dv2->DevBufferRO()) + dv2->DevOffset();
          cudaMemcpy (DevData(), src, sizeof(double)*size, cudaMemcpyDeviceToDevice);
          dev_uptodate = true;
          InvalidateHost();
          return *this;
        }

    FVDouble() = v2.FVDouble();
    UpdateDevice();
    return *this;
  }

  UnifiedVector & UnifiedVector :: operator= (const UnifiedVector & v2)
  {
    operator= (static_cast<const BaseVector&> (v2));
    return *this;
  }


  const double & UnifiedVector :: operator [] (const int ind) const
  {
    UpdateHost();
    return host_data[ind];
  }

  double & UnifiedVector :: operator [] (const int ind)
  {
    UpdateHost();
    InvalidateDevice();
    return host_data[ind];
  }

  AutoVector UnifiedVector :: Range (T_Range<size_t> range) const
  {
    return make_unique<UnifiedVectorWrapper>(*this, range);
  }


  double UnifiedVector :: InnerProductD (const BaseVector & v2) const
  {
    static Timer tdot("CUDA InnerProduct");
    RegionTimer reg(tdot);

    if (auto uv2 = dynamic_cast<const UnifiedVector*> (&v2))
      {
        static Timer tdot("CUDA InnerProduct devdev");
        RegionTimer reg(tdot);
        UpdateDevice();
        uv2->UpdateDevice();

        double res;
        ngla::SetCuBlasStream(ngs_cuda_stream);
        cublasDdot (Get_CuBlas_Handle(),
                    size, (double*)DevData(), 1, (double*)uv2->DevData(), 1, &res);
        return res;
      }

    FlatVector<> fv = FVDouble();
    FlatVector<> fv2 = v2.FVDouble();
    return ngbla::InnerProduct (fv, fv2);
  }

  double UnifiedVector :: L2Norm() const
  {
    UpdateDevice();
    double res;
    ngla::SetCuBlasStream(ngs_cuda_stream);
    cublasDnrm2(Get_CuBlas_Handle(), size, (double*)DevData(), 1, &res);
    return res;
  }


  ostream & UnifiedVector :: Print (ostream & ost) const
  {
    ost << "output unified vector of size " << size;
    ost << ", host = " << IsHostUptodate() << ", dev = " << IsDevUptodate() << endl;
    if (!IsHostUptodate())
      {
        if (IsDevUptodate())
          {
            ost << "host not up-to-date. printing device data" << endl;
            Vector<double> tmp(size);
            cudaMemcpy(tmp.Data(), DevData(), size * sizeof(double), cudaMemcpyDeviceToHost);
            ost << tmp << endl;
          }
        else
          ost << "undefined vector" << endl;
      }
    else
      {
        ost << FlatVector<double>(size, (double*)host_data);
      }
    return ost;
  }

  ostream & UnifiedVector :: PrintStatus (ostream & ost) const
  {
    ost << "output unified vector of size " << size;
    ost << ", host = " << IsHostUptodate() << ", dev = " << IsDevUptodate() << endl;
    return ost;
  }


  AutoVector UnifiedVector :: CreateVector () const
  {
    return make_unique<UnifiedVector> (size);
  }


  UnifiedVectorWrapper :: UnifiedVectorWrapper(const BaseVector & vec_, optional<IntRange> opt_range)
    : vec(vec_)
  {
    IntRange range = opt_range.value_or (IntRange(0, vec.Size()));
    this->size = range.Size();

    if (auto p = dynamic_cast<const DeviceVector<double>*>(&vec_))
      {
        // alias the storage (possibly a sub-range), hand the flags back later
        alias_of = p;
        subrange = (range.Size() != vec.Size());
        // a device write to a sub-range must not leave the rest of the
        // vector stale, so bring the whole vector to the device first
        if (subrange)
          p->UpdateDevice();
        this->devbuffer = p->devbuffer;
        this->queue = p->queue;
        this->devoffset = p->devoffset + range.First();
        this->memtype = p->memtype;
        this->host_data = p->host_data ? p->host_data + range.First() : nullptr;
        this->host_uptodate = p->host_uptodate;
        this->dev_uptodate = p->dev_uptodate;
      }
    else
      {
        // a device buffer of our own, the host side stays with vec
        this->memtype = MemType::Device;
        this->AllocBuffer (this->size);
        this->hostmem.SetSize(0);
        this->host_data = vec.FVDouble().Data() + range.First();
        this->host_uptodate = true;
        this->dev_uptodate = false;
      }
    initial_host_uptodate = this->host_uptodate;
    initial_dev_uptodate = this->dev_uptodate;
  }

  UnifiedVectorWrapper :: ~UnifiedVectorWrapper()
  {
    if (alias_of)
      {
        if (subrange)
          {
            // an update of the sub-range cannot validate the full
            // vector, but staleness propagates
            alias_of->host_uptodate &= this->host_uptodate;
            alias_of->dev_uptodate &= this->dev_uptodate;
          }
        else
          {
            // full alias: same storage, the flags simply carry over
            alias_of->host_uptodate = this->host_uptodate;
            alias_of->dev_uptodate = this->dev_uptodate;
          }
        this->host_data = nullptr;
        return;
      }

    // only pay for the copy back if a kernel actually wrote
    if (initial_host_uptodate && !this->host_uptodate)
      this->UpdateHost();
    this->host_data = nullptr;
  }

  // -------------------------------------------------------
  // UnifiedScalar: only the typed cuda pointer is left here
  // -------------------------------------------------------
  double* UnifiedScalar :: DevPtr() const
  {
    return (double*)ngs_cuda::BufferDevPtr(*DevBuffer());
  }

  // -------------------------------------------------------
  // UnifiedVector overrides for BaseScalar
  // -------------------------------------------------------

  shared_ptr<BaseScalar> UnifiedVector :: CreateScalar() const
  {
    return make_shared<UnifiedScalar>();
  }

  void UnifiedVector :: InnerProduct(const BaseVector& v2,
                                      BaseScalar& scal,
                                      bool conjugate) const
  {
    if (auto uscal = dynamic_cast<UnifiedScalar*>(&scal))
      {
        // check if we are inside graph capture:
        cudaStreamCaptureStatus capture_status;
        cudaStreamGetCaptureInfo(ngs_cuda_stream, &capture_status, nullptr);
        bool capturing = (capture_status == cudaStreamCaptureStatusActive);

        // only call UpdateDevice when NOT capturing
        // (during capture/replay the data is already on GPU)
        if (!capturing) {
            UpdateDevice();
        }

        auto uv2 = dynamic_cast<const UnifiedVector*>(&v2);
        if (!uv2) throw Exception("InnerProduct(BaseScalar): v2 must be UnifiedVector");
        if (!capturing) {
            uv2->UpdateDevice();
        }

        SetCuBlasStream(ngs_cuda_stream);
        cublasSetPointerMode(Get_CuBlas_Handle(), CUBLAS_POINTER_MODE_DEVICE);
        cublasDdot(Get_CuBlas_Handle(),
                   size,
                   (double*)DevData(), 1,
                   (double*)uv2->DevData(), 1,
                   uscal->DevPtr());
        // only restore pointer mode when NOT capturing
        // during capture/replay, caller manages pointer mode
        if (!capturing) {
            cublasSetPointerMode(Get_CuBlas_Handle(), CUBLAS_POINTER_MODE_HOST);
        }
      }
    else
      {
        // fallback — CPU path
        BaseVector::InnerProduct(v2, scal, conjugate);
      }
  }



}
