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
    return (Dev<double>*)ngs_cuda::BufferDevPtr(*devbuffer.Raw()) + devoffset;
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




}
