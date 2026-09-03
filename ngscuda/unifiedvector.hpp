#ifndef UNIFIED_VECTOR_HPP
#define UNIFIED_VECTOR_HPP

/*
  UnifiedVector: cuda-specific extension of the backend-independent
  DeviceVector<double> (linalg/devicevector.hpp).

  Storage, host/device state tracking and the transfers live in the
  base class. This class adds what needs raw CUDA: typed device
  pointers for cuBLAS/cuSPARSE and hand-written kernels, on-device
  reductions, and graph-capturable operations with device-resident
  scalars. It shrinks as those capabilities move into DeviceVector.
*/

namespace ngla
{

  class UnifiedVector : public DeviceVector<double>
  {
  protected:
    UnifiedVector () = default;

  public:
    UnifiedVector (int asize);
    UnifiedVector (const BaseVector & vec);
    UnifiedVector (const UnifiedVector & vec);
    UnifiedVector (UnifiedVector && vec);

    BaseVector & operator= (double d);
    BaseVector & operator= (const BaseVector & v2);
    UnifiedVector & operator= (const UnifiedVector & v2);

    template <typename T2>
    UnifiedVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }

    const double & operator [] (const int ind) const;
    double & operator [] (const int ind);

    virtual shared_ptr<BaseScalar> CreateScalar() const override;

    virtual ostream & Print (ostream & ost) const override;
    virtual ostream & PrintStatus (ostream & ost) const;
    virtual AutoVector CreateVector () const override;

    virtual FlatVector<Complex> FVComplex () const override;

    // typed device pointer into the buffer, no transfer implied
    virtual Dev<double> * DevData() const;

    // update device and invalidate host
    FlatVector<Dev<double>> FVDev() const;

    // update device
    FlatVector<Dev<double>> FVDevRO() const;
  };





  // -------------------------------------------------------
  // UnifiedScalar: cuda view of the backend-independent
  // DeviceScalar (linalg/devicevector.hpp). Adds the typed
  // device pointer and the device expression assignment
  // used for CUDA graph capture.
  // -------------------------------------------------------
  class UnifiedScalar : public DeviceScalar<double>
  {
  public:
    UnifiedScalar (double d = 0.0) : DeviceScalar<double>(d) { }

    double* DevPtr() const;

    using DeviceScalar<double>::operator=;
    UnifiedScalar & operator= (const UnifiedScalar & s2)
    { DeviceScalar<double>::operator= (static_cast<const DeviceScalar<double>&>(s2)); return *this; }
  };

}
#endif
