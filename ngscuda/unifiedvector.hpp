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

    virtual AutoVector Range (T_Range<size_t> range) const override;

    using BaseVector::InnerProduct;
    virtual double InnerProductD (const BaseVector & v2) const override;
    virtual double L2Norm () const override;

    // BaseScalar override for the cublas dot with device-resident result:
    virtual void InnerProduct (const BaseVector & v2, BaseScalar & scal, bool conjugate = false) const override;
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


  class UnifiedVectorWrapper : public UnifiedVector
  {
    const BaseVector & vec;
    const DeviceVector<double> * alias_of = nullptr;
    bool subrange = false;
    bool initial_host_uptodate = false;
    bool initial_dev_uptodate = false;

  public:
    UnifiedVectorWrapper(const BaseVector & vec_, optional<IntRange> opt_range = nullopt);
    ~UnifiedVectorWrapper();

    using UnifiedVector::operator=;
  };



  // -------------------------------------------------------
  // UnifiedScalar: cuda view of the backend-independent
  // DeviceScalar (linalg/devicevector.hpp). Adds the typed
  // device pointer and the device expression assignment
  // used for CUDA graph capture.
  // -------------------------------------------------------
  class UnifiedScalar : public DeviceScalar
  {
  public:
    UnifiedScalar (double d = 0.0) : DeviceScalar(d) { }

    double* DevPtr() const;

    template<typename Expr>
    UnifiedScalar& operator=(const Expr& expr)
    {
      double* d = DevPtr();
      auto kernel = [d, expr] __device__ (size_t) {
          *d = expr.Evaluate();
      };
      DeviceParallelFor(1, kernel);
      return *this;
    }
  };

}
#endif
