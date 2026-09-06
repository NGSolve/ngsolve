#ifndef CUDA_LINALG_HPP
#define CUDA_LINALG_HPP

// partial override of overloaded function (MultAdd)
#pragma nv_diag_suppress 611
#pragma nv_diag_suppress 20013

#include <la.hpp>

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse.h>

#include "cuda_ngstd.hpp"

namespace ngla
{
  cublasHandle_t Get_CuBlas_Handle ();
  cusparseHandle_t Get_CuSparse_Handle ();
}

  
#include "cuda_ngbla.hpp"
#include "linalg_kernels.hpp"
#include "cuda_device.hpp"
#include "unifiedvector.hpp"


namespace ngla
{
  using namespace ngs_cuda;

  void InitCuLinalg();


  // typed device access to any DeviceVector<double> (UnifiedVector,
  // DeviceVectorWrapper, ...), with the transfers the access implies

  inline Dev<double> * DevPtr (const DeviceVector<double> & v)
  {
    return (Dev<double>*)ngs_cuda::BufferDevPtr(*v.DevBufferRO()) + v.DevOffset();
  }

  // kernel reads and writes
  inline FlatVector<Dev<double>> FVDev (const DeviceVector<double> & v)
  {
    auto ptr = (Dev<double>*)ngs_cuda::BufferDevPtr(*v.DevBufferRW()) + v.DevOffset();
    return { v.Size(), ptr };
  }

  // kernel only reads
  inline FlatVector<Dev<double>> FVDevRO (const DeviceVector<double> & v)
  {
    return { v.Size(), DevPtr(v) };
  }


  /* AutoVector CreateUnifiedVector(size_t size); */

  class DevMatrix : public BaseMatrix
  {
  public:
    DevMatrix() { }

    AutoVector CreateRowVector() const override { return make_unique<UnifiedVector>(Width()); }
    AutoVector CreateColVector() const override { return make_unique<UnifiedVector>(Height()); }
  };

  shared_ptr<BaseMatrix> CreateDevMatrix (BaseMatrix &mat);


  class DevSparseMatrix : public DevMatrix
  {
  protected:
    //cusparseMatDescr_t * descr;
    cusparseSpMatDescr_t descr;
    int * dev_ind;
    int * dev_col;
    double * dev_val;
    int height, width, nze;
    size_t spmv_bufferSize = 0;
    void*  spmv_buffer = nullptr;
  public:
    DevSparseMatrix () { }
    DevSparseMatrix (const SparseMatrix<double> & mat);
    virtual ~DevSparseMatrix ();

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;

    virtual int VHeight() const { return height; }
    virtual int VWidth() const { return width; }
  };


  class DevDiagonalMatrix : public DevMatrix
  {
  protected:
    const UnifiedVector diag;

  public:
    DevDiagonalMatrix (const UnifiedVector _diag) : diag(_diag) { }

    virtual xbool IsSymmetric() const { return true; }

    virtual void Mult (const BaseVector & x, BaseVector & y) const;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;

    virtual int VHeight() const { return diag.Size(); }
    virtual int VWidth() const { return diag.Size(); }
  };


  


  
  






  // DevCGSolver — preconditioned CG with GPU-resident scalars (UnifiedScalar)
  // Supports CUDA graph capture of CG iteration body (convergence check via DtoH remains outside graph).
  class DevCGSolver : public KrylovSpaceSolver
  {
    shared_ptr<BaseMatrix> a_dev;  // DevSparseMatrix for graph capture
    shared_ptr<BaseMatrix> c_dev;  // preconditioner for graph capture

  public:
    DevCGSolver() : KrylovSpaceSolver() { }

    DevCGSolver(shared_ptr<BaseMatrix> mat,
                shared_ptr<BaseMatrix> pre)
      : KrylovSpaceSolver(mat, pre) { }

    DevCGSolver(shared_ptr<BaseMatrix> mat,
                shared_ptr<BaseMatrix> pre,
                shared_ptr<BaseMatrix> adev_raw,
                shared_ptr<BaseMatrix> cdev_raw)
      : KrylovSpaceSolver(mat, pre),
        a_dev(adev_raw), c_dev(cdev_raw) { }

    void Mult(const BaseVector& rhs,
              BaseVector& sol) const override;
  };

  // DevTFQMRSolver — preconditioned TFQMR for non-symmetric systems
  // Per-iteration CUDA graph capture: two graphs (even/odd), alternated each step.
  class DevTFQMRSolver : public KrylovSpaceSolver
  {
    shared_ptr<BaseMatrix> a_dev;
    shared_ptr<BaseMatrix> c_dev;

  public:
    DevTFQMRSolver() : KrylovSpaceSolver() { }

    DevTFQMRSolver(shared_ptr<BaseMatrix> mat,
                   shared_ptr<BaseMatrix> pre,
                   shared_ptr<BaseMatrix> adev_raw,
                   shared_ptr<BaseMatrix> cdev_raw)
      : KrylovSpaceSolver(mat, pre),
        a_dev(adev_raw), c_dev(cdev_raw) { }

    void Mult(const BaseVector& rhs,
              BaseVector& sol) const override;
  };
}


#endif
