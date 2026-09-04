#ifndef FILE_DEVICE_SPARSEMATRIX_HPP
#define FILE_DEVICE_SPARSEMATRIX_HPP

/*********************************************************************/
/* File:   device_sparsematrix.hpp                                   */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   3. Sep. 2026                                              */
/*********************************************************************/

/*
  Backend-independent csr matrix on the gpu, built on ngs_gpu like
  DeviceVector. Created from a host SparseMatrix by
  SparseMatrix::CreateDeviceMatrix, which picks T = double where the
  device has fp64 and float otherwise.

  Storage is three device buffers (row starts, column indices, values),
  the spmv kernels are written in the common syntax of ngstd/gpukernel.hpp
  and compiled at runtime. The transposed structure is built on the
  device's first MultTransAdd only.
*/

#include "devicevector.hpp"
#include "sparsematrix.hpp"

namespace ngla
{

  template <typename T>
  class NGS_DLL_HEADER DeviceSparseMatrix : public BaseMatrix
  {
  protected:
    size_t height, width, nze;
    MemType memtype;
    shared_ptr<ngs_gpu::Device> device;
    shared_ptr<ngs_gpu::Queue> queue;

    ngs_gpu::TypedBuffer<int> dev_firsti, dev_colnr;   // int32 indices
    ngs_gpu::TypedBuffer<T> dev_values;

    // kernel and work-items per row, chosen by timing on this matrix
    struct SpMVChoice { shared_ptr<ngs_gpu::Kernel> kernel; int lanes = 1, rows_per_group = 1; };
    SpMVChoice choice;

    // transposed csr, built on demand
    mutable ngs_gpu::TypedBuffer<int> devt_firsti, devt_colnr;
    mutable ngs_gpu::TypedBuffer<T> devt_values;
    mutable SpMVChoice choice_trans;
    mutable std::mutex trans_mutex;

    void BuildTranspose() const;
    SpMVChoice AutoTune (const ngs_gpu::TypedBuffer<int> & firsti,
                         const ngs_gpu::TypedBuffer<int> & colnr,
                         const ngs_gpu::TypedBuffer<T> & values,
                         size_t rows, size_t cols) const;
    void LaunchSpMV (const ngs_gpu::TypedBuffer<int> & firsti,
                     const ngs_gpu::TypedBuffer<int> & colnr,
                     const ngs_gpu::TypedBuffer<T> & values,
                     const SpMVChoice & ch, size_t rows,
                     ngs_gpu::KernelArg x, ngs_gpu::KernelArg y, T s) const;

  public:
    // values are converted to T
    template <typename TM>
    DeviceSparseMatrix (const SparseMatrixTM<TM> & mat);
    virtual ~DeviceSparseMatrix () { }

    virtual int VHeight() const override { return height; }
    virtual int VWidth() const override { return width; }
    virtual size_t NZE () const override { return nze; }
    virtual bool IsComplex() const override { return false; }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;

    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    virtual ostream & Print (ostream & ost) const override;

    // for kernels and backend libraries
    const ngs_gpu::TypedBuffer<int> & DevFirstI() const { return dev_firsti; }
    const ngs_gpu::TypedBuffer<int> & DevColNr() const { return dev_colnr; }
    const ngs_gpu::TypedBuffer<T> & DevValues() const { return dev_values; }
    MemType GetMemType() const { return memtype; }
    const shared_ptr<ngs_gpu::Queue> & GetQueue() const { return queue; }
  };


#if !defined(FILE_DEVICE_SPARSEMATRIX_CPP)
  extern template class DeviceSparseMatrix<double>;
  extern template class DeviceSparseMatrix<float>;
#endif
}

#endif
