#ifndef FILE_DEVICE_SPARSECHOLESKY_HPP
#define FILE_DEVICE_SPARSECHOLESKY_HPP

/*********************************************************************/
/* File:   device_sparsecholesky.hpp                                 */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   4. Sep. 2026                                              */
/*********************************************************************/

/*
  Backend-independent solver with the sparse Cholesky factorization
  computed on the host (SparseCholeskyTM). The factor, the block
  structure and the micro-task dependency graph are uploaded once; a
  solve is: reorder, forward solve, diagonal scale, backward solve,
  reorder-add - each a kernel in the common syntax.

  The two triangular solves are persistent kernels: every simd group
  takes micro-tasks from a global counter in topological order and
  spin-waits on the task's dependency counter, as in ngscuda.
*/

#include "devicevector.hpp"
#include "sparsecholesky.hpp"

namespace ngla
{

  template <typename T>
  class NGS_DLL_HEADER DeviceSparseCholesky : public BaseMatrix
  {
  protected:
    size_t height, width, nused, njobs;
    MemType memtype;
    shared_ptr<ngs_gpu::Device> device;
    shared_ptr<ngs_gpu::Queue> queue;

    ngs_gpu::TypedBuffer<int> dev_microtasks;    // 4 ints per task: blocknr, type, bblock, nbblocks
    ngs_gpu::TypedBuffer<int> dev_depidx, dev_depdata;          // dependency table (csr)
    ngs_gpu::TypedBuffer<int> dev_depidx_trans, dev_depdata_trans;
    ngs_gpu::TypedBuffer<int> dev_incomingdep0, dev_incomingdep0_trans;   // initial counters
    ngs_gpu::TypedBuffer<int> dev_incomingdep, dev_incomingdep_trans;     // working counters
    ngs_gpu::TypedBuffer<int> dev_cnt;           // job counters of the two solves

    ngs_gpu::TypedBuffer<int> dev_blocks;        // block i has dofs [blocks[i], blocks[i+1])
    ngs_gpu::TypedBuffer<int> dev_rowindex2;     // row indices, one set per block
    ngs_gpu::TypedBuffer<int> dev_firstinrow_ri; // into rowindex2
    ngs_gpu::TypedBuffer<int> dev_firstinrow;    // into lfact
    ngs_gpu::TypedBuffer<T> dev_lfact;           // L factor, compressed
    ngs_gpu::TypedBuffer<T> dev_diag;            // inverse diagonal
    ngs_gpu::TypedBuffer<int> dev_order;         // dof i -> order[i], -1 if unused
    ngs_gpu::TypedBuffer<T> dev_hx;              // work vector in factor numbering

  public:
    template <typename TM>
    DeviceSparseCholesky (const SparseCholeskyTM<TM> & mat);
    virtual ~DeviceSparseCholesky () { }

    virtual int VHeight() const override { return height; }
    virtual int VWidth() const override { return width; }
    virtual bool IsComplex() const override { return false; }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    { MultAdd (s, x, y); }

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;

    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    virtual ostream & Print (ostream & ost) const override;
  };


#if !defined(FILE_DEVICE_SPARSECHOLESKY_CPP)
  extern template class DeviceSparseCholesky<double>;
  extern template class DeviceSparseCholesky<float>;
#endif
}

#endif
