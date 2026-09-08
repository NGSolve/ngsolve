#ifndef FILE_DEVICE_EBECONST_HPP
#define FILE_DEVICE_EBECONST_HPP

/*********************************************************************/
/* File:   device_ebeconst.hpp                                       */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   7. Sep. 2026                                              */
/*********************************************************************/

/*
  Backend-independent ConstantElementByElementMatrix on the gpu:
  y += s * sum_blocks P_out^T M P_in x with one small matrix M for all
  blocks. One fused kernel per product: gather the block dofs of x,
  multiply with M (or M^T), scatter into y. Overlapping output dofs
  accumulate with atomic adds, disjoint ones with plain stores.

  On simd-32 devices a group of blocks is staged in group memory and
  multiplied in 8x8 warp tiles (tinybla gemm, compiled per shape); the
  host reference backend uses a plain per-lane kernel.

  Created by ConstantElementByElementMatrix::CreateDeviceMatrix.
*/

#include "devicevector.hpp"
#include "elementbyelement.hpp"

namespace ngla
{
  template <typename T> class DeviceEBEGemmKernels;

  template <typename T>
  class NGS_DLL_HEADER DeviceConstantEBEMatrix : public BaseMatrix
  {
  protected:
    size_t height, width, nblocks;
    int hm, wm;                                 // element matrix shape
    MemType memtype;
    shared_ptr<ngs_gpu::Device> device;
    shared_ptr<ngs_gpu::Queue> queue;

    ngs_gpu::TypedBuffer<T> dev_mat;            // M^T, row-major, zero padded to multiples of 8
    ngs_gpu::TypedBuffer<T> dev_mat_trans;      // M, padded likewise
    ngs_gpu::TypedBuffer<int> dev_rowdnums;     // nblocks x wm, input dofs
    ngs_gpu::TypedBuffer<int> dev_coldnums;     // nblocks x hm, output dofs
    bool disjoint_rows, disjoint_cols;
    bool onto_cols, onto_rows;                  // disjoint and covering: Mult may store
    int lanes, lanes_trans;                     // work-items per block, lane kernel
    shared_ptr<const DeviceEBEGemmKernels<T>> gemm, gemm_trans;   // null: lane kernel

    void Launch (const BaseVector & x, BaseVector & y, T s, T beta, bool trans) const;

  public:
    template <typename TM>
    DeviceConstantEBEMatrix (const ConstantElementByElementMatrix<TM> & mat);
    virtual ~DeviceConstantEBEMatrix () { }

    virtual int VHeight() const override { return height; }
    virtual int VWidth() const override { return width; }
    VecFormat RowFormat () const override { return DeviceVectorFormat<T> (width, memtype); }
    VecFormat ColFormat () const override { return DeviceVectorFormat<T> (height, memtype); }
    virtual bool IsComplex() const override { return false; }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;


    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    virtual ostream & Print (ostream & ost) const override;
  };


#if !defined(FILE_DEVICE_EBECONST_CPP)
  extern template class DeviceConstantEBEMatrix<double>;
  extern template class DeviceConstantEBEMatrix<float>;
#endif
}

#endif
