#ifndef FILE_DEVICE_BLOCKJACOBI_HPP
#define FILE_DEVICE_BLOCKJACOBI_HPP

/*********************************************************************/
/* File:   device_blockjacobi.hpp                                    */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   3. Sep. 2026                                              */
/*********************************************************************/

/*
  Backend-independent block-Jacobi preconditioner on the gpu. The block
  inverses are computed on the host by BlockJacobiPrecond, this class
  only applies them: y += s * sum_blocks P_b^T inv_b P_b x.

  Created by BlockJacobiPrecond::CreateDeviceMatrix. Overlapping blocks
  accumulate with atomic adds, disjoint blocks with plain stores.

  Empty blocks are dropped. Small blocks share a few lanes per block,
  large ones take a whole warp with the block's x staged in group
  memory. Symmetric inverses make the transpose free.
*/

#include "devicevector.hpp"
#include "blockjacobi.hpp"

namespace ngla
{
  template <typename T> class DeviceBlockJacobiWarpKernels;

  template <typename T>
  class NGS_DLL_HEADER DeviceBlockJacobi : public BaseMatrix
  {
  protected:
    size_t height, width, nblocks;
    MemType memtype;
    shared_ptr<ngs_gpu::Device> device;
    shared_ptr<ngs_gpu::Queue> queue;

    ngs_gpu::TypedBuffer<int> dev_blockfirst;   // nblocks+1, into dev_indices
    ngs_gpu::TypedBuffer<int> dev_indices;      // dofs of all blocks
    ngs_gpu::TypedBuffer<int> dev_matfirst;     // nblocks+1, into dev_mats
    ngs_gpu::TypedBuffer<T> dev_mats;           // inverses, block after block
    ngs_gpu::TypedBuffer<int> dev_small, dev_large;   // block numbers per size class
    size_t nsmall, nlarge;
    int maxbs;
    int lanes;                                  // work-items per small block
    bool overlapping;                           // a dof in more than one block
    bool symmetric, colmajor;
    shared_ptr<const DeviceBlockJacobiWarpKernels<T>> warpkern;

    void Launch (const BaseVector & x, BaseVector & y, T s, bool trans) const;

  public:
    template <typename TM>
    DeviceBlockJacobi (const BlockJacobiPrecond<TM> & pre);
    virtual ~DeviceBlockJacobi () { }

    virtual int VHeight() const override { return height; }
    virtual int VWidth() const override { return width; }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    VecFormat RowFormat () const override { return DeviceVectorFormat<T> (width, memtype); }
    VecFormat ColFormat () const override { return DeviceVectorFormat<T> (height, memtype); }

    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    virtual ostream & Print (ostream & ost) const override;
  };


#if !defined(FILE_DEVICE_BLOCKJACOBI_CPP)
  extern template class DeviceBlockJacobi<double>;
  extern template class DeviceBlockJacobi<float>;
#endif
}

#endif
