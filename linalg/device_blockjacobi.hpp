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
*/

#include "devicevector.hpp"
#include "blockjacobi.hpp"

namespace ngla
{

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
    ngs_gpu::TypedBuffer<T> dev_mats;           // inverses, row-major, block after block
    int lanes;                                  // work-items per block
    bool overlapping;                           // a dof in more than one block

    void Launch (const BaseVector & x, BaseVector & y, T s, bool trans) const;

  public:
    template <typename TM>
    DeviceBlockJacobi (const BlockJacobiPrecond<TM> & pre);
    virtual ~DeviceBlockJacobi () { }

    virtual int VHeight() const override { return height; }
    virtual int VWidth() const override { return width; }
    virtual bool IsComplex() const override { return false; }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;

    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    virtual ostream & Print (ostream & ost) const override;
  };


#if !defined(FILE_DEVICE_BLOCKJACOBI_CPP)
  extern template class DeviceBlockJacobi<double>;
  extern template class DeviceBlockJacobi<float>;
#endif
}

#endif
