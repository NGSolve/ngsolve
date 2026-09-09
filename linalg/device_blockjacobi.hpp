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
  only applies them: y += s * sum_blocks P_b^T inv_b P_b x, by the
  batched block gemv of DeviceBlockGemv. Symmetric inverses make the
  transpose the operator itself.

  Created by BlockJacobiPrecond::CreateDeviceMatrix.
*/

#include "device_blockgemv.hpp"
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
    shared_ptr<DeviceBlockGemv<T>> gemv, gemv_trans;   // same object if symmetric
    bool symmetric;

  public:
    template <typename TM>
    DeviceBlockJacobi (const BlockJacobiPrecond<TM> & pre);
    template <typename TM, typename TV>
    DeviceBlockJacobi (const BlockJacobiPrecondSymmetric<TM,TV> & pre);
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
