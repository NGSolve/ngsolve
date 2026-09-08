#ifndef FILE_DEVICE_BLOCKGS_HPP
#define FILE_DEVICE_BLOCKGS_HPP

/*********************************************************************/
/* File:   device_blockgs.hpp                                        */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   8. Sep. 2026                                              */
/*********************************************************************/

/*
  Backend-independent block Gauss-Seidel smoother on the gpu, using the
  block inverses and the block coloring of a BlockJacobiPrecond. One
  kernel launch per color: a warp per block computes the block residual
  b - A x from the csr rows, applies the inverse, and updates x. The
  backward sweep runs the colors in reverse.

  Created by SymmetricBlockGaussSeidelPrecond::CreateDeviceMatrix (Mult
  is a forward and a backward sweep from zero) and usable as a smoother
  through Smooth / SmoothBack.
*/

#include "devicevector.hpp"
#include "device_sparsematrix.hpp"
#include "blockjacobi.hpp"

namespace ngla
{
  template <typename T> class DeviceBlockGSKernels;

  template <typename T>
  class NGS_DLL_HEADER DeviceBlockGaussSeidel : public BaseMSMPrecond
  {
  protected:
    size_t height, nblocks;
    MemType memtype;
    shared_ptr<ngs_gpu::Device> device;
    shared_ptr<ngs_gpu::Queue> queue;

    shared_ptr<DeviceSparseMatrix<T>> devmat;   // csr of A
    ngs_gpu::TypedBuffer<int> dev_blockfirst;   // nblocks+1, into dev_indices
    ngs_gpu::TypedBuffer<int> dev_indices;      // dofs of all blocks
    ngs_gpu::TypedBuffer<int> dev_matfirst;     // nblocks+1, into dev_mats
    ngs_gpu::TypedBuffer<T> dev_mats;           // inverses, column-major, block after block
    ngs_gpu::TypedBuffer<int> dev_colorblocks;  // non-empty block numbers, color after color
    Array<int> colorfirst;                      // ncolors+1, into dev_colorblocks
    Array<int> colorsplit;                      // ncolors, first large block of the color
    int maxbs;
    shared_ptr<const DeviceBlockGSKernels<T>> kern_small, kern_large;   // one warp / one group per block

    void Sweep (BaseVector & x, const BaseVector & b, bool backward) const;

  public:
    template <typename TM>
    DeviceBlockGaussSeidel (const BlockJacobiPrecond<TM> & pre);
    virtual ~DeviceBlockGaussSeidel () { }

    int VHeight() const override { return height; }
    int VWidth() const override { return height; }
    VecFormat RowFormat () const override { return DeviceVectorFormat<T> (height, memtype); }
    VecFormat ColFormat () const override { return DeviceVectorFormat<T> (height, memtype); }
    bool IsComplex() const override { return false; }

    void Smooth (BaseVector & x, const BaseVector & b, int steps = 1) const override;
    void SmoothBack (BaseVector & x, const BaseVector & b, int steps = 1) const override;

    // symmetric Gauss-Seidel as preconditioner
    void Mult (const BaseVector & x, BaseVector & y) const override;
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    { MultAdd (s, x, y); }

    BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    ostream & Print (ostream & ost) const override;
  };


#if !defined(FILE_DEVICE_BLOCKGS_CPP)
  extern template class DeviceBlockGaussSeidel<double>;
  extern template class DeviceBlockGaussSeidel<float>;
#endif
}

#endif
