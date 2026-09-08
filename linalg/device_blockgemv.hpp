#ifndef FILE_DEVICE_BLOCKGEMV_HPP
#define FILE_DEVICE_BLOCKGEMV_HPP

/*********************************************************************/
/* File:   device_blockgemv.hpp                                      */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   9. Sep. 2026                                              */
/*********************************************************************/

/*
  Batched dense block gemv with gathered input and scattered output,

      y[out_b[r]] += s * sum_c M_b(r,c) x[in_b[c]]     for all blocks b,

  the common kernel of DeviceBlockJacobi and DeviceEBEMatrix. Blocks
  have individual sizes; empty ones are dropped. Blocks with up to 16
  output rows share a few lanes per block, larger ones take a group with
  the block's x staged in group memory and up to four rows per lane.
  Output dofs shared between blocks accumulate atomically, disjoint ones
  with plain adds. The matrices are stored column-major w.r.t. the
  output index; Transpose() is a view with in/out swapped.
*/

#include "devicevector.hpp"

namespace ngla
{
  template <typename T> class DeviceBlockGemvKernels;

  // host-side assembly of the blocks
  template <typename T>
  struct BlockGemvBuilder
  {
    Array<int> infirst = { 0 }, outfirst = { 0 }, matfirst = { 0 };
    Array<int> inidx, outidx;
    Array<T> mats;

    // mat(r,c) with r over out, c over in
    template <typename FUNC>
    void AddBlock (FlatArray<int> in, FlatArray<int> out, FUNC mat)
    {
      if (in.Size() == 0 || out.Size() == 0) return;
      inidx += in;
      outidx += out;
      size_t nin = in.Size(), nout = out.Size();
      size_t base = mats.Size();
      mats.SetSize (base + nin*nout);
      for (size_t c = 0; c < nin; c++)
        for (size_t r = 0; r < nout; r++)
          mats[base + c*nout + r] = T(mat(r,c));
      infirst.Append (inidx.Size());
      outfirst.Append (outidx.Size());
      matfirst.Append (mats.Size());
    }
    size_t NBlocks() const { return infirst.Size()-1; }
  };


  template <typename T>
  class NGS_DLL_HEADER DeviceBlockGemv
  {
    shared_ptr<ngs_gpu::Device> device;
    shared_ptr<ngs_gpu::Queue> queue;
    size_t nblocks, nsmall, nlarge;
    int lanes, maxin;
    bool strided, atomic;
    Array<int> nin, nout;                       // host copies of the block shapes

    ngs_gpu::TypedBuffer<int> dev_infirst, dev_outfirst, dev_matfirst;
    ngs_gpu::TypedBuffer<int> dev_inidx, dev_outidx;
    ngs_gpu::TypedBuffer<T> dev_mats;
    ngs_gpu::TypedBuffer<int> dev_small, dev_large;    // block numbers per size class
    shared_ptr<const DeviceBlockGemvKernels<T>> kern_large;

    void Classify ();
    DeviceBlockGemv () = default;

  public:
    DeviceBlockGemv (shared_ptr<ngs_gpu::Device> device, const BlockGemvBuilder<T> & b, size_t width, size_t height);

    // the transposed operator, sharing the data
    shared_ptr<DeviceBlockGemv> Transpose () const;

    size_t NBlocks() const { return nblocks; }
    bool IsAtomic() const { return atomic; }
    string Info() const;

    // y += s * A x, both device-side arguments
    void MultAdd (T s, ngs_gpu::KernelArg x, ngs_gpu::KernelArg y) const;
  };


#if !defined(FILE_DEVICE_BLOCKGEMV_CPP)
  extern template class DeviceBlockGemv<double>;
  extern template class DeviceBlockGemv<float>;
#endif
}

#endif
