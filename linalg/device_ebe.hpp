#ifndef FILE_DEVICE_EBE_HPP
#define FILE_DEVICE_EBE_HPP

/*
  Backend-independent ElementByElementMatrix on the gpu: the element
  matrices with their individual dof lists, applied by the batched block
  gemv of DeviceBlockGemv (the transpose is a view of the same data).
*/

#include "device_blockgemv.hpp"
#include "elementbyelement.hpp"

namespace ngla
{
  template <typename T>
  class NGS_DLL_HEADER DeviceEBEMatrix : public BaseMatrix
  {
  protected:
    size_t height, width;
    MemType memtype;
    shared_ptr<ngs_gpu::Device> device;
    shared_ptr<DeviceBlockGemv<T>> gemv, gemv_trans;

  public:
    template <typename TM>
    DeviceEBEMatrix (const ElementByElementMatrix<TM> & mat);
    virtual ~DeviceEBEMatrix () { }

    virtual int VHeight() const override { return height; }
    virtual int VWidth() const override { return width; }
    virtual bool IsComplex() const override { return false; }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    VecFormat RowFormat () const override { return DeviceVectorFormat<T> (width, memtype); }
    VecFormat ColFormat () const override { return DeviceVectorFormat<T> (height, memtype); }

    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    virtual ostream & Print (ostream & ost) const override;
  };

#if !defined(FILE_DEVICE_EBE_CPP)
  extern template class DeviceEBEMatrix<double>;
  extern template class DeviceEBEMatrix<float>;
#endif
}

#endif
