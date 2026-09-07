#ifndef FILE_DEVICE_EBE_HPP
#define FILE_DEVICE_EBE_HPP

#include "devicevector.hpp"
#include "elementbyelement.hpp"

namespace ngla
{
  template <typename T>
  class NGS_DLL_HEADER DeviceEBEMatrix : public BaseMatrix
  {
  protected:
    size_t height, width, nel;
    MemType memtype;
    shared_ptr<ngs_gpu::Device> device;
    shared_ptr<ngs_gpu::Queue> queue;
    ngs_gpu::TypedBuffer<int> dev_rowfirst, dev_colfirst, dev_matfirst;   // nel+1 each
    ngs_gpu::TypedBuffer<int> dev_rowidx, dev_colidx;                     // dofs of all elements
    ngs_gpu::TypedBuffer<T> dev_mats;                                     // row-major, element after element
    int lanes, lanes_trans;                     // work-items per element

    void Launch (const BaseVector & x, BaseVector & y, T s, bool trans) const;

  public:
    template <typename TM>
    DeviceEBEMatrix (const ElementByElementMatrix<TM> & mat);
    virtual ~DeviceEBEMatrix () { }

    virtual int VHeight() const override { return height; }
    virtual int VWidth() const override { return width; }
    virtual bool IsComplex() const override { return false; }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
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
