#define FILE_DEVICE_EBE_CPP
#include <la.hpp>

namespace ngla
{
  using namespace ngs_gpu;

  namespace
  {
    // reserved but unused element slots carry a -1 entry, as in the host MultAdd
    template <typename TM>
    bool ElementUsed (const ElementByElementMatrix<TM> & mat, size_t i)
    {
      auto rdi = mat.GetElementRowDNums(i);
      auto cdi = mat.GetElementColumnDNums(i);
      return rdi.Size() && cdi.Size() && rdi[0] != -1 && cdi[0] != -1;
    }
  }


  template <typename T>
  template <typename TM>
  DeviceEBEMatrix<T> :: DeviceEBEMatrix (const ElementByElementMatrix<TM> & mat)
    : memtype (PreferredMemType())
  {
    device = GetGpuDevice();
    height = mat.Height();
    width = mat.Width();

    // input dofs are the columns, output dofs the rows of the element matrix
    BlockGemvBuilder<T> builder;
    for (size_t i = 0; i < mat.GetNumElMats(); i++)
      {
        if (!ElementUsed (mat, i)) continue;
        auto rdi = mat.GetElementRowDNums(i);
        auto cdi = mat.GetElementColumnDNums(i);
        auto m = mat.GetElementMatrix(i);
        if (m.Height() != rdi.Size() || m.Width() != cdi.Size())
          throw Exception("DeviceEBEMatrix: element matrix does not match its dof lists");
        builder.AddBlock (cdi, rdi, [&] (size_t r, size_t c) { return m(r,c); });
      }

    gemv = make_shared<DeviceBlockGemv<T>> (device, builder, width, height);
    gemv_trans = gemv->Transpose();

    cout << IM(7) << "DeviceEBEMatrix<" << (is_same_v<T,double> ? "double" : "float")
         << "> " << gemv->Info() << endl;
  }

  template <typename T>
  void DeviceEBEMatrix<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceEBEMatrix::MultAdd"); RegionTimer reg(t);
    if (x.Size() != width || y.Size() != height)
      throw Exception("DeviceEBEMatrix::MultAdd - size mismatch");
    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    gemv->MultAdd (T(s), ux.DevArgRO(), uy.DevArgRW());
  }

  template <typename T>
  void DeviceEBEMatrix<T> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceEBEMatrix::MultTransAdd"); RegionTimer reg(t);
    if (x.Size() != height || y.Size() != width)
      throw Exception("DeviceEBEMatrix::MultTransAdd - size mismatch");
    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    gemv_trans->MultAdd (T(s), ux.DevArgRO(), uy.DevArgRW());
  }

  template <typename T>
  BaseMatrix::OperatorInfo DeviceEBEMatrix<T> :: GetOperatorInfo () const
  {
    return { string("DeviceEBEMatrix<") + (is_same_v<T,double> ? "double" : "float")
             + "> (elements=" + ToString(gemv->NBlocks()) + ")", height, width };
  }

  template <typename T>
  ostream & DeviceEBEMatrix<T> :: Print (ostream & ost) const
  {
    ost << "DeviceEBEMatrix<" << (is_same_v<T,double> ? "double" : "float")
        << ">, " << height << " x " << width << ", " << gemv->Info()
        << ", on " << device->Name() << endl;
    return ost;
  }


  template <typename SCAL>
  shared_ptr<BaseMatrix> ElementByElementMatrix<SCAL> :: CreateDeviceMatrix () const
  {
    if constexpr (is_same_v<SCAL,double> || is_same_v<SCAL,float>)
      if (ngs_gpu::HasDevice())
        {
          if constexpr (is_same_v<SCAL,double>)
            if (GetGpuDevice()->HasFloat64())
              return make_shared<DeviceEBEMatrix<double>> (*this);
          return make_shared<DeviceEBEMatrix<float>> (*this);
        }
    return BaseMatrix::CreateDeviceMatrix();
  }

  template shared_ptr<BaseMatrix> ElementByElementMatrix<double>::CreateDeviceMatrix () const;
  template shared_ptr<BaseMatrix> ElementByElementMatrix<float>::CreateDeviceMatrix () const;
  template shared_ptr<BaseMatrix> ElementByElementMatrix<Complex>::CreateDeviceMatrix () const;

  template class DeviceEBEMatrix<double>;
  template class DeviceEBEMatrix<float>;
  template DeviceEBEMatrix<double>::DeviceEBEMatrix (const ElementByElementMatrix<double>&);
  template DeviceEBEMatrix<double>::DeviceEBEMatrix (const ElementByElementMatrix<float>&);
  template DeviceEBEMatrix<float>::DeviceEBEMatrix (const ElementByElementMatrix<double>&);
  template DeviceEBEMatrix<float>::DeviceEBEMatrix (const ElementByElementMatrix<float>&);
}
