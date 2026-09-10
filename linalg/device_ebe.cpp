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

    static Timer tscan("DeviceEBEMatrix ctor scan"), talloc("DeviceEBEMatrix ctor alloc"),
      tpack("DeviceEBEMatrix ctor pack"), tgemv("DeviceEBEMatrix ctor gemv"),
      ttrans("DeviceEBEMatrix ctor transpose");
    // input dofs are the columns, output dofs the rows of the element matrix
    BlockGemvBuilder<T> builder;
    {
      tscan.Start();
      size_t n = mat.GetNumElMats();
      Array<size_t> used(n);
      builder.infirst.SetSize (n+1);
      builder.outfirst.SetSize (n+1);
      builder.matfirst.SetSize (n+1);
      builder.infirst[0] = builder.outfirst[0] = builder.matfirst[0] = 0;
      size_t nb = 0;
      for (size_t i = 0; i < n; i++)
        {
          if (!ElementUsed (mat, i)) continue;
          size_t nin = mat.GetElementColumnDNums(i).Size();
          size_t nout = mat.GetElementRowDNums(i).Size();
          auto m = mat.GetElementMatrix(i);
          if (size_t(m.Height()) != nout || size_t(m.Width()) != nin)
            throw Exception("DeviceEBEMatrix: element matrix does not match its dof lists");
          used[nb] = i;
          builder.infirst[nb+1] = builder.infirst[nb] + int(nin);
          builder.outfirst[nb+1] = builder.outfirst[nb] + int(nout);
          builder.matfirst[nb+1] = builder.matfirst[nb] + int(nin*nout);
          nb++;
        }
      used.SetSize (nb);
      builder.infirst.SetSize (nb+1);
      builder.outfirst.SetSize (nb+1);
      builder.matfirst.SetSize (nb+1);
      tscan.Stop();

      talloc.Start();
      builder.inidx.SetSize (builder.infirst[nb]);
      builder.outidx.SetSize (builder.outfirst[nb]);
      builder.mats.SetSize (builder.matfirst[nb]);
      talloc.Stop();

      RegionTimer rpack(tpack);
      ParallelFor (nb, [&] (size_t j)
        {
          size_t i = used[j];
          auto rdi = mat.GetElementRowDNums(i);
          auto cdi = mat.GetElementColumnDNums(i);
          auto m = mat.GetElementMatrix(i);
          size_t nin = cdi.Size(), nout = rdi.Size();
          for (size_t k = 0; k < nin; k++) builder.inidx[builder.infirst[j]+k] = cdi[k];
          for (size_t k = 0; k < nout; k++) builder.outidx[builder.outfirst[j]+k] = rdi[k];
          T * dst = builder.mats.Data() + builder.matfirst[j];
          for (size_t c = 0; c < nin; c++)          // builder stores column-major
            for (size_t r = 0; r < nout; r++)
              dst[c*nout + r] = T(m(r,c));
        });
    }

    {
      RegionTimer r(tgemv);
      gemv = make_shared<DeviceBlockGemv<T>> (device, builder, width, height);
    }
    {
      RegionTimer r(ttrans);
      gemv_trans = gemv->Transpose();
    }

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
