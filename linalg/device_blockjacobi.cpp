/*********************************************************************/
/* File:   device_blockjacobi.cpp                                    */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   3. Sep. 2026                                              */
/*********************************************************************/

#define FILE_DEVICE_BLOCKJACOBI_CPP

#include <la.hpp>

namespace ngla
{
  using namespace ngs_gpu;

  template <typename T>
  template <typename TM>
  DeviceBlockJacobi<T> :: DeviceBlockJacobi (const BlockJacobiPrecond<TM> & pre)
  {
    height = pre.Height();
    width = pre.Width();
    const Table<int> & blocktable = *pre.GetBlockTable();
    nblocks = blocktable.Size();
    const auto & inverses = pre.GetInverses();

    device = GetGpuDevice();
    memtype = PreferredMemType();

    // symmetric inverses: the transpose is the matrix itself
    symmetric = true;
    BlockGemvBuilder<T> builder;
    for (size_t i = 0; i < nblocks; i++)
      {
        const auto & inv = inverses[i];
        size_t bs = blocktable[i].Size();
        double scale = 0, asym = 0;
        for (size_t r = 0; r < bs; r++)
          for (size_t c = 0; c < bs; c++)
            {
              scale = max (scale, double(fabs(inv(r,c))));
              asym = max (asym, double(fabs(inv(r,c)-inv(c,r))));
            }
        if (asym > 1e-10*scale) symmetric = false;
        builder.AddBlock (blocktable[i], blocktable[i], [&] (size_t r, size_t c) { return inv(r,c); });
      }

    gemv = make_shared<DeviceBlockGemv<T>> (device, builder, width, height);
    gemv_trans = symmetric ? gemv : gemv->Transpose();

    cout << IM(7) << "DeviceBlockJacobi<" << (is_same_v<T,double> ? "double" : "float")
         << "> " << gemv->Info() << (symmetric ? ", symmetric" : "") << endl;
  }


  template <typename T>
  void DeviceBlockJacobi<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceBlockJacobi::MultAdd"); RegionTimer reg(t);
    if (x.Size() != width || y.Size() != height)
      throw Exception("DeviceBlockJacobi::MultAdd - size mismatch");
    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    gemv->MultAdd (T(s), ux.DevArgRO(), uy.DevArgRW());
  }

  template <typename T>
  void DeviceBlockJacobi<T> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceBlockJacobi::MultTransAdd"); RegionTimer reg(t);
    if (x.Size() != height || y.Size() != width)
      throw Exception("DeviceBlockJacobi::MultTransAdd - size mismatch");
    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    gemv_trans->MultAdd (T(s), ux.DevArgRO(), uy.DevArgRW());
  }


  template <typename T>
  BaseMatrix::OperatorInfo DeviceBlockJacobi<T> :: GetOperatorInfo () const
  {
    return { string("DeviceBlockJacobi<") + (is_same_v<T,double> ? "double" : "float")
             + "> (blocks=" + ToString(nblocks) + ")", height, width };
  }

  template <typename T>
  ostream & DeviceBlockJacobi<T> :: Print (ostream & ost) const
  {
    ost << "DeviceBlockJacobi<" << (is_same_v<T,double> ? "double" : "float")
        << ">, height = " << height << ", " << gemv->Info() << ", on " << device->Name() << endl;
    return ost;
  }


  template class DeviceBlockJacobi<double>;
  template class DeviceBlockJacobi<float>;
  template DeviceBlockJacobi<double>::DeviceBlockJacobi (const BlockJacobiPrecond<double> &);
  template DeviceBlockJacobi<double>::DeviceBlockJacobi (const BlockJacobiPrecond<float> &);
  template DeviceBlockJacobi<float>::DeviceBlockJacobi (const BlockJacobiPrecond<double> &);
  template DeviceBlockJacobi<float>::DeviceBlockJacobi (const BlockJacobiPrecond<float> &);
}
