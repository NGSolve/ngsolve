#ifndef FILE_FMMOPERATOR
#define FILE_FMMOPERATOR
#include <mptools.hpp>
#include <bla.hpp>

namespace ngsbem
{
  
  inline double __Kernel(Vec<3> px, Vec<3> py)
  {
    if (L2Norm2(px-py)==0) return 0.0;
    return 1.0 / (4*M_PI) / L2Norm(px-py);
  }

  inline std::tuple<Vec<3>, double> GetCenterAndRadius(const Array<Vec<3>>& xpts)
  {
    Vec<3> xmax(-1e99, -1e99, -1e99);
    Vec<3> xmin(1e99, 1e99, 1e99);
    for (auto xi : xpts)
      {
        for (int j = 0; j < 3; j++)
          {
            xmin(j) = min(xmin(j), xi(j));
            xmax(j) = max(xmax(j), xi(j));
          }
      }

    Vec<3> cx = 0.5*(xmin+xmax);
    double rx = 0;
    for (int j = 0; j < 3; j++)
      rx = max(rx, xmax(j)-xmin(j));

    return { cx, rx };
  }
  
  
  template <typename TSCAL>
  class Base_FMM_Operator : public BaseMatrix
  {
  protected:
    Array<Vec<3>> xpts, ypts, xnv, ynv;
    Vec<3> cx, cy;
    double rx, ry;
    IVec<2> kernelshape;
  public:
    Base_FMM_Operator(Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                      Array<Vec<3>> _xnv, Array<Vec<3>> _ynv, IVec<2> _kernelshape)
      : xpts(std::move(_xpts)), ypts(std::move(_ypts)),
        xnv(std::move(_xnv)), ynv(std::move(_ynv)), kernelshape(_kernelshape)
    {
      std::tie(cx, rx) = GetCenterAndRadius(xpts);      
      std::tie(cy, ry) = GetCenterAndRadius(ypts);      
    }

    int VHeight() const override { return  ypts.Size()*kernelshape[0]; }
    int VWidth() const override { return  xpts.Size()*kernelshape[1]; }
      
    AutoVector CreateRowVector () const override
    {
      return make_unique<VVector<TSCAL>>(xpts.Size() * kernelshape[1]);
    }
    AutoVector CreateColVector () const override
    {
      return make_unique<VVector<TSCAL>>(ypts.Size() * kernelshape[0]);
    }
  };


  
  template <typename KERNEL>
  class FMM_Operator : public Base_FMM_Operator<typename KERNEL::value_type>
  {
    KERNEL kernel;
    typedef Base_FMM_Operator<typename KERNEL::value_type> BASE;
    using BASE::xpts, BASE::ypts, BASE::xnv, BASE::ynv, BASE::cx, BASE::cy, BASE::rx, BASE::ry;

  public:
    FMM_Operator(KERNEL _kernel, Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                 Array<Vec<3>> _xnv, Array<Vec<3>> _ynv)
      : BASE(std::move(_xpts), std::move( _ypts), std::move(_xnv), std::move(_ynv), KERNEL::Shape()),
      kernel(_kernel)
    { }

    void Mult(const BaseVector & x, BaseVector & y) const override
    {
      static Timer tall("ngbem fmm apply "+KERNEL::Name()); RegionTimer reg(tall);

      auto shape = KERNEL::Shape();
      
      // auto fx = x.FV<typename KERNEL::value_type>();
      // auto fy = y.FV<typename KERNEL::value_type>();
      auto matx = x.FV<typename KERNEL::value_type>().AsMatrix(xpts.Size(), shape[0]);
      auto maty = y.FV<typename KERNEL::value_type>().AsMatrix(ypts.Size(), shape[1]);      
      
      maty = 0;
      auto singmp = kernel.CreateMultipoleExpansion (cx, rx);
      ParallelFor (xpts.Size(), [&](int i){
        // kernel.AddSource(*singmp, xpts[i], xnv[i], make_BareSliceVector(fx.Range(shape[0]*i,shape[0]*(i+1))));
        kernel.AddSource(*singmp, xpts[i], xnv[i], matx.Row(i));
      });
      singmp->CalcMP();
      auto regmp = kernel.CreateLocalExpansion (cy, ry);
      ParallelFor (ypts.Size(), [&](int i){
        regmp->AddTarget(ypts[i]);
      });
      regmp->CalcMP(singmp);
      ParallelFor (ypts.Size(), [&](int i) {
        //kernel.EvaluateMP(*regmp, ypts[i], ynv[i], make_BareSliceVector(fy.Range(shape[1]*i,shape[1]*(i+1))));
        kernel.EvaluateMP(*regmp, ypts[i], ynv[i], maty.Row(i)); 
      });
    }

    void MultTrans(const BaseVector & x, BaseVector & y) const override
    {
        static Timer tall("ngbem fmm apply Trans "+KERNEL::Name()); RegionTimer reg(tall);

        auto shape = KERNEL::Shape();
        auto matx = x.FV<typename KERNEL::value_type>().AsMatrix(xpts.Size(), shape[0]);
        auto maty = y.FV<typename KERNEL::value_type>().AsMatrix(ypts.Size(), shape[1]);

        maty = 0;
        auto singmp = kernel.CreateMultipoleExpansion (cy, ry);
        ParallelFor (ypts.Size(), [&](int i){
          kernel.AddSourceTrans(*singmp, ypts[i], ynv[i], matx.Row(i));
        });
        singmp->CalcMP();
        auto regmp = kernel.CreateLocalExpansion (cx, rx);
        ParallelFor (xpts.Size(), [&](int i){
          regmp->AddTarget(xpts[i]);
        });
        regmp->CalcMP(singmp);
        ParallelFor (xpts.Size(), [&](int i) {
          kernel.EvaluateMPTrans(*regmp, xpts[i], xnv[i], maty.Row(i));
        });
    };
  };

}


#endif

