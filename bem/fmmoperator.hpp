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
      auto fx = x.FV<typename KERNEL::value_type>();
      auto fy = y.FV<typename KERNEL::value_type>();
      
      fy = 0;
      
      if constexpr (std::is_same<KERNEL, class LaplaceDLKernel<3>>())
        {
          // This is slow, but works for now:
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
              {
                double norm = L2Norm(xpts[ix]-ypts[iy]);
                if (norm > 0)
                  {
                    // double nxy = InnerProduct(ynv[iy], xpts[ix]-ypts[iy]);                 
                    double nxy = InnerProduct(xnv[ix], ypts[iy]-xpts[ix]);
                    auto kern = nxy / (4 * M_PI * norm*norm*norm);
                    fy(ix) += kern * fx(iy);
                  }
              }
        }
      else
        throw Exception("KERNEL " + std::string(typeid(KERNEL).name()) + " doesn't have MultTrans overloaded");
    };

  };

  // ********************** operators for Laplace ********************
  // using Helmholtz, with small kappa


  template <>
  void FMM_Operator<LaplaceDLKernel<3>> :: MultTrans(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply LaplaceDL MultTrans (ngfmm)"); RegionTimer reg(tall);
    auto fx = x.FV<double>();
    auto fy = y.FV<double>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    
    double kappa = 1e-16;
    auto singmp = make_shared<SingularMLMultiPole<Complex>>(cy, ry, kappa);

    ParallelFor (ypts.Size(), [&](int i){
      singmp->AddCharge(ypts[i], fx(i));
    });
    singmp->CalcMP();

    RegularMLMultiPole<Complex> regmp (cx, rx, kappa);
    ParallelFor (xpts.Size(), [&](int i){
      regmp.AddTarget(xpts[i]);
    });
    regmp.CalcMP(singmp);

    ParallelFor (xpts.Size(), [&](int i) {
      fy(i) = Real (regmp.EvaluateDirectionalDerivative(xpts[i], xnv[i]));
    });
  }
  


  // ********************** operators for Helmholtz ********************  
  

  template <>
  void FMM_Operator<HelmholtzDLKernel<3>> :: MultTrans(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply HelmholtzDL MultTrans (ngfmm)"); RegionTimer reg(tall);
    auto fx = x.FV<Complex>();
    auto fy = y.FV<Complex>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    double kappa = kernel.GetKappa();
    auto singmp = make_shared<SingularMLMultiPole<Complex>>(cy, ry, kappa);

    ParallelFor (ypts.Size(), [&](int i){
      singmp->AddCharge(ypts[i], fx(i));
    });
    singmp->CalcMP();

    RegularMLMultiPole<Complex> regmp (cx, rx, kappa);
    ParallelFor (xpts.Size(), [&](int i){
      regmp.AddTarget(xpts[i]);
    });
    regmp.CalcMP(singmp);

    ParallelFor (xpts.Size(), [&](int i) {
      fy(i) = regmp.EvaluateDirectionalDerivative(xpts[i], xnv[i]);
    });
  }

  // https://weggler.github.io/ngbem/short_and_sweet/Maxwell_Formulations.html
  template <>
  void FMM_Operator<MaxwellDLKernel<3>> :: MultTrans(const BaseVector & x, BaseVector & y) const
  {
    static Timer tall("ngbem fmm apply MaxwellDL MultTrans (ngfmm)"); RegionTimer reg(tall);
    auto fx = x.FV<Complex>();
    auto fy = y.FV<Complex>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    double kappa = kernel.GetKappa();

    auto singmp = make_shared<SingularMLMultiPole<Vec<3,Complex>>>(cy, ry, kappa);

    ParallelFor (ypts.Size(), [&](int i){
      Vec<3,Complex> current = fx.Range(3*i, 3*i+3);

      for (int k = 0; k < 3; k++)
      {
        Vec<3> ek{0.0}; ek(k) = 1;
        Vec<3> current_real = Real(current);
        Vec<3> current_imag = Imag(current);

        singmp->AddDipole(xpts[i], Cross(current_real, ek), ek);
        singmp->AddDipole(xpts[i], Cross(current_imag, ek), Complex(0,1)*ek);
      }
    });
    singmp->CalcMP();

    RegularMLMultiPole<Vec<3,Complex>> regmp (cx, rx, kappa);
    ParallelFor (xpts.Size(), [&](int i){
      regmp.AddTarget(xpts[i]);
    });
    regmp.CalcMP(singmp);

    ParallelFor (xpts.Size(), [&](int i) {
      fy.Range(3*i, 3*i+3) = regmp.Evaluate(xpts[i]);
    });
  }
}


#endif

