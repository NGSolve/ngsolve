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

  inline std::tuple<Vec<3>, double, Vec<3>, double> GetMidAndRadius(const Array<Vec<3>>& xpts, const Array<Vec<3>>& ypts)
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

    Vec<3> ymax(-1e99, -1e99, -1e99);
    Vec<3> ymin(1e99, 1e99, 1e99);
    for (auto yi : ypts)
      {
        for (int j = 0; j < 3; j++)
          {
            ymin(j) = min(ymin(j), yi(j));
            ymax(j) = max(ymax(j), yi(j));
          }
      }

    Vec<3> cy = 0.5*(ymin+ymax);
    double ry = 0;
    for (int j = 0; j < 3; j++)
      ry = max(ry, ymax(j)-ymin(j));

    return std::make_tuple(cx, rx, cy, ry);
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
      std::tie(cx, rx, cy, ry) = GetMidAndRadius(xpts, ypts);
    }
      
    
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
    using BASE::xpts, BASE::ypts, BASE::xnv, BASE::ynv;

  public:
    FMM_Operator(KERNEL _kernel, Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                 Array<Vec<3>> _xnv, Array<Vec<3>> _ynv)
      : BASE(std::move(_xpts), std::move( _ypts), std::move(_xnv), std::move(_ynv), KERNEL::Shape()),
      kernel(_kernel)
    { }

    void Mult(const BaseVector & x, BaseVector & y) const override
    {
      static Timer tall("ngbem fmm apply"+KERNEL::Name()); RegionTimer reg(tall);
      
      auto fx = x.FV<typename KERNEL::value_type>();
      auto fy = y.FV<typename KERNEL::value_type>();
      
      fy = 0;
      // cout << "slow FMM, kernel = " << typeid(KERNEL).name() << endl;
      if constexpr (std::is_same<KERNEL, class LaplaceSLKernel<3>>())
        {
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
              fy(iy) += __Kernel(xpts[ix], ypts[iy]) * fx(ix);
        }

      else if constexpr (std::is_same<KERNEL, class LaplaceDLKernel<3>>())
        {
          // This is slow, but works for now:
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
            {
              double norm = L2Norm(xpts[ix]-ypts[iy]);
              if (norm > 0)
              {
                // double nxy = InnerProduct(ynv[iy], xpts[ix]-ypts[iy]);
                double nxy = InnerProduct(xnv[ix], -xpts[ix]+ypts[iy]);
                auto kern = nxy / (4 * M_PI * norm*norm*norm);
                fy(iy) += kern * fx(ix);
              }
            }
        }

      else if constexpr (std::is_same<KERNEL, class LaplaceHSKernel<3>>())
        {
          if (fx.Size() != 3*xpts.Size() || fy.Size() != 3*ypts.Size())
            throw Exception ("hs, wrong size");
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
              // fy(iy) += __Kernel(xpts[ix], ypts[iy]) * fx(ix);
              fy.Range(3*iy,3*iy+3) += __Kernel(xpts[ix], ypts[iy]) * fx.Range(3*ix,3*ix+3);
        }
      
      else if constexpr (std::is_same<KERNEL, class HelmholtzSLKernel<3>>())
        {
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
              {
                double norm = L2Norm(xpts[ix]-ypts[iy]);
                if (norm > 0)
                  {
                    auto kern = exp(Complex(0,kernel.GetKappa())*norm) / (4 * M_PI * norm);
                    fy(iy) += kern * fx(ix);
                  }
              }
        }
      else if constexpr (std::is_same<KERNEL, class HelmholtzHSKernel<3>>())
       {
          if (fx.Size() != 4*xpts.Size() || fy.Size() != 4*ypts.Size())
            throw Exception ("HelmholtzHSKernel, wrong size");

          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
            {
              double norm = L2Norm(xpts[ix]-ypts[iy]);
              if (norm > 0)
                {
                  auto kappa = kernel.GetKappa();
                  auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
                  fy.Range(4*iy,4*iy+3) += kern * fx.Range(4*ix,4*ix+3);
                  fy(4*iy+3) += -kappa * kappa * InnerProduct(xnv[ix], ynv[iy])* kern * fx(4*ix+3);
                }
            }
        }
      else if constexpr (std::is_same<KERNEL, class CombinedFieldKernel<3>>())
        {
          /*
            T norm = L2Norm(x-y);
            T nxy = InnerProduct(ny, (x-y));
            auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
            * ( nxy * (Complex(1,0)*T(1.) - Complex(0,kappa)*norm)  - Complex(0,kappa)*norm*norm);
            // return kern;
            return Vec<1,decltype(kern)> (kern);
          */
          double kappa = kernel.GetKappa();
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
              {
                double norm = L2Norm(xpts[ix]-ypts[iy]);
                if (norm > 0)
                  {
                    // auto kern = exp(Complex(0,kernel.GetKappa())*norm) / (4 * M_PI * norm);
                    // fy(iy) += kern * fx(ix);
                    double nxy = InnerProduct(ynv[iy], xpts[ix]-ypts[iy]);
                    auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
                      * ( nxy * (1.0 - Complex(0,kappa)*norm) - Complex(0,kappa)*norm*norm);
                    fy(iy) += kern * fx(ix);                    
                  }
              }
        }
      else if constexpr (std::is_same<KERNEL, class MaxwellSLKernel<3>>())
        {
          if (fx.Size() != 4*xpts.Size() || fy.Size() != 4*ypts.Size())
            throw Exception ("MaxwellSL, wrong size");

          auto kappa = kernel.GetKappa();
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
            {
              double norm = L2Norm(xpts[ix]-ypts[iy]);
              if (norm > 0)
              {
                auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
                fy.Range(4*iy,4*iy+3) += kappa * kern * fx.Range(4*ix,4*ix+3);
                fy(4*iy+3) -= 1.0/kappa * kern * fx(4*ix+3);
              }
            }
        }

      else if constexpr (std::is_same<KERNEL, class MaxwellDLKernel<3>>())
        {
          if (fx.Size() != 3*xpts.Size() || fy.Size() != 3*ypts.Size())
            throw Exception ("MaxwellDL, wrong size");

          auto kappa = kernel.GetKappa();
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
            {
              double norm = L2Norm(xpts[ix]-ypts[iy]);
              if (norm > 0)
              {
                auto G = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm * norm * norm)
                  * (Complex(1,0) - Complex(0,kappa)*norm);
                Vec<3, Complex> kern = G * (xpts[ix]-ypts[iy]);
                fy.Range(3*iy,3*iy+3) += Cross(kern, Vec<3, Complex> (fx.Range(3*ix,3*ix+3)));
              }
            }
        }
      else
        throw Exception(string("fmm not available")+typeid(KERNEL).name());
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
  void FMM_Operator<LaplaceSLKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply LaplaceSL (ngfmm)"); RegionTimer reg(tall);
    
    auto fx = x.FV<double>();
    auto fy = y.FV<double>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    double kappa = 1e-16;

    auto singmp = make_shared<SingularMLMultiPole<Complex>>(cx, rx, kappa);

    for (int i = 0; i < xpts.Size(); i++)
      singmp->AddCharge(xpts[i], fx(i));


    static Timer tsing("call sing");
    {
      RegionTimer regsing(tsing);    
      singmp->CalcMP();
    }
    

    static Timer taddtar("call reg add target");    
    
    RegularMLMultiPole<Complex> regmp (cy, ry, kappa);
    taddtar.Start();
    for (int i = 0; i < ypts.Size(); i++)
      regmp.AddTarget(ypts[i]);
    taddtar.Stop();
    
    
    static Timer treg("call reg");
    {
      RegionTimer regreg(treg);    
      regmp.CalcMP(singmp);
    }

    static Timer teval("ngbem fmm apply LaplaceSL eval"); RegionTimer regeval(teval);    
    ParallelFor (ypts.Size(), [&](int i) {
      fy(i) = Real(regmp.Evaluate(ypts[i]));
    });

    /*
    Vector<double> fyslow(fy.Size());
    fyslow = 0.0;
    for (size_t ix = 0; ix < xpts.Size(); ix++)
      for (size_t iy = 0; iy < ypts.Size(); iy++)
        fyslow(iy) += __Kernel(xpts[ix], ypts[iy]) * fx(ix);

    cout << "norm(fy) = " << L2Norm(fy) << endl;
    cout << "SL diff to slow Laplace: " << L2Norm(fy-fyslow) << endl;

    Vector<double> fyslowkappa(fy.Size());
    fyslowkappa = 0.0;
    for (size_t ix = 0; ix < xpts.Size(); ix++)
      for (size_t iy = 0; iy < ypts.Size(); iy++)
        {
          double norm = L2Norm(xpts[ix]-ypts[iy]);
          if (norm > 0)
            {
              auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
              fyslowkappa(iy) += Real(kern * fx(ix));
            }
        }
    cout << "SL diff to slow Helmholtz: " << L2Norm(fy-fyslowkappa) << endl;
    cout << "slow Helmholtz diff to slow Laplace: " << L2Norm(fyslow-fyslowkappa) << endl;
    */
  }

  template <>
  void FMM_Operator<LaplaceDLKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply LaplaceDL (ngfmm)"); RegionTimer reg(tall);
    auto fx = x.FV<double>();
    auto fy = y.FV<double>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    double kappa = 1e-16;

    auto singmp = make_shared<SingularMLMultiPole<Complex>>(cx, rx, kappa);

    for (int i = 0; i < xpts.Size(); i++)
      singmp->AddDipole(xpts[i], -xnv[i], fx(i));
    singmp->CalcMP();

    RegularMLMultiPole<Complex> regmp (cy, ry, kappa);
    for (int i = 0; i < ypts.Size(); i++)
      regmp.AddTarget(ypts[i]);
    regmp.CalcMP(singmp);

    ParallelFor (ypts.Size(), [&](int i) {
      fy(i) = Real(regmp.Evaluate(ypts[i]));
    });

    /*
    Vector<double> fyslow(fy.Size());
    fyslow = 0.0;
    for (size_t ix = 0; ix < xpts.Size(); ix++)
      for (size_t iy = 0; iy < ypts.Size(); iy++)
        {
          // fyslow(iy) += __Kernel(xpts[ix], ypts[iy]) * fx(ix);
          double norm = L2Norm(xpts[ix]-ypts[iy]);
          if (norm > 0)
            {
              double nxy = InnerProduct(xnv[ix], -xpts[ix]+ypts[iy]);
              auto kern = nxy / (4 * M_PI * norm*norm*norm);
              fyslow(iy) += kern * fx(ix);
            }
        }
    
    cout << "norm(fy) = " << L2Norm(fy) << endl;
    cout << "DL diff to slow Laplace: " << L2Norm(fy-fyslow) << endl;

    Vector<double> fyslowkappa(fy.Size());
    fyslowkappa = 0.0;
    for (size_t ix = 0; ix < xpts.Size(); ix++)
      for (size_t iy = 0; iy < ypts.Size(); iy++)
        {
          double norm = L2Norm(xpts[ix]-ypts[iy]);
          if (norm > 0)
            {
              auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
              fyslowkappa(iy) += Real(kern * fx(ix));
            }
        }
    cout << "DL diff to slow Helmholtz: " << L2Norm(fy-fyslowkappa) << endl;
    cout << "slow Helmholtz diff to slow Laplace: " << L2Norm(fyslow-fyslowkappa) << endl;
    */
  }

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

    for (int i = 0; i < ypts.Size(); i++)
      singmp->AddCharge(ypts[i], fx(i));
    singmp->CalcMP();

    RegularMLMultiPole<Complex> regmp (cx, rx, kappa);
    for (int i = 0; i < xpts.Size(); i++)
      regmp.AddTarget(xpts[i]);
    regmp.CalcMP(singmp);

    ParallelFor (xpts.Size(), [&](int i) {
      fy(i) = Real (regmp.EvaluateDirectionalDerivative(xpts[i], xnv[i]));
    });
  }
  

  template <>
  void FMM_Operator<LaplaceHSKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply LaplaceHS (ngfmm)"); RegionTimer reg(tall);
    auto fx = x.FV<double>();
    auto fy = y.FV<double>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    double kappa = 1e-16;

    auto singmp = make_shared<SingularMLMultiPole<Vec<3,Complex>>>(cx, rx, kappa);

    for (int i = 0; i < xpts.Size(); i++)
      singmp->AddCharge(xpts[i], Vec<3> (fx.Range(3*i,3*i+3)));

    singmp->CalcMP();

    RegularMLMultiPole<Vec<3,Complex>> regmp (cy, ry, kappa);
    for (int i = 0; i < ypts.Size(); i++)
      regmp.AddTarget(ypts[i]);
    regmp.CalcMP(singmp);

    ParallelFor (ypts.Size(), [&](int i) {
      fy.Range(3*i,3*i+3) = Real(regmp.Evaluate(ypts[i]));
    });
  }

  // ********************** operators for Helmholtz ********************  
  

  template <>
  void FMM_Operator<HelmholtzSLKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply HelmholtzSL (ngfmm)"); RegionTimer reg(tall);
    auto fx = x.FV<Complex>();
    auto fy = y.FV<Complex>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    double kappa = kernel.GetKappa();

    auto singmp = make_shared<SingularMLMultiPole<Complex>>(cx, rx, kappa);

    for (int i = 0; i < xpts.Size(); i++)
      singmp->AddCharge(xpts[i], fx(i));

    singmp->CalcMP();

    /*
    RegularMLMultiPole regmp (singmp, cy, ry, int(3*kappa*ry));
    */
    RegularMLMultiPole<Complex> regmp (cy, ry, kappa);
    for (int i = 0; i < ypts.Size(); i++)
      regmp.AddTarget(ypts[i]);
    regmp.CalcMP(singmp);

    
    // for (int i = 0; i < ypts.Size(); i++)
    ParallelFor (ypts.Size(), [&](int i) {
      fy(i) = regmp.Evaluate(ypts[i]);
    });
  }

  template <>
  void FMM_Operator<HelmholtzDLKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply HelmholtzDL (ngfmm)"); RegionTimer reg(tall);
    auto fx = x.FV<Complex>();
    auto fy = y.FV<Complex>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    double kappa = kernel.GetKappa();
    auto singmp = make_shared<SingularMLMultiPole<Complex>>(cx, rx, kappa);

    for (int i = 0; i < xpts.Size(); i++)
      singmp->AddDipole(xpts[i], -xnv[i], fx(i));
    singmp->CalcMP();

    RegularMLMultiPole<Complex> regmp (cy, ry, kappa);
    for (int i = 0; i < ypts.Size(); i++)
      regmp.AddTarget(ypts[i]);
    regmp.CalcMP(singmp);

    ParallelFor (ypts.Size(), [&](int i) {
      fy(i) = regmp.Evaluate(ypts[i]);
    });
  }

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

    for (int i = 0; i < ypts.Size(); i++)
      singmp->AddCharge(ypts[i], fx(i));
    singmp->CalcMP();

    RegularMLMultiPole<Complex> regmp (cx, rx, kappa);
    for (int i = 0; i < xpts.Size(); i++)
      regmp.AddTarget(xpts[i]);
    regmp.CalcMP(singmp);

    ParallelFor (xpts.Size(), [&](int i) {
      fy(i) = regmp.EvaluateDirectionalDerivative(xpts[i], xnv[i]);
    });
  }

  template <>
  void FMM_Operator<HelmholtzHSKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply HelmholtzHS (ngfmm)"); RegionTimer reg(tall);
    auto fx = x.FV<Complex>();
    auto fy = y.FV<Complex>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    double kappa = kernel.GetKappa();

    auto singmp = make_shared<SingularMLMultiPole<Vec<6,Complex>>>(cx, rx, kappa);

    for (int i = 0; i < xpts.Size(); i++){
      Vec<6,Complex> charge;
      charge.Range(0,3) = fx.Range(4*i,4*i+3);
      charge.Range(3,6) = -kappa * kappa * fx(4*i+3) * xnv[i];
      singmp->AddCharge(xpts[i], charge);
    }

    singmp->CalcMP();

    RegularMLMultiPole<Vec<6,Complex>> regmp (cy, ry, kappa);
    for (int i = 0; i < ypts.Size(); i++){
      regmp.AddTarget(ypts[i]);
    }
    regmp.CalcMP(singmp);

    ParallelFor (ypts.Size(), [&](int i) {
      Vec<6,Complex> eval = regmp.Evaluate(ypts[i]);
      fy.Range(4*i,4*i+3) = eval.Range(0,3);
      fy(4*i+3) = InnerProduct(eval.Range(3,6), ynv[i]);
    });
  }


  template <>
  void FMM_Operator<CombinedFieldKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply CombinedField (ngfmm)"); RegionTimer reg(tall);
    static Timer t1("ngbem fmm apply CombinedField (ngfmm) find center/rad + fill");
    static Timer t2("ngbem fmm apply CombinedField (ngfmm) sing");
    static Timer t3("ngbem fmm apply CombinedField (ngfmm) reg");
    static Timer t4("ngbem fmm apply CombinedField (ngfmm) eval");        
    auto fx = x.FV<Complex>();
    auto fy = y.FV<Complex>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    double kappa = kernel.GetKappa();
    
    t1.Start();
    auto singmp = make_shared<SingularMLMultiPole<Complex>>(cx, rx, kappa);
    for (int i = 0; i < xpts.Size(); i++)
      {
        singmp->AddCharge(xpts[i], Complex(0,-kappa)*fx(i));
        singmp->AddDipole(xpts[i], -xnv[i], fx(i));        
      }

    t1.Stop();
    t2.Start();
    
    singmp->CalcMP();
    // cout << "sing norm = " << singmp->Norm() << ", coefs = " << singmp->NumCoefficients() << endl;
    t2.Stop();
    t3.Start();
    // RegularMLMultiPole regmp (singmp, cy, ry, int(2*kappa*ry));

    RegularMLMultiPole<Complex> regmp (cy, ry, kappa);
    for (int i = 0; i < ypts.Size(); i++)
      regmp.AddTarget(ypts[i]);

    regmp.CalcMP(singmp);
    // cout << "reg norm = " << regmp.Norm() << ", coefs = " << regmp.NumCoefficients() << endl;
    t3.Stop();
    t4.Start();
    // for (int i = 0; i < ypts.Size(); i++)
    // fy(i) = regmp.Evaluate(ypts[i]);

    ParallelFor (ypts.Size(), [&](int i) {
      fy(i) = regmp.Evaluate(ypts[i]);
    });
    t4.Stop();
  }

  template <>
  void FMM_Operator<MaxwellSLKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply MaxwellSL (ngfmm)"); RegionTimer reg(tall);
    auto fx = x.FV<Complex>();
    auto fy = y.FV<Complex>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    double kappa = kernel.GetKappa();

    auto singmp = make_shared<SingularMLMultiPole<Vec<4,Complex>>>(cx, rx, kappa);

    for (int i = 0; i < xpts.Size(); i++){
      Vec<4,Complex> charge;
      charge.Range(0,3) = kappa * fx.Range(4*i,4*i+3);
      charge[3] = -1.0 / kappa * fx(4*i+3);
      singmp->AddCharge(xpts[i], charge);
    }
    singmp->CalcMP();

    RegularMLMultiPole<Vec<4,Complex>> regmp (cy, ry, kappa);
    for (int i = 0; i < ypts.Size(); i++){
      regmp.AddTarget(ypts[i]);
    }
    regmp.CalcMP(singmp);

    ParallelFor (ypts.Size(), [&](int i) {
      fy.Range(4*i, 4*i+4) = regmp.Evaluate(ypts[i]);
    });
  }

  template <>
  void FMM_Operator<MaxwellDLKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const
  {
    static Timer tall("ngbem fmm apply MaxwellDL (ngfmm)"); RegionTimer reg(tall);
    auto fx = x.FV<Complex>();
    auto fy = y.FV<Complex>();

    fy = 0;
    if (L2Norm(x) == 0) return;
    double kappa = kernel.GetKappa();

    auto singmp = make_shared<SingularMLMultiPole<Vec<3,Complex>>>(cx, rx, kappa);

    for (int i = 0; i < xpts.Size(); i++){
      Vec<3,Complex> n_cross_m = fx.Range(3*i, 3*i+3);

      // Add dipoles similar to AddCurrentDensity
      for (int k = 0; k < 3; k++)
      {
        Vec<3> ek{0.0}; ek(k) = 1;
        Vec<3> n_cross_m_real = Real(n_cross_m);
        Vec<3> n_cross_m_imag = Imag(n_cross_m);

        singmp->AddDipole(xpts[i], Cross(n_cross_m_real, ek), ek);
        singmp->AddDipole(xpts[i], Cross(n_cross_m_imag, ek), Complex(0,1)*ek);
      }
    }
    singmp->CalcMP();

    RegularMLMultiPole<Vec<3,Complex>> regmp (cy, ry, kappa);
    for (int i = 0; i < ypts.Size(); i++){
      regmp.AddTarget(ypts[i]);
    }
    regmp.CalcMP(singmp);

    ParallelFor (ypts.Size(), [&](int i) {
      fy.Range(3*i, 3*i+3) = regmp.Evaluate(ypts[i]);
    });
  }
}


#endif

