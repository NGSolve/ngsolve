#ifndef FILE_FMMOPERATOR
#define FILE_FMMOPERATOR
#include <mptools.hpp>

namespace ngsbem
{
  
  inline double __Kernel(Vec<3> px, Vec<3> py)
  {
    if (L2Norm2(px-py)==0) return 0.0;
    return 1.0 / (4*M_PI) / L2Norm(px-py);
  }


  template <typename TSCAL>
  class Base_FMM_Operator : public BaseMatrix
  {
  protected:
    Array<Vec<3>> xpts, ypts, xnv, ynv;
  public:
    Base_FMM_Operator(Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                   Array<Vec<3>> _xnv, Array<Vec<3>> _ynv)
      : xpts(std::move(_xpts)), ypts(std::move(_ypts)),
        xnv(std::move(_xnv)), ynv(std::move(_ynv))
    { }
    
    AutoVector CreateRowVector () const override
    {
      return make_unique<VVector<TSCAL>>(xpts.Size());
    }
    AutoVector CreateColVector () const override
    {
      return make_unique<VVector<TSCAL>>(ypts.Size());
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
      : BASE(std::move(_xpts), std::move( _ypts), std::move(_xnv), std::move(_ynv)),
      kernel(_kernel)
    { }

    void Mult(const BaseVector & x, BaseVector & y) const override
    {
      static Timer tall("ngbem fmm apply"); RegionTimer reg(tall);
      
      auto fx = x.FV<typename KERNEL::value_type>();
      auto fy = y.FV<typename KERNEL::value_type>();
      
      fy = 0;
      if constexpr (std::is_same<KERNEL, class LaplaceSLKernel<3>>())
        {
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
              fy(iy) += __Kernel(xpts[ix], ypts[iy]) * fx(ix);
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
      else
        throw Exception("fmm not available");
    }

  };

#define USE_NGFMM
#ifdef USE_NGFMM

  template <>
  void FMM_Operator<HelmholtzSLKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply CombinedField (ngfmm)"); RegionTimer reg(tall);
    auto fx = x.FV<Complex>();
    auto fy = y.FV<Complex>();

    fy = 0;
    
    double kappa = kernel.GetKappa();

    
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


    
    auto singmp = make_shared<SingularMLMultiPole<Complex>>(cx, rx, int(3*kappa*rx), kappa);

    for (int i = 0; i < xpts.Size(); i++)
      singmp->AddCharge(xpts[i], fx(i));

    singmp->CalcMP();

    /*
    RegularMLMultiPole regmp (singmp, cy, ry, int(3*kappa*ry));
    */
    RegularMLMultiPole regmp (cy, ry, int(3*kappa*ry), kappa);
    for (int i = 0; i < ypts.Size(); i++)
      regmp.AddTarget(ypts[i]);
    regmp.CalcMP(singmp);

    
    // for (int i = 0; i < ypts.Size(); i++)
    ParallelFor (ypts.Size(), [&](int i) {
      fy(i) = regmp.Evaluate(ypts[i]);
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
    
    t1.Start();
    
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

    double kappa = kernel.GetKappa();
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

    auto singmp = make_shared<SingularMLMultiPole<Complex>>(cx, rx, int(2*kappa*rx), kappa);
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

    RegularMLMultiPole regmp (cy, ry, int(3*kappa*ry), kappa);
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
  
}


#endif

