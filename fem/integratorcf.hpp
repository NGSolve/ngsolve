#ifndef FILE_INTEGRATORCF
#define FILE_INTEGRATORCF


// WIP, need a bit more design ...

#include "coefficient.hpp"
#include "integrator.hpp"
#include "symbolicintegrator.hpp"

namespace ngfem
{
  class Integral;
  class DifferentialSymbol
  {
  public:
    VorB vb;
    VorB element_vb = VOL;
    bool skeleton = false;
    optional<variant<BitArray,string>> definedon;
    int bonus_intorder = 0;
    shared_ptr<ngcomp::GridFunction> deformation;
    std::map<ELEMENT_TYPE,shared_ptr<IntegrationRule>> userdefined_intrules;
    shared_ptr<BitArray> definedonelements;
    
    DifferentialSymbol (VorB _vb) : vb(_vb) { ; }
    DifferentialSymbol (VorB _vb, VorB _element_vb, bool _skeleton, // const BitArray & _definedon,
                        int _bonus_intorder)
      : vb(_vb), element_vb(_element_vb), skeleton(_skeleton), /* definedon(_definedon), */ bonus_intorder(_bonus_intorder) { ; }

    virtual ~DifferentialSymbol() { }
    virtual shared_ptr<Integral> MakeIntegral(shared_ptr<CoefficientFunction> cf) const
    {
      return make_shared<Integral> (cf, *this);
    }
  };
  

  class Integral
  {
  public:
    shared_ptr<CoefficientFunction> cf;
    DifferentialSymbol dx;
    shared_ptr<Integral> linearization;
    Integral (shared_ptr<CoefficientFunction> _cf,
              DifferentialSymbol _dx)
      : cf(_cf), dx(_dx) { ; }
    virtual ~Integral() { }

    template <typename TSCAL>
    TSCAL T_Integrate (const ngcomp::MeshAccess & ma,
                       FlatVector<TSCAL> element_wise);

    NGS_DLL_HEADER virtual double Integrate (const ngcomp::MeshAccess & ma,
                              FlatVector<double> element_wise);

    NGS_DLL_HEADER virtual Complex Integrate (const ngcomp::MeshAccess & ma,
                               FlatVector<Complex> element_wise);
    
    NGS_DLL_HEADER virtual shared_ptr<BilinearFormIntegrator> MakeBilinearFormIntegrator() const;
    NGS_DLL_HEADER virtual shared_ptr<LinearFormIntegrator> MakeLinearFormIntegrator() const;

    NGS_DLL_HEADER virtual shared_ptr<Integral> CreateSameIntegralType (shared_ptr<CoefficientFunction> _cf)
    {
      return make_shared<Integral> (_cf, dx);
    }    
  };
  
  inline Integral operator* (double fac, const Integral & cf)
  {
    return Integral (fac * cf.cf, cf.dx);
  }

  inline Integral operator* (Complex fac, const Integral & cf)
  {
    return Integral (fac * cf.cf, cf.dx);
  }

  
  inline ostream & operator<< (ostream & ost, const Integral & igl)
  {
    ost << *igl.cf << " " << igl.dx.vb << endl;
    return ost;
  }
  
  class SumOfIntegrals
  {
  public:
    Array<shared_ptr<Integral>> icfs;
    Array<shared_ptr<Integral>> linearization_icfs;
    
    SumOfIntegrals() = default;
    SumOfIntegrals (shared_ptr<Integral> icf)
    { icfs += icf; }

    auto begin() const { return icfs.begin(); }
    auto end() const { return icfs.end(); }

    shared_ptr<SumOfIntegrals>
    Replace (std::map<shared_ptr<CoefficientFunction>, shared_ptr<CoefficientFunction>> replace)

    {
      auto repl = make_shared<SumOfIntegrals>();
      CoefficientFunction::T_Transform transform;
      transform.replace = replace;
      for (auto & icf : icfs)
        repl->icfs += icf->CreateSameIntegralType (icf->cf->Transform(transform));
      return repl;
    }

    Array<shared_ptr<ProxyFunction>> GetProxies (bool trialproxies)
    {
      Array<shared_ptr<ProxyFunction>> proxies;

      for (auto & icf : icfs)
        icf->cf->TraverseTree
          ( [&] (CoefficientFunction & nodecf)
          {
            auto proxy = dynamic_pointer_cast<ProxyFunction> ((&nodecf)->shared_from_this());
            if (proxy) 
              {
                if (proxy->IsTrialFunction() == trialproxies)
                  {
                    if (!proxies.Contains(proxy))
                      proxies.Append(proxy);
                  }
              }
          });

      return proxies;
    }
    
    shared_ptr<SumOfIntegrals>
    Diff (shared_ptr<CoefficientFunction> var,
          shared_ptr<CoefficientFunction> dir) const
    {
      auto deriv = make_shared<SumOfIntegrals>();
      for (auto & icf : icfs)
        deriv->icfs += icf->CreateSameIntegralType (icf->cf->Diff(var.get(), dir));
      return deriv;
    }

    shared_ptr<SumOfIntegrals>
    DiffShape (shared_ptr<CoefficientFunction> dir) const
    {
      auto deriv = make_shared<SumOfIntegrals>();
      auto grad = dir->Operator("grad");
      if (!grad)
        throw Exception("In SumOfIntegrals::DiffShape: dir does not have a grad operator");
      auto divdir = TraceCF(grad);
      auto sgrad = dynamic_pointer_cast<ProxyFunction>(grad)->Trace();
      auto sdivdir = TraceCF(sgrad);

      auto tang = TangentialVectorCF(dir->Dimension(), false) -> Reshape(dir->Dimension(), 1);
      auto bsdivdir = InnerProduct(sgrad*tang,tang);

      DiffShapeCF shape;
      //cout << "should add Eulerian here" << endl;
      
      for (auto & icf : icfs)
        {
          switch (icf->dx.vb)
            {
            case VOL:
              if (icf->dx.element_vb == VOL)
                deriv->icfs += icf->CreateSameIntegralType ( icf->cf->Diff(&shape, dir) + divdir*icf->cf);
              else
                throw Exception("In DiffShape: for vb=VOL only element_vb=VOL implemented!");
              break;
            case BND:
              if (icf->dx.element_vb == VOL)
                deriv->icfs += icf->CreateSameIntegralType ( icf->cf->Diff(&shape, dir) + sdivdir*icf->cf);
              else if (icf->dx.element_vb == BND && dir->Dimension() == 3)
                deriv->icfs += icf->CreateSameIntegralType ( icf->cf->Diff(&shape, dir) + bsdivdir*icf->cf);
              else if (icf->dx.element_vb == BND && dir->Dimension() == 2)
                deriv->icfs += icf->CreateSameIntegralType ( icf->cf->Diff(&shape, dir));
              else
                throw Exception("In DiffShape: for vb=BND something went wrong!");
              break;
            default:
              throw Exception("In DiffShape: for vb="+ToString(icf->dx.vb)+" and element_vb="+ToString(icf->dx.element_vb) + " not implemented!");
            }
        }
      return deriv;
    }

    
    shared_ptr<SumOfIntegrals>
    Compile (bool realcompile, bool wait, bool keep_files) const
    {
      auto compiled = make_shared<SumOfIntegrals>();
      for (auto & icf : icfs)
        compiled->icfs += icf->CreateSameIntegralType (::ngfem::Compile (icf->cf, realcompile, 2, wait, keep_files));
      return compiled;
    }

    void SetDefinedOnElements(shared_ptr<BitArray> defon)
    {
      for(auto& icf : icfs)
        icf->dx.definedonelements = defon;
    }
  };

  inline auto operator+ (const SumOfIntegrals & c1, const SumOfIntegrals & c2)
  {
    SumOfIntegrals sum;
    for (auto & ci : c1.icfs) sum.icfs += ci;
    for (auto & ci : c2.icfs) sum.icfs += ci;
    return sum;
  }

  inline auto operator* (double fac, SumOfIntegrals c1)
  {
    SumOfIntegrals faccf;
    for (auto & ci : c1.icfs) faccf.icfs += ci->CreateSameIntegralType(fac*(ci->cf));
    return faccf;
  }

  inline auto operator* (Complex fac, SumOfIntegrals c1)
  {
    SumOfIntegrals faccf;
    for (auto & ci : c1.icfs) faccf.icfs += ci->CreateSameIntegralType(fac*(ci->cf));
    return faccf;
  }

  inline auto operator- (const SumOfIntegrals & c1, const SumOfIntegrals & c2)
  {
    return c1 + (-1)*c2;
  }
  
  inline ostream & operator<< (ostream & ost, const SumOfIntegrals & igls)
  {
    for (auto & igl : igls.icfs)
      ost << *igl;
    return ost;
  }
  
  class Variation
  {
  public:
    shared_ptr<SumOfIntegrals> igls;
    Variation (shared_ptr<SumOfIntegrals> _igls) : igls(_igls) { ; }
    
    auto Compile (bool realcompile, bool wait, bool keep_files) const
    {
      return Variation(igls->Compile(realcompile, wait, keep_files));
    }
  };
  
}


#endif
