// WIP, need a bit more design ...

namespace ngfem
{
  
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
    
    DifferentialSymbol (VorB _vb) : vb(_vb) { ; }
    DifferentialSymbol (VorB _vb, VorB _element_vb, bool _skeleton, // const BitArray & _definedon,
                        int _bonus_intorder)
      : vb(_vb), element_vb(_element_vb), skeleton(_skeleton), /* definedon(_definedon), */ bonus_intorder(_bonus_intorder) { ; } 
  };
  

  class Integral
  {
  public:
    shared_ptr<CoefficientFunction> cf;
    DifferentialSymbol dx;
    Integral (shared_ptr<CoefficientFunction> _cf,
              DifferentialSymbol _dx)
      : cf(_cf), dx(_dx) { ; }

    template <typename TSCAL>
    TSCAL Integrate (const ngcomp::MeshAccess & ma,
                     FlatVector<TSCAL> element_wise);
  };
  
  inline Integral operator* (double fac, const Integral & cf)
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
    
    SumOfIntegrals() = default;
    SumOfIntegrals (shared_ptr<Integral> icf)
    { icfs += icf; }

    shared_ptr<SumOfIntegrals>
    Diff (shared_ptr<CoefficientFunction> var,
          shared_ptr<CoefficientFunction> dir) const
    {
      auto deriv = make_shared<SumOfIntegrals>();
      for (auto & icf : icfs)
        deriv->icfs += make_shared<Integral> (icf->cf->Diff(var.get(), dir), icf->dx);
      return deriv;
    }

    shared_ptr<SumOfIntegrals>
    DiffShape (shared_ptr<CoefficientFunction> dir) const
    {
      auto deriv = make_shared<SumOfIntegrals>();
      auto grad = dir->Operator("grad");
      auto divdir = TraceCF(grad);
      auto sgrad = dynamic_pointer_cast<ProxyFunction>(grad)->Trace();
      auto sdivdir = TraceCF(sgrad);
      
      for (auto & icf : icfs)
        {
          if (icf->dx.vb == VOL)
            deriv->icfs += make_shared<Integral> ( icf->cf->Diff(shape.get(), dir) + divdir*icf->cf, icf->dx);
          else
            {
              auto n = NormalVectorCF(sgrad->Dimensions()[0]);

              deriv->icfs += make_shared<Integral> ( icf->cf->Diff(shape.get(), dir) +
                                                     (sdivdir  /* - InnerProduct(sgrad*n, n) */ ) * icf->cf, icf->dx);
            }
        }
      return deriv;
    }

    
    shared_ptr<SumOfIntegrals>
    Compile (bool realcompile, bool wait) const
    {
      auto compiled = make_shared<SumOfIntegrals>();
      for (auto & icf : icfs)
        compiled->icfs += make_shared<Integral> (::ngfem::Compile (icf->cf, realcompile, 2, wait), icf->dx);
      return compiled;
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
    for (auto & ci : c1.icfs) faccf.icfs += make_shared<Integral>(fac*(*ci));
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
    
    auto Compile (bool realcompile, bool wait) const
    {
      return Variation(igls->Compile(realcompile, wait));
    }
  };
  
}
