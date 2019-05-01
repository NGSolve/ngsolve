// WIP, need a bit more design ...

namespace ngfem
{
  class DifferentialSymbol
  {
  public:
    VorB vb;
    DifferentialSymbol (VorB _vb) : vb(_vb) { ; }
  };

  class IntegratorCF
  {
  public:
    shared_ptr<CoefficientFunction> cf;
    DifferentialSymbol dx;
    IntegratorCF (shared_ptr<CoefficientFunction> _cf,
                  DifferentialSymbol _dx)
      : cf(_cf), dx(_dx) { ; } 
  };

  class SumOfIntegratorCF
  {
  public:
    Array<shared_ptr<IntegratorCF>> icfs;
    
    SumOfIntegratorCF() = default;
    SumOfIntegratorCF (shared_ptr<IntegratorCF> icf)
    { icfs += icf; }

    shared_ptr<SumOfIntegratorCF>
    Derive (shared_ptr<CoefficientFunction> var,
            shared_ptr<CoefficientFunction> dir) const
    {
      auto deriv = make_shared<SumOfIntegratorCF>();
      for (auto & icf : icfs)
        deriv->icfs += make_shared<IntegratorCF> (icf->cf->Derive(var.get(), dir), icf->dx);
      return deriv;
    }
  };
}
