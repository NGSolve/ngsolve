#ifndef MP_COEFFICIENT
#define MP_COEFFICIENT


#include "mptools.hpp"

namespace ngsbem
{

  
  // ******************** Coefficient Functions *********************

  
  class SphericalHarmonicsCF : public CoefficientFunction
  {
    SphericalHarmonics<Complex> sh;
  public:
    SphericalHarmonicsCF (int order)
      : CoefficientFunction(1, true), sh(order) { }
    Complex & Coef(int n, int m) { return sh.Coef(n,m); } 
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("real eval not available"); }

    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const override
    {
      values(0) = sh.Eval(mip.GetPoint());
    }
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const override
    {
      for (int i = 0; i < ir.Size(); i++)
        {
          auto & mip = ir[i];
          values(i,0) = sh.Eval(mip.GetPoint());
        }
    }

    auto & SH() { return sh; }
  };



  template <typename entry_type> class RegularMLExpansionCF;
  

  template <typename RADIAL, typename entry_type=Complex>
  class SphericalExpansionCF : public CoefficientFunction
  {
    SphericalExpansion<RADIAL, entry_type> mp;
    Vec<3> center;
  public:
    SphericalExpansionCF (int order, double kappa, Vec<3> acenter, double rtyp = 1)
      : CoefficientFunction(sizeof(entry_type)/sizeof(Complex), true), mp(order, kappa, rtyp), center(acenter) { }

    entry_type & Coef(int n, int m) { return mp.Coef(n,m); } 
    auto & SH() { return mp.SH(); }
    auto & MP() { return mp; }
    Vec<3> Center() const { return center; }
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("real eval not available"); }

    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const override
    {
      if constexpr (std::is_same<entry_type, Complex>())
        values(0) = mp.Eval(mip.GetPoint()-center);
      else
        values = mp.Eval(mip.GetPoint()-center);
    }

    template <typename TARGET>
    void ShiftZ (double z, SphericalExpansion<TARGET, entry_type> & target) { mp.ShiftZ(z, target); }

    using CoefficientFunction::Transform;        
    template <typename TARGET>
    void Transform (SphericalExpansionCF<TARGET, entry_type> & target)
    {
      mp.Transform (target.MP(), target.Center()-center);
    }
    template <typename TARGET>
    void TransformAdd (SphericalExpansionCF<TARGET, entry_type> & target)
    {
      mp.TransformAdd (target.MP(), target.Center()-center);
    }
  };

  template <typename entry_type>
  class SingularMLExpansionCF : public CoefficientFunction
  {
    shared_ptr<SingularMLExpansion<entry_type>> mlmp;
  public:
    SingularMLExpansionCF (Vec<3> center, double r, double kappa)
      : CoefficientFunction(sizeof(entry_type)/sizeof(Complex), true), mlmp{make_shared<SingularMLExpansion<entry_type>>(center, r, kappa)} { }
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("real eval not available"); }
    
    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const override
    {
      // values(0) = mlmp->Evaluate(mip.GetPoint());

      if constexpr (std::is_same<entry_type, Complex>())
        values(0) = mlmp->Evaluate(mip.GetPoint());
      else
        values = mlmp->Evaluate(mip.GetPoint());

      
    }
    
    shared_ptr<SingularMLExpansion<entry_type>> MLExpansion() const { return mlmp; }
    shared_ptr<RegularMLExpansionCF<entry_type>> CreateRegularExpansion(Vec<3> center, double r) const;
  };
  

  template <typename entry_type>
  class RegularMLExpansionCF : public CoefficientFunction
  {
    shared_ptr<RegularMLExpansion<entry_type>> mlmp;
  public:
    RegularMLExpansionCF (shared_ptr<SingularMLExpansionCF<entry_type>> asingmp, Vec<3> center, double r)
      : CoefficientFunction(sizeof(entry_type)/sizeof(Complex), true), mlmp{make_shared<RegularMLExpansion<entry_type>>(asingmp->MLExpansion(), center, r, FMM_Parameters())} { } 
    RegularMLExpansionCF (shared_ptr<SingularMLExpansion<entry_type>> asingmp, Vec<3> center, double r)
      : CoefficientFunction(sizeof(entry_type)/sizeof(Complex), true), mlmp{make_shared<RegularMLExpansion<entry_type>>(asingmp, center, r, FMM_Parameters())} { } 
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("real eval not available"); }

    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const override
    {
      // values(0) = mlmp->Evaluate(mip.GetPoint());

      if constexpr (std::is_same<entry_type, Complex>())
        values(0) = mlmp->Evaluate(mip.GetPoint());
      else
        values = mlmp->Evaluate(mip.GetPoint());
    }

    shared_ptr<RegularMLExpansion<entry_type>> MLExpansion() { return mlmp; }
  };

  
}

#endif
