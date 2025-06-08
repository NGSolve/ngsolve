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



  template <typename entry_type> class RegularMLMultiPoleCF;
  

  template <typename RADIAL, typename entry_type=Complex>
  class MultiPoleCF : public CoefficientFunction
  {
    MultiPole<RADIAL, entry_type> mp;
    Vec<3> center;
  public:
    MultiPoleCF (int order, double kappa, Vec<3> acenter, double rtyp = 1)
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
    void ShiftZ (double z, MultiPole<TARGET, entry_type> & target) { mp.ShiftZ(z, target); }

    using CoefficientFunction::Transform;        
    template <typename TARGET>
    void Transform (MultiPoleCF<TARGET, entry_type> & target)
    {
      mp.Transform (target.MP(), target.Center()-center);
    }
  };

  template <typename entry_type>
  class SingularMLMultiPoleCF : public CoefficientFunction
  {
    shared_ptr<SingularMLMultiPole<entry_type>> mlmp;
  public:
    SingularMLMultiPoleCF (Vec<3> center, double r, double kappa)
      : CoefficientFunction(sizeof(entry_type)/sizeof(Complex), true), mlmp{make_shared<SingularMLMultiPole<entry_type>>(center, r, kappa)} { }
    
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
    
    shared_ptr<SingularMLMultiPole<entry_type>> MLMP() const { return mlmp; }
    shared_ptr<RegularMLMultiPoleCF<entry_type>> CreateRegularExpansion(Vec<3> center, double r) const;
  };
  

  template <typename entry_type>
  class RegularMLMultiPoleCF : public CoefficientFunction
  {
    shared_ptr<RegularMLMultiPole<entry_type>> mlmp;
  public:
    RegularMLMultiPoleCF (shared_ptr<SingularMLMultiPoleCF<entry_type>> asingmp, Vec<3> center, double r)
      : CoefficientFunction(sizeof(entry_type)/sizeof(Complex), true), mlmp{make_shared<RegularMLMultiPole<entry_type>>(asingmp->MLMP(), center, r)} { } 
    RegularMLMultiPoleCF (shared_ptr<SingularMLMultiPole<entry_type>> asingmp, Vec<3> center, double r)
      : CoefficientFunction(sizeof(entry_type)/sizeof(Complex), true), mlmp{make_shared<RegularMLMultiPole<entry_type>>(asingmp, center, r)} { } 
    
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

    shared_ptr<RegularMLMultiPole<entry_type>> MLMP() { return mlmp; }
  };

  
}

#endif
