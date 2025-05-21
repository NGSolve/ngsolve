#ifndef NGSOLVE_CONTACT_HPP
#define NGSOLVE_CONTACT_HPP

// #include <comp.hpp>
#include "gridfunction.hpp"
#include "bilinearform.hpp"

namespace ngcomp
{
  template<int DIM>
  struct ContactPair
  {
    ElementId primary_el, secondary_el;
    IntegrationPoint primary_ip, secondary_ip;
  };

  class GapFunction : public CoefficientFunctionNoDerivative
  {
  protected:
    shared_ptr<GridFunction> displacement;
    shared_ptr<MeshAccess> ma;
    Region master;
    Region other;
    double h;
    bool both_sides;

  public:
    GapFunction( shared_ptr<MeshAccess> ma_, Region primary_, Region secondary_)
      : CoefficientFunctionNoDerivative(ma_->GetDimension()),
        ma(ma_), master(primary_), other(secondary_)
    { }

    virtual void Update(shared_ptr<GridFunction> gf, int intorder_, double h_,
                        bool both_sides) = 0;
    void Draw();
  };

  template <int DIM>
  class T_GapFunction : public GapFunction
  {
    unique_ptr<netgen::BoxTree<DIM, int>> searchtree;
  public:
    T_GapFunction( shared_ptr<MeshAccess> mesh_, Region primary_, Region secondary_)
      : GapFunction(mesh_, primary_, secondary_)
    { }

    void Update(shared_ptr<GridFunction> gf, int intorder_, double h, bool both_sides) override;

    const netgen::BoxTree<DIM, int>& GetSearchTree() { return *searchtree; }

    using GapFunction::Evaluate;
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception("Scalar evaluate of GapFunction called");
    }

    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<> result) const override;

    void Evaluate(const BaseMappedIntegrationRule & mir,
                  BareSliceMatrix<> result) const override;

    optional<ContactPair<DIM>> CreateContactPair(const MappedIntegrationPoint<DIM-1, DIM>& mip, LocalHeap& lh, bool both_sides) const;
  };

  template<int DIM>
  class DisplacedNormal : public CoefficientFunctionNoDerivative
  {
    shared_ptr<GridFunction> displacement;
  public:
    DisplacedNormal()
      : CoefficientFunctionNoDerivative(DIM, false),
        displacement(nullptr) {}

    void Update(shared_ptr<GridFunction> _displacement)
    { displacement = _displacement; }

    using CoefficientFunctionNoDerivative::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint& ip) const override
    {
      throw Exception("1dim eval called for Normal");
    }

    void Evaluate(const BaseMappedIntegrationPoint& ir, FlatVector<> values) const override;
  };

  class ContactEnergy
  {
    shared_ptr<CoefficientFunction> cf;
    shared_ptr<FESpace> fes;
    Array<ProxyFunction*> trial_proxies;
    bool deformed;

  public:
    ContactEnergy(shared_ptr<CoefficientFunction> _cf,
                  bool _deformed=false);

    bool IsDeformed() const { return deformed; }

    double CalcEnergy(const FiniteElement& m_fel,
                      const FiniteElement& s_fel,
                      const BaseMappedIntegrationRule& m_mir,
                      FlatVector<double> elx,
                      LocalHeap& lh);

    void ApplyAdd(const FiniteElement& m_fel,
                  const FiniteElement& s_fel,
                  const BaseMappedIntegrationRule& m_mir,
                  FlatVector<double> elx,
                  FlatVector<double> ely,
                  LocalHeap& lh);

    void CalcLinearizedAdd(const FiniteElement& m_fel,
                           const FiniteElement& s_fel,
                           const BaseMappedIntegrationRule& m_mir,
                           FlatVector<double> elx,
                           FlatMatrix<double> elmat,
                           LocalHeap& lh);
  };

  class ContactIntegrator
  {
    shared_ptr<CoefficientFunction> cf;
    shared_ptr<FESpace> fes;
    Array<ProxyFunction*> trial_proxies, test_proxies;
    bool deformed;
  public:
    ContactIntegrator(shared_ptr<CoefficientFunction> _cf,
                      bool _deformed);

    bool IsDeformed() const { return deformed; }

    void ApplyAdd(const FiniteElement& m_fel,
                  const FiniteElement& s_fel,
                  const BaseMappedIntegrationRule& m_mir,
                  FlatVector<double> elx,
                  FlatVector<double> ely,
                  LocalHeap& lh);

    void CalcLinearizedAdd(const FiniteElement& m_fel,
                           const FiniteElement& s_fel,
                           const BaseMappedIntegrationRule& m_mir,
                           FlatVector<double> elx,
                           FlatMatrix<double> elmat,
                           LocalHeap& lh);
  };


  class ContactBoundary : public std::enable_shared_from_this<ContactBoundary>,
                          netgen::UserVisualizationObject
  {
    shared_ptr<GapFunction> gap;
    shared_ptr<CoefficientFunction> normal;
    Region master, other;
    Array<shared_ptr<ContactEnergy>> energies, undeformed_energies, deformed_energies;
    Array<shared_ptr<ContactIntegrator>> integrators, undeformed_integrators, deformed_integrators;
    shared_ptr<FESpace> fes_displacement;
    shared_ptr<FESpace> fes;

    // For visualization only
    bool draw_pairs = false;
    Array<Vec<3>> primary_points;
    Array<Vec<3>> secondary_points;
    bool volume, element_boundary;
  public:
    void Draw();
    ContactBoundary(Region _master, Region _other, bool draw_pairs = false, bool _volume=false, bool element_boundary=false);

    ~ContactBoundary();

    void AddEnergy(shared_ptr<CoefficientFunction> form,
                   bool deformed=false);
    void AddIntegrator(shared_ptr<CoefficientFunction> form,
                       bool deformed=false);

    // Update search tree for gap function, if bf is not
    // nullptr, update SpecialElements of bf
    void Update(shared_ptr<GridFunction> gf,
                shared_ptr<BilinearForm> bf,
                int intorder, double h, bool both_sides);

    shared_ptr<CoefficientFunction> Gap() const { return gap; }
    shared_ptr<CoefficientFunction> Normal() const { return normal; }
    const auto& GetEnergies() const { return energies; }
    const auto& GetEnergies(bool def) const { return def ? deformed_energies : undeformed_energies; }    
    const auto& GetIntegrators() const { return integrators; }
    const auto& GetIntegrators(bool def) const { return def ? deformed_integrators : undeformed_integrators; }    
    shared_ptr<FESpace> GetFESpace() const { return fes; }
    tuple<FlatArray<Vec<3>>, FlatArray<Vec<3>>> GetDrawingPairs() { return {primary_points, secondary_points}; }
  };

  template<int DIM>
  class MPContactElement : public SpecialElement
  {
    // ContactPair<DIM> pair;
    ElementId primary_ei, secondary_ei;
    IntegrationRule primary_ir, secondary_ir;
    shared_ptr<ContactBoundary> cb;
    FESpace* fes;
    GridFunction* deformation;
  public:
    MPContactElement(ElementId primary_ei, ElementId secondary_ei,
                     IntegrationRule primary_ir, IntegrationRule secondary_ir,
                     shared_ptr<ContactBoundary> _cb,
                     GridFunction* deformation);

    void GetDofNrs(Array<DofId>& dnums) const override;

    double Energy(FlatVector<double> elx,
                  LocalHeap& lh) const override;

    void Apply(FlatVector<double> elx,
               FlatVector<double> ely,
               LocalHeap& lh) const override;

    void CalcElementMatrix(FlatMatrix<double> elmat,
                           LocalHeap& lh) const override;

    
    void CalcLinearizedElementMatrix(FlatVector<double> elx,
                                     FlatMatrix<double> elmat,
                                     LocalHeap& lh) const override;

    shared_ptr<ContactBoundary> GetContactBoundary() const
    { return cb; }
  };

  
} // namespace ngcomp

#endif // NGSOLVE_CONTACT_HPP
