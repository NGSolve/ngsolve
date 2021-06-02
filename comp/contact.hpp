#ifndef NGSOLVE_CONTACT_HPP
#define NGSOLVE_CONTACT_HPP

#include <comp.hpp>

namespace ngcomp
{
  template<int DIM>
  struct ContactPair
  {
    ElementId master_el, other_el;
    IntegrationPoint master_ip, other_ip;
  };

  class GapFunction : public CoefficientFunctionNoDerivative
  {
  protected:
    shared_ptr<GridFunction> displacement;
    shared_ptr<MeshAccess> ma;
    Region master;
    Region other;
    double h;

  public:
    GapFunction( shared_ptr<MeshAccess> ma_, Region master_, Region other_)
      : CoefficientFunctionNoDerivative(ma_->GetDimension()),
        ma(ma_), master(master_), other(other_)
    { }

    virtual void Update(shared_ptr<GridFunction> gf, int intorder_, double h_) = 0;
    void Draw();
  };

  template <int DIM>
  class T_GapFunction : public GapFunction
  {
    unique_ptr<netgen::BoxTree<DIM, int>> searchtree;
  public:
    T_GapFunction( shared_ptr<MeshAccess> mesh_, Region master_, Region other_)
      : GapFunction(mesh_, master_, other_)
    { }

    void Update(shared_ptr<GridFunction> gf, int intorder_, double h) override;

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

    optional<ContactPair<DIM>> CreateContactPair(const MappedIntegrationPoint<DIM-1, DIM>& mip, LocalHeap& lh) const;
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


  class ContactBoundary : public netgen::UserVisualizationObject
  {
    shared_ptr<GapFunction> gap;
    shared_ptr<CoefficientFunction> normal;
    Region master, other;
    Array<shared_ptr<ContactEnergy>> energies;
    Array<shared_ptr<ContactIntegrator>> integrators;
    shared_ptr<FESpace> fes_displacement;
    shared_ptr<FESpace> fes;

    // For visualization only
    bool draw_pairs = false;
    Array<Vec<3>> master_points;
    Array<Vec<3>> other_points;
  public:
    void Draw();
    ContactBoundary(Region _master, Region _other, bool draw_pairs = false);

    ~ContactBoundary();

    void AddEnergy(shared_ptr<CoefficientFunction> form,
                   bool deformed=false);
    void AddIntegrator(shared_ptr<CoefficientFunction> form,
                       bool deformed=false);

    // Update search tree for gap function, if bf is not
    // nullptr, update SpecialElements of bf
    void Update(shared_ptr<GridFunction> gf,
                shared_ptr<BilinearForm> bf,
                int intorder, double h);

    shared_ptr<CoefficientFunction> Gap() const { return gap; }
    shared_ptr<CoefficientFunction> Normal() const { return normal; }
    const auto& GetEnergies() const { return energies; }
    const auto& GetIntegrators() const { return integrators; }
    shared_ptr<FESpace> GetFESpace() const { return fes; }
  };

  template<int DIM>
  class ContactElement : public SpecialElement
  {
    ContactPair<DIM> pair;
    ContactBoundary* cb;
    FESpace* fes;
    GridFunction* deformation;
  public:
    ContactElement(const ContactPair<DIM>& _pair,
                   ContactBoundary* _cb,
                   GridFunction* deformation);

    void GetDofNrs(Array<DofId>& dnums) const override;

    double Energy(FlatVector<double> elx,
                  LocalHeap& lh) const override;

    void Apply(FlatVector<double> elx,
               FlatVector<double> ely,
               LocalHeap& lh) const override;

    void CalcLinearizedElementMatrix(FlatVector<double> elx,
                                     FlatMatrix<double> elmat,
                                     LocalHeap& lh) const override;

    ContactBoundary* GetContactBoundary() const
    { return cb; }
  };
} // namespace ngcomp

#endif // NGSOLVE_CONTACT_HPP
