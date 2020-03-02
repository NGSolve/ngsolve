#ifndef NGSOLVE_CONTACT_HPP
#define NGSOLVE_CONTACT_HPP

#include <comp.hpp>

namespace ngcomp
{
  template<int DIM>
  struct ContactPair
  {
    ElementId master_el, slave_el;
    IntegrationPoint master_ip, slave_ip;
  };

  class GapFunction : public CoefficientFunction
  {
  protected:
    shared_ptr<GridFunction> displacement;
    shared_ptr<MeshAccess> ma;
    Region master;
    Region slave;
    double h;

  public:
    GapFunction( shared_ptr<MeshAccess> ma_, Region master_, Region slave_)
      : CoefficientFunction(ma_->GetDimension()),
        ma(ma_), master(master_), slave(slave_)
    { }

    virtual void Update(shared_ptr<GridFunction> gf, int intorder_, double h_) = 0;
  };

  template <int DIM>
  class T_GapFunction : public GapFunction
  {
    unique_ptr<netgen::BoxTree<DIM, int>> searchtree;
  public:
    T_GapFunction( shared_ptr<MeshAccess> mesh_, Region master_, Region slave_)
      : GapFunction(mesh_, master_, slave_)
    { }

    void Update(shared_ptr<GridFunction> gf, int intorder_, double h) override;

    const netgen::BoxTree<DIM, int>& GetSearchTree() { return *searchtree; }

    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception("Scalar evaluate of GapFunction called");
    }

    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<> result) const override;

    void Evaluate(const BaseMappedIntegrationRule & mir,
                  BareSliceMatrix<> result) const override;

    optional<ContactPair<DIM>> CreateContactPair(const MappedIntegrationPoint<DIM-1, DIM>& mip) const;
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
    CoefficientFunction * gap_function;

  public:
    ContactEnergy(shared_ptr<CoefficientFunction> _cf,
                  shared_ptr<FESpace> _fes);

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

  class ContactBoundary
  {
    shared_ptr<GapFunction> gap;
    shared_ptr<CoefficientFunction> normal;
    Region master, slave;
    Array<shared_ptr<ContactEnergy>> energies;
    shared_ptr<FESpace> fes;
  public:
    ContactBoundary(shared_ptr<FESpace> _fes, Region _master,
                    Region _slave);

    void AddEnergy(shared_ptr<CoefficientFunction> form);
    // void AddIntegrator(shared_ptr<CoefficientFunction> form);

    // Update search tree for gap function, if bf is not
    // nullptr, update SpecialElements of bf
    void Update(shared_ptr<GridFunction> gf,
                shared_ptr<BilinearForm> bf,
                int intorder, double h);

    shared_ptr<CoefficientFunction> Gap() const { return gap; }
    shared_ptr<CoefficientFunction> Normal() const { return normal; }
    const auto& GetEnergies() const { return energies; }
    shared_ptr<FESpace> GetFESpace() const { return fes; }
  };

  template<int DIM>
  class ContactElement : public SpecialElement
  {
    ContactPair<DIM> pair;
    ContactBoundary* cb;
    FESpace* fes;
  public:
    ContactElement(const ContactPair<DIM>& _pair,
                   ContactBoundary* _cb);

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
