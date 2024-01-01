// #include <comp.hpp>
#include "globalinterfacespace.hpp"
#include <recursive_pol.hpp>

namespace ngcomp
{
  template<typename VOLFE>
  void GlobalInterfaceSpace::VolDiffOp<VOLFE>::CalcMatrix
    (const FiniteElement & bfel,
     const BaseMappedIntegrationPoint & mip,
     BareSliceMatrix<double,ColMajor> mat,
     LocalHeap & lh) const
  {
    auto & ip = mip.IP();
    auto & fel = dynamic_cast<const VOLFE&> (bfel);
    mat.AddSize(Dim(), bfel.GetNDof()) = 0;

    int fnr = ip.FacetNr();

    if (fnr != -1 && fel.GetFacet(fnr))
      {
        fel.GetFacet(fnr)->CalcShape(mip, mat.Row(0));
        return;
      }
    if(fnr == -1)
      {
        fel.CalcShape(mip, mat.Row(0));
      }
  }

  template<typename VOLFE>
  void GlobalInterfaceSpace::ParameterGradDiffOp<VOLFE>::CalcMatrix
    (const FiniteElement & bfel,
     const BaseMappedIntegrationPoint & mip,
     BareSliceMatrix<double,ColMajor> mat,
     LocalHeap & lh) const
  {
    auto & ip = mip.IP();
    auto & fel = dynamic_cast<const VOLFE&> (bfel);
    mat.AddSize(Dim(), bfel.GetNDof()) = 0;

    int fnr = ip.FacetNr();

    if (fnr != -1 && fel.GetFacet(fnr))
      {
        throw Exception("Grad diffop not yet implemented for boundary integrals!");
        return;
      }
    if(fnr == -1)
      {
        fel.CalcDShape(mip, mat.Row(0));
      }
  }

  template<typename INTERFACEFE>
  void GlobalInterfaceSpace::InterfaceDiffOp<INTERFACEFE>::CalcMatrix
    (const FiniteElement & fel,
     const BaseMappedIntegrationPoint & mip,
     BareSliceMatrix<double,ColMajor> mat,
     LocalHeap & lh) const
  {
    if (fel.GetNDof() > 0)
      dynamic_cast<const INTERFACEFE&> (fel).CalcShape(mip, mat.Row(0));
  }

  GlobalInterfaceSpace::GlobalInterfaceSpace(shared_ptr<MeshAccess> ama,
                                             const Flags& flags)
    : FESpace(ama, flags)
  {
    order = int(flags.GetNumFlag("order", 3));
    periodic[0] = periodic[1] = false;
    polar = flags.GetDefineFlag("polar");
    if(flags.GetDefineFlag("periodic"))
      periodic[0] = periodic[1] = true;
    if(flags.GetDefineFlag("periodicu"))
      periodic[0] = true;
    if(flags.GetDefineFlag("periodicv"))
      periodic[1] = true;
    try
      {
        mapping = std::any_cast<shared_ptr<CoefficientFunction>>(flags.GetAnyFlag("mapping"));
      }
    catch(const std::bad_any_cast& ex)
      {
        throw Exception("No mapping or wrong mapping given!\nGlobalInterfacespace needs kwarg: mapping=CoefficientFunction");
      }
  }

  template<int DIM>
  class GlobalInterfaceSpaceD : public GlobalInterfaceSpace
  {
    Array<bool> nitsche_facet;

    class InterfaceFE : public FiniteElement
    {
    protected:
      const GlobalInterfaceSpaceD<DIM> * fes;
      ELEMENT_TYPE et;
    public:
      InterfaceFE (const GlobalInterfaceSpaceD<DIM> * afes,
                   ELEMENT_TYPE _et)
        : FiniteElement(afes->GetNDof(), afes->order), fes(afes),
          et(_et)
      { ; }

      ELEMENT_TYPE ElementType() const { return et; }
      void CalcShape(const BaseMappedIntegrationPoint& mip,
                     BareSliceVector<double> shapes) const
      {
        int order = fes->GetOrder();
        if constexpr(DIM ==1)
          {
            if(fes->polar)
              throw Exception("Polar coordinates need 2 dimensional mapping!");
            double phi = fes->mapping -> Evaluate(mip);
            if(fes->periodic[0])
              {
                shapes(0) = 1;
                for (int i = 1; i <= order; i++)
                  {
                    shapes(2*i-1) = cos(i*phi);
                    shapes(2*i  ) = sin(i*phi);
                  }
              }
            else
              {
                LegendrePolynomial (order, -1 + 2 * phi, shapes);
              }
          }
        else // DIM == 2
          {
            Vec<2, double> phi;
            fes->mapping->Evaluate(mip, FlatVector<double>(2, &phi[0]));
            if(fes->polar)
              {
                ArrayMem<double, 20> shapeu(order/2+1), shapev(2*order);
                LegendrePolynomial(order/2, -1 + 2 * phi[0]*phi[0],
                                   shapeu);
                for(auto i : Range(1, order+1))
                  {
                    shapev[2*(i-1)] = cos(i*phi[1]);
                    shapev[2*(i-1)+1] = sin(i*phi[1]);
                  }
                int j = 0;
                for(auto k : Range(order/2 + 1))
                  shapes(j++) = shapeu[k];
                for(auto k : Range(1, order+1))
                  {
                    JacobiPolynomialAlpha jac(k);
                    jac.Eval(order/2, 1 - 2 * phi[0]*phi[0], shapeu);

                    for(auto p : Range((order-k)/2+1))
                      {
                        shapes(j++) = shapev[2*(k-1)] * pow(phi[0], k) * shapeu[p];
                        shapes(j++) = shapev[2*(k-1)+1] * pow(phi[0], k) * shapeu[p];
                      }
                  }
              }
            else
              {
                int ndof1 = fes->periodic[0] ? 2*order+1 : order+1;
                int ndof2 = fes->periodic[1] ? 2*order+1 : order+1;
                ArrayMem<double, 20> shapeu(ndof1), shapev(ndof2);

                if(fes->periodic[0])
                  {
                    shapeu[0] = 1;
                    for(auto i : Range(1, order+1))
                      {
                        shapeu[2*i-1] = cos(i*phi[0]);
                        shapeu[2*i  ] = sin(i*phi[0]);
                      }
                  }
                else
                  {
                    LegendrePolynomial(order, -1 + 2 * phi[0], shapeu);
                  }

                if(fes->periodic[1])
                  {
                    shapev[0] = 1;
                    for(auto i : Range(1, order+1))
                      {
                        shapev[2*i-1] = cos(i*phi[1]);
                        shapev[2*i  ] = sin(i*phi[1]);
                      }
                  }
                else
                  {
                    LegendrePolynomial(order, -1 + 2 * phi[1], shapev);
                  }

                for(auto i : Range(shapeu.Size()))
                  for(auto j : Range(shapev.Size()))
                    shapes(i * shapev.Size() + j) = shapeu[i] * shapev[j];
              }
          }
      }
    };

    class VolFE : public InterfaceFE
    {
      ArrayMem<InterfaceFE *, 4> facets;
    public:

      // Why do we need these usings? CL
      using FiniteElement::order;
      using FiniteElement::ndof;
      using InterfaceFE::fes;

      VolFE (const GlobalInterfaceSpaceD<DIM> * afes, int _order,
             ELEMENT_TYPE et)
        : InterfaceFE(afes, et)
      {
        order = _order;
      }
      using INTERFACEFE = InterfaceFE;

      void SetFacet(int nr, InterfaceFE* fe)
      {
        int old_size = facets.Size();
        if (nr+1 > old_size)
          {
            facets.SetSize(nr+1);
            for(auto i : Range(old_size, nr))
              facets[i] = nullptr;
          }
        facets[nr] = fe;
      }
      InterfaceFE* GetFacet(int nr) const
      {
        if(facets.Size() <= nr) return nullptr;
        return facets[nr];
      }

      void ComputeNDof()
      {
        if(order > 0)
          {
            if constexpr(DIM ==1)
              ndof = fes->periodic[0] ? 2*order+1 : order+1;
            else // DIM == 2
              {
                ndof = (fes->periodic[0] ? 2*order+1 : order+1) * (fes->periodic[1] ? 2*order+1 : order+1);
                if(fes->polar)
                  {
                    ndof = order/2+1;
                    for(auto k : Range(1, order+1))
                      ndof += 2*((order-k)/2 + 1);
                  }
              }
          }
        else
          ndof = 0;
        for (auto f : facets)
          if(f)
            {
              ndof = f->GetNDof();
              order = f->Order();
            }
      }

      void CalcDShape(const BaseMappedIntegrationPoint& mip,
                      BareSliceVector<double> dshapes) const
      {
        int order = fes->order;
        double phi = fes->mapping->Evaluate(mip);
        if constexpr(DIM == 1)
          {
            if(fes->periodic[0])
              {
                throw Exception("CalcDShape not implemented for periodic!");
              }
            else
              {
                AutoDiff<1> adphi(phi, 0);
                LegendrePolynomial (order, 2*adphi - 1,
                                    SBLambda([&](int i, auto val)
                                    {
                                      dshapes[i] = val.DValue(0);
                                    }));
              }
          }
        else // DIM == 2
          {
            throw Exception("CalcDShape not implemented for 2d space!");
          }
      }
    };

  public:
    GlobalInterfaceSpaceD(shared_ptr<MeshAccess> ama, const Flags & flags)
      : GlobalInterfaceSpace(ama, flags)
    {
      if(DIM == 1)
        SetNDof(periodic[0] ? 2*order+1 : order+1);
      else // DIM == 2
        {
          if(polar)
            {
              int ndof = order/2+1;
              for(auto k : Range(1, order+1))
                ndof += 2*((order-k)/2 + 1);
              SetNDof(ndof);
            }
          else
            SetNDof((periodic[0] ? 2*order+1 : order+1) * (periodic[1] ? 2*order+1 : order+1));
        }
    
      evaluator[VOL] = make_shared<VolDiffOp<VolFE>>();
      evaluator[BND] = make_shared<InterfaceDiffOp<InterfaceFE>>();
      additional_evaluators.Set("ParameterGrad", make_shared<ParameterGradDiffOp<VolFE>>());
    }

    void Update() override
    {
      GlobalInterfaceSpace::Update();

      nitsche_facet.SetSize(ma->GetNNodes(NT_FACET));
      nitsche_facet = false;

      if(definedon[BND].Size() == 0)
        nitsche_facet = true;
      else
        {
          for (auto el : ma->Elements(BND))
            if(definedon[BND][el.GetIndex()])
              nitsche_facet[el.Facets()] = true;
        }
    }
  
    void GetDofNrs (ElementId ei, Array<DofId> & dofs) const override
    {
      dofs.SetSize0();
      switch (ei.VB())
        {
        case VOL:
          {
            auto ngel = ma->GetElement(ei);
            if(definedon[VOL].Size() == 0 || definedon[VOL][ngel.GetIndex()])
              {
                dofs += IntRange(GetNDof());
              }
            else
              {
                for (auto f : ngel.Facets())
                  if (nitsche_facet[f])
                    {
                      dofs += IntRange(GetNDof());
                      break;
                    }
              }
            break;
          }

        case BND:
          {
            //if (ma->GetElement(ei).GetIndex() == bc)
            if (nitsche_facet[ma->GetElement(ei).Facets()[0]])
              dofs += IntRange(GetNDof());
            break;
          }
        default:
          ;
        }
    }

    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override
    {
      switch (ei.VB())
        {
        case VOL:
          {
            auto ngel = ma->GetElement(ei);
            auto & fe = * new (alloc) VolFE(this, DefinedOn(ei) ? order : 0, ngel.GetType());
            const auto& facets = ngel.Facets();
            for(auto i : Range(facets))
              {
                auto f = facets[i];
                if (nitsche_facet[f])
                  fe.SetFacet(i, new (alloc) InterfaceFE(this, ngel.GetType()));
              }
            fe.ComputeNDof();
            return fe;
          }
        case BND:
          {
            // if (ma->GetElement(ei).GetIndex() == bc)
            if (nitsche_facet[ma->GetElement(ei).Facets()[0]])
              {
                return * new (alloc) InterfaceFE(this, ma->GetElement(ei).GetType());
              }
            return * new (alloc) DummyFE<ET_SEGM>();
          }
        default:
          throw Exception ("Nitsche::GetFE(): no other elements");
        }
    
    }

  };

  shared_ptr<GlobalInterfaceSpace> CreateGlobalInterfaceSpace
    (shared_ptr<MeshAccess> ma, shared_ptr<CoefficientFunction> mapping,
     optional<Region> definedon, bool periodic, bool periodicu,
     bool periodicv, int order, bool complex, bool polar, bool autoupdate)
  {
    Flags flags;
    if(complex)
      flags.SetFlag("complex");
    flags.SetFlag("mapping", mapping);
    if(periodic)
      flags.SetFlag("periodic");
    if(periodicu)
      flags.SetFlag("periodicu");
    if(periodicv)
      flags.SetFlag("periodicv");
    if(definedon.has_value())
        flags.SetFlag("definedon", definedon.value());
    if(polar)
      {
        flags.SetFlag("polar");
        flags.SetFlag("periodicv");
      }
    if(autoupdate)
      flags.SetFlag("autoupdate");
    flags.SetFlag("order", order);
    if(mapping->Dimension() == 1)
      return make_shared<GlobalInterfaceSpaceD<1>>(ma, flags);
    else if(mapping->Dimension() == 2)
      return make_shared<GlobalInterfaceSpaceD<2>>(ma, flags);
    throw Exception("Mapping must be 1 or 2 dimensional!");
  }
} // namespace ngcomp
