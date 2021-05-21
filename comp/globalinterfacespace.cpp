#include <comp.hpp>
#include "globalinterfacespace.hpp"

namespace ngcomp
{
  template<typename VOLFE>
  void GlobalInterfaceSpace::VolDiffOp<VOLFE>::CalcMatrix
    (const FiniteElement & bfel,
     const BaseMappedIntegrationPoint & mip,
     SliceMatrix<double,ColMajor> mat,
     LocalHeap & lh) const
  {
    auto & ip = mip.IP();
    auto & fel = dynamic_cast<const VOLFE&> (bfel);
    mat = 0;

    int fnr = ip.FacetNr();

    if (fnr != -1 && fel.bndels[fnr])
      {
        dynamic_cast<const typename VOLFE::INTERFACEFE&> (*fel.bndels[fnr]).CalcShape(mip, mat.Row(0));
        return;
      }
    if(fnr == -1)
      {
        fel.CalcShape(mip, mat.Row(0));
      }
  }

  template<typename INTERFACEFE>
  void GlobalInterfaceSpace::InterfaceDiffOp<INTERFACEFE>::CalcMatrix
    (const FiniteElement & fel,
     const BaseMappedIntegrationPoint & mip,
     SliceMatrix<double,ColMajor> mat,
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
    if(flags.GetDefineFlag("periodic"))
      periodic[0] = periodic[1] = true;
    if(flags.GetDefineFlag("periodicx"))
      periodic[0] = true;
    if(flags.GetDefineFlag("periodicy"))
      periodic[1] = true;
    try
      {
        mapping = std::any_cast<shared_ptr<CoefficientFunction>>(flags.GetAnyFlag("mapping"));
      }
    catch(std::bad_any_cast ex)
      {
        throw Exception("No mapping or wrong mapping given!\nGlobalInterfacespace needs kwarg: mapping=CoefficientFunction");
      }
  }

  class GlobalInterfaceSpace1D : public GlobalInterfaceSpace
  {
    Array<bool> nitsche_facet;

    class InterfaceFE : public FiniteElement
    {
      const GlobalInterfaceSpace1D * fes;
    public:
      InterfaceFE (const GlobalInterfaceSpace1D * afes)
        : FiniteElement(afes->GetNDof(), afes->order), fes(afes) { ; }

      virtual ELEMENT_TYPE ElementType() const override { return ET_SEGM; }

      virtual void CalcShape (const BaseMappedIntegrationPoint & mip,
                              SliceVector<double> shapes) const
      {
        int order = fes->order;
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
    };

    class VolFE : public FiniteElement
    {
      const GlobalInterfaceSpace1D * fes;
    public:
      VolFE (const GlobalInterfaceSpace1D * afes, int order)
        : FiniteElement(afes->GetNDof(), order), fes(afes) { ; }
      using INTERFACEFE = InterfaceFE;
      InterfaceFE *bndels[3] = { nullptr, nullptr, nullptr };

      void ComputeNDof()
      {
        if(order > 0)
          {
            if(fes->periodic[0])
              ndof = 2*order+1;
            else
              ndof = order+1;
          }
        else
          ndof = 0;
        for (int i = 0; i < 3; i++)
          if (bndels[i])
            {
              ndof = bndels[i] -> GetNDof();
              order = bndels[i]->Order();
            }
      }

      void CalcShape(const BaseMappedIntegrationPoint& mip,
                     SliceVector<double> shapes) const
      {
        int order = fes->order;
        double phi = fes->mapping->Evaluate(mip);
        if(fes->periodic[0])
          {
            shapes(0) = 1.;
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

      virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }
    };

  public:
    GlobalInterfaceSpace1D (shared_ptr<MeshAccess> ama, const Flags & flags)
      : GlobalInterfaceSpace(ama, flags)
    {
      if(periodic[0])
        SetNDof(2*order+1);
      else
        SetNDof(order+1);
    
      evaluator[VOL] = make_shared<VolDiffOp<VolFE>>();
      evaluator[BND] = make_shared<InterfaceDiffOp<InterfaceFE>>();
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
            auto & fe = * new (alloc) VolFE(this, DefinedOn(ei) ? order : 0);
            for (size_t i = 0; i < 3; i++)
              {
                auto f = ngel.Facets()[i];
                if (nitsche_facet[f])
                  fe.bndels[i] =  new (alloc) InterfaceFE(this);
              }
            fe.ComputeNDof();
            return fe;
          }
        case BND:
          {
            // if (ma->GetElement(ei).GetIndex() == bc)
            if (nitsche_facet[ma->GetElement(ei).Facets()[0]])
              {
                return * new (alloc) InterfaceFE(this);
              }
            return * new (alloc) DummyFE<ET_SEGM>();
          }
        default:
          throw Exception ("Nitsche::GetFE(): no other elements");
        }
    
    }

    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.Arg("mapping") = "Mapping for global interface space.";
      docu.Arg("periodic") = "Periodic global interface space (in 2d in x and y direction).";
      docu.Arg("periodicx") = "Periodic x-dir (local coordinate system) global interface space.";
      docu.Arg("periodicy") = "Periodic y-dir (local coordinate system) global interface space.";
      return docu;
    }
  };

  /*

  // old code for cylinder 
 
  // scalar cylinder 
  class NitscheSpaceCylinder : public FESpace
  {
  int order;
  int bc;
  Vec<3> pa, pb;
  double r;
  Array<bool> nitsche_face;
  
  class NitscheFE : public FiniteElement
  {
  const NitscheSpaceCylinder * fes;
  public:
  NitscheFE (const NitscheSpaceCylinder * afes)
  : FiniteElement(afes->GetNDof(), afes->order), fes(afes) { ; }
    
  virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }
    
  virtual void CalcShape (const BaseMappedIntegrationPoint & mip,
  SliceVector<double> shapes) const
  {
  int order = fes->order;
  double x = mip.GetPoint()(0);
  double y = mip.GetPoint()(1);
  double z = mip.GetPoint()(2);
  double zr = -1 + 2*z;
  double phi = atan2 (x,y);
  Vector<> shapez(order+1), shapephi(2*order+1);
  LegendrePolynomial (order, zr, shapez);
  shapephi(0) = 1;
  for (int i = 1; i <= order; i++)
  {
  shapephi(2*i-1) = cos(i*phi);
  shapephi(2*i  ) = sin(i*phi);
  }
      
  for (int i = 0, ii = 0; i < shapez.Size(); i++)
  for (int j = 0; j < shapephi.Size(); j++, ii++)
  shapes(ii) = shapez(i)*shapephi(j);
  }
  };  

  class NitscheDiffOp : public DifferentialOperator
  {
  public:
  NitscheDiffOp () 
  : DifferentialOperator(1, 1, BND, 0) { ; }

  virtual void
  CalcMatrix (const FiniteElement & fel,
  const BaseMappedIntegrationPoint & mip,
  SliceMatrix<double,ColMajor> mat,   
  LocalHeap & lh) const
  {
  if (fel.GetNDof() > 0)
  dynamic_cast<const NitscheFE&> (fel).CalcShape(mip, mat.Row(0));
  }
  };


  class NitscheVolFE : public FiniteElement
  {
  public:
  NitscheFE *bndels[4] = { nullptr, nullptr, nullptr, nullptr };

  void ComputeNDof()
  {
  ndof = 0;
  order = 0;
  for (int i = 0; i < 4; i++)
  if (bndels[i])
  {
  ndof += bndels[i] -> GetNDof();
  order = max2 (order, bndels[i]->Order());
  }
  }
    
  virtual ELEMENT_TYPE ElementType() const { return ET_TET; }
  };


  class NitscheVolDiffOp : public DifferentialOperator
  {
  public:
  NitscheVolDiffOp () 
  : DifferentialOperator(1, 1, BND, 0) { ; }

  virtual void
  CalcMatrix (const FiniteElement & bfel,
  const BaseMappedIntegrationPoint & mip,
  SliceMatrix<double,ColMajor> mat,   
  LocalHeap & lh) const
  {
  auto & ip = mip.IP();

  auto & fel = dynamic_cast<const NitscheVolFE&> (bfel);
  mat = 0;

  int fnr = ip.FacetNr();
  if (fnr != -1 && fel.bndels[fnr])
  dynamic_cast<const NitscheFE&> (*fel.bndels[fnr]).CalcShape(mip, mat.Row(0));
  }
  };

  
  
  public:
  NitscheSpaceCylinder (shared_ptr<MeshAccess> ama, const Flags & flags)
  : FESpace(ama, flags)
  {
  cout << "nitsche flags = " << endl << flags << endl;
  order = int(flags.GetNumFlag("order", 3));
  bc = int(flags.GetNumFlag("bc", 1)) - 1;
    
  SetNDof( (order+1) * (2*order+1) );
    
  evaluator[VOL] = make_shared<NitscheVolDiffOp>();
  evaluator[BND] = make_shared<NitscheDiffOp>();
  }

  
  virtual void Update() 
  {
  FESpace::Update();
  nitsche_face.SetSize(ma->GetNNodes(NT_FACE));
  nitsche_face = false;
  for (auto el : ma->Elements(BND))
  if (el.GetIndex() == bc)
  nitsche_face[el.Faces()] = true;
  }
  
  virtual void GetDofNrs (ElementId ei, Array<DofId> & dofs) const
  {
  dofs.SetSize0();
  switch (ei.VB())
  {
  case VOL:
  {
  auto ngel = ma->GetElement(ei);
  for (size_t i = 0; i < 4; i++)
  {
  auto f = ngel.Faces()[i];
  if (nitsche_face[f])
  dofs += IntRange(0, GetNDof());
  }
  break;
  }

  case BND:
  {
  if (ma->GetElement(ei).GetIndex() == bc)
  dofs += IntRange(0, GetNDof());
  break;
  }
  default:
  ;
  }
  }

  virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const
  {
  switch (ei.VB())
  {
  case VOL:
  {
  auto ngel = ma->GetElement(ei);
  auto & fe = * new (alloc) NitscheVolFE();
  for (size_t i = 0; i < 4; i++)
  {
  auto f = ngel.Faces()[i];
  if (nitsche_face[f])
  fe.bndels[i] =  new (alloc) NitscheFE(this);
  }
  fe.ComputeNDof();
  return fe;
  }
  case BND:
  {
  if (ma->GetElement(ei).GetIndex() == bc)
  return * new (alloc) NitscheFE(this);
  return * new (alloc) DummyFE<ET_TRIG>();
  }
  default:
  throw Exception ("Nitsche::GetFE(): no other elements");
  }
    
  }

  
  
  };

  static RegisterFESpace<NitscheSpaceCylinder> initifes ("NitscheSpaceCylinder");
  */

  shared_ptr<GlobalInterfaceSpace> CreateGlobalInterfaceSpace
    (shared_ptr<MeshAccess> ma, shared_ptr<CoefficientFunction> mapping,
     optional<Region> definedon, bool periodic, bool periodicx,
     bool periodicy)
  {
    Flags flags;
    flags.SetFlag("mapping", mapping);
    if(periodic)
      flags.SetFlag("periodic");
    if(periodicx)
      flags.SetFlag("periodicx");
    if(periodicy)
      flags.SetFlag("periodicy");
    if(definedon.has_value())
        flags.SetFlag("definedon", definedon.value());
    if(ma->GetDimension() == 2)
      return make_shared<GlobalInterfaceSpace1D>(ma, flags);
    throw Exception("this space is not yet implemented!");
  }
} // namespace ngcomp
