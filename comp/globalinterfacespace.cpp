#include <comp.hpp>
#include <python_comp.hpp>


namespace ngcomp
{


  // starting with 1D periodic

  class GlobalInterfaceSpace1DPeriodic : public FESpace
  {
    shared_ptr<CoefficientFunction> mapping;
    int order;
    Array<bool> nitsche_facet;

  
    class InterfaceFE : public FiniteElement
    {
      const GlobalInterfaceSpace1DPeriodic * fes;
    public:
      InterfaceFE (const GlobalInterfaceSpace1DPeriodic * afes)
        : FiniteElement(afes->GetNDof(), afes->order), fes(afes) { ; }
    
      virtual ELEMENT_TYPE ElementType() const { return ET_SEGM; }
    
      virtual void CalcShape (const BaseMappedIntegrationPoint & mip,
                              SliceVector<double> shapes) const
      {
        int order = fes->order;
        double phi = fes->mapping -> Evaluate(mip);
      
        shapes(0) = 1;
        for (int i = 1; i <= order; i++)
          {
            shapes(2*i-1) = cos(i*phi);
            shapes(2*i  ) = sin(i*phi);
          }
      }
    };  

    class InterfaceDiffOp : public DifferentialOperator
    {
    public:
      InterfaceDiffOp () 
        : DifferentialOperator(1, 1, BND, 0) { ; }
    
      virtual void
      CalcMatrix (const FiniteElement & fel,
                  const BaseMappedIntegrationPoint & mip,
                  SliceMatrix<double,ColMajor> mat,   
                  LocalHeap & lh) const
      {
        if (fel.GetNDof() > 0)
          dynamic_cast<const InterfaceFE&> (fel).CalcShape(mip, mat.Row(0));
      }
    };


    class VolFE : public FiniteElement
    {
    public:
      InterfaceFE *bndels[3] = { nullptr, nullptr, nullptr };

      void ComputeNDof()
      {
        ndof = 0;
        order = 0;
        for (int i = 0; i < 3; i++)
          if (bndels[i])
            {
              ndof += bndels[i] -> GetNDof();
              order = max2 (order, bndels[i]->Order());
            }
      }
    
      virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }
    };


    class VolDiffOp : public DifferentialOperator
    {
    public:
      VolDiffOp () 
        : DifferentialOperator(1, 1, BND, 0) { ; }

      virtual void
      CalcMatrix (const FiniteElement & bfel,
                  const BaseMappedIntegrationPoint & mip,
                  SliceMatrix<double,ColMajor> mat,   
                  LocalHeap & lh) const
      {
        auto & ip = mip.IP();
        auto & fel = dynamic_cast<const VolFE&> (bfel);
        mat = 0;

        int fnr = ip.FacetNr();

        if (fnr != -1 && fel.bndels[fnr])
          dynamic_cast<const InterfaceFE&> (*fel.bndels[fnr]).CalcShape(mip, mat.Row(0));
      }
    };

  
  
  public:
    GlobalInterfaceSpace1DPeriodic (shared_ptr<MeshAccess> ama, const Flags & flags)
      : FESpace(ama, flags)
    {
      order = int(flags.GetNumFlag("order", 3));
      try
        {
          mapping = std::any_cast<shared_ptr<CoefficientFunction>>(flags.GetAnyFlag("mapping"));
        }
      catch(std::bad_any_cast ex)
        {
          cout << "No mapping or wrong mapping given!" << endl;
          cout << "GlobalInterfaceSpace needs kwarg: mapping=(interface, mapping_func)" << endl;
        }
    
      SetNDof( (2*order+1) );
    
      evaluator[VOL] = make_shared<VolDiffOp>();
      evaluator[BND] = make_shared<InterfaceDiffOp>();
    }

    virtual void Update() 
    {
      FESpace::Update();

      nitsche_facet.SetSize(ma->GetNNodes(NT_FACET));
      nitsche_facet = false;

      for (auto el : ma->Elements(BND))
        if(definedon[BND][el.GetIndex()])
          nitsche_facet[el.Facets()] = true;
    }
  
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dofs) const
    {
      // cout << "dofid = " << ei << ", dofnrs = " << endl;
      dofs.SetSize0();
      switch (ei.VB())
        {
        case VOL:
          {
            auto ngel = ma->GetElement(ei);
            for (auto f : ngel.Facets())
              if (nitsche_facet[f])
                {
                  dofs += IntRange(0, GetNDof());
                  break;
                }
            break;
          }

        case BND:
          {
            //if (ma->GetElement(ei).GetIndex() == bc)
            if (nitsche_facet[ma->GetElement(ei).Facets()[0]])
              dofs += IntRange(0, GetNDof());
            break;
          }
        default:
          ;
        }
      // cout << dofs << endl;
    }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const
    {
      switch (ei.VB())
        {
        case VOL:
          {
            auto ngel = ma->GetElement(ei);
            auto & fe = * new (alloc) VolFE();
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
                // cout << "get interfacefe" << endl;
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

  
  void ExportGlobalInterfaceSpaces(py::module & m)
  {
    ExportFESpace<GlobalInterfaceSpace1DPeriodic>
      (m, "GlobalInterfaceSpace1DPeriodic")
      .def_static("__special_treated_flags__", [] ()
      {
        py::dict special(**py::module::import("ngsolve").attr("FESpace").attr("__special_treated_flags__")(),
        py::arg("mapping") = py::cpp_function
          ([] (shared_ptr<CoefficientFunction> mapping,
               Flags* flags, py::list info)
          {
            flags->SetFlag("mapping", mapping);
          }));
        return special;
      })
      ;
  }

}
  
