#include <comp.hpp>
#include <fem.hpp> 

#include "../fem/l2hofe.hpp"
#include "../fem/diffop_impl.hpp"
#include "../fem/facethofe.hpp"


namespace ngcomp
{ 

  template <int D>
  class DiffOpIdFacetSurface;

  /// Identity
  template <int D>
  class DiffOpIdFacet_ : public DiffOp<DiffOpIdFacet_<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    static bool SupportsVB (VorB checkvb) { return true; }

    typedef DiffOpIdFacetSurface<D> DIFFOP_TRACE;
    
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      int facetnr = mip.IP().FacetNr();
      if (facetnr >= 0)
        {
          mat = 0.0;
          const FacetVolumeFiniteElement<D> & fel_facet = static_cast<const FacetVolumeFiniteElement<D>&> (bfel);
	  
          fel_facet.Facet(facetnr).CalcShape(mip.IP(), 
                                             mat.Row(0).Range(fel_facet.GetFacetDofs(facetnr)));
        }
      else
        {
          // if (mip.BaseMappedIntegrationPoint::VB() == BND)
          if (mip.IP().VB() == BND) 
            {
              const BaseScalarFiniteElement & fel = static_cast<const BaseScalarFiniteElement&> (bfel);
              fel.CalcShape (mip.IP(), mat.Row(0));
            }
          else
            throw Exception("cannot evaluate facet-fe inside element");
        }
    }

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian)
    {
      if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpIdFacet_");      
      return ZeroCF(Array<int>());
    }

  };
  
  /// Identity
  template <int D>
  class DiffOpIdFacetSurface : public DiffOp<DiffOpIdFacetSurface<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      int facetnr = mip.IP().FacetNr();
      if (facetnr >= 0)
        {
          mat = 0.0;
          const FacetVolumeFiniteElement<D-1> & fel_facet = static_cast<const FacetVolumeFiniteElement<D-1>&> (bfel);
          fel_facet.Facet(facetnr).CalcShape(mip.IP(), 
                                             mat.Row(0).Range(fel_facet.GetFacetDofs(facetnr)));
        }
      else
        {
            throw Exception("cannot evaluate facet-fe inside element");
        }
    }


   static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      int facetnr = mir.IR()[0].FacetNr();
      if (facetnr >= 0)
        {
          mat.AddSize(fel.GetNDof(), mir.Size()) = 0.0;
          const FacetVolumeFiniteElement<D-1> & fel_facet = static_cast<const FacetVolumeFiniteElement<D-1>&> (fel);
          fel_facet.Facet(facetnr).CalcShape(mir.IR(), 
                                             mat.Rows(fel_facet.GetFacetDofs(facetnr)));
        }
      else
        {
          throw ExceptionNOSIMD("facet-simd-bnd not ready");
        }
    }

    
    using DiffOp<DiffOpIdFacetSurface<D>>::ApplySIMDIR;          
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      const FacetVolumeFiniteElement<D-1> & fel_facet = static_cast<const FacetVolumeFiniteElement<D-1>&> (bfel);

      int facetnr = mir.IR()[0].FacetNr();
      if (facetnr < 0)
        throw Exception("cannot evaluate facet-fe inside element, apply simd");
      else
        fel_facet.Facet(facetnr).Evaluate(mir.IR(),
                                          x.Range(fel_facet.GetFacetDofs(facetnr)),
                                          y.Row(0));
    }

    using DiffOp<DiffOpIdFacetSurface<D>>::AddTransSIMDIR;          
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      const FacetVolumeFiniteElement<D-1> & fel_facet = static_cast<const FacetVolumeFiniteElement<D-1>&> (bfel);

      int facetnr = mir.IR()[0].FacetNr();
      if (facetnr < 0)
        throw Exception("cannot evaluate facet-fe inside element, add trans simd");
      else
        fel_facet.Facet(facetnr).AddTrans(mir.IR(),
                                          y.Row(0),
                                          x.Range(fel_facet.GetFacetDofs(facetnr)));
    }

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian)
    {
      if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpIdFacetSurface");      
      return ZeroCF(Array<int>());
    }

    
  };


    /// Identity on boundary
  template <int D, typename FEL = ScalarFiniteElement<D-2> >
  class DiffOpIdFacetSurfaceBoundary : public DiffOp<DiffOpIdFacetSurfaceBoundary<D, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-2 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    static const FEL & Cast (const FiniteElement & fel) 
    { return static_cast<const FEL&> (fel); }

    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcShape (mip.IP(), mat.Row(0));
    }

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian)
    {
      if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpIdFacetSurface");      
      return ZeroCF(Array<int>());
    }


  };





  template <int D>
  class DiffOpGradFacetSurface : public DiffOp<DiffOpGradFacetSurface<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };

    static string Name() { return "grad"; }

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      int facetnr = mip.IP().FacetNr();
      if (facetnr >= 0)
        {
          HeapReset hr(lh);
          
          const FacetVolumeFiniteElement<D-1> & fel_facet = static_cast<const FacetVolumeFiniteElement<D-1>&> (bfel);
          auto r = fel_facet.GetFacetDofs(facetnr);
          
          FlatMatrix<> dshaperef(r.Size(), DIM_ELEMENT, lh);
          mat = 0.0;
          /*
          fel_facet.Facet(facetnr).CalcDShape(mip.IP(), 
                                              mat.Row(0).Range(fel_facet.GetFacetDofs(facetnr)));
          */
          fel_facet.Facet(facetnr).CalcDShape(mip.IP(), dshaperef);
          mat.Cols(r) = Trans(mip.GetJacobianInverse()) * Trans(dshaperef);
        }
      else
        {
            throw Exception("cannot evaluate facet-fe inside element");
        }
    }
  };





  

  FacetSurfaceFESpace ::  FacetSurfaceFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags)
    : FESpace(ama, flags)
  {
    name="FacetSurfaceFESpace(facet)";
    type = "facetsurface";
    // defined flags
    DefineNumFlag("relorder");
    DefineDefineFlag("variableorder"); 


    if(checkflags) CheckFlags(flags);
    
    
    // ndlevel.SetSize(0);
    Flags loflags;
    loflags.SetFlag("order",0.0);
    if ( this->IsComplex() )
      loflags.SetFlag("complex");

    // Variable order space: 
    //      in case of (var_order && order) or (relorder) 
    var_order = flags.GetDefineFlag("variableorder");  
    order =  int (flags.GetNumFlag ("order",0)); 

    if(flags.NumFlagDefined("relorder") && !flags.NumFlagDefined("order")) 
      var_order = 1; 
    
    rel_order=int(flags.GetNumFlag("relorder",order-1)); 


    if(flags.NumFlagDefined("order") && flags.NumFlagDefined("relorder")) 
      {
	if(var_order)
	  cerr << " WARNING: FacetSurfaceFESpace: inconsistent flags: variableorder, order and relorder "
	       << "-> variable order space with rel_order " << rel_order << "is used, but order is ignored " << endl; 
	else 
	  cerr << " WARNING: FacetSurfaceFESpace: inconsistent flags: order and rel_order "
	       << "-> uniform order space with order " << order << " is used " << endl; 
      }

    if (flags.NumFlagDefined("order")) 
      { 
	if(var_order) 
	  { 
	    rel_order = int(flags.GetNumFlag("relorder",order-1)); 
	    order = rel_order + 1;
	  }
	else 
	  order =  int (flags.GetNumFlag ("order",0));
      }
    else if(flags.NumFlagDefined("relorder"))
      {
	var_order=1; 
	rel_order = int (flags.GetNumFlag ("relorder",-1));
	order=1+rel_order; 
      }
    else // neither rel_order nor order is given  
      {
	rel_order = -1;  
	order = 0;  
      }

    
    nowirebasket = flags.GetDefineFlag ("nowirebasket");
    
    auto one = make_shared<ConstantCoefficientFunction>(1);
    if (ma->GetDimension() == 2)
      {
        //throw Exception("FacetSurfaceFESpace only implemented for 3d!");
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdFacet_<2>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdFacetSurface<2>>>();
        evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdFacetSurfaceBoundary<2>>>();
	
        integrator[BND] = make_shared<RobinIntegrator<3>> (one);
      }
    else if (ma->GetDimension() == 3)
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdFacet_<3>>>();
      	evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdFacetSurface<3>>>();
      	evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdFacetSurfaceBoundary<3>>>();

      	flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpGradFacetSurface<3>>>();
        integrator[BND] = make_shared<RobinIntegrator<3>> (one);
      }
    else
    {
      throw Exception("FacetSurfaceFESpace only implemented for 2d and 3d");
    }

    additional_evaluators.Set ("dual", evaluator[VOL]);
    //additional_evaluators.Set ("dual", evaluator[BND]);
  }
  

  FacetSurfaceFESpace :: ~FacetSurfaceFESpace ()
  { ; }


  void FacetSurfaceFESpace :: Update()
  {
    FESpace::Update();
    
    if(print) 
      *testout << " FacetSurfaceFEspace with order " << order << " rel_order " << rel_order << " var_order " << var_order << endl; 

    nel = ma->GetNSE();
    nfa = ma->GetNEdges(); 
    first_edge_dof.SetSize(nfa+1);
    first_edge_dof = 0;
    
    // first_edge_dof.SetSize(nfa+1); 
    // first_edge_dof = nfa;

    if(ma->GetDimension() == 3)
      {
        /*
        for(int i = 0; i < nfa; i++)
          {
            first_edge_dof[i] = ndof;
            ndof += order+1;
          }
        first_edge_dof[nfa] = ndof;              
        */
        for (auto el : ma->Elements(BND))
          for (auto edge : el.Edges())
            first_edge_dof[edge] = order+1;
      }
    else if(ma->GetDimension() == 2)
      {
        for (auto el : ma->Elements(BND))
          for (auto vertex : el.Vertices())
            first_edge_dof[vertex] = 1;
      }
    else
      {
	throw Exception("Only implemented for 3d and 2d!");
      }

    size_t ndof = 0;
    for (size_t i = 0; i < nfa; i++)
      {
        size_t tmp = first_edge_dof[i];
        first_edge_dof[i] = ndof;
        ndof += tmp;
      }
    first_edge_dof[nfa] = ndof;    

    
    SetNDof (ndof);
    /*
    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;
    */
    UpdateCouplingDofArray();
    
    if(print)
      {
	*testout << "*** Update FacetSurfaceFESpace: General Information" << endl;
	*testout << " order edge (edge) " << order << endl; 
	*testout << " first_edge_dof (edge)  " << first_edge_dof << endl; 
      } 
  }

   void FacetSurfaceFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(ndof);
    ctofdof = UNUSED_DOF;

    for (ElementId ei : ma->Elements(BND))
      if (DefinedOn(ei))
        {
          if (ma->GetDimension() == 3)
            for (auto ed : ma->GetElEdges (ei))
              ctofdof[GetEdgeDofs(ed)] = WIREBASKET_DOF;
          else if (ma->GetDimension() == 2)
            for (auto ed : ma->GetElVertices (ei))
              ctofdof[GetEdgeDofs(ed)] = WIREBASKET_DOF;
        }
    
    if (print)
      *testout << "couplingtypes = " << endl << ctofdof << endl;
  }

  FlatArray<VorB> FacetSurfaceFESpace :: GetDualShapeNodes (VorB vb) const
  {
    static VorB nodes[] = { BND };
    return FlatArray<VorB> (1, &nodes[0]);
  }

    template <ELEMENT_TYPE ET>
  FiniteElement & FacetSurfaceFESpace :: T_GetFE (int elnr, Allocator & alloc) const
  {
    Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL> (elnr);

    FacetFE<ET> * fe =  new (alloc) FacetFE<ET> ();
    fe -> SetVertexNumbers (ngel.Vertices());
    fe -> SetOrder (order);
    fe -> ComputeNDof();
    
    return *fe;
  }
  
  // ------------------------------------------------------------------------
  FiniteElement & FacetSurfaceFESpace :: GetFE (ElementId ei, Allocator  & lh) const
  {
    auto vnums = ma->GetElVertices(ei);
	  
    switch(ei.VB())
      {
      case VOL:
        {
          throw Exception("Volume elements not available for FacetSurfaceSpace");
	  break;
        }            
      case BND:
        {
	  //FacetFE<ET_TRIG>* fet = 0;
	  //FacetFE<ET_QUAD>* feq = 0;

    if (!DefinedOn (ei))
      {
        return
          SwitchET (ma->GetElType(ei), [&] (auto et) -> FiniteElement&
                      {
                        return *new (lh) DummyFE<et.ElementType()> ();
                      });
      }
	  switch (ma->GetElType(ei))
	    {
	    case ET_TRIG: return T_GetFE<ET_TRIG>(ei.Nr(), lh);//fe = new (lh) FacetFE<ET_TRIG> (); break;
	    case ET_SEGM: return T_GetFE<ET_SEGM>(ei.Nr(), lh);//fe = new (lh) FacetFE<ET_TRIG> (); break;
	    case ET_QUAD: return T_GetFE<ET_QUAD>(ei.Nr(), lh);//fe = new (lh) FacetFE<ET_QUAD> (); break;
	    default:
	      throw Exception (string("FacetSurfaceFESpace::GetFE: unsupported element ")+
			       ElementTopology::GetElementName(ma->GetElType(ei)));
	    }
	  
	  /*switch (ma->GetElType(ei))
	    {
	    case ET_TRIG:
	    case ET_QUAD:
	      {
		fe -> SetVertexNumbers (vnums);
		fe -> SetOrder (order); 
		fe -> ComputeNDof();
		return *fe;
		break;
	      }
	    default:
	      break;
	    }
	    return *fe;*/
        }
      case BBND:
	{
	  switch (ma->GetElType(ei))
	    {
	    case ET_SEGM:
              {
                auto fe1d = new (lh) L2HighOrderFE<ET_SEGM> ();
                fe1d -> SetVertexNumbers (vnums);
                fe1d -> SetOrder (order); 
                fe1d -> ComputeNDof();
                return *fe1d;
              }
	    default:
	      throw Exception (string("FacetSurfaceFESpace::GetFE: unsupported element ")+
			       ElementTopology::GetElementName(ma->GetElType(ei)));
	    }
	  break;
	}
      case BBBND:
        //throw Exception("No BBBND GetFE implemented for FacetSurfaceFESpace");
	return * new (lh) DummyFE<ET_POINT>();
	break;

      default:
        __assume(false);
      }
  }


  /*
  // ------------------------------------------------------------------------
  size_t FacetSurfaceFESpace :: GetNDof () const throw()
  {
    return ndof;
  }

  // ------------------------------------------------------------------------
  size_t FacetSurfaceFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }
  */
  

  // ------------------------------------------------------------------------
  void FacetSurfaceFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    dnums.SetSize0();

    if (!DefinedOn (ei))
      return;

    switch (ei.VB())
      {
      case VOL:
	break;
     
      case BND:
        {
          if (ma->GetDimension() == 3)
            for (auto ed : ma->GetElEdges (ei))
              dnums += GetEdgeDofs(ed);
          else if (ma->GetDimension() == 2)
            for (auto ed : ma->GetElVertices (ei))
              dnums += GetEdgeDofs(ed);
          break;
        }
        
      case BBND:
	{
	  dnums += GetEdgeDofs(ma->GetElEdges(ei)[0]);
          break;
	}

      case BBBND:
	break;
      }
  }
  
  // ------------------------------------------------------------------------

  static RegisterFESpace<FacetSurfaceFESpace> init_facet ("facetsurface");
}
