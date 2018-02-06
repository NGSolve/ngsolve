#include <comp.hpp>
#include <fem.hpp> 

#include "../fem/l2hofe.hpp"
#include "../fem/diffop_impl.hpp"
#include "../fem/facethofe.hpp"


namespace ngcomp
{ 

  /// Identity
  template <int D>
  class DiffOpIdFacet : public DiffOp<DiffOpIdFacet<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
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
          const FacetVolumeFiniteElement<D> & fel_facet = static_cast<const FacetVolumeFiniteElement<D>&> (bfel);
	  
          fel_facet.Facet(facetnr).CalcShape(mip.IP(), 
                                             mat.Row(0).Range(fel_facet.GetFacetDofs(facetnr)));
        }
      else
        {
          if (mip.BaseMappedIntegrationPoint::VB() == BND) 
            {
              const BaseScalarFiniteElement & fel = static_cast<const BaseScalarFiniteElement&> (bfel);
              fel.CalcShape (mip.IP(), mat.Row(0));
            }
          else
            throw Exception("cannot evaluate facet-fe inside element");
        }
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
    
  };

  template <int D>
  class FacetSurfaceMassIntegrator 
    : public T_BDBIntegrator<DiffOpIdFacetSurface<D>, DiagDMat<DiffOpIdFacetSurface<D>::DIM_DMAT>, FiniteElement>
  {
    typedef T_BDBIntegrator<DiffOpIdFacetSurface<D>, DiagDMat<DiffOpIdFacetSurface<D>::DIM_DMAT>, FiniteElement> BASE;
  public:
    using  T_BDBIntegrator<DiffOpIdFacetSurface<D>, DiagDMat<DiffOpIdFacetSurface<D>::DIM_DMAT>, FiniteElement>::T_BDBIntegrator;

    virtual string Name () const { return "FacetSurface-Mass"; }
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
    
    
    ndlevel.SetSize(0);
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
        throw Exception("FacetSurfaceFESpace only implemented for 3d!");
      }
    else
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdFacet<3>>>();
	evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdFacetSurface<3>>>();
	evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdFacetSurfaceBoundary<3>>>();
	
        integrator[BND] = make_shared<RobinIntegrator<3>> (one);
      }

  }
  

  FacetSurfaceFESpace :: ~FacetSurfaceFESpace ()
  { ; }


  void FacetSurfaceFESpace :: Update(LocalHeap & lh)
  {
    FESpace :: Update (lh);
    
    if(print) 
      *testout << " FacetSurfaceFEspace with order " << order << " rel_order " << rel_order << " var_order " << var_order << endl; 

    nel = ma->GetNSE();
    nfa = ma->GetNEdges(); 
    ndof = 0;
    first_edge_dof.SetSize(nfa +1);
    first_edge_dof = ndof;
        
    
    first_edge_dof.SetSize(nfa+1); 
    first_edge_dof = nfa;

    if(ma->GetDimension() == 3)
      {
        for(int i = 0; i < nfa; i++)
          {
            first_edge_dof[i] = ndof;
            ndof += order+1;
        
          }
        first_edge_dof[nfa] = ndof;              
      }
    else
      {
	throw Exception("Only implemented for 3d!");
      }
    
    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;
    
    if(print)
      {
	*testout << "*** Update FacetSurfaceFESpace: General Information" << endl;
	*testout << " order edge (edge) " << order << endl; 
	*testout << " first_edge_dof (edge)  " << first_edge_dof << endl; 
      } 
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

	  switch (ma->GetElType(ei))
	    {
	    case ET_TRIG: return T_GetFE<ET_TRIG>(ei.Nr(), lh);//fe = new (lh) FacetFE<ET_TRIG> (); break;
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
	  DGFiniteElement<1> * fe1d = 0;
	  
	  switch (ma->GetElType(ei))
	    {
	    case ET_SEGM: fe1d = new (lh) L2HighOrderFE<ET_SEGM> (); break;
	    default:
	      throw Exception (string("FacetSurfaceFESpace::GetFE: unsupported element ")+
			       ElementTopology::GetElementName(ma->GetElType(ei)));
	    }
	  
	  fe1d -> SetVertexNumbers (vnums);
	  fe1d -> SetOrder (order); 
	  fe1d -> ComputeNDof();
	  return *fe1d;
	  break;
	}
      case BBBND:
        //throw Exception("No BBBND GetFE implemented for FacetSurfaceFESpace");
	return * new (lh) DummyFE<ET_POINT>();
	break;
      }
  }


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


  // ------------------------------------------------------------------------
  void FacetSurfaceFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    dnums.SetSize0();
    switch (ei.VB())
      {
      case VOL:
	break;
     
      case BND:
	{
	  int fnum = 0;
	  if (ma->GetDimension() == 3)
	    {
	      auto ednums = ma->GetElEdges (ei);
	      for( auto ed : ednums)
		{
		  dnums += GetEdgeDofs(ed);
		}
	    }
	}
	break;
      case BBND:
	{
	  dnums += GetEdgeDofs(ma->GetElEdges(ei)[0]);
	}
	break;

      case BBBND:
	break;
      }
  }
  
  // ------------------------------------------------------------------------

  static RegisterFESpace<FacetSurfaceFESpace> init_facet ("facetsurface");
}
