/*********************************************************************/
/* File:   normalfacetsurfacefespace.cpp                             */
/* Author: Michael Neunteufel                                        */
/* Date:   2020                                                      */
/*********************************************************************/

#include <comp.hpp>
#include <fem.hpp>
#include <../fem/normalfacetfe.hpp>
#include <normalfacetsurfacefespace.hpp>
#include <../fem/hcurlhdiv_dshape.hpp> 

namespace ngcomp
{

  NormalFacetSurfaceFESpace :: NormalFacetSurfaceFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
                                                            bool parseflags )
    : FESpace(ama, flags )
  {
    type = "normalfacetsurface";
    name = "NormalFacetSurfaceFESpace";
    DefineNumFlag("relorder");
    DefineDefineFlag("variableorder");

    if (parseflags) CheckFlags(flags);

    
    order = int (flags.GetNumFlag ("order",0)); 

    if(flags.NumFlagDefined("relorder") && !flags.NumFlagDefined("order")) 
      var_order = 1; 
    else 
      var_order = 0;
    
    rel_order=int(flags.GetNumFlag("relorder",order-1)); 


    if(flags.NumFlagDefined("order") && flags.NumFlagDefined("relorder")) 
      {
	if(var_order)
	  cerr << " WARNING: NormalFacetSurfaceFESpace: inconsistent flags: variableorder, order and relorder "
	       << "-> variable order space with rel_order " << rel_order << "is used, but order is ignored " << endl; 
	else 
	  cerr << " WARNING: NormalFacetSurfaceFESpace: inconsistent flags: order and rel_order "
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


    if (ma->GetDimension() != 3)
      {
        throw Exception ("only 2D manifolds supported"); 
      }
    else
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdHDivSurface<3>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<3>>>();
      }
  }

  
  DocInfo NormalFacetSurfaceFESpace :: GetDocu()
  {
    DocInfo docu = FESpace::GetDocu();
    return docu;
  }

  void NormalFacetSurfaceFESpace :: Update()
  {
    FESpace::Update();
    if ( print ) 
      *testout << "NormalFacetSurfaceFESpace, order " << order << endl 
	       << "rel_order " << rel_order << ", var_order " << var_order << endl;

    bool first_update = GetTimeStamp() < ma->GetTimeStamp();
    if (first_update) timestamp = NGS_Object::GetNextTimeStamp();
    

    size_t nfacets = ma->GetNEdges();

    if (first_update)
      {
	int p = 0; 
	if (!var_order) p = order; 
    
	order_facet.SetSize(nfacets);
	fine_facet.SetSize(nfacets);

	order_facet = p;
	fine_facet = 0; 
    
	for (Ngs_Element el : ma->Elements<BND>())
	  if (DefinedOn(el))
	    fine_facet[el.Edges()] = true;

	for (int i = 0; i < fine_facet.Size(); i++)
	  if (!fine_facet[i]) order_facet[i] = 0;

      }

    // dof tables
    int ndof_lo = nfacets;
    ndof = ndof_lo;
    
    first_facet_dof.SetSize(nfacets+1);
    first_facet_dof = ndof_lo;

    for (size_t i = 0; i < nfacets; i++ )
      {
        first_facet_dof[i] = ndof;
        int inc = order_facet[i][0];
        if (inc > 0) ndof += inc;
      }
    first_facet_dof[nfacets] = ndof;
    
    if(print)
      {
	*testout << "*** Update NormalFacetSurfaceFESpace: General Information" << endl;
	*testout << " order facet (facet) " << order_facet << endl;
	*testout << " first_facet_dof (facet)  " << first_facet_dof << endl; 
      }

    UpdateCouplingDofArray();
  }

  
  void  NormalFacetSurfaceFESpace :: UpdateCouplingDofArray ()
  {
    ctofdof.SetSize(ndof);
    ctofdof = WIREBASKET_DOF;
    int first,next;
    for(int facet=0; facet<ma->GetNEdges(); facet++)
      {
        COUPLING_TYPE ct = fine_facet[facet] ? WIREBASKET_DOF : UNUSED_DOF;
        ctofdof[facet] = ct; // low_order

	first = first_facet_dof[facet];
	next = first_facet_dof[facet+1];
	for(int j=first ; j<next; j++)
	  ctofdof[j] = INTERFACE_DOF;
      }
  }

  FiniteElement & NormalFacetSurfaceFESpace :: GetFE ( ElementId ei, Allocator & lh ) const
  {
    if (!DefinedOn (ei))
      {
        return
          SwitchET (ma->GetElType(ei), [&] (auto et) -> FiniteElement&
                      {
                        return *new (lh) DummyFE<et.ElementType()> ();
                      });
      }
    switch(ei.VB())
      {
      case VOL:
        {
          throw Exception("No volume elements available for NormalFacetSurfaceFESpace. Trace is missing");
        }
      case BND:
        {
          return * SwitchET<ET_TRIG,ET_QUAD>
            (ma->GetElType(ei), [&] (auto et) -> FiniteElement*
            {
              using ET_T = ET_trait<et.ElementType()>;
              auto fe = new (lh) NormalFacetVolumeFE<et.ElementType()>();
              ArrayMem<int,ET_T::N_EDGE> fanums, order_fa;
              auto vnums = ma->GetElVertices(ei);
              fanums = ma->GetElEdges(ei);
         
              order_fa.SetSize(fanums.Size());
              for (auto i : Range(fanums.Size()))
                order_fa[i] = order_facet[fanums[i]][0];
              fe->SetVertexNumbers(vnums);
              fe->SetOrder(order_fa);
              return fe;
            });
        }
      case BBND: case BBBND: default:
        throw Exception ("NormalFacetSurfaceFESpace::GetFE does not support BBND or BBBND");
      }
  }

  void NormalFacetSurfaceFESpace :: GetDofNrs(ElementId ei, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if (!DefinedOn (ei)) return;
    
    if(ei.VB()==BND)
      {
        for (auto facet : ma->GetElEdges(ei))
          {
            dnums += facet;
            dnums += GetFacetDofs(facet);
          }
      }          
    
      
    // *testout << "dnums = " << endl << dnums << endl;
  }

  void NormalFacetSurfaceFESpace :: SetOrder (NodeId ni, int order) 
  {
    if (order_policy == CONSTANT_ORDER || order_policy == NODE_TYPE_ORDER)
      throw Exception("In NormalFacetSurfaceFESpace::SetOrder. Order policy is constant or node-type!");
    else if (order_policy == OLDSTYLE_ORDER)
      order_policy = VARIABLE_ORDER;
      
    if (order < 0)
      order = 0;
    
    if (CoDimension(ni.GetType(), ma->GetDimension()) == 1)
      if (ni.GetNr() < order_facet.Size())
	order_facet[ni.GetNr()] = fine_facet[ni.GetNr()] ? order : 0;
  }
  
  int NormalFacetSurfaceFESpace :: GetOrder (NodeId ni) const
  {
    if (CoDimension(ni.GetType(), ma->GetDimension()) == 1)
      if (ni.GetNr() < order_facet.Size())
	return order_facet[ni.GetNr()][0];
     
    return 0;
  }

  

  void NormalFacetSurfaceFESpace :: GetInnerDofNrs ( int felnr, Array<int> & dnums ) const
  {
    dnums.SetSize0();
  }

  int NormalFacetSurfaceFESpace :: GetNFacetDofs ( int felnr ) const
  {
    // number of low_order_dofs = dimension - 1
    return ( first_facet_dof[felnr+1] - first_facet_dof[felnr] + 1);
  }

  
  ///
  INT<2> NormalFacetSurfaceFESpace :: GetFacetOrder(int fnr) const
  { return order_facet[fnr]; };
  
  int NormalFacetSurfaceFESpace :: GetFirstFacetDof(int fanr) const {return (first_facet_dof[fanr]);}; 


  void NormalFacetSurfaceFESpace :: GetVertexDofNrs ( int elnum, Array<int> & dnums ) const
  {
    dnums.SetSize0();
  }

  void NormalFacetSurfaceFESpace :: GetEdgeDofNrs ( int elnum, Array<int> & dnums ) const
  {
    dnums.SetSize0();

    dnums.Append(elnum);
    for (int j=first_facet_dof[elnum]; j<first_facet_dof[elnum+1]; j++)
      dnums.Append(j);
  }

  void NormalFacetSurfaceFESpace :: GetFaceDofNrs (int felnr, Array<int> & dnums) const
  {
    dnums.SetSize0();
  }


  static RegisterFESpace<NormalFacetSurfaceFESpace> init_nfacetsurf ("normalfacetsurface");

}
