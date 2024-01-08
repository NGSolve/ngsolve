/*********************************************************************/
/* File:   tangentialfacetfespace.cpp                                */
/* Author: A. Sinwel, Joachim Schoeberl                              */
/* Date:   2008                                                      */
/*********************************************************************/

#include "tangentialfacetfespace.hpp"

#include "../fem/tangentialfacetfe.hpp"
#include "../fem/hcurllofe.hpp"
#include <../fem/hcurl_equations.hpp>


namespace ngcomp
{

  TangentialFacetFESpace :: TangentialFacetFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
					    bool parseflags )
    : FESpace(ama, flags )
  {
    type = "tangentialfacet";
    name = "TangentialFacetFESpace";
    DefineNumFlag("relorder");
    DefineDefineFlag("variableorder");

    if (parseflags) CheckFlags(flags);

    print = flags.GetDefineFlag("print");

    ndlevel.SetSize(0);
    Flags loflags;
    loflags.SetFlag("order", 0.0);
    if (IsComplex()) loflags.SetFlag("complex");

    loflags.SetFlag ("low_order");
    if (!flags.GetDefineFlag ("low_order"))
      low_order_space = make_shared<TangentialFacetFESpace>(ma, loflags);

    order = int (flags.GetNumFlag ("order",0)); 

    if(flags.NumFlagDefined("relorder") && !flags.NumFlagDefined("order")) 
      var_order = 1; 
    else 
      var_order = 0;
    
    rel_order=int(flags.GetNumFlag("relorder",order-1)); 


    if(flags.NumFlagDefined("order") && flags.NumFlagDefined("relorder")) 
      {
	if(var_order)
	  cerr << " WARNING: TangentialFacetFESpace: inconsistent flags: variableorder, order and relorder "
	       << "-> variable order space with rel_order " << rel_order << "is used, but order is ignored " << endl; 
	else 
	  cerr << " WARNING: TangentialFacetFESpace: inconsistent flags: order and rel_order "
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


    if (ma->GetDimension() == 2)
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<2>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<2>>>();
      }
    else
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<3>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<3>>>();
      }

    static ConstantCoefficientFunction one(1);
    // Array<CoefficientFunction*> coeffs(1);
    // coeffs[0] = &one;
    // evaluator = GetIntegrators().CreateBFI("massvectorfacet", 2, coeffs);
    integrator[BND] = GetIntegrators().CreateBFI("robinvectorfacet", ma->GetDimension(), &one); 

    highest_order_dc = flags.GetDefineFlag("highest_order_dc");
    if (highest_order_dc) {
      *testout << "highest_order_dc is active!" << endl;
    }
    hide_highest_order_dc = flags.GetDefineFlag("hide_highest_order_dc");
    
    switch (ma->GetDimension())
      {
      case 1:
        break;
      case 2:
        additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHCurlDual<2>>> ());
        break;
      case 3:
        additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHCurlDual<3>>> ());
        break;
      default:
        break;
      }
    
    // Update();
  }
  
  DocInfo TangentialFacetFESpace :: GetDocu()
  {
    DocInfo docu = FESpace::GetDocu();
    docu.Arg("highest_order_dc")  = "bool = False\n"
      "  Splits highest order facet functions into two which are associated with\n  the corresponding neighbors and are local dofs on the corresponding element\n (used to realize projected jumps)";
    docu.Arg("hide_highest_order_dc") = "bool = False\n"
      "  if highest_order_dc is used this flag marks the corresponding local dofs\n  as hidden dofs (reduces number of non-zero entries in a matrix). These dofs\n  can also be compressed.";
    return docu;
  }

  void TangentialFacetFESpace :: Update()
  {
    FESpace::Update();
    if ( print ) 
      *testout << "TangentialFacetFESpace, order " << order << endl 
	       << "rel_order " << rel_order << ", var_order " << var_order << endl;

    bool first_update = GetTimeStamp() < ma->GetTimeStamp();
    if (first_update) timestamp = NGS_Object::GetNextTimeStamp();
    
    if ( low_order_space ) 
      low_order_space -> Update();

    size_t nel = ma->GetNE();
    size_t nfacets = ma->GetNFacets();

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
	    fine_facet[el.Facets()] = true;
#ifdef PARALLEL
	if(var_order)
	  throw Exception("MPI + variable order for TangentialFacetFESpace is not implemented.");
#endif


	for (size_t i = 0; i < nel; i++)
	  {
	    ElementId ei(VOL,i);
	    if (!DefinedOn (ei))
	      continue;
	    INT<3> el_orders = ma->GetElOrders(i); 

	    ELEMENT_TYPE eltype=ma->GetElType(ei); 
	    const POINT3D * points = ElementTopology :: GetVertices (eltype);	
	
	    if (ma->GetDimension() == 2)
	      {
		auto fanums = ma->GetElEdges (ei);
		for (int j=0;j<fanums.Size();j++) 
		  fine_facet[fanums[j]] = 1; 
	    
		if(var_order)
		  {
		    const EDGE * edges = ElementTopology::GetEdges (eltype);
		    for(int j=0; j<fanums.Size(); j++)
		      for(int k=0;k<2;k++)
			if(points[edges[j][0]][k] != points[edges[j][1]][k])
			  { 
			    order_facet[fanums[j]] = INT<2>(max2(el_orders[k]+rel_order, order_facet[fanums[j]][0]),0);
			    break; 
			  }
		  }
	      }
	    else
	      {
		auto elfaces = ma->GetElFaces(ei);
		for (int j=0;j<elfaces.Size();j++) fine_facet[elfaces[j]] = 1; 
	    
		if(var_order) 
		  {
		    auto vnums = ma->GetElVertices (ei);
		    const FACE * faces = ElementTopology::GetFaces (eltype);
		    for(int j=0;j<elfaces.Size();j++)
		      {
			if(faces[j][3]==-1) // trig  
			  {
			    order_facet[elfaces[j]][0] = max2(order_facet[elfaces[j]][0],el_orders[0]+rel_order);
			    order_facet[elfaces[j]][1] = order_facet[elfaces[j]][0]; 
			  }
			else //quad_face
			  {
			    int fmax = 0;
			    for(int k = 1; k < 4; k++) 
			      if(vnums[faces[j][k]] > vnums[faces[j][fmax]]) fmax = k;   
					
			    INT<2> f((fmax+3)%4,(fmax+1)%4); 
			    if(vnums[faces[j][f[1]]] > vnums[faces[j][f[0]]]) swap(f[0],f[1]);
			
			    // fmax > f[0] > f[1]
			    // p[0] for direction fmax,f[0] 
			    // p[1] for direction fmax,f[1] 
			    for(int l=0;l<2;l++)
			      for(int k=0;k<3;k++)
				if(points[faces[j][fmax]][k] != points[faces[j][f[l] ]][k])
				  {
				    order_facet[elfaces[j]][l] = max2(order_facet[elfaces[j]][l], rel_order + el_orders[k]);
				    break; 
				  } 
			
			  }
		      }
		  }
	    
	      }
	
	  }

	for (int i = 0; i < fine_facet.Size(); i++)
	  if (!fine_facet[i]) order_facet[i] = 0;

	ma->AllReduceNodalData ((ma->GetDimension()==2) ? NT_EDGE : NT_FACE, 
				fine_facet, MPI_LOR);
      }

    // dof tables
    // ncfacets = 0;
    int ndof_lo = (ma->GetDimension() == 2) ? nfacets : 2*nfacets;
    ndof = ndof_lo;
    
    first_facet_dof.SetSize(nfacets+1);
    first_facet_dof = ndof_lo;

    if ( ma->GetDimension() == 2 )
      {
	for (size_t i = 0; i < nfacets; i++ )
	  {
	    first_facet_dof[i] = ndof;
            int inc = order_facet[i][0];
	    if (highest_order_dc) inc--;
            if (inc > 0) ndof += inc;
	  }
	first_facet_dof[nfacets] = ndof;
	
	if (highest_order_dc)
	  {
	    size_t ne = ma->GetNE();
	    first_inner_dof.SetSize(ne+1);
	    for (int i = 0; i < ne; i++)
	      {
		first_inner_dof[i] = ndof;
		
		// only trigs supported:
		auto fnums = ma->GetElFacets(ElementId(VOL,i));
		ndof += fnums.Size();
	      }
	    first_inner_dof[ne] = ndof;
	  }
      }

    else // 3D
      {
	int inci = 0;
	for (size_t i = 0; i < nfacets; i++)
	  {
	    INT<2> p = order_facet[i];
	    if (highest_order_dc)
	      { p[0]--; p[1]--; }

	    auto pnums = ma->GetFacePNums(i);

	    inci = ( pnums.Size() == 3 ) ?
	      ( ((p[0]+1)*(p[0]+2)) - 2) :  ( 2 * (p[0]+1) * (p[1]+1) - 2 );

	    first_facet_dof[i] = ndof;
            if (fine_facet[i])
              ndof += inci;
	  }
	first_facet_dof[nfacets] = ndof;


	if (highest_order_dc)
	  {
	    size_t ne = ma->GetNE();
	    first_inner_dof.SetSize(ne+1);
	    for (size_t i = 0; i < ne; i++)
	      {
		first_inner_dof[i] = ndof;
		switch (ma->GetElType(ElementId(VOL,i)))
		  {
		  case ET_TET: ndof += 4*(order+1)*2; break;
		  case ET_PRISM: ndof += 2 * (2*(order+1)+3*(2*order+1)); break;
		  case ET_HEX: ndof += 6*(2*order+1)*2; break;
		  default: throw Exception (string("TangentialFacetFESpace: Element type not implemented"));		  
		  }
	      }
	    first_inner_dof[ne] = ndof;
	  }
      }

    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;
      
    //no prolongation so far       
    //prol->Update();

    if(print)
      {
	*testout << "*** Update TangentialFacetFESpace: General Information" << endl;
	*testout << " order facet (facet) " << order_facet << endl;
	*testout << " first_facet_dof (facet)  " << first_facet_dof << endl; 
	*testout << " first_inner_dof (facet)  " << first_inner_dof << endl; 
      }
    
    UpdateCouplingDofArray();
  }

  
  void  TangentialFacetFESpace :: UpdateCouplingDofArray ()
  {
    ctofdof.SetSize(ndof);
    ctofdof = WIREBASKET_DOF;
    int first,next;
    for(int facet=0; facet<ma->GetNFacets(); facet++)
      {
        COUPLING_TYPE ct = fine_facet[facet] ? WIREBASKET_DOF : UNUSED_DOF;
	if ( ma->GetDimension() == 2 )
	  ctofdof[facet] = ct; // low_order
	else
	  {
	    ctofdof[2*facet] = ct;
	    ctofdof[2*facet+1] = ct;
	  }
	
	first = first_facet_dof[facet];
	next = first_facet_dof[facet+1];
	for(int j=first ; j<next; j++)
	  ctofdof[j] = INTERFACE_DOF;
      }
    if (highest_order_dc)	  
      for(int el=0; el<ma->GetNE(); el++)	      
	{
	  for (int k = first_inner_dof[el]; k < first_inner_dof[el+1]; k++)
	    ctofdof[k] = hide_highest_order_dc ? HIDDEN_DOF : LOCAL_DOF;
	}	  
    *testout << " VECTORFACETFESPACE - ctofdof = \n" << ctofdof << endl;
  }

  FiniteElement & TangentialFacetFESpace :: GetFE ( ElementId ei, Allocator & lh ) const
  {
    if (!DefinedOn (ei))
      {
        return
          SwitchET (ma->GetElType(ei), [&] (auto et) -> FiniteElement&
                      {
                        return *new (lh) HCurlDummyFE<et.ElementType()> ();
                      });
      }
    
    switch(ei.VB())
      {
      case VOL:
        {
          return * SwitchET<ET_TRIG,ET_QUAD,ET_TET,ET_PRISM,ET_PYRAMID,ET_HEX>
            (ma->GetElType(ei), [&] (auto et) -> FiniteElement*
            {
              using ET_T = ET_trait<et.ElementType()>;
              auto fe = new (lh) TangentialFacetVolumeFE<et.ElementType()>();
              // ArrayMem<int,ET_T::N_VERTEX> vnums;
              ArrayMem<int,ET_T::N_FACET> fanums, order_fa;
              auto vnums = ma->GetElVertices(ei);
              if(ET_T::DIM==2)
                fanums = ma->GetElEdges(ei);
              else
                fanums = ma->GetElFaces(ei);
              assert(vnums.Size() == ET_T::N_VERTEX);
              assert(fanums.Size() == ET_T::N_FACET);
              order_fa.SetSize(fanums.Size());
              for (auto i : Range(fanums.Size()))
                order_fa[i] = order_facet[fanums[i]][0];
              fe->SetVertexNumbers(vnums);
              fe->SetOrder(order_fa);
              fe->SetHighestOrderDC(highest_order_dc);
              return fe;
            });
        }

      case BND:
        {
          int reduceorder = highest_order_dc ? 1 : 0;
          return *SwitchET<ET_SEGM,ET_TRIG,ET_QUAD>
            (ma->GetElType(ei), [&] (auto et) -> FiniteElement*
            {
              using ET_T = ET_trait<et.ElementType()>;
              auto fe = new (lh) TangentialFacetFacetFE<et.ElementType()>();
              //ArrayMem<int,ET_T::N_VERTEX> vnums;
              ArrayMem<int,ET_T::N_EDGE> ednums;
              auto vnums = ma->GetElVertices(ei);
              assert(vnums.Size() == ET_T::N_VERTEX);
              fe->SetVertexNumbers(vnums);
              if (et.ElementType() == ET_SEGM)
                {
                  // ArrayMem<int,ET_T::N_EDGE> ednums;
                  auto ednums = ma->GetElEdges(ei);
                  fe->SetOrder(order_facet[ednums[0]][0]-reduceorder);
                }
              else
                {
                  fe->SetOrder(order_facet[ma->GetSElFace(ei.Nr())][0]-reduceorder);
                }
              return fe;
            });
        }
      case BBND: case BBBND: default:
        throw Exception ("TangentialFacetFESpace::GetFE does not support BBND or BBBND");
      }
  }

  void TangentialFacetFESpace :: GetDofNrs(ElementId ei, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if (!DefinedOn (ei)) return;
    
    if(ei.VB()==VOL)
      {
    if (!highest_order_dc)
      {
	Array<int> fanums; // facet numbers
	int first,next;
	
	fanums.SetSize(0);
	
	
	if(ma->GetDimension() == 3)
	  fanums = ma->GetElFaces (ei);
	else // dim=2
	  fanums = ma->GetElEdges (ei);
	
	for(int i=0; i<fanums.Size(); i++)
	  {
	    if ( ma->GetDimension() == 2 )
	      dnums.Append(fanums[i]); // low_order
	    else
	      {
		dnums.Append(2*fanums[i]);
		dnums.Append(2*fanums[i]+1);
	      }
	    
	    first = first_facet_dof[fanums[i]];
	    next = first_facet_dof[fanums[i]+1];
	    for(int j=first ; j<next; j++)
	      dnums.Append(j);
	  }
      }
    else
      {
	if (ma->GetDimension() == 2)
	  {
	    // Array<int> fanums; // facet numbers
	    
	    // fanums.SetSize(0);
            
	    auto fanums = ma->GetElEdges (ei);
	    
	    int innerdof = first_inner_dof[ei.Nr()];
	    for(int i=0; i<fanums.Size(); i++)
	      {
		int facetdof = first_facet_dof[fanums[i]];
		for (int j = 0; j <= order; j++)
		  {
		    if (j == 0)
		      {
			dnums.Append(fanums[i]);
		      }
		    else if (j < order)
		      {
			dnums.Append (facetdof++);
		      }
		    else
		      {
			dnums.Append (innerdof++);
		      }
		  }
	      }
	  }
	else
	  {
	    // Array<int> fanums; // facet numbers
	    
	    // fanums.SetSize(0);
	    
	    ELEMENT_TYPE et = ma->GetElType (ei);
	    auto fanums = ma->GetElFaces (ei);
	    
	    int innerdof = first_inner_dof[ei.Nr()];
	    for(int i=0; i<fanums.Size(); i++)
	      {
		ELEMENT_TYPE ft = ElementTopology::GetFacetType (et, i);
		
		int facetdof = first_facet_dof[fanums[i]];
		
		if (ft == ET_TRIG)
		  {
		    for (int j = 0; j <= order; j++)
		      for (int k = 0; k <= order-j; k++)
			{
			  if (j+k == 0)
			    {
			      dnums.Append(2*fanums[i]);
			      dnums.Append(2*fanums[i]+1);
			    }
			  else if (j+k < order)
			    {
			      dnums.Append (facetdof++);
			      dnums.Append (facetdof++);
			    }
			  else
			    {
			      dnums.Append (innerdof++);
			      dnums.Append (innerdof++);
			    }
			}
		  }
		else
		  {
		    for (int j = 0; j <= order; j++)
		      for (int k = 0; k <= order; k++)
			{
			  if (j+k == 0)
			    {
			      dnums.Append(2*fanums[i]);
			      dnums.Append(2*fanums[i]+1);
			    }
			  else if ( (j < order) && (k < order) )
			    {
			      dnums.Append (facetdof++);
			      dnums.Append (facetdof++);
			    }
			  else
			    {
			      dnums.Append (innerdof++);
			      dnums.Append (innerdof++);
			    }
			}
		  }
	      }
	  }
      }
      
    // *testout << "dnums = " << endl << dnums << endl;
      }
    if(ei.VB()==BND)
      {
    ArrayMem<int, 1> fanums(1);
    int first, next;


    if ( ma->GetDimension() == 2 )
      {
	auto fanums = ma->GetElEdges (ei);
	dnums.Append(fanums[0]);

	first = first_facet_dof[fanums[0]];
	next = first_facet_dof[fanums[0]+1];
	for ( int j = first; j < next; j++ )
	  dnums.Append(j);
      } 
    else // 3D
      {
	fanums[0] = ma->GetSElFace(ei.Nr());
	dnums.Append( 2*fanums[0] );
	dnums.Append( 2*fanums[0]+1 );

	first = first_facet_dof[fanums[0]];
	next = first_facet_dof[fanums[0]+1];
	for ( int j = first; j < next; j++ )
	  dnums.Append(j);

      }

      }
    //if(ei.VB()==BBND)
    //  dnums.SetSize(0);
  }

  void TangentialFacetFESpace :: SetOrder (NodeId ni, int order) 
  {
    if (order_policy == CONSTANT_ORDER || order_policy == NODE_TYPE_ORDER)
      throw Exception("In TangentialFacetFESpace::SetOrder. Order policy is constant or node-type!");
    else if (order_policy == OLDSTYLE_ORDER)
      order_policy = VARIABLE_ORDER;
      
    if (order < 0)
      order = 0;
    
    if (CoDimension(ni.GetType(), ma->GetDimension()) == 1)
      if (ni.GetNr() < order_facet.Size())
	order_facet[ni.GetNr()] = fine_facet[ni.GetNr()] ? order : 0;
  }
  
  int TangentialFacetFESpace :: GetOrder (NodeId ni) const
  {
    if (CoDimension(ni.GetType(), ma->GetDimension()) == 1)
      if (ni.GetNr() < order_facet.Size())
	return order_facet[ni.GetNr()][0];
     
    return 0;
  }


  FlatArray<VorB> TangentialFacetFESpace :: GetDualShapeNodes (VorB vb) const
  {
    static VorB nodes[] = { BND, VOL };
    if (int(vb) > 1)
      { return FlatArray<VorB> (0, nullptr); }
    else
      { return FlatArray<VorB> (1, &nodes[int(vb)]); }
  }


  shared_ptr<Table<int>> TangentialFacetFESpace :: CreateSmoothingBlocks (const Flags & precflags) const
  { 
    return nullptr;
  }

  shared_ptr<Array<int>> TangentialFacetFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    return nullptr;
  }
  
  void TangentialFacetFESpace :: GetFacetDofNrs ( int felnr, Array<int> & dnums ) const
  {
    dnums.SetSize0();
    if ( ma->GetDimension() == 3 )
      {
	dnums.Append( 2*felnr );
	dnums.Append (2*felnr+1);
      }
    else
      dnums.Append ( felnr );

    int first = first_facet_dof[felnr];
    int next = first_facet_dof[felnr+1];
    for (int j = first; j < next; j++ )
      dnums.Append(j);
  }


  void TangentialFacetFESpace :: GetInnerDofNrs ( int felnr, Array<int> & dnums ) const
  {
    dnums.SetSize0();
  }

  int TangentialFacetFESpace :: GetNFacetDofs ( int felnr ) const
  {
    // number of low_order_dofs = dimension - 1
    return ( first_facet_dof[felnr+1] - first_facet_dof[felnr] + dimension - 1);
  }

  
  void TangentialFacetFESpace :: GetVertexNumbers(int elnr, Array<int>& vnums) const
  { 
    vnums = ma->GetElVertices(ElementId(VOL,elnr));
  };

  ///
  INT<2> TangentialFacetFESpace :: GetFacetOrder(int fnr) const
  { return order_facet[fnr]; };
  
  int TangentialFacetFESpace :: GetFirstFacetDof(int fanr) const {return (first_facet_dof[fanr]);}; 


  void TangentialFacetFESpace :: GetVertexDofNrs ( int elnum, Array<int> & dnums ) const
  {
    dnums.SetSize0();
  }

  void TangentialFacetFESpace :: GetEdgeDofNrs ( int elnum, Array<int> & dnums ) const
  {
    dnums.SetSize0();
    if ( ma->GetDimension() == 3 )
      return;

    dnums.Append(elnum);
    for (int j=first_facet_dof[elnum]; j<first_facet_dof[elnum+1]; j++)
      dnums.Append(j);
  }

  void TangentialFacetFESpace :: GetFaceDofNrs (int felnr, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if ( ma->GetDimension() == 2 ) return;

    dnums.Append( 2*felnr);
    dnums.Append (2*felnr+1);
    for (int j=first_facet_dof[felnr]; j<first_facet_dof[felnr+1]; j++)
      dnums.Append(j);
  }

  static RegisterFESpace<TangentialFacetFESpace> init_vfacet ("tangentialfacet");
}



