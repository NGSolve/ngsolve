/*********************************************************************/
/* File:   vectorfacetfespace.cpp                                    */
/* Author: A. Sinwel, Joachim Schoeberl                              */
/* Date:   2008                                                      */
/*********************************************************************/

#include <comp.hpp>
#include <fem.hpp>
#include <cassert>

namespace ngcomp
{

  VectorFacetFESpace :: VectorFacetFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
					    bool parseflags )
    : FESpace(ama, flags )
  {
    name = "VectorFESpace";
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
      low_order_space = make_shared<VectorFacetFESpace>(ma, loflags);

    order = int (flags.GetNumFlag ("order",0)); 

    if(flags.NumFlagDefined("relorder") && !flags.NumFlagDefined("order")) 
      var_order = 1; 
    else 
      var_order = 0;
    
    rel_order=int(flags.GetNumFlag("relorder",order-1)); 


    if(flags.NumFlagDefined("order") && flags.NumFlagDefined("relorder")) 
      {
	if(var_order)
	  cerr << " WARNING: VectorFacetFESpace: inconsistent flags: variableorder, order and relorder "
	       << "-> variable order space with rel_order " << rel_order << "is used, but order is ignored " << endl; 
	else 
	  cerr << " WARNING: VectorFacetFESpace: inconsistent flags: order and rel_order "
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

    // Update();
  }
  


  void VectorFacetFESpace :: Update(LocalHeap& lh)
  {
    FESpace::Update(lh);
    if ( print ) 
      *testout << "VectorFacetFESpace, order " << order << endl 
	       << "rel_order " << rel_order << ", var_order " << var_order << endl;

    if ( low_order_space ) 
      low_order_space -> Update(lh);

    int nel = ma->GetNE();
    int nfacets = ma->GetNFacets();
 
    int p = 0; 
    if (!var_order) p = order; 
    
    order_facet.SetSize(nfacets);
    fine_facet.SetSize(nfacets);

    order_facet = p;
    fine_facet = 0; 
    
    Array<int> fanums;
        
    for (int i = 0; i < nel; i++)
      {
	INT<3> el_orders = ma->GetElOrders(i); 

	ELEMENT_TYPE eltype=ma->GetElType(i); 
	const POINT3D * points = ElementTopology :: GetVertices (eltype);	
	
	if (ma->GetDimension() == 2)
	  {
	    ma->GetElEdges (i, fanums);
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
	    Array<int> elfaces,vnums;
	    ma->GetElFaces(i,elfaces);
	    for (int j=0;j<elfaces.Size();j++) fine_facet[elfaces[j]] = 1; 
	    
	    if(var_order) 
	      {
		ma->GetElVertices (i, vnums);
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


    // update dofs

    // ncfacets = 0;
    int ndof_lo = (ma->GetDimension() == 2) ? nfacets : 2*nfacets;
    ndof = ndof_lo;

    first_facet_dof.SetSize(nfacets+1);
    first_facet_dof = ndof_lo;

    if ( ma->GetDimension() == 2 )
      {
	for ( int i = 0; i < nfacets; i++ )
	  {
	    first_facet_dof[i] = ndof;
            int inc = order_facet[i][0];
	    if (highest_order_dc) inc--;
            if (inc > 0) ndof += inc;
	  }
	first_facet_dof[nfacets] = ndof;
	
	if (highest_order_dc)
	  {
	    int ne = ma->GetNE();
	    first_inner_dof.SetSize(ne+1);
	    Array<int> fnums;
	    for (int i = 0; i < ne; i++)
	      {
		first_inner_dof[i] = ndof;
		
		// only trigs supported:
		ma->GetElFacets(i,fnums);
		ndof += fnums.Size();
	      }
	    first_inner_dof[ne] = ndof;
	  }
      }

    else // 3D
      {
	int inci = 0;
	Array<int> pnums;

	for (int i = 0; i < nfacets; i++)
	  {
	    INT<2> p = order_facet[i];
	    if (highest_order_dc)
	      { p[0]--; p[1]--; }

	    ma->GetFacePNums(i, pnums);

	    inci = ( pnums.Size() == 3 ) ?
	      ( ((p[0]+1)*(p[0]+2)) - 2) :  ( 2 * (p[0]+1) * (p[1]+1) - 2 );

	    first_facet_dof[i] = ndof;
	    ndof += inci;
	  }
	first_facet_dof[nfacets] = ndof;


	if (highest_order_dc)
	  {
	    int ne = ma->GetNE();
	    first_inner_dof.SetSize(ne+1);
	    for (int i = 0; i < ne; i++)
	      {
		first_inner_dof[i] = ndof;
		switch (ma->GetElType(i))
		  {
		  case ET_TET: ndof += 4*(order+1)*2; break;
		  case ET_PRISM: ndof += 2 * (2*(order+1)+3*(2*order+1)); break;
		  default: throw Exception (string("VectorFacetFESpace: Element type not implemented"));
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
	*testout << "*** Update VectorFacetFESpace: General Information" << endl;
	*testout << " order facet (facet) " << order_facet << endl; 
	*testout << " first_facet_dof (facet)  " << first_facet_dof << endl; 
	*testout << " first_inner_dof (facet)  " << first_inner_dof << endl; 
      } 
    
    UpdateCouplingDofArray();
  }


  void  VectorFacetFESpace :: UpdateCouplingDofArray ()
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
	    ctofdof[k] = LOCAL_DOF;
	}	  
    *testout << " VECTORFACETFESPACE - ctofdof = \n" << ctofdof << endl;
  }

  FiniteElement & VectorFacetFESpace :: GetFE ( ElementId ei, Allocator & lh ) const
  {
    switch(ei.VB())
      {
      case VOL:
        {
          return * SwitchET<ET_TRIG,ET_QUAD,ET_TET,ET_PRISM,ET_PYRAMID,ET_HEX>
            (ma->GetElType(ei), [&] (auto et) -> FiniteElement*
            {
              using ET_T = ET_trait<et.ElementType()>;
              auto fe = new (lh) VectorFacetVolumeFE<et.ElementType()>();
              ArrayMem<int,ET_T::N_VERTEX> vnums;
              ArrayMem<int,ET_T::N_FACET> fanums, order_fa;
              ma->GetElVertices(ei,vnums);
              if(ET_T::DIM==2)
                ma->GetElEdges(ei, fanums);
              else
                ma->GetElFaces(ei.Nr(), fanums);
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
            (ma->GetElType(ei), [&] (auto ET) -> FiniteElement*
            {
              using ET_T = ET_trait<ET>;
              auto fe = new (lh) VectorFacetFacetFE<ET>();
              ArrayMem<int,ET_T::N_VERTEX> vnums;
              ArrayMem<int,ET_T::N_EDGE> ednums;
              ma->GetElVertices(ei,vnums);
              assert(vnums.Size() == ET_T::N_VERTEX);
              fe->SetVertexNumbers(vnums);
              if (ET == ET_SEGM)
                {
                  ArrayMem<int,ET_T::N_EDGE> ednums;
                  ma->GetElEdges(ei,ednums);
                  fe->SetOrder(order_facet[ednums[0]][0]-reduceorder);
                }
              else
                {
                  fe->SetOrder(order_facet[ma->GetSElFace(ei.Nr())][0]-reduceorder);
                }
              return fe;
            });
        }
      case BBND:
        throw Exception ("VectorFacetFESpace::GetFE doesnt support BBND");
      }
  }
  // const FiniteElement & VectorFacetFESpace :: GetFE ( int elnr, LocalHeap & lh ) const
  // {
  //   VectorFacetVolumeFiniteElement<2> * fe2d = 0;
  //   VectorFacetVolumeFiniteElement<3> * fe3d = 0;
  //   switch (ma->GetElType(elnr))
  //     {
  //     case ET_TRIG: fe2d = new (lh)  VectorFacetVolumeTrig (); break;
  //     case ET_QUAD: fe2d = new (lh) VectorFacetVolumeQuad (); break;
  //     case ET_TET:  fe3d = new (lh) VectorFacetVolumeTet (); break;
  //     case ET_PYRAMID: fe3d = new (lh)  VectorFacetVolumePyramid (); break;
  //     case ET_PRISM: fe3d = new (lh)  VectorFacetVolumePrism (); break;
  //     case ET_HEX:   fe3d = new (lh) VectorFacetVolumeHex (); break;
  //     default:
  //       throw Exception (string("VectorFacetFESpace::GetFE: unsupported element ")+
  //                        ElementTopology::GetElementName(ma->GetElType(elnr)));
  //     }

  //   Array<int> vnums;
  //   ArrayMem<int, 6> fanums, order_fa;
    
  //   ma->GetElVertices(elnr, vnums);

  //   if (fe2d)
  //     {
  //       ma->GetElEdges(elnr, fanums);
    
  //       order_fa.SetSize(fanums.Size());
  //       for (int j = 0; j < fanums.Size(); j++)
  //         order_fa[j] = order_facet[fanums[j]][0]; 
	
  //       fe2d -> SetVertexNumbers (vnums);
  //       fe2d -> SetOrder (order_fa);
  //       fe2d -> ComputeNDof();
  //       fe2d -> SetHighestOrderDC(highest_order_dc);
  //       return *fe2d;
  //     }    
  //   else
  //     {
  //       ma->GetElFaces(elnr, fanums);
    
  //       order_fa.SetSize(fanums.Size());
  //       for (int j = 0; j < fanums.Size(); j++)
  //         order_fa[j] = order_facet[fanums[j]][0]; 
	
  //       fe3d -> SetVertexNumbers (vnums);
  //       fe3d -> SetOrder (order_fa);
  //       fe3d -> ComputeNDof();
  //       fe3d -> SetHighestOrderDC(highest_order_dc);
    
  //       return *fe3d;
  //     }
  // }


  // const FiniteElement & VectorFacetFESpace :: GetSFE ( int selnr, LocalHeap & lh ) const
  // {
  //   VectorFacetFacetFiniteElement<1> * fe1d = 0;
  //   VectorFacetFacetFiniteElement<2> * fe2d = 0;

  //   switch (ma->GetSElType(selnr))
  //     {
  //     case ET_SEGM: fe1d = new (lh) VectorFacetFacetSegm (); break;
  //     case ET_TRIG: fe2d = new (lh) VectorFacetFacetTrig (); break;
  //     case ET_QUAD: fe2d = new (lh) VectorFacetFacetQuad (); break;
  //     default:
  //       ;
  //     }
     
  //   if (!fe1d && !fe2d)
  //     {
  //       stringstream str;
  //       str << "VectorFacetFESpace " << GetClassName()
  //           << ", undefined eltype "
  //           << ElementTopology::GetElementName(ma->GetSElType(selnr))
  //           << ", order = " << order << endl;
  //       throw Exception (str.str());
  //     }
     
  //   ArrayMem<int,4> vnums;
  //   ArrayMem<int, 4> ednums;
    
  //   ma->GetSElVertices(selnr, vnums);
  //   int reduceorder = highest_order_dc ? 1 : 0;
  //   switch (ma->GetSElType(selnr))
  //     {
  //     case ET_SEGM:
  //       fe1d -> SetVertexNumbers (vnums);
  //       ma->GetSElEdges(selnr, ednums);
  //       fe1d -> SetOrder (order_facet[ednums[0]][0]-reduceorder); 
  //       fe1d -> ComputeNDof();
  //       return *fe1d;
  //       break;
  //     case ET_TRIG: 
  //     case ET_QUAD:
  //       fe2d -> SetVertexNumbers (vnums);
  //       fe2d -> SetOrder (order_facet[ma->GetSElFace(selnr)][0]-reduceorder);
  //       fe2d -> ComputeNDof();
  //       return *fe2d;
  //       break;
  //     default:
  //       throw Exception ("VectorFacetFESpace::GetSFE: unsupported element");
  //     }
  // }

  void VectorFacetFESpace :: GetDofNrs(ElementId ei, Array<int> & dnums) const
  {
    if(ei.VB()==VOL)
      {
    if (!highest_order_dc)
      {
	Array<int> fanums; // facet numbers
	int first,next;
	
	fanums.SetSize(0);
	dnums.SetSize(0);
	
	
	if(ma->GetDimension() == 3)
	  ma->GetElFaces (ei.Nr(), fanums);
	else // dim=2
	  ma->GetElEdges (ei.Nr(), fanums);
	
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
	    Array<int> fanums; // facet numbers
	    
	    fanums.SetSize(0);
	    dnums.SetSize(0);
	    
	    ma->GetElEdges (ei.Nr(), fanums);
	    
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
	    Array<int> fanums; // facet numbers
	    
	    fanums.SetSize(0);
	    dnums.SetSize(0);
	    
	    ELEMENT_TYPE et = ma->GetElType (ei.Nr());
	    ma->GetElFaces (ei.Nr(), fanums);
	    
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
      
    if (!DefinedOn (VOL, ma->GetElIndex (ei.Nr())))
      dnums = -1;
    // *testout << "dnums = " << endl << dnums << endl;
      }
    if(ei.VB()==BND)
      {
    dnums.SetSize(0);
    ArrayMem<int, 1> fanums(1);
    int first, next;


    if ( ma->GetDimension() == 2 )
      {
	ma->GetSElEdges ( ei.Nr(), fanums );
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
    if(ei.VB()==BBND)
      dnums.SetSize(0);
  }

  shared_ptr<Table<int>> VectorFacetFESpace :: CreateSmoothingBlocks (const Flags & precflags) const
  { 
    return nullptr;
  }

  Array<int> * VectorFacetFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    return NULL;
  }
  
  void VectorFacetFESpace :: GetFacetDofNrs ( int felnr, Array<int> & dnums ) const
  {
    dnums.SetSize(0);
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


  void VectorFacetFESpace :: GetInnerDofNrs ( int felnr, Array<int> & dnums ) const
  {
    dnums.SetSize(0);
  }

  int VectorFacetFESpace :: GetNFacetDofs ( int felnr ) const
  {
    // number of low_order_dofs = dimension - 1
    return ( first_facet_dof[felnr+1] - first_facet_dof[felnr] + dimension - 1);
  }

  
  void VectorFacetFESpace :: GetVertexNumbers(int elnr, Array<int>& vnums) const
  { 
    ma->GetElVertices(elnr, vnums); 
  };

  ///
  INT<2> VectorFacetFESpace :: GetFacetOrder(int fnr) const
  { return order_facet[fnr]; };
  
  int VectorFacetFESpace :: GetFirstFacetDof(int fanr) const {return (first_facet_dof[fanr]);}; 


  void VectorFacetFESpace :: GetVertexDofNrs ( int elnum, Array<int> & dnums ) const
  {
    dnums.SetSize(0);
  }

  void VectorFacetFESpace :: GetEdgeDofNrs ( int elnum, Array<int> & dnums ) const
  {
    dnums.SetSize(0);
    if ( ma->GetDimension() == 3 )
      return;

    dnums.Append(elnum);
    for (int j=first_facet_dof[elnum]; j<first_facet_dof[elnum+1]; j++)
      dnums.Append(j);
  }

  void VectorFacetFESpace :: GetFaceDofNrs (int felnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    if ( ma->GetDimension() == 2 ) return;

    dnums.Append( 2*felnr);
    dnums.Append (2*felnr+1);
    for (int j=first_facet_dof[felnr]; j<first_facet_dof[felnr+1]; j++)
      dnums.Append(j);
  }


  static RegisterFESpace<VectorFacetFESpace> init_vfacet ("vectorfacet");

  /*
  // register FESpaces
  namespace vectorfacetfespace_cpp
  {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("vectorfacet", VectorFacetFESpace::Create);
    }
    
    Init init;
  }
  */
}



