/*********************************************************************/
/* File:   tangentialfacetfespace.cpp                                */
/* Author: A. Sinwel, Joachim Schoeberl                              */
/* Date:   2008                                                      */
/*********************************************************************/

#include "tangentialfacetfespace.hpp"

#include "../fem/tangentialfacetfe.hpp"
#include "../fem/hcurllofe.hpp"
#include <../fem/hcurl_equations.hpp>
#include <../fem/diffop_impl.hpp>
#include <prolongation.hpp>


namespace ngcomp
{
  // shape functions live on the facets of the volume element only
  template <int D>
  class DiffOpIdTangentialFacet : public DiffOpIdEdge<D>
  {
  public:
    static constexpr bool ONLY_ON_FACETS = true;
  };




  using ngmg::FacetProlongation;
  
  // 1D facets of 2D elements
  // prolongate only on coarsegrid edges
  class TangentialFacetHOProlongation2D : public FacetProlongation<2>
  {
    int order;

  public:
    TangentialFacetHOProlongation2D(shared_ptr<MeshAccess> ama, int aorder)
      : FacetProlongation<2>(ama), order(aorder) { }

    void CalcMatrices()
    {
      if (haveprols) return;
      haveprols=true;
      for (int classnr = 0; classnr < 2; classnr++)
        {
          std::array<size_t,3> verts{0,1,2};
          if (classnr == 1) swap(verts[0], verts[1]);

          size_t vertsc[2] = { verts[0], verts[1] };
          size_t vertsf[2] = { verts[0], verts[2] };

          TangentialFacetFacetFE<ET_SEGM> felc(order);
          felc.SetVertexNumbers (vertsc);

          TangentialFacetFacetFE<ET_SEGM> felf(order);
          felf.SetVertexNumbers (vertsf);
          
          IntegrationRule ir(ET_SEGM, 2*order);
          size_t ndof = felf.GetNDof();
          dofs_per_facet = ndof;

          Matrix massf(ndof, ndof), massfc(ndof, ndof);          
          Matrix shapef(ndof,1), shapec(ndof,1), shapeftrans(ndof, 1);

          massf = 0.;
          massfc = 0.;
          Mat<1,1> trafo;
          for (IntegrationPoint ip : ir)
            {
              IntegrationPoint ipc(0.5*(1+ip(0)));
              trafo(0,0) = 2;

              felc.CalcShape (ipc, shapec);
              felf.CalcShape (ip, shapef);
              shapeftrans = shapef * trafo;
              massf += ip.Weight() * shapeftrans * Trans(shapeftrans);
              massfc += ip.Weight() * shapeftrans * Trans(shapec);
            }

          CalcInverse (massf);
          facetprol[classnr].SetSize(ndof, ndof);
          facetprol[classnr] = massf * massfc;
          // *testout << "prol " << classnr << " = " << endl << facetprol[classnr] << endl;
        }
    }

    virtual void Update (const FESpace & bfes) override
    {
      auto & fes = dynamic_cast<const TangentialFacetFESpace&> (bfes);
      
      size_t oldnfacet = prol_class.Size();
      size_t nfacet = ma->GetNFacets();
      cout << IM(3) << "update prol, level = " << ma->GetNLevels() <<  ", nfacet = " << nfacet << endl;

      while (facets_on_level.Size() < ma->GetNLevels())
        facets_on_level.Append(oldnfacet);
      facets_on_level[ma->GetNLevels()-1] = nfacet;

      // *testout << "update first_dofs, first_dofs.Size() == " << first_dofs.Size() << ", levels = " << ma->GetNLevels() << endl;
      if (first_dofs.Size() < ma->GetNLevels())
        first_dofs += fes.GetFirstFacetDof();
      
      // *testout << /* IM(3) << */ "facets_on_level = " << endl << facets_on_level << endl;

      if (ma->GetNLevels() == 1)
        {
          prol_class.SetSize(nfacet);
          prol_class = 0;
          return;
        }
      
      CalcMatrices();

      prol_class.SetSize(nfacet);

      // Array<size_t> verts(4);
      for (size_t i = oldnfacet; i < nfacet; i++)
        if (auto parents = get<1>(ma->GetParentEdges(i)); parents[1] == -1)
          {
            auto myverts = ma->GetEdgePNums(i);
            auto paverts = ma->GetEdgePNums(parents[0]);
            auto minpa = Min(paverts);
            prol_class[i] = myverts.Contains(minpa) ? 0 : 1;
          }
    }
  };



  // 1D facets of 2D elements
  // prolongate only on coarsegrid edges
  class TangentialFacetHOProlongation3D : public FacetProlongation<3>
  {
    int order;

  private:
    int GetClassNr (FlatArray<size_t> verts)
    {
      Array<size_t> hverts{verts};
      int classnr = 0;
      // sorting network for 4 elements, 5 comparisons:
      int pairs[][2] = { {0,1}, {2,3}, {0,2}, {1,3}, {1,2} };
      for (int i = 0; i < 5; i++)
         if (hverts[pairs[i][0]] > hverts[pairs[i][1]])
           {
             classnr += 1 << i;
             swap(hverts[pairs[i][0]], hverts[pairs[i][1]]);
           }
      return classnr;
    }
    Array<size_t> GetClassRealization (int classnr)
    {
      Array<size_t> verts{1,2,3,4};
      int pairs[][2] = { {0,1}, {2,3}, {0,2}, {1,3}, {1,2} };
      for (int i = 4; i >= 0; i--)
        if (classnr & (1 << i))
          swap(verts[pairs[i][0]], verts[pairs[i][1]]);
      return verts;
    }

    
  public:
    TangentialFacetHOProlongation3D(shared_ptr<MeshAccess> ama, int aorder)
      : FacetProlongation<3>(ama), order(aorder) { }

    void CalcMatrices()
    {
      if (haveprols) return;
      haveprols=true;

      for (int classnr = 0; classnr < 32; classnr++)
        {
          Array<size_t> verts{GetClassRealization(classnr)};
          
          size_t vertsc[3] = { verts[0], verts[1], verts[2] };
          size_t vertsf[3] = { verts[3], verts[1], verts[2] };
          
          TangentialFacetFacetFE<ET_TRIG> felc(order);
          felc.SetVertexNumbers (vertsc);

          TangentialFacetFacetFE<ET_TRIG> felf(order);
          felf.SetVertexNumbers (vertsf);
          
          IntegrationRule ir(ET_TRIG, 2*order);
          size_t ndof = felf.GetNDof();
          dofs_per_facet = ndof;
          
          Matrix massf(ndof, ndof), massfc(ndof, ndof);
          Matrix shapef(ndof,2), shapeftrans(ndof, 2), shapec(ndof,2);
          
          massf = 0.;
          massfc = 0.;
          Mat<2,2> trafo { { 2, 0 }, { 0, 1 } };
          for (IntegrationPoint ip : ir)
            {
              IntegrationPoint ipc(0.5*ip(0), ip(1));

              felc.CalcShape (ipc, shapec);
              felf.CalcShape (ip, shapef);
              shapeftrans = shapef * trafo;              
              massf += ip.Weight() * shapeftrans * Trans(shapeftrans);
              massfc += ip.Weight() * shapeftrans * Trans(shapec);
            }
          CalcInverse (massf);
          facetprol[classnr].SetSize(ndof, ndof);
          facetprol[classnr] = massf * massfc;
        }
    }

    virtual void Update (const FESpace & bfes) override
    {
      auto & fes = dynamic_cast<const TangentialFacetFESpace&> (bfes);
      
      size_t oldnfacet = prol_class.Size();
      size_t nfacet = ma->GetNFacets();
      cout << IM(3) << "update prol, level = " << ma->GetNLevels() <<  ", nfacet = " << nfacet << endl;

      while (facets_on_level.Size() < ma->GetNLevels())
        facets_on_level.Append(oldnfacet);
      facets_on_level[ma->GetNLevels()-1] = nfacet;

      if (first_dofs.Size() < ma->GetNLevels())
        first_dofs += fes.GetFirstFacetDof();
      
      cout << IM(3) << "facets_on_level = " << endl << facets_on_level << endl;

      if (ma->GetNLevels() == 1)
        {
          prol_class.SetSize(nfacet);
          prol_class = 0;
          return;
        }
      
      CalcMatrices();

      prol_class.SetSize(nfacet);
      
      for (size_t i = oldnfacet; i < nfacet; i++)
        if (auto parents = get<1>(ma->GetParentFaces(i)); parents[1] == -1)
          {
            auto myverts = ma->GetFacePNums(i);
            auto paverts = ma->GetFacePNums(parents[0]);
            
            std::array<size_t,4> verts;
            
            for (int i = 0; i < 3; i++)
              if (!myverts.Contains(paverts[i]))
                verts[0] = paverts[i];
            
            for (int i = 0; i < 3; i++)
              if (!paverts.Contains(myverts[i]))
                verts[3] = myverts[i];
            
            auto between = ma->GetParentNodes(verts[3]);
            verts[2] = between[0]+between[1]-verts[0];
            verts[1] = paverts[0]+paverts[1]+paverts[2] - verts[0]-verts[2];
            
            prol_class[i] = GetClassNr(verts);
          }
    }
  };
  
  

  TangentialFacetFESpace :: TangentialFacetFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
					    bool parseflags )
    : FESpace(ama, flags )
  {
    type = "tangentialfacet";
    name = "TangentialFacetFESpace";
    if (ma->GetDimension() == 1)
      throw Exception("TangentialFacetFESpace is not defined on 1D meshes (a point has no tangent space)");
    DefineNumFlag("relorder");
    DefineNumFlag("hpslope");
    DefineDefineFlag("variableorder");

    if (parseflags) CheckFlags(flags);

    print = flags.GetDefineFlag("print");
    needs_transform_vec = false;

    // ndlevel.SetSize(0);

    all_dofs_together = !flags.GetDefineFlagX ("all_dofs_together").IsFalse();
    if (flags.GetDefineFlag("hoprolongation"))
      all_dofs_together = true;
    
    
    Flags loflags;
    loflags.SetFlag("order", 0.0);
    if (IsComplex()) loflags.SetFlag("complex");

    loflags.SetFlag ("low_order");
    if (!flags.GetDefineFlag ("low_order") && !all_dofs_together)
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

    if (InitBaseOrders (flags, var_order, rel_order))
      {
        var_order = false;
        rel_order = order-1;
      }


    if (ma->GetDimension() == 2)
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<2>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdTangentialFacet<2>>>();
      }
    else
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<3>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdTangentialFacet<3>>>();
      }

    static ConstantCoefficientFunction one(1);
    // Array<CoefficientFunction*> coeffs(1);
    // coeffs[0] = &one;
    // evaluator = GetIntegrators().CreateBFI("massvectorfacet", 2, coeffs);
    // integrator[BND] = GetIntegrators().CreateBFI("robinvectorfacet", ma->GetDimension(), &one); 

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
        additional_evaluators.Set ("tangentialcomponent", make_shared<T_DifferentialOperator<DiffOpTangentialComponentHCurl<2>>> ());        
        break;
      case 3:
        additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHCurlDual<3>>> ());
        additional_evaluators.Set ("tangentialcomponent", make_shared<T_DifferentialOperator<DiffOpTangentialComponentHCurl<3>>> ());                
        break;
      default:
        break;
      }


    if (flags.GetDefineFlag("hoprolongation"))
      {
        int order_continuous = order;
        if (highest_order_dc) order_continuous--;
        if (ma->GetDimension()==2)
          {
            prol = make_shared<TangentialFacetHOProlongation2D> (GetMeshAccess(), order_continuous);
            ma->GetNetgenMesh()->GetTopology().EnableTable("parentedges", true);
          }
        if (ma->GetDimension()==3)
          {
            prol = make_shared<TangentialFacetHOProlongation3D> (GetMeshAccess(), order_continuous);
            ma->GetNetgenMesh()->GetTopology().EnableTable("parentfaces", true);
          }
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
    docu.Arg("orderlink") = "FESpace\n"
      "  Take the orders of regions and elements from this space, offset by relorder.\n"
      "  Only the given space sets these orders.";
    docu.Arg("relorder") = "int\n"
      "  With orderlink: order relative to the linked space.\n"
      "  Without: order relative to the element orders of the mesh (as variableorder).";
    docu.Arg("hpslope") = "double = 0\n"
      "  On meshes refined by RefineHP: the order decreases by hpslope per level\n"
      "  of geometric refinement.";
    return docu;
  }

  void TangentialFacetFESpace :: Update()
  {
    bool base_changed = order_spec.base->timestamp > order_base_resolved;
    FESpace::Update();
    if ( print ) 
      *testout << "TangentialFacetFESpace, order " << order << endl 
	       << "rel_order " << rel_order << ", var_order " << var_order << endl;

    bool first_update = base_changed || GetTimeStamp() < ma->GetTimeStamp();
    if (first_update) timestamp = NGS_Object::GetNextTimeStamp();
    
    if ( low_order_space ) 
      low_order_space -> Update();

    size_t nel = ma->GetNE();
    size_t nfacets = ma->GetNFacets();

    if (first_update)
      {
	order_facet.SetSize(nfacets);
	fine_facet.SetSize(nfacets);
	fine_facet = false;

	for (Ngs_Element el : ma->Elements<BND>())
	  if (DefinedOn(el))
	    fine_facet[el.Facets()] = true;
	for (Ngs_Element el : ma->Elements<VOL>())
	  if (DefinedOn(el))
	    fine_facet[el.Facets()] = true;
	ma->AllReduceNodalData (StdNodeType(NT_FACET, ma->GetDimension()),
				fine_facet, NG_MPI_LOR);

	ResolveFacetOrders (fine_facet, order_facet, 0);

	if (highest_order_dc)
	  for (auto i : Range(order_facet))
	    if (fine_facet[i] && order_facet[i][0] != order)
	      throw Exception ("TangentialFacetFESpace: highest_order_dc needs the same order on all facets");
      }

    // dof tables
    // ncfacets = 0;
    int ndof_lo = (ma->GetDimension() == 2) ? nfacets : 2*nfacets;
    if (all_dofs_together) ndof_lo = 0;
    size_t ndof = ndof_lo;

    
    first_facet_dof.SetSize(nfacets+1);
    first_facet_dof = ndof_lo;

    if ( ma->GetDimension() == 2 )
      {
	for (size_t i = 0; i < nfacets; i++ )
	  {
	    first_facet_dof[i] = ndof;
            int inc = order_facet[i][0];
            if (all_dofs_together) inc++;
	    if (highest_order_dc) inc--;
            if (!fine_facet[i]) inc = 0;
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
	    IVec<2> p = order_facet[i];
	    if (highest_order_dc)
	      { p[0]--; p[1]--; }

	    auto pnums = ma->GetFacePNums(i);

	    inci = ( pnums.Size() == 3 ) ?
	      ( ((p[0]+1)*(p[0]+2)) - 2) :  ( 2 * (p[0]+1) * (p[1]+1) - 2 );
            if (all_dofs_together) inci+=2;
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

    SetNDof (ndof);

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

    for(int facet=0; facet < ma->GetNFacets(); facet++)
      {
        COUPLING_TYPE ct = fine_facet[facet] ? WIREBASKET_DOF : UNUSED_DOF;

        if (!all_dofs_together)
          {
            if ( ma->GetDimension() == 2 )
              ctofdof[facet] = ct; // low_order
            else
              {
                ctofdof[2*facet] = ct;
                ctofdof[2*facet+1] = ct;
              }

            IntRange r = GetFacetDofs(facet);
            ctofdof[r] = INTERFACE_DOF;
          }
        else
          {
            IntRange r = GetFacetDofs(facet);
            ctofdof[r] = INTERFACE_DOF;
            
            if (r.Size())
              {
                if ( ma->GetDimension() == 2 )
                  ctofdof[r.First()] = ct; // low_order
                else
                  {
                    ctofdof[r.First()] = ct;
                    ctofdof[r.First()+1] = ct;
                  }
              }
          }
            
      }
    if (highest_order_dc)	  
      for(int el=0; el<ma->GetNE(); el++)	      
	{
	  for (int k = first_inner_dof[el]; k < first_inner_dof[el+1]; k++)
	    ctofdof[k] = hide_highest_order_dc ? HIDDEN_DOF : LOCAL_DOF;
	}	  
    // *testout << " VECTORFACETFESPACE - ctofdof = \n" << ctofdof << endl;
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
              /*
              if(ET_T::DIM==2)
                fanums = ma->GetElEdges(ei);
              else
                fanums = ma->GetElFaces(ei);
              */
              fanums = ma->GetElFacets(ei);
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
              fe->SetVertexNumbers(vnums);
              if (et.ElementType() == ET_SEGM)
                {
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
            for (auto f : ma->GetElFacets(ei))
              {
                if (!all_dofs_together)
                  {
                    if ( ma->GetDimension() == 2 )
                      dnums.Append(f); // low_order
                    else
                      {
                        dnums.Append(2*f);
                        dnums.Append(2*f+1);
                      }
                  }
                
                dnums += GetFacetDofs(f);
              }
              
          }
        else
          {
            if (ma->GetDimension() == 2)
              {
                auto fanums = ma->GetElEdges (ei);
	    
                int innerdof = first_inner_dof[ei.Nr()];
                for(int i=0; i<fanums.Size(); i++)
                  {
                    int facetdof = first_facet_dof[fanums[i]];
                    for (int j = 0; j <= order; j++)
                      {
                        if (j == 0 && !all_dofs_together)
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
                              if (j+k == 0 && !all_dofs_together)
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
                              if (j+k == 0 && !all_dofs_together)
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

        if ( ma->GetDimension() == 2 )
          {
            auto fanums = ma->GetElEdges (ei);
            if (!all_dofs_together)
              dnums.Append(fanums[0]);
            
            /*
              first = first_facet_dof[fanums[0]];
              next = first_facet_dof[fanums[0]+1];
              for ( int j = first; j < next; j++ )
              dnums.Append(j);
            */
            dnums += GetFacetDofs(fanums[0]);        
          } 
        else // 3D
          {
            fanums[0] = ma->GetSElFace(ei.Nr());
            if (!all_dofs_together)
              {
                dnums.Append( 2*fanums[0] );
                dnums.Append( 2*fanums[0]+1 );
              }
            
            /*
              first = first_facet_dof[fanums[0]];
              next = first_facet_dof[fanums[0]+1];
              for ( int j = first; j < next; j++ )
              dnums.Append(j);
            */
            dnums += GetFacetDofs(fanums[0]);        
          }
      }
    //if(ei.VB()==BBND)
    //  dnums.SetSize(0);
  }

  void TangentialFacetFESpace :: SetOrder (NodeId ni, int order) 
  {
    order_policy = VARIABLE_ORDER;
      
    if (order < 0)
      order = 0;
    
    // kept in the specification, and set in the resolved orders
    int dim = ma->GetDimension();
    if (CoDimension(ni.GetType(), dim) == 1)
      if (ni.GetNr() < order_facet.Size())
        {
          if (dim >= 2 && fine_facet[ni.GetNr()])
            SetNodeOrder (dim-2, ni.GetNr(), order);
          order_facet[ni.GetNr()] = fine_facet[ni.GetNr()] ? order : 0;
        }
    if (!order_stale) order_stale = 1;
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


  shared_ptr<Array<int>> TangentialFacetFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    return nullptr;
  }
  
  void TangentialFacetFESpace :: GetFacetDofNrs ( int felnr, Array<int> & dnums ) const
  {
    dnums.SetSize0();
    if (!all_dofs_together)
      {
        if ( ma->GetDimension() == 3 )
          {
            dnums.Append (2*felnr );
            dnums.Append (2*felnr+1);
          }
        else
          dnums.Append ( felnr );
      }

    dnums += GetFacetDofs(felnr);
    /*
    int first = first_facet_dof[felnr];
    int next = first_facet_dof[felnr+1];
    for (int j = first; j < next; j++ )
      dnums.Append(j);
    */
  }


  void TangentialFacetFESpace :: GetInnerDofNrs ( int felnr, Array<int> & dnums ) const
  {
    dnums.SetSize0();
    if (highest_order_dc)
      dnums += IntRange(first_inner_dof[felnr], first_inner_dof[felnr+1]);
  }

  int TangentialFacetFESpace :: GetNFacetDofs ( int felnr ) const
  {
    // number of low_order_dofs = dimension - 1
    if (all_dofs_together)
      return first_facet_dof[felnr+1] - first_facet_dof[felnr];
    else
      return ( first_facet_dof[felnr+1] - first_facet_dof[felnr] + dimension - 1);
  }

  /*
  void TangentialFacetFESpace :: GetVertexNumbers(int elnr, Array<int>& vnums) const
  { 
    vnums = ma->GetElVertices(ElementId(VOL,elnr));
  };
  */
  ///
  IVec<2> TangentialFacetFESpace :: GetFacetOrder(int fnr) const
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

    if (!all_dofs_together)
      dnums.Append(elnum);
    //for (int j=first_facet_dof[elnum]; j<first_facet_dof[elnum+1]; j++)
    // dnums.Append(j);
    dnums += GetFacetDofs(elnum);
  }

  void TangentialFacetFESpace :: GetFaceDofNrs (int felnr, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if ( ma->GetDimension() == 2 ) return;

    if (!all_dofs_together)
      {
        dnums.Append(2*felnr);
        dnums.Append(2*felnr+1);
      }
    // for (int j=first_facet_dof[felnr]; j<first_facet_dof[felnr+1]; j++)
    // dnums.Append(j);
    dnums += GetFacetDofs(felnr);    
  }

  static RegisterFESpace<TangentialFacetFESpace> init_vfacet ("tangentialfacet");
}



