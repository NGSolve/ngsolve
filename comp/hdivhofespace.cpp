/**
   High Order Finite Element Space for H(Div)
*/
/* Continuous and discontinuous version : 
   Disont space has no extra low-order block(!), all dofs internal! 
   Flags:
   -order  : uniform order 
   -relorder : variable order relative to mesh-order 
   (on facets maximum order of facet-patch elements, also for discont)
   -curl_order,relcurlorder: orders of curl fields in uniform/variable version 
   -disontinuous (DefineFlag) : discontinuous space 
   -orderinner, orderfacet: for variable and uniform space uniform order for inner or facets 
	 
*/ 



#include <comp.hpp>
#include <../fem/hdivlofe.hpp>  
#include <../fem/hdivhofe.hpp>  
#include <../fem/hdivhofefo.hpp>  


namespace ngcomp
{

  HDivHighOrderFESpace ::  
  HDivHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    type = "hdivho";
    name="HDivHighOrderFESpace(hdivho)";
    // allowed flags
    DefineNumFlag("relorder");
    // DefineNumFlag("relcurlorder");
    DefineDefineFlag("discontinuous");
    // DefineNumFlag("curlorder");
    DefineNumFlag("orderinner");
    DefineNumFlag("orderedge");
    DefineNumFlag("orderface");
    DefineNumFlag("orderfacet");
    DefineDefineFlag("hdiv");
    DefineDefineFlag("hdivho");
    // DefineDefineFlag("print");
    DefineDefineFlag("noprint");
    DefineDefineFlag("variableorder"); 
    DefineDefineFlag("hodivfree"); 
    
    if(parseflags) CheckFlags(flags);

    discont = flags.GetDefineFlag("discontinuous"); 


    Flags loflags;
    loflags.SetFlag ("order", 1);
    loflags.SetFlag ("dim", dimension);
    if (iscomplex) loflags.SetFlag ("complex");
    if (discont) loflags.SetFlag ("discontinuous"); // supported ?

    // low_order_space = new RaviartThomasFESpace (ma, loflags);
    low_order_space = 0; 



    // Variable order space: 
    //      in case of (var_order && order) or (relorder) 
    var_order = flags.GetDefineFlag("variableorder");  
    order =  int (flags.GetNumFlag ("order",0)); 
    curl_order =  int (flags.GetNumFlag ("curlorder",order)); 

    
    if(flags.NumFlagDefined("relorder") && !flags.NumFlagDefined("order")) 
      var_order = 1; 
    
    rel_order=int(flags.GetNumFlag("relorder",order-1)); 
    

    if(flags.NumFlagDefined("order") && flags.NumFlagDefined("relorder")) 
      {
	if(var_order)
	  cerr << " WARNING: HDivHoFeSpace: inconsistent flags: variableorder, order and relorder "
	       << "-> variable order space with rel_order " << rel_order << "is used, but order is ignored " << endl; 
	else 
	  cerr << " WARNING: HDivHoFeSpace: inconsistent flags: order and rel_order "
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

    // SZ hack since is not supported(tested) yet 
    rel_curl_order= rel_order; 
    curl_order = order;
        
        
    print = flags.GetDefineFlag("print"); 

    ho_div_free = flags.GetDefineFlag("hodivfree"); 
    fixed_order = flags.GetDefineFlag ("fixedorder");

    uniform_order_inner = int (flags.GetNumFlag ("orderinner", -1));
    
    if(flags.NumFlagDefined("orderedge") || flags.NumFlagDefined ("orderface")) 
      throw Exception("Flag 'orderface' and 'orderedge' for hdivho are obsolete. Use flag 'orderfacet' instead!"); 
	
    uniform_order_facet = int (flags.GetNumFlag ("orderfacet", -1));


    auto one = make_shared<ConstantCoefficientFunction> (1);
    // integrator[VOL]= GetIntegrators().CreateBFI("masshdiv", ma->GetDimension(), one);
    // integrator[BND] = GetIntegrators().CreateBFI("robinhdiv", ma->GetDimension(), one);
    
    if (ma->GetDimension() == 2)
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<2>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdVecHDivBoundary<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<2>>>();
      }
    else
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<3>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdVecHDivBoundary<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<3>>>();
      }

    /*
    if (dimension > 1)
      {
        integrator[VOL]= make_shared<BlockBilinearFormIntegrator> (integrator[VOL], dimension);
        integrator[BND] = make_shared<BlockBilinearFormIntegrator> (integrator[BND], dimension);
      }
    */
    
    highest_order_dc = flags.GetDefineFlag("highest_order_dc");
    if (highest_order_dc) {
      *testout << "highest_order_dc is active!" << endl;
    }
  }
  
  HDivHighOrderFESpace:: ~HDivHighOrderFESpace () 
  {
    ;
  }


  void HDivHighOrderFESpace :: Update(LocalHeap & lh)
  {
    FESpace::Update(lh);
    // In discontinuous spaces the order on edges and faces  
    // are also set to the maximal order of adjacent elements  
    // and not to the element order (Motivation: Equilibrated_EE) 

    // SZ hack since is not supported(tested) yet 
    rel_curl_order= rel_order; 
    curl_order = order; 

    if (low_order_space)
      low_order_space -> Update(lh);
    
    // int nv = ma->GetNV();
    size_t nel = ma->GetNE();
    size_t nfa = ma->GetNFacets();
    size_t dim = ma->GetDimension();
       
    order_facet.SetSize(nfa);
    order_inner.SetSize(nel); 
    order_inner_curl.SetSize(nel); 
    fine_facet.SetSize(nfa); 
    boundary_facet.SetSize(nfa); 
   
    boundary_facet = false;
    /*
    for (int i = 0; i < ma->GetNSE(); i++)
      {
	Array<int> elfacets; 
	ma->GetSElFacets (i,elfacets); 
	boundary_facet[elfacets[0]] = true;
      }
    */
    for (auto el : ma->Elements(BND))
      boundary_facet[el.Facets()] = true;

    // cout << " order hdiv " << order << endl; 
    // cout << " curl_order hdiv " << curl_order << endl; 

    int p = 0, pc = 0; 
    
    if(!var_order)
      {
	p = order; 
	pc = curl_order; 
      } 
    
    order_facet = pc;
    order_inner = p;
    order_inner_curl = pc;
    fine_facet = 0; //!!!! 


    for (auto el : ma->Elements(VOL))
      {
        if (!DefinedOn (el))
          {
            order_inner[el.Nr()] = 0;
            order_inner_curl[el.Nr()] = 0;
            continue;
          }
          
	ELEMENT_TYPE eltype = el.GetType();
	const POINT3D * points = ElementTopology :: GetVertices (eltype);
	
	// Array<int> elfacets; 
	// ma->GetElFacets (el.Nr(), elfacets); 
	// auto elfacets = ma->GetElFacets (el);
        auto elfacets = el.Facets();
        
        fine_facet[elfacets] = true;
	
	if(!var_order) continue; 
	
	INT<3> el_orders = ma->GetElOrders(el.Nr());  
	
        int i = el.Nr();
	for(int k=0;k<dim;k++)
	  {
	    order_inner_curl[i][k]= max2(el_orders[k] + rel_curl_order,0);
	    order_inner[i][k] = max2(el_orders[k]+rel_order,0);
	  }

	if(dim==2)
	  {
	    const EDGE * edges = ElementTopology::GetEdges (eltype);
	    for(int j=0; j<elfacets.Size(); j++)
	      for(int k=0;k<2;k++)
		if(points[edges[j][0]][k] != points[edges[j][1]][k])
		  { 
		    order_facet[elfacets[j]][0] = max2(el_orders[k]+rel_curl_order, order_facet[elfacets[j]][0]);
		    break; 
		  }
	  }
	else
	  {
	    // Array<int> vnums (el.Vertices());
            auto vnums = el.Vertices();
	    const FACE * faces = ElementTopology::GetFaces (eltype);

	    for(int j = 0; j < elfacets.Size(); j++)
	      {
		if(faces[j][3]==-1) // trig  
		  {
		    order_facet[elfacets[j]][0] = max2(order_facet[elfacets[j]][0],el_orders[0]+rel_curl_order);
		    order_facet[elfacets[j]][1] = order_facet[elfacets[j]][0]; 
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
                            order_facet[elfacets[j]][l] = max2(order_facet[elfacets[j]][l], rel_curl_order + 
                                                              el_orders[k]);
                            break; 
                          } 
		    
		  }
	      }
	  }
      }

    ma->AllReduceNodalData ((ma->GetDimension()==2) ? NT_EDGE : NT_FACE,
			   fine_facet, MPI_LOR);

    if(uniform_order_inner > -1) order_inner = uniform_order_inner;
    if(uniform_order_facet > -1) order_facet = uniform_order_facet;

    for (auto i : Range(nfa)) if (!fine_facet[i]) order_facet[i] = INT<2> (0,0); 

    // by SZ ... since order_inner_curl is not working yet for hdivhofe
    order_inner_curl = order_inner; 
    
    if(print) 
      {
	*testout << " discont " << discont << endl;
	*testout << " fine_facet[i] (hdiv) " << fine_facet << endl; 
	
	*testout << " order_facet (hdivho) " << order_facet << endl; 
	*testout << " order_inner (hdivho) " << order_inner << endl; 
	*testout << " order_inner_curl (hdivho) " << order_inner_curl << endl; 
      }
    
    UpdateDofTables(); 
    UpdateCouplingDofArray();
  }

  void HDivHighOrderFESpace :: UpdateDofTables()
  {
    size_t nfa = ma->GetNFacets();
    size_t nel = ma->GetNE();
    size_t dim = ma->GetDimension();
    // Array<int> pnums; 
     
    first_facet_dof.SetSize(nfa+1); 
    first_inner_dof.SetSize(nel+1); 

    ndof = nfa; 
    first_facet_dof = ndof; 

    if(dim==2)
      {
	// int dec_hodc = highest_order_dc ? 1 : 0;
        for (auto i : Range(nfa))
          {
            first_facet_dof[i] = ndof;
            int inc = fine_facet[i] ? order_facet[i][0] : 0;
	    if (highest_order_dc && !boundary_facet[i]) inc--;
            if (inc > 0) ndof += inc;

            /*
            if(fine_facet[i])
	      ndof += order_facet[i][0];
	    if (highest_order_dc && !boundary_facet[i]) ndof--;
            */
          }

        first_facet_dof[nfa] = ndof;
      
	// Array<int> fnums;
        // for (auto i : Range(nel))
        for (size_t i = 0; i < nel; i++)
          {
            ElementId ei(VOL, i);
            INT<3> pc = order_inner_curl[i];
            INT<3> p = order_inner[i];
            int inci = 0;
            switch(ma->GetElType(ei))
              {
              case ET_TRIG:
                if (!ho_div_free)
                  inci = pc[0]*(pc[0]-1)/2 + p[0]*(p[0]-1)/2 + p[0]-1;
                else
                  inci = pc[0]*(pc[0]-1)/2;
                break;
              case ET_QUAD:
                if (!ho_div_free)
                  inci = pc[0]*pc[1] + p[0]*p[1]+p[0]+p[1];
                else
                  inci = pc[0]*pc[1];
                break;
              default: // for the compiler
                break;  
              }

	    if (highest_order_dc)
	      {
                /*
		ma->GetElFacets (ei, fnums);
                for (auto f : fnums)
		  if (!boundary_facet[f]) inci++;
                */
                for (auto f : ma->GetElFacets(ei))
		  if (!boundary_facet[f]) inci++;
	      }

            first_inner_dof[i] = ndof;
            if (inci > 0) ndof+=inci;
          }
        first_inner_dof[nel] = ndof;


        if (highest_order_dc)
          {
            dc_pairs.SetSize (ma->GetNFacets());
            dc_pairs = INT<2> (-1,-1);
            
            // Array<int> fnums;
            for (auto ei : ma->Elements(VOL))
              {
                // auto i = ei.Nr();
                // ma->GetElFacets (ei, fnums);
                // auto fnums = ma->GetElFacets(ei);
		int fid = first_inner_dof[ei.Nr()];
                for (auto f : ma->GetElFacets(ei))
		  if (!boundary_facet[f])
		    {
		      int di = fid++; // first_inner_dof[i]+k;
		      dc_pairs[f][1] = dc_pairs[f][0];
		      dc_pairs[f][0] = di;
		    }
              }
          }
        else
          dc_pairs.SetSize0 ();
      }
    else 
      {
        int inci = 0;
        for (size_t i = 0; i < nfa; i++) 
          {
            inci = 0; 
            if(fine_facet[i])
              {
                INT<2> p = order_facet[i]; 
                // ma->GetFacePNums(i,pnums);
                auto pnums = ma->GetFacePNums(i);
                switch(pnums.Size())
                  {
                  case 3: //Triangle
                    inci= (p[0]*p[0]+3*p[0])/2; 
                    if (highest_order_dc && !boundary_facet[i]) inci -= p[0]+1;
                    break;
                  case 4: //Quad 
                    inci= p[0]*p[1] + p[0] + p[1];
                    if (highest_order_dc && !boundary_facet[i]) inci -= p[0]+p[1]+1;
                    break;
                  }
              }
            else
              inci = 0;

            if (inci < 0) inci = 0;
            first_facet_dof[i] = ndof;
            ndof+= inci;

          }
        first_facet_dof[nfa] = ndof;
	 
	// Array<int> fnums;
        for (size_t i = 0; i < nel; i++)
          {
            INT<3> p = order_inner[i];
            INT<3> pc = order_inner_curl[i];
            int inci = 0;
	     
            switch(ma->GetElType(ElementId(VOL,i)))
              {
              case ET_TET:
                if(p[0]>1 && !ho_div_free)
                  inci = p[0]*(p[0]+1)*(p[0]-1)/6 + p[0]*(p[0]-1)/2 + p[0]-1;
                if(pc[0]>1)
                  inci += pc[0]*(pc[0]+1)*(pc[0]-1)/3 + pc[0]*(pc[0]-1)/2; ;
                if (highest_order_dc) 
		  {
		    // ma->GetElFacets (i, fnums);
                    auto fnums = ma->GetElFacets(i);
		    for (int j = 0; j < fnums.Size(); j++)
		      if (!boundary_facet[fnums[j]]) inci += p[0]+1;
		  }
		// inci += 4*(p[0]+1);
                break;
              case ET_PRISM:
                // inci = (p[0]+1)*(3*(p[0]-1)+(p[0]-1)*(p[0]-2))+ (p[0]-1)*(p[0]+1)*(p[0]+2)/2;
		inci = (p[0]+2)*p[0]*(p[2]+1) + (p[0]+1)*(p[0]+2)*p[2]/2;
		if (ho_div_free)
		  inci -= (p[0]+1)*(p[0]+2)*(p[2]+1)/2 - 1;

                if (highest_order_dc) inci += 2*(p[0]+1)+3*(p[0]+p[2]+1);
                break;
              case ET_HEX:
                inci =  2*pc[0]*pc[1]*pc[2] + pc[0]*pc[1] + pc[1]*pc[2] + pc[0]* pc[2]
                  +  p[0]*(p[1]*p[2] + p[1] + p[2] + 1)  + p[1]*p[2] + p[1] + p[2]; 
                break; 
              case ET_PYRAMID: 
                inci=0; 
                cout << "WARNING: there are hdiv-pyramids (not implemented yet) !! " << endl;
                break; 
              default:
                inci = 0;
                break;
              }
            if (inci < 0) inci = 0;
            first_inner_dof[i] = ndof;
            ndof+= inci;
          }
        first_inner_dof[nel] = ndof;


        
        if (highest_order_dc)
          {
            dc_pairs.SetSize ((order+1)*ma->GetNFacets());
            dc_pairs = INT<2> (-1,-1);
            
            for (int i = 0; i < ma->GetNE(); i++)
              {
                ElementId ei(VOL,i);
                auto fnums = ma->GetElFacets (ei);
		int di = first_inner_dof[i]; // +k*(order+1);
		    
                for (int k = 0; k < fnums.Size(); k++)
		  if (!boundary_facet[fnums[k]])
		    {
		      int base = fnums[k]*(order+1);
		      
		      for (int l = 0; l < order+1; l++)
			{
			  dc_pairs[base+l][1] = dc_pairs[base+l][0];
			  dc_pairs[base+l][0] = di++;
			}
		    }
              }
          }
        else
          dc_pairs.SetSize (0);

      }
   

    if(discont) 
      { 
        ndof = 0; 
        for(int i=0;i<nel;i++)
          {
            ElementId ei(VOL,i);
            int incii = first_inner_dof[i+1]-first_inner_dof[i]; 
	     	     
            auto elfacets = ma->GetElFacets (ei);
	     
            for(int j=0; j<elfacets.Size(); j++) 
              incii+=first_facet_dof[elfacets[j]+1]-first_facet_dof[elfacets[j]]+1; // incl. lowest-order 
	     
            first_inner_dof[i] = ndof; 
            ndof += incii;
	    
          } 
        first_inner_dof[nel] = ndof; 
        first_facet_dof = 0; 
      } 
	
    
    if(print) 
      {
        (*testout) << "ndof hdivhofespace update = " << endl << ndof << endl;
        (*testout) << "first_facet_dof (hdiv)  = " << endl << first_facet_dof << endl;
        (*testout) << "first_facet_dof (hdiv) = " << endl << first_facet_dof << endl;
        (*testout) << "first_inner_dof (hdiv) = " << endl << first_inner_dof << endl;
      }
     
    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;
    //    prol->Update();
  }

  void HDivHighOrderFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(ndof);
    if(discont) 
      {
        ctofdof = LOCAL_DOF;
        return;
      } 
    
    ctofdof = WIREBASKET_DOF;
    
    for (auto facet : Range(ma->GetNFacets()))
      {
        ctofdof[facet] = fine_facet[facet] ? WIREBASKET_DOF : UNUSED_DOF;
        ctofdof[GetFacetDofs(facet)] = INTERFACE_DOF;
      }
    
    for (auto el : Range (ma->GetNE()))
      ctofdof[GetElementDofs(el)] = LOCAL_DOF;
  }


  template <ELEMENT_TYPE ET>
  FiniteElement & HDivHighOrderFESpace :: T_GetFE (int elnr, bool onlyhdiv, Allocator & lh) const
  {
    Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL> (elnr);
    if (!DefinedOn(ngel)) return * new (lh) HDivDummyFE<ET>();
    
    HDivHighOrderFE<ET> * hofe =  new (lh) HDivHighOrderFE<ET> ();

    hofe -> SetVertexNumbers (ngel.Vertices());
    hofe -> SetHODivFree (ho_div_free);
    hofe -> SetOnlyHODiv (onlyhdiv);

    hofe -> SetOrderInner (order_inner[elnr]);
        
    switch (int(ET_trait<ET>::DIM))
      {
      case 2:
        hofe -> SetOrderFacet (order_facet[ngel.Edges()]);
        break;
      case 3:
        hofe -> SetOrderFacet (order_facet[ngel.Faces()]);
        break;
      }
    hofe -> ComputeNDof();
    return *hofe;
  }

  FiniteElement & HDivHighOrderFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    if (ei.IsVolume())
      {
        int elnr = ei.Nr();
        Ngs_Element ngel = ma->GetElement(ei);
        ELEMENT_TYPE eltype = ngel.GetType();
        
        switch (eltype)
          {
            // case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, false, alloc);
            
          case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, false, alloc);
          case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, false, alloc);
            
          case ET_TET:     return T_GetFE<ET_TET> (elnr, false, alloc);
          case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, false, alloc);
            // case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, false, alloc);
          case ET_HEX:     return T_GetFE<ET_HEX> (elnr, false, alloc);
            
          default:
            throw Exception ("illegal element in HDivHOFESpace::GetFE");
          }
      }
    else if (ei.VB()  == BND)
      {
        if (!DefinedOn(ei))
          return SwitchET(ma->GetElType(ei), [&] (auto et) -> FiniteElement&
                          {
                            return * new (alloc) HDivNormalDummyFE<et.ElementType()>();
                          });
        
        int selnr = ei.Nr();
        FiniteElement * fe = 0;
        
        int porder; 
        if (discont) porder = -1; 
        else porder = order; 
        
        // if (highest_order_dc) porder--;
        
        auto vnums = ma->GetElVertices(ei);
        
        switch (ma->GetElType(ei))
          {
          case ET_SEGM:
            {
              auto fe1 = new (alloc) HDivHighOrderNormalSegm<TrigExtensionMonomial> (porder);
              fe1 -> SetVertexNumbers(vnums);
              fe = fe1;
              break;
            }
          case ET_TRIG:
            {
              auto fe1 = new (alloc) HDivHighOrderNormalTrig<TrigExtensionMonomial> (porder);
              fe1 -> SetVertexNumbers(vnums);              
              fe = fe1;
              break;
            }
          case ET_QUAD:
            {
              auto fe1 = new (alloc) HDivHighOrderNormalQuad<TrigExtensionMonomial> (porder);
              fe1 -> SetVertexNumbers(vnums);              
              fe = fe1;
              break;
            }
          default:
            throw Exception (string("HDivHighOrderFESpace::GetSFE: unsupported element ")+
                             ElementTopology::GetElementName(ma->GetElType(ei)));
          }
        
        if (discont) return *fe; 
        
        // ArrayMem<int, 4> ednums, order_ed;
        INT<3> order_fa;
        
        if(ma->GetElType(ei) == ET_SEGM)
          {
            HDivHighOrderNormalFiniteElement<1> * hofe =
              dynamic_cast<HDivHighOrderNormalFiniteElement<1>*> (fe);
            
            // hofe -> SetVertexNumbers (vnums);
            auto ednums = ma->GetElEdges(ei);
            // int dec = (!boundary_facet[ednums[0]] && highest_order_dc) ? 1 : 0;
            hofe -> SetOrderInner (order_facet[ednums[0]][0] /* -dec */);
            hofe -> ComputeNDof();
          }
        else
          {
            HDivHighOrderNormalFiniteElement<2> * hofe =
              dynamic_cast<HDivHighOrderNormalFiniteElement<2>*> (fe);
            
            // hofe -> SetVertexNumbers (vnums);
            
#ifdef NEW_HDIVFE
            INT<3> order_fa = INT<3>(order_facet[ma->GetSElFace(selnr)][0],
                                     order_facet[ma->GetSElFace(selnr)][1],0);
            if (highest_order_dc)
              {
                order_fa[0]--;
                order_fa[1]--;
                order_fa[2]--;
              }
            hofe -> SetOrderInner (order_fa);
#else 
            int order_fa = order_facet[ma->GetSElFace(selnr)][0];
            // if (highest_order_dc) order_fa--;
            hofe -> SetOrderInner (order_fa);
#endif
            hofe -> ComputeNDof();
          }
        
        return *fe;
      }
    else
      switch (ma->GetElement(ei).GetType())
        {
        case ET_POINT: return * new (alloc) DummyFE<ET_POINT>();
        case ET_SEGM: return * new (alloc) DummyFE<ET_SEGM>();
        default:
          __assume(false);
          throw Exception("HDiv - impossible element");
        }
  }
  
  // const FiniteElement & HDivHighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  // {
  //   Ngs_Element ngel = ma->GetElement(elnr);
  //   ELEMENT_TYPE eltype = ngel.GetType();

  //   /*
  //   if (ma->GetElType(elnr) == ET_TRIG && order <= 6 && fixed_order)
  //     {
  //       HDivHighOrderFiniteElementFO<2> * hofe2d = 0;
  //       switch (order)
  //         {
  //         case 1: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,1> (); break;
  //         case 2: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,2> (); break;
  //         case 3: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,3> (); break;
  //         case 4: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,4> (); break;
  //         case 5: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,5> (); break;
  //         case 6: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,6> (); break;
  //         }
	
  //       Ngs_Element ngel = ma->GetElement<2> (elnr);
  //       for (int j = 0; j < 3; j++)
  //         hofe2d->SetVertexNumber (j, ngel.vertices[j]);

  //       hofe2d -> SetHODivFree (ho_div_free);
  //       hofe2d -> SetOnlyHODiv (false);
  //       hofe2d -> ComputeNDof();
	
  //       return *hofe2d;
  //     }  
  //   */

  //   switch (eltype)
  //     {
  //       // case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, false, lh);
        
  //     case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, false, lh);
  //     case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, false, lh);
        
  //     case ET_TET:     return T_GetFE<ET_TET> (elnr, false, lh);
  //     case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, false, lh);
  //           // case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, false, lh);
  //           // case ET_HEX:     return T_GetFE<ET_HEX> (elnr, false, lh);
        
  //     default:
  //       throw Exception ("illegal element in HDivHOFESpace::GetFE");
  //     }
  // }

  const FiniteElement & HDivHighOrderFESpace :: GetHODivFE (int elnr, LocalHeap & lh) const
  {
    Ngs_Element ngel = ma->GetElement(ElementId(VOL,elnr));
    ELEMENT_TYPE eltype = ngel.GetType();
    
    if (!ho_div_free) throw Exception ("You don't have hodivfree active. You are not allow to call GetHODivFE");
    
    /*
    if (ma->GetElType(elnr) == ET_TRIG && order <= 6)
      {
	HDivHighOrderFiniteElementFO<2> * hofe2d = 0; 
	switch (order)
	  {
	  case 1: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,1> (); break;
	  case 2: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,2> (); break;
	  case 3: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,3> (); break;
	  case 4: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,4> (); break;
	  case 5: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,5> (); break;
	  case 6: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,6> (); break;
	  }
	
	Ngs_Element ngel = ma->GetElement<2> (elnr);
	for (int j = 0; j < 3; j++)
	  hofe2d->SetVertexNumber (j, ngel.vertices[j]);

        hofe2d -> SetOnlyHODiv (true);
        hofe2d -> ComputeNDof();
	
	return *hofe2d;
      }  
    */

    switch (eltype)
      {
        // case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, true, lh);
        
      case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, true, lh);
      case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, true, lh);
        
      case ET_TET:     return T_GetFE<ET_TET> (elnr, true, lh);
      case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, true, lh);
        // case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, false, lh);
      case ET_HEX:     return T_GetFE<ET_HEX> (elnr, true, lh);
        
      default:
        throw Exception ("illegal element in HDivHOFeSpace::GetDivFE");
      }
  }



//   const FiniteElement & HDivHighOrderFESpace :: GetSFE (int selnr, LocalHeap & lh) const
//   {
//     FiniteElement * fe = 0;

//     int porder; 
//     if (discont) porder = -1; 
//     else porder = order; 

//     // if (highest_order_dc) porder--;

//     switch (ma->GetSElType(selnr))
//       {
//       case ET_SEGM:
//         fe = new (lh) HDivHighOrderNormalSegm<TrigExtensionMonomial> (porder); 
//         break;
//       case ET_TRIG: 
//         fe = new (lh) HDivHighOrderNormalTrig<TrigExtensionMonomial> (porder); 
//         break; 
//       case ET_QUAD: 
//         fe = new (lh) HDivHighOrderNormalQuad<TrigExtensionMonomial> (porder); 
//         break; 
//       default:
//         throw Exception (string("HDivHighOrderFESpace::GetSFE: unsupported element ")+
//                          ElementTopology::GetElementName(ma->GetSElType(selnr)));
//       }

//     if (discont) return *fe; 

//     ArrayMem<int,4> vnums;
//     ArrayMem<int, 4> ednums, order_ed;
//     INT<3> order_fa;
//     ma->GetSElVertices(selnr, vnums);
    
//     if(ma->GetSElType(selnr) == ET_SEGM)
//       {
// 	HDivHighOrderNormalFiniteElement<1> * hofe =
// 	  dynamic_cast<HDivHighOrderNormalFiniteElement<1>*> (fe);

// 	hofe -> SetVertexNumbers (vnums);
// 	ma->GetSElEdges(selnr, ednums);
// 	// int dec = (!boundary_facet[ednums[0]] && highest_order_dc) ? 1 : 0;
// 	hofe -> SetOrderInner (order_facet[ednums[0]][0] /* -dec */);
// 	hofe -> ComputeNDof();
//       }
//     else
//       {
// 	HDivHighOrderNormalFiniteElement<2> * hofe =
// 	  dynamic_cast<HDivHighOrderNormalFiniteElement<2>*> (fe);

// 	hofe -> SetVertexNumbers (vnums);
	
// #ifdef NEW_HDIVFE
// 	INT<3> order_fa = INT<3>(order_facet[ma->GetSElFace(selnr)][0],
//                                  order_facet[ma->GetSElFace(selnr)][1],0);
//         if (highest_order_dc)
//           {
//             order_fa[0]--;
//             order_fa[1]--;
//             order_fa[2]--;
//           }
// 	hofe -> SetOrderInner (order_fa);
// #else 
// 	int order_fa = order_facet[ma->GetSElFace(selnr)][0];
//         // if (highest_order_dc) order_fa--;
// 	hofe -> SetOrderInner (order_fa);
// #endif
// 	hofe -> ComputeNDof();
//       }
    
//     return *fe;
//   }
  
  size_t HDivHighOrderFESpace :: GetNDof () const throw()
  {
    return ndof;
  }

  size_t HDivHighOrderFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }



  void HDivHighOrderFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    dnums.SetSize0();

    if (order_policy == NODE_TYPE_ORDER)
      {
        auto et = ma->GetElType(ei);
        cout << "new order policy, et = " << et
             << ", ol = " << et_order_left[et]
             << ", or = " << et_order_right[et]
             << endl;
      }

    
    if(discont) 
      {
	// lowest_order included in inner
	if(ei.VB()==VOL)
	  dnums += GetElementDofs (ei.Nr());
	return;
      } 
    if(ei.VB()==VOL)
      {
	// ArrayMem<int,6> fanums;
	// ma->GetElFacets (ei, fanums);
        auto fanums = ma->GetElFacets(ei);
	if (highest_order_dc)
	  {
	    if (ma->GetDimension() == 2)
	      {
		IntRange eldofs = GetElementDofs (ei.Nr());
		
		dnums += fanums;
		
		int first_el_dof = eldofs.First();
                for (auto f : fanums)
		  {
		    dnums += GetFacetDofs (f);
		    if (!boundary_facet[f])
		      dnums += first_el_dof++;
		  }
		dnums += IntRange (first_el_dof, eldofs.Next());
	      }
	    else // dim not 2
	      {
		IntRange eldofs = GetElementDofs (ei.Nr());
		
		// for (int i = 0; i < fanums.Size(); i++)
		// dnums.Append (fanums[i]);
		dnums += fanums;
		
		int firstel = eldofs.First();
		
		for(int i = 0; i < fanums.Size(); i++)
		  {
		    int firstfa = GetFacetDofs(fanums[i]).First();
		    
		    if (!boundary_facet[fanums[i]])
		      {
			for (int i = 0; i <= order-1; i++)
			  {
			    for(int j = 0; j < order-i-1; j++)
			      dnums.Append(firstfa++);
			    dnums.Append(firstel++);
			  }
			for (int i = 0; i < order-1; i++)
			  dnums.Append(firstfa++);
			dnums.Append(firstel++);
		      }
		    else
		      dnums += GetFacetDofs (fanums[i]);
		  }
		dnums += IntRange (firstel, eldofs.Next());
	      }
	  }
	else // not highest order dc
	  {
	    //Raviart-Thomas
            dnums += fanums;
	    // facets
	    for (auto f : fanums)
	      dnums += GetFacetDofs (f);
	    
	    //inner
	    dnums += GetElementDofs (ei.Nr());
	  }
	
	if (!DefinedOn (ei))
          dnums.SetSize0();
      }
    else if (ei.VB()==BND)
      {
        size_t fanum = ma->GetElFacets(ei)[0];
	// lowest-order
        dnums += fanum;
	
	// high order
        dnums += GetFacetDofs (fanum);
        
	if (!DefinedOn (ei))
          dnums.SetSize0();
        // dnums = -1;
      }
  }



  // ****************************
  // 
  //    smoothing blocks
  //
  //  0) Jacobi
  //  1) 2d-Vertex / 3d-Edge blocks + F + I  --- default
  //  2) 2d: edge by edge,  3d: face by face

  shared_ptr<Table<int>> HDivHighOrderFESpace :: CreateSmoothingBlocks (const Flags & precflags) const
  {
    int first;
    int ncnt = 0;
    // int ni = ma->GetNE(); //nel;
  
    int dim = ma->GetDimension();
  
  
    int SmoothingType = int(precflags.GetNumFlag("blocktype",1)); 
  
  
    Array<int> vnums,elnums; 
    Array<int> orient; 
    Array<int> ednums, fanums;
    Array<int> edges(3);
  
    int ned = ma->GetNEdges();
    int nfa = ma->GetNFaces();
    int nnv = ma->GetNV();
    int nel = ma->GetNE();
  
    cout << "SmoothingType " << SmoothingType << endl; 
    switch(SmoothingType) 
      {
      case 0: 
	cout << "Local Preconditioner" << endl;
	ncnt = ndof; //_used;
	break;  
      case 1:
	if( dim == 2 )
	  {
	    cout << "Vertex blocks + E + I" << endl;
	    ncnt = nnv + ned + nel ;
	  }
	else
	  {
	    cout << "Edge blocks + F + I" << endl;
	    ncnt = ned + nfa + nel;
	  }
	break;

      case 2:
	if( dim == 2 )
	  {
	    cout << "Edge blocks" << endl;
	    ncnt = ned ;
	  }
	else
	  {
	    cout << "Face blocks" << endl;
	    ncnt = nfa;
	  }
	break;

      }


    Array<int> cnt(ncnt); 
    cnt = 0;
    // ii=0; 

    int offset = 0;




    switch (SmoothingType)
      {
	//      0 ..... jacobi
      case 0:   // diagonal
	for(int i=0; i<ndof; i++ )
          cnt[i] = 1;
	break;
	
        //      1 ..... vertex/edge blocks -- E -- I
      case 1:   
	if( dim == 2 )
	  {
	    // vertex blocks
	    for(int i=0; i<ned; i++)
	      if(fine_facet[i])
                {
                  /*
                  int pn1, pn2;
                  ma->GetEdgePNums ( i, pn1, pn2);
                  cnt[offset + pn1] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
                  cnt[offset + pn2] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
                  */
                  auto pn = ma->GetEdgePNums(i);
                  cnt[offset+pn[0]] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
                  cnt[offset+pn[1]] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
                }

	    offset += nnv;
	    // edges
	    for (auto i : Range(ned))
	      if (fine_facet[i])
                cnt[offset+i] += first_facet_dof[i+1] - first_facet_dof[i];;
	    offset += ned;

	    // cells
	    for (auto i : Range(nel))
              cnt[offset+i] += first_inner_dof[i+1] - first_inner_dof[i];;
	  }

	else
	  {
	    // vertex blocks
	    for(auto i : Range(nfa))
	      if(fine_facet[i])
		{
		  // Array<int> edges;
		  // ma->GetFaceEdges ( i, edges);
                  /*
		  for ( int j = 0; j < edges.Size(); j++ )
		    cnt[offset + edges[j]] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
                  */
                  for (auto e : ma->GetFaceEdges(i))
                    cnt[offset+e] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
		}
	    
	    offset += ned;
	    // edges
	    for (auto i : Range(nfa))
	      if (fine_facet[i])
                cnt[offset+i] += first_facet_dof[i+1] - first_facet_dof[i];;
	    offset += nfa;

	    // cells
	    for (auto i : Range(nel))
              cnt[offset+i] += first_inner_dof[i+1] - first_inner_dof[i];;
	  }
	break;
      case 2:
	if( dim == 2 )
          cerr << "not implemented" << endl;
	else
	  {
	    for (auto i : Range(nfa))
	      if (fine_facet[i])
		cnt[i] += first_facet_dof[i+1] - first_facet_dof[i];
	  }
      }


    Table<int> table(cnt); 
  
    // ii = 0; 
    cnt = 0;
  
    offset =0;
  
    switch(SmoothingType) 
      {
	//      0 ..... jacobi
      case 0:
	for(int i=0; i<ndof; i++)
	  table[i][cnt[i]] = i;
	break;
	
	//      1 ..... vertex/edge blocks -- E -- I
      case 1:

	if ( dim == 2 )
	  {
	    for (int i = 0; i < ned; i++)
	      {
		auto pts = ma->GetEdgePNums (i);
		int pn1 = pts[0], pn2 = pts[1];
		if( fine_facet[i] )
		  {
		    table[offset + pn1][cnt[offset + pn1]++] = i;
		    table[offset + pn2][cnt[offset + pn2]++] = i;
		    
		    first = first_facet_dof[i];
		    int last = first_facet_dof[i+1];
		    for( int l=first; l<last; l++)
		      {
			table[offset + pn1][cnt[offset + pn1]++] = l;
			table[offset + pn2][cnt[offset + pn2]++] = l;
		      }
		  }
	      }

	    offset += nnv;

	    for (int i = 0; i < ned; i++ )
	      {
		first = first_facet_dof[i];
		int last = first_facet_dof[i+1];
		for ( int l = first; l < last; l++ )
		  table[offset + i ] [cnt[offset+i]++] = l;
	      }

	    for (int i = 0; i < nel; i++ )
	      {
		first = first_inner_dof[i];
		int last = first_inner_dof[i+1];
		for ( int l = first; l < last; l++ )
		  table[offset + i ] [cnt[offset+i]++] = l;
	      }
	    //(*testout) << table << endl;
	  }

	else // 3d

	  {
	    for (int i = 0; i < nfa; i++)
	      {
		if ( ! fine_facet[i] ) continue;
		Array<int> faces;
		ma->GetFaceEdges (i,faces);	      
		
	
		for ( int j = 0; j < faces.Size(); j++ )
		  {
		    table[offset + faces[j]][cnt[offset + faces[j]]++] = i;
		    
		    first = first_facet_dof[i];
		    int last = first_facet_dof[i+1];
		    for( int l=first; l<last; l++)
		      {
			table[offset + faces[j]][cnt[offset + faces[j]]++] = l;
		      }
		  }
	      }

	    offset += ned;

	    for (int i = 0; i < nfa; i++ )
	      {
		first = first_facet_dof[i];
		int last = first_facet_dof[i+1];
		for ( int l = first; l < last; l++ )
		  table[offset + i ] [cnt[offset+i]++] = l;
	      }

	    for (int i = 0; i < nel; i++ )
	      {
		first = first_inner_dof[i];
		int last = first_inner_dof[i+1];
		for ( int l = first; l < last; l++ )
		  table[offset + i ] [cnt[offset+i]++] = l;
	      }
	    //(*testout) << table << endl;
	    
	  }
	break;

      case 2:
	
	if ( dim == 2 )
	  {
	    cout << "not implemented" << endl;
	  }
	
	else // 3d

	  {
	    for (int i = 0; i < nfa; i++ )
	      {
		first = first_facet_dof[i];
		int last = first_facet_dof[i+1];
		for ( int l = first; l < last; l++ )
		  table[i] [cnt[i]++] = l;
	      }
	  }
	break;
      }
    


    //  *testout << table << endl;
    // cout << "success " << endl;
    return make_shared<Table<int>> (table);

  }







  // ******************** //
  // Direct Solver Clusters
  // 0) none
  // 1) low order dofs  --  default


  shared_ptr<Array<int>> HDivHighOrderFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    auto spclusters = make_shared<Array<int>> (ndof);
    Array<int> & clusters = *spclusters;

    int clustertype = int(precflags.GetNumFlag("ds_cluster",1)); 
    cout << " DirectSolverCluster Clustertype " << clustertype << endl; 
  
    Array<int> vnums,elnums; 
    Array<int> orient; 
  
    Array<int> edges(3);
        
    int nfa = ma->GetNFaces();

    Array<int> ednums, fnums, pnums;
  
    switch (clustertype)
      {
        // 0) none
      case 0:
        clusters = 0;

        // 1) low-order dofs
      case 1: 
        clusters = 0;

        for(int i=0; i<nfa; i++ )
          if( fine_facet[i] )
            clusters[i] = 1;
        break;

      }
    return spclusters;

  }


/** calculates [du1/dx1 du2/dx1 (du3/dx1) du1/dx2 du2/dx2 (du3/dx2) (du1/dx3 du2/dx3 du3/dx3)] */
    template<int D, int BMATDIM>
    void CalcDShapeOfHDivFE(const HDivFiniteElement<D>& fel_u, const MappedIntegrationPoint<D,D>& sip, SliceMatrix<> bmatu, LocalHeap& lh){
      HeapReset hr(lh);
      // bmatu = 0;
      // evaluate dshape by numerical diff
      //fel_u, eltrans, sip, returnval, lh
      int nd_u = fel_u.GetNDof();
      const IntegrationPoint& ip = sip.IP();//volume_ir[i];
      const ElementTransformation & eltrans = sip.GetTransformation();
      FlatMatrixFixWidth<D> shape_ul(nd_u, lh);
      FlatMatrixFixWidth<D> shape_ur(nd_u, lh);
      FlatMatrixFixWidth<D> shape_ull(nd_u, lh);
      FlatMatrixFixWidth<D> shape_urr(nd_u, lh);
      FlatMatrixFixWidth<D> dshape_u_ref(nd_u, lh);//(shape_ur); ///saves "reserved lh-memory"
      FlatMatrixFixWidth<D> dshape_u(nd_u, lh);//(shape_ul);///saves "reserved lh-memory"

      double eps = 1e-4;
      for (int j = 0; j < D; j++)   // d / dxj
      {
        IntegrationPoint ipl(ip);
        ipl(j) -= eps;
        MappedIntegrationPoint<D,D> sipl(ipl, eltrans);

        IntegrationPoint ipr(ip);
        ipr(j) += eps;
        MappedIntegrationPoint<D,D> sipr(ipr, eltrans);

        IntegrationPoint ipll(ip);
        ipll(j) -= 2*eps;
        MappedIntegrationPoint<D,D> sipll(ipll, eltrans);

        IntegrationPoint iprr(ip);
        iprr(j) += 2*eps;
        MappedIntegrationPoint<D,D> siprr(iprr, eltrans);

        fel_u.CalcMappedShape (sipl, shape_ul);
        fel_u.CalcMappedShape (sipr, shape_ur);
        fel_u.CalcMappedShape (sipll, shape_ull);
        fel_u.CalcMappedShape (siprr, shape_urr);

        dshape_u_ref = (1.0/(12.0*eps)) * (8.0*shape_ur-8.0*shape_ul-shape_urr+shape_ull);

        // dshape_u_ref = (1.0/(2*eps)) * (shape_ur-shape_ul);
        // dshape_u_ref = (1.0/(4*eps)) * (shape_urr-shape_ull);

        /*
	  for (int k = 0; k < nd_u; k++)
          for (int l = 0; l < D; l++)
          bmatu(k, j*D+l) = dshape_u_ref(k,l);
        */
        for (int l = 0; l < D; l++)
          bmatu.Col(j*D+l) = dshape_u_ref.Col(l);
      }
      
      for (int j = 0; j < D; j++)
	{
	  for (int k = 0; k < nd_u; k++)
	    for (int l = 0; l < D; l++)
	      dshape_u_ref(k,l) = bmatu(k, l*D+j);
	  
	  dshape_u = dshape_u_ref * sip.GetJacobianInverse();

	  for (int k = 0; k < nd_u; k++)
	    for (int l = 0; l < D; l++)
	      bmatu(k, l*D+j) = dshape_u(k,l);
	}
    }
  

  /// Gradient operator for HDiv
  template <int D, typename FEL = HDivFiniteElement<D> >
  class DiffOpGradientHdiv : public DiffOp<DiffOpGradientHdiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 1 };
    static Array<int> GetDimensions() { return Array<int> ( { D, D } ); };
    
    static constexpr double eps() { return 1e-4; } 
    ///
    template <typename AFEL, typename SIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
      static void GenerateMatrix (const AFEL & fel, const SIP & sip,
                                  MAT & mat, LocalHeap & lh)
    {
      cout << "nicht gut" << endl;
      cout << "type(fel) = " << typeid(fel).name() << ", sip = " << typeid(sip).name()
           << ", mat = " << typeid(mat).name() << endl;
    }
    
    // template <typename AFEL, typename SIP>
    // static void GenerateMatrix (const AFEL & fel, const SIP & sip,
    // SliceMatrix<double,ColMajor> mat, LocalHeap & lh)
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
                                                  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                                                              MAT mat, LocalHeap & lh)
    {
      CalcDShapeOfHDivFE<D,D*D>(static_cast<const FEL&>(fel), mip, Trans(mat), lh);
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> mat)
    {
      // throw ExceptionNOSIMD("generategrad not simded");
      auto & fel = static_cast<const FEL&>(bfel);
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D,D>&> (bmir);
      
      size_t nd_u = fel.GetNDof();

      STACK_ARRAY(SIMD<double>, mem1, 6*D*nd_u);
      FlatMatrix<SIMD<double>> shape_ul(nd_u*D, 1, &mem1[0]);
      FlatMatrix<SIMD<double>> shape_ur(nd_u*D, 1, &mem1[D*nd_u]);
      FlatMatrix<SIMD<double>> shape_ull(nd_u*D, 1, &mem1[2*D*nd_u]);
      FlatMatrix<SIMD<double>> shape_urr(nd_u*D, 1, &mem1[3*D*nd_u]);

      FlatMatrix<SIMD<double>> dshape_u_ref(nd_u*D, 1, &mem1[4*D*nd_u]);
      FlatMatrix<SIMD<double>> dshape_u(nd_u*D, 1, &mem1[5*D*nd_u]);

      LocalHeapMem<10000> lh("diffopgrad-lh");

      auto & ir = mir.IR();
      for (size_t i = 0; i < mir.Size(); i++)
        {
          const SIMD<IntegrationPoint> & ip = ir[i];
          const ElementTransformation & eltrans = mir[i].GetTransformation();

          // double eps = 1e-4;
          for (int j = 0; j < D; j++)   // d / dxj
            {
              HeapReset hr(lh);
              SIMD<IntegrationPoint> ipts[4];
              ipts[0] = ip;
              ipts[0](j) -= eps();
              ipts[1] = ip;
              ipts[1](j) += eps();              
              ipts[2] = ip;
              ipts[2](j) -= 2*eps();
              ipts[3] = ip;
              ipts[3](j) += 2*eps();

              SIMD_IntegrationRule ir(4, ipts);
              SIMD_MappedIntegrationRule<D,D> mirl(ir, eltrans, lh);

              fel.CalcMappedShape (mirl[0], shape_ul);
              fel.CalcMappedShape (mirl[1], shape_ur);
              fel.CalcMappedShape (mirl[2], shape_ull);
              fel.CalcMappedShape (mirl[3], shape_urr);

              dshape_u_ref = (1.0/(12.0*eps())) * (8.0*shape_ur-8.0*shape_ul-shape_urr+shape_ull);
              for (size_t l = 0; l < D; l++)
                for (size_t k = 0; k < nd_u; k++)
                  mat(k*D*D+j*D+l, i) = dshape_u_ref(k*D+l, 0);
            }
          
          for (size_t j = 0; j < D; j++)
            for (size_t k = 0; k < nd_u; k++)
              {
                Vec<D,SIMD<double>> dshape_u_ref, dshape_u;
                for (size_t l = 0; l < D; l++)
                  dshape_u_ref(l) = mat(k*D*D+l*D+j, i);
                
                dshape_u = Trans(mir[i].GetJacobianInverse()) * dshape_u_ref;
                
                for (size_t l = 0; l < D; l++)
                  mat(k*D*D+l*D+j, i) = dshape_u(l);
              }
        }
    }
    
    /*
    template <typename AFEL>
    static void GenerateMatrix (const AFEL & fel, 
                                const MappedIntegrationPoint<D,D> & sip,
                                FlatMatrixFixHeight<D> & mat, LocalHeap & lh)
    {
      FlatMatrixFixWidth<D*D> hm(fel.GetNDof(), &mat(0,0));
      CalcDShapeOfHDivFE<D,D*D>(static_cast<const FEL&>(fel), sip, hm, lh);
    }
    */
    ///
    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL & fel, const MIP & mip,
                       const TVX & x, TVY & y,
                       LocalHeap & lh) 
    {
      // typedef typename TVX::TSCAL TSCAL;
      HeapReset hr(lh);
      FlatMatrixFixWidth<D*D> hm(fel.GetNDof(),lh);
      CalcDShapeOfHDivFE<D,D*D>(static_cast<const FEL&>(fel), mip, hm, lh);
      y = Trans(hm)*x;
    }


    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void ApplyTrans (const AFEL & fel, const MIP & mip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCALX;      
      // typedef typename TVY::TSCAL TSCAL;
      // typedef typename MIP::TSCAL TSCAL;

      HeapReset hr(lh);
      FlatMatrixFixWidth<D*D> bmatu(fel.GetNDof(),lh);
      
      auto & fel_u = static_cast<const FEL&>(fel);
      int nd_u = fel.GetNDof();
      const IntegrationPoint& ip = mip.IP();
      const ElementTransformation & eltrans = mip.GetTransformation();
      FlatMatrixFixWidth<D> shape_ul(nd_u, lh);
      FlatMatrixFixWidth<D> shape_ur(nd_u, lh);
      FlatMatrixFixWidth<D> dshape_u_ref(nd_u, lh);
      FlatMatrixFixWidth<D> dshape_u(nd_u, lh);
      
      FlatMatrix<TSCALX> hx(D,D,&x(0));
      Mat<D,D,TSCALX> tx = mip.GetJacobianInverse() * hx;

      y = 0;
      for (int j = 0; j < D; j++)   // d / dxj
	{
	  IntegrationPoint ipl(ip);
	  ipl(j) -= eps();
	  MappedIntegrationPoint<D,D> sipl(ipl, eltrans);

	  IntegrationPoint ipr(ip);
	  ipr(j) += eps();
	  MappedIntegrationPoint<D,D> sipr(ipr, eltrans);

	  fel_u.CalcMappedShape (sipl, shape_ul);
	  fel_u.CalcMappedShape (sipr, shape_ur);
	  dshape_u_ref = (1.0/(2*eps())) * (shape_ur-shape_ul);
          y += dshape_u_ref * tx.Row(j);
	}
    }

    using DiffOp<DiffOpGradientHdiv<D>>::ApplySIMDIR;
    /*
    template <typename AFEL, class MIR, class TVX, class TVY>
    static void ApplySIMDIR (const AFEL & fel, const MIR & bmir,
                             const TVX & x, TVY & y)
    */
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      int size = (bmir.Size()+1)*2000;
      STACK_ARRAY(char, data, size);
      LocalHeap lh(data, size);

      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D,D>&> (bmir);
      auto & ir = mir.IR();
      const ElementTransformation & trafo = mir.GetTransformation();
      auto & fel_u = static_cast<const FEL&>(fel);
      FlatMatrix<SIMD<double>> hxl(D, mir.IR().Size(), lh);
      FlatMatrix<SIMD<double>> hxr(D, mir.IR().Size(), lh);
      FlatMatrix<SIMD<double>> hxll(D, mir.IR().Size(), lh);
      FlatMatrix<SIMD<double>> hxrr(D, mir.IR().Size(), lh);
      FlatMatrix<SIMD<double>> hx(D, mir.IR().Size(), lh);

      for (int k = 0; k < mir.Size(); k++)
        for (int m = 0; m < D*D; m++)
          y(m, k) = SIMD<double> (0.0);
      
      for (int j = 0; j < D; j++)
        {
          // hx = (F^-1 * x).Row(j)
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irl(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irl.Size(); k++)
              {
                irl[k] = ir[k];
                irl[k](j) -= eps();
              }
            SIMD_MappedIntegrationRule<D,D> mirl(irl, trafo, lh);
            fel_u.Evaluate (mirl, x, hxl);
          }
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irr(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irr.Size(); k++)
              {
                irr[k] = ir[k];              
                irr[k](j) += eps();
              }
            SIMD_MappedIntegrationRule<D,D> mirr(irr, trafo, lh);
            fel_u.Evaluate (mirr, x, hxr);
          }
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irll(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irll.Size(); k++)
              {
                irll[k] = ir[k];
                irll[k](j) -= 2*eps();
              }
            SIMD_MappedIntegrationRule<D,D> mirll(irll, trafo, lh);
            fel_u.Evaluate (mirll, x, hxll);
          }
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irrr(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irrr.Size(); k++)
              {
                irrr[k] = ir[k];              
                irrr[k](j) += 2*eps();
              }
            SIMD_MappedIntegrationRule<D,D> mirrr(irrr, trafo, lh);
            fel_u.Evaluate (mirrr, x, hxrr);
          }
          // hx = 1.0/(2*eps()) * (hxr-hxl);
          // dshape_u_ref = (1.0/(12.0*eps)) * (8.0*shape_ur-8.0*shape_ul-shape_urr+shape_ull);
          hx = 1.0/(12*eps()) * (8*hxr-8*hxl-hxrr+hxll);
          for (int k = 0; k < mir.Size(); k++)
            {
              auto jacinv = mir[k].GetJacobianInverse();
              for (int l = 0; l < D; l++)
                {
                  for (int m = 0; m < D; m++)
                    y(m*D+l, k) += jacinv(j,m) * hx(l, k);
                }
            }
        }
    }

    using DiffOp<DiffOpGradientHdiv<D>>::AddTransSIMDIR;    
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
    {
      size_t size = (bmir.Size()+1)*2000;
      STACK_ARRAY(char, data, size);
      LocalHeap lh(data, size);

      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D,D>&> (bmir);
      auto & ir = mir.IR();
      const ElementTransformation & trafo = mir.GetTransformation();
      auto & fel_u = static_cast<const FEL&>(fel);
      // AFlatMatrix<double> hx(D, mir.IR().GetNIP(), lh);
      FlatMatrix<SIMD<double>> hx1(D, mir.Size(), lh);
      FlatMatrix<SIMD<double>> hx2(D, mir.Size(), lh);

      for (size_t j = 0; j < D; j++)
        {
          // hx = (F^-1 * x).Row(j)
          for (size_t k = 0; k < mir.Size(); k++)
            {
              auto jacinv = mir[k].GetJacobianInverse();
              for (int l = 0; l < D; l++)
                {
                  SIMD<double> sum = 0;
                  for (int m = 0; m < D; m++)
                    sum += jacinv(j,m) * x(m*D+l, k);
                  // hx.Get(l,k) = (-(0.5/eps()) * sum).Data();
                  hx1(l,k) = (-(8/(12*eps())) * sum).Data();
                  hx2(l,k) = ( (1/(12*eps())) * sum).Data();
                }
            }
          
          /*
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irl(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irl.Size(); k++)
              {
                irl[k] = ir[k];
                irl[k](j) -= eps();
              }
            SIMD_MappedIntegrationRule<D,D> mirl(irl, trafo, lh);
            fel_u.AddTrans (mirl, hx, y);
          }
          {
            HeapReset hr(lh);
            hx *= -1;
            SIMD_IntegrationRule irr(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irr.Size(); k++)
              {
                irr[k] = ir[k];              
                irr[k](j) += eps();
              }
            SIMD_MappedIntegrationRule<D,D> mirr(irr, trafo, lh);
            fel_u.AddTrans (mirr, hx, y);
          }
          */
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irl(mir.IR().GetNIP(), lh);
            for (size_t k = 0; k < irl.Size(); k++)
              {
                irl[k] = ir[k];
                irl[k](j) -= eps();
              }
            SIMD_MappedIntegrationRule<D,D> mirl(irl, trafo, lh);
            fel_u.AddTrans (mirl, hx1, y);
            irl.NothingToDelete();
          }
          {
            HeapReset hr(lh);
            hx1 *= -1;
            SIMD_IntegrationRule irr(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irr.Size(); k++)
              {
                irr[k] = ir[k];              
                irr[k](j) += eps();
              }
            SIMD_MappedIntegrationRule<D,D> mirr(irr, trafo, lh);
            fel_u.AddTrans (mirr, hx1, y);
          }
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irl(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irl.Size(); k++)
              {
                irl[k] = ir[k];
                irl[k](j) -= 2*eps();
              }
            SIMD_MappedIntegrationRule<D,D> mirl(irl, trafo, lh);
            fel_u.AddTrans (mirl, hx2, y);
          }
          {
            HeapReset hr(lh);
            hx2 *= -1;
            SIMD_IntegrationRule irr(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irr.Size(); k++)
              {
                irr[k] = ir[k];              
                irr[k](j) += 2*eps();
              }
            SIMD_MappedIntegrationRule<D,D> mirr(irr, trafo, lh);
            fel_u.AddTrans (mirr, hx2, y);
          }
        }
    }
    
  };


  SymbolTable<shared_ptr<DifferentialOperator>>
  HDivHighOrderFESpace :: GetAdditionalEvaluators () const
  {
    SymbolTable<shared_ptr<DifferentialOperator>> additional;
    switch (ma->GetDimension())
      {
      case 1:
        additional.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHdiv<1>>> ()); break;
      case 2:
        additional.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHdiv<2>>> ()); break;
      case 3:
        additional.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHdiv<3>>> ()); break;
      default:
        ;
      }
    return additional;
  }
  
  

  
  void HDivHighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize0(); 
  }
  
  void HDivHighOrderFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  { 
    dnums.SetSize0();
    if(ma->GetDimension() == 3 || discont) return; 

    dnums += ednr;
    dnums += GetFacetDofs (ednr);
  }

  void HDivHighOrderFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 2 || discont) return; 
   
    dnums += fanr;
    dnums += GetFacetDofs (fanr);
  }

  void HDivHighOrderFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums = GetElementDofs (elnr);
  }


  static RegisterFESpace<HDivHighOrderFESpace> init ("hdivho");
}



