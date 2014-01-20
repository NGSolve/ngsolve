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
#include <../fem/hdivhofe.hpp>  
#include <../fem/hdivhofefo.hpp>  


namespace ngcomp
{

  HDivHighOrderFESpace ::  
  HDivHighOrderFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
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

    uniform_order_inner = int (flags.GetNumFlag ("orderinner", -1));
    
    if(flags.NumFlagDefined("orderedge") || flags.NumFlagDefined ("orderface")) 
      throw Exception("Flag 'orderface' and 'orderedge' for hdivho are obsolete. Use flag 'orderfacet' instead!"); 
	
    uniform_order_facet = int (flags.GetNumFlag ("orderfacet", -1));
    
    // Evaluator for shape tester 
    if (ma.GetDimension() == 2)
      {
        Array<CoefficientFunction*> coeffs(1);
        coeffs[0] = new ConstantCoefficientFunction(1);
        integrator = GetIntegrators().CreateBFI("masshdiv", 2, coeffs);
        boundary_integrator = GetIntegrators().CreateBFI("robinhdiv", 2, coeffs);
      }
    else
      {
        Array<CoefficientFunction*> coeffs(1);
        coeffs[0] = new ConstantCoefficientFunction(1);
        integrator = GetIntegrators().CreateBFI("masshdiv", 3, coeffs);
        boundary_integrator = GetIntegrators().CreateBFI("robinhdiv", 3, coeffs);
      }
    if (dimension > 1)
      {
        integrator = new BlockBilinearFormIntegrator (*integrator, dimension);
        boundary_integrator =
          new BlockBilinearFormIntegrator (*boundary_integrator, dimension);
      }

    vefc_dofblocks = Vec<4,int> (0,0,2,1);
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
    
    // int nv = ma.GetNV();
    int nel = ma.GetNE();
    int nfa = ma.GetNFacets();
       
    order_facet.SetSize(nfa);
    order_inner.SetSize(nel); 
    order_inner_curl.SetSize(nel); 
    fine_facet.SetSize(nfa); 
   
     
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

    int dim = ma.GetDimension();
    
    for(int i = 0; i < nel; i++)
      {
	ELEMENT_TYPE eltype = ma.GetElType(i); 
	const POINT3D * points = ElementTopology :: GetVertices (eltype);
	
	Array<int> elfacets; 
	ma.GetElFacets (i,elfacets); 
	
	for (int j=0;j<elfacets.Size();j++) fine_facet[elfacets[j]] = 1; 
	
	if(!var_order) continue; 
	
	INT<3> el_orders = ma.GetElOrders(i);  
	
	for(int k=0;k<dim;k++)
	  {
	    order_inner_curl[i][k]= max(el_orders[k] + rel_curl_order,0);
	    order_inner[i][k] = max(el_orders[k]+rel_order,0);
	  }
	  
	if(dim==2)
	  {
	    const EDGE * edges = ElementTopology::GetEdges (eltype);
	    for(int j=0; j<elfacets.Size(); j++)
	      for(int k=0;k<2;k++)
		if(points[edges[j][0]][k] != points[edges[j][1]][k])
		  { 
		    order_facet[elfacets[j]][0] = max(el_orders[k]+rel_curl_order, order_facet[elfacets[j]][0]);
		    break; 
		  }
	  }
	else
	  {
	    Array<int> vnums; 
	    ma.GetElVertices (i, vnums);
	    const FACE * faces = ElementTopology::GetFaces (eltype);
	    for(int j=0;j<elfacets.Size();j++)
	      {
		if(faces[j][3]==-1) // trig  
		  {
		    order_facet[elfacets[j]][0] = max(order_facet[elfacets[j]][0],el_orders[0]+rel_curl_order);
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
                            order_facet[elfacets[j]][l] = max(order_facet[elfacets[j]][l], rel_curl_order + 
                                                              el_orders[k]);
                            break; 
                          } 
		    
		  }
	      }
	  }
      }

    ma.AllReduceNodalData ((ma.GetDimension()==2) ? NT_EDGE : NT_FACE,
			   fine_facet, MPI_LOR);

    if(uniform_order_inner > -1) order_inner = uniform_order_inner;

    if(uniform_order_facet > -1) order_facet = uniform_order_facet;

    for(int i=0;i<nfa;i++) if(!fine_facet[i]) order_facet[i] = INT<2> (0,0); 

    // by SZ ... since order_inner_curl is not working yet for hdivhofe 	 
    for(int i=0; i<order_inner_curl.Size(); i++) 
      order_inner_curl[i] = order_inner[i]; 

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
    int nfa = ma.GetNFacets();
    int nel = ma.GetNE();
    int dim = ma.GetDimension();
    Array<int> pnums; 
     
    first_facet_dof.SetSize(nfa+1); 
    first_inner_dof.SetSize(nel+1); 

    ndof = nfa; 
    first_facet_dof = ndof; 

    if(dim==2)
      {
        for (int i = 0; i < nfa; i++)
          {
            first_facet_dof[i] = ndof;
            if(fine_facet[i])
              ndof += order_facet[i][0];
          }

        first_facet_dof[nfa] = ndof;
      
        int inci = 0;
        for(int i=0; i< nel; i++)
          {
            INT<3> pc = order_inner_curl[i];
            INT<3> p = order_inner[i];
            switch(ma.GetElType(i))
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
            if (inci < 0) inci = 0;
            {
              first_inner_dof[i] = ndof;
              ndof+=inci;
            }
          }
        first_inner_dof[nel] = ndof;
      }
    else 
      {
        int inci = 0;
        for (int i=0; i< nfa; i++) 
          {
            inci =0; 
            if(fine_facet[i])
              {
                INT<2> p = order_facet[i]; 
                ma.GetFacePNums(i,pnums);
                switch(pnums.Size())
                  {
                  case 3: //Triangle
                    inci= (p[0]*p[0]+3*p[0])/2; 
                    break;
                  case 4: //Quad 
                    inci= p[0]*p[1] + p[0] + p[1];
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
	 
        for (int i=0; i< nel; i++)
          {
            INT<3> p = order_inner[i];
            INT<3> pc = order_inner_curl[i];
            int inci = 0;
	     
            switch(ma.GetElType(i))
              {
              case ET_TET:
                if(p[0]>1 && !ho_div_free)
                  inci = p[0]*(p[0]+1)*(p[0]-1)/6 + p[0]*(p[0]-1)/2 + p[0]-1;
                if(pc[0]>1)
                  inci += pc[0]*(pc[0]+1)*(pc[0]-1)/3 + pc[0]*(pc[0]-1)/2; ;
                break;
              case ET_PRISM:
                // inci = (p[0]+1)*(3*(p[0]-1)+(p[0]-1)*(p[0]-2))+ (p[0]-1)*(p[0]+1)*(p[0]+2)/2;
		inci = (p[0]+2)*p[0]*(p[2]+1) + (p[0]+1)*(p[0]+2)*p[2]/2;
		if (ho_div_free)
		  inci -= (p[0]+1)*(p[0]+2)*(p[2]+1)/2 - 1;
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
      }
   

    if(discont) 
      { 
        ndof = 0; 
        for(int i=0;i<nel;i++)
          {
            int incii = first_inner_dof[i+1]-first_inner_dof[i]; 
	     	     
            Array<int> elfacets; 
            if(dim==2) 
              ma.GetElEdges(i,elfacets);
            else 
              ma.GetElFaces(i,elfacets);
	     
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
     
    while (ma.GetNLevels() > ndlevel.Size())
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
    
    for (int facet = 0; facet < ma.GetNFacets(); facet++)
      {
        ctofdof[facet] = WIREBASKET_DOF;
        ctofdof.Range (first_facet_dof[facet], 
                       first_facet_dof[facet+1]) = INTERFACE_DOF;
      }
    
    for (int el = 0; el < ma.GetNE(); el++)
      ctofdof.Range (first_inner_dof[el], 
                     first_inner_dof[el+1]) = LOCAL_DOF;
  }

  const FiniteElement & HDivHighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    /*
    if (ma.GetElType(elnr) == ET_TRIG && order <= 6)
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
	
	Ngs_Element ngel = ma.GetElement<2> (elnr);
	for (int j = 0; j < 3; j++)
	  hofe2d->SetVertexNumber (j, ngel.vertices[j]);

        hofe2d -> SetHODivFree (ho_div_free);
        hofe2d -> SetOnlyHODiv (false);
        hofe2d -> ComputeNDof();
	
	return *hofe2d;
      }  
    */


    FiniteElement * fe;
    
    // typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
    // typedef TrigExtensionMonomial T_ORTHOPOL;
    

    
    switch (ma.GetElType(elnr))
      {
      case ET_TET: fe = new (lh)  HDivHighOrderFE<ET_TET> (); break;
        // case ET_PYRAMID: fe = new (lh)  HDivHighOrderPyramid<ET_PYRAMID> (order);  break;
      case ET_PRISM: fe = new (lh)  HDivHighOrderFE<ET_PRISM> (); break;
      case ET_HEX:   fe = new (lh)  HDivHighOrderFE<ET_HEX> (order); break;
      case ET_TRIG:  fe = new (lh)  HDivHighOrderFE<ET_TRIG> (order); break;
      case ET_QUAD:  fe = new (lh)  HDivHighOrderFE<ET_QUAD> (order); break;
      default:
	fe = 0; 
      }
  
    if (!fe)
      {
	stringstream str;
	str << "HDivHighOrderFESpace " << GetClassName() 
	    << ", undefined eltype " 
	    << ElementTopology::GetElementName(ma.GetElType(elnr))
	    << ", order = " << order << endl;
	throw Exception (str.str());
      }
    
    Array<int> vnums;
    ma.GetElVertices(elnr, vnums);
    if ( ma.GetDimension() == 2)
      {
	
	HDivHighOrderFiniteElement<2> * hofe =
	  dynamic_cast<HDivHighOrderFiniteElement<2>*> (fe);
	ArrayMem<int, 12> ednums, order_ed;
	
	ma.GetElEdges(elnr, ednums);
	order_ed.SetSize (ednums.Size());
	
	for (int j = 0; j < ednums.Size(); j++)
	  order_ed[j] = order_facet[ednums[j]][0];

	hofe -> SetVertexNumbers (vnums);
	hofe -> SetOrderEdge (order_ed);
	
	// #ifndef NEW_HDIVFE // old version 
	//  hofe -> SetOrderInner (order_inner[elnr][0]);
	// hofe -> SetDiscontinuous(discont); 
	// #else
        // new anisotropic FE
	hofe -> SetOrderInner (order_inner[elnr]); 
	//hofe -> SetOrderInnerCurl (order_inner_curl[elnr]);
	hofe -> SetDiscontinuous(discont); 
        hofe -> SetHODivFree (ho_div_free);
        hofe -> SetOnlyHODiv (false);
	hofe -> ComputeNDof();
      }
    else // dim=3
      {
	HDivHighOrderFiniteElement<3> * hofe =
	  dynamic_cast<HDivHighOrderFiniteElement<3>*> (fe);
	
	ArrayMem<int, 6> fanums; 
	ArrayMem<INT<2>, 6> order_fa;
	ma.GetElFaces(elnr, fanums);
	order_fa.SetSize (fanums.Size());
	
	for (int j = 0; j < fanums.Size(); j++)
	  order_fa[j] = order_facet[fanums[j]];

	hofe -> SetVertexNumbers (vnums);
	
	hofe -> SetOrderFace (order_fa);
	hofe -> SetOrderInner (order_inner[elnr]);
	//hofe -> SetOrderInnerCurl (order_inner_curl[elnr]); // under construction
	hofe -> SetDiscontinuous(discont); 
        hofe -> SetHODivFree (ho_div_free);
        hofe -> SetOnlyHODiv (false);
	hofe -> ComputeNDof();
      }
    return *fe;
  }

  const FiniteElement & HDivHighOrderFESpace :: GetHODivFE (int elnr, LocalHeap & lh) const
  {
    if (!ho_div_free) throw Exception ("You don't have hodivfree active. You are not allow to call GetHODivFE");
    if (ma.GetElType(elnr) == ET_TRIG && order <= 6 && 0)
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
	
	Ngs_Element ngel = ma.GetElement<2> (elnr);
	for (int j = 0; j < 3; j++)
	  hofe2d->SetVertexNumber (j, ngel.vertices[j]);

        hofe2d -> SetOnlyHODiv (true);
        hofe2d -> ComputeNDof();
	
	return *hofe2d;
      }  



    FiniteElement * fe;
    
    switch (ma.GetElType(elnr))
      {
      case ET_TET: fe = new (lh)  HDivHighOrderFE<ET_TET> (); break;
        // case ET_PYRAMID: fe = new (lh)  HDivHighOrderPyramid<ET_PYRAMID> (order);  break;
      case ET_PRISM: fe = new (lh)  HDivHighOrderFE<ET_PRISM> (); break;
      case ET_HEX:   fe = new (lh)  HDivHighOrderFE<ET_HEX> (order); break;
      case ET_TRIG:  fe = new (lh)  HDivHighOrderFE<ET_TRIG> (order); break;
      case ET_QUAD:  fe = new (lh)  HDivHighOrderFE<ET_QUAD> (order); break;
      default:
	fe = 0; 
      }
  
    if (!fe)
      {
	stringstream str;
	str << "HDivHighOrderFESpace " << GetClassName() 
	    << ", undefined eltype " 
	    << ElementTopology::GetElementName(ma.GetElType(elnr))
	    << ", order = " << order << endl;
	throw Exception (str.str());
      }
    
    Array<int> vnums;
    ma.GetElVertices(elnr, vnums);
    if ( ma.GetDimension() == 2)
      {
	
	HDivHighOrderFiniteElement<2> * hofe =
	  dynamic_cast<HDivHighOrderFiniteElement<2>*> (fe);
	ArrayMem<int, 12> ednums, order_ed;
	
	ma.GetElEdges(elnr, ednums);
	order_ed.SetSize (ednums.Size());
	
	for (int j = 0; j < ednums.Size(); j++)
	  order_ed[j] = order_facet[ednums[j]][0];
      
	hofe -> SetVertexNumbers (vnums);
	hofe -> SetOrderEdge (order_ed);
	
	// #ifndef NEW_HDIVFE // old version 
	//  hofe -> SetOrderInner (order_inner[elnr][0]);
	// hofe -> SetDiscontinuous(discont); 
	// #else
        // new anisotropic FE
	hofe -> SetOrderInner (order_inner[elnr]); 
	//hofe -> SetOrderInnerCurl (order_inner_curl[elnr]);
	hofe -> SetDiscontinuous(discont); 
        hofe -> SetOnlyHODiv (true);
	hofe -> ComputeNDof();
      }
    else // dim=3
      {
	HDivHighOrderFiniteElement<3> * hofe =
	  dynamic_cast<HDivHighOrderFiniteElement<3>*> (fe);
	
	ArrayMem<int, 6> fanums; 
	ArrayMem<INT<2>, 6> order_fa;
	ma.GetElFaces(elnr, fanums);
	order_fa.SetSize (fanums.Size());
	
	for (int j = 0; j < fanums.Size(); j++)
	  order_fa[j] = order_facet[fanums[j]];

	hofe -> SetVertexNumbers (vnums);
	
	hofe -> SetOrderFace (order_fa);
	hofe -> SetOrderInner (order_inner[elnr]);
	//hofe -> SetOrderInnerCurl (order_inner_curl[elnr]); // under construction
	hofe -> SetDiscontinuous(discont); 
        hofe -> SetOnlyHODiv (ho_div_free);
	hofe -> ComputeNDof();
      }
    return *fe;
  }





  const FiniteElement & HDivHighOrderFESpace :: GetSFE (int selnr, LocalHeap & lh) const
  {
    FiniteElement * fe = 0;

    int porder; 
    if(discont) porder = -1; 
    else porder = order; 

    switch (ma.GetSElType(selnr))
      {
      case ET_SEGM:
        fe = new (lh) HDivHighOrderNormalSegm<TrigExtensionMonomial> (porder); 
        break;
      case ET_TRIG: 
        fe = new (lh) HDivHighOrderNormalTrig<TrigExtensionMonomial> (porder); 
        break; 
      case ET_QUAD: 
        fe = new (lh) HDivHighOrderNormalQuad<TrigExtensionMonomial> (porder); 
        break; 
      default:
        throw Exception (string("HDivHighOrderFESpace::GetSFE: unsupported element ")+
                         ElementTopology::GetElementName(ma.GetSElType(selnr)));
      }

    if (!fe)
      {
	stringstream str;
	str << "HDivHighOrderFESpace " << GetClassName()
	    << ", undefined eltype "
	    << ElementTopology::GetElementName(ma.GetSElType(selnr))
	    << ", order = " << order << endl;
	throw Exception (str.str());
      }

    if(discont) return *fe; 

    ArrayMem<int,4> vnums;
    ArrayMem<int, 4> ednums, order_ed;
    INT<3> order_fa;
    ma.GetSElVertices(selnr, vnums);
    
    if(ma.GetSElType(selnr) == ET_SEGM)
      {
	HDivHighOrderNormalFiniteElement<1> * hofe =
	  dynamic_cast<HDivHighOrderNormalFiniteElement<1>*> (fe);

	hofe -> SetVertexNumbers (vnums);
	ma.GetSElEdges(selnr, ednums);
	
	hofe -> SetOrderInner (order_facet[ednums[0]][0]);
	hofe -> ComputeNDof();
      }
    else
      {
	HDivHighOrderNormalFiniteElement<2> * hofe =
	  dynamic_cast<HDivHighOrderNormalFiniteElement<2>*> (fe);

	hofe -> SetVertexNumbers (vnums);
	
#ifdef NEW_HDIVFE
	INT<3> order_fa = INT<3>(order_facet[ma.GetSElFace(selnr)][0],order_facet[ma.GetSElFace(selnr)][1],0);
	hofe -> SetOrderInner (order_fa);
#else 
	int order_fa = order_facet[ma.GetSElFace(selnr)][0];
	hofe -> SetOrderInner (order_fa);
#endif
	hofe -> ComputeNDof();
      }
    
    return *fe;
  }
  
  int HDivHighOrderFESpace :: GetNDof () const
  {
    return ndof;
  }

  int HDivHighOrderFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }


  void HDivHighOrderFESpace :: 
  GetDofRanges (ElementId ei, Array<IntRange> & dranges) const
  {
    dranges.SetSize(0);
    if (!DefinedOn (ei)) return;

    // Ngs_Element ngel = ma.GetElement(ei);

    if (ei.IsVolume())
      {
        if(discont) 
          {
            dranges += GetElementDofs(ei.Nr());
            return;
          } 
        
        ArrayMem<int,6> fanums;
        ma.GetElFacets (ei.Nr(), fanums);
        
        //Raviart-Thomas
        for (int i = 0; i < fanums.Size(); i++)
          dranges.Append (fanums[i]);
        // facets
        for(int i=0; i<fanums.Size(); i++)
          dranges += GetFacetDofs (fanums[i]);
        
        //inner
        dranges += GetElementDofs(ei.Nr());
      }
    else
      {
        ArrayMem<int,12> fanums; 
        ma.GetSElFacets (ei.Nr(), fanums);

        dranges.Append (fanums[0]);
        dranges += GetFacetDofs (fanums[0]);
      }
  }



  void HDivHighOrderFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    if(discont) 
      {
	// lowest_order included in inner 
        dnums += GetElementDofs (elnr);
	return;
      } 

    ArrayMem<int,6> fanums;
    ma.GetElFacets (elnr, fanums);

    //Raviart-Thomas
    for (int i = 0; i < fanums.Size(); i++)
      dnums.Append (fanums[i]);
    // facets
    for(int i=0; i<fanums.Size(); i++)
      dnums += GetFacetDofs (fanums[i]);

    //inner
    dnums += GetElementDofs (elnr);
    
    if (!DefinedOn (ma.GetElIndex (elnr)))
      dnums = -1;
  }


  void HDivHighOrderFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    if (discont) return; 

    Array<int> fanums; 
       
    ma.GetSElFacets (selnr, fanums);

    // lowest-order
    for(int i=0;i<fanums.Size();i++) 
      dnums.Append (fanums[i]);

    // high order
    for (int i = 0; i < fanums.Size(); i++)
      dnums += IntRange (first_facet_dof[fanums[i]],
                         first_facet_dof[fanums[i]+1]);
  }


  // ****************************
  // 
  //    smoothing blocks
  //
  //  0) Jacobi
  //  1) 2d-Vertex / 3d-Edge blocks + F + I  --- default
  //  2) 2d: edge by edge,  3d: face by face

  Table<int> * HDivHighOrderFESpace :: CreateSmoothingBlocks (const Flags & precflags) const
  {
    int first;
    int ncnt = 0;
    // int ni = ma.GetNE(); //nel;
  
    int dim = ma.GetDimension();
  
  
    int SmoothingType = int(precflags.GetNumFlag("blocktype",1)); 
  
  
    Array<int> vnums,elnums; 
    Array<int> orient; 
    Array<int> ednums, fanums;
    Array<int> edges(3);
  
    int ned = ma.GetNEdges();
    int nfa = ma.GetNFaces();
    int nnv = ma.GetNV();
    int nel = ma.GetNE();
  
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
                  int pn1, pn2;
                  ma.GetEdgePNums ( i, pn1, pn2);
                  cnt[offset + pn1] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
                  cnt[offset + pn2] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
		
                }

	    offset += nnv;
	    // edges
	    for(int i=0; i<ned; i++)
	      if( fine_facet[i] )
                {
                  cnt[offset + i] += first_facet_dof[i+1] - first_facet_dof[i];;
                }
	    offset += ned;

	    // cells
	    for(int i=0; i<nel; i++)
	      {
		cnt[offset + i] += first_inner_dof[i+1] - first_inner_dof[i];;
	      }

	  }

	else
	  {
	    // vertex blocks
	    for(int i=0; i<nfa; i++)
	      if(fine_facet[i])
		{
		  Array<int> edges;
		  ma.GetFaceEdges ( i, edges);
		  for ( int j = 0; j < edges.Size(); j++ )
		    cnt[offset + edges[j]] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
		  
		}
	    
	    offset += ned;
	    // edges
	    for(int i=0; i<nfa; i++)
	      if( fine_facet[i] )
		{
		  cnt[offset + i] += first_facet_dof[i+1] - first_facet_dof[i];;
		}
	    offset += nfa;

	    // cells
	    for(int i=0; i<nel; i++)
	      {
		cnt[offset + i] += first_inner_dof[i+1] - first_inner_dof[i];;
	      }
	    
	  }
	break;
      case 2:
	if( dim == 2 )
	  {
	    cerr << "not implemented" << endl;
	  }

	else
	  {
	    for(int i=0; i<nfa; i++)
	      if( fine_facet[i] )
		cnt[i] += first_facet_dof[i+1] - first_facet_dof[i];
	  }
      }


    Table<int> & table = *new Table<int> (cnt); 
  
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
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);	      
		
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
		ma.GetFaceEdges (i,faces);	      
		
	
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
    // cout << "sucess " << endl;
    return & table;

  }







  // ******************** //
  // Direct Solver Clusters
  // 0) none
  // 1) low order dofs  --  default


  Array<int> * HDivHighOrderFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    Array<int> & clusters = *new Array<int> (ndof);

    int clustertype = int(precflags.GetNumFlag("ds_cluster",1)); 
    cout << " DirectSolverCluster Clustertype " << clustertype << endl; 
  
    Array<int> vnums,elnums; 
    Array<int> orient; 
  
    Array<int> edges(3);
        
    int nfa = ma.GetNFaces();

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
    return &clusters;

  }

  void HDivHighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize(0); 
  }
  
  void HDivHighOrderFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  { 
    dnums.SetSize(0);
    if(ma.GetDimension() == 3 || discont) return; 

    dnums += ednr;
    dnums += GetFacetDofs (ednr);
  }

  void HDivHighOrderFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    if(ma.GetDimension() == 2 || discont) return; 
   
    dnums += fanr;
    dnums += GetFacetDofs (fanr);
  }

  void HDivHighOrderFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    dnums += GetElementDofs (elnr);
  }

  // h(div) order=0 --->  RT0 was created
  static RegisterFESpace<HDivHighOrderFESpace> init ("hdivho");
}



