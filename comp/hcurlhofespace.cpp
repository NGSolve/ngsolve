
/**
   High Order Finite Element Space for H(Curl) 
*/
#include <comp.hpp>
#include <../fem/hcurlhofe.hpp> 
#include <multigrid.hpp>

#ifdef PARALLEL
#include <parallelngs.hpp>
#endif

#ifdef PARALLEL
extern MPI_Comm MPI_HIGHORDER_COMM;
#endif

namespace ngcomp 
{
  using namespace ngcomp; 

#ifdef PARALLEL
  using namespace ngparallel;
#endif

  HCurlHighOrderFESpace ::  
  HCurlHighOrderFESpace (const MeshAccess & ama, const Flags & aflags, bool parseflags)
    : FESpace (ama, aflags),
      flags(aflags)
  {
    name="HCurlHighOrderFESpace(hcurlho)";
    
    // define flags
    DefineDefineFlag("hcurlho");
    DefineNumFlag("face");
    DefineNumListFlag("gradientdomains");
    DefineNumListFlag("gradientboundaries");
    DefineDefineFlag("nograds");
    DefineDefineFlag("variableorder");
    DefineNumFlag("orderinner");
    DefineNumFlag("orderface");
    DefineNumFlag("orderedge");
    DefineNumFlag("relorder");
    DefineDefineFlag("print"); 
    DefineDefineFlag("noprint"); 
    DefineNumFlag("augmented");
    DefineDefineFlag("fast"); 
    DefineDefineFlag ("discontinuous");
    
    if(parseflags) CheckFlags(flags);
    
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
	  cerr << " WARNING: HCurlHoFeSpace: inconsistent flags: variableorder, order and relorder "
	       << "-> variable order space with rel_order, but order is ignored " << endl; 
	else 
	  { 
	    var_order = 0; 
	    cerr << " WARNING: HCurlHoFeSpace: inconsistent flags: order and rel_order "
		 << "-> uniform order space with order " << endl; 
	  }
      }
      
    int ndom = ma.GetNDomains();
    gradientdomains.SetSize (ndom);
    gradientdomains.Set();
    
    gradientboundaries.SetSize (ma.GetNBoundaries());
    gradientboundaries.Set();

    fn=int(flags.GetNumFlag("face",1));

    if (flags.NumListFlagDefined("gradientdomains"))
      {
	const Array<double> & graddomains = flags.GetNumListFlag ("gradientdomains");
	for (int i = 0; i < gradientdomains.Size(); i++)
	  if (!graddomains[i])
	    gradientdomains.Clear(i);
	cout << endl << " ******* gradientdomains" << gradientdomains << endl;
	
      }
    if (flags.NumListFlagDefined("gradientboundaries"))
      {
	const Array<double> & gradbounds = flags.GetNumListFlag ("gradientboundaries");
	for (int i = 0; i < gradientboundaries.Size(); i++)
	  if (!gradbounds[i])
	    gradientboundaries.Clear(i);
	cout << endl << " ******* gradientboundaries" << gradientboundaries << endl;
      }
    
    if(flags.GetDefineFlag("nograds"))
      {
	gradientdomains.Clear();   
	gradientboundaries.Clear();   
      } 
       
    fast_pfem = flags.GetDefineFlag ("fast");
    print=flags.GetDefineFlag("print"); 
    discontinuous = flags.GetDefineFlag ("discontinuous");


    Flags loflags;
    loflags.SetFlag ("order", 1);
    loflags.SetFlag ("dim", dimension);
    if (iscomplex) loflags.SetFlag ("complex");
    if (discontinuous) loflags.SetFlag ("disontinuous");
    if (flags.NumListFlagDefined ("dirichlet")) 
      loflags.SetFlag ("dirichlet", flags.GetNumListFlag ("dirichlet"));

    
#ifndef PARALLEL
    low_order_space = new  NedelecFESpace (ma, loflags);
#else
    low_order_space = new  ParallelNedelecFESpace (ma, loflags);
#endif

    if (low_order_space)
      prol = new ngmg::EdgeProlongation (*static_cast<NedelecFESpace*> (low_order_space));
    
 
   
    uniform_order_inner = int (flags.GetNumFlag ("orderinner", -1));
    uniform_order_face = int (flags.GetNumFlag ("orderface", -1));
    uniform_order_edge = int (flags.GetNumFlag ("orderedge", -1));

    
        
    if (flags.NumFlagDefined("smoothing")) 
      throw Exception ("Flag 'smoothing' for fespace is obsolete \n Please use flag 'blocktype' in preconditioner instead");
    if (flags.NumFlagDefined("cluster")) 
      throw Exception ("Flag 'cluster' for fespace is obsolete \n Please use flag 'ds_cluster' in preconditioner instead");
       
    augmented = int (flags.GetNumFlag ("augmented", 0));

    /*    if (flags.GetDefineFlag ("optext"))
	  { 
	  trig = new HCurlHighOrderTrig<TrigExtensionMin> (order);
	  tet =  new HCurlHighOrderTet<TrigExtensionMin> (order);
	  prism =  new HCurlHighOrderPrism<TrigExtensionMin> (order);
	  }
	  else
    */


    // Evaluator 
    static ConstantCoefficientFunction one(1);
    if (ma.GetDimension() == 2)
      {
	Array<CoefficientFunction*> coeffs(1);
	coeffs[0] = &one;
	evaluator = GetIntegrators().CreateBFI("massedge", 2, coeffs);
	if ( !discontinuous )
	  boundary_evaluator = GetIntegrators().CreateBFI("robinedge",2,coeffs); 
      }
    else if(ma.GetDimension() == 3) 
      {
	Array<CoefficientFunction*> coeffs(1); 
	coeffs[0] = &one;
	evaluator = GetIntegrators().CreateBFI("massedge",3,coeffs); 
	if ( !discontinuous )
	  boundary_evaluator = GetIntegrators().CreateBFI("robinedge",3,coeffs); 
	
      }

  }
  
  HCurlHighOrderFESpace :: ~HCurlHighOrderFESpace ()
  {
    ;
  }
  
  FESpace * HCurlHighOrderFESpace :: 
  Create (const MeshAccess & ma, const Flags & flags)
  {
    return new HCurlHighOrderFESpace (ma, flags, true);
  }
  
  void HCurlHighOrderFESpace :: Update(LocalHeap & lh)
  {
    FESpace :: Update (lh);


    const int dim = ma.GetDimension(); 

    if (order < 0) 
      throw Exception("HCurlHighOrderFESpace::Update() order < 0 !" ) ;
    
    if (low_order_space)
      low_order_space -> Update(lh);

    nv = ma.GetNV();
    nel = ma.GetNE();
    int nse = ma.GetNSE();
    maxorder = -1; 
    minorder = 99; 

    for(int i=0; i<specialelements.Size(); i++)
      delete specialelements[i];
    specialelements.DeleteAll();

    switch(ma.GetDimension())
      {
      case 2: 
	ned = ma.GetNEdges();
	nfa = 0;  
	break; 
      case 3: 
	ned = ma.GetNEdges();
	nfa = ma.GetNFaces();
	break; 
      }

    order_edge.SetSize (ned);   
    order_face.SetSize (nfa); 
    order_inner.SetSize (nel);
    usegrad_edge.SetSize (ned);                                
    usegrad_face.SetSize (nfa); 
    usegrad_cell.SetSize (nel);
    fine_edge.SetSize(ned); 
    fine_face.SetSize(nfa); 

    int p = var_order ? 0 : order; 
    order_edge = p; 
    order_inner = INT<3> (p,p,p); 

 
   
    fine_edge = 0; 
    if(dim==3) 
      { 
	fine_face = 0; 
	order_face = INT<2> (p,p);
	usegrad_face = 0; 
      } 
     
    usegrad_edge = 0;                                
    usegrad_cell = 0; 

    Array<int> eledges, elfaces, vnums;


    for(int i=0; i< nel; i++) 
      if(gradientdomains[ma.GetElIndex(i)]) 
	{
	  ma.GetElEdges(i,eledges);
	  
	  for(int j=0;j<eledges.Size();j++)
	    usegrad_edge[eledges[j]]=1;
	  
	  if(ma.GetDimension()==3)
	    {
	      ma.GetElFaces(i,elfaces);
	      
	      for(int j=0;j<elfaces.Size();j++)
		usegrad_face[elfaces[j]]=1; 
	    }
	  usegrad_cell[i] = 1; 
	} 

    for(int i=0; i<nse && gradientboundaries.Size(); i++)
      if(gradientboundaries[ma.GetSElIndex(i)])
	{
	  ma.GetSElEdges(i,eledges);
	  for(int j=0; j<eledges.Size();j++)
	    usegrad_edge[eledges[j]] = 1;
	
	  if(ma.GetDimension()==3)
	    usegrad_face[ma.GetSElFace(i)] = 1;
	}

      
	
    for (int i = 0; i < nel; i++)
      {
	int index = ma.GetElIndex(i);
	if (!DefinedOn (index)) continue;
	
	INT<3> el_orders = ma.GetElOrders(i);
	
	ELEMENT_TYPE eltype=ma.GetElType(i); 
	const FACE * faces = ElementTopology::GetFaces (eltype);
	const EDGE * edges = ElementTopology::GetEdges (eltype);
	const POINT3D * points = ElementTopology :: GetVertices (eltype);
	ma.GetElVertices (i, vnums);
	
	ma.GetElEdges (i, eledges);		
	for(int j=0;j<eledges.Size();j++) fine_edge[eledges[j]] = 1; 
	
	if (dim == 3)
	  {
	    ma.GetElFaces(i, elfaces);
	    for(int j=0;j<elfaces.Size();j++) fine_face[elfaces[j]] = 1; 
	  }

	for(int j=0;j<3;j++)
	  {
	    el_orders[j] = el_orders[j]+rel_order;
	    if(el_orders[j] > maxorder) maxorder=el_orders[j];
	    if(el_orders[j] < minorder) minorder=el_orders[j];
	  }

	
	if(!var_order) continue; // for fixed order find only fine-edges-faces! 
	
	     

	for(int j=0;j<eledges.Size();j++)
	  for(int k=0;k<dim;k++)
	    if(points[edges[j][0]][k] != points[edges[j][1]][k]) // find edge-dir x,y,z
	      { 
		order_edge[eledges[j]] = max2(el_orders[k],order_edge[eledges[j]]);
		break; 
	      }
	
	for(int j=0;j<dim;j++)
	  order_inner[i][j] = int(max2(order_inner[i][j],el_orders[j]));
	
	if(dim==2) continue; 
	
	for(int j=0;j<elfaces.Size();j++)
	  {
	    if(faces[j][3]==-1) 
	      {
		
		order_face[elfaces[j]][0] = int(max2(order_face[elfaces[j]][0], 
						     el_orders[0]));
		order_face[elfaces[j]][1] = order_face[elfaces[j]][0]; 
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
			order_face[elfaces[j]][l] = int(max2(order_face[elfaces[j]][l], 
							     el_orders[k]));
			break; 
		      } 
	      }  
	  }
      }
      	
    if(!var_order) { maxorder = order; minorder = order;} 
    
    if(uniform_order_inner>-1) 
      order_inner = INT<3> (uniform_order_inner,uniform_order_inner,uniform_order_inner); 
    if(uniform_order_edge>-1) 
      order_edge = uniform_order_edge; 
    if(uniform_order_face>-1 && dim == 3) 
      order_face = INT<2> (uniform_order_face, uniform_order_face);
	
    // order of FINE FACES and EDGES for savety reasons set to 0 
    for(int i=0;i<order_edge.Size();i++) 
      if(!fine_edge[i]) order_edge[i] = 0;  
    
    for(int i=0;i<order_face.Size();i++) 
      if(!fine_face[i]) order_face[i] = INT<2> (0,0);  



    UpdateDofTables(); 
    UpdateCouplingDofArray();
    FinalizeUpdate (lh);

    if (timing) Timing();    

#ifdef PARALLEL
    UpdateParallelDofs();
#endif
  }
		
  void HCurlHighOrderFESpace :: UpdateDofTables()
  {
    nv = ma.GetNV();
    ned = ma.GetNEdges();
    
    nfa = (ma.GetDimension()==3) ? ma.GetNFaces() : 0;
    nel = ma.GetNE();
    Array<int> pnums; 
    int i; 
    
    nedfine = 0; 
    for(int i=0; i<ned; i++) 
      if(fine_edge[i] == 1) nedfine++; 

    //ndof = nedfine; 
    ndof = ned; // Nedelec (order = 0) !!   
    if(augmented == 1 ) ndof+=nv; 
       
    first_edge_dof.SetSize(ned+1); 
    for (i = 0; i < ned; i++)
      {
	first_edge_dof[i] = ndof;
	if(order_edge[i]>0)
	  ndof += usegrad_edge[i]*order_edge[i];
      }
    first_edge_dof[ned] = ndof;
    // cout << " after edges ndof is " << ndof << endl; 
    
    first_face_dof.SetSize(nfa+1);
    face_ngrad.SetSize(nfa); 
    face_ngrad = 0; 
    for (i=0; i< nfa; i++) 
      { 
	first_face_dof[i] = ndof; 
	ma.GetFacePNums(i,pnums);  
	INT<2> p = order_face[i]; 
	switch(pnums.Size())
	  {
	  case 3: //Triangle   
	    if(p[0]>1)
	      {
		ndof += ((usegrad_face[i]+1)*p[0] + 2)*(p[0]-1)/2;
		face_ngrad[i] = usegrad_face[i] *p[0]*(p[0]-1)/2;
	      }
	    break; 
	  case 4: //Quad 
	    //if(p[0]>0 && p[1]>0)
	    {
	      ndof += (usegrad_face[i]+1)*p[0]*p[1] + p[0] + p[1]; 
	      face_ngrad[i] = usegrad_face[i]*p[0]*p[1];; 
	      break; 
	    }
	  }
      } 
    first_face_dof[nfa] = ndof;   
    // cout << " after faces ndof is " << ndof << endl; 
    
    cell_ngrad.SetSize(nel); 
    cell_ngrad = 0; 
    first_inner_dof.SetSize(nel+1);
    for (i=0; i< nel; i++)
      {
	first_inner_dof[i] = ndof;
	INT<3> p = order_inner[i];
	switch(ma.GetElType(i)) 
	  {
	  case ET_TRIG:
	    if(p[0]>1)
	      {
		ndof += ((usegrad_cell[i]+1)*p[0] + 2) * (p[0]-1) /2;
		cell_ngrad[i] = ((usegrad_cell[i])*p[0]) * (p[0]-1) /2;
	      }
	    break; 
	  case ET_QUAD: 
	    if(p[0]>=0 &&  p[1]>=0) 
	      {
		ndof += (usegrad_cell[i]+1) * p[0] * p[1] + p[0] + p[1]; 
		cell_ngrad[i] = (usegrad_cell[i]) * p[0] * p[1];
	      }
	    break; 
	  case ET_TET: 
	    if(p[0]>2)
	      { 
		ndof += ((usegrad_cell[i] + 2) *  p[0] + 3) * (p[0]-2) * (p[0]-1) / 6; 
		cell_ngrad[i] = ((usegrad_cell[i] ) *  p[0]) * (p[0]-2) * (p[0]-1) / 6; 
	      }
	    break; 
	  case ET_PRISM:
	    if(p[0]>1 &&  p[2]>=0) 
	      {
		ndof += (usegrad_cell[i]+2)*p[2] * p[0]*(p[0]-1)/2
		  + (p[0]-1)*p[2] + p[0]*(p[0]-1)/2; 
		cell_ngrad[i] = (usegrad_cell[i])*p[2] * p[0]*(p[0]-1)/2;
	      } 
	    break; 
	  case ET_HEX:
	    if(p[0]>=0 && p[1]>=0 && p[2] >=0)
	      {
		ndof += (usegrad_cell[i] + 2)* p[0]*p[1]*p[2] +  p[0]*p[1] 
		  + p[0]*p[2] + p[1]*p[2];
		cell_ngrad[i] = (usegrad_cell[i])* p[0]*p[1]*p[2]; 
	      }
	    break;  
	  case ET_PYRAMID:
	    if (p[0] > 1)
	      {
		ndof += usegrad_cell[i]*(p[0]-1)*p[0]*(2*p[0]-1)/6 + p[0]*(2*p[0]*p[0]+3*p[0]-2)/3; 
		cell_ngrad[i] = usegrad_cell[i]*(p[0]-1)*p[0]*(2*p[0]-1)/6; 
	      }
	    break; 
          default:  // for the compiler
            break; 
	  }
      }
    first_inner_dof[nel] = ndof;    
    // cout << " after elements ndof is " << ndof << endl; 

    if ( discontinuous )
      {
	Array<int> edges, faces;

	ndof = 0;
	for ( int el = 0; el < nel; el++ )
	  {
	    int ndof_inner = first_inner_dof[el+1] - first_inner_dof[el];
	    first_inner_dof[el] = ndof;
	    ma.GetElEdges(el, edges);
	    if ( ma.GetDimension() == 3 )
	      ma.GetElFaces(el, faces);

	    for ( int ed = 0; ed < edges.Size(); ed++ )
	      {
		int ndof_edge = 1 + first_edge_dof[edges[ed]+1] - first_edge_dof[edges[ed]];
		ndof += ndof_edge;
	      }

	    if ( ma.GetDimension() == 3 )
	      for ( int fa = 0; fa < faces.Size(); fa++ )
		{
		  int ndof_face = first_face_dof[faces[fa]+1] - first_face_dof[faces[fa]];
		  ndof += ndof_face;
		}
	    ndof += ndof_inner;
	  }
	first_inner_dof[nel] = ndof;

	first_edge_dof = 0;
	first_face_dof = 0;

      }


    if(print)
      {
	(*testout) << " HCURLHO " << endl; 
	(*testout) << "augmented " << augmented << endl; 
	(*testout) << "ndof  = " << ndof << endl; 
	(*testout) << "variableorder  = " << var_order << endl; 
	(*testout) << "order = " << order << " , relorder = " << rel_order << endl;

	// if(var_order)
	{
	  (*testout) << "order_edge  (hcurlho) = " << order_edge << endl;
	  (*testout) << "order_face  (hcurlho) = " << order_face << endl;
	  (*testout) << "order_inner (hcurlho) = " << order_inner << endl;	
	} 
	
	(*testout) << "hcurlho first_edge   = " << first_edge_dof << endl;
	(*testout) << "hcurlho first_face   = " << first_face_dof << endl;
	(*testout) << "hcurlho first_inner  = " << first_inner_dof << endl;

	*testout << " hcurlho fine_edge = " << fine_edge << endl; 
	*testout << " hcurlho fine_face = " << fine_face << endl; 
	
	
      } 
  }

  void  HCurlHighOrderFESpace :: UpdateCouplingDofArray ()
  {
    
    ctofdof.SetSize(ndof);
    ctofdof = WIREBASKET_DOF;
    
    if(discontinuous) 
    {
      ctofdof = LOCAL_DOF;
      return;
    } 
    int first, next;
    for (int edge = 0; edge < ma.GetNEdges(); edge++) {
      ctofdof[edge] = WIREBASKET_DOF; //Nedelec0
      IntRange range = GetEdgeDofs (edge);
      first = range.First();
      next = range.Next();
      for (int j=first; j<next;j++)
	ctofdof[j] = WIREBASKET_DOF;
    }
      
    // faces 
    if (ma.GetDimension() == 3)
      for (int face = 0; face < ma.GetNFaces(); face++){
	IntRange range = GetFaceDofs (face);
	first = range.First();
	next = range.Next();
	for (int j=first;j<next;j++)
	  ctofdof[j] = INTERFACE_DOF; //NOGRAD / GRAD
      }      
    
    for (int el = 0; el < ma.GetNE(); el++){
      IntRange range = GetElementDofs (el);
      for (int j=range.First();j<range.Next();j++)
	ctofdof[j] = LOCAL_DOF;
    }
//     *testout << "ctofdof = \n" << ctofdof << endl;
  }

  
  const FiniteElement & HCurlHighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    FiniteElement * fe = 0;
    
    switch (ma.GetElType(elnr))
      {
      case ET_TET:     fe = new (lh) HCurlHighOrderFE<ET_TET> (); break;
      case ET_PYRAMID: fe = new (lh) HCurlHighOrderFE<ET_PYRAMID> (); break;
      case ET_PRISM:   fe = new (lh) HCurlHighOrderFE<ET_PRISM> (); break;
      case ET_TRIG:    fe = new (lh) HCurlHighOrderFE<ET_TRIG> (); break;
      case ET_QUAD:    fe = new (lh) HCurlHighOrderFE<ET_QUAD> (); break;
      case ET_HEX:     fe = new (lh) HCurlHighOrderFE<ET_HEX> (); break; 

      default:
	stringstream str;
	str << "HCurlHighOrderFESpace " << GetClassName() 
	    << ", undefined eltype " 
	    << ElementTopology::GetElementName(ma.GetElType(elnr))
	    << ", order = " << order << endl;
	throw Exception (str.str());
      }
    

    ArrayMem<int,12> vnums;
    ma.GetElVertices(elnr, vnums);

    if(ma.GetDimension() == 2) 
      {	
	HCurlHighOrderFiniteElement<2> * hofe = 
	  static_cast<HCurlHighOrderFiniteElement<2>*> (fe);

               
    	ArrayMem<int, 4> ednums, ord_edge, ug_edge;
	
	ma.GetElEdges(elnr, ednums);
	
	ord_edge.SetSize (ednums.Size());
	ug_edge.SetSize (ednums.Size()); 
	
	for (int j = 0; j < ednums.Size(); j++)
	  {
	    ord_edge[j] = order_edge[ednums[j]];
	    ug_edge[j]  = usegrad_edge[ednums[j]]; 
	  }

	hofe -> SetOrderEdge (ord_edge);
        hofe -> SetOrderCell (order_inner[elnr]);   // old style
        INT<2> p(order_inner[elnr][0], order_inner[elnr][1]);
        FlatArray<INT<2> > of(1, &p);
        hofe -> SetOrderFace (of);

	hofe -> SetUsegradEdge (ug_edge); 
	hofe -> SetUsegradCell (usegrad_cell[elnr]);  // old style
        FlatArray<int> augf(1,&usegrad_cell[elnr]);
	hofe -> SetUsegradFace (augf); 
	hofe -> ComputeNDof();

#ifdef PARALLEL
	if ( ntasks > 1 )
	  for ( int i = 0; i < vnums.Size(); i++ )
	    vnums[i] = parallelma->GetDistantPNum(0, vnums[i]);
#endif
	hofe -> SetVertexNumbers (vnums);
      }   

    else if (ma.GetDimension() == 3) 

      {
        Ng_Element ngel = ma.GetElement<3> (elnr);

	HCurlHighOrderFiniteElement<3> * hofe = 
	  static_cast<HCurlHighOrderFiniteElement<3>*> (fe);

#ifdef PARALLEL
	if ( ntasks > 1 )
	for ( int i = 0; i < vnums.Size(); i++ )
	  vnums[i] = parallelma->GetDistantPNum(0, vnums[i]);
#endif

        for (int i = 0; i < ngel.vertices.Size(); i++)
          hofe -> SetVertexNumber (i, ngel.vertices[i]);
        
        for (int i = 0; i < ngel.edges.Size(); i++)
          {
            hofe -> SetOrderEdge (i, order_edge[ngel.edges[i]]);
            hofe -> SetUseGradEdge (i, usegrad_edge[ngel.edges[i]]);
          }

        for (int i = 0; i < ngel.faces.Size(); i++)
          {
            hofe -> SetOrderFace (i, order_face[ngel.faces[i]]);
            hofe -> SetUseGradFace (i, usegrad_face[ngel.faces[i]]);
          }

	hofe -> SetOrderCell (order_inner[elnr]);
	hofe -> SetUsegradCell (usegrad_cell[elnr]); 

	hofe -> ComputeNDof();
	hofe -> SetDiscontinuous(discontinuous);
      }

    return *fe;
  }
 
  const FiniteElement & HCurlHighOrderFESpace :: GetSFE (int selnr, LocalHeap & lh) const
  {
    FiniteElement * fe = 0;
    
    if ( discontinuous )
      {
	switch (ma.GetSElType(selnr))
	  {
	  case ET_SEGM: fe = new (lh) DummyFE<ET_SEGM>; break; 
	  case ET_TRIG: fe = new (lh) DummyFE<ET_TRIG>; break;
	  case ET_QUAD: fe = new (lh) DummyFE<ET_QUAD>; break;
	  default:
	    fe = 0;
	  }
	return *fe;
      }

    switch (ma.GetSElType(selnr))
      {
      case ET_SEGM: fe = new (lh) HCurlHighOrderFE<ET_SEGM> (); break; 
      case ET_TRIG: fe = new (lh) HCurlHighOrderFE<ET_TRIG> (); break;
      case ET_QUAD: fe = new (lh) HCurlHighOrderFE<ET_QUAD> (); break;
      default:
	fe = 0;
      }

    if (!fe)
      {
	stringstream str;
	str << "HCurlHighOrderFESpace " << GetClassName() 
	    << ", undefined eltype " 
	    << ElementTopology::GetElementName(ma.GetSElType(selnr))
	    << ", order = " << order << endl;
	throw Exception (str.str());
      }

    ArrayMem<int, 8> vnums;
    ArrayMem<int, 4> ednums, ord_edge, ug_edge; 
       
    ma.GetSElVertices(selnr, vnums);

    if(ma.GetSElType(selnr) == ET_SEGM)
      {
	HCurlHighOrderFiniteElement<1> * hofe =
	  dynamic_cast<HCurlHighOrderFiniteElement<1>*> (fe);

#ifdef PARALLEL
	if ( ntasks > 1 )
	  for ( int i = 0; i < vnums.Size(); i++ )
	    vnums[i] = parallelma->GetDistantPNum(0, vnums[i]);
#endif
	hofe -> SetVertexNumbers (vnums);
	ma.GetSElEdges(selnr, ednums);
	hofe -> SetOrderCell (order_edge[ednums[0]]);  // old style
        FlatArray<int> aoe(1, &order_edge[ednums[0]]);
	// hofe -> SetOrderEdge (order_edge[ednums[0]]);
        hofe -> SetOrderEdge (aoe);
	hofe -> SetUsegradCell (usegrad_edge[ednums[0]]);  // old style
        FlatArray<int> auge(1, &usegrad_edge[ednums[0]]);
	hofe -> SetUsegradEdge (auge);
	hofe -> ComputeNDof();
      } 
    else 
      {     
	HCurlHighOrderFiniteElement<2> * hofe =
	  dynamic_cast<HCurlHighOrderFiniteElement<2>*> (fe);
      
	ma.GetSElEdges(selnr, ednums);
	
	ord_edge.SetSize (ednums.Size());
	ug_edge.SetSize (ednums.Size()); 
	
	for (int j = 0; j < ednums.Size(); j++)
	  {
	    ord_edge[j] = order_edge[ednums[j]]; 
	    ug_edge[j] = usegrad_edge[ednums[j]];
	  }
	
#ifdef PARALLEL
	if ( ntasks > 1 )
	for ( int i = 0; i < vnums.Size(); i++ )
	  vnums[i] = parallelma->GetDistantPNum(0, vnums[i]);
#endif
	hofe -> SetVertexNumbers (vnums);
	hofe -> SetOrderEdge (ord_edge);

        INT<2> p = order_face[ma.GetSElFace(selnr)];
	hofe -> SetOrderCell (INT<3> (p[0],p[1],0));  // old style
        FlatArray<INT<2> > of(1, &p);
        hofe -> SetOrderFace (of);

	hofe -> SetUsegradEdge(ug_edge); 
        FlatArray<int> augf(1, &usegrad_face[ma.GetSElFace(selnr)]);
	hofe -> SetUsegradFace(augf); 
	hofe -> SetUsegradCell(usegrad_face[ma.GetSElFace(selnr)]);   // old style

	hofe -> ComputeNDof();
      }
    
    return *fe;
  }


  int HCurlHighOrderFESpace :: GetNDof () const
  {
    return ndof;
  }

  void HCurlHighOrderFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    // ordering of shape functions
    // (1*e1),.. (1*e_ne)  
    // (p_t)*tang_e1, ... p_t*tang_ne
    // p_n * el1, p_i * el1, ... , p_n * nel_n , p_i *el_nel 

    Ng_Element ngel = ma.GetElement (elnr);
    dnums.SetSize(0);
     
      //Nedelec0
    if ( !discontinuous )
      for (int i = 0; i < ngel.edges.Size(); i++) 
	dnums.Append (ngel.edges[i]);
        
    //edges
    for (int i = 0; i < ngel.edges.Size(); i++) 
      dnums += GetEdgeDofs(ngel.edges[i]); 
       
    // faces 
    if (ma.GetDimension() == 3)
      for(int i = 0; i < ngel.faces.Size(); i++)     
        dnums += GetFaceDofs(ngel.faces[i]);  
        
    dnums += GetElementDofs (elnr); 
    
    if (!DefinedOn (ma.GetElIndex (elnr)))
      dnums = -1;
  }
 



  void HCurlHighOrderFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    Array<int> vnums, ednums;
    int fnum; 
    int i, j;
    int first,next; 

    ma.GetSElEdges (selnr, ednums);
    fnum = ma.GetSElFace (selnr); 

    dnums.SetSize(0);

    if(order <0) throw (" HCurlHighOrderFESpace :: GetSDofNrs() order < 0 "); 

    if ( !discontinuous )
      //Nedelec
      for (i = 0; i < ednums.Size(); i++) 
	dnums.Append (ednums[i]);

    if(augmented==1) 
      {
	
	ma.GetSElVertices (selnr, vnums);

	for(i=0;i< vnums.Size();i++) 
	  dnums.Append(ned+vnums[i]); 
      } 
    
    for (i = 0; i < ednums.Size(); i++)
      {
	first = first_edge_dof[ednums[i]];
	next = first_edge_dof[ednums[i]+1]; 
	for (j = first; j <next; j++)
	  dnums.Append (j);
      }

    if(first_face_dof.Size()>1)		       
      {
	// inner = normal + bubbles 
	first = first_face_dof[fnum];
	next = first_face_dof[fnum+1]; 
	for(j=first; j<next; j++) 
	  dnums.Append(j);
      }
    
    if (!DefinedOnBoundary (ma.GetSElIndex (selnr)))
      dnums = -1;

    // *testout << " GetSDNums selnr " << selnr << " \t " << dnums << endl; 
  
  }


  void HCurlHighOrderFESpace ::
  SetGradientDomains (const BitArray & adoms)
  {
    gradientdomains = adoms;
  }


  void HCurlHighOrderFESpace ::
  SetGradientBoundaries (const BitArray & abnds)
  {
    gradientboundaries = abnds;
  }

  
  void HCurlHighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
  }

  void HCurlHighOrderFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    if ( discontinuous ) return;

    dnums.Append(ednr);
    int first = first_edge_dof[ednr];
    int next = first_edge_dof[ednr+1]; 
    for (int j = first; j <next; j++)
      dnums.Append (j);
  }

  void HCurlHighOrderFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    if ( discontinuous ) return;

    int first = first_face_dof[fanr];
    int next = first_face_dof[fanr+1]; 
    for (int j = first; j <next; j++)
      dnums.Append (j);
  }

  void HCurlHighOrderFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    int first = first_inner_dof[elnr];
    int next = first_inner_dof[elnr+1]; 
    for (int j = first; j <next; j++)
      dnums.Append (j);  
  }


  Table<int> * HCurlHighOrderFESpace :: 
  CreateSmoothingBlocks (const Flags & precflags) const
  {
    bool eliminate_internal = precflags.GetDefineFlag("eliminate_internal");
    int i, j, k, first,ii;
    int ncnt; 
    int ni = nel;
    if (eliminate_internal) ni = 0; 
    
    int SmoothingType = int(precflags.GetNumFlag("blocktype",2)); 
    bool excl_grads = precflags.GetDefineFlag("exclude_grads"); 
    cout << " EXCLUDE GRADS " << excl_grads << endl; 
    
    Array<int> vnums,elnums; 
    Array<int> orient; 
    Array<int> ednums, fanums, enums, f2ed;
    
    // int augv = augmented; 

        
    if(nfa == 0 && SmoothingType == 1) 
      SmoothingType = 4; 
    if(nfa == 0 && (SmoothingType == 2 || SmoothingType ==3)) 
      SmoothingType = 5;






    if (precflags.GetDefineFlag("subassembled"))
      {
	
	TableCreator<int> creator;
	for ( ; !creator.Done(); creator++)
	  {

	    if (creator.GetMode() == 1)
	      cout << "High order AFW blocks " << endl;
		
	    for (int i = 0; i < ned; i++)
	      if (!IsDirichletEdge(i) && fine_edge[i])
		{
		  Ng_Node<1> edge = ma.GetNode<1> (i);
		  for (int k = 0; k < 2; k++)
		    {
		      creator.Add (edge.vertices[k], i);
		      creator.Add (edge.vertices[k], GetEdgeDofs(i));
		    }
		}
	  }

	
	return creator.GetTable();
      }









    cout << "SmoothingType " << SmoothingType << endl; 
    cout << " Use H(Curl)-Block smoothing " ;
    switch(SmoothingType) 
      {
      case 4: // former 11 
	cout << " AFW(loE) + E + F + I (noClusters)" << endl; 
	ncnt = nv  + ned + nfa + ni + augmented*nv;
	break;  
      case 5: //	
	cout << " AFW(hoE) + E + F + I (noClusters)" << endl; 
	ncnt = nv  + ned + nfa + ni;
	break; 	
      case 1: // former 10 
	cout << " Line-Clustering and  AFW(loE)  " << endl; 
	ncnt = nv + ned + nfa + ni; 
	break;
      case 2: // former 26 
	cout << " Line-Clustering and AFW(hoE)  " << endl; 
	ncnt = nv + ned + nfa + ni;
	break;
      case 3: // former 23 (equal to HCurl-Old-SmoothingBlocks)  
	cout << " Line-Clustering: AFW(hoE) - (horizE)F - (trigF)I " << endl; 
	// horiz-quadface-stack(horiz-edge,top-bot-trig-Faces, inner)" << endl; 
	ncnt = nv + ned + nfa +ni; 
	break;
      case 6: // der hier is nur fuer testzwecke !!!  
	cout << "Jacobi (Diag)" << endl; 
	ncnt = nv  + ned + nfa + ni;
	break;  
      case 7: 
	cout << "EDGE-blocks (noClusters)" << endl; 
	ncnt = nv  + ned;
	break; 
      case 21:
	cout << "wave equation blocks" << endl;
	ncnt = ned + nfa;
	break;
      default: 
	throw Exception("HCurlHoFeSpace:: CreateSmoothingBlocks chosen Blocktype not valid \n Choose blocktype 1-6 "); 
	ncnt =0; 
	return 0; 
      } 

        
    Array<int> cnt(ncnt); 
    cnt = 0;
    ii=0; 

    switch (SmoothingType)
      {
      case 4: // AFW(loE) + E + F + I (noCluster)
	{
	  for (i = 0; i < ned; i++)
	    if(fine_edge[i])
	      {		
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);
		cnt[pn1]++;
		cnt[pn2]++;
	      }
	  int nnv = nv; 
	  if(augmented==1 && !excl_grads)
	    { 
	      nnv += nv; 
	      for (i=0;i<nv;i++)
		cnt[nv+i] = 1; 
	    } 
	 
	  for(i=0;i<ned;i++)
	    cnt[nnv+i] =  (first_edge_dof[i+1] - first_edge_dof[i])*(1-excl_grads);
	  for (i = 0; i < nfa; i++)
	    cnt[nnv+ned+i] = first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i];
	  for (i = 0; i < ni; i++)
	    cnt[nnv+ned+nfa+i] = first_inner_dof[i+1] - first_inner_dof[i] - excl_grads*cell_ngrad[i];
	  break;
	}
      case 5: // AFW(hoE) + E + F + I (noCluster)
	{
	  for (i = 0; i < ned; i++)
	    if(fine_edge[i])
	      {
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);
		cnt[pn1]+= 1 + (first_edge_dof[i+1]-first_edge_dof[i])*(1-excl_grads);
		cnt[pn2]+= 1 + (first_edge_dof[i+1]-first_edge_dof[i])*(1-excl_grads);
	      }
	  int nnv = nv; 
	  /* 
	     if(augmented==1 && wgrads)
	     { 
	     nnv += nv; 
	     for (i=0;i<nv;i++)
	     cnt[nv+i] = 1; 
	     } */ 
	  
	  for (i = 0; i < nfa; i++)
	    cnt[nnv+ned+i] = first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i];
	  for (i = 0; i < ni; i++)
	    cnt[nnv+ned+nfa+i] = first_inner_dof[i+1] - first_inner_dof[i] - excl_grads*cell_ngrad[i];
	  break;
	}

      case 1: // Clustering with AFW(loE) 
	{
	  
	  for (i = 0; i < ned; i++)
	    if(fine_edge[i])
	      {
		int pn1, pn2;
	 
		ma.GetEdgePNums (i,pn1,pn2);
		pn1 = ma.GetClusterRepVertex(pn1);
	  
		pn2 = ma.GetClusterRepVertex(pn2);
		cnt[pn1]++;
		if (pn1 != pn2)
		  cnt[pn2]++;
		// ev. noch lo-edge zu horiz.edges 
		cnt[ma.GetClusterRepEdge(i)] +=
		  (first_edge_dof[i+1] - first_edge_dof[i])*(1-excl_grads);
	      }
	 
	  for (i = 0; i < nfa; i++)
	    if(fine_face[i])
	      cnt[ma.GetClusterRepFace(i)] += 
		first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i];
	 
	  for (i = 0; i < ni; i++)
	    cnt[ma.GetClusterRepElement(i)] += 
	      first_inner_dof[i+1] - first_inner_dof[i] - excl_grads*cell_ngrad[i];
	  
	  break;
	}
      case 2: // Clustering with AFW(hoE)
	{
	  // Stack-AFW(hoEdges) and Stack of horiz edges 
	  for(i=0;i<ned;i++) 
	    if(fine_edge[i] && !IsDirichletEdge(i))
	      {
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);
		pn1 = ma.GetClusterRepVertex(pn1);
		pn2 = ma.GetClusterRepVertex(pn2); 
		int nde = 1+ (1-excl_grads)*(first_edge_dof[i+1] - first_edge_dof[i]); 
		cnt[pn1] += nde; 
		if (pn1 != pn2)
		  { 
		    cnt[pn2]+= nde;
		    cnt[ma.GetClusterRepEdge(i)] +=nde;
		  }
	      }
	  
	  // Stack of horizontal edges: + quadfaces and prism faces 
	  //  (+ inner) 
	  for(i=0;i<nfa;i++) 
	    if(fine_face[i] && !IsDirichletFace(i)) 
	      cnt[ma.GetClusterRepFace(i)] += first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i]; 
	     	 
	  for(i=0;i<ni;i++) 
	    { 
	      int ccl = ma.GetClusterRepElement(i);
	      // Stack: plane-face + inner 
	      cnt[ccl] += first_inner_dof[i+1] - first_inner_dof[i]- excl_grads*cell_ngrad[i]; 
	    }
	  break;
	}

      case 3: 
	{
	  // Stack-AFW(hoEdges) and Stack of horiz edges 
	  for(i=0;i<ned;i++) 
	    if(fine_edge[i])
	      {
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);
		pn1 = ma.GetClusterRepVertex(pn1);
		pn2 = ma.GetClusterRepVertex(pn2); 
		int nde = 1+ (1-excl_grads)*(first_edge_dof[i+1] - first_edge_dof[i]); 
		cnt[pn1] += nde; 
		if (pn1 != pn2)
		  { 
		    cnt[pn2]+= nde;
		    cnt[ma.GetClusterRepEdge(i)] +=nde;
		  }
	      }
	  
	  // Stack of horizontal edges: + quadfaces and prism faces 
	  //  (+ inner) 
	  for(i=0;i<nfa;i++) 
	    if(fine_face[i]) 
	      cnt[ma.GetClusterRepFace(i)] += first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i]; 
	  
	  for(i=0;i<ni;i++) 
	    { 
	      int ccl = ma.GetClusterRepElement(i);
	      int ndi = first_inner_dof[i+1] - first_inner_dof[i] - excl_grads*cell_ngrad[i]; 
	      // Stack: plane-face + inner 
	      cnt[ccl] += ndi; 
	      
	      
	      
	      // inner to quad-face (horizontal edge stack) 
	      // each face different stack 
	      
	      
	      ma.GetElFaces(i,fanums); 
	      for(j=0;j<fanums.Size();j++) 
		{ 
		  int fcl = ma.GetClusterRepFace(fanums[j]);
		  if(fcl != ccl) 
		    cnt[fcl] += ndi; 
		}  
	      // if nocluster -> FI + I 
	    }  
	  
	  for(i=0;i<nfa;i++) // Trig-Faces to horiz-edge stack 
	    if(fine_face[i]) 
	      { 
		ma.GetFacePNums(i,vnums); 
		if(vnums.Size()==4) continue; 
		int fcl = ma.GetClusterRepFace(i); 
		ma.GetFaceEdges(i,ednums); 
		
		int nd = first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i]; 
		for(j=0;j<ednums.Size();j++)
		  {
		    int ecl = ma.GetClusterRepEdge(ednums[j]); 
		    if(ecl==fcl) continue; 
		    cnt[ecl] += nd; 
		  }
		// if nocluster -> AFW(hoE) + EF + FI + I 
	      } 
	  
	  break;
	}
      case 6: // Jacobi(Diag)
	{
	  int ii=0; 
	  for (i = 0; i < ned; i++)
	    if(fine_edge[i])
	      cnt[ii++]=1; 
	  
	  for(i=0;i<ned;i++)
	    if(fine_edge[i])
	      for(int k= first_edge_dof[i]; k< first_edge_dof[i+1]; k++) 
		cnt[ii++] = 1;
	  
	  for (i = 0; i < nfa; i++)
	    if(fine_face[i])
	      for(int k= first_face_dof[i]; k< first_face_dof[i+1]; k++) 
		cnt[ii++] = 1; 
	  
	  for (i = 0; i < ni; i++)
	    for(int k= first_inner_dof[i]; k< first_inner_dof[i+1]; k++) 
	      cnt[ii++] = 1; 
	  
	  break;
	}
      case 7: 
	{
	  // AFW(ho-E) - Edge-Blocks 
	  for(i=0;i<ned;i++) 
	    if(fine_edge[i])
	      {
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);
		pn1 = ma.GetClusterRepVertex(pn1);
		pn2 = ma.GetClusterRepVertex(pn2); 
		int nde = 1+ (1-excl_grads)*(first_edge_dof[i+1] - first_edge_dof[i]); 
		cnt[pn1] += nde; 
		if (pn1 != pn2)
		  { 
		    cnt[pn2]+= nde;
		    cnt[ma.GetClusterRepEdge(i)] +=nde;
		  }
	      }
	  for (i = 0; i < ned; i++)
	    cnt[nv+i]= first_edge_dof[i+1]-first_edge_dof[i];
	  for (i = 0; i < nfa; i++)
	    {
	      ma.GetFaceEdges (i, f2ed);
	      for (j = 0; j < f2ed.Size(); j++)
		cnt[nv+f2ed[j]] +=  first_face_dof[i+1]-first_face_dof[i];
	    }
	  for (i=0; i< ni; i++) 
	    {
	      ma.GetElEdges (i, enums, orient);
	      int ndi = first_inner_dof[i+1] - first_inner_dof[i];
	      for (j = 0; j < enums.Size(); j++)
		cnt[nv+enums[j]] += ndi;
	    }

	  break;
      
	}
     case 21: // wave equation
	{
	  int ds_order = precflags.GetNumFlag ("ds_order", 0);
	  if (ds_order < 0) ds_order = 0;	  

	  for(i=0;i<ned;i++)
	    {
	      cnt[i] =  first_edge_dof[i+1] - first_edge_dof[i] - ds_order;
	      if (cnt[i] < 0) cnt[i] = 0;
	    }
	  
	  for (i = 0; i < nfa; i++)
	    {
	      // cnt[ned+i] = first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i];

	      // int first = first_face_dof[i];
	      int p = order_face[i][0];
	      
	      int ii = 0;
	      for (int j = 0; j <= p-2; j++)
		for (int k = 0; k <= p-2-j; k++, ii++)
		  if (j+k+2 > ds_order)
		    cnt[ned+i]++;
	      //clusters[first+ii] = 1;
	      
	      // other combination
	      for (int j = 0; j <= p-2; j++)
		for (int k = 0; k <= p-2-j; k++, ii++)
		  if (j+k+2 > ds_order)
		    cnt[ned+i]++;
	      // clusters[first+ii] = 1;
	      
	      // type 3
	      for (int j = 0; j <= p-2; j++, ii++)
		if (j+2 > ds_order)
		  cnt[ned+i]++;
	      // clusters[first+ii] = 1;
	    }
	  break;
	}
      }
    
   
        
    Table<int> & table = *new Table<int> (cnt); 
    ii = 0; 
    cnt = 0;
    switch(SmoothingType) 
      {
      case 4: // AFW(loE) + E  + F + I (noLineClusters) 
	{
	  cnt = 0;
	  for (i = 0; i < ned; i++)
	    if(fine_edge[i])
	      {
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);	      
		table[pn1][cnt[pn1]++] = i;
		table[pn2][cnt[pn2]++] = i;
	      }
	  int nnv = nv; 
	  
	 
	  if(augmented==1)
	    { 
	      nnv += nv; 
	      for(i=0;i<nv;i++)
		table[nv+i][cnt[nv+i]++] = ned+i;
	    } 
	  
	  if(!excl_grads)
	    for (i =0; i<ned ; i++)
	      for (j = first_edge_dof[i]; j < first_edge_dof[i+1]; j++)
		table[nnv+i][cnt[nnv+i]++] = j;
	  
	  for (i = 0; i < nfa; i++)
	    { 
	      int first = first_face_dof[i] + excl_grads*face_ngrad[i]; 
	      for (j = first; j < first_face_dof[i+1]; j++)
		table[nnv+ned+i][cnt[nnv+ned+i]++] = j;
	    }
	      
	  for (i = 0; i < ni; i++)
	    { 
	      // int first = first_inner_dof[i] + excl_grads*cell_ngrad[i];
	      for (j = first_inner_dof[i]; j < first_inner_dof[i+1]; j++)
		table[nnv+ned+nfa+i][cnt[nnv+ned+nfa+i]++] = j;
	    }
	  break;
	}
      case 5: // AFW(hoE) + E  + F + I (noLineClusters) 
	{
	  cnt = 0;
	  
	  for (i = 0; i < ned; i++)
	    if(fine_edge[i])
	      {
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);
		table[pn1][cnt[pn1]++]  = i; 
		table[pn2][cnt[pn2]++]  = i; 
		if(!excl_grads)
		  for(j=first_edge_dof[i];j<first_edge_dof[i+1];j++) 
		    {
		      table[pn1][cnt[pn1]++]  = j; 
		      table[pn2][cnt[pn2]++]  = j;
		    }
	      }
	  int nnv = nv; 
	  /* if(augmented==1)
	     { 
	     nnv += nv; 
	     for(i=0;i<nv;i++)
	     table[i][cnt[i]++] = ned+i;
	     } */ 
	  for (i = 0; i < nfa; i++)
	    { 
	      int first = first_face_dof[i] + excl_grads*face_ngrad[i]; 
	      for (j = first; j < first_face_dof[i+1]; j++)
		table[nnv+ned+i][cnt[nnv+ned+i]++] = j;
	    }
	  
	  for (i = 0; i < ni; i++)
	    { 
	      // int first = first_inner_dof[i] + excl_grads*cell_ngrad[i];
	      for (j = first_inner_dof[i]; j < first_inner_dof[i+1]; j++)
		table[nnv+ned+nfa+i][cnt[nnv+ned+nfa+i]++] = j;
	    }
	  break;
	}
	
      case 1: // Clustering + AFW(loE) 
	{
	  cnt = 0; 

	  for (i = 0; i < ned; i++)
	    if(fine_edge[i])
	      {
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);
		pn1 = ma.GetClusterRepVertex(pn1);
		pn2 = ma.GetClusterRepVertex(pn2);
		table[pn1][cnt[pn1]++] = i;
		if (pn1 != pn2)
		  table[pn2][cnt[pn2]++] = i;

		int nr = ma.GetClusterRepEdge (i);
		if(!excl_grads)
		  for (j = first_edge_dof[i]; j < first_edge_dof[i+1]; j++)
		    table[nr][cnt[nr]++] = j;
	      }
	  for (i = 0; i < nfa; i++)
	    if(fine_face[i])
	      {
		int nr = ma.GetClusterRepFace (i);
		int first = first_face_dof[i] + excl_grads*face_ngrad[i];
		for (j = first; j < first_face_dof[i+1]; j++)
		  table[nr][cnt[nr]++] = j;
	      }
	  for (i = 0; i < ni; i++)
	    {
	      int nr = ma.GetClusterRepElement (i);
	      int first = first_inner_dof[i] + excl_grads*cell_ngrad[i];
	      for (j = first; j < first_inner_dof[i+1]; j++)
		table[nr][cnt[nr]++] = j;
	    }
	  break;
	}
	
      case 2: // Clustering with AFW(hoEdges) 
	{
	  for(i=0;i<ned;i++) 
	    if(fine_edge[i] && !IsDirichletEdge(i))
	      {
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);
		pn1 = ma.GetClusterRepVertex(pn1);
		pn2 = ma.GetClusterRepVertex(pn2);
		int ecl = ma.GetClusterRepEdge(i);
		table[pn1][cnt[pn1]++] = i;
		if(!excl_grads)
		  for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
		    table[pn1][cnt[pn1]++] = j;
		if(pn1!=pn2)
		  {
		    table[pn2][cnt[pn2]++] = i; 
		    table[ecl][cnt[ecl]++] = i; 
		    // ho-edges
		    if(!excl_grads)
		      for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
			{		    
			  table[pn2][cnt[pn2]++] = j;
			  table[ecl][cnt[ecl]++] = j;
			}
		  }
	      }
	  
	  // Stack of horizontal edges, quadfaces and prism faces 
	  //  (+ inner) 
	  for(i=0;i<nfa;i++) 
	    if(fine_face[i] && !IsDirichletFace(i)) 
	      { 
		int fcl = ma.GetClusterRepFace(i); 
		int first = first_face_dof[i] + excl_grads*face_ngrad[i];
		//	if(fcl  >= nv + ned) // quad-faces (prism) 
		for(j=first; j< first_face_dof[i+1]; j++) 
		  table[fcl][cnt[fcl]++] = j; 
	      }
	  
	  for(i=0;i<ni;i++) 
	    { 
	      // Stack: plane-face + inner
	      int ccl = ma.GetClusterRepElement(i); 
	      int first = first_inner_dof[i] + excl_grads*cell_ngrad[i];
	      for(j=first;j<first_inner_dof[i+1];j++) 
		table[ccl][cnt[ccl]++] = j; 
	    }
	  break;
	}
      case 3: 
	{
	  // Stack-AFW(hoEdges) 
	  for(i=0;i<ned;i++) 
	    if(fine_edge[i])
	      {
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);
		pn1 = ma.GetClusterRepVertex(pn1);
		pn2 = ma.GetClusterRepVertex(pn2);
		int ecl = ma.GetClusterRepEdge(i);
		table[pn1][cnt[pn1]++] = i;
		if(!excl_grads)
		  for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
		    table[pn1][cnt[pn1]++] = j;
		if(pn1!=pn2)
		  {
		    table[pn2][cnt[pn2]++] = i; 
		    table[ecl][cnt[ecl]++] = i; 
		    // ho-edges
		    if(!excl_grads)
		      for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
			{		    
			  table[pn2][cnt[pn2]++] = j;
			  table[ecl][cnt[ecl]++] = j;
			}
		  }
	      }
	  
	  // Stack of horizontal edges, quadfaces and prism faces 
	  //  (+ inner) 
	  for(i=0;i<nfa;i++) 
	    if(fine_face[i]) 
	      { 
		int fcl = ma.GetClusterRepFace(i); 
		int first = first_face_dof[i] + excl_grads*face_ngrad[i];
		
		for(j=first; j< first_face_dof[i+1]; j++) 
		  table[fcl][cnt[fcl]++] = j; 
	      }
	  
	  for(i=0;i<ni;i++) 
	    { 
	      // Stack: plane-face + inner
	      int ccl = ma.GetClusterRepElement(i);
	      int first = first_inner_dof[i] + excl_grads*cell_ngrad[i]; 
	      for(j=first;j<first_inner_dof[i+1];j++) 
		table[ccl][cnt[ccl]++] = j; 
	         	      
	      // inner to quad-face (horizontal edge stack) 
	      // each face different stack 
	      
	      ma.GetElFaces(i,fanums); 
	      
	      for(j=0;j<fanums.Size();j++) 
		{ 
		  int fcl = ma.GetClusterRepFace(fanums[j]);
		  if(fcl != ccl)
		    for(k=first;k<first_inner_dof[i+1];k++) 
		      table[fcl][cnt[fcl]++] = k; 
		} 
	      
	    }  
	  
	  for(i=0;i<nfa;i++) 
	    if(fine_face[i]) 
	      { 
		ma.GetFacePNums(i,vnums); 
		if(vnums.Size()==4) continue; 
		int fcl = ma.GetClusterRepFace(i); 
		ma.GetFaceEdges(i,ednums); 
		
		for(j=0;j<ednums.Size();j++)
		  {
		    int ecl = ma.GetClusterRepEdge(ednums[j]); 
		    if(ecl==fcl) continue; 
		    int first = first_face_dof[i] + excl_grads*face_ngrad[i];
		    for(k=first; k<first_face_dof[i+1]; k++)
		      table[ecl][cnt[ecl]++] = k; 
		  }
	      }
	  
	  break;
	}
      case 6: // Jacobi(Diag)
	{
	  cnt = 0;
	  int ii=0; 
	  
	  for (i = 0; i < ned; i++)
	    if(fine_edge[i])
	      table[ii++][0] = i;

	  for (i =0; i<ned ; i++)
	    if(fine_edge[i])
	      for (j = first_edge_dof[i]; j < first_edge_dof[i+1]; j++)
		table[ii++][0] = j;
	  
	  for (i = 0; i < nfa; i++)
	    { 
	      if(fine_face[i]) 
		{
		  int first = first_face_dof[i]; 
		  for (j = first; j < first_face_dof[i+1]; j++)
		    table[ii++][0] = j;
		}
	    }
	      
	  for (i = 0; i < ni; i++)
	    { 
	      for (j = first_inner_dof[i]; j < first_inner_dof[i+1]; j++)
		table[ii++][0] = j;
	    }
	  break;
	}
	
      case 7: 
	{
	  // Stack-AFW(hoEdges) 
	  for(i=0;i<ned;i++) 
	    if(fine_edge[i])
	      {
		int pn1, pn2;
		ma.GetEdgePNums (i,pn1,pn2);
		pn1 = ma.GetClusterRepVertex(pn1);
		pn2 = ma.GetClusterRepVertex(pn2);
		//	int ecl = ma.GetClusterRepEdge(i);
		table[pn1][cnt[pn1]++] = i;
		if(!excl_grads)
		  for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
		    table[pn1][cnt[pn1]++] = j;
		if(pn1!=pn2)
		  {
		    table[pn2][cnt[pn2]++] = i; 
		    //  table[ecl][cnt[ecl]++] = i; 
		    // ho-edges
		    if(!excl_grads)
		      for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
			{		    
			  table[pn2][cnt[pn2]++] = j;
			  //  table[ecl][cnt[ecl]++] = j;
			}
		  }
	      }
	  for (i = 0; i < ned; i++)
	    {
	      first = first_edge_dof[i];
	      int ndof = first_edge_dof[i+1]-first;
	      for (j = 0; j < ndof; j++)
		table[nv+i][cnt[nv+i]++] = first+j;
	    }
	  for (i = 0; i < nfa; i++)
	    {
	      first = first_face_dof[i];
	      int ndof = first_face_dof[i+1]-first;
	      ma.GetFaceEdges (i, f2ed);
	      for (k = 0; k < f2ed.Size(); k++)
		for (j = 0; j < ndof; j++)
		  table[nv+f2ed[k]][cnt[nv+f2ed[k]]++] = first+j;
	    }
	 
	  for (i = 0; i < ni; i++)
	    {
	      ma.GetElEdges (i, enums, orient);
	      first = first_inner_dof[i];
	      int ndof = first_inner_dof[i+1] - first_inner_dof[i];
	      for (j = 0; j < enums.Size(); j++)
		for (k = 0; k < ndof; k++)
		  table[nv+enums[j]][cnt[nv+enums[j]]++] = first + k;
	    } 
	  
	  
	  break;
	}
      case 21: // wave equation
	{
	  cnt = 0;
	  int ds_order = precflags.GetNumFlag ("ds_order", 0);
	  if (ds_order < 0) ds_order = 0;	  
	  
	  for (i =0; i<ned ; i++)
	    for (j = first_edge_dof[i]+ds_order; j < first_edge_dof[i+1]; j++)
	      table[i][cnt[i]++] = j;
	  /*
	  for (i = 0; i < nfa; i++)
	    { 
	      int first = first_face_dof[i]; //  + excl_grads*face_ngrad[i]; 
	      for (j = first; j < first_face_dof[i+1]; j++)
		table[ned+i][cnt[ned+i]++] = j;
	    }
	  */

	  for (i = 0; i < nfa; i++)
	    {
	      // cnt[ned+i] = first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i];

	      int first = first_face_dof[i];
	      int p = order_face[i][0];
	      
	      int ii = first;
	      for (int j = 0; j <= p-2; j++)
		for (int k = 0; k <= p-2-j; k++, ii++)
		  if (j+k+2 > ds_order)
		    table[ned+i][cnt[ned+i]++] = ii;
	      // cnt[ned+i]++;
	      //clusters[first+ii] = 1;
	      
	      // other combination
	      for (int j = 0; j <= p-2; j++)
		for (int k = 0; k <= p-2-j; k++, ii++)
		  if (j+k+2 > ds_order)
		    table[ned+i][cnt[ned+i]++] = ii;
	      // cnt[ned+i]++;
	      // clusters[first+ii] = 1;
	      
	      // type 3
	      for (int j = 0; j <= p-2; j++, ii++)
		if (j+2 > ds_order)
		  table[ned+i][cnt[ned+i]++] = ii;
		    // cnt[ned+i]++;
	      // clusters[first+ii] = 1;
	    }
	  break;
	}      
      }
    //(*testout) << "H(Curl)-table = " << table << endl;	
    
    return & table; 
  }

    

  Array<int> *   HCurlHighOrderFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    cout << "called createdirectsolverclusters" << endl;
    // 
    if (precflags.NumFlagDefined ("ds_order"))
      {
	int ds_order = int (precflags.GetNumFlag ("ds_order", 1));

	Array<int> & clusters = *new Array<int> (GetNDof());
	clusters = 0;
	
	int ned = ma.GetNEdges();
	for (int i = 0; i < ned; i++)
	  clusters[i] = 1;

	for (int i = 0; i < ned; i++)
	  {
	    int first = first_edge_dof[i];
	    int next = first_edge_dof[i+1];
	    for (int j = 0; (j < ds_order) && (first+j < next) ; j++)
	      clusters[first+j] = 1;
	  }

	int nfa = ma.GetNFaces();
	for (int i = 0; i < nfa; i++)
	  {
	    int first = first_face_dof[i];
	    // int next = first_face_dof[i+1];
	    int p = order_face[i][0];
	    
	    // if (usegrad_face[i])
	    int ii = 0;
            for (int j = 0; j <= p-2; j++)
              for (int k = 0; k <= p-2-j; k++, ii++)
                if (j+k+2 <= ds_order)
		  clusters[first+ii] = 1;
	    
	    // other combination
	    for (int j = 0; j <= p-2; j++)
	      for (int k = 0; k <= p-2-j; k++, ii++)
		if (j+k+2 <= ds_order)
		  clusters[first+ii] = 1;
	    
	    // type 3
	    for (int j = 0; j <= p-2; j++, ii++)
	      if (j+2 <= ds_order)
	      clusters[first+ii] = 1;
	  }

	return &clusters;
      }


    int clustertype = int(precflags.GetNumFlag("ds_cluster",4));  

    cout << " DirectSolverCluster Clustertype " << clustertype << endl; 
    if(clustertype==0)
      return(0);
	   
    // int nv = ma.GetNV();
    // int nd = GetNDof();
    int ne = ma.GetNE();
    int ned = ma.GetNEdges();

    Array<int> ednums, fnums, pnums;





    if (precflags.GetDefineFlag("subassembled"))
      {

	Array<int> & clusters = *new Array<int> (GetNDof());
	clusters = 0;


	for (int i = 0; i < ned; i++)
	  if (!IsDirichletEdge(i) && fine_edge[i])
	    clusters[i] = 1;

	return &clusters;
      }




    int i, j, k;
    bool hasprism = false;

    for (i = 0; !hasprism && i < ne; i++)
      if (ma.GetElType(i) == ET_PRISM)
	hasprism = true;
    
    if (!hasprism && adddirectsolverdofs.Size() == 0 &&
	directedgeclusters.Size() == 0 && directfaceclusters.Size() == 0 && directelementclusters.Size() == 0) 
      return NULL;
        
    Array<int> & clusters = *new Array<int> (GetNDof());
    clusters = 0;


    if(directedgeclusters.Size() != 0 || 
       directfaceclusters.Size() != 0 ||
       directelementclusters.Size() != 0)
      {
	for(i=0; i<directedgeclusters.Size(); i++)
	  for(j=first_edge_dof[directedgeclusters[i]]; j<first_edge_dof[directedgeclusters[i]+1]; j++)
	    clusters[j] = 6;

	for(i=0; i<directfaceclusters.Size(); i++)
	  {
	    for(j=first_face_dof[directfaceclusters[i]];  j< first_face_dof[directfaceclusters[i]] + face_ngrad[directfaceclusters[i]]; j++)
	      clusters[j] = 6;
	  }

	for(i=0; i<directelementclusters.Size(); i++)
	  {
	    for(j=first_inner_dof[directelementclusters[i]];  j< first_inner_dof[directelementclusters[i]] + cell_ngrad[directelementclusters[i]]; j++)
	      clusters[j] = 6;
	  }	    
      }

    if(hasprism)
      {
	// for long prisms, JS June 2005
	switch(clustertype)
	  {

	  case 0: 
	    //clusters = 0; 
	    break; 
	  case 1:
	    //clusters = 0;
	    for (i = 0; i < ne; i++)
	      {
		if (ma.GetElType(i) == ET_PRISM)
		  {
		    ma.GetElEdges (i, ednums);
		
		    for (j = 6; j < 9; j++)  //vertical Edges 
		      { 
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (k = first; k < next; k++)
			  clusters[k] = 1;      //verthoedge
			clusters[ednums[j]] = 1;
		      }

		    ma.GetElFaces (i, fnums); 

		    for (j=2; j<5 ; j++)  // vertical faces
		      {
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			for(k=first; k<next; k++)
			  clusters[k] = 1;
		    
			// INT<2> p = order_face[fnums[j]]; 
			// for(k=next-p[0]-p[1]; k<next; k++)
			// clusters[k] = 1;
		      }

		    int first = first_inner_dof[i]; 
		    int next = first_inner_dof[i+1]; 

		    for(k=first; k<next; k++)
		      clusters[k] = 1;
		  }
	      }
	    break; 
	  case 2: 
	    //clusters = 0;
	
	    // All Vertical Edges in one Cluster for Hex and Prism + const_x*poly faces + const_y*polx faces (-> 2d Problems !) 
	
	
	    //lo 
	    for(i=0;i<ma.GetNEdges();i++)
	      clusters[i]=1; 

       
	    for (i = 0; i < ne; i++)
	      {
		/*if (ma.GetElType(i) == ET_PYRAMID)
		  {
		  GetDofNrs(i,ednums);
		  for(j=0; j<ednums.Size(); j++) clusters[ednums[j]] = 3;
	    
		  }   
		*/
		if (ma.GetElType(i) == ET_PRISM)
		  {
		    ma.GetElEdges (i, ednums);
		    for (j = 0; j < 6; j++)  //horizontal Edges 
		      { 
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (k = first; k < next; k++)
			  clusters[k] = 1; 
			clusters[ednums[j]]=1;
		      }
		
		
		    for (j = 6; j < 9; j++)  //vertical Edges 
		      { 
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (k = first; k < next; k++)
			  clusters[k] = 1;      //verthoedge
			clusters[ednums[j]]=1; //ned
		      }
		    ma.GetElFaces (i, fnums); 
		
		    for (j=2; j<5 ; j++)  // vertical faces
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			// TEST SZ eigentlich nurin eine richtung konst ausreichend
			for (k=first; k < next; k++) 
			  clusters[k]=1; 
		    
		    
			INT<2> p = order_face[fnums[j]]; 
		    
			for(k=next-p[0]-p[1]; k<next; k++)
			  clusters[k] = 1;
			for(k=0; k<4; k++)
			  clusters[k] = 1;
		    
		    
			// for (k=first; k < next; k++) 
			//   clusters[k]=3; 
		    
		      }
		    for (j=0; j<2 ; j++)  // horizontal faces
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			for (k=first; k < next; k++) 
			  clusters[k]=1; 
		      }
		  }
	    
		else if (ma.GetElType(i) == ET_HEX)
		  {
		    ma.GetElEdges (i, ednums);
		
		    for(j=0;j<8;j++) // horizontal edges
		      {	
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (k = first; k < next; k++)
			  clusters[k] = 0;
		
			clusters[ednums[j]]=0; 
		    
		      }
			
		    for (j = 8; j < 12; j++)  // vertical edges 
		      {	
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (k = first; k < next; k++)
			  clusters[k] = 3;
			clusters[ednums[j]]=0; 
		      }
		
		    ma.GetElFaces(i,fnums); // vertical faces 
		    for(j=2;j<6;j++) 
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			/*	for (k=first; k < next; k++) 
			  clusters[k]=0; 
		      
			  INT<2> p = order_face[fnums[j]];
			  for(k=2*(p[0]+1)*(p[0]+1);k<next;k++) //achtung
			  clusters[k]=3;  
			*/ 
			// TEST SZ  eigentlich in eine richtung konstante nur benoetigt
			for (k=first; k < next; k++) 
			  clusters[k]=3; 
		    
		      }
		    for(j=0;j<2;j++)  //horizontal faces 
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			for (k=first; k < next; k++) 
			  clusters[k]=0; 
		      } 
		  }
	    
		for(k=first_inner_dof[i];k<first_inner_dof[i+1];k++) 
		  clusters[k]=0; 
	    
	      }
    
	    break; 
	  case 3: 
	    //lo 
	    //for(i=0;i<ma.GetNEdges();i++)
	    //  clusters[i]=0; 
	
	    for (i = 0; i < ne; i++)
	      {
		ma.GetElPNums(i,pnums); 
		if (ma.GetElType(i) == ET_PRISM)
		  {
		    ma.GetElEdges (i, ednums);
		    for (j = 0; j < 6; j++)  //horizontal Edges 
		      { 
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (k = first; k < next; k++)
			  clusters[k] = 0;
			clusters[ednums[j]]=0; 
		      }
		    for (j = 6; j < 9; j++)  //vertical Edges 
		      { 
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (k = first; k < next; k++)
			  clusters[k] = 3;      //verthoedge
			clusters[ednums[j]]=0; //ned
		      }
		    ma.GetElFaces (i, fnums); 
		
		    for (j=2; j<5 ; j++)  // vertical faces
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			/* TEST SZ eigentlich nurin eine richtung konst ausreichend
			   for (k=first; k < next; k++) 
			   clusters[k]=0; 
		       
			   int p = order_face[fnums[j]][0]; // achtung !!! 
			   for(k=2*(p+1)*(p+1);k<next;k++)
			   clusters[k]=3;  
		       
			*/
		    
			for (k=first; k < next; k++) 
			  clusters[k]=3; 
		    
		      }
		    for (j=0; j<2 ; j++)  // horizontal faces
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			for (k=first; k < next; k++) 
			  clusters[k]=0; 
		      }
		  }
	    
		else if (ma.GetElType(i) == ET_HEX)
		  {
		    ma.GetElEdges (i, ednums);
		    for(j=0;j<8;j++) // horizontal edges
		      {	
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (k = first; k < next; k++)
			  clusters[k] = 0;
		    
			clusters[ednums[j]]=0; 
		    
		      }
		    for (j = 8; j < 12; j++)  // vertical edges 
		      {	
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (k = first; k < next; k++)
			  clusters[k] = 3;
			clusters[ednums[j]]=0; 
		      }
		
		    ma.GetElFaces(i,fnums); // vertical faces 
		    for(j=2;j<6;j++) 
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			/*	for (k=first; k < next; k++) 
			  clusters[k]=0; 
		      
			  INT<2> p = order_face[fnums[j]];
			  for(k=2*(p[0]+1)*(p[0]+1);k<next;k++) //achtung
			  clusters[k]=3;  
			*/ 
			// TEST SZ  eigentlich in eine richtung konstante nur benoetigt
			for (k=first; k < next; k++) 
			  clusters[k]=3; 
		    
		      }
		    for(j=0;j<2;j++)  //horizontal faces 
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			for (k=first; k < next; k++) 
			  clusters[k]=0; 
		      } 
		  }
	    
		for(k=first_inner_dof[i];k<first_inner_dof[i+1];k++) 
		  clusters[k]=0; 
	    
	      }
    
	    //  (*testout) << "direct clusters = " << endl << clusters << endl;
	
	    for(i=0; directsolverclustered.Size() > 0 && i<ne; i++)
	      {
		if(directsolverclustered[ma.GetElIndex(i)])
		  {
		    GetDofNrs(i,ednums);
		    for(k=0; k<ednums.Size(); k++)
		      {
			clusters[ednums[k]] = 4;
		      }
		  }
	      }
	    break; 

	  case 4:  // just like the old hcurl 
	    //lo 
	    //clusters = 0; 
	    for(i=0;i<ned;i++) 
	      { 
		int pi1,pi2; 
		ma.GetEdgePNums(i,pi1,pi2); 
		pi1 = ma.GetClusterRepVertex (pi1); 
		pi2 = ma.GetClusterRepVertex (pi2);
		if(pi1 == pi2) 
		  // stack of horizontal edge-blocks (without low-order)
		  // -> decoupled 
		  for(j=first_edge_dof[i];j <first_edge_dof[i+1];j++) 
		    clusters[j] = 1; 
	      }

	    for(i=0;i<nfa;i++)
	      {
		// Attenzione das is zuviel !!! 

		ma.GetFacePNums(i,pnums); 
		if(pnums.Size() == 4) // quad face 
		  { 
		    INT<2> p = order_face[i];         
		    int first =first_face_dof[i];  
		    int next = first_face_dof[i+1]; 
		    // Ned_0*pol_z 
		    // int first = next - p[0] - p[1]; 
	
		    for(j= first ; j<next; j++)
		      clusters[j] = 1; 
		  } 
	      }
	    break; 
	  case 5:  // just like the old hcurl horizontal only constant ... 
	    //lo 
	    //clusters = 0; 
	    for(i=0;i<ned;i++) 
	      { 
	    
		int pi1,pi2; 
		ma.GetEdgePNums(i,pi1,pi2); 
		pi1 = ma.GetClusterRepVertex (pi1); 
		pi2 = ma.GetClusterRepVertex (pi2);
		if(pi1 == pi2) 
		  // stack of horizontal edge-blocks (without low-order)
		  // -> decoupled 
		  {
		    //clusters[i] = 1; 
		    for(j=first_edge_dof[i];j <first_edge_dof[i+1];j++) 
		      clusters[j] = 1; 
		  }
	      }

	    for(i=0;i<nfa;i++)
	      {
		// Attenzione das is zuviel !!! 

		ma.GetFacePNums(i,pnums); 
		if(pnums.Size() == 4) // quad face 
		  { 
		    INT<2> p = order_face[i];        
		    // int first =first_face_dof[i];  
		    int next = first_face_dof[i+1]; 
		    // Ned_0*pol_z 
		    int first = next - p[0] - p[1]; 
	
		    for(j= first ; j<next; j++)
		      clusters[j] = 1; 
		  } 
	  
	    
	    
	      }


	    //  (*testout) << "direct clusters = " << endl << clusters << endl;
	
	    for(i=0; directsolverclustered.Size() > 0 && i<ne; i++)
	      {
		if(directsolverclustered[ma.GetElIndex(i)])
		  {
		    ELEMENT_TYPE eltype = ma.GetElType(i); 
		    if(eltype != ET_PRISM) continue; 
		
		    GetDofNrs(i,ednums);
		    for(k=0; k<ednums.Size(); k++)
		      if(ednums[k]>=0) clusters[ednums[k]] = 2;
		  }
	      }
	    break; 
      
      
	  }

      }

    //if(adddirectsolverdofs.Size())
    //  (*testout) << "addirectsolverdofs " << adddirectsolverdofs << endl;

    //int numadd = 0;

    for(i=0; i< adddirectsolverdofs.Size(); i++)
      {
	//if(clusters[adddirectsolverdofs[i]] != 5)
	//  numadd++;
	clusters[adddirectsolverdofs[i]] = 5;
      }

    //(*testout) << "number added dofs " << numadd << endl;

    

    //(*testout) << "clusters = " << clusters << endl;
    return &clusters;
    
  }
  


    
  
  /*  BitArray * HCurlHighOrderFESpace :: 
      CreateIntermediatePlanes (int type) const
      {
      int i;
      Array<int> vnums;
 
      //bool has_cluster = 0;
      for (i = 0; i < ned; i++)
      {
      int pi1, pi2;
      ma.GetEdgePNums (i, pi1, pi2);
	
      if (ma.GetClusterRepVertex (pi1) ==
      ma.GetClusterRepVertex (pi2))
      has_cluster = 1;
      }

      if (!has_cluster) return 0;
    
   
      BitArray & ba = *new BitArray (GetNDof());
      ba.Clear();

      for (i = 0; i < ned; i++)
      {
      int pi1, pi2;
      ma.GetEdgePNums (i, pi1, pi2);
	
      if (ma.GetClusterRepVertex (pi1) ==
      ma.GetClusterRepVertex (pi2))
      {
      ba.Set (i);
      for (int l = first_edge_dof[i];
      l < first_edge_dof[i+1]; l++)
      ba.Set (l);
      }
      }
    
    
      for (i = 0; i < nfa; i++)
      {
      ma.GetFacePNums (i, vnums);
      if (vnums.Size() == 4)
      {
      for (int l = first_face_dof[i];
      l < first_face_dof[i+1]; l++)
      ba.Set (l);
      }	    
      }
      return &ba;
      }*/


  SparseMatrix<double> * 
  HCurlHighOrderFESpace :: CreateGradient() const
  {
    Flags flags2 (flags);

    if(iscomplex)
      flags2.SetFlag("complex");
    flags2.SetFlag("order", order+1);
    flags2.SetFlag("relorder", rel_order+1); 
    flags2.SetFlag("orderinner",uniform_order_inner+1); 
    flags2.SetFlag("orderface",uniform_order_face+1); 
    flags2.SetFlag("orderedge",uniform_order_edge+1); 
    flags2.SetFlag ("orderquad", -1);
    flags2.SetFlag ("ordertrig", -1);  
    flags2.SetFlag("variableorder",var_order); 
    
    flags2.SetFlag("print");
   
    // if the following three flags are set -> relorder is used 
    flags2.SetFlag("relorder",rel_order+1);
    flags2.SetFlag("order",order+1); 
    flags2.SetFlag("variableorder"); 
    
    /* 
       if(var_order)
       flags2.SetFlag("relorder",rel_order+1); 
       else
       flags2.SetFlag("order",order+1); 
    */ 

    if(uniform_order_inner>-1)
      flags2.SetFlag("orderinner",uniform_order_inner+1); 
    if(uniform_order_face>-1)
      flags2.SetFlag("orderface",uniform_order_face+1); 
    if(uniform_order_edge>-1)
      flags2.SetFlag("orderedge",uniform_order_edge+1); 
    
    // flags2.SetFlag ("orderquad", -1);
    // flags2.SetFlag ("ordertrig", -1);  
    // flags2.SetFlag("variableorder",var_order); 
    
    H1HighOrderFESpace  fesh1(ma, flags2); 

    BitArray h1def(ma.GetNDomains());
    if(definedon.Size() == 0)
      h1def.Set();
    else
      h1def.Clear();
    for(int i=0; i<definedon.Size(); i++)
      if(definedon[i] > 0)
	h1def.Set(i);
    fesh1.SetDefinedOn(h1def);

    BitArray h1defb(ma.GetNBoundaries());
    if(definedonbound.Size() == 0)
      h1defb.Set();
    else
      h1defb.Clear();
    for(int i=0; i<definedonbound.Size(); i++)
      if(definedonbound[i] > 0)
	h1defb.Set(i);
    fesh1.SetDefinedOnBoundary(h1defb);

    LocalHeap lh(100008);
    fesh1.Update(lh);
     
    int ned = ma.GetNEdges(); 
    // int nv  = ma.GetNV(); 
    int nel = ma.GetNE(); 
    int nfa = 0;
    if(ma.GetDimension()==3) nfa= ma.GetNFaces();  

    // int dim     = fesh1.GetDimension();
    // int dimcurl = GetDimension();

    Array<int> dnums_h1l; 
    Array<int> dnums_hcl;
    
    // Matrix Graph for gradient matrix , start (end) position for assembling 
    Array<int> cnts(ndof);
    cnts = 0; 

    // *testout << "grad fine edges " << fine_edge << endl ; 
    // *testout << "grad fine faces " << fine_edge << endl ; 
   
    for(int i=0; i<ned; i++)
      {
	if(fine_edge[i])
	  {
	    cnts[i] = 2;  // vertices-nedelec
	    int l = first_edge_dof[i]; 
	    for( int k = fesh1.GetFirstEdgeDof(i); k< fesh1.GetFirstEdgeDof(i+1); 
		 k++, l++)
	      cnts[l] = 1;
	  }
      }
    
    for(int i=0; i< nfa; i++) 
      {
	if(fine_face[i])
	  {
	    int l= first_face_dof[i]; 
	    for (int k = fesh1.GetFirstFaceDof(i); k<fesh1.GetFirstFaceDof(i+1); 
		 k++, l++) 
	      cnts[l] = 1;
	  }
      }
    
    for(int i=0; i<nel; i++)
      {
	int l= first_inner_dof[i]; 
	for (int k = fesh1.GetFirstElementDof(i); k<fesh1.GetFirstElementDof(i+1); 
	     k++, l++) 
	  cnts[l] = 1; 
      }
  
    //sparse matrix with above matrix graph cnts 
    SparseMatrix<double> & grad = *new SparseMatrix<double>(cnts); 
    
    // vertices - nedelec
    // ho-edges
    for (int i = 0; i < ned; i++) 
      {
	if(fine_edge[i])
	  { 
	    int p1,p2; 
	    ma.GetEdgePNums(i,p1,p2); 

	    grad.CreatePosition(i,p1); 
	    grad.CreatePosition(i,p2); 

	    if (p1 < p2) 
	      {	      
		grad(i,p1) = -1.;
		grad(i,p2) = 1.; 
	      }
	    else
	      {
		grad(i,p1) = 1.; 
		grad(i,p2) = -1.; 
	      }
	   
	    int l = first_edge_dof[i]; 

	    for( int k = fesh1.GetFirstEdgeDof(i); k< fesh1.GetFirstEdgeDof(i+1); 
		 k++, l++)
	      {
		grad.CreatePosition(l,k); 
		grad(l,k)=1.; 
	      }
	  }
      }
    
    for(int i=0; i< nfa; i++) 
      {
	if(fine_face[i])
	  {
	    int l= first_face_dof[i]; 
	    for (int k = fesh1.GetFirstFaceDof(i); k<fesh1.GetFirstFaceDof(i+1); 
		 k++, l++) 
	      {
		grad.CreatePosition(l,k); 
		grad(l,k)=1.;
	      }
	  }
      }

    for(int i=0; i<nel; i++)
      {
	int l= first_inner_dof[i]; 
	for (int k = fesh1.GetFirstElementDof(i); k<fesh1.GetFirstElementDof(i+1); 
	     k++, l++) 
	  {
	    grad.CreatePosition(l,k); 
	    grad(l,k)=1.;
	  }	
      }
    //(*testout) << " Global grad " << grad << endl; 

    return &grad;
  }










#ifdef PARALLEL


#ifdef OLD_PARALLEL_UPDATE
  void HCurlHighOrderFESpace :: UpdateParallelDofs_hoproc()
  {


    if ( discontinuous )
      {
	*testout << "HCurlHOFESpace::UpdateParallelDofs_hoproc -- discontinuous" << endl;
	// Find number of exchange dofs
	Array<int> nexdof(ntasks);
	nexdof = 0;
	
	const MeshAccess & ma = (*this). GetMeshAccess();
	
	paralleldofs->SetNExDof(nexdof);
	
// 	paralleldofs->localexchangedof = new Table<int> (nexdof);
// 	paralleldofs->distantexchangedof = new Table<int> (nexdof);
	
	paralleldofs->sorted_exchangedof = new Table<int> (nexdof);

	for ( int i = 0; i < ndof; i++ )
	  {
	    paralleldofs->ClearExchangeDof(i);
	    for ( int dest = 0; dest < ntasks; dest++ )
	      paralleldofs->ClearExchangeDof(dest, i);
	  }
	
	
	int * ndof_dist = new int[ntasks];
	MPI_Allgather ( &ndof, 1, MPI_INT, &ndof_dist[1], 
			1, MPI_INT, MPI_HIGHORDER_COMM);
	for ( int dest = 0; dest < ntasks; dest++ )
	  paralleldofs -> SetDistNDof( dest, ndof_dist[dest]) ;
	
	delete [] ndof_dist;
	
	return;
      }

    // ******************************
    // update exchange dofs 
    // ******************************

    *testout << "HCurlHOFESpace::UpdateParallelDofs_hoproc" << endl;
    // Find number of exchange dofs
    Array<int> nexdof(ntasks);
    nexdof = 0;

    MPI_Status status;
    MPI_Request * sendrequest = new MPI_Request [ntasks];
    MPI_Request * recvrequest = new MPI_Request [ntasks];

    // number of edge exchange dofs
    for ( int edge = 0; edge < ma.GetNEdges(); edge++ )
      {
	if ( !parallelma->IsExchangeEdge ( edge ) ) continue;
	
	Array<int> dnums;
	GetEdgeDofNrs ( edge, dnums );
	nexdof[id] += dnums.Size() ;  
	
	for ( int dest = 1; dest < ntasks; dest ++ )
	  if (  parallelma -> GetDistantEdgeNum ( dest, edge ) >= 0 )
	    nexdof[dest] += dnums.Size() ;  
      }

    // + number of face exchange dofs
    for ( int face = 0; face < ma.GetNFaces(); face++ )
      {
	if ( !parallelma->IsExchangeFace ( face ) ) continue;
	
	Array<int> dnums;
	GetFaceDofNrs ( face, dnums );
	nexdof[id] += dnums.Size() ; 

	for ( int dest = 1; dest < ntasks; dest ++ )
	  if (  parallelma -> GetDistantFaceNum ( dest, face ) >= 0 )
	    nexdof[dest] += dnums.Size() ; 
      }

    // + number of inner exchange dofs
    for ( int el = 0; el < ma.GetNE(); el++ )
      {
	if ( !parallelma->IsExchangeElement ( el ) ) continue;
	
	Array<int> dnums;
	GetInnerDofNrs ( el, dnums );
	nexdof[id] += dnums.Size() ; 

	for ( int dest = 1; dest < ntasks; dest ++ )
	  if (  parallelma -> GetDistantElNum ( dest, el ) >= 0 )
	    nexdof[dest] += dnums.Size() ; 
      }

    nexdof[0] = LowOrderFESpace() . GetNDof();

    paralleldofs->SetNExDof(nexdof);

//     paralleldofs->localexchangedof = new Table<int> (nexdof);
//     paralleldofs->distantexchangedof = new Table<int> (nexdof);
    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);

    Array<int> ** owndofs, ** distantdofs;
    owndofs = new Array<int>* [ntasks];
    distantdofs = new Array<int>* [ntasks];

    for ( int i = 0; i < ntasks; i++ )
      {
	owndofs[i] = new Array<int>(1);
	(*owndofs[i])[0] = ndof;
	distantdofs[i] = new Array<int>(0);
      }

    Array<int> cnt_nexdof(ntasks);
    cnt_nexdof = 0;
    int exdof = 0;
    int ii;


    // *****************
    // Parallel Edge dofs
    // *****************


    for ( int edge = 0; edge < ma.GetNEdges(); edge++ )
      if ( parallelma->IsExchangeEdge ( edge ) )
	{
	  Array<int> dnums;
	  GetEdgeDofNrs ( edge, dnums );
	  if ( dnums.Size() == 0 ) continue;

	  for ( int i=0; i<dnums.Size(); i++ )
	    (*(paralleldofs->sorted_exchangedof))[id][exdof++] = dnums[i];

	  for ( int dest = 1; dest < ntasks; dest++ )
	    {
	      int distedge = parallelma -> GetDistantEdgeNum ( dest, edge );
	      if( distedge < 0 ) continue;

	      owndofs[dest]->Append ( distedge );
	      owndofs[dest]->Append (int(paralleldofs->IsGhostDof(dnums[0])) );
	      for ( int i=0; i<dnums.Size(); i++)
		{
		  paralleldofs->SetExchangeDof ( dest, dnums[i] );
		  paralleldofs->SetExchangeDof ( dnums[i] );
		  owndofs[dest]->Append ( dnums[i] );
		}
	    }
	}   


    for ( int sender = 1; sender < ntasks; sender ++ )
      {
        if ( id == sender )
          for ( int dest = 1; dest < ntasks; dest ++ )
            if ( dest != id )
              {
		MyMPI_ISend ( *owndofs[dest], dest, sendrequest[dest]);
              }
	  
        if ( id != sender )
          {
	    MyMPI_IRecv ( *distantdofs[sender], sender, recvrequest[sender]);
          }
 	  
	  
      }

    for( int dest=1; dest<ntasks; dest++) 
      {
	if ( dest == id ) continue;
	MPI_Wait ( recvrequest + dest, &status );
	paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
	// low order nedelec dofs first
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int ednum = (*distantdofs[dest])[ii++];
	    int isdistghost = (*distantdofs[dest])[ii++];
	    Array<int> dnums;
	    GetEdgeDofNrs (ednum, dnums);
// 	    (*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = dnums[0];
// 	    (*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];
	    (*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = dnums[0];
	    if ( dest < id && !isdistghost )
	      paralleldofs->ismasterdof.Clear ( dnums[0] ) ;
	    ii += dnums.Size();
	    cnt_nexdof[dest]++;
	  }
	// then the high order dofs, without nedelecs
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int ednum = (*distantdofs[dest])[ii++];
	    int isdistghost = (*distantdofs[dest])[ii++];
	    Array<int> dnums;
	    GetEdgeDofNrs (ednum, dnums);
	    ii++; 
	    for ( int i=1; i<dnums.Size(); i++)
	      {
// 		(*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
// 		(*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];
		(*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
		if ( dest < id && !isdistghost )
		  paralleldofs->ismasterdof.Clear ( dnums[i] ) ;
		ii++; cnt_nexdof[dest]++;
	      }
	  }
      }


    // *****************
    // Parallel Face dofs
    // *****************

    for ( int dest = 0; dest < ntasks; dest++)
      {
	owndofs[dest]->SetSize(1);
	distantdofs[dest]->SetSize(0);
      }
    for ( int face = 0; face < ma.GetNFaces(); face++ )
      if ( parallelma->IsExchangeFace ( face ) )
	{
	  Array<int> dnums;
	  GetFaceDofNrs ( face, dnums );
	  if ( dnums.Size() == 0 ) continue;

	  for ( int i=0; i<dnums.Size(); i++ )
	    (*(paralleldofs->sorted_exchangedof))[id][exdof++] = dnums[i];

	  for ( int dest = 1; dest < ntasks; dest++ )
	    {
	      int distface = parallelma -> GetDistantFaceNum ( dest, face );
	      if( distface < 0 ) continue;

	      owndofs[dest]->Append ( distface );
	      owndofs[dest]->Append (int(paralleldofs->IsGhostDof(dnums[0])) );

	      for ( int i=0; i<dnums.Size(); i++)
		{
		  paralleldofs->SetExchangeDof ( dest, dnums[i] );
		  paralleldofs->SetExchangeDof ( dnums[i] );
		  owndofs[dest]->Append ( dnums[i] );
		}
	    }
	}   


    for ( int sender = 1; sender < ntasks; sender ++ )
      {
        if ( id == sender )
          for ( int dest = 1; dest < ntasks; dest ++ )
            if ( dest != id )
              {
		MyMPI_ISend ( *owndofs[dest], dest, sendrequest[dest]);
              }
	  
        if ( id != sender )
          {
	    MyMPI_IRecv ( *distantdofs[sender], sender, recvrequest[sender] );
          }
 	  
	  
      }

    for( int dest=1; dest<ntasks; dest++) 
      {
	if ( dest == id ) continue;
	MPI_Wait ( recvrequest + dest, &status );
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int fanum = (*distantdofs[dest])[ii++];
	    int isdistghost = (*distantdofs[dest])[ii++];
	    Array<int> dnums;
	    GetFaceDofNrs (fanum, dnums);
	    for ( int i=0; i<dnums.Size(); i++)
	      {
// 		(*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
// 		(*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];
		(*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
		if ( dest < id && !isdistghost )
		  paralleldofs->ismasterdof.Clear ( dnums[i] ) ;
		ii++; cnt_nexdof[dest]++;
	      }
	  }
      }



    // *****************
    // Parallel Element dofs
    // *****************

    for ( int dest = 0; dest < ntasks; dest++)
      {
	owndofs[dest]->SetSize(1);
	distantdofs[dest]->SetSize(0);
      }
    for ( int el = 0; el < ma.GetNE(); el++ )
      if ( parallelma->IsExchangeElement ( el ) )
	{
	  Array<int> dnums;
	  GetInnerDofNrs ( el, dnums );
	  if ( dnums.Size() == 0 ) continue;

	  for ( int i=0; i<dnums.Size(); i++ )
	    (*(paralleldofs->sorted_exchangedof))[id][exdof++] = dnums[i];

	  for ( int dest = 1; dest < ntasks; dest++ )
	    {
	      int distel = parallelma -> GetDistantElNum ( dest, el );
	      if( distel < 0 ) continue;

	      owndofs[dest]->Append ( distel );
	      owndofs[dest]->Append (int(paralleldofs->IsGhostDof(dnums[0])) );

	      for ( int i=0; i<dnums.Size(); i++)
		{
		  paralleldofs->SetExchangeDof ( dest, dnums[i] );
		  paralleldofs->SetExchangeDof ( dnums[i] );
		  owndofs[dest]->Append ( dnums[i] );
		}
	    }
	}   


    for ( int sender = 1; sender < ntasks; sender ++ )
      {
        if ( id == sender )
          for ( int dest = 1; dest < ntasks; dest ++ )
            if ( dest != id )
              {
		MyMPI_ISend ( *owndofs[dest], dest, sendrequest[dest]);
              }
	  
        if ( id != sender )
          {
	    MyMPI_IRecv ( *distantdofs[sender], sender, recvrequest[sender] );
          }
 	  
	  
      }

    for( int dest=1; dest<ntasks; dest++) 
      {
	if ( dest == id ) continue;
	MPI_Wait ( recvrequest + dest, &status );
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int elnum = (*distantdofs[dest])[ii++];
	    int isdistghost = (*distantdofs[dest])[ii++];
	    Array<int> dnums;
	    GetInnerDofNrs (elnum, dnums);
	    for ( int i=0; i<dnums.Size(); i++)
	      {
// 		(*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
// 		(*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];
		(*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
		if ( dest < id && !isdistghost )
		  paralleldofs->ismasterdof.Clear ( dnums[i] ) ;
		ii++; cnt_nexdof[dest]++;
	      }
	  }
      }


    /*******************************

         update low order space

    *****************************/

    for ( int i = 1; i < ntasks; i++ )
      {
	delete distantdofs[i];
	delete  owndofs[i];
      }

    int ndof_lo = low_order_space->GetNDof();

    // all dofs are exchange dofs
    nexdof = ndof_lo;
 
    exdof = 0;
    cnt_nexdof = 0;


    // *****************
    // Parallel Edge dofs
    // *****************

    owndofs[0]->SetSize(1);
    (*owndofs[0])[0] = ndof;
    distantdofs[0]->SetSize(0);
    
    // find local and distant dof numbers for vertex exchange dofs
    for ( int edge = 0; edge < ma.GetNEdges(); edge++ )
      {
	int dest = 0;
	int distedge = parallelma -> GetDistantEdgeNum ( dest, edge );
	owndofs[0]->Append ( distedge );
	paralleldofs->SetExchangeDof ( dest, edge );
	owndofs[0]->Append ( edge );
      }   
    
    int dest = 0;
    MyMPI_ISend ( *owndofs[0], dest, sendrequest[dest]);
    MyMPI_IRecv ( *distantdofs[0], dest, recvrequest[dest]);
    MPI_Wait ( recvrequest + dest, &status );
	
    paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
    ii = 1;
    while ( ii < distantdofs[0]->Size() )
      {
	int ednum = (*distantdofs[0])[ii++];
// 	(*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = ednum;
// 	(*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[0])[ii];
	(*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = ednum;
	ii++; cnt_nexdof[dest]++;
      }

    for ( int dest = id+1; dest < ntasks; dest++ )
      QuickSort ( (*(paralleldofs->sorted_exchangedof))[dest] );

    for ( int i = 0; i < 1; i++ )
      {
	delete distantdofs[i];
	delete  owndofs[i];
      }

    delete [] owndofs;
    delete [] distantdofs;
    delete [] sendrequest;
    delete [] recvrequest;

  }
#endif

  void HCurlHighOrderFESpace :: UpdateParallelDofs_loproc()
  {
    *testout << "HCurlHo::UpdateParallelDofs_loproc" << endl;
    // lo-proc should never use high order space
    return;
 
  }
#endif // PARALLEL


  // register FESpaces
  namespace hcurlhofespace_cpp
  {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("hcurlho", HCurlHighOrderFESpace::Create);
    }
    
    Init init;
  }
}
