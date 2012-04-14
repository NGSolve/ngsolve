
#include <comp.hpp>
#include <fem.hpp> 

#include "../fem/facethofe.hpp"

namespace ngcomp
{


  FacetFESpace ::  FacetFESpace (const MeshAccess & ama, const Flags & flags, bool checkflags)
    : FESpace(ama, flags)
  {
    name="FacetFESpace(facet)";
    // defined flags
    DefineNumFlag("relorder");
    DefineDefineFlag("print"); 
    DefineDefineFlag("variableorder"); 


    if(checkflags) CheckFlags(flags);
    
    print = flags.GetDefineFlag("print"); 
    
    ndlevel.SetSize(0);
    Flags loflags;
    loflags.SetFlag("order",0.0);
    if ( this->IsComplex() )
      loflags.SetFlag("complex");

    low_order_space = 0;

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
	  cerr << " WARNING: FacetFESpace: inconsistent flags: variableorder, order and relorder "
	       << "-> variable order space with rel_order " << rel_order << "is used, but order is ignored " << endl; 
	else 
	  cerr << " WARNING: FacetFESpace: inconsistent flags: order and rel_order "
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

    highest_order_dc = flags.GetDefineFlag("highest_order_dc");
    if (order == 0)
      highest_order_dc = false;

    
    // TODO: Evaluator for shape tester 
    static ConstantCoefficientFunction one(1);
    if (ma.GetDimension() == 2)
      {
        evaluator = new MassIntegrator<2> (&one);
        boundary_evaluator = new RobinIntegrator<2> (&one);
      }
    else
      {
        evaluator = new MassIntegrator<3> (&one);
        boundary_evaluator = new RobinIntegrator<3> (&one);
      }

    if (dimension > 1)
      {
        evaluator = new BlockBilinearFormIntegrator (*evaluator, dimension);
        boundary_evaluator =
	  new BlockBilinearFormIntegrator (*boundary_evaluator, dimension);
      }
  }

  

  FacetFESpace :: ~FacetFESpace ()
  { ; }




  void FacetFESpace :: Update(LocalHeap & lh)
  {
    FESpace :: Update (lh);
    
    if(print) 
      *testout << " FacetFEspace with order " << order << " rel_order " << rel_order << " var_order " << var_order << endl; 

    if (low_order_space)
      low_order_space -> Update(lh);

    // number of facets
    nel = ma.GetNE();
    nfa = ma.GetNFacets(); // (ma.GetDimension() == 2 ? ma.GetNEdges() : ma.GetNFaces()); 
    
    int p = 0; 
    if(!var_order) p = order; 
    
    order_facet.SetSize(nfa);
    order_facet = INT<2> (p,p);

    fine_facet.SetSize(nfa);
    fine_facet = 0; 
    
    Array<int> fanums;
        
    for (int i = 0; i < nel; i++)
      {
	INT<3> el_orders = ma.GetElOrders(i); 

	ELEMENT_TYPE eltype=ma.GetElType(i); 
	const POINT3D * points = ElementTopology :: GetVertices (eltype);	
	
	if (ma.GetDimension() == 2)
	  {
	    ma.GetElEdges (i, fanums);
	    for (int j=0;j<fanums.Size();j++) 
	      fine_facet[fanums[j]] = 1; 
	    
	    if(var_order)
	      {
		const EDGE * edges = ElementTopology::GetEdges (eltype);
		for(int j=0; j<fanums.Size(); j++)
		  for(int k=0;k<2;k++)
		    if(points[edges[j][0]][k] != points[edges[j][1]][k])
		      { 
			order_facet[fanums[j]] = INT<2>(max(el_orders[k]+rel_order, order_facet[fanums[j]][0]),0);
			break; 
		      }
	      }
	  }
	else
	  {
	    Array<int> elfaces,vnums;
	    ma.GetElFaces(i,elfaces);
	    for (int j=0;j<elfaces.Size();j++) fine_facet[elfaces[j]] = 1; 
	    
	    if(var_order) 
	      {
		ma.GetElVertices (i, vnums);
		const FACE * faces = ElementTopology::GetFaces (eltype);
		for(int j=0;j<elfaces.Size();j++)
		  {
		    if(faces[j][3]==-1) // trig  
		      {
			order_facet[elfaces[j]][0] = max(order_facet[elfaces[j]][0],el_orders[0]+rel_order);
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
				order_facet[elfaces[j]][l] = max(order_facet[elfaces[j]][l], rel_order + el_orders[k]);
				break; 
			      } 
			
		      }
		  }
	      }
	    
	  }
	
      }
    
    // distribute dofs
    ncfa = 0; 
    ndof = nfa; // low_order space
    
    first_facet_dof.SetSize(nfa+1); 
    first_facet_dof = nfa;
      
    if (ma.GetDimension() == 2)
      {
        for (int i = 0; i < nfa; i++)
          {
            first_facet_dof[i] = ndof;
            ndof += order_facet[i][0];
	    if (highest_order_dc) ndof--;
          }
        first_facet_dof[nfa] = ndof;

	if (highest_order_dc)
	  {
	    int ne = ma.GetNE();
	    first_inner_dof.SetSize(ne+1);
	    for (int i = 0; i < ne; i++)
	      {
		first_inner_dof[i] = ndof;
		ELEMENT_TYPE eltype = ma.GetElType(i);
		if (eltype == ET_TRIG)
		  ndof += 3;
		else
		  ndof += 4;
	      }
	    first_inner_dof[ne] = ndof;
	  }
	
      } // 2D
    else  // 3D
      {
        int inci = 0;
        Array<int> pnums;
        for (int i=0; i< nfa; i++)
          {
            int p = order_facet[i][0];
	    if (highest_order_dc) p--;
            ma.GetFacePNums(i,pnums);

            switch(pnums.Size())
              {
              case 3: inci = ((p+1)*(p+2))/2 - 1; break;
              case 4: inci= (p+1)*(p+1) - 1; break;
              }
            first_facet_dof[i] = ndof;
            ndof+= inci;
          }
        first_facet_dof[nfa] = ndof;

	if (highest_order_dc)
	  {
	    int ne = ma.GetNE();
	    first_inner_dof.SetSize(ne+1);
	    for (int i = 0; i < ne; i++)
	      {
		first_inner_dof[i] = ndof;
		
		ELEMENT_TYPE eltype = ma.GetElType(i);
		for (int k = 0; k < ElementTopology::GetNFacets(eltype); k++)
		  if (ElementTopology::GetFacetType(eltype, k) == ET_TRIG)
		    ndof += order+1;
		  else
		    ndof += 2*order+1;
		
		// ndof += 4*(order+1);
	      }
	    first_inner_dof[ne] = ndof;
	  }
      } // 3D


    while (ma.GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;
      
    if(print)
      {
	*testout << "*** Update FAcetFESpace: General Information" << endl;
	*testout << " order facet (facet) " << order_facet << endl; 
	*testout << " first_facet_dof (facet)  " << first_facet_dof << endl; 
      } 

    UpdateCouplingDofArray();

    if (timing) Timing();
  }


  void FacetFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(ndof);
    ctofdof = INTERFACE_DOF;

    for (int facet = 0; facet < ma.GetNFacets(); facet++)
      {
	ctofdof[facet] = WIREBASKET_DOF;
	ctofdof[GetFacetDofs(facet)] = INTERFACE_DOF;
      }

    if (highest_order_dc)
      ctofdof.Range(first_inner_dof[0], ndof) = LOCAL_DOF;
    
    if (print)
      *testout << "FacetFESpace, ctofdof = " << endl << ctofdof << endl;
  }


  // ------------------------------------------------------------------------
  const FiniteElement & FacetFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    FacetVolumeFiniteElement<2> * fe2d = NULL;
    FacetVolumeFiniteElement<3> * fe3d = NULL;;

    switch (ma.GetElType(elnr))
      {
      case ET_SEGM:    break;
      case ET_TRIG:    fe2d = new (lh) FacetFE<ET_TRIG> (); break;
      case ET_QUAD:    fe2d = new (lh) FacetFE<ET_QUAD> (); break;
      case ET_TET:     fe3d = new (lh) FacetFE<ET_TET> (); break;
      case ET_PYRAMID: fe3d = new (lh) FacetFE<ET_PYRAMID> (); break;
      case ET_PRISM:   fe3d = new (lh) FacetFE<ET_PRISM> (); break;
      case ET_HEX:     fe3d = new (lh) FacetFE<ET_HEX> (); break;
      }

    if (!fe2d && !fe3d)
      {
        stringstream str;
        str << "FacetFESpace " << GetClassName() 
            << ", undefined eltype " 
            << ElementTopology::GetElementName(ma.GetElType(elnr))
            << ", order = " << order << endl;
        throw Exception (str.str());
      }

    ArrayMem<int, 12> vnums;
    ArrayMem<int, 6> fanums, order_fa;
    
    ma.GetElVertices(elnr, vnums);
    ma.GetElFacets (elnr, fanums);

    order_fa.SetSize(fanums.Size());
    for (int j = 0; j < fanums.Size(); j++)
      order_fa[j] = order_facet[fanums[j]][0]; //SZ not yet anisotropric
    
    if (ma.GetDimension() == 2)
      {
        fe2d -> SetVertexNumbers (vnums);
        fe2d -> SetOrder (order_fa);
        fe2d -> ComputeNDof();
        return *fe2d;
      }
    else
      {
        fe3d -> SetVertexNumbers (vnums);
        fe3d -> SetOrder (order_fa);
        fe3d -> ComputeNDof();
        return *fe3d;
      }
  }


  // ------------------------------------------------------------------------
  const FiniteElement & FacetFESpace :: GetSFE (int selnr, LocalHeap & lh) const
  {
    L2HighOrderFiniteElement<1> * fe1d = 0;
    L2HighOrderFiniteElement<2> * fe2d = 0;

    switch (ma.GetSElType(selnr))
      {
      case ET_SEGM: fe1d = new (lh) L2HighOrderFE<ET_SEGM> (); break;
      case ET_TRIG: fe2d = new (lh) L2HighOrderFE<ET_TRIG> (); break;
      case ET_QUAD: fe2d = new (lh) L2HighOrderFE<ET_QUAD> (); break;
      default:
        throw Exception (string("FacetFESpace::GetSFE: unsupported element ")+
                         ElementTopology::GetElementName(ma.GetSElType(selnr)));
      }
     
    ArrayMem<int,4> vnums;
    ArrayMem<int,4> ednums;
    
    ma.GetSElVertices(selnr, vnums);
    switch (ma.GetSElType(selnr))
      {
      case ET_SEGM:
	{
	  fe1d -> SetVertexNumbers (vnums);
	  ma.GetSElEdges(selnr, ednums);
	  int p = order_facet[ednums[0]][0];
	  if (highest_order_dc) p--;
	  fe1d -> SetOrder (p); 
	  fe1d -> ComputeNDof();
	  return *fe1d;
	  break;
	}
      case ET_TRIG: 
      case ET_QUAD: 
	{
	  fe2d -> SetVertexNumbers (vnums);
	  int p = order_facet[ma.GetSElFace(selnr)][0];
	  if (highest_order_dc) p--;
	  fe2d -> SetOrder (p);   // SZ not yet anisotropic order for facet fe !!! 
	  fe2d -> ComputeNDof();
	  return *fe2d;
	  break;
	}
      default:
        break;
      }
    return *fe2d;
  }



  // ------------------------------------------------------------------------
  int FacetFESpace :: GetNDof () const
  {
    return ndof;
  }

  // ------------------------------------------------------------------------
  int FacetFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }



  // ------------------------------------------------------------------------
  void FacetFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    ArrayMem<int,6> fanums;
    ma.GetElFacets (elnr, fanums);

    dnums.SetSize(0);

    if (!highest_order_dc)
      {
	for(int i=0; i<fanums.Size(); i++)
	  {
	    dnums.Append(fanums[i]); // low_order
	    dnums += GetFacetDofs (fanums[i]);
	  }
      }
    else
      {
	int innerdof = first_inner_dof[elnr];
	ELEMENT_TYPE eltype = ma.GetElType (elnr);

	for(int i=0; i<fanums.Size(); i++)
	  {
	    int facetdof = first_facet_dof[fanums[i]];

	    if (ma.GetDimension()==2)
	      {
		for (int j = 0; j <= order; j++)
		  {
		    if (j == 0) dnums.Append (fanums[i]);
		    else if (j == order) dnums.Append (innerdof++);
		    else dnums.Append (facetdof++);
		  }
	      }
	    else
	      {
		if (ElementTopology::GetFacetType(eltype, i) == ET_TRIG)
		  for (int j = 0; j <= order; j++)
		    for (int k = 0; k <= order-j; k++)
		      {
			if (j + k == 0) dnums.Append (fanums[i]);
			else if (j + k == order) dnums.Append (innerdof++);
			else dnums.Append (facetdof++);
		      }
		else
		  for (int j = 0; j <= order; j++)
		    for (int k = 0; k <= order; k++)
		      {
			if (j + k == 0) dnums.Append (fanums[i]);
			else if (j == order || k == order) dnums.Append (innerdof++);
			else dnums.Append (facetdof++);
		      }
	      }
		  
	  }
      }
  }

  // ------------------------------------------------------------------------
  void FacetFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    
    int fnum = 0;
    if (ma.GetDimension() == 2)
      {
        ArrayMem<int, 4> fanums;
        ma.GetSElEdges (selnr, fanums);
	fnum = fanums[0];
      }
    else
      fnum = ma.GetSElFace(selnr);
      
    dnums.Append (fnum);
    dnums += GetFacetDofs(fnum);
  }

  // ------------------------------------------------------------------------
  Table<int> * FacetFESpace :: CreateSmoothingBlocks (const Flags & precflags) const
  {
    int ncnt;

    // 1 x low order + faces/edges
    ncnt = nfa-ncfa;
      
    Array<int> cnt(ncnt);
    cnt = 0;

    // setup table
    for (int i = ncfa; i < nfa; i++)
      cnt[i-ncfa] = 1 + first_facet_dof[i+1]-first_facet_dof[i];


    Table<int> & table = *new Table<int> (cnt);
    
    // face/edges
    int ii;
    for (int i = ncfa; i < nfa; i++)
      {
        table[i-ncfa][0] = i-ncfa;
        ii=1;
        for (int j = first_facet_dof[i]; j < first_facet_dof[i+1]; j++, ii++)
          table[i][ii] = j;
      }
      
    // cout << "smoothingblocks = " << endl << table << endl;
    return &table;

  }


  
  Array<int> * FacetFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    Array<int> & clusters = *new Array<int> (GetNDof());

    clusters.SetSize(ndof);
    clusters = 0;
    
    for (int i = 0; i < nfa; i++)
      clusters[i] = 1;
    
    return &clusters;
    
    //
  
    for (int i=0; i<nfa-ncfa; i++)
      clusters[i]  = 1;
  
    // cout << "direct solver cluster = " << clusters << endl;
    return & clusters;
  }








class HybridDGFESpace : public CompoundFESpace
{
public:
  HybridDGFESpace (const MeshAccess & ama, 
                   // const Array<FESpace*> & aspaces,
                   const Flags & flags)
    : CompoundFESpace (ama, flags)
  { 
    Flags l2flags(flags), facetflags(flags);

    int order = int (flags.GetNumFlag ("order", 1));
    
    l2flags.SetFlag ("orderinner", order);
    if (flags.GetDefineFlag("l2_dofs_together")){
      l2flags.SetFlag ("all_dofs_together");
      cout << "l2_dofs_together active" << endl; 
    }

    facetflags.SetFlag("orderfacet", order);
    if (flags.NumListFlagDefined ("dirichlet"))
	facetflags.SetFlag ("dirichlet", flags.GetNumListFlag ("dirichlet"));

    if (flags.NumFlagDefined ("relorder")) facetflags.SetFlag("variableorder");
    
    const FESpaceClasses::FESpaceInfo * info;
    info = GetFESpaceClasses().GetFESpace("DGhotp");
    if (!info) info = GetFESpaceClasses().GetFESpace("l2hotp");
    if (!info) info = GetFESpaceClasses().GetFESpace("l2ho");
    
    // spaces[0] = info->creator(ma, l2flags);
    AddSpace (info->creator(ma, l2flags));
    
    // spaces[0] = new L2HighOrderFESpace (ma, l2flags);    

    // spaces[1] = new FacetFESpace (ma, facetflags);        
    AddSpace (new FacetFESpace (ma, facetflags));        


    if (flags.GetDefineFlag ("edges"))
      throw Exception ("HDG fespace with edges not supported");

    static ConstantCoefficientFunction one(1);
    if (ma.GetDimension() == 2)
      {
        evaluator = new MassIntegrator<2> (&one);
        boundary_evaluator = new RobinIntegrator<2> (&one);
      }
    else
      {
        evaluator = new MassIntegrator<3> (&one);
        boundary_evaluator = new RobinIntegrator<3> (&one);
      }
    evaluator = new CompoundBilinearFormIntegrator (*evaluator, 0);
    boundary_evaluator = new CompoundBilinearFormIntegrator (*boundary_evaluator, 1);
  }

  virtual ~HybridDGFESpace () { ; }

  /*
  static FESpace * Create (const MeshAccess & ma, const Flags & flags)
  {

    HybridDGFESpace * fes = new HybridDGFESpace (ma, spaces, flags);
    return fes;
  }
  */


  virtual Array<int> * CreateDirectSolverClusters (const Flags & flags) const
  {
    if (flags.GetDefineFlag("subassembled"))
      {
	cout << "creating bddc-coarse grid(vertices)" << endl;
	Array<int> & clusters = *new Array<int> (GetNDof());
	clusters = 0;
	/*
	int nv = ma.GetNV();
	int ned = ma.GetNEdges();
	// int nfa = ma.GetNFaces();
	int basefac = spaces[0]->GetNDof();;
	int baseh1 = basefac + spaces[1]->GetNDof();
	  
	Array<int> dnums;
	//low order: 2D: vertices
	if (ma.GetDimension() == 2 && withedges)
	  for (int i = 0; i < nv; i++)
	    if (!IsDirichletVertex(i) && spaces.Size()>2){
	      spaces[2]->GetVertexDofNrs(i,dnums);
	      clusters[dnums[0]+baseh1]=1;
	    }
	    
	//low order: l.o. edges (2D: from facet-space, 3D: from edge-space)
	if (ma.GetDimension() == 2 || ((ma.GetDimension() == 3) && withedges))
	  for (int i = 0; i < ned; i++)
	    if (!IsDirichletEdge(i))
	    {
	      dnums.SetSize(0);
	      if (ma.GetDimension() == 2){
// 		spaces[1]->GetEdgeDofNrs(i,dnums);
// 		clusters[dnums[0]+basefac]=1;
	      }else{
		spaces[2]->GetEdgeDofNrs(i,dnums);
		clusters[dnums[0]+baseh1]=1;
	      }
	    }
	*/


	/*
	//low order: 3D: l.o. face
	if (ma.GetDimension() == 3)
	  for (int i = 0; i < nfa; i++)
	    if (!IsDirichletFace(i))
	    {
		dnums.SetSize(0);
		spaces[1]->GetFaceDofNrs(i,dnums);
// 		for (int l=0; l<dnums.Size(); l++)
// 		  clusters[dnums[l]+basefac]=1;		
		clusters[dnums[0]+basefac]=1;
	      //end-if isdirichletvertex
	    }
	*/

	return &clusters;	
    }
    else
      {
	Array<int> & clusters = *new Array<int> (GetNDof());
	clusters = 0;

	Array<int> dnums;
	int nfa = ma.GetNFacets();

	for (int i = 0; i < nfa; i++)
	  {
	    if (ma.GetDimension() == 2)
	      GetEdgeDofNrs (i, dnums);
	    else
	      GetFaceDofNrs (i, dnums);

	    clusters[dnums[0]] = 1;
	  }

	const BitArray & freedofs = *GetFreeDofs();
	for (int i = 0; i < freedofs.Size(); i++)
	  if (!freedofs.Test(i)) clusters[i] = 0;
	*testout << "Hybrid-FESpace, dsc = " << endl << clusters << endl;
	return &clusters;
      }
  }

  virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const
  {
    bool eliminate_internal = precflags.GetDefineFlag("eliminate_internal");
    bool subassembled = precflags.GetDefineFlag("subassembled");
    int smoothing_type = int(precflags.GetNumFlag("blocktype",1)); 
    COUPLING_TYPE dof_mode = eliminate_internal? (subassembled? WIREBASKET_DOF : EXTERNAL_DOF) : ANY_DOF;
    BitArray filter;
    GetFilteredDofs(dof_mode, filter, true);
    
    int nv = ma.GetNV();
    int ned = ma.GetNEdges();
    cout << " dof_mode " << dof_mode << endl; 
    cout << " blocktype " << smoothing_type << endl; 
    cout << " Use HDG-Block Smoother:  "; 
    Array<int> dnums;

    FilteredTableCreator creator(&filter);
    for ( ; !creator.Done(); creator++)
      {
	switch (smoothing_type)
	  {
	  case 1: 
	  //for BDDC: we have only the condensated (after subassembling) dofs, 
	  //and build patches around each vertex Vertices + Edges
		
	    if (creator.GetMode() == 1)
	      cout << "BDDC-Edges-around-Vertex-Block" << endl;
		
// 	    int ds_order = precflags.GetNumFlag ("ds_order", 1);
// 	    cout << "ds_order = " << ds_order << endl;
// 		
	    if (ma.GetDimension() == 2)
	      {  
		for (int i = 0; i < nv; i++)
		{
		  dnums.SetSize(0);
		  GetVertexDofNrs(i,dnums);
		  if (dnums.Size()>0)
		    creator.Add (i, dnums[0]);
		}
	      }
	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma.GetNode<1> (i);
		for (int k = 0; k < 2; k++){
		  dnums.SetSize(0);
		  if (ma.GetDimension() == 2){
		    GetEdgeDofNrs(i,dnums);
		    creator.Add (edge.vertices[k], dnums[0]);
		  }
		}
	      }
	    break; 	    
	  case 2: 
	    //for BDDC: we have only the condensated (after subassembling) dofs, 
	    //and build patches around each edge: [Faces?!] + Edges
		
	    if (creator.GetMode() == 1)
	      cout << "BDDC-Faces-around-Edges" << endl;
		
	    if (ma.GetDimension() == 2)
	      {
		;
	      }

	    else
	      
	      {
		Array<int> dnums, ednums;

		for (int i = 0; i < ned; i++)
		  if (!IsDirichletEdge(i))
		    {
		      GetEdgeDofNrs(i,dnums);
		      for (int l=0;l<dnums.Size();l++)
			creator.Add (i, dnums[l]);
		    }
	      }
	      break; 	    


	  case 3: 
	    // facet-by-facet
		
	    if (creator.GetMode() == 1)
	      cout << "Facet-by-facet blocks" << endl;
		
	    
	    Array<int> dnums;
	    int nfa = ma.GetNFacets();
	    for (int i = 0; i < nfa; i++)
	      {
		if (ma.GetDimension() == 2)
		  {
		    if (IsDirichletEdge (i)) continue;
		    GetEdgeDofNrs (i, dnums);
		  }
		else
		  {
		    if (IsDirichletFace (i)) continue;
		    GetFaceDofNrs (i, dnums);
		  }

		for (int l = 0; l < dnums.Size(); l++)
		  creator.Add (i, dnums[l]);
	      }
	    break; 	    
	  }
      }
    return creator.GetTable();
  }




  virtual string GetClassName () const { return "HybridDGFESpace"; }
};


  
// ------------------------------------------------------------------------

  
  static RegisterFESpace<FacetFESpace> init_facet ("facet");
  static RegisterFESpace<HybridDGFESpace> init_hde ("HDG");
}


