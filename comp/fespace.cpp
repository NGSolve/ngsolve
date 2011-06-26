/**********************************************************************/
/* File:   fespace.cpp                                                */
/* Author: Joachim Schoeberl                                          */
/* Date:   25. Mar. 2000                                              */
/**********************************************************************/

/* 
   Finite Element Space
*/

#include <comp.hpp>
#include <multigrid.hpp>

#include "../fem/h1lofe.hpp"
#include <parallelngs.hpp>



using namespace ngmg;

namespace ngcomp
{
  using namespace ngcomp;


  FESpace :: FESpace (const MeshAccess & ama, const Flags & flags, bool checkflags)
    : NGS_Object (ama, "FESpace")
  {
    // register flags
    DefineStringFlag("type");
    DefineNumFlag("order");
    DefineNumFlag("dim");
    DefineDefineFlag("vec");
    DefineDefineFlag("complex");
/*    DefineDefineFlag("eliminate_internal");
    DefineDefineFlag("noeliminate_internal");*/
    DefineDefineFlag("timing");
    DefineNumListFlag("directsolverdomains");
    DefineNumListFlag("dirichlet");
    DefineNumListFlag("definedon");
    DefineNumFlag ("definedon");
    DefineStringListFlag ("definedon");
    DefineNumListFlag("definedonbound");
    DefineNumFlag ("definedonbound");
    DefineStringListFlag ("definedonbound");
    DefineDefineFlag("dgjumps");

    order = int (flags.GetNumFlag ("order", 1));
    dimension = int (flags.GetNumFlag ("dim", 1));
    iscomplex = flags.GetDefineFlag ("complex");
//     eliminate_internal = flags.GetDefineFlag("eliminate_internal");
    timing = flags.GetDefineFlag("timing");
    dgjumps = flags.GetDefineFlag("dgjumps");
    if (dgjumps) 
      *testout << "ATTENTION: flag dgjumps is used!\n This leads to a \
lot of new non-zero entries in the matrix!\n" << endl;
    else *testout << "\n (" << order << ") flag dgjumps is not used!" << endl;
    
    if(flags.NumListFlagDefined("directsolverdomains"))
      {
	directsolverclustered.SetSize(ama.GetNDomains());
	directsolverclustered = false;
	Array<double> clusters(flags.GetNumListFlag("directsolverdomains"));
	for(int i=0; i<clusters.Size(); i++) 
	  directsolverclustered[static_cast<int>(clusters[i])-1] = true; // 1-based!!
      }
    
    if(flags.NumListFlagDefined("dirichlet"))
      {
	dirichlet_boundaries.SetSize (ma.GetNBoundaries());
	dirichlet_boundaries.Clear();
	Array<double> db (flags.GetNumListFlag("dirichlet"));
	for(int i = 0; i< db.Size(); i++) 
	  {
	    int bnd = int(db[i]-1);
	    if (bnd >= 0 && bnd < dirichlet_boundaries.Size())
	      dirichlet_boundaries.Set (int(db[i])-1);
	    /*
	    else
	      cerr << "Illegal Dirichlet boundary index " << bnd+1 << endl;
	    */
	  }
        *testout << "dirichlet_boundaries" << dirichlet_boundaries << endl;
      }
    
    if(flags.NumListFlagDefined("definedon") || 
       flags.NumFlagDefined("definedon") ||
       flags.StringListFlagDefined("definedon"))
      {
	definedon.SetSize (ma.GetNDomains());
	definedon = 0;
	Array<double> defon;
	if (flags.NumListFlagDefined("definedon")) 
	  defon = flags.GetNumListFlag("definedon");
	else if (flags.NumFlagDefined("definedon"))
	  {
	    defon.SetSize(1);
	    defon[0] = flags.GetNumFlag("definedon",0);
	  }
	for(int i=0; i< defon.Size(); i++)
	  if(defon[i] <= ma.GetNDomains() && defon[i] > 0)
	    definedon[int(defon[i])-1] = 1;

	if(flags.StringListFlagDefined("definedon"))
	  {
	    Array<string> dmaterials(flags.GetStringListFlag ("definedon").Size());
	    for(int i=0; i<dmaterials.Size(); i++)
	      dmaterials[i] = flags.GetStringListFlag ("definedon")[i];
	    for(int i = 0; i < ma.GetNDomains(); i++)
	      {
		for(int j = 0; definedon[i] == 0 && j < dmaterials.Size(); j++)
		  if(StringFitsPattern(ma.GetDomainMaterial(i),dmaterials[j]))
		    definedon[i] = 1;
	      }
	  }

	// default:
	// fespace only defined on boundaries matching definedon-domains
	definedonbound.SetSize (ma.GetNBoundaries());
	definedonbound = 0;
	for ( int sel=0; sel<ma.GetNSE(); sel++ )
	  {
	    int index = ma.GetSElIndex(sel);
	    int dom1, dom2;
	    ma.GetSElNeighbouringDomains(sel, dom1, dom2);
	    dom1--; dom2--;
	    if ( dom1 >= 0 )
	      if ( definedon[dom1] )
		definedonbound[index] = 1;

	    if ( dom2 >= 0 )
	      if ( definedon[dom2] )
		definedonbound[index] = 1;
	  }
      }

    // additional boundaries
    if(flags.NumListFlagDefined("definedonbound")|| flags.NumFlagDefined("definedonbound") )
      {
	if ( definedonbound.Size() == 0 )
	  {
	    definedonbound.SetSize (ma.GetNBoundaries());
	    definedonbound = 0;
	  }
	Array<double> defon;
	if ( flags.NumListFlagDefined("definedonbound") )
	  defon = (flags.GetNumListFlag("definedonbound"));
	else
	  {
	    defon.SetSize(1);
	    defon[0] = flags.GetNumFlag("definedonbound",0);
	  }

	for(int i=0; i< defon.Size(); i++) 
	  definedonbound[int(defon[i])-1] = 1;
      }
    

    else if(flags.StringListFlagDefined("definedonbound") || flags.StringFlagDefined("definedonbound"))
      {
	if ( definedonbound.Size() == 0 )
	  {
	    definedonbound.SetSize (ma.GetNBoundaries());
	    definedonbound = 0;
	  }

	Array<string*> defon;

	if(flags.StringFlagDefined("definedonbound"))
	  defon.Append(new string(flags.GetStringFlag("definedonbound","")));
	else
	  for(int i=0; i<flags.GetStringListFlag ("definedonbound").Size(); i++)
	    defon.Append(new string(flags.GetStringListFlag("definedonbound")[i]));
	
	for(int selnum = 0; selnum < ma.GetNSE(); selnum++)
	  {
	    if(definedonbound[ma.GetSElIndex(selnum)]!=1)
	      {
		for(int i=0; i<defon.Size(); i++)
		  {
		    if(StringFitsPattern(ma.GetSElBCName(selnum),*(defon[i])))	
		      {		
		 	definedonbound[ma.GetSElIndex(selnum)] =1;
			continue;
		      }
		  }
	      }
	  }
	for(int i=0; i<defon.Size(); i++)
	  delete defon[i];
	
      }
    

    
    prol = 0;
    low_order_space = 0;
    is_low_order_space = false;

    tet = 0;
    pyramid = 0;
    prism = 0;
    hex = 0;
    trig = 0;
    quad = 0;
    segm = 0;

    evaluator = 0;
    boundary_evaluator = 0;

    // first_lodof[4] = -1;   // indicates that nodes are not used

    element_coloring = NULL;

    paralleldofs = NULL;

    ctofdof.SetSize(0);
  }

  
  FESpace :: ~FESpace ()
  {
    delete low_order_space;
    delete boundary_evaluator;
    delete evaluator;
    delete prol;

    delete tet;
    delete pyramid;
    delete prism;
    delete hex;
    delete trig;
    delete quad;
    delete segm;

    delete element_coloring;
    delete paralleldofs;
  }
  

  void FESpace :: Update(LocalHeap & lh)
  {
    *testout << "Update FESpace, type = " << typeid(*this).name() << endl;
    *testout << "name = " << name << endl;

    for (int i=0; i< specialelements.Size(); i++)
      delete specialelements[i]; 
    specialelements.SetSize(0);


      {
	int dim = ma.GetDimension();

	dirichlet_vertex.SetSize (ma.GetNV());
	dirichlet_edge.SetSize (ma.GetNEdges());
        if (dim == 3)
          dirichlet_face.SetSize (ma.GetNFaces());
	
	dirichlet_vertex = false;
	dirichlet_edge = false;
	dirichlet_face = false;

	if (dirichlet_boundaries.Size())
	  for (int i = 0; i < ma.GetNSE(); i++)
	    {
	      int ind = ma.GetSElIndex (i);
	      if (dirichlet_boundaries.Test(ind))
		{
		  Ng_Element ngel = ma.GetSElement(i);
		  
		  for (int j = 0; j < ngel.vertices.Size(); j++)
		    dirichlet_vertex[ngel.vertices[j]] = true;
		  
		  for (int j = 0; j < ngel.edges.Size(); j++)
		    dirichlet_edge[ngel.edges[j]] = true;
		  
		  if (dim == 3)
		    dirichlet_face[ngel.faces[0]] = true;
		  // ma.GetSElVertices (i, vnums);
		  // ma.GetSElEdges (i, ednums);
		  // if (dim == 3)
		  // fanum = ma.GetSElFace (i);
		  
		  // for (int j = 0; j < vnums.Size(); j++)
		  // dirichlet_vertex[vnums[j]] = true;
		  // for (int j = 0; j < ednums.Size(); j++)
		  // dirichlet_edge[ednums[j]] = true;			
		  
		  // if (dim == 3)
		  // dirichlet_face[fanum] = true;
		}
	    }


	
	// cout << "exchange dirichlet vertices" << endl;
	
	DynamicTable<int> dist_dir_vertex(ntasks);
	Array<int[2]> distnums;
	if (id != 0)
	  {
	    for (int i = 0; i < dirichlet_vertex.Size(); i++)
	      if (dirichlet_vertex[i])
		{
		  ma.GetDistantNodeNums (NT_VERTEX, i, distnums);
		  for (int j = 0; j < distnums.Size(); j++)
		    if (distnums[j][0] > 0)  // not master
		      dist_dir_vertex.Add (distnums[j][0], distnums[j][1]);
		}
		/*
		for (int dist = 1; dist < ntasks; dist++)
		  if (dist != id)
		    {
		      int distnum = NgPar_GetDistantPNum (dist, i);
		      if (distnum >= 0)
			dist_dir_vertex.Add (dist, distnum);
		    }
		*/
	    *testout << "dist_dir_vertex = " << dist_dir_vertex << endl;

	    for (int dist = 1; dist < ntasks; dist++)
	      if (id != dist)
		{
		  Array<int> dirvert;
		  if (dist < id)
		    {
		      MyMPI_Send (dist_dir_vertex[dist], dist);
		      MyMPI_Recv (dirvert, dist);
		    }
		  else
		    {
		      MyMPI_Recv (dirvert, dist);
		      MyMPI_Send (dist_dir_vertex[dist], dist);
		    }
		  // *testout << "got dirvert from " << dist << ": " << dirvert << endl;
		  for (int i = 0; i < dirvert.Size(); i++)
		    dirichlet_vertex[dirvert[i]] = true;
		}	   

	    /*
	    // check ordering
	    for (int i = 1; i < ma.GetNV(); i++)
	      {
		int m0 =  NgPar_GetDistantPNum (0, i-1);
		int m1 =  NgPar_GetDistantPNum (0, i);
		if (m0 > m1) cerr << endl << endl << "Wrong ordering !!!!" << endl << endl;
	      }
	    */
	  }

	// cout << "exchange dirichlet vertices done" << endl;		

	(*testout) << "Dirichlet_vertex = " << endl << dirichlet_vertex << endl;
	(*testout) << "Dirichlet_edge = " << endl << dirichlet_edge << endl;
	(*testout) << "Dirichlet_face = " << endl << dirichlet_face << endl;
      }
  }

  void FESpace :: FinalizeUpdate(LocalHeap & lh)
  {
    static Timer timer ("FESpace::FinalizeUpdate");

    if (low_order_space) low_order_space -> FinalizeUpdate(lh);

    RegionTimer reg (timer);

      {
	dirichlet_dofs.SetSize (GetNDof());
	dirichlet_dofs.Clear();
	Array<int> dnums;

	if (dirichlet_boundaries.Size())
	  for (int i = 0; i < ma.GetNSE(); i++)
	    {
	      if (dirichlet_boundaries[ma.GetSElIndex(i)])
		{
		  GetSDofNrs (i, dnums);
		  for (int j = 0; j < dnums.Size(); j++)
		    if (dnums[j] != -1)
		      dirichlet_dofs.Set (dnums[j]);
		}
	    }


	free_dofs.SetSize (GetNDof());
	free_dofs = dirichlet_dofs;
	free_dofs.Invert();

	*testout << "freedofs = " << endl << free_dofs << endl;
      }
    

    UpdateParallelDofs();


    *testout << "coloring ... " << flush;

    Array<int> col(ma.GetNE());
    col = -1;

    bool found;
    Array<int> dnums;
    int maxcolor = 0;


    int basecol = 0;
    Array<unsigned int> mask(GetNDof());
    
    do
      {
	mask = 0;
	found = false;
	
	for (int i = 0; i < ma.GetNE(); i++)
	  {
	    if (col[i] >= 0) continue;
	    GetDofNrs (i, dnums);	

	    unsigned check = 0;
	    for (int j = 0; j < dnums.Size(); j++)
	      check |= mask[dnums[j]];

	    if (check != UINT_MAX) // 0xFFFFFFFF)
	      {
		found = true;
		unsigned checkbit = 1;
		int color = basecol;
		while (check & checkbit)
		  {
		    color++;
		    checkbit *= 2;
		  }
		// *testout << "set color = " << color << endl;
		col[i] = color;
		if (color > maxcolor) maxcolor = color;
		
		for  (int j = 0; j < dnums.Size(); j++)
		  mask[dnums[j]] |= checkbit;
	      }
	  }
	
	basecol += 8*sizeof(unsigned int); // 32;
      }
    while (found);

    // int color = maxcolor+1;

    Array<int> cntcol(maxcolor+1);
    cntcol = 0;
    for (int i = 0; i < ma.GetNE(); i++)
      cntcol[col[i]]++;

    delete element_coloring;
    element_coloring = new Table<int> (cntcol);
    cntcol = 0;
    for (int i = 0; i < ma.GetNE(); i++)
      (*element_coloring)[col[i]][cntcol[col[i]]++] = i;
    
    *testout << "needed " << maxcolor+1 << " colors" << endl;
  }


  const FiniteElement & FESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    FiniteElement * fe = 0;
    
    switch (ma.GetElType(elnr))
      {
      case ET_TET: fe = tet; break;
      case ET_PYRAMID: fe = pyramid; break;
      case ET_PRISM: fe = prism; break;
      case ET_HEX: fe = hex; break;
      case ET_TRIG: fe = trig; break;
      case ET_QUAD: fe = quad; break;
      default:
        ;
      }

    if (!fe)
      {
        stringstream str;
        str << "FESpace " << GetClassName() 
            << ", undefined eltype " 
            << ElementTopology::GetElementName(ma.GetElType(elnr))
            << ", order = " << order << endl;
        throw Exception (str.str());
      }
    
    return *fe;
  }

  /*
  void FESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    TopologicElement topel;
    ma.GetTopologicElement (elnr, topel);

    LocalHeapMem<10003> lh("FESpace - GetDofNrs");
    ArrayMem<Dof, 100> dofs;

    GetFE (elnr, lh) . GetDofs(dofs);

    dnums.SetSize(dofs.Size());

    for (int i = 0; i < dofs.Size(); i++)
      {
        Node gn = topel.GlobalNode (dofs[i].GetNode());
        int nr_on_node = dofs[i].GetNrOnNode();

        if (nr_on_node < lodofs_per_node[gn.GetType()])
          dnums[i] = 
            first_lodof[gn.GetType()] + 
            lodofs_per_node[gn.GetType()] * gn.GetNodeNr() +
            nr_on_node;

        else

          dnums[i] = first_hodofs[gn.GetType()][gn.GetNodeNr()] + 
            nr_on_node - lodofs_per_node[gn.GetType()];
      }
  }
  */

  /// get coupling type of dof
  COUPLING_TYPE FESpace :: GetDofCouplingType (int dof) const 
  {
    if (ctofdof.Size()==0) //this is the standard case if the FESpace does not specify anything.
      return WIREBASKET_DOF;
    if (dof >= ctofdof.Size()) throw Exception("FESpace::GetDofCouplingType out of range. Did you set up ctofdof?");
      return ctofdof[dof];
  }

  /// get coupling types of dofs
  void  FESpace :: GetDofCouplingTypes (int elnr, Array<COUPLING_TYPE> & ctypes) const
  { 
    Array<int> dnums;
    GetDofNrs(elnr, dnums);
    ctypes.SetSize(dnums.Size());
    if (ctofdof.Size()==0){
      for (int i=0;i<dnums.Size();i++)
	ctypes[i] = INTERFACE_DOF;
    }
    else
    {
      for (int i=0; i<dnums.Size(); i++)
	ctypes[i] = ctofdof[dnums[i]];
    }
  }

  void FESpace :: GetDofNrs (int elnr, Array<int> & dnums,COUPLING_TYPE ctype) const
  {
      Array<int> alldnums; 
      Array<COUPLING_TYPE> couptype;
      GetDofNrs(elnr, alldnums);
      GetDofCouplingTypes(elnr, couptype);
      dnums.SetSize(0);
      for (int i=0;i<alldnums.Size();i++){
	if ( (couptype[i] & ctype) != 0){
	  dnums.Append(alldnums[i]);
	}
      }
  }

  void FESpace :: GetNodeDofNrs (NODE_TYPE nt, int nr, Array<int> & dnums) const
  {
    /*
    if (first_lodof[4] != -1)
      {
	dnums.SetSize(0);
	for (int i = 0; i < lodofs_per_node[nt]; i++)
	  dnums.Append (first_lodof[nt]+nr*lodofs_per_node[nt] + i);
	if (first_hodofs[nt].Size())
	  for (int i = first_hodofs[nt][nr]; i < first_hodofs[nt][nr+1]; i++)
	    dnums.Append (i);
	return;
      }
    */

    switch (nt)
      {
      case NT_VERTEX: GetVertexDofNrs(nr, dnums); break;
      case NT_EDGE:   GetEdgeDofNrs(nr, dnums); break;
      case NT_FACE:   
        if (ma.GetDimension() == 3)
          GetFaceDofNrs(nr, dnums); 
        else
          GetInnerDofNrs(nr, dnums); 
        break;
      case NT_CELL:   GetInnerDofNrs(nr, dnums); break;
      }
  }

  void FESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    // if (first_lodof[4] != -1)
    // GetNodeDofNrs (NT_VERTEX, vnr, dnums);
    throw Exception ("FESpace::GetVertexDofNrs called");
  }

  void FESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    // if (first_lodof[4] != -1)
    // GetNodeDofNrs (NT_EDGE, ednr, dnums);
    throw Exception ("FESpace::GetEdgeDofNrs called");
  }

  void FESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    // if (first_lodof[4] != -1)
    // GetNodeDofNrs (NT_FACE, fanr, dnums);
    throw Exception ("FESpace::GetFaceDofNrs called");
  }

  void FESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    // if (first_lodof[4] != -1)
    // GetNodeDofNrs (NT_CELL, elnr, dnums);
    throw Exception (string("FESpace::GetInnerDofNrs called for class")+
		     typeid(*this).name());
  }

  /*
  int FESpace :: GetNLowOrderNodeDofs ( NODE_TYPE nt ) const
  {
    throw Exception (string("FESpace::GetNLowOrderNodeDofs called for class")+
		     typeid(*this).name());
  }
  */

  const FiniteElement & FESpace :: GetSFE (int selnr, LocalHeap & lh) const
  {
    FiniteElement * fe = 0;
    
    switch (ma.GetSElType(selnr))
      {
      case ET_TRIG:
	fe = trig; break;
      case ET_QUAD:
	fe = quad; break;
      case ET_SEGM:
	fe = segm; break;
      default:
        ;
      }
    
    if (!fe)
      {
	stringstream str;
	str << "FESpace " << GetClassName() 
	    << ", undefined surface eltype " << ma.GetSElType(selnr) 
	    << ", order = " << order << endl;
	throw Exception (str.str());
      }

    return *fe;
  }
  

  const FiniteElement & FESpace :: GetFE (ELEMENT_TYPE type) const
  {
    switch (type)
      {
      case ET_SEGM:
	return *segm;
      case ET_TRIG:
	return *trig;
      case ET_QUAD:
	return *quad;
      case ET_TET:
	return *tet;
      case ET_PYRAMID:
	return *pyramid;
      case ET_PRISM:
	return *prism;
      case ET_HEX:
	return *hex;
      }
    throw Exception ("GetFE::Unknown type");
  }


  
  void FESpace :: PrintReport (ostream & ost)
  {
    ost << "type  = " << GetClassName() << endl
	<< "order = " << order << endl
	<< "dim   = " << dimension << endl
	<< "dgjmps= " << dgjumps << endl
	<< "complex = " << iscomplex << endl;
  }
  


  int FESpace :: GetNDofLevel (int level) const
  {
    return GetNDof();
  } 

  void FESpace :: Timing () const
  {
    clock_t starttime;
    double time;
    
    starttime = clock();
    
    int steps = 0;
    Array<int> dnums;
    do
      {
        for (int i = 0; i < ma.GetNE(); i++)
          GetDofNrs (i, dnums);
        steps++;
        time = double(clock() - starttime) / CLOCKS_PER_SEC;
      }
    while (time < 2.0);
    
    cout << ma.GetNE()*steps / time << " GetDofNrs per second" << endl;



    starttime = clock();
    steps = 0;
    LocalHeap lh (100000, "FESpace - Timing");
    do
      {
        for (int i = 0; i < ma.GetNE(); i++)
          {
            HeapReset hr(lh);
            GetFE (i, lh);
          }
        steps++;
        time = double(clock() - starttime) / CLOCKS_PER_SEC;
      }
    while (time < 2.0);
    
    cout << ma.GetNE()*steps / time << " GetFE per second" << endl;
  }



  void FESpace::GetFilteredDofs(COUPLING_TYPE doffilter, BitArray & output, bool freedofsonly) const{
    int ndof = GetNDof();
    output.SetSize(ndof);
    output.Clear();
    if (ctofdof.Size()>0)
      for (int i = 0; i < ndof; i++)
	if ((ctofdof[i] & doffilter) != 0)
	  output.Set(i);
    if (freedofsonly && free_dofs.Size()) {
      output.And(free_dofs);
    }
  }



  Table<int> * FESpace :: CreateSmoothingBlocks (const Flags & flags) const
  {
    int nd = GetNDof();
    TableCreator<int> creator;

    for ( ; !creator.Done(); creator++)
      {
	for (int i = 0; i < nd; i++)
	  if (!IsDirichletDof(i))
	    creator.Add (i, i);
      }
    return creator.GetTable();
  }

    
  void FESpace :: SetDefinedOn (const BitArray & defon)
  {
    definedon.SetSize(defon.Size());
    for (int i = 0; i < defon.Size(); i++)
      definedon[i] = defon.Test(i) ? 1 : 0;
  }

  void FESpace :: SetDefinedOnBoundary (const BitArray & defon)
  {
    definedonbound.SetSize(defon.Size());
    for (int i = 0; i < defon.Size(); i++)
      definedonbound[i] = defon.Test(i) ? 1 : 0;
  }

  void FESpace :: SetDirichletBoundaries (const BitArray & dirbnds)
  {
    dirichlet_boundaries = dirbnds;
  }


  const BitArray * FESpace :: GetFreeDofs () const
  {
    if (!free_dofs.Size())
      return NULL;
    else
      return &free_dofs;
  }




  // Aendern, Bremse!!!
  template < int S, class T >
  void FESpace :: TransformVec (int elnr, bool boundary,
				const FlatVector< Vec<S,T> >& vec, TRANSFORM_TYPE type) const
  {
    //cout << "Achtung, Bremse" << endl;

    Vector<T> partvec(vec.Size());

    for(int i=0; i<S; i++)
      {
	for(int j=0;j<vec.Size(); j++)
	  partvec(j) = vec[j](i);

	TransformVec(elnr,boundary,partvec,type);

	for(int j=0;j<vec.Size(); j++)
	  vec[j](i) = partvec(j);
      }
  }

  template void FESpace::TransformVec(int elnr, bool boundary,
				      const FlatVector< Vec<2,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
				      const FlatVector< Vec<3,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
				      const FlatVector< Vec<4,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
				      const FlatVector< Vec<5,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
				      const FlatVector< Vec<6,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
				      const FlatVector< Vec<7,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<8,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<9,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<10,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<11,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<12,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<13,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<14,double> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<15,double> >& vec, TRANSFORM_TYPE type) const;

  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<2,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<3,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<4,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<5,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<6,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<7,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<8,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<9,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<10,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<11,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<12,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<13,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<14,Complex> >& vec, TRANSFORM_TYPE type) const;
  template void FESpace::TransformVec(int elnr, bool boundary,
					const FlatVector< Vec<15,Complex> >& vec, TRANSFORM_TYPE type) const;


  void FESpace :: UpdateParallelDofs ( )
  {
    if  ( ntasks == 1 ) 
      return;


    if ( id == 0 )
      {
	Array<int> nexdof(ntasks); 
	nexdof = 0;
	
	Table<int> * exdofs = new Table<int> (nexdof);
	// paralleldofs = new ParallelDofs (GetNDof(), exdofs, this);
	paralleldofs = new ParallelDofs (0, exdofs, this);
      }

    else
      {
	bool print = true;

	if (print)
	  {
	    *testout << "FESpace::UpdateParallelDofs_hoproc" << endl;
      
	    for (NODE_TYPE nt = NT_VERTEX; nt <= NT_CELL; nt++)
	      for (int node = 0; node < ma.GetNNodes (nt); node++)
		{
		  *testout << "distnode[" << nt << "-" << node << "] = " << endl;
		  Array<int[2]> distantnodenums;
		  // parallelma -> 
		  ma.GetDistantNodeNums (nt, node, distantnodenums);
		  for (int j = 0; j < distantnodenums.Size(); j++)
		    *testout << distantnodenums[j][0] << " " << distantnodenums[j][1] << endl;
		}
	  }
    

	// Find number of exchange dofs
	Array<int> nexdof(ntasks);
	nexdof = 0;

	Array<int> dnums;


	for (NODE_TYPE nt = NT_VERTEX; nt <= NT_CELL; nt++)
	  for ( int nr = 0; nr < ma.GetNNodes (nt); nr++ )
	    {
	      GetNodeDofNrs (nt, nr, dnums);

	      Array<int[2]> distantnodenums;
	      ma.GetDistantNodeNums (nt, nr, distantnodenums);

	      for (int j = 0; j < distantnodenums.Size(); j++)
		{
		  int dest = distantnodenums[j][0];
		  if (dest > 0) nexdof[dest] += dnums.Size(); 
		}
	    }
    
	nexdof[0] = 0;

	Table<int> * exdofs = new Table<int> (nexdof);
	nexdof = 0;

	for (NODE_TYPE nt = NT_VERTEX; nt <= NT_CELL; nt++)
	  {
	    for ( int node = 0; node < ma.GetNNodes(nt); node++ )
	      {
		GetNodeDofNrs (nt, node, dnums); 
		if ( dnums.Size() == 0 ) continue;

		Array<int[2]> distantnodenums;
		ma.GetDistantNodeNums ( nt, node, distantnodenums);
		
		for (int j = 0; j < distantnodenums.Size(); j++)
		  {
		    int dest = distantnodenums[j][0];
		    if (dest > 0)
		      for (int i = 0; i < dnums.Size(); i++)
			(*exdofs)[dest][nexdof[dest]++] = dnums[i];
		  }
	      } 
	  }

	for (int dest = 1; dest < ntasks; dest++ )
	  QuickSort ( (*exdofs)[dest] );

	paralleldofs = new ParallelDofs (GetNDof(), exdofs, this);

	*testout << "sorted_exdof = " << endl << *exdofs << endl;
      }
  }







  NodalFESpace :: NodalFESpace (const MeshAccess & ama,
				const Flags & flags,
                                bool parseflags)
    : FESpace (ama, flags)
  {
    name="NodalFESpace";
    
    prol = new LinearProlongation(*this);

    if (order >= 2)
      {
	Flags loflags;
	loflags.SetFlag ("order", 1);
	loflags.SetFlag ("dim", dimension);
	if (dgjumps) loflags.SetFlag ("dgjumps");
	if (iscomplex) loflags.SetFlag ("complex");
	low_order_space = new NodalFESpace (ma, loflags);
      }

    if (order == 1)
      {
	tet     = new FE_Tet1;
	prism   = new FE_Prism1;
	pyramid = new FE_Pyramid1;
	hex     = new FE_Hex1;
	trig    = new FE_Trig1;
	quad    = new FE_Quad1;
	segm    = new FE_Segm1;
      }
    else
      {
	if (flags.GetDefineFlag ("hb"))
	  {
	    tet     = new FE_Tet2HB;
	    prism   = new FE_Prism1;
	    pyramid = new FE_Pyramid1;
	    trig    = new FE_Trig2HB;
	    quad    = new FE_Quad1;
	    segm    = new FE_Segm2;
	  }
	else
	  {
	    tet     = new FE_Tet2;
	    prism   = new FE_Prism1;
	    pyramid = new FE_Pyramid1;
	    trig    = new FE_Trig2;
	    quad    = new FE_Quad1;
	    segm    = new FE_Segm2;
	  }
      }


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
	boundary_evaluator = new BlockBilinearFormIntegrator (*boundary_evaluator, dimension);  
      }
  }

  NodalFESpace :: ~NodalFESpace ()
  {
    ;
  }

  int NodalFESpace :: GetNDof () const
  {
    return ndlevel.Last();
  }

  void NodalFESpace :: Update(LocalHeap & lh)
  {
    FESpace :: Update (lh);

    if (low_order_space) low_order_space -> Update(lh);

    if (ma.GetNLevels() > ndlevel.Size())
      {
	Array<int> dnums;
	int i, j;
	int ne = ma.GetNE();
	int nse = ma.GetNSE();
	int ndof = 0;
	for (i = 0; i < ne; i++)
	  {
	    GetDofNrs (i, dnums);
	    for (j = 0; j < dnums.Size(); j++)
	      if (dnums[j] > ndof)
		ndof = dnums[j];
	  }
	for (i = 0; i < nse; i++)
	  {
	    GetSDofNrs (i, dnums);
	    for (j = 0; j < dnums.Size(); j++)
	      if (dnums[j] > ndof)
		ndof = dnums[j];
	  }

	ndlevel.Append (ndof+1);
      }
      

    prol->Update();

    if (dirichlet_boundaries.Size())
      {
	dirichlet_dofs.SetSize (GetNDof());
	dirichlet_dofs.Clear();
	Array<int> dnums;
	for (int i = 0; i < ma.GetNSE(); i++)
	  {
	    if (dirichlet_boundaries[ma.GetSElIndex(i)])
	      {
		GetSDofNrs (i, dnums);
		for (int j = 0; j < dnums.Size(); j++)
		  if (dnums[j] != -1)
		    dirichlet_dofs.Set (dnums[j]);
	      }
	  }
      }


    if (timing) Timing();
  }

  int NodalFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }







 
  void NodalFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    ma.GetElPNums (elnr, dnums);

    if (order == 1)
      { // Ng-mesh may be second order, but FE space is 1st order
	int np = dnums.Size();
	switch (ma.GetElType(elnr))
	  {
	  case ET_TET: np = 4; break;
	  case ET_TRIG: np = 3; break;
	  case ET_QUAD: np = 4; break;
	  case ET_PRISM: np = 6; break;
          default:
            ;
	  }
	if (dnums.Size() > np) dnums.SetSize (np);
      }

    if (!DefinedOn (ma.GetElIndex (elnr)))
      dnums = -1;
  }


  void NodalFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    ma.GetSElPNums (selnr, dnums);

    if (order == 1)
      { // Ng-mesh may be second order, but FE space is 1st order
	int np = dnums.Size();
	switch (ma.GetSElType(selnr))
	  {
	  case ET_SEGM: np = 2; break;
	  case ET_TRIG: np = 3; break;
	  case ET_QUAD: np = 4; break;
          default:
            ;
	  }
	if (dnums.Size() > np) dnums.SetSize (np);
      }

    if (!DefinedOnBoundary (ma.GetSElIndex (selnr)))
      dnums = -1;
  }
  

  void NodalFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize(1);
    dnums[0] = vnr;
  }

  void NodalFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
  }

  void NodalFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
  }

  void NodalFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
  }

  
  Array<int> * 
  NodalFESpace :: CreateDirectSolverClusters (const Flags & flags) const
  {
    Array<int> & clusters = *new Array<int> (GetNDof());
    clusters = 0;

    const int stdoffset = 1;

    /*
    (*testout) << "directvertexclusters" << directvertexclusters << endl;
    (*testout) << "directedgeclusters" << directedgeclusters << endl;
    (*testout) << "directfaceclusters" << directfaceclusters << endl;
    (*testout) << "directelementclusters" << directelementclusters << endl;
    */

    int i;

    for(i=0; i<directvertexclusters.Size(); i++)
      if(directvertexclusters[i] >= 0)
	clusters[i] = directvertexclusters[i] + stdoffset;


    bool nonzero = false;
    for (i = 0; !nonzero && i < clusters.Size(); i++)
      if (clusters[i]) nonzero = true;
    if (!nonzero)
      {
	delete &clusters;
	return 0;
      }

    return &clusters;
  }








  NonconformingFESpace :: 
  NonconformingFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name="NonconformingFESpace(nonconforming)";
    // defined flags
    DefineDefineFlag("nonconforming");
    if (parseflags) CheckFlags(flags);
    
    // prol = new LinearProlongation(*this);
    

    trig = new FE_NcTrig1;

    if (ma.GetDimension() == 2)
      {
	evaluator = 
	  new MassIntegrator<2> (new ConstantCoefficientFunction(1));
	boundary_evaluator = 0;
      }
    else
      {
	evaluator = 
	  new MassIntegrator<3> (new ConstantCoefficientFunction(1));
	boundary_evaluator = 
	  new RobinIntegrator<3> (new ConstantCoefficientFunction(1));
      }

    if (dimension > 1)
      {
	evaluator = new BlockBilinearFormIntegrator (*evaluator, dimension);
	boundary_evaluator = 
	  new BlockBilinearFormIntegrator (*boundary_evaluator, dimension);  
      }
  }

  NonconformingFESpace :: ~NonconformingFESpace ()
  {
    ;
  }


  int NonconformingFESpace :: GetNDof () const
  {
    return ma.GetNEdges();
  }


  void NonconformingFESpace :: Update(LocalHeap & lh)
  {
    /*
    if (ma.GetNLevels() > ndlevel.Size())
      {
	Array<int> dnums;
	int i, j;
	int ne = ma.GetNE();
	int nse = ma.GetNSE();
	int ndof = 0;
	for (i = 0; i < ne; i++)
	  {
	    GetDofNrs (i, dnums);
	    for (j = 0; j < dnums.Size(); j++)
	      if (dnums[j] > ndof)
		ndof = dnums[j];
	  }
	for (i = 0; i < nse; i++)
	  {
	    GetSDofNrs (i, dnums);
	    for (j = 0; j < dnums.Size(); j++)
	      if (dnums[j] > ndof)
		ndof = dnums[j];
	  }

	ndlevel.Append (ndof+1);
      }

    prol->Update();

    if (dirichlet_boundaries.Size())
      {
	dirichlet_dofs.SetSize (GetNDof());
	dirichlet_dofs.Clear();
	Array<int> dnums;
	for (int i = 0; i < ma.GetNSE(); i++)
	  {
	    if (dirichlet_boundaries[ma.GetSElIndex(i)])
	      {
		GetSDofNrs (i, dnums);
		for (int j = 0; j < dnums.Size(); j++)
		  if (dnums[j] != -1)
		    dirichlet_dofs.Set (dnums[j]);
	      }
	  }
      }
    */
  }

 
  void NonconformingFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    ma.GetElEdges (elnr, dnums);
    if (!DefinedOn (ma.GetElIndex (elnr)))
      dnums = -1;
  }


  void NonconformingFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    ma.GetSElEdges (selnr, dnums);
    if (!DefinedOnBoundary (ma.GetSElIndex (selnr)))
      dnums = -1;
  }
  








  ElementFESpace :: ElementFESpace (const MeshAccess & ama, const Flags& flags, 
                                    bool parseflags)
    : FESpace (ama, flags)
  {
    name="ElementFESpace(l2)";
    if (parseflags) CheckFlags(flags);
    
    order = int(flags.GetNumFlag ("order", 0));

    prol = new ElementProlongation (*this);

    if (order == 0)
    {
      tet     = new FE_Tet0;
      prism   = new FE_Prism0;
      pyramid = new FE_Pyramid0;
      hex     = new FE_Hex0;
      trig    = new FE_Trig0;
      quad    = new FE_Quad0;
      segm    = new FE_Segm0;

      n_el_dofs = 1;
    }
    else
    {
      tet     = new FE_Tet1;
      prism   = new FE_Prism1;
      pyramid = new FE_Pyramid1;
      trig    = new FE_Trig1;
      quad    = new FE_Quad1;
      segm    = new FE_Segm1;

      if (ma.GetDimension() == 2)
        n_el_dofs = 4;
      else
        n_el_dofs = 6;
    }


    if (ma.GetDimension() == 2)
    {
      evaluator = 
          new MassIntegrator<2> (new ConstantCoefficientFunction(1));
      boundary_evaluator = 0;
    }
    else
    {
      evaluator = 
          new MassIntegrator<3> (new ConstantCoefficientFunction(1));
      boundary_evaluator = 0;
    }

    if (dimension > 1)
      evaluator = new BlockBilinearFormIntegrator (*evaluator, dimension);
  }
  
  ElementFESpace :: ~ElementFESpace ()
  {
    ;
  }

  void ElementFESpace :: Update(LocalHeap & lh)
  {
    while (ma.GetNLevels() > ndlevel.Size())
      ndlevel.Append (n_el_dofs * ma.GetNE());
  }

  
  void ElementFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    if (order == 0)
      {
	dnums.SetSize(1);
	dnums[0] = elnr;
      }
    else if (order == 1)
      {
	switch (GetMeshAccess().GetElType(elnr))
	  {
	  case ET_TRIG:
	    dnums.SetSize(3);
	    break;
	  case ET_QUAD:
	    dnums.SetSize(4);
	    break;
	  case ET_TET:
	    dnums.SetSize(4);
	    break;
	  case ET_PRISM:
	    dnums.SetSize(6);
	    break;
	  case ET_PYRAMID:
	    dnums.SetSize(5);
	    break;
	  default:
	    throw Exception ("ElementFESpace::GetDofNrs, unknown element type");
	    break;
	  }

	for (int i = 0; i < dnums.Size(); i++)
	  dnums[i] = n_el_dofs*elnr+i;
      }
  }
  void ElementFESpace :: GetSDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
  }
  int ElementFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }

  /*
  const FiniteElement & ElementFESpace :: GetSFE (int selnr) const
  {
    throw Exception ("ElementFESpace::GetSFE not available");
  }
  */


 
  SurfaceElementFESpace :: 
      SurfaceElementFESpace (const MeshAccess & ama, const Flags& flags, bool parseflags)
  : FESpace (ama, flags)
  {
    name="SurfaceElementFESpace(surfl2)";
    if(parseflags) CheckFlags(flags);
    
    // prol = new SurfaceElementProlongation (GetMeshAccess(), *this);

    if (order == 0)
    {
      trig    = new FE_Trig0;
      quad    = new FE_Quad0;
      segm    = new FE_Segm0;

      n_el_dofs = 1;
    }

    else if (order == 1)
    {
      trig    = new FE_Trig1;
      quad    = new FE_Quad1;
      segm    = new FE_Segm1;
	
      if (ma.GetDimension() == 2)
        n_el_dofs = 2;
      else
        n_el_dofs = 4;
    }

    else if (order == 2)
    {
      trig    = new FE_Trig2HB;
      quad    = new FE_Quad1;
      segm    = new FE_Segm2;

      if (ma.GetDimension() == 2)
        n_el_dofs = 3;
      else
        n_el_dofs = 9;
    }

    boundary_evaluator = 
        new RobinIntegrator<3> (new ConstantCoefficientFunction(1));

    if (dimension > 1)
      boundary_evaluator = new BlockBilinearFormIntegrator (*boundary_evaluator, dimension);
  }

  
  SurfaceElementFESpace :: ~SurfaceElementFESpace ()
  {
    ;
  }

  void SurfaceElementFESpace :: Update(LocalHeap & lh)
  {
    const MeshAccess & ma = GetMeshAccess();

    while (ma.GetNLevels() > ndlevel.Size())
      ndlevel.Append (n_el_dofs * ma.GetNSE());
  }

  const FiniteElement & SurfaceElementFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    throw Exception ("SurfaceElementFESpace::GetFE not available");
  }

  
  void SurfaceElementFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize (0);
  }

  int SurfaceElementFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }


  void SurfaceElementFESpace :: GetSDofNrs (int elnr, Array<int> & dnums) const
  {
    if (order == 0)
      {
	dnums.SetSize(1);
	dnums[0] = elnr;
      }
    else if (order == 1)
      {
	switch (GetMeshAccess().GetSElType(elnr))
	  {
	  case ET_SEGM:
	    dnums.SetSize(2);
	    break;
	  case ET_TRIG:
	    dnums.SetSize(3);
	    break;
	  case ET_QUAD:
	    dnums.SetSize(4);
	    break;
	  default:
	    dnums.SetSize(4);
	    break;
	  }
	for (int i = 0; i < dnums.Size(); i++)
	  dnums[i] = n_el_dofs*elnr+i;
      }
    else if (order == 2)
      {
	switch (GetMeshAccess().GetSElType(elnr))
	  {
	  case ET_SEGM:
	    dnums.SetSize(3);
	    break;
	  case ET_TRIG:
	    dnums.SetSize(6);
	    break;
	  case ET_QUAD:
	    dnums.SetSize(4);
	    break;
	  default:
	    dnums.SetSize(4);
	    break;
	  }
	for (int i = 0; i < dnums.Size(); i++)
	  dnums[i] = n_el_dofs*elnr+i;
      }
  }





 











  CompoundFESpace :: CompoundFESpace (const MeshAccess & ama,
				      const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name="CompoundFESpaces";
    DefineDefineFlag("compound");
    DefineStringListFlag("spaces");
    if(parseflags) CheckFlags(flags);
    
    Array<const Prolongation*> prols;
    prol = new CompoundProlongation (this, prols);
  }







  CompoundFESpace :: CompoundFESpace (const MeshAccess & ama,
				      const Array<FESpace*> & aspaces,
				      const Flags & flags, bool parseflags)
    : FESpace (ama, flags), spaces(aspaces)
  {
    name="CompoundFESpaces";
    DefineDefineFlag("compound");
    DefineStringListFlag("spaces");
    if(parseflags) CheckFlags(flags);
    
    Array<const Prolongation*> prols(spaces.Size());
    for (int i = 0; i < spaces.Size(); i++)
      prols[i] = spaces[i]->GetProlongation();

    prol = new CompoundProlongation (this, prols);
  }


  void CompoundFESpace :: AddSpace (FESpace * fes)
  {
    spaces.Append (fes);
    dynamic_cast<CompoundProlongation*> (prol) -> AddProlongation (fes->GetProlongation());
  }

  CompoundFESpace :: ~CompoundFESpace ()
  {
    ;
  }

  void CompoundFESpace :: Update(LocalHeap & lh)
  {
    FESpace :: Update (lh);

    cummulative_nd.SetSize (spaces.Size()+1);
    cummulative_nd[0] = 0;
    for (int i = 0; i < spaces.Size(); i++)
      {
	const_cast<FESpace*> (spaces[i])->Update(lh);
	cummulative_nd[i+1] = cummulative_nd[i] + spaces[i]->GetNDof();
      }

    while (ma.GetNLevels() > ndlevel.Size())
      ndlevel.Append (cummulative_nd.Last());

    prol -> Update();


    UpdateCouplingDofArray();
    // FinalizeUpdate (lh);

    // dirichlet-dofs from sub-spaces
    // ist das umsonst ? (JS)
    bool has_dirichlet_dofs = false;
    for (int i = 0; i < spaces.Size(); i++)
      if (spaces[i]->GetFreeDofs()) 
	has_dirichlet_dofs = true;

    if (has_dirichlet_dofs)
      {
	free_dofs.SetSize(GetNDof());
	free_dofs.Set();

	for (int i = 0; i < spaces.Size(); i++)
	  {
	    const BitArray * free_dofs_sub = spaces[i]->GetFreeDofs();
	    if (free_dofs_sub)
	      {
		int base = cummulative_nd[i];
		int nd = cummulative_nd[i+1] - base;
		for (int i = 0; i < nd; i++)
		  if (!free_dofs_sub->Test(i))
		    free_dofs.Clear (base+i);
	      }
	  }

	dirichlet_dofs = free_dofs;
	dirichlet_dofs.Invert();
      }


    // (*testout) << "free_dofs = " << endl << free_dofs << endl;

    (*testout) << "Update compound fespace" << endl;
    (*testout) << "cummulative dofs start at " << cummulative_nd << endl;
    (*testout) << "dof coupling types " << ctofdof << endl;
  }

  void CompoundFESpace :: FinalizeUpdate(LocalHeap & lh)
  {
    for (int i = 0; i < spaces.Size(); i++)
      spaces[i] -> FinalizeUpdate(lh);
    FESpace::FinalizeUpdate (lh);
  }



  void CompoundFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(this->GetNDof());

    for (int i = 0; i < spaces.Size(); i++)
      for (int j=0; j< spaces[i]->GetNDof();j++)
	ctofdof[cummulative_nd[i]+j] = spaces[i]->GetDofCouplingType(j);	

    *testout << "CompoundFESpace :: UpdateCouplingDofArray() presents \n ctofdof = \n" << ctofdof << endl;
  }



  const FiniteElement & CompoundFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    ArrayMem<const FiniteElement*, 10> fea(spaces.Size());
    for (int i = 0; i < fea.Size(); i++)
      fea[i] = &spaces[i]->GetFE(elnr, lh);
    return *new (lh) CompoundFiniteElement (fea);
  }

  
  void CompoundFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    ArrayMem<int,500> hdnums;
    dnums.SetSize(0);
    for (int i = 0; i < spaces.Size(); i++)
      {
	spaces[i]->GetDofNrs (elnr, hdnums);
	for (int j = 0; j < hdnums.Size(); j++)
	  if (hdnums[j] != -1)
	    dnums.Append (hdnums[j]+cummulative_nd[i]);
	  else
	    dnums.Append (-1);
      }
  }


  void CompoundFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    ArrayMem<int,500> hdnums;
    dnums.SetSize(0);
    for (int i = 0; i < spaces.Size(); i++)
      {
	spaces[i]->GetVertexDofNrs (vnr, hdnums);
	for (int j = 0; j < hdnums.Size(); j++)
	  if (hdnums[j] != -1)
	    dnums.Append (hdnums[j]+cummulative_nd[i]);
	  else
	    dnums.Append (-1);
      }
  }

  void CompoundFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    ArrayMem<int,500> hdnums;
    dnums.SetSize(0);
    for (int i = 0; i < spaces.Size(); i++)
      {
	spaces[i]->GetEdgeDofNrs (ednr, hdnums);
	for (int j = 0; j < hdnums.Size(); j++)
	  if (hdnums[j] != -1)
	    dnums.Append (hdnums[j]+cummulative_nd[i]);
	  else
	    dnums.Append (-1);
      }

  }
  void CompoundFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    ArrayMem<int,500> hdnums;
    dnums.SetSize(0);
    for (int i = 0; i < spaces.Size(); i++)
      {
	spaces[i]->GetFaceDofNrs (fanr, hdnums);
	for (int j = 0; j < hdnums.Size(); j++)
	  if (hdnums[j] != -1)
	    dnums.Append (hdnums[j]+cummulative_nd[i]);
	  else
	    dnums.Append (-1);
      }
  }

  void CompoundFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    ArrayMem<int,500> hdnums;
    dnums.SetSize(0);

    for (int i = 0; i < spaces.Size(); i++)
      {
	spaces[i]->GetInnerDofNrs (elnr, hdnums);
	for (int j = 0; j < hdnums.Size(); j++)
	  if (hdnums[j] != -1)
	    dnums.Append (hdnums[j]+cummulative_nd[i]);
	  else
	    dnums.Append (-1);
      }
  }
  










  const FiniteElement & CompoundFESpace :: GetSFE (int elnr, LocalHeap & lh) const
  {
    ArrayMem<const FiniteElement*, 10> fea(spaces.Size());
    for (int i = 0; i < fea.Size(); i++)
      fea[i] = &spaces[i]->GetSFE(elnr, lh);
    return *new (lh) CompoundFiniteElement (fea);
  }


  void CompoundFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    ArrayMem<int,500> hdnums;
    dnums.SetSize(0);
    for (int i = 0; i < spaces.Size(); i++)
      {
	spaces[i]->GetSDofNrs (selnr, hdnums);
	for (int j = 0; j < hdnums.Size(); j++)
	  if (hdnums[j] != -1)
	    dnums.Append (hdnums[j]+cummulative_nd[i]);
	  else
	    dnums.Append (-1);
      }
  }




  template <class MAT>
  void CompoundFESpace::TransformMat (int elnr, bool boundary,
				      MAT & mat, TRANSFORM_TYPE tt) const
  {
    int base = 0;
    LocalHeapMem<100005> lh("CompoundFESpace - transformmat");
    for (int i = 0; i < spaces.Size(); i++)
      {
	int nd = boundary  ?
	  spaces[i]->GetSFE(elnr, lh).GetNDof()
	  : spaces[i]->GetFE(elnr, lh).GetNDof();

	lh.CleanUp();

	spaces[i]->TransformMat (elnr, boundary, mat.Rows(base, base+nd), TRANSFORM_MAT_LEFT);
	spaces[i]->TransformMat (elnr, boundary, mat.Cols(base, base+nd), TRANSFORM_MAT_RIGHT);

	base += nd;
      }
  }

  template <class VEC>
  void CompoundFESpace::TransformVec (int elnr, bool boundary,
				      VEC & vec, TRANSFORM_TYPE tt) const
  {
    LocalHeapMem<100006> lh("CompoundFESpace - transformvec");
    for (int i = 0, base = 0; i < spaces.Size(); i++)
      {
	int nd = boundary ? 
	  spaces[i]->GetSFE(elnr, lh).GetNDof() :
	  spaces[i]->GetFE(elnr, lh).GetNDof();
	
	lh.CleanUp();

	// VEC svec (nd, &vec(base));
	spaces[i]->TransformVec (elnr, boundary, vec.Range(base, base+nd), tt);
	base += nd;
      }
  }



template
void CompoundFESpace::TransformVec<FlatVector<double> >
(int elnr, bool boundary, FlatVector<double> & vec, TRANSFORM_TYPE tt) const;
template
void CompoundFESpace::TransformVec<FlatVector<Complex> >
(int elnr, bool boundary, FlatVector<Complex> & vec, TRANSFORM_TYPE tt) const;

template
void CompoundFESpace::TransformMat<FlatMatrix<double> > 
(int elnr, bool boundary, FlatMatrix<double> & mat, TRANSFORM_TYPE tt) const;
template
void CompoundFESpace::TransformMat<FlatMatrix<Complex> > 
(int elnr, bool boundary, FlatMatrix<Complex> & mat, TRANSFORM_TYPE tt) const;

template
void CompoundFESpace::TransformMat<SliceMatrix<double> > 
(int elnr, bool boundary, SliceMatrix<double> & mat, TRANSFORM_TYPE tt) const;
template
void CompoundFESpace::TransformMat<SliceMatrix<Complex> > 
(int elnr, bool boundary, SliceMatrix<Complex> & mat, TRANSFORM_TYPE tt) const;


  void CompoundFESpace::VTransformMR (int elnr, bool boundary,
				      const SliceMatrix<double> & mat, TRANSFORM_TYPE tt) const 
  {
    TransformMat (elnr, boundary, mat, tt);
  }
  
  void CompoundFESpace::VTransformMC (int elnr, bool boundary,
				      const SliceMatrix<Complex> & mat, TRANSFORM_TYPE tt) const
  {
    TransformMat (elnr, boundary, mat, tt);
  }
  

  void CompoundFESpace::VTransformVR (int elnr, bool boundary,
				      const FlatVector<double> & vec, TRANSFORM_TYPE tt) const 
  {
    TransformVec (elnr, boundary, vec, tt);
  }
  
  void CompoundFESpace::VTransformVC (int elnr, bool boundary,
				      const FlatVector<Complex> & vec, TRANSFORM_TYPE tt) const 
  {
    TransformVec (elnr, boundary, vec, tt);
  }












  FESpaceClasses::FESpaceInfo::
  FESpaceInfo (const string & aname,
	       FESpace* (*acreator)(const MeshAccess & ma, const Flags & flags))
    : name(aname), creator(acreator)
  {
    ;
  }
  
  FESpaceClasses :: FESpaceClasses ()
  {
    ;
  }

  FESpaceClasses :: ~FESpaceClasses()
  {
    for (int i = 0; i < fesa.Size(); i++)
      delete fesa[i];
  }
  
  void FESpaceClasses :: 
  AddFESpace (const string & aname,
	      FESpace* (*acreator)(const MeshAccess & ma, const Flags & flags))
  {
    fesa.Append (new FESpaceInfo(aname, acreator));
  }

  const FESpaceClasses::FESpaceInfo * 
  FESpaceClasses::GetFESpace(const string & name)
  {
    for (int i = 0; i < fesa.Size(); i++)
      {
	if (name == fesa[i]->name)
	  return fesa[i];
      }
    return 0;
  }

  void FESpaceClasses :: Print (ostream & ost) const
  {
    ost << endl << "FESpaces:" << endl;
    ost <<         "---------" << endl;
    ost << setw(20) << "Name" << endl;
    for (int i = 0; i < fesa.Size(); i++)
      ost << setw(20) << fesa[i]->name << endl;
  }

 
  FESpaceClasses & GetFESpaceClasses ()
  {
    static FESpaceClasses fecl;
    return fecl;
  }


  // standard fespaces:

  static RegisterFESpace<NodalFESpace> initnodalfes ("nodal");
  static RegisterFESpace<NonconformingFESpace> initncfes ("nonconforming");

}



