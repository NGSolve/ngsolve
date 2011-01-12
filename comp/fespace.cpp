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

#include <parallelngs.hpp>

#include "../fem/h1lofe.hpp"


#ifdef PARALLEL
extern MPI_Group MPI_HIGHORDER_WORLD;
extern MPI_Comm MPI_HIGHORDER_COMM;
#endif

using namespace ngmg;

namespace ngcomp
{
  using namespace ngcomp;
  using namespace ngparallel;




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
      cout << "ATTENTION: flag dgjumps is used!\n This leads to a \
lot of new non-zero entries in the matrix!\n" << endl;
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
	    else
	      cerr << "Illegal Dirichlet boundary index " << bnd+1 << endl;
	  }
        cout << "dirichlet_boundaries" << dirichlet_boundaries << endl;
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

    tet = 0;
    pyramid = 0;
    prism = 0;
    hex = 0;
    trig = 0;
    quad = 0;
    segm = 0;

    evaluator = 0;
    boundary_evaluator = 0;

    first_lodof[4] = -1;   // indicates that nodes are not used

    element_coloring = NULL;

#ifdef PARALLEL
    paralleldofs = 0;
#endif
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

 #ifdef PARALLEL
    delete paralleldofs;
 #endif
  }
  

  void FESpace :: Update(LocalHeap & lh)
  {
    for (int i=0; i< specialelements.Size(); i++)
      delete specialelements[i]; 
    specialelements.SetSize(0);
    if (dirichlet_boundaries.Size())
      {
	int dim = ma.GetDimension();

	dirichlet_vertex.SetSize (ma.GetNV());
	dirichlet_edge.SetSize (ma.GetNEdges());
        if (dim == 3)
          dirichlet_face.SetSize (ma.GetNFaces());
	
	dirichlet_vertex = false;
	dirichlet_edge = false;
	dirichlet_face = false;

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

	(*testout) << "Dirichlet_vertex = " << endl << dirichlet_vertex << endl;
	(*testout) << "Dirichlet_edge = " << endl << dirichlet_edge << endl;
	(*testout) << "Dirichlet_face = " << endl << dirichlet_face << endl;
      }
      else 
	(*testout) << "No Dirichlet boundaries!" << endl;
            
  }

  void FESpace :: FinalizeUpdate(LocalHeap & lh)
  {

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


	free_dofs.SetSize (GetNDof());
	free_dofs = dirichlet_dofs;
	free_dofs.Invert();
      }
    



    *testout << "coloring ... " << flush;
    Array<int> col(ma.GetNE());
    Array<int> dofcol(GetNDof());
    
    dofcol = -1;
    col = -1;
    bool found;
    Array<int> dnums;
    int color = 0;
    do
      {
	// cout << "color = " << color << endl;
	found = false;
	for (int i = 0; i < ma.GetNE(); i++)
	  {
	    if (col[i] >= 0) continue;
	    GetDofNrs (i, dnums);
	    bool ok = true;
	    for (int j = 0; j < dnums.Size(); j++)
	      if (dnums[j] != -1)
		if (dofcol[dnums[j]] == color)
		  {
		    ok = false;
		    break;
		  }
	    if (ok)
	      {
		col[i] = color;
		for  (int j = 0; j < dnums.Size(); j++)
		  if (dnums[j] != -1)
		    dofcol[dnums[j]] = color;
		found = true;
	      }
	    // cout << "element " << i << " dnums = " << dnums << " ok = " << ok << endl;
	  }
	color++;
      }
    while (found);

    // cout << "col = " << endl << col << endl;
    color--;

    Array<int> cntcol(color);
    cntcol = 0;
    for (int i = 0; i < ma.GetNE(); i++)
      cntcol[col[i]]++;

    delete element_coloring;
    element_coloring = new Table<int> (cntcol);
    cntcol = 0;
    for (int i = 0; i < ma.GetNE(); i++)
      (*element_coloring)[col[i]][cntcol[col[i]]++] = i;
    
    *testout << "needed " << color << " colors" << endl;
    // cout << "elcolors = " << endl << (*element_coloring) << endl;
    
  }

  const FiniteElement & FESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    FiniteElement * fe = 0;
    
    switch (ma.GetElType(elnr))
      {
      case ET_TET:
	fe = tet; break;
      case ET_PYRAMID:
	fe = pyramid; break;
      case ET_PRISM:
	fe = prism; break;
      case ET_HEX:
	fe = hex; break;
      case ET_TRIG:
	fe = trig; break;
      case ET_QUAD:
	fe = quad; break;
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
    if (first_lodof[4] != -1)
      {
	dnums.SetSize(0);
	for (int i = 0; i < lodofs_per_node[nt]; i++)
	  dnums.Append (first_lodof[nt] + i);
	if (first_hodofs[nt].Size())
	  for (int i = first_hodofs[nt][nr]; i < first_hodofs[nt][nr+1]; i++)
	    dnums.Append (i);
	return;
      }


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
    if (first_lodof[4] != -1)
      GetNodeDofNrs (NT_VERTEX, vnr, dnums);
    throw Exception ("FESpace::GetVertexDofNrs called");
  }

  void FESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    if (first_lodof[4] != -1)
      GetNodeDofNrs (NT_EDGE, ednr, dnums);
    throw Exception ("FESpace::GetEdgeDofNrs called");
  }

  void FESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    if (first_lodof[4] != -1)
      GetNodeDofNrs (NT_FACE, fanr, dnums);
    throw Exception ("FESpace::GetFaceDofNrs called");
  }

  void FESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    if (first_lodof[4] != -1)
      GetNodeDofNrs (NT_CELL, elnr, dnums);
    throw Exception (string("FESpace::GetInnerDofNrs called for class")+
		     typeid(*this).name());
  }




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
    /*
    int nd = GetNDof();
    
    Array<int> cnt(nd);
    for (int i = 0; i < nd; i++)
      cnt[i] = IsDirichletDof(i) ? 0 : 1;

    Table<int> * it = new Table<int>(cnt);

    for (int i = 0; i < nd; i++)
      if (!IsDirichletDof(i))
	(*it)[i][0] = i;

    return it;
    */

    
    int nd = GetNDof();
    TableCreator<int> creator;

    /*
    for (int mode = 1; mode <= 3; mode++)
      {
	creator.SetMode (mode);
    */
    for ( ; !creator.Done(); creator++)
      {
	for (int i = 0; i < nd; i++)
	  if (!IsDirichletDof(i))
	    creator.Add (i, i);
      }
    return creator.GetTable();
  }

    
  /*
  Table<int> * FESpace :: CreateSmoothingBlocks (const Flags & flags) const
  {
    SmoothingBlocksCreator sbc;
    sbc.SetMode (1);
    CreateSmoothingBlocks2 (sbc, flags);
    sbc.SetMode (2);
    CreateSmoothingBlocks2 (sbc, flags);
    sbc.SetMode (3);
    CreateSmoothingBlocks2 (sbc, flags);

    return sbc.GetTable();
  }

  void FESpace :: CreateSmoothingBlocks2 (SmoothingBlocksCreator & sbc, const Flags & flags) const
  {
  }
  */


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


#ifdef PARALLEL
  void FESpace :: UpdateParallelDofs ( )
  {
    if  ( ntasks == 1 ) 
      {
	if ( paralleldofs ) delete paralleldofs;
	paralleldofs = 0;
	return;
      }

    if ( id == 0 )
      {
	if ( paralleldofs == 0 ) 
	  paralleldofs = new ParallelDofs (*this);
	paralleldofs -> Update();

	(*this).UpdateParallelDofs_loproc();

	this->paralleldofs->UpdateMPIType();
	return;
      }

    // no updating of parallel dofs for the low-order processor 0
    // or in case of only one processor

    if ( paralleldofs == 0 ) 
      paralleldofs = new ParallelDofs (*this);

    paralleldofs -> Update();
    (*this). UpdateParallelDofs_hoproc();
    this->paralleldofs->UpdateMPIType();
  }

  void FESpace :: UpdateParallelDofs ( LocalHeap & lh )
  {
    cout << "UPDATE PARALLELDOF (LH) NOT IMPLEMENTED!!!!" << endl;
  }



  void FESpace :: ResetParallelDofs ()
  {
    paralleldofs = 0;
  }

  void FESpace :: UpdateParallelDofs_loproc()
  {
    cout << "ERROR -- FESpace::UpdateParallelDofs_loproc() called!" << endl;
  }

  void FESpace :: UpdateParallelDofs_hoproc()
  {
    // ******************************
    // update exchange dofs 
    // ******************************
    *testout << "FESpace::UpdateParallelDofs_hoproc" << endl;
    // Find number of exchange dofs
    Array<int> nexdof(ntasks);
    nexdof = 0;

    Array<MPI_Request> sendrequest(ntasks);
    Array<MPI_Request> recvrequest(ntasks);
    MPI_Status status;

    Array<int> dnums;


    for (NODE_TYPE nt = NT_VERTEX; nt <= NT_CELL; nt++)
      {
	if ( eliminate_internal && nt == NT_CELL ) break;
	for ( int nr = 0; nr < ma.GetNNodes(nt); nr++ )
	  {
	    if ( !parallelma->IsExchangeNode ( nt, nr ) ) continue;
	    
	    GetNodeDofNrs ( nt, nr, dnums );
	    nexdof[id] += dnums.Size();
	    
	    for ( int dest = 1; dest < ntasks; dest ++ )
	      if (  parallelma -> GetDistantNodeNum ( dest, nt, nr ) >= 0 )
		nexdof[dest] += dnums.Size(); 
	  }
      }

    nexdof[0] = LowOrderFESpace() . GetNDof();

    paralleldofs->SetNExDof(nexdof);

    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);

    Array<int> ** owndofs, ** distantdofs;
    owndofs = new Array<int>* [ntasks];
    distantdofs = new Array<int>* [ntasks];

    for ( int i = 0; i < ntasks; i++ )
      {
	owndofs[i] = new Array<int>(1);
	(*owndofs[i])[0] = GetNDof();
	distantdofs[i] = new Array<int>(0);
      }



    Array<int> cnt_nexdof(ntasks);
    cnt_nexdof = 0;
    int exdof = 0;
    int ii = 1;

    // *******************
    // Parallel Node dofs
    // *******************



    for (NODE_TYPE nt = NT_VERTEX; nt <= NT_CELL; nt++)
      {
	if ( eliminate_internal && nt == NT_CELL ) break;

	for ( int dest = 0; dest < ntasks; dest++)
	  {
	    owndofs[dest]->SetSize(1);
	    distantdofs[dest]->SetSize(0);
	  }

	for ( int node = 0; node < ma.GetNNodes(nt); node++ )
	  if ( parallelma->IsExchangeNode(nt, node) )
	    {
	      GetNodeDofNrs (nt, node, dnums); 
	      if ( dnums.Size() == 0 ) continue;
	      
	      for ( int i=0; i<dnums.Size(); i++ )
		(*(paralleldofs->sorted_exchangedof))[id][exdof++] = dnums[i];

	      Array<int[2]> distantnodenums;
	      parallelma -> GetDistantNodeNums ( nt, node, distantnodenums);
	      for ( int idest = 1; idest < distantnodenums.Size(); idest++ )
		{
		  int dest = distantnodenums[idest][0];
		  int distnode = distantnodenums[idest][1];

		  owndofs[dest]->Append ( distnode );
		  owndofs[dest]->Append (int(paralleldofs->IsGhostDof(dnums[0])) );

		  for ( int i=0; i<dnums.Size(); i++)
		    {
		      paralleldofs->SetExchangeDof ( dest, dnums[i] );
		      paralleldofs->SetExchangeDof ( dnums[i] );
		      owndofs[dest]->Append ( dnums[i] );
		    }
		}
	    } 


	for ( int dest = 1; dest < ntasks; dest ++ )
	  if ( dest != id )
	    {
	      MyMPI_ISend ( *owndofs[dest], dest, sendrequest[dest]);
	    }
	
	for ( int sender = 1; sender < ntasks; sender ++ )
	  if ( sender != id )
	    {
	      MyMPI_IRecv ( *distantdofs[sender], sender, recvrequest[sender]);
	    }
	
	for( int dest=1; dest<ntasks; dest++) 
	  {
	    ii = 1;
	    if ( dest == id ) continue;
	    MPI_Wait ( &recvrequest[dest], &status );

	    // low order dofs first
	    int n_lo = GetNLowOrderNodeDofs (nt);
	    if ( n_lo )
	      while ( ii < distantdofs[dest]->Size() )
		{
		  int nodenum = (*distantdofs[dest])[ii++];
		  int isdistghost = (*distantdofs[dest])[ii++];
		  Array<int> dnums, lodnums;
		  GetNodeDofNrs (nt, nodenum, dnums);
		  for ( int i = 0; i < n_lo; i++ )
		    {
		      (*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
		      if ( dest < id && !isdistghost )
			paralleldofs->ismasterdof.Clear ( dnums[i] ) ;
		      ii++;
		      cnt_nexdof[dest]++;
		    }
		  ii += dnums.Size() - n_lo;
		}
	    // then the high order dofs
	    ii = 1;
	    while ( ii < distantdofs[dest]->Size() )
	      {
		int nodenum = (*distantdofs[dest])[ii++];
		int isdistghost = (*distantdofs[dest])[ii++];
		Array<int> dnums;
		GetNodeDofNrs (nt, nodenum, dnums);

		ii += n_lo;
		for ( int i=n_lo; i<dnums.Size(); i++)
		  {
		    (*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
		    if ( dest < id && !isdistghost )
		      paralleldofs->ismasterdof.Clear ( dnums[i] ) ;
		    ii++; cnt_nexdof[dest]++;
		  }
	      }
	    
	    
	  }

	for ( int dest = 1; dest < ntasks; dest ++ )
	  if ( dest != id )
	    MPI_Wait ( &sendrequest[dest], &status );
      }



    for ( int dest = id+1; dest < ntasks; dest++ )
      QuickSort ( (*(paralleldofs->sorted_exchangedof))[dest] );









    /*******************************

         update low order space

    *****************************/


    int ndof_lo = low_order_space->GetNDof();

    // all dofs are exchange dofs
    nexdof = ndof_lo;
 
    exdof = 0;
    cnt_nexdof = 0;


    // *****************
    // Parallel Node dofs
    // *****************

    for (NODE_TYPE nt = NT_VERTEX; nt <= NT_CELL; nt++)
      {
	if ( eliminate_internal && nt == NT_CELL ) break;

	int n_lo = GetNLowOrderNodeDofs(nt);
	if ( !n_lo ) continue;

	int dest = 0;

	owndofs[0]->SetSize(1);
	(*owndofs[0])[0] = GetNDof();
	distantdofs[0]->SetSize(0);
    
	// find local and distant dof numbers for vertex exchange dofs
	for ( int node = 0; node < ma.GetNNodes(nt); node++ )
	  {
	    GetNodeDofNrs (nt, node, dnums);

	    Array<int[2]> distantnodenums;
	    parallelma -> GetDistantNodeNums ( nt, node, distantnodenums);
	    owndofs[0]->Append (distantnodenums[0][1] );

	    for ( int i = 0; i < n_lo; i++ )
	      {
		paralleldofs->SetExchangeDof ( dest, dnums[i]  );
		owndofs[0]->Append ( dnums[i] );
	      }
	  }
	
	
	MyMPI_ISend ( *owndofs[0], dest, sendrequest[dest]);
	MyMPI_IRecv ( *distantdofs[0], dest, recvrequest[dest] );
	
	MPI_Wait ( &recvrequest[dest], &status );
	
	ii = 1;
	while ( ii < distantdofs[0]->Size() )
	  {
	    if ( dest == id ) continue;
	    paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
	    int nodenum = (*distantdofs[dest])[ii++];
	    Array<int> dnums, lodnums;
	    GetNodeDofNrs (nt, nodenum, dnums);
	    
	    for ( int i = 0; i < n_lo; i++ )
	      {
		(*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
		ii++; cnt_nexdof[dest]++;
	      }
	  }


      }
    for ( int i = 0; i < ntasks; i++ )
      delete distantdofs[i], owndofs[i];

    delete [] owndofs;
    delete [] distantdofs;
  }
#endif






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
	if (dgjumps){ *testout << "(NodalFES:)setting loflag dgjumps " << endl; loflags.SetFlag ("dgjumps");}
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


#ifdef PARALLEL
    paralleldofs = 0;
#endif

  }

  NodalFESpace :: ~NodalFESpace ()
  {
    ;
  }

  FESpace * NodalFESpace :: 
  Create (const MeshAccess & ma, const Flags & flags)
  {
    return new NodalFESpace (ma, flags, true);
  }


  int NodalFESpace :: GetNDof () const
  {
    return ndlevel.Last();
  }

  void NodalFESpace :: Update(LocalHeap & lh)
  {
    if (low_order_space)
	low_order_space -> Update(lh);

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

    FinalizeUpdate (lh);


    if (timing) Timing();

#ifdef PARALLEL
    
    try
      {
	UpdateParallelDofs();
      }
    catch (exception & e)
      {
	throw Exception (e.what() + 
			 string ("\nthrown by NodalFESpace, UpdateParallelDofs") );
      }
    catch (Exception & e)
      {
	e.Append (string ("\nthrown by NodalFESpace, UpdateParallelDofs, id = "));
	throw;
      }

#endif
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

    /*
    if (dirichlet_dofs.Size())
      for (int j = 0; j < dnums.Size(); j++)
	if (dnums[j] != -1 && dirichlet_dofs[dnums[j]])
	  dnums[j] = -1;
    */
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

    /*
    if (dirichlet_dofs.Size())
      for (int j = 0; j < dnums.Size(); j++)
	if (dnums[j] != -1 && dirichlet_dofs[dnums[j]])
	  dnums[j] = -1;
    */
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




#ifdef PARALLEL
  void NodalFESpace :: UpdateParallelDofs_loproc()
  {
    *testout << "Nodal::UpdateParallelDofs_loproc" << endl;

    const MeshAccess & ma = (*this). GetMeshAccess();

    int ndof = GetNDof();

    // Find number of exchange dofs
    Array<int> nexdof(ntasks); 
    nexdof = 0;


    // number of vertex exchange dofs
    for ( int vert = 0; vert < ma.GetNV(); vert++ )
      {
	nexdof[id] ++;//= dnums.Size() ; 
	for ( int dest = 1; dest < ntasks; dest ++ )
	  if (  parallelma -> GetDistantPNum ( dest, vert ) >= 0 )
	    { 
	      nexdof[dest] ++; 
	    }
      }

    paralleldofs->SetNExDof(nexdof);

//     paralleldofs->localexchangedof = new Table<int> (nexdof);
//     paralleldofs->distantexchangedof = new Table<int> (nexdof);
    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);

    MPI_Status status;
    MPI_Request * sendrequest = new MPI_Request [ntasks];
    MPI_Request * recvrequest = new MPI_Request [ntasks];

    Array<int> ** owndofs,** distantdofs;
    owndofs = new Array<int> * [ntasks];
    distantdofs = new Array<int> * [ntasks];

    for ( int i = 0; i < ntasks; i++)
      {
	owndofs[i] = new Array<int> (1);
	(*owndofs[i])[0] = ndof;
	distantdofs[i] = new Array<int> (0);
      }

    int exdof = 0;
    Array<int> cnt_nexdof(ntasks);
    cnt_nexdof = 0;

    // *****************
    // Parallel Vertex dofs
    // *****************


    // find local and distant dof numbers for vertex exchange dofs
    for ( int vert = 0; vert < ma.GetNV(); vert++ )
      {
	(*(paralleldofs->sorted_exchangedof))[id][exdof++] = vert;

	for ( int dest = 1; dest < ntasks; dest++ )
	  {
	    int distvert = parallelma -> GetDistantPNum ( dest, vert );
	    if( distvert < 0 ) continue;

	    owndofs[dest]->Append ( distvert );
	    paralleldofs->SetExchangeDof ( dest, vert );
	    paralleldofs->SetExchangeDof ( vert );
	    owndofs[dest]->Append ( vert );
  	  }
      }   


    for ( int dest = 1; dest < ntasks; dest ++ )
      {
	MyMPI_ISend ( *owndofs[dest], dest, sendrequest[dest]);
	MyMPI_IRecv ( *distantdofs[dest], dest, recvrequest[dest]);
      }
    
    int ii = 1;
    for( int dest=1; dest<ntasks; dest++) 
      {
	if ( dest == id ) continue;
	MPI_Wait ( recvrequest + dest, &status );
	paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int vnum = (*distantdofs[dest])[ii++];
// 	    (*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = vnum;
// 	    (*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];
	    (*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = vnum;
		ii++; cnt_nexdof[dest]++;
	  }
      }

    for ( int dest = id+1; dest < ntasks; dest++ )
      QuickSort ( (*(paralleldofs->sorted_exchangedof))[dest] );

     for ( int i = 0; i < ntasks; i++ )
       delete distantdofs[i], owndofs[i];

     delete [] owndofs, distantdofs;
     delete [] sendrequest, recvrequest;
  }

  void NodalFESpace :: UpdateParallelDofs_hoproc ()
  {
    *testout << "Nodal::UpdateParallelDofs_hoproc" << endl;

    const MeshAccess & ma = (*this). GetMeshAccess();

    int ndof = GetNDof();

    // Find number of exchange dofs
    Array<int> nexdof(ntasks); 
    nexdof = 0;

    MPI_Status status;
    MPI_Request * sendrequest = new MPI_Request [ntasks];
    MPI_Request * recvrequest = new MPI_Request [ntasks];

    // number of vertex exchange dofs
    for ( int vert = 0; vert < ma.GetNV(); vert++ )
      {
	if ( !parallelma->IsExchangeVert ( vert ) ) continue;

	nexdof[id] ++;
	for ( int dest = 1; dest < ntasks; dest ++ )
	  if (  parallelma -> GetDistantPNum ( dest, vert ) >= 0 )
	    { 
	      nexdof[dest] ++; 
	    }
      }


    nexdof[0] = ndof;

    paralleldofs->SetNExDof(nexdof);

//     paralleldofs->localexchangedof = new Table<int> (nexdof);
//     paralleldofs->distantexchangedof = new Table<int> (nexdof);
    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);
 

    Array<int> ** owndofs,** distantdofs;
    owndofs = new Array<int> * [ntasks];
    distantdofs = new Array<int> * [ntasks];

    for ( int i = 0; i < ntasks; i++)
      {
	owndofs[i] = new Array<int> (1);
	(*owndofs[i])[0] = ndof;
	distantdofs[i] = new Array<int> (0);
      }

    int exdof = 0;
    Array<int> cnt_nexdof(ntasks);
    cnt_nexdof = 0;

    // *****************
    // Parallel Vertex dofs
    // *****************

    // find local and distant dof numbers for vertex exchange dofs
    for ( int vert = 0; vert < ma.GetNV(); vert++ )
      {
	if ( !parallelma->IsExchangeVert ( vert ) ) continue;
	Array<int> dnums;

	(*(paralleldofs->sorted_exchangedof))[id][exdof++] = vert;

	for ( int dest = 1; dest < ntasks; dest++ )
	  {
	    int distvert = parallelma -> GetDistantPNum ( dest, vert );
	    if( distvert < 0 ) continue;

	    owndofs[dest]->Append ( distvert );
	    owndofs[dest]->Append( int (paralleldofs->IsGhostDof(vert)) );
	    paralleldofs->SetExchangeDof ( dest, vert );
	    paralleldofs->SetExchangeDof ( vert );
	    owndofs[dest]->Append ( vert );
  	  }
      }   



     for ( int sender = 1; sender < ntasks; sender ++ )
      {
       if ( id == sender )
          for ( int dest = 1; dest < ntasks; dest ++ )
            if ( dest != id )
              {
		MyMPI_ISend ( (*owndofs[dest]), dest, sendrequest[dest] );
              }
	  
        if ( id != sender )
          {
	    int n2;
	    MyMPI_IRecv ( *distantdofs[sender], sender, recvrequest[sender]);
          }
 	  
	 
      }

    int ii = 1;
    for( int dest=1; dest<ntasks; dest++) 
      {
	if ( dest == id ) continue;
	MPI_Wait ( recvrequest + dest, &status );
	paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int vert = (*distantdofs[dest])[ii++];
	    int isdistghost = (*distantdofs[dest])[ii++];
// 	    (*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = vert;
// 	    (*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];
	    (*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = vert;
	    if ( dest < id && !isdistghost )
		  paralleldofs->ismasterdof.Clear ( vert ) ;
	    ii++; cnt_nexdof[dest]++;
	  }
      }

 
    exdof = 0;
    cnt_nexdof = 0;


    // *****************
    // Parallel Vertex dofs
    // *****************

    owndofs[0]->SetSize(1);
    (*owndofs[0])[0] = ndof;
    distantdofs[0]->SetSize(0);
    
    // find local and distant dof numbers for vertex exchange dofs
    for ( int vert = 0; vert < ma.GetNV(); vert++ )
      {
	int dest = 0;
	
	int distvert = parallelma -> GetDistantPNum ( dest, vert );
	owndofs[0]->Append ( distvert );
	paralleldofs->SetExchangeDof ( dest, vert );
	owndofs[0]->Append ( vert );
      }   
    
    int dest = 0;
    MyMPI_ISend ( *owndofs[0], dest, sendrequest[dest]);
    MyMPI_IRecv ( *distantdofs[0], dest, recvrequest[dest] );
   
    MPI_Wait ( recvrequest + dest, &status );

    ii = 1;
    while ( ii < distantdofs[0]->Size() )
      {
	if ( dest == id ) continue;
	paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
	int vnum = (*distantdofs[0])[ii++];
// 	(*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = vnum;
// 	(*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[0])[ii];
	(*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = vnum;
	ii++; cnt_nexdof[dest]++;
      }


    for ( int dest = id+1; dest < ntasks; dest++ )
      QuickSort ( (*(paralleldofs->sorted_exchangedof))[dest] );

     for ( int i = 0; i < ntasks; i++ )
       delete distantdofs[i], owndofs[i];

     delete [] owndofs, distantdofs;
     delete [] sendrequest, recvrequest;

  }
#endif



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

#ifdef PARALLEL
    paralleldofs = 0;
#endif

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

#ifdef PARALLEL
    paralleldofs = 0;
#endif

  }
  
  ElementFESpace :: ~ElementFESpace ()
  {
    ;
  }

  void ElementFESpace :: Update(LocalHeap & lh)
  {
    while (ma.GetNLevels() > ndlevel.Size())
      ndlevel.Append (n_el_dofs * ma.GetNE());

    FinalizeUpdate (lh);

#ifdef PARALLEL
    UpdateParallelDofs();
#endif
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

#ifdef PARALLEL

void ElementFESpace :: UpdateParallelDofs_hoproc()
  {
    
     // ******************************
     // update exchange dofs 
     // ******************************
    *testout << "Element::UpdateParallelDofs_hoproc" << endl;
    // Find number of exchange dofs
    Array<int> nexdof(ntasks);
    nexdof = 0;

    const MeshAccess & ma = (*this). GetMeshAccess();

    paralleldofs->SetNExDof(nexdof);

//     paralleldofs->localexchangedof = new Table<int> (nexdof);
//     paralleldofs->distantexchangedof = new Table<int> (nexdof);
    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);

    for ( int i = 0; i < ma.GetNE(); i++ )
      {
	paralleldofs->ClearExchangeDof(i);
	for ( int dest = 0; dest < ntasks; dest++ )
	  paralleldofs->ClearExchangeDof(dest, i);
      }


    int * ndof_dist = new int[ntasks];
    int ndof = GetNDof();
    MPI_Allgather ( &ndof, 1, MPI_INT, &ndof_dist[0], 
		    1, MPI_INT, MPI_COMM_WORLD);
    for ( int dest = 0; dest < ntasks; dest++ )
      paralleldofs -> SetDistNDof( dest, ndof_dist[dest]) ;

    delete [] ndof_dist;

  }

  void ElementFESpace :: UpdateParallelDofs_loproc()
  {
    *testout << "Element::UpdateParallelDofs_loproc" << endl;

    const MeshAccess & ma = (*this). GetMeshAccess();

    int ndof = GetNDof();

    // Find number of exchange dofs
    Array<int> nexdof(ntasks); 
    nexdof = 0;

    paralleldofs->SetNExDof(nexdof);

//     paralleldofs->localexchangedof = new Table<int> (nexdof);
//     paralleldofs->distantexchangedof = new Table<int> (nexdof);
    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);

    int * ndof_dist = new int[ntasks];
    MPI_Allgather ( &ndof, 1, MPI_INT, &ndof_dist[0], 
		    1, MPI_INT, MPI_COMM_WORLD);
    for ( int dest = 0; dest < ntasks; dest++ )
      paralleldofs -> SetDistNDof( dest, ndof_dist[dest]) ;

    delete [] ndof_dist;
 
  }
#endif // PARALLEL


 
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

#ifdef PARALLEL
    paralleldofs = 0;
#endif

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

    FinalizeUpdate (lh);

#ifdef PARALLEL
    UpdateParallelDofs();  
#endif
  }

  const FiniteElement & SurfaceElementFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    throw Exception ("SurfaceElementFESpace::GetFE not available");
  }

  
  void SurfaceElementFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize (0);
    // throw Exception ("SurfaceElementFESpace::GetDofNrs not available");
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





 









  







#ifdef ABC

  NonConformingFESpace :: 
      NonConformingFESpace (const MeshAccess & ama, const Flags& flags, bool parseflags)
  : FESpace (ama, flags),
  node2face2d (NULL), node2face3d(NULL)
  {
    name="NonConformingFESpace";
    if(parseflags) CheckFlags(flags);
    
    node2face2d = new HashTable<INT<2>,int>(10000);
    prol = new NonConformingProlongation (GetMeshAccess(), *this);

#ifdef PARALLEL
    paralleldofs = 0;
#endif

    // Update();
  }


  NonConformingFESpace :: ~NonConformingFESpace ()
  {
    ;
  }

  void NonConformingFESpace :: Update(LocalHeap & lh)
  {
    cout << "update non-conforming space currently not implemented" << endl;
#ifdef NONE

    const MeshAccess & ma = GetMeshAccess();
    int i, j, k;

    int level = ma.GetNLevels();
    int ne = ma.GetNE();
    int nse = ma.GetNSE();
    int np = ma.GetNP();

    if (level == nflevel.Size())
      return;

    Array<int> pnums;

    elementfaces.SetSize (ne);
    surfelementfaces.SetSize (nse);

    static int loctrigfaces[3][2] =
      { { 2, 3 },
	{ 1, 3 },
	{ 1, 2 } 
      };



    int nface = 0;
    if (nflevel.Size())
      nface = nflevel.Last();


    // insert faces from split edges:
    for (i = 1; i <= np; i++)
      {
	int parents[2];
	ma.GetParentNodes (i, parents);
	if (parents[0] && parents[1])
	  {
	    INT<2> face(parents[0], parents[1]);
	    face.Sort ();

	    if (!node2face2d->Used (face))
	      {
		nface++;
		node2face2d->Set (face, nface);
		faces.Append (face);
		finelevelofedge.Append (level);
	      }
	  }
      }  

    for (i = 1; i <= ne; i++)
      {
	ma.GetElPNums (i, pnums);
	int nelfaces = 0;
	int (*elfaces)[2];
      
	switch (ma.GetElType (i))
	  {
	  case ET_TRIG:
	    {
	      nelfaces = 3;
	      elfaces = loctrigfaces;
	      break;
	    }
	  }

	for (j = 0; j < nelfaces; j++)
	  {
	    int facenum;
	  
	    INT<2> face(pnums.Get(elfaces[j][0]),
			pnums.Get(elfaces[j][1]));
	  
	    face.Sort();

	    if (node2face2d->Used (face))
	      facenum = node2face2d->Get(face);
	    else
	      {
		nface++;
		node2face2d->Set (face, nface);
		faces.Append (face);
		finelevelofedge.Append (level);
		facenum = nface;
	      }
	  
	    elementfaces.Elem(i)[j] = facenum;
	    finelevelofedge.Elem(facenum) = level;
	  }
      }

    cout << "update surfelements" << endl;

    BitArray surfdof(nface);
    surfdof.Clear();

    for (i = 1; i <= nse; i++)
      {
	ma.GetSElPNums (i, pnums);

	INT<2> face (pnums.Get(1), pnums.Get(2));
	face.Sort();

	surfelementfaces.Elem(i) = 
	  node2face2d->Get (face);

	surfdof.Set (surfelementfaces.Elem(i));
      }

    cout << "Build Hierarchy" << endl;

    parentfaces.SetSize (nface);
    for (i = 1; i <= nface; i++)
      for (j = 0; j < 5; j++)
	parentfaces.Elem(i)[j] = 0;
  
    for (i = 1; i <= nface; i++)
      {
	INT<2> i2 (faces.Get(i).I1(), faces.Get(i).I2());
	int pa1[2], pa2[2];
	ma.GetParentNodes (i2[0], pa1);
	ma.GetParentNodes (i2[1], pa2);

	int issplitface = 0;
	if (pa1[0] == i2.I2() || pa1[1] == i2.I2())
	  issplitface = 1;
	if (pa2[0] == i2.I1() || pa2[1] == i2.I1())
	  issplitface = 2;

	if (issplitface)
	  {
	    // edge is obtained by splitting one edge into two parts:
	    INT_2 paface;
	    if (issplitface == 1)
	      paface = INT_2 (pa1[0], pa1[1]);
	    else
	      paface = INT_2 (pa2[0], pa2[1]);
	    paface.Sort();
	  
	    int pafacenr = node2face2d->Get (paface);
	    //	  (*testout) << "set split pa of " << i << ": " << pafacenr << endl;
	    /*
	      parentfaces.Elem(i).I1() = pafacenr;
	      parentfaces.Elem(i).I2() = pafacenr;
	    */
	  }
	else
	  {
	    // edge is splitting edge in middle of triangle:
	    for (j = 1; j <= 2; j++)
	      {
		INT_2 paface1, paface2;
		if (j == 1)
		  {
		    paface1 = INT_2 (pa1[0], i2.I2());
		    paface2 = INT_2 (pa1[1], i2.I2());
		  }
		else
		  {
		    paface1 = INT_2 (pa2[0], i2.I1());
		    paface2 = INT_2 (pa2[1], i2.I1());
		  }
		paface1.Sort();
		paface2.Sort();
	      
		int pafacenr1 = 0, pafacenr2 = 0;
	      
		if (node2face2d->Used (paface1) && node2face2d->Used (paface2))
		  {
		    pafacenr1 = node2face2d->Get (paface1);
		    pafacenr2 = node2face2d->Get (paface2);
		  
		    parentfaces.Elem(i)[0] = pafacenr1;	      
		    parentfaces.Elem(i)[1] = pafacenr2;	      
		    parentfaces.Elem(i)[2] = 0;	      
		    parentfaces.Elem(i)[3] = 0;	      
		    parentfaces.Elem(i)[4] = 0;	      


		    INT_2 splited, son1, son2;
		    if (j == 1)
		      {
			splited = INT_2 (pa1[0], pa1[1]);
			son1.I1() = son2.I1() = i2.I1();
		      }
		    else
		      {
			splited = INT_2 (pa2[0], pa2[1]);
			son1.I1() = son2.I1() = i2.I2();
		      }
		    son1.I2() = splited.I1();
		    son2.I2() = splited.I2();

		    splited.Sort();
		    son1.Sort();
		    son2.Sort();
		    /*
		      (*testout) << "split = " << splited 
		      << ", sons = " << son1 << son2 << endl;
		    */
		    int son1nr = node2face2d->Get (son1);
		    int son2nr = node2face2d->Get (son2);
		    int splitednr = node2face2d->Get (splited);
		    /*
		      (*testout) << "splitnr = " << splitednr
		      << ", sonsnr = " << son1nr << " " << son2nr << endl;
		    */
		  
		    parentfaces.Elem(son1nr)[0] = splitednr;

		    parentfaces.Elem(son1nr)[3] = 
		      parentfaces.Elem(son1nr)[1];
		    parentfaces.Elem(son1nr)[4] = 
		      parentfaces.Elem(son1nr)[2];
		    parentfaces.Elem(son1nr)[1] = pafacenr1;
		    parentfaces.Elem(son1nr)[2] = pafacenr2;

		    parentfaces.Elem(son2nr)[0] = splitednr;
		    parentfaces.Elem(son2nr)[3] = 
		      parentfaces.Elem(son2nr)[1];
		    parentfaces.Elem(son2nr)[4] = 
		      parentfaces.Elem(son2nr)[2];
		    parentfaces.Elem(son2nr)[1] = pafacenr2;
		    parentfaces.Elem(son2nr)[2] = pafacenr1;

		    /*
		      if (surfdof.Test(son1nr))
		      for (j = 1; j < 5; j++)
		      {
		      parentfaces.Elem(son1nr)[j] = splitednr;
		      parentfaces.Elem(son2nr)[j] = splitednr;
		      }
		    */
		  }
	      }
	  }
      }


    cout << "done" << endl;

    (*testout) << "non-conf elements: " << endl;
    for (i = 1; i <= ne; i++)
      (*testout) << elementfaces.Get(i)[0] << " - "
		 << elementfaces.Get(i)[1] << " - "
		 << elementfaces.Get(i)[2] << endl;
  
    (*testout) << "parents:" << endl;
    for (i = 1; i <= nface; i++)
      {
	for (j = 0; j < 5; j++)
	  (*testout) << parentfaces.Elem(i)[j] << " - " ;
	(*testout) << endl;
      }

    nflevel.Append (nface);

#ifdef PARALLEL
    UpdateParallelDofs();
#endif

#endif
  }





  int NonConformingFESpace :: GetNDof () const
  {
    return nflevel.Last();
  }

  int NonConformingFESpace :: GetNDofLevel (int level) const
  {
    return nflevel[level];
  }

  const FiniteElement & NonConformingFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    const MeshAccess & ma = GetMeshAccess();

    switch (ma.GetElType(elnr))
      {
      case ET_TRIG:
	return trig1;
      case ET_TET:
	return tet1;
      default:
        cerr << "NonConformingFESpace, GetFE: unknown type" << endl;
      }
    return tet1;
  }
  
  void NonConformingFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    const MeshAccess & ma = GetMeshAccess();

    switch (ma.GetElType(elnr))
      {
      case ET_TRIG:
	{
	  dnums.SetSize(3);
	  break;
	}
      case ET_TET:
	{
	  dnums.SetSize(4);
	  break;
	}
      default:
        cerr << "NonConformingFESpace, GetFE: unknown type" << endl;
      }
    for (int j = 0; j < dnums.Size(); j++)
      dnums[j] = abs (elementfaces[elnr][j]);
  }

  const FiniteElement & NonConformingFESpace :: GetSFE (int selnr, LocalHeap & lh) const
  {
    const MeshAccess & ma = GetMeshAccess();

    switch (ma.GetSElType(selnr))
      {
      case ET_TRIG:
	return trig1;
      case ET_SEGM:
	return segm1;
      default:
        throw Exception ("NonconformingFESpace::GetSFE unsupported element");
      }  
    return trig1;
  }

  void NonConformingFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    dnums.SetSize(1);
    dnums[0] = surfelementfaces[selnr];
  }

#endif




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

#ifdef PARALLEL
    if ( ntasks == 1 ) return;
    Flags loflags;
    loflags.SetFlag("order",0.0);
    paralleldofs = 0;
    if (dgjumps){ *testout << "(CompFES:)setting loflag dgjumps " << endl; loflags.SetFlag ("dgjumps");}
    /*
    Array<const FESpace*> lospaces (0);
    for ( int i = 0; i < spaces.Size(); i++)
      if ( &spaces[i]->LowOrderFESpace() )
	lospaces . Append ( &spaces[i]->LowOrderFESpace() );
      else
	{ low_order_space = 0; return; }
    */
    low_order_space = new CompoundFESpace(ma, lospaces, flags, parseflags);
#endif
  }







  CompoundFESpace :: CompoundFESpace (const MeshAccess & ama,
				      const Array<const FESpace*> & aspaces,
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

#ifdef PARALLEL
    if ( ntasks == 1 ) return;
    Flags loflags;
    loflags.SetFlag("order",0.0);
    paralleldofs = 0;
    if (dgjumps){ *testout << "(CompFES:)setting loflag dgjumps " << endl; loflags.SetFlag ("dgjumps");}
    Array<const FESpace*> lospaces (0);
    for ( int i = 0; i < spaces.Size(); i++)
      if ( &spaces[i]->LowOrderFESpace() )
	lospaces . Append ( &spaces[i]->LowOrderFESpace() );
      else
	{ low_order_space = 0; return; }
    low_order_space = new CompoundFESpace(ma, lospaces, flags, parseflags);
#endif
  }


  void CompoundFESpace :: AddSpace (const FESpace * fes)
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
    FinalizeUpdate (lh);

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

#ifdef PARALLEL
    if (low_order_space)
      low_order_space -> Update(lh);
    UpdateParallelDofs();
#endif

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
    
    // void * mem = lh.Alloc (sizeof(CompoundFiniteElement));
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
    // int j, k; 
    LocalHeapMem<100005> lh("CompoundFESpace - transformmat");
    for (int i = 0; i < spaces.Size(); i++)
      {
	int nd = boundary  ?
	  spaces[i]->GetSFE(elnr, lh).GetNDof()
	  : spaces[i]->GetFE(elnr, lh).GetNDof();

	lh.CleanUp();

	/*
	typedef typename MAT::TSCAL TSCAL;
	Matrix<TSCAL> rmat(nd, mat.Width());
	Matrix<TSCAL> cmat(mat.Height(), nd);
	*/ 

	/*
	for (j = 0; j < nd; j++)
	  for (k = 0; k < mat.Width(); k++)
	    rmat(j,k) = mat(base+j, k);

	spaces[i]->TransformMat (elnr, boundary, rmat, TRANSFORM_MAT_LEFT);

	for (j = 0; j < nd; j++)
	  for (k = 0; k < mat.Width(); k++)
	    mat(base+j, k) = rmat(j,k);
	*/

	spaces[i]->TransformMat (elnr, boundary, mat.Rows(base, base+nd), TRANSFORM_MAT_LEFT);

	/*
	for (j = 0; j < nd; j++)
	  for (k = 0; k < mat.Height(); k++)
	    cmat(k,j) = mat(k,base+j);

	spaces[i]->TransformMat (elnr, boundary, cmat, TRANSFORM_MAT_RIGHT);

	for (j = 0; j < nd; j++)
	  for (k = 0; k < mat.Height(); k++)
	    mat(k,base+j) = cmat(k,j);
	*/
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

	VEC svec (nd, &vec(base));
	spaces[i]->TransformVec (elnr, boundary, svec, tt);
	base += nd;
      }
  }


#ifdef PARALLEL
  void CompoundFESpace :: UpdateParallelDofs_loproc()
  {
    // gather the paralleldofs of all the fespaces to cumulatedparalleldof
    paralleldofs -> UpdateCompound ();   
  }

  void CompoundFESpace :: UpdateParallelDofs_hoproc()
  {
    paralleldofs -> UpdateCompound ();   
  }

#endif


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



#ifdef PARALLEL

  /*
  ParallelElementFESpace :: ParallelElementFESpace (const MeshAccess & ama,
				    int aorder, int adim, bool acomplex)
    : ElementFESpace (ama, aorder, adim, acomplex)
  {
    ;
  }
  */

  ParallelElementFESpace :: ParallelElementFESpace (const MeshAccess & ama, const Flags& flags, bool parseflags)
    : ElementFESpace (ama, flags, parseflags)
  {
    ;
  }

void ParallelElementFESpace :: UpdateParallelDofs_hoproc()
  {

     // ******************************
     // update exchange dofs 
     // ******************************
    *testout << "Element::UpdateParallelDofs_hoproc" << endl;
    // Find number of exchange dofs
    Array<int> nexdof(ntasks);
    nexdof = 0;

    const MeshAccess & ma = (*this). GetMeshAccess();
    
    paralleldofs->SetNExDof(nexdof);

//     paralleldofs->localexchangedof = new Table<int> (nexdof);
//     paralleldofs->distantexchangedof = new Table<int> (nexdof);
    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);

    for ( int i = 0; i < ma.GetNE(); i++ )
      {
	paralleldofs->ClearExchangeDof(i);
	for ( int dest = 0; dest < ntasks; dest++ )
	  paralleldofs->ClearExchangeDof(dest, i);
      }


    int * ndof_dist = new int[ntasks];
    int ndof = GetNDof();
    MPI_Allgather ( &ndof, 1, MPI_INT, &ndof_dist[1], 
		    1, MPI_INT, MPI_HIGHORDER_COMM);
    for ( int dest = 0; dest < ntasks; dest++ )
      paralleldofs -> SetDistNDof( dest, ndof_dist[dest]) ;

    delete [] ndof_dist;

  }

  void ParallelElementFESpace :: UpdateParallelDofs_loproc()
  {
    *testout << "Element::UpdateParallelDofs_loproc" << endl;

    const MeshAccess & ma = (*this). GetMeshAccess();

    int ndof = GetNDof();

    // Find number of exchange dofs
    Array<int> nexdof(ntasks); 
    nexdof = 0;

    MPI_Status status;
    MPI_Request * sendrequest = new MPI_Request [ntasks];
    MPI_Request * recvrequest = new MPI_Request [ntasks];

    // number of element exchange dofs
    for ( int el = 0; el < ma.GetNE(); el++ )
      {
	nexdof[id] ++;//= dnums.Size() ; 
	for ( int dest = 1; dest < ntasks; dest ++ )
	  if (  parallelma -> GetDistantElNum ( dest, el ) >= 0 )
	    { 
	      nexdof[dest] ++; 
	    }
      }

    paralleldofs->SetNExDof(nexdof);

//     paralleldofs->localexchangedof = new Table<int> (nexdof);
//     paralleldofs->distantexchangedof = new Table<int> (nexdof);
    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);



    Array<int> ** owndofs,** distantdofs;
    owndofs = new Array<int> * [ntasks];
    distantdofs = new Array<int> * [ntasks];

    for ( int i = 0; i < ntasks; i++)
      {
	owndofs[i] = new Array<int> (1);
	(*owndofs[i])[0] = ndof;
	distantdofs[i] = new Array<int> (0);
      }

    int exdof = 0;
    Array<int> cnt_nexdof(ntasks);
    cnt_nexdof = 0;


    // *****************
    // Parallel Element dofs
    // *****************


    // find local and distant dof numbers for edge exchange dofs
    for ( int el = 0; el < ma.GetNE(); el++ )
      {
	(*(paralleldofs->sorted_exchangedof))[id][exdof++] = el;

	for ( int dest = 1; dest < ntasks; dest++ )
	  {
	    int distel = parallelma -> GetDistantElNum ( dest, el );
	    if( distel < 0 ) continue;

	    owndofs[dest]->Append ( distel );
	    paralleldofs->SetExchangeDof ( dest, el );
	    paralleldofs->SetExchangeDof ( el );
	    owndofs[dest]->Append ( el );
  	  }
      }   


    for ( int dest = 1; dest < ntasks; dest ++ )
      {
	MyMPI_ISend ( *owndofs[dest], dest, sendrequest[dest]);
	MyMPI_IRecv ( *distantdofs[dest], dest, recvrequest[dest]);
      }
    


    int ii = 1;
    for( int dest=1; dest<ntasks; dest++) 
      {
	if ( dest == id ) continue;
	MPI_Wait(recvrequest+dest, &status);
	paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int elnum = (*distantdofs[dest])[ii++];
// 	    (*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = elnum;
// 	    (*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];
	    (*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = elnum;
	    ii++; cnt_nexdof[dest]++;
	     
	  }
      }

    for ( int dest = id+1; dest < ntasks; dest++ )
      QuickSort ( (*(paralleldofs->sorted_exchangedof))[dest] );

    for ( int i = 0; i < ntasks; i++ )
      {
	delete distantdofs[i];
	delete owndofs[i];
      }

    delete [] owndofs;
    delete [] distantdofs;
    delete [] sendrequest;
    delete [] recvrequest;

 
  }


  
  ParallelNodalFESpace :: ParallelNodalFESpace (const MeshAccess & ama,
				const Flags & flags,
                                bool parseflags)
    : NodalFESpace (ama, flags, parseflags)
  {
    name="ParallelNodalFESpace";
  }


  void ParallelNodalFESpace :: UpdateParallelDofs_loproc()
  {
    *testout << "ParallelNodal::UpdateParallelDofs_loproc" << endl;

    const MeshAccess & ma = (*this). GetMeshAccess();

    int ndof = GetNDof();

    // Find number of exchange dofs
    Array<int> nexdof(ntasks); 
    nexdof = 0;


    // number of vertex exchange dofs
    for ( int vert = 0; vert < ma.GetNV(); vert++ )
      {
	nexdof[id] ++;//= dnums.Size() ; 
	for ( int dest = 1; dest < ntasks; dest ++ )
	  if (  parallelma -> GetDistantPNum ( dest, vert ) >= 0 )
	    { 
	      nexdof[dest] ++; 
	    }
      }

    paralleldofs->SetNExDof(nexdof);

//     paralleldofs->localexchangedof = new Table<int> (nexdof);
//     paralleldofs->distantexchangedof = new Table<int> (nexdof);
    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);

    MPI_Status status;
    MPI_Request * sendrequest = new MPI_Request [ntasks];
    MPI_Request * recvrequest = new MPI_Request [ntasks];

    Array<int> ** owndofs,** distantdofs;
    owndofs = new Array<int> * [ntasks];
    distantdofs = new Array<int> * [ntasks];

    for ( int i = 0; i < ntasks; i++)
      {
	owndofs[i] = new Array<int> (1);
	(*owndofs[i])[0] = ndof;
	distantdofs[i] = new Array<int> (0);
      }

    int exdof = 0;
    Array<int> cnt_nexdof(ntasks);
    cnt_nexdof = 0;

    // *****************
    // Parallel Vertex dofs
    // *****************


    // find local and distant dof numbers for vertex exchange dofs
    for ( int vert = 0; vert < ma.GetNV(); vert++ )
      {
	(*(paralleldofs->sorted_exchangedof))[id][exdof++] = vert;

	for ( int dest = 1; dest < ntasks; dest++ )
	  {
	    int distvert = parallelma -> GetDistantPNum ( dest, vert );
	    if( distvert < 0 ) continue;

	    owndofs[dest]->Append ( distvert );
	    paralleldofs->SetExchangeDof ( dest, vert );
	    paralleldofs->SetExchangeDof ( vert );
	    owndofs[dest]->Append ( vert );
  	  }
      }   


    for ( int dest = 1; dest < ntasks; dest ++ )
      {
	MyMPI_ISend ( *owndofs[dest], dest, sendrequest[dest]);
	MyMPI_IRecv ( *distantdofs[dest], dest, recvrequest[dest]);
      }
    
    int ii = 1;
    for( int dest=1; dest<ntasks; dest++) 
      {
	if ( dest == id ) continue;
	MPI_Wait ( recvrequest + dest, &status );
	paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int vnum = (*distantdofs[dest])[ii++];
// 	    (*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = vnum;
// 	    (*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];

	    (*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = vnum;
		ii++; cnt_nexdof[dest]++;
	     
	  }
      }

    for ( int dest = id+1; dest < ntasks; dest++ )
      QuickSort ( (*(paralleldofs->sorted_exchangedof))[dest] );

     for ( int i = 0; i < ntasks; i++ )
       delete distantdofs[i], owndofs[i];

     delete [] owndofs, distantdofs;
     delete [] sendrequest, recvrequest;
  }

  void ParallelNodalFESpace :: UpdateParallelDofs_hoproc ()
  {
    *testout << "ParallelNodal::UpdateParallelDofs_hoproc" << endl;

    const MeshAccess & ma = (*this). GetMeshAccess();

    int ndof = GetNDof();

    // Find number of exchange dofs
    Array<int> nexdof(ntasks); 
    nexdof = 0;

    MPI_Status status;
    MPI_Request * sendrequest = new MPI_Request [ntasks];
    MPI_Request * recvrequest = new MPI_Request [ntasks];

    // number of vertex exchange dofs
    for ( int vert = 0; vert < ma.GetNV(); vert++ )
      {
	if ( !parallelma->IsExchangeVert ( vert ) ) continue;

	nexdof[id] ++;
	for ( int dest = 1; dest < ntasks; dest ++ )
	  if (  parallelma -> GetDistantPNum ( dest, vert ) >= 0 )
	    { 
	      nexdof[dest] ++; 
	    }
      }


    nexdof[0] = 0;

    paralleldofs->SetNExDof(nexdof);

//     paralleldofs->localexchangedof = new Table<int> (nexdof);
//     paralleldofs->distantexchangedof = new Table<int> (nexdof);
 
    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);

    Array<int> ** owndofs,** distantdofs;
    owndofs = new Array<int> * [ntasks];
    distantdofs = new Array<int> * [ntasks];

    for ( int i = 0; i < ntasks; i++)
      {
	owndofs[i] = new Array<int> (1);
	(*owndofs[i])[0] = ndof;
	distantdofs[i] = new Array<int> (0);
      }

    int exdof = 0;
    Array<int> cnt_nexdof(ntasks);
    cnt_nexdof = 0;

    // *****************
    // Parallel Vertex dofs
    // *****************
    

    // find local and distant dof numbers for vertex exchange dofs
    for ( int vert = 0; vert < ma.GetNV(); vert++ )
      {
	if ( !parallelma->IsExchangeVert ( vert ) ) continue;
	Array<int> dnums;

	(*(paralleldofs->sorted_exchangedof))[id][exdof++] = vert;

	for ( int dest = 1; dest < ntasks; dest++ )
	  {
	    int distvert = parallelma -> GetDistantPNum ( dest, vert );
	    if( distvert < 0 ) continue;

	    owndofs[dest]->Append ( distvert );
	    owndofs[dest]->Append( int (paralleldofs->IsGhostDof(vert)) );
	    paralleldofs->SetExchangeDof ( dest, vert );
	    paralleldofs->SetExchangeDof ( vert );
	    owndofs[dest]->Append ( vert );
  	  }
      }   



     for ( int sender = 1; sender < ntasks; sender ++ )
      {
       if ( id == sender )
          for ( int dest = 1; dest < ntasks; dest ++ )
            if ( dest != id )
              {
		MyMPI_ISend ( (*owndofs[dest]), dest, sendrequest[dest] );
              }
	  
        if ( id != sender )
          {
	    int n2;
	    MyMPI_IRecv ( *distantdofs[sender], sender, recvrequest[sender]);
          }
 	  
	 
      }

    int ii = 1;
    for( int dest=1; dest<ntasks; dest++) 
      {
	if ( dest == id ) continue;
	MPI_Wait ( recvrequest + dest, &status );
	paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int vert = (*distantdofs[dest])[ii++];
	    int isdistghost = (*distantdofs[dest])[ii++];
// 	    (*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = vert;
// 	    (*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];

	    (*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = vert;
	    if ( dest < id && !isdistghost )
		  paralleldofs->ismasterdof.Clear ( vert ) ;
	    ii++; cnt_nexdof[dest]++;
	  }
      }

    for ( int dest = id+1; dest < ntasks; dest++ )
      QuickSort ( (*(paralleldofs->sorted_exchangedof))[dest] );

     for ( int i = 0; i < ntasks; i++ )
       delete distantdofs[i], owndofs[i];

     delete [] owndofs, distantdofs;
     delete [] sendrequest, recvrequest;
  }

#endif // PARALLEL










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





    
  // register FESpaces
  namespace fespace_cpp
  {
     
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("nonconforming", NonconformingFESpace::Create);
      GetFESpaceClasses().AddFESpace ("nodal", NodalFESpace::Create);
    }
    
    Init init;
  }




}



