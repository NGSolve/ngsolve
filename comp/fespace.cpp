/**********************************************************************/
/* File:   fespace.cpp                                                */
/* Author: Joachim Schoeberl                                          */
/* Date:   25. Mar. 2000                                              */
/**********************************************************************/

/* 
   Finite Element Space
*/

#pragma implementation "fespace.hpp"

#include <comp.hpp>
#include <multigrid.hpp>

#include "../fem/h1lofe.hpp"
#include <parallelngs.hpp>


using namespace ngmg;

using FE_Quad1 = ScalarFE<ET_QUAD,1>;

namespace ngcomp
{

  FESpace :: FESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags)
    : NGS_Object (ama, "FESpace")
  {
    // register flags
    DefineStringFlag("type");
    DefineNumFlag("order");
    DefineNumFlag("dim");
    DefineDefineFlag("vec");
    DefineDefineFlag("complex");
    DefineDefineFlag("timing");
    DefineDefineFlag("print");
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
    print = flags.GetDefineFlag("print");
    dgjumps = flags.GetDefineFlag("dgjumps");
    no_low_order_space = flags.GetDefineFlag("no_low_order_space");
    if (dgjumps) 
      *testout << "ATTENTION: flag dgjumps is used!\n This leads to a \
lot of new non-zero entries in the matrix!\n" << endl;
    // else *testout << "\n (" << order << ") flag dgjumps is not used!" << endl;
    
    if(flags.NumListFlagDefined("directsolverdomains"))
      {
	directsolverclustered.SetSize(ama->GetNDomains());
	directsolverclustered = false;
	Array<double> clusters(flags.GetNumListFlag("directsolverdomains"));
	for(int i=0; i<clusters.Size(); i++) 
	  directsolverclustered[static_cast<int>(clusters[i])-1] = true; // 1-based!!
      }
    
    if(flags.NumListFlagDefined("dirichlet"))
      {
	dirichlet_boundaries.SetSize (ma->GetNBoundaries());
	dirichlet_boundaries.Clear();
        for (double dbi : flags.GetNumListFlag("dirichlet"))
          {
	    int bnd = int(dbi-1);
	    if (bnd >= 0 && bnd < dirichlet_boundaries.Size())
	      dirichlet_boundaries.Set (bnd);
	    // else
            //   cerr << "Illegal Dirichlet boundary index " << bnd+1 << endl;
          }
	if (print)
	  *testout << "dirichlet_boundaries:" << endl << dirichlet_boundaries << endl;
      }
    
    if (flags.NumListFlagDefined("definedon") || 
        flags.NumFlagDefined("definedon") ||
        flags.StringListFlagDefined("definedon"))
      {
	definedon.SetSize (ma->GetNDomains());
	definedon = false;
	Array<double> defon;
	if (flags.NumListFlagDefined("definedon")) 
	  defon = flags.GetNumListFlag("definedon");
	else if (flags.NumFlagDefined("definedon"))
	  {
	    defon.SetSize(1);
	    defon[0] = flags.GetNumFlag("definedon",0);
	  }
        /*
	for(int i = 0; i< defon.Size(); i++)
	  if (defon[i] <= ma->GetNDomains() && defon[i] > 0)
	    definedon[int(defon[i])-1] = true;
        */
        for (int di : defon)
          if (di > 0 && di <= ma->GetNDomains())
            definedon[di-1] = true;
          
	if(flags.StringListFlagDefined("definedon"))
	  {
	    Array<string> dmaterials(flags.GetStringListFlag ("definedon").Size());
	    for(int i=0; i<dmaterials.Size(); i++)
	      dmaterials[i] = flags.GetStringListFlag ("definedon")[i];
	    for(int i = 0; i < ma->GetNDomains(); i++)
	      {
		for(int j = 0; j < dmaterials.Size(); j++)
		  if(StringFitsPattern(ma->GetDomainMaterial(i),dmaterials[j]))
		    {
		      definedon[i] = true;
		      break;
		    }
	      }
	  }

	// default:
	// fespace only defined on boundaries matching definedon-domains
	definedonbound.SetSize (ma->GetNBoundaries());
	definedonbound = false;
	for (int sel = 0; sel < ma->GetNSE(); sel++)
	  {
	    int index = ma->GetSElIndex(sel);
	    int dom1, dom2;
	    ma->GetSElNeighbouringDomains(sel, dom1, dom2);
	    dom1--; dom2--;
	    if ( dom1 >= 0 )
	      if ( definedon[dom1] )
		definedonbound[index] = true;

	    if ( dom2 >= 0 )
	      if ( definedon[dom2] )
		definedonbound[index] = true;
	  }
      }

    // additional boundaries
    if(flags.NumListFlagDefined("definedonbound")|| flags.NumFlagDefined("definedonbound") )
      {
	if ( definedonbound.Size() == 0 )
	  {
	    definedonbound.SetSize (ma->GetNBoundaries());
	    definedonbound = false;
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
	  if(defon[i] <= ma->GetNBoundaries() && defon[i] > 0)
	    definedonbound[int(defon[i])-1] = true;
      }
    

    else if(flags.StringListFlagDefined("definedonbound") || flags.StringFlagDefined("definedonbound"))
      {
	if ( definedonbound.Size() == 0 )
	  {
	    definedonbound.SetSize (ma->GetNBoundaries());
	    definedonbound = false;
	  }

	Array<string*> defon;

	if(flags.StringFlagDefined("definedonbound"))
	  defon.Append(new string(flags.GetStringFlag("definedonbound","")));
	else
	  for(int i=0; i<flags.GetStringListFlag ("definedonbound").Size(); i++)
	    defon.Append(new string(flags.GetStringListFlag("definedonbound")[i]));
	
	for(int selnum = 0; selnum < ma->GetNSE(); selnum++)
	  {
	    if(definedonbound[ma->GetSElIndex(selnum)] == false)
	      {
		for(int i=0; i<defon.Size(); i++)
		  {
		    if(StringFitsPattern(ma->GetSElBCName(selnum),*(defon[i])))	
		      {		
		 	definedonbound[ma->GetSElIndex(selnum)] = true;
			continue;
		      }
		  }
	      }
	  }
	for(int i=0; i<defon.Size(); i++)
	  delete defon[i];
      }
    
    level_updated = -1;


    point = NULL;
    segm = NULL;
    trig = NULL;
    quad = NULL;
    tet = NULL;
    prism = NULL;
    pyramid = NULL;
    hex = NULL;
    
    dummy_tet = new DummyFE<ET_TET>();
    dummy_pyramid = new DummyFE<ET_PYRAMID>();
    dummy_prism = new DummyFE<ET_PRISM>();
    dummy_hex = new DummyFE<ET_HEX>();
    dummy_trig = new DummyFE<ET_TRIG>();
    dummy_quad = new DummyFE<ET_QUAD>();
    dummy_segm = new DummyFE<ET_SEGM>();
    dummy_point = new DummyFE<ET_POINT>();

    evaluator = NULL; 
    boundary_evaluator = NULL;
    flux_evaluator = NULL;

    integrator = NULL;
    boundary_integrator = NULL;
    low_order_space = NULL;
    prol = NULL;


    // element_coloring = NULL;
    // selement_coloring = NULL;
    paralleldofs = NULL;

    ctofdof.SetSize(0);
  }

  
  FESpace :: ~FESpace ()
  {
    delete tet;
    delete pyramid;
    delete prism;
    delete hex;
    delete trig;
    delete quad;
    delete segm;
    delete point;

    delete dummy_tet;
    delete dummy_pyramid;
    delete dummy_prism;
    delete dummy_hex;
    delete dummy_trig;
    delete dummy_quad;
    delete dummy_segm;
    delete dummy_point;

    delete paralleldofs;
  }
  

  void FESpace :: Update(LocalHeap & lh)
  {
    if (print)
      {
 	*testout << "Update FESpace, type = " << typeid(*this).name() << endl;
	*testout << "name = " << name << endl;
      }

    for (int i = 0; i < specialelements.Size(); i++)
      delete specialelements[i]; 
    specialelements.SetSize(0);


    int dim = ma->GetDimension();
    
    dirichlet_vertex.SetSize (ma->GetNV());
    dirichlet_edge.SetSize (ma->GetNEdges());
    if (dim == 3)
      dirichlet_face.SetSize (ma->GetNFaces());
    
    dirichlet_vertex = false;
    dirichlet_edge = false;
    dirichlet_face = false;

#pragma omp parallel
    {
      if (dirichlet_boundaries.Size())
        for (Ngs_Element ngel : ma->Elements(BND).OmpSplit())
          if (dirichlet_boundaries[ngel.GetIndex()])
            {
              dirichlet_vertex[ngel.Vertices()] = true;
              dirichlet_edge[ngel.Edges()] = true;
              if (dim == 3)
                dirichlet_face[ngel.Faces()[0]] = true;
            }
    }

    if (print)
      {
	(*testout) << "Dirichlet_vertex,1 = " << endl << dirichlet_vertex << endl;
	(*testout) << "Dirichlet_edge,1 = " << endl << dirichlet_edge << endl;
	(*testout) << "Dirichlet_face,1 = " << endl << dirichlet_face << endl;
      }


    ma->AllReduceNodalData (NT_VERTEX, dirichlet_vertex, MPI_LOR);
    ma->AllReduceNodalData (NT_EDGE, dirichlet_edge, MPI_LOR);
    ma->AllReduceNodalData (NT_FACE, dirichlet_face, MPI_LOR);
    
    if (print)
      {
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

    dirichlet_dofs.SetSize (GetNDof());
    dirichlet_dofs.Clear();

    /*
    if (dirichlet_boundaries.Size())
      for (ElementId ei : ma->Elements<BND>())
	if (dirichlet_boundaries[ma->GetElIndex(ei)])
	  {
	    GetDofNrs (ei, dnums);
	    for (int d : dnums)
	      if (d != -1) dirichlet_dofs.Set (d);
	  }
    */

    if (dirichlet_boundaries.Size())
      for (FESpace::Element el : Elements(BND))
        if (dirichlet_boundaries[el.GetIndex()])
          for (int d : el.GetDofs())
            if (d != -1) dirichlet_dofs.Set (d);

    Array<int> dnums;
    for (int i : Range(dirichlet_vertex))
      if (dirichlet_vertex[i])
	{
	  GetVertexDofNrs (i, dnums);
	  for (int d : dnums)
	    if (d != -1) dirichlet_dofs.Set (d);
	}

    for (int i : Range(dirichlet_edge))
      if (dirichlet_edge[i])
	{
	  GetEdgeDofNrs (i, dnums);
	  for (int d : dnums)
	    if (d != -1) dirichlet_dofs.Set (d);
	}

    for (int i : Range(dirichlet_face))
      if (dirichlet_face[i])
	{
	  GetFaceDofNrs (i, dnums);
	  for (int d : dnums)
	    if (d != -1) dirichlet_dofs.Set (d);
	}


    free_dofs.SetSize (GetNDof());
    free_dofs = dirichlet_dofs;
    free_dofs.Invert();
    
    for (auto i : Range(ctofdof))
      if (ctofdof[i] == UNUSED_DOF)
	free_dofs.Clear(i);

    external_free_dofs.SetSize (GetNDof());
    external_free_dofs = free_dofs;
    for (auto i : Range(ctofdof))
      if (ctofdof[i] & LOCAL_DOF)
	external_free_dofs.Clear(i);

    if (print)
      *testout << "freedofs = " << endl << free_dofs << endl;
    
    UpdateParallelDofs();

    if (print)
      *testout << "coloring ... " << flush;


    for (auto vb = VOL; vb <= BND; vb++)
      {
        Array<int> col(ma->GetNE(vb));
        col = -1;
        // bool found;
        int maxcolor = 0;
        
        int basecol = 0;
        Array<unsigned int> mask(GetNDof());

        int cnt = 0, found = 0;
        for (ElementId el : Elements(vb)) { cnt++; (void)el; } // no warning 

        do
          {
            mask = 0;

            for (auto el : Elements(vb))
              {
                if (col[el.Nr()] >= 0) continue;

                unsigned check = 0;
                for (auto d : el.GetDofs())
                  if (d != -1) check |= mask[d];

                if (check != UINT_MAX) // 0xFFFFFFFF)
                  {
                    found++;
                    unsigned checkbit = 1;
                    int color = basecol;
                    while (check & checkbit)
                      {
                        color++;
                        checkbit *= 2;
                      }

                    col[el.Nr()] = color;
                    if (color > maxcolor) maxcolor = color;
		
                    for (auto d : el.GetDofs())
                      if (d != -1) mask[d] |= checkbit;
                  }
              }
            
            basecol += 8*sizeof(unsigned int); // 32;
          }
        while (found < cnt);


	/*
	  // didn't help anything ...
	cout << "have " << maxcolor+1 << " colors, try to reduce ***************** " << endl;

	DynamicTable<int> dof2color(GetNDof());
	for (auto el : Elements(vb))
	  for (auto d : el.GetDofs())
	    dof2color.Add(d, col[el.Nr()]);

	*testout << "dof2color = " << endl << dof2color << endl;

	for (int c = 0; c <= maxcolor; c++)
	  {
	    *testout << "try to reschedule color " << c << endl;
	    
	    bool remove_col = true;
	    int cntmove = 0, cnttot = 0;
	    
	    for (auto el : Elements(vb))
	      if (col[el.Nr()] == c)
		{
		  bool resched = false;

		  for (int newcol = c+1; newcol <= maxcolor; newcol++)
		    {
		      bool possible = true;
		      for (auto d : el.GetDofs())
			if (dof2color[d].Contains(newcol))
			  {
			    possible = false; 
			    break;
			  }
		      if (possible) resched = true;
		    }

		  if (resched) cntmove++;
		  cnttot++;
		  if (!resched) 
		    {
		      *testout << "cannot move el with dofs " << el.GetDofs() << endl;
		      remove_col = false;
		    }
		}

	    *testout << "could move " << cntmove << " out of " << cnttot << endl;
	    if (remove_col) cout << "could remove color" << endl;
	  }

	*/



        Array<int> cntcol(maxcolor+1);
        cntcol = 0;
        for (ElementId el : Elements(vb))
          cntcol[col[el.Nr()]]++;

        Table<int> & coloring = (vb == VOL) ? element_coloring : selement_coloring;
        coloring = Table<int> (cntcol);

	cntcol = 0;
        for (ElementId el : Elements(vb))
          coloring[col[el.Nr()]][cntcol[col[el.Nr()]]++] = el.Nr();

        if (print)
          *testout << "needed " << maxcolor+1 << " colors" 
                   << " for " << ((vb == VOL) ? "vol" : "bnd") << endl;
      }


    level_updated = ma->GetNLevels();
    if (timing) Timing();

    // CheckCouplingTypes();
  }


  const FiniteElement & FESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    FiniteElement * fe = NULL;
    
    if (DefinedOn (ElementId (VOL, elnr)))
      {
        switch (ma->GetElType(elnr))
          {
          case ET_TET: fe = tet; break;
          case ET_PYRAMID: fe = pyramid; break;
          case ET_PRISM: fe = prism; break;
          case ET_HEX: fe = hex; break;
          case ET_TRIG: fe = trig; break;
          case ET_QUAD: fe = quad; break;
          case ET_SEGM: fe = segm; break;
          case ET_POINT: fe = point; break;
          }
      }
    else
      {
        switch (ma->GetElType(elnr))
          {
          case ET_TET: fe = dummy_tet; break;
          case ET_PYRAMID: fe = dummy_pyramid; break;
          case ET_PRISM: fe = dummy_prism; break;
          case ET_HEX: fe = dummy_hex; break;
          case ET_TRIG: fe = dummy_trig; break;
          case ET_QUAD: fe = dummy_quad; break;
          case ET_SEGM: fe = dummy_segm; break;
          case ET_POINT: fe = dummy_point; break;
          }
      }

    if (!fe)
      {
        Exception ex;
        ex << "FESpace" << GetClassName() << ", undefined eltype " 
           << ElementTopology::GetElementName(ma->GetElType(elnr))
           << ", order = " << ToString (order) << "\n";
        throw ex;
      }
    
    return *fe;
  }

  /*
    // not such a great idea ..
  void FESpace :: GetDofRanges (ElementId ei, Array<IntRange> & dranges) const
  {
    // cout << "getdofrangs called for fespace " << GetClassName() << endl;
    Array<int> dnums;
    GetDofNrs (ei, dnums);
    dranges.SetSize(0);
    for (int i = 0; i < dnums.Size(); i++)
      dranges.Append (IntRange (dnums[i], dnums[i]+1));
  }
  */

  /*
  FlatArray<int> FESpace :: GetDofNrs (ElementId ei, LocalHeap & lh) const
  {
    Vec<4,int> nnodes = ElementTopology::GetNNodes (ma->GetElType (ei));
    Array<IntRange> dranges(InnerProduct (nnodes, vefc_dofblocks), lh);
    dranges.SetSize(0);
    GetDofRanges (ei, dranges);
    int nd = 0;
    for (int i = 0; i < dranges.Size(); i++) nd += dranges[i].Size();
    
    FlatArray<int> dnums(nd, lh);
    for (int i = 0, cnt = 0; i < dranges.Size(); i++)
      {
        if (dranges[i].First() != -1)
          for (int j = dranges[i].First(); j < dranges[i].Next(); j++)
            dnums[cnt++] = j;
        else
          for (int j = dranges[i].First(); j < dranges[i].Next(); j++)
            dnums[cnt++] = -1;
      }

    return dnums;
  }
  */


  Table<int> FESpace :: CreateDofTable (VorB vorb) const
  {
    TableCreator<int> creator;
    
    for ( ; !creator.Done(); creator++)
      for (FESpace::Element el : Elements(vorb))
        creator.Add(el.Nr(), el.GetDofs());

    return creator.MoveTable();
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
    ArrayMem<int,100> dnums;
    GetDofNrs(elnr, dnums);
    ctypes.SetSize(dnums.Size());

    if (ctofdof.Size()==0)
      ctypes = INTERFACE_DOF;
    else
      {
        for (int i = 0; i < dnums.Size(); i++)
          if (dnums[i] != -1)
            ctypes[i] = ctofdof[dnums[i]];
          else
            ctypes[i] = UNUSED_DOF;
      }
  }

  void FESpace :: CheckCouplingTypes() const
  {
    cout << "checking coupling-types, type = " << typeid(*this).name() << endl;
    int ndof = GetNDof();
    if (ctofdof.Size() != ndof) 
      cout << "ndof = " << ndof
           << ", but couplingtype.size = " << ctofdof.Size() << endl;

    Array<int> cnt(ndof);
    cnt = 0;

    Array<int> dnums;
    for (ElementId id : ma->Elements<VOL>())
      {
        GetDofNrs(id, dnums);
        for (auto d : dnums) cnt[d]++;
      }
    for (int i : IntRange(0,ndof))
      {
        if (cnt[i] == 0 && ctofdof[i] != UNUSED_DOF)
          cout << "dof " << i << " not used, but coupling-type = " << ctofdof[i] << endl;
      }


    cout << "check dofs" << endl;
    for (ElementId id : ma->Elements<VOL>())
      {
        GetDofNrs (id, dnums);
        for (auto d : dnums)
          if (d < 0 || d >= ndof)
            cout << "dof out of range: " << d << endl;
      }
    for (ElementId id : ma->Elements<BND>())
      {
        GetDofNrs (id, dnums);
        for (auto d : dnums)
          if (d < 0 || d >= ndof)
            cout << "dof out of range: " << d << endl;
      }
  }


  void FESpace :: GetDofNrs (int elnr, Array<int> & dnums, COUPLING_TYPE ctype) const
  {
    ArrayMem<int,100> alldnums; 
    GetDofNrs(elnr, alldnums);

    dnums.SetSize(0);
    if (ctofdof.Size() == 0)
      {
	if ( (INTERFACE_DOF & ctype) != 0)
          dnums = alldnums;
      }
    else
      {
        /*
	for (int i = 0; i < alldnums.Size(); i++)
	  if ( (ctofdof[alldnums[i]] & ctype) != 0)
	    dnums.Append(alldnums[i]);
        */
        for (auto d : alldnums)
	  if ( (d != -1) && ((ctofdof[d] & ctype) != 0) )
            dnums.Append(d);
      }
  }

  void FESpace :: GetNodeDofNrs (NODE_TYPE nt, int nr, Array<int> & dnums) const
  {
    switch (nt)
      {
      case NT_VERTEX: GetVertexDofNrs(nr, dnums); break;
      case NT_EDGE:   GetEdgeDofNrs(nr, dnums); break;
      case NT_FACE:   
        if (ma->GetDimension() == 3)
          GetFaceDofNrs(nr, dnums); 
        else
          GetInnerDofNrs(nr, dnums); 
        break;
      case NT_CELL:   GetInnerDofNrs(nr, dnums); break;
      }
  }

  void FESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize0 ();
  }

  void FESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    dnums.SetSize0 ();
  }

  void FESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize0 ();
  }

  void FESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize0 ();
  }





  const FiniteElement & FESpace :: GetSFE (int selnr, LocalHeap & lh) const
  {
    switch (ma->GetSElType(selnr))
      {
      case ET_TRIG:  return *trig; 
      case ET_QUAD:  return *quad; 
      case ET_SEGM:  return *segm; 
      case ET_POINT: return *point; 
      case ET_TET: case ET_PYRAMID:
      case ET_PRISM: case ET_HEX: 
        ;
      }
    throw Exception ("GetSFE: unknown type");
  }

  const FiniteElement & FESpace :: GetFE (ELEMENT_TYPE type) const
  {
    switch (type)
      {
      case ET_SEGM:  return *segm;
      case ET_TRIG:  return *trig;
      case ET_QUAD:  return *quad;
      case ET_TET:   return *tet;
      case ET_PYRAMID: return *pyramid;
      case ET_PRISM: return *prism;
      case ET_HEX:   return *hex;
      case ET_POINT: return *point;
      }
    throw Exception ("GetFE: unknown type");
  }

  
  void FESpace :: PrintReport (ostream & ost) const
  {
    ost << "type  = " << GetClassName() << endl
	<< "order = " << order << endl
	<< "dim   = " << dimension << endl
	<< "dgjmps= " << dgjumps << endl
	<< "complex = " << iscomplex << endl;

    if (!free_dofs.Size()) return;

    ost << "ndof = " << GetNDof() << endl;
    int ntype[8] = { 0 };
    // for (int i = 0; i < ctofdof.Size(); i++)
    // ntype[ctofdof[i]]++;
    for (auto ct : ctofdof) ntype[ct]++;
    if (ntype[UNUSED_DOF]) ost << "unused = " << ntype[UNUSED_DOF] << endl;
    if (ntype[LOCAL_DOF])  ost << "local  = " << ntype[LOCAL_DOF] << endl;

    int nfree = 0;
    for (int i = 0; i < free_dofs.Size(); i++)
      if (free_dofs[i])
	nfree++;
  }
  
  void FESpace :: DoArchive (Archive & archive)
  {
    archive & order & dimension & iscomplex & dgjumps & print & level_updated;
    archive & definedon & definedonbound;
    archive & dirichlet_boundaries & dirichlet_dofs & free_dofs & external_free_dofs;
    archive & dirichlet_vertex & dirichlet_edge & dirichlet_face;
  }


  int FESpace :: GetNDofLevel (int level) const
  {
    return GetNDof();
  } 

  void FESpace :: Timing () const
  {
    double starttime;
    double time;
    int steps;
    LocalHeap lh (100000, "FESpace - Timing");

    cout << endl << "timing fespace " << GetName() 
         << (low_order_space ? "" : " low-order")
         << " ..." << endl;
    
    starttime = WallTime();
    steps = 0;
    do
      {
#pragma omp parallel
        {
	  LocalHeap &clh = lh, lh = clh.Split();
          Array<int> dnums;
#pragma omp for
	  for (int i = 0; i < ma->GetNE(); i++)
            GetDofNrs (i, dnums);
	}
	steps++;
	time = WallTime()-starttime;
      }
    while (time < 2.0);
    
    cout << 1e9*time / (ma->GetNE()*steps) << " ns per GetDofNrs (parallel)" << endl;

    /*
    starttime = WallTime();
    steps = 0;
    do
      {
#pragma omp parallel
        {
	  LocalHeap &clh = lh, lh = clh.Split();
#pragma omp for
	  for (int i = 0; i < ma->GetNE(); i++)
	    {
              HeapReset hr(lh);
              // FlatArray<int> dnums = 
              GetDofNrs (ElementId (VOL, i), lh);
	    }
	}
	steps++;
	time = WallTime()-starttime;
      }
    while (time < 2.0);
    
    cout << 1e9*time / (ma->GetNE()*steps) << " ns per GetDofNrs(lh) (parallel)" << endl;
    */




    starttime = WallTime();
    steps = 0;
    do
      {
#pragma omp parallel
        {
	  LocalHeap &clh = lh, lh = clh.Split();

#pragma omp for
	  for (int i = 0; i < ma->GetNE(); i++)
	    {
	      HeapReset hr(lh);
	      GetFE (i, lh);
	    }
	}
        steps++;
        time = WallTime()-starttime;
      }
    while (time < 2.0);
    
    cout << 1e9 * time / (ma->GetNE()*steps) << " ns per GetFE (parallel)" << endl;




    starttime = WallTime();
    steps = 0;
    do
      {
#pragma omp parallel for
        for (int i = 0; i < ma->GetNE(); i++)
          {
	    /* Ng_Element ngel = */ ma->GetElement(i);
          }
        steps++;
        time = WallTime()-starttime;
      }
    while (time < 2.0);
    
    cout << 1e9 * time / (ma->GetNE()*steps) << " ns per Get - Ng_Element (parallel)" << endl;


    starttime = WallTime();
    steps = 0;
    do
      {
#pragma omp parallel
        {
	  LocalHeap &clh = lh, lh = clh.Split();
#pragma omp for
          for (int i = 0; i < ma->GetNE(); i++)
            {
              HeapReset hr(lh);
              /* ElementTransformation & trafo = */ ma->GetTrafo(i, VOL, lh);
            }
        }
        steps++;
        time = WallTime()-starttime;
      }
    while (time < 2.0);
    
    cout << 1e9 * time / (ma->GetNE()*steps) << " ns per GetTrafo(parallel)" << endl;




#ifdef TIMINGSEQUENTIAL


    starttime = WallTime();
    steps = 0;
    do
      {
        for (int i = 0; i < ma->GetNE(); i++)
	  {
	    ArrayMem<int,100> dnums;
	    GetDofNrs (i, dnums);
	  }
        steps++;
        time = WallTime()-starttime;
      }
    while (time < 2.0);
    
    cout << 1e9*time / (ma->GetNE()*steps) << " ns per GetDofNrs" << endl;

    // ***************************

    starttime = WallTime();
    steps = 0;
    do
      {
        for (int i = 0; i < ma->GetNE(); i++)
	  {
            HeapReset hr(lh);
	    /* FlatArray<int> dnums = */ GetDofNrs (ElementId(VOL, i), lh);
	  }
        steps++;
        time = WallTime()-starttime;
      }
    while (time < 2.0);
    
    cout << 1e9*time / (ma->GetNE()*steps) << " ns per GetDofNrs(lh)" << endl;

    // ***************************


    starttime = WallTime();
    steps = 0;
    do
      {
        for (int i = 0; i < ma->GetNE(); i++)
          {
            HeapReset hr(lh);
            GetFE (i, lh);
          }
        steps++;
        time = WallTime()-starttime;
      }
    while (time < 2.0);
    
    cout << 1e9 * time / (ma->GetNE()*steps) << " ns per GetFE" << endl;


    starttime = WallTime();
    steps = 0;
    do
      {
        for (int i = 0; i < ma->GetNE(); i++)
          {
	    /* Ng_Element ngel = */ ma->GetElement(i);
          }
        steps++;
        time = WallTime()-starttime;
      }
    while (time < 2.0);
    
    cout << 1e9 * time / (ma->GetNE()*steps) << " ns per Get - Ng_Element" << endl;


    starttime = WallTime();
    steps = 0;
    do
      {
        for (int i = 0; i < ma->GetNE(); i++)
          {
            HeapReset hr(lh);
            /* ElementTransformation & trafo = */ ma->GetTrafo(i, VOL, lh);
          }
        steps++;
        time = WallTime()-starttime;
      }
    while (time < 2.0);
    
    cout << 1e9 * time / (ma->GetNE()*steps) << " ns per GetTrafo" << endl;

#endif


  }



  void FESpace :: GetFilteredDofs (COUPLING_TYPE doffilter, BitArray & output, 
                                   bool freedofsonly) const
  {
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
      definedon[i] = defon.Test(i);

    if (low_order_space)
      low_order_space -> SetDefinedOn (defon);
  }

  void FESpace :: SetDefinedOnBoundary (const BitArray & defon)
  {
    definedonbound.SetSize(defon.Size());
    for (int i = 0; i < defon.Size(); i++)
      definedonbound[i] = defon.Test(i);

    if (low_order_space)
      low_order_space -> SetDefinedOnBoundary (defon);
  }

  void FESpace :: SetDirichletBoundaries (const BitArray & dirbnds)
  {
    dirichlet_boundaries = dirbnds;
    if (low_order_space)
      low_order_space -> SetDirichletBoundaries (dirbnds);
  }


  const BitArray * FESpace :: GetFreeDofs (bool external) const
  {
    if (external)
      return &external_free_dofs;
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


  ostream & operator<< (ostream & ost, COUPLING_TYPE ct)
  {
    switch (ct)
      {
      case UNUSED_DOF: ost << "unused"; break;
      case LOCAL_DOF:  ost << "local"; break;
      case INTERFACE_DOF: ost << "interface"; break;
      case NONWIREBASKET_DOF: ost << "non-wirebasket"; break;
      case WIREBASKET_DOF: ost << "wirebasket"; break;
      case EXTERNAL_DOF: ost << "external"; break;
      case ANY_DOF: ost << "any"; break;
      };
    return ost;
  }



  void FESpace :: UpdateParallelDofs ( )
  {
    if (MyMPI_GetNTasks() == 1) return;

    Array<Node> dofnodes (GetNDof());
    dofnodes = Node (NT_VERTEX, -1);
    Array<int> dnums;

    for (NODE_TYPE nt = NT_VERTEX; nt <= NT_CELL; nt++)
      for ( int nr = 0; nr < ma->GetNNodes (nt); nr++ )
	{
	  GetNodeDofNrs (nt, nr, dnums);
	  for (int d : dnums)
	    dofnodes[d] = Node (nt, nr);
	} 

    paralleldofs = new ParallelMeshDofs (ma, dofnodes, dimension, iscomplex);

    if (MyMPI_AllReduce (ctofdof.Size(), MPI_SUM))
      AllReduceDofData (ctofdof, MPI_MAX, GetParallelDofs());
  }


  bool FESpace :: IsParallel() const
  { 
    return paralleldofs != NULL; 
  }

  int FESpace :: GetNDofGlobal() const 
  { 
    return paralleldofs ?
      paralleldofs -> GetNDofGlobal() : GetNDof(); 
  }





  NodalFESpace :: NodalFESpace (shared_ptr<MeshAccess> ama,
				const Flags & flags,
                                bool parseflags)
    : FESpace (ama, flags)
  {
    name="NodalFESpace";
    
    prol = make_shared<LinearProlongation> (*this);

    if (order >= 2)
      {
	Flags loflags;
	loflags.SetFlag ("order", 1);
	loflags.SetFlag ("dim", dimension);
	if (dgjumps) loflags.SetFlag ("dgjumps");
	if (iscomplex) loflags.SetFlag ("complex");
	low_order_space = make_shared<NodalFESpace> (ma, loflags);
      }

    if (order == 1)
      {
	tet     = new ScalarFE<ET_TET,1>;
	prism   = new FE_Prism1;
	pyramid = new FE_Pyramid1;
	hex     = new FE_Hex1;
	trig    = new ScalarFE<ET_TRIG,1>;
	quad    = new ScalarFE<ET_QUAD,1>;
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
	    quad    = new ScalarFE<ET_QUAD,1>;
	    segm    = new FE_Segm2;
	  }
	else
	  {
	    tet     = new FE_Tet2;
	    prism   = new FE_Prism1;
	    pyramid = new FE_Pyramid1;
	    trig    = new FE_Trig2;
	    quad    = new ScalarFE<ET_QUAD,1>;
	    segm    = new FE_Segm2;
	  }
      }
    point = new FE_Point;

    SetDummyFE<ScalarDummyFE> ();

    auto one = make_shared<ConstantCoefficientFunction> (1);
    if (ma->GetDimension() == 2)
      {
	integrator = make_shared<MassIntegrator<2>> (one);
        boundary_integrator = make_shared<RobinIntegrator<2>> (one);
      }
    else
      {
	integrator = make_shared<MassIntegrator<3>> (one);
	boundary_integrator = make_shared<RobinIntegrator<3>> (one);
      }
    
    if (dimension > 1)
      {
	integrator = make_shared<BlockBilinearFormIntegrator> (integrator, dimension);
	boundary_integrator = make_shared<BlockBilinearFormIntegrator> (boundary_integrator, dimension);  
      }
  }

  NodalFESpace :: ~NodalFESpace ()
  {
    ;
  }

  int NodalFESpace :: GetNDof () const throw()
  {
    return ndlevel.Last();
  }

  void NodalFESpace :: Update(LocalHeap & lh)
  {
    FESpace :: Update (lh);
    if (low_order_space) low_order_space -> Update(lh);

    if (ma->GetNLevels() > ndlevel.Size())
      {
	int ndof = ma->GetNV();

        for (auto el : Elements(VOL))
          for (int d : el.GetDofs()) ndof = max2(ndof, d+1);              

        for (auto el : Elements(BND))
          for (int d : el.GetDofs()) ndof = max2(ndof, d+1);           

	ndlevel.Append (ndof);
      }

    prol->Update();

    if (dirichlet_boundaries.Size())
      {
	dirichlet_dofs.SetSize (GetNDof());
	dirichlet_dofs.Clear();
	for (auto el : Elements(BND))
          if (dirichlet_boundaries[el.GetIndex()])
            for (int d : el.GetDofs())
              if (d != -1) dirichlet_dofs.Set (d);
      }
  }

  void NodalFESpace :: DoArchive (Archive & archive)
  {
    if (archive.Input())
      {
	ndlevel.SetSize(1);
	ndlevel[0] = ma->GetNV();
      }
  }

  int NodalFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }


  void NodalFESpace :: GetDofRanges (ElementId ei, Array<IntRange> & dranges) const
  {
    if (!DefinedOn (ei)) 
      {
        dranges.SetSize(0);
        return;
      }

    ArrayMem<int,27> pnums;

    if (order == 1)
      ma->GetElVertices(ei, pnums);
    else
      ma->GetElPNums (ei, pnums);

    dranges.SetSize(pnums.Size());
    for (int i = 0; i < pnums.Size(); i++)
      dranges[i] = IntRange (pnums[i], pnums[i]+1);
  }



 
  void NodalFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    if (order == 1)
      ma->GetElVertices (elnr, dnums);
    else
      ma->GetElPNums (elnr, dnums);

    if (!DefinedOn (ElementId (VOL, elnr))) dnums = -1;
  }


  void NodalFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    ma->GetSElPNums (selnr, dnums);

    if (order == 1)
      { // Ng-mesh may be second order, but FE space is 1st order
	int np = dnums.Size();
	switch (ma->GetSElType(selnr))
	  {
	  case ET_SEGM: np = 2; break;
	  case ET_TRIG: np = 3; break;
	  case ET_QUAD: np = 4; break;
          default:
            ;
	  }
	if (dnums.Size() > np) dnums.SetSize (np);
      }

    if (!DefinedOnBoundary (ma->GetSElIndex (selnr)))
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
  NonconformingFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name="NonconformingFESpace(nonconforming)";
    // defined flags
    DefineDefineFlag("nonconforming");
    if (parseflags) CheckFlags(flags);
    
    // prol = new LinearProlongation(*this);
    

    trig = new FE_NcTrig1;

    if (ma->GetDimension() == 2)
      {
	integrator.reset (new MassIntegrator<2> (new ConstantCoefficientFunction(1)));
	boundary_integrator = 0;
      }
    else
      {
	integrator.reset (new MassIntegrator<3> (new ConstantCoefficientFunction(1)));
	boundary_integrator.reset (new RobinIntegrator<3> (new ConstantCoefficientFunction(1)));
      }

    if (dimension > 1)
      {
	integrator = make_shared<BlockBilinearFormIntegrator> (integrator, dimension);
	boundary_integrator = make_shared<BlockBilinearFormIntegrator> (boundary_integrator, dimension);
      }
  }

  NonconformingFESpace :: ~NonconformingFESpace ()
  {
    ;
  }


  int NonconformingFESpace :: GetNDof () const throw()
  {
    return ma->GetNEdges();
  }


  void NonconformingFESpace :: Update(LocalHeap & lh)
  {
    /*
    if (ma->GetNLevels() > ndlevel.Size())
      {
	Array<int> dnums;
	int i, j;
	int ne = ma->GetNE();
	int nse = ma->GetNSE();
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
	for (int i = 0; i < ma->GetNSE(); i++)
	  {
	    if (dirichlet_boundaries[ma->GetSElIndex(i)])
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
    ma->GetElEdges (elnr, dnums);
    if (!DefinedOn (ma->GetElIndex (elnr)))
      dnums = -1;
  }


  void NonconformingFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    ma->GetSElEdges (selnr, dnums);
    if (!DefinedOnBoundary (ma->GetSElIndex (selnr)))
      dnums = -1;
  }
  








  ElementFESpace :: ElementFESpace (shared_ptr<MeshAccess> ama, const Flags& flags, 
                                    bool parseflags)
    : FESpace (ama, flags)
  {
    name="ElementFESpace(l2)";
    if (parseflags) CheckFlags(flags);
    
    order = int(flags.GetNumFlag ("order", 0));

    prol = make_shared<ElementProlongation> (*this);

    if (order == 0)
    {
      tet     = new ScalarFE<ET_TET,0>;
      prism   = new FE_Prism0;
      pyramid = new FE_Pyramid0;
      hex     = new FE_Hex0;
      trig    = new ScalarFE<ET_TRIG,0>;
      quad    = new ScalarFE<ET_QUAD,0>;
      segm    = new FE_Segm0;

      n_el_dofs = 1;
    }
    else
    {
      tet     = new ScalarFE<ET_TET,1>;
      prism   = new FE_Prism1;
      pyramid = new FE_Pyramid1;
      trig    = new ScalarFE<ET_TRIG,1>;
      quad    = new ScalarFE<ET_QUAD,1>;
      segm    = new FE_Segm1;

      if (ma->GetDimension() == 2)
        n_el_dofs = 4;
      else
        n_el_dofs = 6;
    }

    SetDummyFE<ScalarDummyFE> ();
    static ConstantCoefficientFunction one(1);

    if (ma->GetDimension() == 2)
      {
        integrator.reset (new MassIntegrator<2> (&one));
        boundary_integrator = 0;
      }
    else
      {
        integrator.reset (new MassIntegrator<3> (&one));
        boundary_integrator = 0;
      }
    
    if (dimension > 1)
      integrator = make_shared<BlockBilinearFormIntegrator> (integrator, dimension);
  }
  
  ElementFESpace :: ~ElementFESpace ()
  {
    ;
  }

  void ElementFESpace :: Update(LocalHeap & lh)
  {
    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (n_el_dofs * ma->GetNE());
  }

  void ElementFESpace :: DoArchive (Archive & archive)
  {
    FESpace :: DoArchive (archive);
    archive & ndlevel & n_el_dofs;
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
	switch (ma->GetElType(elnr))
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
      SurfaceElementFESpace (shared_ptr<MeshAccess> ama, const Flags& flags, bool parseflags)
  : FESpace (ama, flags)
  {
    name="SurfaceElementFESpace(surfl2)";
    if(parseflags) CheckFlags(flags);
    
    // prol = new SurfaceElementProlongation (GetMeshAccess(), *this);

    if (order == 0)
    {
      trig    = new ScalarFE<ET_TRIG,0>;
      quad    = new ScalarFE<ET_QUAD,0>;
      segm    = new FE_Segm0;

      n_el_dofs = 1;
    }

    else if (order == 1)
    {
      trig    = new ScalarFE<ET_TRIG,1>;
      quad    = new ScalarFE<ET_QUAD,1>;
      segm    = new FE_Segm1;
	
      if (ma->GetDimension() == 2)
        n_el_dofs = 2;
      else
        n_el_dofs = 4;
    }

    else if (order == 2)
    {
      trig    = new FE_Trig2HB;
      quad    = new ScalarFE<ET_QUAD,1>;
      segm    = new FE_Segm2;

      if (ma->GetDimension() == 2)
        n_el_dofs = 3;
      else
        n_el_dofs = 9;
    }

    boundary_integrator.reset (new RobinIntegrator<3> (new ConstantCoefficientFunction(1)));

    if (dimension > 1)
      boundary_integrator = make_shared<BlockBilinearFormIntegrator> (boundary_integrator, dimension);
  }

  
  SurfaceElementFESpace :: ~SurfaceElementFESpace ()
  {
    ;
  }

  void SurfaceElementFESpace :: Update(LocalHeap & lh)
  {
    // const MeshAccess & ma = GetMeshAccess();

    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (n_el_dofs * ma->GetNSE());
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
	switch (ma->GetSElType(elnr))
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
	switch (ma->GetSElType(elnr))
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





 











  CompoundFESpace :: CompoundFESpace (shared_ptr<MeshAccess> ama,
				      const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name="CompoundFESpaces";
    DefineDefineFlag("compound");
    DefineStringListFlag("spaces");
    if (parseflags) CheckFlags(flags);
    
    prol = make_shared<CompoundProlongation> (this);
  }



  CompoundFESpace :: CompoundFESpace (shared_ptr<MeshAccess> ama,
				      const Array<shared_ptr<FESpace>> & aspaces,
				      const Flags & flags, bool parseflags)
    : FESpace (ama, flags), spaces(aspaces)
  {
    name="CompoundFESpaces";
    DefineDefineFlag("compound");
    DefineStringListFlag("spaces");
    if(parseflags) CheckFlags(flags);
    
    auto hprol = make_shared<CompoundProlongation> (this);
    /*
    for (int i = 0; i < spaces.Size(); i++)
      hprol -> AddProlongation (spaces[i]->GetProlongation());
    */
    for (auto space : spaces)
      hprol -> AddProlongation (space->GetProlongation());      
    prol = hprol;
  }


  void CompoundFESpace :: AddSpace (shared_ptr<FESpace> fes)
  {
    spaces.Append (fes);
    dynamic_cast<CompoundProlongation*> (prol.get()) -> AddProlongation (fes->GetProlongation());
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
	spaces[i] -> Update(lh);
	cummulative_nd[i+1] = cummulative_nd[i] + spaces[i]->GetNDof();
      }

    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (cummulative_nd.Last());


    
    free_dofs.SetSize (GetNDof());
    free_dofs.Clear();
    for (int i = 0; i < spaces.Size(); i++)
      {
	const BitArray & freedofsi = *spaces[i]->GetFreeDofs(false);
	for (int j = 0; j < freedofsi.Size();j++)
	  if (freedofsi.Test(j)) 
	    free_dofs.Set(cummulative_nd[i]+j);
      }
    external_free_dofs.SetSize (GetNDof());
    external_free_dofs.Clear();
    for (int i = 0; i < spaces.Size(); i++)
      {
	const BitArray & freedofsi = *spaces[i]->GetFreeDofs(true);
	for (int j = 0; j < freedofsi.Size();j++)
	  if (freedofsi.Test(j)) 
	    external_free_dofs.Set(cummulative_nd[i]+j);
      }
    
    



    prol -> Update();

    UpdateCouplingDofArray();

    if (print)
      {
	(*testout) << "Update compound fespace" << endl;
	(*testout) << "cummulative dofs start at " << cummulative_nd << endl;
      }
  }

  void CompoundFESpace :: FinalizeUpdate(LocalHeap & lh)
  {
    for (int i = 0; i < spaces.Size(); i++)
      spaces[i] -> FinalizeUpdate(lh);

    FESpace::FinalizeUpdate (lh);


    // dirichlet-dofs from sub-spaces
    // ist das umsonst ? (JS)
    // haben jetzt ja immer dirichlet-dofs
    bool has_dirichlet_dofs = false;
    for (int i = 0; i < spaces.Size(); i++)
      if (spaces[i]->GetFreeDofs()) 
	has_dirichlet_dofs = true;

    has_dirichlet_dofs = MyMPI_AllReduce (has_dirichlet_dofs, MPI_LOR);

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

        for (int i = 0; i < ctofdof.Size(); i++)
          if (ctofdof[i] == UNUSED_DOF)
            free_dofs.Clear(i);

	dirichlet_dofs = free_dofs;
	dirichlet_dofs.Invert();

        external_free_dofs.SetSize (GetNDof());
        external_free_dofs = free_dofs;
        for (int i = 0; i < ctofdof.Size(); i++)
          if (ctofdof[i] & LOCAL_DOF)
            external_free_dofs.Clear(i);


        if (print)
          (*testout) << "compound fespace freedofs:" << endl
                     << free_dofs << endl;

      }
  }



  void CompoundFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(this->GetNDof());

    for (int i = 0; i < spaces.Size(); i++)
      for (int j=0; j< spaces[i]->GetNDof();j++)
	ctofdof[cummulative_nd[i]+j] = spaces[i]->GetDofCouplingType(j);	

    // *testout << "CompoundFESpace :: UpdateCouplingDofArray() presents \n ctofdof = \n" << ctofdof << endl;
  }



  const FiniteElement & CompoundFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    FlatArray<const FiniteElement*> fea(spaces.Size(), lh);
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
  

  void CompoundFESpace :: GetDofRanges (ElementId ei, Array<IntRange> & dranges) const
  {
  /*
    dranges.SetSize (0);
    for (int i = 0; i < spaces.Size(); i++)
      {
        int osize = dranges.Size();

        Array<IntRange> hdranges(dranges.AllocSize()-osize, &dranges[osize]);
        hdranges.SetSize(0);
	spaces[i]->GetDofRanges (ei, hdranges);
        int nsize = osize + hdranges.Size();

        dranges.SetSize (nsize);
        for (int j = 0; j < hdranges.Size(); j++)
          if (hdranges[j].First() != -1)
            dranges[osize+j] = hdranges[j]+cummulative_nd[i];
          else
            dranges[osize+j] = hdranges[j];
      }
  */
  }




  const FiniteElement & CompoundFESpace :: GetSFE (int elnr, LocalHeap & lh) const
  {
    FlatArray<const FiniteElement*> fea(spaces.Size(), lh);
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



  template NGS_DLL_HEADER
  void CompoundFESpace::TransformVec<FlatVector<double> >
  (int elnr, bool boundary, FlatVector<double> & vec, TRANSFORM_TYPE tt) const;
  template NGS_DLL_HEADER
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






  Table<int> * Nodes2Table (const MeshAccess & ma,
			    const Array<Node> & dofnodes)
  {
    int ndof = dofnodes.Size();

    Array<int> distprocs;
    Array<int> ndistprocs(ndof);
    ndistprocs = 0;
    for (int i = 0; i < ndof; i++)
      {
	if (dofnodes[i].GetNr() == -1) continue;
	ma.GetDistantProcs (dofnodes[i], distprocs);
	ndistprocs[i] = distprocs.Size();
      }

    Table<int> * dist_procs = new Table<int> (ndistprocs);

    for (int i = 0; i < ndof; i++)
      {
	if (dofnodes[i].GetNr() == -1) continue;
	ma.GetDistantProcs (dofnodes[i], distprocs);
	(*dist_procs)[i] = distprocs;
      }

    return dist_procs;
  }


#ifdef PARALLEL
  ParallelMeshDofs :: ParallelMeshDofs (shared_ptr<MeshAccess> ama, 
					const Array<Node> & adofnodes, 
					int dim, bool iscomplex)
    : ParallelDofs (ma->GetCommunicator(),
		    Nodes2Table (*ama, adofnodes), dim, iscomplex),		    
      ma(ama), dofnodes(adofnodes)
  { ; }
#endif







  FESpaceClasses :: ~FESpaceClasses() { ; }
  
  void FESpaceClasses :: 
  AddFESpace (const string & aname,
	      shared_ptr<FESpace> (*acreator)(shared_ptr<MeshAccess> ma, const Flags & flags))
  {
    fesa.Append (make_shared<FESpaceInfo> (aname, acreator));
  }

  const shared_ptr<FESpaceClasses::FESpaceInfo> 
  FESpaceClasses::GetFESpace(const string & name)
  {
    for (auto & fes : fesa)
      if (name == fes->name) return fes;
    return NULL;
  }

  void FESpaceClasses :: Print (ostream & ost) const
  {
    ost << endl << "FESpaces:" << endl;
    ost <<         "---------" << endl;
    ost << setw(20) << "Name" << endl;
    for (auto & fes : fesa)
      ost << setw(20) << fes->name << endl;
  }

 
  FESpaceClasses & GetFESpaceClasses ()
  {
    static FESpaceClasses fecl;
    return fecl;
  }

  extern NGS_DLL_HEADER shared_ptr<FESpace> CreateFESpace (const string & type,
                                                           shared_ptr<MeshAccess> ma,
                                                           const Flags & flags)
  {
    shared_ptr<FESpace> space;
    for (int i = 0; i < GetFESpaceClasses().GetFESpaces().Size(); i++)
      if (type == GetFESpaceClasses().GetFESpaces()[i]->name ||
	  flags.GetDefineFlag (GetFESpaceClasses().GetFESpaces()[i]->name) )
	{
	  space = GetFESpaceClasses().GetFESpaces()[i]->creator (ma, flags);
          space -> type = type;
	}
    if (!space)
      throw Exception (string ("undefined fespace '") + type + '\'');
    return space;
  }



  // standard fespaces:

  static RegisterFESpace<NodalFESpace> initnodalfes ("nodal");
  static RegisterFESpace<NonconformingFESpace> initncfes ("nonconforming");

}



