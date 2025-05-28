/**********************************************************************/
/* File:   fespace.cpp                                                */
/* Author: Joachim Schoeberl                                          */
/* Date:   25. Mar. 2000                                              */
/**********************************************************************/


/* 
   Finite Element Space
*/

// #include <comp.hpp>
#include "fespace.hpp"
#include <multigrid.hpp>
#include <diagonalmatrix.hpp>
#include <special_matrix.hpp>

#include "../fem/h1lofe.hpp"
#include <bdbequations.hpp>
#include <parallelngs.hpp>
#include "compressedfespace.hpp"
#include "periodic.hpp"


namespace ngcomp
{
  using namespace ngmg;

  FESpace :: FESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags)
    : NGS_Object (ama, flags, "FESpace")
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
    DefineDefineFlag("autoupdate");

    order = int (flags.GetNumFlag ("order", 1));

    if (flags.NumFlagDefined("order_policy"))
      SetOrderPolicy(ORDER_POLICY(int(flags.GetNumFlag("order_policy",1))));

    /*
    if (flags.NumFlagDefined("order_left"))
      {
        auto order_left = int(flags.GetNumFlag("order_left", 1));
        order = max(order, order_left);
        for (auto et : element_types)
          SetOrderLeft (et, order_left);
      }
    if (flags.NumFlagDefined("order_right"))
      {
        auto order_right = int(flags.GetNumFlag("order_right", 1));
        order = max(order, order_right);    
        for (auto et : element_types)
          SetOrderRight (et, order_right);
      }
    */
    
    dimension = int (flags.GetNumFlag ("dim", 1));

    if (flags.GetDefineFlag ("vec"))
      dimension = ma->GetDimension();
    if (flags.GetDefineFlag ("tensor")) 
      dimension = sqr (ma->GetDimension());
    if (flags.GetDefineFlag ("symtensor")) 
      dimension = ma->GetDimension()*(ma->GetDimension()+1) / 2;
    
    iscomplex = flags.GetDefineFlag ("complex");
//     eliminate_internal = flags.GetDefineFlag("eliminate_internal");
    timing = flags.GetDefineFlag("timing");
    print = flags.GetDefineFlag("print");
    dgjumps = flags.GetDefineFlag("dgjumps");
    
    autoupdate = flags.GetDefineFlagX("autoupdate").IsMaybeTrue();
    ConnectAutoUpdate();

    
    no_low_order_space = flags.GetDefineFlagX("low_order_space").IsFalse() ||
      flags.GetDefineFlag("no_low_order_space");
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
    
    dirichlet_constraints[BND].SetSize (ma->GetNBoundaries());
    dirichlet_constraints[BND].Clear();

    if(flags.AnyFlagDefined("dirichlet"))
      {
        dirichlet_constraints[BND] = std::any_cast<const Region&>(flags.GetAnyFlag("dirichlet")).Mask();
      }

    if(flags.NumListFlagDefined("dirichlet"))
      {
        for (double dbi : flags.GetNumListFlag("dirichlet"))
          {
	    int bnd = int(dbi-1);
	    if (bnd >= 0 && bnd < dirichlet_boundaries.Size())
	      dirichlet_constraints[BND].SetBit (bnd);
	    // else
            //   cerr << "Illegal Dirichlet boundary index " << bnd+1 << endl;
          }
	if (print)
	  *testout << "dirichlet_boundaries:" << endl << dirichlet_boundaries << endl;
      }

    if (flags.StringFlagDefined("dirichlet"))
      {
        Region dir(ma, BND, flags.GetStringFlag("dirichlet"));
        dirichlet_constraints[BND] |= dir.Mask(); // that's what it was before, can we just assign ? 
      }

    dirichlet_constraints[BBND].SetSize (ma->GetNRegions(BBND));
    dirichlet_constraints[BBND].Clear();

    if(flags.AnyFlagDefined("dirichlet_bbnd"))
      {
        dirichlet_constraints[BBND] = std::any_cast<const Region&>(flags.GetAnyFlag("dirichlet_bbnd")).Mask();
      }

    if (flags.StringFlagDefined("dirichlet_bbnd"))
      {
        Region dir(ma, BBND, flags.GetStringFlag("dirichlet_bbnd"));
        dirichlet_constraints[BBND] |= dir.Mask(); // that's what it was before, can we just assign ? 
      }

    dirichlet_constraints[BBBND].SetSize (ma->GetNRegions(BBBND));
    dirichlet_constraints[BBBND].Clear();

    if(flags.AnyFlagDefined("dirichlet_bbbnd"))
      {
        dirichlet_constraints[BBBND] = std::any_cast<const Region&>(flags.GetAnyFlag("dirichlet_bbbnd")).Mask();
      }

    if (flags.StringFlagDefined("dirichlet_bbbnd"))
      {
        Region dir(ma, BBBND, flags.GetStringFlag("dirichlet_bbbnd"));
        dirichlet_constraints[BBBND] |= dir.Mask(); // that's what it was before, can we just assign ? 
      }
    
    if (flags.NumListFlagDefined("definedon") || 
        flags.NumFlagDefined("definedon") ||
        flags.StringListFlagDefined("definedon") ||
        flags.AnyFlagDefined("definedon"))
      {
	definedon[VOL].SetSize (ma->GetNDomains());
	definedon[VOL] = false;
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
            definedon[VOL][di-1] = true;
          
	if(flags.StringListFlagDefined("definedon"))
	  {
	    Array<string> dmaterials(flags.GetStringListFlag ("definedon").Size());
	    for(int i=0; i<dmaterials.Size(); i++)
	      dmaterials[i] = flags.GetStringListFlag ("definedon")[i];
	    for(int i = 0; i < ma->GetNDomains(); i++)
	      {
		for(int j = 0; j < dmaterials.Size(); j++)
		  if(StringFitsPattern(ma->GetMaterial(VOL,i),dmaterials[j]))
		    {
		      definedon[VOL][i] = true;
		      break;
		    }
	      }
	  }

	// default:
	// fespace only defined on boundaries matching definedon-domains
	definedon[BND].SetSize (ma->GetNBoundaries());
	definedon[BND] = false;

        if(flags.AnyFlagDefined("definedon"))
          {
            std::any anydefon = flags.GetAnyFlag("definedon");
            if(anydefon.type() == typeid(std::map<VorB, Region>))
              {
                auto map = std::any_cast<std::map<VorB, Region>>(anydefon);
                for(const auto& [vb, reg] : map)
                  {
                    if(vb == BBND || vb == BBBND)
                      definedon[vb].SetSize(ma->GetNRegions(vb));
                    for(auto i : Range(definedon[vb]))
                      definedon[vb][i] = reg.Mask().Test(i);
                  }
                // if vol is not in map then definedon vol = empty
                if(map.find(VOL) == map.end())
                  definedon[VOL] = false;
              }
            else if(anydefon.type() == typeid(Region))
              {
                auto reg = std::any_cast<Region>(anydefon);
                definedon[reg.VB()].SetSize(reg.Mask().Size());
                for(auto i : Range(definedon[reg.VB()]))
                  definedon[reg.VB()][i] = reg.Mask().Test(i);
                for(auto vb = VOL; vb < reg.VB(); vb++)
                  definedon[vb] = false;
              }
            else
              throw Exception("definedon with wrong argument type!");
          }

	for (int sel = 0; sel < ma->GetNSE(); sel++)
	  {
            ElementId sei(BND, sel);
	    int index = ma->GetElIndex(sei);
	    int dom1, dom2;
	    ma->GetSElNeighbouringDomains(sel, dom1, dom2);
	    dom1--; dom2--;
	    if ( dom1 >= 0 )
	      if ( definedon[VOL][dom1] )
		definedon[BND][index] = true;

	    if ( dom2 >= 0 )
	      if ( definedon[VOL][dom2] )
		definedon[BND][index] = true;
	  }
      }

    // additional boundaries
    if(flags.NumListFlagDefined("definedonbound")|| flags.NumFlagDefined("definedonbound") )
      {
	if ( definedon[BND].Size() == 0 )
	  {
	    definedon[BND].SetSize (ma->GetNBoundaries());
	    definedon[BND] = false;
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
	    definedon[BND][int(defon[i])-1] = true;
      }
    

    else if(flags.StringListFlagDefined("definedonbound") || flags.StringFlagDefined("definedonbound"))
      {
	if ( definedon[BND].Size() == 0 )
	  {
	    definedon[BND].SetSize (ma->GetNBoundaries());
	    definedon[BND] = false;
	  }

	Array<string*> defon;

	if(flags.StringFlagDefined("definedonbound"))
	  defon.Append(new string(flags.GetStringFlag("definedonbound","")));
	else
	  for(int i=0; i<flags.GetStringListFlag ("definedonbound").Size(); i++)
	    defon.Append(new string(flags.GetStringListFlag("definedonbound")[i]));
	
	for(int selnum = 0; selnum < ma->GetNSE(); selnum++)
	  {
            ElementId sei(BND, selnum);
	    if(definedon[BND][ma->GetElIndex(sei)] == false)
	      {
		for(int i=0; i<defon.Size(); i++)
		  {
		    // if(StringFitsPattern(ma->GetSElBCName(selnum),*(defon[i])))
                    if(StringFitsPattern(ma->GetMaterial(ElementId(BND, selnum)),*(defon[i])))	
		      {		
		 	definedon[BND][ma->GetElIndex(sei)] = true;
			continue;
		      }
		  }
	      }
	  }
	for(int i=0; i<defon.Size(); i++)
	  delete defon[i];
      }
    
    level_updated = -1;

    for(auto vb : {VOL,BND,BBND})
      {
	evaluator[vb] = nullptr;
	flux_evaluator[vb] = nullptr;
	integrator[vb] = nullptr;
      }
    low_order_space = NULL;
    prol = NULL;


    // element_coloring = NULL;
    // selement_coloring = NULL;
    // paralleldofs = NULL;

    ctofdof.SetSize(0);

    for (auto & et : et_bonus_order) et = 0;
    et_bonus_order[ET_QUAD] = int (flags.GetNumFlag("quadbonus",0));


    for (auto et : element_types)
      {
        string flagname = string(ElementTopology::GetElementName(et)) + "order";
        std::transform(flagname.begin(), flagname.end(), flagname.begin(), ::tolower);
        if (flags.NumFlagDefined(flagname))
          SetOrder (et, int(flags.GetNumFlag(flagname, 0)));
      }


    
    this->GetMemoryTracer().Track(
        dirichlet_dofs, "dirichlet_dofs",
        dirichlet_vertex, "dirichlet_vertex",
        dirichlet_edge, "dirichlet_edge",
        dirichlet_face, "dirichlet_face",
        ctofdof, "ctofdof"
    );
  }

  
  FESpace :: ~FESpace ()
  {
    this->GetMeshAccess()->updateSignal.Remove(this);
  }


  DocInfo FESpace :: GetDocu ()
  {
    DocInfo docu;
    docu.Arg("order") = "int = 1\n"
      "  order of finite element space";
    docu.Arg("complex") = "bool = False\n"
      "  Set if FESpace should be complex";
    docu.Arg("dirichlet") = "regexpr\n"
      "  Regular expression string defining the dirichlet boundary.\n"
      "  More than one boundary can be combined by the | operator,\n"
      "  i.e.: dirichlet = 'top|right'";
    docu.Arg("dirichlet_bbnd") = "regexpr\n"
      "  Regular expression string defining the dirichlet bboundary,\n"
      "  i.e. points in 2D and edges in 3D.\n"
      "  More than one boundary can be combined by the | operator,\n"
      "  i.e.: dirichlet_bbnd = 'top|right'";
    docu.Arg("dirichlet_bbbnd") = "regexpr\n"
      "  Regular expression string defining the dirichlet bbboundary,\n"
      "  i.e. points in 3D.\n"
      "  More than one boundary can be combined by the | operator,\n"
      "  i.e.: dirichlet_bbbnd = 'top|right'";
    docu.Arg("definedon") = "Region or regexpr\n"
      "  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n"
      "  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n"
      "  mesh.Materials('regexpr').";
    docu.Arg("dim") = "int = 1\n"
      "  Create multi dimensional FESpace (i.e. [H1]^3)";
    docu.Arg("dgjumps") = "bool = False\n"
      "  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n"
      "  since the dofs have a different coupling then and this changes the sparsity\n"
      "  pattern of matrices.";
    docu.Arg("autoupdate") = "bool = False\n"
      "  Automatically update on a change to the mesh.";
    docu.Arg("low_order_space") = "bool = True\n"
      "  Generate a lowest order space together with the high-order space,\n"
      "  needed for some preconditioners.";
    docu.Arg("hoprolongation") 
      = "bool = False\n"
        "  Create high order prolongation operators,\n"
        "  only available for H1 and L2 on simplicial meshes";
    docu.Arg("order_policy") = "ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n"
                               "  CONSTANT .. use the same fixed order for all elements,\n"
                               "  NODAL ..... use the same order for nodes of same shape,\n"
                               "  VARIABLE ... use an individual order for each edge, face and cell,\n"
                               "  OLDSTYLE .. as it used to be for the last decade";
    docu.Arg("print") = "bool = False\n"
      "  (historic) print some output into file set by 'SetTestoutFile'";
    
    for (auto et : { ET_SEGM, ET_TRIG, ET_QUAD, ET_TET, ET_HEX, ET_PRISM, ET_PYRAMID })
      {
        string flagname = string(ElementTopology::GetElementName(et)) + "order";
        std::transform(flagname.begin(), flagname.end(), flagname.begin(), ::tolower);        
        docu.Arg(flagname);
      }
    return docu;
  }

  
  void FESpace :: SetNDof (size_t _ndof)
  {
    ndof = _ndof;
    while (ma->GetNLevels() > ndof_level.Size())
      ndof_level.Append (ndof);
    ndof_level.Last() = ndof;
  }

  void FESpace :: SetHarmonicProlongation (shared_ptr<BilinearForm> bfa)
  {
    prol = make_shared<HarmonicProlongation> (prol, // dynamic_pointer_cast<FESpace>(this->shared_from_this()),
                                              bfa);
    if (prol) prol->Update(*this);
  }
  
  void FESpace :: Update()
  {
    if (print)
      {
 	*testout << "Update FESpace, type = " << typeid(*this).name() << endl;
	*testout << "name = " << name << endl;
      }

    /*
    for (int i = 0; i < specialelements.Size(); i++)
      delete specialelements[i]; 
    specialelements.SetSize(0);
    */
    
    ma->UpdateBuffers();  // is free if netgen-mesh did not change
    int dim = ma->GetDimension();
    
    dirichlet_vertex.SetSize (ma->GetNV());
    dirichlet_edge.SetSize (ma->GetNEdges());
    if (dim == 3)
      dirichlet_face.SetSize (ma->GetNFaces());
    
    dirichlet_vertex = false;
    dirichlet_edge = false;
    dirichlet_face = false;

    /*
    {
      if (dirichlet_boundaries.Size())
        for (Ngs_Element ngel : ma->Elements(BND))  // .OmpSplit())
          if (dirichlet_boundaries[ngel.GetIndex()])
            {
              dirichlet_vertex[ngel.Vertices()] = true;
              if (dim >= 2)
                dirichlet_edge[ngel.Edges()] = true;
              if (dim == 3)
                dirichlet_face[ngel.Faces()[0]] = true;
            }
    }
    */
    for (auto vb : { BND, BBND, BBBND })
      {
        auto & dc = dirichlet_constraints[int(vb)];
        if (dc.Size())
          for (Ngs_Element ngel : ma->Elements(vb)) 
            if (dc[ngel.GetIndex()])
              {
                dirichlet_vertex[ngel.Vertices()] = true;
                // if (dim >= 2)
                if (dim-vb >= 1)
                  dirichlet_edge[ngel.Edges()] = true;
                // if (dim == 3)
                if (dim-vb >= 2)
                  dirichlet_face[ngel.Faces()[0]] = true;
              }
      }

    

    if (print)
      {
	(*testout) << "Dirichlet_vertex,1 = " << endl << dirichlet_vertex << endl;
	(*testout) << "Dirichlet_edge,1 = " << endl << dirichlet_edge << endl;
	(*testout) << "Dirichlet_face,1 = " << endl << dirichlet_face << endl;
      }


    ma->AllReduceNodalData (NT_VERTEX, dirichlet_vertex, NG_MPI_LOR);
    ma->AllReduceNodalData (NT_EDGE, dirichlet_edge, NG_MPI_LOR);
    ma->AllReduceNodalData (NT_FACE, dirichlet_face, NG_MPI_LOR);
    
    if (print)
      {
	(*testout) << "Dirichlet_vertex = " << endl << dirichlet_vertex << endl;
	(*testout) << "Dirichlet_edge = " << endl << dirichlet_edge << endl;
	(*testout) << "Dirichlet_face = " << endl << dirichlet_face << endl;
      }
  }


  void FESpace :: FinalizeUpdate()
  {
    static Timer timer ("FESpace::FinalizeUpdate");
    /*
    static Timer timer1 ("FESpace::FinalizeUpdate 1");
    static Timer timer2 ("FESpace::FinalizeUpdate 2");
    static Timer timer3 ("FESpace::FinalizeUpdate 3");
    static Timer tcol ("FESpace::FinalizeUpdate - coloring");
    static Timer tcolbits ("FESpace::FinalizeUpdate - bitarrays");
    static Timer tcolmutex ("FESpace::FinalizeUpdate - coloring, init mutex");
    */
    if (low_order_space) low_order_space -> FinalizeUpdate();

    RegionTimer reg (timer);
    // timer1.Start();
    dirichlet_dofs.SetSize (GetNDof());
    dirichlet_dofs.Clear();

    /*
    if (dirichlet_boundaries.Size())
      for (FESpace::Element el : Elements(BND))
        if (dirichlet_boundaries[el.GetIndex()])
          for (int d : el.GetDofs())
            if (IsRegularDof(d)) dirichlet_dofs.SetBit (d);
    */

    for (auto vb : { BND, BBND, BBBND })
      if (auto & dc = dirichlet_constraints[int(vb)]; dc.Size())
        for (FESpace::Element el : Elements(vb))
          if (dc[el.GetIndex()])
            for (int d : el.GetDofs())
              if (IsRegularDof(d)) dirichlet_dofs.SetBit (d);
    
    /*
    Array<DofId> dnums;
    for (auto i : Range(dirichlet_vertex))
      if (dirichlet_vertex[i])
	{
	  GetDofNrs (NodeId(NT_VERTEX,i), dnums);
	  for (DofId d : dnums)
	    if (d != -1) dirichlet_dofs.Set (d);
	}
    */
    // timer1.Stop();
    // timer2.Start();
    ParallelForRange
      (dirichlet_vertex.Size(),
       [&] (IntRange r)
       {
         Array<DofId> dnums;
         for (auto i : r)
           if (dirichlet_vertex[i])
             {
               GetDofNrs (NodeId(NT_VERTEX,i), dnums);
               for (DofId d : dnums)
                 if (IsRegularDof(d)) dirichlet_dofs.SetBitAtomic (d);
             }
       });
    // timer2.Stop();
    /*
    for (auto i : Range(dirichlet_edge))
      if (dirichlet_edge[i])
	{
	  GetDofNrs (NodeId(NT_EDGE,i), dnums);
	  for (DofId d : dnums)
	    if (d != -1) dirichlet_dofs.Set (d);
	}
    */
    ParallelForRange
      (dirichlet_edge.Size(),
       [&] (IntRange r)
       {
         Array<DofId> dnums;         
         for (auto i : r)
           if (dirichlet_edge[i])
             {
               GetDofNrs (NodeId(NT_EDGE,i), dnums);
               for (DofId d : dnums)
                 if (IsRegularDof(d)) dirichlet_dofs.SetBitAtomic (d);
             }
       });

    Array<DofId> dnums;             
    for (int i : Range(dirichlet_face))
      if (dirichlet_face[i])
	{
	  GetFaceDofNrs (i, dnums);
	  for (DofId d : dnums)
	    if (IsRegularDof(d)) dirichlet_dofs.SetBit (d);
	}
    
    // tcolbits.Start();
    free_dofs = make_shared<BitArray>(GetNDof());
    *free_dofs = dirichlet_dofs;
    free_dofs->Invert();
    
    for (auto i : Range(ctofdof))
      if (!(ctofdof[i] & VISIBLE_DOF)) //hidden or unused
	free_dofs->Clear(i);

    external_free_dofs = make_shared<BitArray>(GetNDof());
    *external_free_dofs = *free_dofs;
    for (auto i : Range(ctofdof))
      if (ctofdof[i] & CONDENSABLE_DOF)
	external_free_dofs->Clear(i);

    if (print)
      *testout << "freedofs = " << endl << *free_dofs << endl;
    // tcolbits.Stop();
    
    UpdateParallelDofs();

    if (print)
      *testout << "coloring ... " << flush;

    if (low_order_space)
      {
	for(auto vb : {VOL, BND, BBND, BBBND})
	  element_coloring[vb] = Table<int>(low_order_space->element_coloring[vb]);
      }
    else
      {
        // tcolmutex.Start();
      Array<MyMutex> locks(GetNDof());
      // tcolmutex.Stop();
      
      for (auto vb : { VOL, BND, BBND, BBBND })
      {
        /*
        tcol.Start();
        Array<int> col(ma->GetNE(vb));
        col = -1;

        int maxcolor = 0;
        
        int basecol = 0;
        Array<unsigned int> mask(GetNDof());

        size_t cnt = 0, found = 0;
        for (ElementId el : Elements(vb)) { cnt++; (void)el; } // no warning 

        do
          {
            mask = 0;

            Array<DofId> dofs;
            // for (auto el : Elements(vb))
            for (auto el : ma->Elements(vb))
              {
                if (!DefinedOn(el)) continue;
                if (col[el.Nr()] >= 0) continue;

                unsigned check = 0;
                GetDofNrs(el, dofs);
                for (auto d : dofs) // el.GetDofs())
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

                    if (HasAtomicDofs())
                      {
                        for (auto d : dofs) // el.GetDofs())
                          if (d != -1 && !IsAtomicDof(d)) mask[d] |= checkbit;
                      }
                    else
                      {
                        for (auto d : dofs) // el.GetDofs())
                          if (d != -1) mask[d] |= checkbit;
                      }
                  }
              }
            
            basecol += 8*sizeof(unsigned int); // 32;
          }
        while (found < cnt);

        tcol.Stop();

        Array<int> cntcol(maxcolor+1);
        cntcol = 0;

        for (ElementId el : Elements(vb))
          cntcol[col[el.Nr()]]++;
        
        Table<int> & coloring = element_coloring[vb];
        coloring = Table<int> (cntcol);

	cntcol = 0;
        for (ElementId el : Elements(vb))
          coloring[col[el.Nr()]][cntcol[col[el.Nr()]]++] = el.Nr();
        */



        // tcol.Start();
        Array<int> col(ma->GetNE(vb));
        col = -1;

        int maxcolor = 0;
        
        int basecol = 0;
        Array<unsigned int> mask(GetNDof());

        atomic<int> found(0);
        size_t cnt = 0;
        for (ElementId el : Elements(vb)) { cnt++; (void)el; } // no warning 

        while (found < cnt)
          {
            // mask = 0   | tasks;

            ParallelForRange
              (mask.Size(),
               [&] (IntRange myrange) { mask[myrange] = 0; });

            size_t ne = ma->GetNE(vb);

            ParallelForRange
              (ne, [&] (IntRange myrange)
               {
                 Array<DofId> dofs;
                 size_t myfound = 0;
                 
                 for (size_t nr : myrange)
                   {
                     ElementId el = { vb, nr };
                     if (!DefinedOn(el)) continue;
                     if (col[el.Nr()] >= 0) continue;
                     
                     unsigned check = 0;
                     GetDofNrs(el, dofs);
                     
                     if (HasAtomicDofs())
                       {
                         for (int i = dofs.Size()-1; i >= 0; i--)
                           if (!IsRegularDof(dofs[i]) || IsAtomicDof(dofs[i])) dofs.DeleteElement(i);
                       }
                     else
                       for (int i = dofs.Size()-1; i >= 0; i--)
                         if (!IsRegularDof(dofs[i])) dofs.DeleteElement(i);

                     QuickSort (dofs);   // sort to avoid dead-locks
                     // remove duplicates
                     size_t prev = 0;
                     for (size_t i = 1; i < dofs.Size(); i++)
                       if (dofs[i] != dofs[prev])
                         dofs[++prev] = dofs[i];
                     if (dofs.Size() > 0)                         
                       dofs.SetSize(prev+1);

                     /*
                     for (int i = 0; i+1 < dofs.Size(); )
                       if (dofs[i] == dofs[i+1])
                         {
                           for (int j = i+1; j+1 < dofs.Size(); j++)
                             dofs[j] = dofs[j+1];
                           dofs.SetSize(dofs.Size()-1);
                         }
                       else
                         i++;
                     */
                     
                     
                     for (auto d : dofs) 
                       locks[d].lock();
                     
                     for (auto d : dofs) 
                       check |= mask[d];
                     
                     if (check != UINT_MAX) // 0xFFFFFFFF)
                       {
                         myfound++;
                         unsigned checkbit = 1;
                         int color = basecol;
                         while (check & checkbit)
                           {
                             color++;
                             checkbit *= 2;
                           }
                         
                         col[el.Nr()] = color;
                         if (color > maxcolor) maxcolor = color;
                         
                         for (auto d : dofs) // el.GetDofs())
                           mask[d] |= checkbit;
                       }
                     
                     for (auto d : dofs) 
                       locks[d].unlock();
                   }
                 found += myfound;
               });
                 
            basecol += 8*sizeof(unsigned int); // 32;
          }

        // tcol.Stop();

        Array<int> cntcol(maxcolor+1);
        cntcol = 0;

        for (ElementId el : Elements(vb))
          cntcol[col[el.Nr()]]++;
        
        Table<int> & coloring = element_coloring[vb];
        coloring = Table<int> (cntcol);

	cntcol = 0;
        for (ElementId el : Elements(vb))
          coloring[col[el.Nr()]][cntcol[col[el.Nr()]]++] = el.Nr();
        
        if (print)
          *testout << "needed " << maxcolor+1 << " colors" 
                   << " for " << ((vb == VOL) ? "vol" : "bnd") << endl;
      }
      }
    // invalidate facet_coloring
    facet_coloring = Table<int>();
       
    level_updated = ma->GetNLevels();
    if (prol)
        prol->Update(*this);

    if (timing) Timing();
    updateSignal.Emit();
    // CheckCouplingTypes();
  }

  const Table<int> & FESpace :: FacetColoring() const
  {
    if (facet_coloring.Size()) return facet_coloring;

    size_t nf = ma->GetNFacets();
    Array<int> col(nf);
    col = -1;
    
    int maxcolor = 0;
    int basecol = 0;
    Array<unsigned int> mask(GetNDof());

    size_t cnt = nf, found = 0;
    Array<DofId> dofs, dofs1;
    Array<int> elnums, elnums_per;
    
    do
      {
        mask = 0;
        for (size_t f = 0; f < nf; f++)
          {
            if (col[f] >= 0) continue;
            
            ma->GetFacetElements(f,elnums);
            dofs.SetSize0();

            if (elnums.Size() == 1)
              {
                size_t f2 = ma->GetPeriodicFacet(f);
                // if (f2 < f) continue;
                if (f2 != f) // color both, left and right facet
                  {
                    ma->GetFacetElements (f2, elnums_per);
		    // if the facet is identified across subdomain
		    // boundary, we only have the surface element
		    // and not the other volume element!
		    // that case does not impact coloring
		    if (elnums_per.Size())
		      elnums.Append(elnums_per[0]);
                  }
              }
            for (auto el : elnums)
              {
                GetDofNrs(ElementId(VOL, el), dofs1);
                dofs += dofs1;
              }
            
            unsigned check = 0;
            for (auto d : dofs)
              if (IsRegularDof(d)) check |= mask[d];
            
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
                
                col[f] = color;
                if (color > maxcolor) maxcolor = color;
		
                for (auto d : dofs)
                  if (IsRegularDof(d)) mask[d] |= checkbit;
              }
          }
        
        basecol += 8*sizeof(unsigned int); // 32;
      }
    while (found < cnt);

    Array<int> cntcol(maxcolor+1);
    cntcol = 0;
    
    for (auto f : Range(nf))
      cntcol[col[f]]++;
    
    const_cast<Table<int>&> (facet_coloring) = Table<int> (cntcol);
    
    cntcol = 0;
    for (auto f : Range(nf))          
      facet_coloring[col[f]][cntcol[col[f]]++] = f;

    if (print)
      *testout << "needed " << maxcolor+1 << " colors for facet-coloring" << endl;

    return facet_coloring;
  }
  

  // FiniteElement & FESpace :: GetFE (ElementId ei, Allocator & alloc) const
  // {
  //   LocalHeap & lh = dynamic_cast<LocalHeap&> (alloc);
  //   switch(ei.VB())
  //     {
  //     case VOL:
  //       {
  //       FiniteElement * fe = NULL;
    
  //       if (DefinedOn (ei))
  //         {
  //           switch (ma->GetElType(ei))
  //             {
  //             case ET_TET: fe = tet; break;
  //             case ET_PYRAMID: fe = pyramid; break;
  //             case ET_PRISM: fe = prism; break;
  //             case ET_HEX: fe = hex; break;
  //             case ET_TRIG: fe = trig; break;
  //             case ET_QUAD: fe = quad; break;
  //             case ET_SEGM: fe = segm; break;
  //             case ET_POINT: fe = point; break;
  //             }
  //         }
  //       else
  //         {
  //           switch (ma->GetElType(ei))
  //             {
  //             case ET_TET: fe = dummy_tet; break;
  //             case ET_PYRAMID: fe = dummy_pyramid; break;
  //             case ET_PRISM: fe = dummy_prism; break;
  //             case ET_HEX: fe = dummy_hex; break;
  //             case ET_TRIG: fe = dummy_trig; break;
  //             case ET_QUAD: fe = dummy_quad; break;
  //             case ET_SEGM: fe = dummy_segm; break;
  //             case ET_POINT: fe = dummy_point; break;
  //             }
  //         }

  //       if (!fe)
  //         {
  //           /*
  //             Exception ex;
  //             ex << "FESpace" << GetClassName() << ", undefined eltype " 
  //             << ElementTopology::GetElementName(ma->GetElType(elnr))
  //             << ", order = " << ToString (order) << "\n";
  //             throw ex;
  //           */
  //           stringstream str;
  //           str << "FESpace" << GetClassName() << ", undefined eltype " 
  //               << ElementTopology::GetElementName(ma->GetElType(ei.Nr()))
  //               << ", order = " << ToString (order) << "\n";
  //           throw Exception (str.str());
  //         }
    
  //       return *fe;
  //       }
  //     case BND:
  //       switch (ma->GetElement(ei).GetType())
  //         {
  //         case ET_TRIG:  return *trig; 
  //         case ET_QUAD:  return *quad; 
  //         case ET_SEGM:  return *segm; 
  //         case ET_POINT: return *point; 
  //         case ET_TET: case ET_PYRAMID:
  //         case ET_PRISM: case ET_HEX: 
  //           ;
  //         }
  //       throw Exception ("GetFE BND: unknown type");
  //     case BBND:
  //       switch (ma->GetElement(ei).GetType())
  //         {
  //         case ET_SEGM: return *segm;
  //         case ET_POINT: return *point;
  //         default: ;
  //         }
  //       throw Exception("GetFE BBND: unknown type");
  //     default:
  //       __assume(false);
  //     }
  // }



  void FESpace :: SetOrder (ELEMENT_TYPE et, TORDER order)
  {
    if (order_policy == CONSTANT_ORDER || order_policy == OLDSTYLE_ORDER)
      order_policy = NODE_TYPE_ORDER;
    et_bonus_order[et] = order - this->order;
    
    string flagname = string(ElementTopology::GetElementName(et)) + "order";
    std::transform(flagname.begin(), flagname.end(), flagname.begin(), ::tolower);        
    flags.SetFlag(flagname, order);
    
    timestamp = 0;  // rerun first_update
  }


  
  void FESpace :: SetOrder (NodeId ni, int order)
  {
    throw Exception (string("FESpace::SetOrder not overloaded for space") + typeid(*this).name());
  }

  int FESpace :: GetOrder (NodeId ni) const
  {
    throw Exception (string("FESpace::GetOrder not overloaded for space") + typeid(*this).name());
  }
    


  /*
  const FiniteElement & FESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    return const_cast<FiniteElement&>(GetFE(ElementId(VOL,elnr),lh));
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

  /*
  void FESpace :: CheckCouplingTypeArray() const
  {
    if (ctofdof.Size() < GetNDof())
      throw Exception ("don't have a proper coupling type array, needed for static condensation");
  }
  */
  
  /*
  /// get coupling type of dof
  COUPLING_TYPE FESpace :: GetDofCouplingType (DofId dof) const 
  {
    if (ctofdof.Size()==0) //this is the standard case if the FESpace does not specify anything.
      return WIREBASKET_DOF;
    if (dof >= ctofdof.Size()) throw Exception("FESpace::GetDofCouplingType out of range. Did you set up ctofdof?");
    return ctofdof[dof];
  }
  */
  
  void FESpace :: SetDofCouplingType (DofId dof, COUPLING_TYPE ct) const
  {
    if (dof >= ctofdof.Size()) throw Exception("FESpace::SetDofCouplingType out of range");
    ctofdof[dof] = ct;
  }
  
  
  /// get coupling types of dofs
  void  FESpace :: GetDofCouplingTypes (int elnr, Array<COUPLING_TYPE> & ctypes) const
  { 
    ArrayMem<int,100> dnums;
    GetDofNrs(ElementId(VOL,elnr), dnums);
    ctypes.SetSize(dnums.Size());

    if (ctofdof.Size()==0)
      ctypes = INTERFACE_DOF;
    else
      {
        for (int i = 0; i < dnums.Size(); i++)
          if (IsRegularDof(dnums[i]))
            ctypes[i] = ctofdof[dnums[i]];
          else
            {
              if (dnums[i] == NO_DOF_NR)
                ctypes[i] = UNUSED_DOF;
              else
                ctypes[i] = HIDDEN_DOF;
            }
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
        for (auto d : dnums)
          if (IsRegularDof(d)) cnt[d]++;
      }
    for (int i : IntRange(0,ndof))
      {
        if (cnt[i] == 0 && ctofdof[i] != UNUSED_DOF)
          cout << "dof " << i << " not used, but coupling-type = " << ctofdof[i] << endl;
      }

    
    // cout << "check dofs" << endl;
    for(auto vb : {VOL,BND,BBND})
      for (ElementId id : ma->Elements(vb))
	{
	  GetDofNrs (id, dnums);
	  for (auto d : dnums)
	    if (IsRegularDof(d) && d >= ndof)
	      cout << "dof out of range: " << d << endl;
	}
  }
  /*
  void FESpace :: SetIrregularDofNrs (Array<DofId> & dnums) const
  {
    for (DofId & d : dnums)
      if (IsRegularDof(d))
        {
          auto ct = GetDofCouplingType(d);
          if (ct == UNUSED_DOF) d = NO_DOF_NR;
          if (ct == HIDDEN_DOF) d = NO_DOF_NR_CONDENSE;
        }
  }
  */


  FlatArray<VorB> FESpace :: GetDualShapeNodes (VorB vb) const
  {
    static VorB nodes[] = { VOL, BND, BBND, BBBND };
    return FlatArray<VorB> (ma->GetDimension()-int(vb)+1, &nodes[0]); 
  }

  
  void FESpace :: GetDofNrs (int elnr, Array<int> & dnums, COUPLING_TYPE ctype) const
  {
    GetDofNrs(ElementId(VOL,elnr),dnums,ctype);
  }

  void FESpace :: GetDofNrs (ElementId ei, Array<int> & dnums, COUPLING_TYPE ctype) const
  {
    ArrayMem<int,100> alldnums; 
    GetDofNrs(ei, alldnums);
    dnums.SetSize(0);

    if (ctofdof.Size() == 0)
      {
	if ( (INTERFACE_DOF & ctype) != 0)
          dnums = alldnums;
      }
    else
        for (auto d : alldnums)
	  if ( IsRegularDof(d) && ((ctofdof[d] & ctype) != 0) )
            dnums.Append(d);
      }

  void FESpace :: GetElementDofsOfType (ElementId ei, Array<DofId> & dnums, COUPLING_TYPE ctype) const
  {
    ArrayMem<int,100> alldnums; 
    GetDofNrs(ei, alldnums);
    dnums.SetSize(0);
    
    if (ctofdof.Size() == 0)
      {
        if ( (INTERFACE_DOF & ctype) != 0)
        {
          dnums = alldnums;
          for (int i : Range(alldnums.Size()))
            dnums[i] = i;
      }
      }
    else
      for (int i : Range(alldnums.Size()))
      {
        auto d = alldnums[i];
        if ( (d != -1) && ((ctofdof[d] & ctype) != 0) )
          dnums.Append(i);
      }
  }



  void FESpace :: GetNodeDofNrs (NODE_TYPE nt, int nr, Array<int> & dnums) const
  {
    GetDofNrs(NodeId(nt,nr),dnums);
  }

  void FESpace :: GetDofNrs (NodeId ni, Array<int> & dnums) const
  {
    switch (ni.GetType())
      {
      case NT_VERTEX: GetVertexDofNrs(ni.GetNr(), dnums); break;
      case NT_EDGE:   GetEdgeDofNrs(ni.GetNr(), dnums); break;
      case NT_FACE:   
        if (ma->GetDimension() == 3)
          GetFaceDofNrs(ni.GetNr(), dnums); 
        else
          {
            // have to convert from face# to element# on refined meshes
            auto surfel = ma->GetNode<2>(ni.GetNr()).surface_el;
            if (surfel >= 0)
              GetInnerDofNrs(surfel, dnums);
            else
              dnums.SetSize0();
          }
        break;
      case NT_CELL:   GetInnerDofNrs(ni.GetNr(), dnums); break;
      case NT_ELEMENT: case NT_FACET:
        GetDofNrs (NodeId(StdNodeType(ni.GetType(), ma->GetDimension()), ni.GetNr()), dnums);
        break;
      case NT_GLOBAL:
	GetGlobalDofNrs(ni.GetNr(), dnums);
	break;
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

  void FESpace :: GetGlobalDofNrs (int gnr, Array<int> & dnums) const
  {
    dnums.SetSize0 ();
  }

  shared_ptr<Prolongation> FESpace :: GetProlongation () const
  {
    return prol;
  }
  
  SymbolTable<shared_ptr<DifferentialOperator>> FESpace :: GetAdditionalEvaluators () const
  {
    return additional_evaluators;
  } 

  shared_ptr<BilinearFormIntegrator> FESpace :: GetIntegrator (VorB vb) const
  {
    if (integrator[vb])
      return integrator[vb];

    /*
        auto single_evaluator =  flux_evaluator;
        if (dynamic_pointer_cast<BlockDifferentialOperator>(single_evaluator))
          single_evaluator = dynamic_pointer_cast<BlockDifferentialOperator>(single_evaluator)->BaseDiffOp();
        
        auto trial = make_shared<ProxyFunction>(false, false, single_evaluator,
                                                nullptr, nullptr, nullptr, nullptr, nullptr);
        auto test  = make_shared<ProxyFunction>(true, false, single_evaluator,
                                                nullptr, nullptr, nullptr, nullptr, nullptr);
        fluxbli = make_shared<SymbolicBilinearFormIntegrator> (InnerProduct(trial,test), vb, false);
        single_fluxbli = fluxbli;
    */
    auto evaluator = GetEvaluator(vb);
    if (!evaluator) return nullptr;
    
    bool is_block = false;
    int block_dim;
    auto block_evaluator = dynamic_pointer_cast<BlockDifferentialOperator> (evaluator);
    if (block_evaluator)
      {
        is_block = true;
        block_dim = block_evaluator->BlockDim();
        evaluator = block_evaluator->BaseDiffOp();
      }
    auto trial = make_shared<ProxyFunction>(nullptr, // dynamic_pointer_cast<FESpace>(const_cast<FESpace*>(this)->shared_from_this()),
                                            false, false, evaluator,
                                            nullptr, nullptr, nullptr, nullptr, nullptr);
    auto test  = make_shared<ProxyFunction>(nullptr, // dynamic_pointer_cast<FESpace>(const_cast<FESpace*>(this)->shared_from_this()),
                                            true, false, evaluator,
                                            nullptr, nullptr, nullptr, nullptr, nullptr);
    shared_ptr<BilinearFormIntegrator> bli =
      make_shared<SymbolicBilinearFormIntegrator> (InnerProduct(trial,test), vb, VOL);
    
    if (is_block)
      bli = make_shared<BlockBilinearFormIntegrator> (bli, block_dim);
    const_cast<shared_ptr<BilinearFormIntegrator>&> (integrator[vb]) = bli;
    return bli;
  }


  
  void FESpace :: PrintReport (ostream & ost) const
  {
    ost << "type  = " << GetClassName() << endl
	<< "order = " << order << endl
	<< "dim   = " << dimension << endl
	<< "dgjmps= " << dgjumps << endl
    << "autoupdate= " << autoupdate << endl
	<< "complex = " << iscomplex << endl;
    ost << "definedon = " << definedon[VOL] << endl;
    ost << "definedon boundary = " << definedon[BND] << endl;
    ost << "definedon codim 2 = " << definedon[BBND] << endl;
    if (!free_dofs) return;

    ost << "ndof = " << GetNDof() << endl;
    int ntype[16] = { 0 };

    for (auto ct : ctofdof) ntype[ct]++;
    if (ntype[UNUSED_DOF]) ost << "unused = " << ntype[UNUSED_DOF] << endl;
    if (ntype[HIDDEN_DOF]) ost << "hidden = " << ntype[HIDDEN_DOF] << endl;
    if (ntype[LOCAL_DOF])  ost << "local  = " << ntype[LOCAL_DOF] << endl;

    /*
    int nfree = 0;
    for (int i = 0; i < free_dofs->Size(); i++)
      if ((*free_dofs)[i])
	nfree++;
    */
  }
  
  void FESpace :: DoArchive (Archive & archive)
  {
    archive & order & dimension & iscomplex & dgjumps & autoupdate & print & level_updated;
    archive & definedon[VOL] & definedon[BND] & definedon[BBND];
    archive & dirichlet_boundaries & dirichlet_dofs & free_dofs & external_free_dofs;
    archive & dirichlet_vertex & dirichlet_edge & dirichlet_face;
  }

  std::tuple<Shallow<shared_ptr<MeshAccess>>, Flags> FESpace :: GetCArgs ()
  {
    return std::make_tuple(Shallow(GetMeshAccess()), GetFlags());
  }

  Array<MemoryUsage> FESpace :: GetMemoryUsage () const
  {
    Array<MemoryUsage> mu;
    mu += { "coupling types", ctofdof.Size()*sizeof(COUPLING_TYPE), 1 };
    return mu;
  }

  
  std::list<std::tuple<std::string,double>> FESpace :: Timing () const
  {
    double time;
    std::list<std::tuple<std::string,double>> results;
    LocalHeap lh (100000, "FESpace - Timing");

    time = RunTiming([&]() {
        ParallelForRange( IntRange(ma->GetNE()), [&] ( IntRange r )
        {
	  LocalHeap &clh = lh, lh = clh.Split();
          Array<int> dnums;
	  for (int i : r)
            GetDofNrs (ElementId(VOL,i), dnums);
	});
    });
    results.push_back(std::make_tuple<std::string,double>("GetDofNrs",1e9*time / (ma->GetNE())));

    time = RunTiming([&]() {
        ParallelForRange( IntRange(ma->GetNE()), [&] ( IntRange r )
        {
	  LocalHeap &clh = lh, lh = clh.Split();

	  for (int i : r)
	    {
	      HeapReset hr(lh);
	      GetFE (ElementId(VOL,i), lh);
	    }
	});
    });
    results.push_back(std::make_tuple<std::string,double>("GetFE",1e9 * time / (ma->GetNE()))); 

    time = RunTiming([&]() {
        ParallelFor( IntRange(ma->GetNE()), [&] (size_t i)
          {
	    /* Ng_Element ngel = */ ma->GetElement(ElementId(VOL,i));
          });
    });
    results.push_back(std::make_tuple<std::string,double>("Get Ng_Element", 1e9 * time / (ma->GetNE())));


    time = RunTiming([&]() {
        ParallelForRange( IntRange(ma->GetNE()), [&] ( IntRange r )
        {
	  LocalHeap &clh = lh, lh = clh.Split();
          for (int i : r)
            {
              HeapReset hr(lh);
              /* ElementTransformation & trafo = */ ma->GetTrafo(ElementId(VOL, i), lh);
            }
        });
    });
    results.push_back(std::make_tuple<std::string,double>("GetTrafo", 1e9 * time / (ma->GetNE())));


    Array<int> global(GetNDof());
    global = 0;
    time = RunTiming([&]() {
        ParallelForRange( IntRange(ma->GetNE()), [&] ( IntRange r )
        {
          Array<DofId> dnums;
          for (size_t i : r)
            {
              GetDofNrs( { VOL, i }, dnums);
              for (auto d : dnums)
                {
                // global[d]++;
                  AsAtomic (global[d])++;
                  AsAtomic (global[d])++;
                  // int oldval = 0;
                  // AsAtomic (global[d]).compare_exchange_strong (oldval, 1);
                  // AsAtomic (global[d]).compare_exchange_strong (oldval, 0);
                }
                                                                
            }
        });
    });
    results.push_back(std::make_tuple<std::string,double>("Count els of dof", 1e9 * time / (ma->GetNE())));


    
#ifdef TIMINGSEQUENTIAL

    double starttime = WallTime();
    int steps = 0;
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

    return results;

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
	  output.SetBit(i);
    if (freedofsonly && free_dofs->Size()) {
      output.And(*free_dofs);
    }
  }



  shared_ptr<Table<int>> FESpace :: CreateSmoothingBlocks (const Flags & flags) const
  {
    bool eliminate_internal =
      flags.GetDefineFlag("eliminate_internal") ||
      flags.GetDefineFlag("condense");
      
    auto freedofs = GetFreeDofs(eliminate_internal);

    if (flags.GetDefineFlag("subassembled"))
      {
        freedofs = make_shared<BitArray>(*freedofs);
        for (auto i : Range(GetNDof()))
          if (GetDofCouplingType(i) != WIREBASKET_DOF)
            freedofs->Clear(i);
      }

    FilteredTableCreator creator(freedofs.get());

    /*
    for ( ; !creator.Done(); creator++)
      {
	for (size_t i = 0; i < nd; i++)
	  if (freedofs->Test(i))
	    creator.Add (i, i);
      }
    */


    // new version for ngs24:

    Array<string> blocktypes { flags.GetStringListFlag("blocktype") };

    
    if (flags.GetStringFlag("blocktype")=="vertex" ||
        flags.GetStringFlag("blocktype")=="edge" ||
        flags.GetStringFlag("blocktype")=="face" ||
        flags.GetStringFlag("blocktype")=="facet" ||
        flags.GetStringFlag("blocktype")=="element" ||        
        flags.GetStringFlag("blocktype")=="vertexpatch" ||        
        flags.GetStringFlag("blocktype")=="vertexedge")
      blocktypes += flags.GetStringFlag("blocktype");

    
    if (blocktypes.Size())
      {
        Array<DofId> dofs;
        
        for ( ; !creator.Done(); creator++)
          {
            size_t base = 0;

            if (blocktypes.Contains("vertex"))
              {
                for (size_t i : Range(ma->GetNV()))        
                  {
                    GetDofNrs (NodeId(NT_VERTEX, i), dofs);
                    for (auto d : dofs)
                      if (IsRegularDof(d))
                        creator.Add (base+i, d);
                  }
                base += ma->GetNV();
              }

            if (blocktypes.Contains("edge") ||
                (ma->GetDimension()==2 && blocktypes.Contains("facet")) )
              {
                for (size_t i : Range(ma->GetNEdges()))        
                  {
                    GetDofNrs (NodeId(NT_EDGE, i), dofs);
                    for (auto d : dofs)
                      if (IsRegularDof(d))
                        creator.Add (base+i, d);
                  }
                base += ma->GetNEdges();
              }

            if (blocktypes.Contains("face") ||
                (ma->GetDimension()==3 && blocktypes.Contains("facet")) )              
              {
                for (size_t i : Range(ma->GetNFaces()))        
                  {
                    GetDofNrs (NodeId(NT_FACE, i), dofs);
                    for (auto d : dofs)
                      if (IsRegularDof(d))
                        creator.Add (base+i, d);
                  }
                base += ma->GetNFaces();
              }

            if (blocktypes.Contains("element"))
              {
                for (int i : Range(ma->GetNE(VOL)))
                  {
                    GetDofNrs ( ElementId(VOL,i), dofs);
                    for (auto d : dofs)
                      if (IsRegularDof(d))
                        creator.Add (base+i, d);
                  }
                base += ma->GetNE(VOL);
              }
            
            if (blocktypes.Contains("vertexedge"))
              {
                for (size_t i : Range(ma->GetNV()))
                  {
                    GetDofNrs (NodeId(NT_VERTEX, i), dofs);
                    for (auto d : dofs)
                      if (IsRegularDof(d))              
                        creator.Add (base+i, d);
                  }
                for (size_t i : Range(ma->GetNEdges()))        
                  {
                    Ng_Node<1> edge = ma->GetNode<1> (i);
                    GetDofNrs (NodeId(NT_EDGE, i), dofs);
                    for (auto d : dofs)
                      if (IsRegularDof(d))
                        for (int k = 0; k < 2; k++)
                          creator.Add (base+edge.vertices[k], d);
                  }
                base += ma->GetNV();                
              }

            if (blocktypes.Contains("vertexpatch"))
              {
                for (size_t i : Range(ma->GetNV()))
                  {
                    GetDofNrs (NodeId(NT_VERTEX, i), dofs);
                    for (auto d : dofs)
                      if (IsRegularDof(d))              
                        creator.Add (base+i, d);
                  }
                for (size_t i : Range(ma->GetNEdges()))        
                  {
                    Ng_Node<1> edge = ma->GetNode<1> (i);
                    GetDofNrs (NodeId(NT_EDGE, i), dofs);
                    for (auto d : dofs)
                      if (IsRegularDof(d))
                        for (int k = 0; k < 2; k++)
                          creator.Add (base+edge.vertices[k], d);
                  }
                for (size_t i : Range(ma->GetNFaces()))        
                  {
                    auto vnums = ma->GetFacePNums(i);
                    GetDofNrs (NodeId(NT_FACE, i), dofs);
                    for (auto d : dofs)
                      if (IsRegularDof(d))
                        for (auto v : vnums)
                          creator.Add (base+v, d);
                  }
                // 3D only
                for (size_t i : Range(ma->GetNElements(3)))
                  {
                    auto vnums = ma->GetElement({VOL,i}).Vertices();
                    GetDofNrs (NodeId(NT_CELL, i), dofs);
                    for (auto d : dofs)
                      if (IsRegularDof(d))
                        for (auto v : vnums)
                          creator.Add (base+v, d);
                  }
                base += ma->GetNV();                
              }

            if (creator.GetMode()==1)
              creator.SetSize(base);
          }

        Table<int> table = creator.MoveTable();
        return make_shared<Table<int>> (std::move(table));
      }
    


    
    
    if (flags.GetStringFlag("blocktype") == "edgepatch")
      {
        
        Array<DofId> dofs;
        for ( ; !creator.Done(); creator++)
          {
            // EFI

            /*
            for (size_t i : Range(ma->GetNV()))
              {
                GetDofNrs (NodeId(NT_VERTEX, i), dofs);
                for (auto d : dofs)
                  if (IsRegularDof(d))              
                    creator.Add (i, d);
              }
            */
            
            for (size_t i : Range(ma->GetNEdges()))        
              {
                // Ng_Node<1> edge = ma->GetNode<1> (i);
                
                GetDofNrs (NodeId(NT_EDGE, i), dofs);
                for (auto d : dofs)
                  if (IsRegularDof(d))
                    creator.Add (i, d);
              }
            
            for (size_t i : Range(ma->GetNFaces()))        
              {
                // Ng_Node<2> face = ma->GetNode<2> (i);
                
                GetDofNrs (NodeId(NT_FACE, i), dofs);
                for (auto d : dofs)
                  if (IsRegularDof(d))
                    for (auto e : ma->GetFaceEdges(i))
                      creator.Add (e, d);
              }
            
            if(ma->GetDimension() == 3)
              for(size_t i : Range(ma->GetNE()))
                {
                  GetDofNrs(NodeId(NT_CELL, i), dofs);
                  auto eledges = ma->GetElEdges(ElementId(VOL, i));
                  for(auto d : dofs)
                    if(IsRegularDof(d))
                      for(auto v : eledges)
                        creator.Add(v, d);
                }
        
          }
      }

    if (flags.GetStringFlag("blocktype") == "facepatch")
      {
        
        Array<DofId> dofs;
        for ( ; !creator.Done(); creator++)
          {
            // EFI
            
            for (size_t i : Range(ma->GetNFaces()))        
              {
                GetDofNrs (NodeId(NT_FACE, i), dofs);
                for (auto d : dofs)
                  if (IsRegularDof(d))
                    creator.Add (i, d);
              }
            
            if(ma->GetDimension() == 3)
              for(size_t i : Range(ma->GetNE()))
                {
                  GetDofNrs(NodeId(NT_CELL, i), dofs);
                  auto elfaces = ma->GetElFaces(ElementId(VOL, i));
                  for(auto d : dofs)
                    if(IsRegularDof(d))
                      for(auto f : elfaces)
                        creator.Add(f, d);
                }
        
          }
      }

    
    else
      
      { // default is vertexpatch
        Array<DofId> dofs;
        for ( ; !creator.Done(); creator++)
          {
            // VEFI

            for (size_t i : Range(ma->GetNV()))
              {
                GetDofNrs (NodeId(NT_VERTEX, i), dofs);
                for (auto d : dofs)
                  if (IsRegularDof(d))              
                    creator.Add (i, d);
              }
            for (size_t i : Range(ma->GetNEdges()))        
              {
                Ng_Node<1> edge = ma->GetNode<1> (i);
                
                GetDofNrs (NodeId(NT_EDGE, i), dofs);
                for (auto d : dofs)
                  if (IsRegularDof(d))
                    for (int k = 0; k < 2; k++)
                      creator.Add (edge.vertices[k], d);
              }
            
            for (size_t i : Range(ma->GetNFaces()))        
              {
                Ng_Node<2> face = ma->GetNode<2> (i);
                
                GetDofNrs (NodeId(NT_FACE, i), dofs);
                for (auto d : dofs)
                  if (IsRegularDof(d))
                    for (int k = 0; k < face.vertices.Size(); k++)
                      creator.Add (face.vertices[k], d);
              }
            
            if(ma->GetDimension() == 3)
              for(size_t i : Range(ma->GetNE()))
                {
                  GetDofNrs(NodeId(NT_CELL, i), dofs);
                  auto elverts = ma->GetElVertices(ElementId(VOL, i));
                  for(auto d : dofs)
                    if(IsRegularDof(d))
                      for(auto v : elverts)
                        creator.Add(v, d);
                }
          }
        /*
          
          for (int i = 0; i < nfa; i++)
          {
          Ng_Node<2> face = ma->GetNode<2> (i);
          for (int k = 0; k < face.vertices.Size(); k++)
          creator.Add (face.vertices[k], GetFaceDofs(i));
          }
          
          for (int i = 0; i < ni; i++)
          for (auto v : ma->GetElement(ElementId(VOL,i)).Vertices())
          creator.Add (v, GetElementDofs(i));
          
          break; 
        */
      }

    
        
    // return make_shared<Table<int>> (creator.MoveTable());
    Table<int> table = creator.MoveTable();
    if (print)
      *testout << "smoothing blocks = " << endl << table << endl;
    return make_shared<Table<int>> (std::move(table));
  }


  void FESpace :: ConnectAutoUpdate() 
  {
    /*
    if (this->weak_from_this().expired())
      throw Exception("Given pointer is not managed by a shared ptr.");
    if (this->DoesAutoUpdate())
      this->GetMeshAccess()->updateSignal.Connect(this, [this]()
      {
        this->Update();
        this->FinalizeUpdate();
      });
    */
    if (this->DoesAutoUpdate())
      this->GetMeshAccess()->updateSignal.Connect(this, [this]()
      {
        this->Update();
        this->FinalizeUpdate();
      });
  }
  
  void FESpace :: SetDefinedOn (VorB vb, const BitArray & defon)
  {
    definedon[vb].SetSize(defon.Size());
    
    for (int i = 0; i < defon.Size(); i++)
      definedon[vb][i] = defon.Test(i);

    if (low_order_space)
      low_order_space -> SetDefinedOn (vb, defon);

    timestamp = 0;
  }

  void FESpace :: SetDirichletBoundaries (const BitArray & dirbnds)
  {
    dirichlet_boundaries = dirbnds;
    if (low_order_space)
      low_order_space -> SetDirichletBoundaries (dirbnds);
  }

  shared_ptr<BitArray> FESpace :: GetFreeDofs (bool external) const
  {
    if (external)
      return external_free_dofs;
    else
      return free_dofs;
  }

  void IterateElements (const FESpace & fes, 
			VorB vb, 
			LocalHeap & clh, 
			const function<void(FESpace::Element,LocalHeap&)> & func)
  {
    static mutex copyex_mutex;
    const Table<int> & element_coloring = fes.ElementColoring(vb);
    
    if (task_manager)
      {
        for (FlatArray<int> els_of_col : element_coloring)
          {
            SharedLoop2 sl(els_of_col.Range());

            task_manager -> CreateJob
              ( [&] (const TaskInfo & ti) 
                {
                  LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
                  ArrayMem<int,100> temp_dnums;

                  for (int mynr : sl)
                    {
                      HeapReset hr(lh);
                      FESpace::Element el(fes, 
                                          ElementId (vb, els_of_col[mynr]), 
                                          temp_dnums, lh);
                      
                      func (std::move(el), lh);
                    }

                  ProgressOutput::SumUpLocal();
                } );
          }
        return;
      }


    Exception * ex = nullptr;

    for (FlatArray<int> els_of_col : element_coloring)
      ParallelForRange( IntRange(els_of_col.Size()), [&] ( IntRange r )
      {
        LocalHeap lh = clh.Split();
        Array<int> temp_dnums;

        // lh.ClearValues();
        
        for (int i : r)
          {
            try
              {
                HeapReset hr(lh);
                FESpace::Element el(fes, ElementId (vb, els_of_col[i]), temp_dnums, lh);
                
                func (std::move(el), lh);
              }
	    
            catch (const Exception & e)
              {
                {
                  lock_guard<mutex> guard(copyex_mutex);
                  if (!ex)
		    ex = new Exception (e);
                }
              }
	    catch (...)
	      { ; }
          }
      // cout << "lh, used size = " << lh.UsedSize() << endl;
    });
    
    if (ex)
      {
        throw Exception (*ex);
      }
  }
  

  ostream & operator<< (ostream & ost, COUPLING_TYPE ct)
  {
    switch (ct)
      {
      case UNUSED_DOF: ost << "unused"; break;
      case HIDDEN_DOF:  ost << "hidden"; break;
      case LOCAL_DOF:  ost << "local"; break;
      case CONDENSABLE_DOF:  ost << "condensable"; break;
      case INTERFACE_DOF: ost << "interface"; break;
      case NONWIREBASKET_DOF: ost << "non-wirebasket"; break;
      case WIREBASKET_DOF: ost << "wirebasket"; break;
      case EXTERNAL_DOF: ost << "external"; break;
      case VISIBLE_DOF: ost << "visible"; break;
      case ANY_DOF: ost << "any"; break;
      };
    return ost;
  }


  shared_ptr<BaseMatrix> FESpace ::
  GetMassOperator (shared_ptr<CoefficientFunction> rho,
                   shared_ptr<Region> defon,
                   LocalHeap & lh) const
  {
    return make_shared<ApplyMass> (dynamic_pointer_cast<FESpace>(const_cast<FESpace*>(this)->shared_from_this()),
                                   rho, false, defon, lh);    
  }

  shared_ptr<BaseMatrix> FESpace :: CreateMassOperator (shared_ptr<CoefficientFunction> rho,
                                                        shared_ptr<Region> defon,
                                                        bool inverse,
                                                        LocalHeap & lh) const
  {
    throw Exception(string("CreateMassOperator not overloaded, space = ") + typeid(*this).name());
  }
  
  
  void FESpace :: SolveM (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                          LocalHeap & lh) const
  {
    cout << "SolveM is only available for L2-space, not for " << typeid(*this).name() << endl;
  }

  void FESpace :: ApplyM (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                          LocalHeap & lh) const
  {
    cout << "ApplyM is only available for L2-space, not for " << typeid(*this).name() << endl;
  }


  shared_ptr<BaseMatrix> FESpace :: GetTraceOperator (shared_ptr<FESpace> tracespace, bool avg) const
  {
    throw Exception("GetTraceOperator not overloaded");
  }

  void FESpace :: GetTrace (const FESpace & tracespace, const BaseVector & in, BaseVector & out, bool avg,
                            LocalHeap & lh) const
  {
    throw Exception("GetTrace not overloaded");
  }

  void FESpace :: GetTraceTrans (const FESpace & tracespace, const BaseVector & in, BaseVector & out, bool avg,
                                 LocalHeap & lh) const
  {
    throw Exception("GetTrace not overloaded");
  }

  

  shared_ptr<BaseMatrix> FESpace :: ConvertL2Operator (shared_ptr<FESpace> l2space) const
  {
    LocalHeap lh(100000000);
    
    Array<int> classnr(ma->GetNE());
    
    FlatArray<int> vertex_map;
    if (const PeriodicFESpace * periodic = dynamic_cast<const PeriodicFESpace*> (this))
      vertex_map.Assign (periodic->GetVertexMap());
    
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         classnr[el.Nr()] = 
           SwitchET<ET_TRIG,ET_QUAD,ET_TET,ET_HEX>
           (el.GetType(),
            [el, vertex_map] (auto et) {
             if (vertex_map.Size())
               return ET_trait<et.ElementType()>::GetClassNr(vertex_map[el.Vertices()]);
             else
               return ET_trait<et.ElementType()>::GetClassNr(el.Vertices());
            });
       });

    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);
    Table<size_t> table = creator.MoveTable();

    shared_ptr<BaseMatrix> sum;
  
    // size_t ne = ma->GetNE();

    // scaling for multiple dofs (if actually not an L2 space)
    Array<int> cnt(l2space->GetNDof());
    cnt = 0;
    
    for (auto elclass_inds : table)
      {
        HeapReset hr(lh);
        if (elclass_inds.Size() == 0) continue;
        
        ElementId ei(VOL,elclass_inds[0]);
        auto & fel = GetFE (ei, lh);
        auto & fel_l2 = l2space->GetFE (ei, lh);

        auto & trafo = GetFEElementTransformation(fel.ElementType());
        MixedFiniteElement fel_mixed(fel, fel_l2);
        auto evaluator = GetEvaluator(VOL);
        auto l2evaluator = l2space->GetEvaluator(VOL);

        auto trial = make_shared<ProxyFunction>(dynamic_pointer_cast<FESpace>(const_cast<FESpace*>(this)->shared_from_this()),
                                                false, false, evaluator,
                                                nullptr, nullptr, nullptr, nullptr, nullptr);
        auto trial_l2 = make_shared<ProxyFunction>(l2space, 
                                                   false, false, l2evaluator,
                                                   nullptr, nullptr, nullptr, nullptr, nullptr);
        auto test_l2  = make_shared<ProxyFunction>(l2space,
                                                   true, false, l2evaluator,
                                                   nullptr, nullptr, nullptr, nullptr, nullptr);
        shared_ptr<BilinearFormIntegrator> bfi_mass_mixed =
          make_shared<SymbolicBilinearFormIntegrator> (InnerProduct(trial,test_l2), VOL, VOL);
        shared_ptr<BilinearFormIntegrator> bfi_mass_l2 =
          make_shared<SymbolicBilinearFormIntegrator> (InnerProduct(trial_l2,test_l2), VOL, VOL);

        Matrix<> mass_l2(fel_l2.GetNDof(), fel_l2.GetNDof());
        Matrix<> mass_mixed(fel_l2.GetNDof(), fel.GetNDof());

        bfi_mass_l2->CalcElementMatrix (fel_l2, trafo, mass_l2, lh);
        bfi_mass_mixed->CalcElementMatrix (fel_mixed, trafo, mass_mixed, lh);

        CalcInverse (mass_l2);
        Matrix<> trans = mass_l2 * mass_mixed;

        Table<DofId> xdofs(elclass_inds.Size(), fel.GetNDof()),
          ydofs(elclass_inds.Size(), fel_l2.GetNDof());
        
        Array<DofId> dnumsx, dnumsy;
        for (auto i : Range(elclass_inds))
          {
            ElementId ei(VOL, elclass_inds[i]);
            GetDofNrs(ei, dnumsx);
            l2space->GetDofNrs(ei, dnumsy);
            for (auto d : dnumsy)
              cnt[d]++;
            xdofs[i] = dnumsx;
            ydofs[i] = dnumsy;
          }

        auto mat = make_shared<ConstantElementByElementMatrix<>>
          (l2space->GetNDof(), this->GetNDof(),
           trans, std::move(ydofs), std::move(xdofs));

        if (sum)
          sum = make_shared<SumMatrix>(sum, mat);
        else
          sum = mat;
      }

    if (sum == nullptr) {
      // return dummy op for empty FESpace
      Table<DofId> xdofs(0, 0), ydofs(0, 0);
      Matrix<> mat(0,0);
      sum = make_shared<ConstantElementByElementMatrix<>>
	(l2space->GetNDof(), this->GetNDof(),
	 mat, std::move(ydofs), std::move(xdofs));
    }

    bool multiple = false;
    for (auto c : cnt)
      if (c > 1) multiple = true;

    if (multiple)
      {
        VVector<> scaling(cnt.Size());
        for (auto i : Range(cnt))
          scaling(i) = cnt[i] ? 1.0/cnt[i] : 0;
        auto diagmat = make_shared<DiagonalMatrix<>> (scaling);
        sum = make_shared<ProductMatrix> (diagmat, sum);
      }

    if (IsParallel())
      {
        if (multiple)
          throw Exception("parallel convertL2 operator only available for non-coupling spaces");
        sum = make_shared<ParallelMatrix> (sum, GetParallelDofs(), l2space->GetParallelDofs(), C2C);
      }
    
    return sum;    
  }

  
  
  void FESpace :: UpdateParallelDofs ( )
  {
    if (!ma->GetCommunicator().ValidCommunicator()) return;
    // if (ma->GetCommunicator().Size() == 1) return;

    static Timer timer ("FESpace::UpdateParallelDofs"); RegionTimer reg(timer);

    size_t ndof = GetNDof();
    Array<NodeId> dofnodes (ndof);
    dofnodes = NodeId (NT_VERTEX, -1);

    Array<int> dnums;
    for (NODE_TYPE nt : { NT_VERTEX, NT_EDGE, NT_FACE, NT_CELL, NT_GLOBAL })
      for (NodeId ni : ma->Nodes(nt))
	{
	  GetDofNrs (ni, dnums);
	  for (auto d : dnums)
	    if (IsRegularDof(d)) dofnodes[d] = ni;
	} 

    TableCreator<int> creator(ndof);
    for ( ; !creator.Done(); creator++)
      for (size_t i = 0; i < ndof; i++)
        if (dofnodes[i].GetNr() != -1) 
          creator.Add (i, ma->GetDistantProcs (dofnodes[i]));
    
    paralleldofs = make_shared<ParallelDofs>
      (ma->GetCommunicator(), creator.MoveTable(), dimension, iscomplex);

    if (ma->GetCommunicator().AllReduce (ctofdof.Size(), NG_MPI_SUM))
      paralleldofs -> AllReduceDofData (ctofdof, NG_MPI_MAX);
  }


  bool FESpace :: IsParallel() const
  {
    return ma->GetCommunicator().Size() > 1;
  }

  size_t FESpace :: GetNDofGlobal() const 
  { 
    return paralleldofs ?
      paralleldofs -> GetNDofGlobal() : GetNDof(); 
  }

  BitArray FESpace :: GetDofs (const Region & reg) const
  {
    BitArray ba(GetNDof());
    ba.Clear();
    for (auto el : Elements(reg.VB()))
      if (reg.Mask().Test(el.GetIndex()))
        for (auto d : el.GetDofs())
          if (IsRegularDof(d))
            ba.SetBit(d);
    return ba;
  }



  ProxyNode FESpace ::
  MakeProxyFunction (bool testfunction,
                     const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const
  {
    shared_ptr<FESpace> fes = dynamic_pointer_cast<FESpace> (const_cast<FESpace*>(this)->shared_from_this());

    auto proxy = make_shared<ProxyFunction>  (fes, testfunction, fes->IsComplex(),
                                              fes->GetEvaluator(),
                                              fes->GetFluxEvaluator(),
                                              fes->GetEvaluator(BND),
                                              fes->GetFluxEvaluator(BND),
                                              fes->GetEvaluator(BBND),
                                              fes->GetFluxEvaluator(BBND));
    
    auto add_diffops = fes->GetAdditionalEvaluators();
    for (int i = 0; i < add_diffops.Size(); i++)
      proxy->SetAdditionalEvaluator (add_diffops.GetName(i), add_diffops[i]);
    
    proxy = addblock(proxy);
    return ProxyNode (proxy);
  }
  

  ProxyNode FESpace :: GetProxyFunction(bool testfunction) const
  {
    return MakeProxyFunction
      ( testfunction, [&] (shared_ptr<ProxyFunction> proxy) { return proxy; });
  }
  

  

  



  NodalFESpace :: NodalFESpace (shared_ptr<MeshAccess> ama,
				const Flags & flags,
                                bool parseflags)
    : FESpace (ama, flags)
  {
    name="NodalFESpace";
    type = "nodal";
    
    prol = make_shared<LinearProlongation> (GetMeshAccess());

    if (order >= 2)
      {
	Flags loflags;
	loflags.SetFlag ("order", 1);
	loflags.SetFlag ("dim", dimension);
	if (dgjumps) loflags.SetFlag ("dgjumps");
	if (iscomplex) loflags.SetFlag ("complex");
	low_order_space = make_shared<NodalFESpace> (ma, loflags);
      }
    hb_defined = flags.GetDefineFlag("hb");

    // SetDummyFE<ScalarDummyFE> ();

    auto one = make_shared<ConstantCoefficientFunction> (1);
    if (ma->GetDimension() == 2)
      {
	integrator[VOL] = make_shared<MassIntegrator<2>> (one);
        integrator[BND] = make_shared<RobinIntegrator<2>> (one);
      }
    else
      {
	integrator[VOL] = make_shared<MassIntegrator<3>> (one);
	integrator[BND] = make_shared<RobinIntegrator<3>> (one);
      }
    
    if (dimension > 1)
      {
	integrator[VOL] = make_shared<BlockBilinearFormIntegrator> (integrator[VOL], dimension);
	integrator[BND] = make_shared<BlockBilinearFormIntegrator> (integrator[BND], dimension);  
      }


    switch (ma->GetDimension())
      {
      case 1:
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<1>>>();
	  evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<1>>>();
	  flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<1>>>();
          break;
        }
      case 2:
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
	  evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
	  flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpGradientBoundary<2>>>();
          break;
        }
      case 3:
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
	  evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
	  flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpGradientBoundary<3>>>();
          break;
        }
      }
    if (dimension > 1)
      {
	evaluator[VOL] = make_shared<BlockDifferentialOperator> (evaluator[VOL], dimension);
	evaluator[BND] = make_shared<BlockDifferentialOperator> (evaluator[BND], dimension);
	flux_evaluator[VOL] = make_shared<BlockDifferentialOperator> (flux_evaluator[VOL], dimension);
	flux_evaluator[BND] = make_shared<BlockDifferentialOperator> (flux_evaluator[BND], dimension);
      }

    
  }

  NodalFESpace :: ~NodalFESpace ()
  {
    ;
  }

  FiniteElement & NodalFESpace :: GetFE(ElementId ei, Allocator & lh) const
  {
    ELEMENT_TYPE et = ma->GetElType(ei);
    if(order == 1)
      {
        return SwitchET(et, [&lh] (auto type) -> FiniteElement&
                        { return * new (lh) ScalarFE<type.ElementType(),1>(); });

      }
    else
      {
        if(hb_defined)
          {
            switch (ma->GetElType(ei))
              {
              case ET_TET:     return *(new (lh) FE_Tet2HB);
              case ET_PRISM:   return *(new (lh) FE_Prism1);
              case ET_PYRAMID: return *(new (lh) FE_Pyramid1);
              case ET_TRIG:    return *(new (lh) FE_Trig2HB);
              case ET_QUAD:    return *(new (lh) ScalarFE<ET_QUAD,1>);
              case ET_SEGM:    return *(new (lh) FE_Segm2);
              case ET_POINT:   return *(new (lh) FE_Point);
              default:
                throw Exception ("Inconsistent element type in NodalFESpace::GetFE, hb defined");
              }
          }
        else
          {
            switch (ma->GetElType(ei))
              {
                // case ET_TET:     return * new (lh) FE_Tet2;
              case ET_TET:     return * new (lh) ScalarFE<ET_TET,2>();
              case ET_PRISM:   return * new (lh) FE_Prism1;
              case ET_PYRAMID: return * new (lh) FE_Pyramid1;
              case ET_HEX:     return * new (lh) FE_Hex20;                
              case ET_TRIG:    return * new (lh) FE_Trig2;
                // case ET_QUAD:    return *(new (lh) ScalarFE<ET_QUAD,1>;
              case ET_QUAD:    return * new (lh) FE_Quad2Serendipity;
              case ET_SEGM:    return * new (lh) FE_Segm2;
              case ET_POINT:   return * new (lh) FE_Point;
              default:
                throw Exception ("Inconsistent element type in NodalFESpace::GetFE, no hb defined");
              }
          }
      }
  }

  
  void NodalFESpace :: Update()
  {
    FESpace::Update();
    if (low_order_space) low_order_space -> Update();

    bool first_update = GetTimeStamp() < ma->GetTimeStamp();
    if (first_update) timestamp = NGS_Object::GetNextTimeStamp();

    if (first_update)
      {
	used_vertex.SetSize(ma->GetNV());
	used_vertex = false;
        used_edge.SetSize(ma->GetNEdges());
        used_edge = false;

	for (auto vb : { VOL, BND })
	  ParallelFor
	    (ma->GetNE(vb), [&] (size_t nr)
	     {
	       ElementId ei(vb, nr);
	       Ngs_Element el = (*ma)[ei];

	       if (!DefinedOn (el)) return;

	       used_vertex[el.Vertices()] = true;
               if (order > 1)
                 used_edge[el.Edges()] = true;
	     });
      }

    // if (ma->GetNLevels() > ndlevel.Size())
      {
	size_t ndof = ma->GetNV();
        if (order > 1)
          {
            for (auto el : Elements(VOL))
              for (DofId d : el.GetDofs()) ndof = max2(ndof, size_t(d+1));              
            
            for (auto el : Elements(BND))
              for (DofId d : el.GetDofs()) ndof = max2(ndof, size_t(d+1));           
          }
	// ndlevel.Append (ndof);
        SetNDof(ndof);

        for (int i = 0; i < ma->GetNLevels()-1; i++)
          ndof_level[i] = ma->GetNVLevel(i);
      }

      prol->Update(*this);
      
    if (dirichlet_boundaries.Size())
      {
	dirichlet_dofs.SetSize (GetNDof());
	dirichlet_dofs.Clear();
	for (auto el : Elements(BND))
          if (dirichlet_boundaries[el.GetIndex()])
            for (DofId d : el.GetDofs())
              if (IsRegularDof(d)) dirichlet_dofs.SetBit (d);
      }

    UpdateCouplingDofArray();
  }

  void NodalFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(GetNDof());
    ParallelFor
      (ma->GetNV(), [&] (size_t i)
       {
         ctofdof[i] = used_vertex[i] ? WIREBASKET_DOF : UNUSED_DOF;
       });

    if(order > 1 && GetNDof() > ma->GetNV())
      ParallelFor
        (ma->GetNEdges(), [&] (size_t i)
         {
           ctofdof[ma->GetNV() + i] = used_edge[i] ? INTERFACE_DOF : UNUSED_DOF;
         });
  }

  void NodalFESpace :: DoArchive (Archive & archive)
  {
    ;
    /*
    if (archive.Input())
      {
	ndlevel.SetSize(1);
	ndlevel[0] = ma->GetNV();
      }
    */
  }

 
  void NodalFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    if (order == 1)
      dnums = ma->GetElement(ei).Vertices();
    else
      dnums = ma->GetElement(ei).Points();
      
    if (!DefinedOn (ei)) dnums = -1;
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

  
  shared_ptr<Array<int>>
  NodalFESpace :: CreateDirectSolverClusters (const Flags & flags) const
  {
    auto spclusters = make_shared<Array<int>> (GetNDof());
    Array<int> & clusters = *spclusters;
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
	return nullptr;
      }

    return spclusters;
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
    

    // trig = new FE_NcTrig1;
    // segm = new FE_Segm0;
      
    auto one = make_shared<ConstantCoefficientFunction> (1);
    if (ma->GetDimension() == 2)
      {
	integrator[VOL] = make_shared<MassIntegrator<2>> (one);
        integrator[BND] = make_shared<RobinIntegrator<2>> (one);
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
      }
    else
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
        
	integrator[VOL].reset (new MassIntegrator<3> (new ConstantCoefficientFunction(1)));
	integrator[BND].reset (new RobinIntegrator<3> (new ConstantCoefficientFunction(1)));
      }

    if (dimension > 1)
      {
	integrator[VOL] = make_shared<BlockBilinearFormIntegrator> (integrator[VOL], dimension);
	integrator[BND] = make_shared<BlockBilinearFormIntegrator> (integrator[BND], dimension);
      }
  }

  NonconformingFESpace :: ~NonconformingFESpace ()
  {
    ;
  }

  FiniteElement & NonconformingFESpace :: GetFE (ElementId ei, Allocator & lh) const
  {
    if (ei.IsVolume())
      {
        switch (ma->GetElType(ei))
          {
          case ET_TRIG: return *(new (lh) FE_NcTrig1);
          case ET_TET: return *(new (lh) FE_NcTet1);
          default: throw Exception ("Element type not available in NonconformingFESpace::GetFE, vol");
          }
      }
        
    if (ei.IsBoundary())
      {
        switch (ma->GetElType(ei))
          {
          case ET_SEGM: return *(new (lh) FE_Segm0);
          case ET_TRIG: return *(new (lh) FE_Trig0);
          default: throw Exception ("Element type not available in NonconformingFESpace::GetFE, bnd");
          }
      }
    throw Exception ("NonconormingFE: only vol or bnd");
    /*
    switch (ma->GetElType(ei))
      {
      case ET_TRIG: return *(new (lh) FE_NcTrig1);
      case ET_SEGM: return *(new (lh) FE_Segm0);
      default: throw Exception ("Element type not available in NonconformingFESpace::GetFE");
      }
    */
  }

  /*
  size_t NonconformingFESpace :: GetNDof () const throw()
  {
    // return ma->GetNEdges();
    return ma->GetNFacets();
  }
  */

  void NonconformingFESpace :: Update()
  {
    ctofdof.SetSize (ma->GetNFacets());
    ctofdof = UNUSED_DOF;
    for (auto el : ma->Elements(VOL))
      ctofdof[el.Facets()] = WIREBASKET_DOF;

    SetNDof (ma->GetNFacets());
    
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

  void NonconformingFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    // dnums = ma->GetElEdges(ei);
    dnums = ma->GetElFacets(ei);
    if (!DefinedOn(ei))
      dnums = -1;
  }


  NonconformingSurfaceFESpace :: 
  NonconformingSurfaceFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name="NonconformingSurfaceFESpace(nonconforming)";
    // defined flags
    DefineDefineFlag("nonconforming");
    if (parseflags) CheckFlags(flags);
    
    // prol = new LinearProlongation(*this);
    

    // trig = new FE_NcTrig1;
    // segm = new FE_Segm0;
      
    auto one = make_shared<ConstantCoefficientFunction> (1);
    if (ma->GetDimension() < 3)
      {
        throw Exception("NonconformingSurfaceFESpace only for dim=3");
      }
    else
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
        flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpGradientBoundary<3>>>();
        
	integrator[VOL].reset (new MassIntegrator<3> (new ConstantCoefficientFunction(1)));
	integrator[BND].reset (new RobinIntegrator<3> (new ConstantCoefficientFunction(1)));
      }

    if (dimension > 1)
      {
	integrator[VOL] = make_shared<BlockBilinearFormIntegrator> (integrator[VOL], dimension);
	integrator[BND] = make_shared<BlockBilinearFormIntegrator> (integrator[BND], dimension);
      }
  }

  NonconformingSurfaceFESpace :: ~NonconformingSurfaceFESpace ()
  {
    ;
  }

  FiniteElement & NonconformingSurfaceFESpace :: GetFE (ElementId ei, Allocator & lh) const
  {
    if (ei.IsBoundary())
      {
        switch (ma->GetElType(ei))
          {
          case ET_TRIG: return *(new (lh) FE_NcTrig1);
          default: throw Exception ("Element type not available in NonconformingSurfaceFESpace::GetFE, bnd");
          }
      }
    throw Exception ("NonconormingFE: only bnd");
  }

  /*
  size_t NonconformingSurfaceFESpace :: GetNDof () const throw()
  {
    return ma->GetNEdges();
  }
  */

  void NonconformingSurfaceFESpace :: Update()
  {
    ctofdof.SetSize (ma->GetNEdges());
    ctofdof = UNUSED_DOF;
    for (auto el : ma->Elements(BND))
      ctofdof[el.Edges()] = WIREBASKET_DOF;
    SetNDof (ma->GetNEdges()); 
  }

  void NonconformingSurfaceFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    dnums = ma->GetElEdges(ei);
    if (!DefinedOn(ei))
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
      n_el_dofs = 1;
    }
    else
    {
      if (ma->GetDimension() == 2)
        n_el_dofs = 4;
      else
        n_el_dofs = 6;
    }

    // SetDummyFE<ScalarDummyFE> ();
    static ConstantCoefficientFunction one(1);

    if (ma->GetDimension() == 2)
      {
        integrator[VOL].reset (new MassIntegrator<2> (&one));
        integrator[BND] = 0;
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
      }
    else
      {
        integrator[VOL].reset (new MassIntegrator<3> (&one));
	integrator[BND] = 0;
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
      }
    
    if (dimension > 1)
      integrator[VOL] = make_shared<BlockBilinearFormIntegrator> (integrator[VOL], dimension);
  }
  
  ElementFESpace :: ~ElementFESpace ()
  {
    ;
  }

  void ElementFESpace :: Update()
  {
    /*
    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (n_el_dofs * ma->GetNE());
    */
    SetNDof (n_el_dofs * ma->GetNE());
  }

  void ElementFESpace :: DoArchive (Archive & archive)
  {
    FESpace :: DoArchive (archive);
    archive /* & ndlevel */ & n_el_dofs;
  }

  FiniteElement & ElementFESpace :: GetFE (ElementId ei, Allocator & lh) const
  {
    ELEMENT_TYPE et = ma->GetElType(ei);
    if (ei.VB() != VOL)
      return SwitchET(et, [&lh] (auto type) -> FiniteElement&
                      { return * new (lh) DummyFE<type.ElementType()>(); });
    
    if (order == 0)
      {
        return SwitchET(et, [&lh] (auto type) -> FiniteElement&
                        { return * new (lh) ScalarFE<type.ElementType(),0>(); });
      }
    else
      {
        return SwitchET(et, [&lh] (auto type) -> FiniteElement&
                        { return * new (lh) ScalarFE<type.ElementType(),1>(); });
      }
    throw Exception ("Illegal element type in ElementFESpace::GetFE");
  }
  
  void ElementFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    if(ei.VB()!=VOL) { dnums.SetSize(0); return; }
    if (order == 0)
      {
	dnums.SetSize(1);
	dnums[0] = ei.Nr();
      }
    else if (order == 1)
      {
	switch (ma->GetElType(ei))
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
	  dnums[i] = n_el_dofs*ei.Nr()+i;
      }
  }


 
  SurfaceElementFESpace :: 
      SurfaceElementFESpace (shared_ptr<MeshAccess> ama, const Flags& flags, bool parseflags)
  : FESpace (ama, flags)
  {
    name="SurfaceElementFESpace(surfl2)";
    if(parseflags) CheckFlags(flags);
    
    // prol = new SurfaceElementProlongation (GetMeshAccess(), *this);

    if (order == 0)
    {
      /*
      trig    = new ScalarFE<ET_TRIG,0>;
      quad    = new ScalarFE<ET_QUAD,0>;
      segm    = new FE_Segm0;
      */
      n_el_dofs = 1;
    }

    else if (order == 1)
    {
      /*
      trig    = new ScalarFE<ET_TRIG,1>;
      quad    = new ScalarFE<ET_QUAD,1>;
      segm    = new FE_Segm1;
      */
      if (ma->GetDimension() == 2)
        n_el_dofs = 2;
      else
        n_el_dofs = 4;
    }

    else if (order == 2)
    {
      /*
      trig    = new FE_Trig2HB;
      quad    = new ScalarFE<ET_QUAD,1>;
      segm    = new FE_Segm2;
      */
      if (ma->GetDimension() == 2)
        n_el_dofs = 3;
      else
        n_el_dofs = 9;
    }

    integrator[BND].reset (new RobinIntegrator<3> (new ConstantCoefficientFunction(1)));

    if (dimension > 1)
      integrator[BND] = make_shared<BlockBilinearFormIntegrator> (integrator[BND], dimension);
  }

  
  SurfaceElementFESpace :: ~SurfaceElementFESpace ()
  {
    ;
  }

  void SurfaceElementFESpace :: Update()
  {
    // const MeshAccess & ma = GetMeshAccess();
    /*
      while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (n_el_dofs * ma->GetNSE());
    */
    SetNDof (n_el_dofs*ma->GetNE(BND));
  }
  
  FiniteElement & SurfaceElementFESpace :: GetFE (ElementId ei, Allocator & lh) const
  {
    auto et = ma->GetElement(ei).GetType();
    if (ei.IsBoundary())
      {
        if (order == 0)
          switch (et)
            {
            case ET_SEGM: return * new (lh) ScalarFE<ET_SEGM,0>;
            case ET_TRIG: return * new (lh) ScalarFE<ET_TRIG,0>;
            case ET_QUAD: return * new (lh) ScalarFE<ET_QUAD,0>;
            default:
              ;
            }
        else if (order == 1)
          switch (et)
            {
            case ET_SEGM: return * new (lh) FE_Segm1;
            case ET_TRIG: return * new (lh) ScalarFE<ET_TRIG,1>;
            case ET_QUAD: return * new (lh) ScalarFE<ET_QUAD,1>;
            default:
              ;
            }
        else if (order == 2)
          switch (et)
            {
            case ET_SEGM: return * new (lh) FE_Segm2;
            case ET_TRIG: return * new (lh) FE_Trig2HB;
            case ET_QUAD: throw Exception("SurfaceFESpace, second order quad not here");
              // return * new (lh) ScalarFE<ET_QUAD,1>;
            default:
              ;
            }
      }
    return SwitchET(et, [&lh] (auto type) -> FiniteElement&
                    { return * new (lh) DummyFE<type.ElementType()>(); });
  }
  
  
  void SurfaceElementFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    if(ei.VB()!=BND) {dnums.SetSize0 (); return; }
    if (order == 0)
      {
	dnums.SetSize(1);
	dnums[0] = ei.Nr();
      }
    else if (order == 1)
      {
	switch (ma->GetElType(ei))
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
	  dnums[i] = n_el_dofs*ei.Nr()+i;
      }
    else if (order == 2)
      {
	switch (ma->GetElType(ei))
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
	  dnums[i] = n_el_dofs*ei.Nr()+i;
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

    if (flags.GetDefineFlag("low_order_space"))
      {
        Flags loflags = flags;
        loflags.SetFlag("low_order_space", false);
        low_order_space = make_shared<CompoundFESpace> (ama, loflags, parseflags);
      }
    prol = make_shared<CompoundProlongation> (this);

    needs_transform_vec = false;
    all_the_same = true;
  }

  
  CompoundFESpace :: CompoundFESpace (shared_ptr<MeshAccess> ama,
				      const Array<shared_ptr<FESpace>> & aspaces,
				      const Flags & flags, bool parseflags)
    : CompoundFESpace (ama, flags, parseflags)
  {
    for (auto space : aspaces)
      AddSpace (space);
  }

  
  void CompoundFESpace :: AddSpace (shared_ptr<FESpace> fes)
  {
    if (fes->IsComplex() != this->IsComplex())
      throw Exception("Product space of spaces with complex and real spaces is not allowed");
    if (fes->GetDimension() != this->GetDimension())
      throw Exception("Product space of spaces with different dimensions is not allowed");
    
    
    spaces.Append (fes);
    dynamic_pointer_cast<CompoundProlongation> (prol) -> AddProlongation (fes->GetProlongation());
    if (fes->NeedsTransformVec())      
      needs_transform_vec = true;
    if (fes != spaces[0])
      all_the_same = false;
    
    if (low_order_space)
      {
        if (fes->LowOrderFESpacePtr())
          dynamic_pointer_cast<CompoundFESpace>(low_order_space)
            -> AddSpace (fes->LowOrderFESpacePtr());
        else
          low_order_space.reset();
      }
  }

  CompoundFESpace :: ~CompoundFESpace ()
  {
    ;
  }
  
  void CompoundFESpace :: SetDefinedOn (VorB vb, const BitArray & defon)
  {
    for (auto & space : spaces)
      space -> SetDefinedOn (vb, defon);
    FESpace::SetDefinedOn (vb, defon);
  }
    
  void CompoundFESpace :: Update()
  {
    FESpace::Update();
    if (low_order_space)
      low_order_space->Update();
    
    cummulative_nd.SetSize (spaces.Size()+1);
    cummulative_nd[0] = 0;
    for (int i = 0; i < spaces.Size(); i++)
      {
        if(do_subspace_update)
          spaces[i] -> Update();
	cummulative_nd[i+1] = cummulative_nd[i] + spaces[i]->GetNDof();
      }

    SetNDof (cummulative_nd.Last());

    /*
    free_dofs = make_shared<BitArray> (GetNDof());
    free_dofs->Clear();
    for (int i = 0; i < spaces.Size(); i++)
      {
	shared_ptr<BitArray> freedofsi = spaces[i]->GetFreeDofs(false);
	for (int j = 0; j < freedofsi->Size();j++)
	  if (freedofsi->Test(j)) 
	    free_dofs->Set(cummulative_nd[i]+j);
      }
    external_free_dofs = make_shared<BitArray> (GetNDof());
    external_free_dofs->Clear();
    for (int i = 0; i < spaces.Size(); i++)
      {
	shared_ptr<BitArray> freedofsi = spaces[i]->GetFreeDofs(true);
	for (int j = 0; j < freedofsi->Size();j++)
	  if (freedofsi->Test(j)) 
	    external_free_dofs->Set(cummulative_nd[i]+j);
      }
    */
    
    
    bool has_atomic = false;
    for (auto & space : spaces)
      if (space->HasAtomicDofs())
        has_atomic = true;
    if (has_atomic)
      {
        is_atomic_dof = BitArray(GetNDof());
        is_atomic_dof = false;
        for (int i = 0; i < spaces.Size(); i++)
          {
            FESpace & spacei = *spaces[i];
            IntRange r(cummulative_nd[i], cummulative_nd[i+1]);
            if (spacei.HasAtomicDofs())
              {
                for (size_t j = 0; j < r.Size(); j++)
                  if (spacei.IsAtomicDof(j))
                    is_atomic_dof.SetBit(r.begin()+j);
              }
          }       
      }
    // cout << "AtomicDofs = " << endl << is_atomic_dof << endl;

    // prol -> Update(*this);  // do we need that ? 

    UpdateCouplingDofArray();


    if (low_order_space)
      {
        shared_ptr<BaseMatrix> sum_emb;
        for (size_t i = 0; i < spaces.Size(); i++)
          {
            auto emb_i = spaces[i]->LowOrderEmbedding();
            auto hi_range = GetRange(i);
            auto lo_range = dynamic_pointer_cast<CompoundFESpace>(low_order_space)->GetRange(i);
            emb_i = make_shared<EmbeddedMatrix> (GetNDof(), hi_range, emb_i);
            emb_i = make_shared<EmbeddedTransposeMatrix> (low_order_space->GetNDof(), lo_range, emb_i);
            if (sum_emb)
              sum_emb = make_shared<SumMatrix> (sum_emb, emb_i);
            else
              sum_emb = emb_i;
          }
        low_order_embedding = sum_emb;
        // cout << "embedding = " << *low_order_embedding << endl;
      }

    
    if (print)
      {
	(*testout) << "Update compound fespace" << endl;
	(*testout) << "cumulative dofs start at " << cummulative_nd << endl;
      }
  }

  void CompoundFESpace :: FinalizeUpdate()
  {
    for (int i = 0; i < spaces.Size(); i++)
      if(do_subspace_update)
        spaces[i] -> FinalizeUpdate();

    FESpace::FinalizeUpdate();


    // dirichlet-dofs from sub-spaces
    // ist das umsonst ? (JS)
    // haben jetzt ja immer dirichlet-dofs
    bool has_dirichlet_dofs = false;
    for (int i = 0; i < spaces.Size(); i++)
      if (spaces[i]->GetFreeDofs()) 
	has_dirichlet_dofs = true;

    auto comm = ma->GetCommunicator();
    has_dirichlet_dofs = comm.AllReduce (has_dirichlet_dofs, NG_MPI_LOR);

    if (has_dirichlet_dofs)
      {
	free_dofs = make_shared<BitArray> (GetNDof());
	free_dofs->Set();

	for (int i = 0; i < spaces.Size(); i++)
	  {
	    shared_ptr<BitArray> free_dofs_sub = spaces[i]->GetFreeDofs();
	    if (free_dofs_sub)
	      {
		int base = cummulative_nd[i];
		int nd = cummulative_nd[i+1] - base;
		for (int i = 0; i < nd; i++)
		  if (!free_dofs_sub->Test(i))
		    free_dofs->Clear (base+i);
	      }
	  }

        for (int i = 0; i < ctofdof.Size(); i++)
          if (ctofdof[i] == UNUSED_DOF)
            free_dofs->Clear(i);

	dirichlet_dofs = *free_dofs;
	dirichlet_dofs.Invert();

        external_free_dofs = make_shared<BitArray>(GetNDof());
        *external_free_dofs = *free_dofs;
        for (int i = 0; i < ctofdof.Size(); i++)
          if (ctofdof[i] & CONDENSABLE_DOF)
            external_free_dofs->Clear(i);


        if (print)
          (*testout) << "compound fespace freedofs:" << endl
                     << free_dofs << endl;

      }
  }



  void CompoundFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(this->GetNDof());

    for (int i = 0; i < spaces.Size(); i++)
      {
        if (spaces[i] -> CouplingTypeArrayAvailable())
          for (int j=0; j< spaces[i]->GetNDof();j++)
            ctofdof[cummulative_nd[i]+j] = spaces[i]->GetDofCouplingType(j);
        else
          ctofdof.Range(cummulative_nd[i], cummulative_nd[i+1]) = WIREBASKET_DOF;
      }
    // *testout << "CompoundFESpace :: UpdateCouplingDofArray() presents \n ctofdof = \n" << ctofdof << endl;
  }

  ProxyNode CompoundFESpace ::
  MakeProxyFunction (bool testfunction,
                     const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const
  {
    if (this->GetEvaluator())
      return FESpace::MakeProxyFunction (testfunction, addblock);

    auto fes = dynamic_pointer_cast<CompoundFESpace> (const_cast<CompoundFESpace*>(this)->shared_from_this());    
    
    std::vector<ProxyNode> l;
    int nspace = fes->GetNSpaces();
    for (int i = 0; i < nspace; i++)
      {
        l.push_back ((*fes)[i]->MakeProxyFunction
                     (testfunction,
                      [&] (shared_ptr<ProxyFunction> proxy)
                      {
                        // auto addcomp = [] (shared_ptr<DifferentialOperator> diffop, int comp) -> shared_ptr<CompoundDifferentialOperator>
                        // { return diffop ? make_shared<CompoundDifferentialOperator> (diffop, comp) : nullptr; };
                        
                        auto block_eval = make_shared<CompoundDifferentialOperator> (proxy->Evaluator(), i);

                        // auto block_deriv_eval = addcomp (proxy->DerivEvaluator(), i);
                        shared_ptr <CompoundDifferentialOperator> block_deriv_eval = nullptr;
                        if (proxy->DerivEvaluator() != nullptr)
                          block_deriv_eval = make_shared<CompoundDifferentialOperator> (proxy->DerivEvaluator(), i);
                        shared_ptr <CompoundDifferentialOperator> block_trace_eval = nullptr;
                        if (proxy->TraceEvaluator() != nullptr)
                          block_trace_eval = make_shared<CompoundDifferentialOperator> (proxy->TraceEvaluator(), i);
                        shared_ptr <CompoundDifferentialOperator> block_ttrace_eval = nullptr;
                        if (proxy->TTraceEvaluator() != nullptr)
                          block_ttrace_eval = make_shared<CompoundDifferentialOperator> (proxy->TTraceEvaluator(),i);
                        shared_ptr <CompoundDifferentialOperator> block_trace_deriv_eval = nullptr;
                        if (proxy->TraceDerivEvaluator() != nullptr)
                          block_trace_deriv_eval = make_shared<CompoundDifferentialOperator> (proxy->TraceDerivEvaluator(), i);
                        shared_ptr <CompoundDifferentialOperator> block_ttrace_deriv_eval = nullptr;
                        if (proxy->TTraceDerivEvaluator() != nullptr)
                          block_ttrace_deriv_eval = make_shared<CompoundDifferentialOperator> (proxy->TTraceDerivEvaluator(),i);
                        
                        auto block_proxy = make_shared<ProxyFunction> (fes, testfunction, fes->IsComplex(),                                                                                                                  block_eval, block_deriv_eval, block_trace_eval, block_trace_deriv_eval,
                                                                       block_ttrace_eval, block_ttrace_deriv_eval);
                        
                        SymbolTable<shared_ptr<DifferentialOperator>> add = proxy->GetAdditionalEvaluators();
                        for (int j = 0; j < add.Size(); j++)
                          block_proxy->SetAdditionalEvaluator(add.GetName(j),
                                                              make_shared<CompoundDifferentialOperator> (add[j], i));
                        
                        block_proxy = addblock(block_proxy);
                        return block_proxy;
                      }));
      }
    return ProxyNode (std::move(l));
  }
  
  shared_ptr<BaseMatrix> CompoundFESpace :: EmbeddingOperator (int spacenr) const
  {
    shared_ptr<BaseMatrix> emb = make_shared<Embedding> (GetNDof(), GetRange(spacenr), IsComplex());
    if (IsParallel())
      emb = make_shared<ParallelMatrix> (emb, spaces[spacenr] -> GetParallelDofs(),
                                         GetParallelDofs(), C2C);
    return emb;
  }

  
  shared_ptr<BaseMatrix> CompoundFESpace :: RestrictionOperator (int spacenr) const
  {
    shared_ptr<BaseMatrix> emb = make_shared<EmbeddingTranspose> (GetNDof(), GetRange(spacenr), IsComplex());
    if (IsParallel())
      emb = make_shared<ParallelMatrix> (emb,
                                         GetParallelDofs(),
                                         spaces[spacenr] -> GetParallelDofs(),
                                         C2C);
    return emb;
  }
  
  

  FiniteElement & CompoundFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    FlatArray<FiniteElement*> fea(spaces.Size(), alloc);
    if (!all_the_same)
      {
        for (int i = 0; i < fea.Size(); i++)
          fea[i] = &spaces[i]->GetFE(ei, alloc);
      }
    else
      {
        if (fea.Size() > 0)
          {
            fea[0] = &spaces[0]->GetFE(ei, alloc);
            for (int i = 1; i < fea.Size(); i++)
              fea[i] = fea[0];
          }
      }
    return *new (alloc) CompoundFiniteElement (fea);
  }


  void CompoundFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    if (spaces.Size() == 0)
      {
        dnums.SetSize0();
        return;
      }
    
    spaces[0]->GetDofNrs(ei, dnums);
      
    ArrayMem<DofId,500> hdnums;
    for (int i = 1; i < spaces.Size(); i++)
      {
	spaces[i]->GetDofNrs (ei, hdnums);
        int base = dnums.Size();
        int base_cum = cummulative_nd[i];
        dnums.SetSize(base+hdnums.Size());

        for (auto j : Range(hdnums))
          {
            DofId val = hdnums[j];
            if (IsRegularDof(val))
              val += base_cum;
            dnums[base+j] = val;
          }
      }
  }

  void CompoundFESpace :: GetDofNrs (NodeId ni, Array<DofId> & dnums) const
  {
    ArrayMem<DofId,500> hdnums;
    dnums.SetSize0();
    for (int i = 0; i < spaces.Size(); i++)
      {
	spaces[i]->GetDofNrs (ni, hdnums);
        int base = dnums.Size();
        int base_cum = cummulative_nd[i];
        dnums.SetSize(base+hdnums.Size());

        for (auto j : Range(hdnums))
          {
            DofId val = hdnums[j];
            if (IsRegularDof(val))
              val += base_cum;
            dnums[base+j] = val;
          }
      }
  }
  

  void CompoundFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    ArrayMem<DofId,500> hdnums;
    dnums.SetSize(0);
    for (int i = 0; i < spaces.Size(); i++)
      {
	spaces[i]->GetVertexDofNrs (vnr, hdnums);
	for (int j = 0; j < hdnums.Size(); j++)
	  if (IsRegularDof(hdnums[j]))
	    dnums.Append (hdnums[j]+cummulative_nd[i]);
	  else
	    dnums.Append (hdnums[j]);
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
	  if (IsRegularDof(hdnums[j]))
	    dnums.Append (hdnums[j]+cummulative_nd[i]);
	  else
	    dnums.Append (hdnums[j]);
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
	  if (IsRegularDof(hdnums[j]))
	    dnums.Append (hdnums[j]+cummulative_nd[i]);
	  else
	    dnums.Append (hdnums[j]);
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
	  if (IsRegularDof(hdnums[j]))
	    dnums.Append (hdnums[j]+cummulative_nd[i]);
	  else
	    dnums.Append (hdnums[j]);
      }
  }
  
  void CompoundFESpace :: GetElementDofsOfType (ElementId ei, Array<DofId> & dnums, 
                                                COUPLING_TYPE ctype) const
  {
    ArrayMem<int,500> hdnums;
    dnums.SetSize(0);
    int offset = 0;
    for (int i = 0; i < spaces.Size(); i++)
      {
        spaces[i]->GetElementDofsOfType (ei, hdnums, ctype);
        for (int j = 0; j < hdnums.Size(); j++)
          dnums.Append (hdnums[j]+offset);

        if (i+1 < spaces.Size())
        {
          spaces[i]->GetDofNrs (ei, hdnums);
          offset += hdnums.Size();
        }
      }      
  }

  

  void CompoundFESpace :: SolveM(CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                                 LocalHeap & lh) const
  {
    for (size_t i = 0; i < spaces.Size(); i++)
      {
        auto veci = vec.Range (GetRange(i));
        spaces[i] -> SolveM (rho, veci, definedon, lh);
      }
  }
    
  void CompoundFESpace :: ApplyM(CoefficientFunction * rho, BaseVector & vec,
                                 Region * definedon,
                                 LocalHeap & lh) const
  {
    for (size_t i = 0; i < spaces.Size(); i++)
      {
        auto veci = vec.Range (GetRange(i));
        spaces[i] -> ApplyM (rho, veci, definedon, lh);
      }
  }
    




  template <class T>
  void CompoundFESpace::T_TransformMat (ElementId ei, 
                                        SliceMatrix<T> mat, TRANSFORM_TYPE tt) const
  {
    size_t base = 0;
    LocalHeapMem<100005> lh("CompoundFESpace - transformmat");
    for (int i = 0; i < spaces.Size(); i++)
      {
        HeapReset hr(lh);
        size_t nd = spaces[i]->GetFE(ei, lh).GetNDof();

        if (tt & TRANSFORM_MAT_LEFT)
          spaces[i]->TransformMat (ei, mat.Rows(base, base+nd), TRANSFORM_MAT_LEFT);
        if (tt & TRANSFORM_MAT_RIGHT)
          spaces[i]->TransformMat (ei, mat.Cols(base, base+nd), TRANSFORM_MAT_RIGHT);

	base += nd;
      }
  }

  template <class T>
  void CompoundFESpace::T_TransformVec (ElementId ei, 
                                        SliceVector<T> vec, TRANSFORM_TYPE tt) const
  {
    if (!needs_transform_vec) return;
    
    LocalHeapMem<100006> lh("CompoundFESpace - transformvec");
    for (int i = 0, base = 0; i < spaces.Size(); i++)
      {
        HeapReset hr(lh);
        int nd = spaces[i]->GetFE(ei, lh).GetNDof();

	spaces[i]->TransformVec (ei, vec.Range(base, base+nd), tt);
	base += nd;
      }
  }



  template NGS_DLL_HEADER
  void CompoundFESpace::T_TransformVec<double> 
  (ElementId ei, SliceVector<double> vec, TRANSFORM_TYPE tt) const;
  template NGS_DLL_HEADER
  void CompoundFESpace::T_TransformVec<Complex>
  (ElementId ei, SliceVector<Complex> vec, TRANSFORM_TYPE tt) const;
  
  template
  void CompoundFESpace::T_TransformMat<double> 
  (ElementId ei, SliceMatrix<double> mat, TRANSFORM_TYPE tt) const;
  template
  void CompoundFESpace::T_TransformMat<Complex> 
  (ElementId ei, SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const;
  
  
  void CompoundFESpace::VTransformMR (ElementId ei, 
				      SliceMatrix<double> mat, TRANSFORM_TYPE tt) const 
  {
    T_TransformMat (ei, mat, tt);
  }
  
  void CompoundFESpace::VTransformMC (ElementId ei, 
				      SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const
  {
    T_TransformMat (ei, mat, tt);
  }
  

  void CompoundFESpace::VTransformVR (ElementId ei, 
				      SliceVector<double> vec, TRANSFORM_TYPE tt) const 
  {
    T_TransformVec (ei, vec, tt);
  }
  
  void CompoundFESpace::VTransformVC (ElementId ei, 
				      SliceVector<Complex> vec, TRANSFORM_TYPE tt) const 
  {
    T_TransformVec (ei, vec, tt);
  }









  CompoundFESpaceAllSame ::
  CompoundFESpaceAllSame (shared_ptr<FESpace> space, int dim, const Flags & flags,
                          bool checkflags)
    : CompoundFESpace (space->GetMeshAccess(), flags)
  {
    order = space->GetOrder();
    for (int i = 0; i < dim; i++)
      AddSpace (space);

    // interleaved: not yet implemented, will copy from VectorH1
    interleaved = flags.GetDefineFlag("interleaved");
    if (interleaved)
      throw Exception("CompoundFESpaceAllSame: Interleaved not yet supported");
    
    for (auto vb : { VOL, BND, BBND, BBBND })
      {
        if (auto eval = spaces[0] -> GetEvaluator(vb))
          evaluator[vb] = make_shared<VectorDifferentialOperator> (eval, dim);
        if (auto fluxeval = spaces[0] -> GetFluxEvaluator(vb))
            flux_evaluator[vb] = make_shared<VectorDifferentialOperator> (fluxeval, dim);
      }
    
    auto additional = spaces[0]->GetAdditionalEvaluators();
    for (int i = 0; i < additional.Size(); i++)
      additional_evaluators.Set (additional.GetName(i),
                                 make_shared<VectorDifferentialOperator>(additional[i], dim));
    
    type = "Vector"+(*this)[0]->type;
  }
  
  string CompoundFESpaceAllSame :: GetClassName () const
  {
    return "Vector"+ (*this)[0]->GetClassName();
  }


  shared_ptr<BaseMatrix> CompoundFESpaceAllSame ::
  GetTraceOperator (shared_ptr<FESpace> tracespace, bool avg) const
  {
    shared_ptr<BaseMatrix> sum;
    auto compoundtrace = dynamic_pointer_cast<CompoundFESpaceAllSame>(tracespace);
    for (int i = 0; i < Spaces().Size(); i++)
      {
        auto tracei = ComposeOperators( compoundtrace->EmbeddingOperator(i),
                                        ComposeOperators ( Spaces()[i]->GetTraceOperator(compoundtrace->Spaces()[i], avg),
                                                           this->RestrictionOperator(i)));
        if (sum)
          sum = AddOperators(sum, tracei, 1, 1);
        else
          sum = tracei;
      }
    return sum;
  }

  

  MatrixFESpace ::
  MatrixFESpace (shared_ptr<FESpace> space, int avdim, const Flags & flags,
                 bool checkflags)
    : CompoundFESpace (space->GetMeshAccess(), flags), vdim(avdim)
  {
    order = space->GetOrder();    
    symmetric = flags.GetDefineFlag("symmetric");
    deviatoric = flags.GetDefineFlag("deviatoric");
    skewsymmetric = flags.GetDefineFlag("skewsymmetric");
    dgjumps = space->UsesDGCoupling();
    
    if (deviatoric && !symmetric) throw Exception ("non-symmetric and deviatoric not supported");
    if (symmetric && skewsymmetric) throw Exception ("symmetric and skewsymmetric - not much would remain"); 
    int dim = symmetric ? vdim*(vdim+1)/2 : sqr(vdim);
    if (skewsymmetric) dim = vdim*(vdim-1)/2;
    
    if (deviatoric && !skewsymmetric) dim--;
    for (int i = 0; i < dim; i++)
      AddSpace (space);
    
    for (auto vb : { VOL, BND, BBND, BBBND })
      {
        if (auto eval = spaces[0] -> GetEvaluator(vb))
          {
            if (symmetric)
              {
                if (deviatoric)
                  evaluator[vb] = make_shared<SymDevMatrixDifferentialOperator> (eval, vdim);
                else
                  evaluator[vb] = make_shared<SymMatrixDifferentialOperator> (eval, vdim);
              }
            else if (skewsymmetric)
              evaluator[vb] = make_shared<SkewMatrixDifferentialOperator> (eval, vdim);
            else
              evaluator[vb] = make_shared<MatrixDifferentialOperator> (eval, vdim);
          }
        
        if (auto fluxeval = spaces[0] -> GetFluxEvaluator(vb))
          {
            if (!symmetric && !skewsymmetric && !deviatoric)
              flux_evaluator[vb] = make_shared<MatrixDifferentialOperator> (fluxeval, vdim);
          }
      }

    auto additional = spaces[0]->GetAdditionalEvaluators();
    if (!symmetric && !skewsymmetric && !deviatoric)
      for (int i = 0; i < additional.Size(); i++)
        additional_evaluators.Set (additional.GetName(i),
                                   make_shared<MatrixDifferentialOperator>(additional[i], vdim));
    
    type = string((symmetric) ? "Sym" : "") + string((skewsymmetric) ? "Skew" : "") + (deviatoric ? "Dev" : "") + "Matrix" + (*this)[0]->type;

    for (VorB vb : {VOL, BND, BBND, BBBND})
      definedon[vb] = space->DefinedOn(vb);
  }
  
  string MatrixFESpace :: GetClassName () const
  {
    return string(symmetric ? "Sym" : "") + string(skewsymmetric ? "Skew" : "") + (deviatoric ? "Dev" : "") + "Matrix" + (*this)[0]->GetClassName();
  }


  FiniteElement & MatrixFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    const FiniteElement &fe0 = spaces[0]->GetFE(ei, alloc);
    if (symmetric)
      return *new (alloc) SymMatrixFiniteElement(fe0, vdim, deviatoric);
    else if (skewsymmetric)
      return *new (alloc) SkewMatrixFiniteElement(fe0, vdim);
    else
      return *new (alloc) VectorFiniteElement(fe0, static_cast<int>(spaces.Size()));
  }


  


  
  ApplyMass :: ApplyMass (shared_ptr<FESpace> afes,
                          shared_ptr<CoefficientFunction> arho,
                          bool ainverse,
                          shared_ptr<Region> adefinedon,
                          LocalHeap & alh)
    : fes(afes), rho(arho), inverse(ainverse), definedon(adefinedon), lh(alh) { ; }
  
  ApplyMass :: ~ApplyMass()
  { ; } 

  shared_ptr<BaseMatrix> ApplyMass :: InverseMatrix (shared_ptr<BitArray> subset) const
  {
    return make_shared<ApplyMass> (fes, rho, !inverse, definedon, lh);
  }

 
  void ApplyMass :: Mult (const BaseVector & v, BaseVector & prod) const 
  {
    prod = v;
    if (inverse)
      fes->SolveM(rho.get(), prod, definedon.get(), lh);
    else
      fes->ApplyM(rho.get(), prod, definedon.get(), lh);
  }

  void ApplyMass :: MultAdd (double val, const BaseVector & v, BaseVector & prod) const 
  {
    auto hv = prod.CreateVector();
    hv = v;
    if (inverse)    
      fes->SolveM(rho.get(), hv, definedon.get(), lh);
    else
      fes->ApplyM(rho.get(), hv, definedon.get(), lh);
    prod += val * hv;
  }
  
  void ApplyMass :: MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const
  {
    auto hv = prod.CreateVector();
    hv = v;
    if (inverse)
      fes->SolveM(rho.get(), hv, definedon.get(), lh);
    else
      fes->ApplyM(rho.get(), hv, definedon.get(), lh);
    prod += val * hv;
  }
  
  void ApplyMass :: MultTransAdd (double val, const BaseVector & v, BaseVector & prod) const
  {
    MultAdd (val, v, prod);
  }
    
  AutoVector ApplyMass :: CreateVector () const
  {
    // should go to fespace
    return CreateBaseVector(fes->GetNDof(), fes->IsComplex(), fes->GetDimension());
  }
  
  AutoVector ApplyMass :: CreateRowVector () const
  {
    // should go to fespace    
    return CreateBaseVector(fes->GetNDof(), fes->IsComplex(), fes->GetDimension());    
  }
  
  AutoVector ApplyMass :: CreateColVector () const
  {
    // should go to fespace
    return CreateBaseVector(fes->GetNDof(), fes->IsComplex(), fes->GetDimension());    
  }



  
  ApplyTrace :: ApplyTrace (shared_ptr<FESpace> afes,
                            shared_ptr<FESpace> afestrace,
                            bool aaverage,
                            LocalHeap & alh)
    : fes(afes), festrace(afestrace), average(aaverage), lh(alh) { ; }
  
  ApplyTrace :: ~ApplyTrace() { ; } 


  void ApplyTrace :: Mult (const BaseVector & v, BaseVector & prod) const 
  {
    fes->GetTrace(*festrace, v, prod, average, lh);
  }

  void ApplyTrace :: MultAdd (double val, const BaseVector & v, BaseVector & prod) const 
  {
    auto hv = prod.CreateVector();
    fes->GetTrace(*festrace, v, hv, average, lh);
    prod += val * hv;
  }
  
  void ApplyTrace :: MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const
  {
    auto hv = prod.CreateVector();
    fes->GetTrace(*festrace, v, hv, average, lh);
    prod += val * hv;
  }
  
  void ApplyTrace :: MultTransAdd (double val, const BaseVector & v, BaseVector & prod) const
  {
    auto hv = prod.CreateVector();
    fes->GetTraceTrans(*festrace, v, hv, average, lh);
    prod += val * hv;
  }
    
  AutoVector ApplyTrace :: CreateVector () const
  {
    // should go to fespace
    return CreateBaseVector(fes->GetNDof(), fes->IsComplex(), fes->GetDimension());
  }
  
  AutoVector ApplyTrace :: CreateRowVector () const
  {
    // should go to fespace    
    return CreateBaseVector(fes->GetNDof(), fes->IsComplex(), fes->GetDimension());    
  }
  
  AutoVector ApplyTrace :: CreateColVector () const
  {
    // should go to fespace
    return CreateBaseVector(festrace->GetNDof(), fes->IsComplex(), fes->GetDimension());    
  }









  FESpaceClasses :: ~FESpaceClasses() { ; }
  
  void FESpaceClasses :: 
  AddFESpace (const string & aname,
	      shared_ptr<FESpace> (*acreator)(shared_ptr<MeshAccess> ma, const Flags & flags),
              DocInfo (*getdocu)())
  {
    fesa.Append (make_shared<FESpaceInfo> (aname, acreator, getdocu));
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

  static RegisterClassForArchive<FESpace> regfes;
  static RegisterFESpace<NodalFESpace> initnodalfes ("nodal");
  static RegisterFESpace<NonconformingFESpace> initncfes ("nonconforming");
  static RegisterFESpace<NonconformingSurfaceFESpace> initncsurfes ("nonconformingsurface");

}



