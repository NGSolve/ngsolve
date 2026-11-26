/**********************************************************************/
/* File:   fespace.cpp                                                */
/* Author: Joachim Schoeberl                                          */
/* Date:   25. Mar. 2000                                              */
/**********************************************************************/


/* 
   Finite Element Space
*/

#include <bit>
#include "fespace.hpp"
#include <multigrid.hpp>
#include <diagonalmatrix.hpp>
#include <special_matrix.hpp>

#include "../fem/h1lofe.hpp"
#include <bdbequations.hpp>
#include <parallelngs.hpp>
#include "compressedfespace.hpp"
#include "periodic.hpp"
#include "gridfunction.hpp"
#include "tpfes.hpp"


namespace ngcomp
{
  using namespace ngmg;

  // regions in flags index the region arrays of the space's own mesh
  static const Region & OwnRegion (const Region & reg, const MeshAccess & ma, const string & what)
  {
    if (reg.Mesh().get() != &ma)
      throw Exception (what + ": the region belongs to another mesh");
    return reg;
  }

  static const Region & OwnRegion (const Region & reg, const MeshAccess & ma, VorB vb, const string & what)
  {
    OwnRegion (reg, ma, what);
    if (reg.VB() != vb)
      throw Exception (what + ": needs a " + ToString(vb) + " region, got a " + ToString(reg.VB()) + " region");
    return reg;
  }

  static const MixedRegion & OwnRegion (const MixedRegion & reg, const MeshAccess & ma, const string & what)
  {
    if (reg.Mesh().get() != &ma)
      throw Exception (what + ": the region belongs to another mesh");
    for (auto vb : { VOL, BND, BBND, BBBND })
      if (reg.Mask(vb).Size() != ma.GetNRegions(vb))
        throw Exception (what + ": " + ToString(vb) + " mask does not fit the mesh, mesh changed after creating the region?");
    return reg;
  }

  // a MixedRegion flag, or a region descriptor flag resolved on the mesh of the space
  static optional<MixedRegion> MixedRegionFlag (const Flags & flags, const string & name,
                                                const shared_ptr<MeshAccess> & ma)
  {
    if (!flags.AnyFlagDefined(name)) return nullopt;
    const std::any & val = flags.GetAnyFlag(name);
    if (auto reg = std::any_cast<MixedRegion>(&val))
      return OwnRegion(*reg, *ma, name);
    if (auto reg = std::any_cast<Region>(&val))
      return OwnRegion(MixedRegion(*reg), *ma, name);
    if (auto rd = std::any_cast<MixedRegionDescriptor>(&val))
      return MixedRegion(ma, *rd);
    if (auto rd = std::any_cast<RegionDescriptor>(&val))
      return MixedRegion(ma, MixedRegionDescriptor(*rd));
    return nullopt;
  }

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

    if (flags.GetNumFlag ("order", 1) < 0)
      throw Exception ("FESpace: order must be non-negative, got " + ToString(flags.GetNumFlag ("order", 1)));
    CheckOrder (flags.GetNumFlag ("order", 1));
    order = int (flags.GetNumFlag ("order", 1));
    order_spec.base->order = order;
    order_spec.base->creator = this;

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
    
    if (!(flags.GetNumFlag ("dim", 1) >= 1 && flags.GetNumFlag ("dim", 1) <= max_dim))
      throw Exception ("dim must be between 1 and " + ToString(max_dim) + ", got " + ToString(flags.GetNumFlag ("dim", 1)));
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
    
    dirichlet_region.Mask(BND).SetSize (ma->GetNBoundaries());
    dirichlet_region.Mask(BND).Clear();

    if(flags.NumListFlagDefined("dirichlet"))
      {
        for (double dbi : flags.GetNumListFlag("dirichlet"))
          {
	    int bnd = int(dbi-1);
	    if (bnd >= 0 && bnd < dirichlet_boundaries.Size())
	      dirichlet_region.Mask(BND).SetBit (bnd);
	    // else
            //   cerr << "Illegal Dirichlet boundary index " << bnd+1 << endl;
          }
	if (print)
	  *testout << "dirichlet_boundaries:" << endl << dirichlet_boundaries << endl;
      }

    if (flags.StringFlagDefined("dirichlet"))
      {
        Region dir(ma, BND, flags.GetStringFlag("dirichlet"));
        dirichlet_region.Mask(BND) |= dir.Mask(); // that's what it was before, can we just assign ? 
      }

    dirichlet_region.Mask(BBND).SetSize (ma->GetNRegions(BBND));
    dirichlet_region.Mask(BBND).Clear();

    if(flags.AnyFlagDefined("dirichlet_bbnd"))
      {
        dirichlet_region.Mask(BBND) = OwnRegion(std::any_cast<const Region&>(flags.GetAnyFlag("dirichlet_bbnd")), *ma, BBND, "dirichlet_bbnd").Mask();
      }

    if (flags.StringFlagDefined("dirichlet_bbnd"))
      {
        Region dir(ma, BBND, flags.GetStringFlag("dirichlet_bbnd"));
        dirichlet_region.Mask(BBND) |= dir.Mask(); // that's what it was before, can we just assign ? 
      }

    dirichlet_region.Mask(BBBND).SetSize (ma->GetNRegions(BBBND));
    dirichlet_region.Mask(BBBND).Clear();

    if(flags.AnyFlagDefined("dirichlet_bbbnd"))
      {
        dirichlet_region.Mask(BBBND) = OwnRegion(std::any_cast<const Region&>(flags.GetAnyFlag("dirichlet_bbbnd")), *ma, BBBND, "dirichlet_bbbnd").Mask();
      }

    if (flags.StringFlagDefined("dirichlet_bbbnd"))
      {
        Region dir(ma, BBBND, flags.GetStringFlag("dirichlet_bbbnd"));
        dirichlet_region.Mask(BBBND) |= dir.Mask(); // that's what it was before, can we just assign ? 
      }

    // Region, MixedRegion or region descriptor: ORed with the per-VorB flags above
    if(auto dir = MixedRegionFlag(flags, "dirichlet", ma))
      {
        if (dir->Mask(VOL).NumSet())
          throw Exception ("dirichlet: VOL regions are not supported, use BND, BBND or BBBND regions");
        for (auto vb : { BND, BBND, BBBND })
          dirichlet_region.Mask(vb).Or (dir->Mask(vb));
      }
    else if(flags.AnyFlagDefined("dirichlet"))
      throw Exception ("dirichlet with wrong argument type");
    
    Array<bool> defon_given[4];
    if (flags.NumListFlagDefined("definedon") || 
        flags.NumFlagDefined("definedon") ||
        flags.StringListFlagDefined("definedon") ||
        flags.AnyFlagDefined("definedon"))
      {
	defon_given[VOL].SetSize (ma->GetNDomains());
	defon_given[VOL] = false;
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
	    defon_given[int(defon[i])-1] = true;
        */
        for (int di : defon)
          if (di > 0 && di <= ma->GetNDomains())
            defon_given[VOL][di-1] = true;
          
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
		      defon_given[VOL][i] = true;
		      break;
		    }
	      }
	  }

	// default:
	// fespace only defined on boundaries matching definedon-domains
	defon_given[BND].SetSize (ma->GetNBoundaries());
	defon_given[BND] = false;

        if(flags.AnyFlagDefined("definedon"))
          {
            std::any anydefon = flags.GetAnyFlag("definedon");
            if(anydefon.type() == typeid(std::map<VorB, Region>))
              {
                auto map = std::any_cast<std::map<VorB, Region>>(anydefon);
                for(const auto& [vb, reg] : map)
                  {
                    OwnRegion (reg, *ma, vb, "definedon");
                    if(vb == BBND || vb == BBBND)
                      defon_given[vb].SetSize(ma->GetNRegions(vb));
                    for(auto i : Range(defon_given[vb]))
                      defon_given[vb][i] = reg.Mask().Test(i);
                  }
                // if vol is not in map then definedon vol = empty
                if(map.find(VOL) == map.end())
                  defon_given[VOL] = false;
              }
            else if(anydefon.type() == typeid(Region))
              {
                auto reg = OwnRegion(std::any_cast<Region>(anydefon), *ma, "definedon");
                defon_given[reg.VB()].SetSize(reg.Mask().Size());
                for(auto i : Range(defon_given[reg.VB()]))
                  defon_given[reg.VB()][i] = reg.Mask().Test(i);
                for(auto vb = VOL; vb < reg.VB(); vb++)
                  defon_given[vb] = false;
              }
            else if(auto reg = MixedRegionFlag(flags, "definedon", ma))
              {
                // like the {VOL: ..., BND: ...} dict; empty BBND/BBBND parts do not restrict
                for (auto vb : { VOL, BND, BBND, BBBND })
                  {
                    if (vb >= BBND && !reg->Mask(vb).NumSet()) continue;
                    defon_given[vb].SetSize(ma->GetNRegions(vb));
                    for (auto i : Range(defon_given[vb]))
                      defon_given[vb][i] = reg->Mask(vb).Test(i);
                  }
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
	      if ( defon_given[VOL][dom1] )
		defon_given[BND][index] = true;

	    if ( dom2 >= 0 )
	      if ( defon_given[VOL][dom2] )
		defon_given[BND][index] = true;
	  }

        // closed: defined on the sub-entities of the given regions
        MixedRegion given(ma);
        for (auto vb : { VOL, BND, BBND, BBBND })
          for (auto i : Range(min(defon_given[vb].Size(), given.Mask(vb).Size())))
            if (defon_given[vb][i])
              given.Mask(vb).SetBit(i);
        MixedRegion closed = given.Closure();
        for (auto vb : { BND, BBND, BBBND })
          {
            Array<int> on(closed.Mask(vb).Size());
            for (auto i : Range(on))
              on[i] = closed.Mask(vb).Test(i);
            // region incidence is known per rank only
            if (on.Size())
              ma->GetCommunicator().AllReduce (on, NG_MPI_MAX);
            defon_given[vb].SetSize(on.Size());
            for (auto i : Range(on))
              defon_given[vb][i] = on[i];
          }
      }

    if(flags.AnyFlagDefined("definedonbound"))
        throw Exception("definedonbound with Region is not supported, use definedon= {VOL: mats, BND: bnds}");

    // additional boundaries
    if(flags.NumListFlagDefined("definedonbound")|| flags.NumFlagDefined("definedonbound") )
      {
	if ( defon_given[BND].Size() == 0 )
	  {
	    defon_given[BND].SetSize (ma->GetNBoundaries());
	    defon_given[BND] = false;
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
	    defon_given[BND][int(defon[i])-1] = true;
      }
    

    else if(flags.StringListFlagDefined("definedonbound") || flags.StringFlagDefined("definedonbound"))
      {
	if ( defon_given[BND].Size() == 0 )
	  {
	    defon_given[BND].SetSize (ma->GetNBoundaries());
	    defon_given[BND] = false;
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
	    if(defon_given[BND][ma->GetElIndex(sei)] == false)
	      {
		for(int i=0; i<defon.Size(); i++)
		  {
		    // if(StringFitsPattern(ma->GetSElBCName(selnum),*(defon[i])))
                    if(StringFitsPattern(ma->GetMaterial(ElementId(BND, selnum)),*(defon[i])))	
		      {		
		 	defon_given[BND][ma->GetElIndex(sei)] = true;
			continue;
		      }
		  }
	      }
	  }
	for(int i=0; i<defon.Size(); i++)
	  delete defon[i];
      }

    for (auto vb : { VOL, BND, BBND, BBBND })
      if (defon_given[vb].Size())
        {
          auto & mask = definedon_region.Mask(vb);
          mask.Clear();
          for (auto i : Range(min(mask.Size(), defon_given[vb].Size())))
            if (defon_given[vb][i]) mask.SetBit(i);
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
          {
            CheckOrder (flags.GetNumFlag(flagname, 0));
            SetOrder (et, int(flags.GetNumFlag(flagname, 0)));
          }
      }


    
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
    docu.Arg("dirichlet") = "regexpr or Region\n"
      "  Regular expression string defining the dirichlet boundary.\n"
      "  More than one boundary can be combined by the | operator,\n"
      "  i.e.: dirichlet = 'top|right'\n"
      "  A region of BND, BBND or BBBND, also of several codimensions:\n"
      "  dirichlet = mesh.Boundaries('top') | mesh.BBoundaries('edge'),\n"
      "  or without mesh: dirichlet = BND('top') | BBND('edge')";
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
      "  mesh.Materials('regexpr'). Regions of several codimensions can be combined:\n"
      "  definedon = mesh.Materials('a') | mesh.Boundaries('b'). The space is also\n"
      "  defined on all sub-entities (faces, edges, vertices) of the given regions.";
    docu.Arg("definedonbound") = "list of int or list of str\n"
      "  Additional boundaries the FESpace is defined on, besides the boundaries\n"
      "  adjacent to the definedon domains. Given as 1-based boundary numbers or\n"
      "  wildcard patterns ('*', '?') of boundary names. Region is not supported,\n"
      "  use definedon = {VOL: mats, BND: bnds} instead.";
    docu.Arg("dim") = "int = 1\n"
      "  Create multi dimensional FESpace (i.e. [H1]^3)";
    docu.Arg("dgjumps") = "bool = False\n"
      "  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n"
      "  since the dofs have a different coupling then and this changes the sparsity\n"
      "  pattern of matrices.";
    docu.Arg("autoupdate") = "bool = True\n"
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

  void FESpace :: SetHarmonicProlongation (shared_ptr<BilinearForm> bfa, string inverse)
  {
    if (!bfa)
      throw Exception ("SetHarmonicProlongation: bilinear form is None");
    if (bfa->GetFESpace().get() != this)
      throw Exception ("SetHarmonicProlongation: the bilinear form must be defined on this space");
    prol = make_shared<HarmonicProlongation> (prol, bfa, inverse);
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
    order_stale = 0;
    order_base_resolved = order_spec.base->timestamp;

    for (auto vb : { VOL, BND, BBND, BBBND })
      if (auto & mask = definedon_region.Mask(vb); mask.Size() != ma->GetNRegions(vb))
        {
          if (!DefinedOnAll(vb))
            throw Exception("FESpace: number of regions of the mesh changed, definedon is not valid anymore");
          mask.SetSize(ma->GetNRegions(vb));
          mask.Set();
        }
    
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
        auto & dc = dirichlet_region.Mask(vb);
        if (dc.Size())
          for (Ngs_Element ngel : ma->Elements(vb)) 
            if (DefinedAt(dc, ngel.GetIndex()))
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


  // blocks of B elements sorted by smallest vertex number, greedy coloring of blocks
  static void BuildBlockColoring (const FESpace & fes, const MeshAccess & ma, VorB vb,
                                  Table<int> & blocks, Table<int> & coloring)
  {
    static Timer t("FESpace - block coloring"); RegionTimer reg(t);
    Array<int> els, key;
    for (size_t nr = 0; nr < ma.GetNE(vb); nr++)
      {
        ElementId ei(vb, nr);
        if (!fes.DefinedOn(ei)) continue;
        int minv = numeric_limits<int>::max();
        for (auto v : ma.GetElVertices(ei)) minv = min(minv, int(v));
        els.Append(nr);
        key.Append(minv);
      }
    Array<int> index(els.Size());
    for (int i = 0; i < index.Size(); i++) index[i] = i;
    QuickSortI (key, index);

    // smaller blocks on small meshes, to keep threads busy
    size_t B = clamp<size_t> (els.Size() / (64*max(1, TaskManager::GetMaxThreads())), 16, 256);
    size_t nb = (els.Size()+B-1) / B;
    Array<int> bsize(nb);
    for (size_t b = 0; b < nb; b++) bsize[b] = min(B, els.Size()-b*B);
    blocks = Table<int>(bsize);
    for (size_t i = 0; i < els.Size(); i++)
      blocks[i/B][i%B] = els[index[i]];

    Array<int> bcol(nb);
    bcol = -1;
    Array<uint64_t> mask(fes.GetNDof());
    Array<DofId> dofs;
    size_t found = 0;
    int basecol = 0, maxcolor = -1;
    while (found < nb)
      {
        mask = 0;
        for (size_t b = 0; b < nb; b++)
          {
            if (bcol[b] >= 0) continue;
            auto block_dofs = [&] (auto func)
            {
              for (int el : blocks[b])
                {
                  fes.GetDofNrs(ElementId(vb, el), dofs);
                  for (auto d : dofs)
                    if (IsRegularDof(d) && !fes.IsAtomicDof(d)) func(d);
                }
            };
            uint64_t check = 0;
            block_dofs ([&] (DofId d) { check |= mask[d]; });
            if (check == ~uint64_t(0)) continue;
            int bit = std::countr_zero(~check);
            block_dofs ([&] (DofId d) { mask[d] |= uint64_t(1) << bit; });
            bcol[b] = basecol + bit;
            maxcolor = max(maxcolor, bcol[b]);
            found++;
          }
        basecol += 64;
      }

    Array<int> cntcol(maxcolor+1);
    cntcol = 0;
    for (int c : bcol) cntcol[c]++;
    coloring = Table<int>(cntcol);
    cntcol = 0;
    for (size_t b = 0; b < nb; b++)
      coloring[bcol[b]][cntcol[bcol[b]]++] = b;
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

    UpdateFreeDofs();
    
    
    UpdateParallelDofs();

    if (print)
      *testout << "coloring ... " << flush;

    if (low_order_space)
      {
	for(auto vb : {VOL, BND, BBND, BBBND})
	  element_coloring[vb] = Table<int>(low_order_space->element_coloring[vb]);
        element_blocks[VOL] = Table<int>(low_order_space->element_blocks[VOL]);
        block_coloring[VOL] = Table<int>(low_order_space->block_coloring[VOL]);
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
      BuildBlockColoring (*this, *ma, VOL, element_blocks[VOL], block_coloring[VOL]);
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


  void FESpace :: UpdateFreeDofs()
  {
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
      if (auto & dc = dirichlet_region.Mask(vb); dc.Size())
        for (FESpace::Element el : Elements(vb))
          if (DefinedAt(dc, el.GetIndex()))
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



  void FESpace :: CheckOrder (double order)
  {
    if (!(order <= max_order))
      throw Exception ("FESpace: order " + ToString(order) + " exceeds the maximal order " + ToString(max_order));
  }

  void FESpace :: SetOrder (ELEMENT_TYPE et, TORDER order)
  {
    if (order_policy == CONSTANT_ORDER || order_policy == OLDSTYLE_ORDER)
      order_policy = NODE_TYPE_ORDER;
    et_bonus_order[et] = order - this->order;
    
    string flagname = string(ElementTopology::GetElementName(et)) + "order";
    std::transform(flagname.begin(), flagname.end(), flagname.begin(), ::tolower);        
    flags.SetFlag(flagname, order);
    
    timestamp = 0;  // rerun first_update
    order_stale = 2;
  }


  
  void FESpace :: SetOrder (NodeId ni, int order)
  {
    throw Exception (string("FESpace::SetOrder not overloaded for space") + typeid(*this).name());
  }

  int FESpace :: GetOrder (NodeId ni) const
  {
    throw Exception (string("FESpace::GetOrder not overloaded for space") + typeid(*this).name());
  }
    


  template <typename TO>
  void FESpace :: ResolveOrders (int minorder,
                                 FlatArray<bool> used_edge, FlatArray<bool> used_face,
                                 FlatArray<TO> order_edge,
                                 FlatArray<IVec<2,TO>> order_face,
                                 FlatArray<IVec<3,TO>> order_cell,
                                 int edge_offset, bool increments_everywhere, bool use_increments) const
  {
    int dim = ma->GetDimension();
    static const int no_bonus[30] = { 0 };
    const int * bonus = use_increments ? et_bonus_order : no_bonus;
    int p = order_spec.base->mesh_orders ? minorder : order;

    for (auto & ro : order_spec.base->region_orders)
      if (ro.mask.Size() != ma->GetNRegions(ro.vb))
        throw Exception ("FESpace: number of regions of the mesh changed, region orders are not valid anymore");

    // minimum rule: an edge or face takes the lowest order of its cells
    const TO unset = std::numeric_limits<TO>::max();
    order_edge = unset;
    order_face = IVec<2,TO> (unset);
    order_cell = IVec<3,TO> (p);

    ParallelFor (ma->GetNE(VOL), [&] (size_t nr)
                 {
                   Ngs_Element el = (*ma)[ElementId(VOL, nr)];
                   if (!DefinedOn (el)) return;

                   bool elementwise;
                   IVec<3> pc = CellBaseOrder (el, el.GetIndex(), elementwise);
                   // derived orders are not below the lower bound of the space
                   int lo = (elementwise || order_spec.link_offset) ? minorder : 0;
                   ELEMENT_TYPE eltype = el.GetType();

                   for (int k = 0; k < 3; k++)
                     order_cell[nr][k] = max2 (lo, pc[k] + bonus[eltype]);

                   if (pc[0] == pc[1] && pc[1] == pc[2])
                     {
                       if (order_edge.Size())
                         for (auto e : el.Edges())
                           AtomicMin (order_edge[e], TO (max2 (lo, pc[0] + edge_offset + bonus[ET_SEGM])));

                       if (order_face.Size())
                         for (auto f : el.Faces())
                           {
                             TO pf = max2 (lo, pc[0] + bonus[ma->GetFaceType(f)]);
                             AtomicMin (order_face[f][0], pf);
                             AtomicMin (order_face[f][1], pf);
                           }
                       return;
                     }

                   // different orders in the local directions
                   auto faces = ElementTopology::GetFaces (eltype);
                   auto edges = ElementTopology::GetEdges (eltype);
                   const POINT3D * points = ElementTopology::GetVertices (eltype);
                   auto vnums = el.Vertices();

                   if (order_edge.Size())
                     {
                       auto eledges = el.Edges();
                       for (int j = 0; j < eledges.Size(); j++)
                         for (int k = 0; k < dim; k++)
                           if (points[edges[j][0]][k] != points[edges[j][1]][k]) // find edge-dir x,y,z
                             {
                               AtomicMin (order_edge[eledges[j]], TO (max2 (lo, pc[k] + edge_offset + bonus[ET_SEGM])));
                               break;
                             }
                     }

                   if (order_face.Size())
                     {
                       auto elfaces = el.Faces();
                       for (int j = 0; j < elfaces.Size(); j++)
                         {
                           auto & of = order_face[elfaces[j]];
                           if (faces[j][3] == -1)
                             {
                               TO pf = max2 (lo, pc[0] + bonus[ET_TRIG]);
                               AtomicMin (of[0], pf);
                               AtomicMin (of[1], pf);
                             }
                           else
                             {
                               // face directions by global vertex numbers
                               int fmax = 0;
                               for (int k = 1; k < 4; k++)
                                 if (vnums[faces[j][k]] > vnums[faces[j][fmax]]) fmax = k;

                               IVec<2> f((fmax+3)%4, (fmax+1)%4);
                               if (vnums[faces[j][f[1]]] > vnums[faces[j][f[0]]]) swap (f[0], f[1]);

                               // fmax > f[0] > f[1]
                               // p[0] for direction fmax,f[0]
                               // p[1] for direction fmax,f[1]
                               for (int l = 0; l < 2; l++)
                                 for (int k = 0; k < 3; k++)
                                   if (points[faces[j][fmax]][k] != points[faces[j][f[l]]][k])
                                     {
                                       AtomicMin (of[l], TO (max2 (lo, pc[k] + bonus[ET_QUAD])));
                                       break;
                                     }
                             }
                         }
                     }
                 });

    if ((order_spec.base->region_orders.Size() || order_spec.base->mesh_orders || HasElementOrders() || order_spec.base->hp_slope > 0)
        && ma->GetCommunicator().Size() > 1)
      {
        // cells of other ranks
        Array<int> hedge (order_edge.Size()), hface (order_face.Size());
        for (auto i : Range(hedge)) hedge[i] = order_edge[i];
        if (hedge.Size())
          ma->AllReduceNodalData (NT_EDGE, hedge, NG_MPI_MIN);
        for (auto i : Range(hedge)) order_edge[i] = TO (hedge[i]);
        for (int k = 0; k < 2 && hface.Size(); k++)
          {
            for (auto i : Range(hface)) hface[i] = order_face[i][k];
            ma->AllReduceNodalData (NT_FACE, hface, NG_MPI_MIN);
            for (auto i : Range(hface)) order_face[i][k] = TO (hface[i]);
          }
      }

    // edges and faces without a cell of the space
    for (auto & oe : order_edge)
      if (oe == unset)
        oe = increments_everywhere ? TO (max2 (0, p + edge_offset + bonus[ET_SEGM])) : TO (p);
    for (auto i : Range(order_face))
      if (order_face[i][0] == unset)
        order_face[i] = IVec<2,TO> (increments_everywhere ? max2 (0, p + bonus[ma->GetFaceType(i)]) : p);

    // faces and edges of boundary regions, the last region wins
    bool bnd_regions = false;
    for (auto & ro : order_spec.base->region_orders)
      if (ro.vb != VOL) bnd_regions = true;
    if (bnd_regions)
      {
        auto & regs = order_spec.base->region_orders;
        Array<int> edge_reg (order_edge.Size()), face_reg (order_face.Size());
        edge_reg = -1;
        face_reg = -1;
        for (int k : Range(regs))
          if (regs[k].vb == BND || regs[k].vb == BBND)
            for (Ngs_Element el : ma->Elements(regs[k].vb))
              {
                if (el.GetType() == ET_POINT || !regs[k].mask.Test(el.GetIndex())) continue;
                if (edge_reg.Size())
                  for (auto e : el.Edges())
                    edge_reg[e] = k;
                if (face_reg.Size() && regs[k].vb == BND)
                  face_reg[ma->GetSElFace(el.Nr())] = k;
              }
        if (edge_reg.Size())
          ma->AllReduceNodalData (NT_EDGE, edge_reg, NG_MPI_MAX);
        int lobnd = order_spec.link_offset ? minorder : 0;

        for (auto i : Range(edge_reg))
          if (edge_reg[i] >= 0)
            order_edge[i] = TO (max2 (lobnd, regs[edge_reg[i]].order + order_spec.link_offset + edge_offset + bonus[ET_SEGM]));
        for (auto i : Range(face_reg))
          if (face_reg[i] >= 0)
            order_face[i] = IVec<2,TO> (max2 (lobnd, regs[face_reg[i]].order + order_spec.link_offset + bonus[ma->GetFaceType(i)]));
      }

    if (order_spec.order_cell > -1)
      order_cell = IVec<3,TO> (order_spec.order_cell);
    if (order_spec.order_face > -1 && dim == 3)
      order_face = IVec<2,TO> (order_spec.order_face);
    if (order_spec.order_edge > -1)
      order_edge = TO (order_spec.order_edge);

    for (auto i : Range(order_edge))
      if (!used_edge[i]) order_edge[i] = minorder;

    for (auto i : Range(order_face))
      if (!used_face[i]) order_face[i] = IVec<2,TO> (minorder);

    for (ElementId ei : ma->Elements<VOL>())
      if (!DefinedOn(ei)) order_cell[ei.Nr()] = IVec<3,TO> (minorder);

    if (HasNodeOrders())
      {
        auto & no = order_spec.node_order;
        if (no[0].Size() == order_edge.Size())
          for (auto i : Range(order_edge))
            if (no[0][i] >= 0) order_edge[i] = TO (no[0][i]);
        if (no[1].Size() == order_face.Size())
          for (auto i : Range(order_face))
            if (no[1][i] >= 0) order_face[i] = IVec<2,TO> (no[1][i]);
        if (no[2].Size() == order_cell.Size())
          for (auto i : Range(order_cell))
            if (no[2][i] >= 0) order_cell[i] = IVec<3,TO> (no[2][i]);
      }
  }

  template void FESpace :: ResolveOrders<unsigned char> (int, FlatArray<bool>, FlatArray<bool>, FlatArray<unsigned char>,
                                                       FlatArray<IVec<2,unsigned char>>, FlatArray<IVec<3,unsigned char>>,
                                                       int, bool, bool) const;
  template void FESpace :: ResolveOrders<short> (int, FlatArray<bool>, FlatArray<bool>, FlatArray<short>,
                                               FlatArray<IVec<2,short>>, FlatArray<IVec<3,short>>,
                                               int, bool, bool) const;
  template void FESpace :: ResolveOrders<int> (int, FlatArray<bool>, FlatArray<bool>, FlatArray<int>,
                                             FlatArray<IVec<2,int>>, FlatArray<IVec<3,int>>,
                                             int, bool, bool) const;

  void FESpace :: ResolveFacetOrders (FlatArray<bool> fine_facet, FlatArray<IVec<2>> order_facet, int order_1d)
  {
    int dim = ma->GetDimension();
    UpdateOrderSpecification();
    Array<IVec<3>> order_cell (ma->GetNE());
    FlatArray<bool> no_used;
    FlatArray<int> no_edges;
    FlatArray<IVec<2>> no_faces;
    if (dim == 3)
      {
        ResolveOrders<int> (0, no_used, fine_facet, no_edges, order_facet, order_cell, 0, false, false);
        for (auto & of : order_facet)
          of = IVec<2> (min2 (of[0], of[1]));
      }
    else if (dim == 2)
      {
        Array<int> order_edge (order_facet.Size());
        ResolveOrders<int> (0, fine_facet, no_used, order_edge, no_faces, order_cell, 0, false, false);
        for (auto i : Range(order_facet))
          order_facet[i] = IVec<2> (order_edge[i]);
      }
    else
      for (auto i : Range(order_facet))
        order_facet[i] = IVec<2> (fine_facet[i] ? order_1d : 0);
  }

  bool FESpace :: InitBaseOrders (const Flags & flags, bool mesh_orders, int rel_order)
  {
    if (LinkOrders (flags)) return true;
    auto & base = *order_spec.base;
    base.order = order;
    base.mesh_orders = mesh_orders;
    base.rel_order = mesh_orders ? rel_order : 0;
    base.hp_slope = flags.GetNumFlag ("hpslope", 0);
    base.hp_minorder = 0;
    if (base.hp_slope < 0)
      throw Exception (type + ": hpslope must not be negative");
    if (base.hp_slope > 0 && mesh_orders)
      throw Exception (type + ": hpslope not together with variableorder / relorder");
    return false;
  }

  void FESpace :: SetLinkFlags (Flags & flags, int offset) const
  {
    auto & base = *order_spec.base;
    base.has_links = true;
    if (base.creator_name.empty())
      base.creator_name = GetClassName();
    flags.SetFlag ("orderlink", std::any (order_spec.base));
    flags.SetFlag ("relorder", double (order_spec.link_offset + offset));
  }

  void FESpace :: SetNodeOrder (int type, size_t nr, int order)
  {
    auto & spec = order_spec;
    if (spec.node_order_level != ma->GetNLevels())
      {
        for (auto & a : spec.node_order) a.SetSize(0);
        spec.node_order_level = ma->GetNLevels();
      }
    size_t n = (type == 0) ? ma->GetNEdges() : (type == 1) ? ma->GetNFaces() : ma->GetNE();
    auto & a = spec.node_order[type];
    if (a.Size() != n)
      {
        a.SetSize (n);
        a = -1;
      }
    if (nr < n) a[nr] = order;
  }

  void FESpace :: SetRegionOrder (const Region & region, int order, int order_dir,
                                  shared_ptr<CoefficientFunction> direction)
  {
    CheckOrderCreator();
    if (region.Mesh() != ma)
      throw Exception ("SetOrder: region belongs to a different mesh");
    if (region.VB() == BBBND)
      throw Exception ("SetOrder: no order for regions of points");
    if (region.VB() == VOL && order_spec.base->mesh_orders)
      throw Exception ("SetOrder: no order for volume regions together with variableorder / relorder");
    if (order_dir >= 0)
      {
        if (region.VB() != VOL)
          throw Exception ("SetOrder: different orders in different directions only for volume regions");
        if (direction && direction->Dimension() != ma->GetDimension())
          throw Exception ("SetOrder: the direction needs one component per space dimension");
        for (Ngs_Element el : ma->Elements<VOL>())
          if (region.Mask().Test(el.GetIndex()))
            switch (el.GetType())
              {
              case ET_PRISM:
                break;
              case ET_QUAD: case ET_HEX:
                if (!direction)
                  throw Exception ("SetOrder: quadrilaterals and hexahedra need a direction for the second order");
                break;
              default:
                throw Exception ("SetOrder: different orders in different directions only for quadrilaterals, hexahedra and prisms");
              }
      }

    auto & regs = order_spec.base->region_orders;
    // a repeated call for the same region replaces its entry
    for (size_t i = 0; i < regs.Size(); i++)
      if (regs[i].vb == region.VB() && regs[i].mask == region.Mask())
        {
          regs.RemoveElement (i);
          break;
        }
    regs.Append (RegionOrder { region.VB(), BitArray(region.Mask()), order, order_dir, direction });

    order_spec.base->timestamp = GetNextTimeStamp();
    timestamp = 0;
    order_stale = 2;
  }

  IVec<3> FESpace :: CellBaseOrder (ElementId ei, int index, bool & elementwise) const
  {
    IVec<3> o = UnlinkedCellBaseOrder (ei, index, elementwise);
    for (int k = 0; k < 3; k++)
      o[k] += order_spec.link_offset;
    return o;
  }

  IVec<3> FESpace :: UnlinkedCellBaseOrder (ElementId ei, int index, bool & elementwise) const
  {
    elementwise = true;
    auto & eo = order_spec.base->element_order;
    if (ei.Nr() < eo.Size() && eo[ei.Nr()][0] >= 0)
      return IVec<3> (int(eo[ei.Nr()][0]), int(eo[ei.Nr()][1]), int(eo[ei.Nr()][2]));

    if (order_spec.base->mesh_orders)
      {
        // SZ: fuer EE damit ich mit rel_order auch p-ref machen kann :)
        // no element orders in 1D meshes
        IVec<3> o = (ma->GetDimension() >= 2) ? ma->GetElOrders(ei.Nr()) : IVec<3>(1);
        for (int k = 0; k < 3; k++)
          o[k] += order_spec.base->rel_order;
        return o;
      }

    elementwise = false;
    auto & rco = order_spec.base->region_cell_order;
    if (ei.Nr() < rco.Size() && rco[ei.Nr()][0] >= 0)
      return IVec<3> (int(rco[ei.Nr()][0]), int(rco[ei.Nr()][1]), int(rco[ei.Nr()][2]));

    int p = order_spec.base->order;
    for (auto & ro : order_spec.base->region_orders)
      if (ro.vb == VOL && ro.mask.Test(index))
        p = ro.order;

    if (order_spec.base->hp_slope > 0 && ma->GetDimension() >= 2)
      {
        // the mesh stores maxorder - level in each local direction
        elementwise = true;
        IVec<3> mo = ma->GetElOrders(ei.Nr());
        if (ma->GetDimension() == 2) mo[2] = mo[1];
        IVec<3> o;
        for (int k = 0; k < 3; k++)
          {
            int level = max2 (0, order_spec.base->hp_maxmeshorder - mo[k]);
            o[k] = max2 (order_spec.base->hp_minorder, p - int (floor (order_spec.base->hp_slope*level + 1e-9)));
          }
        return o;
      }

    return IVec<3> (p);
  }

  template <int D>
  static int AlignedLocalDirection (const MeshAccess & ma, ElementId ei,
                                    const CoefficientFunction & direction, LocalHeap & lh)
  {
    IntegrationPoint ip (0.5, 0.5, (D == 3) ? 0.5 : 0);
    ElementTransformation & trafo = ma.GetTrafo (ei, lh);
    MappedIntegrationPoint<D,D> mip (ip, trafo);
    Vec<D> dir;
    direction.Evaluate (mip, FlatVector<> (D, &dir(0)));
    Mat<D,D> jac = mip.GetJacobian();
    int best = 0;
    double bestval = -1;
    for (int k = 0; k < D; k++)
      {
        double sum = 0, len = 0;
        for (int i = 0; i < D; i++)
          {
            sum += dir(i) * jac(i,k);
            len += sqr (jac(i,k));
          }
        double val = fabs(sum) / sqrt(len);
        if (val > bestval) { bestval = val; best = k; }
      }
    return best;
  }

  void FESpace :: UpdateOrderSpecification ()
  {
    auto & spec = order_spec;
    // linked spaces: done by the space which updates first
    if (spec.base->prepared_mesh == ma->GetTimeStamp() && spec.base->prepared_spec == spec.base->timestamp)
      return;
    spec.base->prepared_mesh = ma->GetTimeStamp();
    spec.base->prepared_spec = spec.base->timestamp;

    InheritElementOrders();
    size_t ne = ma->GetNE();

    if (spec.base->hp_slope > 0 && ma->GetDimension() >= 2)
      {
        int maxorder = 1;
        for (size_t i = 0; i < ne; i++)
          {
            IVec<3> mo = ma->GetElOrders(i);
            for (int k = 0; k < ma->GetDimension(); k++)
              maxorder = max2 (maxorder, mo[k]);
          }
        spec.base->hp_maxmeshorder = ma->GetCommunicator().AllReduce (maxorder, NG_MPI_MAX);
      }

    bool directions = false;
    for (auto & ro : spec.base->region_orders)
      if (ro.order_dir >= 0) directions = true;
    spec.base->region_cell_order.SetSize (directions ? ne : 0);
    if (!directions) return;

    LocalHeap lh (100000, "orderspec");
    for (Ngs_Element el : ma->Elements<VOL>())
      {
        HeapReset hr(lh);
        auto & rco = spec.base->region_cell_order[el.Nr()];
        rco = IVec<3,short> (-1);

        const RegionOrder * last = nullptr;
        for (auto & ro : spec.base->region_orders)
          if (ro.vb == VOL && ro.mask.Test(el.GetIndex()))
            last = &ro;
        if (!last || last->order_dir < 0) continue;

        short pt = last->order, pn = last->order_dir;
        switch (el.GetType())
          {
          case ET_PRISM:
            rco = IVec<3,short> (pt, pt, pn);
            break;
          case ET_QUAD:
            {
              int k = AlignedLocalDirection<2> (*ma, el, *last->direction, lh);
              rco = (k == 0) ? IVec<3,short> (pn, pt, pt) : IVec<3,short> (pt, pn, pn);
              break;
            }
          case ET_HEX:
            {
              int k = AlignedLocalDirection<3> (*ma, el, *last->direction, lh);
              rco = IVec<3,short> (pt);
              rco[k] = pn;
              break;
            }
          default:
            throw Exception ("order with direction only for quadrilaterals, hexahedra and prisms");
          }
      }
  }

  void FESpace :: InheritElementOrders ()
  {
    auto & eo = order_spec.base->element_order;
    size_t n0 = eo.Size(), n1 = ma->GetNE();
    if (n0 == 0 || n0 == n1) return;

    // a child with a new number has its parent, a refined element keeps its number
    Array<IVec<3,short>> neo(n1);
    bool ok = n1 > n0;
    for (size_t i = 0; ok && i < n1; i++)
      {
        size_t a = i;
        while (a >= n0)
          {
            a = ma->GetParentElement (ElementId(VOL, a)).Nr();
            if (a >= n1) { ok = false; break; }
          }
        if (ok) neo[i] = eo[a];
      }

    if (ok)
      eo = std::move(neo);
    else
      {
        eo.SetSize(0);
        cerr << "WARNING: " << GetClassName() << ": element orders are dropped, "
             << "this refinement provides no parent elements" << endl;
      }
  }

  void FESpace :: SetCellOrder (ElementId ei, IVec<3> o)
  {
    CheckOrderCreator();
    if (ei.VB() != VOL)
      throw Exception ("SetOrder: element orders are for volume elements, use a boundary region or NodeId");
    size_t ne = ma->GetNE();
    if (ei.Nr() >= ne)
      throw Exception ("SetOrder: element number out of range");

    bool isotropic = (o[0] == o[1] && o[1] == o[2]);
    switch (ma->GetElType(ei))
      {
      case ET_SEGM:
        o = IVec<3> (o[0]); break;
      case ET_QUAD:
        o[2] = o[1]; break;
      case ET_PRISM:
        if (o[0] != o[1])
          throw Exception ("SetOrder: a prism has one order in its triangle and one along its axis");
        break;
      case ET_HEX:
        break;
      default:
        if (!isotropic)
          throw Exception ("SetOrder: different orders in different directions only for quadrilaterals, hexahedra and prisms");
      }

    auto & eo = order_spec.base->element_order;
    if (eo.Size() != ne) InheritElementOrders();
    if (eo.Size() != ne)
      {
        eo.SetSize (ne);
        eo = IVec<3,short> (-1);
      }
    eo[ei.Nr()] = IVec<3,short> (short(o[0]), short(o[1]), short(o[2]));

    order_spec.base->timestamp = GetNextTimeStamp();
    timestamp = 0;
    order_stale = 2;
  }

  void FESpace :: SetOrder (ElementId ei, IVec<3> order)
  {
    throw Exception (string("FESpace::SetOrder for elements not available for space ") + GetClassName());
  }

  IVec<3> FESpace :: GetOrder (ElementId ei) const
  {
    throw Exception (string("FESpace::GetOrder for elements not available for space ") + GetClassName());
  }

  void FESpace :: SetOrderSpecification (const OrderSpecification & spec, bool creator)
  {
    order_spec.base = spec.base;
    if (creator)
      order_spec.base->creator = this;
    for (int i = 0; i < 3; i++)
      order_spec.node_order[i] = spec.node_order[i];
    order_spec.node_order_level = spec.node_order_level;

    timestamp = 0;
    order_stale = 2;
  }

  bool FESpace :: LinkOrders (const Flags & flags)
  {
    if (!flags.AnyFlagDefined ("orderlink")) return false;
    order_spec.base = std::any_cast<shared_ptr<BaseOrders>> (flags.GetAnyFlag ("orderlink"));
    order_spec.link_offset = int (flags.GetNumFlag ("relorder", 0));
    order = max2 (0, order_spec.base->order + order_spec.link_offset);
    order_linked = true;
    return true;
  }

  void FESpace :: CheckOrderCreator () const
  {
    auto & base = *order_spec.base;
    if (base.creator != this)
      throw Exception ("SetOrder: this space is linked, the orders belong to "
                       + (base.creator_name.empty() ? string("another space") : base.creator_name));
  }

  void FESpace :: SetOrder (const Region & region, int order)
  {
    throw Exception (string("FESpace::SetOrder for regions not available for space ") + GetClassName());
  }

  void FESpace :: SetOrder (const Region & region, IVec<2> order, shared_ptr<CoefficientFunction> direction)
  {
    throw Exception (string("FESpace::SetOrder for regions not available for space ") + GetClassName());
  }

  ORDER_POLICY FESpace :: DerivedOrderPolicy () const
  {
    if (order_spec.base->mesh_orders || HasNodeOrders() || HasElementOrders() || order_spec.base->region_orders.Size()
        || order_spec.base->hp_slope > 0)
      return VARIABLE_ORDER;
    for (int bonus : et_bonus_order)
      if (bonus) return NODE_TYPE_ORDER;
    if (order_spec.order_edge > -1 || order_spec.order_face > -1 || order_spec.order_cell > -1)
      return NODE_TYPE_ORDER;
    return CONSTANT_ORDER;
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

  void FESpace :: SelectDofs (const string & name, BitArray & dofs) const
  {
    dofs.Clear();
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
      case NT_EDGE:
        GetEdgeDofNrs(ni.GetNr(), dnums);
        if (ma->GetDimension() == 1)
          {
            // the edge is an element: add its inner dofs
            Array<int> elnums, inner;
            ma->GetEdgeElements(ni.GetNr(), elnums);
            for (auto el : elnums)
              {
                GetInnerDofNrs(el, inner);
                dnums += inner;
              }
          }
        break;
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

  bool FESpace :: OwnsOperator (const DifferentialOperator * op) const
  {
    if (!op) return true;
    for (VorB vb : { VOL, BND, BBND, BBBND })
      if (GetEvaluator(vb).get() == op || GetFluxEvaluator(vb).get() == op)
        return true;
    auto additional = GetAdditionalEvaluators();
    for (size_t i = 0; i < additional.Size(); i++)
      if (additional[i].get() == op)
        return true;
    return false;
  }

  // spaces handing out the same elements: wrappers and low-order spaces
  static Array<const FESpace*> ElementSpaces (const FESpace & fes)
  {
    Array<const FESpace*> list { &fes };
    for (size_t i = 0; i < list.Size(); i++)
      {
        auto add = [&] (const FESpace * s) { if (s && !list.Contains(s)) list.Append(s); };
        add (list[i]->ElementSpace());
        add (list[i]->LowOrderFESpacePtr().get());
      }
    return list;
  }

  // one element per element type of vb
  Array<ElementId> RepresentativeElements (const MeshAccess & ma, VorB vb)
  {
    Array<ElementId> reps;
    Array<ELEMENT_TYPE> types;
    if (ma.GetNE(vb) == 0) return reps;
    auto & classes = ma.GetElementsOfClass(vb);
    for (size_t i = 0; i < classes.Size(); i++)
      if (classes[i].Size())
        {
          ElementId ei(vb, classes[i][0]);
          auto et = ma.GetElType(ei);
          if (!types.Contains(et))
            {
              types.Append(et);
              reps.Append(ei);
            }
        }
    return reps;
  }

  void CheckOperatorSpace (const DifferentialOperator & op, const FESpace * opspace,
                           const FESpace & space, VorB elvb, bool strict)
  {
    // level 1
    if (opspace)
      {
        auto a = ElementSpaces(*opspace);
        for (auto s : ElementSpaces(space))
          if (a.Contains(s)) return;
      }
    else if (space.OwnsOperator(&op))
      return;

    auto describe = [&] ()
    {
      return "operator '" + op.Name() + "'" + (opspace ? string(" of ") + opspace->GetClassName() : string("")) +
        " applied to elements of " + space.GetClassName();
    };
    
    // level 3
    if (!op.DeclaresElementType())
      {
        if (!strict) return;
        throw Exception (describe() + ": the operator does not declare its element type, use trial/test-functions of the form's space");
      }

    // level 2: the operator acting on elements of elvb
    shared_ptr<DifferentialOperator> hold;
    const DifferentialOperator * eop = &op;
    for (int i = 0; i < 4 && eop && eop->VB() != elvb; i++)
      {
        hold = eop->GetTrace();
        eop = hold.get();
      }
    if (!eop || eop->VB() != elvb)
      eop = &op;

    // the operator must fit the geometry: space dimension and dimension of the evaluated elements
    auto & ma = *space.GetMeshAccess();
    int dimspace = eop->DimSpace(), dimel = eop->DimElement();
    if (dimspace < 0 || dimel < 0)
      throw Exception (describe() + ": the operator does not declare its dimensions, use trial/test-functions of the form's space");
    if (dimspace > 0)   // 0: geometry free (NumberSpace)
      {
        int eldim = ma.GetDimension() - int(elvb);
        if (dimspace != ma.GetDimension())
          throw Exception (describe() + ": the operator acts in " + ToString(dimspace) + "D space, the mesh of " +
                           space.GetClassName() + " has dimension " + ToString(ma.GetDimension()));
        if (dimel != eldim && !eop->SupportsVB(elvb))
          throw Exception (describe() + ": the operator acts on elements of dimension " + ToString(dimel) +
                           ", evaluated on elements of dimension " + ToString(eldim));
      }

    LocalHeap lh(1000*1000, "CheckOperatorSpace");
    for (ElementId ei : RepresentativeElements (ma, elvb))
      {
        HeapReset hr(lh);
        const FiniteElement & fel = space.GetFE (ei, lh);
        try
          {
            eop->CheckElement (fel);
          }
        catch (const Exception & e)
          {
            throw Exception (describe() + " on " + ToString(ma.GetElType(ei)) + " (" + ToString(elvb) + ") elements: " + e.What());
          }
      }
  }

  CFContext MakeCFContext (const MeshAccess & ma, CFContext::Purpose purpose,
                           VorB vb, VorB element_vb, bool skeleton, bool allow_complex)
  {
    // every element transformation evaluates the deformation
    if (auto def = ma.GetDeformation())
      if (def->GetLevelUpdated() < ma.GetNLevels())
        throw Exception ("the deformation of the mesh is not up to date with its mesh (refined without Update)");
    CFContext ctx;
    ctx.purpose = purpose;
    ctx.mesh = &ma;
    ctx.dim = ma.GetDimension();
    ctx.vb = vb;
    ctx.element_vb = element_vb;
    ctx.skeleton = skeleton;
    ctx.allow_complex = allow_complex;
    for (VorB v : { VOL, BND, BBND, BBBND })
      ctx.nregions[v] = ma.GetNRegions(v);
    ctx.check_proxy = [] (const ProxyFunction & proxy, const FESpace & space, VorB elvb)
      {
        CheckOperatorSpace (*proxy.Evaluator(), proxy.GetFESpace().get(), space, elvb, true);
      };
    return ctx;
  }

  void CheckPointContext (const CoefficientFunction & cf, const MeshAccess * ma,
                          VorB vb, int dim, VorB element_vb)
  {
    int key = int(vb) | (dim << 2) | (int(element_vb) << 6);
    // a refined mesh is a new context
    size_t stamp = ma ? ma->GetTimeStamp() : 0;
    if (cf.gate_key == key && cf.gate_mesh == ma && cf.gate_stamp == stamp) return;
    CFContext ctx;
    if (ma)
      ctx = MakeCFContext (*ma, CFContext::EVALUATE, vb, element_vb);
    ctx.purpose = CFContext::EVALUATE;
    ctx.vb = vb;
    ctx.element_vb = element_vb;
    ctx.dim = dim;
    ngfem::CheckContext (cf, ctx);
    cf.gate_mesh = ma;
    cf.gate_key = key;
    cf.gate_stamp = stamp;
  }

  void CheckPointContext (const CoefficientFunction & cf, const BaseMappedIntegrationPoint & mip)
  {
    auto & trafo = mip.GetTransformation();
    VorB eb = VOL;
    if (mip.IP().FacetNr() != -1)
      eb = mip.IP().VB() == VOL ? BND : mip.IP().VB();
    CheckPointContext (cf, static_cast<const MeshAccess*> (trafo.GetMesh()), trafo.VB(),
                       mip.DimSpace(), eb);
  }

  // sub-entities the element loops integrate over, boundary facets for ds(skeleton=True)
  static void CheckIntegrationDomain (const MeshAccess & ma, VorB vb, VorB element_vb, bool skeleton)
  {
    int eldim = ma.GetDimension() - int(skeleton ? VOL : vb);
    if (element_vb != VOL && int(element_vb) > eldim)
      throw Exception ("element_vb = " + ToString(element_vb) + " on elements of dimension " + ToString(eldim) +
                       ": the elements have no such sub-entities");
    if (skeleton && vb == BND && ma.GetNE(BND))
      {
        // the loop over boundary elements continues on the volume element of the facet
        if (ma.GetNE(VOL) == 0)
          throw Exception ("ds(skeleton=True) integrates over boundary facets of volume elements, the mesh has no volume elements");
        Array<int> elnums;
        for (size_t i = 0; i < ma.GetNE(BND); i++)
          {
            ma.GetFacetElements (ma.GetElFacets(ElementId(BND, i))[0], elnums);
            if (elnums.Size() == 0)
              throw Exception ("ds(skeleton=True): boundary element " + ToString(i) + " is not the facet of a volume element");
          }
      }
  }
  
  void CheckIntegratorContext (const Integrator & igt, CFContext::Purpose purpose,
                               shared_ptr<FESpace> trialspace, shared_ptr<FESpace> testspace,
                               bool allow_complex)
  {
    if (auto cbfi = dynamic_cast<const CompoundBilinearFormIntegrator*> (&igt))
      {
        auto ctrial = dynamic_pointer_cast<CompoundFESpace> (trialspace);
        auto ctest = dynamic_pointer_cast<CompoundFESpace> (testspace);
        int comp = cbfi->GetComponent();
        if (!ctrial || !ctest || comp < 0 || comp >= ctrial->GetNSpaces() || comp >= ctest->GetNSpaces())
          throw Exception ("component integrator " + ToString(comp) + " needs a compound space with that component");
        CheckIntegratorContext (*cbfi->GetBFI(), purpose, (*ctrial)[comp], (*ctest)[comp], allow_complex);
        return;
      }
    if (auto clfi = dynamic_cast<const CompoundLinearFormIntegrator*> (&igt))
      {
        auto ctest = dynamic_pointer_cast<CompoundFESpace> (testspace);
        int comp = clfi->GetComponent();
        if (!ctest || comp < 0 || comp >= ctest->GetNSpaces())
          throw Exception ("component integrator " + ToString(comp) + " needs a compound space with that component");
        CheckIntegratorContext (*clfi->GetLFI(), purpose, nullptr, (*ctest)[comp], allow_complex);
        return;
      }

    // tensor-product spaces have their own geometry
    if (dynamic_pointer_cast<TPHighOrderFESpace> (trialspace) ||
        dynamic_pointer_cast<TPHighOrderFESpace> (testspace))
      return;
    
    auto space = trialspace ? trialspace : testspace;
    if (!space) return;
    auto & ma = *space->GetMeshAccess();

    // loops test definedonelements by element number, facet loops by facet number
    // (boundary skeleton terms: boundary element in assembly, facet in apply)
    if (auto doe = igt.GetDefinedOnElements())
      {
        size_t needed = ma.GetNE(igt.VB());
        string what = ToString(igt.VB()) + " elements";
        if (igt.SkeletonForm())
          {
            if (igt.GetDGFormulation().element_boundary)
              { needed = ma.GetNE(VOL); what = "VOL elements"; }
            else if (igt.VB() == VOL)
              { needed = ma.GetNFacets(); what = "facets"; }
            else
              { needed = max(ma.GetNE(BND), ma.GetNFacets()); what = "BND elements/facets"; }
          }
        if (doe->Size() < needed)
          throw Exception ("definedonelements has size " + ToString(doe->Size()) + ", mesh has " +
                           ToString(needed) + " " + what);
      }
    
    if (igt.GetDeformation())
      CheckDeformation (*igt.GetDeformation(), ma);

    // coefficients of classic integrators are evaluated on the elements of igt.VB()
    if (igt.context_coefficients.Size())
      {
        CFContext cctx = MakeCFContext (ma, CFContext::EVALUATE, igt.VB(), VOL, false, allow_complex);
        for (auto & c : igt.context_coefficients)
          if (c) CheckContext (*c, cctx);
      }
    
    shared_ptr<CoefficientFunction> cf;
    VorB element_vb = VOL;
    bool skeleton = false;
    if (auto s = dynamic_cast<const SymbolicBilinearFormIntegrator*> (&igt))
      {
        cf = s->GetIntegrandCF();
        element_vb = s->ElementVB();
      }
    else if (auto s = dynamic_cast<const SymbolicFacetBilinearFormIntegrator*> (&igt))
      {
        cf = s->GetIntegrandCF();
        bool eb = s->GetDGFormulation().element_boundary;
        skeleton = !eb;
        element_vb = eb ? BND : VOL;
      }
    else if (auto s = dynamic_cast<const SymbolicLinearFormIntegrator*> (&igt))
      {
        cf = s->GetIntegrandCF();
        element_vb = s->ElementVB();
      }
    else if (auto s = dynamic_cast<const SymbolicFacetLinearFormIntegrator*> (&igt))
      {
        cf = s->GetIntegrandCF();
        skeleton = true;
      }
    else if (auto s = dynamic_cast<const SymbolicEnergy*> (&igt))
      {
        cf = s->GetIntegrandCF();
        element_vb = s->ElementVB();
        purpose = CFContext::ENERGY;
      }
    if (!cf)
      {
        CheckIntegrationDomain (ma, igt.VB(), VOL, igt.SkeletonForm());
        // classic integrators: element family, one element per element type
        if (trialspace && testspace && trialspace != testspace) return;
        VorB vb = igt.SkeletonForm() ? VOL : igt.VB();
        if (int(vb) > ma.GetDimension()) return;
        LocalHeap lh(1000*1000, "CheckIntegratorContext");
        Array<ELEMENT_TYPE> types;
        auto & classes = ma.GetElementsOfClass(vb);
        for (size_t i = 0; i < classes.Size(); i++)
          for (auto nr : classes[i])
            {
              ElementId ei(vb, nr);
              if (!space->DefinedOn(ei) || !igt.DefinedOn(ma.GetElIndex(ei))) continue;
              auto et = ma.GetElType(ei);
              if (!types.Contains(et))
                {
                  types.Append (et);
                  HeapReset hr(lh);
                  try
                    {
                      igt.CheckElement (space->GetFE(ei, lh));
                    }
                  catch (const Exception & e)
                    {
                      throw Exception ("integrator '" + igt.Name() + "' used with " + space->GetClassName() +
                                       " on " + ToString(et) + " (" + ToString(vb) + ") elements: " + e.What());
                    }
                }
              break;
            }

        // classic integrators compute element matrices of size dimension*ndof
        const Integrator * plain = &igt;
        shared_ptr<Integrator> hold;
        if (auto c = dynamic_cast<const ComplexBilinearFormIntegrator*> (plain))
          plain = (hold = c->GetBFI()).get();
        if (auto a = dynamic_cast<const BilinearFormIntegratorAnyDim*> (plain))
          plain = (hold = a->GetBFI(ma.GetDimension())).get();
        else if (auto a = dynamic_cast<const LinearFormIntegratorAnyDim*> (plain))
          plain = (hold = a->GetLFI(ma.GetDimension())).get();
        if (!dynamic_cast<const ComplexLinearFormIntegrator*> (plain) && !dynamic_pointer_cast<CompoundFESpace>(space)
            && plain->GetDimension() != space->GetDimension())
          throw Exception ("integrator '" + igt.Name() + "' needs a space of dimension " + ToString(plain->GetDimension()) +
                           ", " + space->GetClassName() + " has dimension " + ToString(space->GetDimension()));
        return;
      }

    CheckIntegrationDomain (ma, igt.VB(), element_vb, skeleton);
    CFContext ctx = MakeCFContext (ma, purpose, igt.VB(), element_vb, skeleton, allow_complex);
    ctx.trialspace = trialspace;
    ctx.testspace = testspace;
    CheckContext (*cf, ctx);
  }

  
  void FESpace :: CheckUpToDate () const
  {
    if (level_updated >= 0 && level_updated < ma->GetNLevels())
      throw Exception (GetClassName() + " is not up to date with its mesh (refined without Update)");
    if (level_updated >= 0 && (order_stale || order_spec.base->timestamp > order_base_resolved))
      throw Exception (GetClassName() + " is not up to date, its order was changed without Update");
  }
  
  void FESpace :: CheckDofConsistency (VorB vb, bool all_elements) const
  {
    CheckUpToDate();
    size_t stamp = ma->GetTimeStamp();
    if (!all_elements && dofcheck_meshstamp[vb] == stamp && dofcheck_ndof[vb] == GetNDof()) return;
    if (!NeedsDofTransformation() && int(vb) <= ma->GetDimension() && ma->GetNE(vb))
      {
        LocalHeap lh(1000*1000, "CheckDofConsistency");
        Array<DofId> dnums;
        auto & classes = ma->GetElementsOfClass(vb);
        // one element per element class and region
        Array<int> seen;
        for (size_t i = 0; i < classes.Size(); i++)
          {
            seen.SetSize0();
            for (auto nr : classes[i])
              {
                // assembly loops skip elements the space is not defined on
                ElementId ei(vb, nr);
                if (!DefinedOn(ei)) continue;
                int index = ma->GetElIndex(ei);
                if (!all_elements)
                  {
                    if (seen.Contains(index)) continue;
                    seen.Append(index);
                  }
                HeapReset hr(lh);
                const FiniteElement & fel = GetFE (ei, lh);
                GetDofNrs (ei, dnums);
                if (fel.GetNDof() != dnums.Size())
                  throw Exception (string(GetClassName()) + ": inconsistent number of degrees of freedom on " +
                                   ToString(ma->GetElType(ei)) + " (" + ToString(vb) + ") elements, the element has " +
                                   ToString(fel.GetNDof()) + ", the space provides " + ToString(dnums.Size()));
              }
          }
      }
    dofcheck_meshstamp[vb] = stamp;
    dofcheck_ndof[vb] = GetNDof();
  }
  
  void FESpace :: CheckElementFamilies (LocalHeap & lh) const
  {
    auto additional = GetAdditionalEvaluators();
    for (VorB vb : { VOL, BND, BBND, BBBND })
      {
        if (int(vb) > ma->GetDimension()) continue;

        Array<pair<string, shared_ptr<DifferentialOperator>>> ops;
        if (evaluator[vb]) ops.Append ({ "evaluator", evaluator[vb] });
        if (flux_evaluator[vb]) ops.Append ({ "flux_evaluator", flux_evaluator[vb] });
        for (size_t i = 0; i < additional.Size(); i++)
          if (additional[i] && additional[i]->SupportsVB(vb))
            {
              // an operator of higher dimension acts on lower-dimensional elements through its trace
              auto op = additional[i];
              while (op && op->VB() != vb)
                op = op->GetTrace();
              ops.Append ({ "additional evaluator '" + additional.GetName(i) + "'", op ? op : additional[i] });
            }
        if (!ops.Size()) continue;

        for (ElementId ei : ma->Elements(vb))
          {
            HeapReset hr(lh);
            const FiniteElement & fel = GetFE (ei, lh);
            for (auto & [name, op] : ops)
              try
                {
                  op->CheckElement (fel);
                }
              catch (Exception & e)
                {
                  throw Exception (string(GetClassName()) + ": " + name + " " + op->Name() + ", element " + ToString(ei.Nr()) +
                                   ", VorB " + ToString(vb) + ", element type " + ToString(ma->GetElType(ei)) +
                                   ": " + e.What());
                }
          }
      }
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
                                            false, IsComplex(), evaluator);
    auto test  = make_shared<ProxyFunction>(nullptr, // dynamic_pointer_cast<FESpace>(const_cast<FESpace*>(this)->shared_from_this()),
                                            true, false, evaluator);
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
    ost << "definedon = " << definedon_region.Mask(VOL) << endl;
    ost << "definedon boundary = " << definedon_region.Mask(BND) << endl;
    ost << "definedon codim 2 = " << definedon_region.Mask(BBND) << endl;
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
    Array<bool> defon_ar[3];
    if (archive.Output())
      for (int i = 0; i < 3; i++)
        {
          auto & mask = definedon_region.Mask(VorB(i));
          defon_ar[i].SetSize(mask.Size());
          for (auto j : Range(mask.Size()))
            defon_ar[i][j] = mask.Test(j);
        }
    archive & defon_ar[0] & defon_ar[1] & defon_ar[2];
    if (archive.Input())
      for (int i = 0; i < 3; i++)
        {
          auto & mask = definedon_region.Mask(VorB(i));
          if (!defon_ar[i].Size())
            mask.Set();
          else
            {
              mask.SetSize(defon_ar[i].Size());
              mask.Clear();
              for (auto j : Range(defon_ar[i]))
                if (defon_ar[i][j]) mask.SetBit(j);
            }
        }
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



  shared_ptr<BitArray> FESpace :: GetSmoothingFreeDofs (const Flags & flags) const
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

    if (flags.AnyFlagDefined("additional_dirichlet_constraints"))
      {
        Region reg = OwnRegion(std::any_cast<Region>(flags.GetAnyFlag("additional_dirichlet_constraints")), *ma, "additional_dirichlet_constraints");
        BitArray dofs = GetDofs(reg);
        dofs.Invert();
        freedofs = make_shared<BitArray>(*freedofs);
        freedofs->And(dofs);
      }
    if (flags.AnyFlagDefined("additional_dirichlet_boundaries"))
      {
        auto dirbcs = std::any_cast<Array<DirichletBoundary>>(flags.GetAnyFlag("additional_dirichlet_boundaries"));
        freedofs = make_shared<BitArray>(*freedofs);
        for (auto dbc : dirbcs)
          {
            MixedRegion reg = ResolveRegion(dbc, ma);
            BitArray dofs = GetDofs(reg, dbc.proxy->Evaluator().get());
            dofs.Invert();
            freedofs->And(dofs);
          }
      }
    return freedofs;
  }


  MixedRegion ResolveRegion (const DirichletBoundary & dbc, const shared_ptr<MeshAccess> & ma)
  {
    if (auto reg = std::any_cast<MixedRegion>(&dbc.region))
      {
        if (reg->Mesh() != ma)
          throw Exception("Dirichlet condition: the region belongs to another mesh");
        return *reg;
      }
    return MixedRegion(ma, dbc.vbn);
  }


  shared_ptr<Table<int>> FESpace :: CreateSmoothingBlocks (const Flags & aflags) const
  {
    static Timer t("FESpace::CreateSmoothingBlocks");
    RegionTimer reg(t);

    Flags flags(aflags);
    if (flags.GetDefineFlag("condense"))
      flags.SetFlag("eliminate_internal");

    auto freedofs = GetSmoothingFreeDofs(flags);

    Array<string> blocktypes { flags.GetStringListFlag("blocktype") };
    if (flags.StringFlagDefined("blocktype"))
      blocktypes += flags.GetStringFlag("blocktype");

    shared_ptr<Table<int>> blocks;
    if (blocktypes.Size())
      blocks = CreateGenericSmoothingBlocks (blocktypes, freedofs);
    else
      {
        blocks = CreateSpecificSmoothingBlocks (flags);
        if (!blocks)
          {
            blocktypes += string("vertexpatch");
            blocks = CreateGenericSmoothingBlocks (blocktypes, freedofs);
          }
        else
          {
            // space-specific tables may contain constrained dofs
            Array<int> cnt(blocks->Size());
            size_t removed = 0;
            for (auto i : Range(blocks->Size()))
              {
                int c = 0;
                for (auto d : (*blocks)[i])
                  {
                    if (d < 0 || d >= GetNDof())
                      throw Exception("CreateSmoothingBlocks: dof "+ToString(d)+" in block "+ToString(i)+" out of range");
                    if (freedofs->Test(d)) c++;
                  }
                cnt[i] = c;
                removed += (*blocks)[i].Size()-c;
              }
            if (removed)
              {
                auto filtered = make_shared<Table<int>> (cnt);
                for (auto i : Range(blocks->Size()))
                  {
                    int j = 0;
                    for (auto d : (*blocks)[i])
                      if (freedofs->Test(d))
                        (*filtered)[i][j++] = d;
                  }
                blocks = filtered;
              }
          }
      }

    if (print)
      *testout << "smoothing blocks = " << endl << *blocks << endl;
    return blocks;
  }


  shared_ptr<Table<int>> FESpace :: CreateSpecificSmoothingBlocks (const Flags & flags) const
  {
    return nullptr;
  }


  shared_ptr<Table<int>> FESpace ::
  CreateGenericSmoothingBlocks (FlatArray<string> blocktypes, shared_ptr<BitArray> freedofs) const
  {
    int dim = ma->GetDimension();
    Array<SmoothingBlockType> types;
    
    for (auto fulltype : blocktypes)
      {
        auto pos = fulltype.find(':');
        string name = fulltype.substr(0, pos);

        SmoothingBlockType type;
        if (name == "vertex")           type = { NT_VERTEX, NT_VERTEX };
        else if (name == "edge")        type = { NT_EDGE, NT_EDGE };
        else if (name == "face")        type = { NT_FACE, NT_FACE };
        else if (name == "facet")       type = { dim-1, dim-1 };
        else if (name == "cell")        type = { NT_CELL, NT_CELL };
        else if (name == "inner")       type = { NT_ELEMENT, NT_ELEMENT };
        else if (name == "element")     type = { -1, -1 };
        else if (name == "vertexedge")  type = { NT_VERTEX, NT_EDGE };
        else if (name == "vertexpatch") type = { NT_VERTEX, NT_CELL };
        else if (name == "edgepatch")   type = { NT_EDGE, NT_CELL };
        else if (name == "facepatch")   type = { NT_FACE, NT_CELL };
        else
          throw Exception("FESpace::CreateSmoothingBlocks type '"+name+"' not available, available types are:\n"
                          "vertex ... dofs of a vertex\n"
                          "edge ..... dofs of an edge\n"
                          "face ..... dofs of a face\n"
                          "cell ..... dofs of a cell (3D)\n"
                          "facet .... dofs of a facet (codimension 1)\n"
                          "inner .... inner dofs of an element (codimension 0)\n"
                          "element .. all dofs of an element (including its boundary)\n"
                          "vertexpatch .. dofs of the vertex and of all edges, faces, cells containing it\n"
                          "edgepatch .... dofs of the edge and of all faces, cells containing it\n"
                          "facepatch .... dofs of the face and of the cells containing it\n"
                          "vertexedge ... dofs of the vertex and of the edges containing it\n"
                          "options: 'type:periodic', 'type:dofset', 'type:dofset,periodic'\n"
                          );

        if (pos != string::npos)
          {
            string opts = fulltype.substr(pos+1);
            while (opts.size())
              {
                auto end = opts.find_first_of(",: ");
                string opt = opts.substr(0, end);
                opts = (end == string::npos) ? string() : opts.substr(end+1);
                if (opt.empty()) continue;
                if (opt == "periodic")
                  type.periodic = true;
                else if (!type.dofset)
                  type.dofset = opt;
                else
                  throw Exception("FESpace::CreateSmoothingBlocks: more than one dof set in '"+fulltype+"'");
              }
          }

        // blocks in the order of the list, each type once
        // inner is the node type of the elements
        auto nodetype = [dim] (int nt) { return nt == NT_ELEMENT ? dim : nt; };
        bool found = false;
        for (auto & other : types)
          if (nodetype(other.anchor) == nodetype(type.anchor) &&
              nodetype(other.maxnt) == nodetype(type.maxnt) &&
              other.periodic == type.periodic && other.dofset == type.dofset)
            found = true;
        if (!found)
          types += type;
      }

    return CreateGenericSmoothingBlocks (types, freedofs);
  }


  shared_ptr<Table<int>> FESpace ::
  CreateGenericSmoothingBlocks (FlatArray<SmoothingBlockType> types, shared_ptr<BitArray> freedofs) const
  {
    static Timer t("FESpace::CreateGenericSmoothingBlocks");
    RegionTimer reg(t);

    int dim = ma->GetDimension();
    size_t ntypes = types.Size();
    Array<size_t> base(ntypes);
    Array<int> anchors(ntypes), maxnts(ntypes);
    Array<Array<int>> entry(ntypes);    // anchor node -> block, if not the identity
    Array<BitArray> mask(ntypes);

    size_t nblocks = 0;
    bool neednodes[4] = { false, false, false, false };
    bool needelements = false;
    
    for (auto ti : Range(ntypes))
      {
        auto & type = types[ti];
        anchors[ti] = type.anchor;
        maxnts[ti] = type.maxnt;
        size_t nb;
        if (type.anchor < 0)
          {
            if (type.periodic)
              throw Exception("FESpace::CreateSmoothingBlocks: periodic not available for element blocks");
            nb = ma->GetNE(VOL);
            needelements = true;
          }
        else if (type.anchor == NT_ELEMENT)
          {
            // inner dofs, blocks numbered like the elements
            if (type.periodic)
              throw Exception("FESpace::CreateSmoothingBlocks: periodic not available for element blocks");
            anchors[ti] = maxnts[ti] = dim;
            neednodes[dim] = true;
            nb = ma->GetNNodes(NODE_TYPE(dim));
            if (dim == 2)
              {
                entry[ti].SetSize(nb);
                for (auto i : Range(nb))
                  entry[ti][i] = ma->GetNode<2>(i).surface_el;
                nb = ma->GetNE(VOL);
              }
          }
        else if (type.anchor > dim)
          nb = 0;   // no such nodes on this mesh
        else
          {
            maxnts[ti] = min(type.maxnt, dim);
            auto nt = NODE_TYPE(type.anchor);
            nb = ma->GetNNodes(nt);
            for (int i = type.anchor; i <= maxnts[ti]; i++)
              neednodes[i] = true;

            if (type.periodic)
              {
                // blocks of identified anchors are merged
                Array<int> rep(nb);
                for (auto i : Range(rep))
                  rep[i] = i;
                if (nt != NT_CELL)
                  for (auto idnr : Range(ma->GetNPeriodicIdentifications()))
                    for (auto & pair : ma->GetPeriodicNodes(nt, idnr))
                      rep[pair[1]] = pair[0];
                for (auto i : Range(rep))
                  while (rep[rep[i]] != rep[i])
                    rep[i] = rep[rep[i]];

                entry[ti].SetSize(rep.Size());
                nb = 0;
                for (auto i : Range(rep))
                  if (rep[i] == i)
                    entry[ti][i] = nb++;
                for (auto i : Range(rep))
                  if (rep[i] != i)
                    entry[ti][i] = entry[ti][rep[i]];
              }
          }

        if (type.dofset)
          {
            mask[ti].SetSize (GetNDof());
            SelectDofs (*type.dofset, mask[ti]);
          }
        
        base[ti] = nblocks;
        nblocks += nb;
      }


    // free dofs per node / per element
    auto CollectDofs = [&] (size_t n, auto getdofs)
    {
      Array<int> cnt(n);
      ParallelForRange (n, [&] (IntRange r)
      {
        Array<DofId> dofs;
        for (auto i : r)
          {
            getdofs (i, dofs);
            int c = 0;
            for (auto d : dofs)
              if (IsRegularDof(d) && freedofs->Test(d)) c++;
            cnt[i] = c;
          }
      });
      Table<DofId> table(cnt);
      ParallelForRange (n, [&] (IntRange r)
      {
        Array<DofId> dofs;
        for (auto i : r)
          {
            getdofs (i, dofs);
            int c = 0;
            for (auto d : dofs)
              if (IsRegularDof(d) && freedofs->Test(d))
                table[i][c++] = d;
          }
      });
      return table;
    };
    
    Table<DofId> nodedofs[4];
    for (int nt = 0; nt < 4; nt++)
      if (neednodes[nt])
        nodedofs[nt] = CollectDofs (ma->GetNNodes(NODE_TYPE(nt)), [&] (size_t i, Array<DofId> & dofs)
                                    { GetDofNrs (NodeId(NODE_TYPE(nt), i), dofs); });
    Table<DofId> eldofs;
    if (needelements)
      eldofs = CollectDofs (ma->GetNE(VOL), [&] (size_t i, Array<DofId> & dofs)
                            { GetDofNrs (ElementId(VOL,i), dofs); });


    auto ForAnchors = [&] (int nt, size_t nr, int anchor, auto func)
    {
      if (nt == anchor)
        {
          func(nr);
          return;
        }
      switch (nt)
        {
        case NT_EDGE:
          for (auto v : ma->GetEdgePNums(nr)) func(v);
          break;
        case NT_FACE:
          if (anchor == NT_VERTEX)
            for (auto v : ma->GetFacePNums(nr)) func(v);
          else
            for (auto e : ma->GetFaceEdges(nr)) func(e);
          break;
        case NT_CELL:
          {
            auto el = ma->GetElement(ElementId(VOL,nr));
            if (anchor == NT_VERTEX)
              for (auto v : el.Vertices()) func(v);
            else if (anchor == NT_EDGE)
              for (auto e : el.Edges()) func(e);
            else
              for (auto f : el.Faces()) func(f);
            break;
          }
        default:
          break;
        }
    };

    // func (block, dofs, mask)
    auto ForContributions = [&] (auto func)
    {
      for (auto ti : Range(ntypes))
        {
          const BitArray * pmask = types[ti].dofset ? &mask[ti] : nullptr;
          if (anchors[ti] < 0)
            {
              ParallelFor (eldofs.Size(), [&] (size_t i)
              { func (base[ti]+i, eldofs[i], pmask); });
              continue;
            }
          FlatArray<int> map = entry[ti];
          for (int nt = anchors[ti]; nt <= maxnts[ti]; nt++)
            ParallelFor (nodedofs[nt].Size(), [&] (size_t i)
            {
              auto dofs = nodedofs[nt][i];
              if (dofs.Size() == 0) return;
              ForAnchors (nt, i, anchors[ti], [&] (size_t a)
              {
                if (map.Size() == 0)
                  func (base[ti]+a, dofs, pmask);
                else if (map[a] >= 0)
                  func (base[ti]+map[a], dofs, pmask);
              });
            });
        }
    };

    Array<int> cnt(nblocks);
    cnt = 0;
    ForContributions ([&] (size_t block, FlatArray<DofId> dofs, const BitArray * pmask)
    {
      int c = dofs.Size();
      if (pmask)
        {
          c = 0;
          for (auto d : dofs)
            if (pmask->Test(d)) c++;
        }
      AsAtomic(cnt[block]) += c;
    });

    Table<int> table(cnt);
    cnt = 0;
    ForContributions ([&] (size_t block, FlatArray<DofId> dofs, const BitArray * pmask)
    {
      int c = dofs.Size();
      if (pmask)
        {
          c = 0;
          for (auto d : dofs)
            if (pmask->Test(d)) c++;
        }
      int pos = AsAtomic(cnt[block]).fetch_add(c);
      auto row = table[block];
      for (auto d : dofs)
        if (!pmask || pmask->Test(d))
          row[pos++] = d;
    });

    // sorted rows without duplicates, independent of the number of threads
    atomic<size_t> nduplicates{0};
    ParallelFor (table.Size(), [&] (size_t i)
    {
      auto row = table[i];
      QuickSort (row);
      int n = 0;
      for (auto j : Range(row))
        if (j == 0 || row[j] != row[j-1])
          row[n++] = row[j];
      cnt[i] = n;
      if (n != row.Size())
        nduplicates += row.Size()-n;
    });

    if (nduplicates)
      {
        Table<int> compressed(cnt);
        ParallelFor (table.Size(), [&] (size_t i)
        { compressed[i] = table[i].Range(0, cnt[i]); });
        table = std::move(compressed);
      }
    
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
    definedon_region.Mask(vb) = defon;

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
    const Table<int> & element_coloring = fes.ElementColoring(vb);

    const Table<int> & blocks = fes.ElementBlocks(vb);
    const Table<int> & block_coloring = fes.BlockColoring(vb);
    // estimated parallel time (in elements) of both colorings
    size_t nt = TaskManager::GetNumThreads();
    size_t tblocks = 0, tels = 0;
    for (auto c : block_coloring)
      if (c.Size()) tblocks += (c.Size()+nt-1)/nt * blocks[c[0]].Size();
    for (auto c : element_coloring) tels += (c.Size()+nt-1)/nt;
    // blocks gain locality, but lose parallelism on small meshes
    bool use_blocks = blocks.Size() && (!GetTaskManager() || tblocks <= 1.5 * tels);
    if (use_blocks)
      {
        auto do_block = [&] (int b, LocalHeap & lh)
        {
          for (int nr : blocks[b])
            {
              HeapReset hr(lh);
              func (FESpace::Element(fes, ElementId (vb, nr), lh), lh);
            }
        };
        if (auto *task_manager = GetTaskManager(); task_manager)
          for (FlatArray<int> blocks_of_col : block_coloring)
            {
              SharedLoop2 sl(blocks_of_col.Range());
              task_manager -> CreateJob
                ( [&] (const TaskInfo & ti)
                  {
                    LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
                    for (int mynr : sl)
                      do_block (blocks_of_col[mynr], lh);
                    ProgressOutput::SumUpLocal();
                  } );
            }
        else
          for (FlatArray<int> blocks_of_col : block_coloring)
            for (int b : blocks_of_col)
              do_block (b, clh);
        return;
      }

    if (auto *task_manager = GetTaskManager(); task_manager)
      {
        for (FlatArray<int> els_of_col : element_coloring)
          {
            SharedLoop2 sl(els_of_col.Range());

            task_manager -> CreateJob
              ( [&] (const TaskInfo & ti) 
                {
                  LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);

                  for (int mynr : sl)
                    {
                      HeapReset hr(lh);
                      FESpace::Element el(fes, 
                                          ElementId (vb, els_of_col[mynr]), 
                                          lh);
                      
                      func (std::move(el), lh);
                    }

                  ProgressOutput::SumUpLocal();
                } );
          }
        return;
      }


    // no task manager: serial, the first exception stops the loop and is passed on unchanged
    for (FlatArray<int> els_of_col : element_coloring)
      for (int nr : els_of_col)
        {
          HeapReset hr(clh);
          func (FESpace::Element(fes, ElementId (vb, nr), clh), clh);
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


  void FESpace :: Interpolate (const CoefficientFunction & cf, BaseVector & vec1,
                               const Region * reg, LocalHeap & clh, int bonus_intorder)
  {
    static Timer t("FESpace::Interpolate"); RegionTimer r(t);

    shared_ptr<MeshAccess> ma = GetMeshAccess(); 
    int dim   = GetDimension();
    // dual shapes evaluate cf also on element facets
    ngfem::CheckContext (cf, MakeCFContext (*ma, CFContext::EVALUATE, reg ? reg->VB() : VOL, BND, false, IsComplex()));

    Array<int> cnti(GetNDof());
    cnti = 0;


    auto tempvec = vec1.CreateVector();
    
    if (reg)
      {
        auto regdofs = make_shared<BitArray> (GetDofs(*reg));
        Projector proj(regdofs, false);
        tempvec = proj * vec1;
      }
    else
      tempvec = 0.0;

   
    // complex: interpolate real and imaginary part separately
    shared_ptr<CoefficientFunction> cfre, cfim;
    if (vec1.IsComplex() && cf.IsComplex())
      {
        auto spcf = const_cast<CoefficientFunction&>(cf).shared_from_this();
        cfre = Real(spcf);
        cfim = Imag(spcf);
      }

    auto vb = reg ? reg->VB() : VOL;
    if (auto eval = GetEvaluator(vb))
      if (cf.Dimension() != eval->Dim())
        throw Exception ("Interpolate: coefficient function has dimension "+ToString(cf.Dimension())
                         +", the space "+GetClassName()+" expects "+ToString(eval->Dim()));
    auto iterate = [&] (auto tscal)
    {
      using SCAL = decltype(tscal);
      IterateElements
        (*this, vb, clh,
         [&] (FESpace::Element el, LocalHeap & lh)
         {
           if (reg)
             if (!reg->Mask().Test(el.GetIndex())) return;

           if (!DefinedOn(el)) return;

           const FiniteElement & fel = GetFE (el, lh);
           int ndof = fel.GetNDof();
           int fesdim = GetDimension();

           const ElementTransformation & eltrans = ma->GetTrafo (el, lh);
           FlatVector<SCAL> elvec(ndof*fesdim, lh), elvec1(ndof*fesdim, lh);

           if constexpr (is_same_v<SCAL,double>)
             fel.Interpolate (eltrans, cf, elvec.AsMatrix(ndof, fesdim), lh, bonus_intorder);
           else
             {
               FlatVector<> re(ndof*fesdim, lh), im(ndof*fesdim, lh);
               fel.Interpolate (eltrans, cfre ? *cfre : cf, re.AsMatrix(ndof, fesdim), lh, bonus_intorder);
               im = 0.0;
               if (cfim)
                 fel.Interpolate (eltrans, *cfim, im.AsMatrix(ndof, fesdim), lh, bonus_intorder);
               elvec = re + Complex(0,1) * im;
             }

           TransformVec (el, elvec, TRANSFORM_SOL_INVERSE);

           tempvec.GetIndirect (el.GetDofs(), elvec1);
           elvec1 += elvec;
           tempvec.SetIndirect (el.GetDofs(), elvec1);

           for (auto d : el.GetDofs())
             if (IsRegularDof(d)) cnti[d]++;
         });
    };
    if (vec1.IsComplex())
      iterate (Complex());
    else
      iterate (double());

#ifdef PARALLEL
    if (GetParallelDofs())
      {
        GetParallelDofs()->AllReduceDofData (cnti, NG_MPI_SUM);
        tempvec.SetParallelStatus(DISTRIBUTED);
        tempvec.Cumulate();
      }
#endif

    ParallelForRange
      (cnti.Size(), [&] (IntRange r)
       {
         auto average = [&] (auto & fluxi)
         {
           ArrayMem<int,1> dnums(1);
           for (auto i : r)
             if (cnti[i])
               {
                 dnums[0] = i;
                 tempvec.GetIndirect (dnums, fluxi);
                 fluxi /= double (cnti[i]);
                 tempvec.SetIndirect (dnums, fluxi);
               }
         };
         if (vec1.IsComplex())
           {
             VectorMem<10,Complex> fluxi(dim);
             average (fluxi);
           }
         else
           {
             VectorMem<10> fluxi(dim);
             average (fluxi);
           }
       });

    vec1 = tempvec;
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
    throw Exception ("SolveM is only available for L2-type spaces, not for " + GetClassName());
  }

  void FESpace :: ApplyM (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                          LocalHeap & lh) const
  {
    throw Exception ("ApplyM is only available for L2-type spaces, not for " + GetClassName());
  }


  shared_ptr<FESpace> FESpace :: CreateMultiplierSpace (Flags flags) const
  {
    throw Exception("CreateMultiplierSpace not available for "+GetClassName());
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
    if (l2space->GetMeshAccess() != ma)
      throw Exception ("ConvertL2Operator: target space is defined on another mesh");
    for (const FESpace * fes : { this, (const FESpace*)l2space.get() })
      if (!fes->GetEvaluator(VOL) || fes->GetEvaluator(VOL)->OnlyOnFacets())
        throw Exception ("ConvertL2Operator: " + fes->GetClassName() + " has no volume evaluator");
    if (GetEvaluator(VOL)->Dim() != l2space->GetEvaluator(VOL)->Dim())
      throw Exception ("ConvertL2Operator: spaces have different value dimensions");
    LocalHeap lh(100000000);
    
    Array<int> classnr(ma->GetNE());
    
    FlatArray<int> vertex_map;
    if (const PeriodicFESpace * periodic = dynamic_cast<const PeriodicFESpace*> (this))
      vertex_map.Assign (periodic->GetVertexOrientation());
    
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         if (!DefinedOn (el) || !l2space->DefinedOn (el)) { classnr[el.Nr()] = -1; return; }
         if (vertex_map.Size())
           classnr[el.Nr()] = ElementClassNr (el.GetType(), vertex_map[el.Vertices()]);
         else
           classnr[el.Nr()] = ElementClassNr (el.GetType(), el.Vertices());
       });
    CompressClassNrs (classnr);

    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        if(classnr[i] >= 0)
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
                                                false, false, evaluator);
        auto trial_l2 = make_shared<ProxyFunction>(l2space, 
                                                   false, false, l2evaluator);
        auto test_l2  = make_shared<ProxyFunction>(l2space,
                                                   true, false, l2evaluator);
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

  BitArray FESpace :: GetDofs (const Region & reg, const DifferentialOperator * diffop) const
  {
    if (reg.Mesh() != ma)
      throw Exception ("GetDofs: region is defined on another mesh");
    BitArray ba(GetNDof());
    ba.Clear();
    for (auto el : Elements(reg.VB()))
      if (reg.Mask().Test(el.GetIndex()))
        for (auto d : el.GetDofs())
          if (IsRegularDof(d))
            ba.SetBit(d);

    if (diffop)
      {
        IntRange r = GetRange(*diffop);
        for (size_t i = 0; i < r.First(); i++)
          ba.Clear(i);
        for (size_t i = r.Next(); i < ba.Size(); i++)
          ba.Clear(i);
      }
    return ba;
  }

  BitArray FESpace :: GetDofs (const MixedRegion & reg, const DifferentialOperator * diffop) const
  {
    if (reg.Mesh() != ma)
      throw Exception ("GetDofs: region is defined on another mesh");
    BitArray ba(GetNDof());
    ba.Clear();
    for (auto vb : { VOL, BND, BBND, BBBND })
      if (int(vb) <= ma->GetDimension() && reg.Mask(vb).NumSet())
        ba.Or (GetDofs(reg[vb], diffop));
    return ba;
  }

  IntRange FESpace :: GetRange (const DifferentialOperator & diffop) const
  {
    return IntRange(0, GetNDof());
  }

  ProxyNode FESpace ::
  MakeProxyFunction (bool testfunction,
                     const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const
  {
    shared_ptr<FESpace> fes = dynamic_pointer_cast<FESpace> (const_cast<FESpace*>(this)->shared_from_this());

    EvaluatorTable evaluators, deriv_evaluators;
    for (VorB vb : { VOL, BND, BBND, BBBND })
      {
        evaluators[vb] = fes->GetEvaluator(vb);
        deriv_evaluators[vb] = fes->GetFluxEvaluator(vb);
      }
    auto proxy = make_shared<ProxyFunction>  (fes, testfunction, testfunction ? false : fes->IsComplex(),
                                              evaluators, deriv_evaluators);
    
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
    needs_transform_vec = false;
    name="NodalFESpace";
    type = "nodal";
    if (order != 1 && order != 2)
      throw Exception ("NodalFESpace supports order 1 and 2, got order " + ToString(order));

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
        if (flux_evaluator[BND])
          flux_evaluator[BND] = make_shared<BlockDifferentialOperator> (flux_evaluator[BND], dimension);
      }


  }

  NodalFESpace :: ~NodalFESpace ()
  {
    ;
  }

  // local edge of each midpoint shape of the second-order element; empty: linear element
  static FlatArray<int> NodalMidpointEdges (ELEMENT_TYPE et, bool hb)
  {
    static int segm[] = { 0 };
    static int trig[] = { 1, 0, 2 };
    static int quad[] = { 0, 1, 2, 3 };
    static int tet[] = { 3, 4, 0, 5, 1, 2 };
    static int hex[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
    switch (et)
      {
      case ET_SEGM: return FlatArray<int>(1, segm);
      case ET_TRIG: return FlatArray<int>(3, trig);
      case ET_TET: return FlatArray<int>(6, tet);
      case ET_QUAD: return hb ? FlatArray<int>(0, nullptr) : FlatArray<int>(4, quad);
      case ET_HEX: return hb ? FlatArray<int>(0, nullptr) : FlatArray<int>(12, hex);
      default: return FlatArray<int>(0, nullptr);
      }
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
              case ET_HEX:     return *(new (lh) ScalarFE<ET_HEX,1>);
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
                 for (int k : NodalMidpointEdges(el.GetType(), hb_defined))
                   used_edge[el.Edges()[k]] = true;
	     });
      }

    // if (ma->GetNLevels() > ndlevel.Size())
      {
	size_t ndof = ma->GetNV();
        if (order > 1)
          ndof += ma->GetNEdges();
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

    if(order > 1)
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
    auto el = ma->GetElement(ei);
    dnums = el.Vertices();
    if (order > 1)
      {
        auto edges = el.Edges();
        size_t nv = ma->GetNV();
        for (int k : NodalMidpointEdges(el.GetType(), hb_defined))
          dnums.Append (nv + edges[k]);
      }
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
    if (order > 1 && used_edge[ednr])
      dnums.Append (ma->GetNV() + ednr);
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
    needs_transform_vec = false;
    name="NonconformingFESpace(nonconforming)";
    // defined flags
    DefineDefineFlag("nonconforming");
    if (parseflags) CheckFlags(flags);
    
    // prol = new LinearProlongation(*this);
    

    // trig = new FE_NcTrig1;
    // segm = new FE_Segm0;
      
    if (ma->GetDimension() == 2)
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
      }
    else
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
      }

    if (dimension > 1)
      for (auto vb : { VOL, BND })
        {
          if (evaluator[vb])
            evaluator[vb] = make_shared<BlockDifferentialOperator> (evaluator[vb], dimension);
          if (flux_evaluator[vb])
            flux_evaluator[vb] = make_shared<BlockDifferentialOperator> (flux_evaluator[vb], dimension);
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
    needs_transform_vec = false;
    name="NonconformingSurfaceFESpace(nonconforming)";
    // defined flags
    DefineDefineFlag("nonconforming");
    if (parseflags) CheckFlags(flags);
    
    // prol = new LinearProlongation(*this);
    

    // trig = new FE_NcTrig1;
    // segm = new FE_Segm0;
      
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
      }

    if (dimension > 1)
      for (auto vb : { VOL, BND })
        {
          if (evaluator[vb])
            evaluator[vb] = make_shared<BlockDifferentialOperator> (evaluator[vb], dimension);
          if (flux_evaluator[vb])
            flux_evaluator[vb] = make_shared<BlockDifferentialOperator> (flux_evaluator[vb], dimension);
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
    needs_transform_vec = false;
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
    if (ma->GetDimension() == 2)
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    else
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
    
    if (dimension > 1)
      evaluator[VOL] = make_shared<BlockDifferentialOperator> (evaluator[VOL], dimension);
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
    needs_transform_vec = false;
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

    switch (ma->GetDimension())
      {
      case 1: evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<1>>>(); break;
      case 2: evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>(); break;
      default: evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
      }

    if (dimension > 1)
      evaluator[BND] = make_shared<BlockDifferentialOperator> (evaluator[BND], dimension);
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
                    { return * new (lh) ScalarDummyFE<type.ElementType()>(); });
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

    needs_dof_transformation = false;
    for (auto & space : spaces)
      if (space->NeedsDofTransformation())
        needs_dof_transformation = true;

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
            if (!emb_i)  // component space has no embedding
              { sum_emb.reset(); break; }
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


  }

  void CompoundFESpace :: UpdateFreeDofs()
  {
    // dirichlet-dofs from sub-spaces
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
    if (GetEvaluator(VOL) || GetEvaluator(BND))
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
                        auto component = [i] (EvaluatorTable table)
                        {
                          for (int j = 0; j < 4; j++)
                            if (table[j])
                              table[j] = make_shared<CompoundDifferentialOperator> (table[j], i);
                          return table;
                        };
                        auto block_proxy = make_shared<ProxyFunction> (fes, testfunction, fes->IsComplex(),
                                                                       component(proxy->Evaluators()),
                                                                       component(proxy->DerivEvaluators()));
                        
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

  void CompoundFESpace :: SelectDofs (const string & name, BitArray & dofs) const
  {
    dofs.Clear();
    for (int i = 0; i < spaces.Size(); i++)
      {
        IntRange r = GetRange(i);
        BitArray dofsi(r.Size());
        spaces[i]->SelectDofs(name, dofsi);

        for (size_t j = 0; j < r.Size(); j++)
          if (dofsi.Test(j))
            dofs.SetBit(r.First()+j);
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

  optional<FlatMatrix<double>> CompoundFESpace :: VGetDofTransformation (ElementId ei, LocalHeap & lh) const
  {
    size_t ns = spaces.Size();
    ArrayMem<optional<FlatMatrix<double>>,8> Ti(ns);
    bool any = false;
    for (size_t i = 0; i < ns; i++)
      {
        Ti[i] = spaces[i]->GetDofTransformation (ei, lh);
        if (Ti[i]) any = true;
      }
    if (!any) return nullopt;

    // block diagonal, components without rectangular map give their square transformation
    ArrayMem<size_t,8> nloc(ns), nglob(ns);
    size_t sumloc = 0, sumglob = 0;
    Array<DofId> dnums;
    for (size_t i = 0; i < ns; i++)
      {
        nloc[i] = spaces[i]->GetFE(ei, lh).GetNDof();
        if (Ti[i])
          nglob[i] = Ti[i]->Width();
        else
          {
            spaces[i]->GetDofNrs (ei, dnums);
            nglob[i] = dnums.Size();
          }
        sumloc += nloc[i];
        sumglob += nglob[i];
      }

    FlatMatrix<double> T(sumloc, sumglob, lh);
    T = 0.0;
    for (size_t i = 0, r = 0, c = 0; i < ns; r += nloc[i], c += nglob[i], i++)
      {
        auto block = T.Rows(r, r+nloc[i]).Cols(c, c+nglob[i]);
        if (Ti[i])
          block = *Ti[i];
        else
          {
            if (nloc[i] != nglob[i])
              throw Exception ("CompoundFESpace::VGetDofTransformation: components with dimension > 1 not supported");
            FlatMatrix<double> Ii(nloc[i], nloc[i], lh);
            Ii = Identity(nloc[i]);
            spaces[i]->TransformMat (ei, Ii, TRANSFORM_MAT_LEFT);
            block = Trans(Ii);
          }
      }
    return T;
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
    if (!compoundtrace || compoundtrace->Spaces().Size() != Spaces().Size())
      throw Exception ("TraceOperator: trace space must be a vector space with " + ToString(Spaces().Size()) + " components");
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
      definedon_region.Mask(vb) = space->DefinedOn(vb);
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
    : fes(afes), rho(arho), inverse(ainverse), definedon(adefinedon), lh(alh)
  {
    CheckMassDensity (*fes, rho.get(), definedon.get());
  }

  void CheckMassDensity (const FESpace & fes, const CoefficientFunction * rho, const Region * definedon)
  {
    if (!rho) return;
    CheckContext (*rho, MakeCFContext (*fes.GetMeshAccess(), CFContext::EVALUATE,
                                       definedon ? definedon->VB() : VOL, VOL, false, fes.IsComplex()));
  }
  
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
  
  VecFormat FESpace :: GetVectorFormat () const
  {
    VecFormat f (GetNDof(), iscomplex ? Scalar(Complex(0)) : Scalar(double(0)), dimension);
    if (IsParallel())
      f.pardofs = GetParallelDofs();
    return f;
  }

  VecFormat ApplyMass :: RowFormat () const { return fes->GetVectorFormat(); }
  VecFormat ApplyMass :: ColFormat () const { return fes->GetVectorFormat(); }
  VecFormat ApplyTrace :: RowFormat () const { return fes->GetVectorFormat(); }
  VecFormat ApplyTrace :: ColFormat () const
  { return fes->GetVectorFormat().WithSize (festrace->GetNDof()); }

  



  
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

  shared_ptr<FESpace> CreateCompanionSpace (const FESpace & base, const string & type, Flags flags)
  {
    if (base.IsComplex()) flags.SetFlag("complex");
    if (base.DoesAutoUpdate()) flags.SetFlag("autoupdate");
    auto fes = CreateFESpace (type, base.GetMeshAccess(), flags);
    if (!flags.AnyFlagDefined("definedon"))
      for (VorB vb : { VOL, BND })
        if (!base.DefinedOnAll(vb))
          fes->SetDefinedOn(vb, base.DefinedOn(vb));
    fes->Update();
    fes->FinalizeUpdate();
    return fes;
  }



  // standard fespaces:

  static RegisterClassForArchive<FESpace> regfes;
  static RegisterFESpace<NodalFESpace> initnodalfes ("nodal");
  static RegisterFESpace<NonconformingFESpace> initncfes ("nonconforming");
  static RegisterFESpace<NonconformingSurfaceFESpace> initncsurfes ("nonconformingsurface");

}



