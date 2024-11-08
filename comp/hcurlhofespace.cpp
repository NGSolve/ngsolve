/**
   High Order Finite Element Space for H(Curl) 
*/
// #include <comp.hpp>
#include "hcurlhofespace.hpp"
#include "hcurlhdivfes.hpp"
#include <../fem/hcurlhofe.hpp> 
#include <../fem/hcurllofe.hpp>
#include <../fem/hcurl_equations.hpp> 
#include <diffop_impl.hpp>

// #include <../fem/hcurlhdiv_dshape.hpp> 
#include <multigrid.hpp>

// extern template class ngla::VFlatVector<double>;

namespace ngcomp 
{
  

  HCurlHighOrderFESpace ::  
  HCurlHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & aflags, bool parseflags)
    : FESpace (ama, aflags),
      flags(aflags)
  {
    type = "hcurlho";
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
    DefineNumFlag("augmented");
    DefineDefineFlag("fast"); 
    DefineDefineFlag ("discontinuous");
    DefineDefineFlag ("type1");
    
    if(parseflags) CheckFlags(flags);

    type1 = flags.GetDefineFlag("type1");
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
      
    int ndom = ma->GetNDomains();
    gradientdomains.SetSize (ndom);
    gradientdomains.Set();
    
    gradientboundaries.SetSize (ma->GetNBoundaries());
    gradientboundaries.Set();

    fn=int(flags.GetNumFlag("face",1));

    if (flags.NumListFlagDefined("gradientdomains"))
      {
	const Array<double> & graddomains = flags.GetNumListFlag ("gradientdomains");
	for (int i = 0; i < gradientdomains.Size(); i++)
	  if (!graddomains[i])
	    gradientdomains.Clear(i);
      }

    if (flags.StringFlagDefined("gradientdomains"))
      {
        Region gd(ma, VOL, flags.GetStringFlag("gradientdomains"));
        gradientdomains = gd.Mask(); 
      }
    
    if (flags.NumListFlagDefined("gradientboundaries"))
      {
	const Array<double> & gradbounds = flags.GetNumListFlag ("gradientboundaries");
	for (int i = 0; i < gradientboundaries.Size(); i++)
	  if (!gradbounds[i])
	    gradientboundaries.Clear(i);
      }

    if (flags.StringFlagDefined("gradientboundaries"))
      {
        Region gd(ma, BND, flags.GetStringFlag("gradientboundaries"));
        gradientboundaries = gd.Mask(); 
      }

    
    if(flags.GetDefineFlag("nograds"))
      {
	gradientdomains.Clear();   
	gradientboundaries.Clear();   
      } 
       
    fast_pfem = flags.GetDefineFlag ("fast");
    discontinuous = flags.GetDefineFlag ("discontinuous");
    highest_order_dc = flags.GetDefineFlag ("highest_order_dc");    
    if (discontinuous)
      SetDefinedOn(BND, BitArray(ma->GetNRegions(BND)).Clear());      

    
    if (flags.GetDefineFlag ("no_couplingtype_upgrade"))
      ctupgrade = false;
    Flags loflags = flags;
    loflags.SetFlag ("order", 1);
    /*
    loflags.SetFlag ("dim", dimension);
    if (iscomplex) loflags.SetFlag ("complex");
    if (discontinuous) loflags.SetFlag ("disontinuous");
    if (flags.NumListFlagDefined ("dirichlet")) 
      loflags.SetFlag ("dirichlet", flags.GetNumListFlag ("dirichlet"));
    */
    low_order_space = make_shared<NedelecFESpace> (ma, loflags);
    prol = make_shared<ngmg::EdgeProlongation> 
      (*static_cast<NedelecFESpace*> (low_order_space.get()));
   
    uniform_order_inner = int (flags.GetNumFlag ("orderinner", -1));
    uniform_order_face = int (flags.GetNumFlag ("orderface", -1));
    uniform_order_edge = int (flags.GetNumFlag ("orderedge", -1));

    wb_loedge = flags.GetDefineFlag("wb_loedge");
        
    if (flags.NumFlagDefined("smoothing")) 
      throw Exception ("Flag 'smoothing' for fespace is obsolete \n Please use flag 'blocktype' in preconditioner instead");
    if (flags.NumFlagDefined("cluster")) 
      throw Exception ("Flag 'cluster' for fespace is obsolete \n Please use flag 'ds_cluster' in preconditioner instead");
       
    augmented = int (flags.GetNumFlag ("augmented", 0));


    // Evaluator
    /*
    static ConstantCoefficientFunction one(1);
    integrator[VOL] = GetIntegrators().CreateBFI("massedge", ma->GetDimension(), &one);
    if (!discontinuous)
      integrator[BND] = GetIntegrators().CreateBFI("robinedge", ma->GetDimension(), &one); 
    */
    
    if (ma->GetDimension() == 2)
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<2>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<2>>>();
      }
    else
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<3>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<3>>>();
        flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpCurlBoundaryEdgeVec<>>>();
	evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdBBoundaryEdge<3>>>();
      }

    switch (ma->GetDimension())
      {
      case 1:
        additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHCurl<1>>> ()); break;
      case 2:
        additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHCurl<2>>> ());
        additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHCurlDual<2>>> ());
        break;
      case 3:
        additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHCurl<3>>> ());
        additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHCurlDual<3>>> ());
        break;
      default:
        ;
      }

    this->GetMemoryTracer().Track(
        order_edge, "order_edge",
        fine_edge, "fine_edge",
        fine_face, "fine_face",
        cell_ngrad, "cell_ngrad",
        face_ngrad, "face_ngrad",
        order_face, "order_face",
        order_inner, "order_inner",
        order_avertex, "order_avertex",
        usegrad_edge, "usegrad_edge",
        usegrad_face, "usegrad_face",
        usegrad_cell, "usegrad_cell",
        dom_order_min, "dom_order_min",
        dom_order_max, "dom_order_max"
        );
  }
  
  HCurlHighOrderFESpace :: ~HCurlHighOrderFESpace () { ; }
  
  DocInfo HCurlHighOrderFESpace :: GetDocu ()
  {
    auto docu = FESpace::GetDocu();
    docu.Arg ("nograds") = "bool = False\n"
      "  Remove higher order gradients of H1 basis functions from HCurl FESpace";
    docu.Arg("type1") = "bool = False\n"
      "  Use type 1 Nedelec elements";
    docu.Arg("discontinuous") = "bool = False\n"
      "  Create discontinuous HCurl space";
    docu.Arg("gradientdomains") = "List[int] = None\n"
      "  Remove high order gradients from domains where the value is 0.\n"
      "  This list can be generated for example like this:\n"
      "  graddoms = [1 if mat == 'iron' else 0 for mat in mesh.GetMaterials()]";
    docu.Arg("gradientboundaries") = "List[int] = None\n"
      "  Remove high order gradients from boundaries where the value is 0.\n"
      "  This list can be generated for example like this:\n"
      "  gradbnds = [1 if bnd == 'iron_bnd' else 0 for bnd in mesh.GetBoundaries()]";
    docu.Arg("highest_order_dc") = "bool = False\n"
      "  Activates relaxed H(curl)-conformity. Allows tangential discontinuity of highest order edge basis functions";
    return docu;
  }

  
  void HCurlHighOrderFESpace :: Update()
  {
    FESpace::Update();

    int dim = ma->GetDimension(); 

    if (order < 0) 
      throw Exception("HCurlHighOrderFESpace::Update() order < 0 !" ) ;
    
    if (low_order_space)
      low_order_space -> Update();

    bool first_update = GetTimeStamp() < ma->GetTimeStamp();
    if (first_update) timestamp = NGS_Object::GetNextTimeStamp();

    int ne = ma->GetNE();
    int nse = ma->GetNSE();
    int ned = ma->GetNEdges();
    int nfa = (ma->GetDimension() == 3 || order_policy == VARIABLE_ORDER) ? ma->GetNFaces() : 0;

    maxorder = -1; 
    minorder = 99; 

    /*
    for(int i = 0; i < specialelements.Size(); i++)
      delete specialelements[i];
    specialelements.DeleteAll();
    */
    
    /*
    if (order_policy == VARIABLE_ORDER &&
        ma->GetTimeStamp() > order_timestamp)
      {
        FESpace::order_edge.SetSize(ned);
        FESpace::order_face_left.SetSize(ma->GetNNodes(NT_FACE));
        FESpace::order_face_right.SetSize(ma->GetNNodes(NT_FACE));
        FESpace::order_cell_left.SetSize(ma->GetNNodes(NT_CELL));
        FESpace::order_cell_right.SetSize(ma->GetNNodes(NT_CELL));
        FESpace::order_edge = order;
        FESpace::order_face_left = order;
        FESpace::order_face_right = order;
        FESpace::order_cell_left = order;
        FESpace::order_cell_right = order;
        order_timestamp = ma->GetTimeStamp();
      }
    */

    if (first_update)
      {
        order_edge.SetSize (ned);   
        order_face.SetSize (nfa); 
        order_inner.SetSize (ne);
        usegrad_edge.SetSize (ned);                                
        usegrad_face.SetSize (nfa); 
        usegrad_cell.SetSize (ne);
        fine_edge.SetSize (ned); 
        fine_face.SetSize (nfa); 
        
        int p = var_order ? 0 : order; 
        order_edge = max(0, p - (type1 ? 1 : 0) + et_bonus_order[ET_SEGM]);

        // order_inner = IVec<3> (p,p,p); 
        order_inner = IVec<3> (0,0,0); 
        
        fine_edge = 0; 
        // if (nfa > 0)
        { 
          fine_face = 0;
          if (et_bonus_order[ET_TRIG] || et_bonus_order[ET_QUAD])
            for (auto f : Range(nfa))
              order_face[f] = p+et_bonus_order[ma->GetFaceType(f)];
          else
            order_face = IVec<2> (p,p);
          usegrad_face = 0; 
        } 
        
        usegrad_edge = false;                                
        usegrad_cell = false; 
        
        Array<int> elfaces;
        
        /*
          for(int i = 0; i < ne; i++) 
          if(gradientdomains[ma->GetElIndex(i)]) 
          {
	  ma->GetElEdges(i,eledges);
	  
	  for(int j=0;j<eledges.Size();j++)
          usegrad_edge[eledges[j]]=1;
	  
	  if(ma->GetDimension()==3)
          {
          ma->GetElFaces(i, elfaces);
          for(int j = 0; j < elfaces.Size(); j++)
          usegrad_face[elfaces[j]]=1; 
          }
	  usegrad_cell[i] = 1; 
          } 
        */
        
        /*
          for (Ngs_Element el : ma->Elements(VOL))
          if (gradientdomains[el.GetIndex()]) 
          {
          usegrad_edge[el.Edges()] = true;
          if (ma->GetDimension() == 3)
          usegrad_face[el.Faces()] = true;
          usegrad_cell[el.Nr()] = true;
          }
        */
        ma -> IterateElements
          (VOL,
           [&](auto el)
           {
             if (gradientdomains[el.GetIndex()]) 
               {
                 usegrad_edge[el.Edges()] = true;
                 // if (ma->GetDimension() == 3)
                 if (nfa)
                   usegrad_face[el.Faces()] = true;
                 usegrad_cell[el.Nr()] = true;
               }
           });
        
        if (gradientboundaries.Size())
          // for (int i = 0; i < nse; i++)
          for (ElementId ei : ma->Elements(BND))
            if (gradientboundaries[ma->GetElIndex(ei)])
              {
                auto eledges = ma->GetElEdges(ei);
                for(int j=0; j<eledges.Size();j++)
                  usegrad_edge[eledges[j]] = 1;
                
                if(ma->GetDimension()==3)
                  usegrad_face[ma->GetSElFace(ei.Nr())] = 1;
              }
        
        ma->AllReduceNodalData (NT_EDGE, usegrad_edge, NG_MPI_LOR);
        ma->AllReduceNodalData (NT_FACE, usegrad_face, NG_MPI_LOR);
        
	
        for (int i = 0; i < ne; i++)
          {
            ElementId ei(VOL, i);
            int index = ma->GetElIndex(ei);
            auto et = ma->GetElType(ei);
            if (!DefinedOn (VOL, index)) continue;
            auto pi = p + et_bonus_order[et];
            if (type1 && et==ET_QUAD && pi >= 1) pi--;
            order_inner[i] = IVec<3> (pi,pi,pi);
            IVec<3,TORDER> el_orders = ma->GetElOrders(i);
            
            ELEMENT_TYPE eltype=ma->GetElType(ei); 
            const FACE * faces = ElementTopology::GetFaces (eltype);
            const EDGE * edges = ElementTopology::GetEdges (eltype);
            const POINT3D * points = ElementTopology :: GetVertices (eltype);
            auto vnums = ma->GetElVertices (ei);
            
            auto eledges = ma->GetElEdges (ei);		
            for(int j=0;j<eledges.Size();j++) fine_edge[eledges[j]] = 1; 
            
            if (dim == 3)
              {
                auto elfaces = ma->GetElFaces(ei);
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
                    
                    
                    IVec<2> f((fmax+3)%4,(fmax+1)%4); 
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
        
        for (int i = 0; i < nse; i++)
          {
            ElementId sei(BND, i);
            if (!DefinedOn (BND, ma->GetElIndex (sei))) continue;
            
            // auto eledges = ma->GetElEdges (sei);		
            // for (int j=0;j<eledges.Size();j++) fine_edge[eledges[j]] = 1;
            fine_edge[ma->GetElEdges(sei)] = true;
            if(dim==3) 
              fine_face[ma->GetSElFace(i)] = true; 
          }
        
        ma->AllReduceNodalData (NT_EDGE, fine_edge, NG_MPI_LOR);
        ma->AllReduceNodalData (NT_FACE, fine_face, NG_MPI_LOR);
        
      	
        if(!var_order) { maxorder = order; minorder = order;} 
    
        if(uniform_order_inner>-1) 
          order_inner = IVec<3> (uniform_order_inner,uniform_order_inner,uniform_order_inner); 
        if(uniform_order_edge>-1) 
          order_edge = uniform_order_edge; 
        if(uniform_order_face>-1 && dim == 3) 
          order_face = IVec<2> (uniform_order_face, uniform_order_face);
        /*
          Array<int> pnums;
          for (int i = 0; i < order_face.Size(); i++)
          {
          ma->GetFacePNums (i,pnums);  
          if (pnums.Size()==4)
          order_face[i] = uniform_order_face+1;
          }
        */
        // order of FINE FACES and EDGES for safety reasons set to 0 
        for(int i=0;i<order_edge.Size();i++) 
          if(!fine_edge[i]) order_edge[i] = 0;  
        
        for(int i=0;i<order_face.Size();i++) 
          if(!fine_face[i]) order_face[i] = IVec<2> (0,0);  
      }

    UpdateDofTables(); 
    UpdateCouplingDofArray();
    if (low_order_space)
      low_order_embedding =
        make_shared<Embedding> (GetNDof(),
                                IntRange(low_order_space->GetNDof()),
                                IsComplex());

  }
		
  void HCurlHighOrderFESpace :: DoArchive(Archive & archive)
  {
    low_order_space -> DoArchive (archive);
    FESpace::DoArchive(archive);
    archive & level;
    archive & first_edge_dof & first_inner_dof & first_face_dof;
    archive & fn & rel_order & rel_orders;
    archive & order_edge & fine_edge & fine_face;
    archive & cell_ngrad & face_ngrad & order_face & order_inner & order_avertex;
    archive & usegrad_edge & usegrad_face & usegrad_cell;
    archive & dom_order_min & dom_order_max;
    archive & maxorder & minorder;
    archive & gradientdomains & gradientboundaries;
    archive & usegrad & var_order;
    archive & ndof & nedfine & uniform_order_inner & 
      uniform_order_face & uniform_order_edge & augmented;
    archive & flags;
    archive & smoother & level_adapted_order & nograds;
    archive & fast_pfem & discontinuous;
  }


  void HCurlHighOrderFESpace :: UpdateDofTables()
  {
    // if (order_policy == VARIABLE_ORDER)
    if(false)
      {
        int dim = ma->GetDimension();
        size_t nedge = ma->GetNEdges(); 
        size_t nface = ma->GetNFaces();
        size_t ncell = ma->GetNE();

        ndof = nedge;
        
        first_edge_dof.SetSize (nedge+1); 
        for (auto i : Range(nedge))
          {
            first_edge_dof[i] = ndof;
            if(order_edge[i] > 0)
              ndof += order_edge[i];
          }
        first_edge_dof[nedge] = ndof;

        first_face_dof.SetSize (nface+1);
        for (auto i : Range(nface))
          {
            first_face_dof[i] = ndof;
            /*
            IVec<2> pl = FESpace::order_face_left[i];
            IVec<2> pr = FESpace::order_face_right[i];
            */
            IVec<2> pl = order_face[i];
            IVec<2> pr = order_face[i];
            int ngrad = 0, ncurl = 0;
            switch (ma->GetFaceType(i))
              {
              case ET_TRIG: 
                {
                  ngrad = pl[0]*(pl[0]-1)/2;
                  ncurl = (pr[0]+2)*(pr[0]-1)/2;
                  break;
                }
              case ET_QUAD:
                {
                  /*
                  ndof += (usegrad_face[i]+1)*p[0]*p[1] + p[0] + p[1]; 
                  face_ngrad[i] = usegrad_face[i]*p[0]*p[1];; 
                  */
                  break; 
                }
              default:
                __assume(false);
              }
            
            if (ngrad < 0) ngrad = 0;
            if (ncurl < 0) ncurl = 0;
            ndof += ngrad + ncurl;
          } 
        first_face_dof[nface] = ndof;


        if (dim == 3)
          {
            first_inner_dof.SetSize(ncell + 1);
            for (auto i : Range(ncell))
              {
                first_inner_dof[i] = ndof;
                /*
                  IVec<2> pl = FESpace::order_face_left[i];
                  IVec<2> pr = FESpace::order_face_right[i];
                */
                IVec<2> pl = order_inner[i];
                switch (ma->GetElType(ElementId(VOL,i)))
                  {
                  case ET_TET:
                    {
                      if (pl[0] > 2)
                        ndof += (pl[0] * pl[0] - 1)*(pl[0] - 2) / 2;
                      break;
                    }
                  case ET_HEX:
                    {
                      ndof += 3 * pl[0] * pl[0] * (pl[0] + 1);
                      /*
                        ndof += (usegrad_face[i]+1)*p[0]*p[1] + p[0] + p[1];
                        face_ngrad[i] = usegrad_face[i]*p[0]*p[1];;
                      */
                      break;
                    }
                  default:
                    __assume(false);
                  }
                
              }
            first_inner_dof[ncell] = ndof;
          }
        return;
      }



    
    int ne = ma->GetNE();

    int ned = ma->GetNEdges();
    int nfa = (ma->GetDimension() == 2) ? 0 : ma->GetNFaces();


    nedfine = 0; 
    for(int i = 0; i < ned; i++) 
      if (fine_edge[i] == 1) nedfine++; 

    ndof = ned; // Nedelec (order = 0) !!   
       
    first_edge_dof.SetSize (ned+1); 
    for (int i = 0; i < ned; i++)
      {
	first_edge_dof[i] = ndof;
        if (usegrad_edge[i])
          {
            int oe = order_edge[i];
            if (highest_order_dc) oe--;
            if (oe > 0) ndof += oe;
          }
      }
    first_edge_dof[ned] = ndof;
    
    first_face_dof.SetSize (nfa+1);
    face_ngrad.SetSize (nfa); 
    face_ngrad = 0; 
    for (int i = 0; i < nfa; i++) 
      { 
	first_face_dof[i] = ndof; 
	auto pnums = ma->GetFacePNums (i);  
	IVec<2> p = order_face[i]; 
	switch(pnums.Size())
	  {
	  case 3: //Triangle   
	    if (p[0]>1)
	      {
                /*
		ndof += ((usegrad_face[i]+1)*p[0] + 2)*(p[0]-1)/2;
		face_ngrad[i] = usegrad_face[i]*p[0]*(p[0]-1)/2;
                */
                int pg = p[0] - (type1 ? 1 : 0);
		face_ngrad[i] = usegrad_face[i]*pg*(pg-1)/2;
                ndof += face_ngrad[i] + (p[0] + 2)*(p[0]-1)/2;
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
    
    cell_ngrad.SetSize(ne); 
    cell_ngrad = 0; 
    first_inner_dof.SetSize(ne+1);
    for (int i = 0; i < ne; i++)
      {
        ElementId ei(VOL, i);
	first_inner_dof[i] = ndof;
	IVec<3> p = order_inner[i];
	switch(ma->GetElType(ei)) 
	  {
	  case ET_TRIG:
	    if(p[0]>1)
	      {
                int pg = p[0] - (type1 ? 1 : 0);
		cell_ngrad[i] = usegrad_cell[i]*pg*(pg-1)/2;
                ndof += cell_ngrad[i] + (p[0] + 2)*(p[0]-1)/2;
                /*
		ndof += ((usegrad_cell[i]+1)*p[0] + 2) * (p[0]-1) /2;
		cell_ngrad[i] = ((usegrad_cell[i])*p[0]) * (p[0]-1) /2;
                */
	      }
	    break; 
	  case ET_QUAD: 
	    if(p[0]>=0 && p[1]>=0) 
	      {
		ndof += (usegrad_cell[i]+1) * p[0] * p[1] + p[0] + p[1]; 
		cell_ngrad[i] = (usegrad_cell[i]) * p[0] * p[1];
	      }
	    break; 
	  case ET_TET: 
	    if(p[0]>2)
	      {
		if (type1) {
		  cell_ngrad[i] = usegrad_cell[i] * (p[0]-3)*(p[0]-2)*(p[0]-1)/6;
		  ndof += (p[0]-2)*(p[0]-1)*(2*p[0]+3)/6 + cell_ngrad[i];
				  
		}
		else {
		  ndof += ((usegrad_cell[i] + 2) *  p[0] + 3) * (p[0]-2) * (p[0]-1) / 6; 
		  cell_ngrad[i] = ((usegrad_cell[i] ) *  p[0]) * (p[0]-2) * (p[0]-1) / 6;
		}
	      }
	    break; 
	  case ET_PRISM:
	    if(p[0]>1 && p[2]>=0) 
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
        if (highest_order_dc)
          ndof += ElementTopology::GetNEdges(ma->GetElType(ei));
      }
    first_inner_dof[ne] = ndof;    

    if (discontinuous)
      {
	Array<int> faces;

	ndof = 0;
	for (int el = 0; el < ne; el++ )
	  {
            ElementId ei(VOL,el);
	    int ndof_inner = first_inner_dof[el+1] - first_inner_dof[el];
	    first_inner_dof[el] = ndof;
	    auto edges = ma->GetElEdges(ei);
	    if (ma->GetDimension() == 3)
	      faces = ma->GetElFaces(ei);

	    for ( int ed = 0; ed < edges.Size(); ed++ )
	      {
		int ndof_edge = 1 + first_edge_dof[edges[ed]+1] - first_edge_dof[edges[ed]];
		ndof += ndof_edge;
	      }

	    if ( ma->GetDimension() == 3 )
	      for ( int fa = 0; fa < faces.Size(); fa++ )
		{
		  int ndof_face = first_face_dof[faces[fa]+1] - first_face_dof[faces[fa]];
		  ndof += ndof_face;
		}
	    ndof += ndof_inner;
	  }
	first_inner_dof[ne] = ndof;

	first_edge_dof = 0;
	first_face_dof = 0;
      }

    SetNDof(ndof);
    *testout << "Hcurlho edge dofs: " << first_edge_dof[0] << "-" << first_edge_dof[ned] << endl;
    *testout << "Hcurlho face dofs: " << first_face_dof[0] << "-" << first_face_dof[nfa] << endl;
    *testout << "Hcurlho inner dofs: " << first_inner_dof[0] << "-" << first_inner_dof[ne] << endl;


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

    for (int edge = 0; edge < ma->GetNEdges(); edge++) 
      {
	ctofdof[edge] = 
	  fine_edge[edge] ? WIREBASKET_DOF : UNUSED_DOF; //Nedelec0

	IntRange range = GetEdgeDofs (edge);
	ctofdof.Range (range) = INTERFACE_DOF;
        if (wb_loedge && range.Size() >= 1)
	  ctofdof[range.First()] = WIREBASKET_DOF;
      }

    if (order_policy == VARIABLE_ORDER && ma->GetDimension() == 2)
      { 
        for (auto face : Range(ma->GetNFaces()))
          ctofdof.Range(GetFaceDofs(face)) = LOCAL_DOF;
        return;
      }
    
    // faces 
    if (ma->GetDimension() == 3)
      for (int face = 0; face < ma->GetNFaces(); face++)
	{
	  IntRange range = GetFaceDofs (face);
	  ctofdof.Range (range) = INTERFACE_DOF;
	  /*
	  if (ma->GetFacetType (face) == ET_QUAD)
	    {
	      IVec<2> p = order_face[face];
	      int hnext = range.Next();
	      int hfirst = hnext-p[0]-p[1];
	      ctofdof.Range (hfirst, hnext) = WIREBASKET_DOF;
	    }
	  */
	}

    
    LocalHeap lh(1000000, "HCurlHighOrderFESpace::UpdateCouplingDofArray");
    Array<int> dnums;
    for (int el = 0; el < ma->GetNE(); el++)
      {
        ElementId ei(VOL, el);
	if (!DefinedOn (VOL, ma->GetElIndex (ei))) continue;
	HeapReset hr(lh);
	IntRange range = GetElementDofs (el);
	ctofdof.Range (range) = LOCAL_DOF;

	bool upgrade = false;
	ELEMENT_TYPE eltype = ma->GetElType (ei);

	if (eltype == ET_PRISM) 
	  {
	    ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
	    IntegrationPoint ip(0.3333, 0.3333, 0.5);
	    MappedIntegrationPoint<3,3> mip(ip, eltrans);

	    Mat<3> jac = mip.GetJacobian();
	    double jaclong = L2Norm (jac.Col(2));
	    double jacplane = L2Norm (jac.Col(0)) + L2Norm (jac.Col(1));
	      
	    auto edge_nums = ma->GetElEdges (ei);
	    auto face_nums = ma->GetElFaces (ei);

	    // vertical edges
	    if (jaclong > 3 * jacplane)
	      {
		for (int j = 6; j < 9; j++)
		  {
		    int enr = edge_nums[j];
		    ctofdof.Range (GetEdgeDofs(enr)) = WIREBASKET_DOF;
		  }
		for (int j = 2; j < 5; j++)
		  {
		    int fnr = face_nums[j];
		    ctofdof.Range (GetFaceDofs(fnr)) = WIREBASKET_DOF;
		    /*
		    range = GetFaceDofs (fnr);
		    IVec<2> p = order_face[fnr];

		    int hnext = range.Next();
		    int hfirst = hnext-p[0]-p[1];
		    ctofdof.Range (hfirst, hnext) = WIREBASKET_DOF;
		    */
		  }
	      }



	    if (jaclong < 0.33 * jacplane)
	      {
		for (int j = 0; j < 6; j++)
		  {
		    int enr = edge_nums[j];
		    ctofdof.Range (GetEdgeDofs(enr)) = WIREBASKET_DOF;
		  }
		for (int j = 0; j < 2; j++)
		  {
		    int fnr = face_nums[j];
		    ctofdof.Range (GetFaceDofs(fnr)) = WIREBASKET_DOF;
		  }

		for (int j = 2; j < 5; j++)
		  {
		    int fnr = face_nums[j];
		    // ctofdof.Range (GetFaceDofs(fnr)) = WIREBASKET_DOF;
		    range = GetFaceDofs (fnr);
		    IVec<2> p = order_face[fnr];

		    int hnext = range.Next();
		    int hfirst = hnext-p[0]-p[1];
		    ctofdof.Range (hfirst, hnext) = WIREBASKET_DOF;
		  }
	      }
	  }

	if (eltype == ET_TET)
	  {
	    ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL, el), lh);
	    IntegrationPoint ip(0.25, 0.25, 0.25);
	    MappedIntegrationPoint<3,3> mip(ip, eltrans);

	    double cond = L2Norm (mip.GetJacobian()) * L2Norm (mip.GetJacobianInverse());
	    if (cond > 10) upgrade = true;
	  }
	
	if (eltype == ET_PYRAMID) upgrade = true;
	
	if (upgrade && ctupgrade)
	  {
	    GetDofNrs (el, dnums);
	    for (int j = 0; j < dnums.Size(); j++)
	      if (dnums[j] != -1 && ctofdof[dnums[j]] == INTERFACE_DOF)
		ctofdof[dnums[j]] = WIREBASKET_DOF;
	  }
      }
    

  }

  void HCurlHighOrderFESpace :: DoCouplingDofUpgrade(bool actupgrade) {
      ctupgrade = actupgrade;
      UpdateCouplingDofArray();
    }

  void HCurlHighOrderFESpace :: SetOrder (NodeId ni, int order) 
  {
    if (order_policy == CONSTANT_ORDER || order_policy == NODE_TYPE_ORDER)
      throw Exception("In HCurlHighOrderFESpace::SetOrder. Order policy is constant or node-type!");
    else if (order_policy == OLDSTYLE_ORDER)
      order_policy = VARIABLE_ORDER;
      
    if (order < 0)
      order = 0;
    
    switch (ni.GetType())
      {
      case NT_VERTEX: case NT_GLOBAL:
        break;
      case NT_EDGE:
        if (ni.GetNr() < order_edge.Size())
          order_edge[ni.GetNr()] = order;
        break;
      case NT_FACE:
        if (ni.GetNr() < order_face.Size())
          order_face[ni.GetNr()] = order;
        break;
      case NT_CELL: case NT_ELEMENT:
        if (ni.GetNr() < order_inner.Size())
          order_inner[ni.GetNr()] = order;
        break;
      case NT_FACET:
	//TODO
        break;
      }
  }
  
  int HCurlHighOrderFESpace :: GetOrder (NodeId ni) const
  {
    switch (ni.GetType())
      {
      case NT_VERTEX: case NT_GLOBAL:
        return 0;
      case NT_EDGE:
        if (ni.GetNr() < order_edge.Size())
          return order_edge[ni.GetNr()];
        break;
      case NT_FACE:
        if (ni.GetNr() < order_face.Size())
          return order_face[ni.GetNr()][0];
        break;
      case NT_CELL: case NT_ELEMENT:
        if (ni.GetNr() < order_inner.Size())
          return order_inner[ni.GetNr()][0];
        break;
      case NT_FACET:
        break;
      }
    return 0;
  }

  

  FiniteElement & HCurlHighOrderFESpace :: GetFE (ElementId ei, Allocator & lh) const
  {
    switch(ma->GetElType(ei))
      {
      case ET_POINT:   return * new (lh) DummyFE<ET_POINT>; 
      case ET_SEGM:    return T_GetFE<ET_SEGM> (ei, lh);
        
      case ET_TRIG:    return T_GetFE<ET_TRIG> (ei, lh);
      case ET_QUAD:    return T_GetFE<ET_QUAD> (ei, lh);
        
      case ET_TET:     return T_GetFE<ET_TET> (ei, lh);
      case ET_PRISM:   return T_GetFE<ET_PRISM> (ei, lh);
      case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (ei, lh);
      case ET_HEX:     return T_GetFE<ET_HEX> (ei, lh);

      default:
        throw Exception ("illegal element in HCurlHoFeSpace::GetFE");
      }
      
  }
  
  template <ELEMENT_TYPE ET>
  FiniteElement & HCurlHighOrderFESpace :: T_GetFE (ElementId ei, Allocator & lh) const
  {
    // if (order_policy == VARIABLE_ORDER)
    if(false)
      {
        // cout << "ei = " << ei << endl;
        Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL> (ei.Nr());
        if (!DefinedOn (ngel))
          return * new (lh) HCurlDummyFE<ET>();
        
        HCurlHighOrderFE<ET> * hofe =  new (lh) HCurlHighOrderFE<ET> ();
        
        hofe -> SetVertexNumbers (ngel.vertices);
        bool ta[12] = { true, true, true, true, true, true, true, true, true, true, true, true };
        hofe -> SetUseGradEdge (ta);
        hofe -> SetUseGradFace (ta);
        /*
        cout << "order-edge = " << int(FESpace::order_edge[ngel.Edges()[0]])
             << int(FESpace::order_edge[ngel.Edges()[1]])
             << int(FESpace::order_edge[ngel.Edges()[2]]) << endl;
        cout << "order-face = " << FESpace::order_face_right[ngel.Faces()] << endl;
        */
        hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
        hofe -> SetOrderFace (order_face[ngel.Faces()]);
        if (ma->GetDimension() == 3) {
          hofe->SetUseGradCell(true);
          hofe->SetOrderCell(order_inner[ngel.Nr()]);
        }
        hofe -> SetType1 (false);
        hofe -> ComputeNDof();
        // cout << "                                neldof = " << hofe->GetNDof() << ", order = " << hofe->Order() << endl;
        return *hofe;
      }

    
    switch(ei.VB())
      {
      case VOL:
        {
          Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL> (ei.Nr());
          if (!DefinedOn (ngel))
            return * new (lh) HCurlDummyFE<ET>();
          
          HCurlHighOrderFE<ET> * hofe =  new (lh) HCurlHighOrderFE<ET> ();
          
          hofe -> SetVertexNumbers (ngel.vertices);
          
          hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
          hofe -> SetUseGradEdge (usegrad_edge[ngel.Edges()]);
          
          switch (int(ET_trait<ET>::DIM))
            {
            case 1:
              throw Exception("no 1D elements in H(curl)");
            case 2:
              {
                hofe -> SetOrderCell (order_inner[ei.Nr()]);   // old style
                IVec<2,TORDER> p(order_inner[ei.Nr()][0], order_inner[ei.Nr()][1]);
                FlatArray<IVec<2,TORDER> > of(1, &p);
                hofe -> SetOrderFace (of);
                
                hofe -> SetUseGradCell (usegrad_cell[ei.Nr()]);  // old style
                FlatArray<bool> augf(1,&usegrad_cell[ei.Nr()]);
                hofe -> SetUseGradFace (augf); 
                break;
              }
            case 3:
              {
                hofe -> SetOrderFace (order_face[ngel.Faces()]);
                hofe -> SetUseGradFace (usegrad_face[ngel.Faces()]);
                
                hofe -> SetOrderCell (order_inner[ei.Nr()]);
                hofe -> SetUseGradCell (usegrad_cell[ei.Nr()]); 
                break;
              }
            }
          hofe -> SetType1 (type1);          
          hofe -> ComputeNDof();
          // hofe -> SetDiscontinuous(discontinuous);
          
          return *hofe;
        }
      case BND:
        {
          if ( discontinuous )
            return * new (lh) DummyFE<ET>; 

          if (!DefinedOn (ei))
            return * new (lh) HCurlDummyFE<ET> ();
          
          Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,BND> (ei.Nr());
          
          HCurlHighOrderFE<ET> * hofe =  new (lh) HCurlHighOrderFE<ET> ();
          hofe -> SetVertexNumbers (ngel.vertices);
          
          hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
          hofe -> SetUseGradEdge (usegrad_edge[ngel.Edges()]);
          
          
          if(ma->GetElType(ei) == ET_SEGM)
            {
              hofe -> SetOrderCell (order_edge[ngel.edges[0]]);  // old style
              FlatArray<TORDER> aoe(1, &order_edge[ngel.edges[0]]);
              hofe -> SetOrderEdge (aoe);
              if (highest_order_dc)
                hofe->SetOrderEdge (0, aoe[0]-1);
              hofe -> SetUseGradCell (usegrad_edge[ngel.edges[0]]);  // old style
            } 
          else 
            {     
              IVec<2> p = order_face[ma->GetSElFace(ei.Nr())];
              hofe -> SetOrderCell (IVec<3> (p[0],p[1],0));  // old style
              FlatArray<IVec<2> > of(1, &p);
              hofe -> SetOrderFace (of);
              
              FlatArray<bool> augf(1, &usegrad_face[ma->GetSElFace(ei.Nr())]);
              hofe -> SetUseGradFace(augf); 
              hofe -> SetUseGradCell(usegrad_face[ma->GetSElFace(ei.Nr())]);   // old style
            }
          hofe -> SetType1 (type1);              
          hofe -> ComputeNDof();
          return *hofe;
        }
      case BBND:
        {
          if (!DefinedOn (ei))
            return * new (lh) DummyFE<ET_SEGM>; 
          
          Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,BBND> (ei.Nr());
          
          HCurlHighOrderFE<ET> * hofe =  new (lh) HCurlHighOrderFE<ET> ();
          hofe -> SetVertexNumbers (ngel.vertices);
          
          hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
          hofe -> SetUseGradEdge (usegrad_edge[ngel.Edges()]);
          
          
          if(ma->GetElement(ei).GetType() == ET_SEGM)
            {
              hofe -> SetOrderCell (order_edge[ngel.edges[0]]);  // old style
              FlatArray<TORDER> aoe(1, &order_edge[ngel.edges[0]]);
              hofe -> SetOrderEdge (aoe);
              hofe -> SetUseGradCell (usegrad_edge[ngel.edges[0]]);  // old style
            } 
          else 
            {
              throw Exception("Only SEGM possible for codim 2 element of hcurlhofe space");
            }
          hofe -> SetType1 (type1);              
          hofe -> ComputeNDof();
          return *hofe;    
        }
      case BBBND: default:
        return * new (lh) DummyFE<ET_POINT>; 
      }
  }
  
  // const FiniteElement & HCurlHighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  // {
  //   Ngs_Element ngel = ma->GetElement(elnr);
  //   ELEMENT_TYPE eltype = ngel.GetType();

  //   /*
  //   if (!DefinedOn (ma->GetElIndex (elnr)))
  //     {
  //       switch (eltype)
  //         {
  //         case ET_TRIG:    return * new (lh) HCurlDummyFE<ET_TRIG> (); 
  //         case ET_QUAD:    return * new (lh) HCurlDummyFE<ET_QUAD> (); 
  //         case ET_TET:     return * new (lh) HCurlDummyFE<ET_TET> (); 
  //         case ET_PYRAMID: return * new (lh) HCurlDummyFE<ET_PYRAMID> (); 
  //         case ET_PRISM:   return * new (lh) HCurlDummyFE<ET_PRISM> (); 
  //         case ET_HEX:     return * new (lh) HCurlDummyFE<ET_HEX> (); 
  //         case ET_SEGM:    break;
  //         case ET_POINT:   break;
  //         }
  //     }
  //   */

    
  //   switch (eltype)
  //     {
  //     case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, lh);
        
  //     case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, lh);
  //     case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, lh);
        
  //     case ET_TET:     return T_GetFE<ET_TET> (elnr, lh);
  //     case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, lh);
  //     case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, lh);
  //     case ET_HEX:     return T_GetFE<ET_HEX> (elnr, lh);

  //     default:
  //       throw Exception ("illegal element in HCurlHoFeSpace::GetFE");
  //     }

  //       /*
  //   FiniteElement * fe = 0;
  //   switch (ma->GetElType(elnr))
  //     {
  //     case ET_TET:     fe = new (lh) HCurlHighOrderFE<ET_TET> (); break;
  //     case ET_PYRAMID: fe = new (lh) HCurlHighOrderFE<ET_PYRAMID> (); break;
  //     case ET_PRISM:   fe = new (lh) HCurlHighOrderFE<ET_PRISM> (); break;
  //     case ET_TRIG:    fe = new (lh) HCurlHighOrderFE<ET_TRIG> (); break;
  //     case ET_QUAD:    fe = new (lh) HCurlHighOrderFE<ET_QUAD> (); break;
  //     case ET_HEX:     fe = new (lh) HCurlHighOrderFE<ET_HEX> (); break; 

  //     default:
  //       stringstream str;
  //       str << "HCurlHighOrderFESpace " << GetClassName() 
  //           << ", undefined eltype " 
  //           << ElementTopology::GetElementName(ma->GetElType(elnr))
  //           << ", order = " << order << endl;
  //       throw Exception (str.str());
  //     }
  //       */
    
  //   /*
  //   ArrayMem<int,12> vnums;
  //   ma->GetElVertices(elnr, vnums);

  //   if(ma->GetDimension() == 2) 
  //     {	
  //       HCurlHighOrderFiniteElement<2> * hofe = 
  //         static_cast<HCurlHighOrderFiniteElement<2>*> (fe);

               
  //   	ArrayMem<int, 4> ednums, ord_edge;
  //       ArrayMem<bool, 4> ug_edge;
	
  //       ma->GetElEdges(elnr, ednums);
	
  //       ord_edge.SetSize (ednums.Size());
  //       ug_edge.SetSize (ednums.Size()); 
	
  //       for (int j = 0; j < ednums.Size(); j++)
  //         {
  //           ord_edge[j] = order_edge[ednums[j]];
  //           ug_edge[j]  = usegrad_edge[ednums[j]]; 
  //         }

  //       hofe -> SetOrderEdge (ord_edge);
  //       hofe -> SetOrderCell (order_inner[elnr]);   // old style
  //       IVec<2> p(order_inner[elnr][0], order_inner[elnr][1]);
  //       FlatArray<IVec<2> > of(1, &p);
  //       hofe -> SetOrderFace (of);

  //       hofe -> SetUseGradEdge (ug_edge); 
  //       hofe -> SetUseGradCell (usegrad_cell[elnr]);  // old style
  //       FlatArray<bool> augf(1,&usegrad_cell[elnr]);
  //       hofe -> SetUseGradFace (augf); 
  //       hofe -> ComputeNDof();

  //       hofe -> SetVertexNumbers (vnums);
  //     }   

  //   else if (ma->GetDimension() == 3) 

  //     {
  //       Ngs_Element ngel = ma->GetElement<3> (elnr);

  //       HCurlHighOrderFiniteElement<3> * hofe = 
  //         static_cast<HCurlHighOrderFiniteElement<3>*> (fe);

  //       hofe -> SetVertexNumbers (ngel.vertices);

  //       hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
  //       hofe -> SetUseGradEdge (usegrad_edge[ngel.Edges()]);

  //       hofe -> SetOrderFace (order_face[ngel.Faces()]);
  //       hofe -> SetUseGradFace (usegrad_face[ngel.Faces()]);

  //       hofe -> SetOrderCell (order_inner[elnr]);
  //       hofe -> SetUseGradCell (usegrad_cell[elnr]); 

  //       hofe -> ComputeNDof();
  //       // hofe -> SetDiscontinuous(discontinuous);
  //     }

  //   return *fe;
  //   */
  // }

 
  // const FiniteElement & HCurlHighOrderFESpace :: GetSFE (int selnr, LocalHeap & lh) const
  // {

  //   switch (ma->GetSElType(selnr))
  //     {
  //     case ET_SEGM: return T_GetSFE<ET_SEGM> (selnr, lh);
  //     case ET_TRIG: return T_GetSFE<ET_TRIG> (selnr, lh);
  //     case ET_QUAD: return T_GetSFE<ET_QUAD> (selnr, lh);

  //     default:
  //       throw Exception ("illegal element in HCurlHoFeSpace::GetSFE");
  //     }

    
  //   /*
  //   if (!DefinedOnBoundary (ma->GetSElIndex (selnr)))
  //     {
  //       switch (ma->GetSElType(selnr))
  //         {
  //         case ET_TRIG:    return * new (lh) HCurlDummyFE<ET_TRIG> ();
  //         case ET_QUAD:    return * new (lh) HCurlDummyFE<ET_QUAD> ();
  //         default: throw Exception ("not all case treated in HCurlHighOrderFESpace::GetSFE");
  //         }
  //     }

  //   FiniteElement * fe = 0;
    
  //   if ( discontinuous )
  //     {
  //       switch (ma->GetSElType(selnr))
  //         {
  //         case ET_SEGM: fe = new (lh) DummyFE<ET_SEGM>; break; 
  //         case ET_TRIG: fe = new (lh) DummyFE<ET_TRIG>; break;
  //         case ET_QUAD: fe = new (lh) DummyFE<ET_QUAD>; break;
  //         default:
  //           fe = 0;
  //         }
  //       return *fe;
  //     }

  //   switch (ma->GetSElType(selnr))
  //     {
  //     case ET_SEGM: fe = new (lh) HCurlHighOrderFE<ET_SEGM> (); break; 
  //     case ET_TRIG: fe = new (lh) HCurlHighOrderFE<ET_TRIG> (); break;
  //     case ET_QUAD: fe = new (lh) HCurlHighOrderFE<ET_QUAD> (); break;
  //     default:
  //       fe = 0;
  //     }

  //   if (!fe)
  //     {
  //       stringstream str;
  //       str << "HCurlHighOrderFESpace " << GetClassName() 
  //           << ", undefined eltype " 
  //           << ElementTopology::GetElementName(ma->GetSElType(selnr))
  //           << ", order = " << order << endl;
  //       throw Exception (str.str());
  //     }

  //   ArrayMem<int, 8> vnums;
  //   ArrayMem<int, 4> ednums, ord_edge;
  //   ArrayMem<bool, 4> ug_edge; 
       
  //   ma->GetSElVertices(selnr, vnums);

  //   if(ma->GetSElType(selnr) == ET_SEGM)
  //     {
  //       HCurlHighOrderFiniteElement<1> * hofe =
  //         dynamic_cast<HCurlHighOrderFiniteElement<1>*> (fe);

  //       hofe -> SetVertexNumbers (vnums);
  //       ma->GetSElEdges(selnr, ednums);
  //       hofe -> SetOrderCell (order_edge[ednums[0]]);  // old style
  //       FlatArray<int> aoe(1, &order_edge[ednums[0]]);
  //       // hofe -> SetOrderEdge (order_edge[ednums[0]]);
  //       hofe -> SetOrderEdge (aoe);
  //       hofe -> SetUseGradCell (usegrad_edge[ednums[0]]);  // old style
  //       FlatArray<bool> auge(1, &usegrad_edge[ednums[0]]);
  //       hofe -> SetUseGradEdge (auge);
  //       hofe -> ComputeNDof();
  //     } 
  //   else 
  //     {     
  //       HCurlHighOrderFiniteElement<2> * hofe =
  //         dynamic_cast<HCurlHighOrderFiniteElement<2>*> (fe);
      
  //       ma->GetSElEdges(selnr, ednums);
	
  //       ord_edge.SetSize (ednums.Size());
  //       ug_edge.SetSize (ednums.Size()); 
	
  //       for (int j = 0; j < ednums.Size(); j++)
  //         {
  //           ord_edge[j] = order_edge[ednums[j]]; 
  //           ug_edge[j] = usegrad_edge[ednums[j]];
  //         }
	
  //       hofe -> SetVertexNumbers (vnums);
  //       hofe -> SetOrderEdge (ord_edge);

  //       IVec<2> p = order_face[ma->GetSElFace(selnr)];
  //       hofe -> SetOrderCell (IVec<3> (p[0],p[1],0));  // old style
  //       FlatArray<IVec<2> > of(1, &p);
  //       hofe -> SetOrderFace (of);

  //       hofe -> SetUseGradEdge(ug_edge); 
  //       FlatArray<bool> augf(1, &usegrad_face[ma->GetSElFace(selnr)]);
  //       hofe -> SetUseGradFace(augf); 
  //       hofe -> SetUseGradCell(usegrad_face[ma->GetSElFace(selnr)]);   // old style

  //       hofe -> ComputeNDof();
  //     }
    
  //   return *fe;
  //   */
  // }

  // template <ELEMENT_TYPE ET>
  // const FiniteElement & HCurlHighOrderFESpace :: T_GetSFE (int selnr, LocalHeap & lh) const
  // {
  // }


  // const FiniteElement & HCurlHighOrderFESpace :: GetCD2FE(int cd2elnr, LocalHeap & lh) const
  // {
  //   switch (ma->GetElement(ElementId(BBND,cd2elnr)).GetType())
  //     {
  //     case ET_SEGM: return T_GetCD2FE<ET_SEGM> (cd2elnr, lh);

  //     default:
  //       throw Exception ("illegal element in HCurlHoFeSpace::GetCD2FE");
  //     }
  // }

  // template <ELEMENT_TYPE ET>
  // const FiniteElement & HCurlHighOrderFESpace :: T_GetCD2FE(int cd2elnr, LocalHeap & lh) const
  // {
    
  // }
    


  size_t HCurlHighOrderFESpace :: GetNDof () const throw()
  {
    return ndof;
  }


  void HCurlHighOrderFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    if (!DefinedOn (ei)) return;

    // ordering of shape functions
    // (1*e1),.. (1*e_ne)  
    // (p_t)*tang_e1, ... p_t*tang_ne
    // p_n * el1, p_i * el1, ... , p_n * ne_n , p_i *el_ne 

    Ngs_Element ngel = ma->GetElement (ei);
     
    //Nedelec0
    if ( !discontinuous )
      for (auto edge : ngel.Edges())
        dnums.Append (edge);
    
    // new style
    // if (order_policy == VARIABLE_ORDER && ma->GetDimension() >= 2)
    if(false)
      {
        for (auto edge : ngel.Edges())
          dnums += GetEdgeDofs(edge);
        for (auto face : ngel.Faces())
          dnums += GetFaceDofs(face);

        if (ma->GetDimension() == 3 && ei.IsVolume())
          {
            dnums += GetElementDofs(ngel.Nr());
          }
        return;
      }


    IntRange eldofs;
    if (ei.IsVolume())
      eldofs = GetElementDofs (ei.Nr());
    
    for (auto edge : ngel.Edges())
      {
        dnums += GetEdgeDofs(edge);
        if (ei.IsVolume() && highest_order_dc)
          {
            dnums += eldofs.First();
            eldofs.First()++;
          }
      }
    
    // faces 
    if (ma->GetDimension() == 3)
      for (auto face : ngel.Faces())
        dnums += GetFaceDofs(face);

    if (ei.IsVolume())
      dnums += eldofs;
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
    if (discontinuous) return;

    dnums.Append(ednr);
    dnums += GetEdgeDofs(ednr);
  }

  void HCurlHighOrderFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    if (order_policy == VARIABLE_ORDER)
      {
        dnums = GetFaceDofs(fanr);
        return;
      }
   
    if (discontinuous) 
      {
        dnums.SetSize(0);
        return;
      }

    dnums = GetFaceDofs (fanr);
  }

  void HCurlHighOrderFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    if (order_policy == VARIABLE_ORDER && ma->GetDimension() == 2)
      {
        dnums.SetSize0();
        return;
      }
    dnums = GetElementDofs (elnr);
  }  


  
  shared_ptr<Table<int>> HCurlHighOrderFESpace :: 
  CreateSmoothingBlocks (const Flags & precflags) const
  {
    if (precflags.StringFlagDefined("blocktype") || precflags.StringListFlagDefined("blocktype"))
      return FESpace::CreateSmoothingBlocks(precflags);

    
    size_t nv = ma->GetNV();
    size_t ne = ma->GetNE();
    // int nse = ma->GetNSE();

    size_t ned = ma->GetNEdges();
    size_t nfa = (ma->GetDimension() == 2) ? 0 : ma->GetNFaces();


    bool eliminate_internal = precflags.GetDefineFlag("eliminate_internal");
    int i, j, k, first; 
    int ni = ne;
    if (eliminate_internal) ni = 0; 
    
    int SmoothingType = int(precflags.GetNumFlag("blocktype",2)); 
    bool excl_grads = precflags.GetDefineFlag("exclude_grads"); 
    cout << IM(5) << " EXCLUDE GRADS " << excl_grads << endl;
    
    Array<int> vnums; 
    // Array<int> orient; 
    Array<int> fanums;
    
    // int augv = augmented; 

        
    if(nfa == 0 && SmoothingType == 1) 
      SmoothingType = 4; 
    if(nfa == 0 && (SmoothingType == 2 || SmoothingType ==3)) 
      SmoothingType = 5;




    cout << "hcurl smoothingblocks, SmoothingType = " << SmoothingType << endl;

    if (precflags.GetDefineFlag("subassembled"))
      {
	
	TableCreator<int> creator;
	for ( ; !creator.Done(); creator++)
	  {

	    if (creator.GetMode() == 1)
	      cout << IM(5) << "High order AFW blocks " << endl;
		
	    for (int i = 0; i < ned; i++)
	      if (!IsDirichletEdge(i) && fine_edge[i])
		{
		  Ng_Node<1> edge = ma->GetNode<1> (i);
		  for (int k = 0; k < 2; k++)
		    {
		      creator.Add (edge.vertices[k], i);
		      creator.Add (edge.vertices[k], GetEdgeDofs(i));
		    }
		}
	  }

	
	// return shared_ptr<Table<int>> (creator.GetTable());
        return make_shared<Table<int>> (creator.MoveTable());
      }


    cout << IM(5) << "SmoothingType " << SmoothingType << endl;
    cout << IM(5) << " Use H(Curl)-Block smoothing " ;



    FilteredTableCreator creator(GetFreeDofs().get());

    for ( ; !creator.Done(); creator++)
      {
        switch(SmoothingType)           
	  {
          case 1: // Clustering + AFW(loE)   (former 10)
            {
              if (creator.GetMode() == 1)
             	cout << IM(5) << " Line-Clustering and  AFW(loE)  " << endl;

              for (size_t i = 0; i < ned; i++)
                if (fine_edge[i])
                  {
                    IVec<2> vts = ma->GetEdgePNums (i);
                    int pn1 = vts[0], pn2 = vts[1];
                    
                    pn1 = ma->GetClusterRepVertex(pn1);
                    pn2 = ma->GetClusterRepVertex(pn2);
                    // table[pn1][cnt[pn1]++] = i;
                    creator.Add (pn1, i);
                    if (pn1 != pn2)
                      // table[pn2][cnt[pn2]++] = i;
                      creator.Add (pn2, i);
                    
                    int nr = ma->GetClusterRepEdge (i);
                    if(!excl_grads)
                      for (j = first_edge_dof[i]; j < first_edge_dof[i+1]; j++)
                        // table[nr][cnt[nr]++] = j;
                        creator.Add (nr, j);
                  }

              for (size_t i = 0; i < nfa; i++)
                if(fine_face[i])
                  {
                    int nr = ma->GetClusterRepFace (i);
                    int first = first_face_dof[i] + excl_grads*face_ngrad[i];
                    for (j = first; j < first_face_dof[i+1]; j++)
                      // table[nr][cnt[nr]++] = j;
                      creator.Add(nr,j);
                  }
              
              for (size_t i = 0; i < ni; i++)
                {
                  int nr = ma->GetClusterRepElement (i);
                  int first = first_inner_dof[i] + excl_grads*cell_ngrad[i];
                  for (j = first; j < first_inner_dof[i+1]; j++)
                    // table[nr][cnt[nr]++] = j;
                    creator.Add(nr, j);
                }
               
              break;
            }

          case 2: // Clustering with AFW(hoEdges) 
            {
              if (creator.GetMode() == 1)
                cout << IM(5) << " Line-Clustering and AFW(hoE)  " << endl;

              for (size_t i=0;i<ned;i++) 
                if(fine_edge[i] && !IsDirichletEdge(i))
                  {
                    IVec<2> vts = ma->GetEdgePNums (i);
                    int pn1 = vts[0], pn2 = vts[1];
                    
                    pn1 = ma->GetClusterRepVertex(pn1);
                    pn2 = ma->GetClusterRepVertex(pn2);
                    int ecl = ma->GetClusterRepEdge(i);
                    // table[pn1][cnt[pn1]++] = i;
                    creator.Add (pn1, i);
                    if(!excl_grads)
                      for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
                        // table[pn1][cnt[pn1]++] = j;
                        creator.Add (pn1, j);
                    if(pn1!=pn2)
                      {
                        // table[pn2][cnt[pn2]++] = i; 
                        // table[ecl][cnt[ecl]++] = i;
                        creator.Add (pn2, i);
                        creator.Add (ecl, i);
                        
                        // ho-edges
                        if(!excl_grads)
                          for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
                            {		    
                              // table[pn2][cnt[pn2]++] = j;
                              // table[ecl][cnt[ecl]++] = j;
                              creator.Add (pn2, j);
                              creator.Add (ecl, j);
                            }
                      }
                  }
              
              // Stack of horizontal edges, quadfaces and prism faces 
              //  (+ inner) 
              for(size_t i=0;i<nfa;i++) 
                if(fine_face[i] && !IsDirichletFace(i)) 
                  { 
                    int fcl = ma->GetClusterRepFace(i); 
                    int first = first_face_dof[i] + excl_grads*face_ngrad[i];
                    //	if(fcl  >= nv + ned) // quad-faces (prism) 
                    for(j=first; j< first_face_dof[i+1]; j++) 
                      // table[fcl][cnt[fcl]++] = j;
                      creator.Add (fcl, j);
                  }
              
              for(size_t i=0;i<ni;i++) 
                { 
                  // Stack: plane-face + inner
                  int ccl = ma->GetClusterRepElement(i); 
                  int first = first_inner_dof[i] + excl_grads*cell_ngrad[i];
                  for(j=first;j<first_inner_dof[i+1];j++) 
                    // table[ccl][cnt[ccl]++] = j;
                    creator.Add (ccl, j);
                }
              
              break;
            }

          case 3:  // Stack-AFW(hoEdges) 
            {
              if (creator.GetMode() == 1)              
                cout << IM(5) << " Line-Clustering: AFW(hoE) - (horizE)F - (trigF)I " << endl;

              	  // Stack-AFW(hoEdges) 
              for(size_t i=0;i<ned;i++) 
                if(fine_edge[i])
                  {
                    IVec<2> vts = ma->GetEdgePNums (i);
                    int pn1 = vts[0], pn2 = vts[1];
                    
                    pn1 = ma->GetClusterRepVertex(pn1);
                    pn2 = ma->GetClusterRepVertex(pn2);
                    int ecl = ma->GetClusterRepEdge(i);
                    // table[pn1][cnt[pn1]++] = i;
                    creator.Add (pn1, i);
                    if(!excl_grads)
                      for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
                        // table[pn1][cnt[pn1]++] = j;
                        creator.Add (pn1, j);
                    
                    if(pn1!=pn2)
                      {
                        // table[pn2][cnt[pn2]++] = i; 
                        // table[ecl][cnt[ecl]++] = i;
                        creator.Add (pn2, i);
                        creator.Add (ecl, i);
                        // ho-edges
                        if(!excl_grads)
                          for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
                            {		    
                              // table[pn2][cnt[pn2]++] = j;
                              // table[ecl][cnt[ecl]++] = j;
                              creator.Add (pn2, j);
                              creator.Add (ecl, j);
                            }
                      }
                  }
	  
              // Stack of horizontal edges, quadfaces and prism faces 
              //  (+ inner) 
              for (size_t i=0;i<nfa;i++) 
                if(fine_face[i]) 
                  { 
                    int fcl = ma->GetClusterRepFace(i); 
                    int first = first_face_dof[i] + excl_grads*face_ngrad[i];
                    
                    for(j=first; j< first_face_dof[i+1]; j++) 
                      // table[fcl][cnt[fcl]++] = j;
                      creator.Add (fcl, j);
                  }
              
              for(size_t i=0;i<ni;i++) 
                { 
                  // Stack: plane-face + inner
                  int ccl = ma->GetClusterRepElement(i);
                  int first = first_inner_dof[i] + excl_grads*cell_ngrad[i]; 
                  for(j=first;j<first_inner_dof[i+1];j++) 
                    // table[ccl][cnt[ccl]++] = j;
                    creator.Add (ccl, j);
                  
                  // inner to quad-face (horizontal edge stack) 
                  // each face different stack 
                  
                  auto fanums = ma->GetElFaces(ElementId(VOL,i)); 
                  
                  for(j=0;j<fanums.Size();j++) 
                    { 
                      int fcl = ma->GetClusterRepFace(fanums[j]);
                      if(fcl != ccl)
                        for(k=first;k<first_inner_dof[i+1];k++) 
                          // table[fcl][cnt[fcl]++] = k;
                          creator.Add (fcl, k);
                    } 
                  
                }  
              
              for(size_t i=0;i<nfa;i++) 
                if(fine_face[i]) 
                  { 
                    vnums = ma->GetFacePNums(i); 
                    if(vnums.Size()==4) continue; 
                    int fcl = ma->GetClusterRepFace(i); 
                    auto ednums = ma->GetFaceEdges(i); 
                    
                    for(j=0;j<ednums.Size();j++)
                      {
                        int ecl = ma->GetClusterRepEdge(ednums[j]); 
                        if(ecl==fcl) continue; 
                        int first = first_face_dof[i] + excl_grads*face_ngrad[i];
                        for(k=first; k<first_face_dof[i+1]; k++)
                          // table[ecl][cnt[ecl]++] = k;
                          creator.Add (ecl, k);
                      }
                  }

              
              break;
            }
            
          case 4: // AFW(loE) + E  + F + I (noLineClusters) 
            {
              if (creator.GetMode() == 1)              
              	cout << IM(5) << " AFW(loE) + E + F + I (noClusters)" << endl;

              for (size_t i = 0; i < ned; i++)
                if(fine_edge[i])
                  {
                    IVec<2> vts = ma->GetEdgePNums (i);
                    int pn1 = vts[0], pn2 = vts[1];
                    // table[pn1][cnt[pn1]++] = i;
                    // table[pn2][cnt[pn2]++] = i;
                    creator.Add (pn1, i);
                    creator.Add (pn2, i);
                  }
              int nnv = nv; 
              
              
              if(augmented==1)
                { 
                  nnv += nv; 
                  for(i=0;i<nv;i++)
                    // table[nv+i][cnt[nv+i]++] = ned+i;
                    creator.Add (nv+i, ned+i);
                } 
              
              if(!excl_grads)
                for (i =0; i<ned ; i++)
                  for (j = first_edge_dof[i]; j < first_edge_dof[i+1]; j++)
                    // table[nnv+i][cnt[nnv+i]++] = j;
                    creator.Add(nnv+i, j);
              
              for (size_t i = 0; i < nfa; i++)
                { 
                  int first = first_face_dof[i] + excl_grads*face_ngrad[i]; 
                  for (j = first; j < first_face_dof[i+1]; j++)
                    // table[nnv+ned+i][cnt[nnv+ned+i]++] = j;
                    creator.Add (nnv+ned+i, j);
                }
	      
              for (size_t i = 0; i < ni; i++)
                { 
                  // int first = first_inner_dof[i] + excl_grads*cell_ngrad[i];
                  for (j = first_inner_dof[i]; j < first_inner_dof[i+1]; j++)
                    // table[nnv+ned+nfa+i][cnt[nnv+ned+nfa+i]++] = j;
                    creator.Add (nnv+ned+nfa+i, j);
                }
              
              break;
            }
            
          case 5: // AFW(hoE) + E  + F + I (noLineClusters) 
            {
              if (creator.GetMode() == 1)              
                cout << IM(5) << " AFW(hoE) + E + F + I (noClusters)" << endl;

              for (size_t i = 0; i < ned; i++)
                if(fine_edge[i])
                  {
                    IVec<2> vts = ma->GetEdgePNums (i);
                    int pn1 = vts[0], pn2 = vts[1];
                    
                    // table[pn1][cnt[pn1]++]  = i; 
                    // table[pn2][cnt[pn2]++]  = i;
                    creator.Add (pn1, i);
                    creator.Add (pn2, i);
                    
                    if(!excl_grads)
                      for(j=first_edge_dof[i];j<first_edge_dof[i+1];j++) 
                        {
                          // table[pn1][cnt[pn1]++]  = j; 
                          // table[pn2][cnt[pn2]++]  = j;
                          creator.Add (pn1, j);
                          creator.Add (pn2, j);
                        }
                  }
              
              int nnv = nv; 
              /* if(augmented==1)
                 { 
                 nnv += nv; 
                 for(i=0;i<nv;i++)
                 table[i][cnt[i]++] = ned+i;
                 } */ 
              for (size_t i = 0; i < nfa; i++)
                { 
                  int first = first_face_dof[i] + excl_grads*face_ngrad[i]; 
                  for (j = first; j < first_face_dof[i+1]; j++)
                    // table[nnv+ned+i][cnt[nnv+ned+i]++] = j;
                    creator.Add (nnv+ned+i, j);
                }
              
              for (size_t i = 0; i < ni; i++)
                { 
                  // int first = first_inner_dof[i] + excl_grads*cell_ngrad[i];
                  for (j = first_inner_dof[i]; j < first_inner_dof[i+1]; j++)
                    // table[nnv+ned+nfa+i][cnt[nnv+ned+nfa+i]++] = j;
                    creator.Add (nnv+ned+nfa+i, j);
                }
              
              break;
            }

          case 6: // der hier is nur fuer testzwecke !!!
            {
              if (creator.GetMode() == 1)
              	cout << IM(5) << "Jacobi (Diag)" << endl;

              int ii=0; 
              
              for (size_t i = 0; i < ned; i++)
                if(fine_edge[i])
                  //table[ii++][0] = i;
                  creator.Add (ii++, i);

              for (size_t i =0; i<ned ; i++)
                if(fine_edge[i])
                  for (j = first_edge_dof[i]; j < first_edge_dof[i+1]; j++)
                    // table[ii++][0] = j;
                    creator.Add (ii++, j);
	  
              for (size_t i = 0; i < nfa; i++)
                { 
                  if(fine_face[i]) 
                    {
                      int first = first_face_dof[i]; 
                      for (j = first; j < first_face_dof[i+1]; j++)
                        // table[ii++][0] = j;
                        creator.Add (ii++, j);
                    }
                }
	      
              for (size_t i = 0; i < ni; i++)
                { 
                  for (j = first_inner_dof[i]; j < first_inner_dof[i+1]; j++)
                    // table[ii++][0] = j;
                    creator.Add (ii++,j);
                }

              
              break;
            }
            
          case 7:
            {
              // AFW(ho-E) - Edge-Blocks ???
              // Stack-AFW(hoEdges) ???
              
              if (creator.GetMode() == 1)
              	cout << IM(5) << "EDGE-blocks (noClusters)" << endl;
              
              for (size_t i=0;i<ned;i++) 
                if(fine_edge[i])
                  {
                    // int pn1, pn2;
                    // ma->GetEdgePNums (i,pn1,pn2);
                    IVec<2> vts = ma->GetEdgePNums (i);
                    int pn1 = vts[0], pn2 = vts[1];
                    
                    pn1 = ma->GetClusterRepVertex(pn1);
                    pn2 = ma->GetClusterRepVertex(pn2);
                    //	int ecl = ma->GetClusterRepEdge(i);
                    // table[pn1][cnt[pn1]++] = i;
                    creator.Add (pn1, i);
                    
                    if(!excl_grads)
                      for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
                        // table[pn1][cnt[pn1]++] = j;
                        creator.Add (pn1, j);
                    if(pn1!=pn2)
                      {
                        // table[pn2][cnt[pn2]++] = i;
                        creator.Add (pn2, i);
                        // //  table[ecl][cnt[ecl]++] = i;  was already commented
                        // ho-edges
                        if(!excl_grads)
                          for(j=first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
                            {		    
                              // table[pn2][cnt[pn2]++] = j;
                              creator.Add (pn2, j);
                              
                              // //  table[ecl][cnt[ecl]++] = j; was already commendted
                            }
                      }
                  }
              
              for (size_t i = 0; i < ned; i++)
                {
                  first = first_edge_dof[i];
                  int ndof = first_edge_dof[i+1]-first;
                  for (j = 0; j < ndof; j++)
                    // table[nv+i][cnt[nv+i]++] = first+j;
                    creator.Add (nv+i, first+j);
                }
              
              for (size_t i = 0; i < nfa; i++)
                {
                  first = first_face_dof[i];
                  int ndof = first_face_dof[i+1]-first;
                  auto f2ed = ma->GetFaceEdges (i);
                  for (k = 0; k < f2ed.Size(); k++)
                    for (j = 0; j < ndof; j++)
                      // table[nv+f2ed[k]][cnt[nv+f2ed[k]]++] = first+j;
                      creator.Add (nv+f2ed[k], first+j);
                }
              
              for (size_t i = 0; i < ni; i++)
                {
                  auto enums = ma->GetElEdges (ElementId(VOL,i));
                  first = first_inner_dof[i];
                  int ndof = first_inner_dof[i+1] - first_inner_dof[i];
                  for (j = 0; j < enums.Size(); j++)
                    for (k = 0; k < ndof; k++)
                      // table[nv+enums[j]][cnt[nv+enums[j]]++] = first + k;
                      creator.Add (nv+enums[j], first+k);
                } 
              
              
              break;
            }
            
          case 21:
            {
              if (creator.GetMode() == 1)
                cout << IM(5) << "wave equation blocks" << endl;
              
              
              int ds_order = int(precflags.GetNumFlag ("ds_order", 0));
              if (ds_order < 0) ds_order = 0;	  
              
              for (size_t i =0; i<ned ; i++)
                for (j = first_edge_dof[i]+ds_order; j < first_edge_dof[i+1]; j++)
                  // table[i][cnt[i]++] = j;
                  creator.Add(i,j);
              
              for (size_t i = 0; i < nfa; i++)
                {
                  // cnt[ned+i] = first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i];
                  
                  int first = first_face_dof[i];
                  int p = order_face[i][0];
                  
                  int ii = first;
                  for (int j = 0; j <= p-2; j++)
                    for (int k = 0; k <= p-2-j; k++, ii++)
                      if (j+k+2 > ds_order)
                        // table[ned+i][cnt[ned+i]++] = ii;
                        creator.Add (ned+i, ii);
                  // cnt[ned+i]++;
                  //clusters[first+ii] = 1;
	      
                  // other combination
                  for (int j = 0; j <= p-2; j++)
                    for (int k = 0; k <= p-2-j; k++, ii++)
                      if (j+k+2 > ds_order)
                        // table[ned+i][cnt[ned+i]++] = ii;
                        creator.Add (ned+i, ii);
                  // cnt[ned+i]++;
                  // clusters[first+ii] = 1;
                  
                  // type 3
                  for (int j = 0; j <= p-2; j++, ii++)
                    if (j+2 > ds_order)
                      // table[ned+i][cnt[ned+i]++] = ii;
                      creator.Add(ned+i, ii);
                  // cnt[ned+i]++;
                  // clusters[first+ii] = 1;
                }

              break;
            }

          default: 
            throw Exception("HCurlHoFeSpace:: CreateSmoothingBlocks chosen Blocktype not valid \n Choose blocktype 1-6 "); 
          }
      }
    
    return make_shared<Table<int>> (creator.MoveTable());
 
    
#ifdef OLD
    int ncnt; 

    switch(SmoothingType) 
      {
      case 4: // former 11 
	cout << IM(5) << " AFW(loE) + E + F + I (noClusters)" << endl;
	ncnt = nv  + ned + nfa + ni + augmented*nv;
	break;  
      case 5: //	
	cout << IM(5) << " AFW(hoE) + E + F + I (noClusters)" << endl;
	ncnt = nv  + ned + nfa + ni;
	break; 	
      case 1: // former 10 
	cout << IM(5) << " Line-Clustering and  AFW(loE)  " << endl;
	ncnt = nv + ned + nfa + ni; 
	break;
      case 2: // former 26 
	cout << IM(5) << " Line-Clustering and AFW(hoE)  " << endl;
	ncnt = nv + ned + nfa + ni;
	break;
      case 3: // former 23 (equal to HCurl-Old-SmoothingBlocks)  
	cout << IM(5) << " Line-Clustering: AFW(hoE) - (horizE)F - (trigF)I " << endl;
	// horiz-quadface-stack(horiz-edge,top-bot-trig-Faces, inner)" << endl; 
	ncnt = nv + ned + nfa +ni; 
	break;
      case 6: // der hier is nur fuer testzwecke !!!  
	cout << IM(5) << "Jacobi (Diag)" << endl;
	ncnt = nv  + ned + nfa + ni;
	break;  
      case 7: 
	cout << IM(5) << "EDGE-blocks (noClusters)" << endl;
	ncnt = nv  + ned;
	break; 
      case 21:
	cout << IM(5) << "wave equation blocks" << endl;
	ncnt = ned + nfa;
	break;
      default: 
	throw Exception("HCurlHoFeSpace:: CreateSmoothingBlocks chosen Blocktype not valid \n Choose blocktype 1-6 "); 
	ncnt =0; 
	return 0; 
      } 

        
    Array<int> cnt(ncnt); 
    cnt = 0;
    // ii=0; 

    switch (SmoothingType)
      {
      case 4: // AFW(loE) + E + F + I (noCluster)
	{
	  for (i = 0; i < ned; i++)
	    if(fine_edge[i])
	      {		
		auto pts = ma->GetEdgePNums (i);
		cnt[pts[0]]++;
		cnt[pts[1]]++;
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
	  for (size_t i = 0; i < ned; i++)
	    if(fine_edge[i])
	      {
		auto pts = ma->GetEdgePNums (i);
		cnt[pts[0]]+= 1 + (first_edge_dof[i+1]-first_edge_dof[i])*(1-excl_grads);
		cnt[pts[1]]+= 1 + (first_edge_dof[i+1]-first_edge_dof[i])*(1-excl_grads);
	      }
	  int nnv = nv; 
	  /* 
	     if(augmented==1 && wgrads)
	     { 
	     nnv += nv; 
	     for (i=0;i<nv;i++)
	     cnt[nv+i] = 1; 
	     } */ 
	  
	  for (size_t i = 0; i < nfa; i++)
	    cnt[nnv+ned+i] = first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i];
	  for (size_t i = 0; i < ni; i++)
	    cnt[nnv+ned+nfa+i] = first_inner_dof[i+1] - first_inner_dof[i] - excl_grads*cell_ngrad[i];
	  break;
	}

      case 1: // Clustering with AFW(loE) 
	{
	  
	  for (i = 0; i < ned; i++)
	    if(fine_edge[i])
	      {
		auto pts = ma->GetEdgePNums (i);
                int pn1 = pts[0], pn2 = pts[1];
		pn1 = ma->GetClusterRepVertex(pn1);
	  
		pn2 = ma->GetClusterRepVertex(pn2);
		cnt[pn1]++;
		if (pn1 != pn2)
		  cnt[pn2]++;
		// ev. noch lo-edge zu horiz.edges 
		cnt[ma->GetClusterRepEdge(i)] +=
		  (first_edge_dof[i+1] - first_edge_dof[i])*(1-excl_grads);
	      }
	 
	  for (i = 0; i < nfa; i++)
	    if(fine_face[i])
	      cnt[ma->GetClusterRepFace(i)] += 
		first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i];
	 
	  for (i = 0; i < ni; i++)
	    cnt[ma->GetClusterRepElement(i)] += 
	      first_inner_dof[i+1] - first_inner_dof[i] - excl_grads*cell_ngrad[i];
	  
	  break;
	}
      case 2: // Clustering with AFW(hoE)
	{
	  // Stack-AFW(hoEdges) and Stack of horiz edges 
	  for(i=0;i<ned;i++) 
	    if(fine_edge[i] && !IsDirichletEdge(i))
	      {
		auto pts = ma->GetEdgePNums (i);
                int pn1 = pts[0], pn2 = pts[1];
		pn1 = ma->GetClusterRepVertex(pn1);
		pn2 = ma->GetClusterRepVertex(pn2); 
		int nde = 1+ (1-excl_grads)*(first_edge_dof[i+1] - first_edge_dof[i]); 
		cnt[pn1] += nde; 
		if (pn1 != pn2)
		  { 
		    cnt[pn2]+= nde;
		    cnt[ma->GetClusterRepEdge(i)] +=nde;
		  }
	      }
	  
	  // Stack of horizontal edges: + quadfaces and prism faces 
	  //  (+ inner) 
	  for(i=0;i<nfa;i++) 
	    if(fine_face[i] && !IsDirichletFace(i)) 
	      cnt[ma->GetClusterRepFace(i)] += first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i]; 
	     	 
	  for(i=0;i<ni;i++) 
	    { 
	      int ccl = ma->GetClusterRepElement(i);
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
		auto pts = ma->GetEdgePNums (i);
                int pn1 = pts[0], pn2 = pts[1];

		pn1 = ma->GetClusterRepVertex(pn1);
		pn2 = ma->GetClusterRepVertex(pn2); 
		int nde = 1+ (1-excl_grads)*(first_edge_dof[i+1] - first_edge_dof[i]); 
		cnt[pn1] += nde; 
		if (pn1 != pn2)
		  { 
		    cnt[pn2]+= nde;
		    cnt[ma->GetClusterRepEdge(i)] +=nde;
		  }
	      }
	  
	  // Stack of horizontal edges: + quadfaces and prism faces 
	  //  (+ inner) 
	  for(i=0;i<nfa;i++) 
	    if(fine_face[i]) 
	      cnt[ma->GetClusterRepFace(i)] += first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i]; 
	  
	  for(i=0;i<ni;i++) 
	    { 
	      int ccl = ma->GetClusterRepElement(i);
	      int ndi = first_inner_dof[i+1] - first_inner_dof[i] - excl_grads*cell_ngrad[i]; 
	      // Stack: plane-face + inner 
	      cnt[ccl] += ndi; 
	      
	      
	      
	      // inner to quad-face (horizontal edge stack) 
	      // each face different stack 
	      
	      
	      auto fanums = ma->GetElFaces(ElementId(VOL,i)); 
	      for(j=0;j<fanums.Size();j++) 
		{ 
		  int fcl = ma->GetClusterRepFace(fanums[j]);
		  if(fcl != ccl) 
		    cnt[fcl] += ndi; 
		}  
	      // if nocluster -> FI + I 
	    }  
	  
	  for(i=0;i<nfa;i++) // Trig-Faces to horiz-edge stack 
	    if(fine_face[i]) 
	      { 
		auto vnums = ma->GetFacePNums(i); 
		if(vnums.Size()==4) continue; 
		int fcl = ma->GetClusterRepFace(i); 
		auto ednums = ma->GetFaceEdges(i); 
		
		int nd = first_face_dof[i+1] - first_face_dof[i] - excl_grads*face_ngrad[i]; 
		for(j=0;j<ednums.Size();j++)
		  {
		    int ecl = ma->GetClusterRepEdge(ednums[j]); 
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
		IVec<2> vts = ma->GetEdgePNums (i);
                int pn1 = vts[0], pn2 = vts[1];
		pn1 = ma->GetClusterRepVertex(pn1);
		pn2 = ma->GetClusterRepVertex(pn2); 
		int nde = 1+ (1-excl_grads)*(first_edge_dof[i+1] - first_edge_dof[i]); 
		cnt[pn1] += nde; 
		if (pn1 != pn2)
		  { 
		    cnt[pn2]+= nde;
		    cnt[ma->GetClusterRepEdge(i)] +=nde;
		  }
	      }
	  for (i = 0; i < ned; i++)
	    cnt[nv+i]= first_edge_dof[i+1]-first_edge_dof[i];
	  for (i = 0; i < nfa; i++)
	    {
	      auto f2ed = ma->GetFaceEdges (i);
	      for (j = 0; j < f2ed.Size(); j++)
		cnt[nv+f2ed[j]] +=  first_face_dof[i+1]-first_face_dof[i];
	    }
	  for (i=0; i< ni; i++) 
	    {
	      auto enums = ma->GetElEdges (ElementId(VOL,i));
	      int ndi = first_inner_dof[i+1] - first_inner_dof[i];
	      for (j = 0; j < enums.Size(); j++)
		cnt[nv+enums[j]] += ndi;
	    }

	  break;
      
	}
     case 21: // wave equation
	{
	  int ds_order = int(precflags.GetNumFlag ("ds_order", 0));
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
    
   
        
    Table<int> table(cnt); 
    // ii = 0; 
    cnt = 0;
    switch(SmoothingType) 
      {
      case 4: // AFW(loE) + E  + F + I (noLineClusters) 
	{
	  cnt = 0;
	  for (i = 0; i < ned; i++)
	    if(fine_edge[i])
	      {
		IVec<2> vts = ma->GetEdgePNums (i);
                int pn1 = vts[0], pn2 = vts[1];
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
		IVec<2> vts = ma->GetEdgePNums (i);
                int pn1 = vts[0], pn2 = vts[1];
                
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
		IVec<2> vts = ma->GetEdgePNums (i);
                int pn1 = vts[0], pn2 = vts[1];
                
		pn1 = ma->GetClusterRepVertex(pn1);
		pn2 = ma->GetClusterRepVertex(pn2);
		table[pn1][cnt[pn1]++] = i;
		if (pn1 != pn2)
		  table[pn2][cnt[pn2]++] = i;

		int nr = ma->GetClusterRepEdge (i);
		if(!excl_grads)
		  for (j = first_edge_dof[i]; j < first_edge_dof[i+1]; j++)
		    table[nr][cnt[nr]++] = j;
	      }
	  for (i = 0; i < nfa; i++)
	    if(fine_face[i])
	      {
		int nr = ma->GetClusterRepFace (i);
		int first = first_face_dof[i] + excl_grads*face_ngrad[i];
		for (j = first; j < first_face_dof[i+1]; j++)
		  table[nr][cnt[nr]++] = j;
	      }
	  for (i = 0; i < ni; i++)
	    {
	      int nr = ma->GetClusterRepElement (i);
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
		IVec<2> vts = ma->GetEdgePNums (i);
                int pn1 = vts[0], pn2 = vts[1];
                
		pn1 = ma->GetClusterRepVertex(pn1);
		pn2 = ma->GetClusterRepVertex(pn2);
		int ecl = ma->GetClusterRepEdge(i);
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
		int fcl = ma->GetClusterRepFace(i); 
		int first = first_face_dof[i] + excl_grads*face_ngrad[i];
		//	if(fcl  >= nv + ned) // quad-faces (prism) 
		for(j=first; j< first_face_dof[i+1]; j++) 
		  table[fcl][cnt[fcl]++] = j; 
	      }
	  
	  for(i=0;i<ni;i++) 
	    { 
	      // Stack: plane-face + inner
	      int ccl = ma->GetClusterRepElement(i); 
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
		IVec<2> vts = ma->GetEdgePNums (i);
                int pn1 = vts[0], pn2 = vts[1];
                
		pn1 = ma->GetClusterRepVertex(pn1);
		pn2 = ma->GetClusterRepVertex(pn2);
		int ecl = ma->GetClusterRepEdge(i);
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
		int fcl = ma->GetClusterRepFace(i); 
		int first = first_face_dof[i] + excl_grads*face_ngrad[i];
		
		for(j=first; j< first_face_dof[i+1]; j++) 
		  table[fcl][cnt[fcl]++] = j; 
	      }
	  
	  for(i=0;i<ni;i++) 
	    { 
	      // Stack: plane-face + inner
	      int ccl = ma->GetClusterRepElement(i);
	      int first = first_inner_dof[i] + excl_grads*cell_ngrad[i]; 
	      for(j=first;j<first_inner_dof[i+1];j++) 
		table[ccl][cnt[ccl]++] = j; 
	         	      
	      // inner to quad-face (horizontal edge stack) 
	      // each face different stack 
	      
	      auto fanums = ma->GetElFaces(ElementId(VOL,i)); 
	      
	      for(j=0;j<fanums.Size();j++) 
		{ 
		  int fcl = ma->GetClusterRepFace(fanums[j]);
		  if(fcl != ccl)
		    for(k=first;k<first_inner_dof[i+1];k++) 
		      table[fcl][cnt[fcl]++] = k; 
		} 
	      
	    }  
	  
	  for(i=0;i<nfa;i++) 
	    if(fine_face[i]) 
	      { 
		vnums = ma->GetFacePNums(i); 
		if(vnums.Size()==4) continue; 
		int fcl = ma->GetClusterRepFace(i); 
		auto ednums = ma->GetFaceEdges(i); 
		
		for(j=0;j<ednums.Size();j++)
		  {
		    int ecl = ma->GetClusterRepEdge(ednums[j]); 
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
		// int pn1, pn2;
		// ma->GetEdgePNums (i,pn1,pn2);
		IVec<2> vts = ma->GetEdgePNums (i);
                int pn1 = vts[0], pn2 = vts[1];
                
		pn1 = ma->GetClusterRepVertex(pn1);
		pn2 = ma->GetClusterRepVertex(pn2);
		//	int ecl = ma->GetClusterRepEdge(i);
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
	      auto f2ed = ma->GetFaceEdges (i);
	      for (k = 0; k < f2ed.Size(); k++)
		for (j = 0; j < ndof; j++)
		  table[nv+f2ed[k]][cnt[nv+f2ed[k]]++] = first+j;
	    }
	 
	  for (i = 0; i < ni; i++)
	    {
	      auto enums = ma->GetElEdges (ElementId(VOL,i));
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
	  int ds_order = int(precflags.GetNumFlag ("ds_order", 0));
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
    
    return make_shared<Table<int>> (table);
#endif
    
  }

    

  shared_ptr<Array<int>> HCurlHighOrderFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    // int nv = ma->GetNV();
    size_t ne = ma->GetNE();
    // int nse = ma->GetNSE();

    size_t ned = ma->GetNEdges();
    size_t nfa = (ma->GetDimension() == 2) ? 0 : ma->GetNFaces();


    cout << IM(3) << "called createdirectsolverclusters" << endl;
    // 
    if (precflags.NumFlagDefined ("ds_order"))
      {
	int ds_order = int (precflags.GetNumFlag ("ds_order", 1));

        auto spclusters = make_shared<Array<int>> (GetNDof());
	Array<int> & clusters = *spclusters;
	clusters = 0;
	
	size_t ned = ma->GetNEdges();
	for (int i = 0; i < ned; i++)
	  clusters[i] = 1;

	for (size_t i = 0; i < ned; i++)
	  {
	    int first = first_edge_dof[i];
	    int next = first_edge_dof[i+1];
	    for (int j = 0; (j < ds_order) && (first+j < next) ; j++)
	      clusters[first+j] = 1;
	  }

	size_t nfa = ma->GetNFaces();
	for (size_t i = 0; i < nfa; i++)
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

	return spclusters;
      }


    int clustertype = int(precflags.GetNumFlag("ds_cluster",4));  

    cout << IM(3) << " DirectSolverCluster Clustertype " << clustertype << endl;
    if(clustertype==0)
      return(0);
	   
    // int nv = ma->GetNV();
    // int nd = GetNDof();
    // int ne = ma->GetNE();
    // int ned = ma->GetNEdges();

    Array<int> ednums, fnums, pnums;





    if (precflags.GetDefineFlag("subassembled"))
      {
        auto spclusters = make_shared<Array<int>> (GetNDof());
	Array<int> & clusters = *spclusters;
	clusters = 0;


	for (int i = 0; i < ned; i++)
	  if (!IsDirichletEdge(i) && fine_edge[i])
	    clusters[i] = 1;

	return spclusters;
      }




    // int i, j, k;
    bool hasprism = false;

    for (size_t i = 0; !hasprism && i < ne; i++)
      if (ma->GetElType(ElementId(VOL, i)) == ET_PRISM)
	hasprism = true;
    
    if (!hasprism && adddirectsolverdofs.Size() == 0 &&
	directedgeclusters.Size() == 0 && directfaceclusters.Size() == 0 && directelementclusters.Size() == 0) 
      return NULL;

    auto spclusters = make_shared<Array<int>> (GetNDof());
    Array<int> & clusters = *spclusters;
    clusters = 0;


    if(directedgeclusters.Size() != 0 || 
       directfaceclusters.Size() != 0 ||
       directelementclusters.Size() != 0)
      {
	for(size_t i=0; i<directedgeclusters.Size(); i++)
	  for(size_t j=first_edge_dof[directedgeclusters[i]]; j<first_edge_dof[directedgeclusters[i]+1]; j++)
	    clusters[j] = 6;

	for(size_t i=0; i<directfaceclusters.Size(); i++)
	  {
	    for(size_t j=first_face_dof[directfaceclusters[i]];  j< first_face_dof[directfaceclusters[i]] + face_ngrad[directfaceclusters[i]]; j++)
	      clusters[j] = 6;
	  }

	for(size_t i=0; i<directelementclusters.Size(); i++)
	  {
	    for(size_t j=first_inner_dof[directelementclusters[i]];  j< first_inner_dof[directelementclusters[i]] + cell_ngrad[directelementclusters[i]]; j++)
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
	    for (size_t i = 0; i < ne; i++)
	      {
                ElementId ei(VOL, i);
		if (ma->GetElType(ei) == ET_PRISM)
		  {
		    auto ednums = ma->GetElEdges (ei);
		
		    for (int j = 6; j < 9; j++)  //vertical Edges 
		      { 
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (int k = first; k < next; k++)
			  clusters[k] = 1;      //verthoedge
			clusters[ednums[j]] = 1;
		      }

		    fnums = ma->GetElFaces (ei); 

		    for (int j=2; j<5 ; j++)  // vertical faces
		      {
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			for(int k=first; k<next; k++)
			  clusters[k] = 1;
		    
			// IVec<2> p = order_face[fnums[j]]; 
			// for(k=next-p[0]-p[1]; k<next; k++)
			// clusters[k] = 1;
		      }

		    int first = first_inner_dof[i]; 
		    int next = first_inner_dof[i+1]; 

		    for(int k=first; k<next; k++)
		      clusters[k] = 1;
		  }
	      }
	    break; 
	  case 2: 
	    //clusters = 0;
	
	    // All Vertical Edges in one Cluster for Hex and Prism + const_x*poly faces + const_y*polx faces (-> 2d Problems !) 
	
	
	    //lo 
	    for(size_t i=0;i<ma->GetNEdges();i++)
	      clusters[i]=1; 

       
	    for (size_t i = 0; i < ne; i++)
	      {
                ElementId ei(VOL, i);
		/*if (ma->GetElType(i) == ET_PYRAMID)
		  {
		  GetDofNrs(i,ednums);
		  for(j=0; j<ednums.Size(); j++) clusters[ednums[j]] = 3;
	    
		  }   
		*/
		if (ma->GetElType(ei) == ET_PRISM)
		  {
		    auto ednums = ma->GetElEdges (ei);
		    for (int j = 0; j < 6; j++)  //horizontal Edges 
		      { 
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (int k = first; k < next; k++)
			  clusters[k] = 1; 
			clusters[ednums[j]]=1;
		      }
		
		
		    for (int j = 6; j < 9; j++)  //vertical Edges 
		      { 
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (int k = first; k < next; k++)
			  clusters[k] = 1;      //verthoedge
			clusters[ednums[j]]=1; //ned
		      }
		    fnums = ma->GetElFaces (ei); 
		
		    for (int j=2; j<5 ; j++)  // vertical faces
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			// TEST SZ eigentlich nurin eine richtung konst ausreichend
			for (int k=first; k < next; k++) 
			  clusters[k]=1; 
		    
		    
			IVec<2> p = order_face[fnums[j]]; 
		    
			for(int k=next-p[0]-p[1]; k<next; k++)
			  clusters[k] = 1;
			for(int k=0; k<4; k++)
			  clusters[k] = 1;
		    
		    
			// for (k=first; k < next; k++) 
			//   clusters[k]=3; 
		    
		      }
		    for (int j=0; j<2 ; j++)  // horizontal faces
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			for (int k=first; k < next; k++) 
			  clusters[k]=1; 
		      }
		  }
	    
		else if (ma->GetElType(ei) == ET_HEX)
		  {
		    auto ednums = ma->GetElEdges (ei);
		
		    for(int j=0;j<8;j++) // horizontal edges
		      {	
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (int k = first; k < next; k++)
			  clusters[k] = 0;
		
			clusters[ednums[j]]=0; 
		    
		      }
			
		    for (int j = 8; j < 12; j++)  // vertical edges 
		      {	
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (int k = first; k < next; k++)
			  clusters[k] = 3;
			clusters[ednums[j]]=0; 
		      }
		
		    fnums = ma->GetElFaces(ei); // vertical faces 
		    for(int j=2;j<6;j++) 
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			/*	for (k=first; k < next; k++) 
			  clusters[k]=0; 
		      
			  IVec<2> p = order_face[fnums[j]];
			  for(k=2*(p[0]+1)*(p[0]+1);k<next;k++) //achtung
			  clusters[k]=3;  
			*/ 
			// TEST SZ  eigentlich in eine richtung konstante nur benoetigt
			for (int k=first; k < next; k++) 
			  clusters[k]=3; 
		    
		      }
		    for(int j=0;j<2;j++)  //horizontal faces 
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			for (int k=first; k < next; k++) 
			  clusters[k]=0; 
		      } 
		  }
	    
		for(int k=first_inner_dof[i];k<first_inner_dof[i+1];k++) 
		  clusters[k]=0; 
	    
	      }
    
	    break; 
	  case 3: 
	    //lo 
	    //for(i=0;i<ma->GetNEdges();i++)
	    //  clusters[i]=0; 
	
	    for (size_t i = 0; i < ne; i++)
	      {
                ElementId ei(VOL, i);
		// auto pnums = ma->GetElPNums(ei); 
		if (ma->GetElType(ei) == ET_PRISM)
		  {
		    auto ednums = ma->GetElEdges (ei);
		    for (int j = 0; j < 6; j++)  //horizontal Edges 
		      { 
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (int k = first; k < next; k++)
			  clusters[k] = 0;
			clusters[ednums[j]]=0; 
		      }
		    for (int j = 6; j < 9; j++)  //vertical Edges 
		      { 
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (int k = first; k < next; k++)
			  clusters[k] = 3;      //verthoedge
			clusters[ednums[j]]=0; //ned
		      }
		    fnums = ma->GetElFaces (ei); 
		
		    for (int j=2; j<5 ; j++)  // vertical faces
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
		    
			for (int k=first; k < next; k++) 
			  clusters[k]=3; 
		    
		      }
		    for (int j=0; j<2 ; j++)  // horizontal faces
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			for (int k=first; k < next; k++) 
			  clusters[k]=0; 
		      }
		  }
	    
		else if (ma->GetElType(ei) == ET_HEX)
		  {
		    auto ednums = ma->GetElEdges (ei);
		    for(int j=0;j<8;j++) // horizontal edges
		      {	
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (int k = first; k < next; k++)
			  clusters[k] = 0;
		    
			clusters[ednums[j]]=0; 
		    
		      }
		    for (int j = 8; j < 12; j++)  // vertical edges 
		      {	
			int first = first_edge_dof[ednums[j]];
			int next = first_edge_dof[ednums[j]+1];
			for (int k = first; k < next; k++)
			  clusters[k] = 3;
			clusters[ednums[j]]=0; 
		      }
		
		    fnums = ma->GetElFaces(ei); // vertical faces 
		    for(int j=2;j<6;j++) 
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			/*	for (k=first; k < next; k++) 
			  clusters[k]=0; 
		      
			  IVec<2> p = order_face[fnums[j]];
			  for(k=2*(p[0]+1)*(p[0]+1);k<next;k++) //achtung
			  clusters[k]=3;  
			*/ 
			// TEST SZ  eigentlich in eine richtung konstante nur benoetigt
			for (int k=first; k < next; k++) 
			  clusters[k]=3; 
		    
		      }
		    for(int j=0;j<2;j++)  //horizontal faces 
		      {
		    
			int first = first_face_dof[fnums[j]]; 
			int next = first_face_dof[fnums[j]+1]; 
		    
			for (int k=first; k < next; k++) 
			  clusters[k]=0; 
		      } 
		  }
	    
		for(int k=first_inner_dof[i];k<first_inner_dof[i+1];k++) 
		  clusters[k]=0; 
	    
	      }
    
	    //  (*testout) << "direct clusters = " << endl << clusters << endl;
	
	    for(size_t i=0; directsolverclustered.Size() > 0 && i<ne; i++)
	      {
                ElementId ei(VOL,i);
		if(directsolverclustered[ma->GetElIndex(ei)])
		  {
		    GetDofNrs(i,ednums);
		    for(int k=0; k<ednums.Size(); k++)
		      {
			clusters[ednums[k]] = 4;
		      }
		  }
	      }
	    break; 

	  case 4:  // just like the old hcurl 
	    //lo 
	    //clusters = 0; 
	    for(size_t i=0;i<ned;i++) 
	      { 
		// int pi1,pi2; 
		// ma->GetEdgePNums(i,pi1,pi2);
		IVec<2> vts = ma->GetEdgePNums (i);
                int pi1 = vts[0], pi2 = vts[1];
                
		pi1 = ma->GetClusterRepVertex (pi1); 
		pi2 = ma->GetClusterRepVertex (pi2);
		if(pi1 == pi2) 
		  // stack of horizontal edge-blocks (without low-order)
		  // -> decoupled 
		  for(int j=first_edge_dof[i];j <first_edge_dof[i+1];j++) 
		    clusters[j] = 1; 
	      }

	    for(size_t i=0;i<nfa;i++)
	      {
		// Attenzione das is zuviel !!! 

		auto pnums = ma->GetFacePNums(i); 
		if(pnums.Size() == 4) // quad face 
		  { 
		    //IVec<2> p = order_face[i];         
		    int first =first_face_dof[i];  
		    int next = first_face_dof[i+1]; 
		    // Ned_0*pol_z 
		    // int first = next - p[0] - p[1]; 
	
		    for(int j= first ; j<next; j++)
		      clusters[j] = 1; 
		  } 
	      }
	    break; 
	  case 5:  // just like the old hcurl horizontal only constant ... 
	    //lo 
	    //clusters = 0; 
	    for(size_t i=0;i<ned;i++) 
	      { 
		// int pi1,pi2; 
		// ma->GetEdgePNums(i,pi1,pi2);
		IVec<2> vts = ma->GetEdgePNums (i);
                int pi1 = vts[0], pi2 = vts[1];
                
		pi1 = ma->GetClusterRepVertex (pi1); 
		pi2 = ma->GetClusterRepVertex (pi2);
		if(pi1 == pi2) 
		  // stack of horizontal edge-blocks (without low-order)
		  // -> decoupled 
		  {
		    //clusters[i] = 1; 
		    for(int j=first_edge_dof[i];j <first_edge_dof[i+1];j++) 
		      clusters[j] = 1; 
		  }
	      }

	    for(size_t i=0;i<nfa;i++)
	      {
		// Attenzione das is zuviel !!! 

		auto pnums = ma->GetFacePNums(i); 
		if(pnums.Size() == 4) // quad face 
		  { 
		    IVec<2> p = order_face[i];        
		    // int first =first_face_dof[i];  
		    int next = first_face_dof[i+1]; 
		    // Ned_0*pol_z 
		    int first = next - p[0] - p[1]; 
	
		    for(int j= first ; j<next; j++)
		      clusters[j] = 1; 
		  } 
	  
	    
	    
	      }


	    //  (*testout) << "direct clusters = " << endl << clusters << endl;
	
	    for(size_t i=0; directsolverclustered.Size() > 0 && i<ne; i++)
	      {
                ElementId ei(VOL, i);
		if(directsolverclustered[ma->GetElIndex(ei)])
		  {
		    ELEMENT_TYPE eltype = ma->GetElType(ei); 
		    if(eltype != ET_PRISM) continue; 
		
		    GetDofNrs(i,ednums);
		    for(int k=0; k<ednums.Size(); k++)
		      if(ednums[k]>=0) clusters[ednums[k]] = 2;
		  }
	      }
	    break; 
      
      
	  }

      }

    //if(adddirectsolverdofs.Size())
    //  (*testout) << "addirectsolverdofs " << adddirectsolverdofs << endl;

    //int numadd = 0;

    for(size_t i=0; i< adddirectsolverdofs.Size(); i++)
      {
	//if(clusters[adddirectsolverdofs[i]] != 5)
	//  numadd++;
	clusters[adddirectsolverdofs[i]] = 5;
      }

    //(*testout) << "number added dofs " << numadd << endl;

    

    //(*testout) << "clusters = " << clusters << endl;
    return spclusters;
    
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
      ma->GetEdgePNums (i, pi1, pi2);
	
      if (ma->GetClusterRepVertex (pi1) ==
      ma->GetClusterRepVertex (pi2))
      has_cluster = 1;
      }

      if (!has_cluster) return 0;
    
   
      BitArray & ba = *new BitArray (GetNDof());
      ba.Clear();

      for (i = 0; i < ned; i++)
      {
      int pi1, pi2;
      ma->GetEdgePNums (i, pi1, pi2);
	
      if (ma->GetClusterRepVertex (pi1) ==
      ma->GetClusterRepVertex (pi2))
      {
      ba.Set (i);
      for (int l = first_edge_dof[i];
      l < first_edge_dof[i+1]; l++)
      ba.Set (l);
      }
      }
    
    
      for (i = 0; i < nfa; i++)
      {
      ma->GetFacePNums (i, vnums);
      if (vnums.Size() == 4)
      {
      for (int l = first_face_dof[i];
      l < first_face_dof[i+1]; l++)
      ba.Set (l);
      }	    
      }
      return &ba;
      }*/

  shared_ptr<H1HighOrderFESpace> HCurlHighOrderFESpace::CreateGradientSpace() const
  {
    Flags flags2(flags);
    if(iscomplex)
      flags2.SetFlag("complex");

    flags2.SetFlag("order", order+1 - (type1? 1: 0));
    
    flags2.SetFlag("print");
    if(uniform_order_inner>-1)
      flags2.SetFlag("orderinner",uniform_order_inner+1);
    if(uniform_order_face>-1)
      flags2.SetFlag("orderface",uniform_order_face+1);
    if(uniform_order_edge>-1)
      flags2.SetFlag("orderedge",uniform_order_edge+1);
    auto fesh1 = make_shared<H1HighOrderFESpace>(ma, flags2);

    BitArray h1def(ma->GetNDomains());
    if(definedon[VOL].Size() == 0)
      h1def.Set();
    else
      h1def.Clear();
    for(int i=0; i<definedon[VOL].Size(); i++)
      if(definedon[VOL][i])
	h1def.SetBit(i);
    fesh1->SetDefinedOn(VOL,h1def);

    BitArray h1defb(ma->GetNBoundaries());
    if(definedon[BND].Size() == 0)
      h1defb.Set();
    else
      h1defb.Clear();
    for(int i=0; i<definedon[BND].Size(); i++)
      if(definedon[BND][i])
	h1defb.SetBit(i);
    fesh1->SetDefinedOn(BND,h1defb);

    int ne = ma -> GetNE();
    // int ned = ma->GetNEdges();
    // int nfa = 0;
    int nse = ma->GetNSE();
    // if(ma->GetDimension()==3) nfa = ma->GetNFaces();

    LocalHeap lh(100008, "HCurlHOFeSpace::CreateGradient");
    fesh1->Update();
    fesh1->FinalizeUpdate();

    for(int i=0; i < ne; i++){
      ElementId ei(VOL,i);
      if(!gradientdomains[ma->GetElIndex(ei)]){
	fesh1->SetElementOrder(i, 1);
	for(auto edge : ma->GetElement(ei).Edges()) {
	  fesh1->SetEdgeOrder(edge,1);
	}
	for(auto face : ma->GetElement(ei).Faces()) {
	  fesh1->SetFaceOrder(face,1);
	}
      }
    }

    int value;
    for(int i=0; i<nse; i++) {
      ElementId sei(BND, i);
      //cout << "SElIndex: " << ma->GetSElIndex(i) << endl;
      if(!gradientboundaries[ma->GetElIndex(sei)]){
	value = 1;
      }
      else{
	value = order+1 - (type1? 1: 0);
      }
        auto eledges = ma->GetElEdges(sei);
	for(int j=0;j<eledges.Size();j++){
	  fesh1->SetEdgeOrder(eledges[j],value);
	}
	if(ma->GetDimension()==3) {
	  fesh1->SetFaceOrder(ma->GetSElFace(i),value);
	}
    }
    fesh1->UpdateDofTables();
    fesh1->UpdateCouplingDofArray();
    fesh1->FinalizeUpdate();

    return fesh1;
  }

  shared_ptr<SparseMatrix<double>>
  HCurlHighOrderFESpace::CreateGradient(const H1HighOrderFESpace & fesh1) const
  {
    int ne = ma->GetNE();
    int ned = ma->GetNEdges();
    int nfa = 0;
    if(ma->GetDimension()==3) nfa = ma->GetNFaces();

    Array<int> dnums_h1l;
    Array<int> dnums_hcl;

    Array<int> cnts(ndof);
    cnts = 0;

    for(int i=0; i<ned; i++)
      {
	if(fine_edge[i])
	  {
	    cnts[i] = 2;  // vertices-nedelec
	    int l = first_edge_dof[i];
            IntRange edofs = fesh1.GetEdgeDofs(i);
	    for( int k = edofs.First(); k < edofs.Next(); k++, l++)
	      cnts[l] = 1;
	  }
      }
    
    for(int i=0; i< nfa; i++)
      {
	if(fine_face[i])
	  {
	    int l= first_face_dof[i]; 
            IntRange fdofs = fesh1.GetFaceDofs(i);
	    for (int k = fdofs.First(); k < fdofs.Next(); k++, l++) 
	      cnts[l] = 1;
	  }
      }
    
    for(int i=0; i<ne; i++)
      {
	int l= first_inner_dof[i]; 
        IntRange edofs = fesh1.GetElementDofs(i);
	for (int k = edofs.First(); k < edofs.Next(); k++, l++) 
	  cnts[l] = 1; 
      }
  
    //sparse matrix with above matrix graph cnts 
    auto grad = make_shared<SparseMatrix<double>>(cnts, fesh1.GetNDof()); 
    
    // vertices - nedelec
    // ho-edges
    for (int i = 0; i < ned; i++) 
      {
	if(fine_edge[i])
	  { 
            IVec<2> vts = ma->GetEdgePNums (i);
            int p1 = vts[0], p2 = vts[1];

	    grad->CreatePosition(i,p1); 
	    grad->CreatePosition(i,p2); 

	    if (p1 < p2) 
	      {	      
		(*grad)(i,p1) = -1.;
		(*grad)(i,p2) = 1.; 
	      }
	    else
	      {
		(*grad)(i,p1) = 1.; 
		(*grad)(i,p2) = -1.; 
	      }
	   
	    int l = first_edge_dof[i]; 
            
            IntRange edofs = fesh1.GetEdgeDofs(i);

	    for( int k = edofs.First(); k < edofs.Next(); k++, l++)
	      {
		grad->CreatePosition(l,k); 
		(*grad)(l,k)=1.; 
	      }
	  }
      }
    
    for(int i=0; i< nfa; i++) 
      {
	if(fine_face[i])
	  {
	    int l= first_face_dof[i]; 
            IntRange fdofs = fesh1.GetFaceDofs(i);
	    for (int k = fdofs.First(); k < fdofs.Next(); k++, l++) 
	      {
		grad->CreatePosition(l,k); 
		(*grad)(l,k)=1.;
	      }
	  }
      }

    for(int i=0; i<ne; i++)
      {
	int l= first_inner_dof[i]; 
        IntRange edofs = fesh1.GetElementDofs(i);
	for (int k = edofs.First(); k < edofs.Next(); k++, l++) 
	  {
	    grad->CreatePosition(l,k); 
	    (*grad)(l,k)=1.;
	  }	
      }
    //(*testout) << " Global grad " << grad << endl; 

    return grad;
    
  }
  /*
  SparseMatrix<double> * 
  HCurlHighOrderFESpace :: CreateGradient() const
  {
    Flags flags2 (flags);

    if(iscomplex)
      flags2.SetFlag("complex");
    flags2.SetFlag("order", order+1);
    flags2.SetFlag("relorder", rel_order+1); 
    // flags2.SetFlag("orderinner",uniform_order_inner+1); 
    // flags2.SetFlag("orderface",uniform_order_face+1); 
    // flags2.SetFlag("orderedge",uniform_order_edge+1); 
    flags2.SetFlag ("orderquad", -1);
    flags2.SetFlag ("ordertrig", -1);  
    flags2.SetFlag("variableorder",var_order); 
    
    flags2.SetFlag("print");
   
    // if the following three flags are set -> relorder is used 
    flags2.SetFlag("relorder",rel_order+1);
    flags2.SetFlag("order",order+1); 
    // flags2.SetFlag("variableorder"); 
    
    
       if(var_order)
       flags2.SetFlag("relorder",rel_order+1); 
       else
       flags2.SetFlag("order",order+1); 
    

    if(uniform_order_inner>-1)
      flags2.SetFlag("orderinner",uniform_order_inner+1); 
    if(uniform_order_face>-1)
      flags2.SetFlag("orderface",uniform_order_face+1); 
    if(uniform_order_edge>-1)
      flags2.SetFlag("orderedge",uniform_order_edge+1); 
    
    // flags2.SetFlag ("orderquad", -1);
    // flags2.SetFlag ("ordertrig", -1);  
    // flags2.SetFlag("variableorder",var_order); 

    *testout << "Flags for h1space: " << endl << flags2 << endl;
    H1HighOrderFESpace  fesh1(ma, flags2); 

    BitArray h1def(ma->GetNDomains());
    if(definedon.Size() == 0)
      h1def.Set();
    else
      h1def.Clear();
    for(int i=0; i<definedon.Size(); i++)
      if(definedon[i] > 0)
	h1def.Set(i);
    fesh1.SetDefinedOn(h1def);

    BitArray h1defb(ma->GetNBoundaries());
    if(definedonbound.Size() == 0)
      h1defb.Set();
    else
      h1defb.Clear();
    for(int i=0; i<definedonbound.Size(); i++)
      if(definedonbound[i] > 0)
	h1defb.Set(i);
    fesh1.SetDefinedOnBoundary(h1defb);

    LocalHeap lh(100008, "HCurlHOFeSpace::CreateGradient");
    fesh1.Update(lh);
    fesh1.FinalizeUpdate(lh);
     
    int ned = ma->GetNEdges(); 
    // int nv  = ma->GetNV(); 
    int ne = ma->GetNE(); 
    int nfa = 0;
    if(ma->GetDimension()==3) nfa= ma->GetNFaces();  

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
            IntRange edofs = fesh1.GetEdgeDofs(i);
	    for( int k = edofs.First(); k < edofs.Next(); k++, l++)
	      cnts[l] = 1;
	  }
      }
    
    for(int i=0; i< nfa; i++) 
      {
	if(fine_face[i])
	  {
	    int l= first_face_dof[i]; 
            IntRange fdofs = fesh1.GetFaceDofs(i);
	    for (int k = fdofs.First(); k < fdofs.Next(); k++, l++) 
	      cnts[l] = 1;
	  }
      }
    
    for(int i=0; i<ne; i++)
      {
	int l= first_inner_dof[i]; 
        IntRange edofs = fesh1.GetElementDofs(i);
	for (int k = edofs.First(); k < edofs.Next(); k++, l++) 
	  cnts[l] = 1; 
      }
  
    //sparse matrix with above matrix graph cnts 
    SparseMatrix<double> & grad = *new SparseMatrix<double>(cnts, fesh1.GetNDof()); 
    
    // vertices - nedelec
    // ho-edges
    for (int i = 0; i < ned; i++) 
      {
	if(fine_edge[i])
	  { 
	    int p1,p2; 
	    ma->GetEdgePNums(i,p1,p2); 

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
            
            IntRange edofs = fesh1.GetEdgeDofs(i);

	    for( int k = edofs.First(); k < edofs.Next(); k++, l++)
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
            IntRange fdofs = fesh1.GetFaceDofs(i);
	    for (int k = fdofs.First(); k < fdofs.Next(); k++, l++) 
	      {
		grad.CreatePosition(l,k); 
		grad(l,k)=1.;
	      }
	  }
      }

    for(int i=0; i<ne; i++)
      {
	int l= first_inner_dof[i]; 
        IntRange edofs = fesh1.GetElementDofs(i);
	for (int k = edofs.First(); k < edofs.Next(); k++, l++) 
	  {
	    grad.CreatePosition(l,k); 
	    grad(l,k)=1.;
	  }	
      }
    //(*testout) << " Global grad " << grad << endl; 

    return &grad;
  }

*/





  static RegisterFESpace<HCurlHighOrderFESpace> init ("hcurlho");
}
