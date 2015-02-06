/*********************************************************************/
/* File:   h1hofespace.cpp                                           */
/* Author: Start                                                     */
/* Date:   10. Feb. 2003                                             */
/*********************************************************************/

/**
   High Order Finite Element Space
*/

#include <comp.hpp>
#include <multigrid.hpp> 
#include "../fem/h1hofe.hpp"
#include "../fem/h1hofefo.hpp"

using namespace ngmg; 


#ifdef PARALLEL

#include "../parallel/dump.hpp"


template <NODE_TYPE NT, typename TELEM, typename TSIZE>
class NodalArray
{
  const MeshAccess & ma;
  Array<TELEM,TSIZE> & a;
public:
  NodalArray (const MeshAccess & ama, Array<TELEM,TSIZE> & aa) : ma(ama), a(aa) { ; }
  const MeshAccess & GetMeshAccess() const { return ma; }
  Array<TELEM,TSIZE> & A() { return a; }
};

template <NODE_TYPE NT, typename TELEM, typename TSIZE>
auto NodalData (const MeshAccess & ama, Array<TELEM,TSIZE> & a) -> NodalArray<NT,TELEM,TSIZE> 
{ return NodalArray<NT,TELEM,TSIZE> (ama, a); }



template <NODE_TYPE NT, typename T, typename TSIZE> 
Archive & operator & (Archive & archive, NodalArray<NT,T,TSIZE> && a)
{
  if (MyMPI_GetNTasks() == 1) return archive & a.A();
  
  auto g = [&] (int size) { archive & size; };    

  typedef typename key_trait<NT>::TKEY TKEY;
  auto f = [&] (TKEY key, T val) { archive & val; };
      
  GatherNodalData<NT> (a.GetMeshAccess(), a.A(), g, f);

  return archive;
}



#else
template <NODE_TYPE NT, typename TELEM>
auto NodalData (MeshAccess & ma, Array<TELEM> & a) -> Array<TELEM> & { return a; }
#endif


namespace ngcomp
{

  H1HighOrderFESpace ::  
  H1HighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name = "H1HighOrderFESpace(h1ho)";
    // define h1ho flags
    DefineDefineFlag("h1ho");
    DefineNumFlag("relorder");
    DefineNumFlag("orderinner");
    DefineNumFlag("orderface");
    DefineNumFlag("orderedge");
    DefineNumFlag("orderquad");
    DefineNumFlag("ordertrig");
    DefineNumFlag("variableorder"); 
    //  DefineNumListFlag("dom_order_min_x");
    //  DefineNumListFlag("dom_order_max_x");
    //  DefineNumListFlag("dom_order_min_y");
    //  DefineNumListFlag("dom_order_max_y");
    //  DefineNumListFlag("dom_order_min_z");
    //  DefineNumListFlag("dom_order_max_z");
    DefineNumFlag("smoothing");
    DefineDefineFlag("wb_withedges");
    if (parseflags) CheckFlags(flags);

    wb_loedge = ma->GetDimension() == 3;
    if (flags.GetDefineFlag("wb_withedges")) wb_loedge = true;
    if (flags.GetDefineFlag("wb_withoutedges")) wb_loedge = false;

    
    // Variable order space: 
    //      in case of (var_order && order) or (relorder) 
    var_order = flags.GetDefineFlag("variableorder");  
    fixed_order = flags.GetDefineFlag("fixedorder");  
    order = int (flags.GetNumFlag ("order",1)); 
    if (order < 1) order = 1;

    if(flags.NumFlagDefined("relorder") && !flags.NumFlagDefined("order")) 
      var_order = 1; 
    
    rel_order=int(flags.GetNumFlag("relorder",order-1)); 

    print = flags.GetDefineFlag("print");

    if(flags.NumFlagDefined("order") && flags.NumFlagDefined("relorder")) 
      {
	if(var_order)
	  cerr << " WARNING: H1HoFeSpace: inconsistent flags: variableorder, order and relorder "
	       << "-> variable order space with rel_order " << rel_order << "is used, but order is ignored " << endl; 
	else 
	  cerr << " WARNING: H1HoFeSpace: inconsistent flags: order and rel_order "
	       << "-> uniform order space with order " << order << " is used " << endl; 
      }
    
    uniform_order_inner = int (flags.GetNumFlag ("orderinner", -1));
    uniform_order_face = int (flags.GetNumFlag ("orderface", -1));
    uniform_order_edge = int (flags.GetNumFlag ("orderedge", -1));
    uniform_order_quad = int (flags.GetNumFlag ("orderquad", -1));
    uniform_order_trig = int (flags.GetNumFlag ("ordertrig", -1));
    
    if (flags.NumFlagDefined("smoothing")) 
      throw Exception ("Flag 'smoothing' for fespace is obsolete \n Please use flag 'blocktype' in preconditioner instead");
    nodalp2 = flags.GetDefineFlag ("nodalp2");
          
    Flags loflags;
    loflags.SetFlag ("order", 1);
    loflags.SetFlag ("dim", dimension);
    if (iscomplex) loflags.SetFlag ("complex");
    // if (timing) loflags.SetFlag ("timing");
    if (flags.NumListFlagDefined ("dirichlet")) 
      loflags.SetFlag ("dirichlet", flags.GetNumListFlag ("dirichlet"));
    if (dgjumps){ *testout << "(L2HOFES:)setting loflag dgjumps " << endl; loflags.SetFlag ("dgjumps");}

    if (!no_low_order_space)
      low_order_space = make_shared<NodalFESpace> (ma, loflags);
    switch (ma->GetDimension())
      {
      case 1:
        {
          evaluator = make_shared<T_DifferentialOperator<DiffOpId<1>>>();
          flux_evaluator = make_shared<T_DifferentialOperator<DiffOpGradient<1>>>();
          boundary_evaluator = make_shared<T_DifferentialOperator<DiffOpIdBoundary<1>>>();
          break;
        }
      case 2:
        {
          evaluator = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
          flux_evaluator = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
          boundary_evaluator = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
          break;
        }
      case 3:
        {
          evaluator = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
          flux_evaluator = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
          boundary_evaluator = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
          break;
        }
      }

    if (dimension > 1)
      {
	evaluator = make_shared<BlockDifferentialOperator> (evaluator, dimension);
	boundary_evaluator = 
	  make_shared<BlockDifferentialOperator> (boundary_evaluator, dimension);
      }

    auto one = make_shared<ConstantCoefficientFunction> (1);
    integrator = CreateBFI("mass", ma->GetDimension(), one);
    boundary_integrator = CreateBFI("robin", ma->GetDimension(), one);
    if (dimension > 1)
      {
	integrator = make_shared<BlockBilinearFormIntegrator> (integrator, dimension);
        boundary_integrator = make_shared<BlockBilinearFormIntegrator> (boundary_integrator, dimension);
      }

    prol = make_shared<LinearProlongation> (*this);
  }


  H1HighOrderFESpace :: ~H1HighOrderFESpace ()
  {
    ;
  }


  void H1HighOrderFESpace :: Update(LocalHeap & lh)
  {
    static Timer timer ("H1HighOrderFESpace::Update");
    RegionTimer reg(timer);

    FESpace :: Update (lh);

    TORDER maxorder = 0;
    TORDER minorder = 99; 

    if (low_order_space) low_order_space -> Update(lh);
    
    int dim = ma->GetDimension();
    int nv = ma->GetNV();
    int ned = (dim <= 1) ? 0 : ma->GetNEdges();
    int nfa = (dim <= 2) ? 0 : ma->GetNFaces();
    int ne = ma->GetNE();

    used_edge.SetSize(ned); 
    used_face.SetSize(nfa); 
    used_vertex.SetSize(nv); 

    used_edge = false; 
    used_face = false; 
    used_vertex = false; 

    for (FESpace::Element el : Elements (VOL))
      {
        used_vertex[el.Vertices()] = true;
        if (dim >= 2) used_edge[el.Edges()] = true;
        if (dim == 3) used_face[el.Faces()] = true;
      }
    
    for (FESpace::Element el : Elements (BND))
      {
        used_vertex[el.Vertices()] = true;
        if (dim >= 2) used_edge[el.Edges()] = true;
        if (dim == 3) used_face[el.Faces()] = true;
      }
    
    ma->AllReduceNodalData (NT_VERTEX, used_vertex, MPI_LOR);
    ma->AllReduceNodalData (NT_EDGE, used_edge, MPI_LOR);
    ma->AllReduceNodalData (NT_FACE, used_face, MPI_LOR);

    
    order_edge.SetSize (ned);
    order_face.SetSize (nfa);
    order_inner.SetSize (ne);

    int p = var_order ?  1 : order; 
    
    order_edge = p; 
    order_face = p; 
    order_inner = p;
	
    Array<int> eledges, elfaces, vnums;
    if(var_order) 
      for (Ngs_Element el : ma->Elements<VOL>())
        {	
          if (!DefinedOn (el)) continue;
          int i = el.Nr();
          
          ELEMENT_TYPE eltype = el.GetType(); 
          const FACE * faces = ElementTopology::GetFaces (eltype);
          const EDGE * edges = ElementTopology::GetEdges (eltype);
          const POINT3D * points = ElementTopology :: GetVertices (eltype);

          vnums = el.Vertices();
          eledges = el.Edges();
	
          INT<3,TORDER> el_orders = ma->GetElOrders(i) + INT<3> (rel_order); 

          maxorder = max2 (maxorder, Max(el_orders));
          minorder = min2 (minorder, Min(el_orders));
          // for(int l=0;l<3;l++) maxorder = max2(el_orders[l],maxorder); 
          // for(int l=0;l<3;l++) minorder = min2(el_orders[l],minorder); 
          
          order_inner[i] = Max (order_inner[i], el_orders);
          // for(int j=0;j<dim;j++) order_inner[i][j] = max2(order_inner[i][j],el_orders[j]);

          for(int j=0;j<eledges.Size();j++)
            {
              for(int k=0;k<dim;k++)
                if(points[edges[j][0]][k] != points[edges[j][1]][k])
                  { 
                    order_edge[eledges[j]] = max2(order_edge[eledges[j]],el_orders[k]);
                    k=dim; 
                  }
            }
          
          if(dim==3)
            {
              for(int j=0;j<elfaces.Size();j++)
                {
                  // trig_face
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
        }
    
    /* 
       if (ma->GetDimension() == 2 && uniform_order_trig != -1 && uniform_order_quad != -1)
       {
       for (int i = 0; i < nel; i++)
       {
       if (ma->GetElType(i) == ET_TRIG)
       order_inner = INT<3> (uniform_order_trig, uniform_order_trig, uniform_order_trig);
       else
       order_inner = INT<3> (uniform_order_quad, uniform_order_quad, uniform_order_quad);
       }
       }
    */ 
    
    if(uniform_order_inner > -1)  
      order_inner = uniform_order_inner;
    if(uniform_order_face > -1 && dim == 3) 
      order_face = uniform_order_face;
    if(uniform_order_edge > -1)   
      order_edge = uniform_order_edge; 

    for (auto i : Range(used_edge))
      if (!used_edge[i]) order_edge[i] = 1; 

    for (auto i : used_face.Range())
      if (!used_face[i]) order_face[i] = 1; 

    for (ElementId ei : ma->Elements<VOL>())
      if (!DefinedOn(ei)) order_inner[ei.Nr()] = 1;

    if(print) 
      {
	*testout << " H1HoFESpace order " << order << " , var_order " << var_order << " , relorder " << rel_order << endl;  
	(*testout) << "used_vertex (h1): " << used_vertex << endl;
	(*testout) << "used_edge (h1): " << used_edge << endl;
	(*testout) << "used_face (h1): " << used_face << endl;
	
	
	(*testout) << "order_edge (h1): " << order_edge << endl;
	
	(*testout) << "order_face (h1): " << order_face <<  endl;
	(*testout) << "order_inner (h1): " << order_inner << endl;
      }
    
#ifdef SABINE 
    // fuer EE damit ich mit rel_order auch p-ref machen kann :) 
    order=order_inner[0][0]; 

    if (!var_order) { maxorder = order; minorder = order; }; 
    order = maxorder;

    cout << " H1FESPACE : " << minorder << " <= order <= " << maxorder << endl;  
#endif 

    UpdateDofTables ();
    UpdateCouplingDofArray ();
  }


  void H1HighOrderFESpace :: UpdateDofTables ()
  {
    int dim = ma->GetDimension();
    int nv = ma->GetNV();
    int ned = (dim <= 1) ? 0 : ma->GetNEdges();
    int nfa = (dim <= 2) ? 0 : ma->GetNFaces();
    int ne = ma->GetNE();

    int hndof = nv;

    first_edge_dof.SetSize (ned+1);
    for (auto i : Range (ned))
      {
	first_edge_dof[i] = hndof;
	if (order_edge[i] > 1)
	  hndof += order_edge[i] - 1;
      }
    first_edge_dof[ned] = hndof;

    first_face_dof.SetSize (nfa+1);
    for (auto i : Range (nfa))
      {
	first_face_dof[i] = hndof;
	INT<2> p = order_face[i];
	switch(ma->GetFacetType(i))
	  {
	  case ET_TRIG:
            if (p[0] > 2)
              hndof += (p[0]-1)*(p[0]-2)/2;
	    break;
	  case ET_QUAD:
	    if (p[0] > 1 && p[1] > 1)
	      hndof += (p[0]-1)*(p[1]-1);
	    break; 
	  default:
            ;
	  }
      }
    first_face_dof[nfa] = hndof;
 
    first_element_dof.SetSize(ne+1);
    for (auto i : Range(ne))
      {
	first_element_dof[i] = hndof;
	INT<3> p = order_inner[i];	
	switch (ma->GetElType(i))
	  {
	  case ET_TRIG:
	    if(p[0] > 2)
	      hndof += (p[0]-1)*(p[0]-2)/2;
	    break;
	  case ET_QUAD:
	    if(p[0] > 1 && p[1] > 1)
	      hndof += (p[0]-1)*(p[1]-1);
	    break;
	  case ET_TET:
	    if(p[0] > 3)
	      hndof += (p[0]-1)*(p[0]-2)*(p[0]-3)/6;
	    break;
	  case ET_PRISM:
	    if(p[0] > 2 && p[2] > 1)
	      hndof += (p[0]-1)*(p[0]-2)*(p[2]-1)/2;
	    break;
	  case ET_PYRAMID:
	    if(p[0] > 2)
	      hndof += (p[0]-1)*(p[0]-2)*(2*p[0]-3)/6;
	    break;
	  case ET_HEX:
	    if(p[0] > 1 && p[1] > 1 && p[2] > 1) 
	      hndof += (p[0]-1)*(p[1]-1)*(p[2]-1);
	    break;
          case ET_SEGM:
            if (p[0] > 1)
	      hndof += p[0]-1;
            break;
          case ET_POINT:
	    break;
	  }
      } 
    first_element_dof[ne] = hndof;
    ndof = hndof;
   

    if (print)
      {
        (*testout) << "h1 first edge = " << first_edge_dof << endl;
        (*testout) << "h1 first face = " << first_face_dof << endl;
        (*testout) << "h1 first inner = " << first_element_dof << endl;
      }

    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;

    prol->Update();
  }


  void H1HighOrderFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(ndof);

    for (auto i : Range (ma->GetNV()))
      ctofdof[i] = used_vertex[i] ? WIREBASKET_DOF : UNUSED_DOF;

    for (auto edge : Range (ma->GetNEdges()))
      {
	IntRange range = GetEdgeDofs (edge);
	ctofdof[range] = INTERFACE_DOF;
	if (wb_loedge && (range.Size() > 0))
	  ctofdof[range.First()] = WIREBASKET_DOF;
      }

    if (ma->GetDimension() == 3)
      for (auto face : Range (ma->GetNFaces()))
	ctofdof[GetFaceDofs(face)] = INTERFACE_DOF;

    for (int el = 0; el < ma->GetNE(); el++)
      ctofdof[GetElementDofs(el)] = LOCAL_DOF;
    
    if (print)
      *testout << "ctofdof: " << endl << ctofdof << endl;
  }


  void H1HighOrderFESpace :: DoArchive (Archive & archive)
  {
    low_order_space -> DoArchive (archive);
    FESpace::DoArchive(archive);
    archive & level;

    archive & NodalData<NT_EDGE> (*ma, order_edge);
    archive & NodalData<NT_FACE> (*ma, order_face);
    archive & NodalData<NT_CELL> (*ma, order_inner);

    if (archive.Input())
      UpdateDofTables();

    archive & rel_order & var_order & fixed_order & wb_loedge;
    archive & used_vertex & used_edge & used_face;
    archive & uniform_order_inner & uniform_order_face 
      & uniform_order_edge & uniform_order_quad & uniform_order_trig;
    archive & dom_order_min & dom_order_max;
    // archive & smoother;
    archive & ndlevel;
    archive & level_adapted_order & nodalp2;
  }

  const FiniteElement & H1HighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    Ngs_Element ngel = ma->GetElement(elnr);
    ELEMENT_TYPE eltype = ngel.GetType();

    if (!DefinedOn (ElementId (VOL, elnr)))
      {
        switch (eltype)
          {
          case ET_SEGM:    return * new (lh) ScalarDummyFE<ET_SEGM> (); break;
          case ET_TRIG:    return * new (lh) ScalarDummyFE<ET_TRIG> (); break;
          case ET_QUAD:    return * new (lh) ScalarDummyFE<ET_QUAD> (); break;
          case ET_TET:     return * new (lh) ScalarDummyFE<ET_TET> (); break;
          case ET_PYRAMID: return * new (lh) ScalarDummyFE<ET_PYRAMID> (); break;
          case ET_PRISM:   return * new (lh) ScalarDummyFE<ET_PRISM> (); break;
          case ET_HEX:     return * new (lh) ScalarDummyFE<ET_HEX> (); break;
	  case ET_POINT:   break;
	  }
      }

    if (fixed_order && eltype == ET_TRIG)
      {
        switch (order)
          {
          case 1: return *(new (lh) H1HighOrderFEFO<ET_TRIG,1> ()) -> SetVertexNumbers(ngel.vertices);
          case 2: return *(new (lh) H1HighOrderFEFO<ET_TRIG,2> ()) -> SetVertexNumbers(ngel.vertices);
          case 3: return *(new (lh) H1HighOrderFEFO<ET_TRIG,3> ()) -> SetVertexNumbers(ngel.vertices);
          default:
            ; 
          }
      }

    if (fixed_order && eltype == ET_TET)
      {
        switch (order)
          {
          case 1: return *(new (lh) H1HighOrderFEFO<ET_TET,1> ()) -> SetVertexNumbers(ngel.vertices);
          case 2: return *(new (lh) H1HighOrderFEFO<ET_TET,2> ()) -> SetVertexNumbers(ngel.vertices);
          case 3: return *(new (lh) H1HighOrderFEFO<ET_TET,3> ()) -> SetVertexNumbers(ngel.vertices);
          default:
            ; 
          }
      }


    try
      {
        switch (eltype)
          {
          case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, lh);

          case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, lh);
          case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, lh);

          case ET_TET:     return T_GetFE<ET_TET> (elnr, lh);
          case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, lh);
          case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, lh);
          case ET_HEX:     return T_GetFE<ET_HEX> (elnr, lh);

          default:
            throw Exception ("illegal element in H1HoFeSpace::GetFE");
          }
      }
    catch (Exception & e)
      {
        e.Append ("in H1HoFESpace::GetElement\n");
        throw;
      }
  }
 



  template <ELEMENT_TYPE ET>
  const FiniteElement & H1HighOrderFESpace :: T_GetFE (int elnr, LocalHeap & lh) const
  {
    Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL> (elnr);

    H1HighOrderFE<ET> * hofe =  new (lh) H1HighOrderFE<ET> ();
    
    hofe -> SetVertexNumbers (ngel.Vertices());

    switch (int(ET_trait<ET>::DIM))
      {
      case 1:
        {
          hofe -> SetOrderEdge (0, order_inner[elnr][0]);
          break;
        }

      case 2:
        {
          hofe -> SetOrderEdge (order_edge[ngel.Edges()] );
          hofe -> SetOrderFace (0, order_inner[elnr]);
          break;
        }

      case 3: default:  
        {
          hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
          hofe -> SetOrderFace (order_face[ngel.Faces()]);
          hofe -> SetOrderCell (order_inner[elnr]);
          break;
        }
      }

    hofe -> ComputeNDof();
    return *hofe;
  }








  const FiniteElement & H1HighOrderFESpace :: GetSFE (int elnr, LocalHeap & lh) const
  {
    if (!DefinedOnBoundary (ma->GetSElIndex (elnr)))
      {
        switch (ma->GetSElType(elnr))
          {
          case ET_POINT:   return * new (lh) ScalarDummyFE<ET_POINT> (); 
          case ET_SEGM:    return * new (lh) ScalarDummyFE<ET_SEGM> (); 
          case ET_TRIG:    return * new (lh) ScalarDummyFE<ET_TRIG> (); 
          case ET_QUAD:    return * new (lh) ScalarDummyFE<ET_QUAD> (); 
	  default: ;
	  }
      }

    try
      {
        switch (ma->GetSElType(elnr))
          {
          case ET_POINT:   return T_GetSFE<ET_POINT> (elnr, lh);
          case ET_SEGM:    return T_GetSFE<ET_SEGM> (elnr, lh);

          case ET_TRIG:    return T_GetSFE<ET_TRIG> (elnr, lh);
          case ET_QUAD:    return T_GetSFE<ET_QUAD> (elnr, lh);

          default:
            throw Exception ("illegal element in H1HoFeSpace::GetSFE");
          }
      }
    catch (Exception & e)
      {
        e.Append ("in H1HoFESpace::GetSElement\n");
        throw;
      }
  }


  template <ELEMENT_TYPE ET>
  const FiniteElement & H1HighOrderFESpace :: T_GetSFE (int elnr, LocalHeap & lh) const
  {
    Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,BND> (elnr);

    H1HighOrderFE<ET> * hofe =  new (lh) H1HighOrderFE<ET> ();
    
    hofe -> SetVertexNumbers (ngel.vertices);
    
    switch (int (ET_trait<ET>::DIM))
      {
      case 0:
        {
          break;
        }

      case 1:
        {
          hofe -> SetOrderEdge ( order_edge[ngel.Edges()] );
          break;
        }

      case 2: default:  
        {
          hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
	  hofe -> SetOrderFace (0, order_face[ma->GetSElFace(elnr)]);
          break;
        }
      }

    hofe -> ComputeNDof();
    return *hofe;
  }
 

  int H1HighOrderFESpace :: GetNDofLevel (int alevel) const
  {
    return ndlevel[alevel];
  }


  void H1HighOrderFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    Ngs_Element ngel = ma->GetElement(elnr);

    if (!DefinedOn (ngel.GetIndex()))
      {
	dnums.SetSize0();
	return;
      }
    
    dnums = ngel.Vertices();

    if (ma->GetDimension() >= 2)
      for (auto edge : ngel.Edges())
        dnums += GetEdgeDofs (edge);

    if (ma->GetDimension() == 3)
      for (auto face : ngel.Faces())
        dnums += GetFaceDofs (face);

    dnums += GetElementDofs (elnr);
  }


  void H1HighOrderFESpace :: 
  GetDofRanges (ElementId ei, Array<IntRange> & dranges) const
  {
    dranges.SetSize(0);

    if (!DefinedOn (ei)) return;
    Ngs_Element ngel = ma->GetElement(ei);

    for (int i = 0; i < ngel.vertices.Size(); i++)
      dranges.Append (IntRange (ngel.vertices[i], ngel.vertices[i]+1));
         
    if (ma->GetDimension() >= 2)
      for (int i = 0; i < ngel.edges.Size(); i++)
        if (GetEdgeDofs (ngel.edges[i]).Size())
          dranges += GetEdgeDofs (ngel.edges[i]);

    if (ma->GetDimension() == 3)
      for (int i = 0; i < ngel.faces.Size(); i++)
        if (GetFaceDofs (ngel.faces[i]).Size())
          dranges += GetFaceDofs (ngel.faces[i]);

    if (ei.IsVolume())
      if (GetElementDofs (ei.Nr()).Size())
        dranges += GetElementDofs (ei.Nr());
  }



  
  void H1HighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize(1);
    dnums[0] = vnr;
  }
  
  
  void H1HighOrderFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    dnums = GetEdgeDofs (ednr);
  }

  void H1HighOrderFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if (ma->GetDimension() < 3) return;
    dnums = GetFaceDofs (fanr);
  }

  void H1HighOrderFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums = GetElementDofs (elnr);
  }

  void H1HighOrderFESpace :: 
  GetSDofNrs (int elnr, Array<int> & dnums) const
  {
    if (!DefinedOnBoundary (ma->GetSElIndex (elnr)))
      {
	dnums.SetSize0();
	return;
      }

    Ngs_Element ngel = ma->GetSElement(elnr);

    dnums = ngel.Vertices();

    for (int i = 0; i < ngel.edges.Size(); i++)
      dnums += GetEdgeDofs (ngel.edges[i]);

    if (ma->GetDimension() == 3)
      dnums += GetFaceDofs (ngel.faces[0]);
  }
  
  Table<int> * H1HighOrderFESpace :: 
  CreateSmoothingBlocks (const Flags & precflags) const
  {
    bool eliminate_internal = precflags.GetDefineFlag("eliminate_internal");
    bool subassembled = precflags.GetDefineFlag("subassembled");
    // smoothing_types: 
    // 1: 2d V + E + I 
    // 2: 2d VE + I 
    // 3: V + E + F + I 
    // 4: VE + F + I 
    // 5: VE + FI 
    // 6: VEF + I 
    // 7: VEFI 
    // 8: V + E + FI 
    // 9: V + EF + I 
    // 10:V + EF + FI 

    int smoothing_type = int(precflags.GetNumFlag("blocktype",4)); 
    if (subassembled) smoothing_type = 51;

    int nv = ma->GetNV();
    int ned = ma->GetNEdges();
    int nfa = (ma->GetDimension() == 2) ? 0 : ma->GetNFaces();
    int ni = (eliminate_internal) ? 0 : ma->GetNE(); 
   
    cout << " blocktype " << smoothing_type << endl; 
    // cout << " Use H1-Block Smoother:  "; 

    FilteredTableCreator creator(GetFreeDofs());
    for ( ; !creator.Done(); creator++)
      {
	switch (smoothing_type)
	  {

	  case 1:  // 2d: V + E + I
		
	    if (creator.GetMode() == 1)
	      cout << " V + E + I " << endl;
		
	    for (int i = 0; i < nv; i++)
	      creator.Add (i, i);
		
	    for (int i = 0; i < ned; i++)
	      creator.Add (nv+i, GetEdgeDofs(i));

	    for (int i = 0; i < ni; i++)
	      creator.Add (nv+ned+i, GetElementDofs(i));
		
	    break; 
		
	  case 2: // 2d VE + I

	    if (creator.GetMode() == 1)
	      cout << " 2d VE + I " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add(i, i);

	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		for (int k = 0; k < 2; k++)
		  creator.Add (edge.vertices[k], GetEdgeDofs(i));
	      }
		
	    for (int i = 0; i < ni; i++)
	      creator.Add (nv+ned+i, GetElementDofs(i));
		
	    break;


	  case 3: // V + E + F + I
		
	    if (creator.GetMode() == 1)
	      cout << " V + E + F + I " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add(i, i);

	    for (int i = 0; i < ned; i++)
	      creator.Add (nv+i, GetEdgeDofs(i));

	    for (int i = 0; i < nfa; i++)
	      creator.Add(nv+ned+i, GetFaceDofs(i));

	    for (int i = 0; i < ni; i++)
	      creator.Add (nv+ned+nfa+i, GetElementDofs(i));
	    
	    break; 

	  case 4: // VE + F + I
		
	    if (creator.GetMode() == 1)
	      cout << " VE + F + I " << endl;

	    for (int i = 0; i < nv; i++)
	      creator.Add(i, i);
		
	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		for (int k = 0; k < 2; k++)
		  creator.Add (edge.vertices[k], GetEdgeDofs(i));
	      }
	    
	    for (int i = 0; i < nfa; i++)
	      creator.Add(nv+i, GetFaceDofs(i));
	    
	    for (int i = 0; i < ni; i++)
	      creator.Add (nv+nfa+i, GetElementDofs(i));
		
	    break; 

	  case 5: // VE + FI

	    if (creator.GetMode() == 1)
	      cout << " VE + FI " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add(i, i);
		
	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
                for (int k = 0; k < 2; k++)
                  creator.Add (edge.vertices[k], GetEdgeDofs(i));
	      }
	    
	    for (int i = 0; i < nfa; i++)
	      creator.Add(nv+i, GetFaceDofs(i));

	    for (int i = 0; i < ni; i++)
	      {
		const Ngs_Element & ngel = ma->GetElement(i);
		for (int j = 0; j < ngel.faces.Size(); j++)
		  creator.Add (nv+ngel.faces[j], GetElementDofs(i));
	      }
	    break; 


	  case 6: // VEF + I

	    if (creator.GetMode() == 1)
	      cout << " VEF + I " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add (i, i);
		
	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		for (int k = 0; k < 2; k++)
		  creator.Add (edge.vertices[k], GetEdgeDofs(i));
	      }
	    
		
	    for (int i = 0; i < nfa; i++)
	      {
		Ng_Node<2> face = ma->GetNode<2> (i);
		for (int k = 0; k < face.vertices.Size(); k++)
		  creator.Add (face.vertices[k], GetFaceDofs(i));
	      }
	    
	    for (int i = 0; i < ni; i++)
	      creator.Add (nv+i, GetElementDofs(i));
		
	    break; 

	  case 7: // VEFI

	    if (creator.GetMode() == 1)
	      cout << " VEFI " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add (i, i);
		
	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		for (int k = 0; k < 2; k++)
		  creator.Add (edge.vertices[k], GetEdgeDofs(i));
	      }

	    for (int i = 0; i < nfa; i++)
	      {
		Ng_Node<2> face = ma->GetNode<2> (i);
		for (int k = 0; k < face.vertices.Size(); k++)
		  creator.Add (face.vertices[k], GetFaceDofs(i));
	      }
	    
	    for (int i = 0; i < ni; i++)
              for (auto v : ma->GetElement(i).Vertices())
                creator.Add (v, GetElementDofs(i));

	    break; 
	    
	  case 8: // V + E + FI
	    if (creator.GetMode() == 1)
	      cout << " V + E + FI " << endl; 
		
	    for (int i = 0; i < nv; i++)
	      creator.Add (i, i);
		
	    for (int i = 0; i < ned; i++)
	      creator.Add (nv+i, GetEdgeDofs(i));
		
	    for (int i = 0; i < nfa; i++)
	      creator.Add(nv+ned+i, GetFaceDofs(i));
		
	    for (int i = 0; i < ni; i++)
	      for (int f : ma->GetElement(i).Faces())
		creator.Add (nv+ned+f, GetElementDofs(i));                  

	    break;


	  case 9: // V + EF + I
	    {
	      if (creator.GetMode() == 1)
                {
                  cout << " V + EF + I " << endl; 
		  creator.SetSize(nv+ned+ni);
                }
              else
                {
                  for (int i = 0; i < nv; i++)
                    creator.Add (i, i);
		  
                  for (int i = 0; i < ned; i++)
                    creator.Add (nv+i, GetEdgeDofs(i));
		  
                 
                  Array<int> f2ed;
                  for (int i = 0; i < nfa; i++)
                    {
                      /*
                        Ng_Node<2> face = ma->GetNode<2> (i);
                        for (int k = 0; k < face.edges.Size(); k++)
                        creator.Add (face.edges[k], GetFaceDofs(i));
                      */
                      
                      ma->GetFaceEdges (i, f2ed);
                      for (int k = 0; k < f2ed.Size(); k++)
                        creator.Add (nv+f2ed[k], GetFaceDofs(i));
                    }
                  
                  for (int i = 0; i < ni; i++)
                    creator.Add (nv+ned+i, GetElementDofs(i));
                }
	      break; 
	    }


	  case 10: 
	    if (creator.GetMode() == 1)
	      cout << " V + E + F +I (Cluster)  " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add(ma->GetClusterRepVertex(i), i);

	    for (int i = 0; i < ned; i++)
	      creator.Add (ma->GetClusterRepEdge(i), GetEdgeDofs(i));

	    for (int i = 0; i < nfa; i++)
	      creator.Add(ma->GetClusterRepFace(i), GetFaceDofs(i));

	    for (int i = 0; i < ni; i++)
	      creator.Add (ma->GetClusterRepElement(i), GetElementDofs(i));

	    break; 

	  case 11: 
	    if (creator.GetMode() == 1)
	      cout << " 2d VEI " << endl; 

	    for (int i = 0; i < nv; i++)
		creator.Add (i, i);
		
	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		for (int k = 0; k < 2; k++)
		  creator.Add (edge.vertices[k], GetEdgeDofs(i));
	      }

	    for (int i = 0; i < ni; i++)
	      {
		const Ngs_Element & ngel = ma->GetElement(i);
		for (int j = 0; j < ngel.vertices.Size(); j++)
		  creator.Add (ngel.vertices[j], GetElementDofs(i));
	      }

	    break;


	  case 12: 
	    if (creator.GetMode() == 1)
	      cout << " VEFI Cluster " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add(ma->GetClusterRepVertex(i), i);

	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		int rep[2];
		for (int k = 0; k < 2; k++)
		  rep[k] = ma->GetClusterRepVertex(edge.vertices[k]);
		
		creator.Add (rep[0], GetEdgeDofs(i));
		if (rep[0] != rep[1])
		  creator.Add (rep[1], GetEdgeDofs(i));
	      }

	    for (int i = 0; i < nfa; i++)
	      {
		Ng_Node<2> face = ma->GetNode<2> (i);
		int rep[4];
		
		for (int k = 0; k < face.vertices.Size(); k++)
		  {
		    rep[k] = ma->GetClusterRepVertex(face.vertices[k]);
		    
		    bool ok = true;
		    for (int j = 0; j < k; j++)
		      if (rep[j] == rep[k]) ok = false;
		    if (ok) creator.Add (rep[k], GetFaceDofs(i));
		  }
	      }

	    for (int i = 0; i < ni; i++)
	      {
		Ngs_Element ngel = ma->GetElement (i);
		int rep[8];
		      
		for (int k = 0; k < ngel.vertices.Size(); k++)
		  {
		    rep[k] = ma->GetClusterRepVertex(ngel.vertices[k]);
			
		    bool ok = true;
		    for (int j = 0; j < k; j++)
		      if (rep[j] == rep[k]) ok = false;
		    if (ok) creator.Add (rep[k], GetElementDofs(i));
		  }
	      }
	    break; 

	  case 20: // VE + EF + I
	    {
	      if (creator.GetMode() == 1)
		cout << "VE + EF + I" << endl;
		  
	      for (int i = 0; i < nv; i++)
		creator.Add (i, i);
		  
	      for (int i = 0; i < ned; i++)
		{
		  creator.Add (nv+i, GetEdgeDofs(i));
		  Ng_Node<1> edge = ma->GetNode<1> (i);
		  for (int k = 0; k < 2; k++)
		    creator.Add (edge.vertices[k], GetEdgeDofs(i));
		}
	      
	      Array<int> f2ed;
	      for (int i = 0; i < nfa; i++)
		{
		  ma->GetFaceEdges (i, f2ed);
		  for (int k = 0; k < f2ed.Size(); k++)
		    creator.Add (nv+f2ed[k], GetFaceDofs(i));
		}
	      
	      for (int i = 0; i < ni; i++)
		creator.Add (nv+ned+i, GetElementDofs(i));
	      break;
	    }
	    
	  case 21: // E + F 
	  {
	    if (creator.GetMode() == 1)
	      cout << "Helmholtz" << endl;
		
	    int ds_order = int(precflags.GetNumFlag ("ds_order", 1));
	    cout << "ds_order = " << ds_order << endl;
		
	    for (int i = 0; i < ned; i++)
	      {
		int first = first_edge_dof[i] + ds_order - 1;
		int ndof = first_edge_dof[i+1]-first;
		for (int j = 0; j < ndof; j++)
		  creator.Add (i, first+j);
	      }
	    for (int i = 0; i < nfa; i++)
	      creator.Add(ned+i, GetFaceDofs(i));
		
	    break; 
	  }

	  case 51: 
	  //for BDDC: we have only the condensated (after subassembling) dofs, 
	  //and build patches around each vertex Vertices + Edges
	    
	    if (creator.GetMode() == 1)
	      cout << "BDDC-Edges-around-Vertex-Block" << endl;

	    for (int i = 0; i < nv; i++)
	      creator.Add (i, i);
	    
	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		IntRange edgedofs = GetEdgeDofs(i);

		for (int k = 0; k < 2; k++)
		  for (int l = edgedofs.First(); l < edgedofs.Next(); l++)
		    if (ctofdof[l] == WIREBASKET_DOF)
		      creator.Add (edge.vertices[k], l);
	      }
	    
	    break; 	    
	  }
      }
    return creator.GetTable();
  }

    


  Array<int> * 
  H1HighOrderFESpace :: CreateDirectSolverClusters (const Flags & flags) const
  {
    if (flags.GetDefineFlag("subassembled"))
    {
	cout << "creating bddc-coarse grid(vertices)" << endl;
	Array<int> & clusters = *new Array<int> (GetNDof());
	clusters = 0;
	int nv = ma->GetNV();
	for (int i = 0; i < nv; i++)
	  if (!IsDirichletVertex(i))
	    clusters[i] = 1;		
	return &clusters;	
    }
    
    if (flags.NumFlagDefined ("ds_order"))
      {
	int ds_order = int (flags.GetNumFlag ("ds_order", 1));

	Array<int> & clusters = *new Array<int> (GetNDof());
	clusters = 0;
	
	int ned = ma->GetNEdges();
	int nv = ma->GetNV();

	// for (int i = 0; i < nv; i++)
        //   clusters[i] = 1;
        clusters.Range(0,nv) = 1;

	for (int i = 0; i < ned; i++)
	  {
	    int first = first_edge_dof[i];
	    int next = first_edge_dof[i+1];
	    for (int j = 0; (j+2 <= ds_order) && (first+j < next) ; j++)
	      clusters[first+j] = 1;
	  }

	/*
	  int nfa = ma->GetNFaces();
	  for (int i = 0; i < nfa; i++)
	  {
	  int first = first_face_dof[i];
	  int next = first_face_dof[i+1];
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
	*/

	return &clusters;
      }




    // return 0;

    // int nv = ma->GetNV();
    // int nd = GetNDof();
    int ne = ma->GetNE();
    Array<int> & clusters = *new Array<int> (GetNDof());
    clusters = 0;

    // all vertices in global space
    /*
      for (int i = 0; i < nv; i++)
      clusters[i] = 1;

      // all edges in direct solver cluster

      for (i = first_edge_dof[0]; i < first_edge_dof.Last(); i++)
      clusters[i] = 1;
      return &clusters;

    */
   
    // All Vertical Edges in one Cluster for Hex and Prism (-> 2d Problems !) 
    Array<int> ednums,fnums;

    //Array<int> & clusters = *new Array<int> (nd);
    //clusters = 0;


    
    for (int i = 0; i < ne; i++)
      {
	if (ma->GetElType(i) == ET_PRISM)
	  {
	    ma->GetElEdges (i, ednums);
	    for (int j = 6; j < 9; j++)  //vertical Edges 
              clusters[GetEdgeDofs(ednums[j])] = 2;
	    
	    ma->GetElFaces(i,fnums); // vertical faces 
	    for (int j =2;j<5;j++) 
	      {
                clusters[GetFaceDofs(fnums[j])] = 0;
                /*
		int first = first_face_dof[fnums[j]]; 
		int next = first_face_dof[fnums[j]+1]; 
		
		for (k=first; k < next; k++) 
		  clusters[k]=0; 
                */

		//INT<2> p = order_face[fnums[j]];
		//for(k=first + 2*(p[0]+1)*(p[1]+1);k<next;k++)
		//  clusters[k]=3;  
	      }
	  }

	else if (ma->GetElType(i) == ET_HEX)  
	  {
	    ma->GetElEdges (i, ednums);
	    for (int j = 8; j < 12; j++) //vertical edges
              clusters[GetEdgeDofs(ednums[j])] = 2;

	    ma->GetElFaces(i,fnums); // vertical faces 
	    for (int j = 2; j < 6; j++) 
              clusters[GetFaceDofs(fnums[j])] = 3;
	  } 
      }

   
    for (int i =0; directsolverclustered.Size() > 0 && i<ne; i++)
      {
	if(directsolverclustered[ma->GetElIndex(i)])
	  {
	    GetDofNrs(i,ednums);
	    for (int k = 0; k<ednums.Size(); k++)
	      {
		clusters[ednums[k]] = 4;
	      }
	  }
      }

   
    clusters[adddirectsolverdofs] = 5;
    /*
    for (int i =0; i< adddirectsolverdofs.Size(); i++)
      {
	clusters[adddirectsolverdofs[i]] = 5;
      }
    */
    const int stdoffset = 6;

    /*
      (*testout) << "directvertexclusters" << directvertexclusters << endl;
      (*testout) << "directedgeclusters" << directedgeclusters << endl;
      (*testout) << "directfaceclusters" << directfaceclusters << endl;
      (*testout) << "directelementclusters" << directelementclusters << endl;
    */

    for (int i = 0; i<directvertexclusters.Size(); i++)
      if(directvertexclusters[i] >= 0)
	clusters[i] = directvertexclusters[i] + stdoffset;

    for (int i = 0; i<directedgeclusters.Size(); i++)
      if(directedgeclusters[i] >= 0)
	for(int j = first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
	  clusters[j] = directedgeclusters[i] + stdoffset;

    for (int i = 0; i<directfaceclusters.Size(); i++)
      if(directfaceclusters[i] >= 0)
	for(int j = first_face_dof[i]; j<first_face_dof[i+1]; j++)
	  clusters[j] = directfaceclusters[i] + stdoffset;
	  
    for (int i = 0; i<directelementclusters.Size(); i++)
      if(directelementclusters[i] >= 0)
	for(int j = first_element_dof[i]; j<first_element_dof[i+1]; j++)
	  clusters[j] = directelementclusters[i] + stdoffset;


    //    (*testout) << "clusters " << clusters << endl;
    
    bool nonzero = false;
    for (int i = 0; !nonzero && i < clusters.Size(); i++)
      if (clusters[i]) nonzero = true;
    if (!nonzero)
      {
	delete &clusters;
	return 0;
      }

    return &clusters;
  }

  template<>
  shared_ptr<FESpace> RegisterFESpace<H1HighOrderFESpace> :: 
  Create (shared_ptr<MeshAccess> ma, const Flags & flags)
  {
    // we will check the -periodic flag here
    return make_shared<H1HighOrderFESpace> (ma, flags);
  }

  static RegisterFESpace<H1HighOrderFESpace> init ("h1ho");
}
 
