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
#include <solve.hpp>
#include <parallelngs.hpp>


namespace ngfem
{
#include "../fem/h1hofefo.hpp"
}


using namespace ngmg; 

namespace ngcomp
{
  using namespace ngcomp;
  using namespace ngparallel;

  H1HighOrderFESpace ::  
  H1HighOrderFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name = "H1HighOrderFESpace(h1ho)";
    // define h1ho flags
    DefineDefineFlag("h1ho");
    DefineNumFlag("augmented");
    DefineDefineFlag("plate");
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
    DefineDefineFlag("minext");
    DefineDefineFlag("optext");
    DefineDefineFlag("fast");
    DefineDefineFlag("fastsz");
    DefineDefineFlag("print");
    DefineDefineFlag("noprint");
    if (parseflags) ParseFlags(flags);
    
    augmented = int (flags.GetNumFlag ("augmented", 0));
    plate = int (flags.GetDefineFlag ("plate"));
    print = (flags.GetDefineFlag("print")); 
  
    // Variable order space: 
    //      in case of (var_order && order) or (relorder) 
    var_order = flags.GetDefineFlag("variableorder");  
    order =  int (flags.GetNumFlag ("order",1)); 
    
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
          
    Flags loflags;
    loflags.SetFlag ("order", 1);
    loflags.SetFlag ("dim", dimension);
    if (iscomplex) loflags.SetFlag ("complex");

#ifndef PARALLEL
    low_order_space = new NodalFESpace (ma, loflags);
#else
    low_order_space = new ParallelNodalFESpace (ma, loflags);
#endif
          
    minext = flags.GetDefineFlag ("minext");
    optext = flags.GetDefineFlag ("optext");
    fast_pfem = flags.GetDefineFlag ("fast");
    

    static ConstantCoefficientFunction one(1);
    if (ma.GetDimension() == 2)
      {
	evaluator = new MassIntegrator<2> (&one);
	boundary_evaluator = 0;
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

    prol = new LinearProlongation(*this);
  }


  H1HighOrderFESpace :: ~H1HighOrderFESpace ()
  {
    ;
  }

  FESpace * H1HighOrderFESpace :: 
  Create (const MeshAccess & ma, const Flags & flags)
  {
    return new H1HighOrderFESpace (ma, flags, true); // true: parse flags
  }
    
  void H1HighOrderFESpace :: Update(LocalHeap & lh)
  {
    int maxorder = -1; 
    int minorder = 99; 

    try
      {
	if (low_order_space)
	  low_order_space -> Update(lh);
      }
    catch (exception & e)
      {
	throw Exception (e.what() + 
			 string ("\nthrown by H1HoFESpace, UpdateLowOrderSpace "));
      }
    catch (Exception & e)
      {
	e.Append (string ("\nthrown by allocate matrix ") +
		  string ("\nthrown by H1HoFESpace, UpdateLowOrderSpace "));
	throw;
      }

    
    const int dim = ma.GetDimension();
    nv = ma.GetNV();
    ned = ma.GetNEdges();
    nfa = (dim == 2) ? 0 : ma.GetNFaces();
    // nfa = ma.GetNFaces();
    nel = ma.GetNE();
    
    order_edge.SetSize (ned);
    order_face.SetSize (nfa);
    order_inner.SetSize (nel);
    order_avertex.SetSize (nv); 
    fine_edge.SetSize(ned); 
    fine_face.SetSize(nfa); 

    fine_edge = 0;
    fine_face = 0; 
    int p = var_order ?  1 : order; 
    
    order_edge = p; 
    order_face = INT<2>(p,p);
    order_inner = INT<3>(p,p,p); 
    order_avertex = p; 
	
    Array<int> eledges, elfaces, vnums;
    
    for (int i = 0; i < nel; i++)
      {	
	ELEMENT_TYPE eltype=ma.GetElType(i); 
	const FACE * faces = ElementTopology::GetFaces (eltype);
	const EDGE * edges = ElementTopology::GetEdges (eltype);
	const POINT3D * points = ElementTopology :: GetVertices (eltype);
	ma.GetElVertices (i, vnums);
	ma.GetElEdges (i, eledges);		
	

	if(dim==3) 
	  {
	    ma.GetElFaces(i, elfaces);
	    for (int j=0;j<elfaces.Size();j++) fine_face[elfaces[j]] = 1; 
	  }
	for (int j=0;j<eledges.Size();j++) fine_edge[eledges[j]] = 1; 
	
	if(!var_order) continue;  
	
	INT<3> el_orders = ma.GetElOrders(i); 
	for(int l=0;l<3;l++) el_orders[l] += rel_order; 

	for(int l=0;l<3;l++) 
	  maxorder = int(max2(el_orders[l],maxorder)); 
	for(int l=0;l<3;l++) 
	  minorder = int(min2(el_orders[l],minorder)); 
	

	for(int j=0;j<dim;j++)
	  order_inner[i][j] = int(max2(order_inner[i][j],el_orders[j]));
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
      if(uniform_order_inner > -1) 
      order_inner = INT<3>(uniform_order_inner,uniform_order_inner,uniform_order_inner);
      if(uniform_order_face > -1) 
      order_face = INT<2>(uniform_order_face,uniform_order_face); 
      if(uniform_order_edge > -1) 
      order_edge = uniform_order_edge; 
    */ 

 
 

    /*
    // skip dofs on Dirichlet boundary hack [JS]
    int nsel = ma.GetNE();
    for (i = 0; i < nsel; i++)
    {
    ma.GetSElEdges (i, eledges);
    int face = ma.GetSElFace (i);
	  
    for (j = 0; j < eledges.Size(); j++) 
    order_edge[eledges[j]] = 1;  
    order_face[face] = 1;
    }
    */

    /* 
       if (ma.GetDimension() == 2 && uniform_order_trig != -1 && uniform_order_quad != -1)
       {
       for (int i = 0; i < nel; i++)
       {
       if (ma.GetElType(i) == ET_TRIG)
       order_inner = INT<3> (uniform_order_trig, uniform_order_trig, uniform_order_trig);
       else
       order_inner = INT<3> (uniform_order_quad, uniform_order_quad, uniform_order_quad);
       }
       }
    */ 
    
    if(uniform_order_inner > -1) 
      order_inner = INT<3>(uniform_order_inner,uniform_order_inner,uniform_order_inner);
    if(uniform_order_face > -1 && dim==3) 
      order_face = INT<2>(uniform_order_face,uniform_order_face); 
    if(uniform_order_edge > -1) 
      order_edge = uniform_order_edge; 
    

    for(int i=0;i<ned;i++) if(!fine_edge[i]) {order_edge[i] = 1;}
    if(dim == 3) for(int i=0;i<nfa;i++) if(!fine_face[i]) order_face[i] = INT<2> (1,1); 
    if(dim==2) 
      { 
	nfa = 0; 
	fine_face.SetSize(0); 
	order_face.SetSize(0);
      } 

    if(print) 
      {
	*testout << " H1HoFESpace order " << order << " , var_order " << var_order << " , relorder " << rel_order << endl;  
	(*testout) << "fine_edge (h1): " << fine_edge << endl;
	(*testout) << "fine_face (h1): " << fine_face << endl;
	
	
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

#ifdef PARALLEL
    try
      {
	UpdateParallelDofs();
      }
    catch (exception & e)
      {
	throw Exception (e.what() + 
			 string ("\nthrown by H1HoFESpace, UpdateParallelDofs "));
      }
    catch (Exception & e)
      {
	e.Append (string ("\nthrown by allocate matrix ") +
		  string ("\nthrown by H1HoFESpace, UpdateParallelDofs "));
	throw;
      }

#endif

  }


  void H1HighOrderFESpace :: UpdateDofTables ()
  {
    int dim = ma.GetDimension();

    nv = ma.GetNV();
    ned = ma.GetNEdges();
    nfa = (dim == 2) ? 0 : ma.GetNFaces();
    //     nfa = ma.GetNFaces();
    nel = ma.GetNE();

    ndof = nv;

    first_edge_dof.SetSize (ned+1);
    for (int i = 0; i < ned; i++)
      {
	first_edge_dof[i] = ndof;
	if(order_edge[i]>1)
	  ndof += order_edge[i] - 1;
      }
    first_edge_dof[ned] = ndof;
         
    first_face_dof.SetSize (nfa+1);
    Array<int> fapnums;
    for (int i = 0; i < nfa; i++)
      {
	first_face_dof[i] = ndof;
	ma.GetFacePNums (i, fapnums);
	INT<2> p = order_face[i];
	switch(fapnums.Size())
	  {
	  case 3: 
	    if(p[0]>2)
	      ndof += (p[0]-1)*(p[0]-2)/2;
	    break;
	  case 4:
	    if(p[0] > 1 && p[1]>1)
	      ndof += (p[0]-1)*(p[1]-1);
	    break; 
	  }
      }
    first_face_dof[nfa] = ndof;
 
    first_element_dof.SetSize(nel+1);
    for (int i = 0; i < nel; i++)
      {
	first_element_dof[i] = ndof;
	INT<3> p = order_inner[i];	
	switch (ma.GetElType(i))
	  {
	  case ET_TRIG:
	    if(p[0]>2)
	      ndof += (p[0]-1)*(p[0]-2)/2;
	    break;
	  case ET_QUAD:
	    if(p[0]>1 && p[1]>1)
	      ndof += (p[0]-1)*(p[1]-1);
	    break;
	  case ET_TET:
	    if(p[0] > 3)
	      ndof += (p[0]-1)*(p[0]-2)*(p[0]-3)/6;
	    break;
	  case ET_PRISM:
	    if(p[0]>2 && p[2]>1)
	      ndof += (p[0]-1)*(p[0]-2)*(p[2]-1)/2;
	    break;
	  case ET_PYRAMID:
	    if(p[0]>2)
	      ndof += (p[0]-1)*(p[0]-2)*(2*p[0]-3)/6;
	    break;
	  case ET_HEX:
	    if(p[0]>1 && p[1] > 1 && p[2]>1) 
	      ndof += (p[0]-1)*(p[1]-1)*(p[2]-1);
	    break;
	  }
      } 
    first_element_dof[nel] = ndof;
   
    /*
      (*testout) << "augmented " << augmented << endl; 
   
      (*testout) << "h1 first edge = " << first_edge_dof << endl;
      (*testout) << "h1 first face = " << first_face_dof << endl;
      (*testout) << "h1 first inner = " << first_element_dof << endl;
    
      (*testout) << "order_edge H1 = " << order_edge << endl;
      (*testout) << "order_face H1 = " << order_face << endl;
      (*testout) << "order_inner H1 = " << order_inner << endl;
      (*testout) << " NDOF H1" << ndof << endl; 
    
      (*testout) << specialelements.Size() << endl;
    */

    if (dirichlet_boundaries.Size())
      {
	dirichlet_vertex.SetSize (ma.GetNV());
	dirichlet_edge.SetSize (ma.GetNEdges());
	dirichlet_face.SetSize (ma.GetNFaces());
	
	Array<int> vnums, ednums;
	int fanum;

	dirichlet_vertex = false;
	dirichlet_edge = false;
	dirichlet_face = false;

	for (int i = 0; i < ma.GetNSE(); i++)
	  {
	    int ind = ma.GetSElIndex (i);
	    if (dirichlet_boundaries.Test(ind))
	      {
		ma.GetSElVertices (i, vnums);
		ma.GetSElEdges (i, ednums);
		fanum = ma.GetSElFace (i);

		for (int j = 0; j < vnums.Size(); j++)
		  dirichlet_vertex[vnums[j]] = true;
		for (int j = 0; j < ednums.Size(); j++)
		  dirichlet_edge[ednums[j]] = true;
		dirichlet_face[fanum] = true;
	      }
	  }

	(*testout) << "Dirichlet_vertex = " << endl << dirichlet_vertex << endl;
	(*testout) << "Dirichlet_edge = " << endl << dirichlet_edge << endl;
	(*testout) << "Dirichlet_face = " << endl << dirichlet_face << endl;
      }



    lodofs_per_node[0] = 1;
    for (int i = 1; i < 4; i++) lodofs_per_node[i] = 0;
    first_lodof[0] = 0;
    for (int i = 1; i < 5; i++) first_lodof[i] = ma.GetNV();
    first_hodofs[0].SetSize(0);
    first_hodofs[1] = first_edge_dof;

    if (dim == 3)
      {
        first_hodofs[2] = first_face_dof;
        first_hodofs[3] = first_element_dof;
      }
    else
      first_hodofs[2] = first_element_dof;
    
    while (ma.GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;
    prol->Update();
  }


  void H1HighOrderFESpace :: PrintReport (ostream & ost)
  {
    FESpace::PrintReport (ost);

    //     (*testout) << "first edge = " << first_edge_dof << endl;
    //     (*testout) << "first face = " << first_face_dof << endl;
    //     (*testout) << "first inner = " << first_element_dof << endl;
  }


  


  const FiniteElement & H1HighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    H1HighOrderFiniteElement<2> * hofe2d = 0;
    H1HighOrderFiniteElement<3> * hofe3d = 0;

    /*
      // order_inner etc. must all be the same ...
    if (!var_order && ma.GetElType(elnr) == ET_TRIG && order <= 6)
      {
        H1HighOrderFiniteElementFO<2> * hofe2d = 0;
        switch (order)
          {
          case 1: hofe2d = new (lh.Alloc (sizeof(H1HighOrderFEFO<ET_TRIG,1>)))  H1HighOrderFEFO<ET_TRIG,1> (); break;
          case 2: hofe2d = new (lh.Alloc (sizeof(H1HighOrderFEFO<ET_TRIG,2>)))  H1HighOrderFEFO<ET_TRIG,2> (); break;
          case 3: hofe2d = new (lh.Alloc (sizeof(H1HighOrderFEFO<ET_TRIG,3>)))  H1HighOrderFEFO<ET_TRIG,3> (); break;
          case 4: hofe2d = new (lh.Alloc (sizeof(H1HighOrderFEFO<ET_TRIG,4>)))  H1HighOrderFEFO<ET_TRIG,4> (); break;
          case 5: hofe2d = new (lh.Alloc (sizeof(H1HighOrderFEFO<ET_TRIG,5>)))  H1HighOrderFEFO<ET_TRIG,5> (); break;
          case 6: hofe2d = new (lh.Alloc (sizeof(H1HighOrderFEFO<ET_TRIG,6>)))  H1HighOrderFEFO<ET_TRIG,6> (); break;
          }

    
        ArrayMem<int,12> vnums; 
        ma.GetElVertices(elnr, vnums);
        hofe2d -> SetVertexNumbers (vnums);
        return *hofe2d;
      }
    */

    try
      {
        typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
    
        switch (ma.GetElType(elnr))
          {
          case ET_TET:
            { 
              hofe3d = new (lh.Alloc (sizeof(H1HighOrderFE<ET_TET>)))  H1HighOrderFE<ET_TET> ();
              break;
            }
          case ET_PYRAMID:
            {
              hofe3d = new (lh.Alloc (sizeof(H1HighOrderFE<ET_PYRAMID>)))  H1HighOrderFE<ET_PYRAMID> ();
              break;
            }
          case ET_PRISM:
            {
              hofe3d = new (lh.Alloc (sizeof(H1HighOrderFE<ET_PRISM>)))  H1HighOrderFE<ET_PRISM> ();
              break;
            }
          case ET_HEX:
            {
              hofe3d = new (lh.Alloc (sizeof(H1HighOrderFE<ET_HEX>)))  H1HighOrderFE<ET_HEX> ();
              break;
            }
          case ET_TRIG:
            { 
              hofe2d = new (lh.Alloc (sizeof(H1HighOrderFE<ET_TRIG>)))  H1HighOrderFE<ET_TRIG> ();
              break;
            }
          case ET_QUAD:
            {
              hofe2d = new (lh.Alloc (sizeof(H1HighOrderFE<ET_QUAD>)))  H1HighOrderFE<ET_QUAD> ();
              break;
            }
          default:
            {
              throw Exception ("GetFE not supported for element");
            }
          }
    
        if (!hofe2d && !hofe3d)
          {
            stringstream str;
            str << "H1HighOrderFESpace " << GetClassName() 
                << ", undefined eltype " 
                << ElementTopology::GetElementName(ma.GetElType(elnr))
                << ", order = " << order << endl;
            throw Exception (str.str());
          }


    
        ArrayMem<int,12> vnums; // calls GetElPNums -> max 12 for PRISM12
        ArrayMem<int, 12> ednums, order_ed;
        ArrayMem<int, 6> fanums;
        ArrayMem<INT<2>, 6> order_fa;


        ma.GetElVertices(elnr, vnums);

        ma.GetElEdges(elnr, ednums);

        order_ed.SetSize (ednums.Size());

        for (int j = 0; j < ednums.Size(); j++)
          order_ed[j] = order_edge[ednums[j]];

        ArrayMem<int, 8> order_vert(vnums.Size());
        for (int j = 0; j < vnums.Size(); j++)
          order_vert[j] = order_avertex[vnums[j]];

#ifdef PARALLEL
        if ( ntasks > 1 )
          for ( int i = 0; i < vnums.Size(); i++ )
            vnums[i] = parallelma->GetDistantPNum(0, vnums[i]);
#endif


        if (ma.GetDimension() == 2)
          {
            hofe2d -> SetVertexNumbers (vnums);

            hofe2d -> SetOrderEdge (order_ed);
            // hofe2d -> SetOrderInner (order_inner[elnr]); // old style
            
            INT<2> p(order_inner[elnr][0], order_inner[elnr][1]);
            FlatArray<INT<2> > of(1, &p);
            hofe2d -> SetOrderFace (of);


            hofe2d -> ComputeNDof();

            return *hofe2d;
          }

        else

          {
            hofe3d -> SetVertexNumbers (vnums);

            hofe3d -> SetOrderEdge (order_ed);
    
            ma.GetElFaces(elnr, fanums);
            order_fa.SetSize (fanums.Size());
            for (int j = 0; j < fanums.Size(); j++)
              order_fa[j] = order_face[fanums[j]];
	
            hofe3d -> SetOrderFace (order_fa);
            hofe3d -> SetOrderCell (order_inner[elnr]);
            hofe3d -> ComputeNDof();

            return *hofe3d;
          }
      }

    catch (Exception & e)
      {
        e.Append ("in H1HoFESpace::GetElement");
        e.Append ("\n");
        throw;
      }
  }
 
 


  const FiniteElement & H1HighOrderFESpace :: GetSFE (int elnr, LocalHeap & lh) const
  {
    /*
      FiniteElement * fe = 0;

      switch (ma.GetSElType(elnr))
      {
      case ET_SEGM:
      fe = segm; break;
      case ET_TRIG:
      fe = trig; break;
      case ET_QUAD:
      fe = quad; break;
      default:
      fe = 0;
      }

      H1HighOrderFiniteElement * hofe =
      dynamic_cast<H1HighOrderFiniteElement*> (fe);
    */

    H1HighOrderFiniteElement<1> * hofe1d;
    H1HighOrderFiniteElement<2> * hofe2d;

    {
      typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
	
      switch (ma.GetSElType(elnr))
        {
        case ET_TRIG:
          {
            hofe2d = new (lh.Alloc (sizeof(H1HighOrderFE<ET_TRIG>)))  H1HighOrderFE<ET_TRIG> ();
            break;
          }
        case ET_QUAD:
          {
            hofe2d = new (lh.Alloc (sizeof(H1HighOrderFE<ET_QUAD>)))  H1HighOrderFE<ET_QUAD> ();
            break;
          }
        case ET_SEGM:
          {
            hofe1d = new (lh.Alloc (sizeof(H1HighOrderFE<ET_SEGM>)))  H1HighOrderFE<ET_SEGM> ();
            break;
          }
        default:
          {
            throw Exception ("GetFE not supported for element");
          }
        }
    }
      
    if (!hofe1d && !hofe2d)
      {
        stringstream str;
        str << "FESpace " << GetClassName() 
            << ", undefined surface eltype " << ma.GetSElType(elnr) 
            << ", order = " << order << endl;
        throw Exception (str.str());
      }

    
    ArrayMem<int, 8> vnums, order_vert(8);
    ArrayMem<int, 4> ednums, order_ed;
    
    ma.GetSElVertices(elnr, vnums);
    ma.GetSElEdges(elnr, ednums);
    
    order_ed.SetSize (ednums.Size());
    
    for (int j = 0; j < ednums.Size(); j++)
      order_ed[j] = order_edge[ednums[j]];
    
    for (int j = 0; j < 8; j++)
      order_vert[j] = order;
    
    
    if (ma.GetDimension() == 2)
      {
        // hofe1d -> SetOrderInner (order_ed[0]);  // old style
        hofe1d -> SetOrderEdge (order_ed);

        hofe1d -> ComputeNDof();
        
#ifdef PARALLEL
        if ( ntasks > 1 )
          for ( int i = 0; i < vnums.Size(); i++ )
            vnums[i] = parallelma->GetDistantPNum(0, vnums[i]);
#endif 
        hofe1d -> SetVertexNumbers (vnums);  
        
        return *hofe1d;
      }
    else
      {
        hofe2d -> SetOrderEdge (order_ed);
        INT<2> p = order_face[ma.GetSElFace(elnr)];
        // hofe2d -> SetOrderInner (INT<3> (p[0], p[1], 0));  // old style

        FlatArray<INT<2> > of(1, &p);
        hofe2d -> SetOrderFace (of);

        hofe2d  -> ComputeNDof();
#ifdef PARALLEL
        if ( ntasks > 1 )
          for ( int i = 0; i < vnums.Size(); i++ )
            vnums[i] = parallelma->GetDistantPNum(0, vnums[i]);
#endif 
        hofe2d -> SetVertexNumbers (vnums);  
        
        return *hofe2d;
      }
  }
 






  int H1HighOrderFESpace :: GetNDof () const
  {
    return ndof;
  }



  int H1HighOrderFESpace :: GetNDofLevel (int alevel) const
  {
    return ndlevel[alevel];
  }



  void H1HighOrderFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    // FESpace :: GetDofNrs(elnr, dnums);
    // return;

    ArrayMem<int,12> vnums, ednums, fanums; 
    int i, j;
    int first, next; 

    ma.GetElVertices (elnr, vnums);
    ma.GetElEdges (elnr, ednums);
    if (ma.GetDimension() == 3)
      ma.GetElFaces (elnr, fanums);
    else 
      fanums.SetSize(0);
    dnums.SetSize(0); 

    
    if(order < 1 && !var_order) 
      throw Exception(" H1HighOrderFESpace :: GetDofNrs() order < 1 "); 

    for (i = 0; i < vnums.Size(); i++)
      dnums.Append (vnums[i]);

    //(*testout) << "after verts dnums.Size() " << dnums.Size() << endl;

    for (i = 0; i < ednums.Size(); i++)
      {
        first = first_edge_dof[ednums[i]];
        next = first_edge_dof[ednums[i]+1];
        for (j = first; j < next; j++)
          dnums.Append (j);
      }

    //(*testout) << "after edges dnums.Size() " << dnums.Size() << endl;

    for (i = 0; i < fanums.Size(); i++)
      {
        first = first_face_dof[fanums[i]];
        next = first_face_dof[fanums[i]+1];
	  
        for (j = first; j < next; j++)
          dnums.Append (j);
      }

    //(*testout) << "after faces dnums.Size() " << dnums.Size() << endl;

    first = first_element_dof[elnr];
    next = first_element_dof[elnr+1];
 
    
    for (j = first; j < next; j++)
      dnums.Append (j);
     
    //(*testout) << "after inner dnums.Size() " << dnums.Size() << endl;

    if (!DefinedOn (ma.GetElIndex (elnr)))
      dnums = -1;

    // cout << "orig dnums = " << dnums << endl;
  }


  void  H1HighOrderFESpace :: GetExternalDofNrs (int elnr, Array<int> & dnums) const
  {
    if (!eliminate_internal) 
      {
        GetDofNrs (elnr, dnums);
        return;
      }

    ArrayMem<int,12> vnums, ednums, fanums; 
    int i, j;
    int first,next; 

    ma.GetElVertices (elnr, vnums);
    ma.GetElEdges (elnr, ednums);
    if (ma.GetDimension() == 3)
      ma.GetElFaces (elnr, fanums);
    else 
      fanums.SetSize(0); 
    dnums.SetSize(0);

    for (i = 0; i < vnums.Size(); i++)
      dnums.Append (vnums[i]);

    if(order < 1) 
      throw Exception(" H1HighOrderFESpace :: GetDofNrs() order < 1 "); 

    for (i = 0; i < ednums.Size(); i++)
      {
        first = first_edge_dof[ednums[i]];
        next = first_edge_dof[ednums[i]+1];
	
        for (j = first; j < next; j++)
          dnums.Append (j);
      }

    if (ma.GetDimension() == 3)
      {
        int first_face = 0;
        if (plate && vnums.Size() == 6) first_face = 2;
        for (i = first_face; i < fanums.Size(); i++)
          {
            first = first_face_dof[fanums[i]];
            next = first_face_dof[fanums[i]+1]; 
            for (j = first; j < next; j++)
              dnums.Append (j);
          }
      }
    if (!DefinedOn (ma.GetElIndex (elnr)))
      dnums = -1;
  }





  void H1HighOrderFESpace :: GetWireBasketDofNrs (int elnr, Array<int> & dnums) const
  {
    ArrayMem<int,12> vnums, ednums;

    dnums.SetSize(0);

    ma.GetElVertices (elnr, vnums);
    ma.GetElEdges (elnr, ednums);
    
    for (int i = 0; i < vnums.Size(); i++)
      dnums.Append (vnums[i]);

    for (int i = 0; i < ednums.Size(); i++)
      if (order_edge[ednums[i]] > 1)
        dnums.Append (first_edge_dof[ednums[i]]);
  }





  
  void H1HighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    dnums.Append (vnr);
  }
  
  
  void H1HighOrderFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    int first = first_edge_dof[ednr];
    int neddofs = first_edge_dof[ednr+1] - first;
    for (int j = 0; j < neddofs; j++)
      dnums.Append (first+j);
  }

  void H1HighOrderFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    if (ma.GetDimension() < 3) return;
    
    int first = first_face_dof[fanr];
    int nfadofs = first_face_dof[fanr+1] - first;
    for (int j = 0; j < nfadofs; j++)
      dnums.Append (first+j);
  }


  void H1HighOrderFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);

    int first = first_element_dof[elnr];
    int neldofs = first_element_dof[elnr+1] - first;
     
    for (int j = 0; j < neldofs; j++)
      dnums.Append (first+j);
  }
  
  void H1HighOrderFESpace :: 
  GetSDofNrs (int elnr, Array<int> & dnums) const
  {
    static bool getall = false;

    Array<int> vnums, ednums;
    int fanum;
    int i, j;

    ma.GetSElVertices (elnr, vnums);
    ma.GetSElEdges (elnr, ednums);
    if (ma.GetDimension() == 3)
      fanum = ma.GetSElFace (elnr);
    
    dnums.SetSize(0);
    for (i = 0; i < vnums.Size(); i++)
      dnums.Append (vnums[i]);
    
    Array<int> edge_start;

    for (i = 0; i < ednums.Size(); i++)
      {
        edge_start.Append(dnums.Size());
        int first = first_edge_dof[ednums[i]];
        int neddofs = first_edge_dof[ednums[i]+1] - first;
        for (j = 0; j < neddofs; j++)
          dnums.Append (first+j);
      }
    edge_start.Append(dnums.Size());

    int face_start  = dnums.Size();

    if (ma.GetDimension() == 3)
      {
        int first = first_face_dof[fanum];
        int nfadofs = first_face_dof[fanum+1] - first;
        for (j = 0; j < nfadofs; j++)
          dnums.Append (first+j);
      }

    if (!DefinedOnBoundary (ma.GetSElIndex (elnr)) && !getall)
      dnums = -1;

    if(defined_on_one_side_of_bounding_curve.Size() > 0 && !getall)
      {
        Array<bool> keep_dnum(dnums.Size());
        keep_dnum = true;

        Array<int> vnums;
        Array<int> neighbours;

        ma.GetSElVertices(elnr,vnums);
        for(int i=0; i<vnums.Size(); i++)
          {
            Array<int> auxn;
            ma.GetVertexSurfaceElements(vnums[i],auxn);
            for(int j=0; j<auxn.Size(); j++)
              if(auxn[j] != elnr)
                neighbours.Append(auxn[j]);
          }
	

        bool isfirst = true;
        for(int i=0; i<defined_on_one_side_of_bounding_curve.Size(); i++)
          {
            if(defined_on_one_side_of_bounding_curve[i][0] == ma.GetSElIndex (elnr))
              {
                if(isfirst)
                  {
                    keep_dnum = false;
                    isfirst = false;
                  }

                for(int j=0; j<neighbours.Size(); j++)
                  {
                    if(ma.GetSElIndex(neighbours[j]) == 
                       defined_on_one_side_of_bounding_curve[i][1])
                      {
                        //(*testout) << "sel " << elnr << " neighbour " << neighbours[j] << endl;
                        getall = true;
                        Array<int> neighbour_dnums;
                        GetSDofNrs(neighbours[j],neighbour_dnums);
                        getall = false;
			
                        for(int k=0; k<neighbour_dnums.Size(); k++)
                          {
                            if(neighbour_dnums[k] == -1)
                              continue;

                            int pos = dnums.Pos(neighbour_dnums[k]);

                            if(pos >= 0)
                              keep_dnum[pos] = true;
                          }
                      }
                  }		
              }
          }

        for(int i=0; i<dnums.Size(); i++)
          if(!keep_dnum[i])
            dnums[i] = -1;

        //if(!isfirst)
        //  (*testout) << "elnr " << elnr << ", keep_dnum " << keep_dnum << endl;
      }

  }
  

  
  void H1HighOrderFESpace :: RestrictToOneSideOfBoundingCurve(int index1, int index2)
  {
    defined_on_one_side_of_bounding_curve.Append(INT<2>(index1,index2));
  }
  void H1HighOrderFESpace :: DeleteOneSideRestrictions(void)
  {
    defined_on_one_side_of_bounding_curve.DeleteAll();
  }

  
  Table<int> * H1HighOrderFESpace :: 
  CreateSmoothingBlocks (const Flags & precflags ) const
  {
    int i, j, k, first;
    int SmoothingType = int(precflags.GetNumFlag("blocktype",0)); 
    

    Array<int> orient, ednums, fanums, vnums,f2ed; 
    


    int ni = nel, ncnt; 
    if (eliminate_internal) ni = 0; 
   
    // SmoothingTypes: 
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
    
    // default smoother
    
    if (SmoothingType == 0) 
      if(augmented) SmoothingType=3; 
      else SmoothingType = 4; 
    


    cout << " blocktype " << SmoothingType << endl; 
    // *testout << " h1ho-Smoother with order " << order << " level " << level << endl; 
    cout << " Use H1-Block Smoother:  "; 
    switch(SmoothingType) 
      {
      case 1 : 
        cout << " 2d V + E + I " << endl; 
        ncnt = nv + ned + ni;
        break;
      case 2 : 
        cout << " 2d VE + I " << endl; 
        ncnt = nv + ni;
        break; 
      case 3: 
        cout << " V + E + F + I " << endl; 
        ncnt = nv + ned + nfa + ni;
        break; 
      case 4:
        cout << " VE + F + I " << endl;
        ncnt = nv + nfa + ni;
        break; 
      case 5: 
        cout << " VE + FI " << endl; 
        ncnt = nv + nfa; 
        break; 
      case 6: 
        cout << " VEF + I " << endl; 
        ncnt = nv + ni; 
        break; 
      case 7: 
        cout << " VEFI " << endl; 
        ncnt = nv; 
        break; 
      case 8: 
        cout << " V + E + FI " << endl; 
        ncnt = nv + ned + nfa; 
        break; 
      case 9: 
        cout << " V + EF + I " << endl; 
        ncnt = nv + ned + ni; 
        break; 
      case 10: 
        cout << " V + E + F +I (Cluster)  " << endl; 
        ncnt = nv + ned + nfa; 
        break; 
      case 11 : 
        cout << " 2d VEI " << endl; 
        ncnt = nv + ned + ni;
        break;
      case 12: 
        cout << " VEFI Cluster " << endl; 
        ncnt = nv; 
        break; 

      case 20:
        cout << "Joachim's special" << endl;
        ncnt = nv+ned;
        break;

      default: 
        cout << " Error in H1HOFESpace :: CreateSmoothingBlocks no SmoothingType " << endl; 
        return 0; 
      }
    
    Array<int> cnt(ncnt); 
    cnt = 0; 
    int nvdof = 1; 
    if (augmented == 1) nvdof = 2; 
    if (augmented == 2) nvdof = order; //SZ: hier noch v_order(i) einsetzen 
    for (i=0; i<nv; i++)
      cnt[i] = nvdof; 
      
    int ii =0; 
    switch(SmoothingType)
      { 
      case 1: // 2d V + E + I 
        for (i = 0; i < ned; i++, ii++)
          if(fine_edge[i])
            cnt[nv+ii] = first_edge_dof[i+1]-first_edge_dof[i];
        for (i = 0; i < ni; i++)
          cnt[nv+ned+i] = first_element_dof[i+1]-first_element_dof[i];
        break;
      case 2: // 2d VE + I 
        for (i = 0; i < ned; i++ )
          if(fine_edge[i])
            {
              int v1, v2;
              ma.GetEdgePNums (i, v1, v2);
              cnt[v1] += first_edge_dof[i+1]-first_edge_dof[i];
              cnt[v2] += first_edge_dof[i+1]-first_edge_dof[i];
            }
        for (i = 0; i < ni; i++)
          cnt[nv+i] = first_element_dof[i+1]-first_element_dof[i];
        break; 
      case 3: // V + E + F + I 
        for (i = 0; i < ned; i++)
          cnt[nv+i] = first_edge_dof[i+1]-first_edge_dof[i];
        for (i = 0; i < nfa; i++)
          cnt[nv+ned+i] = first_face_dof[i+1]-first_face_dof[i];
        for (i = 0; i < ni; i++)
          cnt[nv+ned+nfa+i] = first_element_dof[i+1]-first_element_dof[i];
        break; 
      case 4: // VE + F + I 
        for (i = 0; i < ned; i++)
          {
            int v1, v2;
            int ndof = first_edge_dof[i+1]-first_edge_dof[i];
            ma.GetEdgePNums (i, v1, v2);
            cnt[v1] += ndof;
            cnt[v2] += ndof;
          }
        for (i = 0; i < nfa; i++)
          cnt[nv+i] = first_face_dof[i+1]-first_face_dof[i];
        for (i = 0; i < ni; i++)
          cnt[nv+nfa+i] = first_element_dof[i+1]-first_element_dof[i];
        break; 
      case 5: // VE + FI 
        for (i = 0; i < ned; i++)
          {
            int v1, v2;
            int ndof = first_edge_dof[i+1]-first_edge_dof[i];
            ma.GetEdgePNums (i, v1, v2);
            cnt[v1] += ndof;
            cnt[v2] += ndof;
          }
        for (i = 0; i < nfa; i++)
          cnt[nv+i] = first_face_dof[i+1]-first_face_dof[i];
        for (i = 0; i < ni; i++)
          {
            ma.GetElFaces (i, fanums, orient);
            int ndof = first_element_dof[i+1] - first_element_dof[i];
            for (j = 0; j < fanums.Size(); j++)
              cnt[nv+fanums[j]] += ndof;
          }
        break; 
      case 6: // VEF + I 
        for (i = 0; i < ned; i++)
          {
            int v1, v2;
            ma.GetEdgePNums (i, v1, v2);
            cnt[v1] += first_edge_dof[i+1]-first_edge_dof[i];
            cnt[v2] += first_edge_dof[i+1]-first_edge_dof[i];
          }
        for (i = 0; i < nfa; i++)
          { 
            Array<int>  pnums; 
            ma.GetFacePNums(i,pnums); 
            for(j=0;j<pnums.Size();j++) 
              cnt[pnums[j]] +=  first_face_dof[i+1] - first_face_dof[i];
          }
        for (i = 0; i < ni; i++)
          cnt[nv + i] +=  first_element_dof[i+1] - first_element_dof[i];
        break;
      case 7: // VEFI 
        for(i=0; i<ned; i++)
          {
            int v1, v2;
            int ndof = first_edge_dof[i+1]-first_edge_dof[i];
            ma.GetEdgePNums (i, v1, v2);
            cnt[v1] += ndof;
            cnt[v2] += ndof;
          }
        for (i = 0; i < nfa; i++)
          { 
            Array<int>  pnums; 
            ma.GetFacePNums(i,pnums); 
            int ndof =  first_face_dof[i+1] - first_face_dof[i];
            for(j=0;j<pnums.Size();j++) 
              cnt[pnums[j]] += ndof ;
          }
        for (i = 0; i < ni; i++)
          {
            Array<int>  pnums; 
            ma.GetElPNums(i,pnums); 
            int ndof = first_element_dof[i+1] - first_element_dof[i];
            for (j = 0; j < pnums.Size(); j++)
              cnt[pnums[j]] += ndof;
          }
        break;
      case 12: // VEFI-Cluster 
        cnt =0; 
        for(i=0;i<nv;i++)
          {
            cnt[ma.GetClusterRepVertex(i)]++;
          }
        for(i=0; i<ned; i++)
          {
            int v1, v2;
            int ndof = first_edge_dof[i+1]-first_edge_dof[i];
            ma.GetEdgePNums (i, v1, v2);
            cnt[ma.GetClusterRepVertex(v1)] += ndof;
            if(ma.GetClusterRepVertex(v1)!=ma.GetClusterRepVertex(v2))
              cnt[ma.GetClusterRepVertex(v2)] += ndof;
          }
        for (i = 0; i < nfa; i++)
          { 
            Array<int>  pnums; 
            ma.GetFacePNums(i,pnums); 
            int ndof =  first_face_dof[i+1] - first_face_dof[i];
            Array<int> repv; 
            for(j=0;j<pnums.Size();j++) 
              repv.Append(ma.GetClusterRepVertex(pnums[j]));
            for(j=0;j<pnums.Size();j++) 
              { 
                bool ok=1; 
                for(k=0;k<j;k++) if(repv[j] == repv[k]) ok = 0;  
                if(ok) cnt[repv[j]] += ndof ;
              }
          }
        for (i = 0; i < ni; i++)
          {
            Array<int>  pnums; 
            ma.GetElPNums(i,pnums); 
            Array<int> repv; 
	   
            for(j=0;j<pnums.Size();j++) 
              repv.Append(ma.GetClusterRepVertex(pnums[j]));
            int ndof = first_element_dof[i+1] - first_element_dof[i];
            for (j = 0; j < pnums.Size(); j++)
              { 
                bool ok=1; 
                for(k=0;k<j;k++) if(repv[j] == repv[k]) ok = 0;  
                if(ok) cnt[repv[j]] += ndof ;
              }
	     
          }
        break;

      case 8: // V + E + FI 
        for (i = 0; i < ned; i++)
          cnt[nv+i] = first_edge_dof[i+1]-first_edge_dof[i];
        for (i = 0; i < nfa; i++)
          cnt[nv+ned+i] = first_face_dof[i+1]-first_face_dof[i];
        for (i = 0; i < ni; i++)
          {
            ma.GetElFaces (i,fanums,orient);
            for (k = 0; k < fanums.Size(); k++)
              cnt[nv+ned+fanums[k]] += first_element_dof[i+1] - first_element_dof[i];
          }
        break; 

      case 9: // V + EF + I 
        for (i = 0; i < ned; i++)
          cnt[nv+i]= first_edge_dof[i+1]-first_edge_dof[i];
        for (i = 0; i < nfa; i++)
          {
            ma.GetFaceEdges (i, f2ed);
            int fdof = first_face_dof[i+1]-first_face_dof[i];
            for (j = 0; j < f2ed.Size(); j++)
              cnt[nv+f2ed[j]] +=  fdof; 
          }
        for (i = 0; i < ni; i++)
          cnt[nv+ned+i] = first_element_dof[i+1]-first_element_dof[i];
        break;     

      case 10: // V + EI + FI 
        cnt = 0;  
        for(i=0; i< nv;i++)
          cnt[ma.GetClusterRepVertex(i)]++; 
        for (i = 0; i < ned; i++)
          if(fine_edge[i])
            cnt[ma.GetClusterRepEdge(i)] += first_edge_dof[i+1]-first_edge_dof[i];
       
        for (i = 0; i < nfa; i++)
          if(fine_face[i])
            cnt[ma.GetClusterRepFace(i)] += first_face_dof[i+1]-first_face_dof[i];
        for (i = 0; i < ni; i++)
          {
            int ccl = ma.GetClusterRepElement(i);
            cnt[ccl] +=first_element_dof[i+1] - first_element_dof[i];
          }
        /*
          for (i = 0; i < ni; i++)
          {
          int ccl = ma.GetClusterRepElement(i);
          ma.GetElEdges (i,ednums,orient);
	     
          for (k = 0; k < ednums.Size(); k++)
          if(ccl!= ma.GetClusterRepFace(ednums[k]))
          cnt[ma.GetClusterRepEdge(ednums[k])] += first_element_dof[i+1] - first_element_dof[i];
	      
          ma.GetElFaces (i,fanums,orient);
          for (k = 0; k < fanums.Size(); k++)
          if(ccl != ma.GetClusterRepFace(fanums[k]))
          cnt[ma.GetClusterRepFace(fanums[k])] += first_element_dof[i+1] - first_element_dof[i];
          }*/
        //cout << "cnt " << cnt << endl; 
        //*testout  << "cnt " << cnt << endl; 
        break; 

      case 11: // 2d VEI 
        for(i=0; i<ned; i++)
          if(fine_edge[i])
            {
              int v1, v2;
              int ndof = first_edge_dof[i+1]-first_edge_dof[i];
              ma.GetEdgePNums (i, v1, v2);
              cnt[v1] += ndof;
              cnt[v2] += ndof;
            }
        for (i = 0; i < ni; i++)
          {
            Array<int>  pnums; 
            ma.GetElPNums(i,pnums); 
            int ndof = first_element_dof[i+1] - first_element_dof[i];
            for (j = 0; j < pnums.Size(); j++)
              cnt[pnums[j]] += ndof;
          }
        break;


      case 20:
        {
          cnt = 0;
          for (int i = 0; i < nv; i++)
            cnt[i] = 1;
          for (int i = 0; i < ned; i++)
            if (fine_edge[i])
              {
                int v1, v2;
                int ndof = first_edge_dof[i+1]-first_edge_dof[i];
                ma.GetEdgePNums (i, v1, v2);
                cnt[v1] += ndof;
                cnt[v2] += ndof;
                cnt[nv+i] += ndof;
              }
          for (int i = 0; i < nfa; i++)
            {
              ma.GetFaceEdges (i, f2ed);
              int fdof = first_face_dof[i+1]-first_face_dof[i];
              for (j = 0; j < f2ed.Size(); j++)
                cnt[nv+f2ed[j]] +=  fdof; 
            }
          break;
        }
      }
     
    //    *testout << " cnt " << cnt << endl; 
    
    Table<int> & table = *new Table<int> (cnt); 
    cnt = 0; 
    if(SmoothingType != 10) 
      {
        for (i = 0; i < nv; i++)
          table[i][0] = i;
        if (augmented == 1)
          for (i = 0; i < nv; i++)
            table[i][1] = nv+i;
        else if (augmented == 2)
          for (i = 0; i < nv; i++)
            for (j = 0; j < order-1; j++)
              table[i][j+1] = nv+i*(order-1)+j; 
        for (i = 0; i < nv; i++)
          cnt[i] = nvdof;
      }
 

    ii=0; 
    switch(SmoothingType)
      {
      case 1:  // 2d: V + E + I
        for (i = 0; i < ned; i++,ii++)
          if(fine_edge[i])
            {
              first = first_edge_dof[i];
              int ndof = first_edge_dof[i+1]-first_edge_dof[i]; 
              for (j = 0; j < ndof; j++)
                table[nv+i][j] = first+j;
            }
        for (i = 0; i < nel; i++)
          {
            first = first_element_dof[i];
            int ndof =  first_element_dof[i+1]-first_element_dof[i]; 
            for (j = 0; j < ndof; j++)
              table[nv+ned+i][j] = first+j;
          }
        break; 
      case 2: // 2d VE + I
        for (i = 0; i < ned; i++)
          if(fine_edge[i])
            {
              int v1, v2;
              first = first_edge_dof[i];
              int ndof = first_edge_dof[i+1]-first;
              ma.GetEdgePNums (i, v1, v2);
              for (j = 0; j < ndof; j++)
                {
                  table[v1][cnt[v1]++] = first+j;
                  table[v2][cnt[v2]++] = first+j;
                }
            }
        for (i = 0; i < ni; i++)
          {
            first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first;
            for (j = 0; j < ndof; j++)
              table[nv+i][j] = first+j;
          }
      case 3 : // V + E + F + I 
        for (i = 0; i < ned; i++)
          {
            int ndof = first_edge_dof[i+1]-first_edge_dof[i]; 
            first = first_edge_dof[i];
            for (j = 0; j < ndof; j++)
              table[nv+i][j] = first+j;
          }
        for (i = 0; i < nfa; i++)
          {
            int ndof = first_face_dof[i+1]-first_face_dof[i]; 
            first = first_face_dof[i];
            for (j = 0; j < ndof; j++)
              table[nv+ned+i][j] = first+j;
          }
        for (i = 0; i < ni; i++)
          {
            first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i]; 
            for (j = 0; j < ndof; j++)
              table[nv+ned+nfa+i][j] = first+j;
          }
        break; 
      case 4: // VE + F + I
        for (i = 0; i < ned; i++)
          {
            int v1, v2;
            first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first_edge_dof[i];
            ma.GetEdgePNums (i, v1, v2);
            for (j = 0; j < ndof; j++)
              {
                table[v1][cnt[v1]++] = first+j;
                table[v2][cnt[v2]++] = first+j;
              }
          }
        for (i = 0; i < nfa; i++)
          {
            first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first_face_dof[i];
            for (j = 0; j < ndof; j++)
              table[nv+i][j] = first+j;
            cnt[nv+i] = ndof;
          }
        for (i = 0; i < ni; i++)
          {
            first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i]; 
            for (j = 0; j < ndof; j++)
              table[nv+nfa+i][j] = first+j;
          }
        break; 
      case 5: // VE + FI
        for (i = 0; i < ned; i++)
          {
            int v1, v2;
            first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first;
            ma.GetEdgePNums (i, v1, v2);
            for (j = 0; j < ndof; j++)
              {
                table[v1][cnt[v1]++] = first+j;
                table[v2][cnt[v2]++] = first+j;
              }
          }
        for (i = 0; i < nfa; i++)
          {
            first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first_face_dof[i];
            for (j = 0; j < ndof; j++)
              table[nv+i][j] = first+j;
            cnt[nv+i] = ndof;
          }
	
        for (i = 0; i < ni; i++)
          {
            ma.GetElFaces (i, fanums, orient);
            first = first_element_dof[i];
            int ndof = first_element_dof[i+1] - first_element_dof[i];
            for (j = 0; j < fanums.Size(); j++)
              for (k = 0; k < ndof; k++)
                table[nv+fanums[j]][cnt[nv+fanums[j]]++] = first + k;
          }
        break; 
      case 6: // VEF + I
        for (i = 0; i < ned; i++)
          {
            int v1, v2;
            first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first;
            ma.GetEdgePNums (i, v1, v2);
            for (j = 0; j < ndof; j++)
              {
                table[v1][cnt[v1]++] = first+j;
                table[v2][cnt[v2]++] = first+j;
              }
          }
        for (i = 0; i < nfa; i++)
          {
            Array<int> pnums; 
            ma.GetFacePNums(i,pnums); 
            first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first_face_dof[i];
            for(k=0;k<pnums.Size();k++)
              for (j = 0; j < ndof; j++)
                table[pnums[k]][cnt[pnums[k]]++] = first+j;
          }
        for (i = 0; i < ni; i++)
          {
            first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first;
            for (j = 0; j < ndof; j++)
              table[nv+i][j] = first+j;
          }
        break; 
      case 7: // VEFI
        for (i = 0; i < ned; i++)
          {
            int v1, v2;
            first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first;
            ma.GetEdgePNums (i, v1, v2);
            for (j = 0; j < ndof; j++)
              {
                table[v1][cnt[v1]++] = first+j;
                table[v2][cnt[v2]++] = first+j;
              }
          }
        for (i = 0; i < nfa; i++)
          {
            Array<int> pnums; 
            ma.GetFacePNums(i,pnums); 
            first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first_face_dof[i];
            for(k=0;k<pnums.Size();k++)
              for (j = 0; j < ndof; j++)
                table[pnums[k]][cnt[pnums[k]]++] = first+j;
          }
        for (i = 0; i < ni; i++)
          {
            Array<int> pnums; 
            ma.GetElPNums(i,pnums); 
            first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i];
            for(k=0;k<pnums.Size();k++)
              for (j = 0; j < ndof; j++)
                table[pnums[k]][cnt[pnums[k]]++] = first+j;
          }
        break;
      case 12: // VEFI
        //cnt =0; 
        for(i=0;i<nv;i++)
          table[ma.GetClusterRepVertex(i)][cnt[ma.GetClusterRepVertex(i)]++]=i;
	     
        for (i = 0; i < ned; i++)
          {
            int v1, v2;
            first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first;
            ma.GetEdgePNums (i, v1, v2);
            int r1 = ma.GetClusterRepVertex(v1); 
            int r2 = ma.GetClusterRepVertex(v2); 
	    
            for (j = 0; j < ndof; j++)
              {
                table[r1][cnt[r1]++] = first+j;
                if(r1!=r2)
                  table[v2][cnt[v2]++] = first+j;
              }
          }
        for (i = 0; i < nfa; i++)
          {
            Array<int> pnums; 
            ma.GetFacePNums(i,pnums); 
	    
            first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first_face_dof[i];
            Array<int> repv; 
	  
            for(j=0;j<pnums.Size();j++) 
              repv.Append(ma.GetClusterRepVertex(pnums[j]));
            for(j=0;j<pnums.Size();j++) 
              { 
                bool ok=1; 
                for(k=0;k<j;k++) if(repv[j] == repv[k]) ok = 0;  
		
                if(ok)
                  for (k = 0; k < ndof; k++)
                    table[repv[j]][cnt[repv[j]]++] += first+k ;
              }
          }
        for (i = 0; i < ni; i++)
          {
            Array<int> pnums; 
            ma.GetElPNums(i,pnums); 
            first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i];
	    
	    
            Array<int> repv; 
            for(j=0;j<pnums.Size();j++) 
              repv.Append( ma.GetClusterRepVertex(pnums[j]));
            for(j=0;j<pnums.Size();j++) 
              { 
                bool ok=1; 
                for(k=0;k<j;k++) if(repv[j] == repv[k]) ok = 0;  
		
                if(ok)
                  for (k = 0; k < ndof; k++)
                    table[repv[j]][cnt[repv[j]]++] += first+k ;
              }
          }
	   
        break;
 
        for (i = 0; i < ned; i++)
          {
            int v1, v2;
            first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first;
            ma.GetEdgePNums (i, v1, v2);
            for (j = 0; j < ndof; j++)
              {
                table[v1][cnt[v1]++] = first+j;
                table[v2][cnt[v2]++] = first+j;
              }
          }
        for (i = 0; i < nfa; i++)
          {
            Array<int> pnums; 
            ma.GetFacePNums(i,pnums); 
            first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first_face_dof[i];
            for(k=0;k<pnums.Size();k++)
              for (j = 0; j < ndof; j++)
                table[pnums[k]][cnt[pnums[k]]++] = first+j;
          }
        for (i = 0; i < ni; i++)
          {
            Array<int> pnums; 
            ma.GetElPNums(i,pnums); 
            first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i];
            for(k=0;k<pnums.Size();k++)
              for (j = 0; j < ndof; j++)
                table[pnums[k]][cnt[pnums[k]]++] = first+j;
          }
        break;
 
      case 8: // V + E + FI 
        for (i = 0; i < ned; i++)
          {
            first = first_edge_dof[i];
            for (j = first; j < first_edge_dof[i+1]; j++)
              table[nv+i][cnt[nv+i]++] = j;
          }
        for (i = 0; i < nfa; i++)
          {
            first = first_face_dof[i];
            for (j = first; j < first_face_dof[i+1]; j++)
              table[nv+ned+i][cnt[nv+ned+i]++] = j;
          }
        for (i = 0; i < nel; i++)
          {
            ma.GetElFaces (i,fanums,orient);
            for (k = 0; k < fanums.Size(); k++)
              for (j = first_element_dof[i]; j < first_element_dof[i+1]; j++)
                table[nv+ned+fanums[k]][cnt[nv+ned+fanums[k]]++] = j;
          }
        break;
      case 9: // V + EF + I
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
            first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first;
            for (j = 0; j < ndof; j++)
              table[nv+ned+i][j] = first+j;
          }
        break; 
      case 10: // V + EI + FI
        cnt =0; 

        for (i=0;i<nv ; i++) 
          table[ma.GetClusterRepVertex(i)][cnt[ma.GetClusterRepVertex(i)]++] = i; 	

	
        for (i = 0; i < ned; i++)
          if(fine_edge[i])
            for (j = first_edge_dof[i]; j < first_edge_dof[i+1]; j++)
              table[ma.GetClusterRepEdge(i)][cnt[ma.GetClusterRepEdge(i)]++] = j;
	

        for (i = 0; i < nfa; i++)
          if(fine_face[i])
            for (j = first_face_dof[i]; j < first_face_dof[i+1]; j++)
              table[ma.GetClusterRepFace(i)][cnt[ma.GetClusterRepFace(i)]++] = j;

        for(i=0;i<ni;i++)
          {
            int ccl = ma.GetClusterRepElement(i);
            for (j = first_element_dof[i]; j < first_element_dof[i+1]; j++)
              table[ma.GetClusterRepElement(i)][cnt[ccl]++] = j;
          }
        /*	for (i = 0; i < ni; i++)
          {
          int ccl = ma.GetClusterRepElement(i);
          ma.GetElEdges (i,ednums,orient);
          for (k = 0; k < ednums.Size(); k++)
          //  if(ccl!=ma.GetClusterRepEdge(ednums[k]))
          for (j = first_element_dof[i]; j < first_element_dof[i+1]; j++)
          table[ma.GetClusterRepEdge(ednums[k])][cnt[ma.GetClusterRepEdge(ednums[k])]++] = j;
	    
          ma.GetElFaces (i,fanums,orient);
	   
          for (k = 0; k < fanums.Size(); k++)
          //  if(ccl!=ma.GetClusterRepFace(fanums[k]))
          for (j = first_element_dof[i]; j < first_element_dof[i+1]; j++)
          table[ma.GetClusterRepFace(fanums[k])][cnt[ma.GetClusterRepFace(fanums[k])]++] = j;
          }	*/
        //cout << " test 5 " << endl; 
        //cout << "table " << table << endl; 
        //*testout  << "table " << table << endl; 
        break;
      case 11: //2d VEF
        for (i = 0; i < ned; i++)
          if(fine_edge[i])
            {
              int v1, v2;
              first = first_edge_dof[i];
              int ndof = first_edge_dof[i+1]-first;
              ma.GetEdgePNums (i, v1, v2);
              for (j = 0; j < ndof; j++)
                {
                  table[v1][cnt[v1]++] = first+j;
                  table[v2][cnt[v2]++] = first+j;
                }
            }
        for (i = 0; i < ni; i++)
          {
            Array<int> pnums; 
            ma.GetElPNums(i,pnums); 
            first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i];
            for(k=0;k<pnums.Size();k++)
              for (j = 0; j < ndof; j++)
                table[pnums[k]][cnt[pnums[k]]++] = first+j;
          }
        break;


      case 20:
        {
          cnt = 0;
          for (int i = 0; i < nv; i++)
            cnt[i] = 1;
          for (int i = 0; i < ned; i++)
            if (fine_edge[i])
              {
                int v1, v2;
                int first = first_edge_dof[i];
                int ndof = first_edge_dof[i+1]-first;
                ma.GetEdgePNums (i, v1, v2);
                for (int j = 0; j < ndof; j++)
                  {
                    table[v1][cnt[v1]++] = first+j;
                    table[v2][cnt[v2]++] = first+j;
                    table[nv+i][cnt[nv+i]++] = first+j;
                  }
              }
          for (int i = 0; i < nfa; i++)
            {
              first = first_face_dof[i];
              int ndof = first_face_dof[i+1]-first;
              ma.GetFaceEdges (i, f2ed);
              for (k = 0; k < f2ed.Size(); k++)
                for (j = 0; j < ndof; j++)
                  table[nv+f2ed[k]][cnt[nv+f2ed[k]]++] = first+j;
            }
          break;
        }

      }
    
    /*
      (*testout) << "H1HO-table = " << table << endl;
    */
        return &table;
  }


  Array<int> * 
  H1HighOrderFESpace :: CreateDirectSolverClusters (const Flags & flags) const
  {
    // return 0;

    int i, j, k;

    int nv = ma.GetNV();
    int nd = GetNDof();
    int ne = ma.GetNE();
    Array<int> & clusters = *new Array<int> (GetNDof());
    clusters = 0;

    // all vertices in global space
    /*
      for (i = 0; i < nv; i++)
      clusters[i] = 1;

      // all edges in direct solver cluster

      for (i = first_edge_dof[0]; i < first_edge_dof.Last(); i++)
      clusters[i] = 1;
      return &clusters;

    */
   
    // All Vertical Edges in one Cluster for Hex and Prism (-> 2d Problems !) 
    Array<int> ednums, edorient,fnums, forient;

    //Array<int> & clusters = *new Array<int> (nd);
    //clusters = 0;


    
    for (i = 0; i < ne; i++)
      {
        if (ma.GetElType(i) == ET_PRISM)
          {
            ma.GetElEdges (i, ednums, edorient);
            for (j = 6; j < 9; j++)  //vertical Edges 
              { 
                int first = first_edge_dof[ednums[j]];
                int next = first_edge_dof[ednums[j]+1];
                for (k = first; k < next; k++)
                  clusters[k] = 2;
              }
	    
            ma.GetElFaces(i,fnums,forient); // vertical faces 
            for(j=2;j<5;j++) 
              {
	
                int first = first_face_dof[fnums[j]]; 
                int next = first_face_dof[fnums[j]+1]; 
		
                for (k=first; k < next; k++) 
                  clusters[k]=0; 
		
                //INT<2> p = order_face[fnums[j]];
                //for(k=first + 2*(p[0]+1)*(p[1]+1);k<next;k++)
                //  clusters[k]=3;  
              }
          }

        else if (ma.GetElType(i) == ET_HEX)  
          {
            ma.GetElEdges (i, ednums, edorient);
            for (j = 8; j < 12; j++) //vertical edges
              {
                int first = first_edge_dof[ednums[j]];
                int next = first_edge_dof[ednums[j]+1];
                for (k = first; k < next; k++)
                  clusters[k] = 2;
              }
            ma.GetElFaces(i,fnums,forient); // vertical faces 
            for(j=2;j<6;j++) 
              {
                
                int first = first_face_dof[fnums[j]]; 
                int next = first_face_dof[fnums[j]+1]; 
		
                for (k=first; k < next; k++) 
                  clusters[k]=3; 
              }
          } 
      }


   
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

   

    for(i=0; i< adddirectsolverdofs.Size(); i++)
      {
        clusters[adddirectsolverdofs[i]] = 5;
      }

    const int stdoffset = 6;

    /*
      (*testout) << "directvertexclusters" << directvertexclusters << endl;
      (*testout) << "directedgeclusters" << directedgeclusters << endl;
      (*testout) << "directfaceclusters" << directfaceclusters << endl;
      (*testout) << "directelementclusters" << directelementclusters << endl;
    */

    for(i=0; i<directvertexclusters.Size(); i++)
      if(directvertexclusters[i] >= 0)
        clusters[i] = directvertexclusters[i] + stdoffset;

    for(i=0; i<directedgeclusters.Size(); i++)
      if(directedgeclusters[i] >= 0)
        for(j = first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
          clusters[j] = directedgeclusters[i] + stdoffset;

    for(i=0; i<directfaceclusters.Size(); i++)
      if(directfaceclusters[i] >= 0)
        for(j = first_face_dof[i]; j<first_face_dof[i+1]; j++)
          clusters[j] = directfaceclusters[i] + stdoffset;
	  
    for(i=0; i<directelementclusters.Size(); i++)
      if(directelementclusters[i] >= 0)
        for(j = first_element_dof[i]; j<first_element_dof[i+1]; j++)
          clusters[j] = directelementclusters[i] + stdoffset;


    //    (*testout) << "clusters " << clusters << endl;
    
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


#ifdef OLD_PARALLEL_UPDATE

  das sollte nicht definiert sein

  void H1HighOrderFESpace :: UpdateParallelDofs_hoproc()
  {
    // ******************************
    // update exchange dofs 
    // ******************************
    *testout << "H1Ho::UpdateParallelDofs_hoproc" << endl;
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

    //     // brauchen wir eigentlich nicht mehr
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

            while ( ii < distantdofs[dest]->Size() )
              {
                int nodenum = (*distantdofs[dest])[ii++];
                int isdistghost = (*distantdofs[dest])[ii++];
                Array<int> dnums;
                GetNodeDofNrs (nt, nodenum, dnums);
                for ( int i=0; i<dnums.Size(); i++)
                  {
                    // 		    (*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
                    // 		    (*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];


                    (*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];

                    if ( dest < id && !isdistghost )
                      paralleldofs->ismasterdof.Clear ( dnums[i] ) ;
                    ii++; 
                    cnt_nexdof[dest]++;
                  }
              }
          }

        // das ist abgegangen (JS -> AS) !
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
   
    MPI_Wait ( &recvrequest[dest], &status );

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

    for ( int i = 0; i < ntasks; i++ )
      delete distantdofs[i], owndofs[i];

    delete [] owndofs, distantdofs;
    // delete [] sendrequest, recvrequest;


    //     *testout << "localexchangedof = " << endl 
    //              << *paralleldofs->localexchangedof << endl;
    //     *testout << "distantexchangedof = " << endl 
    //              << *paralleldofs->distantexchangedof << endl;
  }

#endif

  void H1HighOrderFESpace :: UpdateParallelDofs_loproc()
  {
    *testout << "H1Ho::UpdateParallelDofs_loproc" << endl;
    // lo-proc should never use high order space
    return;
 
  }
#endif

  // register FESpaces
  namespace h1hofespace_cpp
  {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("h1hotp", H1HighOrderFESpace::Create);
    }
    
    Init init_h1hofespace;
  }
  
}
