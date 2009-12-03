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
#include <parallelngs.hpp>
#include "../fem/h1hofe.hpp"
#include "../fem/h1hofefo.hpp"


using namespace ngmg; 

namespace ngfem
{
  extern int link_it_h1hofefo;
}

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
    // DefineNumFlag("augmented");
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
    DefineDefineFlag("print");
    DefineDefineFlag("noprint");
    if (parseflags) CheckFlags(flags);
    
    // augmented = int (flags.GetNumFlag ("augmented", 0));
    // plate = int (flags.GetDefineFlag ("plate"));
    print = (flags.GetDefineFlag("print")); 
  
    // Variable order space: 
    //      in case of (var_order && order) or (relorder) 
    var_order = flags.GetDefineFlag("variableorder");  
    fixed_order = flags.GetDefineFlag("fixedorder");  
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
    if (flags.NumListFlagDefined ("dirichlet")) 
      loflags.SetFlag ("dirichlet", flags.GetNumListFlag ("dirichlet"));

#ifndef PARALLEL
    low_order_space = new NodalFESpace (ma, loflags);
#else
    low_order_space = new ParallelNodalFESpace (ma, loflags);
#endif
          
    minext = flags.GetDefineFlag ("minext");
    optext = flags.GetDefineFlag ("optext");
    // fast_pfem = flags.GetDefineFlag ("fast");
    

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
    FESpace :: Update (lh);

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
    fine_edge.SetSize(ned); 
    fine_face.SetSize(nfa); 

    fine_edge = 0;
    fine_face = 0; 
    int p = var_order ?  1 : order; 
    
    order_edge = p; 
    order_face = INT<2>(p,p);
    order_inner = INT<3>(p,p,p); 
	
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
          case ET_SEGM:
            cerr << "should not be" << endl;
            break;
	  }
      } 
    first_element_dof[nel] = ndof;
   

    if (print)
      {
        (*testout) << "h1 first edge = " << first_edge_dof << endl;
        (*testout) << "h1 first face = " << first_face_dof << endl;
        (*testout) << "h1 first inner = " << first_element_dof << endl;
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

    if (fixed_order && ma.GetElType(elnr) == ET_TRIG && order <= 6)
      {
        H1HighOrderFiniteElementFO<2> * hofe2d = 0;
        switch (order)
          {
          case 1: hofe2d = new (lh)  H1HighOrderFEFO<ET_TRIG,1> (); break;
          case 2: hofe2d = new (lh)  H1HighOrderFEFO<ET_TRIG,2> (); break;
          case 3: hofe2d = new (lh)  H1HighOrderFEFO<ET_TRIG,3> (); break;
          case 4: hofe2d = new (lh)  H1HighOrderFEFO<ET_TRIG,4> (); break;
          case 5: hofe2d = new (lh)  H1HighOrderFEFO<ET_TRIG,5> (); break;
          case 6: hofe2d = new (lh)  H1HighOrderFEFO<ET_TRIG,6> (); break;
          }
    
        Ng_Element ngel = ma.GetElement<2> (elnr);
        for (int j = 0; j < 3; j++)
          hofe2d->SetVertexNumber (j, ngel.vertices[j]);
        return *hofe2d;
      }

    try
      {
        switch (ma.GetElType(elnr))
          {
          case ET_TET:     hofe3d = new (lh) H1HighOrderFE<ET_TET> (); break;
          case ET_PYRAMID: hofe3d = new (lh) H1HighOrderFE<ET_PYRAMID> (); break;
          case ET_PRISM:   hofe3d = new (lh) H1HighOrderFE<ET_PRISM> (); break;
          case ET_HEX:     hofe3d = new (lh) H1HighOrderFE<ET_HEX> ();  break;
          case ET_TRIG:    hofe2d = new (lh) H1HighOrderFE<ET_TRIG> (); break;
          case ET_QUAD:    hofe2d = new (lh) H1HighOrderFE<ET_QUAD> (); break;

          default:
            {
              throw Exception (string ("GetFE not supported for element") + 
                               ElementTopology::GetElementName(ma.GetElType(elnr)));
            }
          }


        if (ma.GetDimension() == 2)
          {
            Ng_Element ngel = ma.GetElement<2> (elnr);
        
#ifdef PARALLEL
            for (int j = 0; j < ngel.vertices.Size())
              if (ntasks > 1)
                for (int j = 0; j < ngel.vertices.Size(); j++)
                  hofe2d -> SetVertexNumber (j, parallelma->GetDistantPNum(0, ngel.vertices[j]));
              else
#endif
                for (int j = 0; j < ngel.vertices.Size(); j++)
                  hofe2d -> SetVertexNumber (j, ngel.vertices[j]);
            
            for (int j = 0; j < ngel.edges.Size(); j++)
              hofe2d -> SetOrderEdge (j, order_edge[ngel.edges[j]]);
            
            INT<2> p(order_inner[elnr][0], order_inner[elnr][1]);
            hofe2d -> SetOrderFace (0, p);

            hofe2d -> ComputeNDof();

            return *hofe2d;
          }

        else

          {
            Ng_Element ngel = ma.GetElement<3> (elnr);

#ifdef PARALLEL
            for (int j = 0; j < ngel.vertices.Size())
              if (ntasks > 1)
                for (int j = 0; j < ngel.vertices.Size(); j++)
                  hofe3d -> SetVertexNumber (j, parallelma->GetDistantPNum(0, ngel.vertices[i]));
              else
#endif
                for (int j = 0; j < ngel.vertices.Size(); j++)
                  hofe3d -> SetVertexNumber (j, ngel.vertices[j]);

            for (int j = 0; j < ngel.edges.Size(); j++)
              hofe3d -> SetOrderEdge (j, order_edge[ngel.edges[j]]);
            
            for (int j = 0; j < ngel.faces.Size(); j++)
              hofe3d -> SetOrderFace (j, order_face[ngel.faces[j]]);

            hofe3d -> SetOrderCell (order_inner[elnr]);

            hofe3d -> ComputeNDof();
            return *hofe3d;
          }
      }

    catch (Exception & e)
      {
        e.Append ("in H1HoFESpace::GetElement\n");
        throw;
      }
  }
 
 


  const FiniteElement & H1HighOrderFESpace :: GetSFE (int elnr, LocalHeap & lh) const
  {
    H1HighOrderFiniteElement<1> * hofe1d = NULL;
    H1HighOrderFiniteElement<2> * hofe2d = NULL;

    switch (ma.GetSElType(elnr))
      {
      case ET_TRIG: hofe2d = new (lh) H1HighOrderFE<ET_TRIG> (); break;
      case ET_QUAD: hofe2d = new (lh) H1HighOrderFE<ET_QUAD> (); break;
      case ET_SEGM: hofe1d = new (lh) H1HighOrderFE<ET_SEGM> (); break;
      default:
        {
          throw Exception ("GetFE not supported for element");
        }
      }
  
    Ng_Element ngel = ma.GetSElement(elnr);
    
#ifdef PARALLEL
    if ( ntasks > 1 )
      for ( int i = 0; i < vnums.Size(); i++ )
        vnums[i] = parallelma->GetDistantPNum(0, vnums[i]);
    cout << "have to set vertex numbers!!!" << endl;
#endif 
    
    if (ma.GetDimension() == 2)
      {
        for (int j = 0; j < ngel.vertices.Size(); j++)
          hofe1d -> SetVertexNumber (j, ngel.vertices[j]);

        for (int j = 0; j < ngel.edges.Size(); j++)
          hofe1d -> SetOrderEdge (j, order_edge[ngel.edges[j]]);

        hofe1d -> ComputeNDof();
        return *hofe1d;
      }
    else
      {
        for (int j = 0; j < ngel.vertices.Size(); j++)
          hofe2d -> SetVertexNumber (j, ngel.vertices[j]);

        for (int j = 0; j < ngel.edges.Size(); j++)
          hofe2d -> SetOrderEdge (j, order_edge[ngel.edges[j]]);

        INT<2> p = order_face[ma.GetSElFace(elnr)];
        FlatArray<INT<2> > of(1, &p);
        hofe2d -> SetOrderFace (of);

        hofe2d  -> ComputeNDof();
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
    Ng_Element ngel = ma.GetElement(elnr);

    dnums.SetSize(ngel.vertices.Size()); 
    for (int i = 0; i < ngel.vertices.Size(); i++)
      dnums[i] = ngel.vertices[i];

    for (int i = 0; i < ngel.edges.Size(); i++)
      dnums += GetEdgeDofs (ngel.edges[i]);

    if (ma.GetDimension() == 3)
      for (int i = 0; i < ngel.faces.Size(); i++)
        dnums += GetFaceDofs (ngel.faces[i]);

    dnums += GetElementDofs (elnr);
    if (!DefinedOn (ma.GetElIndex (elnr)))
      dnums = -1;
  }


  void  H1HighOrderFESpace :: GetExternalDofNrs (int elnr, Array<int> & dnums) const
  {
    if (!eliminate_internal) 
      {
        GetDofNrs (elnr, dnums);
        return;
      }

    Ng_Element ngel = ma.GetElement(elnr);

    dnums.SetSize(ngel.vertices.Size()); 
    for (int i = 0; i < ngel.vertices.Size(); i++)
      dnums[i] = ngel.vertices[i];

    for (int i = 0; i < ngel.edges.Size(); i++)
      dnums += GetEdgeDofs (ngel.edges[i]);

    if (ma.GetDimension() == 3)
      for (int i = 0; i < ngel.faces.Size(); i++)
        dnums += GetFaceDofs (ngel.faces[i]);

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

    if (ma.GetDimension() == 2)
      {
        for (int i = 0; i < ednums.Size(); i++)
          if (order_edge[ednums[i]] > 1)
            dnums.Append (first_edge_dof[ednums[i]]);
      }

    if (ma.GetDimension() == 3)
      {
        for (int i = 0; i < ednums.Size(); i++)
          for (int j = first_edge_dof[ednums[i]]; j < first_edge_dof[ednums[i]+1]; j++)
            dnums.Append (j);

	/*
        ArrayMem<int,12> fanums;
        ma.GetElFaces (elnr, fanums);
        for (int i = 0; i < fanums.Size(); i++)
          if (order_face[fanums[i]][0] > 1)
            dnums.Append (first_face_dof[fanums[i]]);
	*/
      }
  }





  
  void H1HighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    dnums.Append (vnr);
  }
  
  
  void H1HighOrderFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    dnums += GetEdgeDofs (ednr);
    /*
    int first = first_edge_dof[ednr];
    int next = first_edge_dof[ednr+1];
    for (int j = first; j < next; j++)
      dnums.Append (j);
    */
  }

  void H1HighOrderFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    if (ma.GetDimension() < 3) return;

    dnums += GetFaceDofs (fanr);

    /*
    int first = first_face_dof[fanr];
    int next = first_face_dof[fanr+1];
    for (int j = first; j < next; j++)
      dnums.Append (j);
    */
  }


  void H1HighOrderFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    dnums += GetElementDofs (elnr);

    /*
    int first = first_element_dof[elnr];
    int next = first_element_dof[elnr+1];
    for (int j = first; j < next; j++)
      dnums.Append (j);
    */
  }
  
  void H1HighOrderFESpace :: 
  GetSDofNrs (int elnr, Array<int> & dnums) const
  {
    Ng_Element ngel = ma.GetSElement(elnr);

    dnums.SetSize(ngel.vertices.Size()); 

    for (int i = 0; i < ngel.vertices.Size(); i++)
      dnums[i] = ngel.vertices[i];

    for (int i = 0; i < ngel.edges.Size(); i++)
      dnums += GetEdgeDofs (ngel.edges[i]);

    if (ma.GetDimension() == 3)
      dnums += GetFaceDofs (ngel.faces[0]);
    

    /*
    // what's that (JS) ?
    static bool getall = false;

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
        for (int i = 0; i<defined_on_one_side_of_bounding_curve.Size(); i++)
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
    */

    // what's that (JS) ?
    // I still don't understand it .... but it's now threadsafe
    if (!DefinedOnBoundary (ma.GetSElIndex (elnr)))
      dnums = -1;


    if(defined_on_one_side_of_bounding_curve.Size() > 0)
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
        for (int i = 0; i<defined_on_one_side_of_bounding_curve.Size(); i++)
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
                        Array<int> neighbour_dnums;
                        // GetSDofNrs(neighbours[j],neighbour_dnums);

			Ng_Element nbel = ma.GetSElement(neighbours[j]);
			
			neighbour_dnums.SetSize(nbel.vertices.Size()); 
			
			for (int i = 0; i < nbel.vertices.Size(); i++)
			  neighbour_dnums[i] = nbel.vertices[i];
			for (int i = 0; i < nbel.edges.Size(); i++)
			  neighbour_dnums += GetEdgeDofs (nbel.edges[i]);
			if (ma.GetDimension() == 3)
			  neighbour_dnums += GetFaceDofs (nbel.faces[0]);

			
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
    int SmoothingType = int(precflags.GetNumFlag("blocktype",0)); 

    Array<int> ednums, fanums, vnums,f2ed; 
    
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
      SmoothingType = 4; 

    cout << " blocktype " << SmoothingType << endl; 
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
        cout << "VE + EF + I" << endl;
        ncnt = nv+ned;
        break;

      case 21:
        cout << "Helmholtz" << endl;
        ncnt = ned+nfa;
        break;

      default: 
        cout << " Error in H1HOFESpace :: CreateSmoothingBlocks no SmoothingType " << endl; 
        return 0; 
      }
    
    Array<int> cnt(ncnt); 
    cnt = 0; 
    int nvdof = 1; 
    if (SmoothingType <= 20)
      for (int i = 0; i < nv; i++)
	cnt[i] = nvdof; 
      
    int ii =0; 
    switch(SmoothingType)
      { 
      case 1: // 2d V + E + I 
        for (int i = 0; i < ned; i++, ii++)
          if(fine_edge[i])
            cnt[nv+ii] = GetEdgeDofs(i).Size(); 
        for (int i = 0; i < ni; i++)
          cnt[nv+ned+i] = GetElementDofs(i).Size(); 
        break;
      case 2: // 2d VE + I 
        for (int i = 0; i < ned; i++ )
          if(fine_edge[i])
            {
              Ng_Node<1> edge = ma.GetNode<1> (i);
              cnt[edge.vertices[0]] += GetEdgeDofs(i).Size();
              cnt[edge.vertices[1]] += GetEdgeDofs(i).Size();
            }
        for (int i = 0; i < ni; i++)
          cnt[nv+i] = first_element_dof[i+1]-first_element_dof[i];
        break; 
      case 3: // V + E + F + I 
        for (int i = 0; i < ned; i++)
          cnt[nv+i] = first_edge_dof[i+1]-first_edge_dof[i];
        for (int i = 0; i < nfa; i++)
          cnt[nv+ned+i] = first_face_dof[i+1]-first_face_dof[i];
        for (int i = 0; i < ni; i++)
          cnt[nv+ned+nfa+i] = first_element_dof[i+1]-first_element_dof[i];
        break; 
      case 4: // VE + F + I 
	for (int i = 0; i < nv; i++)
          cnt[i] = IsDirichletVertex(i) ? 0 : 1;
        for (int i = 0; i < ned; i++)
          if (!IsDirichletEdge(i))
            {
              Ng_Node<1> enode = ma.GetNode<1> (i);
              cnt[enode.vertices[0]] += GetEdgeDofs(i).Size();
              cnt[enode.vertices[1]] += GetEdgeDofs(i).Size();
            }
        for (int i = 0; i < nfa; i++)
	  if (!IsDirichletFace(i))
	    cnt[nv+i] = first_face_dof[i+1]-first_face_dof[i];
        for (int i = 0; i < ni; i++)
          cnt[nv+nfa+i] = first_element_dof[i+1]-first_element_dof[i];
        break; 
      case 5: // VE + FI 
        for (int i = 0; i < ned; i++)
          {
            int v1, v2;
            int ndof = first_edge_dof[i+1]-first_edge_dof[i];
            ma.GetEdgePNums (i, v1, v2);
            cnt[v1] += ndof;
            cnt[v2] += ndof;
          }
        for (int i = 0; i < nfa; i++)
          cnt[nv+i] = first_face_dof[i+1]-first_face_dof[i];
        for (int i = 0; i < ni; i++)
          {
            ma.GetElFaces (i, fanums);
            int ndof = first_element_dof[i+1] - first_element_dof[i];
            for (int j = 0; j < fanums.Size(); j++)
              cnt[nv+fanums[j]] += ndof;
          }
        break; 
      case 6: // VEF + I 
        for (int i = 0; i < ned; i++)
          {
            int v1, v2;
            ma.GetEdgePNums (i, v1, v2);
            cnt[v1] += first_edge_dof[i+1]-first_edge_dof[i];
            cnt[v2] += first_edge_dof[i+1]-first_edge_dof[i];
          }
        for (int i = 0; i < nfa; i++)
          { 
            Array<int>  pnums; 
            ma.GetFacePNums(i,pnums); 
            for(int j=0;j<pnums.Size();j++) 
              cnt[pnums[j]] +=  first_face_dof[i+1] - first_face_dof[i];
          }
        for (int i = 0; i < ni; i++)
          cnt[nv + i] +=  first_element_dof[i+1] - first_element_dof[i];
        break;
      case 7: // VEFI 
        for (int i=0; i<ned; i++)
          {
            int v1, v2;
            int ndof = first_edge_dof[i+1]-first_edge_dof[i];
            ma.GetEdgePNums (i, v1, v2);
            cnt[v1] += ndof;
            cnt[v2] += ndof;
          }
        for (int i = 0; i < nfa; i++)
          { 
            Array<int>  pnums; 
            ma.GetFacePNums(i,pnums); 
            int ndof =  first_face_dof[i+1] - first_face_dof[i];
            for (int j =0;j<pnums.Size();j++) 
              cnt[pnums[j]] += ndof ;
          }
        for (int i = 0; i < ni; i++)
          {
            Array<int>  pnums; 
            ma.GetElPNums(i,pnums); 
            int ndof = first_element_dof[i+1] - first_element_dof[i];
            for (int j = 0; j < pnums.Size(); j++)
              cnt[pnums[j]] += ndof;
          }
        break;
      case 12: // VEFI-Cluster 
        cnt =0; 
        for(int i=0;i<nv;i++)
          {
            cnt[ma.GetClusterRepVertex(i)]++;
          }
        for (int i =0; i<ned; i++)
          {
            int v1, v2;
            int ndof = first_edge_dof[i+1]-first_edge_dof[i];
            ma.GetEdgePNums (i, v1, v2);
            cnt[ma.GetClusterRepVertex(v1)] += ndof;
            if(ma.GetClusterRepVertex(v1)!=ma.GetClusterRepVertex(v2))
              cnt[ma.GetClusterRepVertex(v2)] += ndof;
          }
        for (int i = 0; i < nfa; i++)
          { 
            Array<int>  pnums; 
            ma.GetFacePNums(i,pnums); 
            int ndof =  first_face_dof[i+1] - first_face_dof[i];
            Array<int> repv; 
            for (int j =0;j<pnums.Size();j++) 
              repv.Append(ma.GetClusterRepVertex(pnums[j]));
            for (int j =0;j<pnums.Size();j++) 
              { 
                bool ok=1; 
                for (int k = 0;k<j;k++) if(repv[j] == repv[k]) ok = 0;  
                if(ok) cnt[repv[j]] += ndof ;
              }
          }
        for (int i = 0; i < ni; i++)
          {
            Array<int>  pnums; 
            ma.GetElPNums(i,pnums); 
            Array<int> repv; 
	   
            for (int j =0;j<pnums.Size();j++) 
              repv.Append(ma.GetClusterRepVertex(pnums[j]));
            int ndof = first_element_dof[i+1] - first_element_dof[i];
            for (int j = 0; j < pnums.Size(); j++)
              { 
                bool ok=1; 
                for (int k = 0;k<j;k++) if(repv[j] == repv[k]) ok = 0;  
                if(ok) cnt[repv[j]] += ndof ;
              }
	     
          }
        break;

      case 8: // V + E + FI 
        for (int i = 0; i < ned; i++)
          cnt[nv+i] = first_edge_dof[i+1]-first_edge_dof[i];
        for (int i = 0; i < nfa; i++)
          cnt[nv+ned+i] = first_face_dof[i+1]-first_face_dof[i];
        for (int i = 0; i < ni; i++)
          {
            ma.GetElFaces (i,fanums);
            for (int k = 0; k < fanums.Size(); k++)
              cnt[nv+ned+fanums[k]] += first_element_dof[i+1] - first_element_dof[i];
          }
        break; 

      case 9: // V + EF + I 
        for (int i = 0; i < ned; i++)
          cnt[nv+i]= first_edge_dof[i+1]-first_edge_dof[i];
        for (int i = 0; i < nfa; i++)
          {
            ma.GetFaceEdges (i, f2ed);
            int fdof = first_face_dof[i+1]-first_face_dof[i];
            for (int j = 0; j < f2ed.Size(); j++)
              cnt[nv+f2ed[j]] +=  fdof; 
          }
        for (int i = 0; i < ni; i++)
          cnt[nv+ned+i] = first_element_dof[i+1]-first_element_dof[i];
        break;     

      case 10: // V + EI + FI 
        cnt = 0;  
        for (int i =0; i< nv;i++)
          cnt[ma.GetClusterRepVertex(i)]++; 
        for (int i = 0; i < ned; i++)
          if(fine_edge[i])
            cnt[ma.GetClusterRepEdge(i)] += first_edge_dof[i+1]-first_edge_dof[i];
       
        for (int i = 0; i < nfa; i++)
          if(fine_face[i])
            cnt[ma.GetClusterRepFace(i)] += first_face_dof[i+1]-first_face_dof[i];
        for (int i = 0; i < ni; i++)
          {
            int ccl = ma.GetClusterRepElement(i);
            cnt[ccl] +=first_element_dof[i+1] - first_element_dof[i];
          }
        /*
          for (int i = 0; i < ni; i++)
          {
          int ccl = ma.GetClusterRepElement(i);
          ma.GetElEdges (i,ednums,orient);
	     
          for (int k = 0; k < ednums.Size(); k++)
          if(ccl!= ma.GetClusterRepFace(ednums[k]))
          cnt[ma.GetClusterRepEdge(ednums[k])] += first_element_dof[i+1] - first_element_dof[i];
	      
          ma.GetElFaces (i,fanums,orient);
          for (int k = 0; k < fanums.Size(); k++)
          if(ccl != ma.GetClusterRepFace(fanums[k]))
          cnt[ma.GetClusterRepFace(fanums[k])] += first_element_dof[i+1] - first_element_dof[i];
          }*/
        //cout << "cnt " << cnt << endl; 
        //*testout  << "cnt " << cnt << endl; 
        break; 

      case 11: // 2d VEI 
        for (int i =0; i<ned; i++)
          if(fine_edge[i])
            {
              int v1, v2;
              int ndof = first_edge_dof[i+1]-first_edge_dof[i];
              ma.GetEdgePNums (i, v1, v2);
              cnt[v1] += ndof;
              cnt[v2] += ndof;
            }
        for (int i = 0; i < ni; i++)
          {
            Array<int>  pnums; 
            ma.GetElPNums(i,pnums); 
            int ndof = first_element_dof[i+1] - first_element_dof[i];
            for (int j = 0; j < pnums.Size(); j++)
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
              for (int j = 0; j < f2ed.Size(); j++)
                cnt[nv+f2ed[j]] +=  fdof; 
            }
          break;
        }
      case 21: // V + E + F + I 
	int ds_order = precflags.GetNumFlag ("ds_order", 1);
	if (ds_order < 1) ds_order = 1;
        for (int i = 0; i < ned; i++)
	  {
	    cnt[i] = first_edge_dof[i+1]-first_edge_dof[i]-ds_order+1;
	    if (cnt[i] < 0) cnt[i] = 0;
	  }
        for (int i = 0; i < nfa; i++)
          cnt[ned+i] = first_face_dof[i+1]-first_face_dof[i];
        break; 
      }
     
    //    *testout << " cnt " << cnt << endl; 
    
    Table<int> & table = *new Table<int> (cnt); 
    cnt = 0; 
    if(SmoothingType != 10 && SmoothingType != 4 && SmoothingType <= 20) 
      {
        for (int i = 0; i < nv; i++)
          table[i][0] = i;
        /*
        if (augmented == 1)
          for (int i = 0; i < nv; i++)
            table[i][1] = nv+i;
        else if (augmented == 2)
          for (int i = 0; i < nv; i++)
            for (int j = 0; j < order-1; j++)
              table[i][j+1] = nv+i*(order-1)+j; 
        */
        for (int i = 0; i < nv; i++)
          cnt[i] = nvdof;
      }
 

    ii=0; 
    switch(SmoothingType)
      {
      case 1:  // 2d: V + E + I
        for (int i = 0; i < ned; i++,ii++)
          if(fine_edge[i])
            {
              int first = first_edge_dof[i];
              int ndof = first_edge_dof[i+1]-first_edge_dof[i]; 
              for (int j = 0; j < ndof; j++)
                table[nv+i][j] = first+j;
            }
        for (int i = 0; i < nel; i++)
          {
            int first = first_element_dof[i];
            int ndof =  first_element_dof[i+1]-first_element_dof[i]; 
            for (int j = 0; j < ndof; j++)
              table[nv+ned+i][j] = first+j;
          }
        break; 
      case 2: // 2d VE + I
        for (int i = 0; i < ned; i++)
          if(fine_edge[i])
            {
              int v1, v2;
              int first = first_edge_dof[i];
              int ndof = first_edge_dof[i+1]-first;
              ma.GetEdgePNums (i, v1, v2);
              for (int j = 0; j < ndof; j++)
                {
                  table[v1][cnt[v1]++] = first+j;
                  table[v2][cnt[v2]++] = first+j;
                }
            }
        for (int i = 0; i < ni; i++)
          {
            int first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first;
            for (int j = 0; j < ndof; j++)
              table[nv+i][j] = first+j;
          }
      case 3 : // V + E + F + I 
        for (int i = 0; i < ned; i++)
          {
            int ndof = first_edge_dof[i+1]-first_edge_dof[i]; 
            int first = first_edge_dof[i];
            for (int j = 0; j < ndof; j++)
              table[nv+i][j] = first+j;
          }
        for (int i = 0; i < nfa; i++)
          {
            int ndof = first_face_dof[i+1]-first_face_dof[i]; 
            int first = first_face_dof[i];
            for (int j = 0; j < ndof; j++)
              table[nv+ned+i][j] = first+j;
          }
        for (int i = 0; i < ni; i++)
          {
            int first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i]; 
            for (int j = 0; j < ndof; j++)
              table[nv+ned+nfa+i][j] = first+j;
          }
        break; 
      case 4: // VE + F + I
        for (int i = 0; i < nv; i++)
	  if (!IsDirichletVertex(i))
	    {
	      table[i][0] = i;
	      cnt[i] = 1;
	    }

        for (int i = 0; i < ned; i++)
          if (!IsDirichletEdge(i))
            {
              IntRange range = GetEdgeDofs(i);
              Ng_Node<1> edge = ma.GetNode<1> (i);
              int v[2] = { edge.vertices[0], edge.vertices[1] };
              
              for (int j = 0; j < range.Size(); j ++)
                for (int k = 0; k < 2; k++)
                  table[v[k]][cnt[v[k]]++] = range.First()+j;
            }
        for (int i = 0; i < nfa; i++)
	  if (!IsDirichletFace(i))
	    {
	      IntRange range = GetFaceDofs(i);
	      for (int j = 0; j < range.Size(); j++)
		table[nv+i][j] = range.First()+j;
	      cnt[nv+i] = range.Size();
	    }
        for (int i = 0; i < ni; i++)
          {
            int first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i]; 
            for (int j = 0; j < ndof; j++)
              table[nv+nfa+i][j] = first+j;
          }
        break; 
      case 5: // VE + FI
        for (int i = 0; i < ned; i++)
          {
            int v1, v2;
            int first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first;
            ma.GetEdgePNums (i, v1, v2);
            for (int j = 0; j < ndof; j++)
              {
                table[v1][cnt[v1]++] = first+j;
                table[v2][cnt[v2]++] = first+j;
              }
          }
        for (int i = 0; i < nfa; i++)
          {
            int first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first_face_dof[i];
            for (int j = 0; j < ndof; j++)
              table[nv+i][j] = first+j;
            cnt[nv+i] = ndof;
          }
	
        for (int i = 0; i < ni; i++)
          {
            ma.GetElFaces (i, fanums);
            int first = first_element_dof[i];
            int ndof = first_element_dof[i+1] - first_element_dof[i];
            for (int j = 0; j < fanums.Size(); j++)
              for (int k = 0; k < ndof; k++)
                table[nv+fanums[j]][cnt[nv+fanums[j]]++] = first + k;
          }
        break; 
      case 6: // VEF + I
        for (int i = 0; i < ned; i++)
          {
            int v1, v2;
            int first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first;
            ma.GetEdgePNums (i, v1, v2);
            for (int j = 0; j < ndof; j++)
              {
                table[v1][cnt[v1]++] = first+j;
                table[v2][cnt[v2]++] = first+j;
              }
          }
        for (int i = 0; i < nfa; i++)
          {
            Array<int> pnums; 
            ma.GetFacePNums(i,pnums); 
            int first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first_face_dof[i];
            for (int k = 0;k<pnums.Size();k++)
              for (int j = 0; j < ndof; j++)
                table[pnums[k]][cnt[pnums[k]]++] = first+j;
          }
        for (int i = 0; i < ni; i++)
          {
            int first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first;
            for (int j = 0; j < ndof; j++)
              table[nv+i][j] = first+j;
          }
        break; 
      case 7: // VEFI
        for (int i = 0; i < ned; i++)
          {
            int v1, v2;
            int first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first;
            ma.GetEdgePNums (i, v1, v2);
            for (int j = 0; j < ndof; j++)
              {
                table[v1][cnt[v1]++] = first+j;
                table[v2][cnt[v2]++] = first+j;
              }
          }
        for (int i = 0; i < nfa; i++)
          {
            Array<int> pnums; 
            ma.GetFacePNums(i,pnums); 
            int first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first_face_dof[i];
            for (int k = 0;k<pnums.Size();k++)
              for (int j = 0; j < ndof; j++)
                table[pnums[k]][cnt[pnums[k]]++] = first+j;
          }
        for (int i = 0; i < ni; i++)
          {
            Array<int> pnums; 
            ma.GetElPNums(i,pnums); 
            int first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i];
            for (int k = 0;k<pnums.Size();k++)
              for (int j = 0; j < ndof; j++)
                table[pnums[k]][cnt[pnums[k]]++] = first+j;
          }
        break;
      case 12: // VEFI
        //cnt =0; 
        for (int i =0;i<nv;i++)
          table[ma.GetClusterRepVertex(i)][cnt[ma.GetClusterRepVertex(i)]++]=i;
	     
        for (int i = 0; i < ned; i++)
          {
            int v1, v2;
            int first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first;
            ma.GetEdgePNums (i, v1, v2);
            int r1 = ma.GetClusterRepVertex(v1); 
            int r2 = ma.GetClusterRepVertex(v2); 
	    
            for (int j = 0; j < ndof; j++)
              {
                table[r1][cnt[r1]++] = first+j;
                if(r1!=r2)
                  table[v2][cnt[v2]++] = first+j;
              }
          }
        for (int i = 0; i < nfa; i++)
          {
            Array<int> pnums; 
            ma.GetFacePNums(i,pnums); 
	    
            int first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first_face_dof[i];
            Array<int> repv; 
	  
            for (int j =0;j<pnums.Size();j++) 
              repv.Append(ma.GetClusterRepVertex(pnums[j]));
            for (int j =0;j<pnums.Size();j++) 
              { 
                bool ok=1; 
                for (int k = 0;k<j;k++) if(repv[j] == repv[k]) ok = 0;  
		
                if(ok)
                  for (int k = 0; k < ndof; k++)
                    table[repv[j]][cnt[repv[j]]++] += first+k ;
              }
          }
        for (int i = 0; i < ni; i++)
          {
            Array<int> pnums; 
            ma.GetElPNums(i,pnums); 
            int first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i];
	    
	    
            Array<int> repv; 
            for (int j =0;j<pnums.Size();j++) 
              repv.Append( ma.GetClusterRepVertex(pnums[j]));
            for (int j =0;j<pnums.Size();j++) 
              { 
                bool ok=1; 
                for (int k = 0;k<j;k++) if(repv[j] == repv[k]) ok = 0;  
		
                if(ok)
                  for (int k = 0; k < ndof; k++)
                    table[repv[j]][cnt[repv[j]]++] += first+k ;
              }
          }
	   
        break;
 
        for (int i = 0; i < ned; i++)
          {
            int v1, v2;
            int first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first;
            ma.GetEdgePNums (i, v1, v2);
            for (int j = 0; j < ndof; j++)
              {
                table[v1][cnt[v1]++] = first+j;
                table[v2][cnt[v2]++] = first+j;
              }
          }
        for (int i = 0; i < nfa; i++)
          {
            Array<int> pnums; 
            ma.GetFacePNums(i,pnums); 
            int first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first_face_dof[i];
            for (int k = 0;k<pnums.Size();k++)
              for (int j = 0; j < ndof; j++)
                table[pnums[k]][cnt[pnums[k]]++] = first+j;
          }
        for (int i = 0; i < ni; i++)
          {
            Array<int> pnums; 
            ma.GetElPNums(i,pnums); 
            int first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i];
            for (int k = 0;k<pnums.Size();k++)
              for (int j = 0; j < ndof; j++)
                table[pnums[k]][cnt[pnums[k]]++] = first+j;
          }
        break;
 
      case 8: // V + E + FI 
        for (int i = 0; i < ned; i++)
          {
            int first = first_edge_dof[i];
            for (int j = first; j < first_edge_dof[i+1]; j++)
              table[nv+i][cnt[nv+i]++] = j;
          }
        for (int i = 0; i < nfa; i++)
          {
            int first = first_face_dof[i];
            for (int j = first; j < first_face_dof[i+1]; j++)
              table[nv+ned+i][cnt[nv+ned+i]++] = j;
          }
        for (int i = 0; i < nel; i++)
          {
            ma.GetElFaces (i,fanums);
            for (int k = 0; k < fanums.Size(); k++)
              for (int j = first_element_dof[i]; j < first_element_dof[i+1]; j++)
                table[nv+ned+fanums[k]][cnt[nv+ned+fanums[k]]++] = j;
          }
        break;
      case 9: // V + EF + I
        for (int i = 0; i < ned; i++)
          {
            int first = first_edge_dof[i];
            int ndof = first_edge_dof[i+1]-first;
            for (int j = 0; j < ndof; j++)
              table[nv+i][cnt[nv+i]++] = first+j;
          }
        for (int i = 0; i < nfa; i++)
          {
            int first = first_face_dof[i];
            int ndof = first_face_dof[i+1]-first;
            ma.GetFaceEdges (i, f2ed);
            for (int k = 0; k < f2ed.Size(); k++)
              for (int j = 0; j < ndof; j++)
                table[nv+f2ed[k]][cnt[nv+f2ed[k]]++] = first+j;
          }
        for (int i = 0; i < ni; i++)
          {
            int first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first;
            for (int j = 0; j < ndof; j++)
              table[nv+ned+i][j] = first+j;
          }
        break; 
      case 10: // V + EI + FI
        cnt =0; 

        for (int i = 0; i < nv ; i++) 
          table[ma.GetClusterRepVertex(i)][cnt[ma.GetClusterRepVertex(i)]++] = i; 	

	
        for (int i = 0; i < ned; i++)
          if(fine_edge[i])
            for (int j = first_edge_dof[i]; j < first_edge_dof[i+1]; j++)
              table[ma.GetClusterRepEdge(i)][cnt[ma.GetClusterRepEdge(i)]++] = j;
	

        for (int i = 0; i < nfa; i++)
          if(fine_face[i])
            for (int j = first_face_dof[i]; j < first_face_dof[i+1]; j++)
              table[ma.GetClusterRepFace(i)][cnt[ma.GetClusterRepFace(i)]++] = j;

        for (int i =0;i<ni;i++)
          {
            int ccl = ma.GetClusterRepElement(i);
            for (int j = first_element_dof[i]; j < first_element_dof[i+1]; j++)
              table[ma.GetClusterRepElement(i)][cnt[ccl]++] = j;
          }
        /*	for (int i = 0; i < ni; i++)
          {
          int ccl = ma.GetClusterRepElement(i);
          ma.GetElEdges (i,ednums,orient);
          for (int k = 0; k < ednums.Size(); k++)
          //  if(ccl!=ma.GetClusterRepEdge(ednums[k]))
          for (int j = first_element_dof[i]; j < first_element_dof[i+1]; j++)
          table[ma.GetClusterRepEdge(ednums[k])][cnt[ma.GetClusterRepEdge(ednums[k])]++] = j;
	    
          ma.GetElFaces (i,fanums,orient);
	   
          for (int k = 0; k < fanums.Size(); k++)
          //  if(ccl!=ma.GetClusterRepFace(fanums[k]))
          for (int j = first_element_dof[i]; j < first_element_dof[i+1]; j++)
          table[ma.GetClusterRepFace(fanums[k])][cnt[ma.GetClusterRepFace(fanums[k])]++] = j;
          }	*/
        //cout << " test 5 " << endl; 
        //cout << "table " << table << endl; 
        //*testout  << "table " << table << endl; 
        break;
      case 11: //2d VEF
        for (int i = 0; i < ned; i++)
          if(fine_edge[i])
            {
              int v1, v2;
              int first = first_edge_dof[i];
              int ndof = first_edge_dof[i+1]-first;
              ma.GetEdgePNums (i, v1, v2);
              for (int j = 0; j < ndof; j++)
                {
                  table[v1][cnt[v1]++] = first+j;
                  table[v2][cnt[v2]++] = first+j;
                }
            }
        for (int i = 0; i < ni; i++)
          {
            Array<int> pnums; 
            ma.GetElPNums(i,pnums); 
            int first = first_element_dof[i];
            int ndof = first_element_dof[i+1]-first_element_dof[i];
            for (int k=0;k<pnums.Size();k++)
              for (int j = 0; j < ndof; j++)
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
              int first = first_face_dof[i];
              int ndof = first_face_dof[i+1]-first;
              ma.GetFaceEdges (i, f2ed);
              for (int k = 0; k < f2ed.Size(); k++)
                for (int j = 0; j < ndof; j++)
                  table[nv+f2ed[k]][cnt[nv+f2ed[k]]++] = first+j;
            }
          break;
        }
      case 21: // E + F 
	
	int ds_order = precflags.GetNumFlag ("ds_order", 1);
	cout << "ds_order = " << ds_order << endl;
        for (int i = 0; i < ned; i++)
          {
            int first = first_edge_dof[i] + ds_order - 1;
            int ndof = first_edge_dof[i+1]-first;
            for (int j = 0; j < ndof; j++)
              table[i][j] = first+j;
          }
        for (int i = 0; i < nfa; i++)
          {
            int ndof = first_face_dof[i+1]-first_face_dof[i]; 
            int first = first_face_dof[i];
            for (int j = 0; j < ndof; j++)
              table[ned+i][j] = first+j;
          }
        break; 

      }
    
    /*
      (*testout) << "H1HO-table = " << table << endl;
    */
        return &table;
  }


  Array<int> * 
  H1HighOrderFESpace :: CreateDirectSolverClusters (const Flags & flags) const
  {
    if (flags.NumFlagDefined ("ds_order"))
      {
	int ds_order = int (flags.GetNumFlag ("ds_order", 1));

	Array<int> & clusters = *new Array<int> (GetNDof());
	clusters = 0;
	
	int ned = ma.GetNEdges();
	int nv = ma.GetNV();

	for (int i = 0; i < nv; i++)
	  clusters[i] = 1;

	for (int i = 0; i < ned; i++)
	  {
	    int first = first_edge_dof[i];
	    int next = first_edge_dof[i+1];
	    for (int j = 0; (j+2 <= ds_order) && (first+j < next) ; j++)
	      clusters[first+j] = 1;
	  }

	/*
	int nfa = ma.GetNFaces();
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

    int i, j, k;
    
    // int nv = ma.GetNV();
    // int nd = GetNDof();
    int ne = ma.GetNE();
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
                  clusters[k] = 2;
              }
	    
            ma.GetElFaces(i,fnums); // vertical faces 
            for (int j =2;j<5;j++) 
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
            ma.GetElEdges (i, ednums);
            for (j = 8; j < 12; j++) //vertical edges
              {
                int first = first_edge_dof[ednums[j]];
                int next = first_edge_dof[ednums[j]+1];
                for (k = first; k < next; k++)
                  clusters[k] = 2;
              }
            ma.GetElFaces(i,fnums); // vertical faces 
            for (int j =2;j<6;j++) 
              {
                
                int first = first_face_dof[fnums[j]]; 
                int next = first_face_dof[fnums[j]+1]; 
		
                for (k=first; k < next; k++) 
                  clusters[k]=3; 
              }
          } 
      }


   
    for (int i =0; directsolverclustered.Size() > 0 && i<ne; i++)
      {
        if(directsolverclustered[ma.GetElIndex(i)])
          {
            GetDofNrs(i,ednums);
            for (int k = 0; k<ednums.Size(); k++)
              {
                clusters[ednums[k]] = 4;
              }
          }
      }

   

    for (int i =0; i< adddirectsolverdofs.Size(); i++)
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

    for (int i = 0; i<directvertexclusters.Size(); i++)
      if(directvertexclusters[i] >= 0)
        clusters[i] = directvertexclusters[i] + stdoffset;

    for (int i = 0; i<directedgeclusters.Size(); i++)
      if(directedgeclusters[i] >= 0)
        for(j = first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
          clusters[j] = directedgeclusters[i] + stdoffset;

    for (int i = 0; i<directfaceclusters.Size(); i++)
      if(directfaceclusters[i] >= 0)
        for(j = first_face_dof[i]; j<first_face_dof[i+1]; j++)
          clusters[j] = directfaceclusters[i] + stdoffset;
	  
    for (int i = 0; i<directelementclusters.Size(); i++)
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
