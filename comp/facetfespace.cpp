
#include <comp.hpp>
#include <fem.hpp> 

#include <parallelngs.hpp>


namespace ngcomp
{
  using namespace ngfem ;
  using namespace ngparallel;

// ------------------------------------------------------------------------
  FacetFESpace ::  FacetFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags)
  : FESpace(ama, flags)
  {
    name="FacetFESpace(facet)";
    // defined flags
    DefineNumFlag("relorder");
    DefineDefineFlag("print"); 
    DefineDefineFlag("variableorder"); 


    if(parseflags) ParseFlags(flags);
    
    print = flags.GetDefineFlag("print"); 
    
    ndlevel.SetSize(0);
    Flags loflags;
    loflags.SetFlag("order",0.0);
    if ( this->IsComplex() )
      loflags.SetFlag("complex");
    if (order != 0)
      low_order_space = new FacetFESpace(ma, loflags);
    else
      low_order_space = 0;
        
// #ifdef PARALLEL
//     low_order_space = 0; // new FacesFESpace
// #endif

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
	  cerr << " WARNING: H1HoFeSpace: inconsistent flags: variableorder, order and relorder "
	       << "-> variable order space with rel_order " << rel_order << "is used, but order is ignored " << endl; 
	else 
	  cerr << " WARNING: H1HoFeSpace: inconsistent flags: order and rel_order "
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
    
    // TODO: Evaluator for shape tester 
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
    // Update();
  }
  
// ------------------------------------------------------------------------
  FacetFESpace :: ~FacetFESpace ()
  {
    ;
  }

// ------------------------------------------------------------------------
  FESpace * FacetFESpace :: Create (const MeshAccess & ma, const Flags & flags)
  {
    return new FacetFESpace (ma, flags, true);
  }

// ------------------------------------------------------------------------
  void FacetFESpace :: Update(LocalHeap & lh)
  {
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
    
    ARRAY<int> fanums;
        
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
	    ARRAY<int> elfaces,vnums;
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
      }
      first_facet_dof[nfa] = ndof;
    } // 2D
    else  // 3D
    {
      int inci = 0;
      INT<2> p; 
      ARRAY<int> pnums;
      for (int i=0; i< nfa; i++)
      {
        p = order_facet[i];
        ma.GetFacePNums(i,pnums);

        switch(pnums.Size())
        {
          case 3: 
            inci = ((p[0]+1)*(p[0]+2))/2 - 1;
            break;
          case 4: 
            inci= (p[0]+1)*(p[1]+1) - 1;
            break;
        }
        first_facet_dof[i] = ndof;
        ndof+= inci;
      }
      first_facet_dof[nfa] = ndof;
    } // 3D

      //HERBERT: h1hofe code. is that correct ??
    while (ma.GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;
      
      //HERBERT: no prolongation so far       
      //prol->Update();
    if(print)
      {
	*testout << "*** Update FAcetFESpace: General Information" << endl;
	*testout << " order facet (facet) " << order_facet << endl; 
	*testout << " first_facet_dof (facet)  " << first_facet_dof << endl; 
      } 
   /* 
    cout << "    ne = " << nel << ", nfa = " << nfa << ", ncfa=" << ncfa << ", ndof=" << ndof << endl;
    cout << "    order_facet = " << order_facet << endl;*/




    for (int i = 0; i < 4; i++) 
      {
	lodofs_per_node[i] = 0;
	first_lodof[i] = 0;
      }

    lodofs_per_node[ma.GetDimension()-1] = 1;
    for (int i = ma.GetDimension(); i < 5; i++) 
      first_lodof[i] = nfa;

    for (int i = 0; i < 4; i++)
      first_hodofs[i].SetSize(0);
    first_hodofs[ma.GetDimension()-1] = first_facet_dof;


#ifdef PARALLEL
    *testout << "update parallel dofs in facet-fespace, ndof " << ndof << endl;
    UpdateParallelDofs();
#endif

  }



// ------------------------------------------------------------------------
  const FiniteElement & FacetFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    FacetVolumeFiniteElement<2> * fe2d;
    FacetVolumeFiniteElement<3> * fe3d;

    switch (ma.GetElType(elnr))
    {
      case ET_TRIG:
        fe2d = new (lh.Alloc (sizeof(FacetVolumeTrig)))  FacetVolumeTrig ();
        break;
      case ET_QUAD:
        fe2d = new (lh.Alloc (sizeof(FacetVolumeQuad)))  FacetVolumeQuad ();
        break;
      case ET_TET:
        fe3d = new (lh.Alloc (sizeof(FacetVolumeTet)))  FacetVolumeTet ();
        break;
      case ET_PYRAMID:
        fe3d = new (lh.Alloc (sizeof(FacetVolumePyramid)))  FacetVolumePyramid ();
        break;
      case ET_PRISM:
        fe3d = new (lh.Alloc (sizeof(FacetVolumePrism)))  FacetVolumePrism ();
        break;
      case ET_HEX:
        fe3d = new (lh.Alloc (sizeof(FacetVolumeHex)))  FacetVolumeHex ();
        break;
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

    ARRAY<int> vnums;
    ArrayMem<int, 6> fanums, order_fa;
    
    ma.GetElVertices(elnr, vnums);

    if (ma.GetDimension() == 2)
      ma.GetElEdges(elnr, fanums);
    else
      ma.GetElFaces(elnr, fanums);
    
    order_fa.SetSize(fanums.Size());
    for (int j = 0; j < fanums.Size(); j++)
      order_fa[j] = order_facet[fanums[j]][0]; //SZ not yet anisotropric
    
#ifdef PARALLEL
    if ( ntasks > 1 )
      for ( int i = 0; i < vnums.Size(); i++ )
        vnums[i] = parallelma->GetDistantPNum(0, vnums[i]);
#endif
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

    FacetFacetFiniteElement<1> * fe1d = 0;
    FacetFacetFiniteElement<2> * fe2d = 0;

    switch (ma.GetSElType(selnr))
    {
      case ET_SEGM:
        fe1d = new (lh.Alloc (sizeof(FacetFacetSegm)))  FacetFacetSegm ();
        break;
      case ET_TRIG:
        fe2d = new (lh.Alloc (sizeof(FacetFacetTrig)))  FacetFacetTrig ();
        break;
      case ET_QUAD:
        fe2d = new (lh.Alloc (sizeof(FacetFacetQuad)))  FacetFacetQuad ();
        break;
      default:
        ;
    }
     
    if (!fe1d && !fe2d)
    {
      stringstream str;
      str << "FacetFESpace " << GetClassName()
          << ", undefined eltype "
          << ElementTopology::GetElementName(ma.GetSElType(selnr))
          << ", order = " << order << endl;
      throw Exception (str.str());
    }
     
    ArrayMem<int,4> vnums;
    ArrayMem<int, 4> ednums;
    
    ma.GetSElVertices(selnr, vnums);
    switch (ma.GetSElType(selnr))
    {
      case ET_SEGM:
#ifdef PARALLEL
	if ( ntasks > 1 )
	  for ( int i = 0; i < vnums.Size(); i++ )
	    vnums[i] = parallelma->GetDistantPNum(0, vnums[i]);
#endif
        fe1d -> SetVertexNumbers (vnums);
        ma.GetSElEdges(selnr, ednums);
        fe1d -> SetOrder (order_facet[ednums[0]][0]); 
        fe1d -> ComputeNDof();
        return *fe1d;
        break;
      case ET_TRIG: 
      case ET_QUAD:
#ifdef PARALLEL
	if ( ntasks > 1 )
	  for ( int i = 0; i < vnums.Size(); i++ )
	    vnums[i] = parallelma->GetDistantPNum(0, vnums[i]);
#endif
        fe2d -> SetVertexNumbers (vnums);
        fe2d -> SetOrder (order_facet[ma.GetSElFace(selnr)][0]);// SZ not yet anisotropic order for facet fe !!! 
        fe2d -> ComputeNDof();
        return *fe2d;
        break;
    }
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
  void FacetFESpace :: GetDofNrs (int elnr, ARRAY<int> & dnums) const
  {

    ARRAY<int> fanums; // facet numbers
    int first,next;

    fanums.SetSize(0);
    dnums.SetSize(0);
    

    if(ma.GetDimension() == 3)
      ma.GetElFaces (elnr, fanums);
    else // dim=2
      ma.GetElEdges (elnr, fanums);

    for(int i=0; i<fanums.Size(); i++)
      {
        dnums.Append(fanums[i]); // low_order
        first = first_facet_dof[fanums[i]];
        next = first_facet_dof[fanums[i]+1];
        for(int j=first ; j<next; j++)
        dnums.Append(j);
      }

    
    // old 
    /*
    if(ma.GetDimension() == 3)
    {
      ma.GetElFaces (elnr, fdnums);
      for(int i=0; i<fdnums.Size(); i++)
      {
        for (int j=0; j<fdnums.Size(); j++
        first = first_facet_dof[fdnums[i]];
        next = first_facet_dof[fdnums[i]+1];
        for(int j=first ; j<next; j++)
          dnums.Append(j);
      }
    }
    else // 2D
    {
      ma.GetElEdges (elnr, ednums); //neu
      for (int i = 0; i < ednums.Size(); i++)
      {
        first = first_facet_dof[ednums[i]];
        next = first_facet_dof[ednums[i]+1];

        for (int j = first; j <next; j++)
          dnums.Append (j);
      }
    }
    */
//     cout << "** GetDofNrs(" << elnr << "): " << dnums << endl;
  }

  void FacetFESpace :: GetWireBasketDofNrs (int elnr, ARRAY<int> & dnums) const
  {
    ArrayMem<int,12> facets;

    dnums.SetSize(0);

    ma.GetElFacets (elnr, facets);

    for (int i = 0; i < facets.Size(); i++)
      {
	dnums.Append (facets[i]);
      }
  }

// ------------------------------------------------------------------------
//   void FacetFESpace :: GetExternalDofNrs (int elnr, ARRAY<int> & dnums) const
//   {
//     GetDofNrs(elnr, dnums);
//   }


// ------------------------------------------------------------------------
  void FacetFESpace :: GetSDofNrs (int selnr, ARRAY<int> & dnums) const
  {
    int first,next;

    dnums.SetSize(0);
    if (ma.GetDimension() == 2)
    {
      ArrayMem<int, 4> fanums;
      ma.GetSElEdges (selnr, fanums);

      dnums.Append(fanums[0]);
      
      first = first_facet_dof[fanums[0]];
      next = first_facet_dof[fanums[0]+1];
      for (int j = first; j <next; j++)
        dnums.Append (j);
      
/*      for (int i = 0; i < fanums.Size(); i++)
      {
        dnums.Append(fanums[i]);
      
        first = first_facet_dof[ednums[i]];
        next = first_facet_dof[ednums[i]+1];
        for (int j = first; j <next; j++)
          dnums.Append (j);
      }*/
    }
    else// 3D
    {
        int fnum = ma.GetSElFace(selnr);
        dnums.Append(fnum);
        first = first_facet_dof[fnum];
        next = first_facet_dof[fnum+1];
        for(int j=first; j<next; j++)
          dnums.Append(j);
      }

/*      if(first_facet_dof.Size()>1)
      {
        int fnum = ma.GetSElFace(selnr);
        first = first_facet_dof[fnum];
        next = first_facet_dof[fnum+1];
        for(int j=first; j<next; j++)
          dnums.Append(j);
      }*/
    }



// ------------------------------------------------------------------------
  Table<int> * FacetFESpace :: CreateSmoothingBlocks (const Flags & precflags) const
  {
    int ncnt;

    // 1 x low order + faces/edges
    ncnt = nfa-ncfa;
      
    ARRAY<int> cnt(ncnt);
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
      
    cout << "smoothingblocks = " << endl << table << endl;
    return &table;

  }


  
  ARRAY<int> * FacetFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    ARRAY<int> & clusters = *new ARRAY<int> (GetNDof());

    clusters.SetSize(ndof);
    clusters = 0;
    
    return &clusters;
    
    //
  
    for (int i=0; i<nfa-ncfa; i++)
      clusters[i]  = 1;
  
    cout << "direct solver cluster = " << clusters << endl;
    return & clusters;
  }



#ifdef PARALLEL_NOT_JS

void FacetFESpace :: UpdateParallelDofs_hoproc()
  {
     // ******************************
     // update exchange dofs 
     // ******************************
    *testout << "Facet::UpdateParallelDofs_hoproc" << endl;
    // Find number of exchange dofs
    ARRAY<int> nexdof(ntasks);
    nexdof = 0;

    MPI_Status status;
    MPI_Request * sendrequest = new MPI_Request[ntasks];
    MPI_Request * recvrequest = new MPI_Request[ntasks];

    // number of face exchange dofs
    for ( int face = 0; face < ma.GetNFaces(); face++ )
      {
	if ( !parallelma->IsExchangeFace ( face ) ) continue;
	
	ARRAY<int> dnums;
	GetFaceDofNrs ( face, dnums );
	nexdof[id] += dnums.Size() ; 

	for ( int dest = 1; dest < ntasks; dest ++ )
	  if (  parallelma -> GetDistantFaceNum ( dest, face ) >= 0 )
	    nexdof[dest] += dnums.Size() ; 
      }

    // + number of inner exchange dofs --> set eliminate internal

    nexdof[0] = ma.GetNFaces();

    paralleldofs->SetNExDof(nexdof);

//     paralleldofs->localexchangedof = new Table<int> (nexdof);
//     paralleldofs->distantexchangedof = new Table<int> (nexdof);
    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);

    ARRAY<int> ** owndofs, ** distantdofs;
    owndofs = new ARRAY<int>* [ntasks];
    distantdofs = new ARRAY<int>* [ntasks];

    for ( int i = 0; i < ntasks; i++ )
      {
	owndofs[i] = new ARRAY<int>(1);
	(*owndofs[i])[0] = ndof;
	distantdofs[i] = new ARRAY<int>(0);
      }

    ARRAY<int> cnt_nexdof(ntasks);
    cnt_nexdof = 0;
    int exdof = 0;
    int ii;

    // *****************
    // Parallel Face dofs
    // *****************


   for ( int face = 0; face < ma.GetNFaces(); face++ )
      if ( parallelma->IsExchangeFace ( face ) )
      {
	ARRAY<int> dnums;
	GetFaceDofNrs ( face, dnums );
	if ( dnums.Size() == 0 ) continue;

 	for ( int i=0; i<dnums.Size(); i++ )
 	  (*(paralleldofs->sorted_exchangedof))[id][exdof++] = dnums[i];

	for ( int dest = 1; dest < ntasks; dest++ )
	  {
	    int distface = parallelma -> GetDistantFaceNum ( dest, face );
	    if( distface < 0 ) continue;

	    owndofs[dest]->Append ( distface );
	    owndofs[dest]->Append (int(paralleldofs->IsGhostDof(dnums[0])) );
	    for ( int i=0; i<dnums.Size(); i++)
	      {
		paralleldofs->SetExchangeDof ( dest, dnums[i] );
		paralleldofs->SetExchangeDof ( dnums[i] );
		owndofs[dest]->Append ( dnums[i] );
	      }
  	  }
      }   


    for ( int sender = 1; sender < ntasks; sender ++ )
      {
        if ( id == sender )
          for ( int dest = 1; dest < ntasks; dest ++ )
            if ( dest != id )
              {
		MyMPI_ISend ( *owndofs[dest], dest, sendrequest[dest]);
              }
	  
        if ( id != sender )
          {
	    MyMPI_IRecv ( *distantdofs[sender], sender, recvrequest[sender]);
          }
 	  
	  
      }

     for( int dest=1; dest<ntasks; dest++) 
      {
	if ( dest == id ) continue;
	MPI_Wait(recvrequest+dest, &status);
	paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
	// low order raviart thomas dofs first
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int fanum = (*distantdofs[dest])[ii++];
	    int isdistghost = (*distantdofs[dest])[ii++];
	    ARRAY<int> dnums;
	    GetFaceDofNrs (fanum, dnums);
// 	    (*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = dnums[0];
// 	    (*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];
	    (*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = dnums[0];
	    if ( dest < id && !isdistghost )
	      paralleldofs->ismasterdof.Clear ( dnums[0] ) ;
	    ii += dnums.Size();
	    cnt_nexdof[dest]++;
	  }
	// then the high order dofs, without raviart thomas
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int fanum = (*distantdofs[dest])[ii++];
	    int isdistghost = (*distantdofs[dest])[ii++];
	    ARRAY<int> dnums;
	    GetFaceDofNrs (fanum, dnums);
	    ii++; 
	    for ( int i=1; i<dnums.Size(); i++)
	      {
// 		(*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
// 		(*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];
		(*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = dnums[i];
		if ( dest < id && !isdistghost )
		  paralleldofs->ismasterdof.Clear ( dnums[i] ) ;
		ii++; cnt_nexdof[dest]++;
	      }
	  }
      }



     /*******************************

         update low order space

     *****************************/

     if ( order == 0 ) return;

     int ndof_lo = ma.GetNFaces();

     // all dofs are exchange dofs
     nexdof = ndof_lo;
     
     exdof = 0;
     cnt_nexdof = 0;
     

     // *****************
     // Parallel Face dofs
     // *****************
     
     owndofs[0]->SetSize(1);
     (*owndofs[0])[0] = ndof;
     distantdofs[0]->SetSize(0);
     
     // find local and distant dof numbers for face exchange dofs
     for ( int face = 0; face < ma.GetNFaces(); face++ )
       {
	 int dest = 0;
	 
	 int distface = parallelma -> GetDistantFaceNum ( dest, face );
	owndofs[0]->Append ( distface );
	paralleldofs->SetExchangeDof ( dest, face );
	owndofs[0]->Append ( face );
      }   
    
    int dest = 0;
    MyMPI_ISend ( *owndofs[0], dest, sendrequest[dest]);
    MyMPI_IRecv ( *distantdofs[0], dest, recvrequest[dest]);
   
    MPI_Wait ( recvrequest+dest, &status);

    paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
    ii = 1;
    while ( ii < distantdofs[0]->Size() )
      {
	int fanum = (*distantdofs[0])[ii++];
// 	(*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = fanum;
// 	(*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[0])[ii];
	(*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = fanum;
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

#ifdef PARALLEL

  void FacetFESpace :: UpdateParallelDofs_loproc()
  {
    *testout << "FacetFESpace::UpdateParallelDofs_loproc" << endl;

    const MeshAccess & ma = (*this). GetMeshAccess();

    int ndof = GetNDof();

    // Find number of exchange dofs
    ARRAY<int> nexdof(ntasks); 
    nexdof = 0;

    MPI_Status status;
    MPI_Request * sendrequest = new MPI_Request[ntasks];
    MPI_Request * recvrequest = new MPI_Request[ntasks];

    // number of face exchange dofs
    for ( int face = 0; face < ma.GetNFaces(); face++ )
      {
	nexdof[id] ++;//= dnums.Size() ; 
	for ( int dest = 1; dest < ntasks; dest ++ )
	  if (  parallelma -> GetDistantFaceNum ( dest, face ) >= 0 )
	    { 
	      nexdof[dest] ++; 
	    }
      }

    paralleldofs->SetNExDof(nexdof);

//     paralleldofs->localexchangedof = new Table<int> (nexdof);
//     paralleldofs->distantexchangedof = new Table<int> (nexdof);
    paralleldofs->sorted_exchangedof = new Table<int> (nexdof);



    ARRAY<int> ** owndofs,** distantdofs;
    owndofs = new ARRAY<int> * [ntasks];
    distantdofs = new ARRAY<int> * [ntasks];

    for ( int i = 0; i < ntasks; i++)
      {
	owndofs[i] = new ARRAY<int> (1);
	(*owndofs[i])[0] = ndof;
	distantdofs[i] = new ARRAY<int> (0);
      }

    int exdof = 0;
    ARRAY<int> cnt_nexdof(ntasks);
    cnt_nexdof = 0;

    // *****************
    // Parallel Face dofs
    // *****************


    // find local and distant dof numbers for face exchange dofs
    for ( int face = 0; face < ma.GetNFaces(); face++ )
      {
	(*(paralleldofs->sorted_exchangedof))[id][exdof++] = face;

	for ( int dest = 1; dest < ntasks; dest++ )
	  {
	    int distface = parallelma -> GetDistantFaceNum ( dest, face );
	    if( distface < 0 ) continue;

	    owndofs[dest]->Append ( distface );
	    paralleldofs->SetExchangeDof ( dest, face );
	    paralleldofs->SetExchangeDof ( face );
	    owndofs[dest]->Append ( face );
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
	MPI_Wait ( recvrequest+dest, &status );
	paralleldofs -> SetDistNDof( dest, (*distantdofs[dest])[0]) ;
	ii = 1;
	while ( ii < distantdofs[dest]->Size() )
	  {
	    int fanum = (*distantdofs[dest])[ii++];
// 	    (*(paralleldofs->localexchangedof))[dest][ cnt_nexdof[dest] ] = fanum;
// 	    (*(paralleldofs->distantexchangedof))[dest][ cnt_nexdof[dest] ] = (*distantdofs[dest])[ii];
	    (*(paralleldofs->sorted_exchangedof))[dest][ cnt_nexdof[dest] ] = fanum;
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



  
// ------------------------------------------------------------------------
   // register FESpaces
  namespace facefespace_cpp
  {
    class Init
    {
      public:
        Init ();
    };

    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("facet", FacetFESpace::Create);
    }

    Init init;
  }

}


