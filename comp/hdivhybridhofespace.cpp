/**
   High Order Finite Element Space for H(Div)
*/

#include <comp.hpp>
#include <fem.hpp> 

namespace ngcomp
{
  using namespace ngcomp; 

  HDivHybridHighOrderFESpace ::  
  HDivHybridHighOrderFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags)
  : FESpace (ama, flags)
  {
    name="HDivHybridHighOrderFESpace(hdivhybridho)";
    // defined flags
    DefineDefineFlag("hdivhybridho");
    DefineNumFlag("relorder");
    
    if(parseflags) ParseFlags(flags);
    
    //dimension = 1;
    low_order_space = 0; //new RaviartThomasFESpace (ma,dimension, iscomplex);


    rel_order = int (flags.GetNumFlag ("relorder", 0));

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
    // typedef TrigExtensionMonomial T_ORTHOPOL;

    segm = new HDivHighOrderNormalSegm<T_ORTHOPOL>(order);
    if (ma.GetDimension() == 2)
      {
	trig = new HDivHighOrderFE<ET_TRIG>(order);
	quad = new HDivHighOrderFE<ET_QUAD>(order);
      }
    else
      {
	trig = new HDivHighOrderNormalTrig<T_ORTHOPOL>(order);
	quad = new HDivHighOrderNormalQuad<T_ORTHOPOL>(order);
      }   
    //quad = new HDivHighOrderQuad (order);
    tet = new HDivHighOrderFE<ET_TET>(order);
    hex = new HDivHighOrderFE<ET_HEX>(order);
    prism = new HDivHighOrderFE<ET_PRISM>(order);
    
    // Evaluator for shape tester 
    if (ma.GetDimension() == 2)
      {
	Array<CoefficientFunction*> coeffs(1);
	coeffs[0] = new ConstantCoefficientFunction(1);
	evaluator = GetIntegrators().CreateBFI("masshdiv", 2, coeffs);

      }
    else
      {
        Array<CoefficientFunction*> coeffs(1);
	coeffs[0] = new ConstantCoefficientFunction(1);
	evaluator = GetIntegrators().CreateBFI("masshdiv", 3, coeffs);
	boundary_evaluator = GetIntegrators().CreateBFI("robinhdiv", 3, coeffs);

      }
     if (dimension > 1)
      {
	evaluator = new BlockBilinearFormIntegrator (*evaluator, dimension);
	boundary_evaluator =
	  new BlockBilinearFormIntegrator (*boundary_evaluator, dimension);
      }
  }
  
  HDivHybridHighOrderFESpace :: ~HDivHybridHighOrderFESpace ()
  {
    ;
  }

  FESpace * HDivHybridHighOrderFESpace ::
  Create (const MeshAccess & ma, const Flags & flags)
  {
    int order = int(flags.GetNumFlag ("order", 0));

    if (order < 0)
      return new RaviartThomasFESpace (ma, flags, true);
    // Space with MG
    else
      return new   HDivHybridHighOrderFESpace (ma, flags, true);

  }

  void HDivHybridHighOrderFESpace :: Update()
  {
    cout<<"HDivHybridFeSpace"<<endl;
    int i, j;

    if (low_order_space)
      low_order_space -> Update();

    nv = ma.GetNV();
    ned = ma.GetNEdges();
    nel = ma.GetNE();
    nfa = ma.GetNFaces();
   // (*testout)<<"nel= "<<nel<<" nv= "<<nv<<" ned= "<<ned<<" nfa= "<<nfa<<endl;
	
    int nedel; // no of edges per element

    if (ma.GetDimension() == 2)
      {
        order_edge.SetSize (ned);
        order_inner.SetSize (nel);

        order_edge = order;
        order_inner = order;

        Array<int> eledges;
  
        for (i = 0; i < nel; i++)

	  {
	    int elorder = ma.GetElOrder(i) + rel_order; // start with 0

	    ma.GetElEdges (i, eledges);

	    for (j = 0; j < eledges.Size(); j++)
	      if (elorder > order_edge[eledges[j]])
		order_edge[eledges[j]] = elorder;

	    if (elorder > order_inner[i])
	      order_inner[i] = elorder;
	  }


	(*testout) << "H(div)ho: order_edge = " << endl << order_edge << endl;

        ndof = ned; // Thomas-Raviart
	order_edge = order;
	//	order_inner = 1;

        if(1 || order>0)
	  {
	    first_edge_dof.SetSize(ned+1); // last element = ndof;
	    first_edge_dof = 0;

	    first_inner_dof.SetSize(nel+1);
	    int inci = 0;
	    for(i=0; i< nel; i++)
	      {
		switch(ma.GetElType(i))
		  {
		  case ET_TRIG:
		    inci = (order_inner[i]+1)*(order_inner[i]+2);
		    break;
		  case ET_QUAD:
		    inci = 2*(order_inner[i]*order_inner[i]+ order_inner[i]);
                    break;
		  default:
		    ;
		  }
		if (inci < 0) inci = 0;

		first_inner_dof[i] = ndof;
		ndof+=inci;
	      }
	    first_inner_dof[nel] = ndof;

	  }
        else
	  {
	    first_inner_dof.SetSize(1);
	    first_inner_dof[0] = ndof;
	    first_edge_dof.SetSize(1);
	    first_edge_dof[0] =ndof;
	  }
      }
    else
      {
	int p;
	order_face.SetSize (nfa);
	order_inner.SetSize (nel);

	order_face = order;
	order_inner = order;

        Array<int> pnums;
	Array<int> elfaces;

	ndof = nfa;

	for (i = 0; i < nel; i++)

	  {
	    int elorder = ma.GetElOrder(i) + rel_order; // start with 0

	    (*testout)<<"GetElOrder = "<<i<<"-> "<<ma.GetElOrder(i)<<endl;
	    ma.GetElFaces (i, elfaces);

	    for (j = 0; j < elfaces.Size(); j++)
	      if (elorder > order_face[elfaces[j]])
		order_face[elfaces[j]] = elorder;

	    if (elorder > order_inner[i])
	      order_inner[i] = elorder;

	  }

        if(1 || order>0)
	  {
            // face_dof = edge_face_normal_dof + face_normal_dof;
            first_face_dof.SetSize(nfa+1);
	    int inci = 0;
	    for (i=0; i< nfa; i++)
	      {
	        p = order_face[i];
		ma.GetFacePNums(i,pnums);
		//(*testout)<<"Facenummer "<<i<<"="<<pnums<<endl;
		switch(pnums.Size())
		  {
		  case 3: //Triangle
		    inci= (p*p+3*p)/2;//3*p+ (p-1)*(p-2)/2*p;
		    break;
		  case 4: //Quad
		    inci= p*p+2*p;
		    break;
		  }
		if (inci < 0) inci = 0;
		first_face_dof[i] = ndof;

		ndof+= inci;


	      }
	    first_face_dof[nfa] = ndof;



	    first_inner_dof.SetSize(nel+1);
	    for (i=0; i< nel; i++)
	      {
		p = order_inner[i];
		//(*testout)<<"p update"<<p<<endl;
		//(*testout)<<"ElType="<<ma.GetElType(i)<<endl;
		switch(ma.GetElType(i))
		  {
		  case ET_TET:
		    if (p >= 2)
		     inci = 6*(p-1) + 4*(p-1)*(p-2) + (p-1)*(p-2)*(p-3)/2;

		   else inci = 0;
		    break;
		  case ET_PRISM:
		    inci = (p+1)*(3*(p-1)+(p-1)*(p-2))+ (p-1)*(p+1)*(p+2)/2;
		    break;
		  case ET_HEX:
		    inci = 3*p*(p+1)*(p+1);
		    break;
		  default:
		    inci = 0;
		    break;
		  }
		if (inci < 0) inci = 0;

		first_inner_dof[i] = ndof;

		ndof+= inci;

	      }



	    first_inner_dof[first_inner_dof.Size()-1] = ndof;

	  }
	else
	  {
	    first_face_dof.SetSize(1);
	    first_face_dof[0] = ndof;
	    first_inner_dof.SetSize(1);
	    first_inner_dof[0] = ndof;
	  }


      }

  /*  (*testout) << "ndof hdivhofespace update = " << endl << ndof << endl;
    (*testout) << "first_edge_dof = " << endl << first_edge_dof << endl;
    (*testout) << "first_face_dof = " << endl << first_face_dof << endl;
    (*testout) << "first_inner_dof = " << endl << first_inner_dof << endl;*/


        
    while (ma.GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;
    //    prol->Update();
  }











  const FiniteElement & HDivHybridHighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    FiniteElement * fe;
    
    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
    // typedef TrigExtensionMonomial T_ORTHOPOL;
    
    switch (ma.GetElType(elnr))
      {
      case ET_TET:
	{ 
	  fe = new (lh.Alloc (sizeof(HDivHighOrderFE<ET_TET>)))  HDivHighOrderFE<ET_TET> (order);
	  break;
	}
	/*
      case ET_PYRAMID:
	{
	  fe = new (lh.Alloc (sizeof(HDivHighOrderPyramid<T_ORTHOPOL>)))  HDivHighOrderPyramid<T_ORTHOPOL> (order);
	  break;
	}
	*/
      case ET_PRISM:
	{
	  fe = new (lh.Alloc (sizeof(HDivHighOrderFE<ET_PRISM>)))  HDivHighOrderFE<ET_PRISM> (order);
	  break;
	}
      case ET_HEX:
	{
	  fe = new (lh.Alloc (sizeof(HDivHighOrderFE<ET_HEX>)))  HDivHighOrderFE<ET_HEX> (order);
	  break;
	}
      case ET_TRIG:
	{ 
	  fe = new (lh.Alloc (sizeof(HDivHighOrderFE<ET_TRIG>)))  HDivHighOrderFE<ET_TRIG> (order);
	  break;
	}
      case ET_QUAD:
	{
	  fe = new (lh.Alloc (sizeof(HDivHighOrderFE<ET_QUAD>)))  HDivHighOrderFE<ET_QUAD> (order);
	  break;
	}
      default:
	fe = 0; 
      }
  
    if (!fe)
      {
	stringstream str;
	str << "HDivHybridHighOrderFESpace " << GetClassName() 
	    << ", undefined eltype " 
	    << ElementTopology::GetElementName(ma.GetElType(elnr))
	    << ", order = " << order << endl;
	throw Exception (str.str());
      }
    


    Array<int> vnums;
    ma.GetElVertices(elnr, vnums);


    if ( ma.GetDimension() == 2)
      {
	HDivHighOrderFiniteElement<2> * hofe =
	  dynamic_cast<HDivHighOrderFiniteElement<2>*> (fe);
	ArrayMem<int, 12> ednums, order_ed;
	
	ma.GetElEdges(elnr, ednums);
	order_ed.SetSize (ednums.Size());
	
	for (int j = 0; j < ednums.Size(); j++)
	  order_ed[j] = order_edge[ednums[j]];
      
	hofe -> SetVertexNumbers (vnums);
	hofe -> SetOrderEdge (order_ed);
	hofe -> SetOrderInner (order_inner[elnr]);
	hofe -> ComputeNDof();
      }
    else
      {
	
	HDivHighOrderFiniteElement<3> * hofe =
	  dynamic_cast<HDivHighOrderFiniteElement<3>*> (fe);
	
	ArrayMem<int, 6> fanums, order_fa;
	
	ma.GetElFaces(elnr, fanums);
	
	order_fa.SetSize (fanums.Size());
	
	for (int j = 0; j < fanums.Size(); j++)
	  order_fa[j] = order_face[fanums[j]];
	
	
	hofe -> SetVertexNumbers (vnums);
	hofe -> SetOrderFace (order_fa);
	hofe -> SetOrderInner (order_inner[elnr]);
	hofe -> ComputeNDof();
	
      }
    return *fe;
  }


  const FiniteElement & HDivHybridHighOrderFESpace :: GetSFE (int selnr, LocalHeap & lh) const
  {
    int i, j;

    FiniteElement * fe = 0;

    switch (ma.GetSElType(selnr))
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

    if (!fe)
      {
	stringstream str;
	str << "HDivHighOrderFESpace " << GetClassName()
	    << ", undefined eltype "
	    << ElementTopology::GetElementName(ma.GetSElType(selnr))
	    << ", order = " << order << endl;
	throw Exception (str.str());
      }


    ArrayMem<int,4> vnums;
    ArrayMem<int, 4> ednums, order_ed;

    int order_fa;

    ma.GetSElVertices(selnr, vnums);


    if(fe == segm)
      {
	HDivHighOrderNormalFiniteElement<1> * hofe =
	  dynamic_cast<HDivHighOrderNormalFiniteElement<1>*> (fe);

	hofe -> SetVertexNumbers (vnums);


	ma.GetSElEdges(selnr, ednums);

	/*(*testout)<<"selnr="<<selnr<<endl;
       (*testout)<<"ednums="<<ednums<<endl;
       (*testout)<<"vnums="<<vnums<<endl;*/

	hofe -> SetOrderInner (order_edge[ednums[0]]);
	hofe -> ComputeNDof();
      }
      
    else
      {

     // ma.GetSElEdges(selnr, ednums);
      order_fa = order_face[ma.GetSElFace(selnr)];

      HDivHighOrderNormalFiniteElement<2> * hofe =
      dynamic_cast<HDivHighOrderNormalFiniteElement<2>*> (fe);

      hofe -> SetVertexNumbers (vnums);
      hofe -> SetOrderInner (order_fa);
      hofe -> ComputeNDof();


      }

    return *fe;
  }



  int HDivHybridHighOrderFESpace :: GetNDof () const
  {
    
    return ndof;
  }

  int HDivHybridHighOrderFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }




  void HDivHybridHighOrderFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {

    Array<int> ednums, fnums;
    int i, j;
    int first,next;

   // ma.GetElEdges (elnr, ednums);
    if(ma.GetDimension() == 3)
    {
      ma.GetElFaces (elnr, fnums);
      ednums.SetSize(0);//neu
      dnums.SetSize(0);//neu
    }
    else
    {
      ma.GetElEdges (elnr, ednums); //neu
      fnums.SetSize(0);
      dnums.SetSize(0);
    }  



    if(order < 0)
      throw Exception(" HDivHybridHighOrderFESpace :: GetDofNrs() order < 0 ");


    if(ma.GetDimension() == 3)
      {

	//Ravier-Thomas
	for (i = 0; i < fnums.Size(); i++)
	  dnums.Append (fnums[i]);

	
	// faces
	for(i=0; i<fnums.Size(); i++)
	  {
	    first = first_face_dof[fnums[i]];
	    next = first_face_dof[fnums[i]+1];
	    for(j=first ; j<next; j++)
	      dnums.Append(j);
	  }
	
	//inner
	first = first_inner_dof[elnr];
	next = first_inner_dof[elnr+1];
	for(j=first; j<next; j++)
	  dnums.Append(j);

	
      }
    else
      {
	//Ravier-Thomas
	for (i = 0; i < ednums.Size(); i++)
	  dnums.Append (ednums[i]);
	
	//edges
	for (i = 0; i < ednums.Size(); i++)
	  {
	    first = first_edge_dof[ednums[i]];
	    next = first_edge_dof[ednums[i]+1];
	    
	    for (j = first; j <next; j++)
	      dnums.Append (j);
	  }
       //inner
	first = first_inner_dof[elnr];
	next = first_inner_dof[elnr+1];
	for(j=first; j<next; j++)
	  dnums.Append(j);
	
	
      }



       // (*testout) << "el " << elnr << " has dofs " << dnums << endl;
  }











  void HDivHybridHighOrderFESpace :: GetExternalDofNrs (int elnr, Array<int> & dnums) const
  {
    if (!eliminate_internal) 
      {
	GetDofNrs (elnr, dnums);
	return;
      }

    Array<int> vnums, ednums, fnums;
    int i, j;
    int first,next;


    if(order < 0)
      throw Exception(" HDivHybridHighOrderFESpace :: GetDofNrs() order < 0 ");

    dnums.SetSize(0);
    
    if(ma.GetDimension() == 3)
      {
	ma.GetElFaces (elnr, fnums);

	//Ravier-Thomas
	for (i = 0; i < fnums.Size(); i++)
	  dnums.Append (fnums[i]);

	
	// faces
	for(i=0; i<fnums.Size(); i++)
	  {
	    first = first_face_dof[fnums[i]];
	    next = first_face_dof[fnums[i]+1];
	    for(j=first ; j<next; j++)
	      dnums.Append(j);
	  }

	/*	
	//inner
	first = first_inner_dof[elnr];
	next = first_inner_dof[elnr+1];
	for(j=first; j<next; j++)
	  dnums.Append(j);
	*/
	
      }
    else
      {
	ma.GetElEdges (elnr, ednums);

	
	//Ravier-Thomas
	for (i = 0; i < ednums.Size(); i++)
	  dnums.Append (ednums[i]);
	
	//Edges_normal
	for (i = 0; i < ednums.Size(); i++)
	  {
	    first = first_edge_dof[ednums[i]];
	    next = first_edge_dof[ednums[i]+1];
	    
	    for (j = first; j <next; j++)
	      dnums.Append (j);
	  }
	
	
      }


   // (*testout) << "el " << elnr << " external has dofs " << dnums << endl;



  }






  void HDivHybridHighOrderFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {


    Array<int> vnums, ednums, eorient;
    int fnum,forient;
    int i, j;
    int first,next;

    ma.GetSElEdges (selnr, ednums, eorient);
    ma.GetSElFace (selnr, fnum, forient);

    dnums.SetSize(0);

    if(order <0) throw (" HDivHybridHighOrderFESpace :: GetSDofNrs() order < 0 ");

    if (ma.GetDimension() == 2)
      {
	// RT_0
	for (i = 0; i < ednums.Size(); i++)
          dnums.Append (ednums[i]);
	//edges_normal
	for (i = 0; i < ednums.Size(); i++)
	  {
	    first = first_edge_dof[ednums[i]];
	    next = first_edge_dof[ednums[i]+1];
	    for (j = first; j <next; j++)
	      dnums.Append (j);
	  }
      }
    else
      {
        dnums.Append (fnum);

        if(first_face_dof.Size()>1)
	  {
	    first = first_face_dof[fnum];
	    next = first_face_dof[fnum+1];
	    for(j=first; j<next; j++)
	      dnums.Append(j);
	  }



      }
    // (*testout)<<"SElNr= "<<selnr<<endl;
    //(*testout)<<"SDofNr= "<<dnums<<endl;



  }


  Table<int> * HDivHybridHighOrderFESpace ::
  CreateSmoothingBlocks (const Flags & precflags) const
  {
    int i, j, first;
    int ncnt;

    int ii;
    int SmoothingType=1;

    cout << "SmoothingType " << SmoothingType << endl;

     ArrayMem<int,12> vnums, ednums, fnums, orient;


/*
    for(i=0; i< nel ; i++)
      ma.GetElVertices (i, vnums);


    for(i=0; i< nel; i++)
      ma.GetElEdges (i,vnums,orient);

    int pn1,pn2;
    for(i=0; i< ned; i++)
      ma.GetEdgePNums (i,pn1,pn2);*/

    switch(SmoothingType)
      {
      case 1:  // RT_0 - Faeod - Iaecd - block  (RT_0 - Eaeod - Iaecd - block for 2D)
        if (ma.GetDimension() == 3)
        {
	   ncnt = 1 + nfa;
	   if (!eliminate_internal)
	      ncnt += nel;
	}
	else
	{
           ncnt = 1 + ned;
	   if (!eliminate_internal)
	      ncnt += nel;
	}
	break;
	case 2:  // RT_0 - ElFaces -I  -block ( - block for 2D)
        if (ma.GetDimension() == 3)
        {
	   ncnt = 1 + nel;
	   if (!eliminate_internal)
	      ncnt += nel;
	}
	else
	{
           ncnt = 1 + ned;
	   if (!eliminate_internal)
	      ncnt += nel;
	}
	break;
      }

      
    Array<int> cnt(ncnt);
    cnt = 0;

    cout << " ncnt " << ncnt << endl;

    switch(SmoothingType)
    {
      case 1:
      // RT_0 - Faeod - Iaecd - block  (RT_0 - Eaeod - Iaecd - block for 2D)
      if (ma.GetDimension() == 3)
      {
         int i, j, first;
         cnt[0] = nfa;
         for (i = 0; i < nfa; i++)
	   cnt[i+1] = first_face_dof[i+1]-first_face_dof[i];
           //cnt[i+1] = 1+first_face_dof[i+1]-first_face_dof[i];


	 if (!eliminate_internal)
	   for (i = 0; i < nel; i++)
	      cnt[nfa+1+i] = first_inner_dof[i+1]-first_inner_dof[i];


	 Table<int> & table = *new Table<int> (cnt);
         cnt = 0;
         for (i = 0; i < nfa; i++)
	     table[0][i] = i; //RT_0
         for (i = 0; i < nfa; i++)
	 {
	     //table[i+1][0] = i;
	     //cnt[i+1] = 1;
	     for (j = first_face_dof[i]; j < first_face_dof[i+1]; j++)
	        table[i+1][cnt[i+1]++] = j;
	 }
         if (!eliminate_internal)
	    for (i = 0; i < nel; i++)
	    {
               for (j = first_inner_dof[i]; j < first_inner_dof[i+1]; j++)
                  table[nfa+1+i][cnt[nfa+1+i]++] = j;

	    }
	(*testout) << "smoothingblocks = " << endl << table << endl;
       return &table;

     }
     else
     {
        int i, j, first;
         cnt[0] = ned;
         for (i = 0; i < ned; i++)
	   cnt[i+1] = 1+first_edge_dof[i+1]-first_edge_dof[i];
	 if (!eliminate_internal)
	   for (i = 0; i < nel; i++)
	      cnt[nfa+1+i] = first_inner_dof[i+1]-first_inner_dof[i];

	 Table<int> & table = *new Table<int> (cnt);
         cnt = 0;
         for (i = 0; i < ned; i++)
	     table[0][i] = i; //RT_0
         for (i = 0; i < ned; i++)
	 {
	     table[i+1][0] = i;
	     cnt[i+1] = 1;
	     for (j = first_edge_dof[i]; j < first_edge_dof[i+1]; j++)
	        table[i+1][cnt[i+1]++] = j;
	 }
         if (!eliminate_internal)
	    for (i = 0; i < nel; i++)
	    {
               for (j = first_inner_dof[i]; j < first_inner_dof[i+1]; j++)
	           table[ned+1+i][cnt[ned+1+i]++] = j;
	    }
	    (*testout) << "smoothingblocks = " << endl << table << endl;
	  return &table;
    }

      break;
    case 2:
      // // RT_0 - ElFaces -I  -block ( - block for 2D)
      if (ma.GetDimension() == 3)
      {
         //cout << "cnt = " <<endl;
         int i, j, k, first;
         cnt[0] = nfa;
         for (i = 0; i < nel; i++)
	 {
	    ma.GetElFaces(i,fnums,orient);
	    //cout << "fnums.Size = " <<fnums.Size()<<endl;
	    for (j = 0; j < fnums.Size(); j++)
	    {
               // cout << "firstfacedof = " <<first_face_dof[fnums[j]]<<endl;
	       cnt[i+1] += first_face_dof[fnums[j]+1]-first_face_dof[fnums[j]];
	    }
	 }   
	    //cout << "cnt = " <<endl;
           //cnt[i+1] = 1+first_face_dof[i+1]-first_face_dof[i];

         if (!eliminate_internal)
	   for (i = 0; i < nel; i++)
	      cnt[nel+1+i] = first_inner_dof[i+1]-first_inner_dof[i];

         // cout << "cnt = " << endl << cnt << endl;
	 Table<int> & table = *new Table<int> (cnt);
         cnt = 0;
         for (i = 0; i < nfa; i++)
	     table[0][i] = i; //RT_0
         for (i = 0; i < nel; i++)
	 {
             ma.GetElFaces(i,fnums,orient);
	     for (j = 0; j < fnums.Size(); j++)
	        for (k = first_face_dof[fnums[j]]; k < first_face_dof[fnums[j]+1]; k++)
	         table[i+1][cnt[i+1]++] = k;
         }
         if (!eliminate_internal)
	    for (i = 0; i < nel; i++)
	    {
               for (j = first_inner_dof[i]; j < first_inner_dof[i+1]; j++)
                  table[nel+1+i][cnt[nel+1+i]++] = j;

	    }
	(*testout) << "smoothingblocks = " << endl << table << endl;
       return &table;

     }
     else
     {
       ;
     }

      break;

   }
    return 0;
  }


  /// 
  void HDivHybridHighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  { dnums.SetSize(0); return; }
  /// 
  void HDivHybridHighOrderFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  { 
    dnums.SetSize(0);
    if ( ma.GetDimension() == 2 )
      {
	dnums.Append (ednr);
	
	//Edges_normal
	int first = first_edge_dof[ednr];
	int next = first_edge_dof[ednr+1];
	for (int j = first; j <next; j++)
	  dnums.Append (j);
      }
  }
  /// 
  void HDivHybridHighOrderFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    if(ma.GetDimension() == 3)
      {
	//Ravier-Thomas
	dnums.Append (fanr);
	
	// faces
	int first = first_face_dof[fanr];
	int next = first_face_dof[fanr+1];
	for(int j=first ; j<next; j++)
	  dnums.Append(j);
      }

  }  
  /// 
  void HDivHybridHighOrderFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    int first = first_inner_dof[elnr];
    int next = first_inner_dof[elnr+1];
    for(int j=first; j<next; j++)
      dnums.Append(j);

  }





  // register FESpaces
  namespace
#ifdef MACOS
  hdivhybridhofespace_cpp
#endif
  {

    class Init
    {
    public:
      Init ();
    };

    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("hdivhybridho", HDivHybridHighOrderFESpace::Create);
    }

    Init init;
  }

}


