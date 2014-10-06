/*********************************************************************/
/* File:   l2hofespace.cpp                                         */
/* Author: Start                                                     */
/* Date:   24. Feb. 2003                                             */
/*********************************************************************/

/**
   High Order Finite Element Space for L2
*/

/* ***********************************************
To do: *Internal External Dofs (eliminate internal) 
       *Flag for low-order dofs eliminated ...   
************************* */ 

#include <comp.hpp>
#include <multigrid.hpp>

using namespace ngmg;

namespace ngcomp
{

  L2HighOrderFESpace ::  
  L2HighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name="L2HighOrderFESpace(l2ho)";
    
    // defined flags
    DefineNumFlag("relorder");
    DefineDefineFlag("l2ho");
    DefineDefineFlag("all_dofs_together");

    if (parseflags) CheckFlags(flags);

    var_order = 0; 
    
    if (flags.NumFlagDefined("order"))
      order =  int (flags.GetNumFlag("order",0));
    else 
      {
	if(flags.NumFlagDefined("relorder"))
	  {
	    order=0; 
	    var_order = 1; 
	    rel_order = int (flags.GetNumFlag("relorder",0));
	  }
	else 
	  order = 0; 
      }
    
    if(flags.GetDefineFlag("variableorder") )
      {
	throw Exception ("Flag 'variableorder' for l2ho is obsolete. \n  Either choose uniform order by -order= .. \n -relorder=.. for relative mesh order "); 
      }
    
    integrator = CreateBFI("mass", ma->GetDimension(),
                           make_shared<ConstantCoefficientFunction>(1));

    if (dimension > 1)
      integrator = make_shared<BlockBilinearFormIntegrator> (integrator, dimension);

    

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





    all_dofs_together = flags.GetDefineFlag ("all_dofs_together");

    Flags loflags;
    loflags.SetFlag ("order", 0.0);
    loflags.SetFlag ("dim", dimension);
    if (dgjumps){ *testout << "(L2HOFES:)setting loflag dgjumps " << endl; loflags.SetFlag ("dgjumps");}
    if (iscomplex) loflags.SetFlag ("complex");


    if(all_dofs_together)
      prol = make_shared<L2HoProlongation>(ma,first_element_dof);
    else
      {
        low_order_space = make_shared<ElementFESpace> (ma, loflags);
        prol = make_shared<ElementProlongation> (*static_cast<ElementFESpace*> (low_order_space.get()));
      }
  }

  L2HighOrderFESpace :: ~L2HighOrderFESpace ()
  { ; }

  shared_ptr<FESpace> L2HighOrderFESpace :: 
  Create (shared_ptr<MeshAccess> ma, const Flags & flags)
  {
    int order = int(flags.GetNumFlag ("order", 0));
    if (order == 0)
      return make_shared<ElementFESpace> (ma, flags);
    else
      return make_shared<L2HighOrderFESpace> (ma, flags, true);
  }  

  
  void L2HighOrderFESpace :: Update(LocalHeap & lh)
  {
    FESpace::Update(lh);
    if(low_order_space) low_order_space -> Update(lh);

    nel = ma->GetNE();
    order_inner.SetSize(nel); 
 
    order_inner = INT<3>(order);

    if(var_order) 
      for(int i = 0; i < nel; i++) 
        order_inner[i] = Max(ma->GetElOrders(i)+INT<3>(rel_order), INT<3>(0));

    for(int i = 0; i < nel; i++) 
      if (!DefinedOn (ma->GetElIndex (i))) 
        order_inner[i] = 0;

    if(print) 
      *testout << " order_inner (l2ho) " << order_inner << endl; 

    UpdateDofTables();
    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;

    if(low_order_space) prol->Update();

    UpdateCouplingDofArray();
  } 
  
  void L2HighOrderFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(ndof);
    ctofdof = WIREBASKET_DOF;

    if (!all_dofs_together)
      for (int i=0; i<ma->GetNE(); i++)
	{
          if (!DefinedOn (ma->GetElIndex (i)))
            {
              ctofdof[i] = UNUSED_DOF;
              continue;
            }

	  ctofdof[i] = LOCAL_DOF; //lowest order (constants)
	  int first = first_element_dof[i];
	  int next = first_element_dof[i+1];
	  for (int j = first; j < next; j++)
	    ctofdof[j] = LOCAL_DOF; //higher order
	}
    else
      for (int i=0; i<ma->GetNE(); i++)
	{
	  int first = first_element_dof[i];
	  int next = first_element_dof[i+1];
	  if (next > first)
	    ctofdof[first] = LOCAL_DOF;  //lowest order (constants)
	  for (int j = first+1; j < next; j++)
	    ctofdof[j] = LOCAL_DOF;
	}
  }

  void L2HighOrderFESpace :: UpdateDofTables()
  {
    ndof = all_dofs_together ? 0 : nel;
    first_element_dof.SetSize(nel+1);
    for (int i = 0; i < nel; i++)
      {
	first_element_dof[i] = ndof;
	INT<3> pi = order_inner[i]; 
	switch (ma->GetElType(i))
	  {
	  case ET_SEGM:
	    ndof += pi[0]+1;
	    break;
	  case ET_TRIG:
	    ndof += (pi[0]+1)*(pi[0]+2)/2 ;
	    break;
	  case ET_QUAD:
	    ndof += (pi[0]+1)*(pi[1]+1);
	    break;
	  case ET_TET:
	    ndof += (pi[0]+1)*(pi[0]+2)*(pi[0]+3)/6;
	    break;
	  case ET_PRISM:
	    ndof += (pi[0]+1)*(pi[0]+2)*(pi[2]+1)/2;
	    break;
	  case ET_PYRAMID:
	    ndof += 5 + 8*(pi[0]-1) + 2*(pi[0]-1)*(pi[0]-2) + (pi[0]-1)*(pi[0]-1) 
	      + (pi[0]-1)*(pi[0]-2)*(2*pi[0]-3)/6;
	    break;
	  case ET_HEX:
	    ndof += (pi[0]+1)*(pi[1]+1)*(pi[2]+1);
	    break;
          default:  // for the compiler
            break;
	  }
	if (!all_dofs_together)
	  ndof--; // subtract constant 
      }
    first_element_dof[nel] = ndof;
    
    if(print) 
      *testout << " first_element dof (l2hofe) " << first_element_dof << endl;  

    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;

    prol->Update();
  }


  const FiniteElement & L2HighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    try
      { 
        Ngs_Element ngel = ma->GetElement(elnr);
        ELEMENT_TYPE eltype = ngel.GetType();
        
        if (!DefinedOn (ma->GetElIndex (elnr)))
          {
            switch (eltype)
              {
              case ET_POINT:   return * new (lh) ScalarDummyFE<ET_POINT> (); break;
              case ET_SEGM:    return * new (lh) ScalarDummyFE<ET_SEGM> (); break;
              case ET_TRIG:    return * new (lh) ScalarDummyFE<ET_TRIG> (); break;
              case ET_QUAD:    return * new (lh) ScalarDummyFE<ET_QUAD> (); break;
              case ET_TET:     return * new (lh) ScalarDummyFE<ET_TET> (); break;
              case ET_PYRAMID: return * new (lh) ScalarDummyFE<ET_PYRAMID> (); break;
              case ET_PRISM:   return * new (lh) ScalarDummyFE<ET_PRISM> (); break;
              case ET_HEX:     return * new (lh) ScalarDummyFE<ET_HEX> (); break;
              }
          }

	if (ma->GetElType(elnr) == ET_TRIG) 
	  {
	    ArrayMem<int,3> vnums(3);
	    ma->GetElVertices (elnr, vnums);
	    return *CreateL2HighOrderFE<ET_TRIG> (order, vnums, lh);
	  }

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
            throw Exception ("illegal element in L2HoFeSpace::GetFE");
          }
      } 
    catch (Exception & e)
      {
	e.Append ("in L2HoFESpace::GetElement");
	e.Append ("\n");
	throw; 
      }
  }

  
  template <ELEMENT_TYPE ET>
  const FiniteElement & L2HighOrderFESpace :: T_GetFE (int elnr, LocalHeap & lh) const
  {
    Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM>(elnr);
    L2HighOrderFE<ET> * hofe =  new (lh) L2HighOrderFE<ET> ();
    
    hofe -> SetVertexNumbers (ngel.vertices);
    hofe -> SetOrder (order_inner[elnr]);
    hofe -> ComputeNDof();
    
    return *hofe;
  }




  const FiniteElement & L2HighOrderFESpace :: GetFacetFE (int fnr, LocalHeap & lh) const
  {
    DGFiniteElement<2> * fe2d = NULL;

    ArrayMem<int,4> vnums;
    ma->GetFacetPNums (fnr, vnums);

    switch (vnums.Size())
      {
      case 1: return *new (lh) L2HighOrderFE<ET_POINT> (0); 
      case 2: return *CreateL2HighOrderFE<ET_SEGM> (order, vnums, lh);
      case 3: fe2d = new (lh) L2HighOrderFE<ET_TRIG> (); break;
      case 4: fe2d = new (lh) L2HighOrderFE<ET_QUAD> (); break;
      default:
	{
	  stringstream str;
	  str << "L2HighOrderFESpace " << GetClassName() 
	      << ", undefined facet-eltype" << endl;
	  throw Exception (str.str());
	}
      }
    
    fe2d-> SetVertexNumbers (vnums); 
    fe2d-> SetOrder(order);
    fe2d-> ComputeNDof(); 
    return *fe2d;
  } 







 
  const FiniteElement & L2HighOrderFESpace :: GetSFE (int elnr, LocalHeap & lh) const
  {
    switch (ma->GetSElType(elnr))
      {
      case ET_POINT: return *new (lh) DummyFE<ET_POINT>; 
      case ET_SEGM:  return *new (lh) DummyFE<ET_SEGM>; break;
      case ET_TRIG:  return *new (lh) DummyFE<ET_TRIG>; break;
      case ET_QUAD:  return *new (lh) DummyFE<ET_QUAD>; break;

      default:
	stringstream str;
	str << "FESpace " << GetClassName() 
	    << ", undefined surface eltype " << ma->GetSElType(elnr) 
	    << ", order = " << order << endl;
	throw Exception (str.str());
      }
  }

  int L2HighOrderFESpace :: GetNDof () const
  {
    return ndof;
  }

  int L2HighOrderFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }


  void L2HighOrderFESpace :: GetDofRanges (ElementId ei, Array<IntRange> & dranges) const
  {
    dranges.SetSize(0);

    if (!ei.IsVolume()) return;
    if (!DefinedOn (ma->GetElIndex (ei))) return;
    
    if (!all_dofs_together)
      dranges.Append (ei.Nr());
    dranges.Append (GetElementDofs(ei.Nr()));
  }


  void L2HighOrderFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if (!DefinedOn (ma->GetElIndex (elnr))) return;

    if (!all_dofs_together)
      dnums.Append (elnr); // lowest_order 
    dnums += GetElementDofs(elnr);
  }
  
  void L2HighOrderFESpace :: 
  GetSDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize0 ();
  }
  
  Table<int> * L2HighOrderFESpace :: 
  CreateSmoothingBlocks (const Flags & precflags) const
  {
    int i, j, first;
    Array<int> cnt(nel);
    cnt = 0;
    for (i = 0; i < nel; i++)
      cnt[i] = first_element_dof[i+1]-first_element_dof[i];
	
    Table<int> & table = *new Table<int> (cnt);
    
    for (i = 0; i < nel; i++)
      {
	first = first_element_dof[i];
	for (j = 0; j < cnt[i]; j++)
	  table[i][j] = first+j;
      }
    return &table;
  }

  void  L2HighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  { 
    dnums.SetSize0();
  }
  
  void  L2HighOrderFESpace ::GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  { 
    dnums.SetSize0(); 
  }
  
  void  L2HighOrderFESpace ::GetFaceDofNrs (int fanr, Array<int> & dnums) const
  { 
    dnums.SetSize0(); 
  }
  
  void  L2HighOrderFESpace ::GetInnerDofNrs (int elnr, Array<int> & dnums) const
  { 
    GetDofNrs (elnr, dnums); 
  }










  L2SurfaceHighOrderFESpace ::  
  L2SurfaceHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name="L2SurfaceHighOrderFESpace(l2surf)";
    // defined flags 
    DefineDefineFlag("l2surf");

    if(parseflags) CheckFlags(flags);
    
    
    if(flags.NumFlagDefined("relorder"))
      throw Exception("Variable order not implemented for L2SurfaceHighOrderFESpace"); 
    
    segm = new L2HighOrderFE<ET_SEGM> (order);
    trig = new L2HighOrderFE<ET_TRIG> (order);
    quad = new L2HighOrderFE<ET_QUAD> (order);

    if (ma->GetDimension() == 2)
      {
	// boundary_integrator.reset (new RobinIntegrator<2> (new ConstantCoefficientFunction(1)));
        boundary_integrator = 
          make_shared<RobinIntegrator<2>>(make_shared<ConstantCoefficientFunction>(1));
      }
    else
      {
	boundary_integrator = 
          make_shared<RobinIntegrator<3>> (make_shared<ConstantCoefficientFunction>(1));
      }

    if (dimension > 1)
      boundary_integrator = make_shared<BlockBilinearFormIntegrator> (boundary_integrator, dimension);
  }

  L2SurfaceHighOrderFESpace :: ~L2SurfaceHighOrderFESpace ()
  {
    ;
  }

  shared_ptr<FESpace> L2SurfaceHighOrderFESpace :: 
  Create (shared_ptr<MeshAccess> ma, const Flags & flags)
  {
    return make_shared<L2SurfaceHighOrderFESpace> (ma, flags, true);
  }

  void L2SurfaceHighOrderFESpace :: Update(LocalHeap & lh)
  {
    nel = ma->GetNSE();

    ndof = 0;
    first_element_dof.SetSize(nel+1);
    for (int i = 0; i < nel; i++)
      {
	first_element_dof[i] = ndof;
	switch (ma->GetSElType(i))
	  {
	  case ET_SEGM:
	    ndof += order+1;
	    break;
	  case ET_TRIG:
	    ndof += (order+1)*(order+2)/2;
	    break;
	  case ET_QUAD:
	    ndof += (order+1)*(order+1);
	    break;
	  default:
	    ;
	  }
      }
    first_element_dof[nel] = ndof;
  }

  const FiniteElement & L2SurfaceHighOrderFESpace :: GetSFE (int elnr, LocalHeap & lh) const
  {
    if (ma->GetDimension() == 2)
      {
	DGFiniteElement<1> * fe1d = 0;
	
	Ngs_Element ngel = ma->GetElement<1> (elnr);

	switch (ngel.GetType())
	  {
	  case ET_SEGM: fe1d = new (lh) L2HighOrderFE<ET_SEGM> (); break;
	  default:
	    ;
	  }

	fe1d -> SetVertexNumbers (ngel.vertices);
	fe1d -> SetOrder (INT<1> (order));
	fe1d -> ComputeNDof(); 
	return *fe1d;
      }
    else
      {
	DGFiniteElement<2> * fe2d = 0;
	
	Ngs_Element ngel = ma->GetElement<2> (elnr);
	
	switch (ngel.GetType())
	  {
	  case ET_TRIG: fe2d = new (lh) L2HighOrderFE<ET_TRIG> (); break;
	  case ET_QUAD: fe2d = new (lh) L2HighOrderFE<ET_QUAD> (); break;
	  default:
	    ;
	  }
	
	fe2d -> SetVertexNumbers (ngel.vertices);
	fe2d -> SetOrder (INT<2> (order));
	fe2d -> ComputeNDof(); 
	return *fe2d;
      }
	/*

    FiniteElement * fe = 0;
    
    switch (ma->GetSElType(elnr))
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
    
    ArrayMem<int,12> vnums; // calls GetElPNums -> max 12 for PRISM12
    ma->GetSElVertices(elnr, vnums);

    if (!fe)
      {
	stringstream str;
	str << "L2SurfaceHighOrderFESpace " << GetClassName() 
	    << ", undefined eltype " 
	    << ElementTopology::GetElementName(ma->GetSElType(elnr))
	    << ", order = " << order << endl;
	throw Exception (str.str());
      }

    dynamic_cast<L2HighOrderFiniteElement<2>*> (fe) -> SetVertexNumbers (vnums);

    return *fe;
	*/
  }
 
  const FiniteElement & L2SurfaceHighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    throw Exception ("Volume elements not available for L2SurfaceHighOrderFESpace");
  }
 
  int L2SurfaceHighOrderFESpace :: GetNDof () const
  {
    return ndof;
  }

  void L2SurfaceHighOrderFESpace :: GetSDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    int first = first_element_dof[elnr];
    int neldofs = first_element_dof[elnr+1] - first;
    for (int j = 0; j < neldofs; j++)
      dnums.Append (first+j);

    if (!DefinedOn (ma->GetSElIndex (elnr)))
      dnums = -1;
  }
  
  void L2SurfaceHighOrderFESpace :: 
  GetDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize (0);
  }
  
  Table<int> * L2SurfaceHighOrderFESpace :: 
  CreateSmoothingBlocks ( int type) const
  {
    int i, j, first;
    Array<int> cnt(nel);
    cnt = 0;
    for (i = 0; i < nel; i++)
      cnt[i] = first_element_dof[i+1]-first_element_dof[i];
	
    Table<int> & table = *new Table<int> (cnt);
    
    for (i = 0; i < nel; i++)
      {
	first = first_element_dof[i];
	for (j = 0; j < cnt[i]; j++)
	  table[i][j] = first+j;
      }
    return &table;
  }


  void  L2SurfaceHighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  { dnums.SetSize(0); return; }
  
  void  L2SurfaceHighOrderFESpace ::GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  { dnums.SetSize(0); return; }
  
  void  L2SurfaceHighOrderFESpace ::GetFaceDofNrs (int fanr, Array<int> & dnums) const
  { GetDofNrs ( fanr, dnums ); return; }
  
  void  L2SurfaceHighOrderFESpace ::GetInnerDofNrs (int elnr, Array<int> & dnums) const
  { GetDofNrs ( elnr, dnums ); return; }
    




  // register FESpaces
  namespace l2hofespace_cpp
  {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("l2", L2HighOrderFESpace::Create);
      GetFESpaceClasses().AddFESpace ("l2ho", L2HighOrderFESpace::CreateHO);
      GetFESpaceClasses().AddFESpace ("l2surf", L2SurfaceHighOrderFESpace::Create);
    }
    
    Init init;
  }
}
