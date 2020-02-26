/*********************************************************************/
/* File:   l2hofespace.cpp                                         */
/* Author: Start                                                     */
/* Date:   24. Feb. 2003                                             */
/*********************************************************************/

/**
   High Order Finite Element Space for L2
*/

#include <comp.hpp>
#include <multigrid.hpp>
#include "../fem/l2hofetp.hpp"

using namespace ngmg;

namespace ngcomp
{

  class BlockDifferentialOperatorId : public BlockDifferentialOperator
  {
  public:
    using BlockDifferentialOperator::BlockDifferentialOperator;

    virtual void Apply (const FiniteElement & fel,
                        const SIMD_BaseMappedIntegrationRule & mir,
                        BareSliceVector<double> x, 
                        BareSliceMatrix<SIMD<double>> flux) const
    {
      if (comp == -1)
        static_cast<const BaseScalarFiniteElement&> (fel).
          Evaluate(mir.IR(), SliceMatrix<double> (fel.GetNDof(), dim, dim, &x(0)), flux);
      else
        diffop->Apply(fel, mir, x.Slice(comp, dim), flux.RowSlice(comp,dim));
    }
    virtual void
    AddTrans (const FiniteElement & fel,
              const SIMD_BaseMappedIntegrationRule & mir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const
    {
    if (comp == -1)
      static_cast<const BaseScalarFiniteElement&> (fel).      
        AddTrans(mir.IR(), flux, SliceMatrix<double> (fel.GetNDof(), dim, dim, &x(0)));
    else
      diffop->AddTrans(fel, mir, flux.RowSlice(comp,dim), x.Slice(comp,dim));
    }
    
  };


  
  L2HighOrderFESpace ::  
  L2HighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name="L2HighOrderFESpace(l2ho)";
    type = "l2ho";
   
    // defined flags
    DefineNumFlag("relorder");
    DefineDefineFlag("l2ho");
    DefineDefineFlag("all_dofs_together");
    DefineDefineFlag("hide_all_dofs");

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

    /*
    integrator[VOL] = CreateBFI("mass", ma->GetDimension(),
                           make_shared<ConstantCoefficientFunction>(1));
    if (dimension > 1)
      integrator[VOL] = make_shared<BlockBilinearFormIntegrator> (integrator[VOL], dimension);
    */

    SetDefinedOn(BND, BitArray(ma->GetNRegions(BND)).Clear());

    switch (ma->GetDimension())
      {
      case 1:
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<1>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<1>>>();
          // evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<1>>>();
          break;
        }
      case 2:
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
          // evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
          break;
        }
      case 3:
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
          // evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
          break;
        }
      }
    if (dimension > 1) 
      {
        evaluator[VOL] = make_shared<BlockDifferentialOperatorId> (evaluator[VOL], dimension);
        // evaluator[VOL] = make_shared<BlockDifferentialOperator> (evaluator[VOL], dimension);
	flux_evaluator[VOL] = make_shared<BlockDifferentialOperator> (flux_evaluator[VOL], dimension);
	// evaluator[BND] = make_shared<BlockDifferentialOperator> (evaluator[BND], dimension);
        /*
	boundary_flux_evaluator = 
	  make_shared<BlockDifferentialOperator> (boundary_flux_evaluator, dimension);
        */
      }

    additional_evaluators.Set ( "dual", evaluator[VOL] );


    tensorproduct = flags.GetDefineFlag ("tp");
    all_dofs_together = flags.GetDefineFlag ("all_dofs_together");
    hide_all_dofs = flags.GetDefineFlag ("hide_all_dofs");

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

    lowest_order_ct =
      flags.GetDefineFlagX ("lowest_order_wb").IsTrue() ? WIREBASKET_DOF : LOCAL_DOF;

    switch (ma->GetDimension())
      {
      case 1:
        additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpHesse<1>>> ());
        break;
      case 2:
        additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpHesse<2>>> ()); 
        break;
      case 3:
        additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpHesse<3>>> ());
	break;
      default:
        ;
      }

    needs_transform_vec = false;    
  }

  L2HighOrderFESpace :: ~L2HighOrderFESpace ()
  { ; }

  DocInfo L2HighOrderFESpace :: GetDocu ()
  {
    DocInfo docu = FESpace::GetDocu(); 
    docu.short_docu = "An L2-conforming finite element space.";
    docu.long_docu =
      R"raw_string(The L2 finite element space consists of element-wise polynomials,
which are discontinuous from element to element. It uses an
L2-orthogonal hierarchical basis which leads to orthogonal
mass-matrices on non-curved elements.

Boundary values are not meaningful for an L2 function space.

The L2 space supports element-wise variable order, which can be set
for ELEMENT-nodes.

Per default, all dofs are local dofs and are condensed if static
condensation is performed. The lowest order can be kept in the
WIRE_BASKET via the flag 'lowest_order_wb=True'.

All dofs can be hidden. Then the basis functions don't show up in the
global system.
)raw_string";
    
    docu.Arg("all_dofs_together") = "bool = False\n"
      "  Change ordering of dofs. If this flag ist set,\n"
      "  all dofs of an element are ordered successively.\n"
      "  Otherwise, the lowest order dofs (the constants)\n"
      "  of all elements are ordered first.";

    docu.Arg("lowest_order_wb") = "bool = False\n"
      "  Keep lowest order dof in WIRE_BASKET";

    docu.Arg("hide_all_dofs") = "bool = False\n"
      "  Set all used dofs to HIDDEN_DOFs";
    return docu;
  }

  
  shared_ptr<FESpace> L2HighOrderFESpace :: 
  Create (shared_ptr<MeshAccess> ma, const Flags & flags)
  {
    int order = int(flags.GetNumFlag ("order", 0));
    if (order == 0)
      return make_shared<ElementFESpace> (ma, flags);
    else
      return make_shared<L2HighOrderFESpace> (ma, flags, true);
  }  

  
  void L2HighOrderFESpace :: Update()
  {
    FESpace::Update();
    if(low_order_space) low_order_space -> Update();

    nel = ma->GetNE();

    bool first_update = GetTimeStamp() < ma->GetTimeStamp();
    if (first_update) timestamp = NGS_Object::GetNextTimeStamp();

    if (first_update)
      {
	order_inner.SetSize(nel); 
	
	order_inner = INT<3>(order);
	
	if(var_order) 
	  for(int i = 0; i < nel; i++) 
	    order_inner[i] = ma->GetElOrders(i)+INT<3>(rel_order);
    
	for(int i = 0; i < nel; i++) 
	  {
	    ElementId ei(VOL,i);
	    order_inner[i] = order_inner[i] + INT<3> (et_bonus_order[ma->GetElType(ei)]);
	    order_inner[i] = Max(order_inner[i], INT<3>(0));
	    if (!DefinedOn (ei))
	      order_inner[i] = 0;
	  }
	if(print) 
	  *testout << " order_inner (l2ho) " << order_inner << endl;
      }

    UpdateDofTables();
    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;

    if(low_order_space) prol->Update(*this);

    UpdateCouplingDofArray();
  } 
  
  void L2HighOrderFESpace :: UpdateCouplingDofArray()
  {
    auto ct_local = hide_all_dofs ? HIDDEN_DOF : LOCAL_DOF;
    auto ct_lowest_order = hide_all_dofs ? HIDDEN_DOF : lowest_order_ct;
    ctofdof.SetSize(ndof);
    for (auto i : Range(ma->GetNE()))
      {
        bool definedon = DefinedOn(ElementId(VOL,i));
        auto r = GetElementDofs(i);
        ctofdof[r] = definedon ? ct_local : UNUSED_DOF;
        
        if (!all_dofs_together)
	  ctofdof[i] = definedon ? ct_lowest_order : UNUSED_DOF;
        else
          if (r.Size() != 0)
            ctofdof[r.First()] = definedon ? ct_lowest_order : UNUSED_DOF;            
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
	switch (ma->GetElType(ElementId(VOL,i)))
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

    prol->Update(*this);
  }


  

  FiniteElement & L2HighOrderFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    ELEMENT_TYPE eltype = ngel.GetType();

    if (ei.IsVolume())
      {
        int elnr = ei.Nr();
        
        if (!DefinedOn (ngel))
          {
            /*
            switch (eltype)
              {
              case ET_POINT:   return * new (alloc) ScalarDummyFE<ET_POINT> (); break;
              case ET_SEGM:    return * new (alloc) ScalarDummyFE<ET_SEGM> (); break;
              case ET_TRIG:    return * new (alloc) ScalarDummyFE<ET_TRIG> (); break;
              case ET_QUAD:    return * new (alloc) ScalarDummyFE<ET_QUAD> (); break;
              case ET_TET:     return * new (alloc) ScalarDummyFE<ET_TET> (); break;
              case ET_PYRAMID: return * new (alloc) ScalarDummyFE<ET_PYRAMID> (); break;
              case ET_PRISM:   return * new (alloc) ScalarDummyFE<ET_PRISM> (); break;
              case ET_HEX:     return * new (alloc) ScalarDummyFE<ET_HEX> (); break;
              }
            */
            return SwitchET(eltype,
                            [&alloc] (auto et) -> FiniteElement&
                            { return *new(alloc) ScalarDummyFE<et.ElementType()>(); });
          }

	if (eltype == ET_TRIG && order_policy == CONSTANT_ORDER) 
          return *CreateL2HighOrderFE<ET_TRIG> (order, INT<3>(ngel.Vertices()), alloc);

        if (tensorproduct)
          {
            if (eltype == ET_QUAD)
              return * new (alloc) L2HighOrderFETP<ET_QUAD> (order, ngel.Vertices(), alloc);
            if (eltype == ET_TET)
              return * new (alloc) L2HighOrderFETP<ET_TET> (order, ngel.Vertices(), alloc);
            if (eltype == ET_HEX)
              return * new (alloc) L2HighOrderFETP<ET_HEX> (order, ngel.Vertices(), alloc);
          }
        if (eltype == ET_TET && order_policy == CONSTANT_ORDER) 
          return *CreateL2HighOrderFE<ET_TET> (order, INT<4>(ngel.Vertices()), alloc);

        /*
        return SwitchET(eltype,
                        [this,elnr,&alloc] (auto et) -> FiniteElement&
                        {
                        // return T_GetFE<et.ElementType()>(elnr, alloc);
                          // Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL>(elnr);
                          auto * hofe =  new (alloc) L2HighOrderFE<et.ElementType()> ();
                          hofe -> SetVertexNumbers (ngel.vertices);
                          hofe -> L2HighOrderFE<et.ElementType()>::SetOrder (order_inner[ngel.Nr()]);
                          hofe -> L2HighOrderFE<et.ElementType()>::ComputeNDof();
                          return *hofe;
                        });
        */

        switch (eltype)
          {
          case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, alloc);

          case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, alloc);
          case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, alloc);
            
          case ET_TET:     return T_GetFE<ET_TET> (elnr, alloc);
          case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, alloc);
          case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, alloc);
          case ET_HEX:     return T_GetFE<ET_HEX> (elnr, alloc);
            
          default:
            throw Exception ("illegal element in L2HoFeSpace::GetFE");
          }
      }
    else
      {
        try
          {
            return SwitchET<ET_POINT,ET_SEGM,ET_TRIG,ET_QUAD>
              (eltype,
               [&alloc] (auto et) -> FiniteElement&
               { return * new (alloc) DummyFE<et.ElementType()>; });
          }
        catch (Exception e)
          {
            throw Exception("illegal element type in L2::GetSurfaceFE");
          }
        /*
        switch (eltype)
          {
          case ET_POINT: return *new (alloc) DummyFE<ET_POINT>; 
          case ET_SEGM:  return *new (alloc) DummyFE<ET_SEGM>;
          case ET_TRIG:  return *new (alloc) DummyFE<ET_TRIG>;
          case ET_QUAD:  return *new (alloc) DummyFE<ET_QUAD>;
            
          default:
            stringstream str;
            str << "FESpace " << GetClassName() 
                << ", undefined surface eltype " << ma->GetElType(ei) 
                << ", order = " << order << endl;
            throw Exception (str.str());
          }
        */
      }
  }

  // const FiniteElement & L2HighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  // {
  //   try
  //     { 
  //       Ngs_Element ngel = ma->GetElement(elnr);
  //       ELEMENT_TYPE eltype = ngel.GetType();
        
  //       if (!DefinedOn (ngel.GetIndex()))
  //         {
  //           switch (eltype)
  //             {
  //             case ET_POINT:   return * new (lh) ScalarDummyFE<ET_POINT> (); break;
  //             case ET_SEGM:    return * new (lh) ScalarDummyFE<ET_SEGM> (); break;
  //             case ET_TRIG:    return * new (lh) ScalarDummyFE<ET_TRIG> (); break;
  //             case ET_QUAD:    return * new (lh) ScalarDummyFE<ET_QUAD> (); break;
  //             case ET_TET:     return * new (lh) ScalarDummyFE<ET_TET> (); break;
  //             case ET_PYRAMID: return * new (lh) ScalarDummyFE<ET_PYRAMID> (); break;
  //             case ET_PRISM:   return * new (lh) ScalarDummyFE<ET_PRISM> (); break;
  //             case ET_HEX:     return * new (lh) ScalarDummyFE<ET_HEX> (); break;
  //             }
  //         }

  //       if (ngel.GetType() == ET_TRIG) 
  //         {
  //           int ia[3];
  //           FlatArray<int> vnums(3, &ia[0]);
  //           vnums = ngel.Vertices();
  //           return *CreateL2HighOrderFE<ET_TRIG> (order, vnums, lh);
  //         }

  //       if (eltype == ET_TET)         
  //         return *CreateL2HighOrderFE<ET_TET> (order, INT<4>(ngel.Vertices()), lh);

  //       switch (eltype)
  //         {
  //         case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, lh);

  //         case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, lh);
  //         case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, lh);
            
  //         case ET_TET:     return T_GetFE<ET_TET> (elnr, lh);
  //         case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, lh);
  //         case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, lh);
  //         case ET_HEX:     return T_GetFE<ET_HEX> (elnr, lh);
            
  //         default:
  //           throw Exception ("illegal element in L2HoFeSpace::GetFE");
  //         }
  //     } 
  //   catch (Exception & e)
  //     {
  //       e.Append ("in L2HoFESpace::GetElement");
  //       e.Append ("\n");
  //       throw; 
  //     }
  // }

  
  template <ELEMENT_TYPE ET>
  FiniteElement & L2HighOrderFESpace :: T_GetFE (int elnr, Allocator & lh) const
  {
    Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL>(elnr);
    L2HighOrderFE<ET> * hofe =  new (lh) L2HighOrderFE<ET> ();

    hofe -> SetVertexNumbers (ngel.vertices);
    hofe -> L2HighOrderFE<ET>::SetOrder (order_inner[elnr]);
    hofe -> L2HighOrderFE<ET>::ComputeNDof();

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







 
  // const FiniteElement & L2HighOrderFESpace :: GetSFE (int elnr, LocalHeap & lh) const
  // {
  //   switch (ma->GetSElType(elnr))
  //     {
  //     case ET_POINT: return *new (lh) DummyFE<ET_POINT>; 
  //     case ET_SEGM:  return *new (lh) DummyFE<ET_SEGM>; break;
  //     case ET_TRIG:  return *new (lh) DummyFE<ET_TRIG>; break;
  //     case ET_QUAD:  return *new (lh) DummyFE<ET_QUAD>; break;

  //     default:
  //       stringstream str;
  //       str << "FESpace " << GetClassName() 
  //           << ", undefined surface eltype " << ma->GetSElType(elnr) 
  //           << ", order = " << order << endl;
  //       throw Exception (str.str());
  //     }
  // }

  size_t L2HighOrderFESpace :: GetNDof () const throw()
  {
    return ndof;
  }

  size_t L2HighOrderFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }


  void L2HighOrderFESpace :: GetDofRanges (ElementId ei, Array<IntRange> & dranges) const
  {
    dranges.SetSize(0);

    if (!ei.IsVolume()) return;
    if (!DefinedOn (ei)) return;
    
    if (!all_dofs_together)
      dranges.Append (IntRange (ei.Nr(), ei.Nr()+1));
    dranges.Append (GetElementDofs(ei.Nr()));
  }


  void L2HighOrderFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if (!DefinedOn (ei) || ei.VB() != VOL) return;

    auto eldofs = GetElementDofs(ei.Nr());
    size_t size = eldofs.Size();
    size_t base = all_dofs_together ? 0 : 1;
    size += base;
    dnums.SetSize(size);
    if (!all_dofs_together) dnums[0] = ei.Nr();
    dnums.Range(base, size) = eldofs;
    /*
    if (!all_dofs_together)
      dnums.Append (elnr); // lowest_order 
    dnums += GetElementDofs(elnr);
    */
  }

  void L2HighOrderFESpace :: SetOrder (NodeId ni, int order)
  {
    if (order_policy == CONSTANT_ORDER || order_policy == NODE_TYPE_ORDER)
      throw Exception("In L2HighOrderFESpace::SetOrder. Order policy is constant or node-type!");
    else if (order_policy == OLDSTYLE_ORDER)
      order_policy = VARIABLE_ORDER;
      
    if (order < 0)
      order = 0;

    if (CoDimension(ni.GetType(), ma->GetDimension()) == 0)
      {
	if (ma->GetDimension() == 2 && ni.GetType() == NT_FACE)
	  {
	    Array<int> elnr;
	    ma->GetFacetSurfaceElements(ni.GetNr(),elnr);
	    if (elnr[0] < order_inner.Size())
	      order_inner[elnr[0]] = order;
	  }
        else if (ni.GetNr() < order_inner.Size())
          order_inner[ni.GetNr()] = order;
      }
    else
      throw Exception ("L2HighOrderFESpace::SetOrder requires NodeType of codimension 0!");
  }
  
  int L2HighOrderFESpace :: GetOrder (NodeId ni) const
  {
    if (CoDimension(ni.GetType(), ma->GetDimension()) == 0)
      {
	if (ma->GetDimension() == 2 && ni.GetType() == NT_FACE)
	  {
	    Array<int> elnr;
	    ma->GetFacetSurfaceElements(ni.GetNr(),elnr);
	    if (elnr[0] < order_inner.Size())
	      return order_inner[elnr[0]][0];
	  }
        else if (ni.GetNr() < order_inner.Size())
          return order_inner[ni.GetNr()][0];
      }
    
    return 0;
  }
  
  shared_ptr<Table<int>> L2HighOrderFESpace :: 
  CreateSmoothingBlocks (const Flags & precflags) const
  {
    int i, j, first;
    Array<int> cnt(nel);
    cnt = 0;
    for (i = 0; i < nel; i++)
      cnt[i] = first_element_dof[i+1]-first_element_dof[i];
	
    Table<int> table(cnt);
    
    for (i = 0; i < nel; i++)
      {
	first = first_element_dof[i];
	for (j = 0; j < cnt[i]; j++)
	  table[i][j] = first+j;
      }
    return make_shared<Table<int>> (table);
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


#include "applymassl2.hpp"

  shared_ptr<BaseMatrix> L2HighOrderFESpace ::
  GetMassOperator (shared_ptr<CoefficientFunction> rho,
                   shared_ptr<Region> defon,
                   LocalHeap & lh) const
  {
    if (rho->ElementwiseConstant() && all_dofs_together && order_policy == CONSTANT_ORDER)
      {
        return make_shared<ApplyMassL2Const>
          (dynamic_pointer_cast<FESpace>(const_cast<L2HighOrderFESpace*>(this)->shared_from_this()),
           rho, defon, lh);    
      }
    
    return FESpace::GetMassOperator(rho, defon, lh);
  }

  /*
  shared_ptr<BaseMatrix> L2HighOrderFESpace ::
  GetMassOperator (shared_ptr<CoefficientFunction> rho,
                   shared_ptr<Region> defon,
                   LocalHeap & lh) const
  {
    return FESpace::GetMassOperator(rho, defon, lh);
  }
  */

  
  
  void L2HighOrderFESpace :: SolveM (CoefficientFunction * rho, BaseVector & vec, Region * def,
                                     LocalHeap & lh) const
  {
    static Timer t("SolveM"); RegionTimer reg(t);
    if (rho && rho->Dimension() != 1)
      throw Exception("L2HighOrderFESpace::SolveM needs a scalar density");
    IterateElements (*this, VOL, lh,
                     [&rho, &vec, def, this] (FESpace::Element el, LocalHeap & lh)
                     {
                       auto & fel = static_cast<const BaseScalarFiniteElement&>(el.GetFE());
                       const ElementTransformation & trafo = el.GetTrafo();
                       
                       Array<int> dnums(fel.GetNDof(), lh);
                       GetDofNrs (el.Nr(), dnums);

                       FlatVector<double> elx(fel.GetNDof()*dimension, lh);

                       if (def && !def->Mask()[ma->GetElIndex(el)])
                         {
                           elx = 0.0;
                           vec.SetIndirect (dnums, elx);
                           return;
                         }
                       
                       vec.GetIndirect(dnums, elx);
		       auto melx = elx.AsMatrix(fel.GetNDof(),dimension);

                       FlatVector<double> diag_mass(fel.GetNDof(), lh);
                       fel.GetDiagMassMatrix (diag_mass);

                       bool curved = trafo.IsCurvedElement();
                       if (rho && !rho->ElementwiseConstant()) curved = true;

                       if (!curved)
                         {
                           IntegrationRule ir(fel.ElementType(), 0);
                           BaseMappedIntegrationRule & mir = trafo(ir, lh);
                           double jac = mir[0].GetMeasure();
                           if (rho) jac *= rho->Evaluate(mir[0]);
                           diag_mass *= jac;
                           for (int i = 0; i < melx.Height(); i++)
                             melx.Row(i) /= diag_mass(i);
                         }
                       else
                         {
                           SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
                           auto & mir = trafo(ir, lh);
                           FlatVector<SIMD<double>> pntvals(ir.Size(), lh);
                           FlatMatrix<SIMD<double>> rhovals(1, ir.Size(), lh);
                           if (rho) rho->Evaluate (mir, rhovals);
                           
                           for (int i = 0; i < melx.Height(); i++)
                             melx.Row(i) /= diag_mass(i);
                           for (int comp = 0; comp < dimension; comp++)
                             {
                               fel.Evaluate (ir, melx.Col(comp), pntvals);
                               if (rho)
                                 for (size_t i = 0; i < ir.Size(); i++)
                                   pntvals(i) *= ir[i].Weight() / (mir[i].GetMeasure() * rhovals(0,i));
                               else
                                 for (size_t i = 0; i < ir.Size(); i++)
                                   pntvals(i) *= ir[i].Weight() / mir[i].GetMeasure();
                               
                               melx.Col(comp) = 0.0;
                               fel.AddTrans (ir, pntvals, melx.Col(comp));
                             }
                           for (int i = 0; i < melx.Height(); i++)
                             melx.Row(i) /= diag_mass(i);
                         }
                       vec.SetIndirect(dnums, elx);
                     });
  }
  

  void L2HighOrderFESpace :: ApplyM (CoefficientFunction * rho, BaseVector & vec, Region * def,
                                     LocalHeap & lh) const
  {
    static Timer t("ApplyM"); RegionTimer reg(t);
    static Timer tall("ApplyM - all");
    static Timer tel("ApplyM - el");    
    static Timer ttrafo("ApplyM - trafo");    
    static Timer tdofs("ApplyM - getdofs");
    static Timer tgetx("ApplyM - getx");
    static Timer tsety("ApplyM - sety");
    static Timer tcalc("ApplyM - calc");
    static Timer tcalc1("ApplyM - calc1");
    static Timer tcalc2("ApplyM - calc2");
    static Timer tsetup("ApplyM - setup");
    if (rho && rho->Dimension() != 1)
      throw Exception("L2HighOrderFESpace::ApplyM needs a scalar density");

    auto fv = vec.FV<double>();
    
    IterateElements (*this, VOL, lh,
                     [&rho, &vec, fv, def, this] (FESpace::Element el, LocalHeap & lh)
                     {
                       auto tid = TaskManager::GetThreadId();
                       NgProfiler::StartThreadTimer(tall, tid);                       
                       NgProfiler::StartThreadTimer(tel, tid);
                       
                       auto & fel = static_cast<const BaseScalarFiniteElement&>(el.GetFE());
                       NgProfiler::StopThreadTimer(tel, tid);
                       NgProfiler::AddThreadFlops(tel, tid, 1);                       
                       NgProfiler::StartThreadTimer(ttrafo, tid);                       
                       const ElementTransformation & trafo = el.GetTrafo();
                       NgProfiler::StopThreadTimer(ttrafo, tid);
                       NgProfiler::StartThreadTimer(tdofs, tid);
                       
                       Array<int> dnums(fel.GetNDof(), lh);
                       auto dofrange = GetElementDofs(el.Nr());
                       FlatVector<double> elx(fel.GetNDof()*dimension, lh);

                       bool lindofs = all_dofs_together && dimension==1;

                       if (def && !def->Mask()[ma->GetElIndex(el)])
                         {
                           if (lindofs)
                             fv.Range(dofrange) = 0.0;
                           else
                             {
                               elx = 0.0;
                               GetDofNrs (el, dnums);                               
                               vec.SetIndirect(dnums, elx);
                             }
                           return;
                         }
                       
                       if (!lindofs)
                         {
                           GetDofNrs (el, dnums);
                           vec.GetIndirect(dnums, elx);
                         }
                       else
                         elx = fv.Range(dofrange);
                       
                       NgProfiler::StopThreadTimer(tdofs, tid);
                       NgProfiler::StartThreadTimer(tgetx, tid);
                       
		       auto melx = elx.AsMatrix(fel.GetNDof(),dimension);
                       
                       NgProfiler::StopThreadTimer(tgetx, tid);
                       
                       NgProfiler::StartThreadTimer(tsetup, tid);                       
                       FlatVector<double> diag_mass(fel.GetNDof(), lh);
                       fel.GetDiagMassMatrix (diag_mass);

                       bool curved = trafo.IsCurvedElement();
                       if (rho && !rho->ElementwiseConstant()) curved = true;
                       NgProfiler::StopThreadTimer(tsetup, tid);

                       NgProfiler::StartThreadTimer(tcalc, tid);                       
                       if (!curved)
                         {
                           NgProfiler::StartThreadTimer(tcalc1, tid);                                                  
                           IntegrationRule ir(fel.ElementType(), 0);
                           BaseMappedIntegrationRule & mir = trafo(ir, lh);
                           double jac = mir[0].GetMeasure();
                           if (rho) jac *= rho->Evaluate(mir[0]);
                           
                           
                           NgProfiler::StopThreadTimer(tcalc1, tid);
                           NgProfiler::StartThreadTimer(tcalc2, tid);

                           if (dimension == 1)
                             for (size_t i = 0; i < elx.Size(); i++)
                               elx(i) *= jac*diag_mass(i);
                           else
                             for (size_t i = 0; i < melx.Height(); i++)
                               melx.Row(i) *= jac*diag_mass(i);
                           NgProfiler::StopThreadTimer(tcalc2, tid);                           
                         }
                       else
                         {
                           throw Exception ("L2HighOrderFESpace::ApplyM for curved not available");

                           SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
                           auto & mir = trafo(ir, lh);
                           FlatVector<SIMD<double>> pntvals(ir.Size(), lh);
                           FlatMatrix<SIMD<double>> rhovals(1, ir.Size(), lh);
                           if (rho) rho->Evaluate (mir, rhovals);
                           
                           for (int i = 0; i < melx.Height(); i++)
                             melx.Row(i) /= diag_mass(i);
                           for (int comp = 0; comp < dimension; comp++)
                             {
                               fel.Evaluate (ir, melx.Col(comp), pntvals);
                               if (rho)
                                 for (size_t i = 0; i < ir.Size(); i++)
                                   pntvals(i) *= ir[i].Weight() / (mir[i].GetMeasure() * rhovals(0,i));
                               else
                                 for (size_t i = 0; i < ir.Size(); i++)
                                   pntvals(i) *= ir[i].Weight() / mir[i].GetMeasure();
                               
                               melx.Col(comp) = 0.0;
                               fel.AddTrans (ir, pntvals, melx.Col(comp));
                             }
                           for (int i = 0; i < melx.Height(); i++)
                             melx.Row(i) /= diag_mass(i);
                         }
                       NgProfiler::StopThreadTimer(tcalc, tid);
                       
                       NgProfiler::StartThreadTimer(tsety, tid);
                       
                       if (!lindofs)
                         vec.SetIndirect(dnums, elx);
                       else
                         fv.Range(dofrange) = elx;
                       
                       NgProfiler::StopThreadTimer(tsety, tid);
                       NgProfiler::StopThreadTimer(tall, tid);                                              
                     });
  }


  Matrix<> GetTraceMatrix (const FiniteElement & fel)
  {
    auto * dgfel2 = dynamic_cast<const DGFiniteElement<2>*> (&fel);
    if (dgfel2)
      {
        int order = fel.Order();
        Matrix<> trace(3*(order+1), fel.GetNDof());
        for (int j = 0; j < 3; j++)
          dgfel2->CalcTraceMatrix(j, trace.Rows(j*(order+1), (j+1)*(order+1)));
        return trace;
      }
    auto * dgfel3 = dynamic_cast<const DGFiniteElement<3>*> (&fel);
    if (dgfel3)
      {
        int order = fel.Order();
        int nd2d = (order+1)*(order+2)/2;
        Matrix<> trace(4*nd2d, fel.GetNDof());
        for (int j = 0; j < 4; j++)
          dgfel3->CalcTraceMatrix(j, trace.Rows(j*nd2d, (j+1)*nd2d));
        return trace;
      }
    throw Exception("no trace");
  }


  shared_ptr<BaseMatrix> L2HighOrderFESpace ::
  GetTraceOperator (shared_ptr<FESpace> tracespace) const
  {
    LocalHeap lh(1000000);
    Array<short> classnr(ma->GetNE());
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         classnr[el.Nr()] = 
           SwitchET<ET_TRIG,ET_TET>
           (el.GetType(),
            [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
       });
    
    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);
    Table<size_t> table = creator.MoveTable();

    shared_ptr<BaseMatrix> sum;
  
    // size_t ne = ma->GetNE();
  
    for (auto elclass_inds : table)
      {
        if (elclass_inds.Size() == 0) continue;
        
        ElementId ei(VOL,elclass_inds[0]);
        auto & felx = GetFE (ei, lh);
        //auto & trafo = GetMeshAccess()->GetTrafo(ei, lh);
        
        Matrix<> trace_op_x = GetTraceMatrix (felx);


        Table<DofId> xdofs(elclass_inds.Size(), felx.GetNDof()),
          ydofs(elclass_inds.Size(), trace_op_x.Height());

        Array<DofId> dnumsx, dnumsy;
        for (auto i : Range(elclass_inds))
          {
            ElementId ei(VOL, elclass_inds[i]);
            GetDofNrs(ei, dnumsx);
            tracespace->GetDofNrs(ei, dnumsy);
            xdofs[i] = dnumsx;
            ydofs[i] = dnumsy;
          }

        auto mat = make_shared<ConstantElementByElementMatrix>
          (tracespace->GetNDof(), this->GetNDof(),
           trace_op_x, std::move(ydofs), std::move(xdofs));

        if (sum)
          sum = make_shared<SumMatrix>(sum, mat);
        else
          sum = mat;
      }
    return sum;
  }

  
  void L2HighOrderFESpace ::
  GetTrace (const FESpace & tracespace, const BaseVector & in, BaseVector & out, bool avg,
            LocalHeap & lh) const 
  {
    static Timer t("GetTrace"); RegionTimer reg(t);

    out = 0.0;
    Array<short> classnr(ma->GetNE());
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         classnr[el.Nr()] = 
           SwitchET<ET_TRIG,ET_TET>
           (el.GetType(),
            [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
       });
    
    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);
    Table<size_t> table = creator.MoveTable();


    // size_t ne = ma->GetNE();
  
    for (auto elclass_inds : table)
      {
        if (elclass_inds.Size() == 0) continue;
        
        ElementId ei(VOL,elclass_inds[0]);
        auto & felx = GetFE (ei, lh);
        // auto & trafo = GetMeshAccess()->GetTrafo(ei, lh);
        
        Matrix<> trace_op_x = GetTraceMatrix (felx);
        
        Matrix<> temp_x(elclass_inds.Size(), trace_op_x.Width());
        Matrix<> temp_trace(elclass_inds.Size(), trace_op_x.Height());

        ParallelForRange
          (Range(elclass_inds),
           [&] (IntRange myrange)
           {
             Array<DofId> dofs;
             
             // int tid = TaskManager::GetThreadId();
             {
               for (auto i : myrange)
                 {
                   GetDofNrs(ElementId(VOL,elclass_inds[i]), dofs);
                   in.GetIndirect(dofs, temp_x.Row(i));
                 }
             }
             
             {
               temp_trace.Rows(myrange) = 0;      
               // RegionTracer t(tid, tmulttracex, 0, temp_x.Width());
               AddABt(temp_x.Rows(myrange), trace_op_x, temp_trace.Rows(myrange));
             }
             
             {
               for (auto i : myrange)
                 {
                   tracespace.GetDofNrs(ElementId(VOL,elclass_inds[i]), dofs);
                   out.AddIndirect (dofs, temp_trace.Row(i), true);
                 }
             }
           });
      }
  }
  
  void L2HighOrderFESpace ::
  GetTraceTrans (const FESpace & tracespace, const BaseVector & in, BaseVector & out, bool avg,
                 LocalHeap & lh) const
  {
    static Timer t("GetTraceTrans"); RegionTimer reg(t);

    out = 0.0;
    Array<short> classnr(ma->GetNE());
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         classnr[el.Nr()] = 
           SwitchET<ET_TRIG,ET_TET>
           (el.GetType(),
            [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
       });
    
    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);
    Table<size_t> table = creator.MoveTable();


    // size_t ne = ma->GetNE();
  
    for (auto elclass_inds : table)
      {
        if (elclass_inds.Size() == 0) continue;
        
        ElementId ei(VOL,elclass_inds[0]);
        auto & felx = GetFE (ei, lh);
        // auto & trafo = GetMeshAccess()->GetTrafo(ei, lh);
        
        Matrix<> trace_op_x = GetTraceMatrix (felx);
        
        Matrix<> temp_x(elclass_inds.Size(), trace_op_x.Width());
        Matrix<> temp_trace(elclass_inds.Size(), trace_op_x.Height());

        ParallelForRange
          (Range(elclass_inds),
           [&] (IntRange myrange)
           {
             Array<DofId> dofs;
             
             // int tid = TaskManager::GetThreadId();
             {
               for (auto i : myrange)
                 {
                   tracespace.GetDofNrs(ElementId(VOL,elclass_inds[i]), dofs);
                   in.GetIndirect(dofs, temp_trace.Row(i));
                 }
             }
             
             {
               // temp_trace.Rows(myrange) = 0;      
               // AddABt(temp_x.Rows(myrange), trace_op_x, temp_trace.Rows(myrange));
               temp_x.Rows(myrange) = temp_trace.Rows(myrange) * trace_op_x;
             }
             
             {
               for (auto i : myrange)
                 {
                   GetDofNrs(ElementId(VOL,elclass_inds[i]), dofs);
                   out.AddIndirect (dofs, temp_x.Row(i));
                 }
             }
           });
      }
  }


  template <int D, typename FEL = ScalarFiniteElement<D-1> >
  class DiffOpSurfaceGradient : public DiffOp<DiffOpSurfaceGradient<D, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };

    ///
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      mat = Trans (mip.GetJacobianInverse ()) * 
	Trans (static_cast<const FEL&>(fel).GetDShape(mip.IP(),lh));
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      static_cast<const FEL&>(fel).CalcMappedDShape (mir, mat);      
    }

    using DiffOp<DiffOpSurfaceGradient<D, FEL> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      static_cast<const FEL&>(fel).EvaluateGrad (mir, x, y);
    }

    using DiffOp<DiffOpSurfaceGradient<D, FEL> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      static_cast<const FEL&>(fel).AddGradTrans (mir, y, x);
    }

    static string Name() { return "grad"; }
  };


  L2SurfaceHighOrderFESpace ::  
  L2SurfaceHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    type = "l2surf";
    name="L2SurfaceHighOrderFESpace(l2surf)";
    // defined flags 
    DefineDefineFlag("l2surf");

    if(parseflags) CheckFlags(flags);
    
    if(flags.NumFlagDefined("relorder"))
      throw Exception("Variable order not implemented for L2SurfaceHighOrderFESpace"); 

    /*
    segm = new L2HighOrderFE<ET_SEGM> (order);
    trig = new L2HighOrderFE<ET_TRIG> (order);
    quad = new L2HighOrderFE<ET_QUAD> (order);
    */

    lowest_order_wb = flags.GetDefineFlagX ("lowest_order_wb").IsTrue();

    discontinuous = flags.GetDefineFlagX ("discontinuous").IsTrue();

    if (lowest_order_wb && discontinuous)
      throw Exception("In L2SurfaceFESpace: lowest_order_wb and discontinuous flag are set!");
    
    if (ma->GetDimension() == 2)
      {
        integrator[BND] = 
          make_shared<RobinIntegrator<2>>(make_shared<ConstantCoefficientFunction>(1));
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>(); // for dimension
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>(); // to avoid exception "grad does not exist"
	      flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpSurfaceGradient<2>>>();
      }
    else
      {
	integrator[BND] = 
          make_shared<RobinIntegrator<3>> (make_shared<ConstantCoefficientFunction>(1));
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>(); // to avoid exception "grad does not exist"
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>(); // for dimension
	flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpSurfaceGradient<3>>>();
          
      }

    if (dimension > 1)
    {
      integrator[BND] = make_shared<BlockBilinearFormIntegrator> (integrator[BND], dimension);
      for (auto vb : { VOL,BND, BBND, BBBND })
      {
        if (evaluator[vb])
          evaluator[vb] = make_shared<BlockDifferentialOperator> (evaluator[vb], dimension);
        if (flux_evaluator[vb])
          flux_evaluator[vb] = make_shared<BlockDifferentialOperator> (flux_evaluator[vb], dimension);            
      }
    }
  }

  L2SurfaceHighOrderFESpace :: ~L2SurfaceHighOrderFESpace ()
  {
    ;
  }

  DocInfo L2SurfaceHighOrderFESpace :: GetDocu ()
  {
    DocInfo docu = FESpace::GetDocu(); 
    docu.short_docu = "An L2-conforming finite element space.";
    docu.long_docu =
      R"raw_string(The L2 finite element space on surfaces consists of element-wise polynomials,
which are discontinuous from element to element. It uses an
L2-orthogonal hierarchical basis which leads to orthogonal
mass-matrices on non-curved elements.

The L2 space supports element-wise variable order, which can be set
for ELEMENT-nodes.

Per default, all dofs are local dofs and are condensed if static
condensation is performed. The lowest order can be kept in the
WIRE_BASKET via the flag 'lowest_order_wb=True'.

)raw_string";

    docu.Arg("lowest_order_wb") = "bool = False\n"
      "  Keep lowest order dof in WIRE_BASKET and make other dofs LOCAL";
    docu.Arg("discontinuous") = "bool = False\n"
      "  Make all dofs LOCAL";

    return docu;
  }

  shared_ptr<FESpace> L2SurfaceHighOrderFESpace :: 
  Create (shared_ptr<MeshAccess> ma, const Flags & flags)
  {
    return make_shared<L2SurfaceHighOrderFESpace> (ma, flags, true);
  }

  void L2SurfaceHighOrderFESpace :: Update()
  {
    nel = ma->GetNSE();

    bool first_update = GetTimeStamp() < ma->GetTimeStamp();
    if (first_update) timestamp = NGS_Object::GetNextTimeStamp();
    
    
    if (first_update)
      {
	order_inner.SetSize(nel);
	order_inner = INT<3>(order);

	for(int i = 0; i < nel; i++) 
	  {
	    ElementId ei(BND,i);
	    order_inner[i] = order_inner[i] + INT<3> (et_bonus_order[ma->GetElType(ei)]);
	    order_inner[i] = Max(order_inner[i], INT<3>(0));
	    if (!DefinedOn (ei))
	      order_inner[i] = 0;
	  }
    
	if(print) 
	  *testout << " order_inner (l2surf) " << order_inner << endl;
      }
    
    ndof = 0;
    first_element_dof.SetSize(nel+1);
    for (int i = 0; i < nel; i++)
      {
	first_element_dof[i] = ndof;
	INT<3> pi = order_inner[i];
	switch (ma->GetElType(ElementId(BND, i)))
	  {
	  case ET_SEGM:
	    ndof += pi[0]+1;
	    break;
	  case ET_TRIG:
	    ndof += (pi[0]+1)*(pi[1]+2)/2;
	    break;
	  case ET_QUAD:
	    ndof += (pi[0]+1)*(pi[1]+1);
	    break;
	  default:
	    ;
	  }
      }
    first_element_dof[nel] = ndof;

    if(print) 
      *testout << " first_element dof (l2surf) " << first_element_dof << endl;
    
    UpdateCouplingDofArray();
  }

  void L2SurfaceHighOrderFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(ndof);
    ctofdof = UNUSED_DOF;
    
    for (auto i : Range(ma->GetNSE()))
      if (DefinedOn({BND,i}))
        {
          auto r = GetElementDofs(i);
          ctofdof[r] =  (discontinuous || lowest_order_wb) ? LOCAL_DOF : WIREBASKET_DOF;
          
          if (lowest_order_wb && r.Size() != 0)
            ctofdof[r.First()] = WIREBASKET_DOF;
        }
  }


  void L2SurfaceHighOrderFESpace ::SetOrder (NodeId ni, int order)
  {
    if (order_policy == CONSTANT_ORDER || order_policy == NODE_TYPE_ORDER)
      throw Exception("In L2SurfaceHighOrderFESpace::SetOrder. Order policy is constant or node-type!");
    else if (order_policy == OLDSTYLE_ORDER)
      order_policy = VARIABLE_ORDER;
      
    if (order < 0)
      order = 0;

    if (CoDimension(ni.GetType(), ma->GetDimension()) == 1)
      {
	if (ma->GetDimension() == 3 && ni.GetType() == NT_FACE)
	  {
	    Array<int> elnr;
	    ma->GetFacetSurfaceElements(ni.GetNr(),elnr);
	    if (elnr[0] < order_inner.Size())
	      order_inner[elnr[0]] = order;
	  }
        else if (ni.GetNr() < order_inner.Size())
          order_inner[ni.GetNr()] = order;
      }
    else
      throw Exception ("L2SurfaceHighOrderFESpace::SetOrder requires NodeType of codimension 1!");
  }
  
  int L2SurfaceHighOrderFESpace ::GetOrder (NodeId ni) const
  {
    if (CoDimension(ni.GetType(), ma->GetDimension()) == 1)
      {
	if (ma->GetDimension() == 3 && ni.GetType() == NT_FACE)
	  {
	    Array<int> elnr;
	    ma->GetFacetSurfaceElements(ni.GetNr(),elnr);
	    if (elnr[0] < order_inner.Size())
	      return order_inner[elnr[0]][0];
	  }
        else if (ni.GetNr() < order_inner.Size())
          return order_inner[ni.GetNr()][0];
      }
        
    return 0;
  }
  
  // const FiniteElement & L2SurfaceHighOrderFESpace :: GetSFE (int elnr, LocalHeap & lh) const
  // {
  //   if (ma->GetDimension() == 2)
  //     {
  //       DGFiniteElement<1> * fe1d = 0;
	
  //       Ngs_Element ngel = ma->GetElement<1,BND> (elnr);

  //       switch (ngel.GetType())
  //         {
  //         case ET_SEGM: fe1d = new (lh) L2HighOrderFE<ET_SEGM> (); break;
  //         default:
  //           ;
  //         }

  //       fe1d -> SetVertexNumbers (ngel.vertices);
  //       fe1d -> SetOrder (INT<1> (order));
  //       fe1d -> ComputeNDof(); 
  //       return *fe1d;
  //     }
  //   else
  //     {
  //       DGFiniteElement<2> * fe2d = 0;
	
  //       Ngs_Element ngel = ma->GetElement<2,BND> (elnr);
	
  //       switch (ngel.GetType())
  //         {
  //         case ET_TRIG: fe2d = new (lh) L2HighOrderFE<ET_TRIG> (); break;
  //         case ET_QUAD: fe2d = new (lh) L2HighOrderFE<ET_QUAD> (); break;
  //         default:
  //           ;
  //         }
	
  //       fe2d -> SetVertexNumbers (ngel.vertices);
  //       fe2d -> SetOrder (INT<2> (order));
  //       fe2d -> ComputeNDof(); 
  //       return *fe2d;
  //     }
  //       /*

  //   FiniteElement * fe = 0;
    
  //   switch (ma->GetSElType(elnr))
  //     {
  //     case ET_SEGM:
  //       fe = segm; break;
  //     case ET_TRIG:
  //       fe = trig; break;
  //     case ET_QUAD:
  //       fe = quad; break;
  //     default:
  //       fe = 0;
  //     }
    
  //   ArrayMem<int,12> vnums; // calls GetElPNums -> max 12 for PRISM12
  //   ma->GetSElVertices(elnr, vnums);

  //   if (!fe)
  //     {
  //       stringstream str;
  //       str << "L2SurfaceHighOrderFESpace " << GetClassName() 
  //           << ", undefined eltype " 
  //           << ElementTopology::GetElementName(ma->GetSElType(elnr))
  //           << ", order = " << order << endl;
  //       throw Exception (str.str());
  //     }

  //   dynamic_cast<L2HighOrderFiniteElement<2>*> (fe) -> SetVertexNumbers (vnums);

  //   return *fe;
  //       */
  // }

  FiniteElement & L2SurfaceHighOrderFESpace :: GetFE (ElementId ei, Allocator & lh) const
  {
    if (ei.VB() == BND && DefinedOn(ei))
      {
        if (ma->GetDimension() == 2)
          {
            DGFiniteElement<1> * fe1d = 0;
	
            Ngs_Element ngel = ma->GetElement<1,BND> (ei.Nr());

            switch (ngel.GetType())
              {
              case ET_SEGM: fe1d = new (lh) L2HighOrderFE<ET_SEGM> (); break;
              default:
                ;
              }

            fe1d -> SetVertexNumbers (ngel.vertices);
            fe1d -> SetOrder (order_inner[ei.Nr()]);
            fe1d -> ComputeNDof(); 
            return *fe1d;
          }
        else
          {
            DGFiniteElement<2> * fe2d = 0;
	
            Ngs_Element ngel = ma->GetElement<2,BND> (ei.Nr());
	
            switch (ngel.GetType())
              {
              case ET_TRIG: fe2d = new (lh) L2HighOrderFE<ET_TRIG> (); break;
              case ET_QUAD: fe2d = new (lh) L2HighOrderFE<ET_QUAD> (); break;
              default:
                ;
              }
	
            fe2d -> SetVertexNumbers (ngel.vertices);
            fe2d -> SetOrder (order_inner[ei.Nr()]);
            fe2d -> ComputeNDof(); 
            return *fe2d;
          }
      }
    
    else

      return * SwitchET (ma->GetElement(ei).GetType(),
                         [&lh] (auto et) -> FiniteElement*
                         {
                           return new (lh) ScalarDummyFE<et.ElementType()>();
                         });
    
  }
  // const FiniteElement & L2SurfaceHighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  // {
  //   throw Exception ("Volume elements not available for L2SurfaceHighOrderFESpace");
  // }
 
  size_t L2SurfaceHighOrderFESpace :: GetNDof () const throw()
  {
    return ndof;
  }

  void L2SurfaceHighOrderFESpace :: 
  GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if (ei.VB()!=BND || !DefinedOn (ei)) return;
 
    int first = first_element_dof[ei.Nr()];
    int neldofs = first_element_dof[ei.Nr()+1] - first;
    for (int j = 0; j < neldofs; j++)
      dnums.Append (first+j);
  }
  
  void L2SurfaceHighOrderFESpace :: SolveM (CoefficientFunction * rho, BaseVector & vec, Region * def,
                                  LocalHeap & lh) const
  {
    static Timer t("SolveM"); RegionTimer reg(t);
    if (rho && rho->Dimension() != 1)
      throw Exception("L2HighSurfaceOrderFESpace::SolveM needs a scalar density");
    IterateElements (*this, BND, lh,
                     [&rho, &vec, def, this] (FESpace::Element el, LocalHeap & lh)
                     {
                       auto & fel = static_cast<const BaseScalarFiniteElement&>(el.GetFE());
                       const ElementTransformation & trafo = el.GetTrafo();
                       
                       Array<int> dnums(fel.GetNDof(), lh);
                       GetDofNrs (ElementId(BND,el.Nr()), dnums);

                       FlatVector<double> elx(fel.GetNDof()*dimension, lh);
                       if (def && !def->Mask()[ma->GetElIndex(el)])
                         {
                           elx = 0.0;
                           vec.SetIndirect (dnums, elx);
                           return;
                         }
                       
                       vec.GetIndirect(dnums, elx);
                       auto melx = elx.AsMatrix(fel.GetNDof(),dimension);

                       FlatVector<double> diag_mass(fel.GetNDof(), lh);
                       fel.GetDiagMassMatrix (diag_mass);

                       bool curved = trafo.IsCurvedElement();
                       if (rho && !rho->ElementwiseConstant()) curved = true;

                       if (!curved)
                         {
                           IntegrationRule ir(fel.ElementType(), 0);
                           BaseMappedIntegrationRule & mir = trafo(ir, lh);
                           double jac = mir[0].GetMeasure();
                           if (rho) jac *= rho->Evaluate(mir[0]);
                           diag_mass *= jac;
                           for (int i = 0; i < melx.Height(); i++)
                             melx.Row(i) /= diag_mass(i);
                         }
                       else
                         {
                           SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
                           auto & mir = trafo(ir, lh);
                           FlatVector<SIMD<double>> pntvals(ir.Size(), lh);
                           FlatMatrix<SIMD<double>> rhovals(1, ir.Size(), lh);
                           if (rho) rho->Evaluate (mir, rhovals);
                           
                           for (int i = 0; i < melx.Height(); i++)
                             melx.Row(i) /= diag_mass(i);
                           for (int comp = 0; comp < dimension; comp++)
                             {
                               fel.Evaluate (ir, melx.Col(comp), pntvals);
                               if (rho)
                                 for (size_t i = 0; i < ir.Size(); i++)
                                   pntvals(i) *= ir[i].Weight() / (mir[i].GetMeasure() * rhovals(0,i));
                               else
                                 for (size_t i = 0; i < ir.Size(); i++)
                                   pntvals(i) *= ir[i].Weight() / mir[i].GetMeasure();
                               
                               melx.Col(comp) = 0.0;
                               fel.AddTrans (ir, pntvals, melx.Col(comp));
                             }
                           for (int i = 0; i < melx.Height(); i++)
                             melx.Row(i) /= diag_mass(i);
                         }
                       vec.SetIndirect(dnums, elx);
                     });
  }

  void L2SurfaceHighOrderFESpace :: ApplyM (CoefficientFunction * rho, BaseVector & vec, Region * def,
                                     LocalHeap & lh) const
  {
    static Timer t("ApplyM"); RegionTimer reg(t);
    if (rho && rho->Dimension() != 1)
      throw Exception("L2HighOrderFESpace::ApplyM needs a scalar density");

    // auto fv = vec.FV<double>();
    
    IterateElements (*this, BND, lh,
                     [&rho, &vec, /* fv, */ def, this] (FESpace::Element el, LocalHeap & lh)
                     {
                       // auto tid = TaskManager::GetThreadId();
                       auto & fel = static_cast<const BaseScalarFiniteElement&>(el.GetFE());
                       const ElementTransformation & trafo = el.GetTrafo();
                       
                       Array<int> dnums(fel.GetNDof(), lh);
                       // auto dofrange = GetElementDofs(el.Nr());
                       FlatVector<double> elx(fel.GetNDof()*dimension, lh);

                       //bool lindofs = all_dofs_together && dimension==1;

                       if (def && !def->Mask()[ma->GetElIndex(el)])
                         {
                           /*if (lindofs)
                             fv.Range(dofrange) = 0.0;
                           else
                             {*/
                               elx = 0.0;
                               GetDofNrs (el, dnums);                               
                               vec.SetIndirect(dnums, elx);
                             //}
                           return;
                         }
                       
                       /*if (!lindofs)
                         {*/
                           GetDofNrs (el, dnums);
                           vec.GetIndirect(dnums, elx);
                         /*}
                       else
                         elx = fv.Range(dofrange);
                       */
                       
		                   auto melx = elx.AsMatrix(fel.GetNDof(),dimension);
                       
                       bool curved = trafo.IsCurvedElement();
                       if (rho && !rho->ElementwiseConstant()) curved = true;
                       if (!curved)
                         {
                           FlatVector<double> diag_mass(fel.GetNDof(), lh);
                           fel.GetDiagMassMatrix (diag_mass);
                           IntegrationRule ir(fel.ElementType(), 0);
                           BaseMappedIntegrationRule & mir = trafo(ir, lh);
                           double jac = mir[0].GetMeasure();
                           if (rho) jac *= rho->Evaluate(mir[0]);
                           
                           

                           if (dimension == 1)
                             for (size_t i = 0; i < elx.Size(); i++)
                               elx(i) *= jac*diag_mass(i);
                           else
                             for (size_t i = 0; i < melx.Height(); i++)
                               melx.Row(i) *= jac*diag_mass(i);
                         }
                       else
                         {
                           SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
                           auto & mir = trafo(ir, lh);
                           FlatVector<SIMD<double>> pntvals(ir.Size(), lh);
                           FlatMatrix<SIMD<double>> rhovals(1, ir.Size(), lh);
                           if (rho) rho->Evaluate (mir, rhovals);
                           
                           for (int comp = 0; comp < dimension; comp++)
                             {
                               fel.Evaluate (ir, melx.Col(comp), pntvals);
                               if (rho)
                                 for (size_t i = 0; i < ir.Size(); i++)
                                   pntvals(i) *= ir[i].Weight()* rhovals(0,i) * mir[i].GetMeasure() ;
                               else
                                 for (size_t i = 0; i < ir.Size(); i++)
                                   pntvals(i) *= ir[i].Weight() * mir[i].GetMeasure();
                               
                               melx.Col(comp) = 0.0;
                               fel.AddTrans (ir, pntvals, melx.Col(comp));
                             }
                         }
                       
                       //if (!lindofs)
                         vec.SetIndirect(dnums, elx);
                       //else
                        // fv.Range(dofrange) = elx;
                       
                     });
  }
  shared_ptr<Table<int>> L2SurfaceHighOrderFESpace :: 
  // CreateSmoothingBlocks ( int type) const
  CreateSmoothingBlocks (const Flags & precflags) const    
  {
    int i, j, first;
    Array<int> cnt(nel);
    cnt = 0;
    for (i = 0; i < nel; i++)
      cnt[i] = first_element_dof[i+1]-first_element_dof[i];
	
    Table<int> table(cnt);
    
    for (i = 0; i < nel; i++)
      {
	first = first_element_dof[i];
	for (j = 0; j < cnt[i]; j++)
	  table[i][j] = first+j;
      }
    return make_shared<Table<int>> (table);
  }


  void  L2SurfaceHighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  { dnums.SetSize0(); return; }
  
  void  L2SurfaceHighOrderFESpace ::GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  { dnums.SetSize0(); return; }
  
  void  L2SurfaceHighOrderFESpace ::GetFaceDofNrs (int fanr, Array<int> & dnums) const
  { GetDofNrs ( fanr, dnums ); return; }
  
  void  L2SurfaceHighOrderFESpace ::GetInnerDofNrs (int elnr, Array<int> & dnums) const
  { GetDofNrs ( elnr, dnums ); return; }

  template <int DIM_SPC, VorB VB = VOL>
  class DiffOpIdVectorL2Piola : public DiffOp<DiffOpIdVectorL2Piola<DIM_SPC> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = DIM_SPC };
    enum { DIM_ELEMENT = DIM_SPC-VB };
    enum { DIM_DMAT = DIM_SPC };
    enum { DIFFORDER = 0 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      mat = 0.0;
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      for (int i = 0; i < DIM_SPACE; i++)
        feli.CalcShape (mip.IP(), mat.Row(i).Range(fel.GetRange(i)));
      Mat<DIM_SPACE> trafo = (1.0/mip.GetJacobiDet()) * mip.GetJacobian();
      for (int i = 0; i < DIM_SPACE*feli.GetNDof(); i++)
        {
          Vec<DIM_SPACE> hv = mat.Col(i);
          mat.Col(i) = trafo * hv;
        } 
    }
    
    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir, 
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM_SPACE,DIM_SPACE>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);

      size_t ndofi = feli.GetNDof();
      auto scalmat = mat.Rows( (sqr(DIM_SPC)-1)*ndofi, sqr(DIM_SPC)*ndofi);
      feli.CalcShape (mir.IR(), scalmat);
      
      for (auto i_ip : Range(mir))
        {
          auto & mip = mir[i_ip];
          Mat<DIM_SPACE,DIM_SPACE,SIMD<double>> trafo = (1.0/mip.GetJacobiDet()) * mip.GetJacobian();

          auto col = mat.Col(i_ip);
          auto scalcol = scalmat.Col(i_ip);

          size_t base = 0;
          for (size_t k = 0; k < DIM_SPACE; k++)
            for (size_t i = 0; i < feli.GetNDof(); i++, base += DIM_SPACE)
              col.Range(base, base+DIM_SPACE) = scalcol(i) * trafo.Col(k);
        }
    }
    
    
    using DiffOp<DiffOpIdVectorL2Piola<DIM_SPC>>::ApplySIMDIR;        
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      // static Timer t("DiffOpIdVectorL2Piola::ApplySIMDIR");
      // RegionTracer rt(TaskManager::GetThreadId(), t);
      
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPC>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      size_t ndofi = feli.GetNDof();
      
      STACK_ARRAY(double, memx, DIM_SPACE*ndofi);
      FlatMatrixFixWidth<DIM_SPACE, double> matx(ndofi, &memx[0]);
      for (size_t k = 0; k < DIM_SPACE; k++)
        matx.Col(k) = x.Range(k*ndofi, (k+1)*ndofi);

      {
        // RegionTracer rt(TaskManager::GetThreadId(), t);                      
        feli.Evaluate(mir.IR(), matx, y);
      }
      for (size_t i = 0; i < mir.Size(); i++)
        {
          auto jac = mir[i].GetJacobian();
          Vec<DIM_SPACE,SIMD<double>> val = y.Col(i);
          val *= 1/mir[i].GetJacobiDet();
          y.Col(i).Range(0,DIM_SPACE) = jac * val;
        }
    }    
    
    
    using DiffOp<DiffOpIdVectorL2Piola<DIM_SPC>>::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      //static Timer t("DiffpIdVectorL2Piola::AddTransSIMD");
      // static Timer tc("DiffpIdVectorL2Piola::AddTransSIMD calc");
      // RegionTracer rt(TaskManager::GetThreadId(), t);            

      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPC>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      size_t ndofi = feli.GetNDof();
      
      STACK_ARRAY(SIMD<double>, mempt, mir.Size()*DIM_SPACE);
      FlatMatrix<SIMD<double>> hy(DIM_SPACE, mir.Size(), &mempt[0]);

      for (size_t i = 0; i < mir.Size(); i++)
        {
          auto jac = mir[i].GetJacobian();
          Vec<DIM_SPACE,SIMD<double>> val = y.Col(i);
          val *= 1/mir[i].GetJacobiDet();
          hy.Col(i) = Trans(jac) * val;
        }
      
      STACK_ARRAY(double, memx, DIM_SPACE*ndofi);
      FlatMatrixFixWidth<DIM_SPACE, double> matx(ndofi, &memx[0]);

      for (size_t i = 0; i < ndofi; i++)
        matx.Row(i) = x.Slice(i, ndofi);

      feli.AddTrans(mir.IR(), hy, matx);

      for (size_t k = 0; k < DIM_SPACE; k++)
        x.Range(k*ndofi, (k+1)*ndofi) = matx.Col(k);
    }    

    
/*   

  using DiffOp<DiffOpIdVectorH1<DIM_SPC>>::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      for (int i = 0; i < DIM_SPC; i++)
        {
          auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[i]);
          feli.Evaluate (mir.IR(), x.Range(fel.GetRange(i)), y.Row(i));
        }
    }

    using DiffOp<DiffOpIdVectorH1<DIM_SPC>>::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      for (int i = 0; i < DIM_SPC; i++)
        {
          auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[i]);
          feli.AddTrans (mir.IR(), y.Row(i), x.Range(fel.GetRange(i)));
        }
    }    
    */
  };




  template <int DIM_SPC>
  class DiffOpDivVectorL2Piola : public DiffOp<DiffOpDivVectorL2Piola<DIM_SPC>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = DIM_SPC };
    enum { DIM_ELEMENT = DIM_SPC };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };

    static string Name() { return "div"; }
    
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      auto & bfel = static_cast<const CompoundFiniteElement&> (fel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (bfel[0]);
      
      int ndofi = feli.GetNDof();
      FlatMatrix<> grad (ndofi, DIM_SPC, lh);
      feli.CalcDShape(mip.IP(), grad);
      double idet = 1.0/mip.GetJacobiDet();

      for (int k = 0; k < DIM_SPC; k++)
        mat.Row(0).Range(k*ndofi, (k+1)*ndofi) = idet * grad.Col(k);
    }


    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir, 
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM_SPACE,DIM_SPACE>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);

      size_t ndofi = feli.GetNDof();
      feli.CalcMappedDShape (mir, mat);

      STACK_ARRAY (SIMD<double>, mem, ndofi*DIM_SPC);
      FlatMatrix<SIMD<double>> tmp(ndofi, DIM_SPC, &mem[0]);
      for (auto i_ip : Range(mir))
        {
          auto col = mat.Col(i_ip);
          auto & mip = mir[i_ip];
          
          for (size_t i = 0; i < ndofi; i++)
            for (size_t k = 0; k < DIM_SPC; k++)
              tmp(i, k) = col(i*DIM_SPC+k);

          for (size_t k = 0; k < DIM_SPC; k++)
            {
              Vec<DIM_SPC, SIMD<double>> dir = 1.0/mip.GetJacobiDet() * mip.GetJacobian().Col(k);
              col.Range(k*ndofi, (k+1)*ndofi) = tmp * dir;
            }
        }
    }

    

    /*
    using DiffOp<DiffOpCurlVectorL2Covariant>::ApplySIMDIR;        
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<3,3>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      size_t ndofi = feli.GetNDof();
      y.AddSize(3,mir.Size()) = SIMD<double>(0.0);

      STACK_ARRAY(SIMD<double>, mem, 3*mir.Size());
      FlatMatrix<SIMD<double>> grad(3, mir.Size(), &mem[0]);
      for (size_t k = 0; k < 3; k++)
        {
          feli.EvaluateGrad (mir, x.Range(k*ndofi, (k+1)*ndofi), grad);
          for (size_t i = 0; i < mir.Size(); i++)
            {
              auto trafo = Trans(mir[i].GetJacobianInverse());
              Vec<3,SIMD<double>> gi = grad.Col(i);
              Vec<3,SIMD<double>> tek = trafo.Col(k);
              Vec<3,SIMD<double>> hv = Cross(gi, tek);
              y.Col(i).AddSize(3) += hv;
            }
        }             
    }    
      
    using DiffOp<DiffOpCurlVectorL2Covariant>::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<3,3>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      size_t ndofi = feli.GetNDof();

      STACK_ARRAY(SIMD<double>, mem, 3*mir.Size());
      FlatMatrix<SIMD<double>> grad(3, mir.Size(), &mem[0]);
      for (size_t k = 0; k < 3; k++)
        {
          for (size_t i = 0; i < mir.Size(); i++)
            {
              auto trafo = Trans(mir[i].GetJacobianInverse());
              Vec<3,SIMD<double>> cy = y.Col(i);
              Vec<3,SIMD<double>> tek = trafo.Col(k);
              Vec<3,SIMD<double>> hv = Cross(cy, tek);
              grad.Col(i) = -hv;
            }
          
          feli.AddGradTrans (mir, grad, x.Range(k*ndofi, (k+1)*ndofi));
        }             
    }    
    */
  };
  


  template <int DIM_SPC>
  class DiffOpGradVectorL2Piola : public DiffOp<DiffOpGradVectorL2Piola<DIM_SPC>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = DIM_SPC };
    enum { DIM_ELEMENT = DIM_SPC };
    enum { DIM_DMAT = DIM_SPC*DIM_SPC };
    enum { DIFFORDER = 1 };

    static string Name() { return "grad"; }

    static Array<int> GetDimensions() { return Array<int> ( { DIM_SPC, DIM_SPC } ); };
    
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      auto & bfel = static_cast<const CompoundFiniteElement&> (fel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (bfel[0]);
      
      int ndofi = feli.GetNDof();
      FlatMatrix<> grad (ndofi, DIM_SPC, lh);
      feli.CalcMappedDShape(mip, grad);

      Mat<DIM_SPC,DIM_SPC> trans = 1/(mip.GetJacobiDet())*mip.GetJacobian();

      for (int i = 0; i < DIM_SPC; i++)
        for (int j = 0; j < DIM_SPC; j++)
          for (int k = 0; k < DIM_SPC; k++)
            mat.Row(i*DIM_SPC+j).Range(k*ndofi, (k+1)*ndofi) = trans(i,k) * grad.Col(j);

      FlatVector<> val (ndofi, lh);
      feli.CalcShape(mip.IP(), val);

      // 1/J ( H - (F^{-T}:H) F ) 
      if (!mip.GetTransformation().IsCurvedElement())
        return;

      auto jac = mip.GetJacobian();
      auto inv = mip.GetJacobianInverse();
      auto invJ = 1/mip.GetJacobiDet();
      Vec<DIM_SPC, Mat<DIM_SPC,DIM_SPC>> hesse;
      mip.CalcHesse(hesse);

      Vec<DIM_SPC, Mat<DIM_SPC,DIM_SPC>> invjac_hesse;
      for (int i = 0; i < DIM_SPC; i++)
        invjac_hesse(i) = Trans(inv) * hesse(i);
      
      Vec<DIM_SPC> inv_hesse = 0.0;
      for (int i = 0; i < DIM_SPC; i++)
        for (int j = 0; j < DIM_SPC; j++)
          inv_hesse(i) += invjac_hesse(j)(j,i);
      inv_hesse = Trans(inv) * inv_hesse;
      
      for (int i = 0; i < DIM_SPC; i++)
        for (int j = 0; j < DIM_SPC; j++)
          for (int k = 0; k < DIM_SPC; k++)
            mat.Row(i*DIM_SPC+j).Range(k*ndofi, (k+1)*ndofi) +=
              invJ * (invjac_hesse(i)(j,k)-inv_hesse(j)*jac(i,k)) * val;
    }

    /*
    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir, 
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM_SPACE,DIM_SPACE>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);

      size_t ndofi = feli.GetNDof();
      feli.CalcMappedDShape (mir, mat);

      STACK_ARRAY (SIMD<double>, mem, ndofi*DIM_SPC);
      FlatMatrix<SIMD<double>> tmp(ndofi, DIM_SPC, &mem[0]);
      for (auto i_ip : Range(mir))
        {
          auto col = mat.Col(i_ip);
          auto & mip = mir[i_ip];
          
          for (size_t i = 0; i < ndofi; i++)
            for (size_t k = 0; k < DIM_SPC; k++)
              tmp(i, k) = col(i*DIM_SPC+k);

          for (size_t k = 0; k < DIM_SPC; k++)
            {
              Vec<DIM_SPC, SIMD<double>> dir = 1.0/mip.GetJacobiDet() * mip.GetJacobian().Col(k);
              col.Range(k*ndofi, (k+1)*ndofi) = tmp * dir;
            }
        }
    }
    */
    

    /*
    using DiffOp<DiffOpCurlVectorL2Covariant>::ApplySIMDIR;        
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<3,3>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      size_t ndofi = feli.GetNDof();
      y.AddSize(3,mir.Size()) = SIMD<double>(0.0);

      STACK_ARRAY(SIMD<double>, mem, 3*mir.Size());
      FlatMatrix<SIMD<double>> grad(3, mir.Size(), &mem[0]);
      for (size_t k = 0; k < 3; k++)
        {
          feli.EvaluateGrad (mir, x.Range(k*ndofi, (k+1)*ndofi), grad);
          for (size_t i = 0; i < mir.Size(); i++)
            {
              auto trafo = Trans(mir[i].GetJacobianInverse());
              Vec<3,SIMD<double>> gi = grad.Col(i);
              Vec<3,SIMD<double>> tek = trafo.Col(k);
              Vec<3,SIMD<double>> hv = Cross(gi, tek);
              y.Col(i).AddSize(3) += hv;
            }
        }             
    }    
    */
    using DiffOp<DiffOpGradVectorL2Piola>::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM_SPACE,DIM_SPACE>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      size_t ndofi = feli.GetNDof();

      STACK_ARRAY(SIMD<double>, mem, DIM_SPC*DIM_SPC*mir.Size());
      FlatMatrix<SIMD<double>> grad(DIM_SPC*DIM_SPC, mir.Size(), &mem[0]);
      grad = SIMD<double>(0.0);

      for (size_t i = 0; i < mir.Size(); i++)
        {
          auto trans = 1/(mir[i].GetJacobiDet())*mir[i].GetJacobian();
          for (int j = 0; j < DIM_SPC; j++)    
            for (int k = 0; k < DIM_SPC; k++)
              for (int l = 0; l < DIM_SPC; l++)
                grad(j*DIM_SPC+k, i) += trans(l,j)*y(k*DIM_SPC+l, i);
        }
          
      for (size_t k = 0; k < DIM_SPC; k++)
        feli.AddGradTrans (mir, grad.Rows(k*DIM_SPC, (k+1)*DIM_SPC), x.Range(k*ndofi, (k+1)*ndofi));

      if (!mir.GetTransformation().IsCurvedElement())
        return;

      STACK_ARRAY(SIMD<double>, mem2, DIM_SPC*mir.Size());
      FlatMatrix<SIMD<double>> val(DIM_SPC, mir.Size(), &mem2[0]);
      val = SIMD<double>(0.0);

      for (size_t ip = 0; ip < mir.Size(); ip++)
        {
          auto jac = mir[ip].GetJacobian();
          auto inv = mir[ip].GetJacobianInverse();
          auto invJ = 1/mir[ip].GetJacobiDet();
          Vec<DIM_SPC, Mat<DIM_SPC,DIM_SPC,SIMD<double>>> hesse;
          mir[ip].CalcHesse(hesse);
          
          Vec<DIM_SPC, Mat<DIM_SPC,DIM_SPC,SIMD<double>>> invjac_hesse;
          for (int i = 0; i < DIM_SPC; i++)
            invjac_hesse(i) = Trans(inv) * hesse(i);
          
          Vec<DIM_SPC,SIMD<double>> inv_hesse = SIMD<double>(0.0);
          for (int i = 0; i < DIM_SPC; i++)
            for (int j = 0; j < DIM_SPC; j++)
              inv_hesse(i) += invjac_hesse(j)(j,i);
          inv_hesse = Trans(inv) * inv_hesse;
          
          for (int i = 0; i < DIM_SPC; i++)
            for (int j = 0; j < DIM_SPC; j++)
              for (int k = 0; k < DIM_SPC; k++)
                val(k,ip) +=
                  invJ * (invjac_hesse(i)(j,k)-inv_hesse(j)*jac(i,k)) *
                  y(i*DIM_SPC+j,ip);
        }

      for (size_t k = 0; k < DIM_SPC; k++)
        feli.AddTrans (mir.IR(), val.Row(k), x.Range(k*ndofi, (k+1)*ndofi));
    }    
  };




  
  template <int DIM_SPC, VorB VB = VOL>
  class DiffOpIdVectorL2Covariant : public DiffOp<DiffOpIdVectorL2Covariant<DIM_SPC> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = DIM_SPC };
    enum { DIM_ELEMENT = DIM_SPC-VB };
    enum { DIM_DMAT = DIM_SPC };
    enum { DIFFORDER = 0 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      mat = 0.0;
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      for (int i = 0; i < DIM_SPACE; i++)
        feli.CalcShape (mip.IP(), mat.Row(i).Range(fel.GetRange(i)));
      Mat<DIM_SPACE> trafo = mip.GetJacobianInverse();
      for (int i = 0; i < DIM_SPACE*feli.GetNDof(); i++)
        {
          Vec<DIM_SPACE> hv = mat.Col(i);
          mat.Col(i) = Trans(trafo) * hv;
        } 
    }
    
    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir, 
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      size_t ndofi = feli.GetNDof();
      
      STACK_ARRAY(SIMD<double>, mem, ndofi*mir.Size());
      FlatMatrix<SIMD<double>> shapes(ndofi, mir.Size(), &mem[0]);
      feli.CalcShape (mir.IR(), shapes);

      for (auto i_ip : Range(mir))
        {
          auto col = mat.Col(i_ip);
          auto & mip = static_cast<const SIMD<ngfem::MappedIntegrationPoint<DIM_ELEMENT,DIM_SPC>>&>(mir[i_ip]);
          auto trafo = mip.GetJacobianInverse();
          
          for (int k = 0; k < DIM_SPACE; k++)
            {
              size_t offset = DIM_SPACE*k*ndofi;
              for (size_t i = 0; i < feli.GetNDof(); i++, offset += DIM_SPACE)
                col.Range(offset,offset+DIM_SPACE) = shapes(i,i_ip) * trafo.Row(k);
	    }
        }
    }


    using DiffOp<DiffOpIdVectorL2Covariant<DIM_SPC>>::ApplySIMDIR;        
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPC>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      size_t ndofi = feli.GetNDof();
      
      STACK_ARRAY(double, memx, DIM_SPACE*ndofi);
      FlatMatrix<double> matx(ndofi, DIM_SPACE, &memx[0]);
      for (size_t k = 0; k < DIM_SPACE; k++)
        matx.Col(k) = x.Range(k*ndofi, (k+1)*ndofi);
      
      feli.Evaluate(mir.IR(), matx, y);

      for (size_t i = 0; i < mir.Size(); i++)
        {
          auto jacinv = mir[i].GetJacobianInverse();
          Vec<DIM_SPACE,SIMD<double>> val = y.Col(i);
          y.Col(i).Range(0,DIM_SPACE) = Trans(jacinv) * val;
        }
    }    
    

    using DiffOp<DiffOpIdVectorL2Covariant<DIM_SPC>>::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPC>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      size_t ndofi = feli.GetNDof();
      
      STACK_ARRAY(SIMD<double>, mempt, mir.Size()*DIM_SPACE);
      FlatMatrix<SIMD<double>> hy(DIM_SPACE, mir.Size(), &mempt[0]);

      for (size_t i = 0; i < mir.Size(); i++)
        {
          auto jacinv = mir[i].GetJacobianInverse();
          Vec<DIM_SPACE,SIMD<double>> val = y.Col(i);
          hy.Col(i) = jacinv * val;
        }
      
      STACK_ARRAY(double, memx, DIM_SPACE*ndofi);
      FlatMatrix<double> matx(ndofi, DIM_SPACE, &memx[0]);

      for (size_t k = 0; k < DIM_SPACE; k++)
        matx.Col(k) = x.Range(k*ndofi, (k+1)*ndofi);

      feli.AddTrans(mir.IR(), hy, matx);

      for (size_t k = 0; k < DIM_SPACE; k++)
        x.Range(k*ndofi, (k+1)*ndofi) = matx.Col(k);
    }    
  };


  class DiffOpCurlVectorL2Covariant : public DiffOp<DiffOpCurlVectorL2Covariant>
  {
  public:
    enum { DIM = 3 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 3 };
    enum { DIM_DMAT = 3 };
    enum { DIFFORDER = 1 };

    static string Name() { return "curl"; }
    
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      auto & bfel = static_cast<const CompoundFiniteElement&> (fel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (bfel[0]);
      mat = 0;
      
      int ndofi = feli.GetNDof();
      FlatMatrix<> grad (ndofi, 3, lh);
      feli.CalcDShape(mip.IP(), grad);
      Mat<DIM_SPACE> trafo = (1.0/mip.GetJacobiDet()) * mip.GetJacobian();

      for (int k = 0; k < 3; k++)
        for (int i = 0; i < ndofi; i++)
          {
            Vec<3> gi = grad.Row(i);
            Vec<3> ek = 0.0; ek(k) = 1;
            Vec<3> hv = Cross(gi, ek);
            mat.Col(k*ndofi+i) = trafo * hv;
          }
    }

    using DiffOp<DiffOpCurlVectorL2Covariant>::ApplySIMDIR;        
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<3,3>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      size_t ndofi = feli.GetNDof();
      y.AddSize(3,mir.Size()) = SIMD<double>(0.0);

      STACK_ARRAY(SIMD<double>, mem, 3*mir.Size());
      FlatMatrix<SIMD<double>> grad(3, mir.Size(), &mem[0]);
      for (size_t k = 0; k < 3; k++)
        {
          feli.EvaluateGrad (mir, x.Range(k*ndofi, (k+1)*ndofi), grad);
          for (size_t i = 0; i < mir.Size(); i++)
            {
              auto trafo = Trans(mir[i].GetJacobianInverse());
              Vec<3,SIMD<double>> gi = grad.Col(i);
              Vec<3,SIMD<double>> tek = trafo.Col(k);
              Vec<3,SIMD<double>> hv = Cross(gi, tek);
              y.Col(i).Range(0,3) += hv;
            }
        }             
    }    
      
    using DiffOp<DiffOpCurlVectorL2Covariant>::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<3,3>&> (bmir);
      auto & fel = static_cast<const CompoundFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      size_t ndofi = feli.GetNDof();

      STACK_ARRAY(SIMD<double>, mem, 3*mir.Size());
      FlatMatrix<SIMD<double>> grad(3, mir.Size(), &mem[0]);
      for (size_t k = 0; k < 3; k++)
        {
          for (size_t i = 0; i < mir.Size(); i++)
            {
              auto trafo = Trans(mir[i].GetJacobianInverse());
              Vec<3,SIMD<double>> cy = y.Col(i);
              Vec<3,SIMD<double>> tek = trafo.Col(k);
              Vec<3,SIMD<double>> hv = Cross(cy, tek);
              grad.Col(i) = -hv;
            }
          
          feli.AddGradTrans (mir, grad, x.Range(k*ndofi, (k+1)*ndofi));
        }             
    }    
  };
  
  /*class DiffOpCurlVectorL2Covariant : public DiffOp<DiffOpCurlVectorL2Covariant>
  {
  public:
    enum { DIM = 2 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      auto & bfel = static_cast<const CompoundFiniteElement&> (fel);
      mat = 0.0;
      auto & feli = static_cast<const BaseScalarFiniteElement&> (bfel[0]);
      for (int i = 0; i < 2; i++)
        feli.CalcMappedDShape(mip.IP(), Trans(mat.Row(i).Range(fel.GetRange(i))));
      //~ int nd = fel.GetNDof();
    
      //~ mat = 0;
      for (int i = 0; i < nd; i++)
      {
	    mat(0, DIM*i  ) =  grad(1, i);
	    mat(0, DIM*i+1) = -grad(0, i);
	  }
    }
  };*/


  DocInfo VectorL2FESpace :: GetDocu ()
  {
    DocInfo docu = FESpace::GetDocu(); 
    docu.short_docu = "A vector-valued L2-conforming finite element space.";
    docu.long_docu =
      R"raw_string(The Vector-L2 finite element space is a product-space of L2 spaces,
where the number of components coincides with the mesh dimension.

It is implemented by means of a CompoundFESpace, as one could do it at the
user-level. Additionally, some operators are added for convenience and performance:
One can evaluate the vector-valued function, and one can take the gradient.
)raw_string";
    docu.Arg("piola") = "bool = False\n"
      "  Use Piola transform to map to physical element\n"
      "  allows to use the div-differential operator.";
    docu.Arg("covariant") = "bool = False\n"
      "  Use the covariant transform to map to physical element\n"
      "  allows to use the curl-differential operator.";
    docu.Arg("all_dofs_together") = "bool = True\n"
      "  dofs within one scalar component are together.";
    docu.Arg("hide_all_dofs") = "bool = False\n"
      "  all dofs are condensed without a global dofnr";
    
    return docu;
  }

  VectorL2FESpace::VectorL2FESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
                     bool checkflags)
      : CompoundFESpace(ama, flags)
    {
      type = "VectorL2";
      Flags compflags = flags;
      if (!flags.GetDefineFlagX("all_dofs_together").IsFalse())
        compflags.SetFlag("all_dofs_together");
      for (int i = 0; i <  ma->GetDimension(); i++)
        AddSpace (make_shared<L2HighOrderFESpace> (ama, compflags));

      piola = flags.GetDefineFlag ("piola");
      covariant = flags.GetDefineFlag ("covariant");

      if (piola)
        {
          switch (ma->GetDimension())
            {
            case 2:
              evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdVectorL2Piola<2>>>();
              flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivVectorL2Piola<2>>> ();
              additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradVectorL2Piola<2>>> ());
              break;
            case 3:
              evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdVectorL2Piola<3>>>();
              flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivVectorL2Piola<3>>> ();
              additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradVectorL2Piola<3>>> ());
              break;
            }
        }
      else if (covariant)
        {
          switch (ma->GetDimension())
            {
            case 2:
              evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdVectorL2Covariant<2>>>();
              break;
            case 3:
              evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdVectorL2Covariant<3>>>();
              flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlVectorL2Covariant>>();
              break;
            }
        }
      else
        switch (ma->GetDimension())
          {
          case 2:
            evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdVectorH1<2>>>();
            flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradVectorH1<2>>>();
            additional_evaluators.Set ("div", make_shared<T_DifferentialOperator<DiffOpDivVectorH1<2>>> ());
            additional_evaluators.Set ("Grad", make_shared<T_DifferentialOperator<DiffOpGradVectorH1<2>>> ());
	    additional_evaluators.Set ("hesse", make_shared<VectorDifferentialOperator>(make_shared<T_DifferentialOperator<DiffOpHesse<2>>>(), 2));
            break;
          case 3:
            evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdVectorH1<3>>>();
            flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradVectorH1<3>>>();
            additional_evaluators.Set ("div", make_shared<T_DifferentialOperator<DiffOpDivVectorH1<3>>> ());
            additional_evaluators.Set ("Grad", make_shared<T_DifferentialOperator<DiffOpGradVectorH1<3>>> ());
            additional_evaluators.Set ("hesse", make_shared<VectorDifferentialOperator>(make_shared<T_DifferentialOperator<DiffOpHesse<3>>>(), 3));
            break;
          }
    }


  void VectorL2FESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if (!DefinedOn (ei) || ei.VB() != VOL) return;
    auto & fes1 = static_cast<const L2HighOrderFESpace&> (*spaces[0]);
    auto eldofs = fes1.GetElementDofs(ei.Nr());
    size_t ndof1 = fes1.GetNDof();
    size_t locndof1 = eldofs.Size();
    dnums.SetSize (locndof1 * spaces.Size());
    for (size_t j = 0; j < spaces.Size(); j++)
      dnums.Range(j*locndof1, (j+1)*locndof1) = eldofs + j*ndof1;
    /*
    size_t size = eldofs.Size();
    size_t base = all_dofs_together ? 0 : 1;
    size += base;
    dnums.SetSize(size);
    if (!all_dofs_together) dnums[0] = ei.Nr();
    dnums.Range(base, size) = eldofs;
    */
  }


  shared_ptr<BaseMatrix> VectorL2FESpace ::
  GetMassOperator (shared_ptr<CoefficientFunction> rho,
                   shared_ptr<Region> defon,
                   LocalHeap & lh) const
  {
    /*
    cout << "VectorL2, GetMassOp" << endl
         << "rho = " << *rho << endl
         << "elwiseconst = " << rho->ElementwiseConstant() << endl;
    */
    if (rho->ElementwiseConstant() && order_policy == CONSTANT_ORDER)
      // && (covariant || piola || (rho && rho->Dimension() > 1) )
      {
        // cout << "optimized vector-apply-mass" << endl;
        switch (ma->GetDimension())
          {
          case 1:
            return make_shared<ApplyMassVectorL2Const<1>>
              (dynamic_pointer_cast<FESpace>(const_cast<VectorL2FESpace*>(this)->shared_from_this()),
               rho, defon, lh);    
          case 2:
            return make_shared<ApplyMassVectorL2Const<2>>
              (dynamic_pointer_cast<FESpace>(const_cast<VectorL2FESpace*>(this)->shared_from_this()),
               rho, defon, lh);    
          case 3:
            return make_shared<ApplyMassVectorL2Const<3>>
              (dynamic_pointer_cast<FESpace>(const_cast<VectorL2FESpace*>(this)->shared_from_this()),
               rho, defon, lh);    
          }
      }
    return FESpace::GetMassOperator(rho, defon, lh);
  }

  
  /*
  shared_ptr<BaseMatrix> VectorL2FESpace ::
  GetMassOperator (shared_ptr<CoefficientFunction> rho,
                   shared_ptr<Region> defon,
                   LocalHeap & lh) const
  {
    return FESpace::GetMassOperator(rho, defon, lh);
  }
  */
  
  
  void VectorL2FESpace :: SolveM (CoefficientFunction * rho, BaseVector & vec, Region * def,
                                  LocalHeap & lh) const
  {
    if (covariant || piola || (rho && rho->Dimension() > 1) )
      {
        switch (ma->GetDimension())
          {
          case 1: SolveM_Dim<1>(rho, vec, def, lh); break;
          case 2: SolveM_Dim<2>(rho, vec, def, lh); break;
          case 3: SolveM_Dim<3>(rho, vec, def, lh); break;
          default: throw Exception("VectorL2FESpace::SolveM: illegal dimension");
          }
        return;
      }
    
    for (size_t i = 0; i < spaces.Size(); i++)
      {
        auto veci = vec.Range (GetRange(i));
        spaces[i] -> SolveM (rho, veci, def, lh);
      }
  }


  
  
  void VectorL2FESpace :: ApplyM (CoefficientFunction * rho, BaseVector & vec,
                                  Region * defon,
                                  LocalHeap & lh) const
  {
    if (piola)
      {
        switch (ma->GetDimension())
          {
          case 1: ApplyMPiola<1>(rho, vec, defon, lh); break;
          case 2: ApplyMPiola<2>(rho, vec, defon, lh); break;
          case 3: ApplyMPiola<3>(rho, vec, defon, lh); break;
          default: throw Exception("VectorL2FESpace::ApplyM: illegal dimension");
          }
        return;
      }

    if (covariant)
      {
        switch (ma->GetDimension())
          {
          case 1: ApplyMCovariant<1>(rho, vec, defon, lh); break;
          case 2: ApplyMCovariant<2>(rho, vec, defon, lh); break;
          case 3: ApplyMCovariant<3>(rho, vec, defon, lh); break;
          default: throw Exception("VectorL2FESpace::ApplyM: illegal dimension");
          }
        return;
      }

    if (rho && rho->Dimension() > 1)
      {
        switch (ma->GetDimension())
          {
          case 1: ApplyM_Dim<1>(rho, vec, defon, lh); break;
          case 2: ApplyM_Dim<2>(rho, vec, defon, lh); break;
          case 3: ApplyM_Dim<3>(rho, vec, defon, lh); break;
          default: throw Exception("VectorL2FESpace::ApplyM: illegal dimension");
          }
        return;
      }
    for (size_t i = 0; i < spaces.Size(); i++)
      {
        auto veci = vec.Range (GetRange(i));
        spaces[i] -> ApplyM (rho, veci, defon, lh);
      }
  }
    
    
  
  template <int DIM>
  void VectorL2FESpace ::
  SolveM_Dim (CoefficientFunction * rho, BaseVector & vec, Region * def,
              LocalHeap & lh) const
  {
    static Timer t("SolveM - Vec"); RegionTimer reg(t);
    
    IterateElements
      (*this, VOL, lh,
       [&rho, &vec, def, this] (FESpace::Element el, LocalHeap & lh)
       {
         auto & fel = static_cast<const CompoundFiniteElement&>(el.GetFE());
         auto & feli = static_cast<const BaseScalarFiniteElement&>(fel[0]);
         const ElementTransformation & trafo = el.GetTrafo();
         
         Array<int> dnums(fel.GetNDof(), lh);
         GetDofNrs (el.Nr(), dnums);
         
         FlatVector<double> elx(feli.GetNDof()*DIM, lh);

         if (def && !def->Mask()[ma->GetElIndex(el)])
           {
             elx = 0.0;
             vec.SetIndirect (dnums, elx);
             return;
           }

         vec.GetIndirect(dnums, elx);
         auto melx = elx.AsMatrix(DIM, feli.GetNDof());
         
         FlatVector<double> diag_mass(feli.GetNDof(), lh);
         feli.GetDiagMassMatrix (diag_mass);
         
         bool curved = trafo.IsCurvedElement();
         if (rho && !rho->ElementwiseConstant()) curved = true;
         curved = false;  // curved not implemented
         
         if (!curved)
           {
             IntegrationRule ir(fel.ElementType(), 0);
             MappedIntegrationRule<DIM,DIM> mir(ir, trafo, lh);

             Mat<DIM,DIM> rhoi(0.0);
             if (!rho)
               rhoi = Identity(DIM);
             else if (rho->Dimension() == 1)
               rhoi = rho->Evaluate(mir[0]) * Identity(DIM);
             else
               rho -> Evaluate(mir[0], FlatVector<> (DIM*DIM, &rhoi(0,0)));
             
             Mat<DIM> trans(0.0);
             if (piola)
               trans = (1/mir[0].GetMeasure()) * Trans(mir[0].GetJacobian()) * rhoi * mir[0].GetJacobian();
             else if (covariant)
               trans = mir[0].GetMeasure() * mir[0].GetJacobianInverse() * rhoi * Trans(mir[0].GetJacobianInverse());               
             else
               trans = mir[0].GetMeasure() * rhoi;

             Mat<DIM> invtrans = Inv(trans);
             
             for (int i = 0; i < melx.Width(); i++)
               {
                 Vec<DIM> hv = melx.Col(i);
                 hv /=  diag_mass(i);
                 melx.Col(i) = invtrans * hv;
               }
           }
         /*
           else
           {
           SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
           auto & mir = trafo(ir, lh);
           FlatVector<SIMD<double>> pntvals(ir.Size(), lh);
           FlatMatrix<SIMD<double>> rhovals(1, ir.Size(), lh);
           if (rho) rho->Evaluate (mir, rhovals);
                 
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           for (int comp = 0; comp < dimension; comp++)
           {
           fel.Evaluate (ir, melx.Col(comp), pntvals);
           if (rho)
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / (mir[i].GetMeasure() * rhovals(0,i));
           else
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / mir[i].GetMeasure();
           
           melx.Col(comp) = 0.0;
           fel.AddTrans (ir, pntvals, melx.Col(comp));
           }
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           }
         */
         vec.SetIndirect(dnums, elx);
       });
    }

    
#ifdef OLD
  template <int DIM>
  void VectorL2FESpace ::
  SolveMPiola (CoefficientFunction * rho, BaseVector & vec,
               LocalHeap & lh) const
  {
    static Timer t("SolveM - Piola"); RegionTimer reg(t);
        
    IterateElements
      (*this, VOL, lh,
       [&rho, &vec,this] (FESpace::Element el, LocalHeap & lh)
       {
         auto & fel = static_cast<const CompoundFiniteElement&>(el.GetFE());
         auto & feli = static_cast<const BaseScalarFiniteElement&>(fel[0]);
         const ElementTransformation & trafo = el.GetTrafo();
         
         Array<int> dnums(fel.GetNDof(), lh);
         GetDofNrs (el.Nr(), dnums);
         
         FlatVector<double> elx(feli.GetNDof()*DIM, lh);
         vec.GetIndirect(dnums, elx);
         auto melx = elx.AsMatrix(DIM, feli.GetNDof());
         
         FlatVector<double> diag_mass(feli.GetNDof(), lh);
         feli.GetDiagMassMatrix (diag_mass);
         
         bool curved = trafo.IsCurvedElement();
         if (rho && !rho->ElementwiseConstant()) curved = true;
         curved = false;  // curved not implemented
         
         if (!curved)
           {
             IntegrationRule ir(fel.ElementType(), 0);
             BaseMappedIntegrationRule & mir = trafo(ir, lh);
             Mat<DIM> trans = (1/mir[0].GetMeasure()) * Trans(mir[0].GetJacobian()) * mir[0].GetJacobian();
             Mat<DIM> invtrans = Inv(trans);
             
             // double jac = mir[0].GetMeasure();
             // if (rho) jac *= rho->Evaluate(mir[0]);
             // diag_mass *= jac;
             
             for (int i = 0; i < melx.Width(); i++)
               {
                 Vec<DIM> hv = melx.Col(i);
                 hv /=  diag_mass(i);
                 melx.Col(i) = invtrans * hv;
               }
           }
         /*
           else
           {
           SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
           auto & mir = trafo(ir, lh);
           FlatVector<SIMD<double>> pntvals(ir.Size(), lh);
           FlatMatrix<SIMD<double>> rhovals(1, ir.Size(), lh);
           if (rho) rho->Evaluate (mir, rhovals);
                 
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           for (int comp = 0; comp < dimension; comp++)
           {
           fel.Evaluate (ir, melx.Col(comp), pntvals);
           if (rho)
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / (mir[i].GetMeasure() * rhovals(0,i));
           else
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / mir[i].GetMeasure();
           
           melx.Col(comp) = 0.0;
           fel.AddTrans (ir, pntvals, melx.Col(comp));
           }
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           }
         */
         vec.SetIndirect(dnums, elx);
       });
  }


  template <int DIM>
  void VectorL2FESpace ::
  SolveMCovariant (CoefficientFunction * rho, BaseVector & vec,
                   LocalHeap & lh) const
  {
    static Timer t("SolveM - Covariant"); RegionTimer reg(t);
    IterateElements
      (*this, VOL, lh,
       [&rho, &vec,this] (FESpace::Element el, LocalHeap & lh)
       {
         auto & fel = static_cast<const CompoundFiniteElement&>(el.GetFE());
         auto & feli = static_cast<const BaseScalarFiniteElement&>(fel[0]);
         const ElementTransformation & trafo = el.GetTrafo();
         
         Array<int> dnums(fel.GetNDof(), lh);
         GetDofNrs (el.Nr(), dnums);
         
         FlatVector<double> elx(feli.GetNDof()*DIM, lh);
         vec.GetIndirect(dnums, elx);
         auto melx = elx.AsMatrix(DIM, feli.GetNDof());
         
         FlatVector<double> diag_mass(feli.GetNDof(), lh);
         feli.GetDiagMassMatrix (diag_mass);
         
         bool curved = trafo.IsCurvedElement();
         if (rho && !rho->ElementwiseConstant()) curved = true;
         curved = false;  // curved not implemented
         
         if (!curved)
           {
             IntegrationRule ir(fel.ElementType(), 0);
             MappedIntegrationRule<DIM,DIM> mir(ir, trafo, lh);
             
             Mat<DIM,DIM> rhoi(0.0);
             if (!rho)
               rhoi = Identity(DIM);
             else if (rho->Dimension() == 1)
               rhoi = rho->Evaluate(mir[0]) * Identity(DIM);
             else
               rho -> Evaluate(mir[0], FlatVector<> (DIM*DIM, &rhoi(0,0)));
             
             // Mat<DIM> trans = (1/mir[0].GetMeasure()) * Trans(mir[0].GetJacobian()) * mir[0].GetJacobian();
             Mat<DIM> trans = mir[0].GetMeasure() * mir[0].GetJacobianInverse() * rhoi * Trans(mir[0].GetJacobianInverse());
             Mat<DIM> invtrans = Inv(trans);
             
             // double jac = mir[0].GetMeasure();
             // if (rho) jac *= rho->Evaluate(mir[0]);
             // diag_mass *= jac;
             
             for (int i = 0; i < melx.Width(); i++)
               {
                 Vec<DIM> hv = melx.Col(i);
                 hv /=  diag_mass(i);
                 melx.Col(i) = invtrans * hv;
               }
           }
         /*
           else
           {
           SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
           auto & mir = trafo(ir, lh);
           FlatVector<SIMD<double>> pntvals(ir.Size(), lh);
           FlatMatrix<SIMD<double>> rhovals(1, ir.Size(), lh);
           if (rho) rho->Evaluate (mir, rhovals);
                 
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           for (int comp = 0; comp < dimension; comp++)
           {
           fel.Evaluate (ir, melx.Col(comp), pntvals);
           if (rho)
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / (mir[i].GetMeasure() * rhovals(0,i));
           else
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / mir[i].GetMeasure();
           
           melx.Col(comp) = 0.0;
           fel.AddTrans (ir, pntvals, melx.Col(comp));
           }
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           }
         */
         vec.SetIndirect(dnums, elx);
       });
  }
#endif
  
  template <int DIM>
  void VectorL2FESpace ::
  ApplyM_Dim (CoefficientFunction * rho, BaseVector & vec, Region * def,
              LocalHeap & lh) const
  {
    static Timer t("ApplyM - Piola"); RegionTimer reg(t);
    IterateElements
      (*this, VOL, lh,
       [&rho, &vec, def,this] (FESpace::Element el, LocalHeap & lh)
       {
         auto & fel = static_cast<const CompoundFiniteElement&>(el.GetFE());
         auto & feli = static_cast<const BaseScalarFiniteElement&>(fel[0]);
         const ElementTransformation & trafo = el.GetTrafo();
         
         Array<int> dnums(fel.GetNDof(), lh);
         GetDofNrs (el.Nr(), dnums);

         
         FlatVector<double> elx(feli.GetNDof()*DIM, lh);
         vec.GetIndirect(dnums, elx);
         auto melx = elx.AsMatrix(DIM, feli.GetNDof());
         
         FlatVector<double> diag_mass(feli.GetNDof(), lh);
         feli.GetDiagMassMatrix (diag_mass);
         
         bool curved = trafo.IsCurvedElement();
         if (rho && !rho->ElementwiseConstant()) curved = true;
         curved = false;  // curved not implemented
         
         if (!curved)
           {
             IntegrationRule ir(fel.ElementType(), 0);
             MappedIntegrationRule<DIM,DIM> mir(ir, trafo, lh);
             
             Mat<DIM,DIM> rhoi(0.0);
             if (!rho)
               rhoi = Identity(DIM);
             else if (rho->Dimension() == 1)
               rhoi = rho->Evaluate(mir[0]) * Identity(DIM);
             else
               rho -> Evaluate(mir[0], FlatVector<> (DIM*DIM, &rhoi(0,0)));
             
             Mat<DIM> trans = mir[0].GetMeasure() * rhoi;
             
             for (int i = 0; i < melx.Width(); i++)
               {
                 Vec<DIM> hv = melx.Col(i);
                 hv *=  diag_mass(i);
                 melx.Col(i) = trans * hv;
               }
           }
         /*
           else
           {
           SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
           auto & mir = trafo(ir, lh);
           FlatVector<SIMD<double>> pntvals(ir.Size(), lh);
           FlatMatrix<SIMD<double>> rhovals(1, ir.Size(), lh);
           if (rho) rho->Evaluate (mir, rhovals);
                 
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           for (int comp = 0; comp < dimension; comp++)
           {
           fel.Evaluate (ir, melx.Col(comp), pntvals);
           if (rho)
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / (mir[i].GetMeasure() * rhovals(0,i));
           else
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / mir[i].GetMeasure();
           
           melx.Col(comp) = 0.0;
           fel.AddTrans (ir, pntvals, melx.Col(comp));
           }
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           }
         */
         if (def && !def->Mask()[ma->GetElIndex(el)])
           elx = 0.0;
         vec.SetIndirect(dnums, elx);
       });
    
  }

  

  template <int DIM>
  void VectorL2FESpace ::
  ApplyMPiola (CoefficientFunction * rho, BaseVector & vec, Region * def,
               LocalHeap & lh) const
  {
    static Timer t("ApplyM - Piola"); RegionTimer reg(t);
    IterateElements
      (*this, VOL, lh,
       [&rho, &vec, def,this] (FESpace::Element el, LocalHeap & lh)
       {
         auto & fel = static_cast<const CompoundFiniteElement&>(el.GetFE());
         auto & feli = static_cast<const BaseScalarFiniteElement&>(fel[0]);
         const ElementTransformation & trafo = el.GetTrafo();
         
         Array<int> dnums(fel.GetNDof(), lh);
         GetDofNrs (el.Nr(), dnums);

         
         FlatVector<double> elx(feli.GetNDof()*DIM, lh);
         vec.GetIndirect(dnums, elx);
         auto melx = elx.AsMatrix(DIM, feli.GetNDof());
         
         FlatVector<double> diag_mass(feli.GetNDof(), lh);
         feli.GetDiagMassMatrix (diag_mass);
         
         bool curved = trafo.IsCurvedElement();
         if (rho && !rho->ElementwiseConstant()) curved = true;
         curved = false;  // curved not implemented
         
         if (!curved)
           {
             IntegrationRule ir(fel.ElementType(), 0);
             MappedIntegrationRule<DIM,DIM> mir(ir, trafo, lh);
             
             Mat<DIM,DIM> rhoi(0.0);
             if (!rho)
               rhoi = Identity(DIM);
             else if (rho->Dimension() == 1)
               rhoi = rho->Evaluate(mir[0]) * Identity(DIM);
             else
               rho -> Evaluate(mir[0], FlatVector<> (DIM*DIM, &rhoi(0,0)));
             
             Mat<DIM> trans = 1/mir[0].GetMeasure() * Trans(mir[0].GetJacobian()) * rhoi * mir[0].GetJacobian();
             
             for (int i = 0; i < melx.Width(); i++)
               {
                 Vec<DIM> hv = melx.Col(i);
                 hv *=  diag_mass(i);
                 melx.Col(i) = trans * hv;
               }
           }
         /*
           else
           {
           SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
           auto & mir = trafo(ir, lh);
           FlatVector<SIMD<double>> pntvals(ir.Size(), lh);
           FlatMatrix<SIMD<double>> rhovals(1, ir.Size(), lh);
           if (rho) rho->Evaluate (mir, rhovals);
                 
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           for (int comp = 0; comp < dimension; comp++)
           {
           fel.Evaluate (ir, melx.Col(comp), pntvals);
           if (rho)
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / (mir[i].GetMeasure() * rhovals(0,i));
           else
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / mir[i].GetMeasure();
           
           melx.Col(comp) = 0.0;
           fel.AddTrans (ir, pntvals, melx.Col(comp));
           }
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           }
         */
         if (def && !def->Mask()[ma->GetElIndex(el)])
           elx = 0.0;
         vec.SetIndirect(dnums, elx);
       });
    
  }


  template <int DIM>
  void VectorL2FESpace ::
  ApplyMCovariant (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                   LocalHeap & lh) const
  {
    static Timer t("ApplyM - Covariant"); RegionTimer reg(t);
    IterateElements
      (*this, VOL, lh,
       [&rho, &vec,this] (FESpace::Element el, LocalHeap & lh)
       {
         auto & fel = static_cast<const CompoundFiniteElement&>(el.GetFE());
         auto & feli = static_cast<const BaseScalarFiniteElement&>(fel[0]);
         const ElementTransformation & trafo = el.GetTrafo();
         
         Array<int> dnums(fel.GetNDof(), lh);
         GetDofNrs (el.Nr(), dnums);
         
         FlatVector<double> elx(feli.GetNDof()*DIM, lh);
         vec.GetIndirect(dnums, elx);
         auto melx = elx.AsMatrix(DIM, feli.GetNDof());
         
         FlatVector<double> diag_mass(feli.GetNDof(), lh);
         feli.GetDiagMassMatrix (diag_mass);
         
         bool curved = trafo.IsCurvedElement();
         if (rho && !rho->ElementwiseConstant()) curved = true;
         curved = false;  // curved not implemented
         
         if (!curved)
           {
             IntegrationRule ir(fel.ElementType(), 0);
             MappedIntegrationRule<DIM,DIM> mir(ir, trafo, lh);
             
             Mat<DIM,DIM> rhoi(0.0);
             if (!rho)
               rhoi = Identity(DIM);
             else if (rho->Dimension() == 1)
               rhoi = rho->Evaluate(mir[0]) * Identity(DIM);
             else
               rho -> Evaluate(mir[0], FlatVector<> (DIM*DIM, &rhoi(0,0)));
             
             Mat<DIM> trans = mir[0].GetMeasure() * mir[0].GetJacobianInverse() * rhoi * Trans(mir[0].GetJacobianInverse());
             
             for (int i = 0; i < melx.Width(); i++)
               {
                 Vec<DIM> hv = melx.Col(i);
                 hv *=  diag_mass(i);
                 melx.Col(i) = trans * hv;
               }
           }
         /*
           else
           {
           SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
           auto & mir = trafo(ir, lh);
           FlatVector<SIMD<double>> pntvals(ir.Size(), lh);
           FlatMatrix<SIMD<double>> rhovals(1, ir.Size(), lh);
           if (rho) rho->Evaluate (mir, rhovals);
                 
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           for (int comp = 0; comp < dimension; comp++)
           {
           fel.Evaluate (ir, melx.Col(comp), pntvals);
           if (rho)
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / (mir[i].GetMeasure() * rhovals(0,i));
           else
           for (size_t i = 0; i < ir.Size(); i++)
           pntvals(i) *= ir[i].Weight() / mir[i].GetMeasure();
           
           melx.Col(comp) = 0.0;
           fel.AddTrans (ir, pntvals, melx.Col(comp));
           }
           for (int i = 0; i < melx.Height(); i++)
           melx.Row(i) /= diag_mass(i);
           }
         */
         vec.SetIndirect(dnums, elx);
       });

  }

  
  static RegisterFESpace<VectorL2FESpace> initvecl2 ("VectorL2");
  


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
      GetFESpaceClasses().AddFESpace ("l2", L2HighOrderFESpace::Create, L2HighOrderFESpace::GetDocu);
      GetFESpaceClasses().AddFESpace ("l2ho", L2HighOrderFESpace::CreateHO, L2HighOrderFESpace::GetDocu);
      GetFESpaceClasses().AddFESpace ("l2surf", L2SurfaceHighOrderFESpace::Create, L2SurfaceHighOrderFESpace::GetDocu);
    }
    
    Init init;
  }
}
