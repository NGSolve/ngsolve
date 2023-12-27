#include <comp.hpp>
#include <../fem/hdiv_equations.hpp>  
#include <../fem/hdivlofe.hpp>  
#include <../fem/hdivhofe.hpp>  
#include <../fem/hdivhofefo.hpp>  
#include <../fem/hcurlhdiv_dshape.hpp> 

namespace ngcomp
{


template <int D, typename FEL = HDivNormalFiniteElement<D-2> >
class DiffOpIdVecHDivSurfaceBoundary : public DiffOp<DiffOpIdVecHDivSurfaceBoundary<D, FEL> >
{
public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-2 };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 0 };

    static const FEL & Cast(const FiniteElement & fel)
    {
        return static_cast<const FEL&> (fel);
    }

    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {

      throw Exception("Does not work yet!!!! eltrans 1D -> 3D, no normal vector available");
      // auto normal = Vec<3>(mip.GetNV());      
      // auto tangential = mip.GetTV();      
      // Vec<3> normalel = Cross(normal, tangential);
      getchar();
	
      auto scaled_nv = (1.0/mip.GetJacobiDet()) * mip.GetNV();
      mat = scaled_nv * Trans(Cast(fel).GetShape (mip.IP(), lh));
      //throw Exception("in DiffOpIdVecHDivSurfaceBoundary::GenerateMatrix");      
    }

};

  
  

template <int D, typename FEL = HDivFiniteElement<D-1> >
class DiffOpDivHDivSurface : public DiffOp<DiffOpDivHDivSurface<D, FEL> >
{
public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    static const FEL & Cast(const FiniteElement & fel)
    {
        return static_cast<const FEL&> (fel);
    }

    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      mat = (1.0 / mip.GetJacobiDet()) * 
	Trans (static_cast<const FEL&>(fel).GetDivShape(mip.IP(),lh));
    }

  static string Name() { return "div"; }

  static shared_ptr<CoefficientFunction>
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian)
  {
    if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpDivHDivSurface");    
    return -TraceCF(dir->Operator("Gradboundary"))*proxy;     
  }

};

  
  
  //Gradient operator of HDivSurface
  template <int D, typename FEL = HDivFiniteElement<D-1> >
  class DiffOpGradientHDivSurface : public DiffOp<DiffOpGradientHDivSurface<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 1 };

    static Array<int> GetDimensions() { return Array<int> ( { D, D } ); };
    static constexpr double eps() { return 1e-4; }

    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
                                                  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                                                              MAT mat, LocalHeap & lh)
    {      
      CalcDShapeFE<FEL,D,D-1,D>(static_cast<const FEL&>(fel), mip, Trans(mat), lh, eps());
    }

    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL & fel, const MIP & mip,
                       const TVX & x, TVY && y,
                       LocalHeap & lh) 
    {
      // typedef typename TVX::TSCAL TSCAL;
      HeapReset hr(lh);
      FlatMatrixFixWidth<D*D> hm(fel.GetNDof(),lh);
      CalcDShapeFE<FEL,D,D-1,D>(static_cast<const FEL&>(fel), mip, hm, lh, eps());
      y = Trans(hm)*x;
    }

    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void ApplyTrans (const AFEL & fel, const MIP & mip,
			    const TVX & x, TVY & by,
			    LocalHeap & lh) 
    {
      ApplyTransDShapeFE<FEL,D,D-1,D>(static_cast<const FEL&>(fel), mip, x, by, lh, eps());
    }

    
    /*static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> mat)
    {
      CalcSIMDDShapeFE<FEL,D,D-1,D>(static_cast<const FEL&>(fel), static_cast<const SIMD_MappedIntegrationRule<D-1,D> &>(bmir), mat, eps());
      }*/

    
    using DiffOp<DiffOpGradientHDivSurface<D>>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      ApplySIMDDShapeFE<FEL,D,D-1,D>(static_cast<const FEL&>(fel), bmir, x, y, eps());
    }


    using DiffOp<DiffOpGradientHDivSurface<D>>::AddTransSIMDIR;    
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
    {
      AddTransSIMDDShapeFE<FEL,D,D-1,D>(static_cast<const FEL&>(fel), bmir, x, y, eps());
    }
  };

  HDivHighOrderSurfaceFESpace ::  
  HDivHighOrderSurfaceFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    type = "hdivhosurface";
    name="HDivHighOrderSurfaceFESpace(hdivhosurf)";
       
    DefineDefineFlag("discontinuous");   
    DefineDefineFlag("hodivfree");
    DefineNumFlag("orderinner");
    DefineDefineFlag("RT");

    
    if(parseflags) CheckFlags(flags);

    discont = flags.GetDefineFlag("discontinuous"); 

    order =  int (flags.GetNumFlag ("order",0));       

    if (flags.NumFlagDefined("order")) 
      {
	order =  int (flags.GetNumFlag ("order",0));
      }
    else 
      {       
	order = 0;  
      }

    uniform_order_inner = int (flags.GetNumFlag ("orderinner", -1));
    uniform_order_facet = int (flags.GetNumFlag ("orderfacet", -1));

    *testout << "uniform_order_inner = " << uniform_order_inner << endl;

    ho_div_free = flags.GetDefineFlag("hodivfree");
    RT = flags.GetDefineFlag ("RT");

           
    auto one = make_shared<ConstantCoefficientFunction> (1);
    
    if (ma->GetDimension() < 3)
      {
	throw Exception ("only 2D manifolds supported");       
      }
    else
      {
	evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<3>>>();	
	evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdHDivSurface<3>>>();
	evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdVecHDivSurfaceBoundary<3>>>();

	flux_evaluator[VOL] =  make_shared<T_DifferentialOperator<DiffOpDivHDiv<3>>>();
	flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpDivHDivSurface<3>>>();

        additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHDivSurface<3>>> ());
        additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHDivDualSurface<3>>> ());
      }
    
    highest_order_dc = flags.GetDefineFlag("highest_order_dc");
    if (highest_order_dc) {
      *testout << "highest_order_dc is active!" << endl;
    }
  }
  
  HDivHighOrderSurfaceFESpace:: ~HDivHighOrderSurfaceFESpace () 
  {
    ;
  }


  DocInfo HDivHighOrderSurfaceFESpace :: GetDocu ()
  {
    auto docu = FESpace::GetDocu();
    docu.Arg("discontinuous") = "bool = False\n"
      "  Create discontinuous HDivSurface space";
    docu.Arg("hodivfree") = "bool = False\n"
      "  Remove high order element bubbles with non zero divergence";
    docu.Arg("RT") = "bool = False\n"
      "  RT elements for simplicial elements: P^k subset RT_k subset P^{k+1}";
    docu.Arg("orderinner") = "optional<int> = undefined\n"
      "  order for inner space if defined, otherwise use order for orderinner";
    return docu;
  }
  
void HDivHighOrderSurfaceFESpace :: Average (BaseVector & vec) const
{
  // auto & pairs = GetDCPairs();
  auto fu = vec.FV<double>();
  for (auto pair : dc_pairs)
      {
        auto f1 = pair[0];
        auto f2 = pair[1];
        if (f2 != -1)
          {
            double mean = 0.5 * (fu(f1) + fu(f2));
            fu(f1) = fu(f2) = mean;
          }
        else if (f1 != -1)
          fu(f1) = 0.0;
      }
  }
  void HDivHighOrderSurfaceFESpace :: Update()
  {
    FESpace::Update();

    size_t nel = ma->GetNSE();
    size_t nfa = ma->GetNEdges();
    // size_t dim = ma->GetDimension();
       
    first_facet_dof.SetSize(nfa+1);
    first_inner_dof.SetSize(nel+1);

    order_facet.SetSize(nfa);
    order_inner.SetSize(nel);

    if(uniform_order_inner > -1)
      order_inner = uniform_order_inner;
    else
      order_inner = order;
    if(uniform_order_facet > -1)
      order_facet = uniform_order_facet;
    else
      order_facet = order;
    
    fine_facet.SetSize(nfa); 
    fine_facet = 0; //!!!! 


    for (auto el : ma->Elements(BND))
      {
        if (!DefinedOn (el))
          {
            order_inner[el.Nr()] = 0;
            // order_inner_curl[el.Nr()] = 0;
            continue;
          }
          
	// ELEMENT_TYPE eltype = el.GetType();
	// const POINT3D * points = ElementTopology :: GetVertices (eltype);
        auto elfacets = el.Edges();

        fine_facet[elfacets] = true;
	
	// if(!var_order) continue; 
	
	// INT<3> el_orders = ma->GetElOrders(el.Nr());  
	
        // int i = el.Nr();
	// for(int k=0;k<dim;k++)
	//   {
	//     // order_inner_curl[i][k]= max2(el_orders[k] + rel_curl_order,0);
	//     order_inner[i][k] = max2(el_orders[k]+rel_order,0);
	//   }

	// if(dim==2)
	//   {
	//     const EDGE * edges = ElementTopology::GetEdges (eltype);
	//     for(int j=0; j<elfacets.Size(); j++)
	//       for(int k=0;k<2;k++)
	// 	if(points[edges[j][0]][k] != points[edges[j][1]][k])
	// 	  { 
	// 	    order_facet[elfacets[j]][0] = max2(el_orders[k]+rel_curl_order, order_facet[elfacets[j]][0]);
	// 	    break; 
	// 	  }
	//   }
      }

    
    if(print) 
      {
	*testout << " discont " << discont << endl;
	*testout << " fine_facet[i] (hdiv) " << fine_facet << endl; 
	
	*testout << " order_facet (hdivho) " << order_facet << endl; 
	*testout << " order_inner (hdivho) " << order_inner << endl; 	
      }

    // for (int i = 0; i < nel; i++)
    //   {
    //     order_inner[i] = INT<3> (order,order,order);
    //   }

    // if(uniform_order_inner>-1) 
    //   order_inner = INT<3> (uniform_order_inner,uniform_order_inner,uniform_order_inner);

    UpdateDofTables(); 
    UpdateCouplingDofArray();
  }

  void HDivHighOrderSurfaceFESpace :: UpdateDofTables()
  {
    size_t nel = ma->GetNSE();
    size_t nfa = ma->GetNEdges();
    size_t dim = ma->GetDimension();
    
    ndof = nfa;
    first_facet_dof = ndof;

    if(dim==3)
      {
	for (auto i : Range(nfa))
          {
            first_facet_dof[i] = ndof;
	    // ndof += order;
            int inc = fine_facet[i] ? order_facet[i][0] : 0;
	    //if (highest_order_dc && !boundary_facet[i]) inc--;
            if (inc > 0) ndof += inc;            
          }

        first_facet_dof[nfa] = ndof;
      
        for (size_t i = 0; i < nel; i++)
          {
            ElementId ei(BND, i);
            //INT<3> pc = order_inner_curl[i];
            INT<3> p = order_inner[i];
            int inci = 0;
            switch(ma->GetElType(ei))
              {
              case ET_TRIG:
                if (!ho_div_free)
		  inci = p[0]*(p[0]-1)/2+p[0]*(p[0]-1)/2+p[0]-1;
		//inci = pc[0]*(pc[0]-1)/2 + p[0]*(p[0]-1)/2 + p[0]-1;
                else
		  inci = p[0]*(p[0]-1)/2;
		    //inci = pc[0]*(pc[0]-1)/2;
                if (RT)
                  inci += p[0] + 1;
                break;
              case ET_QUAD:
		if (!ho_div_free)
                  inci = p[0]*p[1] + p[0]*p[1]+p[0]+p[1];
                else
                  inci = p[0]*p[1];
                break;
              default: // for the compiler
                break;  
              }

	    if (highest_order_dc)
	      {
                /*
		ma->GetElFacets (ei, fnums);
                for (auto f : fnums)
		  if (!boundary_facet[f]) inci++;
                */
                for (auto f : ma->GetElFacets(ei))
		  if (!boundary_facet[f]) inci++;
	      }

            first_inner_dof[i] = ndof;
            if (inci > 0) ndof+=inci;
          }
        first_inner_dof[nel] = ndof;


        if (highest_order_dc)
          {
            dc_pairs.SetSize (ma->GetNFacets());
            dc_pairs = INT<2> (-1,-1);
            
            // Array<int> fnums;
            for (auto ei : ma->Elements(VOL))
              {
                // auto i = ei.Nr();
                // ma->GetElFacets (ei, fnums);
                // auto fnums = ma->GetElFacets(ei);
		int fid = first_inner_dof[ei.Nr()];
                for (auto f : ma->GetElFacets(ei))
		  if (!boundary_facet[f])
		    {
		      int di = fid++; // first_inner_dof[i]+k;
		      dc_pairs[f][1] = dc_pairs[f][0];
		      dc_pairs[f][0] = di;
		    }
              }
          }
        else
          dc_pairs.SetSize0 ();

        if (discont)
          {
            ndof = 0;
            for(auto i : Range(nel))
              {
                ElementId ei(BND, i);
                int incii = first_inner_dof[i+1] - first_inner_dof[i];
                auto elfacets = ma->GetElEdges(ei);
                if(DefinedOn(ei))
                  for(const auto& facet : elfacets)
                    incii += first_facet_dof[facet+1] - first_facet_dof[facet] + 1;
                first_inner_dof[i] = ndof;
                ndof += incii;
              }
            first_inner_dof[nel] = ndof;
            first_facet_dof = 0;
          }
      }
    else
      {
	throw Exception ("you should not be here - only 2D manifolds supported");
      }    
        
    if(print) 
      {
        (*testout) << "ndof hdivhofespace update = " << endl << ndof << endl;        
        (*testout) << "first_facet_dof (hdiv) = " << endl << first_facet_dof << endl;
        (*testout) << "first_inner_dof (hdiv) = " << endl << first_inner_dof << endl;
      }
  }

  void HDivHighOrderSurfaceFESpace :: UpdateCouplingDofArray()
  {
    // throw Exception("This is copy paste! Not updated yet!!!");
    ctofdof.SetSize(ndof);
    if(discont) 
      {
        ctofdof = LOCAL_DOF;
        return;
      } 
    
    ctofdof = WIREBASKET_DOF;
    
    for (auto facet : Range(ma->GetNEdges()))
      {
        ctofdof[facet] = fine_facet[facet] ? WIREBASKET_DOF : UNUSED_DOF;
        ctofdof[GetFacetDofs(facet)] = INTERFACE_DOF;
      }
    
    for (auto el : Range (ma->GetNE(BND)))
      ctofdof[GetElementDofs(el)] = LOCAL_DOF;

    // cout << ctofdof << endl;
  }


  FiniteElement & HDivHighOrderSurfaceFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {

    if (!DefinedOn (ei))
      {
        return
          SwitchET (ma->GetElType(ei), [&] (auto et) -> FiniteElement&
                      {
                        return *new (alloc) DummyFE<et.ElementType()> ();
                      });
      }
    
     if (ei.IsVolume())
      {
	throw Exception("No volume elements available");
      }
     else if (ei.VB()  == BND)
       {
	 switch (ma->GetElType(ei))
	   {
	   case ET_TRIG: return T_GetSFE<ET_TRIG>(ei, false, alloc);
	   case ET_QUAD: return T_GetSFE<ET_QUAD>(ei, false, alloc);
      
	   default: throw Exception("illigal element in HDivHighOrderSurfaceFESpace::GetSFE");
	   }   
       }
     else if (ei.VB() == BBND)
       {
	 if(ma->GetElType(ei) == ET_SEGM)
          {
            if(!DefinedOn(ei) || discont) return *new (alloc) HDivNormalDummyFE<ET_SEGM>();
	    Ngs_Element ngel = ma->GetElement(ei);
	    //HDivHighOrderNormalFiniteElement<1> * hofe = new (alloc)HDivHighOrderNormalFiniteElement<1>();
	    auto * fe = new (alloc)HDivHighOrderNormalSegm<TrigExtensionMonomial>(order);
	    fe -> SetVertexNumbers(ngel.Vertices());

	    //I thnk this is not needed if we have a global edge order...
	    HDivHighOrderNormalFiniteElement<1> * hofe =
              dynamic_cast<HDivHighOrderNormalFiniteElement<1>*> (fe);
	    	    
	    hofe -> SetOrderInner(order);
	    hofe -> ComputeNDof();

	    return *hofe;
	  }
	 else
	   throw Exception("illegal element in HDivHighOrderSurfaceFESpace :: GetFE BBND");
       }
     else // BBBND
       return * new (alloc) DummyFE<ET_POINT>();
  }

  template<ELEMENT_TYPE ET>
    FiniteElement & HDivHighOrderSurfaceFESpace::T_GetSFE(ElementId ei, bool onlyhdiv, Allocator & lh) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    if(!DefinedOn(ngel)) return *new (lh) HDivDummyFE<ET>();
    HDivHighOrderFE<ET>* hofe = new (lh)HDivHighOrderFE<ET>();
    hofe->SetOrderInner(order_inner[ei.Nr()][0]);
    hofe->SetVertexNumbers(ngel.Vertices());
    hofe->SetRT(RT);


    Array<int> facet_order(ngel.Edges());
    facet_order = order;
    hofe->SetOrderFacet(facet_order);
    hofe->ComputeNDof();    

    return *hofe;
  }
  
  const FiniteElement & HDivHighOrderSurfaceFESpace :: GetHODivFE (int elnr, LocalHeap & lh) const
  {
    Ngs_Element ngel = ma->GetElement(ElementId(VOL,elnr));
    ELEMENT_TYPE eltype = ngel.GetType();
    
    if (!ho_div_free) throw Exception ("You don't have hodivfree active. You are not allow to call GetHODivFE");
    
    switch (eltype)
      {        
	//case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, true, lh);
	//case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, true, lh);       	
        
      default:
        throw Exception ("illegal element in HDivHOFeSpace::GetDivFE");
      }
  }
  
  size_t HDivHighOrderSurfaceFESpace :: GetNDof () const throw()
  {
    return ndof;
  }

  void HDivHighOrderSurfaceFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {    
    dnums.SetSize0();
    if(!DefinedOn(ei)) return;
    if(discont)
      {
        if(ei.VB() == BND)
          dnums += GetElementDofs(ei.Nr());
        return;
      }
    if(ei.VB()==VOL)
      {
	dnums.SetSize0();
      }
    else if(ei.VB()==BND)
      {
	auto fanums = ma->GetElEdges(ei);
	
	if (highest_order_dc)
	  {
	    IntRange eldofs = GetElementDofs (ei.Nr());
	    
	    dnums += fanums;
		
	    int first_el_dof = eldofs.First();
	    for (auto f : fanums)
	      {
		dnums += GetFacetDofs (f);
		if (!boundary_facet[f])
		  dnums += first_el_dof++;
	      }
	    dnums += IntRange (first_el_dof, eldofs.Next());	    
	  }	         
	else
	  {
	    //lowest order RT0 dofs
            dnums += fanums;
	    
	    // edges
	    for (auto f : fanums)
	      dnums += GetFacetDofs (f);   
	    	    
	    //inner
	    dnums += GetElementDofs (ei.Nr());
	  }
      }
    else if (ei.VB()==BBND)
      {
	auto fanums = ma->GetElEdges(ei);	
	dnums += fanums;	
	dnums += GetFacetDofs (fanums[0]);	
      }    
    
  }

  void HDivHighOrderSurfaceFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize0(); 
  }
  
  void HDivHighOrderSurfaceFESpace :: GetFacetDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 2 || discont) return; 

    dnums += fanr;
    //high order facet dofs
    dnums += GetFacetDofs (fanr);
  }

  void HDivHighOrderSurfaceFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums = GetElementDofs (elnr);
  }


  static RegisterFESpace<HDivHighOrderSurfaceFESpace> init ("hdivhosurface");
}



