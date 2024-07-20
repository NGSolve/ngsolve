/**
   High Order Finite Element Space for H(Div)
*/
/* Continuous and discontinuous version : 
   Disont space has no extra low-order block(!), all dofs internal! 
   Flags:
   -order  : uniform order 
   -relorder : variable order relative to mesh-order 
   (on facets maximum order of facet-patch elements, also for discont)
   -curl_order,relcurlorder: orders of curl fields in uniform/variable version 
   -disontinuous (DefineFlag) : discontinuous space 
   -orderinner, orderfacet: for variable and uniform space uniform order for inner or facets 
	 
*/ 



// #include <comp.hpp>
#include "hdivhofespace.hpp"
#include "../fem/hdiv_equations.hpp"
#include <../fem/hdivlofe.hpp>  
#include <../fem/hdivhofe.hpp>  
#include <../fem/hdivhofefo.hpp>  
#include <../fem/hcurlhdiv_dshape.hpp> 
#include <../fem/diffop_impl.hpp>
#include <multigrid.hpp> 
#include <fesconvert.hpp>

namespace ngcomp
{



  template<int D>
  class DiffOpNormalComponentHDiv: public DiffOp<DiffOpNormalComponentHDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };
    // enum { DIM_STRESS = D };

    // static Array<int> GetDimensions() { return Array<int> ({1}); }
    static Array<int> GetDimensions() { return Array<int>(0); }

    /*
      template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
    }
    */
    
    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT && mat,LocalHeap & lh)
    {
      HeapReset hr(lh);
      auto & fel = dynamic_cast<const HDivFiniteElement<D>&> (bfel);

      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,D,lh);
      
      Vec<D> n = sip.GetNV();
      fel.CalcMappedShape(sip,shape);
      mat.Row(0) = shape * n;
    }

    
    static int DimRef() { return 1; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      FlatMatrix<> tmp(fel.GetNDof(), D, lh);
      static_cast<const BaseHDivFiniteElement&> (fel).CalcShape (ip, tmp);
      Vec<D> nv = ElementTopology::GetNormals<D>(fel.ElementType())[ip.FacetNr()];
      mat.Row(0) = tmp * nv;
    }

  template <typename MIP, typename MAT>
  static void CalcTransformationMatrix (const MIP & bmip,
                                        MAT & mat, LocalHeap & lh)
  {
    mat(0,0) = 1./bmip.GetMeasure();
  }
    

    

    
    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      // static Timer t("HDivDivFE - DiffOpNormalComponent", NoTracing);
      // RegionTracer regtr(TaskManager::GetThreadId(), t);    

      auto & fel = static_cast<const HDivFiniteElement<D>&> (bfel);
      fel.CalcMappedNormalShape (bmir, mat);
      /*
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D,D> &> (bmir);
      LocalHeapMem<100000> lh("normalcomp");
      auto & fel = dynamic_cast<const HDivFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<SIMD<double>> shape(nd*D, mir.Size(), lh);
      fel.CalcMappedShape (mir, shape);
      for (size_t i = 0; i < mir.Size(); i++)
        {
          auto nv = mir[i].GetNV();
          for (size_t j = 0; j < nd; j++)
            {
              SIMD<double> sum = 0.0;
              for (size_t k = 0; k < D; k++)
                sum += shape(j*D+k, i) * nv(k);
              mat(j, i) = sum;
            }
        }
      */
    }
    
    /*
    using DiffOp<DiffOpIdHDivDiv<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).Evaluate_Matrix (mir, x, y);
    }

    using DiffOp<DiffOpIdHDivDiv<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).AddTrans_Matrix (mir, y, x);
    }    
    */
  };



  class HDivHOProlongation : public Prolongation
  {
    HDivHighOrderFESpace * fes;
    shared_ptr<FESpace> fesL2;

    mutable Array<shared_ptr<BaseMatrix>> convL2toH1;
    mutable Array<shared_ptr<BaseMatrix>> convH1toL2;
  public:
    HDivHOProlongation (HDivHighOrderFESpace * afes)
    {
      fes = afes;
      Flags flagsL2;
      flagsL2.SetFlag ("order", fes->GetOrder());
      flagsL2.SetFlag ("piola2", true);
      flagsL2.SetFlag ("hoprolongation", true);
      fesL2 = CreateFESpace("VectorL2", fes->GetMeshAccess(), flagsL2);
    }

    virtual void Update (const FESpace & cfes) override
    {
      FESpace & fes = const_cast<FESpace&>(cfes);
      fesL2->Update();
      fesL2->FinalizeUpdate();
      
      int levels = fes.GetMeshAccess()->GetNLevels();
      if (convL2toH1.Size() < levels)
        {
          convL2toH1.SetSize(levels);
          convH1toL2.SetSize(levels);
          
          LocalHeap lh(10*1000*1000);
          convL2toH1[levels-1] = ConvertOperator(fesL2, dynamic_pointer_cast<FESpace>(fes.shared_from_this()), VOL, lh, 
                                        nullptr, nullptr, NULL, nullptr, false, true, false, 0, 0, true);
          convH1toL2[levels-1] = ConvertOperator(dynamic_pointer_cast<FESpace>(fes.shared_from_this()), fesL2, VOL, lh, 
                                        nullptr, nullptr, NULL, nullptr, false, true, false, 0, 0, true);
        }
    }

    virtual size_t GetNDofLevel (int level) override
    {
      return fes->GetNDofLevel(level);
    }
  
    ///
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    { return NULL; }

    ///
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override
    {
      auto vl2 = convL2toH1[finelevel]->CreateRowVector();

      auto shapec = convH1toL2[finelevel-1]->Shape();
      auto shapef = convL2toH1[finelevel]->Shape();

      vl2.Range(get<0>(shapec)) = *convH1toL2[finelevel-1] * v.Range(get<1>(shapec));      
      fesL2->GetProlongation()->ProlongateInline(finelevel, vl2);
      v.Range(get<0>(shapef)) = *convL2toH1[finelevel] * vl2.Range(get<1>(shapef));
    }    

    
    ///
    virtual void RestrictInline (int finelevel, BaseVector & v) const override
    {
      auto vl2 = convL2toH1[finelevel]->CreateRowVector();

      auto shapec = convH1toL2[finelevel-1]->Shape();
      auto shapef = convL2toH1[finelevel]->Shape();

      vl2.Range(get<1>(shapef)) = Transpose(*convL2toH1[finelevel]) * v.Range(get<0>(shapef));      
      fesL2->GetProlongation()->RestrictInline(finelevel, vl2);
      v.Range(get<1>(shapec)) = Transpose(*convH1toL2[finelevel-1]) * vl2.Range(get<0>(shapec));
    }    
  };



  
  
  HDivHighOrderFESpace ::  
  HDivHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    type = "hdivho";
    name="HDivHighOrderFESpace(hdivho)";
    // allowed flags
    DefineNumFlag("relorder");
    // DefineNumFlag("relcurlorder");
    DefineDefineFlag("discontinuous");
    // DefineNumFlag("curlorder");
    DefineNumFlag("orderinner");
    DefineNumFlag("orderedge");
    DefineNumFlag("orderface");
    DefineNumFlag("orderfacet");
    DefineDefineFlag("hdiv");
    DefineDefineFlag("hdivho");
    // DefineDefineFlag("print");
    DefineDefineFlag("noprint");
    DefineDefineFlag("variableorder"); 
    DefineDefineFlag("hodivfree"); 

    DefineDefineFlag("RT");
    
    if(parseflags) CheckFlags(flags);

    discont = flags.GetDefineFlag("discontinuous"); 


    Flags loflags;
    loflags.SetFlag ("order", 1);
    loflags.SetFlag ("dim", dimension);
    if (iscomplex) loflags.SetFlag ("complex");
    if (discont) loflags.SetFlag ("discontinuous"); // supported ?

    // low_order_space = new RaviartThomasFESpace (ma, loflags);
    low_order_space = nullptr;
    if (flags.GetDefineFlag("loworderp1"))
      low_order_space = CreateFESpace ("BDM1", ma, loflags);
    


    // Variable order space: 
    //      in case of (var_order && order) or (relorder) 
    var_order = flags.GetDefineFlag("variableorder");  
    order =  int (flags.GetNumFlag ("order",0)); 
    curl_order =  int (flags.GetNumFlag ("curlorder",order)); 

    
    if(flags.NumFlagDefined("relorder") && !flags.NumFlagDefined("order")) 
      var_order = 1; 
    
    rel_order=int(flags.GetNumFlag("relorder",order-1)); 


    if(flags.NumFlagDefined("order") && flags.NumFlagDefined("relorder")) 
      {
	if(var_order)
	  cerr << " WARNING: HDivHoFeSpace: inconsistent flags: variableorder, order and relorder "
	       << "-> variable order space with rel_order " << rel_order << "is used, but order is ignored " << endl; 
	else 
	  cerr << " WARNING: HDivHoFeSpace: inconsistent flags: order and rel_order "
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

    // SZ hack since is not supported(tested) yet 
    rel_curl_order= rel_order; 
    curl_order = order;
        
        
    print = flags.GetDefineFlag("print"); 

    ho_div_free = flags.GetDefineFlag("hodivfree"); 
    fixed_order = flags.GetDefineFlag ("fixedorder");
    RT = flags.GetDefineFlag ("RT");

    uniform_order_inner = int (flags.GetNumFlag ("orderinner", -1));
    
    if(flags.NumFlagDefined("orderedge") || flags.NumFlagDefined ("orderface")) 
      throw Exception("Flag 'orderface' and 'orderedge' for hdivho are obsolete. Use flag 'orderfacet' instead!"); 
	
    uniform_order_facet = int (flags.GetNumFlag ("orderfacet", -1));


    auto one = make_shared<ConstantCoefficientFunction> (1);
    // integrator[VOL]= GetIntegrators().CreateBFI("masshdiv", ma->GetDimension(), one);
    // integrator[BND] = GetIntegrators().CreateBFI("robinhdiv", ma->GetDimension(), one);
    
    if (ma->GetDimension() == 2)
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<2>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdVecHDivBoundary<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<2>>>();
      }
    else
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<3>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdVecHDivBoundary<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<3>>>();
      }

    /*
    if (dimension > 1)
      {
        integrator[VOL]= make_shared<BlockBilinearFormIntegrator> (integrator[VOL], dimension);
        integrator[BND] = make_shared<BlockBilinearFormIntegrator> (integrator[BND], dimension);
      }
    */
    
    highest_order_dc = flags.GetDefineFlag("highest_order_dc");
    if (highest_order_dc) {
      *testout << "highest_order_dc is active!" << endl;
    }
    hide_all_dofs = flags.GetDefineFlag("hide_all_dofs");

    switch (ma->GetDimension())
      {
      case 1:
        additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHDiv<1>>> ()); break;
      case 2:
        additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHDiv<2>>> ());
	additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHDivDual<2>>> ());
        additional_evaluators.Set ("normalcomponent",make_shared<T_DifferentialOperator<DiffOpNormalComponentHDiv<2>>> ());
	break;
      case 3:
	additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHDiv<3>>> ());
	additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHDivDual<3>>> ());
        additional_evaluators.Set ("normalcomponent",make_shared<T_DifferentialOperator<DiffOpNormalComponentHDiv<3>>> ());
	break;
      default:
        ;
      }


   if (flags.GetDefineFlag("hoprolongation"))
        prol = make_shared<HDivHOProlongation> (this);
  }
  
  HDivHighOrderFESpace:: ~HDivHighOrderFESpace () 
  {
    ;
  }

  DocInfo HDivHighOrderFESpace :: GetDocu ()
  {
    auto docu = FESpace::GetDocu();
    docu.Arg("RT") = "bool = False\n"
      "  RT elements for simplicial elements: P^k subset RT_k subset P^{k+1}";
    docu.Arg("discontinuous") = "bool = False\n"
      "  Create discontinuous HDiv space";
    docu.Arg("hodivfree") = "bool = False\n"
      "  Remove high order element bubbles with non zero divergence";
    docu.Arg("highest_order_dc") = "bool = False\n"
      "  Activates relaxed H(div)-conformity. Allows normal discontinuity of highest order facet basis functions";
    docu.Arg("hide_all_dofs") = "bool = False\n"
      "  Set all used dofs to HIDDEN_DOFs";
    docu.Arg("orderinner") = "int = unused\n"
      "  Set order of inner shapes (orderinner=0 for no inner shapes)\n";
    return docu;
  }
  

  void HDivHighOrderFESpace :: Average (BaseVector & vec) const
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
  
  void HDivHighOrderFESpace :: Update()
  {
    FESpace::Update();
    // In discontinuous spaces the order on edges and faces  
    // are also set to the maximal order of adjacent elements  
    // and not to the element order (Motivation: Equilibrated_EE) 

    // SZ hack since is not supported(tested) yet 
    rel_curl_order= rel_order; 
    curl_order = order; 

    bool first_update = GetTimeStamp() < ma->GetTimeStamp();
    if (first_update) timestamp = NGS_Object::GetNextTimeStamp();
    
    if (low_order_space)
      low_order_space -> Update();
    
    // int nv = ma->GetNV();
    size_t nel = ma->GetNE();
    size_t nfa = ma->GetNFacets();
    size_t dim = ma->GetDimension();

    if (first_update)
      {
	order_facet.SetSize(nfa);
	order_inner.SetSize(nel); 
	order_inner_curl.SetSize(nel); 
	fine_facet.SetSize(nfa); 
	boundary_facet.SetSize(nfa); 
   
	boundary_facet = false;
	/*
	  for (int i = 0; i < ma->GetNSE(); i++)
	  {
	  Array<int> elfacets; 
	  ma->GetSElFacets (i,elfacets); 
	  boundary_facet[elfacets[0]] = true;
	  }
	*/
	for (auto el : ma->Elements(BND))
	  boundary_facet[el.Facets()] = true;

	// cout << " order hdiv " << order << endl; 
	// cout << " curl_order hdiv " << curl_order << endl; 

	int p = 0, pc = 0; 
    
	if(!var_order)
	  {
	    p = order; 
	    pc = curl_order; 
	  } 
    
	order_facet = pc;
	order_inner = p;
	order_inner_curl = pc;
	fine_facet = 0; //!!!! 


	for (auto el : ma->Elements(VOL))
	  {
	    if (!DefinedOn (el))
	      {
		order_inner[el.Nr()] = 0;
		order_inner_curl[el.Nr()] = 0;
		continue;
	      }
          
	    ELEMENT_TYPE eltype = el.GetType();
	    const POINT3D * points = ElementTopology :: GetVertices (eltype);
	
	    // Array<int> elfacets; 
	    // ma->GetElFacets (el.Nr(), elfacets); 
	    // auto elfacets = ma->GetElFacets (el);
	    auto elfacets = el.Facets();
        
	    fine_facet[elfacets] = true;
	
	    if(!var_order) continue; 
	
	    IVec<3> el_orders = ma->GetElOrders(el.Nr());  
	
	    int i = el.Nr();
	    for(int k=0;k<dim;k++)
	      {
		order_inner_curl[i][k]= max2(el_orders[k] + rel_curl_order,0);
		order_inner[i][k] = max2(el_orders[k]+rel_order,0);
	      }

	    if(dim==2)
	      {
		const EDGE * edges = ElementTopology::GetEdges (eltype);
		for(int j=0; j<elfacets.Size(); j++)
		  for(int k=0;k<2;k++)
		    if(points[edges[j][0]][k] != points[edges[j][1]][k])
		      { 
			order_facet[elfacets[j]][0] = max2(el_orders[k]+rel_curl_order, order_facet[elfacets[j]][0]);
			break; 
		      }
	      }
	    else
	      {
		// Array<int> vnums (el.Vertices());
		auto vnums = el.Vertices();
		const FACE * faces = ElementTopology::GetFaces (eltype);

		for(int j = 0; j < elfacets.Size(); j++)
		  {
		    if(faces[j][3]==-1) // trig  
		      {
			order_facet[elfacets[j]][0] = max2(order_facet[elfacets[j]][0],el_orders[0]+rel_curl_order);
			order_facet[elfacets[j]][1] = order_facet[elfacets[j]][0]; 
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
				order_facet[elfacets[j]][l] = max2(order_facet[elfacets[j]][l], rel_curl_order + 
								   el_orders[k]);
				break; 
			      } 
		    
		      }
		  }
	      }
	  }

	ma->AllReduceNodalData ((ma->GetDimension()==2) ? NT_EDGE : NT_FACE,
				fine_facet, NG_MPI_LOR);

	if(uniform_order_inner > -1) order_inner = uniform_order_inner;
	if(uniform_order_facet > -1) order_facet = uniform_order_facet;

	for (auto i : Range(nfa)) if (!fine_facet[i]) order_facet[i] = IVec<2> (0,0); 

	// by SZ ... since order_inner_curl is not working yet for hdivhofe
	order_inner_curl = order_inner;
      }
    
    if(print) 
      {
	*testout << " discont " << discont << endl;
	*testout << " fine_facet[i] (hdiv) " << fine_facet << endl; 
	
	*testout << " order_facet (hdivho) " << order_facet << endl; 
	*testout << " order_inner (hdivho) " << order_inner << endl; 
	*testout << " order_inner_curl (hdivho) " << order_inner_curl << endl; 
      }
    
    UpdateDofTables(); 
    UpdateCouplingDofArray();

    if (low_order_space)
      {
        if (low_order_space->GetClassName()=="BDM11FESpace")
          {
            cout << "make p1 embedding" << endl;
            Array<int> ia, ja;
            Array<double> va;
            for (auto f : Range(ma->GetNFacets()))
              {
                ia.Append(f);
                ja.Append(3*f);
                va.Append(1);

                auto fd = GetFacetDofs(f);
                ia.Append(fd[0]);
                ja.Append(3*f+1);
                va.Append(1);
                
                ia.Append(fd[0]+order);
                ja.Append(3*f+2);
                va.Append(1);
              }
            low_order_embedding = SparseMatrix<double>::CreateFromCOO (ia, ja, va,
                                                                       GetNDof(), low_order_space->GetNDof());
            // cout << "low_order_embeddingn = " << *low_order_embedding << endl;
          }
      }
  }

  void HDivHighOrderFESpace :: UpdateDofTables()
  {
    size_t nfa = ma->GetNFacets();
    size_t nel = ma->GetNE();
    size_t dim = ma->GetDimension();
    // Array<int> pnums; 
     
    first_facet_dof.SetSize(nfa+1); 
    first_inner_dof.SetSize(nel+1); 

    size_t ndof = nfa;
    first_facet_dof = ndof; 

    if(dim==2)
      {
	// int dec_hodc = highest_order_dc ? 1 : 0;
        for (auto i : Range(nfa))
          {
            first_facet_dof[i] = ndof;
            int inc = fine_facet[i] ? order_facet[i][0] : 0;
	    if (highest_order_dc && !boundary_facet[i]) inc--;
            if (inc > 0) ndof += inc;
	   
            /*
            if(fine_facet[i])
	      ndof += order_facet[i][0];
	    if (highest_order_dc && !boundary_facet[i]) ndof--;
            */
          }

        first_facet_dof[nfa] = ndof;
      
	// Array<int> fnums;
        // for (auto i : Range(nel))
        for (size_t i = 0; i < nel; i++)
          {
            ElementId ei(VOL, i);
            IVec<3> pc = order_inner_curl[i];
            IVec<3> p = order_inner[i];
            int inci = 0;
            switch(ma->GetElType(ei))
              {
              case ET_TRIG:
                if (!ho_div_free)
                  inci = pc[0]*(pc[0]-1)/2 + p[0]*(p[0]-1)/2 + p[0]-1;
                else
                  inci = pc[0]*(pc[0]-1)/2;
                if (RT)
                  inci += pc[0] + 1;
                break;
              case ET_QUAD:
                if (!ho_div_free)
                  inci = pc[0]*pc[1] + p[0]*p[1]+p[0]+p[1];
                else
                  inci = pc[0]*pc[1];
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
            dc_pairs = IVec<2> (-1,-1);
            
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
      }
    else 
      {
        int inci = 0;
        for (size_t i = 0; i < nfa; i++) 
          {
            inci = 0; 
            if(fine_facet[i])
              {
                IVec<2> p = order_facet[i]; 
                // ma->GetFacePNums(i,pnums);
                auto pnums = ma->GetFacePNums(i);
                switch(pnums.Size())
                  {
                  case 3: //Triangle
                    inci= (p[0]*p[0]+3*p[0])/2; 
                    if (highest_order_dc && !boundary_facet[i]) inci -= p[0]+1;
                    break;
                  case 4: //Quad 
                    inci= p[0]*p[1] + p[0] + p[1];
                    if (highest_order_dc && !boundary_facet[i]) inci -= p[0]+p[1]+1;
                    break;
                  }
              }
            else
              inci = 0;

            if (inci < 0) inci = 0;
            first_facet_dof[i] = ndof;
            ndof+= inci;
	    
          }
        first_facet_dof[nfa] = ndof;
	 
	// Array<int> fnums;
        bool have_pyramids = false;
        for (size_t i = 0; i < nel; i++)
          {
            IVec<3> p = order_inner[i];
            IVec<3> pc = order_inner_curl[i];
            int inci = 0;
	     
            switch(ma->GetElType(ElementId(VOL,i)))
              {
              case ET_TET:
                if(p[0]>1 && !ho_div_free)
                  inci = p[0]*(p[0]+1)*(p[0]-1)/6 + p[0]*(p[0]-1)/2 + p[0]-1;
                if(pc[0]>1)
                  inci += pc[0]*(pc[0]+1)*(pc[0]-1)/3 + pc[0]*(pc[0]-1)/2; ;
                if (highest_order_dc) 
		  {
		    // ma->GetElFacets (i, fnums);
                    /*
                    auto fnums = ma->GetElFacets(i);
		    for (int j = 0; j < fnums.Size(); j++)
                      {
                        IVec<2> pf = order_facet[fnums[j]]; 
                        if (!boundary_facet[fnums[j]]) inci += pf[0]+1;
                      }
                    */
                    for (auto f : ma->GetElFacets(i))
                      if (!boundary_facet[f]) inci += order_facet[f][0]+1;
		  }
		if (RT)
		  inci += (p[0]+1) * (p[0]+2)/2;
		// inci += 4*(p[0]+1);
                break;
              case ET_PRISM:
                // inci = (p[0]+1)*(3*(p[0]-1)+(p[0]-1)*(p[0]-2))+ (p[0]-1)*(p[0]+1)*(p[0]+2)/2;
		inci = (p[0]+2)*p[0]*(p[2]+1) + (p[0]+1)*(p[0]+2)*p[2]/2;
		if (ho_div_free)
		  inci -= (p[0]+1)*(p[0]+2)*(p[2]+1)/2 - 1;

                if (highest_order_dc) inci += 2*(p[0]+1)+3*(p[0]+p[2]+1);
                break;
              case ET_HEX:
                inci =  2*pc[0]*pc[1]*pc[2] + pc[0]*pc[1] + pc[1]*pc[2] + pc[0]* pc[2]
                        + p[0]*(p[1]*p[2] + p[1] + p[2] + 1)  + p[1]*p[2] + p[1] + p[2]; 
                if (ho_div_free)
                  inci -= p[0]*(p[1]*p[2] + p[1] + p[2] + 1)  + p[1]*p[2] + p[1] + p[2]; 
                break; 
              case ET_PYRAMID: 
                inci=0;
                have_pyramids=true;
                break; 
              default:
                inci = 0;
                break;
              }
            if (inci < 0) inci = 0;
            first_inner_dof[i] = ndof;
            ndof+= inci;
	  }
        first_inner_dof[nel] = ndof;

        if (have_pyramids)
          cout << "WARNING: there are hdiv-pyramids (not implemented yet) !! " << endl;

        
        if (highest_order_dc)
          {
            dc_pairs.SetSize ((order+1)*ma->GetNFacets());
            dc_pairs = IVec<2> (-1,-1);
            
            for (int i = 0; i < ma->GetNE(); i++)
              {
                ElementId ei(VOL,i);
                auto fnums = ma->GetElFacets (ei);
		int di = first_inner_dof[i]; // +k*(order+1);
		    
                for (int k = 0; k < fnums.Size(); k++)
		  if (!boundary_facet[fnums[k]])
		    {
		      int base = fnums[k]*(order+1);
		      
		      for (int l = 0; l < order+1; l++)
			{
			  dc_pairs[base+l][1] = dc_pairs[base+l][0];
			  dc_pairs[base+l][0] = di++;
			}
		    }
              }
          }
        else
          dc_pairs.SetSize (0);

      }
   

    if(discont) 
      { 
        ndof = 0; 
        for(int i=0;i<nel;i++)
          {
            ElementId ei(VOL,i);
            int incii = first_inner_dof[i+1]-first_inner_dof[i]; 
	     	     
            auto elfacets = ma->GetElFacets (ei);
	     
            for(int j=0; j<elfacets.Size(); j++) 
              incii+=first_facet_dof[elfacets[j]+1]-first_facet_dof[elfacets[j]]+1; // incl. lowest-order 
	     
            first_inner_dof[i] = ndof; 
            ndof += incii;
	    
          } 
        first_inner_dof[nel] = ndof; 
        first_facet_dof = 0; 
      } 
	
    
    if(print) 
      {
        (*testout) << "ndof hdivhofespace update = " << endl << ndof << endl;
        (*testout) << "first_facet_dof (hdiv)  = " << endl << first_facet_dof << endl;
        (*testout) << "first_facet_dof (hdiv) = " << endl << first_facet_dof << endl;
        (*testout) << "first_inner_dof (hdiv) = " << endl << first_inner_dof << endl;
      }

    SetNDof (ndof);
    /*
    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;
    */
    //    prol->Update();
  }

  void HDivHighOrderFESpace :: UpdateCouplingDofArray()
  {
    auto wirebasket_ct = hide_all_dofs ? HIDDEN_DOF : WIREBASKET_DOF;
    auto interface_ct = hide_all_dofs ? HIDDEN_DOF : INTERFACE_DOF;
    auto local_ct = hide_all_dofs ? HIDDEN_DOF : LOCAL_DOF;
    ctofdof.SetSize(GetNDof());
    if(discont) 
      {
        ctofdof = local_ct;
        return;
      } 
    
    ctofdof = wirebasket_ct;
    
    for (auto facet : Range(ma->GetNFacets()))
      {
        ctofdof[facet] = fine_facet[facet] ? wirebasket_ct : UNUSED_DOF;
        ctofdof[GetFacetDofs(facet)] = interface_ct;
      }
    
    for (auto el : Range (ma->GetNE()))
      ctofdof[GetElementDofs(el)] = local_ct;
  }


  void HDivHighOrderFESpace :: SetOrder (NodeId ni, int order) 
  {
    if (order_policy == CONSTANT_ORDER || order_policy == NODE_TYPE_ORDER)
      throw Exception("In HDivHighOrderFESpace::SetOrder. Order policy is constant or node-type!");
    else if (order_policy == OLDSTYLE_ORDER)
      order_policy = VARIABLE_ORDER;
      
    if (order < 0)
      order = 0;
    
    switch( CoDimension(ni.GetType(), ma->GetDimension()) )
      {
      case 1:
	if (ni.GetNr() < order_facet.Size())
	  order_facet[ni.GetNr()] = fine_facet[ni.GetNr()] ? order : 0;
	break;
      case 0:
        if (ma->GetDimension() == 2 && ni.GetType() == NT_FACE)
	  {
	    Array<int> elnr;
	    ma->GetFacetSurfaceElements(ni.GetNr(),elnr);
	    if (elnr[0] < order_inner.Size())
	      {
		order_inner[elnr[0]] = order;
		order_inner_curl[elnr[0]] = order;
	      }
	  }
        else if (ni.GetNr() < order_inner.Size())
	  {
	    order_inner[ni.GetNr()] = order;
	    order_inner_curl[ni.GetNr()] = order;
	  }
	break;
      default:
	break;
      }
  }
  
  int HDivHighOrderFESpace :: GetOrder (NodeId ni) const
  {
    switch( CoDimension(ni.GetType(), ma->GetDimension()) )
      {
      case 1:
	if (ni.GetNr() < order_facet.Size())
	  return order_facet[ni.GetNr()][0];
	break;
      case 0:
	if (ma->GetDimension() == 2 && ni.GetType() == NT_FACE)
	  {
	    Array<int> elnr;
	    ma->GetFacetSurfaceElements(ni.GetNr(),elnr);
	    if (elnr[0] < order_inner.Size())
	      return order_inner[elnr[0]][0];
	  }
        else if (ni.GetNr() < order_inner.Size())
	  return order_inner[ni.GetNr()][0];
	break;
      default:
	break;
      }
    
    return 0;
  }


  template <ELEMENT_TYPE ET>
  FiniteElement & HDivHighOrderFESpace :: T_GetFE (int elnr, bool onlyhdiv, Allocator & lh) const
  {
    Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL> (elnr);
    if (!DefinedOn(ngel)) return * new (lh) HDivDummyFE<ET>();
    
    HDivHighOrderFE<ET> * hofe =  new (lh) HDivHighOrderFE<ET> ();

    hofe -> SetVertexNumbers (ngel.Vertices());
    hofe -> SetHODivFree (ho_div_free);
    hofe -> SetOnlyHODiv (onlyhdiv);
    hofe -> SetRT(RT);
    hofe -> SetOrderInner (order_inner[elnr]);
        
    switch (int(ET_trait<ET>::DIM))
      {
      case 2:
        hofe -> SetOrderFacet (order_facet[ngel.Edges()]);
        break;
      case 3:
        hofe -> SetOrderFacet (order_facet[ngel.Faces()]);
        break;
      }
    hofe -> ComputeNDof();
    return *hofe;
  }

  FiniteElement & HDivHighOrderFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    if (ei.IsVolume())
      {
        int elnr = ei.Nr();
        Ngs_Element ngel = ma->GetElement(ei);
        ELEMENT_TYPE eltype = ngel.GetType();
        
        switch (eltype)
          {
            // case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, false, alloc);
            
          case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, false, alloc);
          case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, false, alloc);
            
          case ET_TET:     return T_GetFE<ET_TET> (elnr, false, alloc);
          case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, false, alloc);
            // case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, false, alloc);
          case ET_HEX:     return T_GetFE<ET_HEX> (elnr, false, alloc);
            
          default:
            throw Exception ("illegal element in HDivHOFESpace::GetFE");
          }
      }
    else if (ei.VB()  == BND)
      {
        if (!DefinedOn(ei))
          return SwitchET(ma->GetElType(ei), [&] (auto et) -> FiniteElement&
                          {
                            return * new (alloc) HDivNormalDummyFE<et.ElementType()>();
                          });
        
        int selnr = ei.Nr();
        FiniteElement * fe = 0;
        
        int porder; 
        if (discont) porder = -1; 
        else porder = order; 
        
        // if (highest_order_dc) porder--;
        
        auto vnums = ma->GetElVertices(ei);
        
        switch (ma->GetElType(ei))
          {
          case ET_SEGM:
            {
              auto fe1 = new (alloc) HDivHighOrderNormalSegm<TrigExtensionMonomial> (porder);
              fe1 -> SetVertexNumbers(vnums);
              fe = fe1;
              break;
            }
          case ET_TRIG:
            {
              auto fe1 = new (alloc) HDivHighOrderNormalTrig<TrigExtensionMonomial> (porder);
              fe1 -> SetVertexNumbers(vnums);              
              fe = fe1;
              break;
            }
          case ET_QUAD:
            {
              auto fe1 = new (alloc) HDivHighOrderNormalQuad<TrigExtensionMonomial> (porder);
              fe1 -> SetVertexNumbers(vnums);              
              fe = fe1;
              break;
            }
          default:
            throw Exception (string("HDivHighOrderFESpace::GetSFE: unsupported element ")+
                             ElementTopology::GetElementName(ma->GetElType(ei)));
          }
        
        if (discont) return *fe; 
        
        // ArrayMem<int, 4> ednums, order_ed;
        // IVec<3> order_fa;
        
        if(ma->GetElType(ei) == ET_SEGM)
          {
            HDivHighOrderNormalFiniteElement<1> * hofe =
              dynamic_cast<HDivHighOrderNormalFiniteElement<1>*> (fe);
            
            // hofe -> SetVertexNumbers (vnums);
            auto ednums = ma->GetElEdges(ei);
            // int dec = (!boundary_facet[ednums[0]] && highest_order_dc) ? 1 : 0;
            hofe -> SetOrderInner (order_facet[ednums[0]][0] /* -dec */);
            hofe -> ComputeNDof();
          }
        else
          {
            HDivHighOrderNormalFiniteElement<2> * hofe =
              dynamic_cast<HDivHighOrderNormalFiniteElement<2>*> (fe);
            
            // hofe -> SetVertexNumbers (vnums);
            
#ifdef NEW_HDIVFE
            IVec<3> order_fa = IVec<3>(order_facet[ma->GetSElFace(selnr)][0],
                                     order_facet[ma->GetSElFace(selnr)][1],0);
            if (highest_order_dc)
              {
                order_fa[0]--;
                order_fa[1]--;
                order_fa[2]--;
              }
            hofe -> SetOrderInner (order_fa);
#else 
            int order_fa = order_facet[ma->GetSElFace(selnr)][0];
            // if (highest_order_dc) order_fa--;
            hofe -> SetOrderInner (order_fa);
#endif
            hofe -> ComputeNDof();
          }
        
        return *fe;
      }
    else
      switch (ma->GetElement(ei).GetType())
        {
        case ET_POINT: return * new (alloc) DummyFE<ET_POINT>();
        case ET_SEGM: return * new (alloc) DummyFE<ET_SEGM>();
        default:
          __assume(false);
          throw Exception("HDiv - impossible element");
        }
  }
  
  // const FiniteElement & HDivHighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  // {
  //   Ngs_Element ngel = ma->GetElement(elnr);
  //   ELEMENT_TYPE eltype = ngel.GetType();

  //   /*
  //   if (ma->GetElType(elnr) == ET_TRIG && order <= 6 && fixed_order)
  //     {
  //       HDivHighOrderFiniteElementFO<2> * hofe2d = 0;
  //       switch (order)
  //         {
  //         case 1: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,1> (); break;
  //         case 2: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,2> (); break;
  //         case 3: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,3> (); break;
  //         case 4: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,4> (); break;
  //         case 5: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,5> (); break;
  //         case 6: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,6> (); break;
  //         }
	
  //       Ngs_Element ngel = ma->GetElement<2> (elnr);
  //       for (int j = 0; j < 3; j++)
  //         hofe2d->SetVertexNumber (j, ngel.vertices[j]);

  //       hofe2d -> SetHODivFree (ho_div_free);
  //       hofe2d -> SetOnlyHODiv (false);
  //       hofe2d -> ComputeNDof();
	
  //       return *hofe2d;
  //     }  
  //   */

  //   switch (eltype)
  //     {
  //       // case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, false, lh);
        
  //     case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, false, lh);
  //     case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, false, lh);
        
  //     case ET_TET:     return T_GetFE<ET_TET> (elnr, false, lh);
  //     case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, false, lh);
  //           // case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, false, lh);
  //           // case ET_HEX:     return T_GetFE<ET_HEX> (elnr, false, lh);
        
  //     default:
  //       throw Exception ("illegal element in HDivHOFESpace::GetFE");
  //     }
  // }

  const FiniteElement & HDivHighOrderFESpace :: GetHODivFE (int elnr, LocalHeap & lh) const
  {
    Ngs_Element ngel = ma->GetElement(ElementId(VOL,elnr));
    ELEMENT_TYPE eltype = ngel.GetType();
    
    if (!ho_div_free) throw Exception ("You don't have hodivfree active. You are not allow to call GetHODivFE");
    
    /*
    if (ma->GetElType(elnr) == ET_TRIG && order <= 6)
      {
	HDivHighOrderFiniteElementFO<2> * hofe2d = 0; 
	switch (order)
	  {
	  case 1: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,1> (); break;
	  case 2: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,2> (); break;
	  case 3: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,3> (); break;
	  case 4: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,4> (); break;
	  case 5: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,5> (); break;
	  case 6: hofe2d = new (lh)  HDivHighOrderFEFO<ET_TRIG,6> (); break;
	  }
	
	Ngs_Element ngel = ma->GetElement<2> (elnr);
	for (int j = 0; j < 3; j++)
	  hofe2d->SetVertexNumber (j, ngel.vertices[j]);

        hofe2d -> SetOnlyHODiv (true);
        hofe2d -> ComputeNDof();
	
	return *hofe2d;
      }  
    */

    switch (eltype)
      {
        // case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, true, lh);
        
      case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, true, lh);
      case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, true, lh);
        
      case ET_TET:     return T_GetFE<ET_TET> (elnr, true, lh);
      case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, true, lh);
        // case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, false, lh);
      case ET_HEX:     return T_GetFE<ET_HEX> (elnr, true, lh);
        
      default:
        throw Exception ("illegal element in HDivHOFeSpace::GetDivFE");
      }
  }



//   const FiniteElement & HDivHighOrderFESpace :: GetSFE (int selnr, LocalHeap & lh) const
//   {
//     FiniteElement * fe = 0;

//     int porder; 
//     if (discont) porder = -1; 
//     else porder = order; 

//     // if (highest_order_dc) porder--;

//     switch (ma->GetSElType(selnr))
//       {
//       case ET_SEGM:
//         fe = new (lh) HDivHighOrderNormalSegm<TrigExtensionMonomial> (porder); 
//         break;
//       case ET_TRIG: 
//         fe = new (lh) HDivHighOrderNormalTrig<TrigExtensionMonomial> (porder); 
//         break; 
//       case ET_QUAD: 
//         fe = new (lh) HDivHighOrderNormalQuad<TrigExtensionMonomial> (porder); 
//         break; 
//       default:
//         throw Exception (string("HDivHighOrderFESpace::GetSFE: unsupported element ")+
//                          ElementTopology::GetElementName(ma->GetSElType(selnr)));
//       }

//     if (discont) return *fe; 

//     ArrayMem<int,4> vnums;
//     ArrayMem<int, 4> ednums, order_ed;
//     IVec<3> order_fa;
//     ma->GetSElVertices(selnr, vnums);
    
//     if(ma->GetSElType(selnr) == ET_SEGM)
//       {
// 	HDivHighOrderNormalFiniteElement<1> * hofe =
// 	  dynamic_cast<HDivHighOrderNormalFiniteElement<1>*> (fe);

// 	hofe -> SetVertexNumbers (vnums);
// 	ma->GetSElEdges(selnr, ednums);
// 	// int dec = (!boundary_facet[ednums[0]] && highest_order_dc) ? 1 : 0;
// 	hofe -> SetOrderInner (order_facet[ednums[0]][0] /* -dec */);
// 	hofe -> ComputeNDof();
//       }
//     else
//       {
// 	HDivHighOrderNormalFiniteElement<2> * hofe =
// 	  dynamic_cast<HDivHighOrderNormalFiniteElement<2>*> (fe);

// 	hofe -> SetVertexNumbers (vnums);
	
// #ifdef NEW_HDIVFE
// 	IVec<3> order_fa = IVec<3>(order_facet[ma->GetSElFace(selnr)][0],
//                                  order_facet[ma->GetSElFace(selnr)][1],0);
//         if (highest_order_dc)
//           {
//             order_fa[0]--;
//             order_fa[1]--;
//             order_fa[2]--;
//           }
// 	hofe -> SetOrderInner (order_fa);
// #else 
// 	int order_fa = order_facet[ma->GetSElFace(selnr)][0];
//         // if (highest_order_dc) order_fa--;
// 	hofe -> SetOrderInner (order_fa);
// #endif
// 	hofe -> ComputeNDof();
//       }
    
//     return *fe;
//   }

  /*
  size_t HDivHighOrderFESpace :: GetNDof () const throw()
  {
    return ndof;
  }

  size_t HDivHighOrderFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }
  */


  void HDivHighOrderFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if (!DefinedOn (ei)) return;

    /*
    if (order_policy == NODE_TYPE_ORDER)
      {
        auto et = ma->GetElType(ei);
        cout << "new order policy, et = " << et
             << ", ol = " << et_order_left[et]
             << ", or = " << et_order_right[et]
             << endl;
      }
    */
    
    if(discont) 
      {
	// lowest_order included in inner
	if(ei.VB()==VOL)
	  dnums += GetElementDofs (ei.Nr());
	return;
      } 
    if(ei.VB()==VOL)
      {
	// ArrayMem<int,6> fanums;
	// ma->GetElFacets (ei, fanums);
        auto fanums = ma->GetElFacets(ei);
	if (highest_order_dc)
	  {
	    if (ma->GetDimension() == 2)
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
	    else // dim not 2
	      {
		IntRange eldofs = GetElementDofs (ei.Nr());
		
		// for (int i = 0; i < fanums.Size(); i++)
		// dnums.Append (fanums[i]);
		dnums += fanums;
		
		int firstel = eldofs.First();
		
		for(int i = 0; i < fanums.Size(); i++)
		  {
		    int firstfa = GetFacetDofs(fanums[i]).First();
		    
		    if (!boundary_facet[fanums[i]])
		      {
			for (int i = 0; i <= order-1; i++)
			  {
			    for(int j = 0; j < order-i-1; j++)
			      dnums.Append(firstfa++);
			    dnums.Append(firstel++);
			  }
			for (int i = 0; i < order-1; i++)
			  dnums.Append(firstfa++);
			dnums.Append(firstel++);
		      }
		    else
		      dnums += GetFacetDofs (fanums[i]);
		  }
		dnums += IntRange (firstel, eldofs.Next());
	      }
	  }
	else // not highest order dc
	  {
	    //Raviart-Thomas
            dnums += fanums;
	    // facets
	    for (auto f : fanums)
	      dnums += GetFacetDofs (f);
	    
	    //inner
	    dnums += GetElementDofs (ei.Nr());
	  }
      }
    else if (ei.VB()==BND)
      {
        size_t fanum = ma->GetElFacets(ei)[0];
	// lowest-order
        dnums += fanum;
	
	// high order
        dnums += GetFacetDofs (fanum);
        
      }
  }



  // ****************************
  // 
  //    smoothing blocks
  //
  //  0) Jacobi
  //  1) 2d-Vertex / 3d-Edge blocks + F + I  --- default
  //  2) 2d: edge by edge,  3d: face by face

  shared_ptr<Table<int>> HDivHighOrderFESpace :: CreateSmoothingBlocks (const Flags & precflags) const
  {
    {
      // auto blocktype = precflags.GetStringFlag ("blocktype", "");
      // if (blocktype == "edgepatch" || blocktype == "facepatch" || blocktype == "vertexpatch")
      if (precflags.StringFlagDefined("blocktype") || precflags.StringListFlagDefined("blocktype"))
        return FESpace::CreateSmoothingBlocks(precflags);
    }
    
    int first;
    int ncnt = 0;
    // int ni = ma->GetNE(); //nel;
  
    int dim = ma->GetDimension();
  
  
    int SmoothingType = int(precflags.GetNumFlag("blocktype",1)); 
  
  
    Array<int> vnums,elnums; 
    Array<int> orient; 
    Array<int> ednums, fanums;
    Array<int> edges(3);
  
    int ned = ma->GetNEdges();
    int nfa = ma->GetNFaces();
    int nnv = ma->GetNV();
    int nel = ma->GetNE();
    size_t ndof = GetNDof();
    
    cout << "SmoothingType " << SmoothingType << endl; 
    switch(SmoothingType) 
      {
      case 0: 
	cout << "Local Preconditioner" << endl;
	ncnt = ndof; //_used;
	break;  
      case 1:
	if( dim == 2 )
	  {
	    cout << "Vertex blocks + E + I" << endl;
	    ncnt = nnv + ned + nel ;
	  }
	else
	  {
	    cout << "Edge blocks + F + I" << endl;
	    ncnt = ned + nfa + nel;
	  }
	break;

      case 2:
	if( dim == 2 )
	  {
	    cout << "Edge blocks" << endl;
	    ncnt = ned ;
	  }
	else
	  {
	    cout << "Face blocks" << endl;
	    ncnt = nfa;
	  }
	break;

      }


    Array<int> cnt(ncnt); 
    cnt = 0;
    // ii=0; 

    int offset = 0;




    switch (SmoothingType)
      {
	//      0 ..... jacobi
      case 0:   // diagonal
	for(int i=0; i<ndof; i++ )
          cnt[i] = 1;
	break;
	
        //      1 ..... vertex/edge blocks -- E -- I
      case 1:   
	if( dim == 2 )
	  {
	    // vertex blocks
	    for(int i=0; i<ned; i++)
	      if(fine_facet[i])
                {
                  /*
                  int pn1, pn2;
                  ma->GetEdgePNums ( i, pn1, pn2);
                  cnt[offset + pn1] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
                  cnt[offset + pn2] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
                  */
                  auto pn = ma->GetEdgePNums(i);
                  cnt[offset+pn[0]] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
                  cnt[offset+pn[1]] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
                }

	    offset += nnv;
	    // edges
	    for (auto i : Range(ned))
	      if (fine_facet[i])
                cnt[offset+i] += first_facet_dof[i+1] - first_facet_dof[i];;
	    offset += ned;

	    // cells
	    for (auto i : Range(nel))
              cnt[offset+i] += first_inner_dof[i+1] - first_inner_dof[i];;
	  }

	else
	  {
	    // vertex blocks
	    for(auto i : Range(nfa))
	      if(fine_facet[i])
		{
		  // Array<int> edges;
		  // ma->GetFaceEdges ( i, edges);
                  /*
		  for ( int j = 0; j < edges.Size(); j++ )
		    cnt[offset + edges[j]] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
                  */
                  for (auto e : ma->GetFaceEdges(i))
                    cnt[offset+e] += 1 + first_facet_dof[i+1] - first_facet_dof[i];
		}
	    
	    offset += ned;
	    // edges
	    for (auto i : Range(nfa))
	      if (fine_facet[i])
                cnt[offset+i] += first_facet_dof[i+1] - first_facet_dof[i];;
	    offset += nfa;

	    // cells
	    for (auto i : Range(nel))
              cnt[offset+i] += first_inner_dof[i+1] - first_inner_dof[i];;
	  }
	break;
      case 2:
	if( dim == 2 )
          cerr << "not implemented" << endl;
	else
	  {
	    for (auto i : Range(nfa))
	      if (fine_facet[i])
		cnt[i] += first_facet_dof[i+1] - first_facet_dof[i];
	  }
      }


    Table<int> table(cnt); 
  
    // ii = 0; 
    cnt = 0;
  
    offset =0;
  
    switch(SmoothingType) 
      {
	//      0 ..... jacobi
      case 0:
	for(int i=0; i<ndof; i++)
	  table[i][cnt[i]] = i;
	break;
	
	//      1 ..... vertex/edge blocks -- E -- I
      case 1:

	if ( dim == 2 )
	  {
	    for (int i = 0; i < ned; i++)
	      {
		auto pts = ma->GetEdgePNums (i);
		int pn1 = pts[0], pn2 = pts[1];
		if( fine_facet[i] )
		  {
		    table[offset + pn1][cnt[offset + pn1]++] = i;
		    table[offset + pn2][cnt[offset + pn2]++] = i;
		    
		    first = first_facet_dof[i];
		    int last = first_facet_dof[i+1];
		    for( int l=first; l<last; l++)
		      {
			table[offset + pn1][cnt[offset + pn1]++] = l;
			table[offset + pn2][cnt[offset + pn2]++] = l;
		      }
		  }
	      }

	    offset += nnv;

	    for (int i = 0; i < ned; i++ )
	      {
		first = first_facet_dof[i];
		int last = first_facet_dof[i+1];
		for ( int l = first; l < last; l++ )
		  table[offset + i ] [cnt[offset+i]++] = l;
	      }

	    for (int i = 0; i < nel; i++ )
	      {
		first = first_inner_dof[i];
		int last = first_inner_dof[i+1];
		for ( int l = first; l < last; l++ )
		  table[offset + i ] [cnt[offset+i]++] = l;
	      }
	    //(*testout) << table << endl;
	  }

	else // 3d

	  {
	    for (int i = 0; i < nfa; i++)
	      {
		if ( ! fine_facet[i] ) continue;
		// Array<int> faces;
		// ma->GetFaceEdges (i,faces);	      
		auto faces = ma->GetFaceEdges(i);
	
		for ( int j = 0; j < faces.Size(); j++ )
		  {
		    table[offset + faces[j]][cnt[offset + faces[j]]++] = i;
		    
		    first = first_facet_dof[i];
		    int last = first_facet_dof[i+1];
		    for( int l=first; l<last; l++)
		      {
			table[offset + faces[j]][cnt[offset + faces[j]]++] = l;
		      }
		  }
	      }

	    offset += ned;

	    for (int i = 0; i < nfa; i++ )
	      {
		first = first_facet_dof[i];
		int last = first_facet_dof[i+1];
		for ( int l = first; l < last; l++ )
		  table[offset + i ] [cnt[offset+i]++] = l;
	      }

            offset += nfa;
            
	    for (int i = 0; i < nel; i++ )
	      {
		first = first_inner_dof[i];
		int last = first_inner_dof[i+1];
		for ( int l = first; l < last; l++ )
		  table[offset + i ] [cnt[offset+i]++] = l;
	      }
	    //(*testout) << table << endl;
	    
	  }
	break;

      case 2:
	
	if ( dim == 2 )
	  {
	    cout << "not implemented" << endl;
	  }
	
	else // 3d

	  {
	    for (int i = 0; i < nfa; i++ )
	      {
		first = first_facet_dof[i];
		int last = first_facet_dof[i+1];
		for ( int l = first; l < last; l++ )
		  table[i] [cnt[i]++] = l;
	      }
	  }
	break;
      }
    


    //  *testout << table << endl;
    // cout << "success " << endl;
    return make_shared<Table<int>> (table);

  }







  // ******************** //
  // Direct Solver Clusters
  // 0) none
  // 1) low order dofs  --  default


  shared_ptr<Array<int>> HDivHighOrderFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    auto spclusters = make_shared<Array<int>> (GetNDof());
    Array<int> & clusters = *spclusters;

    int default_ds = low_order_space ? 0 : 1;
    int clustertype = int(precflags.GetNumFlag("ds_cluster",default_ds)); 
    cout << IM(3) << " DirectSolverCluster Clustertype " << clustertype << endl; 
  
    Array<int> vnums,elnums; 
    Array<int> orient; 
  
    Array<int> edges(3);
        
    int nfa = ma->GetNFaces();

    Array<int> ednums, fnums, pnums;
    auto freedofs = GetFreeDofs();
  
    switch (clustertype)
      {
        // 0) none
      case 0:
        clusters = 0;
        break;
        
        // 1) low-order dofs
      case 1: 
        clusters = 0;

        for(int i=0; i<nfa; i++ )
          if( fine_facet[i] && freedofs->Test(i))
            clusters[i] = 1;
        break;

      }
    return spclusters;

  }  

  
  void HDivHighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize0(); 
  }
  
  void HDivHighOrderFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  { 
    dnums.SetSize0();
    if(ma->GetDimension() == 3 || discont) return; 

    dnums += ednr;
    dnums += GetFacetDofs (ednr);
  }

  void HDivHighOrderFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 2 || discont) return; 
   
    dnums += fanr;
    dnums += GetFacetDofs (fanr);
  }

  void HDivHighOrderFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums = GetElementDofs (elnr);
  }


  static RegisterFESpace<HDivHighOrderFESpace> init ("hdivho");
}



