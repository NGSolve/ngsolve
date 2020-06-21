/**********************************************************************/
/* File:   interpolate.cpp                                            */
/* Author: L Kogler, M Neunteufel, J Schoeberl                        */
/* Date:   June 2020                                                  */
/**********************************************************************/

/* 
   Interpolation of CoefficientFunctions using
   dual shapes
*/

#include <comp.hpp>
#include <variant>
#include "../fem/integratorcf.hpp"


namespace ngcomp
{
  
  class InterpolationCoefficientFunction : public T_CoefficientFunction<InterpolationCoefficientFunction>
  {
  protected:
    shared_ptr<CoefficientFunction> func;
    shared_ptr<FESpace> fes;
    int bonus_intorder;

    Array<shared_ptr<BilinearFormIntegrator>> bli;
    Array<shared_ptr<BilinearFormIntegrator>> single_bli;
    // shared_ptr<CoefficientFunction> dual_diffop;
    shared_ptr<DifferentialOperator> dual_diffop;

    VorB vb;
    
  public:
    InterpolationCoefficientFunction (shared_ptr<CoefficientFunction> f, shared_ptr<FESpace> afes,
                                      int abonus_intorder)
      : T_CoefficientFunction<InterpolationCoefficientFunction>(f->Dimension(), f->IsComplex()),
      func(f), fes(afes), bonus_intorder(abonus_intorder)
    {
      this->SetDimensions (func->Dimensions());
      this->elementwise_constant = func->ElementwiseConstant();

      // copied from Set (dualshapes)
      
      vb = VOL;    // for the moment only 

      /** Trial-Proxy **/
      auto single_evaluator =  fes->GetEvaluator(vb);
      if (dynamic_pointer_cast<BlockDifferentialOperator>(single_evaluator))
        single_evaluator = dynamic_pointer_cast<BlockDifferentialOperator>(single_evaluator)->BaseDiffOp();
      auto trial = make_shared<ProxyFunction>(fes, false, false, single_evaluator,
                                              nullptr, nullptr, nullptr, nullptr, nullptr);

      /** Test-Proxy (dual) **/
      auto dual_evaluator = fes->GetAdditionalEvaluators()["dual"];
      for (VorB avb = VOL; avb < vb; avb++) {
        dual_evaluator = dual_evaluator->GetTrace();
        if ( dual_evaluator == nullptr )
          { throw Exception(fes->GetClassName() + string(" has no dual trace operator for vb = ") + \
                            to_string(avb) + string(" -> ") + to_string(avb + 1) + string("!")); }
      }
      if (dynamic_pointer_cast<BlockDifferentialOperator>(dual_evaluator))
        dual_evaluator = dynamic_pointer_cast<BlockDifferentialOperator>(dual_evaluator)->BaseDiffOp();
      auto dual = make_shared<ProxyFunction>(fes, true, false, dual_evaluator,
                                             nullptr, nullptr, nullptr, nullptr, nullptr);

      dual_diffop = dual_evaluator;

      for (auto element_vb : fes->GetDualShapeNodes(vb))
        {
          shared_ptr<CoefficientFunction> dual_trial;
          if (dual -> Dimension() == 1)
            { dual_trial = dual * trial; }
          else
            { dual_trial = InnerProduct(dual, trial); }
          auto bfi = make_shared<SymbolicBilinearFormIntegrator> (dual_trial, vb, element_vb);
	  bfi->SetSimdEvaluate(false);  // dual are not simded, yet
	  bli.Append(bfi);
	  if (auto block_bfi = dynamic_pointer_cast<BlockBilinearFormIntegrator> (bfi)) {
	    auto sbfi = block_bfi->BlockPtr();
	    sbfi->SetSimdEvaluate(false);
	    single_bli.Append(sbfi);
	  }
	  else
	    { single_bli.Append(bfi); }
	}
    }


    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate_impl (const MIR & ir, BareSliceMatrix<T,ORD> values) const
    {
      // #ifdef FIRSTDRAFT
      LocalHeapMem<100000> lh("interpolate");

      // static Timer t("interpolate");
      // RegionTracer reg(TaskManager::GetThreadId(), t);    

      const ElementTransformation & trafo = ir.GetTransformation();
      // const MeshAccess & ma = *static_cast<const MeshAccess*> (trafo.GetMesh());
      ElementId ei = trafo.GetElementId();
      auto & fel = fes->GetFE(ei, lh);
      // int dim   = fes->GetDimension();
      int dim = Dimension();

      
      // cout << " eval for ei " << ei << endl;
      // cout << " ndof = " << fel.GetNDof() << endl;

      // if (dim != 1)
      // { throw Exception("Dim != 1 porbably does not work (yet)"); }

      /** func * dual_shape **/
      FlatVector<T> elflux(fel.GetNDof(), lh);

      FlatVector<T> elfluxadd(fel.GetNDof(), lh);  elflux = 0; // non-SIMD version
      for (auto el_vb : fes->GetDualShapeNodes(trafo.VB()))
	{
          if (el_vb == VOL)
            {
              IntegrationRule ir(fel.ElementType(), 2*fel.Order()+bonus_intorder);
              auto & mir = trafo(ir, lh);
              FlatMatrix<T> mflux(ir.Size(), dim, lh);
              func->Evaluate (mir, mflux);
              for (size_t j : Range(mir))
                mflux.Row(j) *= mir[j].GetWeight();
              dual_diffop -> ApplyTrans (fel, mir, mflux, elfluxadd, lh);
              elflux += elfluxadd;
            }
          else
            {
              Facet2ElementTrafo f2el (fel.ElementType(), el_vb);
              for (int locfnr : Range(f2el.GetNFacets()))
                {
                  // SIMD does not work yet
                  // SIMD_IntegrationRule irfacet(f2el.FacetType(locfnr), 2 * fel.Order());
                  IntegrationRule irfacet(f2el.FacetType(locfnr), 2*fel.Order()+bonus_intorder);
                  auto & irvol = f2el(locfnr, irfacet, lh);
                  auto & mir = trafo(irvol, lh);
                  mir.ComputeNormalsAndMeasure(fel.ElementType(), locfnr);
                  
                  // FlatMatrix<T,ORD> mflux(dim, irfacet.Size(), lh);
                  // func->Evaluate (mir, mflux);
                  // for (size_t j : Range(mir))
                  // 	mflux.Col(j) *= mir[j].GetWeight();
                  // SIMD only
                  // dual_diffop->AddTrans (fel, mir, mflux, elflux);
                  // NON-simd version
                  
                  FlatMatrix<T> mflux(irfacet.Size(), dim, lh);
                  func->Evaluate (mir, mflux);
                  for (size_t j : Range(mir))
                    mflux.Row(j) *= mir[j].GetWeight();
                  // cout << "mflux = " << mflux << endl;
                  dual_diffop -> ApplyTrans (fel, mir, mflux, elfluxadd, lh);
                  // cout << " elfluxadd = " << endl << elfluxadd << endl;
                  elflux += elfluxadd;
                }
            }
        }
      
      /** Calc Element Matrix - shape * dual_shape **/
      FlatMatrix<double> elmat(fel.GetNDof(), lh);
      elmat = 0.0;
      bool symmetric_so_far = false;

      auto & nonconst_trafo = const_cast<ElementTransformation&>(trafo);
      auto saveud = nonconst_trafo.userdata;
      for (auto sbfi : single_bli)
        sbfi->CalcElementMatrixAdd (fel, trafo, elmat, symmetric_so_far, lh);
      nonconst_trafo.userdata = saveud;

      /** Invert Element Matrix and Solve for RHS **/
      CalcInverse(elmat); // Not Symmetric !
      

      // cout << "interpolation elmat = " << endl << elmat << endl;
      
      /** Calc coefficients of Interpolation **/
      FlatVector<double> coeffs(fel.GetNDof(), lh);
      coeffs = elmat * elflux;

      // cout << " coeffs: " << endl << coeffs << endl;

      // func->Evaluate(ir, values);
      // cout << " un-interp values: " << endl << values.AddSize(Dimension(), ir.Size()) << endl;
        
      if constexpr (ORD==ColMajor) 
        fes->GetEvaluator(vb)->Apply(fel, ir, coeffs, Trans(values), lh);
      else 
	fes->GetEvaluator(vb)->Apply(fel, ir, coeffs, values, lh);

      // cout << " values: " << endl << values.AddSize(Dimension(), ir.Size()) << endl;
    }


    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
    {
      if constexpr (is_same<MIR, SIMD_BaseMappedIntegrationRule>::value)
                     throw ExceptionNOSIMD ("no simd in InterpolateCF");
                    
      if constexpr(std::is_same<T, double>::value) {
	  // if constexpr(ORD == RowMajor) {
	      T_Evaluate_impl (ir, values);
	    // }
	  // else
	    // { throw Exception("Col-major does not compile (yet)"); }
	}
      else
	{ throw Exception("InterpolateCF::T_Evaluate only for double!"); }
      // func->Evaluate(ir, values);
    }


    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
		     FlatArray<BareSliceMatrix<T,ORD>> input,
		     BareSliceMatrix<T,ORD> values) const
    {
      T_Evaluate (ir, values);
    }

    Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ func }); }

    void PrintReport (ostream & ost) const override
    {
      ost << "InterpolationCF(";
      func->PrintReport(ost);
      ost << ")";
    }

    string GetDescription() const override
    {
      return "InterpolationCF";
    }

    void NonZeroPattern (const class ProxyUserData & ud,
			 FlatVector<AutoDiffDiff<1,bool>> nonzero) const override
    {
      func->NonZeroPattern(ud, nonzero);
    }

    void NonZeroPattern (const class ProxyUserData & ud,
			 FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
			 FlatVector<AutoDiffDiff<1,bool>> values) const override
    {
      func->NonZeroPattern(ud, input, values);
    }

    void TraverseTree (const function<void(CoefficientFunction&)> & func_) override
    {
      func->TraverseTree(func_);
      func_(*this);
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override
    {
      if (this == var) return dir;
      return InterpolateCF(func->Diff(var, dir),fes);
    }
  };
  


  class InterpolateDiffOp : public DifferentialOperator
  {
    shared_ptr<CoefficientFunction> func;
    shared_ptr<FESpace> fes;
    shared_ptr<FESpace> fes_func;
    int bonus_intorder;


    Array<shared_ptr<BilinearFormIntegrator>> single_bli;
    Array<shared_ptr<BilinearFormIntegrator>> m3_bli;
    // shared_ptr<CoefficientFunction> dual_diffop;
    shared_ptr<DifferentialOperator> dual_diffop;
    VorB vb;

    bool testfunction;

    shared_ptr<DifferentialOperator> diffop; // the final evaluation
    
  public:
    InterpolateDiffOp (shared_ptr<CoefficientFunction> afunc,
                       shared_ptr<FESpace> afes,
                       shared_ptr<DifferentialOperator> adiffop,
                       int abonus_intorder,
                       bool atestfunction)
      : DifferentialOperator (adiffop->Dim(), 1, VOL, 0), 
        func(afunc), fes(afes), bonus_intorder(abonus_intorder), testfunction(atestfunction), diffop(adiffop)
    { 
      // copied from Set (dualshapes)
      vb = VOL;    // for the moment only 

      /** Trial-Proxy **/
      auto single_evaluator =  fes->GetEvaluator(vb);
      if (dynamic_pointer_cast<BlockDifferentialOperator>(single_evaluator))
        single_evaluator = dynamic_pointer_cast<BlockDifferentialOperator>(single_evaluator)->BaseDiffOp();
      auto trial = make_shared<ProxyFunction>(fes, false, false, single_evaluator,
                                              nullptr, nullptr, nullptr, nullptr, nullptr);

      /** Test-Proxy (dual) **/
      auto dual_evaluator = fes->GetAdditionalEvaluators()["dual"];
      for (VorB avb = VOL; avb < vb; avb++) {
        dual_evaluator = dual_evaluator->GetTrace();
        if ( dual_evaluator == nullptr )
          { throw Exception(fes->GetClassName() + string(" has no dual trace operator for vb = ") + \
                            to_string(avb) + string(" -> ") + to_string(avb + 1) + string("!")); }
      }
      if (dynamic_pointer_cast<BlockDifferentialOperator>(dual_evaluator))
        dual_evaluator = dynamic_pointer_cast<BlockDifferentialOperator>(dual_evaluator)->BaseDiffOp();
      auto dual = make_shared<ProxyFunction>(fes, true, false, dual_evaluator,
                                             nullptr, nullptr, nullptr, nullptr, nullptr);
      auto dual_tr = make_shared<ProxyFunction>(fes, false, false, dual_evaluator,
                                             nullptr, nullptr, nullptr, nullptr, nullptr);

      dual_diffop = dual_evaluator;


      //Get FESpace of proxy to interpolate
      func->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          if (auto proxy = dynamic_cast<ProxyFunction*> (&nodecf))
              fes_func = proxy->GetFESpace();
        });

     

      for (auto element_vb : fes->GetDualShapeNodes(vb))
        {
          shared_ptr<CoefficientFunction> dual_trial;
          shared_ptr<CoefficientFunction> proxy_cf;
          if (dual -> Dimension() == 1)
            {
              dual_trial = dual * trial;
              proxy_cf = atestfunction ? (dual_tr * afunc) : (dual * afunc) ;
            }
          else
            {
              dual_trial = InnerProduct(dual, trial);
              proxy_cf = atestfunction ? InnerProduct(dual_tr, afunc) : InnerProduct(dual, afunc);
            }
          auto bfi = make_shared<SymbolicBilinearFormIntegrator> (dual_trial, vb, element_vb);
          auto bfi2 = make_shared<SymbolicBilinearFormIntegrator> (proxy_cf, vb, element_vb);
	  bfi->SetSimdEvaluate(true);
	  bfi2->SetSimdEvaluate(true);
	  if (auto block_bfi = dynamic_pointer_cast<BlockBilinearFormIntegrator> (bfi))
            {
              auto sbfi = block_bfi->BlockPtr();
              sbfi->SetSimdEvaluate(true);
              single_bli.Append(sbfi);
            }
	  else
	    single_bli.Append(bfi);
          if (auto block_bfi = dynamic_pointer_cast<BlockBilinearFormIntegrator> (bfi2))
            {
              auto sbfi = block_bfi->BlockPtr();
              sbfi->SetSimdEvaluate(true);
              m3_bli.Append(sbfi);
            }
	  else
	    m3_bli.Append(bfi2);
	}
    }

    void
    CalcMatrix (const FiniteElement & inner_fel,
		const BaseMappedIntegrationRule & mir,
		SliceMatrix<double,ColMajor> mat,   
		LocalHeap & lh) const override
    {
      static Timer t1("interpolateDiffOp, CalcMat");
      static Timer tm2("interpolateDiffOp, CalcMat m2");
      static Timer t23("interpolateDiffOp, mult 23");
      static Timer t23t("interpolateDiffOp, mult 23t");
      RegionTracer reg(TaskManager::GetThreadId(), t1);    

      HeapReset hr(lh);
      
      const ElementTransformation & trafo = mir.GetTransformation();
      ElementId ei = trafo.GetElementId();
      auto & interpol_fel = fes->GetFE(ei, lh);
      // int dim = func->Dimension();

      /** Calc Element Matrix - shape * dual_shape **/
      FlatMatrix<double> elmat(interpol_fel.GetNDof(), lh);
      elmat = 0.0;
      bool symmetric_so_far = false;

      // auto saveud = nonconst_trafo.userdata;
      {
      RegionTracer reg(TaskManager::GetThreadId(), tm2);          
      for (auto & sbfi : single_bli)
        sbfi->CalcElementMatrixAdd (interpol_fel, trafo, elmat, symmetric_so_far, lh);

      CalcInverse(elmat); 
      }
      /** func * dual_shape **/
      // FlatVector<> elfluxadd(interpol_fel.GetNDof(), lh);
      size_t nshape = inner_fel.GetNDof();
      
      // FlatVector<> vtrialtest(nshape, lh);
      auto save_ud = trafo.PushUserData();
      /*
      ProxyUserData myud;
      myud.trial_elvec = &vtrialtest;
      myud.test_elvec = &vtrialtest;
      myud.fel = &inner_fel;
      myud.lh = &lh;

      auto & nonconst_trafo = const_cast<ElementTransformation&>(trafo);
      nonconst_trafo.userdata = &myud;
      */

      // auto & func_fel = fes_func->GetFE(ei, lh);
      auto & func_fel = inner_fel;
      
      MixedFiniteElement mfe = (testfunction)
        ? MixedFiniteElement (interpol_fel, func_fel) 
        : MixedFiniteElement (func_fel, interpol_fel); 

      FlatMatrix<> m2m3 (elmat.Height(), nshape, lh);

      if (testfunction)
        {
          FlatMatrix<> m3T(nshape, interpol_fel.GetNDof(), lh);
          m3T = 0.0;
          for (auto & sbfi : m3_bli)
            sbfi->CalcElementMatrixAdd (mfe, trafo, m3T, symmetric_so_far, lh);
          RegionTracer reg(TaskManager::GetThreadId(), t23t);              
          m2m3 = elmat * Trans(m3T);
        }
      else
        {
          FlatMatrix<> m3(interpol_fel.GetNDof(), nshape, lh);
          m3 = 0.0;
          for (auto & sbfi : m3_bli)
            sbfi->CalcElementMatrixAdd (mfe, trafo, m3, symmetric_so_far, lh);
          RegionTracer reg(TaskManager::GetThreadId(), t23);                        
          m2m3 = elmat * m3;
        }
      
      FlatMatrix<double, ColMajor> m1(mat.Height(), interpol_fel.GetNDof(), lh);
      diffop->CalcMatrix(interpol_fel, mir, m1, lh);
      mat = m1*m2m3;
    }
  };
    


  InterpolateProxy :: InterpolateProxy (shared_ptr<CoefficientFunction> afunc,
                                        shared_ptr<FESpace> aspace,
                                        bool atestfunction,
                                        shared_ptr<DifferentialOperator> diffop,
                                        int abonus_intorder)
    : ProxyFunction(aspace, atestfunction, false,
                    make_shared<InterpolateDiffOp> (afunc, aspace, diffop, abonus_intorder,atestfunction), nullptr, nullptr,
                    nullptr, nullptr, nullptr),
      func(afunc), space(aspace), testfunction(atestfunction),
      final_diffop(diffop),
      bonus_intorder(abonus_intorder)
  { 
    this->SetDimensions (diffop->Dimensions());
  } 


  shared_ptr<ProxyFunction>
  InterpolateProxy :: GetAdditionalProxy (string name) const
  {
    shared_ptr<DifferentialOperator> new_diffop;
    new_diffop = space->GetFluxEvaluator();
    if (!new_diffop || new_diffop->Name()!=name)
      {
        auto evaluators = fes->GetAdditionalEvaluators();
        if (evaluators.Used(name))
          new_diffop = evaluators[name];
      }
    
    return make_shared<InterpolateProxy> (func, space, testfunction, new_diffop, bonus_intorder);
  }

  
  
  shared_ptr<CoefficientFunction> InterpolateCF (shared_ptr<CoefficientFunction> func, shared_ptr<FESpace> space,
                                                 int bonus_intorder)
  {
    func->PrintReport(cout);

    if (func->GetDescription() == "ZeroCF")
      return func;

    bool has_trial = false, has_test = false;

    func->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          
          if (auto proxy = dynamic_cast<ProxyFunction*> (&nodecf))
            {
              if (proxy->IsTestFunction())
                has_test = true;
              else
                has_trial = true;
            }
        });

    if (has_trial != has_test)
      return make_shared<InterpolateProxy> (func, space, has_test, space->GetEvaluator(VOL), bonus_intorder);

    return make_shared<InterpolationCoefficientFunction> (func, space, bonus_intorder);
  }

}
