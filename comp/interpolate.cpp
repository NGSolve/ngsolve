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

  class FECoefficientFunction : public T_CoefficientFunction<FECoefficientFunction>
  {
    shared_ptr<DifferentialOperator> diffop;
    // FiniteElement * fe;
    // FlatVector<> * elvec;
    Array<FiniteElement*> fes;
    Array<FlatVector<>*> elvecs;
  public:

    FECoefficientFunction (shared_ptr<DifferentialOperator> adiffop)
      : T_CoefficientFunction<FECoefficientFunction>(adiffop->Dim(), false),
      diffop(adiffop), fes(TaskManager::GetMaxThreads()), elvecs(TaskManager::GetMaxThreads())
    { }

    void Set (FiniteElement * afe, FlatVector<> * aelvec)
    {
      auto tid = TaskManager::GetThreadId();
      fes[tid] = afe;
      elvecs[tid] = aelvec;
    }
    
    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir, BareSliceMatrix<T,ORD> values) const
    {
      // static Timer t(string("FECF - T_Evaluate")+typeid(T).name());
      // RegionTracer reg(TaskManager::GetThreadId(), t);
      
      LocalHeapMem<10000> lh("fecoef::eval");
      auto tid = TaskManager::GetThreadId();
      [[maybe_unused]] auto fe = fes[tid];
      [[maybe_unused]] auto elvec = elvecs[tid];

      if constexpr (is_same<MIR, BaseMappedIntegrationRule>::value &&
                    is_same<T,double>::value)
                     {
                       diffop->Apply(*fe, mir, *elvec, Trans(values), lh);
                     }
      if constexpr (is_same<MIR, BaseMappedIntegrationRule>::value &&
                    is_same<T,AutoDiffDiff<1,double>>::value)
                     {
                       Matrix<> tmp(mir.Size(), Dimension());
                       diffop->Apply(*fe, mir, *elvec, tmp, lh);
                       values.AddSize(mir.Size(), Dimension()) = Trans(tmp);
                     }
      else if constexpr (is_same<MIR, SIMD_BaseMappedIntegrationRule>::value &&
                         is_same<T,SIMD<double>>::value)
                     {
                       diffop->Apply(*fe, mir, *elvec, values);
                     }
      else if constexpr (is_same<MIR, SIMD_BaseMappedIntegrationRule>::value &&
                         is_same<T,AutoDiffDiff<1,SIMD<double>>>::value)
                     {
                       /*
                       BareSliceMatrix<SIMD<double>> hvalues(3*values.Dist(), &values(0).Value(),
                                                             DummySize(Dimension(), mir.Size()));
                       */
                       BareSliceMatrix<SIMD<double>> hvalues(Dimension(), mir.Size(), 3*values.Dist(), &values(0).Value());
                       
                       // Evaluate (ir, hvalues);
                       diffop->Apply(*fe, mir, *elvec, hvalues);                       
                       for (size_t i = 0; i < Dimension(); i++)
                         for (size_t j = mir.Size(); j-- > 0; )
                           values(i,j) = hvalues(i,j);

                       /*
                       Matrix<SIMD<double>> tmp(Dimension(), mir.Size());
                       diffop->Apply(*fe, mir, *elvec, tmp);
                       values.AddSize(Dimension(), mir.Size()) = tmp;
                       */
                     }
      else
        {
          cout << "FECF, unhandled type: " << typeid(T).name() << endl;
        }
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
		     FlatArray<BareSliceMatrix<T,ORD>> input,
		     BareSliceMatrix<T,ORD> values) const
    {
      T_Evaluate (ir, values);
    }
  };


  
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
    optional<string> opname;
    
  public:
    InterpolationCoefficientFunction (shared_ptr<CoefficientFunction> f, shared_ptr<FESpace> afes,
                                      int abonus_intorder,
                                      optional<string> aopname)
      : T_CoefficientFunction<InterpolationCoefficientFunction>(aopname.has_value() ? afes->GetAdditionalEvaluators()[*aopname]->Dim() : f->Dimension(), f->IsComplex()),
        func(f), fes(afes), bonus_intorder(abonus_intorder), opname(aopname)
    {
      if(opname.has_value())
        this->SetDimensions(fes->GetAdditionalEvaluators()[*opname]->Dimensions());
      else
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
      LocalHeapMem<2000000> lh("interpolate");

      // static Timer t("interpolate");
      // RegionTracer reg(TaskManager::GetThreadId(), t);    

      const ElementTransformation & trafo = ir.GetTransformation();
      // const MeshAccess & ma = *static_cast<const MeshAccess*> (trafo.GetMesh());
      ElementId ei = trafo.GetElementId();
      auto & fel = fes->GetFE(ei, lh);
      // int dim   = fes->GetDimension();
      int dim = func->Dimension();

      
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

      if constexpr (ORD == ColMajor)
        {
          if(opname.has_value())
            fes->GetAdditionalEvaluators()[*opname]->Apply(fel, ir, coeffs, Trans(values), lh);
          else
            fes->GetEvaluator(vb)->Apply(fel, ir, coeffs, Trans(values), lh);
        }
      else
        {
          if(opname.has_value())
            fes->GetAdditionalEvaluators()[*opname]->Apply(fel, ir, coeffs, values, lh);
          fes->GetEvaluator(vb)->Apply(fel, ir, coeffs, values, lh);
        }

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
    shared_ptr<DifferentialOperator> dual_diffop;

    shared_ptr<FECoefficientFunction> dual_fecf;
    Array<shared_ptr<SymbolicEnergy>> hessian_bfi;

    bool testfunction;

    shared_ptr<DifferentialOperator> diffop; // the final evaluation
    
  public:
    InterpolateDiffOp (shared_ptr<CoefficientFunction> afunc,
                       shared_ptr<FESpace> afes,
                       shared_ptr<DifferentialOperator> adiffop,
                       int abonus_intorder,
                       bool atestfunction,
                       VorB avb=VOL)
      : DifferentialOperator (adiffop->Dim(), 1, avb, 0), 
        func(afunc), fes(afes), bonus_intorder(abonus_intorder), testfunction(atestfunction), diffop(adiffop)
    { 
      // copied from Set (dualshapes)

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


      dual_fecf = make_shared<FECoefficientFunction> (dual_evaluator);
      
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
          bfi->SetBonusIntegrationOrder(bonus_intorder);
          bfi2->SetBonusIntegrationOrder(bonus_intorder);
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

          hessian_bfi += make_shared<SymbolicEnergy> (InnerProduct(dual_fecf, func),
                                                      vb, element_vb);
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

      size_t nshape = inner_fel.GetNDof();
      FlatMatrix<> m2m3 (elmat.Height(), nshape, lh);

      // auto saveud = nonconst_trafo.userdata;
      try
        {
          RegionTracer reg(TaskManager::GetThreadId(), tm2);          
          for (auto & sbfi : single_bli)
            sbfi->CalcElementMatrixAdd (interpol_fel, trafo, elmat, symmetric_so_far, lh);

          CalcInverse(elmat);
          
          MixedFiniteElement mfe = (testfunction)
            ? MixedFiniteElement (interpol_fel, inner_fel) 
            : MixedFiniteElement (inner_fel, interpol_fel); 
          
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
        }
      catch (const ExceptionNOSIMD& e)
        {
          cout << IM(6) << e.What() << endl
               << "switching to scalar evaluation" << endl;
          for (auto & sbfi : single_bli)
            sbfi->SetSimdEvaluate(false);
          for (auto & sbfi : m3_bli)
            sbfi->SetSimdEvaluate(false);
          CalcMatrix (inner_fel, mir, mat, lh);
          return;
        }
       
      
      FlatMatrix<double, ColMajor> m1(mat.Height(), interpol_fel.GetNDof(), lh);
      diffop->CalcMatrix(interpol_fel, mir, m1, lh);
      mat = m1*m2m3;
    }


    
    
    void CalcLinearizedMatrix (const FiniteElement & inner_fel,
                               const BaseMappedIntegrationRule & mir,
                               BareSliceVector<double> x,
                               SliceMatrix<double,ColMajor> mat,   
                               LocalHeap & lh) const override
    {
      static Timer t("CAlcLinearizedBMatrix");
      RegionTracer reg(TaskManager::GetThreadId(), t);

      
      HeapReset hr(lh);
      
      const ElementTransformation & trafo = mir.GetTransformation();
      ElementId ei = trafo.GetElementId();
      auto & interpol_fel = fes->GetFE(ei, lh);

      FlatMatrix<double> elmat(interpol_fel.GetNDof(), lh);
      elmat = 0.0;
      bool symmetric_so_far = false;

      try
        {
          for (auto & sbfi : single_bli)
            sbfi->CalcElementMatrixAdd (interpol_fel, trafo, elmat, symmetric_so_far, lh);
        }
      catch (const ExceptionNOSIMD& e)
        {
          cout << IM(6) << e.What() << endl
               << "switching to scalar evaluation" << endl;
          for (auto & sbfi : single_bli)
            sbfi->SetSimdEvaluate(false);
          for (auto & sbfi : m3_bli)
            sbfi->SetSimdEvaluate(false);
          CalcLinearizedMatrix (inner_fel, mir, x, mat, lh);
          return;
        }

      CalcInverse(elmat); 
      size_t nshape = inner_fel.GetNDof();
      
      auto save_ud = trafo.PushUserData();

      auto & func_fel = inner_fel;
      
      MixedFiniteElement mfe = (testfunction)
        ? MixedFiniteElement (interpol_fel, func_fel) 
        : MixedFiniteElement (func_fel, interpol_fel); 

      FlatMatrix<> m2m3 (elmat.Height(), nshape, lh);

      if (testfunction)
        {
          FlatMatrix<> m3T(nshape, interpol_fel.GetNDof(), lh);
          FlatMatrix<> m3Ti(nshape, interpol_fel.GetNDof(), lh);
          FlatVector<> elvec(nshape, lh);
          elvec = x;
          m3T = 0.0;
          for (auto & sbfi : m3_bli)
            {
              sbfi->CalcLinearizedElementMatrix (mfe, trafo, elvec, m3Ti, lh);
              m3T += m3Ti;
            }
          m2m3 = elmat * Trans(m3T);
        }
      else
        {
          FlatMatrix<> m3(interpol_fel.GetNDof(), nshape, lh);
          FlatMatrix<> m3i(interpol_fel.GetNDof(), nshape, lh);
          FlatVector<> elvec(nshape, lh);
          elvec = x;
          m3 = 0.0;
          for (auto & sbfi : m3_bli)
            {
              sbfi->CalcLinearizedElementMatrix (mfe, trafo, elvec, m3i, lh);
              m3 += m3i;
            }
          m2m3 = elmat * m3;
        }
      
      FlatMatrix<double, ColMajor> m1(mat.Height(), interpol_fel.GetNDof(), lh);
      diffop->CalcMatrix(interpol_fel, mir, m1, lh);
      mat = m1*m2m3;
    }


    void Apply (const FiniteElement & inner_fel,
                const BaseMappedIntegrationRule & mir,
                BareSliceVector<double> x, 
                BareSliceMatrix<double> flux,
                LocalHeap & lh) const override
    {
      HeapReset hr(lh);

      const ElementTransformation & trafo = mir.GetTransformation();
      ElementId ei = trafo.GetElementId();
      auto & interpol_fel = fes->GetFE(ei, lh);

      FlatMatrix<double> elmat(interpol_fel.GetNDof(), lh);
      elmat = 0.0;
      bool symmetric_so_far = false;

      try
        {
          for (auto & sbfi : single_bli)
            sbfi->CalcElementMatrixAdd (interpol_fel, trafo, elmat, symmetric_so_far, lh);
        }
      catch (const ExceptionNOSIMD& e)
        {
          cout << IM(6) << e.What() << endl
               << "switching to scalar evaluation" << endl;
          for (auto & sbfi : single_bli)
            sbfi->SetSimdEvaluate(false);
          for (auto & sbfi : m3_bli)
            sbfi->SetSimdEvaluate(false);
          Apply (inner_fel, mir, x, flux, lh);
          return;          
        }
      
      CalcInverse(elmat);
      
      auto save_ud = trafo.PushUserData();

      MixedFiniteElement mfe = (testfunction)
        ? MixedFiniteElement (interpol_fel, inner_fel) 
        : MixedFiniteElement (inner_fel, interpol_fel); 

      if (testfunction)
        throw Exception("ApplyInterpolation only makes sense for trialfunctions");

      FlatVector<> rhs(interpol_fel.GetNDof(), lh);
      FlatVector<> rhsi(interpol_fel.GetNDof(), lh);
      rhs = 0;
      FlatVector<> fvx(inner_fel.GetNDof(), lh);
      fvx = x;
      for (auto & sbfi : m3_bli)
        {
          sbfi->ApplyElementMatrix (mfe, trafo, fvx, rhsi, nullptr, lh);
          rhs += rhsi;
        }

      rhsi = elmat * rhs;
      diffop->Apply(interpol_fel, mir, rhsi, flux, lh);
    }


    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		FlatMatrix<double> flux,
		BareSliceVector<double> x, 
		LocalHeap & lh) const override
    {
      HeapReset hr(lh);
      FlatMatrix<double,ColMajor> mat(flux.Height()*flux.Width(), fel.GetNDof(), lh);
      CalcMatrix (fel, mir, mat, lh);
      x.Range(0,fel.GetNDof()) = Trans(mat)*flux.AsVector();
    }

    
    NGS_DLL_HEADER virtual void
    ApplyLinearizedTrans (const FiniteElement & fel,
                          const BaseMappedIntegrationRule & mir,
                          SliceVector<double> elveclin,
                          FlatMatrix<double> flux,
                          BareSliceVector<double> x, 
                          LocalHeap & lh) const override
    {
      HeapReset hr(lh);
      FlatMatrix<double,ColMajor> mat(flux.Height()*flux.Width(), fel.GetNDof(), lh);
      CalcLinearizedMatrix (fel, mir, elveclin, mat, lh);
      x.Range(0,fel.GetNDof()) = Trans(mat)*flux.AsVector();
    }

    bool IsNonlinear() const override
    {
      return true;
    }

    // second derivative of \sum_ipt wprime * B(u) 
    void CalcHessianAdd (const FiniteElement & inner_fel,
                         const BaseMappedIntegrationRule & mir,
                         SliceMatrix<> wprime,
                         BareSliceVector<> elvecu,
                         SliceMatrix<> hessian,   
                         LocalHeap & lh) const override
    {
      // a first simple implementation by numerical differentiation ....
      static Timer t("interpolateDiffOp, Hessian");
      RegionTracer reg(TaskManager::GetThreadId(), t);
      HeapReset hr(lh);

      /*
      size_t ndof = inner_fel.GetNDof();
      double eps = 1e-6;
      FlatVector<> wprimevec(wprime.Height()*wprime.Width(), lh);
      for (int i = 0, ii = 0; i < wprime.Height(); i++)
        for (int j = 0; j < wprime.Width(); j++, ii++)
          wprimevec(ii) = wprime(i,j) * mir[i].GetWeight();
      
      FlatMatrix<double,ColMajor> bmatl(mir.Size()*diffop->Dim(), ndof, lh);
      FlatMatrix<double,ColMajor> bmatr(mir.Size()*diffop->Dim(), ndof, lh);
      FlatMatrix<double,ColMajor> dbmat(mir.Size()*diffop->Dim(), ndof, lh);
                         
      for (size_t i = 0; i < ndof; i++)
        {
          FlatVector<> elvecur(ndof, lh), elvecul(ndof, lh);
          elvecur = elvecu;
          elvecul = elvecu;
          elvecur(i) += eps;
          elvecul(i) -= eps;
          CalcLinearizedMatrix(inner_fel, mir, elvecul, bmatl, lh);
          CalcLinearizedMatrix(inner_fel, mir, elvecur, bmatr, lh);
          dbmat = 1/(2*eps) * (bmatr-bmatl);
          hessian.Row(i) += Trans(dbmat) * wprimevec;
        }
        // cout << "hessian 1 = " << endl << hessian << endl;
      */
      
      const ElementTransformation & trafo = mir.GetTransformation();
      ElementId ei = trafo.GetElementId();
      auto & interpol_fel = fes->GetFE(ei, lh);

      FlatMatrix<double> elmat(interpol_fel.GetNDof(), lh);
      elmat = 0.0;
      bool symmetric_so_far = false;
      try
        {
          for (auto & sbfi : single_bli)
            sbfi->CalcElementMatrixAdd (interpol_fel, trafo, elmat, symmetric_so_far, lh);
        }
      catch (const ExceptionNOSIMD& e)
        {
          cout << IM(6) << e.What() << endl
               << "switching to scalar evaluation" << endl;
          for (auto & sbfi : single_bli)
            sbfi->SetSimdEvaluate(false);
          for (auto & sbfi : m3_bli)
            sbfi->SetSimdEvaluate(false);
          CalcHessianAdd (inner_fel, mir, wprime, elvecu, hessian, lh);
          return;          
        }
      
      CalcInverse(elmat);

      FlatMatrix<> hessian2(inner_fel.GetNDof(), lh);
      FlatVector<> rhs(interpol_fel.GetNDof(), lh);
      FlatVector<> v2(interpol_fel.GetNDof(), lh);
      
      FlatMatrix wprimeflat = wprime | lh;
      for (int i = 0; i < wprime.Height(); i++)
        wprimeflat.Row(i) *= mir[i].GetWeight();
      
      diffop->ApplyTrans(interpol_fel, mir, wprimeflat, rhs, lh);
      v2 = Trans(elmat) * rhs;

      dual_fecf->Set (&interpol_fel, &v2);
      FlatVector<> fvelvecu(inner_fel.GetNDof(), lh);
      fvelvecu = elvecu;
      for (auto & hbfi : hessian_bfi)
        {
          hbfi->CalcLinearizedElementMatrix(inner_fel, mir.GetTransformation(), fvelvecu, hessian2, lh);
          hessian += hessian2;
        }
      // cout << "hessian2 = " << endl << hessian2 << endl;
    } 
  };
    

  shared_ptr<FESpace> FindProxySpace(shared_ptr<CoefficientFunction> func)
  {
    shared_ptr<FESpace> space;
    
    func->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          if (auto proxy = dynamic_cast<ProxyFunction*> (&nodecf))
            space = proxy->GetFESpace();
        } );
    return space;
  }

  InterpolateProxy :: InterpolateProxy (shared_ptr<CoefficientFunction> afunc,
                                        shared_ptr<FESpace> aspace,
                                        bool atestfunction,
                                        shared_ptr<DifferentialOperator> diffop,
                                        int abonus_intorder,
                                        VorB avb)
    : ProxyFunction( /* aspace */
                    FindProxySpace(afunc), atestfunction, false,
                    make_shared<InterpolateDiffOp> (afunc, aspace, diffop, abonus_intorder,atestfunction, avb), nullptr, nullptr,
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

  
  shared_ptr<CoefficientFunction> InterpolateProxy ::
  Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const 
  {
    if (this == var) return dir;
    return make_shared<InterpolateProxy> (func->Diff(var,dir), space, testfunction, final_diffop, bonus_intorder);
  }
  
  shared_ptr<CoefficientFunction> InterpolateCF (shared_ptr<CoefficientFunction> func, shared_ptr<FESpace> space,
                                                 int bonus_intorder,
                                                 optional<string> opname)
  {
    // func->PrintReport(cout);

    if (func->IsZeroCF())
      return func;

    bool has_trial = false, has_test = false;
    VorB vb = VOL;

    func->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          
          if (auto proxy = dynamic_cast<ProxyFunction*> (&nodecf))
            {
              if (proxy->IsTestFunction())
                {
                  // if e.g. H1 without Trace is combined with space where Trace is mandatory
                  if(opname.has_value())
                    vb = max(proxy->GetAdditionalEvaluator(*opname)->VB(),vb);
                  else
                    vb = max(proxy->Evaluator()->VB(),vb);
                  has_test = true;
                }
              else
                {
                  if(opname.has_value())
                    vb = max(proxy->GetAdditionalEvaluator(*opname)->VB(),vb);
                  else
                    vb = max(proxy->Evaluator()->VB(),vb);
                  has_trial = true;
                }
            }
        });

    if (has_trial != has_test)
      {
        if (opname.has_value())
          return make_shared<InterpolateProxy>(
              func, space, has_test, space->GetAdditionalEvaluators()[*opname],
              bonus_intorder, vb);
        else
          return make_shared<InterpolateProxy>(func, space, has_test,
                                               space->GetEvaluator(vb),
                                               bonus_intorder, vb);
      }
    return make_shared<InterpolationCoefficientFunction> (func, space, bonus_intorder, opname);
  }

}
