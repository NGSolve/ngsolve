#ifndef FILE_SYMBOLICINTEGRATOR
#define FILE_SYMBOLICINTEGRATOR

/*********************************************************************/
/* File:   symbolicintegrator.hpp                                    */
/* Author: Joachim Schoeberl                                         */
/* Date:   August 2015                                               */
/*********************************************************************/



namespace ngfem
{

  

class ProxyFunction : public CoefficientFunction
{
  bool testfunction; // true .. test, false .. trial
  // bool is_complex;
  bool is_other;    // neighbour element (DG)

  shared_ptr<DifferentialOperator> evaluator;
  shared_ptr<DifferentialOperator> deriv_evaluator;
  shared_ptr<DifferentialOperator> trace_evaluator;
  shared_ptr<DifferentialOperator> trace_deriv_evaluator;
  shared_ptr<DifferentialOperator> ttrace_evaluator;
  shared_ptr<DifferentialOperator> ttrace_deriv_evaluator;
  shared_ptr<ProxyFunction> deriv_proxy;
  shared_ptr<CoefficientFunction> boundary_values; // for DG - apply

  SymbolTable<shared_ptr<DifferentialOperator>> additional_diffops;
  // int dim;
public:
  NGS_DLL_HEADER ProxyFunction (bool atestfunction, bool ais_complex,
                 shared_ptr<DifferentialOperator> aevaluator, 
                 shared_ptr<DifferentialOperator> aderiv_evaluator,
                 shared_ptr<DifferentialOperator> atrace_evaluator,
                 shared_ptr<DifferentialOperator> atrace_deriv_evaluator,
		 shared_ptr<DifferentialOperator> attrace_evaluator,
		 shared_ptr<DifferentialOperator> attrace_deriv_evaluator);

  bool IsTestFunction () const { return testfunction; }
  bool IsOther() const { return is_other; }

  NGS_DLL_HEADER virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
  
  const shared_ptr<DifferentialOperator> & Evaluator() const { return evaluator; }
  const shared_ptr<DifferentialOperator> & DerivEvaluator() const { return deriv_evaluator; }
  const shared_ptr<DifferentialOperator> & TraceEvaluator() const { return trace_evaluator; }
  const shared_ptr<DifferentialOperator> & TraceDerivEvaluator() const { return trace_deriv_evaluator; }
  const shared_ptr<DifferentialOperator> & TTraceEvaluator() const { return ttrace_evaluator; }
  const shared_ptr<DifferentialOperator> & TTraceDerivEvaluator() const { return ttrace_deriv_evaluator; }

  shared_ptr<ProxyFunction> Deriv() const
  {
    return deriv_proxy;
  }

  NGS_DLL_HEADER shared_ptr<ProxyFunction> Trace() const;

  shared_ptr<ProxyFunction> Other(shared_ptr<CoefficientFunction> _boundary_values) const
  {
    auto other = make_shared<ProxyFunction> (testfunction, is_complex, evaluator, deriv_evaluator, trace_evaluator, trace_deriv_evaluator,ttrace_evaluator, ttrace_deriv_evaluator);
    other->is_other = true;
    if (other->deriv_proxy)
      other->deriv_proxy->is_other = true;
    other->boundary_values = _boundary_values;

    for (int i = 0; i < additional_diffops.Size(); i++)
      other->SetAdditionalEvaluator (additional_diffops.GetName(i), additional_diffops[i]);
    
    return other;
  }
  const shared_ptr<CoefficientFunction> & BoundaryValues() const { return boundary_values; } 

  void SetAdditionalEvaluator (string name, shared_ptr<DifferentialOperator> diffop)
  {
    additional_diffops.Set (name, diffop);
  }
  
  shared_ptr<DifferentialOperator> GetAdditionalEvaluator (string name) const
  {
    if (additional_diffops.Used(name))
      return additional_diffops[name];
    return shared_ptr<DifferentialOperator>();
  }

  SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const
  {
    return additional_diffops;
  }

  shared_ptr<ProxyFunction> GetAdditionalProxy (string name) const
  {
    if (additional_diffops.Used(name))
    {
      auto adddiffop = make_shared<ProxyFunction> (testfunction, is_complex, additional_diffops[name], nullptr, nullptr, nullptr, nullptr, nullptr);
      if (is_other)
        adddiffop->is_other = true;
      return adddiffop;
    }
    return shared_ptr<ProxyFunction>();
  }

  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    // Vector<> tmp(Dimension());
    STACK_ARRAY(double, mem, Dimension());
    FlatVector<> tmp(Dimension(), &mem[0]);
    Evaluate (ip, tmp);
    return tmp(0);
  }

  NGS_DLL_HEADER virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const;

  NGS_DLL_HEADER virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<Complex> result) const;

  NGS_DLL_HEADER virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         FlatMatrix<> result) const;

  NGS_DLL_HEADER virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         FlatMatrix<Complex> result) const;

  // virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
  // AFlatMatrix<double> values) const;

  NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         BareSliceMatrix<SIMD<double>> values) const;
  NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         BareSliceMatrix<SIMD<Complex>> values) const;

  NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const;

  NGS_DLL_HEADER virtual void EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const;

  NGS_DLL_HEADER virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const;

  NGS_DLL_HEADER virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                              AFlatMatrix<double> values, AFlatMatrix<double> deriv) const;
  
  NGS_DLL_HEADER virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                               AFlatMatrix<double> values, AFlatMatrix<double> deriv,
                               AFlatMatrix<double> dderiv) const;
  
  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                              FlatArray<AFlatMatrix<>*> input,
                              FlatArray<AFlatMatrix<>*> dinput,
                              AFlatMatrix<> result,
                              AFlatMatrix<> deriv) const
  {
    EvaluateDeriv (ir, result, deriv);
  }
  
  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                               FlatArray<AFlatMatrix<>*> input,
                               FlatArray<AFlatMatrix<>*> dinput,
                               FlatArray<AFlatMatrix<>*> ddinput,
                               AFlatMatrix<> result,
                               AFlatMatrix<> deriv,
                               AFlatMatrix<> dderiv) const
  {
    EvaluateDDeriv (ir, result, deriv, dderiv);
  }

  
  virtual bool ElementwiseConstant () const { return true; }

  NGS_DLL_HEADER virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const;
};

class ProxyUserData
{
  FlatArray<const ProxyFunction*> remember_first;
  FlatArray<FlatMatrix<double>> remember_second;
  FlatArray<FlatMatrix<SIMD<double>>> remember_asecond;

  FlatArray<const CoefficientFunction*> remember_cf_first;
  FlatArray<FlatMatrix<SIMD<double>>> remember_cf_asecond;
  FlatArray<bool> remember_cf_computed;
public:
  class ProxyFunction * testfunction = nullptr;
  int test_comp;
  class ProxyFunction * trialfunction = nullptr;
  int trial_comp;
  
  const FiniteElement * fel = nullptr;
  const FlatVector<double> * elx;
  LocalHeap * lh;

  ProxyUserData ()
    : remember_first(0,nullptr), remember_second(0,nullptr), remember_asecond(0,nullptr),
      remember_cf_first(0, nullptr), remember_cf_asecond(0,nullptr), remember_cf_computed(0, nullptr)
  { ; }
  ProxyUserData (int ntrial, int ncf, LocalHeap & lh)
    : remember_first(ntrial, lh), remember_second(ntrial, lh),
      remember_asecond(ntrial, lh),
      remember_cf_first(ncf, lh), remember_cf_asecond(ncf, lh),
      remember_cf_computed(ncf, lh)
  { remember_first = nullptr; remember_cf_first = nullptr; }

  ProxyUserData (int ntrial, LocalHeap & lh)
    : ProxyUserData (ntrial, 0, lh) { ; } 
  
  void AssignMemory (const ProxyFunction * proxy, int h, int w, LocalHeap & lh)
  {
    for (int i = 0; i < remember_first.Size(); i++)
      {
        if (remember_first[i] == nullptr)
          {
            remember_first[i] = proxy;
            new (&remember_second[i]) FlatMatrix<> (h, w, lh);
            new (&remember_asecond[i]) FlatMatrix<SIMD<double>> (w, (h+SIMD<double>::Size()-1)/SIMD<double>::Size(), lh);
            return;
          }
      }
    throw Exception ("no space for userdata - memory available");
  }

  void AssignMemory (const CoefficientFunction * cf, int h, int w, LocalHeap & lh)
  {
    for (int i = 0; i < remember_cf_first.Size(); i++)
      {
        if (remember_cf_first[i] == nullptr)
          {
            remember_cf_first[i] = cf;
            new (&remember_cf_asecond[i]) FlatMatrix<SIMD<double>> (w, (h+SIMD<double>::Size()-1)/SIMD<double>::Size(), lh);
            remember_cf_computed[i] = false;
            return;
          }
      }
    throw Exception ("no space for userdata - memory available");
  }


  bool HasMemory (const ProxyFunction * proxy) const
  {
    return remember_first.Contains(proxy);
  }
  bool HasMemory (const CoefficientFunction * cf) const
  {
    return remember_cf_first.Contains(cf);
  }
  FlatMatrix<> GetMemory (const ProxyFunction * proxy) const
  {
    return remember_second[remember_first.Pos(proxy)];
  }
  FlatMatrix<SIMD<double>> GetAMemory (const ProxyFunction * proxy) const
  {
    return remember_asecond[remember_first.Pos(proxy)];
  }
  FlatMatrix<SIMD<double>> GetAMemory (const CoefficientFunction * cf) const
  {
    return remember_cf_asecond[remember_cf_first.Pos(cf)];
  }
  bool Computed (const CoefficientFunction * cf) const
  {
    return remember_cf_computed[remember_cf_first.Pos(cf)];
  }
  void SetComputed (const CoefficientFunction * cf) const
  {
    remember_cf_computed[remember_cf_first.Pos(cf)] = true;
  }
};

  

class CompoundDifferentialOperator : public DifferentialOperator
{
  shared_ptr<DifferentialOperator> diffop;
  int comp;
public:
  CompoundDifferentialOperator (shared_ptr<DifferentialOperator> adiffop, 
                                int acomp)
    : DifferentialOperator(adiffop->Dim(), adiffop->BlockDim(),
                           adiffop->VB(), adiffop->DiffOrder()),
      diffop(adiffop), comp(acomp)
  {
    dimensions = adiffop->Dimensions();
  }
  
  virtual ~CompoundDifferentialOperator () = default;
  shared_ptr<DifferentialOperator> BaseDiffOp() const { return diffop; } 
  int Component () const { return comp; }

  virtual bool operator== (const DifferentialOperator & diffop2) const
  {
    const CompoundDifferentialOperator * do2 =
      dynamic_cast<const CompoundDifferentialOperator*> (&diffop2);
    if (do2 && do2->Component() == Component())
      return *diffop == *(do2->diffop);
    return false;
  }

  
  virtual string Name() const { return diffop->Name(); }

  virtual IntRange UsedDofs(const FiniteElement & bfel) const
  {
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = BlockDim() * fel.GetRange(comp);
    return r;
  }
  
  NGS_DLL_HEADER virtual void
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & mip,
              SliceMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const
  {
    mat = 0;
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = BlockDim() * fel.GetRange(comp);
    diffop->CalcMatrix (fel[comp], mip, mat.Cols(r), lh);
  }

  NGS_DLL_HEADER virtual void
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & mip,
              SliceMatrix<Complex,ColMajor> mat, 
              LocalHeap & lh) const
  {
    mat = 0;
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = BlockDim() * fel.GetRange(comp);
    diffop->CalcMatrix (fel[comp], mip, mat.Cols(r), lh);
  }

  NGS_DLL_HEADER virtual void
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationRule & mir,
              SliceMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const
  {
    mat = 0;
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = BlockDim() * fel.GetRange(comp);
    diffop->CalcMatrix (fel[comp], mir, mat.Cols(r), lh);
  }
  
  NGS_DLL_HEADER virtual void
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationRule & mir,
              SliceMatrix<Complex,ColMajor> mat,   
              LocalHeap & lh) const
  {
    mat = 0;
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = BlockDim() * fel.GetRange(comp);
    diffop->CalcMatrix (fel[comp], mir, mat.Cols(r), lh);
  }

  NGS_DLL_HEADER virtual void
  CalcMatrix (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & mir,
              BareSliceMatrix<SIMD<double>> mat) const
  {
    // mat = 0;   // take care: unused elements not zerod !!!!
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = Dim() * fel.GetRange(comp);
    diffop->CalcMatrix (fel[comp], mir, mat.Rows(r));
  }
  
  NGS_DLL_HEADER virtual void
  Apply (const FiniteElement & bfel,
         const BaseMappedIntegrationPoint & mip,
         FlatVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = BlockDim() * fel.GetRange(comp);
    diffop->Apply (fel[comp], mip, x.Range(r), flux, lh);
  }


  virtual void
  Apply (const FiniteElement & bfel,
         const SIMD_BaseMappedIntegrationRule & bmir,
         BareSliceVector<double> x, 
         BareSliceMatrix<SIMD<double>> flux) const
  {
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = BlockDim() * fel.GetRange(comp);
    diffop->Apply (fel[comp], bmir, x.Range(r), flux);
  }


  
  NGS_DLL_HEADER virtual void
  ApplyTrans (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              FlatVector<double> x, 
              LocalHeap & lh) const
  {
    x = 0;
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = BlockDim() * fel.GetRange(comp);
    diffop->ApplyTrans (fel[comp], mip, flux, x.Range(r), lh);
  }

  NGS_DLL_HEADER virtual void
  ApplyTrans (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<Complex> flux,
              FlatVector<Complex> x, 
              LocalHeap & lh) const
  {
    x = 0;
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = BlockDim() * fel.GetRange(comp);
    diffop->ApplyTrans (fel[comp], mip, flux, x.Range(r), lh);
  }

  NGS_DLL_HEADER virtual void
  AddTrans (const FiniteElement & bfel,
            const SIMD_BaseMappedIntegrationRule & bmir,
            BareSliceMatrix<SIMD<double>> flux,
            BareSliceVector<double> x) const
  {
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = BlockDim() * fel.GetRange(comp);
    diffop->AddTrans (fel[comp], bmir, flux, x.Range(r));
  }
  
  NGS_DLL_HEADER virtual void
  AddTrans (const FiniteElement & bfel,
            const SIMD_BaseMappedIntegrationRule & bmir,
            BareSliceMatrix<SIMD<Complex>> flux,
            BareSliceVector<Complex> x) const
  {
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = BlockDim() * fel.GetRange(comp);
    diffop->AddTrans (fel[comp], bmir, flux, x.Range(r));
  }
};





  class SymbolicLinearFormIntegrator : public LinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> proxies;
    VorB vb;
    bool element_boundary;
    mutable bool simd_evaluate = true;
    IntegrationRule ir;   // if non-empty use this integration-rule
    SIMD_IntegrationRule simd_ir;   // if non-empty use this integration-rule

  public:
    NGS_DLL_HEADER SymbolicLinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb,
                                  bool aelement_boundary);

    virtual VorB VB() const { return vb; }
    virtual string Name () const { return string ("Symbolic LFI"); }
    virtual int GetDimension() const override { return proxies[0]->Evaluator()->BlockDim(); }

    void SetIntegrationRule (const IntegrationRule & _ir);

    virtual void 
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatVector<double> elvec,
		       LocalHeap & lh) const;
      
    virtual void 
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatVector<Complex> elvec,
		       LocalHeap & lh) const;

    template <typename SCAL> 
    void T_CalcElementVector (const FiniteElement & fel,
                              const ElementTransformation & trafo, 
                              FlatVector<SCAL> elvec,
                              LocalHeap & lh) const;
  };



  class SymbolicBilinearFormIntegrator : public BilinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> trial_proxies, test_proxies;
    Array<CoefficientFunction*> gridfunction_cfs;
    Array<int> trial_cum, test_cum;   // cumulated dimension of proxies
    VorB vb;
    bool element_boundary;
    Matrix<bool> nonzeros;    // do components interact ? 
    Matrix<bool> nonzeros_proxies; // do proxies interact ?
    Matrix<bool> diagonal_proxies; // do proxies interact diagonally ?
    Matrix<bool> same_diffops; // are diffops the same ? 
    bool elementwise_constant;
    mutable bool simd_evaluate;
    IntegrationRule ir;   // if non-empty use this integration-rule
    SIMD_IntegrationRule simd_ir;   // if non-empty use this integration-rule
    int trial_difforder, test_difforder;
  public:
    NGS_DLL_HEADER SymbolicBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb,
                                    bool aelement_boundary);

    virtual VorB VB() const { return vb; }
    virtual bool IsSymmetric() const { return true; }  // correct would be: don't know
    virtual string Name () const { return string ("Symbolic BFI"); }

    NGS_DLL_HEADER virtual IntegrationRule GetIntegrationRule (const FiniteElement & fel, LocalHeap & lh) const;
    NGS_DLL_HEADER virtual SIMD_IntegrationRule Get_SIMD_IntegrationRule (const FiniteElement & fel, LocalHeap & lh) const;
    // virtual IntegrationRule GetIntegrationRuleEB (const FiniteElement & fel, int facetnr, LocalHeap & lh) const;
    // virtual SIMD_IntegrationRule Get_SIMD_IntegrationRuleEB (const FiniteElement & fel, int facetnr, LocalHeap & lh) const;
    
    void SetIntegrationRule (const IntegrationRule & _ir);

    virtual int GetDimension() const override { return trial_proxies[0]->Evaluator()->BlockDim(); }
    
    virtual void 
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const;

    virtual void 
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatMatrix<Complex> elmat,
		       LocalHeap & lh) const;    

    virtual void 
    CalcElementMatrixAdd (const FiniteElement & fel,
                          const ElementTransformation & trafo, 
                          FlatMatrix<double> elmat,
                          LocalHeap & lh) const;
    
    virtual void 
    CalcElementMatrixAdd (const FiniteElement & fel,
                          const ElementTransformation & trafo, 
                          FlatMatrix<Complex> elmat,
                          LocalHeap & lh) const;    

    
    template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
    void T_CalcElementMatrixAdd (const FiniteElement & fel,
                                 const ElementTransformation & trafo, 
                                 FlatMatrix<SCAL_RES> elmat,
                                 LocalHeap & lh) const;

    template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
    void T_CalcElementMatrixEBAdd (const FiniteElement & fel,
                                   const ElementTransformation & trafo, 
                                   FlatMatrix<SCAL_RES> elmat,
                                   LocalHeap & lh) const;
    
    virtual void 
    CalcLinearizedElementMatrix (const FiniteElement & fel,
                                 const ElementTransformation & trafo, 
				 FlatVector<double> elveclin,
                                 FlatMatrix<double> elmat,
                                 LocalHeap & lh) const;

    template <typename SCAL, typename SCAL_SHAPES>
    void T_CalcLinearizedElementMatrixEB (const FiniteElement & fel,
                                          const ElementTransformation & trafo, 
                                          FlatVector<double> elveclin,
                                          FlatMatrix<double> elmat,
                                          LocalHeap & lh) const;
    
    virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & trafo, 
			const FlatVector<double> elx, 
			FlatVector<double> ely,
			void * precomputed,
			LocalHeap & lh) const;

    template <typename SCAL, typename SCAL_SHAPES>
    void T_ApplyElementMatrixEB (const FiniteElement & fel, 
                                 const ElementTransformation & trafo, 
                                 const FlatVector<double> elx, 
                                 FlatVector<double> ely,
                                 void * precomputed,
                                 LocalHeap & lh) const;


  };



  class SymbolicFacetLinearFormIntegrator : public FacetLinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> proxies;
    Array<int> test_cum;    // cumulated dimension of proxies
    VorB vb;                // only BND supported by now
    // bool element_boundary;  /// not needed (by now ???)
    IntegrationRule ir;   // if non-empty use this integration-rule
    SIMD_IntegrationRule simd_ir;   // if non-empty use this integration-rule

  public:
    SymbolicFacetLinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb);

    virtual VorB VB() const { return vb; }
    virtual bool BoundaryForm() const { return vb == BND; }

    virtual void
    CalcFacetVector (const FiniteElement & volumefel, int LocalFacetNr,
                     const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                     const ElementTransformation & seltrans,
                     FlatVector<double> elvec,
                     LocalHeap & lh) const;
  };


  

  class SymbolicFacetBilinearFormIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> trial_proxies, test_proxies;
    Array<CoefficientFunction*> gridfunction_cfs;    
    Array<int> trial_cum, test_cum;   // cumulated dimension of proxies
    VorB vb;
    bool element_boundary;
    bool neighbor_testfunction;
    mutable bool simd_evaluate;
  public:
    NGS_DLL_HEADER SymbolicFacetBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb, bool aelement_boundary);

    virtual VorB VB() const { return vb; }
    virtual bool BoundaryForm() const { return vb == BND; }
    virtual bool IsSymmetric() const { return true; }  // correct would be: don't know
    
    virtual DGFormulation GetDGFormulation() const { return DGFormulation(neighbor_testfunction,
                                                                          element_boundary); }
    
    NGS_DLL_HEADER virtual void
    CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                     const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                     const ElementTransformation & seltrans, FlatArray<int> & SElVertices,  
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    CalcLinearizedFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                               const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                               const ElementTransformation & seltrans, FlatArray<int> & SElVertices,  
                               FlatVector<double> vec, FlatMatrix<double> elmat,
                               LocalHeap & lh) const;
    
    virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                      const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                      const FiniteElement & volumefel2, int LocalFacetNr2,
                      const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                      FlatVector<double> elx, FlatVector<double> ely,
                      LocalHeap & lh) const;

    virtual void
    CalcTraceValues (const FiniteElement & volumefel, int LocalFacetNr,
		     const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
		     FlatVector<double> & trace, FlatVector<double> elx, LocalHeap & lh) const;

    virtual void
    ApplyFromTraceValues (const FiniteElement & volumefel, int LocalFacetNr,
			  const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			  FlatVector<double> trace,
			  FlatVector<double> elx, FlatVector<double> ely, 
			  LocalHeap & lh) const;
    
    virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                      const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                      const ElementTransformation & seltrans, FlatArray<int> & SElVertices,
                      FlatVector<double> elx, FlatVector<double> ely,
                      LocalHeap & lh) const;

  };

  class SymbolicEnergy : public BilinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> cf;
    VorB vb;
    Array<ProxyFunction*> trial_proxies;
    mutable bool simd_evaluate;
    
  public:
    SymbolicEnergy (shared_ptr<CoefficientFunction> acf, VorB avb);

    virtual VorB VB() const { return vb; }
    virtual bool IsSymmetric() const { return true; } 
    virtual string Name () const { return string ("Symbolic Energy"); }
    
    virtual void 
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const
    {
      cout << "SymbolicEnergy :: CalcMatrix not implemented" << endl;
    }

    virtual void 
    CalcLinearizedElementMatrix (const FiniteElement & fel,
                                 const ElementTransformation & trafo, 
				 FlatVector<double> elveclin,
                                 FlatMatrix<double> elmat,
                                 LocalHeap & lh) const;


    virtual double Energy (const FiniteElement & fel, 
			   const ElementTransformation & trafo, 
                           FlatVector<double> elx, 
			   LocalHeap & lh) const;

    virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & trafo, 
			const FlatVector<double> elx, 
			FlatVector<double> ely,
			void * precomputed,
			LocalHeap & lh) const;
    

  };
  

}


#endif
