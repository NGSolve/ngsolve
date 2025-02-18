#ifndef FILE_SYMBOLICINTEGRATOR
#define FILE_SYMBOLICINTEGRATOR

/*********************************************************************/
/* File:   symbolicintegrator.hpp                                    */
/* Author: Joachim Schoeberl                                         */
/* Date:   August 2015                                               */
/*********************************************************************/

namespace ngcomp
{
  class FESpace;
}


#include "integrator.hpp"
#include "coefficient.hpp"


namespace ngfem
{

  class NGS_DLL_HEADER ProxyFunction : public CoefficientFunction
{
protected:
  shared_ptr<ngcomp::FESpace> fes;
  bool testfunction; // true .. test, false .. trial
  // bool is_complex;
  bool is_other;    // neighbour element (DG)
  shared_ptr<ProxyFunction> primaryproxy;   // derivatives and traces point to it
  shared_ptr<DifferentialOperator> evaluator;
  shared_ptr<DifferentialOperator> deriv_evaluator;
  shared_ptr<DifferentialOperator> trace_evaluator;
  shared_ptr<DifferentialOperator> trace_deriv_evaluator;
  shared_ptr<DifferentialOperator> ttrace_evaluator;
  shared_ptr<DifferentialOperator> ttrace_deriv_evaluator;
  weak_ptr<ProxyFunction> deriv_proxy;          // weak since we point to the primary proxy
  shared_ptr<CoefficientFunction> boundary_values; // for DG - apply

  SymbolTable<shared_ptr<DifferentialOperator>> additional_diffops;
  mutable SymbolTable<weak_ptr<ProxyFunction>> additional_proxies;

  mutable weak_ptr<ProxyFunction> dt;
  mutable shared_ptr<ProxyFunction> anti_dt;
  // int dim;
public:
  ProxyFunction (shared_ptr<ngcomp::FESpace> afes,
                                bool atestfunction, bool ais_complex,
                                shared_ptr<DifferentialOperator> aevaluator, 
                                shared_ptr<DifferentialOperator> aderiv_evaluator,
                                shared_ptr<DifferentialOperator> atrace_evaluator,
                                shared_ptr<DifferentialOperator> atrace_deriv_evaluator,
                                shared_ptr<DifferentialOperator> attrace_evaluator,
                                shared_ptr<DifferentialOperator> attrace_deriv_evaluator);

  ProxyFunction (const ProxyFunction &) = default;

  bool IsTrialFunction () const { return !testfunction; }
  bool IsTestFunction () const { return testfunction; }
  bool IsOther() const { return is_other; }

  string GetDescription () const override;

  void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override;
  
  const shared_ptr<DifferentialOperator> & Evaluator() const { return evaluator; }
  const shared_ptr<DifferentialOperator> & DerivEvaluator() const { return deriv_evaluator; }
  const shared_ptr<DifferentialOperator> & TraceEvaluator() const { return trace_evaluator; }
  const shared_ptr<DifferentialOperator> & TraceDerivEvaluator() const { return trace_deriv_evaluator; }
  const shared_ptr<DifferentialOperator> & TTraceEvaluator() const { return ttrace_evaluator; }
  const shared_ptr<DifferentialOperator> & TTraceDerivEvaluator() const { return ttrace_deriv_evaluator; }

  shared_ptr<ProxyFunction> Deriv() const;
  shared_ptr<ProxyFunction> Trace() const;
  shared_ptr<CoefficientFunction> Primary() const override
  { return primaryproxy; }
  
  shared_ptr<ProxyFunction> Dt() const;  
  shared_ptr<ProxyFunction> AntiDt() const;
  int OrderDt() const;    

  
  shared_ptr<ProxyFunction> Other(shared_ptr<CoefficientFunction> _boundary_values) const;

  const shared_ptr<CoefficientFunction> & BoundaryValues() const { return boundary_values; } 

  void SetAdditionalEvaluator (string name, shared_ptr<DifferentialOperator> diffop)
  {
    additional_diffops.Set (name, diffop);
  }
  
  virtual shared_ptr<DifferentialOperator> GetAdditionalEvaluator (string name) const
  {
    if (additional_diffops.Used(name))
      return additional_diffops[name];
    return nullptr; // shared_ptr<DifferentialOperator>();
  }

  SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const
  {
    return additional_diffops;
  }

  virtual shared_ptr<ProxyFunction> GetAdditionalProxy (string name) const;
  
  shared_ptr<CoefficientFunction> Operator (const string & name) const override;
  shared_ptr<CoefficientFunction> Operator (shared_ptr<DifferentialOperator> diffop) const override;
    
  const shared_ptr<ngcomp::FESpace> & GetFESpace() const { return fes; }
  void SetFESpace(shared_ptr<ngcomp::FESpace> fespace) { fes = fespace; }
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    // Vector<> tmp(Dimension());
    STACK_ARRAY(double, mem, Dimension());
    FlatVector<> tmp(Dimension(), &mem[0]);
    Evaluate (ip, tmp);
    return tmp(0);
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const override;

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<Complex> result) const override;

  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                                        BareSliceMatrix<> result) const override;

  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                                        BareSliceMatrix<Complex> result) const override;

  // virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
  // AFlatMatrix<double> values) const;

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         BareSliceMatrix<SIMD<double>> values) const override;
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         BareSliceMatrix<SIMD<Complex>> values) const override;

  /*
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const override;

  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const override;

  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const override;

  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                              AFlatMatrix<double> values, AFlatMatrix<double> deriv) const override;
  
  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                               AFlatMatrix<double> values, AFlatMatrix<double> deriv,
                               AFlatMatrix<double> dderiv) const override;
  */
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         FlatArray<BareSliceMatrix<SIMD<double>>> input,
                         BareSliceMatrix<SIMD<double>> values) const override
  {
    ProxyFunction::Evaluate (ir, values);
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
                         BareSliceMatrix<AutoDiff<1,double>> values) const override;
  
  virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
                         BareSliceMatrix<AutoDiffDiff<1,double>> values) const override;
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                         BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const override;

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         FlatArray<BareSliceMatrix<AutoDiff<1,SIMD<double>>>> input,
                         BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const override
  {
    ProxyFunction::Evaluate (ir, values);
  }
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                         BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const override;

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         FlatArray<BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>>> input,
                         BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const override
  {
    ProxyFunction::Evaluate (ir, values);
  }
  
  /*
  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                              FlatArray<AFlatMatrix<>*> input,
                              FlatArray<AFlatMatrix<>*> dinput,
                              AFlatMatrix<> result,
                              AFlatMatrix<> deriv) const override
  {
    EvaluateDeriv (ir, result, deriv);
  }
  
  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                               FlatArray<AFlatMatrix<>*> input,
                               FlatArray<AFlatMatrix<>*> dinput,
                               FlatArray<AFlatMatrix<>*> ddinput,
                               AFlatMatrix<> result,
                               AFlatMatrix<> deriv,
                               AFlatMatrix<> dderiv) const override
  {
    EvaluateDDeriv (ir, result, deriv, dderiv);
  }
  */
  
  // virtual bool ElementwiseConstant () const  override{ return true; }

  // the old one, to be replaced
  void NonZeroPattern (const class ProxyUserData & ud,
                                      FlatVector<bool> nonzero,
                                      FlatVector<bool> nonzero_deriv,
                                      FlatVector<bool> nonzero_dderiv) const;

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,NonZero>> values) const override;
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                               FlatVector<AutoDiffDiff<1,NonZero>> values) const override
  {
    NonZeroPattern (ud, values);
    /*
    Vector<bool> nz(values.Size()), nzd(values.Size()), nzdd(values.Size());
    NonZeroPattern (ud, nz, nzd, nzdd);
    for (size_t i = 0; i < values.Size(); i++)
      {
        values(i).Value() = nz(i);
        values(i).DValue(0) = nzd(i);
        values(i).DDValue(0) = nzdd(i);
      }
    */
  }

  
  virtual shared_ptr<CoefficientFunction>
  Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override;
};


  class SumOfIntegrals;
  class DualProxyFunction : public ProxyFunction
  {
  public:
    DualProxyFunction(const ProxyFunction & proxy)
      : ProxyFunction(proxy) { ; }

    shared_ptr<SumOfIntegrals> operator() (shared_ptr<CoefficientFunction> u) const;
  };


  
class ProxyUserData
{
  FlatArray<const ProxyFunction*> remember_first;
  FlatArray<FlatMatrix<double>> remember_second;
  FlatArray<FlatMatrix<SIMD<double>>> remember_asecond;

  FlatArray<const CoefficientFunction*> remember_cf_first;
  FlatArray<FlatMatrix<double>> remember_cf_second;  
  FlatArray<FlatMatrix<SIMD<double>>> remember_cf_asecond;
  FlatArray<bool> remember_cf_computed;
public:
  class ProxyFunction * testfunction = nullptr;
  int test_comp;
  class ProxyFunction * trialfunction = nullptr;
  int trial_comp;
  int eval_deriv = 0; // 0 .. evaluate bfi, 1 .. deriv, 2 .. second order deriv
  const FiniteElement * fel = nullptr;
  FlatArray<pair<const CoefficientFunction*, void*>> caches;

  FlatVector<double> *trial_elvec = nullptr, *test_elvec = nullptr; // for shape-wise evaluate
  LocalHeap * lh = nullptr;

  
  ProxyUserData ()
    : remember_first(0,nullptr), remember_second(0,nullptr), remember_asecond(0,nullptr),
      remember_cf_first(0, nullptr), remember_cf_second(0,nullptr),remember_cf_asecond(0,nullptr), remember_cf_computed(0, nullptr)
  { ; }
  ProxyUserData (size_t ntrial, size_t ncf, LocalHeap & lh)
    : remember_first(ntrial, lh), remember_second(ntrial, lh),
      remember_asecond(ntrial, lh),
      remember_cf_first(ncf, lh), remember_cf_second(ncf, lh),
      remember_cf_asecond(ncf, lh),      
      remember_cf_computed(ncf, lh)
  { remember_first = nullptr; remember_cf_first = nullptr; }

  ProxyUserData (int ntrial, LocalHeap & lh)
    : ProxyUserData (ntrial, 0, lh) { ; } 
  
  void AssignMemory (const ProxyFunction * proxy, size_t h, size_t w, LocalHeap & lh)
  {
    for (size_t i = 0; i < remember_first.Size(); i++)
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

  void AssignMemory (const ProxyFunction * proxy, FlatMatrix<SIMD<double>> mat)
  {
    for (size_t i = 0; i < remember_first.Size(); i++)
      {
        if (remember_first[i] == nullptr || remember_first[i] == proxy)
          {
            remember_first[i] = proxy;
            new (&remember_asecond[i]) FlatMatrix<SIMD<double>> (mat);
            return;
          }
      }
    throw Exception ("no space for userdata - memory available");
  }

  
  void AssignMemory (const CoefficientFunction * cf, size_t h, size_t w, LocalHeap & lh)
  {
    for (size_t i = 0; i < remember_cf_first.Size(); i++)
      {
        if (remember_cf_first[i] == nullptr)
          {
            remember_cf_first[i] = cf;
            new (&remember_cf_second[i]) FlatMatrix<double> (h, w, lh);            
            new (&remember_cf_asecond[i]) FlatMatrix<SIMD<double>> (w, (h+SIMD<double>::Size()-1)/SIMD<double>::Size(), lh);
            remember_cf_computed[i] = false;
            return;
          }
      }
    throw Exception ("no space for userdata - memory available");
  }

  void AssignMemory (const CoefficientFunction * cf, FlatMatrix<SIMD<double>> mat)
  {
    for (size_t i = 0; i < remember_cf_first.Size(); i++)
      {
        if (remember_cf_first[i] == nullptr || remember_cf_first[i] == cf)
          {
            remember_cf_first[i] = cf;
            new (&remember_cf_asecond[i]) FlatMatrix<SIMD<double>> (mat);
            remember_cf_computed[i] = true;            
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
    return remember_second[remember_first.PosSure(proxy)];
  }
  FlatMatrix<SIMD<double>> GetAMemory (const ProxyFunction * proxy) const
  {
    return remember_asecond[remember_first.PosSure(proxy)];
  }
  FlatMatrix<double> GetMemory (const CoefficientFunction * cf) const
  {
    return remember_cf_second[remember_cf_first.PosSure(cf)];
  }
  FlatMatrix<SIMD<double>> GetAMemory (const CoefficientFunction * cf) const
  {
    return remember_cf_asecond[remember_cf_first.PosSure(cf)];
  }
  bool Computed (const CoefficientFunction * cf) const
  {
    return remember_cf_computed[remember_cf_first.PosSure(cf)];
  }
  void SetComputed (const CoefficientFunction * cf, bool val = true) const
  {
    remember_cf_computed[remember_cf_first.PosSure(cf)] = val;
  }
};

  





  class SymbolicLinearFormIntegrator : public LinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> proxies;
    Array<CoefficientFunction*> gridfunction_cfs;
    Array<CoefficientFunction*> cache_cfs;
    Array<shared_ptr<CoefficientFunction>> dcf_dtest;

    VorB vb;
    // bool element_boundary;
    VorB element_vb;

  public:
    NGS_DLL_HEADER SymbolicLinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb,
                                                 VorB aelement_vb);

    virtual VorB VB() const override { return vb; }
    virtual string Name () const override { return string ("Symbolic LFI"); }
    virtual int GetDimension() const override { return proxies[0]->Evaluator()->BlockDim(); }

    virtual void 
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatVector<double> elvec,
		       LocalHeap & lh) const override;
      
    virtual void 
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatVector<Complex> elvec,
		       LocalHeap & lh) const override;

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
    Array<CoefficientFunction*> cache_cfs;
    Array<int> trial_cum, test_cum;   // cumulated dimension of proxies
    VorB vb;           // on the boundary of the domain ? 
    // bool element_boundary;
    VorB element_vb;   // on the boundary of the element ? 
    Matrix<bool> nonzeros;    // do components interact ?
    Matrix<bool> nonzeros_deriv;   // do components interact ? 
    Matrix<bool> nonzeros_proxies; // do proxies interact ?
    Matrix<bool> diagonal_proxies; // do proxies interact diagonally ?
    Matrix<bool> same_diffops; // are diffops the same ? 
    bool elementwise_constant;

    int trial_difforder, test_difforder;
    bool is_symmetric;
    bool has_interpolate; // is there an interpolate in the expression tree ? 
    shared_ptr<BilinearFormIntegrator> linearization;
    Array<shared_ptr<CoefficientFunction>> dcf_dtest;  // derivatives by test-functions
    Matrix<shared_ptr<CoefficientFunction>> ddcf_dtest_dtrial;  // derivatives by test- and trial-functions
  public:
    NGS_DLL_HEADER SymbolicBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb,
                                                   VorB aelement_boundary);

    virtual VorB VB() const override { return vb; }
    virtual VorB ElementVB() const { return element_vb; }
    virtual xbool IsSymmetric() const override { return is_symmetric ? xbool(true) : xbool(maybe); } 
    virtual string Name () const override { return string ("Symbolic BFI"); }

    using Integrator::GetIntegrationRule;
    NGS_DLL_HEADER virtual const IntegrationRule& GetIntegrationRule (const FiniteElement & fel, LocalHeap & lh) const;
    NGS_DLL_HEADER virtual const SIMD_IntegrationRule& Get_SIMD_IntegrationRule (const FiniteElement & fel, LocalHeap & lh) const;
    // virtual IntegrationRule GetIntegrationRuleEB (const FiniteElement & fel, int facetnr, LocalHeap & lh) const;
    // virtual SIMD_IntegrationRule Get_SIMD_IntegrationRuleEB (const FiniteElement & fel, int facetnr, LocalHeap & lh) const;
    
    virtual int GetDimension() const override { return trial_proxies[0]->Evaluator()->BlockDim(); }
    void SetLinearization(shared_ptr<BilinearFormIntegrator> _lin)
    { linearization = _lin; }

    NGS_DLL_HEADER virtual void 
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void 
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatMatrix<Complex> elmat,
		       LocalHeap & lh) const override;    

    NGS_DLL_HEADER virtual void 
    CalcElementMatrixAdd (const FiniteElement & fel,
                          const ElementTransformation & trafo, 
                          FlatMatrix<double> elmat,
                          bool & symmetric_so_far,                          
                          LocalHeap & lh) const override;
    
    NGS_DLL_HEADER virtual void 
    CalcElementMatrixAdd (const FiniteElement & fel,
                          const ElementTransformation & trafo, 
                          FlatMatrix<Complex> elmat,
                          bool & symmetric_so_far,                          
                          LocalHeap & lh) const override;    

    
    template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
    void T_CalcElementMatrixAdd (const FiniteElement & fel,
                                 const ElementTransformation & trafo, 
                                 FlatMatrix<SCAL_RES> elmat,
                                 bool & symmetric_so_far, 
                                 LocalHeap & lh) const;

    template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
    void T_CalcElementMatrixEBAdd (const FiniteElement & fel,
                                   const ElementTransformation & trafo, 
                                   FlatMatrix<SCAL_RES> elmat,
                                   LocalHeap & lh) const;

    template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
    void T_CalcElementMatrixAddShapeWise (const FiniteElement & fel,
                                          const ElementTransformation & trafo, 
                                          FlatMatrix<SCAL_RES> elmat,
                                          LocalHeap & lh) const;

    
    NGS_DLL_HEADER virtual void 
    CalcLinearizedElementMatrix (const FiniteElement & fel,
                                 const ElementTransformation & trafo, 
				 FlatVector<double> elveclin,
                                 FlatMatrix<double> elmat,
                                 LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void
    CalcLinearizedElementMatrix (const FiniteElement & fel,
                                 const ElementTransformation & trafo,
				 FlatVector<Complex> elveclin,
                                 FlatMatrix<Complex> elmat,
                                 LocalHeap & lh) const override;

    template<typename SCAL> void
    T_CalcLinearizedElementMatrixFrozen (const FiniteElement & fel,
                                         const ElementTransformation & trafo,
                                         FlatVector<SCAL> elveclin,
                                         FlatMatrix<SCAL> elmat,
                                         LocalHeap & lh) const;

    template <typename SCAL, typename SCAL_SHAPES>
    void T_CalcLinearizedElementMatrixEB (const FiniteElement & fel,
                                          const ElementTransformation & trafo, 
                                          FlatVector<double> elveclin,
                                          FlatMatrix<double> elmat,
                                          LocalHeap & lh) const;
    
    NGS_DLL_HEADER virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & trafo, 
			const FlatVector<double> elx, 
			FlatVector<double> ely,
			void * precomputed,
			LocalHeap & lh) const override;

    template <typename SCAL, typename SCAL_SHAPES>
    void T_ApplyElementMatrixEB (const FiniteElement & fel, 
                                 const ElementTransformation & trafo, 
                                 const FlatVector<double> elx, 
                                 FlatVector<double> ely,
                                 void * precomputed,
                                 LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void 
    ApplyElementMatrixTrans (const FiniteElement & fel, 
                             const ElementTransformation & trafo, 
                             const FlatVector<double> elx, 
                             FlatVector<double> ely,
                             void * precomputed,
                             LocalHeap & lh) const override;
    
    template <typename SCAL, typename SCAL_SHAPES>
    void T_ApplyElementMatrixTransEB (const FiniteElement & fel, 
                                      const ElementTransformation & trafo, 
                                      const FlatVector<double> elx, 
                                      FlatVector<double> ely,
                                      void * precomputed,
                                      LocalHeap & lh) const;

    
    const auto & GetCoefficientFunction() { return cf; }
    const auto & TrialProxies() { return trial_proxies; } 
    const auto & TestProxies() { return test_proxies; } 
    const auto & GridFunctionCoefficients() { return gridfunction_cfs; } 
  };



  class SymbolicFacetLinearFormIntegrator : public FacetLinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> proxies;
    Array<CoefficientFunction*> cache_cfs;
    Array<int> test_cum;    // cumulated dimension of proxies
    VorB vb;                // only BND supported by now
    // bool element_boundary;  /// not needed (by now ???)
    IntegrationRule ir;   // if non-empty use this integration-rule
    SIMD_IntegrationRule simd_ir;   // if non-empty use this integration-rule

  public:
    SymbolicFacetLinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb);

    virtual VorB VB() const override { return vb; }
    virtual bool BoundaryForm() const override { return vb == BND; }

    NGS_DLL_HEADER void
    CalcFacetVector (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
			 FlatVector<double> elvec,
			 LocalHeap & lh) const override;   

    NGS_DLL_HEADER void 
    CalcFacetVector (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,	 
			 FlatVector<Complex> elvec,
			 LocalHeap & lh) const override;

    NGS_DLL_HEADER void
    CalcFacetVector(const FiniteElement & volumefel, int LocalFacetNr,
                    const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                    const ElementTransformation & seltrans,
                    FlatVector<double> elvec,
                    LocalHeap & lh) const override;

    NGS_DLL_HEADER void
    CalcFacetVector(const FiniteElement & volumefel, int LocalFacetNr,
                    const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                    const ElementTransformation & seltrans,
                    FlatVector<Complex> elvec,
                    LocalHeap & lh) const override;

  private:
    template<typename TSCAL>
    void T_CalcFacetVector (const FiniteElement & fel1, int LocalFacetNr1,
                       const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                       const FiniteElement & fel2, int LocalFacetNr2,
                       const ElementTransformation & trafo2, FlatArray<int> & ElVertices2, 
                       FlatVector<TSCAL> elvec,
                       LocalHeap &lh) const;

    template<typename TSCAL>
    void T_CalcFacetVector (const FiniteElement & volumefel, int LocalFacetNr,
                            const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                            const ElementTransformation & seltrans,
                            FlatVector<TSCAL> elvec,
                            LocalHeap & lh) const;
  };


  

  class SymbolicFacetBilinearFormIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> trial_proxies, test_proxies;
    Array<CoefficientFunction*> gridfunction_cfs;
    Array<CoefficientFunction*> cache_cfs;
    Array<int> trial_cum, test_cum;   // cumulated dimension of proxies
    VorB vb;
    bool element_boundary;
    bool neighbor_testfunction;
    Array<shared_ptr<CoefficientFunction>> dcf_dtest;  // derivatives by test-functions    
  public:
    NGS_DLL_HEADER SymbolicFacetBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb, bool aelement_boundary);

    virtual VorB VB() const { return vb; }
    virtual bool BoundaryForm() const { return vb == BND; }
    virtual xbool IsSymmetric() const { return maybe; } 
    
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
    CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<Complex> elmat,
                     LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                     const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                     const ElementTransformation & seltrans, FlatArray<int> & SElVertices,  
                     FlatMatrix<Complex> elmat,
                     LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    CalcLinearizedFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                               const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                               const ElementTransformation & seltrans, FlatArray<int> & SElVertices,  
                               FlatVector<double> vec, FlatMatrix<double> elmat,
                               LocalHeap & lh) const;
    
    NGS_DLL_HEADER virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                      const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                      const FiniteElement & volumefel2, int LocalFacetNr2,
                      const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                      FlatVector<double> elx, FlatVector<double> ely,
                      LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    CalcTraceValues (const FiniteElement & volumefel, int LocalFacetNr,
		     const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
		     FlatVector<double> & trace, FlatVector<double> elx, LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    ApplyFromTraceValues (const FiniteElement & volumefel, int LocalFacetNr,
			  const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			  FlatVector<double> trace,
			  FlatVector<double> elx, FlatVector<double> ely, 
			  LocalHeap & lh) const;
    
    NGS_DLL_HEADER virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                      const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                      const ElementTransformation & seltrans, FlatArray<int> & SElVertices,
                      FlatVector<double> elx, FlatVector<double> ely,
                      LocalHeap & lh) const;

  private:
    template<typename TSCAL, typename SCAL_SHAPES = double>
    void T_CalcFacetMatrix(const FiniteElement & volumefel1, int LocalFacetNr1,
                           const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                           const FiniteElement & volumefel2, int LocalFacetNr2,
                           const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                           FlatMatrix<TSCAL> elmat,
                           LocalHeap & lh) const;

    template<typename TSCAL, typename SCAL_SHAPES = double>
    void T_CalcFacetMatrix(const FiniteElement & volumefel, int LocalFacetNr,
                           const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                           const ElementTransformation & seltrans, FlatArray<int> & SElVertices,  
                           FlatMatrix<TSCAL> elmat,
                           LocalHeap & lh) const;
  };

  class SymbolicEnergy : public BilinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> cf;
    VorB vb;
    Array<ProxyFunction*> trial_proxies;
    VorB element_vb;    

    Timer<TNoTracing> timer{"SymbolicEnergy"};
    Array<int> trial_cum;     // cumulated dimension of proxies
    Matrix<bool> nonzeros;    // do components interact ? 
    Matrix<bool> nonzeros_proxies; // do proxies interact ?
    Array<shared_ptr<CoefficientFunction>> dcf, ddcf;
    
  public:
    SymbolicEnergy (shared_ptr<CoefficientFunction> acf, VorB avb, VorB aelement_vb);

    virtual VorB VB() const { return vb; }
    virtual xbool IsSymmetric() const { return true; } 
    virtual string Name () const { return string ("Symbolic Energy"); }


    using Integrator::GetIntegrationRule;
    NGS_DLL_HEADER virtual const IntegrationRule& GetIntegrationRule (const FiniteElement & fel, LocalHeap & lh) const;
    NGS_DLL_HEADER virtual const SIMD_IntegrationRule& Get_SIMD_IntegrationRule (const FiniteElement & fel, LocalHeap & lh) const;
    
    virtual void 
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const;

    virtual void 
    CalcLinearizedElementMatrix (const FiniteElement & fel,
                                 const ElementTransformation & trafo, 
				 FlatVector<double> elveclin,
                                 FlatMatrix<double> elmat,
                                 LocalHeap & lh) const;

    void AddLinearizedElementMatrix (const FiniteElement & fel,
                                     ProxyUserData & ud, 
                                     const BaseMappedIntegrationRule & mir, 
                                     FlatVector<double> elveclin,
                                     FlatMatrix<double> elmat,
                                     LocalHeap & lh) const;

    void AddLinearizedElementMatrix (const FiniteElement & fel,
                                     const ElementTransformation & trafo, 
                                     const SIMD_BaseMappedIntegrationRule & mir, 
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
