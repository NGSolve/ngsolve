#ifndef FILE_SYMBOLICINTEGRATOR
#define FILE_SYMBOLICINTEGRATOR

/*********************************************************************/
/* File:   symbolicintegrator.hpp                                    */
/* Author: Joachim Schoeberl                                         */
/* Date:   August 2015                                               */
/*********************************************************************/


namespace ngfem
{

  
class ProxyUserData
{
public:
  class ProxyFunction * testfunction = nullptr;
  int test_comp;
  class ProxyFunction * trialfunction = nullptr;
  int trial_comp;

  const FiniteElement * fel = nullptr;
  const FlatVector<double> * elx;
  LocalHeap * lh;
};

class ProxyFunction : public CoefficientFunction
{
  // const FESpace * space;
  bool testfunction; // true .. test, false .. trial
  bool is_complex;

  shared_ptr<DifferentialOperator> evaluator;
  shared_ptr<DifferentialOperator> deriv_evaluator;
  shared_ptr<DifferentialOperator> trace_evaluator;
public:
  ProxyFunction (// const FESpace * aspace, 
                 bool atestfunction, bool ais_complex,
                 shared_ptr<DifferentialOperator> aevaluator, 
                 shared_ptr<DifferentialOperator> aderiv_evaluator,
                 shared_ptr<DifferentialOperator> atrace_evaluator)
                 
    : // space(aspace), 
    testfunction(atestfunction), is_complex(ais_complex),
    evaluator(aevaluator), 
    deriv_evaluator(aderiv_evaluator),
    trace_evaluator(atrace_evaluator)
  { ; }

  bool IsTestFunction () const { return testfunction; }
  virtual int Dimension () const { return evaluator->Dim(); }
  virtual bool IsComplex () const { return is_complex; } // space->IsComplex(); }

  const shared_ptr<DifferentialOperator> & Evaluator() const { return evaluator; }
  const shared_ptr<DifferentialOperator> & DerivEvaluator() const { return deriv_evaluator; }
  const shared_ptr<DifferentialOperator> & TraceEvaluator() const { return trace_evaluator; }
  shared_ptr<ProxyFunction> Deriv() const
  {
    return make_shared<ProxyFunction> (testfunction, is_complex, deriv_evaluator, nullptr, nullptr);
  }
  shared_ptr<ProxyFunction> Trace() const
  {
    return make_shared<ProxyFunction> (testfunction, is_complex, trace_evaluator, nullptr, nullptr);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    Vector<> tmp(Dimension());
    Evaluate (ip, tmp);
    return tmp(0);
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const
  {
    ProxyUserData * ud = (ProxyUserData*)ip.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    if (!testfunction && ud->fel)
      {
        evaluator->Apply (*ud->fel, ip, *ud->elx, result, *ud->lh);
        return;
      }

    result = 0;
    if (ud->testfunction == this)
      result (ud->test_comp) = 1;
    if (ud->trialfunction == this)
      result (ud->trial_comp) = 1;
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<Complex> result) const
  {
    Vector<> result_double(result.Size());
    Evaluate (ip, result_double);
    result = result_double;
  }

  virtual void EvaluateDeriv (const BaseMappedIntegrationPoint & ip,
                              FlatVector<> result,
                              FlatVector<> deriv) const
  {
    ProxyUserData * ud = (ProxyUserData*)ip.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    deriv = 0;
    result = 0;

    if (!testfunction && ud->fel)
      evaluator->Apply (*ud->fel, ip, *ud->elx, result, *ud->lh);

    if (ud->testfunction == this)
      result(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv(ud->trial_comp) = 1;
  }


  virtual void EvaluateDDeriv (const BaseMappedIntegrationPoint & ip,
                               FlatVector<> result,
                               FlatVector<> deriv,
                               FlatVector<> dderiv) const
  {
    ProxyUserData * ud = (ProxyUserData*)ip.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    result = 0;
    deriv = 0;
    dderiv = 0;

    if (!testfunction && ud->fel)
      evaluator->Apply (*ud->fel, ip, *ud->elx, result, *ud->lh);

    if (ud->testfunction == this)
      deriv(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv(ud->trial_comp) = 1;
  }

};


class CompoundDifferentialOperator : public DifferentialOperator
{
  shared_ptr<DifferentialOperator> diffop;
  int comp;
public:
  CompoundDifferentialOperator (shared_ptr<DifferentialOperator> adiffop, 
                                int acomp)
    : diffop(adiffop), comp(acomp) { ; }
  
  virtual ~CompoundDifferentialOperator () = default;
  
  /// dimension of range
  virtual int Dim() const { return diffop->Dim(); }
  virtual bool Boundary() const { return diffop->Boundary(); }
  virtual int DiffOrder() const { return diffop->DiffOrder(); }
  
  
  NGS_DLL_HEADER virtual void
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & mip,
              FlatMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const
  {
    mat = 0;
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = fel.GetRange(comp);
    diffop->CalcMatrix (fel[comp], mip, mat.Cols(r), lh);
  }
  
  
  NGS_DLL_HEADER virtual void
  Apply (const FiniteElement & bfel,
         const BaseMappedIntegrationPoint & mip,
         FlatVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = fel.GetRange(comp);
    diffop->Apply (fel[comp], mip, x.Range(r), flux, lh);
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
    IntRange r = fel.GetRange(comp);
    diffop->ApplyTrans (fel[comp], mip, flux, x.Range(r), lh);
  }

};





  class SymbolicLinearFormIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> proxies;
    VorB vb;
    
  public:
    SymbolicLinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb)
      : cf(acf), vb(avb)
    {
      if (cf->Dimension() != 1)
        throw Exception ("SymbolicLFI needs scalar-valued CoefficientFunction");
      cf->TraverseTree
        ([&] (CoefficientFunction & nodecf)
         {
           auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
           if (proxy && !proxies.Contains(proxy))
             proxies.Append (proxy);
         });
    }

    virtual bool BoundaryForm() const { return vb==BND; }

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
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> trial_proxies, test_proxies;
    VorB vb;
    bool element_boundary;

  public:
    SymbolicBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb,
                                    bool aelement_boundary)
      : cf(acf), vb(avb), element_boundary(aelement_boundary)
    {
      if (cf->Dimension() != 1)
        throw Exception ("SymblicBFI needs scalar-valued CoefficientFunction");

      cf->TraverseTree
        ( [&] (CoefficientFunction & nodecf)
          {
            auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
            if (proxy) 
              {
                if (proxy->IsTestFunction())
                  {
                    if (!test_proxies.Contains(proxy))
                      test_proxies.Append (proxy);
                  }
                else
                  {                                         
                    if (!trial_proxies.Contains(proxy))
                      trial_proxies.Append (proxy);
                  }
              }
          });
    }

    virtual bool BoundaryForm() const { return vb == BND; }
    virtual bool IsSymmetric() const { return true; } 

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
      
    template <typename SCAL, typename SCAL_SHAPES = double>
    void T_CalcElementMatrix (const FiniteElement & fel,
                              const ElementTransformation & trafo, 
                              FlatMatrix<SCAL> elmat,
                              LocalHeap & lh) const;




    template <int D, typename SCAL, typename SCAL_SHAPES>
    void T_CalcElementMatrixEB (const FiniteElement & fel,
                                const ElementTransformation & trafo, 
                                FlatMatrix<SCAL> elmat,
                                LocalHeap & lh) const;

  };



  class SymbolicEnergy : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> trial_proxies;

  public:
    SymbolicEnergy (shared_ptr<CoefficientFunction> acf);

    virtual bool BoundaryForm() const { return false; }
    virtual bool IsSymmetric() const { return true; } 

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
