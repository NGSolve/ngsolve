#ifndef NGBEM_hpp
#define NGBEM_hpp

#include "bem_diffops.hpp"
#include "diffopwithfactor.hpp"
#include "kernels.hpp"
#include "../fem/integratorcf.hpp"

namespace ngsbem
{
  using namespace ngcomp;
  class BasePotentialCF;

  class IntOp_Parameters
  {
    bool use_fmm = true;
    int fmm_maxdirect = 100;
    int fmm_minorder = 20;
  public:
    IntOp_Parameters () = default;
    IntOp_Parameters (const Flags & flags);
    
    bool UseFMM() const { return use_fmm; }
    int FMMMaxDirect() const { return fmm_maxdirect; }
    int FMMMinOrder() const { return fmm_minorder; }

    operator FMM_Parameters() const
    {
      FMM_Parameters fmm_params;
      fmm_params.maxdirect = fmm_maxdirect;
      fmm_params.minorder = fmm_minorder;      
      return fmm_params;
    }
      
  };

  inline ostream & operator<< (ostream & ost, const IntOp_Parameters & ioflags)
  {
    ost << "use_fmm = " << ioflags.UseFMM() << endl;
    ost << "fmm_maxdirect = " << ioflags.FMMMaxDirect() << endl;
    ost << "fmm_minorder = " << ioflags.FMMMinOrder() << endl;    
    return ost;
  }

  
  /** The IntegralOperator provides methods for the assembly of and access to a matrix 
      resulting from a variational formulation of a boundary integral equation.*/
  class IntegralOperator 
  {
  protected:
    shared_ptr<FESpace> trial_space; 
    shared_ptr<FESpace> test_space; 

    optional<Region> trial_definedon;
    optional<Region> test_definedon;

    shared_ptr<DifferentialOperator> trial_evaluator;
    // shared_ptr<CoefficientFunction> trial_factor;
    
    shared_ptr<DifferentialOperator> test_evaluator;
    // shared_ptr<CoefficientFunction> test_factor;
    
    
    // integration order
    int intorder;
    IntOp_Parameters io_params;

    // Sauter-Schwab integration rules:
    Array<Vec<2>> identic_panel_x, identic_panel_y;
    Array<double> identic_panel_weight;
    
    Array<Vec<2>> common_vertex_x, common_vertex_y;
    Array<double> common_vertex_weight;

    Array<Vec<2>> common_edge_x, common_edge_y;
    Array<double> common_edge_weight;
    

    shared_ptr<BaseMatrix> matrix;

  public:

    IntegralOperator (shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                      optional<Region> _definedon_trial, optional<Region> _definedon_test,
                      shared_ptr<DifferentialOperator> _trial_evaluator, //  shared_ptr<CoefficientFunction> _trial_factor,
                      shared_ptr<DifferentialOperator> _test_evaluator, // shared_ptr<CoefficientFunction> _test_factor,
                      int _intorder, const IntOp_Parameters & _io_params);
    
    virtual ~IntegralOperator() = default;

    shared_ptr<BaseMatrix> GetMatrix() const { return matrix; }

    virtual shared_ptr<BaseMatrix> CreateMatrixFMM(LocalHeap & lh) const = 0;

    virtual shared_ptr<BasePotentialCF> GetPotential(shared_ptr<GridFunction> gf,
                                                     optional<int> io, bool nearfield_experimental) const = 0;
  };


  
  /** The GenericIntegralOperator is a templated #IntegralOperator, the template type is 
      the kernel the specific potential,i.e. a fundamental solution or its 
      derivative of specific pde. */
  template <typename KERNEL>
  class GenericIntegralOperator : public IntegralOperator // <typename KERNEL::value_type>
  {
    KERNEL kernel;
    typedef typename KERNEL::value_type value_type;
    typedef IntegralOperator BASE;

    
  public:
    /*
    GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                            optional<Region> _definedon_trial, optional<Region> _definedon_test,
                            shared_ptr<DifferentialOperator> _trial_evaluator, shared_ptr<CoefficientFunction> _trial_factor,
                            shared_ptr<DifferentialOperator> _test_evaluator, shared_ptr<CoefficientFunction> _test_factor,
                            KERNEL _kernel,
                            int _intorder);
    */

    GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                            optional<Region> _definedon_trial, optional<Region> _definedon_test,
                            shared_ptr<DifferentialOperator> _trial_evaluator, 
                            shared_ptr<DifferentialOperator> _test_evaluator, 
                            KERNEL _kernel,
                            int _intorder, const IntOp_Parameters & _io_params = IntOp_Parameters());
    /*
      : GenericIntegralOperator (_trial_space, _test_space, _definedon_trial, _definedon_test,
                                 _trial_evaluator, nullptr, _test_evaluator, nullptr,
                                 _kernel, _intorder) { } 
    */


    
    void CalcElementMatrix(FlatMatrix<value_type> matrix,
                           ElementId ei_trial, ElementId ei_test,
                           LocalHeap &lh) const;

    
    shared_ptr<BaseMatrix> CreateMatrixFMM(LocalHeap & lh) const override;
    
    virtual shared_ptr<BasePotentialCF> GetPotential(shared_ptr<GridFunction> gf,
                                                         optional<int> io, bool nearfield_experimental) const override;
  };




  class BasePotentialCF : public CoefficientFunctionNoDerivative
  {
  protected:
    shared_ptr<GridFunction> gf;
    optional<Region> definedon;
    shared_ptr<DifferentialOperator> evaluator;
    IntOp_Parameters io_params;
       
  public:
    BasePotentialCF (shared_ptr<GridFunction> _gf,
                     optional<Region> _definedon,    
                     shared_ptr<DifferentialOperator> _evaluator, bool is_complex)
      : CoefficientFunctionNoDerivative (_evaluator->Dim(), is_complex), 
        gf(_gf), definedon(_definedon), evaluator(_evaluator) { } 

    virtual  ~BasePotentialCF() = default;
    
    virtual void BuildLocalExpansion(const Region & reg) = 0;
  };

  

  template  <typename KERNEL>
  class PotentialCF : public BasePotentialCF
  {
    KERNEL kernel;
    int intorder;
    bool nearfield;

    using LOCAL_EXPANSION = typename std::invoke_result_t<decltype(&KERNEL::CreateLocalExpansion),KERNEL,Vec<3>,double,FMM_Parameters>;

    LOCAL_EXPANSION local_expansion;
    
  public:
    PotentialCF (shared_ptr<GridFunction> _gf,
                 optional<Region> _definedon,    
                 shared_ptr<DifferentialOperator> _evaluator,
                 KERNEL _kernel, int _intorder, bool anearfield);



    void BuildLocalExpansion(const Region & reg) override;
    
    using CoefficientFunctionNoDerivative::Evaluate;
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("eval not implemented"); }

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<> result) const override
    { T_Evaluate(ip, result); }
    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<Complex> result) const override
    { T_Evaluate(ip, result); }

    virtual void Evaluate(const BaseMappedIntegrationRule & ir,
			  BareSliceMatrix<> result) const override
    { T_Evaluate(ir, result); }
    virtual void Evaluate(const BaseMappedIntegrationRule & ir,
			  BareSliceMatrix<Complex> result) const override
    { T_Evaluate(ir, result); }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           BareSliceMatrix<SIMD<double>> result) const override
    { T_Evaluate(ir, result); }
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           BareSliceMatrix<SIMD<Complex>> result) const override
    { T_Evaluate(ir, result); }
    
  private:
    template <typename T>
    void T_Evaluate(const BaseMappedIntegrationPoint & ip,
                    FlatVector<T> result) const;
    template <typename T>
    void T_Evaluate(const BaseMappedIntegrationRule & ir,
                    BareSliceMatrix<T> result) const;
    template <typename T>
    void T_Evaluate(const SIMD_BaseMappedIntegrationRule & ir,
                    BareSliceMatrix<SIMD<T>> result) const;
  };





  
  
  


  

  class BasePotentialOperator
  {
  public:
    shared_ptr<ProxyFunction> proxy;
    // shared_ptr<CoefficientFunction> factor;
    optional<Region> definedon;
    shared_ptr<DifferentialOperator> evaluator;
    IntOp_Parameters io_params;
    int intorder;
    
  public:
    BasePotentialOperator (shared_ptr<ProxyFunction> _proxy, // shared_ptr<CoefficientFunction> _factor,
                           optional<Region> _definedon,    
                           shared_ptr<DifferentialOperator> _evaluator,
                           int _intorder)
      : proxy(_proxy), /* factor(_factor), */ definedon(_definedon), evaluator(_evaluator), intorder(_intorder) { ; }

    BasePotentialOperator (shared_ptr<ProxyFunction> _proxy, // shared_ptr<CoefficientFunction> _factor,
                           optional<Region> _definedon,    
                           shared_ptr<DifferentialOperator> _evaluator,
                           IntOp_Parameters _io_params, 
                           int _intorder)
      : proxy(_proxy), definedon(_definedon), evaluator(_evaluator), io_params(_io_params), intorder(_intorder) { ; } 



    
    virtual ~BasePotentialOperator() { } 
    virtual shared_ptr<IntegralOperator> MakeIntegralOperator(shared_ptr<ProxyFunction> test_proxy,
                                                              // shared_ptr<CoefficientFunction> test_factor,
                                                              DifferentialSymbol dx) = 0;
    virtual shared_ptr<BasePotentialCF> MakePotentialCF(shared_ptr<GridFunction> gf) = 0;

    virtual void Print (ostream & ost) const = 0;
  };

  inline ostream & operator<< (ostream & ost, const BasePotentialOperator & intop)
  {
    intop.Print(ost);
    return ost;
  }
  
  template  <typename KERNEL>
  class PotentialOperator : public BasePotentialOperator
  {
  public:
    KERNEL kernel;
  public:
    PotentialOperator (shared_ptr<ProxyFunction> _proxy,
                       // shared_ptr<CoefficientFunction> _factor,                       
                       optional<Region> _definedon,    
                       shared_ptr<DifferentialOperator> _evaluator,
                       KERNEL _kernel, int _intorder)
      : BasePotentialOperator (_proxy, /* _factor, */ _definedon, _evaluator, _intorder), kernel(_kernel) { ; }

    PotentialOperator (shared_ptr<ProxyFunction> _proxy,
                       // shared_ptr<CoefficientFunction> _factor,                       
                       optional<Region> _definedon,    
                       shared_ptr<DifferentialOperator> _evaluator,
                       KERNEL _kernel,
                       IntOp_Parameters _io_params, 
                       int _intorder)
      : BasePotentialOperator (_proxy, _definedon, _evaluator, _io_params, _intorder), kernel(_kernel) { ; }

    shared_ptr<IntegralOperator> MakeIntegralOperator(shared_ptr<ProxyFunction> test_proxy,
                                                      // shared_ptr<CoefficientFunction> test_factor,
                                                      DifferentialSymbol dx) override
    {
      auto festest = test_proxy->GetFESpace();
      
      auto tmpfes = festest;
      auto tmpeval = proxy->Evaluator();
      while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
        {
          int component = compeval->Component();
          tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[component];
          tmpeval = compeval->BaseDiffOp();
        }
      
      optional<Region> definedon_test;
      if (dx.definedon)
        definedon_test = Region(festest->GetMeshAccess(), dx.vb, get<1> (*(dx.definedon)));
      
      return make_shared<GenericIntegralOperator<KERNEL>> (proxy->GetFESpace(),
                                                           festest, 
                                                           definedon,
                                                           definedon_test,
                                                           proxy->Evaluator(), // nullptr, // factor,
                                                           test_proxy->Evaluator(), // nullptr, // test_factor, 
                                                           kernel,
                                                           2 + intorder + tmpfes->GetOrder()+dx.bonus_intorder,
                                                           io_params);
    }
    
    shared_ptr<BasePotentialCF> MakePotentialCF(shared_ptr<GridFunction> gf) override
    {
      return make_shared<PotentialCF<KERNEL>>(gf, definedon, evaluator, kernel, 2+intorder, true);
    }

    virtual void Print (ostream & ost) const override
    {
      ost << "Potential operator:" << endl;
      ost << "Kernel: " << kernel.Name() << endl;
      ost << io_params << endl;
    }
  };  



  inline Array < tuple <shared_ptr<ProxyFunction>, shared_ptr<CoefficientFunction>  >>
  CreateProxyLinearization (shared_ptr<CoefficientFunction> cf, bool trial)
  {
    Array<tuple<shared_ptr<ProxyFunction>, shared_ptr<CoefficientFunction>>>  proxylin;
    Array<ProxyFunction*> proxies;
    
    cf->TraverseTree
      ([&] (CoefficientFunction & nodecf)
      {
        if (auto proxy = dynamic_cast<ProxyFunction*> (&nodecf))
          if ((proxy->IsTrialFunction() == trial) &&  (!proxies.Contains(proxy)))
            proxies.Append (proxy);
      });
    
    for (auto proxy : proxies)
      if (!proxy->Evaluator()->SupportsVB(BND))
        throw Exception ("Testfunction does not support BND-forms, maybe a Trace() operator is missing");

    for (auto proxy : proxies)
      {
        CoefficientFunction::T_DJC cache;
        shared_ptr<ProxyFunction> spproxy = dynamic_pointer_cast<ProxyFunction> (proxy->shared_from_this());
        proxylin += tuple { spproxy, cf->DiffJacobi(proxy, cache) };
      }

    return proxylin;
  }

  
  inline int GetFESOrder (shared_ptr<ProxyFunction> proxy)
  {
    auto fes = proxy->GetFESpace();
    
    auto tmpfes = fes;
    auto tmpeval = proxy->Evaluator();

    while (true)
      {
        if (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))\
          {
            tmpfes = (*dynamic_pointer_cast<CompoundFESpace>(tmpfes))[compeval->Component()];
            tmpeval = compeval->BaseDiffOp();
          }
        else if (auto diffopfac = dynamic_pointer_cast<DifferentialOperatorWithFactor> (tmpeval))
          {
            tmpeval = diffopfac->BaseDiffOp();
          }
        else
          break;
      }
    
    return tmpfes->GetOrder();
  }
  
  inline tuple < shared_ptr<ProxyFunction>, shared_ptr<CoefficientFunction> >
  GetProxyAndFactor (shared_ptr<CoefficientFunction> cf, bool trial)
  {
    if (auto proxy = dynamic_pointer_cast<ProxyFunction>(cf))
      return { proxy, nullptr };
    
    auto proxylin = CreateProxyLinearization (cf, trial);
    if (proxylin.Size() != 1)
      throw Exception(string("need exactly one")+(trial?"trial":"test") + "-proxy");

    // return proxylin[0];
    auto [proxy,factor] = proxylin[0];

    auto diffopwith = make_shared<DifferentialOperatorWithFactor> (proxy->Evaluator(),
                                                                   factor->Reshape(cf->Dimension(), proxy->Dimension()));
    
    return { make_shared<ProxyFunction> (proxy->GetFESpace(), proxy->IsTestFunction(), proxy->IsComplex(),
                                         diffopwith, nullptr, nullptr, nullptr, nullptr, nullptr),
             nullptr };

    // return { proxy, factor->Reshape(cf->Dimension(), proxy->Dimension()) };
  }

  inline shared_ptr<ProxyFunction>
  GetProxyWithFactor (shared_ptr<CoefficientFunction> cf, bool trial)
  {
    return std::get<0> (GetProxyAndFactor(cf, trial));
  }

  class BasePotentialOperatorAndTest
  {
    shared_ptr<BasePotentialOperator> pot;
    shared_ptr<ProxyFunction> test_proxy_with_factor;    
    /*
    shared_ptr<ProxyFunction> test_proxy;
    shared_ptr<CoefficientFunction> test_factor;
    */
  public:

    BasePotentialOperatorAndTest (shared_ptr<BasePotentialOperator> _pot,
                                  shared_ptr<CoefficientFunction> _test_proxy)
      : pot(_pot)
    {
      // tie(test_proxy,test_factor) = GetProxyAndFactor(_test_proxy, false);
      test_proxy_with_factor = GetProxyWithFactor(_test_proxy, false);      
    }

    
    shared_ptr<IntegralOperator> MakeIntegralOperator (DifferentialSymbol dx)
    {
      // return pot->MakeIntegralOperator(test_proxy, test_factor, dx);
      return pot->MakeIntegralOperator(test_proxy_with_factor, dx);
    }
  };



}



#endif

