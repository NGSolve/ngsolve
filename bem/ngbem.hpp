#ifndef NGBEM_hpp
#define NGBEM_hpp

#include "bem_diffops.hpp"
#include "mptools.hpp"
#include "../fem/integratorcf.hpp"

namespace ngsbem
{
  using namespace ngcomp;
  class BasePotentialCF;

  
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
    shared_ptr<CoefficientFunction> trial_factor;
    
    shared_ptr<DifferentialOperator> test_evaluator;
    shared_ptr<CoefficientFunction> test_factor;
    
    
    // integration order
    int intorder;


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
                      shared_ptr<DifferentialOperator> _trial_evaluator, shared_ptr<CoefficientFunction> _trial_factor,
                      shared_ptr<DifferentialOperator> _test_evaluator, shared_ptr<CoefficientFunction> _test_factor,
                      int _intorder);
    
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
    GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                            optional<Region> _definedon_trial, optional<Region> _definedon_test,
                            shared_ptr<DifferentialOperator> _trial_evaluator, shared_ptr<CoefficientFunction> _trial_factor,
                            shared_ptr<DifferentialOperator> _test_evaluator, shared_ptr<CoefficientFunction> _test_factor,
                            KERNEL _kernel,
                            int _intorder);


    GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                            optional<Region> _definedon_trial, optional<Region> _definedon_test,
                            shared_ptr<DifferentialOperator> _trial_evaluator, 
                            shared_ptr<DifferentialOperator> _test_evaluator, 
                            KERNEL _kernel,
                            int _intorder)
      : GenericIntegralOperator (_trial_space, _test_space, _definedon_trial, _definedon_test,
                                 _trial_evaluator, nullptr, _test_evaluator, nullptr,
                                 _kernel, _intorder) { } 



    
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

    // typedef decltype (defvar<KERNEL>().CreateLocalExpansion(Vec<3>(), double())) LOCAL_EXPANSION;
    using LOCAL_EXPANSION = typename std::invoke_result_t<decltype(&KERNEL::CreateLocalExpansion),KERNEL,Vec<3>,double>;

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
    shared_ptr<CoefficientFunction> factor;
    optional<Region> definedon;
    shared_ptr<DifferentialOperator> evaluator;
    int intorder;
  public:
    BasePotentialOperator (shared_ptr<ProxyFunction> _proxy, shared_ptr<CoefficientFunction> _factor,
                           optional<Region> _definedon,    
                           shared_ptr<DifferentialOperator> _evaluator,
                           int _intorder)
      : proxy(_proxy), factor(_factor), definedon(_definedon), evaluator(_evaluator), intorder(_intorder) { ; } 
    virtual ~BasePotentialOperator() { } 
    virtual shared_ptr<IntegralOperator> MakeIntegralOperator(shared_ptr<ProxyFunction> test_proxy,
                                                              shared_ptr<CoefficientFunction> test_factor,
                                                              DifferentialSymbol dx) = 0;
    virtual shared_ptr<BasePotentialCF> MakePotentialCF(shared_ptr<GridFunction> gf) = 0;
  };
  
  
  template  <typename KERNEL>
  class PotentialOperator : public BasePotentialOperator
  {
  public:
    KERNEL kernel;
  public:
    PotentialOperator (shared_ptr<ProxyFunction> _proxy,
                       shared_ptr<CoefficientFunction> _factor,                       
                       optional<Region> _definedon,    
                       shared_ptr<DifferentialOperator> _evaluator,
                       KERNEL _kernel, int _intorder)
      : BasePotentialOperator (_proxy, _factor, _definedon, _evaluator, _intorder), kernel(_kernel) { ; }

    shared_ptr<IntegralOperator> MakeIntegralOperator(shared_ptr<ProxyFunction> test_proxy,
                                                      shared_ptr<CoefficientFunction> test_factor,
                                                      DifferentialSymbol dx) override
    {
      auto festest = test_proxy->GetFESpace();
      
      auto tmpfes = festest;
      auto tmpeval = proxy->Evaluator();
      while (auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval))
        {
          // auto compeval = dynamic_pointer_cast<CompoundDifferentialOperator>(tmpeval);
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
                                                           proxy->Evaluator(), factor,
                                                           test_proxy->Evaluator(), test_factor, 
                                                           kernel,
                                                           2 + intorder + tmpfes->GetOrder()+dx.bonus_intorder);
    }
    
    shared_ptr<BasePotentialCF> MakePotentialCF(shared_ptr<GridFunction> gf) override
    {
      return make_shared<PotentialCF<KERNEL>>(gf, definedon, evaluator, kernel, 2+intorder, true);
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
    return { proxy, factor->Reshape(cf->Dimension(), proxy->Dimension()) };
  }        

  class BasePotentialOperatorAndTest
  {
    shared_ptr<BasePotentialOperator> pot;
    shared_ptr<ProxyFunction> test_proxy;
    shared_ptr<CoefficientFunction> test_factor;
  public:

    BasePotentialOperatorAndTest (shared_ptr<BasePotentialOperator> _pot,
                                  shared_ptr<CoefficientFunction> _test_proxy)
      : pot(_pot)
    {
      tie(test_proxy,test_factor) = GetProxyAndFactor(_test_proxy, false);
      /*
      test_proxy = dynamic_pointer_cast<ProxyFunction>(_test_proxy);

      if (!test_proxy)
        {
          auto proxylin = CreateProxyLinearization (_test_proxy, false);

          cout << "got an expression for potential-test, linearization is: " << endl;
          for (auto [factor, proxy] : proxylin)
            {
              cout << "proxy = " << *proxy << endl;
              cout << "factor = " << *factor << endl;
              test_proxy = proxy;
              test_factor = factor -> Reshape (_test_proxy->Dimension(), proxy->Dimension() );
            }
        }
      */
    }

    
    shared_ptr<IntegralOperator> MakeIntegralOperator (DifferentialSymbol dx)
    {
      return pot->MakeIntegralOperator(test_proxy, test_factor, dx);
    }
  };



  //  ****************************  The kernels **********************************


  

  struct KernelTerm
  {
    double fac;
    size_t kernel_comp;
    size_t trial_comp;
    size_t test_comp;
  };
  

  class BaseKernel
  {
  public:
    shared_ptr<SingularMLMultiPole<Complex>> CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      throw Exception("Create Multipole Expansion not implemented");
    }

    shared_ptr<RegularMLMultiPole<Complex>> CreateLocalExpansion (Vec<3> c, double r) const
    {
      throw Exception("Create Local Expansion not implemented");      
    }

    template <typename TV>
    void AddSource (SingularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, const TV & val) const
    {
      throw Exception("Addsource not implemented");            
    }

    template <typename entry, typename TV>
    void AddSourceTrans(SingularMLMultiPole<entry> & mp, Vec<3> pnt, Vec<3> nv, const TV & val) const
    {
      throw Exception("AddSourceTrans not implemented");
    }

    template <typename TV>
    void EvaluateMP (RegularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, const TV & val) const
    {
      throw Exception("Evaluate not implemented");            
    }

    template <typename entry, typename TV>
    void EvaluateMPTrans(RegularMLMultiPole<entry> & mp, Vec<3> pnt, Vec<3> nv, const TV & val) const
    {
      throw Exception("EvaluateMPTrans not implemented");
    }
  };



  
  
  /** LaplaceSLkernel is the kernel for the single layer potential of 
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM> class LaplaceSLKernel;

  /** LaplaceSLkernel in 3D reads 
      $$ G(x-y) = \frac{1}{4\,\pi \, | x-y| }, \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class LaplaceSLKernel<3> : public BaseKernel
  {
  public:
    typedef double value_type;
    static string Name() { return "LaplaceSL"; }
    static auto Shape() { return IVec<2>(1,1); }
    
    template <typename T>        
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      // return 1.0 / (4 * M_PI * norm);   
      return Vec<1,T> (1.0 / (4 * M_PI * norm));
    }
    
    Array<KernelTerm> terms = { KernelTerm{1.0, 0, 0, 0}, };
    
    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLMultiPole<Complex>> (c, r, 1e-16);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLMultiPole<Complex>> (c, r, 1e-16);
    }

    void AddSource (SingularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      mp.AddCharge (pnt, val(0));
    }

    void EvaluateMP (RegularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      val(0) = Real(mp.Evaluate (pnt));
    }
  };


  /** LaplaceDLkernel is the kernel for the double layer potential of 
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM> class LaplaceDLKernel;

  /** LaplaceDLkernel in 3D reads 
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, 
          \frac{ \langle n(y), x-y\rangle }{ | x-y|^3 }, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class LaplaceDLKernel<3> : public BaseKernel
  {
  public:
    typedef double value_type;
    static string Name() { return "LaplaceDL"; }
    static auto Shape() { return IVec<2>(1,1); }
    
    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      // return nxy / (4 * M_PI * norm*norm*norm);
      return Vec<1,T> (nxy / (4 * M_PI * norm*norm*norm));
    }

    Array<KernelTerm> terms = { KernelTerm{1.0, 0, 0, 0}, };    

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLMultiPole<Complex>> (c, r, 1e-16);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLMultiPole<Complex>> (c, r, 1e-16);
    }

    void AddSource (SingularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      mp.AddDipole(pnt, -nv, val(0));
    }

    void AddSourceTrans(SingularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      mp.AddCharge(pnt, val(0));
    }

    void EvaluateMP (RegularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      val(0) = Real(mp.Evaluate (pnt));
    }

    void EvaluateMPTrans(RegularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      val(0) = Real(mp.EvaluateDirectionalDerivative(pnt, nv));
    }
  };


  template <int DIM> class LaplaceHSKernel;
  
  template<>
  class LaplaceHSKernel<3> : public BaseKernel 
  {
  public:
    typedef double value_type;
    static string Name() { return "LaplaceHL"; }
    static auto Shape() { return IVec<2>(3,3); }
    
    template <typename T>        
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      // return 1.0 / (4 * M_PI * norm);   
      return Vec<1,T> (1.0 / (4 * M_PI * norm));
    }
    
    Array<KernelTerm> terms =
      {
        KernelTerm{1.0, 0, 0, 0},
        KernelTerm{1.0, 0, 1, 1},
        KernelTerm{1.0, 0, 2, 2},
      };

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLMultiPole<Vec<3,Complex>>> (c, r, 1e-16);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLMultiPole<Vec<3,Complex>>> (c, r, 1e-16);
    }

    void AddSource (SingularMLMultiPole<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      mp.AddCharge(pnt, val);
    }

    void EvaluateMP (RegularMLMultiPole<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<double> val) const
    {
      val = Real(mp.Evaluate (pnt));
    }
  };


  /** HelmholtzSLkernel is the kernel for the double layer potential of the 
      Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,. $ */
  template <int DIM> class HelmholtzSLKernel;

  /** HelmholtzSLkernel in 3D reads 
      $$ G(x-y) = \frac{1 }{4\,\pi} \,\frac{e^{i\,\kappa \, |x-y| }{|x-y|} \, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class HelmholtzSLKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzSL"; }
    static auto Shape() { return IVec<2>(1,1); }
    
    /** Construction of the kernel specifies the wavenumber $\kappa$. */
    HelmholtzSLKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
      // return kern;
      return Vec<1,decltype(kern)> (kern);
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms = { KernelTerm{1.0, 0, 0, 0}, };    

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLMultiPole<Complex>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLMultiPole<Complex>> (c, r, kappa);
    }

    void AddSource (SingularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      mp.AddCharge(pnt, val(0));
    }

    void EvaluateMP (RegularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val(0) = mp.Evaluate (pnt);
    }
  };


  /** HelmholtzDLkernel is the kernel for the double layer potential of 
      the Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,.$ */
  template <int DIM> class HelmholtzDLKernel;

  /** HelmholtzDLkernel in 3D reads
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \langle n(y), x-y\rangle \cdot \left( 1 - i\,\kappa\, | x-y| \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class HelmholtzDLKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzDL"; }
    static auto Shape() { return IVec<2>(1,1); }
    
    HelmholtzDLKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>    
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
        * nxy * (Complex(1,0)*T(1.) - Complex(0,kappa)*norm);
      // return kern;
      return Vec<1,decltype(kern)> (kern);
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms = { KernelTerm{1.0, 0, 0, 0}, };    

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLMultiPole<Complex>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLMultiPole<Complex>> (c, r, kappa);
    }

    void AddSource (SingularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      mp.AddDipole(pnt, -nv, val(0));
    }

    void AddSourceTrans(SingularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      mp.AddCharge(pnt, val(0));
    }

    void EvaluateMP (RegularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val(0) = mp.Evaluate (pnt);
    }

    void EvaluateMPTrans(RegularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val(0) = mp.EvaluateDirectionalDerivative(pnt, nv);
    }
  };


  template <int DIM> class HelmholtzSLVecKernel;

  template<>
  class HelmholtzSLVecKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzSLVec"; }
    static auto Shape() { return IVec<2>(3,3); }

    HelmholtzSLVecKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
      return Vec<1,decltype(kern)> (kern);
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms =
      {
        KernelTerm{1.0, 0, 0, 0},
        KernelTerm{1.0, 0, 1, 1},
        KernelTerm{1.0, 0, 2, 2},
      };

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLMultiPole<Vec<3,Complex>>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLMultiPole<Vec<3,Complex>>> (c, r, kappa);
    }

    void AddSource (SingularMLMultiPole<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      mp.AddCharge(pnt, val);
    }

    void EvaluateMP (RegularMLMultiPole<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val = mp.Evaluate (pnt);
    }
  };


  template <int DIM> class HelmholtzHSKernel;
  
  template<>
  class HelmholtzHSKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzHS"; }
    static auto Shape() { return IVec<2>(4,4); }
    
    HelmholtzHSKernel (double _kappa) : kappa(_kappa) { }
    template <typename T>        
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxny = InnerProduct(nx, ny);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
      auto kernnxny = -kappa * kappa * kern * nxny;
      // return kern;
      return Vec<2,decltype(kern)> ({kern, kernnxny});
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms =
      {
        KernelTerm{1.0, 0, 0, 0},
        KernelTerm{1.0, 0, 1, 1},
        KernelTerm{1.0, 0, 2, 2},
	KernelTerm{1.0, 1, 3, 3},
      };

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLMultiPole<Vec<6,Complex>>> (c, r, kappa);
    }
 
    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLMultiPole<Vec<6,Complex>>> (c, r, kappa);
    }

    void AddSource (SingularMLMultiPole<Vec<6,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      Vec<6,Complex> charge;
      charge.Range(0,3) = val.Range(0,3);
      charge.Range(3,6) = -kappa * kappa * val(3) * nv;
      mp.AddCharge(pnt, charge);
    }

    void EvaluateMP (RegularMLMultiPole<Vec<6,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      Vec<6,Complex> eval = mp.Evaluate (pnt);
      val.Range(0,3) = eval.Range(0,3);
      val(3) = InnerProduct(eval.Range(3,6), nv);
    }
  };




  

  /** CombinedFieldKernel is a kernel for the combined field integral equation 
      is considered for the Helmholtz equation. */
  template <int DIM> class CombinedFieldKernel;

  /** CombinedFieldKernel in 3D reads
      $$ G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \left( \langle n_y, x-y\rangle (1- i\,\kappa\, | x-y|) - i\,\kappa\,|x-y|^2 \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class CombinedFieldKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "Helmholtz Combined Field"; }
    static auto Shape() { return IVec<2>(1,1); }
    
    CombinedFieldKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>    
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
        * ( nxy * (Complex(1,0)*T(1.) - Complex(0,kappa)*norm)  - Complex(0,kappa)*norm*norm);
      // return kern;
      return Vec<1,decltype(kern)> (kern);
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms = { KernelTerm{1.0, 0, 0, 0}, };    

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLMultiPole<Complex>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLMultiPole<Complex>> (c, r, kappa);
    }

    void AddSource (SingularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      mp.AddCharge(pnt, Complex(0, -kappa)*val(0));
      mp.AddDipole(pnt, -nv, val(0));
    }

    void EvaluateMP (RegularMLMultiPole<Complex> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val(0) = mp.Evaluate (pnt);
    }
  };

  

  

  template <int D> class MaxwellSLKernel;

  template<>
  class MaxwellSLKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "MaxwellSL"; }
    static auto Shape() { return IVec<2>(4,4); }
    
    MaxwellSLKernel (const MaxwellSLKernel&) = default;
    MaxwellSLKernel (MaxwellSLKernel&&) = default;
    MaxwellSLKernel (double _kappa) : kappa(_kappa)
    {
      for (size_t i = 0; i < 3; i++)
        terms += KernelTerm { kappa, 0, i, i };
      terms += KernelTerm { -1.0/kappa, 0, 3, 3 };      
    }
    
    template <typename T>    
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
      return Vec<1,decltype(kern)> (kern);
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms;

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLMultiPole<Vec<4,Complex>>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLMultiPole<Vec<4,Complex>>> (c, r, kappa);
    }

    void AddSource (SingularMLMultiPole<Vec<4,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      Vec<4,Complex> charge;
      charge.Range(0,3) = kappa * val.Range(0, 3);
      charge(3) = -1.0/kappa * val(3);
      mp.AddCharge(pnt, charge);
    }

    void EvaluateMP (RegularMLMultiPole<Vec<4,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val = mp.Evaluate (pnt);
    }
  };




  

  template <int D> class MaxwellDLKernel;

  // https://weggler.github.io/ngbem/short_and_sweet/Maxwell_Formulations.html
  /** MaxwellDLkernel for 3D in matrix representation reads
      $$  \left( \begin{array}{ccc} 0 & -\frac{\partial G_\kappa(x-y)}{\partial x_3} & \frac{\partial G_\kappa(x-y)}{\partial x_2} \\ \frac{\partial G_\kappa(x-y)}{\partial x_3} & 0 & -\frac{\partial G_\kappa(x-y)}{\partial x_1} \\ -\frac{\partial G_\kappa(x-y)}{\partial x_2} & \frac{\partial G_\kappa(x-y)}{\partial x_1} & 0 \end{array}\right)\,,$$ with 
   $$ G_\kappa(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \langle n(y), x-y\rangle \cdot \left( 1 - i\,\kappa\, | x-y| \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class MaxwellDLKernel<3> : public BaseKernel
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "MaxwellDL"; }
    static auto Shape() { return IVec<2>(3,3); }
    
    MaxwellDLKernel (double _kappa) : kappa(_kappa) { }
    
    template <typename T>    
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
        * (Complex(0,kappa)*norm - Complex(1,0)*T(1.)) * (x-y);
      return kern;
    }
    double GetKappa() const { return kappa; }
    Array<KernelTerm> terms =
      {
        KernelTerm{ 1.0, 0, 1, 2},  // factor, comp, trial, test
        KernelTerm{-1.0, 0, 2, 1},
        KernelTerm{ 1.0, 1, 2, 0},
        KernelTerm{-1.0, 1, 0, 2},
        KernelTerm{ 1.0, 2, 0, 1},
        KernelTerm{-1.0, 2, 1, 0},
      };    

    auto CreateMultipoleExpansion (Vec<3> c, double r) const
    {
      return make_shared<SingularMLMultiPole<Vec<3,Complex>>> (c, r, kappa);
    }

    auto CreateLocalExpansion (Vec<3> c, double r) const
    {
      return make_shared<RegularMLMultiPole<Vec<3,Complex>>> (c, r, kappa);
    }

    void AddSource (SingularMLMultiPole<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      Vec<3,Complex> n_cross_m = val.Range(0, 3);
      for (int k = 0; k < 3; k++)
      {
        Vec<3> ek{0.0}; ek(k) = 1;
        Vec<3> n_cross_m_real = Real(n_cross_m);
        Vec<3> n_cross_m_imag = Imag(n_cross_m);
        mp.AddDipole(pnt, Cross(n_cross_m_real, ek), ek);
        mp.AddDipole(pnt, Cross(n_cross_m_imag, ek), Complex(0,1)*ek);
      }
    }

    void AddSourceTrans(SingularMLMultiPole<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      Vec<3,Complex> n_cross_m = val.Range(0, 3);
      for (int k = 0; k < 3; k++)
      {
        Vec<3> ek{0.0}; ek(k) = 1;
        Vec<3> n_cross_m_real = Real(n_cross_m);
        Vec<3> n_cross_m_imag = Imag(n_cross_m);
        mp.AddDipole(pnt, Cross(n_cross_m_real, ek), ek);
        mp.AddDipole(pnt, Cross(n_cross_m_imag, ek), Complex(0,1)*ek);
      }
    }

    void EvaluateMP (RegularMLMultiPole<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val = mp.Evaluate (pnt);
    }

    void EvaluateMPTrans(RegularMLMultiPole<Vec<3,Complex>> & mp, Vec<3> pnt, Vec<3> nv, BareSliceVector<Complex> val) const
    {
      val = mp.Evaluate (pnt);
    }
  };

  
}



#endif

