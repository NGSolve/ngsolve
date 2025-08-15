#ifndef NGBEM_hpp
#define NGBEM_hpp

#include "bem_diffops.hpp"

namespace ngsbem
{
  using namespace ngcomp;

  /** The IntegralOperator provides methods for the assembly of and access to a matrix 
      resulting from a variational formulation of a boundary integral equation.*/
  template <typename T = double>
  class IntegralOperator
  {
  protected:
    shared_ptr<FESpace> trial_space; 
    shared_ptr<FESpace> test_space; 

    optional<Region> trial_definedon;
    optional<Region> test_definedon;

    // integration order
    int intorder;

    /* boundary to global fe dofs mappings */
    Array<DofId> mapglob2bnd;
    Array<DofId> mapbnd2glob;
    Array<DofId> mapglob2bnd2;
    Array<DofId> mapbnd2glob2;

    Table<int> elems4dof; // contains list of elems contributing to bnd-dof
    Table<int> elems4dof2; // contains list of elems contributing to bnd-dof

    shared_ptr<BaseMatrix> matrix;


  public:

    /** Constructor. */
    IntegralOperator (shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                      optional<Region> _definedon_trial, optional<Region> _definedon_test,
                      int _intorder);
    virtual ~IntegralOperator() = default;

    shared_ptr<BaseMatrix> GetMatrix() const { return matrix; }

    virtual shared_ptr<BaseMatrix> CreateMatrixFMM(LocalHeap & lh) const = 0;
    
    virtual void CalcElementMatrix(FlatMatrix<T> matrix,
                                   ElementId ei_trial, ElementId ei_test,
                                   LocalHeap &lh) const = 0;
    
    /** CalcFarFieldBlock computes a low-rank approximation of block with trialdofs and testdofs. */
    /*
    virtual unique_ptr<LowRankMatrix<T>>
    CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
                      LocalHeap &lh) const;
    */
    virtual shared_ptr<CoefficientFunction> GetPotential(shared_ptr<GridFunction> gf,
                                                         optional<int> io, bool nearfield_experimental) const = 0;
  };


  
  /** The GenericIntegralOperator is a templated #IntegralOperator, the template type is 
      the kernel the specific potential,i.e. a fundamental solution or its 
      derivative of specific pde. */
  template <typename KERNEL>
  class GenericIntegralOperator : public IntegralOperator<typename KERNEL::value_type>
  {
    KERNEL kernel;
    typedef typename KERNEL::value_type value_type;
    typedef IntegralOperator<typename KERNEL::value_type> BASE;
    
    using BASE::trial_space; 
    using BASE::test_space;

    using BASE::trial_definedon; 
    using BASE::test_definedon;
       
    using BASE::intorder;

    using BASE::mapglob2bnd;
    using BASE::mapbnd2glob;
    using BASE::mapglob2bnd2;
    using BASE::mapbnd2glob2;

    using BASE::elems4dof; 
    using BASE::elems4dof2; 
    
    // using BASE::trial_ct; 
    // using BASE::test_ct;
    
    using BASE::matrix;

    shared_ptr<DifferentialOperator> trial_evaluator;
    shared_ptr<DifferentialOperator> test_evaluator;


    Array<Vec<2>> identic_panel_x, identic_panel_y;
    Array<double> identic_panel_weight;
    
    Array<Vec<2>> common_vertex_x, common_vertex_y;
    Array<double> common_vertex_weight;

    Array<Vec<2>> common_edge_x, common_edge_y;
    Array<double> common_edge_weight;
    
    
  public:
    GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                            optional<Region> _definedon_trial, optional<Region> _definedon_test,
                            shared_ptr<DifferentialOperator> _trial_evaluator, 
                            shared_ptr<DifferentialOperator> _test_evaluator, 
                            KERNEL _kernel,
                            int _intorder);

    GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                            shared_ptr<DifferentialOperator> _trial_evaluator, 
                            shared_ptr<DifferentialOperator> _test_evaluator, 
                            KERNEL _kernel,
                            int _intorder)
      : GenericIntegralOperator (_trial_space, _test_space,
                                 nullopt, nullopt,
                                 _trial_evaluator, _test_evaluator, _kernel, _intorder) { }

    GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                            KERNEL _kernel,
                            int _intorder)
      : GenericIntegralOperator (_trial_space, _test_space,
                                 _trial_space -> GetEvaluator(BND),
                                 _test_space -> GetEvaluator(BND),
                                 _kernel, _intorder) { }


    void CalcElementMatrix(FlatMatrix<value_type> matrix,
                           ElementId ei_trial, ElementId ei_test,
                           LocalHeap &lh) const override;

    /*
    unique_ptr<LowRankMatrix<value_type>>
    CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
                      LocalHeap &lh) const override;
    */

    
    shared_ptr<BaseMatrix> CreateMatrixFMM(LocalHeap & lh) const override;
    
    virtual shared_ptr<CoefficientFunction> GetPotential(shared_ptr<GridFunction> gf,
                                                         optional<int> io, bool nearfield_experimental) const override;
  };



  template  <typename KERNEL>
  class PotentialCF : public CoefficientFunctionNoDerivative
  {
    shared_ptr<GridFunction> gf;
    optional<Region> definedon;
    shared_ptr<DifferentialOperator> evaluator;
    KERNEL kernel;
    int intorder;
    bool nearfield;
  public:
    PotentialCF (shared_ptr<GridFunction> _gf,
                 optional<Region> _definedon,    
                 shared_ptr<DifferentialOperator> _evaluator,
                 KERNEL _kernel, int _intorder, bool anearfield);

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
  

  struct KernelTerm
  {
    double fac;
    size_t kernel_comp;
    size_t trial_comp;
    size_t test_comp;
  };
  

  /** LaplaceSLkernel is the kernel for the single layer potential of 
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM> class LaplaceSLKernel;

  /** LaplaceSLkernel in 3D reads 
      $$ G(x-y) = \frac{1}{4\,\pi \, | x-y| }, \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class LaplaceSLKernel<3> 
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
  };


  /** LaplaceDLkernel is the kernel for the double layer potential of 
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM> class LaplaceDLKernel;

  /** LaplaceDLkernel in 3D reads 
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, 
          \frac{ \langle n(y), x-y\rangle }{ | x-y|^3 }, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class LaplaceDLKernel<3> 
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
  };


  template <int DIM> class LaplaceHSKernel;
  
  template<>
  class LaplaceHSKernel<3> 
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
  };


  /** HelmholtzSLkernel is the kernel for the double layer potential of the 
      Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,. $ */
  template <int DIM> class HelmholtzSLKernel;

  /** HelmholtzSLkernel in 3D reads 
      $$ G(x-y) = \frac{1 }{4\,\pi} \,\frac{e^{i\,\kappa \, |x-y| }{|x-y|} \, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class HelmholtzSLKernel<3> 
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
  };


  /** HelmholtzDLkernel is the kernel for the double layer potential of 
      the Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,.$ */
  template <int DIM> class HelmholtzDLKernel;

  /** HelmholtzDLkernel in 3D reads
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \langle n(y), x-y\rangle \cdot \left( 1 - i\,\kappa\, | x-y| \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class HelmholtzDLKernel<3> 
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
  };




  template <int DIM> class HelmholtzHSKernel;
  
  template<>
  class HelmholtzHSKernel<3> 
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
  };




  

  /** CombinedFieldKernel is a kernel for the combined field integral equation 
      is considered for the Helmholtz equation. */
  template <int DIM> class CombinedFieldKernel;

  /** CombinedFieldKernel in 3D reads
      $$ G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \left( \langle n_y, x-y\rangle (1- i\,\kappa\, | x-y|) - i\,\kappa\,|x-y|^2 \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class CombinedFieldKernel<3> 
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
  };

  

  

  template <int D> class MaxwellSLKernel;

  template<>
  class MaxwellSLKernel<3> 
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
  };




  

  template <int D> class MaxwellDLKernel;

  /** MaxwellDLkernel for 3D in matrix representation reads
      $$  \left( \begin{array}{ccc} 0 & -\frac{\partial G_\kappa(x-y)}{\partial x_3} & \frac{\partial G_\kappa(x-y)}{\partial x_2} \\ \frac{\partial G_\kappa(x-y)}{\partial x_3} & 0 & -\frac{\partial G_\kappa(x-y)}{\partial x_1} \\ -\frac{\partial G_\kappa(x-y)}{\partial x_2} & \frac{\partial G_\kappa(x-y)}{\partial x_1} & 0 \end{array}\right)\,,$$ with 
   $$ G_\kappa(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \langle n(y), x-y\rangle \cdot \left( 1 - i\,\kappa\, | x-y| \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class MaxwellDLKernel<3> 
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
  };

  
}



#endif

