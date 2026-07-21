#ifndef NGBEM_hpp
#define NGBEM_hpp

#include "../fem/integratorcf.hpp"
#include "../linalg/basematrix.hpp"
#include "../linalg/sparsematrix.hpp"
#include "../comp/fespace.hpp"
#include <mutex>


#include "bem_diffops.hpp"
#include "diffopwithfactor.hpp"
#include "kernels.hpp"
#include "fmmoperator.hpp"
#include "potentialcf.hpp"


namespace ngsbem
{
  using namespace ngcomp;
  using namespace ngla;

  
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
    shared_ptr<DifferentialOperator> test_evaluator;
    
    mutable unique_ptr<Table<tuple<ElementId,int>,DofId>> trial_dof2el;
    mutable unique_ptr<Table<tuple<ElementId,int>,DofId>> test_dof2el;
    mutable std::mutex matrix_mutex;
    mutable std::mutex nearfield_mutex;
    
    // integration order
    int intorder;
    IntOp_Parameters io_params;

    // Sauter-Schwab integration rules:
    Array<Vec<2>> identic_panel_x, identic_panel_y;
    Array<double> identic_panel_weight;

    Array<Vec<2>> identic_panel_quad_x, identic_panel_quad_y;
    Array<double> identic_panel_quad_weight;
    
    Array<Vec<2>> common_vertex_x, common_vertex_y;
    Array<double> common_vertex_weight;

    Array<Vec<2>> common_vertex_quad_x, common_vertex_quad_y;
    Array<double> common_vertex_quad_weight;

    Array<Vec<2>> common_vertex_quadtrig_x, common_vertex_quadtrig_y;
    Array<double> common_vertex_quadtrig_weight;

    Array<Vec<2>> common_edge_x, common_edge_y;
    Array<double> common_edge_weight;
    
    Array<Vec<2>> common_edge_quad_x, common_edge_quad_y;
    Array<double> common_edge_quad_weight;
    
    Array<Vec<2>> common_edge_quadtrig_x, common_edge_quadtrig_y;
    Array<double> common_edge_quadtrig_weight;
    


    
    mutable shared_ptr<BaseMatrix> matrix;
    mutable shared_ptr<BaseMatrix> nearfield_matrix;

  public:

    IntegralOperator (shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                      optional<Region> _definedon_trial, optional<Region> _definedon_test,
                      shared_ptr<DifferentialOperator> _trial_evaluator, //  shared_ptr<CoefficientFunction> _trial_factor,
                      shared_ptr<DifferentialOperator> _test_evaluator, // shared_ptr<CoefficientFunction> _test_factor,
                      int _intorder, const IntOp_Parameters & _io_params);
    
    virtual ~IntegralOperator() = default;

    shared_ptr<BaseMatrix> GetMatrix() const
    {
      if (!matrix)
        {
          std::lock_guard<std::mutex> guard(matrix_mutex);
          if (!matrix)
            {
              LocalHeap lh(100*1000*1000);
              matrix = CreateMatrixFMM(lh);
            }
        }
      return matrix;
    }
    shared_ptr<FESpace> GetTrialSpace() const { return trial_space; }
    shared_ptr<FESpace> GetTestSpace() const { return test_space; }
    optional<Region> GetTrialDefinedOn() const { return trial_definedon; }
    optional<Region> GetTestDefinedOn() const { return test_definedon; }
    shared_ptr<DifferentialOperator> GetTrialEvaluator() const { return trial_evaluator; }
    shared_ptr<DifferentialOperator> GetTestEvaluator() const { return test_evaluator; }
    int GetIntOrder() const { return intorder; }
    IntOp_Parameters GetIOParams() const { return io_params; }

    virtual shared_ptr<BaseMatrix> CreateMatrixFMM(LocalHeap & lh) const = 0;
    virtual shared_ptr<BaseMatrix> CreateNearFieldMatrix(LocalHeap & lh) const = 0;

    virtual shared_ptr<BasePotentialCF> GetPotential(shared_ptr<GridFunction> gf,
                                                     optional<int> io, bool nearfield_experimental) const = 0;

    shared_ptr<BaseMatrix> GetNearFieldMatrix() const
    {
      if (!nearfield_matrix)
        {
          std::lock_guard<std::mutex> guard(nearfield_mutex);
          if (!nearfield_matrix)
            {
              LocalHeap lh(100*1000*1000);
              nearfield_matrix = CreateNearFieldMatrix(lh);
            }
        }
      return nearfield_matrix;
    }

    virtual std::variant<Matrix<double>, Matrix<Complex>> CalcSubMatrix (FlatArray<DofId> rowids, FlatArray<DofId> colids, LocalHeap &lh) const = 0;

    virtual FMMOperatorInfo GetFMMInfo() const
    {
      throw Exception("GetFMMInfo not implemented for this IntegralOperator");
    }
  };

  template <typename TSCAL>
  inline std::variant<Matrix<double>, Matrix<Complex>>
  ScaleDenseMatrixVariant(std::variant<Matrix<double>, Matrix<Complex>> mat, TSCAL fac)
  {
    return std::visit([&](auto const & inmat) -> std::variant<Matrix<double>, Matrix<Complex>>
    {
      using MAT = std::decay_t<decltype(inmat)>;
      if constexpr (std::is_same_v<MAT, Matrix<double>> && std::is_same_v<TSCAL, double>)
        {
          Matrix<double> scaled(inmat);
          scaled *= fac;
          return scaled;
        }
      else
        {
          Matrix<Complex> scaled(inmat.Height(), inmat.Width());
          scaled = fac * inmat;
          return scaled;
        }
    }, std::move(mat));
  }

  inline std::variant<Matrix<double>, Matrix<Complex>>
  AddDenseMatrixVariants(std::variant<Matrix<double>, Matrix<Complex>> a,
                         std::variant<Matrix<double>, Matrix<Complex>> b)
  {
    return std::visit([&](auto const & mata, auto const & matb) -> std::variant<Matrix<double>, Matrix<Complex>>
    {
      using MATA = std::decay_t<decltype(mata)>;
      using MATB = std::decay_t<decltype(matb)>;
      if constexpr (std::is_same_v<MATA, Matrix<double>> && std::is_same_v<MATB, Matrix<double>>)
        {
          Matrix<double> sum(mata);
          sum += matb;
          return sum;
        }
      else
        {
          Matrix<Complex> sum(mata.Height(), mata.Width());
          sum = mata + matb;
          return sum;
        }
    }, std::move(a), std::move(b));
  }

  class SumIntegralOperator : public IntegralOperator
  {
    Array<shared_ptr<IntegralOperator>> summands;

    static shared_ptr<IntegralOperator>
    GetReferenceOperator(const Array<shared_ptr<IntegralOperator>> & asummands)
    {
      if (asummands.Size() == 0)
        throw Exception("SumIntegralOperator needs at least one summand");
      return asummands[0];
    }
  public:
    SumIntegralOperator(Array<shared_ptr<IntegralOperator>> asummands)
      : IntegralOperator(GetReferenceOperator(asummands)->GetTrialSpace(),
                         GetReferenceOperator(asummands)->GetTestSpace(),
                         GetReferenceOperator(asummands)->GetTrialDefinedOn(),
                         GetReferenceOperator(asummands)->GetTestDefinedOn(),
                         GetReferenceOperator(asummands)->GetTrialEvaluator(),
                         GetReferenceOperator(asummands)->GetTestEvaluator(),
                         GetReferenceOperator(asummands)->GetIntOrder(),
                         GetReferenceOperator(asummands)->GetIOParams()),
        summands(std::move(asummands))
    { ; }

    auto const & Summands() const { return summands; }

    shared_ptr<BaseMatrix> CreateMatrixFMM(LocalHeap & lh) const override
    {
      shared_ptr<BaseMatrix> sum;
      for (auto const & op : summands)
        {
          auto mat = op->GetMatrix();
          if (sum)
            sum = sum + mat;
          else
            sum = mat;
        }
      return sum;
    }

    shared_ptr<BasePotentialCF> GetPotential(shared_ptr<GridFunction> gf,
                                             optional<int> io, bool nearfield_experimental) const override
    {
      throw Exception("GetPotential not implemented for SumIntegralOperator");
    }

    shared_ptr<BaseMatrix> CreateNearFieldMatrix(LocalHeap & lh) const override
    {
      shared_ptr<BaseMatrix> sum;
      for (auto const & op : summands)
        {
          auto mat = op->GetNearFieldMatrix();
          if (sum)
            sum = sum + mat;
          else
            sum = mat;
        }
      return sum;
    }

    virtual std::variant<Matrix<double>, Matrix<Complex>> CalcSubMatrix (FlatArray<DofId> rowids, FlatArray<DofId> colids, LocalHeap &lh) const override
    {
      optional<std::variant<Matrix<double>, Matrix<Complex>>> sum;
      for (auto const & op : summands)
        {
          auto mat = op->CalcSubMatrix(rowids, colids, lh);
          if (sum)
            sum = AddDenseMatrixVariants(std::move(*sum), std::move(mat));
          else
            sum = std::move(mat);
        }
      if (!sum)
        throw Exception("SumIntegralOperator needs at least one summand");
      return std::move(*sum);
    }

    FMMOperatorInfo GetFMMInfo() const override
    {
      if (summands.Size() != 1)
        throw Exception("GetFMMInfo for SumIntegralOperator with multiple summands is not implemented");
      return summands[0]->GetFMMInfo();
    }
  };

  template <typename TSCAL>
  class ScaledIntegralOperator : public IntegralOperator
  {
    shared_ptr<IntegralOperator> op;
    TSCAL fac;
  public:
    ScaledIntegralOperator(shared_ptr<IntegralOperator> _op, TSCAL _fac)
      : IntegralOperator(_op->GetTrialSpace(), _op->GetTestSpace(),
                         _op->GetTrialDefinedOn(), _op->GetTestDefinedOn(),
                         _op->GetTrialEvaluator(), _op->GetTestEvaluator(),
                         _op->GetIntOrder(), _op->GetIOParams()),
        op(std::move(_op)), fac(_fac)
    { ; }

    shared_ptr<BaseMatrix> CreateMatrixFMM(LocalHeap & lh) const override
    {
      return make_shared<VScaleMatrix<TSCAL>>(op->GetMatrix(), fac);
    }

    shared_ptr<BaseMatrix> CreateNearFieldMatrix(LocalHeap & lh) const override
    {
      return make_shared<VScaleMatrix<TSCAL>>(op->GetNearFieldMatrix(), fac);
    }

    shared_ptr<BasePotentialCF> GetPotential(shared_ptr<GridFunction> gf,
                                             optional<int> io, bool nearfield_experimental) const override
    {
      throw Exception("GetPotential not implemented for ScaledIntegralOperator");
    }

    std::variant<Matrix<double>, Matrix<Complex>> CalcSubMatrix (FlatArray<DofId> rowids, FlatArray<DofId> colids, LocalHeap &lh) const override
    {
      return ScaleDenseMatrixVariant(op->CalcSubMatrix(rowids, colids, lh), fac);
    }

    FMMOperatorInfo GetFMMInfo() const override { return op->GetFMMInfo(); }
  };

  template <typename TSCAL>
  inline shared_ptr<IntegralOperator>
  ScaleIntegralOperator(shared_ptr<IntegralOperator> op, TSCAL fac)
  {
    return make_shared<ScaledIntegralOperator<TSCAL>>(std::move(op), fac);
  }

  inline shared_ptr<IntegralOperator>
  AddIntegralOperators(shared_ptr<IntegralOperator> opa, shared_ptr<IntegralOperator> opb)
  {
    Array<shared_ptr<IntegralOperator>> summands;
    if (auto suma = dynamic_pointer_cast<SumIntegralOperator>(opa))
      summands += suma->Summands();
    else
      summands.Append(opa);
    if (auto sumb = dynamic_pointer_cast<SumIntegralOperator>(opb))
      summands += sumb->Summands();
    else
      summands.Append(opb);
    return make_shared<SumIntegralOperator>(std::move(summands));
  }


  
  /** The GenericIntegralOperator is a templated #IntegralOperator, the template type is 
      the kernel the specific potential,i.e. a fundamental solution or its 
      derivative of specific pde. */
  template <typename KERNEL>
  class GenericIntegralOperator : public IntegralOperator // <typename KERNEL::value_type>
  {
    KERNEL kernel;
    typedef typename KERNEL::value_type value_type;
    typedef IntegralOperator BASE;

    static size_t CountSparseCorrectionNZE(const BaseMatrix * mat)
    {
      if (!mat) return 0;
      if (auto sparse = dynamic_cast<const BaseSparseMatrix*>(mat))
        return sparse->NZE();
      if (auto sum = dynamic_cast<const SumMatrix*>(mat))
        {
          size_t nze = 0;
          nze += CountSparseCorrectionNZE(sum->SPtrA().get());
          nze += CountSparseCorrectionNZE(sum->SPtrB().get());
          return nze;
        }
      return 0;
    }
    
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
    shared_ptr<BaseMatrix> CreateNearFieldMatrix(LocalHeap & lh) const override;
    virtual std::variant<Matrix<double>, Matrix<Complex>> CalcSubMatrix (FlatArray<DofId> rowids, FlatArray<DofId> colids, LocalHeap &lh) const override;

    virtual shared_ptr<BasePotentialCF> GetPotential(shared_ptr<GridFunction> gf,
                                                         optional<int> io, bool nearfield_experimental) const override;

    FMMOperatorInfo GetFMMInfo() const override
    {
      auto mat = this->GetMatrix();
      auto fmm = FindFMMOperator<value_type>(mat.get());
      if (!fmm)
        throw Exception("GetFMMInfo: no FMM operator found inside matrix");
      auto info = fmm->GetFMMInfo();
      info.source_dofs = trial_space->GetNDof();
      info.target_dofs = test_space->GetNDof();
      info.nearfield_nze = CountSparseCorrectionNZE(mat.get());
      double total_dof_pairs = double(info.source_dofs) * double(info.target_dofs);
      info.nearfield_fraction = (total_dof_pairs > 0) ? double(info.nearfield_nze) / total_dof_pairs : 0.0;
      return info;
    }
  };




  class BasePotentialOperator
  {
  public:
    shared_ptr<ProxyFunction> proxy;
    // shared_ptr<CoefficientFunction> factor;
    VorB source_vb = BND;
    optional<Region> definedon;
    shared_ptr<DifferentialOperator> evaluator;
    IntOp_Parameters io_params;
    int intorder;
    
  public:
    BasePotentialOperator (shared_ptr<ProxyFunction> _proxy, // shared_ptr<CoefficientFunction> _factor,
                           VorB _source_vb,
                           optional<Region> _definedon,    
                           shared_ptr<DifferentialOperator> _evaluator,
                           int _intorder)
      : proxy(_proxy), /* factor(_factor), */ source_vb(_source_vb), definedon(_definedon), evaluator(_evaluator), intorder(_intorder) { ; }

    BasePotentialOperator (shared_ptr<ProxyFunction> _proxy, // shared_ptr<CoefficientFunction> _factor,
                           VorB _source_vb,
                           optional<Region> _definedon,    
                           shared_ptr<DifferentialOperator> _evaluator,
                           IntOp_Parameters _io_params, 
                           int _intorder)
      : proxy(_proxy), source_vb(_source_vb), definedon(_definedon), evaluator(_evaluator), io_params(_io_params), intorder(_intorder) { ; }



    
    virtual ~BasePotentialOperator() { } 
    virtual shared_ptr<IntegralOperator> MakeIntegralOperator(shared_ptr<ProxyFunction> test_proxy,
                                                              // shared_ptr<CoefficientFunction> test_factor,
                                                              DifferentialSymbol dx) = 0;
    virtual shared_ptr<BasePotentialCF> MakePotentialCF(shared_ptr<GridFunction> gf) = 0;
    virtual shared_ptr<BasePotentialOperator> MakeDiffBasePotential(string name) = 0;

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
                       VorB _source_vb,
                       optional<Region> _definedon,    
                       shared_ptr<DifferentialOperator> _evaluator,
                       KERNEL _kernel, int _intorder)
      : BasePotentialOperator (_proxy, /* _factor, */ _source_vb, _definedon, _evaluator, _intorder), kernel(_kernel) { ; }

    PotentialOperator (shared_ptr<ProxyFunction> _proxy,
                       // shared_ptr<CoefficientFunction> _factor,                       
                       VorB _source_vb,
                       optional<Region> _definedon,    
                       shared_ptr<DifferentialOperator> _evaluator,
                       KERNEL _kernel,
                       IntOp_Parameters _io_params, 
                       int _intorder)
      : BasePotentialOperator (_proxy, _source_vb, _definedon, _evaluator, _io_params, _intorder), kernel(_kernel) { ; }

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
      return make_shared<PotentialCF<KERNEL>>(gf, source_vb, definedon, evaluator, kernel, 2+intorder, true, io_params);
    }

    shared_ptr<BasePotentialOperator> MakeDiffBasePotential(string name) override
    {
      if constexpr (HasDiffKernelOverload<KERNEL>::value)
      {
        auto diffkernel = kernel.GetDifferentiatedKernel(name);
        if constexpr (IsStdVariant<decltype(diffkernel)>::value)
        {
          return std::visit([&](auto const & dk) -> shared_ptr<BasePotentialOperator>
          {
            using DK = std::decay_t<decltype(dk)>;
            return std::make_shared<PotentialOperator<DK>>(this->proxy, this->source_vb, this->definedon, this->evaluator, dk, this->io_params, this->intorder);
          }, diffkernel);
        }
        else
          return make_shared<PotentialOperator<decltype(diffkernel)>>(this->proxy, this->source_vb, this->definedon, this->evaluator, diffkernel, this->io_params, this->intorder);
      }
      else
        throw Exception("Kernel does not support differentiated kernel '"+name+"'");
    }

    virtual void Print (ostream & ost) const override
    {
      ost << "Potential operator:" << endl;
      ost << "Kernel: " << kernel.Name() << endl;
      ost << io_params << endl;
    }
  };  

  
  class SumOfPotentialOperators
  {
    Array<tuple<Scalar, shared_ptr<BasePotentialOperator>>> summands;
  public:
    SumOfPotentialOperators (Scalar fac, shared_ptr<BasePotentialOperator> op)
    { summands += tuple(fac, op); }
    SumOfPotentialOperators (Array<tuple<Scalar, shared_ptr<BasePotentialOperator>>> asummands)
      : summands(asummands) { }
    auto & Summands() const { return summands; }
    
    shared_ptr<BasePotentialCF> MakePotentialCF(shared_ptr<GridFunction> gf)
    {
      Array<tuple<Scalar, shared_ptr<BasePotentialCF>>> cfs;
      for (auto [scal, potop] : summands)
        cfs.Append(tuple(scal, potop->MakePotentialCF(gf)));
      return make_shared<SumPotentialCF>(std::move(cfs));
    }

    shared_ptr<BasePotentialCF> MakePotentialCF(shared_ptr<GridFunction> gf, const Region & region)
    {
      auto sum = MakePotentialCF(gf);
      sum->BuildLocalExpansion(region);
      return sum;
    }

    
  };

  inline SumOfPotentialOperators operator+ (const SumOfPotentialOperators & sum1,
                                            const SumOfPotentialOperators & sum2)
    
  {
    Array<tuple<Scalar, shared_ptr<BasePotentialOperator>>> res{sum1.Summands()};
    res += sum2.Summands();
    return SumOfPotentialOperators(std::move(res));
  }

  inline SumOfPotentialOperators operator- (const SumOfPotentialOperators & sum)
  {
    Array<tuple<Scalar, shared_ptr<BasePotentialOperator>>> res;
    for (auto [scal,potop] : sum.Summands())
      res += tuple(-scal, potop);
    return SumOfPotentialOperators(std::move(res));
  }

  inline SumOfPotentialOperators operator- (const SumOfPotentialOperators & sum1,
                                            const SumOfPotentialOperators & sum2)
  {
    return sum1 + (-sum2);
  }
  
  

  inline Array < tuple <shared_ptr<ProxyFunction>, shared_ptr<CoefficientFunction>  >>
  CreateProxyLinearization (shared_ptr<CoefficientFunction> cf, bool trial, VorB vb)
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
      if (!proxy->Evaluator()->SupportsVB(vb))
        throw Exception ("Proxy does not support requested integration domain");

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
  GetProxyAndFactor (shared_ptr<CoefficientFunction> cf, bool trial, VorB vb = BND)
  {
    if (auto proxy = dynamic_pointer_cast<ProxyFunction>(cf))
      return { proxy, nullptr };
    
    auto proxylin = CreateProxyLinearization (cf, trial, vb);
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
  GetProxyWithFactor (shared_ptr<CoefficientFunction> cf, bool trial, VorB vb = BND)
  {
    return std::get<0> (GetProxyAndFactor(cf, trial, vb));
  }

  class BasePotentialOperatorAndTest
  {
    shared_ptr<BasePotentialOperator> pot;
    shared_ptr<ProxyFunction> test_proxy_with_factor;    
  public:

    BasePotentialOperatorAndTest (shared_ptr<BasePotentialOperator> _pot,
                                  shared_ptr<CoefficientFunction> _test_proxy)
      : pot(_pot)
    {
      test_proxy_with_factor = GetProxyWithFactor(_test_proxy, false);      
    }

    
    shared_ptr<IntegralOperator> MakeIntegralOperator (DifferentialSymbol dx)
    {
      return pot->MakeIntegralOperator(test_proxy_with_factor, dx);
    }
  };

  class SumOfPotentialOperatorsAndTest
  {
    Array<tuple<Scalar, shared_ptr<BasePotentialOperator>>> pots;
    shared_ptr<ProxyFunction> test_proxy_with_factor;
  public:
    SumOfPotentialOperatorsAndTest (const Array<tuple<Scalar, shared_ptr<BasePotentialOperator>>> & _pots,
                                    shared_ptr<CoefficientFunction> _test_proxy)
      : pots(_pots)
    {
      test_proxy_with_factor = GetProxyWithFactor(_test_proxy, false);
    }

    shared_ptr<IntegralOperator> MakeIntegralOperator (DifferentialSymbol dx)
    {
      Array<shared_ptr<IntegralOperator>> ops;
      for (auto [scal, pot] : pots)
        {
          auto op = pot->MakeIntegralOperator(test_proxy_with_factor, dx);
          if (std::holds_alternative<Complex>(scal))
            ops.Append(ScaleIntegralOperator(op, std::get<Complex>(scal)));
          else
            ops.Append(ScaleIntegralOperator(op, std::get<double>(scal)));
        }
      return make_shared<SumIntegralOperator>(std::move(ops));
    }
  };



}



#endif
