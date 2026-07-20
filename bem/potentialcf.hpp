#ifndef NGBEM_POTENTIALCF_HPP
#define NGBEM_POTENTIALCF_HPP

#include <functional>
#include <type_traits>
#include <variant>

#include "../fem/integratorcf.hpp"
#include "../comp/fespace.hpp"

#include "intop_parameters.hpp"
#include "kernels.hpp"


namespace ngsbem
{
  using namespace ngcomp;
  using namespace ngla;

  template <typename T>
  struct IsStdVariant : std::false_type { };
  template <typename... Ts>
  struct IsStdVariant<std::variant<Ts...>> : std::true_type { };

  // BaseKernel::GetDifferentiatedKernel returns void; a real overload returns a kernel object.
  template <class K>
  struct HasDiffKernelOverload
    : std::integral_constant<bool,
        !std::is_void<decltype(std::declval<const K&>().GetDifferentiatedKernel(std::declval<const string&>()))>::value> { };

  class BasePotentialCF : public CoefficientFunctionNoDerivative
  {
  protected:
    shared_ptr<GridFunction> gf;
    VorB source_vb = BND;
    optional<Region> definedon;
    shared_ptr<DifferentialOperator> evaluator;
    IntOp_Parameters io_params;
  public:
    BasePotentialCF (shared_ptr<GridFunction> _gf,
                     VorB _source_vb,
                     optional<Region> _definedon,
                     shared_ptr<DifferentialOperator> _evaluator, bool is_complex)
      : CoefficientFunctionNoDerivative (_evaluator->Dim(), is_complex),
        gf(_gf), source_vb(_source_vb), definedon(_definedon), evaluator(_evaluator) { }

    BasePotentialCF (int dim, bool is_complex)
      : CoefficientFunctionNoDerivative (dim, is_complex) { }

    virtual ~BasePotentialCF() = default;

    virtual void BuildLocalExpansion(const Region & reg) = 0;
  };


  class SumPotentialCF : public BasePotentialCF
  {
    Array<tuple<Scalar, shared_ptr<BasePotentialCF>>> summands;
    shared_ptr<CoefficientFunction> summed_cf;

    static int GetDimension(const Array<tuple<Scalar, shared_ptr<BasePotentialCF>>> & asummands)
    {
      if (asummands.Size() == 0)
        throw Exception("SumPotentialCF needs at least one summand");
      return get<1>(asummands[0])->Dimension();
    }

    static bool IsComplexSum(const Array<tuple<Scalar, shared_ptr<BasePotentialCF>>> & asummands)
    {
      for (auto const & entry : asummands)
        {
          auto const & [scal, cf] = entry;
          if (cf->IsComplex() || std::holds_alternative<Complex>(scal))
            return true;
        }
      return false;
    }

    static shared_ptr<CoefficientFunction>
    BuildSummedCF(const Array<tuple<Scalar, shared_ptr<BasePotentialCF>>> & asummands)
    {
      shared_ptr<CoefficientFunction> sum;
      for (auto const & entry : asummands)
        {
          auto const & [scal, cf] = entry;
          shared_ptr<CoefficientFunction> scalecf = scal * static_pointer_cast<CoefficientFunction>(cf);
          if (sum)
            sum = sum + scalecf;
          else
            sum = scalecf;
        }
      return sum;
    }

  public:
    SumPotentialCF (Array<tuple<Scalar, shared_ptr<BasePotentialCF>>> asummands)
      : BasePotentialCF(GetDimension(asummands), IsComplexSum(asummands)),
        summands(asummands),
        summed_cf(BuildSummedCF(asummands))
    { ; }

    void BuildLocalExpansion(const Region & reg) override
    {
      for (auto [scal, cf] : summands)
        cf->BuildLocalExpansion(reg);
    }

    using CoefficientFunctionNoDerivative::Evaluate;
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("eval not implemented"); }

    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<> result) const override
    { summed_cf->Evaluate(ip, result); }

    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<Complex> result) const override
    { summed_cf->Evaluate(ip, result); }

    void Evaluate(const BaseMappedIntegrationRule & ir,
                  BareSliceMatrix<> result) const override
    { summed_cf->Evaluate(ir, result); }

    void Evaluate(const BaseMappedIntegrationRule & ir,
                  BareSliceMatrix<Complex> result) const override
    { summed_cf->Evaluate(ir, result); }

    void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   BareSliceMatrix<SIMD<double>> result) const override
    { summed_cf->Evaluate(ir, result); }

    void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   BareSliceMatrix<SIMD<Complex>> result) const override
    { summed_cf->Evaluate(ir, result); }
  };


  template <typename KERNEL>
  class PotentialCF : public BasePotentialCF
  {
    KERNEL kernel;
    int intorder;
    bool nearfield;

    using LOCAL_EXPANSION = typename std::invoke_result_t<decltype(&KERNEL::target_type::CreateLocalExpansion),typename KERNEL::target_type,Vec<3>,double,FMM_Parameters>;

    LOCAL_EXPANSION local_expansion;

  public:
    PotentialCF (shared_ptr<GridFunction> _gf,
                 VorB _source_vb,
                 optional<Region> _definedon,
                 shared_ptr<DifferentialOperator> _evaluator,
                 KERNEL _kernel, int _intorder, bool anearfield,
                 IntOp_Parameters _io_params = IntOp_Parameters());

    virtual shared_ptr<CoefficientFunction> Operator (const string & name) const override
    {
      if constexpr (HasDiffKernelOverload<KERNEL>::value)
      {
        auto diffkernel = kernel.GetDifferentiatedKernel(name);
        if constexpr (IsStdVariant<decltype(diffkernel)>::value)
          {
            return std::visit([&](auto const & dk) -> shared_ptr<CoefficientFunction>
            {
              using DK = std::decay_t<decltype(dk)>;
              return make_shared<PotentialCF<DK>>(this->gf, this->source_vb, this->definedon, this->evaluator, dk, intorder, nearfield, io_params);
            }, diffkernel);
          }
        else
          return make_shared<PotentialCF<decltype(diffkernel)>>(this->gf, this->source_vb, this->definedon, this->evaluator, diffkernel, intorder, nearfield, io_params);
      }
      else
        throw Exception("Kernel does not support differentiated kernel '"+name+"'");
    }

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
    void AddSourceElementContribution(const BaseMappedIntegrationPoint & mip,
                                      ElementId ei,
                                      const IntegrationRule & ir,
                                      FlatVector<T> result,
                                      T scale,
                                      LocalHeap & lh) const;
    template <typename T>
    void AddTangentCorrection(const BaseMappedIntegrationPoint & mip,
                              ElementId ei,
                              FlatVector<T> result,
                              LocalHeap & lh) const;
    template <typename T>
    void AddLocalExpansionNearfieldCorrection(const BaseMappedIntegrationRule & bmir,
                                              BareSliceMatrix<T> result) const;
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
}

#endif
