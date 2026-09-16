#ifndef FILE_INTEGRALKERNEL
#define FILE_INTEGRALKERNEL

#include <fem.hpp>
#include "fmminterface.hpp"

namespace ngla { class BaseMatrix; }

namespace ngsbem
{
  using namespace ngfem;

  struct KernelTerm
  {
    double fac;
    size_t kernel_comp;
    size_t trial_comp;
    size_t test_comp;
  };

  enum class AnalyticTriangleFormula
  {
    none,
    laplace_sl,
    laplace_dl,
    laplace_grad_sl
  };

  // Keep element assembly independent of the kernel and FMM coefficient types.
  // Each virtual evaluation handles a whole quadrature batch, not one point.
  template <typename TSCAL>
  class BaseIntegralKernel
  {
  public:
    virtual ~BaseIntegralKernel() = default;
    virtual string Name() const = 0;
    virtual IVec<2> Shape() const = 0;
    virtual AnalyticTriangleFormula GetAnalyticTriangleFormula() const = 0;
    virtual const BaseFMMInterface & Source() const = 0;
    virtual const BaseFMMInterface & Target() const = 0;
    virtual shared_ptr<const BaseIntegralKernel<TSCAL>> GetDifferentiatedKernel(const string & name) const = 0;
    virtual size_t NumKernelComponents() const = 0;
    virtual FlatArray<const KernelTerm> Terms() const = 0;

    // Component-major values, including the paired quadrature weights.
    virtual void EvaluatePairs(const SIMD_BaseMappedIntegrationRule & mirx,
                               const SIMD_BaseMappedIntegrationRule & miry,
                               FlatMatrix<SIMD<TSCAL>> values) const = 0;

    // Weighted Cartesian product of two ordinary quadrature rules.
    virtual void EvaluateMatrix(const BaseMappedIntegrationRule & mirx,
                                const BaseMappedIntegrationRule & miry,
                                size_t component, double factor,
                                FlatMatrix<TSCAL> values, bool skip_diagonal = false) const = 0;

    // Accumulate one target's weighted source contribution without pointwise dispatch.
    virtual void AddPotential(const BaseMappedIntegrationPoint & mip,
                              const SIMD_BaseMappedIntegrationRule & miry,
                              FlatMatrix<SIMD<TSCAL>> values,
                              FlatVector<SIMD<TSCAL>> result, VorB source_vb) const = 0;

    virtual shared_ptr<ngla::BaseMatrix> CreateFMMOperator(Array<Vec<3>> xpts, Array<Vec<3>> ypts,
                                                         Array<Vec<3>> xnv, Array<Vec<3>> ynv,
                                                         const FMM_Parameters & params) const = 0;
  };

  template <typename KERNEL>
  NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<typename KERNEL::value_type>> MakeIntegralKernel(KERNEL kernel);
}

#endif
