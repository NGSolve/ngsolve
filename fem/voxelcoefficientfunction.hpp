#ifndef NGSOLVE_VOXELCOEFFICIENTFUNCTION_HPP
#define NGSOLVE_VOXELCOEFFICIENTFUNCTION_HPP

#include "coefficient.hpp"

namespace ngfem
{
  template<typename SCAL>
  class VoxelCoefficientFunction : public CoefficientFunctionNoDerivative
  {
    Array<double> start, end;
    Array<size_t> dim_vals;
    Array<SCAL> values;
    bool linear;
    shared_ptr<CoefficientFunction> trafocf;
  public:
    VoxelCoefficientFunction(const Array<double>& _start,
                             const Array<double>& _end,
                             const Array<size_t>& _dim_vals,
                             Array<SCAL>&& _values,
                             bool _linear,
                             shared_ptr<CoefficientFunction> trafo=nullptr)
      : CoefficientFunctionNoDerivative(1, is_same_v<SCAL, Complex>),
        start(_start), end(_end), dim_vals(_dim_vals),
        values(std::move(_values)), linear(_linear), trafocf(trafo)
    { ; }

    using CoefficientFunctionNoDerivative::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint& ip) const override;
    Complex EvaluateComplex(const BaseMappedIntegrationPoint& ip) const override;

    void Evaluate(const BaseMappedIntegrationPoint& mip, FlatVector<Complex> values) const override;
    auto GetCArgs() const
    { return make_tuple(start, end, dim_vals, values, linear); }

  private:
    SCAL T_Evaluate(const BaseMappedIntegrationPoint& ip) const;
  };
} // namespace ngfem

#endif // NGSOLVE_VOXELCOEFFICIENTFUNCTION_HPP
