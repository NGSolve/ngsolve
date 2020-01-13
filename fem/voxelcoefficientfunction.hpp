#ifndef NGSOLVE_VOXELCOEFFICIENTFUNCTION_HPP
#define NGSOLVE_VOXELCOEFFICIENTFUNCTION_HPP

#include "fem.hpp"

namespace ngfem
{
  template<typename SCAL>
  class VoxelCoefficientFunction : public CoefficientFunctionNoDerivative
  {
    Array<double> start, end;
    Array<size_t> dim_vals;
    FlatArray<SCAL> values;
    bool linear;
  public:
    VoxelCoefficientFunction(const Array<double>& _start,
                             const Array<double>& _end,
                             const Array<size_t>& _dim_vals,
                             FlatArray<SCAL> _values,
                             bool _linear)
      : CoefficientFunctionNoDerivative(1, is_same_v<SCAL, Complex>),
        start(_start), end(_end), dim_vals(_dim_vals),
        values(_values), linear(_linear)
    { ; }

    using CoefficientFunctionNoDerivative::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint& ip) const override;
    Complex EvaluateComplex(const BaseMappedIntegrationPoint& ip) const override;

    void Evaluate(const BaseMappedIntegrationPoint& mip, FlatVector<Complex> values) const override;

  private:
    SCAL T_Evaluate(const BaseMappedIntegrationPoint& ip) const;
  };
} // namespace ngfem

#endif // NGSOLVE_VOXELCOEFFICIENTFUNCTION_HPP
