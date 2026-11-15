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
    {
      int dim = start.Size();
      if (dim < 1 || dim > 3)
        throw Exception("VoxelCoefficient: start must have 1, 2 or 3 entries, got " + ToString(dim));
      if (end.Size() != dim)
        throw Exception("VoxelCoefficient: start has " + ToString(dim) + " entries, end has " + ToString(end.Size()));
      if (dim_vals.Size() != dim)
        throw Exception("VoxelCoefficient: values must be a " + ToString(dim) + "-dimensional array, got "
                        + ToString(dim_vals.Size()) + " dimensions");
      size_t prod = 1;
      for (int i = 0; i < dim; i++)
        {
          if (dim_vals[i] == 0)
            throw Exception("VoxelCoefficient: values must not be empty");
          if (!(start[i] < end[i]))
            throw Exception("VoxelCoefficient: need start < end in every direction");
          prod *= dim_vals[i];
        }
      if (prod != values.Size())
        throw Exception("VoxelCoefficient: number of values does not match the shape");
      if (trafocf && trafocf->Dimension() != dim)
        throw Exception("VoxelCoefficient: trafocf must have dimension " + ToString(dim) + ", got "
                        + ToString(trafocf->Dimension()));
    }

    void CheckContext (const CFContext & ctx) const override
    {
      if (!trafocf && ctx.dim >= 0 && ctx.dim < int(start.Size()))
        throw Exception("VoxelCoefficient of dimension " + ToString(start.Size()) + " evaluated on a mesh of dimension "
                        + ToString(ctx.dim) + ", use trafocf");
    }

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
