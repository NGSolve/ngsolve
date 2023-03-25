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

  private:
    SCAL T_Evaluate(const BaseMappedIntegrationPoint& ip) const;
  };

  template<typename SCAL>
  class VectorialVoxelCoefficientFunction : public CoefficientFunctionNoDerivative
  {
	  // Array of VoxelCoefficientFunction
  	  Array<shared_ptr<VoxelCoefficientFunction<SCAL>>> voxel_cfs;
	  Array<size_t> dimi;  // shape of entries

  public:
    VectorialVoxelCoefficientFunction() = default;
    VectorialVoxelCoefficientFunction (Array<shared_ptr<VoxelCoefficientFunction<SCAL>>> avoxel_cfs) : voxel_cfs(avoxel_cfs), dimi(avoxel_cfs.Size())
    {
    	for (auto cf : voxel_cfs)
      	  if (cf && cf->IsComplex())
            is_complex = true;
    }

    using CoefficientFunctionNoDerivative::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint& ip) const override;
    Complex EvaluateComplex(const BaseMappedIntegrationPoint& ip) const override;
    void Evaluate(const BaseMappedIntegrationPoint& mip, FlatVector<> values) const override;
    void Evaluate(const BaseMappedIntegrationPoint& mip, FlatVector<Complex> values) const override;

  };

} // namespace ngfem

#endif // NGSOLVE_VOXELCOEFFICIENTFUNCTION_HPP
