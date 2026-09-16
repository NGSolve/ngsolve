#include "integralkernel.hpp"
#include "kernels.hpp"
#include "fmmoperator.hpp"

namespace ngsbem
{
  template <typename T>
  struct IsStdVariant : std::false_type { };
  template <typename... Ts>
  struct IsStdVariant<std::variant<Ts...>> : std::true_type { };

  // Construct expansions here so their virtual methods are not emitted by every caller.
  template <typename T_MP, typename T_VALUE, typename T_Kappa>
  shared_ptr<BaseSingularMLExpansion> FMMInterface<T_MP,T_VALUE,T_Kappa> ::
  CreateMultipoleExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
  {
    return make_shared<SingularMLExpansion<T_MP,T_Kappa>> (c, r, kappa, fmm_params);
  }

  template <typename T_MP, typename T_VALUE, typename T_Kappa>
  shared_ptr<BaseRegularMLExpansion> FMMInterface<T_MP,T_VALUE,T_Kappa> ::
  CreateLocalExpansion (Vec<3> c, double r, FMM_Parameters fmm_params) const
  {
    return make_shared<RegularMLExpansion<T_MP,T_Kappa>> (c, r, kappa, fmm_params);
  }

  // Export the supported factory combinations for callers outside this library.
  template class FMMInterface<Vec<1,Complex>,double,double>;
  template class FMMInterface<Vec<1,Complex>,Complex,double>;
  template class FMMInterface<Vec<1,Complex>,Complex,Complex>;
  template class FMMInterface<Vec<1,Complex32>,double,double>;
  template class FMMInterface<Vec<1,Complex32>,Complex,double>;
  template class FMMInterface<Vec<1,Complex32>,Complex,Complex>;
  template class FMMInterface<Vec<3,Complex>,double,double>;
  template class FMMInterface<Vec<3,Complex>,Complex,double>;
  template class FMMInterface<Vec<3,Complex>,Complex,Complex>;
  template class FMMInterface<Vec<3,Complex32>,double,double>;
  template class FMMInterface<Vec<3,Complex32>,Complex,double>;
  template class FMMInterface<Vec<3,Complex32>,Complex,Complex>;
  template class FMMInterface<Vec<6,Complex>,double,double>;
  template class FMMInterface<Vec<6,Complex32>,double,double>;

  template <typename KERNEL>
  class IntegralKernel : public BaseIntegralKernel<typename KERNEL::value_type>, public std::enable_shared_from_this<IntegralKernel<KERNEL>>
  {
    KERNEL kernel;
    using value_type = typename KERNEL::value_type;
    using components = decltype(kernel.Evaluate(Vec<3>(), Vec<3>(), Vec<3>(), Vec<3>()));
  public:
    IntegralKernel(KERNEL akernel) : kernel(std::move(akernel)) { }

    string Name() const override { return KERNEL::Name(); }
    IVec<2> Shape() const override { return KERNEL::Shape(); }
    AnalyticTriangleFormula GetAnalyticTriangleFormula() const override { return KERNEL::analytic_triangle_formula; }
    const BaseFMMInterface & Source() const override { return kernel.source; }
    const BaseFMMInterface & Target() const override { return kernel.target; }

    shared_ptr<const BaseIntegralKernel<value_type>> GetDifferentiatedKernel(const string & name) const override
    {
      if constexpr (!std::is_void_v<decltype(kernel.GetDifferentiatedKernel(name))>)
        {
          auto diffkernel = kernel.GetDifferentiatedKernel(name);
          if constexpr (IsStdVariant<decltype(diffkernel)>::value)
            return std::visit([](auto const & dk) -> shared_ptr<const BaseIntegralKernel<value_type>>
            { return MakeIntegralKernel(dk); }, diffkernel);
          else
            return MakeIntegralKernel(diffkernel);
        }
      else
        throw Exception("Kernel does not support differentiated kernel '"+name+"'");
    }

    size_t NumKernelComponents() const override { return components::SIZE; }

    FlatArray<const KernelTerm> Terms() const override { return { kernel.terms.Size(), kernel.terms.Data() }; }

    void EvaluatePairs(const SIMD_BaseMappedIntegrationRule & mirx,
                       const SIMD_BaseMappedIntegrationRule & miry,
                       FlatMatrix<SIMD<value_type>> values) const override
    {
      auto px = mirx.GetPoints();
      auto py = miry.GetPoints();
      auto nx = mirx.GetNormals();
      auto ny = miry.GetNormals();
      bool boundary = mirx.DimElement() == 2;
      for (size_t i = 0; i < mirx.Size(); i++)
        {
          Vec<3,SIMD<double>> normalx(0.0), normaly(0.0);
          if (boundary)
            {
              normalx = nx.Row(i);
              normaly = ny.Row(i);
            }
          auto weight = mirx[i].GetMeasure()*miry[i].GetMeasure()*mirx.IR()[i].Weight();
          Vec<components::SIZE,SIMD<value_type>> result = weight * kernel.Evaluate(Vec<3,SIMD<double>>(px.Row(i)), Vec<3,SIMD<double>>(py.Row(i)), normalx, normaly);
          for (size_t comp = 0; comp < result.Size(); comp++)
            values(comp,i) = result(comp);
        }
    }

    void EvaluateMatrix(const BaseMappedIntegrationRule & mirx,
                        const BaseMappedIntegrationRule & miry,
                        size_t component, double factor,
                        FlatMatrix<value_type> values, bool skip_diagonal) const override
    {
      auto px = mirx.GetPoints();
      auto py = miry.GetPoints();
      auto normalx = mirx.GetNormals();
      auto normaly = miry.GetNormals();
      bool boundary = mirx.DimElement() == 2;
      for (size_t i = 0; i < mirx.Size(); i++)
        for (size_t j = 0; j < miry.Size(); j++)
          {
            Vec<3> x = px.Row(i);
            Vec<3> y = py.Row(j);
            Vec<3> nx(0.0), ny(0.0);
            if (boundary)
              {
                nx = normalx.Row(i);
                ny = normaly.Row(j);
              }
            value_type value = 0.0;
            if (!skip_diagonal || L2Norm2(x-y) > 0)
              value = kernel.Evaluate(x,y,nx,ny)(component);
            double weight = mirx[i].GetWeight()*miry[j].GetWeight();
            values(i,j) = factor*weight*value;
          }
    }

    void AddPotential(const BaseMappedIntegrationPoint & mip,
                      const SIMD_BaseMappedIntegrationRule & miry,
                      FlatMatrix<SIMD<value_type>> values,
                      FlatVector<SIMD<value_type>> result, VorB source_vb) const override
    {
      for (int iy = 0; iy < miry.Size(); iy++)
        {
          Vec<3,SIMD<double>> x = mip.GetPoint();
          Vec<3,SIMD<double>> nx{0.0};
          if constexpr (KERNEL::target_type::needs_normal)
            nx = dynamic_cast<const MappedIntegrationPoint<2,3>&>(mip).GetNV();

          Vec<3,SIMD<double>> y = miry[iy].GetPoint();
          Vec<3,SIMD<double>> ny{0.0};
          if constexpr (KERNEL::source_type::needs_normal)
            {
              if (source_vb != BND)
                throw Exception("kernel requires boundary source normals");
              ny = static_cast<const SIMD<MappedIntegrationPoint<2,3>>&>(miry[iy]).GetNV();
            }

          auto eval = kernel.Evaluate(x, y, nx, ny);
          for (auto term : kernel.terms)
            {
              auto kernel_ = term.fac * eval(term.kernel_comp);
              result(term.test_comp) += miry[iy].GetWeight() * kernel_ * values(term.trial_comp, iy);
            }
        }
    }

    shared_ptr<ngla::BaseMatrix> CreateFMMOperator(Array<Vec<3>> xpts, Array<Vec<3>> ypts,
                                                 Array<Vec<3>> xnv, Array<Vec<3>> ynv,
                                                 const FMM_Parameters & params) const override
    {
      // Aliasing pointers keep the kernel alive without copying its source/target.
      auto owner = this->shared_from_this();
      shared_ptr<const BaseFMMInterface> source(owner, &kernel.source);
      shared_ptr<const BaseFMMInterface> target(owner, &kernel.target);
      return make_shared<FMM_Operator<value_type>>(KERNEL::Name(), KERNEL::Shape(), std::move(source), std::move(target),
                                                 std::move(xpts), std::move(ypts), std::move(xnv), std::move(ynv), params);
    }
  };

  template <typename KERNEL>
  shared_ptr<const BaseIntegralKernel<typename KERNEL::value_type>> MakeIntegralKernel(KERNEL kernel)
  {
    return make_shared<IntegralKernel<KERNEL>>(std::move(kernel));
  }

  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceSLKernel<3>::value_type>> MakeIntegralKernel(LaplaceSLKernel<3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceSLKernel<3,3>::value_type>> MakeIntegralKernel(LaplaceSLKernel<3,3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceSLKernel<3,1,Complex>::value_type>> MakeIntegralKernel(LaplaceSLKernel<3,1,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceSLKernel<3,3,Complex>::value_type>> MakeIntegralKernel(LaplaceSLKernel<3,3,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceSLKernel<3,1,double,float>::value_type>> MakeIntegralKernel(LaplaceSLKernel<3,1,double,float>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceSLKernel<3,3,double,float>::value_type>> MakeIntegralKernel(LaplaceSLKernel<3,3,double,float>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceSLKernel<3,1,Complex,Complex32>::value_type>> MakeIntegralKernel(LaplaceSLKernel<3,1,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceSLKernel<3,3,Complex,Complex32>::value_type>> MakeIntegralKernel(LaplaceSLKernel<3,3,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceDLKernel<3>::value_type>> MakeIntegralKernel(LaplaceDLKernel<3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceDLKernel<3,3>::value_type>> MakeIntegralKernel(LaplaceDLKernel<3,3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceDLKernel<3,1,Complex>::value_type>> MakeIntegralKernel(LaplaceDLKernel<3,1,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceDLKernel<3,3,Complex>::value_type>> MakeIntegralKernel(LaplaceDLKernel<3,3,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceDLKernel<3,1,double,float>::value_type>> MakeIntegralKernel(LaplaceDLKernel<3,1,double,float>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceDLKernel<3,3,double,float>::value_type>> MakeIntegralKernel(LaplaceDLKernel<3,3,double,float>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceDLKernel<3,1,Complex,Complex32>::value_type>> MakeIntegralKernel(LaplaceDLKernel<3,1,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LaplaceDLKernel<3,3,Complex,Complex32>::value_type>> MakeIntegralKernel(LaplaceDLKernel<3,3,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LameSLKernel<3>::value_type>> MakeIntegralKernel(LameSLKernel<3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<LameSLKernel<3,float>::value_type>> MakeIntegralKernel(LameSLKernel<3,float>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzSLKernel<3>::value_type>> MakeIntegralKernel(HelmholtzSLKernel<3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzSLKernel<3,3>::value_type>> MakeIntegralKernel(HelmholtzSLKernel<3,3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzSLKernel<3,1,Complex>::value_type>> MakeIntegralKernel(HelmholtzSLKernel<3,1,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzSLKernel<3,3,Complex>::value_type>> MakeIntegralKernel(HelmholtzSLKernel<3,3,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzSLKernel<3,1,double,Complex32>::value_type>> MakeIntegralKernel(HelmholtzSLKernel<3,1,double,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzSLKernel<3,3,double,Complex32>::value_type>> MakeIntegralKernel(HelmholtzSLKernel<3,3,double,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzSLKernel<3,1,Complex,Complex32>::value_type>> MakeIntegralKernel(HelmholtzSLKernel<3,1,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzSLKernel<3,3,Complex,Complex32>::value_type>> MakeIntegralKernel(HelmholtzSLKernel<3,3,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzDLKernel<3>::value_type>> MakeIntegralKernel(HelmholtzDLKernel<3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzDLKernel<3,3>::value_type>> MakeIntegralKernel(HelmholtzDLKernel<3,3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzDLKernel<3,1,Complex>::value_type>> MakeIntegralKernel(HelmholtzDLKernel<3,1,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzDLKernel<3,3,Complex>::value_type>> MakeIntegralKernel(HelmholtzDLKernel<3,3,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzDLKernel<3,1,double,Complex32>::value_type>> MakeIntegralKernel(HelmholtzDLKernel<3,1,double,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzDLKernel<3,3,double,Complex32>::value_type>> MakeIntegralKernel(HelmholtzDLKernel<3,3,double,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzDLKernel<3,1,Complex,Complex32>::value_type>> MakeIntegralKernel(HelmholtzDLKernel<3,1,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<HelmholtzDLKernel<3,3,Complex,Complex32>::value_type>> MakeIntegralKernel(HelmholtzDLKernel<3,3,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<CombinedFieldKernel<3>::value_type>> MakeIntegralKernel(CombinedFieldKernel<3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<CombinedFieldKernel<3,3>::value_type>> MakeIntegralKernel(CombinedFieldKernel<3,3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<CombinedFieldKernel<3,1,Complex>::value_type>> MakeIntegralKernel(CombinedFieldKernel<3,1,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<CombinedFieldKernel<3,3,Complex>::value_type>> MakeIntegralKernel(CombinedFieldKernel<3,3,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<CombinedFieldKernel<3,1,double,Complex32>::value_type>> MakeIntegralKernel(CombinedFieldKernel<3,1,double,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<CombinedFieldKernel<3,3,double,Complex32>::value_type>> MakeIntegralKernel(CombinedFieldKernel<3,3,double,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<CombinedFieldKernel<3,1,Complex,Complex32>::value_type>> MakeIntegralKernel(CombinedFieldKernel<3,1,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<CombinedFieldKernel<3,3,Complex,Complex32>::value_type>> MakeIntegralKernel(CombinedFieldKernel<3,3,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<MaxwellDLKernel<3>::value_type>> MakeIntegralKernel(MaxwellDLKernel<3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<MaxwellDLKernel<3,Complex>::value_type>> MakeIntegralKernel(MaxwellDLKernel<3,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<MaxwellDLKernel<3,double,Complex32>::value_type>> MakeIntegralKernel(MaxwellDLKernel<3,double,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<MaxwellDLKernel<3,Complex,Complex32>::value_type>> MakeIntegralKernel(MaxwellDLKernel<3,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffLaplaceSLKernel<3>::value_type>> MakeIntegralKernel(DiffLaplaceSLKernel<3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffLaplaceSLKernel<3,3>::value_type>> MakeIntegralKernel(DiffLaplaceSLKernel<3,3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffLaplaceSLKernel<3,1,Complex>::value_type>> MakeIntegralKernel(DiffLaplaceSLKernel<3,1,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffLaplaceSLKernel<3,3,Complex>::value_type>> MakeIntegralKernel(DiffLaplaceSLKernel<3,3,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffLaplaceSLKernel<3,1,double,float>::value_type>> MakeIntegralKernel(DiffLaplaceSLKernel<3,1,double,float>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffLaplaceSLKernel<3,3,double,float>::value_type>> MakeIntegralKernel(DiffLaplaceSLKernel<3,3,double,float>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffLaplaceSLKernel<3,1,Complex,Complex32>::value_type>> MakeIntegralKernel(DiffLaplaceSLKernel<3,1,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffLaplaceSLKernel<3,3,Complex,Complex32>::value_type>> MakeIntegralKernel(DiffLaplaceSLKernel<3,3,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzSLKernel<3>::value_type>> MakeIntegralKernel(DiffHelmholtzSLKernel<3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzSLKernel<3,3>::value_type>> MakeIntegralKernel(DiffHelmholtzSLKernel<3,3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzSLKernel<3,1,Complex>::value_type>> MakeIntegralKernel(DiffHelmholtzSLKernel<3,1,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzSLKernel<3,3,Complex>::value_type>> MakeIntegralKernel(DiffHelmholtzSLKernel<3,3,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzSLKernel<3,1,double,Complex32>::value_type>> MakeIntegralKernel(DiffHelmholtzSLKernel<3,1,double,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzSLKernel<3,3,double,Complex32>::value_type>> MakeIntegralKernel(DiffHelmholtzSLKernel<3,3,double,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzSLKernel<3,1,Complex,Complex32>::value_type>> MakeIntegralKernel(DiffHelmholtzSLKernel<3,1,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzSLKernel<3,3,Complex,Complex32>::value_type>> MakeIntegralKernel(DiffHelmholtzSLKernel<3,3,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzDLKernel<3>::value_type>> MakeIntegralKernel(DiffHelmholtzDLKernel<3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzDLKernel<3,3>::value_type>> MakeIntegralKernel(DiffHelmholtzDLKernel<3,3>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzDLKernel<3,1,Complex>::value_type>> MakeIntegralKernel(DiffHelmholtzDLKernel<3,1,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzDLKernel<3,3,Complex>::value_type>> MakeIntegralKernel(DiffHelmholtzDLKernel<3,3,Complex>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzDLKernel<3,1,double,Complex32>::value_type>> MakeIntegralKernel(DiffHelmholtzDLKernel<3,1,double,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzDLKernel<3,3,double,Complex32>::value_type>> MakeIntegralKernel(DiffHelmholtzDLKernel<3,3,double,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzDLKernel<3,1,Complex,Complex32>::value_type>> MakeIntegralKernel(DiffHelmholtzDLKernel<3,1,Complex,Complex32>);
  template NGS_DLL_HEADER shared_ptr<const BaseIntegralKernel<DiffHelmholtzDLKernel<3,3,Complex,Complex32>::value_type>> MakeIntegralKernel(DiffHelmholtzDLKernel<3,3,Complex,Complex32>);
}
