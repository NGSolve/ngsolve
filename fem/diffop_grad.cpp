#include <fem.hpp>
#include <diffop_impl.hpp>
  
namespace ngfem
{


  template <int D, typename FEL>
  shared_ptr<CoefficientFunction> DiffOpGradient<D,FEL> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir)
  {
    return -TransposeCF(dir->Operator("Grad")) * proxy;
    // return -dir->Operator("grad") * proxy;
  }

  template <int D, typename FEL>
  shared_ptr<CoefficientFunction> DiffOpGradientBoundary<D,FEL> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir)
  {
    // return -TransposeCF(dir->Operator("Gradboundary")) * proxy;

    int dim = dir->Dimension();
    auto n = NormalVectorCF(dim);
    n -> SetDimensions( Array<int> ( { dim, 1 } ) );
    auto Pn = n * TransposeCF(n);
    return  Pn * dir->Operator("Gradboundary") * proxy
      -TransposeCF(dir->Operator("Gradboundary")) * proxy;
  }

  template <int DIM_SPC>
  shared_ptr<CoefficientFunction> DiffOpGradVectorH1<DIM_SPC> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir)
  {
    return -proxy*dir->Operator("Grad");
  }

  
  template <int DIM_SPC>
  shared_ptr<CoefficientFunction> DiffOpGradBoundaryVectorH1<DIM_SPC> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir)
  {
    int dim = dir->Dimension();
    auto n = NormalVectorCF(dim);
    n -> SetDimensions( Array<int> ( { dim, 1 } ) );
    auto Pn = n * TransposeCF(n);
    //return  Pn * dir->Operator("Gradboundary") * proxy-TransposeCF(dir->Operator("Gradboundary")) * proxy;
    return  proxy * TransposeCF(dir->Operator("Gradboundary"))*Pn
      - proxy * dir->Operator("Gradboundary");
  }


  template class NGS_DLL_HEADER DiffOpGradient<1>;
  template class NGS_DLL_HEADER DiffOpGradient<2>;
  template class NGS_DLL_HEADER DiffOpGradient<3>;
  
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradient<1> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradient<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradient<3> >;

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradientBoundary<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradientBoundary<3> >;

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradVectorH1<1> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradVectorH1<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradVectorH1<3> >;

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradBoundaryVectorH1<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradBoundaryVectorH1<3> >;
}
