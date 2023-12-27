// #include <fem.hpp>

#include "bdbequations.hpp"
#include <diffop_impl.hpp>

namespace ngfem
{
  template <int D, typename FEL>
  shared_ptr<CoefficientFunction> DiffOpGradient<D,FEL> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian)
  {
    if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpGradient");
    return -TransposeCF(dir->Operator("Grad")) * proxy;
    // return -dir->Operator("grad") * proxy;
  }

  template <int D, typename FEL>
  shared_ptr<CoefficientFunction> DiffOpGradientBoundary<D,FEL> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian)
  {
    if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpGradientBoundary");    
    int dim = dir->Dimension();
    auto n = NormalVectorCF(dim) -> Reshape(Array<int> ( { dim, 1 } ));
    //n -> SetDimensions( Array<int> ( { dim, 1 } ) );
    auto Pn = n * TransposeCF(n);

    return (2*SymmetricCF(Pn * dir->Operator("Gradboundary"))
            -TransposeCF(dir->Operator("Gradboundary"))) * proxy;

  }

  template <int DIM_SPC>
  shared_ptr<CoefficientFunction> DiffOpGradVectorH1<DIM_SPC> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian)
  {
    if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpGradVectorH1");        
    return -proxy*dir->Operator("Grad");
  }

  
  template <int DIM_SPC>
  shared_ptr<CoefficientFunction> DiffOpGradBoundaryVectorH1<DIM_SPC> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian)
  {
    if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpGradBoundaryVectorH1");            
    int dim = dir->Dimension();
    auto n = NormalVectorCF(dim) -> Reshape(Array<int> ( { dim, 1 } ));
    //n -> SetDimensions( Array<int> ( { dim, 1 } ) );
    auto Pn = n * TransposeCF(n);

    return proxy * (2*SymmetricCF(Pn * dir->Operator("Gradboundary"))
                    - dir->Operator("Gradboundary"));
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
