#include <fem.hpp>
#include <diffop_impl.hpp>
  
namespace ngfem
{

  template <int D, typename FEL>
  shared_ptr<CoefficientFunction> DiffOpId<D,FEL> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian)
  {
    if (Eulerian) 
      // return dir->Deriv() * proxy;
      return proxy->Operator(make_shared<T_DifferentialOperator<DiffOpGradient<D>>>()) * dir;
    else
      return ZeroCF(Array<int>());
  }

  template <int DIM_SPC, VorB VB>
  shared_ptr<CoefficientFunction> DiffOpIdVectorH1<DIM_SPC,VB> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian)
  {
    if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpIdVectorH1");            
    return ZeroCF(Array<int>( {DIM_SPC} ));
  }

  //Why do I need to specify this explicitly? See end of this file...
  //If I neglect it I get /usr/bin/ld: ../comp/libngcomp.so: undefined reference to `ngfem::DiffOpIdVectorH1<3, (ngfem::VorB)2>::DiffShape(std::shared_ptr<ngfem::CoefficientFunction>, std::shared_ptr<ngfem::CoefficientFunction>)'
  template <>
  shared_ptr<CoefficientFunction> DiffOpIdVectorH1<3,BBND> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian)
  {
    if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpIdVectorH1");                
    return ZeroCF(Array<int>( {3} ));
  }


  template class NGS_DLL_HEADER DiffOpId<1>;
  template class NGS_DLL_HEADER DiffOpId<2>;
  template class NGS_DLL_HEADER DiffOpId<3>;
  
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpId<1>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpId<2>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpId<3>>;

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdBoundary<1>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdBoundary<2>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdBoundary<3>>;

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<1,VOL>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<2,VOL>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<3,VOL>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<1,BND>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<2,BND>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<3,BND>>;
  //Here the explicit case should be already handled!!
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<3,BBND>>;
}
