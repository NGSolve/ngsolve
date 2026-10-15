#include "bdbequations.hpp"
#include "diffop_impl.hpp"
  
namespace ngfem
{ 
  template <int D>
  shared_ptr<CoefficientFunction>
  DiffOpHesse<D>::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool eulerian)
  {
    if(eulerian)
      throw Exception("DiffShape Eulerian not implemented for Hesse!");

    auto grad_proxy = proxy->Primary()->Operator("Grad")->Reshape(Array<int>{-1,D});
    auto hesse_grad_dir = proxy->Reshape(Array<int>({-1,D})) * dir->Operator("Grad");
    const auto& shape = proxy->Dimensions();
    return - (grad_proxy * dir->Operator("hesse")->Reshape(Array<int>{D,-1}))->Reshape(shape)
      - hesse_grad_dir->Reshape(Array<int>{-1,D,D})->TensorTranspose(1,2)->Reshape(shape)
  - hesse_grad_dir->Reshape(shape);
  };

  template <int D, typename FEL>
  shared_ptr<CoefficientFunction>
  DiffOpHesseBoundary<D, FEL>::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool eulerian)
  {
    if(eulerian)
      throw Exception("DiffShape Eulerian not implemented for HesseBoundary!");
    // hb = covariant Hessian (tangential-tangential), G = grad_s V, K = G^T Pn - G
    // d hb = K^T hb + hb K + (n.G grad_s u) W - sum_i (grad_s u)_i hb(V_i)
    Array<int> mdd = {-1,D,D};
    Array<int> md = {-1,D};
    auto shape = proxy->Dimensions();
    auto grad_u = proxy->Primary()->Operator("Gradboundary")->Reshape(Array<int>{-1,D});
    auto grad_dir = dir->Operator("Gradboundary");
    auto hesse_dir = dir->Operator("hesseboundary")->Reshape(Array<int>{D,-1});
    auto n = NormalVectorCF(D)->Reshape(Array<int>({D, 1}));
    auto W = NormalVectorCF(D)->Operator("Grad")->Reshape(Array<int>{1,-1});
    auto Pn = n * TransposeCF(n);
    auto K = TransposeCF(grad_dir) * Pn - grad_dir;

    auto hbK = proxy->Reshape(md) * K;
    auto KthbT = proxy->Reshape(mdd)->TensorTranspose(1,2)->Reshape(md) * K;
    auto curv = ((grad_u * (TransposeCF(grad_dir) * n)) * W - grad_u * hesse_dir)->Reshape(shape);
    return hbK->Reshape(shape) + KthbT->Reshape(mdd)->TensorTranspose(1,2)->Reshape(shape) + curv;
  }

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<1> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<3> >;

  // template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<1> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<3> >;
}
