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
    auto hesse_dir = dir->Operator("hesseboundary");
    auto grad_u = proxy->Primary()->Operator("Gradboundary");
    auto grad_dir = dir->Operator("Gradboundary");
    auto grad_dir_t = TransposeCF(grad_dir);
    auto n = NormalVectorCF(D)->Reshape(Array<int>({D, 1}));
    auto grad_n = n->Operator("Grad");
    auto n_t = TransposeCF(n);
    auto grad_n_t = TransposeCF(grad_n);
    auto Pn = n * n_t;
    auto Ps = IdentityCF(D) - Pn;
    Array<int> md = {-1,D};
    auto d_pn = [&](shared_ptr<CoefficientFunction> v_) {
      auto v = v_->Reshape(md);
      // return ((n_t * v) * grad_n + (grad_n_t * v)->Reshape(Array<int>{D,1}) * n_t)->Reshape(md);
      return ((v * n) * grad_n->Reshape(Array<int>{1,-1}))->Reshape(md) + ((v*grad_n)->Reshape(Array<int>{-1,1})*n_t);
    };
    auto shape = proxy->Dimensions();

    // in comments are scalar valued equivalents (for easier reading)

    // auto part1a = - (TransposeCF(hesse_dir->Reshape(Array<int>{D,D*D})) * (grad_u))->Reshape(Array<int>{D,D}) * Ps;
    auto part1a = - ((grad_u->Reshape(Array<int>{-1,D}) * hesse_dir->Reshape(Array<int>{D,-1}))->Reshape(Array<int>{-1,D,D})->TensorTranspose(1,2)->Reshape(md) * Ps)->Reshape(shape);
    // auto part1b = grad_dir_t * d_pn(grad_u)->Reshape(Array<int>{D,D}) * Ps;
    auto part1b = (((d_pn(grad_u) * Ps)->Reshape(Array<int>{-1,D,D})->TensorTranspose(1,2)->Reshape(md) * grad_dir)->Reshape(Array<int>{-1,D,D})->TensorTranspose(1,2))->Reshape(shape);
    // auto part1c = - grad_dir_t * proxy;
    auto part1c = - (proxy->Reshape(md) * grad_dir)->Reshape(Array<int>{-1,D,D})->TensorTranspose(1,2)->Reshape(shape);
    // auto part2a = d_pn(grad_dir * grad_u) * Ps;
    auto part2a = (d_pn(grad_u->Reshape(md) * grad_dir_t) * Ps)->Reshape(shape);
    // auto part2b = Pn * (hesse_dir->Reshape(Array<int>{D,D,D}) * grad_u) * Ps;
    auto part2b = (TransposeCF((Pn * (hesse_dir->Reshape(Array<int>{-1,D}) * TransposeCF(grad_u->Reshape(Array<int>{-1,D})))->Reshape(Array<int>{D,-1}))->Reshape(Array<int>{D*D,-1}))->Reshape(md) * Ps)->Reshape(shape);
    // auto part2c = Pn * grad_dir * proxy;
    auto part2c = (proxy->Reshape(md) * grad_dir_t * Pn)->Reshape(Array<int>{-1,D,D})->TensorTranspose(1,2)->Reshape(shape);
    // auto part3a = d_pn(grad_u) * Ps * grad_dir;
    auto part3a = ((d_pn(grad_u) * Ps)->Reshape(md) * grad_dir)->Reshape(shape);
    auto part3b = - (proxy->Reshape(md) * grad_dir)->Reshape(shape);
    // auto part4a = - d_pn(grad_u) * grad_dir_t * Pn;
    auto part4a = - ((d_pn(grad_u) * grad_dir_t)->Reshape(md) * Pn)->Reshape(shape);
    auto part4b = ((proxy->Reshape(md) * grad_dir_t)->Reshape(md) * Pn)->Reshape(shape);

    return part1a + part1b + part1c + part2a + part2b + part2c + part3a + part3b + part4a + part4b;
  }

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<1> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<3> >;

  // template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<1> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<3> >;
}
