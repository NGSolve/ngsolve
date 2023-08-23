#include <emscripten/bind.h>
#include <vector>

#include <comp.hpp>
#include <meshing.hpp>
#include "coefficient.hpp"
#include "hdivdivfespace.hpp"
#include "hcurldivfespace.hpp"
#include "hcurlcurlfespace.hpp"
#include "normalfacetfespace.hpp"
#include "../fem/hdivdivfe.hpp"
#include "hdivdivsurfacespace.hpp"
#include "numberfespace.hpp"
#include "irspace.hpp"
#include "compressedfespace.hpp"
#include "../fem/integratorcf.hpp"
#include "../fem/h1lofe.hpp"
#include "contact.hpp"
#include "globalinterfacespace.hpp"
#include "globalspace.hpp"
#include "vtkoutput.hpp"
#include "webgui.hpp"

namespace emscripten {
namespace internal {

// Auto-convert std::vector<T> to JS Array
template <typename T, typename Allocator>
struct BindingType<std::vector<T, Allocator>> {
  using ValBinding = BindingType<val>;
  using WireType = ValBinding::WireType;

  static WireType toWireType(const std::vector<T, Allocator> &vec) {
    std::vector<val> valVec(vec.begin(), vec.end());
    return BindingType<val>::toWireType(val::array(valVec));
  }

  static std::vector<T, Allocator> fromWireType(WireType value) {
    return vecFromJSArray<T>(ValBinding::fromWireType(value));
  }
};

template <typename T>
struct TypeID<
    T,
    typename std::enable_if_t<std::is_same<
        typename Canonicalized<T>::type,
        std::vector<
            typename Canonicalized<T>::type::value_type,
            typename Canonicalized<T>::type::allocator_type>>::value>> {
  static constexpr TYPEID get() {
    return TypeID<val>::get();
  }
};

} // namespace internal
} // namespace emscripten

using namespace emscripten;
using namespace ngcomp;
using namespace webgui;

shared_ptr<MeshAccess> LoadMesh(string s, bool binary) {
  return make_shared<MeshAccess>(FromArchive<netgen::Mesh>(s, binary));
}

shared_ptr<CoefficientFunction> LoadCoefficientFunction(string s, bool binary) {
  return FromArchive<CoefficientFunction>(s,binary);
}

auto MyEvaluate(shared_ptr<MeshAccess> ma, shared_ptr<CoefficientFunction> cf) {
  std::vector<double> res;
  LocalHeapMem<10000> lh("CF(MeshPoint)");
  for(auto el : ma->Elements(VOL))
  {
    HeapReset hrt(lh);
    auto & trafo = ma->GetTrafo(ElementId(VOL, el.Nr()), lh);
    auto ir = SelectIntegrationRule(el.GetType(), 3);
    auto & mir = trafo(ir, lh);
    FlatMatrix<double> values(ir.Size(), cf->Dimension(), lh);
    cf->Evaluate(mir, BareSliceMatrix<double>(values));
    for(auto v : values.AsVector())
      res.push_back(v);
  }
  return res;
}

EMSCRIPTEN_BINDINGS(em_ngs) {
  class_<MeshAccess>("Mesh")
  .smart_ptr<std::shared_ptr<MeshAccess>>("Mesh")
  .function("GetDimension", &MeshAccess::GetDimension)
  .function("GetNP", &MeshAccess::GetNP)
  .function("GetNV", &MeshAccess::GetNV)
  .function("GetNE", select_overload<size_t ()const>(&MeshAccess::GetNE))
  .function("GetNSE", &MeshAccess::GetNSE)
  .function("Curve", &MeshAccess::Curve)
  ;

  class_<CoefficientFunction>("CoefficientFunction")
  .smart_ptr<std::shared_ptr<CoefficientFunction>>("CoefficientFunction")
    ;

  emscripten::function("LoadCoefficientFunction", &FromArchive<CoefficientFunction>);
  emscripten::function("LoadMesh", &LoadMesh);
  emscripten::function("Evaluate", &MyEvaluate);
  emscripten::function("GenerateWebguiData", &GenerateWebguiData);

  class_<WebguiData>("WebguiData")
    .property("mesh_dim", &WebguiData::mesh_dim)
    .property("order2d", &WebguiData::order2d)
    .property("order3d", &WebguiData::order3d)
    .property("funcmin", &WebguiData::funcmin)
    .property("funcmax", &WebguiData::funcmax)
    .property("mesh_radius", &WebguiData::mesh_radius)
    .property("mesh_center", &WebguiData::mesh_center)
    .property("draw_vol", &WebguiData::draw_vol)
    .property("draw_surf", &WebguiData::draw_surf)
    .property("show_wireframe", &WebguiData::show_wireframe)
    .property("show_mesh", &WebguiData::show_mesh)
    .property("Bezier_trig_points", &WebguiData::Bezier_trig_points)
    .property("Bezier_points", &WebguiData::Bezier_points)
    .property("objects", &WebguiData::objects)
    .property("names", &WebguiData::names)
    .property("edges", &WebguiData::edges)
    ;

}

int main() {
  cout << "hello from main!" << endl;

  // the following code is just to check the binary output file size

  // auto m = LoadMesh("");
  // Flags flags;
  // new H1HighOrderFESpace(ma, flags));
  // new L2HighOrderFESpace(ma, flags));

  // new FacetFESpace (m, flags);
  // new FacetSurfaceFESpace (m, flags);
  // new H1HighOrderFESpace (m, flags);
  // new H1LumpingFESpace (m, flags);
  // new HCurlCurlFESpace (m, flags);
  // new HCurlDivFESpace (m, flags);
  // new HDivDivFESpace (m, flags);
  // new HDivDivSurfaceSpace (m, flags);
  // new HDivHighOrderSurfaceFESpace (m, flags);
  // new IntegrationRuleSpace (m, flags);
  // new IntegrationRuleSpaceSurface (m, flags);
  // new L2HighOrderFESpace (m, flags);
  // new L2SurfaceHighOrderFESpace (m, flags);
  // new NodalFESpace (m, flags);
  // new NormalFacetFESpace (m, flags);
  // new NormalFacetSurfaceFESpace (m, flags);
  // new NumberFESpace (m, flags);
  // new TangentialSurfaceL2FESpace (m, flags);
  // new VectorFESpace<FacetFESpace> (m, flags);
  // new VectorFESpace<FacetSurfaceFESpace> (m, flags);
  // new VectorFESpace<L2SurfaceHighOrderFESpace> (m, flags);
  // new VectorFESpace<NodalFESpace> (m, flags);
  // new VectorFacetFESpace (m, flags);
  // new VectorH1FESpace (m, flags);
  // new VectorL2FESpace (m, flags);
}

