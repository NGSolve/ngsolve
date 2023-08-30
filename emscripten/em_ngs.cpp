#include <emscripten/bind.h>
#include <vector>

#include <comp.hpp>
#include <meshing.hpp>
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

shared_ptr<WebguiArchiveData> LoadData(string s, bool binary) {
  return FromArchive<WebguiArchiveData>(s, binary);
}

shared_ptr<MeshAccess> LoadMesh(string s, bool binary) {
  return FromArchive<MeshAccess>(s, binary);
}

shared_ptr<CoefficientFunction> LoadCoefficientFunction(string s, bool binary) {
  return FromArchive<CoefficientFunction>(s,binary);
}

EMSCRIPTEN_BINDINGS(em_ngs) {
  class_<WebguiArchiveData>("WebguiArchiveData")
  .smart_ptr<std::shared_ptr<WebguiArchiveData>>("WebguiArchiveData")
  .property("mesh", &WebguiArchiveData::mesh)
  .property("cf", &WebguiArchiveData::cf)
    ;

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
  emscripten::function("LoadData", &LoadData);
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
    .property("points3d", &WebguiData::points3d)
    .property("objects", &WebguiData::objects)
    .property("names", &WebguiData::names)
    .property("edges", &WebguiData::edges)
    ;

}

int main() {
  cout << "hello from main!" << endl;
}

