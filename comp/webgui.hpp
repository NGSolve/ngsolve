#ifndef FILE_WEBGUI_HPP_INCLUDED
#define FILE_WEBGUI_HPP_INCLUDED

#include <core/archive.hpp>

#include "comp.hpp"

namespace webgui {
using namespace ngcomp;

struct WebguiData {
  int mesh_dim, order2d, order3d;
  double funcmin, funcmax, mesh_radius;
  vector<double> mesh_center;
  bool draw_vol, draw_surf, show_wireframe, show_mesh;

  vector<string> Bezier_trig_points, Bezier_points;
  vector<string> objects, names, edges;
};

unique_ptr<WebguiData> GenerateWebguiData(shared_ptr<MeshAccess> ma,
                                          shared_ptr<CoefficientFunction> cf,
                                          int order = 1);

inline string MeshToTextArchive(shared_ptr<MeshAccess> ma) {
  auto ss = make_shared<stringstream>();
  ma->GetNetgenMesh()->SetGeometry(nullptr);
  ma->GetNetgenMesh()->Save(*ss);
  // TextOutArchive ar(ss);
  // ar & (*ma);
  return ss->str();
}
inline string CFToTextArchive(shared_ptr<CoefficientFunction> cf) {
  auto ss = make_shared<stringstream>();
  TextOutArchive ar(ss);
  CoefficientFunction* p = cf.get();
  ar& p;
  return ss->str();
}

}  // namespace webgui

#endif  // FILE_WEBGUI_HPP_INCLUDED
