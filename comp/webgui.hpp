#ifndef FILE_WEBGUI_HPP_INCLUDED
#define FILE_WEBGUI_HPP_INCLUDED

#include <core/archive.hpp>

// #include "comp.hpp"
#include "meshaccess.hpp"


namespace webgui {
using namespace ngcomp;

struct WebguiArchiveData {
  shared_ptr<MeshAccess> mesh;
  shared_ptr<CoefficientFunction> cf;

  WebguiArchiveData() = default;

  WebguiArchiveData(shared_ptr<MeshAccess> mesh, shared_ptr<CoefficientFunction> cf)
    : mesh(mesh), cf(cf) { }

  shared_ptr<MeshAccess> GetMeshAccess() { return mesh; }

  void DoArchive(Archive & ar) {
    shared_ptr<netgen::NetgenGeometry> geo;
    if(ar.Output())
    {
      geo = mesh->GetNetgenMesh()->GetGeometry();
      mesh->GetNetgenMesh()->SetGeometry(nullptr);
    }

    ar & mesh & cf;

    if(ar.Output())
      mesh->GetNetgenMesh()->SetGeometry(geo);
  }
};

struct WebguiData {
  int mesh_dim, order2d, order3d;
  double funcmin, funcmax, mesh_radius;
  vector<double> mesh_center;
  bool draw_vol, draw_surf, show_wireframe, show_mesh;

  vector<string> Bezier_trig_points, Bezier_points, points3d;
  vector<string> objects, names, edges;
};

unique_ptr<WebguiData> GenerateWebguiData(shared_ptr<MeshAccess> ma,
                                          shared_ptr<CoefficientFunction> cf,
                                          int order = 1);

template<typename T>
string ToArchive(shared_ptr<T> obj, bool binary=true) {
  auto ss = make_shared<stringstream>();
  shared_ptr<netgen::NetgenGeometry> geo;
  if constexpr(is_same_v<T, netgen::Mesh>)
  {
    geo = obj->GetGeometry();
    obj->SetGeometry(nullptr);
  }
  if(binary) {
    BinaryOutArchive ar(ss);
    ar & obj;
  }
  else {
    TextOutArchive ar(ss);
    ar & obj;
  }
  if constexpr(is_same_v<T, netgen::Mesh>)
    obj->SetGeometry(geo);
  return ss->str();
}

template<typename T>
shared_ptr<T> FromArchive(string s, bool binary=true) {
  auto ss = make_shared<stringstream>(s);
  shared_ptr<T> obj;
  if(binary) {
    BinaryInArchive ar(ss);
    ar & obj;
  }
  else {
    TextInArchive ar(ss);
    ar & obj;
  }
  return obj;
}

}  // namespace webgui

#endif  // FILE_WEBGUI_HPP_INCLUDED
