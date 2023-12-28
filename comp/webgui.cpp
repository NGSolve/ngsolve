#include "webgui.hpp"

#include <cmath>
// #include <comp.hpp>
#include "gridfunction.hpp"


namespace webgui {

static std::string base64_encode(FlatArray<unsigned char> in) {
  std::string out;

  int val = 0, valb = -6;
  for (unsigned char c : in) {
    val = (val << 8) + c;
    valb += 8;
    while (valb >= 0) {
      out.push_back(
          "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
              [(val >> valb) & 0x3F]);
      valb -= 6;
    }
  }
  if (valb > -6)
    out.push_back(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
            [((val << 8) >> (valb + 8)) & 0x3F]);
  while (out.size() % 4) out.push_back('=');
  return out;
}

double factorial(int n) { return std::tgamma(n + 1); }

double Binomial(int n, int i) {
  return factorial(n) / (factorial(i) * factorial(n - i));
}

double Bernstein(double x, int i, int n) {
  return Binomial(n, i) * pow(x, i) * pow(1 - x, n - i);
}

double BernsteinTrig(double x, double y, int i, int j, int n) {
  return factorial(n) / factorial(i) / factorial(j) / factorial(n - i - j) *
         pow(x, i) * pow(y, j) * pow(1 - x - y, n - i - j);
}

Matrix<double> GetIBernsteinBasis(ngfem::ELEMENT_TYPE etype, int order) {
  std::map<std::pair<ngfem::ELEMENT_TYPE, int>, Matrix<double>> cache;
  if(cache.count({etype, order})) return cache[{etype, order}];
  if (etype == ET_SEGM) {
    Matrix<double> ret(order + 1, order + 1);
    ret = 0.;
    for (auto i : Range(order + 1))
      for (auto j : Range(order + 1))
        ret(i, j) = Bernstein(1. * i / order, j, order);
    CalcInverse(ret);
    cache.insert({{etype, order}, ret});
    return ret;
  }

  if (etype == ET_TRIG) {
    int n = (order + 1) * (order + 2) / 2;
    Matrix<double> ret(n, n);
    int ii = 0;
    for (auto ix : Range(order + 1))
      for (auto iy : Range(order + 1 - ix)) {
        int jj = 0;
        for (auto jx : Range(order + 1))
          for (auto jy : Range(order + 1 - jx))
            ret(ii, jj++) =
                BernsteinTrig(1. * ix / order, 1. * iy / order, jx, jy, order);
        ii++;
      }
    CalcInverse(ret);
    cache.insert({{etype, order}, ret});
    return ret;
  }
  throw Exception("Element type not supported");
}

const IntegrationRule & GetElementPoints(ngfem::ELEMENT_TYPE etype, int order) {
  static std::map<std::pair<ngfem::ELEMENT_TYPE, int>, IntegrationRule> cache;
  if (cache.count({etype, order}))
    return cache[{etype, order}];

  using IP = IntegrationPoint;
  int n = order + 1;
  IntegrationRule ir;
  if (etype == ET_SEGM)
    for (auto i : Range(n))
      ir.Append(IP{1. * i / order, 0., 0.});
  if (etype == ET_TRIG)
    for (auto i : Range(n))
      for (auto j : Range(n - i))
        ir.Append(IP{1. * j / order, 1. * i / order, 0.});
  if (etype == ET_QUAD) {
    for (auto i : Range(n))
      for (auto j : Range(n - i))
        ir.Append(IP{1. * j / order, 1. * i / order, 0.});
    for (auto i : Range(n))
      for (auto j : Range(n - i))
        ir.Append(
            IP{1. - 1. * j / order, 1. - 1. * i / order, 0.});
  }
  if (etype == ET_TET) {
    ir.Append({ IP(1,0,0), IP(0,1,0), IP(0,0,1), IP(0,0,0) });
  }
  if(etype == ET_PYRAMID) {
    ir.Append( {
      IP(1,0,0), IP(0,1,0), IP(0,0,1), IP(0,0,0),
      IP(1,0,0), IP(0,1,0), IP(0,0,1), IP(1,1,0)
    });
  }
  if(etype == ET_PRISM) {
    ir.Append( {
      IP(1,0,0), IP(0,1,0), IP(0,0,1), IP(0,0,0),
      IP(0,0,1), IP(0,1,0), IP(0,1,1), IP(1,0,0),
      IP(1,0,1), IP(0,1,1), IP(1,0,0), IP(0,0,1)
    });
  }
  if(etype == ET_HEX) {
    ir.Append( {
      IP(1,0,0), IP(0,1,0), IP(0,0,1), IP(0,0,0),
      IP(0,1,1), IP(1,1,1), IP(1,1,0), IP(1,0,1),
      IP(1,0,1), IP(0,1,1), IP(1,0,0), IP(0,0,1),
      IP(0,1,1), IP(1,1,0), IP(0,1,0), IP(1,0,0),
      IP(0,0,1), IP(0,1,0), IP(0,1,1), IP(1,0,0),
      IP(1,0,1), IP(1,1,0), IP(0,1,1), IP(1,0,0)
    });
  }
  if(Dim(etype) == 3 && order == 2) {
    // append edge midpoints to each subtet
    auto mid = [](IP a, IP b) {
      return IP((a(0)+b(0))/2, (a(1)+b(1))/2, (a(2)+b(2))/2);
    };
    auto ntets = ir.Size() / 4;
    IntegrationRule ir_copy = std::move(ir);
    ir.SetSize0();
    for(auto i : Range(ntets))
    {
      auto tet = ir_copy.Range(i*4, i*4+4);
      ir.Append(tet);
      ir.Append(mid(tet[0], tet[3]));
      ir.Append(mid(tet[1], tet[3]));
      ir.Append(mid(tet[2], tet[3]));
      ir.Append(mid(tet[0], tet[1]));
      ir.Append(mid(tet[0], tet[2]));
      ir.Append(mid(tet[1], tet[2]));
    }
  }
  cache[{etype, order}] = std::move(ir);
  return cache[{etype, order}];
}

IntegrationRule GetWireframePoints(ngfem::ELEMENT_TYPE etype, int order) {
  int n = order + 1;
  IntegrationRule ir;
  if (etype == ET_TRIG) {
    for (auto i : Range(n)) ir.Append(IntegrationPoint{1. * i / order, 0., 0.});
    for (auto i : Range(n)) ir.Append(IntegrationPoint{0., 1. * i / order, 0.});
    for (auto i : Range(n))
      ir.Append(IntegrationPoint{1. * i / order, 1. - 1. * i / order, 0.});
  }
  if (etype == ET_QUAD) {
    for (auto i : Range(n)) ir.Append(IntegrationPoint{1. * i / order, 0., 0.});
    for (auto i : Range(n)) ir.Append(IntegrationPoint{0., 1. * i / order, 0.});
    for (auto i : Range(n)) ir.Append(IntegrationPoint{1. * i / order, 1., 0.});
    for (auto i : Range(n)) ir.Append(IntegrationPoint{1., 1. * i / order, 0.});
  }

  return ir;
}

vector<string> MapBernstein(FlatTensor<3, double> input, ngfem::ELEMENT_TYPE eltype,
                            int order, int n_components) {
  auto [ne, nip, comps] = input.Shape();
  Tensor<3, double> output(nip, ne, comps);
  Tensor<3, float> output_f(nip, ne, comps);
  auto trafo = GetIBernsteinBasis(eltype, order);

  for (size_t i = 0; i < ne; i++)
    output(STAR, i, STAR) = trafo * input(i, STAR, STAR);

  // convert double to float
  auto n = output.GetTotalSize();
  for (auto i : Range(n)) output_f.Data()[i] = output.Data()[i];

  vector<string> ret;
  for (auto i : Range(nip)) {
    auto vals = output_f(i, STAR, STAR);
    ret.push_back(base64_encode(FlatArray<unsigned char>(
        vals.Height() * vals.Width() * sizeof(vals(0, 0)),
        reinterpret_cast<unsigned char *>(vals.Data()))));
  }
  return ret;
}

vector<string> GenerateWireframeData( shared_ptr<MeshAccess> ma,
                                      shared_ptr<CoefficientFunction> cf,
                                      int order) {
  // generate wireframe data
  LocalHeapMem<100000> lh("webgui_wireframe_data");
  auto vb = ma->GetDimension() == 3 ? BND : VOL;
  Array<double> surface_data;
  auto comps = cf->Dimension();
  for (auto elnr : Range(ma->GetNElements(2))) {
    HeapReset hr(lh);

    auto el = ma->GetElement({vb, elnr});
    auto ir = GetWireframePoints(el.GetType(), order);
    auto &trafo = ma->GetTrafo(el, lh);
    BaseMappedIntegrationRule &mir = trafo(ir, lh);
    FlatMatrix<> values(ir.Size(), comps, lh);
    cf->Evaluate(mir, values);
    surface_data.Append(
        FlatArray<double>(values.AsVector().Size(), values.Data()));
  }
  auto nip = order+1;
  auto nseg = surface_data.Size() / nip / comps;
  FlatTensor<3, double> data(nseg, nip, comps, surface_data.Data());
  return MapBernstein(data, ET_SEGM, order, cf->Dimension());
}

vector<string> GenerateEdgeData( shared_ptr<MeshAccess> ma,
                                      shared_ptr<CoefficientFunction> cf,
                                      int order) {
  // generate wireframe data
  LocalHeapMem<100000> lh("webgui_edge_data");
  auto vb = ma->GetDimension() == 3 ? BBND : BND;
  Array<double> edge_data;
  auto comps = cf->Dimension();
  for (auto elnr : Range(ma->GetNElements(1))) {
    HeapReset hr(lh);

    auto el = ma->GetElement({vb, elnr});
    auto &ir = GetElementPoints(el.GetType(), order);
    auto &trafo = ma->GetTrafo(el, lh);
    BaseMappedIntegrationRule &mir = trafo(ir, lh);
    FlatMatrix<> values(ir.Size(), comps, lh);
    cf->Evaluate(mir, values);
    edge_data.Append(
        FlatArray<double>(values.AsVector().Size(), values.Data()));
  }
  auto nip = order+1;
  auto nseg = edge_data.Size() / nip / comps;
  FlatTensor<3, double> data(nseg, nip, comps, edge_data.Data());
  return MapBernstein(data, ET_SEGM, order, cf->Dimension());
}

vector<string> GenerateSurfaceData( shared_ptr<MeshAccess> ma,
                                    shared_ptr<CoefficientFunction> cf,
                                    int order) {
  // generate surface data
  LocalHeapMem<100000> lh("webgui_surface_data");
  auto vb = ma->GetDimension() == 3 ? BND : VOL;
  Array<double> surface_data;
  auto comps = cf->Dimension();
  for (auto elnr : Range(ma->GetNElements(2))) {
    HeapReset hr(lh);

    auto el = ma->GetElement({vb, elnr});
    auto &ir = GetElementPoints(el.GetType(), order);
    auto &trafo = ma->GetTrafo(el, lh);
    BaseMappedIntegrationRule &mir = trafo(ir, lh);
    FlatMatrix<> values(ir.Size(), comps, lh);
    cf->Evaluate(mir, values);
    surface_data.Append(
        FlatArray<double>(values.AsVector().Size(), values.Data()));
  }
  auto nip = GetElementPoints(ET_TRIG, order).Size();
  auto nel = surface_data.Size() / nip /
             comps;  // might be different from ma->GetNElements(2), because
                     // quads are divided into 2 trigs
  FlatTensor<3, double> data(nel, nip, comps, surface_data.Data());
  return MapBernstein(data, ET_TRIG, order, cf->Dimension());
}

vector<string> GenerateVolumeData( shared_ptr<MeshAccess> ma,
                                    shared_ptr<CoefficientFunction> cf,
                                    int order) {
  LocalHeapMem<100000> lh("webgui_volume_data");
  Array<double> volume_data;
  auto comps = cf->Dimension();
  for (auto elnr : Range(ma->GetNElements(3))) {
    HeapReset hr(lh);

    auto el = ma->GetElement({VOL, elnr});
    auto &ir = GetElementPoints(el.GetType(), order);
    auto &trafo = ma->GetTrafo(el, lh);
    BaseMappedIntegrationRule &mir = trafo(ir, lh);
    FlatMatrix<> values(ir.Size(), comps, lh);
    cf->Evaluate(mir, values);
    volume_data.Append(
        FlatArray<double>(values.AsVector().Size(), values.Data()));
  }
  auto nip = GetElementPoints(ET_TET, order).Size();
  auto nel = volume_data.Size() / nip /
             comps;  // might be different from ma->GetNElements(3), because
                     // hexes/prisms are divided into tets
  FlatTensor<3, double> data(nel, nip, comps, volume_data.Data());
  Tensor<3, double> output(nip, nel, comps);

  for (size_t i = 0; i < nel; i++)
    output(STAR, i, STAR) = data(i, STAR, STAR);

  Tensor<3, float> output_f(nip, nel, comps);
  auto n = output.GetTotalSize();
  for (auto i : Range(n)) output_f.Data()[i] = output.Data()[i];

  std::vector<string> ret;
  for (auto i : Range(nip)) {
    auto vals = output_f(i, STAR, STAR);
    ret.push_back(base64_encode(FlatArray<unsigned char>(
        vals.Height() * vals.Width() * sizeof(vals(0, 0)),
        reinterpret_cast<unsigned char *>(vals.Data()))));
  }
  return ret;
}

unique_ptr<WebguiData> GenerateWebguiData(shared_ptr<MeshAccess> ma,
                                          shared_ptr<CoefficientFunction> cf,
                                          int order) {
  auto d = make_unique<WebguiData>();
  d->mesh_dim = ma->GetDimension();
  d->order2d = min(order, 3);
  d->order3d = min(order, 2);
  d->funcmin = 0.0;
  d->funcmax = 1.0;
  d->draw_vol = true;
  d->draw_surf = true;
  d->show_wireframe = true;
  d->show_mesh = true;

  for (const auto &mat : ma->GetMaterials(VOL)) d->names.push_back(string(mat));

  netgen::Point3d pmin, pmax;
  ma->GetNetgenMesh()->GetBox(pmin, pmax);
  d->mesh_radius = (pmax - pmin).Length() / 2;
  auto c = Center(pmin, pmax);
  d->mesh_center = {c.X(), c.Y(), c.Z()};

  d->edges = GenerateEdgeData(ma, cf, order);
  d->Bezier_points = GenerateWireframeData(ma, cf, d->order2d);
  d->Bezier_trig_points = GenerateSurfaceData(ma, cf, d->order2d);
  d->points3d = GenerateVolumeData(ma, cf, d->order3d);

  return d;
}

}  // namespace webgui
