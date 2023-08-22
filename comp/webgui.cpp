#include "webgui.hpp"

#include <cmath>
#include <comp.hpp>

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
  if (etype == ET_SEGM) {
    Matrix<double> ret(order + 1, order + 1);
    ret = 0.;
    for (auto i : Range(order + 1))
      for (auto j : Range(order + 1))
        ret(i, j) = Bernstein(1. * i / order, j, order);
    CalcInverse(ret);
    return {ret};
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
    return ret;
  }
  throw Exception("Element type not supported");
}

IntegrationRule GetElementPoints(ngfem::ELEMENT_TYPE etype, int order) {
  int n = order + 1;
  IntegrationRule ir;
  if (etype == ET_TRIG)
    for (auto i : Range(n))
      for (auto j : Range(n - i))
        ir.Append(IntegrationPoint{1. * j / order, 1. * i / order, 0.});
  if (etype == ET_QUAD) {
    for (auto i : Range(n))
      for (auto j : Range(n - i))
        ir.Append(IntegrationPoint{1. * j / order, 1. * i / order, 0.});
    for (auto i : Range(n))
      for (auto j : Range(n - i))
        ir.Append(
            IntegrationPoint{1. - 1. * j / order, 1. - 1. * i / order, 0.});
  }

  return ir;
}

IntegrationRule GetWireframePoints(ngfem::ELEMENT_TYPE etype, int order) {
  int n = order + 1;
  IntegrationRule ir;
  if (etype == ET_TRIG) {
    for (auto i : Range(n)) ir.Append({1. * i / order, 0., 0.});
    for (auto i : Range(n)) ir.Append({0., 1. * i / order, 0.});
    for (auto i : Range(n))
      ir.Append({1. * i / order, 1. - 1. * i / order, 0.});
  }
  if (etype == ET_QUAD) {
    for (auto i : Range(n)) ir.Append({1. * i / order, 0., 0.});
    for (auto i : Range(n)) ir.Append({0., 1. * i / order, 0.});
    for (auto i : Range(n)) ir.Append({1. * i / order, 1., 0.});
    for (auto i : Range(n)) ir.Append({1., 1. * i / order, 0.});
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

unique_ptr<WebguiData> GenerateWebguiData(shared_ptr<MeshAccess> ma,
                                          shared_ptr<CoefficientFunction> cf,
                                          int order) {
  auto d = make_unique<WebguiData>();
  d->mesh_dim = ma->GetDimension();
  d->order2d = order;
  d->order3d = order;
  d->funcmin = 0.0;
  d->funcmax = 1.0;
  d->draw_vol = false;
  d->draw_surf = true;
  d->show_wireframe = true;
  d->show_mesh = true;

  for (const auto &mat : ma->GetMaterials(VOL)) d->names.push_back(string(mat));

  netgen::Point3d pmin, pmax;
  ma->GetNetgenMesh()->GetBox(pmin, pmax);
  d->mesh_radius = (pmax - pmin).Length() / 2;
  auto c = Center(pmin, pmax);
  d->mesh_center = {c.X(), c.Y(), c.Z()};

  // generate wireframe data
  LocalHeapMem<100000> lh("webgui");
  auto vb = d->mesh_dim == 3 ? BND : VOL;
  Array<double> wireframe_data;
  auto comps = cf->Dimension();
  for (auto elnr : Range(ma->GetNElements(2))) {
    HeapReset hr(lh);

    auto el = ma->GetElement({vb, elnr});
    auto ir = GetElementPoints(el.GetType(), order);
    auto &trafo = ma->GetTrafo(el, lh);
    BaseMappedIntegrationRule &mir = trafo(ir, lh);
    FlatMatrix<> values(ir.Size(), comps, lh);
    cf->Evaluate(mir, values);
    wireframe_data.Append(
        FlatArray<double>(values.AsVector().Size(), values.Data()));
  }
  auto nip = GetElementPoints(ET_TRIG, order).Size();
  auto nel = wireframe_data.Size() / nip /
             comps;  // might be different from ma->GetNElements(2), because
                     // quads are divided into 2 trigs
  FlatTensor<3, double> data(wireframe_data.Size() / nip / comps, nip, comps,
                             wireframe_data.Data());
  d->Bezier_trig_points = MapBernstein(data, ET_TRIG, order, cf->Dimension());

  return d;
}

}  // namespace webgui
