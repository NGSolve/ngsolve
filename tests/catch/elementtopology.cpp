
#include "catch.hpp"
#include <fem.hpp>

using namespace ngfem;

constexpr double tolerance = 1e-12;

// distorted, non-affine element, flip mirrors it (negative det)
template <int D>
static Matrix<> ElementPoints (ELEMENT_TYPE et, bool flip)
{
  auto ref = ElementTopology::GetVertices<D>(et);
  Matrix<> p(ref.Size(), D);
  for (int i = 0; i < ref.Size(); i++)
    for (int j = 0; j < D; j++)
      p(i,j) = ref[i](j) * (1+0.1*j) + 0.07*sin(3*i+j) + 0.05*ref[i](0)*ref[i](D-1);
  if (flip) p.Col(0) *= -1;
  return p;
}

// old formula: Cof(J) * n_ref, normalized with sign(det)
template <int D, typename T>
static void RefNormal (const Mat<D,D,T> & J, T det, ELEMENT_TYPE et, int f, Vec<D,T> & nv, T & meas)
{
  Vec<D> nref = ElementTopology::GetNormals<D>(et)[f];
  Vec<D,T> n = Cof(J) * nref;
  meas = L2Norm(n);
  nv = (IfPos(det, T(1), T(-1))/meas) * n;
}

template <int D>
static Vec<D> RefEdgeTangent (ELEMENT_TYPE et, int e)
{
  auto pts = ElementTopology::GetVertices<D>(et);
  auto edge = ElementTopology::GetEdges(et)[e];
  return pts[edge[1]]-pts[edge[0]];
}

template <int D>
static void TestNormals (ELEMENT_TYPE et, bool flip)
{
  LocalHeap lh(1000000, "lh");
  FE_ElementTransformation<D,D> trafo(et, ElementPoints<D>(et, flip));

  for (VorB vb : { BND, BBND })
    {
      if (D == 2 && vb == BBND) continue;
      int nf = (vb == BND) ? ElementTopology::GetNFacets(et) : ElementTopology::GetNEdges(et);
      for (int f = 0; f < nf; f++)
        {
          HeapReset hr(lh);
          Facet2ElementTrafo f2el(et, vb);
          ELEMENT_TYPE ft = (vb == BND) ? ElementTopology::GetFacetType(et, f) : ET_SEGM;

          SIMD_IntegrationRule sirf(ft, 4);
          auto & smir = static_cast<SIMD_MappedIntegrationRule<D,D>&> (trafo(f2el(f, sirf, lh), lh));
          for (size_t i = 0; i < smir.Size(); i++)
            for (int l = 0; l < SIMD<double>::Size(); l++)
              {
                auto & mip = smir[i];
                if (vb == BND)
                  {
                    Vec<D,SIMD<double>> n; SIMD<double> m;
                    RefNormal<D> (mip.GetJacobian(), mip.GetJacobiDet(), et, f, n, m);
                    CHECK(fabs(mip.GetMeasure()[l]-m[l]) < tolerance*m[l]);
                    for (int j = 0; j < D; j++)
                      CHECK(fabs(mip.GetNV()(j)[l]-n(j)[l]) < tolerance);
                    if constexpr (D == 2)
                      {
                        CHECK(fabs(mip.GetTV()(0)[l]+n(1)[l]) < tolerance);
                        CHECK(fabs(mip.GetTV()(1)[l]-n(0)[l]) < tolerance);
                      }
                  }
                else
                  {
                    Vec<D,SIMD<double>> tau = mip.GetJacobian() * RefEdgeTangent<D>(et, f);
                    SIMD<double> m = L2Norm(tau);
                    CHECK(fabs(mip.GetMeasure()[l]-m[l]) < tolerance*m[l]);
                  }
              }

          IntegrationRule irf(ft, 4);
          auto & mir = static_cast<MappedIntegrationRule<D,D>&> (trafo(f2el(f, irf, lh), lh));
          for (auto & mip : mir)
            {
              if (vb == BND)
                {
                  Vec<D> n; double m;
                  RefNormal<D,double> (mip.GetJacobian(), mip.GetJacobiDet(), et, f, n, m);
                  CHECK(fabs(mip.GetMeasure()-m) < tolerance*m);
                  CHECK(L2Norm(mip.GetNV()-n) < tolerance);
                  if constexpr (D == 2)
                    CHECK(L2Norm(mip.GetTV()-Vec<2>(-n(1),n(0))) < tolerance);
                }
              else
                {
                  Vec<D> tau = mip.GetJacobian() * RefEdgeTangent<D>(et, f);
                  double m = L2Norm(tau);
                  CHECK(fabs(mip.GetMeasure()-m) < tolerance*m);
                  CHECK(L2Norm(mip.GetTV()-(1/m)*tau) < tolerance);
                }
            }

          if (vb == BND)
            {
              Mat<D,D> J = mir[0].GetJacobian();
              Vec<D> n = ElementTopology::GetPhysicalNormal(et, f, J);
              Vec<D> nref = Cof(J) * ElementTopology::GetNormals<D>(et)[f];
              CHECK(L2Norm(n-nref) < tolerance*L2Norm(nref));
            }
        }
    }
}

template <ELEMENT_TYPE ET>
static void TestTables ()
{
  using T = ET_trait<ET>;
  CHECK(T::vertices.size() == T::N_VERTEX);
  CHECK(T::edges.size() == T::N_EDGE);
  CHECK(T::faces.size() == T::N_FACE);
  CHECK(ElementTopology::GetEdges(ET).Size() == T::N_EDGE);
  CHECK(ElementTopology::GetFaces(ET).Size() == T::N_FACE);
  for (int i = 0; i < T::N_EDGE; i++)
    CHECK(T::GetEdge(i) == ElementTopology::GetEdges(ET)[i]);
  if constexpr (T::DIM >= 2)
    Iterate<T::N_FACET> ([&] (auto f)
      {
        constexpr int F = decltype(f)::value;
        // facet tangents reproduce the reference normal exactly, except at the pyramid apex
        if constexpr (ET == ET_PYRAMID && F < 4)
          CHECK(RefFacetSign<ET,F>() == 0);
        else
          CHECK(abs(RefFacetSign<ET,F>()) == 1);
      });
}

TEST_CASE ("ElementTopology", "[fem][elementtopology]")
{
  SECTION ("tables")
    {
      TestTables<ET_SEGM>(); TestTables<ET_TRIG>(); TestTables<ET_QUAD>();
      TestTables<ET_TET>(); TestTables<ET_PRISM>(); TestTables<ET_PYRAMID>(); TestTables<ET_HEX>();
      CHECK_THROWS(ElementTopology::GetVertices<3>(ET_TRIG));
      CHECK_THROWS(ElementTopology::GetNormals<3>(ET_HEXAMID));
    }

  for (bool flip : { false, true })
    {
      SECTION (std::string("normals and measures, det ") + (flip ? "< 0" : "> 0"))
        {
          for (auto et : { ET_TRIG, ET_QUAD })
            TestNormals<2> (et, flip);
          for (auto et : { ET_TET, ET_PRISM, ET_PYRAMID, ET_HEX })
            TestNormals<3> (et, flip);
        }
    }
}
