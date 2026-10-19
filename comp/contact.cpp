#include <core/register_archive.hpp>
#include <netgen_config.hpp>

#include "contact.hpp"

// disable contact pair visualization due to linking issues with netgen gui lib
#undef NETGEN_USE_GUI
#define NETGEN_USE_GUI 0

#if NETGEN_USE_GUI
#include <visual.hpp>
#include <visualization/mvdraw.hpp>
#endif // NETGEN_USE_GUI


namespace ngcomp
{
  static IntegrationRule CopyIR(const IntegrationRule& ir)
  {
    IntegrationRule res;
    for (const auto& ip : ir)
      res.AddIntegrationPoint(ip);
    return res;
  }

  inline int GetDomIn(const MeshAccess& ma, const Ngs_Element& el)
  {
    if(ma.GetDimension() ==3)
      return ma.GetNetgenMesh()->GetFaceDescriptor(netgen::FaceRegionIndex::FromNr0(el.GetIndex())).DomainIn();
    else
    {
      auto& seg = (*ma.GetNetgenMesh())[netgen::SegmentIndex::FromNr0(el.Nr())];
      return ma.GetNetgenMesh()->GetEdgeDescriptor(seg.GetIndex()).DomainIn();
    }
  }

  // Quadratic approximation of distance to pmaster
  template<int DIM>
  struct T2
  {
    Mat<DIM> a;
    Vec<DIM> b;
    Vec<DIM> x0;
    double d0;

    template<int DIMR>
    T2( MappedIntegrationPoint<DIM, DIMR> & mip, Vec<DIMR> pmaster ) : x0(mip.IP().Point())
    {
      // IntegrationPoint ip {x0};
      auto jac = mip.GetJacobian();
      auto hesse = mip.CalcHesse();
      Vec<DIMR> dist = mip.GetPoint() - pmaster;

      if constexpr (DIM==1)
        b = 2*InnerProduct(jac,dist);
      else
        b = 2*Trans(jac)*dist;

      a = 2*Trans(jac)*jac;
      for (int d : Range(DIMR))
        a += dist[d] * hesse[d];
      d0 = L2Norm2(dist);
    }

    bool CalcMinimumOnSegment( Vec<DIM> x1, Vec<DIM> x2, Vec<DIM> & res )
    {
      auto v = x2-x1;
      Vec<DIM> va = a*v;
      double lam = - InnerProduct(v,b) - InnerProduct(va,x1-x0);
      lam /= InnerProduct(va,v);

      if(lam<0 || lam>1.0)
          return false;

      res = lam*x2 + (1-lam)*x1;
      return true;
    }

    Vec<DIM> CalcMinimum( )
    {
      return Inv(a)*(a*x0-b);
    }

    double operator() (Vec<DIM> x)
    {
      x -= x0;
      return d0 + 0.5*InnerProduct(x, a*x) + InnerProduct(b,x);
    }
  };

  template<int DIMS, int DIMR>
  double FindClosestPoint( Vec<DIMR> pmaster, Vec<DIMR> n, double h, const DeformedBndElement<DIMR> & el2, IntegrationPoint & ip, Vec<DIMR> & p, bool both_sides)
  // input arguments: pmaster, n (oriented primary normal), h (maximum distance), el2
  // output arguments: ip, p
  {
    auto & trafo = el2.trafo;

    if (el2.affine)
      {
        // the squared distance is quadratic: interior minimum or minimum on an edge
        Vec<DIMR> nv;
        if constexpr (DIMR == 3)
          nv = Cross(Vec<3>(el2.jac.Col(0)), Vec<3>(el2.jac.Col(1)));
        else
          nv = Vec<2>(-el2.jac(1,0), el2.jac(0,0));
        if (!both_sides && el2.orient * InnerProduct(n, nv) >= 0)
          return 1e99;

        Mat<DIMS,DIMS> g = Trans(el2.jac) * el2.jac;
        Vec<DIMS> r = Trans(el2.jac) * (pmaster - el2.x0);
        Vec<DIMS> lam = 0;
        double dist2 = 1e99;
        auto check = [&] (Vec<DIMS> l)
        {
          double d = L2Norm2(Vec<DIMR>(el2.x0 + el2.jac * l - pmaster));
          if (d < dist2) { dist2 = d; lam = l; }
        };

        ArrayMem<Vec<DIMS>, 3> verts;
        if constexpr (DIMS == 1)
          verts = { Vec<1>(0.), Vec<1>(1.) };
        else
          {
            verts = { Vec<2>(0., 0.), Vec<2>(1., 0.), Vec<2>(0., 1.) };
            Vec<2> l = Inv(g) * r;
            if (l(0) >= 0 && l(1) >= 0 && l(0) + l(1) <= 1)
              check(l);
          }
        for (size_t i = 0; i < verts.Size(); i++)
          {
            Vec<DIMS> a = verts[i];
            Vec<DIMS> v = verts[(i+1) % verts.Size()] - a;
            Vec<DIMS> gv = g * v;
            double t = InnerProduct(v, Vec<DIMS>(r - g * a)) / InnerProduct(v, gv);
            check(Vec<DIMS>(a + min(1., max(0., t)) * v));
          }

        if (dist2 > h*h)
          return sqrt(dist2);
        ip = lam;
        p = el2.x0 + el2.jac * lam;
        return L2Norm(pmaster - p);
      }
    Vec<DIMS> min_lam = 0;
    double min_dist = 1e99;

    min_lam = 1./(DIMS+1);

    if constexpr (DIMS==1)
    {
      auto trueDist2 = [&] (double l)
      {
        Vec<DIMR> pt;
        trafo.CalcPoint(IntegrationPoint(l), pt);
        return L2Norm2(pt - pmaster);
      };

      double lam = 0.0, d_lam = trueDist2(0.0);
      for (int k = 1; k <= 4; k++)
        {
          double l = k / 4.0, d = trueDist2(l);
          if (d < d_lam) { d_lam = d; lam = l; }
        }
      for (int it = 0; it < 12; it++)
        {
          IntegrationPoint ipl(lam);
          MappedIntegrationPoint<1, DIMR> mip{ipl, trafo};
          T2<1> t2{mip, pmaster};
          double nl = t2.CalcMinimum()(0);
          nl = min(1.0, max(0.0, nl));               // stay inside the segment
          double step = nl - lam;
          double d_new = trueDist2(nl);
          int cnt = 0;
          while (d_new > d_lam && cnt++ < 10)        // safeguard: no ascent
            {
              step *= 0.5;
              nl = lam + step;
              d_new = trueDist2(nl);
            }
          lam = nl;
          d_lam = d_new;
          if (fabs(step) < 1e-10) break;
        }

      bool is_front = both_sides;
      if (!is_front)
        {
          IntegrationPoint ipm(0.5);
          MappedIntegrationPoint<1, DIMR> mipm{ipm, trafo};
          is_front = el2.orient * InnerProduct(n, mipm.GetNV()) < 0;
        }

      for (double l : { lam, 0.0, 1.0 })
        {
          IntegrationPoint ipl(l);
          MappedIntegrationPoint<1, DIMR> mip{ipl, trafo};
          double d = L2Norm2(mip.GetPoint() - pmaster);
          if (is_front && d < min_dist)
            {
              min_dist = d;
              min_lam = l;
            }
        }
    }
    else
    // Todo: line search, stop criterion
    for([[maybe_unused]] auto i : Range(4) )
    {
      ip = min_lam;
      MappedIntegrationPoint<DIMS, DIMR> mip{ip, trafo};
      T2<DIMS> t2{mip, pmaster};
      bool is_front = el2.orient * InnerProduct(n, mip.GetNV()) < 0;
      if(both_sides)
        is_front = true;

      if constexpr (DIMS==2)
      {
        auto getDist = [&] ( auto l )
        {
          Vec<DIMR> p;
          trafo.CalcPoint( {l}, p );
          return L2Norm2( p-pmaster );
        };

        // check corner points and edges separately
        ArrayMem<Vec<DIMS>, 4> points;
        auto eltype = trafo.GetElementType();
        if(eltype==ET_TRIG)
            points = { {0,0}, {0,1}, {1,0} };
        else if(eltype==ET_QUAD)
            points = { {0,0}, {0,1}, {1,1}, {1,0} };

        for (Vec<DIMS> lam : points)
        {
          // auto d = t2(lam);
          auto d = getDist(lam);
          if(is_front && d < min_dist)
          {
            min_dist = d;
            min_lam = lam;
          }
        }

        auto checkLam = [&]( auto lam )
        {
            if(!is_front)
                return;
            auto dist = getDist(lam);

            if(dist < min_dist)
            {
                min_dist = dist;
                min_lam = lam;
            }
        };

        auto l = t2.CalcMinimum();

        if(l[0]>0 && l[1]>0 && l[0]<1 && l[1]<1 && (l[0]+l[1]<1 || eltype==ET_QUAD))
            checkLam(l);

        auto npoints = points.Size();
        for(auto i : Range(npoints))
            if(t2.CalcMinimumOnSegment( points[i], points[(i+1)%npoints], l ))
                checkLam(l);
      }
    }

    if(min_dist > h*h)
      return sqrt(min_dist);

    ip = min_lam;
    trafo.CalcPoint( ip, p );
    return L2Norm(pmaster-p);
  }

  template<int DIM>
  DeformedBndElement<DIM> T_GapFunction<DIM> :: GetDeformedBndElement(const ElementTransformation & trafo, LocalHeap & lh) const
  {
    auto ei = trafo.GetElementId();
    DeformedBndElement<DIM> el { ei, trafo.AddDeformation(displacement.get(), lh),
                                 GetDomIn(*ma, ma->GetElement(ei)) == 0 ? -1. : 1. };
    el.affine = trafo.GetElementType() == (DIM == 3 ? ET_TRIG : ET_SEGM) && !trafo.IsCurvedElement()
      && (!displacement || displacement->GetFESpace()->GetFE(ei, lh).Order() <= 1);
    if (el.affine)
      {
        el.trafo.CalcPoint(IntegrationPoint(0., 0.), el.x0);
        for (int j = 0; j < DIM-1; j++)
          {
            Vec<DIM> pj;
            el.trafo.CalcPoint(IntegrationPoint(j == 0 ? 1. : 0., j == 1 ? 1. : 0.), pj);
            for (int k = 0; k < DIM; k++)
              el.jac(k, j) = pj(k) - el.x0(k);
          }
      }
    return el;
  }

  template<int DIM>
  auto T_GapFunction<DIM> :: FindClosest(const ElementTransformation & trafo1, const IntegrationPoint & ip1,
                                         LocalHeap & lh, bool both_sides) const -> optional<ClosestPoint>
  {
    auto bel1 = GetDeformedBndElement(trafo1, lh);
    const auto & mip1 = static_cast<const MappedIntegrationPoint<DIM-1, DIM>&>(bel1.trafo(ip1, lh));
    Vec<DIM> p1 = mip1.GetPoint();
    Vec<DIM> n1 = bel1.orient * mip1.GetNV();

    double mindist = h;
    optional<ClosestPoint> closest;

    // find all bound-2 elements closer to p1 than h
    netgen::Point<DIM> ngp1;
    for (int j = 0; j < DIM; j++) ngp1(j) = p1(j);

    auto hcurrent = h/(1024.*1024.);
    int found = 2;
    while(found>0 && hcurrent<=h)
      {
        netgen::Box<DIM> box(ngp1, ngp1);
        box.Increase(hcurrent);

        searchtree->GetFirstIntersecting
          (box.PMin(), box.PMax(),
           [&] (int elnr2)
           {
             HeapReset hr(lh);
             auto bel2 = GetDeformedBndElement(ma->GetTrafo(ElementId(BND, elnr2), lh), lh);
             IntegrationPoint ip2;
             Vec<DIM> p2;
             double dist = FindClosestPoint<DIM-1,DIM>(p1, n1, mindist, bel2, ip2, p2, both_sides);
             if(dist < mindist)
               {
                 mindist = dist;
                 closest = ClosestPoint{elnr2, ip2, p1, p2};
               }
             return false;
           });

        if(mindist < hcurrent)
          found--;
        hcurrent *= 2;
      }
    return closest;
  }

  template<int DIM>
  optional<ContactPair<DIM>> T_GapFunction<DIM> :: CreateContactPair(const MappedIntegrationPoint<DIM-1, DIM>& mip1, LocalHeap& lh, bool both_sides) const
  {
    HeapReset hr(lh);
    auto & trafo1 = mip1.GetTransformation();
    auto closest = FindClosest(trafo1, mip1.IP(), lh, both_sides);
    if (!closest)
      return nullopt;
    return ContactPair<DIM>{trafo1.GetElementId(), ElementId(BND, closest->el2),
        mip1.IP(), closest->ip2};
  }

  template<int DIM>
  void T_GapFunction<DIM> :: Update(shared_ptr<GridFunction> displacement_, int intorder2, double h_, bool both_sides)
  {
    h = h_;
    this->both_sides = both_sides;

    displacement = displacement_;

    auto & mask = other.Mask();
    Array<netgen::Box<DIM>> elboxes(ma->GetNE(BND));
    Array<double> maxh_thread(TaskManager::GetMaxThreads());
    maxh_thread = 0.0;
    LocalHeap clh(32*1000*1000, "T_GapFunction::Update");
    ma->IterateElements
      (BND, clh, [&] (Ngs_Element el2, LocalHeap & lh)
       {
         if (!mask.Test(el2.GetIndex())) return;
         auto & trafo2 = ma->GetTrafo (el2, lh);
         auto & trafo2_def = trafo2.AddDeformation(displacement.get(), lh);
         IntegrationRule ir2_(trafo2.GetElementType(), intorder2);
         netgen::Box<DIM> elbox{netgen::Box<DIM>::EMPTY_BOX};
         double & maxh = maxh_thread[TaskManager::GetThreadId()];
         for(auto ir2 : ir2_.Split())
           {
             HeapReset hr(lh);
             MappedIntegrationRule<DIM-1, DIM> mir2_def(ir2, trafo2_def, lh);
             netgen::Box<DIM> partbox{netgen::Box<DIM>::EMPTY_BOX};
             for (auto & mip : mir2_def)
               {
                 netgen::Point<DIM> p;
                 for (int j = 0; j < DIM; j++)
                   p(j) = mip.GetPoint()(j);
                 partbox.Add(p);
               }
             maxh = max(maxh, partbox.Diam());
             elbox.Add(partbox.PMin());
             elbox.Add(partbox.PMax());
           }
         elboxes[el2.Nr()] = elbox;
       });

    netgen::Box<DIM> bbox{netgen::Box<DIM>::EMPTY_BOX};
    double maxh = 0;
    for (auto v : maxh_thread) maxh = max(maxh, v);
    for (Ngs_Element el2 : ma->Elements(BND))
      if (mask.Test(el2.GetIndex()))
        {
          bbox.Add(elboxes[el2.Nr()].PMin());
          bbox.Add(elboxes[el2.Nr()].PMax());
        }

    // Default-value for h is 2 * maximum_element_diameter
    if(h==0.0)
      h = 2*maxh;

    bbox.Scale(2); // make sure we don't add boxes outside of tree bounding box
    searchtree = make_shared<netgen::BoxTree<DIM, int>>(bbox);
    for (Ngs_Element el2 : ma->Elements(BND))
      {
        if (!mask.Test(el2.GetIndex())) continue;
        auto elbox = elboxes[el2.Nr()];
        elbox.Scale(1.1);
        searchtree->Insert(elbox, el2.Nr());
      }
  }

  template<int DIM>
  void T_GapFunction<DIM> :: Evaluate(const BaseMappedIntegrationPoint & ip,
                                      FlatVector<> result) const
  {
    LocalHeapMem<100000> lh("gapfunction");
    auto & trafo1 = ip.GetTransformation();
    result = 0;
    if (!master.Mask().Test(ma->GetElement(trafo1.GetElementId()).GetIndex())) return;

    result = std::numeric_limits<double>::max();
    if (auto closest = FindClosest(trafo1, ip.IP(), lh, both_sides))
      result = closest->p2 - closest->p1;
  }

  template<int DIM>
  void T_GapFunction<DIM> :: Evaluate(const BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<> hresult) const
  {
    auto result = hresult.AddSize(mir.Size(), Dimension());
    for (auto i : Range(mir))
      Evaluate(mir[i], result.Row(i));
  }

  template<int DIM>
  void T_GapFunction<DIM> :: CreateClipGeometry(const ElementTransformation & trafo1, int intorder, LocalHeap & lh, bool both_sides,
                                                ClipGeometry & geo) const
  {
    geo.pieces.SetSize0();
    geo.rim1.SetSize0();
    geo.rim2.SetSize0();
    geo.rimel.SetSize0();
    HeapReset hr(lh);
    auto bel1 = GetDeformedBndElement(trafo1, lh);
    auto & trafo1_def = bel1.trafo;
    IntegrationPoint ipc = DIM == 3 ? IntegrationPoint(1./3, 1./3) : IntegrationPoint(0.5);
    MappedIntegrationPoint<DIM-1, DIM> mip1(ipc, trafo1_def);
    Vec<DIM> nv1 = bel1.orient * mip1.GetNV();

    if constexpr (DIM == 3)
      {
        if (trafo1.GetElementType() != ET_TRIG)
          throw Exception("ContactBoundary quadrature='clip': only triangle boundary elements supported");
        Vec<3> P0, P1, P2;
        trafo1_def.CalcPoint(IntegrationPoint(0., 0.), P0);
        trafo1_def.CalcPoint(IntegrationPoint(1., 0.), P1);
        trafo1_def.CalcPoint(IntegrationPoint(0., 1.), P2);
        Vec<3> e1 = P1 - P0, e2 = P2 - P0;
        Mat<2,2> G = { { InnerProduct(e1, e1), InnerProduct(e1, e2) },
                       { InnerProduct(e1, e2), InnerProduct(e2, e2) } };
        Mat<2,2> Ginv = Inv(G);
        Vec<3> n1 = Cross(e1, e2);
        n1 /= L2Norm(n1);
        auto ref1 = [&] (Vec<3> x)
        {
          Vec<2> r { InnerProduct(e1, x - P0), InnerProduct(e2, x - P0) };
          return Vec<2>(Ginv * r);
        };

        netgen::Box<3> box(netgen::Box<3>::EMPTY_BOX);
        for (auto & P : { P0, P1, P2 })
          box.Add(netgen::Point<3>(P(0), P(1), P(2)));
        box.Increase(h);
        IntegrationRule irtrig(ET_TRIG, intorder);
        double covered = 0;

        searchtree->GetFirstIntersecting
          (box.PMin(), box.PMax(),
           [&] (int elnr2)
           {
             HeapReset hr2(lh);
             auto & trafo2 = ma->GetTrafo(ElementId(BND, elnr2), lh);
             if (trafo2.GetElementType() != ET_TRIG)
               throw Exception("ContactBoundary quadrature='clip': only triangle boundary elements supported");
             auto bel2 = GetDeformedBndElement(trafo2, lh);
             auto & trafo2_def = bel2.trafo;
             MappedIntegrationPoint<2,3> mip2(ipc, trafo2_def);
             if (!both_sides && bel2.orient * InnerProduct(nv1, mip2.GetNV()) >= 0)
               return false;
             Vec<3> Q0, Q1, Q2;
             trafo2_def.CalcPoint(IntegrationPoint(0., 0.), Q0);
             trafo2_def.CalcPoint(IntegrationPoint(1., 0.), Q1);
             trafo2_def.CalcPoint(IntegrationPoint(0., 1.), Q2);
             for (auto & Q : { Q0, Q1, Q2 })
               if (fabs(InnerProduct(n1, Q - P0)) > h)
                 return false;

             // Sutherland-Hodgman against s >= 0, t >= 0, 1-s-t >= 0
             ArrayMem<Vec<2>, 16> poly, next;
             poly.Append(ref1(Q0)); poly.Append(ref1(Q1)); poly.Append(ref1(Q2));
             for (int edge = 0; edge < 3; edge++)
               {
                 auto f = [edge] (Vec<2> v) { return edge == 0 ? v(0) : edge == 1 ? v(1) : 1 - v(0) - v(1); };
                 next.SetSize0();
                 for (size_t k = 0; k < poly.Size(); k++)
                   {
                     Vec<2> a = poly[k], b = poly[(k + 1) % poly.Size()];
                     double fa = f(a), fb = f(b);
                     if (fa >= 0) next.Append(a);
                     if ((fa >= 0) != (fb >= 0))
                       next.Append(a + fa / (fa - fb) * (b - a));
                   }
                 poly.SetSize0();
                 for (auto & v : next)
                   poly.Append(v);
                 if (poly.Size() < 3)
                   return false;
               }

             Vec<3> f1 = Q1 - Q0, f2 = Q2 - Q0;
             Mat<2,2> G2 = { { InnerProduct(f1, f1), InnerProduct(f1, f2) },
                             { InnerProduct(f1, f2), InnerProduct(f2, f2) } };
             Mat<2,2> G2inv = Inv(G2);
             Mat<2,2> F = { { InnerProduct(f1, e1), InnerProduct(f1, e2) },
                            { InnerProduct(f2, e1), InnerProduct(f2, e2) } };
             Vec<2> d { InnerProduct(f1, P0 - Q0), InnerProduct(f2, P0 - Q0) };
             for (size_t k = 1; k + 1 < poly.Size(); k++)
               {
                 Vec<2> B = poly[k] - poly[0], C = poly[k+1] - poly[0];
                 covered += 0.5 * fabs(B(0) * C(1) - B(1) * C(0));
               }
             geo.pieces.Append(ClipPiece{elnr2, Array<Vec<2>>(poly), G2inv * d, G2inv * F});
             return false;
           });

        // primary area no secondary facet covers (rim of non-matching discretized
        // geometry): closest-point pairing of the Gauss points outside all pieces,
        // weighted to the uncovered measure, keeps the primary side fully coupled
        double uncovered = 0.5 - covered;
        if (uncovered > 1e-10)
          {
            auto inside = [] (const Array<Vec<2>> & poly, Vec<2> p)
            {
              int sign = 0;
              for (size_t k = 0; k < poly.Size(); k++)
                {
                  Vec<2> a = poly[k], b = poly[(k + 1) % poly.Size()];
                  double c = (b(0) - a(0)) * (p(1) - a(1)) - (b(1) - a(1)) * (p(0) - a(0));
                  int sc = c > 1e-14 ? 1 : (c < -1e-14 ? -1 : 0);
                  if (sc == 0) continue;
                  if (sign == 0) sign = sc;
                  else if (sc != sign) return false;
                }
              return true;
            };
            // the gap is a thin strip along the facet boundary: sample there too
            Array<Vec<2>> samples;
            for (auto & ip : irtrig)
              samples.Append(Vec<2>(ip(0), ip(1)));
            Vec<2> corners[3] = { Vec<2>(1, 0), Vec<2>(0, 1), Vec<2>(0, 0) };
            Vec<2> center(1./3, 1./3);
            for (int e = 0; e < 3; e++)
              for (int k = 0; k < 16; k++)
                {
                  Vec<2> q = corners[e] + (k + 0.5) / 16 * (corners[(e + 1) % 3] - corners[e]);
                  samples.Append(q + 1e-6 * (center - q));
                }
            ArrayMem<Vec<2>, 64> outside;
            for (auto & p : samples)
              {
                bool in = false;
                for (auto & piece : geo.pieces)
                  if (inside(piece.poly, p)) { in = true; break; }
                if (!in) outside.Append(p);
              }
            for (auto & p : outside)
              {
                HeapReset hr3(lh);
                IntegrationPoint ip(p(0), p(1), 0, uncovered / outside.Size());
                MappedIntegrationPoint<2,3> mip(ip, trafo1);
                auto pair = CreateContactPair(mip, lh, both_sides);
                if (!pair) continue;
                geo.rim1.Append(ip);
                geo.rim2.Append(pair->secondary_ip);
                geo.rimel.Append(pair->secondary_el.Nr());
              }
          }
      }
    else
      {
        Vec<2> P0, P1;
        trafo1_def.CalcPoint(IntegrationPoint(0.), P0);
        trafo1_def.CalcPoint(IntegrationPoint(1.), P1);
        Vec<2> e = P1 - P0;
        double ee = InnerProduct(e, e);
        Vec<2> n1 { -e(1), e(0) };
        n1 /= L2Norm(n1);
        netgen::Box<2> box(netgen::Box<2>::EMPTY_BOX);
        box.Add(netgen::Point<2>(P0(0), P0(1)));
        box.Add(netgen::Point<2>(P1(0), P1(1)));
        box.Increase(h);
        IntegrationRule irseg(ET_SEGM, intorder);
        double covered = 0;

        searchtree->GetFirstIntersecting
          (box.PMin(), box.PMax(),
           [&] (int elnr2)
           {
             HeapReset hr2(lh);
             auto bel2 = GetDeformedBndElement(ma->GetTrafo(ElementId(BND, elnr2), lh), lh);
             auto & trafo2_def = bel2.trafo;
             MappedIntegrationPoint<1,2> mip2(ipc, trafo2_def);
             if (!both_sides && bel2.orient * InnerProduct(nv1, mip2.GetNV()) >= 0)
               return false;
             Vec<2> Q0, Q1;
             trafo2_def.CalcPoint(IntegrationPoint(0.), Q0);
             trafo2_def.CalcPoint(IntegrationPoint(1.), Q1);
             if (fabs(InnerProduct(n1, Q0 - P0)) > h || fabs(InnerProduct(n1, Q1 - P0)) > h)
               return false;
             double a = InnerProduct(e, Q0 - P0) / ee, b = InnerProduct(e, Q1 - P0) / ee;
             double lo = max(0., min(a, b)), hi = min(1., max(a, b));
             if (hi - lo < 1e-14)
               return false;
             covered += hi - lo;
             Vec<2> f = Q1 - Q0;
             double ff = InnerProduct(f, f);
             Mat<2,2> M2 = 0.;
             M2(0,0) = InnerProduct(f, e) / ff;
             Array<Vec<2>> iv { Vec<2>(lo, 0.), Vec<2>(hi, 0.) };
             geo.pieces.Append(ClipPiece{elnr2, std::move(iv), Vec<2>(InnerProduct(f, P0 - Q0) / ff, 0.), M2});
             return false;
           });

        double uncovered = 1 - covered;
        if (uncovered > 1e-10)
          {
            // the gaps sit at the segment ends: sample there too
            Array<double> samples;
            for (auto & ip : irseg)
              samples.Append(ip(0));
            for (int k = 0; k < 16; k++)
              {
                samples.Append(1e-6 + 0.05 * k / 16);
                samples.Append(1 - 1e-6 - 0.05 * k / 16);
              }
            ArrayMem<double, 64> outside;
            for (double t : samples)
              {
                bool in = false;
                for (auto & piece : geo.pieces)
                  if (t >= piece.poly[0](0) && t <= piece.poly[1](0)) { in = true; break; }
                if (!in) outside.Append(t);
              }
            for (double t : outside)
              {
                HeapReset hr3(lh);
                IntegrationPoint ip(t, 0, 0, uncovered / outside.Size());
                MappedIntegrationPoint<1,2> mip(ip, trafo1);
                auto pair = CreateContactPair(mip, lh, both_sides);
                if (!pair) continue;
                geo.rim1.Append(ip);
                geo.rim2.Append(pair->secondary_ip);
                geo.rimel.Append(pair->secondary_el.Nr());
              }
          }
      }
  }

  // Gauss rule of order intorder on every clipped piece; primary points carry reference weights
  static void ClipRules(const ClipGeometry & geo, int dim, int intorder,
                        Array<IntegrationPoint> & ip1s, Array<IntegrationPoint> & ip2s, Array<int> & el2s)
  {
    if (dim == 3)
      {
        IntegrationRule irtrig(ET_TRIG, intorder);
        for (auto & piece : geo.pieces)
          {
            auto & poly = piece.poly;
            for (size_t k = 1; k + 1 < poly.Size(); k++)
              {
                Vec<2> A = poly[0], B = poly[k] - poly[0], C = poly[k+1] - poly[0];
                double det = fabs(B(0) * C(1) - B(1) * C(0));
                if (det < 1e-14)
                  continue;
                for (auto & ip : irtrig)
                  {
                    Vec<2> st = A + ip(0) * B + ip(1) * C;
                    Vec<2> st2 = piece.c2 + piece.M2 * st;
                    ip1s.Append(IntegrationPoint(st(0), st(1), 0, ip.Weight() * det));
                    ip2s.Append(IntegrationPoint(st2(0), st2(1), 0, 0));
                    el2s.Append(piece.el2);
                  }
              }
          }
      }
    else
      {
        IntegrationRule irseg(ET_SEGM, intorder);
        for (auto & piece : geo.pieces)
          {
            double lo = piece.poly[0](0), hi = piece.poly[1](0);
            for (auto & ip : irseg)
              {
                double s1 = lo + ip(0) * (hi - lo);
                double s2 = piece.c2(0) + piece.M2(0,0) * s1;
                ip1s.Append(IntegrationPoint(s1, 0, 0, ip.Weight() * (hi - lo)));
                ip2s.Append(IntegrationPoint(s2, 0, 0, 0));
                el2s.Append(piece.el2);
              }
          }
      }
    for (size_t i = 0; i < geo.rim1.Size(); i++)
      {
        ip1s.Append(geo.rim1[i]);
        ip2s.Append(geo.rim2[i]);
        el2s.Append(geo.rimel[i]);
      }
  }

  template class T_GapFunction<2>;
  template class T_GapFunction<3>;

  template<int DIM>
  void DisplacedNormal<DIM>::Evaluate(const BaseMappedIntegrationPoint& ip, FlatVector<> values) const
  {
    auto ma = displacement->GetFESpace()->GetMeshAccess();
    const auto& el = ma->GetElement(ip.GetTransformation().GetElementId());
    double inv_fac = GetDomIn(*ma, el) == 0 ? -1. : 1.;

    if(!displacement)
      {
        values = inv_fac * static_cast<const MappedIntegrationPoint<DIM-1, DIM>&>(ip).GetNV();
        return;
      }
    LocalHeapMem<10000> lh("deformednormal");
    auto& trafo_def = ip.GetTransformation().AddDeformation(displacement.get(), lh);
    auto& mip_def = static_cast<const MappedIntegrationPoint<DIM-1, DIM>&>(trafo_def(ip.IP(), lh));
    values = inv_fac * mip_def.GetNV();
  }

  template class DisplacedNormal<2>;
  template class DisplacedNormal<3>;

  ContactEnergy::ContactEnergy(shared_ptr<CoefficientFunction> _cf,
                               bool _deformed)
    : cf(_cf), deformed(_deformed)
  {
    cf->TraverseDAG
      ([&](CoefficientFunction& nodecf)
       {
         auto proxy = dynamic_cast<ProxyFunction*>(&nodecf);
         if (proxy && !proxy->IsTestFunction() && !trial_proxies.Contains(proxy))
           trial_proxies.Append(proxy);
         auto gf = dynamic_cast<GridFunction*>(&nodecf);
         if (gf && !cf_gridfunctions.Contains(&nodecf))
           cf_gridfunctions.Append(&nodecf);
       });
    fes = trial_proxies[0]->GetFESpace();
  }

  double ContactEnergy::CalcEnergy(const FiniteElement& primary_fel,
                                   const FiniteElement& secondary_fel,
                                   const BaseMappedIntegrationRule& primary_mir,
                                   FlatVector<double> elx,
                                   LocalHeap& lh)
  {
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(primary_mir.GetTransformation()).userdata = &ud;
    ud.fel = & primary_fel;

    for(auto proxy : trial_proxies)
      {
        IntRange trial_range = proxy->IsOther() ?
          IntRange(proxy->Evaluator()->BlockDim() * primary_fel.GetNDof(), elx.Size())
          : IntRange(0, proxy->Evaluator()->BlockDim() * primary_fel.GetNDof());
        ud.AssignMemory(proxy, primary_mir.Size(), proxy->Dimension(), lh);
        if(proxy->IsOther())
          proxy->Evaluator()->Apply(secondary_fel, *primary_mir.GetOtherMIR(), elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(primary_fel, primary_mir, elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
      }

    FlatMatrix<> values(primary_mir.Size(), 1, lh);
    cf->Evaluate(primary_mir, values);

    double sum = 0.;
    for (int i = 0; i < primary_mir.Size(); i++)
      sum += primary_mir[i].GetWeight() * values(i,0);
    return sum;
  }

  void ContactEnergy::ApplyAdd(const FiniteElement& primary_fel,
                               const FiniteElement& secondary_fel,
                               const BaseMappedIntegrationRule& primary_mir,
                               FlatVector<double> elx,
                               FlatVector<double> ely,
                               LocalHeap& lh)
  {
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(primary_mir.GetTransformation()).userdata = &ud;
    ud.fel = & primary_fel;

    for(auto proxy : trial_proxies)
      {
        IntRange trial_range = proxy->IsOther() ?
          IntRange(proxy->Evaluator()->BlockDim() * primary_fel.GetNDof(), elx.Size())
          : IntRange(0, proxy->Evaluator()->BlockDim() * primary_fel.GetNDof());
        ud.AssignMemory(proxy, primary_mir.Size(), proxy->Dimension(), lh);
        if(proxy->IsOther())
          proxy->Evaluator()->Apply(secondary_fel, *primary_mir.GetOtherMIR(), elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(primary_fel, primary_mir, elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
      }

    FlatVector<> ely1(ely.Size(), lh);
    FlatMatrix<AutoDiff<1>> dval(primary_mir.Size(), 1, lh);

    for (auto proxy : trial_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(primary_mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            cf -> Evaluate (primary_mir, dval);
            for (size_t i = 0; i < primary_mir.Size(); i++)
              proxyvalues(i,k) = dval(i,0).DValue(0);
          }

        for (int i = 0; i < primary_mir.Size(); i++)
          proxyvalues.Row(i) *= primary_mir[i].GetWeight();

        IntRange test_range  = proxy->IsOther() ?
          IntRange(proxy->Evaluator()->BlockDim()*primary_fel.GetNDof(), ely.Size())
          : IntRange(0, proxy->Evaluator()->BlockDim()*primary_fel.GetNDof());
        ely1 = 0.;
        if(proxy->IsOther())
          proxy->Evaluator()->ApplyTrans(secondary_fel, *primary_mir.GetOtherMIR(), proxyvalues, ely1.Range(test_range), lh);
        else
          proxy->Evaluator()->ApplyTrans(primary_fel, primary_mir, proxyvalues, ely1.Range(test_range), lh);
        ely += ely1;
      }
  }

  void ContactEnergy::CalcLinearizedAdd(const FiniteElement& primary_fel,
                                        const FiniteElement& secondary_fel,
                                        const BaseMappedIntegrationRule& primary_mir,
                                        FlatVector<double> elx,
                                        FlatMatrix<double> elmat,
                                        LocalHeap& lh)
  {
    HeapReset hr(lh);
    ProxyUserData ud(trial_proxies.Size(), cf_gridfunctions.Size(), lh);
    const_cast<ElementTransformation&>(primary_mir.GetTransformation()).userdata = &ud;
    ud.fel = & primary_fel;

    for(auto proxy : trial_proxies)
      {
        IntRange trial_range = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim() * primary_fel.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim() * primary_fel.GetNDof());
        ud.AssignMemory(proxy, primary_mir.Size(), proxy->Dimension(), lh);
        if(proxy->IsOther())
          proxy->Evaluator()->Apply(secondary_fel, *primary_mir.GetOtherMIR(), elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(primary_fel, primary_mir, elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
      }

    for(auto gf : cf_gridfunctions)
        ud.AssignMemory(gf, primary_mir.Size(), gf->Dimension(), lh);

    FlatMatrix<> dderiv(primary_mir.Size(), 1,lh);
    FlatMatrix<AutoDiffDiff<1,double>> ddval(primary_mir.Size(), 1, lh);

    FlatArray<FlatMatrix<>> diags(trial_proxies.Size(), lh);
    for (int k1 : Range(trial_proxies))
      {
        auto proxy = trial_proxies[k1];
        diags[k1].AssignMemory(primary_mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (primary_mir, ddval);
            for (size_t i = 0; i < primary_mir.Size(); i++)
              diags[k1](i,k) = ddval(i,0).DDValue(0);
          }
      }

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(trial_proxies))
        {
          HeapReset hr(lh);

          auto proxy1 = trial_proxies[k1];
          auto proxy2 = trial_proxies[l1];

          FlatTensor<3> proxyvalues(lh, primary_mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;

                {
                  cf -> Evaluate (primary_mir, ddval);
                  for (size_t i = 0; i < primary_mir.Size(); i++)
                    dderiv(i,0) = ddval(i,0).DDValue(0);
                }
                proxyvalues(STAR,l,k) = dderiv.Col(0);

                if (proxy1 != proxy2 || k != l)  // computed mixed second derivatives
                  {
                    proxyvalues(STAR,l,k) -= diags[k1].Col(k);
                    proxyvalues(STAR,l,k) -= diags[l1].Col(l);
                    proxyvalues(STAR,l,k) *= 0.5;
                  }
              }

          for (int i = 0; i < primary_mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= primary_mir[i].GetWeight();

          IntRange trial_range  = proxy1->IsOther() ? IntRange(proxy1->Evaluator()->BlockDim()*primary_fel.GetNDof(), elmat.Width()) : IntRange(0, proxy1->Evaluator()->BlockDim()*primary_fel.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(proxy2->Evaluator()->BlockDim()*primary_fel.GetNDof(), elmat.Height()) : IntRange(0, proxy2->Evaluator()->BlockDim()*primary_fel.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);


          int bs = primary_mir.Size();
          FlatMatrix<double,ColMajor> bdbmat1(bs*proxy2->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bbmat2(bs*proxy2->Dimension(), loc_elmat.Height(), lh);

          bdbmat1 = 0.;
          bbmat2 = 0.;

          const auto& s_mir = *primary_mir.GetOtherMIR();

          for (size_t j = 0; j < bs; j++)
            {
              size_t ii = j;
              IntRange r3 = proxy2->Dimension() * IntRange(j,j+1);

              {
                if(proxy1->IsOther())
                  proxy1->Evaluator()->CalcMatrix(secondary_fel, s_mir[ii], bmat1, lh);
                else
                  proxy1->Evaluator()->CalcMatrix(primary_fel, primary_mir[ii], bmat1, lh);
                if(proxy2->IsOther())
                  proxy2->Evaluator()->CalcMatrix(secondary_fel, s_mir[ii], bmat2, lh);
                else
                  proxy2->Evaluator()->CalcMatrix(primary_fel, primary_mir[ii], bmat2, lh);
              }

              bdbmat1.Rows(r3) = proxyvalues(ii,STAR,STAR) * bmat1;
              bbmat2.Rows(r3) = bmat2;
            }
          loc_elmat += Trans(bbmat2) * bdbmat1;
        }
  }

  ContactIntegrator::ContactIntegrator(shared_ptr<CoefficientFunction> _cf,
                                       bool _deformed)
    : cf(_cf), deformed(_deformed)
  {
    cf->TraverseDAG
      ([&](CoefficientFunction& nodecf)
       {
         auto proxy = dynamic_cast<ProxyFunction*>(&nodecf);
         if (proxy && !proxy->IsTestFunction() && !trial_proxies.Contains(proxy))
           trial_proxies.Append(proxy);

         if(proxy && proxy->IsTestFunction() && !test_proxies.Contains(proxy))
           test_proxies.Append(proxy);

         auto gf = dynamic_cast<GridFunction*>(&nodecf);
         if (gf && !cf_gridfunctions.Contains(&nodecf))
           cf_gridfunctions.Append(&nodecf);
       });
    if(trial_proxies.Size())
      fes = trial_proxies[0]->GetFESpace();
    else
      {
        if (test_proxies.Size())
          fes = test_proxies[0]->GetFESpace();
        else
          throw Exception("No trial or test function found in ContactIntegrator");
      }

    dcf_dtest.SetSize(test_proxies.Size());
    ddcf_dtest_dtrial.SetSize(test_proxies.Size(), trial_proxies.Size());
    for (int i = 0; i < test_proxies.Size(); i++)
      {
        try
          {
            CoefficientFunction::T_DJC cache;
            dcf_dtest[i] = cf->DiffJacobi(test_proxies[i], cache);
          }
        catch (const Exception& e)
          {
            cout << IM(5) << "dcf_dtest has thrown exception " << e.What() << endl;
          }
          for (auto j : Range(trial_proxies))
            {
              if (dcf_dtest[i])
                try
                  {
                    CoefficientFunction::T_DJC cache2;
                      ddcf_dtest_dtrial(i, j) = dcf_dtest[i]->DiffJacobi(trial_proxies[j], cache2);
                  }
                catch (const Exception& e)
                  {
                    cout << IM(1) << "ddcf_dtest_dtrial has thrown exception " << e.What() << endl;
                  }
            }
      }

  }

  void ContactIntegrator::ApplyAdd(const FiniteElement& primary_trial_fel,
                                   const FiniteElement& secondary_trial_fel,
                                   const FiniteElement& primary_test_fel,
                                   const FiniteElement& secondary_test_fel,
                                   const BaseMappedIntegrationRule& primary_mir,
                                   FlatVector<double> elx,
                                   FlatVector<double> ely,
                                   LocalHeap& lh)
  {
    HeapReset hr(lh);
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(primary_mir.GetTransformation()).userdata = &ud;
    ud.fel = & primary_trial_fel;

    for(auto proxy : trial_proxies)
      {
        IntRange trial_range = proxy->IsOther()
          ? IntRange(proxy->Evaluator()->BlockDim() * primary_trial_fel.GetNDof(), elx.Size())
          : IntRange(0, proxy->Evaluator()->BlockDim() * primary_trial_fel.GetNDof());
        ud.AssignMemory(proxy, primary_mir.Size(), proxy->Dimension(), lh);
        if(proxy->IsOther())
          proxy->Evaluator()->Apply(secondary_trial_fel, *primary_mir.GetOtherMIR(), elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(primary_trial_fel, primary_mir, elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
      }

    FlatVector<> ely1(ely.Size(), lh);
    FlatMatrix<> val(primary_mir.Size(), 1, lh);

    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(primary_mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (primary_mir, val);
            proxyvalues.Col(k) = val.Col(0);
          }

        for (int i = 0; i < primary_mir.Size(); i++)
          proxyvalues.Row(i) *= primary_mir[i].GetWeight();

        IntRange test_range  = proxy->IsOther()
          ? IntRange(proxy->Evaluator()->BlockDim()*primary_test_fel.GetNDof(), ely.Size())
          : IntRange(0, proxy->Evaluator()->BlockDim()*primary_test_fel.GetNDof());
        ely1 = 0.;
        if(proxy->IsOther())
          proxy->Evaluator()->ApplyTrans(secondary_test_fel, *primary_mir.GetOtherMIR(), proxyvalues, ely1.Range(test_range), lh);
        else
          proxy->Evaluator()->ApplyTrans(primary_test_fel, primary_mir, proxyvalues, ely1.Range(test_range), lh);
        ely += ely1;
      }
  }

  void ContactIntegrator::CalcLinearizedAdd(const FiniteElement& primary_trial_fel,
                                            const FiniteElement& secondary_trial_fel,
                                            const FiniteElement& primary_test_fel,
                                            const FiniteElement& secondary_test_fel,
                                            const BaseMappedIntegrationRule& primary_mir,
                                            FlatVector<double> elx,
                                            FlatMatrix<double> elmat,
                                            LocalHeap& lh)
  {
    static Timer t("ContactIntegrator::CalcLinearizedAdd");
    static Timer teval("ContactIntegrator::CalcLinearizedAdd-Evaluate");
    static Timer tbdb("ContactIntegrator::CalcLinearizedAdd-BDB");
    // RegionTimer rt(t);
    HeapReset hr(lh);
    ProxyUserData ud(trial_proxies.Size(), cf_gridfunctions.Size(), lh);
    ProxyUserData ud2(trial_proxies.Size(), cf_gridfunctions.Size(), lh);
    auto & secondary_mir = *primary_mir.GetOtherMIR();
    const_cast<ElementTransformation&>(primary_mir.GetTransformation()).userdata = &ud;
    const_cast<ElementTransformation&>(secondary_mir.GetTransformation()).userdata = &ud2;
    ud.fel = & primary_trial_fel;
    ud2.fel = & secondary_trial_fel;

    for(auto proxy : trial_proxies)
      {
        IntRange trial_range = proxy->IsOther()
          ? IntRange(proxy->Evaluator()->BlockDim() * primary_trial_fel.GetNDof(), elx.Size())
          : IntRange(0, proxy->Evaluator()->BlockDim() * primary_trial_fel.GetNDof());
        ud.AssignMemory(proxy, primary_mir.Size(), proxy->Dimension(), lh);
        if(proxy->IsOther())
          proxy->Evaluator()->Apply(secondary_trial_fel, secondary_mir, elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(primary_trial_fel, primary_mir, elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
      }

    for(auto gf : cf_gridfunctions)
      {
        ud.AssignMemory(gf, primary_mir.Size(), gf->Dimension(), lh);
        ud2.AssignMemory(gf, secondary_mir.Size(), gf->Dimension(), lh);
      }

    FlatMatrix<> dderiv(primary_mir.Size(), 1,lh);
    FlatMatrix<AutoDiff<1>> dval(primary_mir.Size(), 1, lh);

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies)) // ss
        {
          HeapReset hr(lh);

          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3> proxyvalues(lh, primary_mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          FlatMatrix<> proxyvalues2(primary_mir.Size(), proxy2->Dimension()*proxy1->Dimension(), lh);

          if (ddcf_dtest_dtrial(l1, k1))
            {
              // RegionTimer rteval(teval);
              ddcf_dtest_dtrial(l1, k1)->Evaluate(primary_mir, proxyvalues2);
              for (int k = 0; k < proxy1->Dimension(); k++)
                for (int l = 0; l < proxy2->Dimension(); l++)
                  for (size_t i = 0; i < primary_mir.Size(); i++)
                    proxyvalues(i,l,k) = proxyvalues2(i, l*proxy1->Dimension() + k);
            }
          else
          {
            // RegionTimer rteval(teval);
            for (int k = 0; k < proxy1->Dimension(); k++)
              for (int l = 0; l < proxy2->Dimension(); l++)
                {
                  ud.trialfunction = proxy1;
                  ud.trial_comp = k;
                  ud.testfunction = proxy2;
                  ud.test_comp = l;

                  cf -> Evaluate (primary_mir, dval);
                  for (size_t i = 0; i < primary_mir.Size(); i++)
                    proxyvalues(i,l,k) = dval(i,0).DValue(0);
                }
          }

          for (int i = 0; i < primary_mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= primary_mir[i].GetWeight();

          IntRange trial_range  = proxy1->IsOther()
            ? IntRange(proxy1->Evaluator()->BlockDim()*primary_trial_fel.GetNDof(), elmat.Width())
            : IntRange(0, proxy1->Evaluator()->BlockDim()*primary_trial_fel.GetNDof());
          IntRange test_range  = proxy2->IsOther()
            ? IntRange(proxy2->Evaluator()->BlockDim()*primary_test_fel.GetNDof(), elmat.Height())
            : IntRange(0, proxy2->Evaluator()->BlockDim()*primary_test_fel.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), trial_range.Size(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), test_range.Size(), lh);

          int bs = primary_mir.Size();
          FlatMatrix<double,ColMajor> bdbmat1(bs*proxy2->Dimension(), trial_range.Size(), lh);
          FlatMatrix<double,ColMajor> bbmat2(bs*proxy2->Dimension(), test_range.Size(), lh);

          bdbmat1 = 0.;
          bbmat2 = 0.;

          // RegionTimer rtbdb(tbdb);
          const auto& s_mir = *primary_mir.GetOtherMIR();

          for (size_t j = 0; j < bs; j++)
            {
              size_t ii = j;
              IntRange r3 = proxy2->Dimension() * IntRange(j,j+1);

              {
                if(proxy1->IsOther())
                  proxy1->Evaluator()->CalcMatrix(secondary_trial_fel, s_mir[ii], bmat1, lh);
                else
                  proxy1->Evaluator()->CalcMatrix(primary_trial_fel, primary_mir[ii], bmat1, lh);
                if(proxy2->IsOther())
                  proxy2->Evaluator()->CalcMatrix(secondary_test_fel, s_mir[ii], bmat2, lh);
                else
                  proxy2->Evaluator()->CalcMatrix(primary_test_fel, primary_mir[ii], bmat2, lh);
              }

              bdbmat1.Rows(r3) = proxyvalues(ii,STAR,STAR) * bmat1;
              bbmat2.Rows(r3) = bmat2;
            }
          loc_elmat += Trans(bbmat2) * bdbmat1;
        }
  }

  ContactBoundary::ContactBoundary(Region _master, Region _other,
                                   bool _draw_pairs, bool _volume,
                                   bool _element_boundary)
    : master(_master), other(_other), draw_pairs(_draw_pairs), volume(_volume),
      element_boundary(_element_boundary)
  {
#if NETGEN_USE_GUI
    if(draw_pairs)
      AddUserVisualizationObject (this);
#endif // NETGEN_USE_GUI

    auto mesh = master.Mesh();
    if (!mesh) return;
    
    if(mesh->GetDimension() == 2)
      {
        gap = make_shared<T_GapFunction<2>>(mesh, master, other);
        normal = make_shared<DisplacedNormal<2>>();
      }
    else
      {
        gap = make_shared<T_GapFunction<3>>(mesh, master, other);
        normal = make_shared<DisplacedNormal<3>>();
      }
  }

  ContactBoundary :: ~ContactBoundary()
  {
#if NETGEN_USE_GUI
    DeleteUserVisualizationObject (this);
#endif // NETGEN_USE_GUI
  }

  void ContactBoundary :: Draw()
  {
    if(!draw_pairs)
      return;

#if NETGEN_USE_GUI
    Array<netgen::Point<3>> points(2*primary_points.Size());
    for (auto i : Range(primary_points.Size()))
      {
        auto & mp = primary_points[i];
        auto & sp = secondary_points[i];
        points[2*i] = netgen::Point<3>(mp[0], mp[1], mp[2]);
        points[2*i+1] = netgen::Point<3>(sp[0], sp[1], sp[2]);
      }
    MyOpenGLLines(points);
#endif // NETGEN_USE_GUI
  }

  void ContactBoundary::AddEnergy(shared_ptr<CoefficientFunction> form,
                                  bool deformed)
  {
    energies.Append(make_shared<ContactEnergy>(form, deformed));
    if (deformed)
      deformed_energies.Append (energies.Last());
    else
      undeformed_energies.Append (energies.Last());
  }

  void ContactBoundary::AddIntegrator(shared_ptr<CoefficientFunction> form,
                                      bool deformed)
  {
    integrators.Append(make_shared<ContactIntegrator>(form, deformed));
    if (deformed)
      deformed_integrators.Append (integrators.Last());
    else
      undeformed_integrators.Append (integrators.Last());
  }

  void ContactBoundary::
  Update(shared_ptr<GridFunction> displacement_,
         shared_ptr<BilinearForm> bf,
         int intorder, double h, bool both_sides, bool keep_pairs, shared_ptr<ContactBoundary> other_cb,
         string quadrature)
  {
    if (quadrature != "gauss" && quadrature != "clip")
      throw Exception("ContactBoundary: unknown quadrature '" + quadrature + "', use 'gauss' or 'clip'");
    bool clip = quadrature == "clip";
    if (clip && element_boundary)
      throw Exception("ContactBoundary: quadrature='clip' not available with element_boundary");
    if(!displacement_ && !bf)
      throw Exception("Either displacement or BilinearForm needed in ContactBoundary update!");
    if(!trial_fes)
      {
        if(bf)
          {
            trial_fes = bf->GetTrialSpace();
            test_fes = bf->GetTestSpace();
          }
        else
          {
            trial_fes = displacement_->GetFESpace();
            test_fes = displacement_->GetFESpace();
          }
      }
    if(bf && (bf->GetTrialSpace().get() != trial_fes.get()))
      throw Exception("BilinearForm on different space as given to ContactBoundary!");
    if(displacement_)
      fes_displacement = displacement_->GetFESpace();

    if(draw_pairs)
    {
      secondary_points.SetSize(0);
      primary_points.SetSize(0);
    }

    shared_ptr<GridFunction> displacement = nullptr;
    if(displacement_)
      {
        auto flags = displacement_->GetFlags();
        flags.SetFlag("novisual");
        displacement = CreateGridFunction(displacement_->GetFESpace(), "_cb_displacement", flags);
        displacement->Update();
        displacement->GetVector() = displacement_->GetVector();
      }
    if (other_cb)
      gap->ShareSearchTree(*other_cb->gap);
    else if (displacement)
      gap->Update(displacement, 10*displacement->GetFESpace()->GetOrder(), h,
                  both_sides);
    else
      gap->Update(nullptr, 10, h, both_sides);
      
    auto mesh = trial_fes->GetMeshAccess();
    if(mesh->GetDimension() == 2)
      static_pointer_cast<DisplacedNormal<2>>(normal)->Update(displacement);
    else
      static_pointer_cast<DisplacedNormal<3>>(normal)->Update(displacement);
    if(bf)
      {
        static Timer t1("Build contact pairs");
        RegionTimer regt1(t1);
        LocalHeap lh(1000000, "ContactBoundary-Update", true);
        contact_rules.SetSize(0);
        has_contact_rules = keep_pairs;
        bool keep_geo = clip && keep_pairs;
        const ClipGeometry * shared_geo = nullptr;
        if (clip && other_cb && other_cb.get() != this && other_cb->has_clip_geometry)
          shared_geo = other_cb->clip_geometry.Data();
        if (keep_geo)
          clip_geometry.SetSize(mesh->GetNE(master.VB()));
        else
          clip_geometry.SetSize0();
        has_clip_geometry = keep_geo;
        Iterate<2>
          ([&](auto i)
           {
             constexpr auto DIM = i.value+2;
             if(mesh->GetDimension() == DIM)
               {
                 // Delete special elements created by me
                 auto& specialels = bf->GetSpecialElements();
                 auto last = specialels.Size()-1;
                 for(auto i : Range(specialels))
                   {
                     auto index = last-i;
                     auto cont_el = dynamic_cast<MPContactElement<DIM>*>(specialels[index].get());
                     if(cont_el && cont_el->GetContactBoundary().get() == this)
                       bf->DeleteSpecialElement(index);
                   }

                 // with 'clip' the pairing is cheap: share the search tree, clip at the own order
                 if(other_cb && !clip)
                   {
                     if(!other_cb->has_contact_rules)
                       throw Exception("ContactBoundary passed to Update has no contact pairs (call its Update with keep_pairs=True first)");
                     for(const auto& rule : other_cb->contact_rules)
                       bf->AddSpecialElement(make_unique<MPContactElement<DIM>>
                                             (rule.primary_ei, rule.secondary_ei,
                                              CopyIR(rule.primary_ir), CopyIR(rule.secondary_ir),
                                              shared_from_this(), displacement.get()));
                     return;
                   }

                 auto tgap = static_pointer_cast<T_GapFunction<DIM>>(gap);
                 static mutex add_mutex;
                 static mutex add_draw_mutex;
                 auto& mask = master.Mask();
                 mesh->IterateElements
                   (master.VB(), lh,
                    [&] (Ngs_Element el, LocalHeap& lh)
                    {
		      constexpr auto DIM = i.value+2; // Declare DIM again (MSVC bug)
                      HeapReset hr(lh);
                      if(!mask.Test(el.GetIndex())) return;
                      auto& trafo = mesh->GetTrafo(el, lh);
                      Array<unique_ptr<MappedIntegrationRule<DIM-1, DIM>>> mirs;
                      Array<IntegrationPoint> clip1, clip2;
                      Array<int> clipnr;
                      if (clip)
                        {
                          ClipGeometry own;
                          const ClipGeometry * geo = shared_geo ? &shared_geo[el.Nr()] : &own;
                          if (!shared_geo)
                            {
                              ClipGeometry & target = keep_geo ? clip_geometry[el.Nr()] : own;
                              tgap->CreateClipGeometry(trafo, intorder, lh, both_sides, target);
                              geo = &target;
                            }
                          ClipRules(*geo, DIM, intorder, clip1, clip2, clipnr);
                        }
                      else if (element_boundary)
                        {
                          auto eltype = trafo.GetElementType();
                          Facet2ElementTrafo transform(eltype, BND);
                          auto nfacet = transform.GetNFacets();
                          for(int facnr : Range(nfacet))
                            {
                              auto etfacet = transform.FacetType(facnr);
                              IntegrationRule ir_facet(etfacet, intorder);
                              IntegrationRule & ir_facet_vol = transform(facnr, ir_facet, lh);
                              mirs.Append(make_unique<MappedIntegrationRule<DIM-1, DIM>>(ir_facet_vol, trafo, lh));
                            }
                        }
                      else
                        {
                          IntegrationRule ir(trafo.GetElementType(), intorder);
                          mirs.Append(make_unique<MappedIntegrationRule<DIM-1, DIM>>(ir, trafo, lh));
                        }
                      int npts = clip1.Size();
                      for(const auto& mir : mirs)
                        npts += mir->Size();
                      FlatArray<IntegrationPoint> this_ir(npts, lh);
                      FlatArray<IntegrationPoint> other_ir(npts, lh);
                      FlatArray<size_t> other_nr(npts, lh);

                      int cntpair = 0;
                      for (size_t i = 0; i < clip1.Size(); i++, cntpair++)
                        {
                          this_ir[cntpair] = clip1[i];
                          other_ir[cntpair] = clip2[i];
                          other_nr[cntpair] = clipnr[i];
                        }
                      for(const auto& mir : mirs)
                      for(auto& mip : *mir)
                        {
                          auto pair = tgap->CreateContactPair(mip, lh,
                                                              both_sides);
                          if (pair.has_value())
                            {
                              this_ir[cntpair] =  (*pair).primary_ip;
                              if(element_boundary)
                                this_ir[cntpair].SetFacetNr
                                  (mip.IP().FacetNr(),
                                   mip.IP().VB());
                              other_ir[cntpair] = (*pair).secondary_ip;
                              other_nr[cntpair] = (*pair).secondary_el.Nr();
                              cntpair++;

                              if(draw_pairs)
                                {
                                  HeapReset hr(lh);
                                  auto & t1_def = trafo.AddDeformation(displacement.get(), lh);
                                  Vec<3> p1 = 0;
                                  t1_def.CalcPoint(pair->primary_ip, p1);

                                  auto & t2 = mesh->GetTrafo(pair->secondary_el, lh);
                                  auto & t2_def = t2.AddDeformation(displacement.get(), lh);
                                  Vec<3> p2 = 0;
                                  t2_def.CalcPoint(pair->secondary_ip, p2);

                                  lock_guard<mutex> guard(add_draw_mutex);
                                  primary_points.Append(p1);
                                  secondary_points.Append(p2);
                                }
                            }
                        }
                      // cout << "other_nr = " << other_nr << endl;
                      
                      FlatArray<int> index(cntpair, lh);
                      for (int i : Range(cntpair))
                        index[i] = i;

                      QuickSort (index, [&] (auto i1, auto i2) { return other_nr[i1] < other_nr[i2]; });
                      
                      // cout << "index = " << index << endl;
                      
                      int first = 0;
                      while (first < cntpair)
                        {
                          int next = first;
                          while ((next < cntpair && other_nr[index[next]] == other_nr[index[first]]) && (!element_boundary || this_ir[index[next]].FacetNr() == this_ir[index[first]].FacetNr()))
                            next++;
                          // cout << "interval = [" << first << ", " << next << ")" << endl;
                          
                          IntegrationRule primary_ir;
                          IntegrationRule secondary_ir;

                          for (int i = first; i < next; i++)
                            {
                              primary_ir.AddIntegrationPoint (this_ir[index[i]]);
                              secondary_ir.AddIntegrationPoint (other_ir[index[i]]);
                            }


                          /*
                          cout << "surfmir =  " << trafo(primary_ir, lh) << endl;
                          auto & t2 = mesh->GetTrafo(ElementId(BND, other_nr[index[first]]), lh);
                          cout << "secondarymir =  " << t2(secondary_ir, lh) << endl;                              
                          // they are not matching ??? 
                          */

                          auto this_ei = trafo.GetElementId();
                          auto other_ei = ElementId(BND, other_nr[index[first]]);
                          
                          
                          if (!volume)
                            {
                              lock_guard<mutex> guard(add_mutex);
                              if(keep_pairs)
                                AddSpecialElements(bf, ContactRule{this_ei, other_ei, CopyIR(primary_ir), CopyIR(secondary_ir)}, keep_pairs);
                              bf->AddSpecialElement(make_unique<MPContactElement<DIM>>
                                                    (el, other_ei,
                                                     std::move(primary_ir), std::move(secondary_ir),
                                                     shared_from_this(), displacement.get()));
                            }
                          else
                            {
                              // LocalHeapMem<100000> lh("lhbfv");
                              // const ElementTransformation & trafo = ir.GetTransformation();
                              // auto ei = trafo.GetElementId();
                              // const MeshAccess & ma = *static_cast<const MeshAccess*> (trafo.GetMesh());

                              int facet = mesh->GetElFacets(el)[0];
                              ArrayMem<int,2> elnums;
                              mesh->GetFacetElements (facet, elnums);
                              if (elnums.Size() != 1)
                                throw Exception("surface element does not have exactly one vol element");

                              ElementId vei(VOL, elnums[0]);
                              int locfacnr = mesh->GetElFacets(vei).Pos(facet);
                                  
                              ElementTransformation & vol_trafo = mesh->GetTrafo (vei, lh);
                              // if (!vol_cf->DefinedOn(vol_trafo)) continue;
                                  
                              Facet2ElementTrafo f2el(vol_trafo.GetElementType(), mesh->GetElVertices(vei));
                              Array<int> surfvnums { mesh->GetElVertices(el) };
                              Facet2SurfaceElementTrafo f2sel(trafo.GetElementType(), surfvnums);
                              
                              auto & ir_ref = f2sel.Inverse(primary_ir, lh);
                              auto & ir_vol = f2el(locfacnr, ir_ref, lh);

                              IntegrationRule volir;
                              for (int i = 0; i < ir_vol.Size(); i++)
                                volir.AddIntegrationPoint (ir_vol[i]);

                              /*
                                // vol and primary surf do match (but not secondary)
                              cout << "volmir =  " << vol_trafo(volir, lh) << endl;
                              cout << "surfmir =  " << trafo(primary_ir, lh) << endl;
                              
                              auto & t2 = mesh->GetTrafo(ElementId(BND, other_nr[index[first]]), lh);
                              cout << "secondarymir =  " << t2(secondary_ir, lh) << endl;                              
                              */
                              
                              lock_guard<mutex> guard(add_mutex);
                              if(keep_pairs)
                                AddSpecialElements(bf, ContactRule{vei, ElementId(BND, other_nr[index[first]]), CopyIR(volir), CopyIR(secondary_ir)}, keep_pairs);
                              bf->AddSpecialElement(make_unique<MPContactElement<DIM>>
                                                    (vei, ElementId(BND, other_nr[index[first]]),
                                                     std::move(volir), std::move(secondary_ir),
                                                     shared_from_this(), displacement.get()));
                            }
                          
                          first = next;
                        }
                    });
               }
           });
      }
  }

  void ContactBoundary :: AddSpecialElements(shared_ptr<BilinearForm> bf, ContactRule rule, bool keep_pairs)
  {
      if(keep_pairs) {
          static mutex add_mutex;
          lock_guard<mutex> guard(add_mutex);
          contact_rules.Append(std::move(rule));
      }
  }

  // ****************************** Contact Special Element Group **********************************
  
  ContactIntegrator2 :: ContactIntegrator2(std::variant<Region,string> _primary,
                                           std::variant<Region,string> _secondary,
                                           shared_ptr<GridFunction> _deformation,
                                           bool _volume, bool element_boundary)
  // : cb(make_shared<ContactBoundary>(_primary, _secondary)), primary(_primary), secondary(_secondary)
    : cb{make_shared<ContactBoundary>(Region(), Region())},
      primary(_primary), secondary(_secondary),
      deformation(_deformation)
  {
  }

  void ContactIntegrator2 :: AddIntegrator(shared_ptr<CoefficientFunction> form)
  {
    cb -> AddIntegrator (form);
    integrators.Append(make_shared<ContactIntegrator>(form, true));
    trial_fes = integrators.Last()->GetTrialSpace();
    test_fes = integrators.Last()->GetTestSpace();

    cb->SetTrialFESpace(trial_fes);
    cb->SetTestFESpace(test_fes);
  }

  int ContactIntegrator2 :: GetNElements()
  {
    // cout << "ContactIntegrator: have " << elements.Size() << " elements" << endl;
    return elements.Size();
  }

  
  void ContactIntegrator2 :: Update()
  {
    // copied from ContactBoundary update

    // cout << "ContactIntegrator Update called" << endl;

    /*
      // TODO
    if(displacement_)
      fes_displacement = displacement_->GetFESpace();

    shared_ptr<GridFunction> displacement = nullptr;
    
    if(displacement_)
      {
        auto flags = displacement_->GetFlags();
        flags.SetFlag("novisual");
        displacement = CreateGridFunction(displacement_->GetFESpace(), "_cb_displacement", flags);
        displacement->Update();
        displacement->GetVector() = displacement_->GetVector();
      }
    */
    
    auto mesh = trial_fes->GetMeshAccess();

    Region primary_reg;
    if (auto p = get_if<Region>(&primary))
      primary_reg = *p;
    else
      primary_reg = Region (mesh, BND, get<string>(primary));

    Region secondary_reg;
    if (auto p = get_if<Region>(&secondary))
      secondary_reg = *p;
    else
      secondary_reg = Region (mesh, BND, get<string>(secondary));
    
    bool element_boundary = false;  // TODO
    bool volume = false;    // TODO
    
    LocalHeap lh(1000000, "ContactIntegrator-Update", true);
    
    Iterate<2>
      ([&](auto i)
      {
        constexpr auto DIM = i.value+2;
        if(mesh->GetDimension() == DIM)
          {
            elements.SetSize0();
            auto tgap = make_shared<T_GapFunction<DIM>>(mesh, primary_reg, secondary_reg);
            tgap -> Update (deformation,
                            10, // intorder
                            5,  //  h
                            true  /* bothsides */ );
            
            static mutex add_mutex;
            auto& mask = primary_reg.Mask();

            mesh->IterateElements
              (primary_reg.VB(), lh,
               [&] (Ngs_Element el, LocalHeap& lh)
               {
                 constexpr auto DIM = i.value+2; // Declare DIM again (MSVC bug)
                 HeapReset hr(lh);
                 if(!mask.Test(el.GetIndex())) return;

                 auto& trafo = mesh->GetTrafo(el, lh);
                 Array<unique_ptr<MappedIntegrationRule<DIM-1, DIM>>> mirs;
                 /*
                 // TODO
                 if (element_boundary)
                 {
                 auto eltype = trafo.GetElementType();
                 Facet2ElementTrafo transform(eltype, BND);
                 auto nfacet = transform.GetNFacets();
                 for(int facnr : Range(nfacet))
                 {
                 auto etfacet = transform.FacetType(facnr);
                 IntegrationRule ir_facet(etfacet, intorder);
                 IntegrationRule & ir_facet_vol = transform(facnr, ir_facet, lh);
                 mirs.Append(make_unique<MappedIntegrationRule<DIM-1, DIM>>(ir_facet_vol, trafo, lh));
                 }
                 }
                 else
                 */
                 {
                   IntegrationRule ir(trafo.GetElementType(), intorder);
                   mirs.Append(make_unique<MappedIntegrationRule<DIM-1, DIM>>(ir, trafo, lh));
                 }
                 
                 int npts = 0;
                 for(const auto& mir : mirs)
                   npts += mir->Size();
                 FlatArray<IntegrationPoint> this_ir(npts, lh);
                 FlatArray<IntegrationPoint> other_ir(npts, lh);
                 FlatArray<size_t> other_nr(npts, lh);
                 
                 int cntpair = 0;
                 for(const auto& mir : mirs)
                   for(auto& mip : *mir)
                     {
                       bool both_sides = true; // ???
                       auto pair = tgap->CreateContactPair(mip, lh,
                                                           both_sides);
                       if (pair.has_value())
                         {
                           this_ir[cntpair] =  (*pair).primary_ip;
                           /*
                             // todo
                           if(element_boundary)
                             this_ir[cntpair].SetFacetNr
                               (mip.IP().FacetNr(),
                                mip.IP().VB());
                           */
                           other_ir[cntpair] = (*pair).secondary_ip;
                           other_nr[cntpair] = (*pair).secondary_el.Nr();
                           cntpair++;
                         }
                     }
                 
                 FlatArray<int> index(cntpair, lh);
                 for (int i : Range(cntpair))
                   index[i] = i;
                 
                 QuickSort (index, [&] (auto i1, auto i2) { return other_nr[i1] < other_nr[i2]; });
                 
                 // cout << "index = " << index << endl;
                 
                 int first = 0;
                 while (first < cntpair)
                   {
                     int next = first;
                     while ((next < cntpair && other_nr[index[next]] == other_nr[index[first]]) &&
                            (!element_boundary || this_ir[index[next]].FacetNr() == this_ir[index[first]].FacetNr()))
                       next++;
                     // cout << "interval = [" << first << ", " << next << ")" << endl;
                     
                     IntegrationRule primary_ir;
                     IntegrationRule secondary_ir;
                     
                     for (int i = first; i < next; i++)
                       {
                         primary_ir.AddIntegrationPoint (this_ir[index[i]]);
                         secondary_ir.AddIntegrationPoint (other_ir[index[i]]);
                       }
                     
                     
                     /*
                       cout << "surfmir =  " << trafo(primary_ir, lh) << endl;
                       auto & t2 = mesh->GetTrafo(ElementId(BND, other_nr[index[first]]), lh);
                       cout << "secondarymir =  " << t2(secondary_ir, lh) << endl;                              
                       // they are not matching ??? 
                       */
                     
                     
                     if (!volume)
                       {
                         lock_guard<mutex> guard(add_mutex);
                         elements.Append (make_unique<MPContactElement<DIM>>
                                          (el, ElementId(BND, other_nr[index[first]]),
                                           std::move(primary_ir), std::move(secondary_ir),
                                           cb, 
                                           nullptr // displacement.get()
                                           ));
                       }
                 /*
                   TODO
                   else
                   {
                   // LocalHeapMem<100000> lh("lhbfv");
                   // const ElementTransformation & trafo = ir.GetTransformation();
                   // auto ei = trafo.GetElementId();
                   // const MeshAccess & ma = *static_cast<const MeshAccess*> (trafo.GetMesh());
                   
                   int facet = mesh->GetElFacets(el)[0];
                   ArrayMem<int,2> elnums;
                   mesh->GetFacetElements (facet, elnums);
                   if (elnums.Size() != 1)
                   throw Exception("surface element does not have exactly one vol element");
                   
                   ElementId vei(VOL, elnums[0]);
                   int locfacnr = mesh->GetElFacets(vei).Pos(facet);
                   
                   ElementTransformation & vol_trafo = mesh->GetTrafo (vei, lh);
                   // if (!vol_cf->DefinedOn(vol_trafo)) continue;
                   
                   Facet2ElementTrafo f2el(vol_trafo.GetElementType(), mesh->GetElVertices(vei));
                   Array<int> surfvnums { mesh->GetElVertices(el) };
                   Facet2SurfaceElementTrafo f2sel(trafo.GetElementType(), surfvnums);
                   
                   auto & ir_ref = f2sel.Inverse(primary_ir, lh);
                   auto & ir_vol = f2el(locfacnr, ir_ref, lh);
                   
                   IntegrationRule volir;
                   for (int i = 0; i < ir_vol.Size(); i++)
                   volir.AddIntegrationPoint (ir_vol[i]);
                   
                   
                   lock_guard<mutex> guard(add_mutex);
                   bf->AddSpecialElement(make_unique<MPContactElement<DIM>>
                   (vei, ElementId(BND, other_nr[index[first]]),
                   std::move(volir), std::move(secondary_ir),
                   shared_from_this(), displacement.get()));
                   }
                 */
                 
                     first = next;
                   }   // while(...)
               }); // iterate elements
          }
      }); // iterate dimension
  }

  void ContactIntegrator2 :: GetDofNrs(std::function<void(int,FlatArray<DofId>)> eldofs)
  {
    Array<DofId> dnums;
    for (int i = 0; i < elements.Size(); i++)
      {
        elements[i]->GetDofNrs(dnums);
        eldofs(i, dnums);
      }
  }

  void ContactIntegrator2 :: Assemble(std::function<void(FlatArray<DofId>,FlatArray<DofId>,FlatMatrix<double>,ElementId,LocalHeap&)> addelmat, LocalHeap& lh)
  {
    Array<DofId> dnums;
    for (auto& el : elements)
      {
        HeapReset hr(lh);
        el->GetDofNrs(dnums);
        FlatMatrix<double> elmat(dnums.Size(), dnums.Size(), lh);
        el->CalcElementMatrix(elmat, lh);

        addelmat(dnums, dnums, elmat, 0, lh);
      }
  }



  // ****************************************************************************************************
  
  void ContactBoundary :: DoArchive(Archive& ar)
  {
    ar & primary_points & secondary_points;
  }

  template<int DIM>
  MPContactElement<DIM>::MPContactElement(ElementId aprimary_ei, ElementId asecondary_ei,
                                          IntegrationRule aprimary_ir, IntegrationRule asecondary_ir,
                                          shared_ptr<ContactBoundary> _cb,
                                          GridFunction* _deformation)
    : primary_ei(aprimary_ei), secondary_ei(asecondary_ei),
      primary_ir(std::move(aprimary_ir)), secondary_ir(std::move(asecondary_ir)),
      cb(_cb),
      trial_fes(_cb->GetTrialFESpace().get()),
      test_fes(_cb->GetTestFESpace().get()),
      deformation(_deformation)
  { }

  template<int DIM>
  void MPContactElement<DIM>::GetDofNrs(Array<DofId>& dnums) const
  {
    trial_fes->GetDofNrs(primary_ei, dnums);
    Array<DofId> s_dofs;
    trial_fes->GetDofNrs(secondary_ei, s_dofs);
    dnums += s_dofs;
  }

  template<int DIM>
  void MPContactElement<DIM>::GetDofNrs2(Array<DofId>& dnums) const
  {
    test_fes->GetDofNrs(primary_ei, dnums);
    Array<DofId> s_dofs;
    test_fes->GetDofNrs(secondary_ei, s_dofs);
    dnums += s_dofs;
  }

  
  template<int DIM>
  double MPContactElement<DIM>::Energy(FlatVector<double> elx,
                                     LocalHeap& lh) const
  {
    if(cb->GetIntegrators().Size())
      throw Exception("Energy does not work with contact integrators!");
    
    auto& primary_trafo = trial_fes->GetMeshAccess()->GetTrafo(primary_ei, lh);
    auto& secondary_trafo = trial_fes->GetMeshAccess()->GetTrafo(secondary_ei, lh);
    auto& primary_deformed_trafo = primary_trafo.AddDeformation(deformation, lh);
    auto& secondary_deformed_trafo = secondary_trafo.AddDeformation(deformation, lh);

    // IntegrationRule primary_ir(1, &const_cast<IntegrationPoint&>(pair.primary_ip));
    // IntegrationRule secondary_ir(1, &const_cast<IntegrationPoint&>(pair.secondary_ip));

    auto& primary_fel = trial_fes->GetFE(primary_ei, lh);
    auto& secondary_fel = trial_fes->GetFE(secondary_ei, lh);
    
    double energy = 0.;

    for (bool def : { false, true })
      {
        auto& energies = cb->GetEnergies(def);
        if (energies.Size())
          {
            MappedIntegrationRule<DIM-1,DIM> primary_mir(primary_ir, def ? primary_deformed_trafo :  primary_trafo, lh);
            MappedIntegrationRule<DIM-1,DIM> secondary_mir(secondary_ir, def ? secondary_deformed_trafo : secondary_trafo, lh);
            primary_mir.SetOtherMIR(&secondary_mir);
            for(const auto& ce : energies) 
              energy += ce->CalcEnergy(primary_fel, secondary_fel, primary_mir, elx, lh);
          }
      }

    return energy;
  }

  template<int DIM>
  void MPContactElement<DIM>::Apply(FlatVector<double> elx,
                                  FlatVector<double> ely,
                                  LocalHeap& lh) const
  {
    // NETGEN_TIMER_FROM_HERE("ContactElement::Apply");    
    ely = 0.;

    auto& primary_trafo = trial_fes->GetMeshAccess()->GetTrafo(primary_ei, lh);
    auto& secondary_trafo = trial_fes->GetMeshAccess()->GetTrafo(secondary_ei, lh);
    auto& primary_deformed_trafo = primary_trafo.AddDeformation(deformation, lh);
    auto& secondary_deformed_trafo = secondary_trafo.AddDeformation(deformation, lh);

    // IntegrationRule primary_ir(1, &const_cast<IntegrationPoint&>(pair.primary_ip));
    // IntegrationRule secondary_ir(1, &const_cast<IntegrationPoint&>(pair.secondary_ip));

    auto& primary_trial_fel = trial_fes->GetFE(primary_ei, lh);
    auto& secondary_trial_fel = trial_fes->GetFE(secondary_ei, lh);
    auto& primary_test_fel = test_fes->GetFE(primary_ei, lh);
    auto& secondary_test_fel = test_fes->GetFE(secondary_ei, lh);
    
    if(primary_ei.IsBoundary())
      {
    for (bool def : { false, true })
      {
        auto& energies = cb->GetEnergies(def);
        auto& integrators = cb->GetIntegrators(def);
        if (energies.Size() || integrators.Size())
          {
            MappedIntegrationRule<DIM-1,DIM> primary_mir(primary_ir, def ? primary_deformed_trafo :  primary_trafo, lh);
            MappedIntegrationRule<DIM-1,DIM> secondary_mir(secondary_ir, def ? secondary_deformed_trafo : secondary_trafo, lh);
            primary_mir.SetOtherMIR(&secondary_mir);

            for(const auto& ci : integrators)
              ci->ApplyAdd(primary_trial_fel, secondary_trial_fel, primary_test_fel, secondary_test_fel, primary_mir, elx, ely, lh);
            for(const auto& ce : energies)
              ce->ApplyAdd(primary_trial_fel, secondary_trial_fel, primary_mir, elx, ely, lh);
          }
      }
      }
    else
      {
        for(bool def : { false, true  })
          {
            auto& energies = cb->GetEnergies(def);
            auto& integrators = cb->GetIntegrators(def);
            if (energies.Size() || integrators.Size())
              {
                MappedIntegrationRule<DIM,DIM> primary_mir(primary_ir, def ? primary_deformed_trafo :  primary_trafo, lh);
                primary_mir.ComputeNormalsAndMeasure (primary_trafo.GetElementType(), primary_ir[0].FacetNr());
                MappedIntegrationRule<DIM-1,DIM> secondary_mir(secondary_ir, def ? secondary_deformed_trafo : secondary_trafo, lh);
                primary_mir.SetOtherMIR(&secondary_mir);
                for(const auto& ci : integrators)
                  ci->ApplyAdd(primary_trial_fel, secondary_trial_fel, primary_test_fel, secondary_test_fel, primary_mir, elx, ely, lh);
                for(const auto& ce : energies)
                  ce->ApplyAdd(primary_trial_fel, secondary_trial_fel, primary_mir, elx, ely, lh);

                }
          }
      }
  }


  template<int DIM>
  void MPContactElement<DIM>::CalcElementMatrix(FlatMatrix<double> elmat,
                                                LocalHeap& lh) const
  {
    HeapReset hr(lh);
    FlatVector<> elx(elmat.Width(), lh);
    elx = 0;
    CalcLinearizedElementMatrix (elx, elmat, lh);
  }

  
  template<int DIM>
  void MPContactElement<DIM>::CalcLinearizedElementMatrix(FlatVector<double> elx,
                                                        FlatMatrix<double> elmat,
                                                        LocalHeap& lh) const
  {
    // NETGEN_TIMER_FROM_HERE("MPContactElement::CalcLinearizedElementMatrix");
    
    elmat = 0.;

    auto& primary_trafo = trial_fes->GetMeshAccess()->GetTrafo(primary_ei, lh);
    auto& secondary_trafo = trial_fes->GetMeshAccess()->GetTrafo(secondary_ei, lh);
    auto& primary_deformed_trafo = primary_trafo.AddDeformation(deformation, lh);
    auto& secondary_deformed_trafo = secondary_trafo.AddDeformation(deformation, lh);

    // IntegrationRule primary_ir(1, &const_cast<IntegrationPoint&>(primary_e));
    // IntegrationRule secondary_ir(1, &const_cast<IntegrationPoint&>(pair.secondary_ip));

    auto& primary_trial_fel = trial_fes->GetFE(primary_ei, lh);
    auto& secondary_trial_fel = trial_fes->GetFE(secondary_ei, lh);
    auto& primary_test_fel = test_fes->GetFE(primary_ei, lh);
    auto& secondary_test_fel = test_fes->GetFE(secondary_ei, lh);

    if (primary_ei.IsBoundary())
      for (bool def : { false, true })
        {
          auto& energies = cb->GetEnergies(def);
          auto& integrators = cb->GetIntegrators(def);
          if (energies.Size() || integrators.Size())
            {
              MappedIntegrationRule<DIM-1,DIM> primary_mir(primary_ir, def ? primary_deformed_trafo :  primary_trafo, lh);
              MappedIntegrationRule<DIM-1,DIM> secondary_mir(secondary_ir, def ? secondary_deformed_trafo : secondary_trafo, lh);
              primary_mir.SetOtherMIR(&secondary_mir);
              
              for(const auto& ci : integrators)
                ci->CalcLinearizedAdd(primary_trial_fel, secondary_trial_fel, primary_test_fel, secondary_test_fel, primary_mir, elx, elmat, lh);            
              for(const auto& ce : energies)
                ce->CalcLinearizedAdd(primary_trial_fel, secondary_trial_fel, primary_mir, elx, elmat, lh);            
            }
        }
    else
      {
        // cout << "volume contact integrator" << endl;
      for (bool def : { false, true })
        {
          auto& energies = cb->GetEnergies(def);
          auto& integrators = cb->GetIntegrators(def);
          if (energies.Size() || integrators.Size())
            {
              MappedIntegrationRule<DIM,DIM> primary_mir(primary_ir, def ? primary_deformed_trafo :  primary_trafo, lh);
              primary_mir.ComputeNormalsAndMeasure (primary_trafo.GetElementType(), primary_ir[0].FacetNr());
              
              MappedIntegrationRule<DIM-1,DIM> secondary_mir(secondary_ir, def ? secondary_deformed_trafo : secondary_trafo, lh);

              /*
              cout << "primary mir: " << endl;
              for (const auto & mip : primary_mir)
                cout << mip << endl;
              cout << "secondary mir: " << endl;
              for (const auto & mip : secondary_mir)
                cout << mip << endl;
              */
              primary_mir.SetOtherMIR(&secondary_mir);
              
              for(const auto& ci : integrators)
                ci->CalcLinearizedAdd(primary_trial_fel, secondary_trial_fel, primary_test_fel, secondary_test_fel, primary_mir, elx, elmat, lh);            
              for(const auto& ce : energies)
                ce->CalcLinearizedAdd(primary_trial_fel, secondary_trial_fel, primary_mir, elx, elmat, lh);            
            }
        }
      }
      
  }

  template class MPContactElement<2>;
  template class MPContactElement<3>;

  RegisterClassForArchive<ContactBoundary> regcb;
} // namespace ngcomp
