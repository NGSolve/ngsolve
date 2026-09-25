#include "potentialcf.hpp"

#include "../comp/meshaccess.hpp"
#include "../comp/gridfunction.hpp"

#include "analytic_integrals.hpp"
#include "kernels.hpp"


namespace ngsbem
{
  template <typename TSCAL>
  PotentialCF<TSCAL> ::
  PotentialCF (shared_ptr<GridFunction> _gf,
               VorB _source_vb,
               optional<Region> _definedon,
               shared_ptr<DifferentialOperator> _evaluator,
               shared_ptr<const BaseIntegralKernel<TSCAL>> _kernel, int _intorder,
               IntOp_Parameters _io_params)
    : BasePotentialCF(_gf, _source_vb, _definedon, _evaluator, std::is_same<TSCAL,Complex>()),
      kernel(std::move(_kernel)), intorder(_intorder)
  {
    io_params = _io_params;
    IVec<2> shape = kernel->Shape();
    if (shape[0] > 1)
      this->SetDimensions( Array<int>( { shape[0] } ));
  }


  template <typename TSCAL>
  shared_ptr<CoefficientFunction> PotentialCF<TSCAL> ::
  Operator (const string & name) const
  {
    return make_shared<PotentialCF<TSCAL>>(gf, source_vb, definedon, evaluator,
                                          kernel->GetDifferentiatedKernel(name), intorder, io_params);
  }


  // x is near a source element if L2Norm(x-c) < r
  template <int DIMS>
  tuple<Vec<3>,double> NearElementBall (const ElementTransformation & trafo, const IntegrationPoint & center)
  {
    MappedIntegrationPoint<DIMS,3> mip(center, trafo);
    return { mip.GetPoint(), L2Norm(mip.GetJacobian()) };
  }

  optional<tuple<Vec<3>,double>> PotentialNearfieldBall (const ElementTransformation & trafo)
  {
    switch (trafo.GetElementType())
      {
      case ET_SEGM: return NearElementBall<1> (trafo, IntegrationPoint(1./2, 0, 0));
      case ET_TRIG: return NearElementBall<2> (trafo, IntegrationPoint(1./3, 1./3));
      case ET_QUAD: return NearElementBall<2> (trafo, IntegrationPoint(1./2, 1./2));
      case ET_TET:  return NearElementBall<3> (trafo, IntegrationPoint(1./4, 1./4, 1./4));
      default:      return nullopt;
      }
  }

  // the source elements with a near ball, as boxes slightly larger than the balls: the ones containing x include all
  // elements x is near
  struct PotentialNearSources
  {
    Array<size_t> elnr;
    unique_ptr<netgen::BoxTree<3,int>> tree;
  };


  template <typename TSCAL>
  void PotentialCF<TSCAL> ::
  BuildLocalExpansion(const Region & reg)
  {
    LocalHeapMem<100000> lh("PotentialCF::BuildLocalExpansion");

    auto space = this->gf->GetFESpace();
    auto mesh = space->GetMeshAccess();


    Vec<3> smax(-1e99, -1e99, -1e99);
    Vec<3> smin(1e99, 1e99, 1e99);
    auto near = make_shared<PotentialNearSources>();
    Array<netgen::Box<3>> near_boxes;

    for (size_t i = 0; i < mesh->GetNE(source_vb); i++)
      {
        HeapReset hr(lh);
        ElementId ei(source_vb, i);
        if (!space->DefinedOn(ei)) continue;
        if (definedon && !(*definedon).Mask().Test(mesh->GetElIndex(ei))) continue;

        const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);
        IntegrationRule ir(trafo.GetElementType(), intorder);
        auto & miry = trafo(ir, lh);

        for (int k = 0; k < miry.Size(); k++)
          for (int j = 0; j < 3; j++)
            {
              smin(j) = min(smin(j), miry[k].GetPoint()(j));
              smax(j) = max(smax(j), miry[k].GetPoint()(j));
            }

        if (auto ball = PotentialNearfieldBall (trafo))
          {
            auto [c, r] = *ball;
            netgen::Box<3> box(netgen::Point<3>(c(0), c(1), c(2)));
            box.Increase ((1+1e-8)*r);
            near_boxes.Append (box);
            near->elnr.Append (i);
          }
      }

    if (near_boxes.Size())
      {
        netgen::Box<3> all(netgen::Box<3>::EMPTY_BOX);
        for (auto & box : near_boxes)
          {
            all.Add (box.PMin());
            all.Add (box.PMax());
          }
        near->tree = make_unique<netgen::BoxTree<3,int>> (all);
        for (size_t k = 0; k < near_boxes.Size(); k++)
          near->tree->Insert (near_boxes[k], int(k));
      }
    near_sources = near;

    Vec<3> cs = 0.5*(smin+smax);
    double rs = MaxNorm(smax-smin);

    auto & source = kernel->Source();
    bool source_needs_normal = source.NeedsNormal();
    auto singmp = source.CreateMultipoleExpansion(cs, rs, io_params);

    typedef TSCAL T;
    for (size_t i = 0; i < mesh->GetNE(source_vb); i++)
      {
        HeapReset hr(lh);
        ElementId ei(source_vb, i);

        if (!space->DefinedOn(ei)) continue;
        if (definedon && !(*definedon).Mask().Test(mesh->GetElIndex(ei))) continue;

        const FiniteElement &fel = space->GetFE(ei, lh);
        if (fel.GetNDof() == 0) continue;   // no dofs on this element (e.g. L2 on curves): nothing to integrate
        const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);

        Array<DofId> dnums(fel.GetNDof(), lh);
        space->GetDofNrs(ei, dnums);
        FlatVector<T> elvec(fel.GetNDof(), lh);
        gf->GetElementVector(dnums, elvec);

        IntegrationRule ir(fel.ElementType(), intorder);
        auto & miry = trafo(ir, lh);
        FlatMatrix<T> vals(miry.Size(), evaluator->Dim(), lh);

        evaluator->Apply (fel, miry, elvec, vals, lh);

        // add vals to multipole ...
        for (int j = 0; j < miry.Size(); j++)
          {
            vals.Row(j) *= miry[j].GetWeight();
            Vec<3> ny = 0.0;
            if (source_needs_normal)
              {
                if (source_vb != BND)
                  throw Exception("kernel requires boundary source normals");
                ny = static_cast<const MappedIntegrationPoint<2,3>&>(miry[j]).GetNV();
              }
            source.AddSource (*singmp, miry[j].GetPoint(), ny, make_BareSliceVector(vals.Row(j)));
          }
      }

    singmp->CalcMP();

    Vec<3> tmax(-1e99, -1e99, -1e99);
    Vec<3> tmin(1e99, 1e99, 1e99);
    Array<Vec<3>> tpoints;
    auto tmesh = reg.Mesh();
    for (auto el : reg.GetElements())
      {
        HeapReset hr(lh);

        const ElementTransformation &trafo = tmesh->GetTrafo(el, lh);
        IntegrationRule ir(trafo.GetElementType(), intorder);
        auto & miry = trafo(ir, lh);

        for (int k = 0; k < miry.Size(); k++)
          {
            tpoints.Append (Vec<3>(miry[k].GetPoint()));
            for (int j = 0; j < 3; j++)
              {
                tmin(j) = min(tmin(j), miry[k].GetPoint()(j));
                tmax(j) = max(tmax(j), miry[k].GetPoint()(j));
              }
          }
      }

    Vec<3> ct = 0.5*(tmin+tmax);
    double rt = MaxNorm(tmax-tmin);


    double l2 = ceil (log2 (rt/rs));
    rt = exp2 (l2) * rs;

    local_expansion = kernel->Target().CreateLocalExpansion(ct, rt, io_params);

    for (auto el : reg.GetElements())
      {
        HeapReset hr(lh);

        const ElementTransformation &trafo = tmesh->GetTrafo(el, lh);
        IntegrationRule ir(trafo.GetElementType(), intorder);
        auto & miry = trafo(ir, lh);


        Vec<3> elmax(-1e99, -1e99, -1e99);
        Vec<3> elmin(1e99, 1e99, 1e99);

        for (int k = 0; k < miry.Size(); k++)
          {
            for (int j = 0; j < 3; j++)
              {
                elmin(j) = min(elmin(j), miry[k].GetPoint()(j));
                elmax(j) = max(elmax(j), miry[k].GetPoint()(j));
              }
          }

        Vec<3> el_center = 0.5 * (elmin+elmax);
        double el_rad = 0.5 * L2Norm(elmax-elmin);
        local_expansion -> AddVolumeTarget (el_center, el_rad);
      }

    local_expansion->CalcMP(singmp, true);
  }


  // minimize the function f(x) = 1/2 x a x + b x + c
  // returns (x, f(x))
  tuple<double,double> MinimizeOnSegm (double a, double b, double c)
  {
    if (a > 0)
      {
        double x = -b/a;
        if (x > 0 && x < 1)
          return { x, 0.5*a*x*x + b*x + c };
      }
    double val0 = c;
    double val1 = 0.5*a+b+c;
    if (val0 < val1)
      return { 0, val0 };
    else
      return { 1, val1 };
  }

  // 1/2 x^T A x + b x + c
  Vec<2> MinimizeOnTrig (Mat<2,2> a, Vec<2> b, double c)
  {
    if (a(0,0) > 0 && Det(a) > 0)
      {
        Vec<2> sol = -Inv(a)*b;
        if (sol(0) > 0 && sol(1) > 0 && sol(0)+sol(1) < 1)
          return sol;
      }

    auto [x0,val0] = MinimizeOnSegm(a(0,0), b(0), c);
    auto [x1,val1] = MinimizeOnSegm(a(1,1), b(1), c);
    Vec<2> p{1,0}, d{-1,1};
    Vec<2> g = a*p + b;
    auto [x2,val2] = MinimizeOnSegm(InnerProduct(a*d,d), InnerProduct(g,d), 0.5*a(0,0)+b(0)+c);  // (1,0) + (-1,1)*s

    if (val0 < val1 && val0 < val2)
      return { x0, 0 };
    if (val1 < val2)
      return { 0, x1 };
    return p+x2*d;
  }


  Vec<2> MinimizeOnQuad (Mat<2,2> a, Vec<2> b, double c)
  {
    if (a(0,0) > 0 && Det(a) > 0)
      {
        Vec<2> sol = -Inv(a)*b;
        if (sol(0) > 0 && sol(0) < 1 && sol(1) > 0 && sol(1) < 1)
          return sol;
      }

    auto [y0,val0] = MinimizeOnSegm(a(1,1), b(1), c);
    auto [y1,val1] = MinimizeOnSegm(a(1,1), a(0,1)+b(1), 0.5*a(0,0)+b(0)+c);
    auto [x0,val2] = MinimizeOnSegm(a(0,0), b(0), c);
    auto [x1,val3] = MinimizeOnSegm(a(0,0), a(0,1)+b(0), 0.5*a(1,1)+b(1)+c);

    if (val0 <= val1 && val0 <= val2 && val0 <= val3)
      return { 0, y0 };
    if (val1 <= val2 && val1 <= val3)
      return { 1, y1 };
    if (val2 <= val3)
      return { x0, 0 };
    return { x1, 1 };
  }


  // IntegrationPoint ProjectPointToReference(Vec<3> x, const ElementTransformation & trafo)
  // {
  //   auto et = trafo.GetElementType();
  //   IntegrationPoint ip = et == ET_TRIG ?
  //     IntegrationPoint(1./3, 1./3) : IntegrationPoint(1./2, 1./2);
  //   constexpr double reference_step_tolerance = 1e-12;
  //   for (int j = 0; j < 5; j++) // SQP steps
  //     {
  //       MappedIntegrationPoint<2,3> mip(ip, trafo);
  //       Mat<3,2> jac = mip.GetJacobian();
  //       Vec<2> ipvec { ip(0), ip(1) };
  //       auto Hesse = mip.CalcHesse();
  //       Vec<3> r = mip.GetPoint()-x;

  //       // Newton model of 1/2 ||F(uv)-x||^2, with a*ipvec+b = J^T*r.
  //       Mat<2,2> a = Trans(jac)*jac;
  //       for (int k = 0; k < 3; k++)
  //         a += r(k)*Hesse[k];
  //       // CalcHesse uses finite differences; enforce symmetry for the minimizer.
  //       a(0,1) = a(1,0) = 0.5*(a(0,1)+a(1,0));
  //       Vec<2> b = Trans(jac)*r-a*ipvec;
  //       Vec<2> uv = et == ET_TRIG ? MinimizeOnTrig(a, b, 0) : MinimizeOnQuad(a, b, 0);
  //       ip = IntegrationPoint(uv(0), uv(1));
  //       if (L2Norm(uv-ipvec) <= reference_step_tolerance)
  //         break;
  //     }
  //   return ip;
  // }


  IntegrationPoint ProjectPointToReference(Vec<3> x, const ElementTransformation & trafo)
  {
    auto et = trafo.GetElementType();
    IntegrationPoint ip = et == ET_TRIG ?
      IntegrationPoint(1./3, 1./3) : IntegrationPoint(1./2, 1./2);
    constexpr int max_iterations = 5;
    constexpr double reference_step_tolerance = 1e-12;
    for (int j = 0; j < max_iterations; j++) // Gauss-newton steps
      {
        MappedIntegrationPoint<2,3> mip(ip, trafo);
        Mat<3,2> jac = mip.GetJacobian();
        Vec<2> ipvec { ip(0), ip(1) };

        Mat<2,2> a = Trans(jac)*jac;
        Vec<2> b = -Trans(jac) * (x-mip.GetPoint()+jac*ipvec);
        Vec<2> uv = et == ET_TRIG ? MinimizeOnTrig(a, b, 0) : MinimizeOnQuad(a, b, 0);

        double reference_step = L2Norm(uv-ipvec);
        ip = IntegrationPoint(uv(0), uv(1));
        if (reference_step <= reference_step_tolerance)
          break;
      }
    return ip;
  }


  IntegrationPoint ProjectPointToReferenceTet(Vec<3> x, const ElementTransformation & trafo)
  {
    Vec<3> xi(1./4, 1./4, 1./4);
    constexpr int max_iterations = 10;
    constexpr double reference_step_tolerance = 1e-12;
    for (int j = 0; j < max_iterations; j++) // Newton steps for the inverse map
      {
        IntegrationPoint ip(xi(0), xi(1), xi(2));
        MappedIntegrationPoint<3,3> mip(ip, trafo);
        if (fabs(Det(mip.GetJacobian())) < 1e-14)  // degenerate/inverted element
          break;
        Vec<3> dx = Inv(mip.GetJacobian()) * (x-mip.GetPoint());
        xi += dx;
        if (L2Norm(dx) <= reference_step_tolerance)
          break;
      }
    if (!(isfinite(xi(0)) && isfinite(xi(1)) && isfinite(xi(2))))
      xi = Vec<3>(1./4, 1./4, 1./4);

    // Clamp into the reference tet via barycentric coordinates. For a point
    // outside the element this is only an approximate closest point.
    double lam[4] = { xi(0), xi(1), xi(2), 1-xi(0)-xi(1)-xi(2) };
    bool outside = false;
    for (int k = 0; k < 4; k++)
      if (lam[k] < 0)
        {
          lam[k] = 0;
          outside = true;
        }
    if (outside)
      {
        double sum = lam[0]+lam[1]+lam[2]+lam[3];
        for (int k = 0; k < 4; k++)
          lam[k] /= sum;
      }
    return IntegrationPoint(lam[0], lam[1], lam[2]);
  }


  IntegrationRule GetIntegrationRule(Vec<3> x, const ElementTransformation & trafo, int intorder, LocalHeap & lh)
  {
    auto et = trafo.GetElementType();
    if (et == ET_TET)
      {
        IntegrationPoint ipc(1./4, 1./4, 1./4);
        MappedIntegrationPoint<3,3> mipc(ipc, trafo);
        double elsize = L2Norm(mipc.GetJacobian());
        double dist = L2Norm(x-mipc.GetPoint());

        if (dist < elsize)
          {
            // Duffy-type rule: split the reference tet into 4 sub-tets meeting at
            // the projection of x, and map a prism rule (trig x segm) onto each
            // sub-tet, collapsing the segment end t=1 onto the apex (Jacobian (1-t)^2).
            IntegrationPoint p = ProjectPointToReferenceTet(x, trafo);
            Vec<3> vp(p(0), p(1), p(2));

            int order = intorder + 2;
            IntegrationRule irtrig(ET_TRIG, order), irsegm(ET_SEGM, order);
            IntegrationRule ir(4*irtrig.Size()*irsegm.Size(), lh);
            size_t cnt = 0;

            auto verts = ElementTopology::GetVertices(ET_TET);
            auto faces = ElementTopology::GetFaces(ET_TET);
            for (int f = 0; f < faces.Size(); f++)
              {
                Vec<3> v0(verts[faces[f][0]][0], verts[faces[f][0]][1], verts[faces[f][0]][2]);
                Vec<3> v1(verts[faces[f][1]][0], verts[faces[f][1]][1], verts[faces[f][1]][2]);
                Vec<3> v2(verts[faces[f][2]][0], verts[faces[f][2]][1], verts[faces[f][2]][2]);
                Mat<3,3> sides;
                sides.Col(0) = v1-v0;
                sides.Col(1) = v2-v0;
                sides.Col(2) = vp-v0;
                double factor = Det(sides);
                if (!(fabs(factor) > 1e-12))
                  continue;

                // Weights of one sub-tet sum to factor/6 = its volume
                // (trig weights sum to 1/2, int_0^1 (1-t)^2 dt = 1/3).
                for (auto ips : irtrig)
                  for (auto ipt : irsegm)
                    {
                      Vec<3> F = v0 + ips(0)*(v1-v0) + ips(1)*(v2-v0);
                      double t = ipt(0);
                      Vec<3> y = F + t*(vp-F);
                      ir[cnt++] = IntegrationPoint(y(0), y(1), y(2),
                                                   ips.Weight()*ipt.Weight()*(1-t)*(1-t)*factor);
                    }
              }
            ir.SetSize(cnt);
            return ir;
          }
        return IntegrationRule(et, intorder);
      }
    if (et == ET_SEGM)
      {
        // curve source. Far away: plain Gauss rule. Otherwise: project x onto the segment and use a
        // sinh-graded composite Gauss rule s = rho*sinh(u), which turns the 1/r and 1/r^2 peaks at
        // distance rho into smooth functions of u (panels of unit length in u, ~5 points each).
        MappedIntegrationPoint<1,3> mipmid(IntegrationPoint(0.5, 0, 0), trafo);
        double J = L2Norm(mipmid.GetJacobian());
        if (!(J > 0) || L2Norm(x - mipmid.GetPoint()) > 10*J)    // degenerate segment, or far away
          return IntegrationRule(et, intorder);

        // foot point: Gauss-Newton for (y(t)-x).y'(t) = 0, clamped to [0,1]; one step is exact for straight segments
        double t0 = 0.5;
        for (int it = 0; it < 5; it++)
          {
            MappedIntegrationPoint<1,3> mipt(IntegrationPoint(t0, 0, 0), trafo);
            Vec<3> dy = mipt.GetJacobian().Col(0);
            if (L2Norm2(dy) < 1e-30) break;
            t0 -= InnerProduct(Vec<3>(mipt.GetPoint() - x), dy) / L2Norm2(dy);
            t0 = std::min(1.0, std::max(0.0, t0));
          }
        MappedIntegrationPoint<1,3> mip0(IntegrationPoint(t0, 0, 0), trafo);
        double rho = std::max(L2Norm(x - mip0.GetPoint()), 1e-10*J);   // distance to the curve, floored for x on the wire

        double umin = -asinh(t0*J/rho), umax = asinh((1-t0)*J/rho);    // u-range mapping onto t in [0,1]
        int npan = int(ceil(umax - umin));
        double du = (umax - umin) / npan;
        IntegrationRule irgauss(ET_SEGM, std::max(intorder, 9));
        IntegrationRule ir(npan*irgauss.Size(), lh);
        size_t cnt = 0;
        for (int k = 0; k < npan; k++)
          for (auto & ip : irgauss)
            {
              double u = umin + (k + ip(0)) * du;
              ir[cnt++] = IntegrationPoint(t0 + rho*sinh(u)/J, 0, 0, ip.Weight() * du * rho*cosh(u)/J);   // dt = ds/J
            }
        return ir;
      }
    if (et != ET_TRIG && et != ET_QUAD)
      return IntegrationRule(et, intorder);


    IntegrationPoint ip = et == ET_TRIG ? IntegrationPoint(1./3, 1./3) : IntegrationPoint(1./2, 1./2);
    MappedIntegrationPoint<2,3>  mip(ip, trafo);
    double elsize = L2Norm(mip.GetJacobian());
    double dist = L2Norm(x-mip.GetPoint());

    if (dist < elsize)
      {
        // Find the projection of x onto the curved triangle/quad.
        IntegrationPoint ip = ProjectPointToReference(x, trafo);

        // Split the reference element into triangles meeting at the projection.
        IntegrationRule irsegm(ET_SEGM, intorder);
        int ncorners = et == ET_TRIG ? 3 : 4;
        IntegrationRule ir(ncorners*irsegm.Size()*irsegm.Size(), lh);
        size_t cnt = 0;

        Vec<2> corners[] = {Vec<2>(0,0), Vec<2>(1,0), Vec<2>(1,1), Vec<2>(0,1)};
        if (et == ET_TRIG)
          corners[2] = Vec<2>(0,1);
        for (int j = 0; j < ncorners; j++)
          {
            Vec<2> v0 = corners[j];
            Vec<2> v1 = corners[(j+1)%ncorners];
            Vec<2> v2 { ip(0), ip(1) };
            Mat<2,2> sides;
            sides.Col(0) = v0-v2;
            sides.Col(1) = v1-v2;
            double factor = Det(sides);

            if (fabs(factor) > 1e-12)
              for (auto ip1 : irsegm)
                for (auto ip2 : irsegm)
                  {
                    Vec<2> ipxy = v0 + ip1(0)*(1-ip2(0))*(v1-v0) + ip2(0)*(v2-v0);
                    ir[cnt++] = IntegrationPoint(ipxy(0), ipxy(1), 0,
                                                 ip1.Weight()*ip2.Weight()*(1-ip2(0))*factor);
                  }
          }
        ir.SetSize(cnt);
        return ir;
      }
    return IntegrationRule(et, intorder);
  }


  bool IsPotentialNearfieldSourceElement(Vec<3> x, const ElementTransformation & trafo)
  {
    auto ball = PotentialNearfieldBall (trafo);
    return ball && L2Norm(x-get<0>(*ball)) < get<1>(*ball);
  }


  template <typename TSCAL> template <typename T>
  void PotentialCF<TSCAL> ::
  AddSourceElementContribution(const BaseMappedIntegrationPoint & mip,
                               ElementId ei,
                               const IntegrationRule & ir,
                               FlatVector<T> result,
                               T scale,
                               LocalHeap & lh) const
  {
    auto space = this->gf->GetFESpace();
    auto mesh = space->GetMeshAccess();

    const FiniteElement &fel = space->GetFE(ei, lh);
    if (fel.GetNDof() == 0) return;   // no dofs on this element: nothing to integrate
    const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);

    Array<DofId> dnums(fel.GetNDof(), lh);
    space->GetDofNrs(ei, dnums);
    FlatVector<T> elvec(fel.GetNDof(), lh);
    gf->GetElementVector(dnums, elvec);

    SIMD_IntegrationRule simd_ir(ir, lh);
    Vector<SIMD<T>> simd_result(Dimension());
    simd_result = SIMD<T>(0.0);

    static constexpr int bs = 64;
    for (int k = 0; k < simd_ir.Size(); k += bs)
      {
        HeapReset hr(lh);
        auto simd_ir_range = simd_ir.Range(k, min(simd_ir.Size(), size_t(k+bs)));
        auto & miry = trafo(simd_ir_range, lh);
        FlatMatrix<SIMD<T>> vals(evaluator->Dim(), miry.Size(), lh);

        evaluator->Apply(fel, miry, elvec, vals);
        kernel->AddPotential(mip, miry, vals, simd_result, source_vb);
      }

    for (int i = 0; i < Dimension(); i++)
      result(i) += scale * HSum(simd_result(i));
  }


  template <typename TSCAL> template <typename T>
  void PotentialCF<TSCAL> ::
  AddTangentCorrection(const BaseMappedIntegrationPoint & mip,
                       ElementId ei,
                       const IntegrationRule & ir,
                       FlatVector<T> result,
                       LocalHeap & lh) const
  {
    auto formula = kernel->GetAnalyticTriangleFormula();
    if (formula == AnalyticTriangleFormula::none)
      throw Exception("no analytic triangle formula available for "+kernel->Name());

    auto space = this->gf->GetFESpace();
    auto mesh = space->GetMeshAccess();

    const FiniteElement &fel = space->GetFE(ei, lh);
    if (fel.GetNDof() == 0) return;   // no dofs on this element: nothing to correct
    const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);
    auto et = trafo.GetElementType();

    if (et == ET_TET)
      {
        // Volume analogue of the tangent correction below. Both identities
        // reduce the integral over the flat tet to face integrals of G, which
        // LaplaceSL_Polygon gives exactly:
        //   int_T grad_y G dy = sum_f n_f int_f G dS            (n_f outward)
        //   int_T G dy        = 1/2 sum_f d_f int_f G dS        (from
        //     div_y [(y-x) G] = 3G + (y-x).grad_y G = 2G, and (y-x).n_f being
        //     the constant face distance d_f on a flat face)
        // DiffLaplaceSLKernel evaluates grad_x G = -grad_y G, hence the sign
        // on the vector branch. LaplaceDL cannot reach here: its source needs
        // a normal, so it is rejected for volume sources before this point.
        if (formula != AnalyticTriangleFormula::laplace_sl &&
            formula != AnalyticTriangleFormula::laplace_grad_sl)
          return;

        Array<DofId> dnums(fel.GetNDof(), lh);
        space->GetDofNrs(ei, dnums);
        FlatVector<T> elvec(fel.GetNDof(), lh);
        gf->GetElementVector(dnums, elvec);

        Vec<3> x = mip.GetPoint();
        IntegrationPoint ip0 = ProjectPointToReferenceTet(x, trafo);
        MappedIntegrationPoint<3,3> mip0(ip0, trafo);
        Vec<3> xi0 { ip0(0), ip0(1), ip0(2) };
        Vec<3> p0 = mip0.GetPoint();
        Mat<3,3> jac = mip0.GetJacobian();
        double measure0 = mip0.GetMeasure();

        // affine image of the reference tet, tangent to the element at p0
        auto verts = ElementTopology::GetVertices(ET_TET);
        Vec<3> v[4];
        for (int i = 0; i < 4; i++)
          {
            Vec<3> c(verts[i][0], verts[i][1], verts[i][2]);
            v[i] = p0 + jac * (c-xi0);
          }
        Vec<3> centroid = 0.25*(v[0]+v[1]+v[2]+v[3]);

        double analytic_sl = 0.0;
        Vec<3> analytic_grad { 0.0, 0.0, 0.0 };
        auto faces = ElementTopology::GetFaces(ET_TET);
        for (int f = 0; f < 4; f++)
          {
            Vec<3> fv[3] = { v[faces[f][0]], v[faces[f][1]], v[faces[f][2]] };
            Vec<3> n = Cross(fv[1]-fv[0], fv[2]-fv[0]);
            double nlen = L2Norm(n);
            if (nlen < 1e-30)
              continue;
            n /= nlen;
            if (InnerProduct(fv[0]-centroid, n) < 0)   // make it outward
              n *= -1;
            FlatArray<Vec<3>> polygon(3, fv);
            double face_sl = LaplaceSL_Polygon(polygon, x);
            analytic_grad -= face_sl * n;
            analytic_sl += 0.5 * InnerProduct(fv[0]-x, n) * face_sl;
          }

        double scalar_correction = 0.0;
        Vec<3> grad_correction { 0.0, 0.0, 0.0 };
        Vec<3> nx{0.0}, ny{0.0};

        if (formula == AnalyticTriangleFormula::laplace_sl)
          {
            double flat_numeric = 0.0;
            LaplaceSLKernel<3> singularity;
            for (auto ip : ir)
              {
                Vec<3> xi { ip(0), ip(1), ip(2) };
                Vec<3> y = p0 + jac * (xi-xi0);
                if (L2Norm(x-y) > 0)
                  flat_numeric += ip.Weight() * measure0 *
                    singularity.Evaluate(x, y, nx, ny)(0);
              }
            scalar_correction = analytic_sl - flat_numeric;
          }
        else
          {
            Vec<3> flat_numeric { 0.0, 0.0, 0.0 };
            DiffLaplaceSLKernel<3> singularity;
            for (auto ip : ir)
              {
                Vec<3> xi { ip(0), ip(1), ip(2) };
                Vec<3> y = p0 + jac * (xi-xi0);
                if (L2Norm(x-y) > 0)
                  flat_numeric += ip.Weight() * measure0 *
                    singularity.Evaluate(x, y, nx, ny);
              }
            grad_correction = analytic_grad - flat_numeric;
          }

        FlatVector<T> vals(evaluator->Dim(), lh);
        evaluator->Apply(fel, mip0, elvec, vals, lh);
        for (auto term : kernel->Terms())
          {
            double correction =
              formula == AnalyticTriangleFormula::laplace_grad_sl ?
              grad_correction(term.kernel_comp) : scalar_correction;
            result(term.test_comp) += term.fac * correction * vals(term.trial_comp);
          }
        return;
      }

    if (et != ET_TRIG && et != ET_QUAD)
      return;

    Array<DofId> dnums(fel.GetNDof(), lh);
    space->GetDofNrs(ei, dnums);
    FlatVector<T> elvec(fel.GetNDof(), lh);
    gf->GetElementVector(dnums, elvec);

    Vec<3> x = mip.GetPoint();
    IntegrationPoint ip0 = ProjectPointToReference(x, trafo);
    MappedIntegrationPoint<2,3> mip0(ip0, trafo);
    Vec<2> xi0 { ip0(0), ip0(1) };
    Vec<3> p0 = mip0.GetPoint();
    Mat<3,2> jac = mip0.GetJacobian();

    Vec<2> corners[] = {
      Vec<2>(0,0), Vec<2>(1,0), Vec<2>(1,1), Vec<2>(0,1)
    };
    int ncorners = et == ET_TRIG ? 3 : 4;
    if (et == ET_TRIG)
      corners[2] = Vec<2>(0,1);
    Vec<3> v[4];
    for (int i = 0; i < ncorners; i++)
      v[i] = p0 + jac * (corners[i]-xi0);
    FlatArray<Vec<3>> polygon(ncorners, v);

    double scalar_correction = 0.0;
    Vec<3> grad_correction { 0.0, 0.0, 0.0 };
    Complex cf_correction = 0.0;
    double measure0 = mip0.GetMeasure();
    // Subtract the tangent kernel using the same rule as the curved kernel.
    Vec<3> nx{0.0};
    Vec<3> ny = mip0.GetNV();

    if (formula == AnalyticTriangleFormula::laplace_sl)
      {
        double flat_numeric = 0.0;
        double analytic = LaplaceSL_Polygon(polygon, x);
        LaplaceSLKernel<3> singularity;
        for (auto ip : ir)
          {
            Vec<2> xi { ip(0), ip(1) };
            Vec<3> y = p0 + jac * (xi-xi0);
            double r = L2Norm(x-y);
            if (r > 0)
              flat_numeric += ip.Weight() * measure0 *
                singularity.Evaluate(x, y, nx, ny)(0);
          }
        scalar_correction = analytic - flat_numeric;
      }
    else if (formula == AnalyticTriangleFormula::laplace_dl)
      {
        double flat_numeric = 0.0;
        double analytic = LaplaceDL_Polygon(polygon, x, ny);
        LaplaceDLKernel<3> singularity;
        for (auto ip : ir)
          {
            Vec<2> xi { ip(0), ip(1) };
            Vec<3> y = p0 + jac * (xi-xi0);
            double r = L2Norm(x-y);
            if (r > 0)
              flat_numeric += ip.Weight() * measure0 *
                singularity.Evaluate(x, y, nx, ny)(0);
          }
        scalar_correction = analytic - flat_numeric;
      }
    else if (formula == AnalyticTriangleFormula::laplace_grad_sl)
      {
        Vec<3> flat_numeric { 0.0, 0.0, 0.0 };
        Vec<3> analytic = LaplaceGradSL_Polygon(polygon, x);
        DiffLaplaceSLKernel<3> singularity;
        for (auto ip : ir)
          {
            Vec<2> xi { ip(0), ip(1) };
            Vec<3> y = p0 + jac * (xi-xi0);
            double r = L2Norm(x-y);
            if (r > 0)
              flat_numeric += ip.Weight() * measure0 *
                singularity.Evaluate(x, y, nx, ny);
          }
        grad_correction = analytic - flat_numeric;
      }
    else if (formula == AnalyticTriangleFormula::helmholtz_cf)
      {
        double sl_correction = LaplaceSL_Polygon(polygon, x);
        double dl_correction = LaplaceDL_Polygon(polygon, x, ny);
        LaplaceSLKernel<3> sl_singularity;
        LaplaceDLKernel<3> dl_singularity;
        for (auto ip : ir)
          {
            Vec<2> xi { ip(0), ip(1) };
            Vec<3> y = p0 + jac * (xi-xi0);
            double r = L2Norm(x-y);
            if (r > 0)
              {
                sl_correction -= ip.Weight() * measure0 * sl_singularity.Evaluate(x, y, nx, ny)(0);
                dl_correction -= ip.Weight() * measure0 * dl_singularity.Evaluate(x, y, nx, ny)(0);
              }
          }
        cf_correction = dl_correction - Complex(0,1) * kernel->GetKappa() * sl_correction;
      }

    FlatVector<T> vals(evaluator->Dim(), lh);
    evaluator->Apply(fel, mip0, elvec, vals, lh);
    for (auto term : kernel->Terms())
      {
        if (formula == AnalyticTriangleFormula::helmholtz_cf)
          {
            if constexpr (std::is_same_v<T,Complex>)
              result(term.test_comp) += term.fac * cf_correction * vals(term.trial_comp);
            continue;
          }
        double correction =
          formula == AnalyticTriangleFormula::laplace_grad_sl ?
          grad_correction(term.kernel_comp) : scalar_correction;
        result(term.test_comp) += term.fac * correction * vals(term.trial_comp);
      }
  }


  template <typename TSCAL> template <typename T>
  void PotentialCF<TSCAL> ::
  AddLocalExpansionNearfieldCorrection(const BaseMappedIntegrationRule & bmir,
                                       BareSliceMatrix<T> result) const
  {
    LocalHeap lh(10*1000*1000, "Potential::LocalExpansionNearfieldCorrection");
    auto space = this->gf->GetFESpace();
    auto mesh = space->GetMeshAccess();
    auto formula = kernel->GetAnalyticTriangleFormula();

    Array<int> candidates;
    for (int ix = 0; ix < bmir.Size(); ix++)
      {
        const auto & mip = bmir[ix];
        FlatVector<T> row = result.Row(ix).Range(0, Dimension());
        Vec<3> x = mip.GetPoint();

        // the elements whose box contains x, in element order
        candidates.SetSize0();
        if (near_sources->tree)
          {
            netgen::Point<3> p(x(0), x(1), x(2));
            near_sources->tree->GetFirstIntersecting (p, p, [&] (int k) { candidates.Append (k); return false; });
          }
        QuickSort (candidates);
        for (int k : candidates)
          {
            HeapReset hr(lh);
            ElementId ei(source_vb, near_sources->elnr[k]);

            const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);
            if (!IsPotentialNearfieldSourceElement(x, trafo))
              continue;

            IntegrationRule near_ir = GetIntegrationRule(x, trafo, intorder, lh);
            // Replace the expansion's standard source quadrature by Duffy.
            IntegrationRule standard_ir(trafo.GetElementType(), intorder);
            AddSourceElementContribution(mip, ei, standard_ir, row, T(-1.0), lh);
            AddSourceElementContribution(mip, ei, near_ir, row, T(1.0), lh);

            // Defined for surface elements and, via the divergence theorem,
            // for tets.
            if (formula != AnalyticTriangleFormula::none)
              {
                auto set = trafo.GetElementType();
                if (set == ET_TRIG || set == ET_QUAD || set == ET_TET)
                  AddTangentCorrection(mip, ei, near_ir, row, lh);
              }
          }
      }
  }


  template <typename TSCAL> template <typename T>
  void PotentialCF<TSCAL> :: T_Evaluate(const BaseMappedIntegrationPoint & mip,
                                         FlatVector<T> result) const
  {
    static Timer t("ngbem evaluate potential (ip)"); RegionTimer reg(t);
    LocalHeapMem<1000000> lh("Potential::Eval");
    auto space = this->gf->GetFESpace();
    auto mesh = space->GetMeshAccess();
    auto formula = kernel->GetAnalyticTriangleFormula();

    Vector<SIMD<T>> simd_result(Dimension());
    simd_result = SIMD<T>(0.0);
    Vector<T> correction_result(Dimension());
    correction_result = T(0.0);
    if constexpr (std::is_same<TSCAL,T>())
      for (size_t i = 0; i < mesh->GetNE(source_vb); i++)
        {
          HeapReset hr(lh);
          ElementId ei(source_vb, i);
          if (!space->DefinedOn(ei)) continue;
          if (definedon &&  !(*definedon).Mask().Test(mesh->GetElIndex(ei))) continue;

          const FiniteElement &fel = space->GetFE(ei, lh);
          if (fel.GetNDof() == 0) continue;   // no dofs on this element (e.g. L2 on curves): nothing to integrate
          const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);

          Array<DofId> dnums(fel.GetNDof(), lh);
          space->GetDofNrs(ei, dnums);
          FlatVector<T> elvec(fel.GetNDof(), lh);
          gf->GetElementVector(dnums, elvec);

          // Defined for surface elements and, via the divergence theorem,
          // for tets.
          bool use_tangent_correction = false;
          if (formula != AnalyticTriangleFormula::none)
            {
              auto set = trafo.GetElementType();
              if (set == ET_TRIG || set == ET_QUAD || set == ET_TET)
                use_tangent_correction = IsPotentialNearfieldSourceElement(mip.GetPoint(), trafo);
            }

          IntegrationRule ir = GetIntegrationRule(mip.GetPoint(), trafo, intorder, lh);

          SIMD_IntegrationRule simd_ir(ir, lh);

          static constexpr int bs = 64;
          for (int k = 0; k < simd_ir.Size(); k += bs)
            {
              HeapReset hr(lh);
              auto simd_ir_range = simd_ir.Range(k, min(simd_ir.Size(), size_t(k+bs)));
              auto & miry = trafo(simd_ir_range, lh);
              FlatMatrix<SIMD<T>> vals(evaluator->Dim(), miry.Size(), lh);

              evaluator->Apply (fel, miry, elvec, vals);
              kernel->AddPotential(mip, miry, vals, simd_result, source_vb);
            }
          if (formula != AnalyticTriangleFormula::none)
            if (use_tangent_correction)
              AddTangentCorrection(mip, ei, ir, correction_result, lh);
        }
    for (int i = 0; i < Dimension(); i++)
      result(i) = HSum(simd_result(i)) + correction_result(i);
  }


  template <typename TSCAL> template <typename T>
  void PotentialCF<TSCAL> :: T_Evaluate(const BaseMappedIntegrationRule & bmir,
                                         BareSliceMatrix<T> result) const
  {
    if constexpr (std::is_same<TSCAL,T>())
      if (local_expansion)
        {
          // static Timer t("ngbem evaluate potential, local expansion (bmir)"); RegionTimer reg(t);

          auto & target = kernel->Target();
          bool target_needs_normal = target.NeedsNormal();
          const MappedIntegrationRule<2,3> * mir23 = nullptr;
          if (target_needs_normal)
            mir23 = &dynamic_cast<const MappedIntegrationRule<2,3>&>(bmir);

          for (int j = 0; j < bmir.Size(); j++)
            {
              Vec<3> nx = 0.0;
              if (target_needs_normal)
                nx = (*mir23)[j].GetNV();
              target.EvaluateMP (*local_expansion, Vec<3>(bmir[j].GetPoint()), nx, make_BareSliceVector(result.Row(j)));
            }
          AddLocalExpansionNearfieldCorrection(bmir, result);
          return;
        }

    for (int i = 0; i < bmir.Size(); i++)
      T_Evaluate(bmir[i], result.Row(i).Range(0,Dimension()));
  }


  template <typename TSCAL> template <typename T>
  void PotentialCF<TSCAL> :: T_Evaluate(const SIMD_BaseMappedIntegrationRule & ir,
                                         BareSliceMatrix<SIMD<T>> result) const
  {
    static Timer t("ngbem evaluate potential (ir-simd), throwing"); RegionTimer reg(t);
    throw ExceptionNOSIMD ("PotentialCF::Evaluate (SIMD) not available");

    result.AddSize(Dimension(), ir.Size()) = SIMD<T>(0.0);
    return;
  }


  template class PotentialCF<double>;
  template class PotentialCF<Complex>;
}
