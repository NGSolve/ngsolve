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


  template <typename TSCAL>
  void PotentialCF<TSCAL> ::
  BuildLocalExpansion(const Region & reg)
  {
    LocalHeapMem<100000> lh("PotentialCF::BuildLocalExpansion");

    auto space = this->gf->GetFESpace();
    auto mesh = space->GetMeshAccess();


    Vec<3> smax(-1e99, -1e99, -1e99);
    Vec<3> smin(1e99, 1e99, 1e99);

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
      }

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


  IntegrationRule GetIntegrationRule(Vec<3> x, const ElementTransformation & trafo, int intorder)
  {
    auto et = trafo.GetElementType();
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
        IntegrationRule ir;

        Vec<2> corners[] = {Vec<2>(0,0), Vec<2>(1,0), Vec<2>(1,1), Vec<2>(0,1)};
        int ncorners = et == ET_TRIG ? 3 : 4;
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
                    ir.AddIntegrationPoint (IntegrationPoint(ipxy(0), ipxy(1), 0,
                                                             ip1.Weight()*ip2.Weight()*(1-ip2(0))*factor));
                  }
          }
        return ir;
      }
    return IntegrationRule(et, intorder);
  }


  bool IsPotentialNearfieldSourceElement(Vec<3> x, const ElementTransformation & trafo)
  {
    auto et = trafo.GetElementType();
    if (et != ET_TRIG && et != ET_QUAD)
      return false;

    IntegrationPoint ip = et == ET_TRIG ?
      IntegrationPoint(1./3, 1./3) : IntegrationPoint(1./2, 1./2);
    MappedIntegrationPoint<2,3> mip(ip, trafo);
    double elsize = L2Norm(mip.GetJacobian());
    double dist = L2Norm(x-mip.GetPoint());
    return dist < elsize;
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
    const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);

    Array<DofId> dnums(fel.GetNDof(), lh);
    space->GetDofNrs(ei, dnums);
    FlatVector<T> elvec(fel.GetNDof(), lh);
    gf->GetElementVector(dnums, elvec);

    SIMD_IntegrationRule simd_ir(ir);
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
    const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);
    auto et = trafo.GetElementType();
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

    FlatVector<T> vals(evaluator->Dim(), lh);
    evaluator->Apply(fel, mip0, elvec, vals, lh);
    for (auto term : kernel->Terms())
      {
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

    // TODO: find a better way to identify nearfield source elements.
    // The current path scans all source elements for every target point.
    for (int ix = 0; ix < bmir.Size(); ix++)
      {
        const auto & mip = bmir[ix];
        FlatVector<T> row = result.Row(ix).Range(0, Dimension());
        Vec<3> x = mip.GetPoint();

        for (size_t i = 0; i < mesh->GetNE(source_vb); i++)
          {
            HeapReset hr(lh);
            ElementId ei(source_vb, i);
            if (!space->DefinedOn(ei)) continue;
            if (definedon && !(*definedon).Mask().Test(mesh->GetElIndex(ei))) continue;

            const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);
            if (!IsPotentialNearfieldSourceElement(x, trafo))
              continue;

            IntegrationRule near_ir = GetIntegrationRule(x, trafo, intorder);
            // Replace the expansion's standard source quadrature by Duffy.
            IntegrationRule standard_ir(trafo.GetElementType(), intorder);
            AddSourceElementContribution(mip, ei, standard_ir, row, T(-1.0), lh);
            AddSourceElementContribution(mip, ei, near_ir, row, T(1.0), lh);

            if (formula != AnalyticTriangleFormula::none)
              AddTangentCorrection(mip, ei, near_ir, row, lh);
          }
      }
  }


  template <typename TSCAL> template <typename T>
  void PotentialCF<TSCAL> :: T_Evaluate(const BaseMappedIntegrationPoint & mip,
                                         FlatVector<T> result) const
  {
    static Timer t("ngbem evaluate potential (ip)"); RegionTimer reg(t);
    LocalHeapMem<100000> lh("Potential::Eval");
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
          const ElementTransformation &trafo = mesh->GetTrafo(ei, lh);

          Array<DofId> dnums(fel.GetNDof(), lh);
          space->GetDofNrs(ei, dnums);
          FlatVector<T> elvec(fel.GetNDof(), lh);
          gf->GetElementVector(dnums, elvec);

          bool use_tangent_correction = false;
          if (formula != AnalyticTriangleFormula::none)
            use_tangent_correction = IsPotentialNearfieldSourceElement(mip.GetPoint(), trafo);

          IntegrationRule ir = GetIntegrationRule(mip.GetPoint(), trafo, intorder);

          SIMD_IntegrationRule simd_ir(ir);

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
