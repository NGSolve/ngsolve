#include <netgen_config.hpp>

#include "contact.hpp"

#undef NETGEN_USE_GUI // TODO: implement interface in netgen to draw lines (to avoid linking opengl here)

#if NETGEN_USE_GUI
#include <incopengl.hpp>
#include <visual.hpp>
#include <visualization/vssolution.hpp>
#endif // NETGEN_USE_GUI


namespace ngcomp
{
  inline int GetDomIn(const MeshAccess& ma, const Ngs_Element& el)
  {
    if(ma.GetDimension() ==3)
      return ma.GetNetgenMesh()->GetFaceDescriptor(el.GetIndex()+1).DomainIn();
    else
      return ma.GetNetgenMesh()->LineSegments()[el.Nr()].domin;
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
  double FindClosestPoint( Vec<DIMR> pmaster, Vec<DIMR> n, double h, const ElementTransformation & trafo, IntegrationPoint & ip, Vec<DIMR> & p, bool both_sides)
  // input arguments: pmaster, n, h (maximum distance), trafo
  // output arguments: ip, p
  {
    Vec<DIMS> min_lam = 0;
    double min_dist = 1e99;

    min_lam = 1./(DIMS+1);

    // Todo: line search, stop criterion
    for([[maybe_unused]] auto i : Range(4) )
    {
      ip = min_lam;
      MappedIntegrationPoint<DIMS, DIMR> mip{ip, trafo};
      T2<DIMS> t2{mip, pmaster};
      bool is_front = InnerProduct(n, mip.GetNV()) < 0;
      if(both_sides)
        is_front = true;

      if constexpr (DIMS==1)
      {
        // check end points
        for(double lam : {0.,1.})
        {
          auto dist = t2(lam);
          if(is_front && dist<min_dist)
          {
            min_dist = dist;
            min_lam = lam;
          }
        }

        auto lam = t2.CalcMinimum();
        auto dist = t2(lam);
        if(is_front && lam[0]>0 && lam[0] < 1)
        {
          min_dist = dist;
          min_lam = lam;
        }
      }
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
  optional<ContactPair<DIM>> T_GapFunction<DIM> :: CreateContactPair(const MappedIntegrationPoint<DIM-1, DIM>& mip1, LocalHeap& lh, bool both_sides) const
  {
    HeapReset hr(lh);
    // int intorder2 = 10*displacement->GetFESpace()->GetOrder();
    auto & ip1 = mip1.IP();
    auto & trafo1 = mip1.GetTransformation();
    const auto & el1 = ma->GetElement(trafo1.GetElementId());
    auto & trafo1_def = trafo1.AddDeformation(displacement.get(), lh);
    const auto & mip1_def = static_cast<const MappedIntegrationPoint<DIM-1, DIM>&>(trafo1_def(ip1, lh));
    double inv_fac = GetDomIn(*ma, el1) == 0 ? -1. : 1.;

    const auto & p1 = mip1_def.GetPoint();

    // find closest point

    double mindist = h;
    IntegrationPoint ip2_min;
    bool intersect = false;
    int el2_min(-1);


    // find all bound-2 elements closer to p1 then h

    netgen::Point<DIM> ngp1;
    for (int j = 0; j < DIM; j++) ngp1(j) = p1(j);

    auto hcurrent = h/(1024.*1024.);
    int found = 2;
    while(found>0 && hcurrent<=h)
    {

     netgen::Box<DIM> box(ngp1, ngp1);
     box.Increase(hcurrent);

    Vec<DIM> p2_min;
    searchtree->GetFirstIntersecting
      (box.PMin(), box.PMax(),
       [&] (int elnr2)
       {
         auto el2 = ma->GetElement(ElementId(BND, elnr2));
         HeapReset hr(lh);
         double inv_fac2 = GetDomIn(*ma, el2) == 0 ? -1. : 1.;

         auto & trafo2 = ma->GetTrafo(el2, lh);
         auto & trafo2_def = trafo2.AddDeformation(displacement.get(), lh);

         IntegrationPoint ip2;
         Vec<DIM> p2;
         double dist = FindClosestPoint<DIM-1,DIM>(p1, inv_fac * inv_fac2 * mip1_def.GetNV(), mindist, trafo2_def, ip2, p2, both_sides );
         if(dist<mindist && dist < h)
         {
           mindist = dist;
           el2_min = el2.Nr();
           ip2_min = ip2;
           p2_min = p2;
           intersect = true;
         }
         return false;
       });
    if(mindist < hcurrent)
      found--;

    hcurrent *= 2;
    }

    if(intersect)
    {
      return ContactPair<DIM>{el1, ElementId(BND,el2_min),
          ip1, ip2_min};
    }
    return nullopt;
  }

  template<int DIM>
  void T_GapFunction<DIM> :: Update(shared_ptr<GridFunction> displacement_, int intorder2, double h_, bool both_sides)
  {
    h = h_;
    this->both_sides = both_sides;

    displacement = displacement_;

    LocalHeap lh(1000000, "T_GapFunction::Update");
    netgen::Box<DIM> bbox{netgen::Box<DIM>::EMPTY_BOX};
    double maxh = 0;
    for (Ngs_Element el2 : ma->Elements(BND))
      {
        HeapReset hr(lh);
        auto & mask = other.Mask();
        if (!mask.Test(el2.GetIndex())) continue;

        auto & trafo2 = ma->GetTrafo (el2, lh);
        auto & trafo2_def = trafo2.AddDeformation(displacement.get(), lh);

        IntegrationRule ir2_(trafo2.GetElementType(), intorder2);
        for(auto ir2 : ir2_.Split())
        {
            HeapReset hr(lh);
            MappedIntegrationRule<DIM-1, DIM> mir2(ir2, trafo2, lh);
            MappedIntegrationRule<DIM-1, DIM> mir2_def(ir2, trafo2_def, lh);

            netgen::Box<DIM> el_box{netgen::Box<DIM>::EMPTY_BOX};
            for (auto & mip : mir2_def)
              {
                netgen::Point<DIM> p;
                for (int j = 0; j < DIM; j++)
                  p(j) = mip.GetPoint()(j);
                bbox.Add(p);
                if(h==0.0)
                  el_box.Add(p);
              }
            maxh = max(maxh, el_box.Diam());
        }
      }

    // Default-value for h is 2 * maximum_element_diameter
    if(h==0.0)
      h = 2*maxh;

    bbox.Scale(2); // make sure we don't add boxes outside of tree bounding box
    searchtree = make_unique<netgen::BoxTree<DIM, int>>(bbox);
    for (Ngs_Element el2 : ma->Elements(BND))
      {
        HeapReset hr(lh);
        auto & mask = other.Mask();
        if (!mask.Test(el2.GetIndex())) continue;

        auto & trafo2 = ma->GetTrafo (el2, lh);
        auto & trafo2_def = trafo2.AddDeformation(displacement.get(), lh);

        IntegrationRule ir2_(trafo2.GetElementType(), intorder2);
        netgen::Box<DIM> elbox{netgen::Box<DIM>::EMPTY_BOX};
        for(auto ir2 : ir2_.Split())
          {
            HeapReset hr(lh);
            MappedIntegrationRule<DIM-1, DIM> mir2(ir2, trafo2, lh);
            MappedIntegrationRule<DIM-1, DIM> mir2_def(ir2, trafo2_def, lh);

            for (auto & mip : mir2_def)
              {
                netgen::Point<DIM> p;
                for (int j = 0; j < DIM; j++)
                  p(j) = mip.GetPoint()(j);
                elbox.Add(p);
              }

        }
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
    const auto & el1 = ma->GetElement(trafo1.GetElementId());
    result = 0;
    if (!master.Mask().Test(el1.GetIndex())) return;

    // int intorder2 = 2*displacement->GetFESpace()->GetOrder();

    auto & trafo1_def = trafo1.AddDeformation(displacement.get(), lh);

    double inv_fac = GetDomIn(*ma, el1) == 0 ? -1. : 1.;

    auto & ip1 = ip.IP();
    Vec<DIM> p1;
    trafo1_def.CalcPoint(ip1, p1);

    double mindist = 1e99;
    result = std::numeric_limits<double>::infinity();

    // find all bound-2 elements closer to p1 than h
    netgen::Point<DIM> ngp1;
    for (int j = 0; j < DIM; j++)
      ngp1(j) = p1(j);

    auto & mip = static_cast<const DimMappedIntegrationPoint<DIM>&>(ip);

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
         auto el2 = ma->GetElement( ElementId (BND, elnr2) );
         double inv_fac2 = GetDomIn(*ma, el2) == 0 ? -1. : 1.;
         HeapReset hr(lh);

         bool common_vertex = false;
         for (auto s_v : el1.Vertices() )
           for (auto v : el2.Vertices() )
             if(s_v==v)
               common_vertex = true;
         if (common_vertex) return false;
         auto & trafo2 = ma->GetTrafo (el2, lh);
         auto & trafo2_def = trafo2.AddDeformation(displacement.get(), lh);

         Vec<DIM> p2;
         IntegrationPoint ip2;
         double dist = FindClosestPoint<DIM-1,DIM>(p1, inv_fac * inv_fac2 * mip.GetNV(), mindist, trafo2_def, ip2, p2, both_sides);
         if(dist<mindist && dist < h)
         {
           mindist = dist;
           result = p2-p1;
         }
         return false;
       });

    if(mindist < hcurrent)
      found--;

    hcurrent *= 2;
    }
  }

  template<int DIM>
  void T_GapFunction<DIM> :: Evaluate(const BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<> hresult) const
  {
    auto result = hresult.AddSize(mir.Size(), Dimension());
    for (auto i : Range(mir))
      Evaluate(mir[i], result.Row(i));
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
    cf->TraverseTree
      ([&](CoefficientFunction& nodecf)
       {
         auto proxy = dynamic_cast<ProxyFunction*>(&nodecf);
         if (proxy && !proxy->IsTestFunction() && !trial_proxies.Contains(proxy))
           trial_proxies.Append(proxy);
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
    ProxyUserData ud(trial_proxies.Size(), lh);
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
    cf->TraverseTree
      ([&](CoefficientFunction& nodecf)
       {
         auto proxy = dynamic_cast<ProxyFunction*>(&nodecf);
         if (proxy && !proxy->IsTestFunction() && !trial_proxies.Contains(proxy))
           trial_proxies.Append(proxy);

         if(proxy && proxy->IsTestFunction() && !test_proxies.Contains(proxy))
           test_proxies.Append(proxy);
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
  }

  void ContactIntegrator::ApplyAdd(const FiniteElement& primary_fel,
                                   const FiniteElement& secondary_fel,
                                   const BaseMappedIntegrationRule& primary_mir,
                                   FlatVector<double> elx,
                                   FlatVector<double> ely,
                                   LocalHeap& lh)
  {
    HeapReset hr(lh);
    ProxyUserData ud(trial_proxies.Size(), lh);
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

        IntRange test_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*primary_fel.GetNDof(), ely.Size()) : IntRange(0, proxy->Evaluator()->BlockDim()*primary_fel.GetNDof());
        ely1 = 0.;
        if(proxy->IsOther())
          proxy->Evaluator()->ApplyTrans(secondary_fel, *primary_mir.GetOtherMIR(), proxyvalues, ely1.Range(test_range), lh);
        else
          proxy->Evaluator()->ApplyTrans(primary_fel, primary_mir, proxyvalues, ely1.Range(test_range), lh);
        ely += ely1;
      }
  }

  void ContactIntegrator::CalcLinearizedAdd(const FiniteElement& primary_fel,
                                            const FiniteElement& secondary_fel,
                                            const BaseMappedIntegrationRule& primary_mir,
                                            FlatVector<double> elx,
                                            FlatMatrix<double> elmat,
                                            LocalHeap& lh)
  {
    HeapReset hr(lh);
    ProxyUserData ud(trial_proxies.Size(), lh);
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

    FlatMatrix<> dderiv(primary_mir.Size(), 1,lh);
    FlatMatrix<AutoDiff<1>> dval(primary_mir.Size(), 1, lh);

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies)) // ss
        {
          HeapReset hr(lh);

          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3> proxyvalues(lh, primary_mir.Size(), proxy2->Dimension(), proxy1->Dimension());
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

          for (int i = 0; i < primary_mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= primary_mir[i].GetWeight();

          IntRange trial_range  = proxy1->IsOther() ? IntRange(proxy1->Evaluator()->BlockDim()*primary_fel.GetNDof(), elmat.Width()) : IntRange(0, proxy1->Evaluator()->BlockDim()*primary_fel.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(proxy2->Evaluator()->BlockDim()*primary_fel.GetNDof(), elmat.Height()) : IntRange(0, proxy2->Evaluator()->BlockDim()*primary_fel.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), trial_range.Size(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), test_range.Size(), lh);

          int bs = primary_mir.Size();
          FlatMatrix<double,ColMajor> bdbmat1(bs*proxy2->Dimension(), trial_range.Size(), lh);
          FlatMatrix<double,ColMajor> bbmat2(bs*proxy2->Dimension(), test_range.Size(), lh);

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

  ContactBoundary::ContactBoundary(Region _master, Region _other,
                                   bool _draw_pairs, bool _volume)
    : master(_master), other(_other), draw_pairs(_draw_pairs), volume(_volume)
  {
#if NETGEN_USE_GUI
    if(draw_pairs)
      AddUserVisualizationObject (this);
#endif // NETGEN_USE_GUI

    auto mesh = master.Mesh();
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
    glBegin (GL_LINES);
    for (auto i : Range(primary_points.Size()))
      {
        auto & mp = primary_points[i];
        auto & sp = secondary_points[i];
        glVertex3d (mp(0), mp(1), mp(2));
        glVertex3d (sp(0), sp(1), sp(2));
      }
    glEnd();
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
         int intorder, double h, bool both_sides)
  {
    if(!displacement_ && !bf)
      throw Exception("Either displacement or BilinearForm needed in ContactBoundary update!");
    if(!fes)
      {
        if(bf)
          fes = bf->GetFESpace();
        else
          fes = displacement_->GetFESpace();
      }
    if(bf && (bf->GetFESpace().get() != fes.get()))
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
    if (displacement)
      gap->Update(displacement, 10*displacement->GetFESpace()->GetOrder(), h,
                  both_sides);
    else
      gap->Update(nullptr, 10, h, both_sides);
      
    auto mesh = fes->GetMeshAccess();
    if(mesh->GetDimension() == 2)
      static_pointer_cast<DisplacedNormal<2>>(normal)->Update(displacement);
    else
      static_pointer_cast<DisplacedNormal<3>>(normal)->Update(displacement);
    if(bf)
      {
        static Timer t1("Build contact pairs");
        RegionTimer regt1(t1);
        LocalHeap lh(1000000, "ContactBoundary-Update", true);
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
                      IntegrationRule ir(trafo.GetElementType(), intorder);
                      MappedIntegrationRule<DIM-1, DIM> mir(ir, trafo, lh);

                      FlatArray<IntegrationPoint> this_ir(ir.Size(), lh);
                      FlatArray<IntegrationPoint> other_ir(ir.Size(), lh);
                      FlatArray<size_t> other_nr(ir.Size(), lh);

                      int cntpair = 0;
                      for(const auto& mip : mir)                        
                        {
                          auto pair = tgap->CreateContactPair(mip, lh,
                                                              both_sides);
                          if (pair.has_value())
                            {
                              this_ir[cntpair] =  (*pair).primary_ip;
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
                          while (next < cntpair && other_nr[index[next]] == other_nr[index[first]])
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
                              bf->AddSpecialElement(make_unique<MPContactElement<DIM>>
                                                    (el, ElementId(BND, other_nr[index[first]]),
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


  template<int DIM>
  MPContactElement<DIM>::MPContactElement(ElementId aprimary_ei, ElementId asecondary_ei,
                                          IntegrationRule aprimary_ir, IntegrationRule asecondary_ir,
                                          shared_ptr<ContactBoundary> _cb,
                                          GridFunction* _deformation)
    : primary_ei(aprimary_ei), secondary_ei(asecondary_ei),
      primary_ir(std::move(aprimary_ir)), secondary_ir(std::move(asecondary_ir)),
      cb(_cb), fes(_cb->GetFESpace().get()), deformation(_deformation)
  { }

  template<int DIM>
  void MPContactElement<DIM>::GetDofNrs(Array<DofId>& dnums) const
  {
    fes->GetDofNrs(primary_ei, dnums);
    Array<DofId> s_dofs;
    fes->GetDofNrs(secondary_ei, s_dofs);
    dnums += s_dofs;
  }

  template<int DIM>
  double MPContactElement<DIM>::Energy(FlatVector<double> elx,
                                     LocalHeap& lh) const
  {
    if(cb->GetIntegrators().Size())
      throw Exception("Energy does not work with contact integrators!");
    
    auto& primary_trafo = fes->GetMeshAccess()->GetTrafo(primary_ei, lh);
    auto& secondary_trafo = fes->GetMeshAccess()->GetTrafo(secondary_ei, lh);
    auto& primary_deformed_trafo = primary_trafo.AddDeformation(deformation, lh);
    auto& secondary_deformed_trafo = secondary_trafo.AddDeformation(deformation, lh);

    // IntegrationRule primary_ir(1, &const_cast<IntegrationPoint&>(pair.primary_ip));
    // IntegrationRule secondary_ir(1, &const_cast<IntegrationPoint&>(pair.secondary_ip));

    auto& primary_fel = fes->GetFE(primary_ei, lh);
    auto& secondary_fel = fes->GetFE(secondary_ei, lh);
    
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

    auto& primary_trafo = fes->GetMeshAccess()->GetTrafo(primary_ei, lh);
    auto& secondary_trafo = fes->GetMeshAccess()->GetTrafo(secondary_ei, lh);
    auto& primary_deformed_trafo = primary_trafo.AddDeformation(deformation, lh);
    auto& secondary_deformed_trafo = secondary_trafo.AddDeformation(deformation, lh);

    // IntegrationRule primary_ir(1, &const_cast<IntegrationPoint&>(pair.primary_ip));
    // IntegrationRule secondary_ir(1, &const_cast<IntegrationPoint&>(pair.secondary_ip));

    auto& primary_fel = fes->GetFE(primary_ei, lh);
    auto& secondary_fel = fes->GetFE(secondary_ei, lh);
    
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
              ci->ApplyAdd(primary_fel, secondary_fel, primary_mir, elx, ely, lh);
            for(const auto& ce : energies)
              ce->ApplyAdd(primary_fel, secondary_fel, primary_mir, elx, ely, lh);
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
                  ci->ApplyAdd(primary_fel, secondary_fel, primary_mir, elx, ely, lh);
                for(const auto& ce : energies)
                  ce->ApplyAdd(primary_fel, secondary_fel, primary_mir, elx, ely, lh);

                }
          }
      }
  }


  template<int DIM>
  void MPContactElement<DIM>::CalcElementMatrix(FlatMatrix<double> elmat,
                                                LocalHeap& lh) const
  {
    HeapReset hr(lh);
    FlatVector<> elx(elmat.Height(), lh);
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

    auto& primary_trafo = fes->GetMeshAccess()->GetTrafo(primary_ei, lh);
    auto& secondary_trafo = fes->GetMeshAccess()->GetTrafo(secondary_ei, lh);
    auto& primary_deformed_trafo = primary_trafo.AddDeformation(deformation, lh);
    auto& secondary_deformed_trafo = secondary_trafo.AddDeformation(deformation, lh);

    // IntegrationRule primary_ir(1, &const_cast<IntegrationPoint&>(primary_e));
    // IntegrationRule secondary_ir(1, &const_cast<IntegrationPoint&>(pair.secondary_ip));

    auto& primary_fel = fes->GetFE(primary_ei, lh);
    auto& secondary_fel = fes->GetFE(secondary_ei, lh);

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
                ci->CalcLinearizedAdd(primary_fel, secondary_fel, primary_mir, elx, elmat, lh);            
              for(const auto& ce : energies)
                ce->CalcLinearizedAdd(primary_fel, secondary_fel, primary_mir, elx, elmat, lh);            
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
                ci->CalcLinearizedAdd(primary_fel, secondary_fel, primary_mir, elx, elmat, lh);            
              for(const auto& ce : energies)
                ce->CalcLinearizedAdd(primary_fel, secondary_fel, primary_mir, elx, elmat, lh);            
            }
        }
      }
      
  }

  template class MPContactElement<2>;
  template class MPContactElement<3>;



  
} // namespace ngcomp
