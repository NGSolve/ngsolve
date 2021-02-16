#include <netgen_config.hpp>

#include "contact.hpp"

#if NETGEN_USE_GUI
#include <incopengl.hpp>

namespace netgen
{
  DLL_HEADER void AddUserVisualizationObject (UserVisualizationObject * vis);
  DLL_HEADER void DeleteUserVisualizationObject (UserVisualizationObject * vis);
}
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

    Vec<DIM> CalcMinimumOnSegment( Vec<DIM> x1, Vec<DIM> x2 )
    {
      auto v = x2-x1;
      Vec<DIM> va = a*v;
      double lam = InnerProduct(v,b) - InnerProduct(va,x2+x0);
      lam /= InnerProduct(va,v);
      lam = min(lam,1.0);
      lam = max(lam,0.0);
      return lam*x1 + (1-lam)*x2;
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
  double FindClosestPoint( Vec<DIMR> pmaster, Vec<DIMR> n, double h, const ElementTransformation & trafo, IntegrationPoint & ip, Vec<DIMR> & p)
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
        // TODO: Handle quad elements

        // check corner points
        ArrayMem<Vec<DIMS>, 5> points{ {0,0}, {0,1}, {1,0} };
        for (Vec<DIMS> lam : points)
        {
          auto d = t2(lam);
          if(is_front && d < min_dist)
          {
            min_dist = d;
            min_lam = lam;
          }
        }

        ArrayMem<Vec<DIMS>, 5> lam;
        ArrayMem<double, 5> dist;

        auto l = t2.CalcMinimum();
        if(l[0]>0 && l[1]>0 && l[0]<1 && l[1]<1 && l[0]+l[1]<1)
        {
          lam.Append(l);
          dist.Append(t2(lam.Last()));
        }

//         lam.Append(t2.CalcMinimumOnSegment( {0,0}, {1,0} ));
//         dist.Append(t2(lam.Last()));
// 
//         lam.Append(t2.CalcMinimumOnSegment( {1,0}, {0,1} ));
//         dist.Append(t2(lam.Last()));
// 
//         lam.Append(t2.CalcMinimumOnSegment( {0,1}, {0,0} ));
//         dist.Append(t2(lam.Last()));

        for(auto i : Range(lam.Size()))
        {
          if(is_front && dist[i]<min_dist && InnerProduct(n , mip.GetNV())<0)
          {
            min_dist = dist[i];
            min_lam = lam[i];
          }
        }
      }
    }

    if(min_dist > h)
      return min_dist;

    ip = min_lam;
    trafo.CalcPoint( ip, p );
    return L2Norm(pmaster-p);
  }

  template<int DIM>
  optional<ContactPair<DIM>> T_GapFunction<DIM> :: CreateContactPair(const MappedIntegrationPoint<DIM-1, DIM>& mip1, LocalHeap& lh) const
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
    netgen::Box<DIM> box(ngp1, ngp1);
    box.Increase(h);

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
         double dist = FindClosestPoint<DIM-1,DIM>(p1, inv_fac * inv_fac2 * mip1_def.GetNV(), mindist, trafo2_def, ip2, p2 );
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

    if(intersect)
    {
      return ContactPair<DIM>{el1, ElementId(BND,el2_min),
          ip1, ip2_min};
    }
    return nullopt;
  }

  template<int DIM>
  void T_GapFunction<DIM> :: Update(shared_ptr<GridFunction> displacement_, int intorder2, double h_)
  {
    h = h_;

    displacement = displacement_;
    auto fes = displacement->GetFESpace();
    intorder2 = 10*fes->GetOrder();

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

        IntegrationRule ir2(trafo2.GetElementType(), intorder2);
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

    // Default-value for h is 2 * maximum_element_diameter
    if(h==0.0)
      h = 2*maxh;

    searchtree = make_unique<netgen::BoxTree<DIM, int>>(bbox);
    for (Ngs_Element el2 : ma->Elements(BND))
      {
        HeapReset hr(lh);
        auto & mask = other.Mask();
        if (!mask.Test(el2.GetIndex())) continue;

        auto & trafo2 = ma->GetTrafo (el2, lh);
        auto & trafo2_def = trafo2.AddDeformation(displacement.get(), lh);

        IntegrationRule ir2(trafo2.GetElementType(), intorder2);
        MappedIntegrationRule<DIM-1, DIM> mir2(ir2, trafo2, lh);
        MappedIntegrationRule<DIM-1, DIM> mir2_def(ir2, trafo2_def, lh);

        netgen::Box<DIM> elbox{netgen::Box<DIM>::EMPTY_BOX};
        for (auto & mip : mir2_def)
          {
            netgen::Point<DIM> p;
            for (int j = 0; j < DIM; j++)
              p(j) = mip.GetPoint()(j);
            elbox.Add(p);
          }

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
         double dist = FindClosestPoint<DIM-1,DIM>(p1, inv_fac * inv_fac2 * mip.GetNV(), mindist, trafo2_def, ip2, p2 );
         if(dist<mindist && dist < h)
         {
           mindist = dist;
           result = p2-p1;
         }
         return false;
       });

    if(mindist < h)
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
    auto ma = displacement->GetMeshAccess();
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

  double ContactEnergy::CalcEnergy(const FiniteElement& m_fel,
                                   const FiniteElement& s_fel,
                                   const BaseMappedIntegrationRule& m_mir,
                                   FlatVector<double> elx,
                                   LocalHeap& lh)
  {
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(m_mir.GetTransformation()).userdata = &ud;
    ud.fel = & m_fel;

    for(auto proxy : trial_proxies)
      {
        IntRange trial_range = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim() * m_fel.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim() * m_fel.GetNDof());
        ud.AssignMemory(proxy, 1, proxy->Dimension(), lh);
        if(proxy->IsOther())
          proxy->Evaluator()->Apply(s_fel, *m_mir.GetOtherMIR(), elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(m_fel, m_mir, elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
      }

    FlatMatrix<> values(m_mir.Size(), 1, lh);
    cf->Evaluate(m_mir, values);

    double sum = 0.;
    for (int i = 0; i < m_mir.Size(); i++)
      sum += m_mir[i].GetWeight() * values(i,0);
    return sum;
  }

  void ContactEnergy::ApplyAdd(const FiniteElement& m_fel,
                               const FiniteElement& s_fel,
                               const BaseMappedIntegrationRule& m_mir,
                               FlatVector<double> elx,
                               FlatVector<double> ely,
                               LocalHeap& lh)
  {
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(m_mir.GetTransformation()).userdata = &ud;
    ud.fel = & m_fel;

    for(auto proxy : trial_proxies)
      {
        IntRange trial_range = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim() * m_fel.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim() * m_fel.GetNDof());
        ud.AssignMemory(proxy, 1, proxy->Dimension(), lh);
        if(proxy->IsOther())
          proxy->Evaluator()->Apply(s_fel, *m_mir.GetOtherMIR(), elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(m_fel, m_mir, elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
      }

    FlatVector<> ely1(ely.Size(), lh);
    FlatMatrix<AutoDiff<1>> dval(m_mir.Size(), 1, lh);

    for (auto proxy : trial_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(m_mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            cf -> Evaluate (m_mir, dval);
            for (size_t i = 0; i < m_mir.Size(); i++)
              proxyvalues(i,k) = dval(i,0).DValue(0);
          }

        for (int i = 0; i < m_mir.Size(); i++)
          proxyvalues.Row(i) *= m_mir[i].GetWeight();

        IntRange test_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*m_fel.GetNDof(), ely.Size()) : IntRange(0, proxy->Evaluator()->BlockDim()*m_fel.GetNDof());
        ely1 = 0.;
        if(proxy->IsOther())
          proxy->Evaluator()->ApplyTrans(s_fel, *m_mir.GetOtherMIR(), proxyvalues, ely1.Range(test_range), lh);
        else
          proxy->Evaluator()->ApplyTrans(m_fel, m_mir, proxyvalues, ely1.Range(test_range), lh);
        ely += ely1;
      }
  }

  void ContactEnergy::CalcLinearizedAdd(const FiniteElement& m_fel,
                                        const FiniteElement& s_fel,
                                        const BaseMappedIntegrationRule& m_mir,
                                        FlatVector<double> elx,
                                        FlatMatrix<double> elmat,
                                        LocalHeap& lh)
  {
    HeapReset hr(lh);
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(m_mir.GetTransformation()).userdata = &ud;
    ud.fel = & m_fel;

    for(auto proxy : trial_proxies)
      {
        IntRange trial_range = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim() * m_fel.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim() * m_fel.GetNDof());
        ud.AssignMemory(proxy, 1, proxy->Dimension(), lh);
        if(proxy->IsOther())
          proxy->Evaluator()->Apply(s_fel, *m_mir.GetOtherMIR(), elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(m_fel, m_mir, elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
      }

    FlatMatrix<> dderiv(m_mir.Size(), 1,lh);
    FlatMatrix<AutoDiffDiff<1,double>> ddval(m_mir.Size(), 1, lh);

    FlatArray<FlatMatrix<>> diags(trial_proxies.Size(), lh);
    for (int k1 : Range(trial_proxies))
      {
        auto proxy = trial_proxies[k1];
        diags[k1].AssignMemory(m_mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (m_mir, ddval);
            for (size_t i = 0; i < m_mir.Size(); i++)
              diags[k1](i,k) = ddval(i,0).DDValue(0);
          }
      }

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(trial_proxies))
        {
          HeapReset hr(lh);

          auto proxy1 = trial_proxies[k1];
          auto proxy2 = trial_proxies[l1];

          FlatTensor<3> proxyvalues(lh, m_mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;

                {
                  cf -> Evaluate (m_mir, ddval);
                  for (size_t i = 0; i < m_mir.Size(); i++)
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

          for (int i = 0; i < m_mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= m_mir[i].GetWeight();

          IntRange trial_range  = proxy1->IsOther() ? IntRange(proxy1->Evaluator()->BlockDim()*m_fel.GetNDof(), elmat.Width()) : IntRange(0, proxy1->Evaluator()->BlockDim()*m_fel.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(proxy2->Evaluator()->BlockDim()*m_fel.GetNDof(), elmat.Height()) : IntRange(0, proxy2->Evaluator()->BlockDim()*m_fel.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);


          int bs = m_mir.Size();
          FlatMatrix<double,ColMajor> bdbmat1(bs*proxy2->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bbmat2(bs*proxy2->Dimension(), loc_elmat.Height(), lh);

          bdbmat1 = 0.;
          bbmat2 = 0.;

          const auto& s_mir = *m_mir.GetOtherMIR();

          for (size_t j = 0; j < bs; j++)
            {
              size_t ii = j;
              IntRange r3 = proxy2->Dimension() * IntRange(j,j+1);

              {
                if(proxy1->IsOther())
                  proxy1->Evaluator()->CalcMatrix(s_fel, s_mir[ii], bmat1, lh);
                else
                  proxy1->Evaluator()->CalcMatrix(m_fel, m_mir[ii], bmat1, lh);
                if(proxy2->IsOther())
                  proxy2->Evaluator()->CalcMatrix(s_fel, s_mir[ii], bmat2, lh);
                else
                  proxy2->Evaluator()->CalcMatrix(m_fel, m_mir[ii], bmat2, lh);
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
    fes = trial_proxies[0]->GetFESpace();
  }

  void ContactIntegrator::ApplyAdd(const FiniteElement& m_fel,
                                   const FiniteElement& s_fel,
                                   const BaseMappedIntegrationRule& m_mir,
                                   FlatVector<double> elx,
                                   FlatVector<double> ely,
                                   LocalHeap& lh)
  {
    HeapReset hr(lh);
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(m_mir.GetTransformation()).userdata = &ud;
    ud.fel = & m_fel;

    for(auto proxy : trial_proxies)
      {
        IntRange trial_range = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim() * m_fel.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim() * m_fel.GetNDof());
        ud.AssignMemory(proxy, 1, proxy->Dimension(), lh);
        if(proxy->IsOther())
          proxy->Evaluator()->Apply(s_fel, *m_mir.GetOtherMIR(), elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(m_fel, m_mir, elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
      }

    FlatVector ely1(ely.Size(), lh);
    FlatMatrix<> val(m_mir.Size(), 1, lh);

    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(m_mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (m_mir, val);
            proxyvalues.Col(k) = val.Col(0);
          }

        for (int i = 0; i < m_mir.Size(); i++)
          proxyvalues.Row(i) *= m_mir[i].GetWeight();

        IntRange test_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*m_fel.GetNDof(), ely.Size()) : IntRange(0, proxy->Evaluator()->BlockDim()*m_fel.GetNDof());
        ely1 = 0.;
        if(proxy->IsOther())
          proxy->Evaluator()->ApplyTrans(s_fel, *m_mir.GetOtherMIR(), proxyvalues, ely1.Range(test_range), lh);
        else
          proxy->Evaluator()->ApplyTrans(m_fel, m_mir, proxyvalues, ely1.Range(test_range), lh);
        ely += ely1;
      }
  }

  void ContactIntegrator::CalcLinearizedAdd(const FiniteElement& m_fel,
                                            const FiniteElement& s_fel,
                                            const BaseMappedIntegrationRule& m_mir,
                                            FlatVector<double> elx,
                                            FlatMatrix<double> elmat,
                                            LocalHeap& lh)
  {
    HeapReset hr(lh);
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(m_mir.GetTransformation()).userdata = &ud;
    ud.fel = & m_fel;

    for(auto proxy : trial_proxies)
      {
        IntRange trial_range = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim() * m_fel.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim() * m_fel.GetNDof());
        ud.AssignMemory(proxy, 1, proxy->Dimension(), lh);
        if(proxy->IsOther())
          proxy->Evaluator()->Apply(s_fel, *m_mir.GetOtherMIR(), elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(m_fel, m_mir, elx.Range(trial_range),
                                    ud.GetMemory(proxy), lh);
      }

    FlatMatrix<> dderiv(m_mir.Size(), 1,lh);
    FlatMatrix<AutoDiff<1>> dval(m_mir.Size(), 1, lh);

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies)) // ss
        {
          HeapReset hr(lh);

          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3> proxyvalues(lh, m_mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;

                cf -> Evaluate (m_mir, dval);
                for (size_t i = 0; i < m_mir.Size(); i++)
                  proxyvalues(i,l,k) = dval(i,0).DValue(0);
              }

          for (int i = 0; i < m_mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= m_mir[i].GetWeight();

          IntRange trial_range  = proxy1->IsOther() ? IntRange(proxy1->Evaluator()->BlockDim()*m_fel.GetNDof(), elmat.Width()) : IntRange(0, proxy1->Evaluator()->BlockDim()*m_fel.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(proxy2->Evaluator()->BlockDim()*m_fel.GetNDof(), elmat.Height()) : IntRange(0, proxy2->Evaluator()->BlockDim()*m_fel.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);

          int bs = m_mir.Size();
          FlatMatrix<double,ColMajor> bdbmat1(bs*proxy2->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bbmat2(bs*proxy2->Dimension(), loc_elmat.Height(), lh);

          bdbmat1 = 0.;
          bbmat2 = 0.;

          const auto& s_mir = *m_mir.GetOtherMIR();

          for (size_t j = 0; j < bs; j++)
            {
              size_t ii = j;
              IntRange r3 = proxy2->Dimension() * IntRange(j,j+1);

              {
                if(proxy1->IsOther())
                  proxy1->Evaluator()->CalcMatrix(s_fel, s_mir[ii], bmat1, lh);
                else
                  proxy1->Evaluator()->CalcMatrix(m_fel, m_mir[ii], bmat1, lh);
                if(proxy2->IsOther())
                  proxy2->Evaluator()->CalcMatrix(s_fel, s_mir[ii], bmat2, lh);
                else
                  proxy2->Evaluator()->CalcMatrix(m_fel, m_mir[ii], bmat2, lh);
              }

              bdbmat1.Rows(r3) = proxyvalues(ii,STAR,STAR) * bmat1;
              bbmat2.Rows(r3) = bmat2;
            }
          loc_elmat += Trans(bbmat2) * bdbmat1;
        }
  }

  ContactBoundary::ContactBoundary(Region _master, Region _other,
                                   bool _draw_pairs)
    : master(_master), other(_other), draw_pairs(_draw_pairs)
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
    for (auto i : Range(master_points.Size()))
      {
        auto & mp = master_points[i];
        auto & sp = other_points[i];
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
  }

  void ContactBoundary::AddIntegrator(shared_ptr<CoefficientFunction> form,
                                      bool deformed)
  {
    integrators.Append(make_shared<ContactIntegrator>(form, deformed));
  }

  void ContactBoundary::Update(shared_ptr<GridFunction> displacement_,
                               shared_ptr<BilinearForm> bf,
                               int intorder, double h)
  {
    if(!fes)
      {
        if(bf)
          fes = bf->GetFESpace();
        else
          fes = displacement_->GetFESpace();
      }
    if(bf && (bf->GetFESpace().get() != fes.get()))
      throw Exception("BilinearForm on different space as given to ContactBoundary!");
    fes_displacement = displacement_->GetFESpace();

    if(draw_pairs)
    {
      other_points.SetSize(0);
      master_points.SetSize(0);
    }

    auto displacement = CreateGridFunction(displacement_->GetFESpace(), "_cb_displacement", displacement_->GetFlags());
    displacement->Update();
    displacement->GetVector() = displacement_->GetVector();
    gap->Update(displacement, intorder, h);
    auto mesh = displacement->GetFESpace()->GetMeshAccess();
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
                     auto cont_el = dynamic_cast<ContactElement<DIM>*>(specialels[index].get());
                     if(cont_el && cont_el->GetContactBoundary() == this)
                       bf->DeleteSpecialElement(index);
                   }

                 auto tgap = static_pointer_cast<T_GapFunction<DIM>>(gap);
                 static mutex add_mutex;
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
                      for(const auto& mip : mir)
                        {
                          auto pair = tgap->CreateContactPair(mip, lh);
                          if(pair.has_value())
                            {
                              lock_guard<mutex> guard(add_mutex);
                              bf->AddSpecialElement(make_unique<ContactElement<DIM>>(*pair, this, displacement.get()));
                              if(draw_pairs)
                              {
                                HeapReset hr(lh);
                                auto & t1_def = trafo.AddDeformation(displacement.get(), lh);
                                Vec<3> p1 = 0;
                                t1_def.CalcPoint(pair->master_ip, p1);
                                master_points.Append(p1);

                                auto & t2 = mesh->GetTrafo(pair->other_el, lh);
                                auto & t2_def = t2.AddDeformation(displacement.get(), lh);
                                Vec<3> p2 = 0;
                                t2_def.CalcPoint(pair->other_ip, p2);
                                other_points.Append(p2);
                              }
                            }
                        }
                    });
               }
           });
      }
  }

  template<int DIM>
  ContactElement<DIM>::ContactElement(const ContactPair<DIM>& _pair,
                                      ContactBoundary* _cb,
                                      GridFunction* _deformation)
    : pair(_pair), cb(_cb), fes(_cb->GetFESpace().get()), deformation(_deformation)
  {
  }

  template<int DIM>
  void ContactElement<DIM>::GetDofNrs(Array<DofId>& dnums) const
  {
    fes->GetDofNrs(pair.master_el, dnums);
    Array<DofId> s_dofs;
    fes->GetDofNrs(pair.other_el, s_dofs);
    dnums += s_dofs;
  }

  template<int DIM>
  double ContactElement<DIM>::Energy(FlatVector<double> elx,
                                     LocalHeap& lh) const
  {
    if(cb->GetIntegrators().Size())
      throw Exception("Energy does not work with contact integrators!");
    auto& m_trafo = fes->GetMeshAccess()->GetTrafo(pair.master_el, lh);
    double energy = 0.;
    m_trafo(pair.master_ip, lh).IntegrationRuleFromPoint
      ([&](const BaseMappedIntegrationRule& m_mir)
       {
         auto& s_trafo = fes->GetMeshAccess()->GetTrafo(pair.other_el, lh);
         s_trafo(pair.other_ip, lh).IntegrationRuleFromPoint
           ([&](const BaseMappedIntegrationRule& s_mir)
            {
              auto& fel = fes->GetFE(pair.master_el, lh);
              auto& s_fel = fes->GetFE(pair.other_el, lh);
              const_cast<BaseMappedIntegrationRule&>(m_mir).SetOtherMIR(&s_mir);
              for(const auto& ce : cb->GetEnergies())
                if(!ce->IsDeformed())
                  energy += ce->CalcEnergy(fel, s_fel, m_mir, elx, lh);
            });
       });

    m_trafo.AddDeformation(deformation, lh);
    m_trafo(pair.master_ip, lh).IntegrationRuleFromPoint
      ([&](const BaseMappedIntegrationRule& m_mir)
       {
         auto& s_trafo = fes->GetMeshAccess()->GetTrafo(pair.other_el, lh);
         s_trafo.AddDeformation(deformation, lh);
         s_trafo(pair.other_ip, lh).IntegrationRuleFromPoint
           ([&](const BaseMappedIntegrationRule& s_mir)
            {
              auto& fel = fes->GetFE(pair.master_el, lh);
              auto& s_fel = fes->GetFE(pair.other_el, lh);
              const_cast<BaseMappedIntegrationRule&>(m_mir).SetOtherMIR(&s_mir);
              for(const auto& ce : cb->GetEnergies())
                if(ce->IsDeformed())
                  energy += ce->CalcEnergy(fel, s_fel, m_mir, elx, lh);
            });
       });

    return energy;
  }

  template<int DIM>
  void ContactElement<DIM>::Apply(FlatVector<double> elx,
                                  FlatVector<double> ely,
                                  LocalHeap& lh) const
  {
    auto& m_trafo = fes->GetMeshAccess()->GetTrafo(pair.master_el, lh);
    ely = 0.;
    m_trafo(pair.master_ip, lh).IntegrationRuleFromPoint
      ([&](const BaseMappedIntegrationRule& m_mir)
       {
         auto& s_trafo = fes->GetMeshAccess()->GetTrafo(pair.other_el, lh);
         s_trafo(pair.other_ip, lh).IntegrationRuleFromPoint
           ([&](const BaseMappedIntegrationRule& s_mir)
            {
              auto& fel = fes->GetFE(pair.master_el, lh);
              auto& s_fel = fes->GetFE(pair.other_el, lh);
              const_cast<BaseMappedIntegrationRule&>(m_mir).SetOtherMIR(&s_mir);
              for(const auto& ce : cb->GetEnergies())
                if(!ce->IsDeformed())
                  ce->ApplyAdd(fel, s_fel, m_mir, elx, ely, lh);
              for(const auto& ci : cb->GetIntegrators())
                if(!ci->IsDeformed())
                  ci->ApplyAdd(fel, s_fel, m_mir, elx, ely, lh);
            });
       });

    m_trafo.AddDeformation(deformation, lh);
    m_trafo(pair.master_ip, lh).IntegrationRuleFromPoint
      ([&](const BaseMappedIntegrationRule& m_mir)
       {
         auto& s_trafo = fes->GetMeshAccess()->GetTrafo(pair.other_el, lh);
         s_trafo.AddDeformation(deformation, lh);
         s_trafo(pair.other_ip, lh).IntegrationRuleFromPoint
           ([&](const BaseMappedIntegrationRule& s_mir)
            {
              auto& fel = fes->GetFE(pair.master_el, lh);
              auto& s_fel = fes->GetFE(pair.other_el, lh);
              const_cast<BaseMappedIntegrationRule&>(m_mir).SetOtherMIR(&s_mir);
              for(const auto& ce : cb->GetEnergies())
                if(ce->IsDeformed())
                  ce->ApplyAdd(fel, s_fel, m_mir, elx, ely, lh);
              for(const auto& ci : cb->GetIntegrators())
                if(ci->IsDeformed())
                  ci->ApplyAdd(fel, s_fel, m_mir, elx, ely, lh);
            });
       });
  }

  template<int DIM>
  void ContactElement<DIM>::CalcLinearizedElementMatrix(FlatVector<double> elx,
                                                        FlatMatrix<double> elmat,
                                                        LocalHeap& lh) const
  {
    auto& m_trafo = fes->GetMeshAccess()->GetTrafo(pair.master_el, lh);
    elmat = 0.;
    m_trafo(pair.master_ip, lh).IntegrationRuleFromPoint
      ([&](const BaseMappedIntegrationRule& m_mir)
       {
         auto& s_trafo = fes->GetMeshAccess()->GetTrafo(pair.other_el, lh);
         s_trafo(pair.other_ip, lh).IntegrationRuleFromPoint
           ([&](const BaseMappedIntegrationRule& s_mir)
            {
              auto& fel = fes->GetFE(pair.master_el, lh);
              auto& s_fel = fes->GetFE(pair.other_el, lh);
              const_cast<BaseMappedIntegrationRule&>(m_mir).SetOtherMIR(&s_mir);
              for(const auto& ce : cb->GetEnergies())
                if(!ce->IsDeformed())
                  ce->CalcLinearizedAdd(fel, s_fel, m_mir, elx, elmat, lh);
              for(const auto& ci : cb->GetIntegrators())
                if(!ci->IsDeformed())
                  ci->CalcLinearizedAdd(fel, s_fel, m_mir, elx, elmat, lh);
            });
       });

    m_trafo.AddDeformation(deformation, lh);
    m_trafo(pair.master_ip, lh).IntegrationRuleFromPoint
      ([&](const BaseMappedIntegrationRule& m_mir)
       {
         auto& s_trafo = fes->GetMeshAccess()->GetTrafo(pair.other_el, lh);
         s_trafo.AddDeformation(deformation, lh);
         s_trafo(pair.other_ip, lh).IntegrationRuleFromPoint
           ([&](const BaseMappedIntegrationRule& s_mir)
            {
              auto& fel = fes->GetFE(pair.master_el, lh);
              auto& s_fel = fes->GetFE(pair.other_el, lh);
              const_cast<BaseMappedIntegrationRule&>(m_mir).SetOtherMIR(&s_mir);
              for(const auto& ce : cb->GetEnergies())
                if(ce->IsDeformed())
                  ce->CalcLinearizedAdd(fel, s_fel, m_mir, elx, elmat, lh);
              for(const auto& ci : cb->GetIntegrators())
                if(ci->IsDeformed())
                  ci->CalcLinearizedAdd(fel, s_fel, m_mir, elx, elmat, lh);
            });
       });

  }

  template class ContactElement<2>;
  template class ContactElement<3>;
} // namespace ngcomp
