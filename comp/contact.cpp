
#include "contact.hpp"

namespace ngcomp
{
  template<int DIM>
  optional<ContactPair<DIM>> T_GapFunction<DIM> :: CreateContactPair(const MappedIntegrationPoint<DIM-1, DIM>& mip1) const
  {
    LocalHeapMem<10000> lh("gapfunction");
    int intorder2 = 10*displacement->GetFESpace()->GetOrder();
    auto & ip1 = mip1.IP();
    auto & trafo1 = mip1.GetTransformation();
    const auto & el1 = ma->GetElement(trafo1.GetElementId());
    auto & trafo1_def = trafo1.AddDeformation(displacement.get(), lh);
    const auto & mip1_def = static_cast<const MappedIntegrationPoint<DIM-1, DIM>&>(trafo1_def(ip1, lh));

    const auto & p1 = mip1_def.GetPoint();

    // find closest point

    double mindist = 1e99;
    IntegrationPoint ip2_min;
    bool intersect = false;
    int el2_min(-1);


    // find all bound-2 elements closer to p1 then h

    netgen::Point<DIM> ngp1;
    for (int j = 0; j < DIM; j++) ngp1(j) = p1(j);
    netgen::Box<DIM> box(ngp1, ngp1);
    box.Increase(h);

    searchtree->GetFirstIntersecting
      (box.PMin(), box.PMax(),
       [&] (int elnr2)
       {
         auto el2 = ma->GetElement(ElementId(BND, elnr2));
         HeapReset hr(lh);

         bool common_vertex = false;
         for (auto s_v : el1.Vertices() )
           for (auto v : el2.Vertices() )
             if(s_v==v)
               common_vertex = true;
         if (common_vertex) return false;

         auto & trafo2 = ma->GetTrafo(el2, lh);
         auto & trafo2_def = trafo2.AddDeformation(displacement.get(), lh);

         IntegrationRule ir2(trafo2.GetElementType(), intorder2);
         MappedIntegrationRule<DIM-1, DIM> mir2(ir2, trafo2, lh);
         MappedIntegrationRule<DIM-1, DIM> mir2_def(ir2, trafo2_def, lh);

         for (auto j : Range(mir2_def))
           {
             const auto & mip2 = mir2_def[j];
             const auto & p2 = mip2.GetPoint();
             double dist = L2Norm(p1-p2);
             if (dist<h && dist < mindist &&
                 InnerProduct(mip1_def.GetNV(), mip2.GetNV()) < 0)
               {
                 mindist = dist;
                 el2_min = el2.Nr();
                 ip2_min = mip2.IP();
                 intersect = true;
               }
           }
         return false;
       });
    if(intersect)
      return ContactPair<DIM>{el1, ElementId(BND,el2_min),
          ip1, ip2_min};
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
        auto & mask = slave.Mask();
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
        auto & mask = slave.Mask();
        if (!mask.Test(el2.GetIndex())) continue;

        auto & trafo2 = ma->GetTrafo (el2, lh);
        auto & trafo2_def = trafo2.AddDeformation(displacement.get(), lh);

        IntegrationRule ir2(trafo2.GetElementType(), intorder2);
        MappedIntegrationRule<DIM-1, DIM> mir2(ir2, trafo2, lh);
        MappedIntegrationRule<DIM-1, DIM> mir2_def(ir2, trafo2_def, lh);

        for (auto & mip : mir2_def)
          {
            netgen::Point<DIM> p;
            for (int j = 0; j < DIM; j++)
              p(j) = mip.GetPoint()(j);
            bbox.Add(p);
          }

        searchtree->Insert(bbox, el2.Nr());
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

    int intorder2 = 10*displacement->GetFESpace()->GetOrder();

    auto & trafo1_def = trafo1.AddDeformation(displacement.get(), lh);

    auto & ip1 = ip.IP();
    Vec<DIM> p1;
    trafo1_def.CalcPoint(ip1, p1);

    double mindist = 1e99;
    result = 0;

    // find all bound-2 elements closer to p1 than h
    netgen::Point<DIM> ngp1;
    for (int j = 0; j < DIM; j++)
      ngp1(j) = p1(j);
    netgen::Box<DIM> box(ngp1, ngp1);
    box.Increase(h);

    searchtree->GetFirstIntersecting
      (box.PMin(), box.PMax(),
       [&] (int elnr2)
       {
         auto el2 = ma->GetElement( ElementId (BND, elnr2) );
         HeapReset hr(lh);

         bool common_vertex = false;
         for (auto s_v : el1.Vertices() )
           for (auto v : el2.Vertices() )
             if(s_v==v)
               common_vertex = true;
         if (common_vertex) return false;
         auto & trafo2 = ma->GetTrafo (el2, lh);
         auto & trafo2_def = trafo2.AddDeformation(displacement.get(), lh);

         IntegrationRule ir2(trafo2.GetElementType(), intorder2);
         MappedIntegrationRule<DIM-1, DIM> mir2(ir2, trafo2, lh);
         MappedIntegrationRule<DIM-1, DIM> mir2_def(ir2, trafo2_def, lh);

         for (auto j : Range(mir2_def))
           {
             const auto & mip2 = mir2_def[j];
             const auto & p2 = mip2.GetPoint();
             double dist = L2Norm(p1-p2);
             if (dist<h && dist < mindist)
               {
                 mindist = dist;
                 result = p2-p1;
               }
           }
         return false;
       });
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
    if(!displacement)
      {
        values = static_cast<const MappedIntegrationPoint<DIM-1, DIM>&>(ip).GetNV();
        return;
      }
    LocalHeapMem<10000> lh("deformednormal");
    auto& trafo_def = ip.GetTransformation().AddDeformation(displacement.get(), lh);
    auto& mip_def = static_cast<const MappedIntegrationPoint<DIM-1, DIM>&>(trafo_def(ip.IP(), lh));
    values = mip_def.GetNV();
  }

  template class DisplacedNormal<2>;
  template class DisplacedNormal<3>;

  ContactEnergy::ContactEnergy(shared_ptr<CoefficientFunction> _cf,
                               shared_ptr<FESpace> _fes)
    : cf(_cf), fes(_fes)
  {
    cf->TraverseTree
      ([&](CoefficientFunction& nodecf)
       {
         auto proxy = dynamic_cast<ProxyFunction*>(&nodecf);
         if (proxy && !proxy->IsTestFunction() && !trial_proxies.Contains(proxy))
           trial_proxies.Append(proxy);
       });
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
                                       shared_ptr<FESpace> _fes)
    : cf(_cf), fes(_fes)
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

  ContactBoundary::ContactBoundary(shared_ptr<FESpace> _fes,
                                   Region _master, Region _slave)
    : master(_master), slave(_slave), fes(_fes)
  {
    auto mesh = fes->GetMeshAccess();
    if(mesh->GetDimension() == 2)
      {
        gap = make_shared<T_GapFunction<2>>(mesh, master, slave);
        normal = make_shared<DisplacedNormal<2>>();
      }
    else
      {
        gap = make_shared<T_GapFunction<3>>(mesh, master, slave);
        normal = make_shared<DisplacedNormal<3>>();
      }
  }

  void ContactBoundary::AddEnergy(shared_ptr<CoefficientFunction> form)
  {
    energies.Append(make_shared<ContactEnergy>(form, fes));
  }

  void ContactBoundary::AddIntegrator(shared_ptr<CoefficientFunction> form)
  {
    integrators.Append(make_shared<ContactIntegrator>(form, fes));
  }

  void ContactBoundary::Update(shared_ptr<GridFunction> displacement_,
                               shared_ptr<BilinearForm> bf,
                               int intorder, double h)
  {
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
                          auto pair = tgap->CreateContactPair(mip);
                          if(pair.has_value())
                            {
                              lock_guard<mutex> guard(add_mutex);
                              bf->AddSpecialElement(make_unique<ContactElement<DIM>>(*pair, this));
                            }
                        }
                    });
               }
           });
      }
  }

  template<int DIM>
  ContactElement<DIM>::ContactElement(const ContactPair<DIM>& _pair,
                                      ContactBoundary* _cb)
    : pair(_pair), cb(_cb), fes(_cb->GetFESpace().get())
  {
  }

  template<int DIM>
  void ContactElement<DIM>::GetDofNrs(Array<DofId>& dnums) const
  {
    fes->GetDofNrs(pair.master_el, dnums);
    Array<DofId> s_dofs;
    fes->GetDofNrs(pair.slave_el, s_dofs);
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
         auto& s_trafo = fes->GetMeshAccess()->GetTrafo(pair.slave_el, lh);
         s_trafo(pair.slave_ip, lh).IntegrationRuleFromPoint
           ([&](const BaseMappedIntegrationRule& s_mir)
            {
              auto& fel = fes->GetFE(pair.master_el, lh);
              auto& s_fel = fes->GetFE(pair.slave_el, lh);
              const_cast<BaseMappedIntegrationRule&>(m_mir).SetOtherMIR(&s_mir);
              for(const auto& ce : cb->GetEnergies())
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
         auto& s_trafo = fes->GetMeshAccess()->GetTrafo(pair.slave_el, lh);
         s_trafo(pair.slave_ip, lh).IntegrationRuleFromPoint
           ([&](const BaseMappedIntegrationRule& s_mir)
            {
              auto& fel = fes->GetFE(pair.master_el, lh);
              auto& s_fel = fes->GetFE(pair.slave_el, lh);
              const_cast<BaseMappedIntegrationRule&>(m_mir).SetOtherMIR(&s_mir);
              for(const auto& ce : cb->GetEnergies())
                ce->ApplyAdd(fel, s_fel, m_mir, elx, ely, lh);
              for(const auto& ci : cb->GetIntegrators())
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
         auto& s_trafo = fes->GetMeshAccess()->GetTrafo(pair.slave_el, lh);
         s_trafo(pair.slave_ip, lh).IntegrationRuleFromPoint
           ([&](const BaseMappedIntegrationRule& s_mir)
            {
              auto& fel = fes->GetFE(pair.master_el, lh);
              auto& s_fel = fes->GetFE(pair.slave_el, lh);
              const_cast<BaseMappedIntegrationRule&>(m_mir).SetOtherMIR(&s_mir);
              for(const auto& ce : cb->GetEnergies())
                ce->CalcLinearizedAdd(fel, s_fel, m_mir, elx, elmat, lh);
              for(const auto& ci : cb->GetIntegrators())
                ci->CalcLinearizedAdd(fel, s_fel, m_mir, elx, elmat, lh);
            });
       });
  }

  template class ContactElement<2>;
  template class ContactElement<3>;
} // namespace ngcomp
