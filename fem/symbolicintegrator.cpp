/*********************************************************************/
/* File:   symbolicintegrator.cpp                                    */
/* Author: Joachim Schoeberl                                         */
/* Date:   August 2015                                               */
/*********************************************************************/
/* 
   Symbolic integrators
*/

#include <fem.hpp>
#include <map>
#include "avector.hpp"
  
namespace ngfem
{
  
  class ProxyUserData
  {
    // map<const ProxyFunction*, FlatMatrix<double>> remember;
    FlatArray<const ProxyFunction*> remember_first;
    FlatArray<FlatMatrix<double>> remember_second;
    FlatArray<AFlatMatrix<double>> remember_asecond;
  public:
    class ProxyFunction * testfunction = nullptr;
    int test_comp;
    class ProxyFunction * trialfunction = nullptr;
    int trial_comp;
    
    const FiniteElement * fel = nullptr;
    const FlatVector<double> * elx;
    LocalHeap * lh;

    ProxyUserData ()
      : remember_first(0,nullptr), remember_second(0,nullptr), remember_asecond(0,nullptr) { ; }
    ProxyUserData (int ntrial, LocalHeap & lh)
      : remember_first(ntrial, lh), remember_second(ntrial, lh), remember_asecond(ntrial, lh)
    { remember_first = nullptr; }
    
    void AssignMemory (const ProxyFunction * proxy, int h, int w, LocalHeap & lh)
    {
      // remember[proxy] = Matrix<> (h, w);
      /*
      remember.emplace(piecewise_construct,
                       forward_as_tuple(proxy),
                       forward_as_tuple(h,w,lh));
      */
      for (int i = 0; i < remember_first.Size(); i++)
        {
          if (remember_first[i] == nullptr)
            {
              remember_first[i] = proxy;
              // remember_second[i].AssignMemory (h,w,lh);
              new (&remember_second[i]) FlatMatrix<> (h, w, lh);
              new (&remember_asecond[i]) AFlatMatrix<double> (w, h, lh);
              return;
            }
        }
      throw Exception ("no space for userdata - memory available");
    }
    bool HasMemory (const ProxyFunction * proxy) const
    {
      // return remember.count (proxy);
      return remember_first.Contains(proxy);
    }
    FlatMatrix<> GetMemory (const ProxyFunction * proxy) const
    {
      // return remember.find(proxy) -> second;
      return remember_second[remember_first.Pos(proxy)];
    }
    AFlatMatrix<double> GetAMemory (const ProxyFunction * proxy) const
    {
      return remember_asecond[remember_first.Pos(proxy)];
    }
  };

  Array<int> ProxyFunction ::
  Dimensions() const
  {
    int dim = evaluator->Dim();
    int blockdim = evaluator->BlockDim();
    if (blockdim == 1)
      return Array<int> ({dim});
    else
      return Array<int> ({dim/blockdim, blockdim});
  }
  
  
  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationPoint & mip,
            FlatVector<> result) const
  {
    ProxyUserData * ud = (ProxyUserData*)mip.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");

    if (!testfunction && ud->fel)
      {
        evaluator->Apply (*ud->fel, mip, *ud->elx, result, *ud->lh);
        return;
      }

    result = 0;
    if (ud->testfunction == this)
      result (ud->test_comp) = 1;
    if (ud->trialfunction == this)
      result (ud->trial_comp) = 1;
  }

  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationRule & mir,
            FlatMatrix<> result) const
  {
    // static Timer t("ProxyFunction :: Evaluate", 2);
    // RegionTimer reg(t);
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");
    
    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          result = ud->GetMemory (this);
        else
          evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
        return;
      }

    result = 0;
    if (ud->testfunction == this)
      result.Col(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      result.Col(ud->trial_comp) = 1;
  }

  void ProxyFunction ::
  Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
            AFlatMatrix<double> result) const
  {
    // static Timer t("ProxyFunction::EvalSIMD"); RegionTimer reg(t);
    
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");
    
    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          // result = Trans(ud->GetMemory (this));
          result = ud->GetAMemory (this);
        else
          throw Exception ("ProxyFunction::Evaluate(simd) : need x-values");
          // evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
        return;
      }

    result = 0;
    if (ud->testfunction == this)
      result.Row(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      result.Row(ud->trial_comp) = 1;
  }


  
  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationPoint & ip,
            FlatVector<Complex> result) const
  {
    Vector<> result_double(result.Size());
    Evaluate (ip, result_double);
    result = result_double;
  }

  void ProxyFunction ::
  EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                 FlatMatrix<> result,
                 FlatMatrix<> deriv) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    deriv = 0;
    result = 0;

    static Timer t("ProxyFunction EvaluateDeriv");
    t.Start();
    if (!testfunction && ud->fel)
      {
        // evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
        /*
        if (ud->remember.count (const_cast<ProxyFunction*>(this)))
          result = ud->remember[const_cast<ProxyFunction*>(this)];
        */
        if (ud->HasMemory(this))
          result = ud->GetMemory (this);          
        else
          evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
      }
    t.Stop();
    
    if (ud->testfunction == this)
      result.Col(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv.Col(ud->trial_comp) = 1;
  }


  void ProxyFunction ::
  EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                  FlatMatrix<> result,
                  FlatMatrix<> deriv,
                  FlatMatrix<> dderiv) const
  {
    static Timer t("ProxyFunction :: Evaluate DDeriv", 2);
    static Timer t2("ProxyFunction :: Evaluate DDeriv, calc only", 2);
    RegionTimer reg(t);
    
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    result = 0;
    deriv = 0;
    dderiv = 0;

    if (!testfunction && ud->fel)
      {
        RegionTimer reg2(t2);
        // if (ud->remember.count (const_cast<ProxyFunction*>(this)))
        // result = ud->remember[const_cast<ProxyFunction*>(this)];
        if (ud->HasMemory (this))
          result = ud->GetMemory (this);
        else
          evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
      }
    if (ud->testfunction == this)
      deriv.Col(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv.Col(ud->trial_comp) = 1;
  }

  void ProxyFunction ::  
  NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const
  {
    if (!testfunction && ud.fel)
      {
        nonzero = true;
        return;
      }

    nonzero = false;
    if (ud.testfunction == this)
      nonzero(ud.test_comp) = true;
    if (ud.trialfunction == this)
      nonzero(ud.trial_comp) = true;
  }


  SymbolicLinearFormIntegrator ::
  SymbolicLinearFormIntegrator(shared_ptr<CoefficientFunction> acf, VorB avb)
    : cf(acf), vb(avb)
  {
    if (cf->Dimension() != 1)
      throw Exception ("SymbolicLFI needs scalar-valued CoefficientFunction");
    cf->TraverseTree
      ([&] (CoefficientFunction & nodecf)
       {
         auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
         if (proxy && !proxies.Contains(proxy))
           proxies.Append (proxy);
       });
  }

  
  void 
  SymbolicLinearFormIntegrator ::
  CalcElementVector (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatVector<double> elvec,
                     LocalHeap & lh) const
  {
    T_CalcElementVector (fel, trafo, elvec, lh);
  }
  
  void 
  SymbolicLinearFormIntegrator ::
  CalcElementVector (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatVector<Complex> elvec,
                     LocalHeap & lh) const
  {
    T_CalcElementVector (fel, trafo, elvec, lh);
  }
  
  template <typename SCAL> 
  void SymbolicLinearFormIntegrator ::
  T_CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatVector<SCAL> elvec,
                       LocalHeap & lh) const
  {
    static Timer t("symbolicLFI - CalcElementVector", 2);
    // static Timer td("symboliLFI - CalcElementVector dvecs", 2);
    // static Timer tb("symboliLFI - CalcElementVector diffops", 2);
    RegionTimer reg(t);
    
    HeapReset hr(lh);
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    
    FlatVector<SCAL> elvec1(elvec.Size(), lh);
    
    FlatMatrix<SCAL> values(ir.Size(), 1, lh);
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    
    elvec = 0;
    for (auto proxy : proxies)
      {
        // td.Start();
        FlatMatrix<SCAL> proxyvalues(ir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            
            cf -> Evaluate (mir, values);
            for (int i = 0; i < mir.Size(); i++)
              proxyvalues(i,k) = mir[i].GetWeight() * values(i,0);
          }
        // td.Stop();
        // tb.Start();
        proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
        // tb.Stop();
        elvec += elvec1;
      }
  }
  
  


  SymbolicBilinearFormIntegrator ::
  SymbolicBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb,
                                  bool aelement_boundary)
    : cf(acf), vb(avb), element_boundary(aelement_boundary)
  {
    if (cf->Dimension() != 1)
        throw Exception ("SymblicBFI needs scalar-valued CoefficientFunction");
    trial_cum.Append(0);
    test_cum.Append(0);    
    cf->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
          if (proxy) 
            {
              if (proxy->IsTestFunction())
                {
                  if (!test_proxies.Contains(proxy))
                    {
                      test_proxies.Append (proxy);
                      test_cum.Append(test_cum.Last()+proxy->Dimension());
                    }
                }
              else
                {                                         
                  if (!trial_proxies.Contains(proxy))
                    {
                      trial_proxies.Append (proxy);
                      trial_cum.Append(trial_cum.Last()+proxy->Dimension());
                    }
                }
            }
        });
    cout << IM(3) << "num test_proxies " << test_proxies.Size() << endl;
    cout << IM(3) << "num trial_proxies " << trial_proxies.Size() << endl;
    cout << IM(3) << "cumulated test_proxy dims  " << test_cum << endl;
    cout << IM(3) << "cumulated trial_proxy dims " << trial_cum << endl;

    elementwise_constant = cf -> ElementwiseConstant();
    cout << IM(3) << "element-wise constant = " << elementwise_constant << endl;

    // find non-zeros
    int cnttest = 0, cnttrial = 0;
    for (auto proxy : trial_proxies)
      cnttrial += proxy->Dimension();
    for (auto proxy : test_proxies)
      cnttest += proxy->Dimension();
    nonzeros = Matrix<bool>(cnttest, cnttrial);

    ProxyUserData ud;
    Vector<bool> nzvec(1);
    int k = 0;
    for (int k1 : test_proxies.Range())
      for (int k2 : Range(0,test_proxies[k1]->Dimension()))
        {
          int l = 0;
          for (int l1 : trial_proxies.Range())
            for (int l2 : Range(0,trial_proxies[l1]->Dimension()))
              {
                ud.trialfunction = trial_proxies[l1];
                ud.trial_comp = l2;
                ud.testfunction = test_proxies[k1];
                ud.test_comp = k2;
                cf -> NonZeroPattern (ud, nzvec);
                nonzeros(k,l) = nzvec(0);
                l++;
              }
          k++;
        }
    cout << IM(3) << "nonzeros: " << endl << nonzeros << endl;
  }
  
  
  void 
  SymbolicBilinearFormIntegrator ::
  CalcElementMatrix (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    T_CalcElementMatrix<double> (fel, trafo, elmat, lh);
  }
  
  void 
  SymbolicBilinearFormIntegrator ::
  CalcElementMatrix (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatMatrix<Complex> elmat,
                     LocalHeap & lh) const
  {
    if (fel.ComplexShapes())
      T_CalcElementMatrix<Complex,Complex> (fel, trafo, elmat, lh);
    else
      T_CalcElementMatrix<Complex,double> (fel, trafo, elmat, lh);
  }
  
  
  template <typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatMatrix<SCAL> elmat,
                       LocalHeap & lh) const
    
  { 
    static Timer t("symbolicBFI - CalcElementMatrix", 2);
    // static Timer tstart("symboliBFI - CalcElementMatrix startup", 2);
    // static Timer tstart1("symboliBFI - CalcElementMatrix startup 1", 2);
    // static Timer tmain("symboliBFI - CalcElementMatrix main", 2);
    static Timer td("symboliBFI - CalcElementMatrix dmats", 2);
    static Timer tb("symboliBFI - CalcElementMatrix diffops", 2);
    static Timer tdb("symboliBFI - CalcElementMatrix D * B", 2);
    static Timer tlapack("symboliBFI - CalcElementMatrix lapack", 2);

    // tstart.Start();
    
    if (element_boundary)
      {
        switch (trafo.SpaceDim())
          {
          case 1:
            T_CalcElementMatrixEB<1,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          case 2:
            T_CalcElementMatrixEB<2,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          case 3:
            T_CalcElementMatrixEB<3,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          default:
            throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
          }
      }
    
    RegionTimer reg(t);


    int trial_difforder = 99, test_difforder = 99;
    for (auto proxy : trial_proxies)
      trial_difforder = min2(trial_difforder, proxy->Evaluator()->DiffOrder());
    for (auto proxy : test_proxies)
      test_difforder = min2(test_difforder, proxy->Evaluator()->DiffOrder());

    int intorder = 2*fel.Order();
    auto et = trafo.GetElementType();
    if (et == ET_TRIG || et == ET_TET)
      intorder -= test_difforder+trial_difforder;
    IntegrationRule ir(trafo.GetElementType(), intorder);
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    
    elmat = 0;

    // tstart.Stop();

    int k1 = 0;
    for (auto proxy1 : trial_proxies)
      {
        int l1 = 0;
        for (auto proxy2 : test_proxies)
          {
            bool is_diagonal = proxy1->Dimension() == proxy2->Dimension();
            bool is_nonzero = false;

            for (int k = 0; k < proxy1->Dimension(); k++)
              for (int l = 0; l < proxy2->Dimension(); l++)
                if (nonzeros(l1+l, k1+k))
                  {
                    if (k != l) is_diagonal = false;
                    is_nonzero = true;
                  }


            if (is_nonzero)
              {
                HeapReset hr(lh);
                bool samediffop = *(proxy1->Evaluator()) == *(proxy2->Evaluator());
                td.Start();
                FlatTensor<3,SCAL> proxyvalues(lh, mir.Size(), proxy1->Dimension(), proxy2->Dimension());
                FlatVector<SCAL> diagproxyvalues(mir.Size()*proxy1->Dimension(), lh);
                FlatMatrix<SCAL> val(mir.Size(), 1, lh);
                
                
                if (!is_diagonal)
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    for (int l = 0; l < proxy2->Dimension(); l++)
                      {
                        if (nonzeros(l1+l, k1+k))
                          {
                            if (k != l) is_diagonal = false;
                            is_nonzero = true;
                            ud.trialfunction = proxy1;
                            ud.trial_comp = k;
                            ud.testfunction = proxy2;
                            ud.test_comp = l;
                            
                            cf -> Evaluate (mir, val);
                            proxyvalues(STAR,k,l) = val.Col(0);
                          }
                        else
                          proxyvalues(STAR,k,l) = 0.0;
                      }
                else
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    {
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = k;

                      if (!elementwise_constant)
                        {
                          cf -> Evaluate (mir, val);
                          diagproxyvalues.Slice(k, proxy1->Dimension()) = val.Col(0);
                        }
                      else
                        {
                          cf -> Evaluate (mir[0], val.Row(0));
                          diagproxyvalues.Slice(k, proxy1->Dimension()) = val(0,0);
                        }
                    }
            
                td.Stop();

                if (!is_diagonal)
                  for (int i = 0; i < mir.Size(); i++)
                    proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
                else
                  for (int i = 0; i < mir.Size(); i++)
                    diagproxyvalues.Range(proxy1->Dimension()*IntRange(i,i+1)) *= mir[i].GetWeight();
                  
                IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                SliceMatrix<SCAL> part_elmat = elmat.Rows(r2).Cols(r1);
                FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

                
                enum { BS = 16 };
                for (int i = 0; i < mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    int bs = min2(int(BS), mir.Size()-i);
                    
                    AFlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);
                    AFlatMatrix<SCAL> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                    AFlatMatrix<SCAL_SHAPES> bbmat2 = samediffop ?
                      bbmat1 : AFlatMatrix<SCAL_SHAPES>(elmat.Height(), bs*proxy2->Dimension(), lh);

                    tb.Start();
                    BaseMappedIntegrationRule & bmir = mir.Range(i, i+bs, lh);
                    proxy1->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat1), lh);
                    if (!samediffop)
                      proxy2->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat2), lh);
                    tb.Stop();

                    tdb.Start();
                    if (is_diagonal)
                      {
                        AFlatVector<SCAL> diagd(bs*proxy1->Dimension(), lh);
                        diagd = diagproxyvalues.Range(i*proxy1->Dimension(),
                                                      (i+bs)*proxy1->Dimension());
                        /*
                        for (int i = 0; i < diagd.Size(); i++)
                          bdbmat1.Col(i) = diagd(i) * bbmat1.Col(i);
                        */
                        MultMatDiagMat(bbmat1, diagd, bdbmat1);
                        // tdb.AddFlops (bbmat1.Height()*bbmat1.Width());
                      }
                    else
                      {
                        for (int j = 0; j < bs; j++)
                          {
                            int ii = i+j;
                            IntRange r1 = proxy1->Dimension() * IntRange(j,j+1);
                            IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                            // bdbmat1.Cols(r2) = bbmat1.Cols(r1) * proxyvalues(ii,STAR,STAR);
                            MultMatMat (bbmat1.Cols(r1), proxyvalues(ii,STAR,STAR), bdbmat1.Cols(r2));
                          }
                        // tdb.AddFlops (proxy1->Dimension()*proxy2->Dimension()*bs*bbmat1.Height());
                      }
                    tdb.Stop();
                    
                    tlapack.Start();
                    // elmat.Rows(r2).Cols(r1) += bbmat2.Rows(r2) * Trans(bdbmat1.Rows(r1));
                    // AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), elmat.Rows(r2).Cols(r1));

                    if (samediffop && is_diagonal)
                      AddABtSym (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);
                    else
                      AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);

                    tlapack.Stop();
                    // tlapack.AddFlops (r2.Size()*r1.Size()*bdbmat1.Width());
                  }

                if (samediffop && is_diagonal)
                  for (int i = 0; i < part_elmat.Height(); i++)
                    for (int j = i+1; j < part_elmat.Width(); j++)
                      part_elmat(i,j) = part_elmat(j,i);
              }
            
            l1 += proxy2->Dimension();  
          }
        k1 += proxy1->Dimension();
      }
  }
  
  

  template <int D, typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcElementMatrixEB (const FiniteElement & fel,
                           const ElementTransformation & trafo, 
                           FlatMatrix<SCAL> elmat,
                           LocalHeap & lh) const
      
    {
      static Timer t("symbolicBFI - CalcElementMatrix EB", 2);
      static Timer td("symbolicBFI - CalcElementMatrix EB - dmats", 2);
      static Timer tb("symbolicBFI - CalcElementMatrix EB - bmats", 2);
      static Timer tmult("symbolicBFI - CalcElementMatrix EB - mult", 2);
      RegionTimer reg(t);
      
      elmat = 0;

      auto eltype = trafo.GetElementType();
      int nfacet = ElementTopology::GetNFacets(eltype);

      Facet2ElementTrafo transform(eltype); 

      for (int k = 0; k < nfacet; k++)
        {
          HeapReset hr(lh);
          ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
        
          IntegrationRule ir_facet(etfacet, 2*fel.Order());
          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
          const BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
          
          ProxyUserData ud;
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          FlatVector<> measure(mir.Size(), lh);
          
          for (int i = 0; i < mir.Size(); i++)
            {
              double len;
              if (!trafo.Boundary())
                {
                  FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);
                  Vec<D> normal_ref = normals[k];
                  auto & mip = static_cast<const MappedIntegrationPoint<D,D>&> (mir[i]);
                  Mat<D> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
                  len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  
                  const_cast<MappedIntegrationPoint<D,D>&> (mip).SetNV(normal);
                }
              else
                {
                  if (D != 3)
                    throw Exception ("element boundary for surface elements is only possible in 3D");
                  FlatVector< Vec<D-1> > normals = ElementTopology::GetNormals<D-1>(eltype);
                  Vec<D-1> normal_ref = normals[k];

                  auto & mip = static_cast<const MappedIntegrationPoint<2,3>&> (mir[i]);
                  Mat<2,3> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<3> normal = det * Trans (inv_jac) * normal_ref;       
                  len = L2Norm (normal);    // that's the surface measure
                  normal /= len;                   // normal vector on physical element
                  Vec<3> tang = Cross(normal, mip.GetNV());
                  const_cast<MappedIntegrationPoint<2,3>&> (mip).SetTV(tang);
                }
              measure(i) = len;
            }

          
          for (auto proxy1 : trial_proxies)
            for (auto proxy2 : test_proxies)
              {
                HeapReset hr(lh);
                FlatTensor<3,SCAL> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
                FlatMatrix<SCAL> val(mir.Size(), 1, lh);
                
                td.Start();
                for (int k = 0; k < proxy1->Dimension(); k++)
                  for (int l = 0; l < proxy2->Dimension(); l++)
                    {
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = l;
                      
                      cf->Evaluate (mir, val);
                      for (int i = 0; i < mir.Size(); i++)
                        val(i) *= ir_facet[i].Weight() * measure(i);
                      proxyvalues(STAR,l,k) = val.Col(0);
                    }
                td.Stop();

                for (int i = 0; i < mir.Size(); i++)
                  {
                    tb.Start();
                    FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                    FlatMatrix<SCAL,ColMajor> dbmat1(proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
                    
                    proxy1->Evaluator()->CalcMatrix(fel, mir[i], bmat1, lh);
                    proxy2->Evaluator()->CalcMatrix(fel, mir[i], bmat2, lh);
                    tb.Stop();
                    tmult.Start();
                    IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                    IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                    
                    dbmat1 = proxyvalues(i,STAR,STAR) * bmat1;                    
                    elmat.Rows(r2).Cols(r1) += Trans (bmat2.Cols(r2)) * dbmat1.Cols(r1);
                    tmult.Stop();
                  }
              }
        }
    }

  
  void 
  SymbolicBilinearFormIntegrator ::
  CalcLinearizedElementMatrix (const FiniteElement & fel,
                               const ElementTransformation & trafo, 
                               FlatVector<double> elveclin,
                               FlatMatrix<double> elmat,
                               LocalHeap & lh) const
  {
    /*
      CalcElementMatrix(fel, trafo, elmat, lh);
      return;
    */

      
    if (element_boundary)
      {
        switch (trafo.SpaceDim())
          {
          case 1:
            T_CalcLinearizedElementMatrixEB<1,double,double> (fel, trafo, elveclin, elmat, lh);
            return;
          case 2:
            T_CalcLinearizedElementMatrixEB<2,double,double> (fel, trafo, elveclin, elmat, lh);            
            return;
          case 3:
            T_CalcLinearizedElementMatrixEB<3,double,double> (fel, trafo, elveclin, elmat, lh);
            return;
          default:
            throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
          }
      }

    
    static Timer t("symbolicbfi - calclinearized", 2);
    static Timer td("symbolicbfi - calclinearized dmats", 2);
    RegionTimer reg(t);
    
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elveclin;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      {
        // ud.remember[proxy] = Matrix<> (ir.Size(), proxy->Dimension());
        // proxy->Evaluator()->Apply(fel, mir, elveclin, ud.remember[proxy], lh);
        ud.AssignMemory (proxy, ir.Size(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetMemory(proxy), lh);
      }
    
    FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh);
    
    elmat = 0;
               
    
    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          td.Start(); 
          FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              if (nonzeros(test_cum[l1]+l, trial_cum[k1]+k))
                {
                  ud.trialfunction = proxy1;
                  ud.trial_comp = k;
                  ud.testfunction = proxy2;
                  ud.test_comp = l;
                  
                  cf -> EvaluateDeriv (mir, val, deriv);
                  proxyvalues(STAR,l,k) = deriv.Col(0);
                }
              else
                proxyvalues(STAR,l,k) = 0;
          td.Stop();

          for (int i = 0; i < mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();

          t.AddFlops (double (mir.Size()) * proxy1->Dimension()*elmat.Width()*elmat.Height());

          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

          enum { BS = 16 };
          for (int i = 0; i < mir.Size(); i+=BS)
            {
              int rest = min2(int(BS), mir.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
              elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }




          
        }
  }



  template <int D, typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcLinearizedElementMatrixEB (const FiniteElement & fel,
                                   const ElementTransformation & trafo, 
                                   FlatVector<double> elveclin,
                                   FlatMatrix<double> elmat,
                                   LocalHeap & lh) const
  {
    static Timer t("symbolicbfi - calclinearized EB", 2);
    static Timer td("symbolicbfi - calclinearized EB dmats", 2);
    RegionTimer reg(t);
    
    // IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    // BaseMappedIntegrationRule & mir = trafo(ir, lh);
    /*
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elveclin;
    ud.lh = &lh;
    */
    elmat = 0;



      auto eltype = trafo.GetElementType();
      int nfacet = ElementTopology::GetNFacets(eltype);

      Facet2ElementTrafo transform(eltype); 

      for (int k = 0; k < nfacet; k++)
        {
          HeapReset hr(lh);
          ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
        
          IntegrationRule ir_facet(etfacet, 2*fel.Order());
          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
          const BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
          
          ProxyUserData ud(trial_proxies.Size(), lh);          
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          FlatVector<> measure(mir.Size(), lh);
          
          for (int i = 0; i < mir.Size(); i++)
            {
              double len;
              if (!trafo.Boundary())
                {
                  FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);
                  Vec<D> normal_ref = normals[k];
                  auto & mip = static_cast<const MappedIntegrationPoint<D,D>&> (mir[i]);
                  Mat<D> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
                  len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  
                  const_cast<MappedIntegrationPoint<D,D>&> (mip).SetNV(normal);
                }
              else
                {
                  if (D != 3)
                    throw Exception ("element boundary for surface elements is only possible in 3D");
                  FlatVector< Vec<D-1> > normals = ElementTopology::GetNormals<D-1>(eltype);
                  Vec<D-1> normal_ref = normals[k];

                  auto & mip = static_cast<const MappedIntegrationPoint<2,3>&> (mir[i]);
                  Mat<2,3> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<3> normal = det * Trans (inv_jac) * normal_ref;       
                  len = L2Norm (normal);    // that's the surface measure
                  normal /= len;                   // normal vector on physical element
                  Vec<3> tang = Cross(normal, mip.GetNV());
                  const_cast<MappedIntegrationPoint<2,3>&> (mip).SetTV(tang);
                }
              measure(i) = len;
            }


          for (ProxyFunction * proxy : trial_proxies)
            {
              // ud.remember[proxy] = Matrix<> (mir.Size(), proxy->Dimension());
              // proxy->Evaluator()->Apply(fel, mir, elveclin, ud.remember[proxy], lh);
              ud.AssignMemory (proxy, mir.Size(), proxy->Dimension(), lh);
              proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetMemory(proxy), lh);
            }
    
    
          FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh);
          for (int k1 : Range(trial_proxies))
            for (int l1 : Range(test_proxies))
              {
                HeapReset hr(lh);
                auto proxy1 = trial_proxies[k1];
                auto proxy2 = test_proxies[l1];
                td.Start();
                FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
                
                for (int k = 0; k < proxy1->Dimension(); k++)
                  for (int l = 0; l < proxy2->Dimension(); l++)
                    if (nonzeros(test_cum[l1]+l, trial_cum[k1]+k))
                      {
                        ud.trialfunction = proxy1;
                        ud.trial_comp = k;
                        ud.testfunction = proxy2;
                        ud.test_comp = l;
                        
                        cf -> EvaluateDeriv (mir, val, deriv);
                        proxyvalues(STAR,l,k) = deriv.Col(0);
                      }
                    else
                      proxyvalues(STAR,l,k) = 0.0;
                        
                td.Stop();

                for (int i = 0; i < mir.Size(); i++)
                  proxyvalues(i,STAR,STAR) *= ir_facet[i].Weight() * measure(i);
                
                t.AddFlops (double (mir.Size()) * proxy1->Dimension()*elmat.Width()*elmat.Height());
                
                FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
                
                enum { BS = 16 };
                for (int i = 0; i < mir.Size(); i+=BS)
                  {
                    int rest = min2(int(BS), mir.Size()-i);
                    HeapReset hr(lh);
                    FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);
                    
                    for (int j = 0; j < rest; j++)
                      {
                        int ii = i+j;
                        IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                        proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                        proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                        bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                        bbmat2.Rows(r2) = bmat2;
                      }
                    
                    IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                    IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                    elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
                  }


                /*
                int i = 0;
                
                enum { BS = 16 };
                for ( ; i+BS <= mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    FlatMatrix<double,ColMajor> bdbmat1(BS*proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<double,ColMajor> bbmat2(BS*proxy2->Dimension(), elmat.Height(), lh);
                    
                    for (int j = 0; j < BS; j++)
                      {
                        int ii = i+j;
                        IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                        proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                        proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                        bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                        bbmat2.Rows(r2) = bmat2;
                      }
                    elmat += Trans (bbmat2) * bdbmat1 | Lapack;
                  }
                
                
                if (i < mir.Size())
                  {
                    HeapReset hr(lh);
                    int rest = mir.Size()-i;
                    FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);
                    
                    for (int j = 0; j < rest; j++)
                      {
                        int ii = i+j;
                        IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                        proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                        proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                        bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                        bbmat2.Rows(r2) = bmat2;
                      }
                    elmat += Trans (bbmat2) * bdbmat1 | Lapack;
                  }
                */
              }
        }
  }









  
  
  void
  SymbolicBilinearFormIntegrator :: ApplyElementMatrix (const FiniteElement & fel, 
                                                        const ElementTransformation & trafo, 
                                                        const FlatVector<double> elx, 
                                                        FlatVector<double> ely,
                                                        void * precomputed,
                                                        LocalHeap & lh) const
  {
    /*
    if (element_boundary)
      {
        FlatMatrix<> elmat(elx.Size(), lh);
        CalcElementMatrix(fel, trafo, elmat, lh);
        ely = elmat * elx;
        return;
      }
    */
    
    if (element_boundary)
      {
        switch (trafo.SpaceDim())
          {
          case 1:
            T_ApplyElementMatrixEB<1,double,double> (fel, trafo, elx, ely, precomputed, lh);
            return;
          case 2:
            T_ApplyElementMatrixEB<2,double,double> (fel, trafo, elx, ely, precomputed, lh);            
            return;
          case 3:
            T_ApplyElementMatrixEB<3,double,double> (fel, trafo, elx, ely, precomputed, lh);            
            return;
          default:
            throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
          }
      }


    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;

    HeapReset hr(lh);
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
      
    ely = 0;
    FlatVector<> ely1(ely.Size(), lh);

    FlatMatrix<> val(mir.Size(), 1,lh);
    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir, val);
            proxyvalues.Col(k) = val.Col(0);
          }

        for (int i = 0; i < mir.Size(); i++)
          proxyvalues.Row(i) *= mir[i].GetWeight();
        
        proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, ely1, lh);
        ely += ely1;
      }
  }


  template <int D, typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_ApplyElementMatrixEB (const FiniteElement & fel, 
                          const ElementTransformation & trafo, 
                          const FlatVector<double> elx, 
                          FlatVector<double> ely,
                          void * precomputed,
                          LocalHeap & lh) const
  {
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;

    ely = 0;
    
    auto eltype = trafo.GetElementType();
    int nfacet = ElementTopology::GetNFacets(eltype);

    Facet2ElementTrafo transform(eltype); 

    for (int k = 0; k < nfacet; k++)
      {
        HeapReset hr(lh);
        ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
          
        IntegrationRule ir_facet(etfacet, 2*fel.Order());
        IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
        const BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
        
        // ProxyUserData ud;
        // const_cast<ElementTransformation&>(trafo).userdata = &ud;
        FlatVector<> measure(mir.Size(), lh);
        
        for (int i = 0; i < mir.Size(); i++)
          {
            double len;
            if (!trafo.Boundary())
              {
                FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);
                Vec<D> normal_ref = normals[k];
                auto & mip = static_cast<const MappedIntegrationPoint<D,D>&> (mir[i]);
                Mat<D> inv_jac = mip.GetJacobianInverse();
                double det = mip.GetMeasure();
                Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
                len = L2Norm (normal);    // that's the surface measure 
                normal /= len;                   // normal vector on physical element
                
                const_cast<MappedIntegrationPoint<D,D>&> (mip).SetNV(normal);
              }
            else
              {
                if (D != 3)
                  throw Exception ("element boundary for surface elements is only possible in 3D");
                FlatVector< Vec<D-1> > normals = ElementTopology::GetNormals<D-1>(eltype);
                Vec<D-1> normal_ref = normals[k];
                
                auto & mip = static_cast<const MappedIntegrationPoint<2,3>&> (mir[i]);
                Mat<2,3> inv_jac = mip.GetJacobianInverse();
                double det = mip.GetMeasure();
                Vec<3> normal = det * Trans (inv_jac) * normal_ref;       
                len = L2Norm (normal);    // that's the surface measure
                normal /= len;                   // normal vector on physical element
                Vec<3> tang = Cross(normal, mip.GetNV());
                const_cast<MappedIntegrationPoint<2,3>&> (mip).SetTV(tang);
              }
            measure(i) = len;
          }
        
    
    
        FlatVector<> ely1(ely.Size(), lh);
        FlatMatrix<> val(mir.Size(), 1,lh);
        for (auto proxy : test_proxies)
          {
            HeapReset hr(lh);
            FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
            for (int k = 0; k < proxy->Dimension(); k++)
              {
                ud.testfunction = proxy;
                ud.test_comp = k;
                cf -> Evaluate (mir, val);
                proxyvalues.Col(k) = val.Col(0);
              }
            
            for (int i = 0; i < mir.Size(); i++)
              proxyvalues.Row(i) *= ir_facet[i].Weight() * measure(i);
            
            proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, ely1, lh);
            ely += ely1;
          }
      }
  }
  
  SymbolicFacetBilinearFormIntegrator ::
  SymbolicFacetBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb, bool eb)
    : cf(acf), vb(avb), element_boundary(eb)
  {
    if (cf->Dimension() != 1)
        throw Exception ("SymblicBFI needs scalar-valued CoefficientFunction");
    trial_cum.Append(0);
    test_cum.Append(0);    
    cf->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
          if (proxy) 
            {
              if (proxy->IsTestFunction())
                {
                  if (!test_proxies.Contains(proxy))
                    {
                      test_proxies.Append (proxy);
                      test_cum.Append(test_cum.Last()+proxy->Dimension());
                    }
                }
              else
                {                                         
                  if (!trial_proxies.Contains(proxy))
                    {
                      trial_proxies.Append (proxy);
                      trial_cum.Append(trial_cum.Last()+proxy->Dimension());
                    }
                }
            }
        });

    neighbor_testfunction = false;
    for (auto proxy : test_proxies)
      if (proxy->IsOther())
        neighbor_testfunction = true;
    
    cout << IM(3) << "num test_proxies " << test_proxies.Size() << endl;
    cout << IM(3) << "num trial_proxies " << trial_proxies.Size() << endl;
    cout << IM(3) << "cumulated test_proxy dims  " << test_cum << endl;
    cout << IM(3) << "cumulated trial_proxy dims " << trial_cum << endl;
  }


  void SymbolicFacetBilinearFormIntegrator ::
  CalcFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                   const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                   const FiniteElement & fel2, int LocalFacetNr2,
                   const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                   FlatMatrix<double> & elmat,
                   LocalHeap & lh) const
  {
    elmat = 0.0;

    if (LocalFacetNr2==-1) throw Exception ("SymbolicFacetBFI: LocalFacetNr2==-1");

    int maxorder = max2 (fel1.Order(), fel2.Order());

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    const BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    
    Facet2ElementTrafo transform2(eltype2, ElVertices2); 
    IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, ir_facet, lh);
    const BaseMappedIntegrationRule & mir2 = trafo2(ir_facet_vol2, lh);

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());
          FlatVector<> measure(mir1.Size(), lh);
          switch (trafo1.SpaceDim())
            {
	    case 1:
              {
                Vec<1> normal_ref = ElementTopology::GetNormals<1>(eltype1)[LocalFacetNr1];
                for (int i = 0; i < mir1.Size(); i++)
                  {
                    auto & mip = static_cast<const MappedIntegrationPoint<1,1>&> (mir1[i]);
                    Mat<1> inv_jac = mip.GetJacobianInverse();
                    double det = mip.GetMeasure();
                    Vec<1> normal = det * Trans (inv_jac) * normal_ref;       
                    double len = L2Norm (normal);    // that's the surface measure 
                    normal /= len;                   // normal vector on physical element
                    const_cast<MappedIntegrationPoint<1,1>&> (mip).SetNV(normal);
                    measure(i) = len;
                  }
                break;
              }
            case 2:
              {
                Vec<2> normal_ref = ElementTopology::GetNormals<2>(eltype1)[LocalFacetNr1];
                for (int i = 0; i < mir1.Size(); i++)
                  {
                    auto & mip = static_cast<const MappedIntegrationPoint<2,2>&> (mir1[i]);
                    Mat<2> inv_jac = mip.GetJacobianInverse();
                    double det = mip.GetMeasure();
                    Vec<2> normal = det * Trans (inv_jac) * normal_ref;       
                    double len = L2Norm (normal);    // that's the surface measure 
                    normal /= len;                   // normal vector on physical element
                    const_cast<MappedIntegrationPoint<2,2>&> (mip).SetNV(normal);
                    measure(i) = len;
                  }
                break;
              }
            default:
              cout << "Symbolic DG in " << trafo1.SpaceDim() << " not available" << endl;
            }
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> Evaluate (mir1, val);
                proxyvalues(STAR,l,k) = val.Col(0);
              }

          for (int i = 0; i < mir1.Size(); i++)
            proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();

          IntRange trial_range = proxy1->IsOther() ? IntRange(fel1.GetNDof(), elmat.Width()) : IntRange(0, fel1.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(fel1.GetNDof(), elmat.Height()) : IntRange(0, fel1.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);

          enum { BS = 16 };
          for (int i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(int(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), loc_elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), loc_elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  if (proxy1->IsOther())
                    proxy1->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat1, lh);
                  else
                    proxy1->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat1, lh);
                  
                  if (proxy2->IsOther())
                    proxy2->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat2, lh);
                  else
                    proxy2->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat2, lh);
                  
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(proxy1->IsOther() ? fel2 : fel1);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(proxy2->IsOther() ? fel2 : fel1);
              loc_elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }
  }



  void SymbolicFacetBilinearFormIntegrator ::
  CalcFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                   const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                   const ElementTransformation & strafo,  
                   FlatMatrix<double> & elmat,
                   LocalHeap & lh) const
  {
    // cout << "calc boundary facet matrix (DG)" << endl;
    elmat = 0.0;

    int maxorder = fel1.Order();

    auto eltype1 = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    const BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          if (proxy1->IsOther() || proxy2->IsOther()) continue;

          FlatTensor<3> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());
          FlatVector<> measure(mir1.Size(), lh);
          switch (trafo1.SpaceDim())
            {
	    case 1:
              {
                Vec<1> normal_ref = ElementTopology::GetNormals<1>(eltype1)[LocalFacetNr1];
                for (int i = 0; i < mir1.Size(); i++)
                  {
                    auto & mip = static_cast<const MappedIntegrationPoint<1,1>&> (mir1[i]);
                    Mat<1> inv_jac = mip.GetJacobianInverse();
                    double det = mip.GetMeasure();
                    Vec<1> normal = det * Trans (inv_jac) * normal_ref;       
                    double len = L2Norm (normal);    // that's the surface measure 
                    normal /= len;                   // normal vector on physical element
                    const_cast<MappedIntegrationPoint<1,1>&> (mip).SetNV(normal);
                    measure(i) = len;
                  }
                break;
              }
            case 2:
              {
                Vec<2> normal_ref = ElementTopology::GetNormals<2>(eltype1)[LocalFacetNr1];
                for (int i = 0; i < mir1.Size(); i++)
                  {
                    auto & mip = static_cast<const MappedIntegrationPoint<2,2>&> (mir1[i]);
                    Mat<2> inv_jac = mip.GetJacobianInverse();
                    double det = mip.GetMeasure();
                    Vec<2> normal = det * Trans (inv_jac) * normal_ref;       
                    double len = L2Norm (normal);    // that's the surface measure 
                    normal /= len;                   // normal vector on physical element
                    const_cast<MappedIntegrationPoint<2,2>&> (mip).SetNV(normal);
                    measure(i) = len;
                  }
                break;
              }
            default:
              cout << "Symbolic DG in " << trafo1.SpaceDim() << " not available" << endl;
            }
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> Evaluate (mir1, val);
                proxyvalues(STAR,l,k) = val.Col(0);
              }

          for (int i = 0; i < mir1.Size(); i++)
            proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();

          // IntRange trial_range = proxy1->IsOther() ? IntRange(fel1.GetNDof(), elmat.Width()) : IntRange(0, fel1.GetNDof());
          // IntRange test_range  = proxy2->IsOther() ? IntRange(fel1.GetNDof(), elmat.Height()) : IntRange(0, fel1.GetNDof());

          // auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

          enum { BS = 16 };
          for (int i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(int(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  proxy1->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat2, lh);
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(fel1);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(fel1);
              elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }

  }

  // #define STDAPPYFACET
#ifdef STDAPPYFACET
  void SymbolicFacetBilinearFormIntegrator ::
  ApplyFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                    const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                    const FiniteElement & fel2, int LocalFacetNr2,
                    const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                    FlatVector<double> elx, FlatVector<double> ely,
                    LocalHeap & lh) const
  {
    static Timer tall("SymbolicFacetBFI::Apply - all", 2); RegionTimer rall(tall);
    /*
    static Timer t("SymbolicFacetBFI::Apply", 2);
    static Timer ts1("SymbolicFacetBFI::Apply start 1", 2);
    static Timer ts2("SymbolicFacetBFI::Apply start 2", 2);
    static Timer t1("SymbolicFacetBFI::Apply 1", 2);
    static Timer t2("SymbolicFacetBFI::Apply 2", 2);
    static Timer t3("SymbolicFacetBFI::Apply 3", 2);
    */
    
    HeapReset hr(lh);
    // ts1.Start();
    /*
    Matrix<> elmat(elx.Size());
    CalcFacetMatrix(fel1, LocalFacetNr1, trafo1, ElVertices1,
                    fel2, LocalFacetNr2, trafo2, ElVertices2, elmat, lh);
    ely = elmat * elx;
    return;
    */
    
    ely = 0;
    
    FlatVector<> ely1(ely.Size(), lh);

    int maxorder = max2 (fel1.Order(), fel2.Order());

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    const BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    
    Facet2ElementTrafo transform2(eltype2, ElVertices2); 
    IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, ir_facet, lh);
    const BaseMappedIntegrationRule & mir2 = trafo2(ir_facet_vol2, lh);


    // ts1.Stop();
    // ts2.Start();

    
    // evaluate proxy-values
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    ud.fel = &fel1;   // necessary to check remember-map
    // ud.elx = &elx;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir_facet.Size(), proxy->Dimension(), lh);

    for (ProxyFunction * proxy : trial_proxies)
      {
        IntRange trial_range  = proxy->IsOther() ? IntRange(fel1.GetNDof(), elx.Size()) : IntRange(0, fel1.GetNDof());
        if (proxy->IsOther()) 
          proxy->Evaluator()->Apply(fel2, mir2, elx.Range(trial_range), ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(fel1, mir1, elx.Range(trial_range), ud.GetMemory(proxy), lh);
      }

    // ts2.Stop();
    // RegionTimer reg(t);
    // t.Start();

    FlatMatrix<> val(ir_facet.Size(), 1,lh);

    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        // t1.Start();
        FlatMatrix<> proxyvalues(ir_facet.Size(), proxy->Dimension(), lh);
        
        FlatVector<> measure(mir1.Size(), lh);
        switch (trafo1.SpaceDim())
          {
	  case 1:
            {
              Vec<1> normal_ref = ElementTopology::GetNormals<1>(eltype1)[LocalFacetNr1];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<1,1>&> (mir1[i]);
                  Mat<1> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<1> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<1,1>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
	      break;
            }
          case 2:
            {
              Vec<2> normal_ref = ElementTopology::GetNormals<2>(eltype1)[LocalFacetNr1];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<2,2>&> (mir1[i]);
                  Mat<2> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<2> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<2,2>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
              break;
            }
          default:
            cout << "Symbolic DG in " << trafo1.SpaceDim() << " not available" << endl;
          }
        
        // t1.Stop();
        // t2.Start();
        
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir1, val);
            proxyvalues.Col(k) = val.Col(0);
          }

        // t2.Stop();
        // t3.Start();

        for (int i = 0; i < mir1.Size(); i++)
          proxyvalues.Row(i) *= measure(i) * ir_facet[i].Weight();

        ely1 = 0.0;
        IntRange test_range  = proxy->IsOther() ? IntRange(fel1.GetNDof(), elx.Size()) : IntRange(0, fel1.GetNDof());
        if (proxy->IsOther()) 
          proxy->Evaluator()->ApplyTrans(fel2, mir2, proxyvalues, ely1.Range(test_range), lh);
        else
          proxy->Evaluator()->ApplyTrans(fel1, mir1, proxyvalues, ely1.Range(test_range), lh);
        
        ely += ely1;
        // t3.Stop();
      }
    // t.Stop();
  }

#else

  void SymbolicFacetBilinearFormIntegrator ::
  ApplyFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                    const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                    const FiniteElement & fel2, int LocalFacetNr2,
                    const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                    FlatVector<double> elx, FlatVector<double> ely,
                    LocalHeap & lh) const
  {
    static Timer tall("SymbolicFacetBFI::Apply - all", 2); RegionTimer rall(tall);

    static Timer t("SymbolicFacetBFI::Apply", 2);
    static Timer ts1("SymbolicFacetBFI::Apply start 1", 2);
    static Timer ts2("SymbolicFacetBFI::Apply start 2", 2);

    static Timer t1("SymbolicFacetBFI::Apply 1", 2);
    static Timer t2("SymbolicFacetBFI::Apply 2", 2);
    static Timer t3("SymbolicFacetBFI::Apply 3", 2);
    
    HeapReset hr(lh);
    // ts1.Start();
    /*
    Matrix<> elmat(elx.Size());
    CalcFacetMatrix(fel1, LocalFacetNr1, trafo1, ElVertices1,
                    fel2, LocalFacetNr2, trafo2, ElVertices2, elmat, lh);
    ely = elmat * elx;
    return;
    */
    
    ely = 0;
    
    FlatVector<> ely1(ely.Size(), lh);

    int maxorder = max2 (fel1.Order(), fel2.Order());

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    Facet2ElementTrafo transform2(eltype2, ElVertices2); 

    SIMD_IntegrationRule simd_ir_facet(etfacet, 2*maxorder);
    
    auto & simd_ir_facet_vol1 = transform1(LocalFacetNr1, simd_ir_facet, lh);
    auto & simd_mir1 = trafo1(simd_ir_facet_vol1, lh);
    
    auto & simd_ir_facet_vol2 = transform2(LocalFacetNr2, simd_ir_facet, lh);
    auto & simd_mir2 = trafo2(simd_ir_facet_vol2, lh);

    simd_mir1.ComputeNormalVectors(eltype1, LocalFacetNr1);
    
    
    // ts1.Stop();
    // ts2.Start();
    
    // evaluate proxy-values
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    ud.fel = &fel1;   // necessary to check remember-map
    // ud.elx = &elx;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir_facet.Size(), proxy->Dimension(), lh);

    for (ProxyFunction * proxy : trial_proxies)
      {
        IntRange trial_range  = proxy->IsOther() ? IntRange(fel1.GetNDof(), elx.Size()) : IntRange(0, fel1.GetNDof());
        if (proxy->IsOther())
          {
            // proxy->Evaluator()->Apply(fel2, mir2, elx.Range(trial_range), ud.GetMemory(proxy), lh);
            static_cast<const ScalarFiniteElement<2>&> (fel2).Evaluate(simd_ir_facet_vol2, elx.Range(trial_range),
                                                                       ud.GetAMemory(proxy).Row(0));
          }
        else
          {
            // proxy->Evaluator()->Apply(fel1, mir1, elx.Range(trial_range), ud.GetMemory(proxy), lh);
            static_cast<const ScalarFiniteElement<2>&> (fel1).Evaluate(simd_ir_facet_vol1, elx.Range(trial_range),
                                                                       ud.GetAMemory(proxy).Row(0));
          }
      }

    // ts2.Stop();
    // RegionTimer reg(t);
    // t.Start();

    // FlatMatrix<> val(ir_facet.Size(), 1,lh);
    // AFlatMatrix<double> aval(1, simd_ir_facet.GetNIP(),lh);

    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        // t1.Start();
        FlatMatrix<> proxyvalues(proxy->Dimension(), ir_facet.Size(), lh);
        AFlatMatrix<double> simd_proxyvalues(proxy->Dimension(), ir_facet.Size(), lh);

        // t1.Stop();
        // t2.Start();
        
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (simd_mir1, simd_proxyvalues.Rows(k,k+1));
          }

        // t2.Stop();
        // t3.Start();

        for (int i = 0; i < simd_proxyvalues.Height(); i++)
          {
            auto row = simd_proxyvalues.Row(i);
            for (int j = 0; j < row.VSize(); j++)
              row.Get(j) *= simd_mir1[i].GetMeasure().Data() * simd_ir_facet[i].Weight().Data();
          }
        
        IntRange test_range  = proxy->IsOther() ? IntRange(fel1.GetNDof(), elx.Size()) : IntRange(0, fel1.GetNDof());
        if (proxy->IsOther())
          {
            static_cast<const ScalarFiniteElement<2>&> (fel2).EvaluateTrans(simd_ir_facet_vol2,
                                                                            simd_proxyvalues.Row(0), ely1.Range(test_range));
            // proxy->Evaluator()->ApplyTrans(fel2, mir2, proxyvalues, ely1.Range(test_range), lh);
          }
        else
          {
            static_cast<const ScalarFiniteElement<2>&> (fel1).EvaluateTrans(simd_ir_facet_vol1,
                                                                            simd_proxyvalues.Row(0), ely1.Range(test_range));
            // proxy->Evaluator()->ApplyTrans(fel1, mir1, proxyvalues, ely1.Range(test_range), lh);
          }
        ely += ely1;

        // t3.Stop();
      }
    // t.Stop();
  }

  
#endif



  
  
  void SymbolicFacetBilinearFormIntegrator ::
  ApplyFacetMatrix (const FiniteElement & fel1, int LocalFacetNr,
                    const ElementTransformation & trafo1, FlatArray<int> & ElVertices,
                    const ElementTransformation & strafo,  
                    FlatVector<double> elx, FlatVector<double> ely,
                    LocalHeap & lh) const
  {
    static Timer t("SymbolicFacetBFI::ApplyFacetMatrix - boundary", 2);
    
    HeapReset hr(lh);

    /*
    Matrix<> elmat(elx.Size());
    CalcFacetMatrix(fel1, LocalFacetNr, trafo1, ElVertices, strafo, elmat, lh);
    ely = elmat * elx;
    return;
    */

    ely = 0;
    
    FlatVector<> ely1(ely.Size(), lh);

    int maxorder = fel1.Order();

    auto eltype1 = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr, ir_facet, lh);
    const BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    
    // evaluate proxy-values
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    ud.fel = &fel1;   // necessary to check remember-map
    // ud.elx = &elx;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      {
        // ud.remember[proxy] = Matrix<> (ir_facet.Size(), proxy->Dimension());
        ud.AssignMemory (proxy, ir_facet.Size(), proxy->Dimension(), lh);
        IntRange trial_range  = proxy->IsOther() ? IntRange(fel1.GetNDof(), elx.Size()) : IntRange(0, fel1.GetNDof());
        if (proxy->IsOther()) 
          // proxy->Evaluator()->Apply(fel2, mir2, elx.Range(trial_range), ud.remember[proxy], lh);
          // ud.remember[proxy] = 0.0;
          ud.GetMemory(proxy) = 0.0;
        else
          proxy->Evaluator()->Apply(fel1, mir1, elx.Range(trial_range), ud.GetMemory(proxy), lh);
      }

    RegionTimer reg(t);
    

    FlatMatrix<> val(ir_facet.Size(), 1,lh);
    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(ir_facet.Size(), proxy->Dimension(), lh);
        
        FlatVector<> measure(mir1.Size(), lh);
        switch (trafo1.SpaceDim())
          {
          case 1:
            {
              Vec<1> normal_ref = ElementTopology::GetNormals<1>(eltype1)[LocalFacetNr];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<1,1>&> (mir1[i]);
                  Mat<1> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<1> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<1,1>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
                break;
            }
          case 2:
            {
              Vec<2> normal_ref = ElementTopology::GetNormals<2>(eltype1)[LocalFacetNr];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<2,2>&> (mir1[i]);
                  Mat<2> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<2> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<2,2>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
                break;
            }
          default:
            cout << "Symbolic DG in " << trafo1.SpaceDim() << " not available" << endl;
          }
        
        
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir1, val);
            proxyvalues.Col(k) = val.Col(0);
          }

        for (int i = 0; i < mir1.Size(); i++)
          proxyvalues.Row(i) *= measure(i) * ir_facet[i].Weight();

        ely1 = 0.0;
        IntRange test_range  = proxy->IsOther() ? IntRange(fel1.GetNDof(), elx.Size()) : IntRange(0, fel1.GetNDof());
        if (proxy->IsOther())
          ;
          // nothing to do ????
          // throw Exception ("ApplyFacetBoundary: other testfunction not allowed here");
        else
          proxy->Evaluator()->ApplyTrans(fel1, mir1, proxyvalues, ely1.Range(test_range), lh);
        ely += ely1;
      }
  }
  
  

  
  SymbolicEnergy :: SymbolicEnergy (shared_ptr<CoefficientFunction> acf,
                                    VorB avb)
    : cf(acf), vb(avb)
  {
    if (cf->Dimension() != 1)
      throw Exception ("SymblicEnergy needs scalar-valued CoefficientFunction");
    
    
    cf->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
          if (proxy) 
            {
              if (!proxy->IsTestFunction())
                {                                         
                  if (!trial_proxies.Contains(proxy))
                    trial_proxies.Append (proxy);
                }
            }
        });
  }
  

  void 
  SymbolicEnergy :: CalcLinearizedElementMatrix (const FiniteElement & fel,
                                                 const ElementTransformation & trafo, 
                                                 FlatVector<double> elveclin,
                                                 FlatMatrix<double> elmat,
                                                 LocalHeap & lh) const
  {
    static Timer t("symbolicenergy - calclinearized", 2);
    static Timer td("symbolicenergy - calclinearized dmats", 2);
    RegionTimer reg(t);
    
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elveclin;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      {
        // ud.remember[proxy] = Matrix<> (ir.Size(), proxy->Dimension());
        // proxy->Evaluator()->Apply(fel, mir, elveclin, ud.remember[proxy], lh);
        ud.AssignMemory (proxy, ir.Size(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetMemory(proxy), lh);        
      }
    
    FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh), dderiv(mir.Size(), 1,lh);
    
    elmat = 0;
    


    FlatArray<FlatMatrix<>> diags(trial_proxies.Size(), lh);
    for (int k1 : Range(trial_proxies))
      {
        auto proxy = trial_proxies[k1];
        diags[k1].AssignMemory(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> EvaluateDDeriv (mir, val, deriv, dderiv);

            diags[k1].Col(k) = dderiv.Col(0);
          }
      }
           
    
    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(trial_proxies))
        {
          HeapReset hr(lh);
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = trial_proxies[l1];
          td.Start();
          // Tensor<3> proxyvalues(mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> EvaluateDDeriv (mir, val, deriv, dderiv);
                proxyvalues(STAR,l,k) = dderiv.Col(0);
                
                if (proxy1 != proxy2 || k != l)  // computed mixed second derivatives
                  {
                    proxyvalues(STAR,l,k) -= diags[k1].Col(k);
                    proxyvalues(STAR,l,k) -= diags[l1].Col(l);
                    proxyvalues(STAR,l,k) *= 0.5;
                  }
              }
          td.Stop();

          /*
          for (int i = 0; i < mir.Size(); i++)
            {
              HeapReset hr(lh);
              proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
              
              FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> dbmat1(proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
              
              proxy1->Evaluator()->CalcMatrix(fel, mir[i], bmat1, lh);
              proxy2->Evaluator()->CalcMatrix(fel, mir[i], bmat2, lh);
              dbmat1 = proxyvalues(i,STAR,STAR) * bmat1;
              elmat += Trans (bmat2) * dbmat1 | Lapack;
            }
          */
          
          for (int i = 0; i < mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();

          t.AddFlops (double (mir.Size()) * proxy1->Dimension()*elmat.Width()*elmat.Height());

          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
          int i = 0;

          enum { BS = 16 };
          for ( ; i+BS <= mir.Size(); i+=BS)
            {
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(BS*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(BS*proxy2->Dimension(), elmat.Height(), lh);

              for (int j = 0; j < BS; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }
              elmat += Trans (bbmat2) * bdbmat1 | Lapack;
            }


          if (i < mir.Size())
            {
              HeapReset hr(lh);
              int rest = mir.Size()-i;
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);
              
              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }
              elmat += Trans (bbmat2) * bdbmat1 | Lapack;
            }
          /*
          for ( ; i < mir.Size(); i++)
            {
              HeapReset hr(lh);
              // proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
              
              FlatMatrix<double,ColMajor> dbmat1(proxy2->Dimension(), elmat.Width(), lh);
              
              proxy1->Evaluator()->CalcMatrix(fel, mir[i], bmat1, lh);
              proxy2->Evaluator()->CalcMatrix(fel, mir[i], bmat2, lh);
              dbmat1 = proxyvalues(i,STAR,STAR) * bmat1;
              elmat += Trans (bmat2) * dbmat1 | Lapack;
            }
          */
        }
  }
  
  
  
  double SymbolicEnergy :: Energy (const FiniteElement & fel, 
                                   const ElementTransformation & trafo, 
                                   FlatVector<double> elx, 
                                   LocalHeap & lh) const
  {
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;

    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    FlatMatrix<> values(mir.Size(), 1, lh);
    cf -> Evaluate(mir, values);

    double sum = 0;
    for (int i = 0; i < mir.Size(); i++)
      sum += mir[i].GetWeight() * values(i,0);
    return sum;
  }

  void
  SymbolicEnergy :: ApplyElementMatrix (const FiniteElement & fel, 
                                        const ElementTransformation & trafo, 
                                        const FlatVector<double> elx, 
                                        FlatVector<double> ely,
                                        void * precomputed,
                                        LocalHeap & lh) const
  {
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;

    HeapReset hr(lh);
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
      
    ely = 0;
    FlatVector<> ely1(ely.Size(), lh);

    FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh);
      
    for (auto proxy : trial_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            cf -> EvaluateDeriv (mir, val, deriv);
            proxyvalues.Col(k) = deriv.Col(0);
          }
        
        for (int i = 0; i < mir.Size(); i++)
          proxyvalues.Row(i) *= mir[i].GetWeight();
        
        proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, ely1, lh);
        ely += ely1;
      }
  }
  
  
}

