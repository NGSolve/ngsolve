/*********************************************************************/
/* File:   symbolicintegrator.cpp                                    */
/* Author: Joachim Schoeberl                                         */
/* Date:   August 2015                                               */
/*********************************************************************/
/* 
   Symbolic integrators
*/

#include <fem.hpp>
  
namespace ngfem
{


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
    HeapReset hr(lh);
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    
    FlatVector<SCAL> elvec1(fel.GetNDof(), lh);
    
    FlatMatrix<SCAL> values(ir.Size(), 1, lh);
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    
    elvec = 0;
    for (auto proxy : proxies)
      {
        FlatMatrix<SCAL> proxyvalues(ir.Size(), proxy->Dimension(), lh);
          for (int k = 0; k < proxy->Dimension(); k++)
            {
              ud.testfunction = proxy;
              ud.test_comp = k;
              
              cf -> Evaluate (mir, values);
              
              for (int i = 0; i < mir.Size(); i++)
                values.Row(i) *= mir[i].GetWeight();
              proxyvalues.Col(k) = values.Col(0);
            }
          
          proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
          elvec += elvec1;
      }
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

    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    
    elmat = 0;
    for (int i = 0; i < mir.Size(); i++)
      {
        auto & mip = mir[i];
        
        for (auto proxy1 : trial_proxies)
          for (auto proxy2 : test_proxies)
            {
              HeapReset hr(lh);
              FlatMatrix<SCAL> proxyvalues(proxy2->Dimension(), proxy1->Dimension(), lh);
              for (int k = 0; k < proxy1->Dimension(); k++)
                for (int l = 0; l < proxy2->Dimension(); l++)
                  {
                    ud.trialfunction = proxy1;
                    ud.trial_comp = k;
                    ud.testfunction = proxy2;
                    ud.test_comp = l;
                    
                    Vec<1,SCAL> result;
                    cf->Evaluate (mip, result);
                    proxyvalues(l,k) = mip.GetWeight() * result(0);
                  }

              FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), fel.GetNDof(), lh);
              FlatMatrix<SCAL,ColMajor> dbmat1(proxy2->Dimension(), fel.GetNDof(), lh);
              FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), fel.GetNDof(), lh);
              
              proxy1->Evaluator()->CalcMatrix(fel, mip, bmat1, lh);
              proxy2->Evaluator()->CalcMatrix(fel, mip, bmat2, lh);
              
              dbmat1 = proxyvalues * bmat1;
              elmat += Trans (bmat2) * dbmat1;
            }
      }
  }
  
  
  

  template <int D, typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcElementMatrixEB (const FiniteElement & fel,
                           const ElementTransformation & trafo, 
                           FlatMatrix<SCAL> elmat,
                           LocalHeap & lh) const
      
    {
      elmat = 0;

      auto eltype = trafo.GetElementType();
      int nfacet = ElementTopology::GetNFacets(eltype);

      Facet2ElementTrafo transform(eltype); 
      FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);

      for (int k = 0; k < nfacet; k++)
        {
          HeapReset hr(lh);
          ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
        
          Vec<2> normal_ref = normals[k];
        
          IntegrationRule ir_facet(etfacet, 2*fel.Order());
          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
          MappedIntegrationRule<D,D> mir(ir_facet_vol, trafo, lh);
          
          
          ProxyUserData ud;
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          
          for (int i = 0; i < mir.Size(); i++)
            {
              auto & mip = mir[i];
              Mat<D> inv_jac = mip.GetJacobianInverse();
              double det = mip.GetMeasure();
              Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
              double len = L2Norm (normal);    // that's the surface measure 
              normal /= len;                   // normal vector on physical element

              mip.SetNV(normal);
              
              for (auto proxy1 : trial_proxies)
                for (auto proxy2 : test_proxies)
                  {
                    HeapReset hr(lh);
                    FlatMatrix<SCAL> proxyvalues(proxy2->Dimension(), proxy1->Dimension(), lh);
                    for (int k = 0; k < proxy1->Dimension(); k++)
                      for (int l = 0; l < proxy2->Dimension(); l++)
                        {
                          ud.trialfunction = proxy1;
                          ud.trial_comp = k;
                          ud.testfunction = proxy2;
                          ud.test_comp = l;
                      
                          Vec<1,SCAL> result;
                          cf->Evaluate (mip, result);
                          proxyvalues(l,k) = ir_facet[i].Weight() * len * result(0);
                        }
                    
                    FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), fel.GetNDof(), lh);
                    FlatMatrix<SCAL,ColMajor> dbmat1(proxy2->Dimension(), fel.GetNDof(), lh);
                    FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), fel.GetNDof(), lh);
                    
                    proxy1->Evaluator()->CalcMatrix(fel, mip, bmat1, lh);
                    proxy2->Evaluator()->CalcMatrix(fel, mip, bmat2, lh);
                    
                    dbmat1 = proxyvalues * bmat1;
                    elmat += Trans (bmat2) * dbmat1;
                  }
            }
        }
    }


  



  
  SymbolicEnergy :: SymbolicEnergy (shared_ptr<CoefficientFunction> acf)
    : cf(acf)
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
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elveclin;
    ud.lh = &lh;
    FlatVector<> val(1,lh), deriv(1,lh), dderiv(1,lh);
    
    elmat = 0;
    
    for (int i = 0; i < mir.Size(); i++)
      {
        auto & mip = mir[i];
        
        for (auto proxy1 : trial_proxies)
          for (auto proxy2 : trial_proxies)
            {
              HeapReset hr(lh);
              FlatMatrix<> proxyvalues(proxy2->Dimension(), proxy1->Dimension(), lh);
              for (int k = 0; k < proxy1->Dimension(); k++)
                for (int l = 0; l < proxy2->Dimension(); l++)
                  {
                    ud.trialfunction = proxy1;
                    ud.trial_comp = k;
                    ud.testfunction = proxy2;
                    ud.test_comp = l;
                    
                    cf -> EvaluateDDeriv (mip, val, deriv, dderiv);
                    proxyvalues(l,k) = dderiv(0);
                    
                    if (proxy1 != proxy2 || k != l)
                      {
                        ud.trialfunction = proxy1;
                        ud.trial_comp = k;
                        ud.testfunction = proxy1;
                        ud.test_comp = k;
                        cf -> EvaluateDDeriv (mip, val, deriv, dderiv);
                        proxyvalues(l,k) -= dderiv(0);
                        
                        ud.trialfunction = proxy2;
                        ud.trial_comp = l;
                        ud.testfunction = proxy2;
                        ud.test_comp = l;
                        cf -> EvaluateDDeriv (mip, val, deriv, dderiv);
                        proxyvalues(l,k) -= dderiv(0);
                      }
                  }
              
              proxyvalues *= mip.GetWeight();
              
              FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), fel.GetNDof(), lh);
              FlatMatrix<double,ColMajor> dbmat1(proxy2->Dimension(), fel.GetNDof(), lh);
              FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), fel.GetNDof(), lh);
              
              proxy1->Evaluator()->CalcMatrix(fel, mip, bmat1, lh);
              proxy2->Evaluator()->CalcMatrix(fel, mip, bmat2, lh);
              
              dbmat1 = proxyvalues * bmat1;
              elmat += Trans (bmat2) * dbmat1;
            }
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
    
    double sum = 0;
    for (int i = 0; i < mir.Size(); i++)
      {
        auto & mip = mir[i];
        sum += mip.GetWeight() * cf -> Evaluate (mip);
      }
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
    FlatVector<> val(1,lh), deriv(1,lh);

    for (int i = 0; i < mir.Size(); i++)
      {
        auto & mip = mir[i];
        
        for (auto proxy : trial_proxies)
          {
            HeapReset hr(lh);
            FlatVector<> proxyvalues(proxy->Dimension(), lh);
            for (int k = 0; k < proxy->Dimension(); k++)
              {
                ud.trialfunction = proxy;
                ud.trial_comp = k;
                cf -> EvaluateDeriv (mip, val, deriv);
                proxyvalues(k) = mip.GetWeight() * deriv(0);
              }

            const MappedIntegrationPoint<2,2> & cmip =
              static_cast<const MappedIntegrationPoint<2,2>&> (mip);

            proxy->Evaluator()->ApplyTrans(fel, mip, proxyvalues, ely1, lh);
            ely += ely1;
          }
      }
  }
  
  
}

