/*********************************************************************/
/* File:   tensorproductintegrator.cpp                               */
/* Author: Gerhard Kitzler                                           */
/* Date:   January 2017                                              */
/*********************************************************************/
/* 
   integrators for tensor product spaces
*/

#include <fem.hpp>

namespace ngfem
{
  void TensorProductBilinearFormIntegrator :: ApplyXElementMatrix(const FiniteElement & fel, 
                          const ElementTransformation & trafo, 
                          const FlatMatrix<double> elx, 
                          void * axheap,
                          BaseMappedIntegrationRule * mirx,
                          LocalHeap & lh) const
  {
     // fel..... x space finite element
     // trafo... x space ElementTransformation
     // elx..... coeffs of a whole slice of Y elements, dim(elx) = (ndofx , space_y.GetNDof() )
     // ely..... is currently not needed
     // precomp. is the LocalHeap to store the x evaluations
     // mirx.... is construced in loop over x elements -> shall be recycled in ApplyYElementMatrix
     LocalHeap & xheap = *(static_cast<LocalHeap *>(axheap));
     ProxyUserData * xevaluations = new (xheap) ProxyUserData(trial_proxies.Size()+test_proxies.Size(), xheap);
     const_cast<ElementTransformation&>(trafo).userdata = xevaluations;
     xevaluations->lh = &xheap;
     for (ProxyFunction * proxy : trial_proxies)
     {
         int dimx = dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(0)->Dim();
         xevaluations->AssignMemory (proxy, mirx->IR().Size()*dimx, elx.Width(), xheap);
         FlatMatrix<double, ColMajor> bmatx( mirx->IR().Size()*dimx, fel.GetNDof(),lh ); 
         dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(0)->CalcMatrix(fel,*mirx,bmatx,lh);
         xevaluations->GetMemory(proxy) = bmatx*elx;
     }
     for (ProxyFunction * proxy : test_proxies)
     {
        int dimx = dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(0)->Dim();
        xevaluations->AssignMemory (proxy, mirx->IR().Size()*dimx, elx.Width(), xheap);
     }
  }

  void TensorProductBilinearFormIntegrator :: ApplyYElementMatrix(const FiniteElement & tpfel, 
                          const ElementTransformation & tptrafo, 
                          IntRange dnumsy, 
                          void * axevaluations,
                          BaseMappedIntegrationRule * mirx,
                          LocalHeap & lh) const
  {
     const FiniteElement & fel = *dynamic_cast<const TPHighOrderFE &>(tpfel).elements[1];
     const ElementTransformation & trafo = dynamic_cast<const TPElementTransformation &>(tptrafo).GetTrafo(1);
     ProxyUserData & xevals = *(static_cast<ProxyUserData *>(axevaluations));
     const IntegrationRule & ir = SelectIntegrationRule(fel.ElementType(),2*fel.Order());
     BaseMappedIntegrationRule & miry = trafo(ir, lh);
     ProxyUserData ud(trial_proxies.Size(), lh);
     const_cast<ElementTransformation&>(tptrafo).userdata = &ud;
     ud.lh = &lh;
     ud.fel = &fel;
     int ndofy = fel.GetNDof();
     int nipx = mirx->Size();
     int nipy = miry.Size();
     int niptp = nipx*nipy;
     TPMappedIntegrationRule * tpmir = new (lh) TPMappedIntegrationRule(*mirx, miry, TPIntegrationRule(niptp), tptrafo);
     tpmir->SetFacet(0);
     for (ProxyFunction * proxy : trial_proxies)
     {
         ud.AssignMemory (proxy, (*mirx).IR().Size()*ir.Size(), proxy->Dimension(), lh);
         dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->ApplyY(fel,miry,ud.GetMemory(proxy), xevals.GetMemory(proxy).Cols(dnumsy),lh);
     }
     FlatMatrix<> val(niptp, 1,lh);
     for (auto proxy : test_proxies)
       {
         FlatMatrix<> proxyvalues(niptp, proxy->Dimension(), lh);
         for (int k = 0; k < proxy->Dimension(); k++)
           {
             ud.testfunction = proxy;
             ud.test_comp = k;
             cf -> Evaluate (*tpmir, val);
             proxyvalues.Col(k) = val.Col(0);
           }                    
         for(int i=0,ii=0;i<nipx;i++)
           for(int j=0;j<nipy;j++,ii++)
             proxyvalues.Row(ii) *= (*mirx)[i].GetWeight()*miry[j].GetWeight();

         dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())
         ->ApplyYTrans(fel,miry,proxyvalues,xevals.GetMemory(proxy).Cols(dnumsy),lh);
       }
  }

  void TensorProductBilinearFormIntegrator::ApplyXElementMatrixTrans(const FiniteElement & fel, 
            const ElementTransformation & trafo, 
            FlatMatrix<double> ely,
            void * yapplytrans,
            BaseMappedIntegrationRule * mirx,
            LocalHeap & lh) const
  {
     // fel..... x space finite element
     // trafo... x space ElementTransformation
     // elx..... coeffs of a whole slice of Y elements, dim(elx) = (ndofx , space_y.GetNDof() )
     // ely..... is currently not needed
     // precomp. is the LocalHeap to store the x evaluations
     // mirx.... is construced in loop over x elements -> shall be recycled in ApplyYElementMatrix
     ProxyUserData * ud = (static_cast<ProxyUserData *>(yapplytrans));//new (xheap) ProxyUserData(trial_proxies.Size(), xheap);
     ely = 0.0;
     for (auto proxy : test_proxies)
     {
         int dimx = dynamic_cast<TPDifferentialOperator * >(proxy->Evaluator().get())->GetEvaluators(0)->Dim();
         FlatMatrix<double, ColMajor> bmatx( mirx->IR().Size()*dimx, fel.GetNDof(),lh ); 
         dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(0)->CalcMatrix(fel,*mirx,bmatx,lh);
         ely += Trans(bmatx)*ud->GetMemory(proxy);
     }
  }
 
   void TensorProductFacetBilinearFormIntegrator :: 
   ApplyFacetMatrix (const FiniteElement & fel, int LocalFacetNr,
                       const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                       const ElementTransformation & seltrans, FlatArray<int> & SElVertices, 
                       FlatVector<double> elx, FlatVector<double> ely, 
                       LocalHeap & lh) const
   {
     int xfacet = 0;
     if(LocalFacetNr>=10)
     {
       xfacet = 1;
       LocalFacetNr-=10;
     }
     HeapReset hr(lh);
     ely = 0;
     FlatVector<> ely1(ely.Size(), lh);
     const TPHighOrderFE & volumefel = dynamic_cast<const TPHighOrderFE &>(fel);
     int order_fac = volumefel.elements[xfacet]->Order();
     int order_vol = volumefel.elements[1-xfacet]->Order();
     auto eltype1 = dynamic_cast<const TPElementTransformation &>(eltrans).GetTrafo(xfacet).GetElementType();
     
     auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);
     auto etvol = dynamic_cast<const TPElementTransformation &>(eltrans).GetTrafo(1-xfacet).GetElementType();
     ArrayMem<const IntegrationRule *,2> irs(2);
     const IntegrationRule & ir_facet = SelectIntegrationRule(etfacet, 2*order_fac);
     Facet2ElementTrafo transform(eltype1, ElVertices); 
     const IntegrationRule & ir_facet_vol = transform(LocalFacetNr, ir_facet, lh);
     const IntegrationRule &ir_vol = SelectIntegrationRule(etvol,2*order_vol);
     irs[xfacet] = &ir_facet_vol;
     irs[1-xfacet] = &ir_vol;
     TPIntegrationRule tpir(irs);
     BaseMappedIntegrationRule & tpmir = eltrans(tpir, lh);
     auto & smir = seltrans(ir_facet, lh);
     TPMappedIntegrationRule tpsmir(tpir,eltrans);
     tpsmir.GetIRs()[xfacet] = &smir;
     tpsmir.GetIRs()[1-xfacet] = dynamic_cast<TPMappedIntegrationRule &>(tpmir).GetIRs()[1-xfacet];
     // evaluate proxy-values
     ProxyUserData ud(trial_proxies.Size(), lh);
     dynamic_cast<TPMappedIntegrationRule &>(tpmir).SetFacet(xfacet);
     dynamic_cast<TPMappedIntegrationRule &>(tpsmir).SetFacet(xfacet);
     const_cast<ElementTransformation&>(eltrans).userdata = &ud;
     const_cast<ElementTransformation&>(seltrans).userdata = &ud;
     ud.fel = &volumefel;   // necessary to check remember-map
     ud.lh = &lh;
     for (ProxyFunction * proxy : trial_proxies)
       ud.AssignMemory (proxy, tpmir.Size(), proxy->Dimension(), lh);
     for (ProxyFunction * proxy : trial_proxies)
       if (! (proxy->IsOther() && proxy->BoundaryValues()))
         proxy->Evaluator()->Apply(volumefel, tpmir, elx, ud.GetMemory(proxy), lh);
     for (ProxyFunction * proxy : trial_proxies)    
       if (proxy->IsOther() && proxy->BoundaryValues())
         proxy->BoundaryValues()->Evaluate (tpsmir, ud.GetMemory(proxy));
 
     dynamic_cast<TPMappedIntegrationRule &>(tpmir).GetIRs()[xfacet]->ComputeNormalsAndMeasure (eltype1, LocalFacetNr);
     FlatMatrixFixWidth<1> val(tpmir.Size(),lh);
     for (auto proxy : test_proxies)
       {
         HeapReset hr(lh);
         FlatMatrix<> proxyvalues(tpmir.Size(), proxy->Dimension(), lh);
         for (int k = 0; k < proxy->Dimension(); k++)
           {
             ud.testfunction = proxy;
             ud.test_comp = k;
             cf -> Evaluate (tpmir, val);
             proxyvalues.Col(k) = val.Col(0);
           }
         TPMappedIntegrationRule & ttpmir = dynamic_cast<TPMappedIntegrationRule &>(tpmir);
         int ii=0;
         for(int i=0;i<ttpmir.GetIRs()[0]->Size();i++)
           for(int j=0;j<ttpmir.GetIRs()[1]->Size();j++)
             proxyvalues.Row(ii++) *= ttpmir.GetIRs()[0]->operator[](i).GetWeight()*ttpmir.GetIRs()[1]->operator[](j).GetWeight();
 
         ely1 = 0.0;
         if (proxy->IsOther() && proxy->BoundaryValues())
           ;  // nothing to do 
         else
           proxy->Evaluator()->ApplyTrans(volumefel, tpmir, proxyvalues, ely1, lh);
         ely += ely1;
       }
  }  

    void TensorProductFacetBilinearFormIntegrator :: 
    ApplyXFacetMatrix(const FiniteElement & felx1, 
            const ElementTransformation & trafox1,
            const FiniteElement & felx2,
            const ElementTransformation & trafox2,
            const FlatMatrix<double> elx, 
            void * axheap,
            BaseMappedIntegrationRule * mirx1,
            BaseMappedIntegrationRule * mirx2,
            LocalHeap & lh) const
    {
      // fel..... x space finite element
      // trafo... x space ElementTransformation
      // elx..... coeffs of a whole slice of Y elements, dim(elx) = (ndofx , space_y.GetNDof() )
      // ely..... is currently not needed
      // precomp. is the LocalHeap to store the x evaluations
      // mirx.... is construced in loop over x elements -> shall be recycled in ApplyYElementMatrix
      LocalHeap & xheap = *(static_cast<LocalHeap *>(axheap));
      ProxyUserData * xevaluations = new (xheap) ProxyUserData(trial_proxies.Size()+test_proxies.Size(), xheap);
      const_cast<ElementTransformation&>(trafox1).userdata = xevaluations;
      xevaluations->lh = &xheap;
      for (ProxyFunction * proxy : trial_proxies)
      {
        int dimx = dynamic_cast<TPDifferentialOperator * >(proxy->Evaluator().get())->GetEvaluators(0)->Dim();
        xevaluations->AssignMemory (proxy, mirx1->IR().Size()*dimx, elx.Width(), xheap);
        IntRange trial_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*felx1.GetNDof(), elx.Height()) : IntRange(0, proxy->Evaluator()->BlockDim()*felx1.GetNDof());       
        if(proxy->IsOther())
        {
          FlatMatrix<double, ColMajor> bmatx( mirx1->IR().Size()*dimx, felx2.GetNDof(),lh ); 
          dynamic_cast<TPDifferentialOperator*>(  proxy->Evaluator().get()  )->GetEvaluators(0)->CalcMatrix(felx2,*mirx2,bmatx,lh);
          xevaluations->GetMemory(proxy) = bmatx*elx.Rows(trial_range);
        }
        else
        {
          FlatMatrix<double, ColMajor> bmatx( mirx1->IR().Size()*dimx, felx1.GetNDof(),lh ); 
          dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(0)->CalcMatrix(felx1,*mirx1,bmatx,lh);
          xevaluations->GetMemory(proxy) = bmatx*elx.Rows(trial_range);
        }
        
      }
      for (ProxyFunction * proxy : test_proxies)
      {
        int dimx = dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(0)->Dim();
        xevaluations->AssignMemory (proxy, mirx1->IR().Size()*dimx, elx.Width(), xheap);
      }
    }

    void TensorProductFacetBilinearFormIntegrator :: 
    ApplyYElementMatrix(const FiniteElement & tpfel, 
                          const ElementTransformation & tptrafo, 
                          IntRange dnumsy, 
                          void * axevaluations,
                          BaseMappedIntegrationRule * mirx,
                          LocalHeap & lh) const
    {
      const FiniteElement & fel = *dynamic_cast<const TPHighOrderFE &>(tpfel).elements[1];
      const ElementTransformation & trafo = dynamic_cast<const TPElementTransformation &>(tptrafo).GetTrafo(1);
      ProxyUserData & xevals = *(static_cast<ProxyUserData *>(axevaluations));
      const IntegrationRule & ir = SelectIntegrationRule(fel.ElementType(),2*fel.Order());
      BaseMappedIntegrationRule & miry = trafo(ir, lh);
      ProxyUserData ud(trial_proxies.Size(), lh);
      const_cast<ElementTransformation&>(tptrafo).userdata = &ud;
      ud.lh = &lh;
      ud.fel = &fel;
      int ndofy = fel.GetNDof();
      for (ProxyFunction * proxy : trial_proxies)
      {
          ud.AssignMemory (proxy, (*mirx).IR().Size()*ir.Size(), proxy->Dimension(), lh);
          dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->ApplyY(fel,miry,ud.GetMemory(proxy), xevals.GetMemory(proxy).Cols(dnumsy),lh);
      }
      int nipx = mirx->Size();
      int nipy = miry.Size();
      int niptp = nipx*nipy;
      TPMappedIntegrationRule * tpmir = new (lh) TPMappedIntegrationRule(*mirx, miry, TPIntegrationRule(niptp), tptrafo);
      tpmir->SetFacet(0);
      FlatMatrix<> val(niptp, 1,lh);
      
      for (auto proxy : test_proxies)
      {
        FlatMatrix<> proxyvalues(niptp, proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
        {
          ud.testfunction = proxy;
          ud.test_comp = k;
          cf -> Evaluate (*tpmir, val);
          proxyvalues.Col(k) = val.Col(0);
        }                    
        for(int i=0,ii=0;i<nipx;i++)
          for(int j=0;j<nipy;j++,ii++)
            proxyvalues.Row(ii) *= (*mirx)[i].GetWeight()*miry[j].GetWeight();

        dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())
        ->ApplyYTrans(fel,miry,proxyvalues,xevals.GetMemory(proxy).Cols(dnumsy),lh);
      }
    }

    void TensorProductFacetBilinearFormIntegrator :: 
    ApplyXFacetMatrixTrans(const FiniteElement & felx1, 
            const ElementTransformation & trafox1,
            const FiniteElement & felx2,
            const ElementTransformation & trafox2,
            FlatMatrix<double> ely,
            void * yapplytrans,
            BaseMappedIntegrationRule * mirx1,
            BaseMappedIntegrationRule * mirx2,
            LocalHeap & lh) const
    {
      ProxyUserData * ud = (static_cast<ProxyUserData *>(yapplytrans));
      for (auto proxy : test_proxies)
      {
        int dimx = dynamic_cast<TPDifferentialOperator * >(proxy->Evaluator().get())->GetEvaluators(0)->Dim();
        IntRange test_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*felx1.GetNDof(), ely.Height()) : IntRange(0, proxy->Evaluator()->BlockDim()*felx1.GetNDof());
        if(proxy->IsOther())
        {
          FlatMatrix<double, ColMajor> bmatx( mirx1->IR().Size()*dimx, felx2.GetNDof(),lh ); 
          dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(0)->CalcMatrix(felx2,*mirx2,bmatx,lh);
          ely.Rows(test_range) += Trans(bmatx)*ud->GetMemory(proxy);
        }
        else
        {
          FlatMatrix<double, ColMajor> bmatx( mirx1->IR().Size()*dimx, felx1.GetNDof(),lh );
          dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(0)->CalcMatrix(felx1,*mirx1,bmatx,lh);
          ely.Rows(test_range) += Trans(bmatx)*ud->GetMemory(proxy);
        }
      }
    }

    void TensorProductFacetBilinearFormIntegrator ::
    ApplyXElementMatrix(const FiniteElement & tpfel, 
            const ElementTransformation & tptrafo, 
            IntRange dnumsx, 
            void * ayevaluations,
            BaseMappedIntegrationRule * miry,
            LocalHeap & lh) const
    {
      const FiniteElement & fel = *dynamic_cast<const TPHighOrderFE &>(tpfel).elements[0];
      const ElementTransformation & trafo = dynamic_cast<const TPElementTransformation &>(tptrafo).GetTrafo(0);
      ProxyUserData & yevals = *(static_cast<ProxyUserData *>(ayevaluations));
      const IntegrationRule & ir = SelectIntegrationRule(fel.ElementType(),2*fel.Order());
      BaseMappedIntegrationRule & mirx = trafo(ir, lh);
      ProxyUserData ud(trial_proxies.Size(), lh);
      const_cast<ElementTransformation&>(tptrafo).userdata = &ud;
      ud.lh = &lh;
      ud.fel = &fel;
      int ndofx = fel.GetNDof();
      for (ProxyFunction * proxy : trial_proxies)
      {
          ud.AssignMemory (proxy, (*miry).IR().Size()*ir.Size(), proxy->Dimension(), lh);
          dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->ApplyX(fel,mirx,ud.GetMemory(proxy), yevals.GetMemory(proxy).Rows(dnumsx),lh);
      }
      int nipx = mirx.Size();
      int nipy = miry->Size();
      int niptp = nipx*nipy;
      TPMappedIntegrationRule * tpmir = new (lh) TPMappedIntegrationRule(mirx, *miry, TPIntegrationRule(niptp), tptrafo);
      tpmir->SetFacet(1);
      FlatMatrix<> val(niptp, 1,lh);

      for (auto proxy : test_proxies)
      {
        FlatMatrix<> proxyvalues(niptp, proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
        {
          ud.testfunction = proxy;
          ud.test_comp = k;
          cf -> Evaluate (*tpmir, val);
          proxyvalues.Col(k) = val.Col(0);
        }                    
        for(int i=0,ii=0;i<nipx;i++)
          for(int j=0;j<nipy;j++,ii++)
            proxyvalues.Row(ii) *= (*miry)[j].GetWeight()*mirx[i].GetWeight();

        dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())
        ->ApplyXTrans(fel,mirx,proxyvalues,yevals.GetMemory(proxy).Rows(dnumsx),lh);
      }
    }
    
    void TensorProductFacetBilinearFormIntegrator ::
    ApplyYFacetMatrix(const FiniteElement & fely1, 
        const ElementTransformation & trafoy1,
        const FiniteElement & fely2,
        const ElementTransformation & trafoy2,
        const FlatMatrix<double> elx, 
        void * ayheap,
        BaseMappedIntegrationRule * miry1,
        BaseMappedIntegrationRule * miry2,
        LocalHeap & lh) const
    {
      // fel..... x space finite element
      // trafo... x space ElementTransformation
      // elx..... coeffs of a whole slice of Y elements, dim(elx) = (ndofx , space_y.GetNDof() )
      // ely..... is currently not needed
      // precomp. is the LocalHeap to store the x evaluations
      // mirx.... is construced in loop over x elements -> shall be recycled in ApplyYElementMatrix
      LocalHeap & yheap = *(static_cast<LocalHeap *>(ayheap));
      ProxyUserData * yevaluations = new (yheap) ProxyUserData(trial_proxies.Size()+test_proxies.Size(), yheap);
      const_cast<ElementTransformation&>(trafoy1).userdata = yevaluations;
      yevaluations->lh = &yheap;
      for (ProxyFunction * proxy : trial_proxies)
      {
        int dimy = dynamic_cast<TPDifferentialOperator * >(proxy->Evaluator().get())->GetEvaluators(1)->Dim();
        yevaluations->AssignMemory (proxy, elx.Width(),miry1->IR().Size()*dimy, yheap);
        IntRange trial_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*fely1.GetNDof(), elx.Height()) : IntRange(0, proxy->Evaluator()->BlockDim()*fely1.GetNDof());       
        if(proxy->IsOther())
        {
          FlatMatrix<double, ColMajor> bmaty( miry1->IR().Size()*dimy, fely2.GetNDof(),lh ); 
          dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(1)->CalcMatrix(fely2,*miry2,bmaty,lh);
          yevaluations->GetMemory(proxy) = Trans(bmaty*elx.Rows(trial_range));
        }
        else
        {
          FlatMatrix<double, ColMajor> bmaty( miry1->IR().Size()*dimy, fely1.GetNDof(),lh ); 
          dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(1)->CalcMatrix(fely1,*miry1,bmaty,lh);
          yevaluations->GetMemory(proxy) = Trans(bmaty*elx.Rows(trial_range));
        }
        
      }
      for (ProxyFunction * proxy : test_proxies)
      {
        int dimy = dynamic_cast<TPDifferentialOperator *>(proxy->Evaluator().get())->GetEvaluators(1)->Dim();
        yevaluations->AssignMemory (proxy, elx.Width(), miry1->IR().Size()*dimy, yheap);
      }
    }
    
    void TensorProductFacetBilinearFormIntegrator ::
    ApplyYFacetMatrixTrans(const FiniteElement & fely1, 
        const ElementTransformation & trafoy1,
        const FiniteElement & fely2,
        const ElementTransformation & trafoy2,
        FlatMatrix<double> ely,
        void * xapplytrans,
        BaseMappedIntegrationRule * miry1,
        BaseMappedIntegrationRule * miry2,
        LocalHeap & lh) const
    {
      ProxyUserData * ud = (static_cast<ProxyUserData *>(xapplytrans));
      for (auto proxy : test_proxies)
      {
        int dimy = dynamic_cast<TPDifferentialOperator * >(proxy->Evaluator().get())->GetEvaluators(1)->Dim();
        IntRange test_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*fely1.GetNDof(), ely.Height()) : IntRange(0, proxy->Evaluator()->BlockDim()*fely1.GetNDof());
        if(proxy->IsOther())
        {
          FlatMatrix<double, ColMajor> bmaty( miry1->IR().Size()*dimy, fely2.GetNDof(),lh ); 
          dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(1)->CalcMatrix(fely2,*miry2,bmaty,lh);
          ely.Rows(test_range) += Trans(ud->GetMemory(proxy)*bmaty);
        }
        else
        {
          FlatMatrix<double, ColMajor> bmaty( miry1->IR().Size()*dimy, fely1.GetNDof(),lh );
          dynamic_cast<TPDifferentialOperator*>(proxy->Evaluator().get())->GetEvaluators(1)->CalcMatrix(fely1,*miry1,bmaty,lh);
          ely.Rows(test_range) += Trans(ud->GetMemory(proxy)*bmaty);
        }
      }
    }
    
}
 