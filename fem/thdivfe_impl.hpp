#include <thdivfe.hpp>

namespace ngfem
{

  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {    
    // Vec<DIM,AutoDiff<DIM>> adp = ip;
    TIP<DIM,AutoDiff<DIM>> pt = ip;
    static_cast<const FEL*> (this) -> 
      T_CalcShape (pt, SBLambda( [shape] (size_t nr, THDiv2Shape<DIM> val)  LAMBDA_INLINE
                                 {
                                   // shape.Row(nr) = Vec<DIM> (val);
                                   FlatVec<DIM> (&shape(nr,0)) = Vec<DIM> (val);
                                 }));
    /*
    double * pshape = &shape(0,0);
    size_t dist = shape.Dist();
    static_cast<const FEL*> (this) -> 
      T_CalcShape (&adp(0), SBLambda( [pshape,dist] (size_t nr, THDiv2Shape<DIM> val)  LAMBDA_INLINE
                                      {
                                        // shape.Row(nr) = Vec<DIM> (val);
					// FlatVec<DIM> (&shape(nr,0)) = Vec<DIM> (val);

                                        Vec<DIM> v = val;
                                        for (size_t j = 0; j < DIM; j++)
                                          pshape[nr*dist+j] = v(j);
                                      }));
    */
  }
  
  
  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  CalcDivShape (const IntegrationPoint & ip, 
		SliceVector<> divshape) const
  {  
    // Vec<DIM,AutoDiff<DIM>> adp = ip;
    TIP<DIM,AutoDiff<DIM>> pt = ip;    
    static_cast<const FEL*> (this) -> 
      T_CalcShape (pt, SBLambda( [divshape] (size_t nr, THDiv2DivShape<DIM> val) LAMBDA_INLINE
                                 {
                                   divshape(nr) = val;
                                 }));
  }

#ifndef FASTCOMPILE
  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
		   SliceMatrix<> shape) const
  {   
    // Vec<DIM,AutoDiff<DIM>> adp = mip;
    // TIP<DIM,AutoDiff<DIM>> pt(Vec<DIM, AutoDiff<DIM>> (mip));
      
    static_cast<const FEL*> (this) -> 
      T_CalcShape (TIP<DIM,AutoDiff<DIM>>(mip),
                   SBLambda([shape] (size_t nr, THDiv2Shape<DIM> val) LAMBDA_INLINE
                            {
                              FlatVec<DIM> (&shape(nr,0)) = Vec<DIM> (val);
                            }));
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET>::
  CalcMappedShape (const MappedIntegrationRule<DIM,DIM> & mir, 
                   SliceMatrix<> shape) const
  {
    for (size_t i = 0; i < mir.Size(); i++)
      CalcMappedShape (mir[i], shape.Cols(i*DIM,(i+1)*DIM));
  }
      
  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET>::
  CalcMappedShape (const SIMD<MappedIntegrationPoint<DIM,DIM>> & mip,
                   BareSliceMatrix<SIMD<double>> shapes) const
  {
    //     Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mip;
    TIP<DIM,AutoDiffRec<DIM,SIMD<double>>> adp = GetTIP(mip);    
    static_cast<const FEL*> (this) ->                 
      T_CalcShape (adp, SBLambda ([shapes] (size_t j, THDiv2Shape<DIM,SIMD<double>> shape)
                                      {
                                        Vec<DIM,SIMD<double>> vshape = shape;                                          
                                        for (size_t k = 0; k < DIM; k++)
                                          shapes(j*DIM+k, 0) = vshape(k);
                                      }));
  }
  
  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET>::
  CalcMappedShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                   BareSliceMatrix<SIMD<double>> shapes) const
  {
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        // Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
        TIP<DIM,AutoDiffRec<DIM,SIMD<double>>> adp = GetTIP(mir[i]);    
        static_cast<const FEL*> (this) ->                 
          T_CalcShape (adp, SBLambda ([shapes,i] (size_t j, THDiv2Shape<DIM,SIMD<double>> shape)
                                      {
                                        Vec<DIM,SIMD<double>> vshape = shape;                                          
                                        for (size_t k = 0; k < DIM; k++)
                                          shapes(j*DIM+k, i) = vshape(k);
                                      }));
      }
  }
  
  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET>::
  CalcMappedDivShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                      BareSliceMatrix<SIMD<double>> divshapes) const
  {
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        // Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
        TIP<DIM,AutoDiffRec<DIM,SIMD<double>>> adp = GetTIP(mir[i]);            
        static_cast<const FEL*> (this) ->                 
          T_CalcShape (adp, SBLambda ([divshapes,i] (size_t j, THDiv2DivShape<DIM,SIMD<double>> val)
                                      {
                                        divshapes(j, i) = val;
                                      }));
      }
  }

  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, 
	    FlatMatrixFixWidth<DIM> vals) const  
  {    
    // static Timer t("HDivFE - evaluate IR");
    // t.AddFlops (ir.GetNIP()* this->GetNDof());
    // RegionTimer reg(t);

    for (size_t i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM, AutoDiff<DIM>> adp = ir[i]; 

        Vec<DIM> sum = 0;
        static_cast<const FEL*> (this) -> 
          T_CalcShape (TIP<DIM,AutoDiff<DIM>>(adp),
                       SBLambda([coefs,&sum] (size_t j, THDiv2Shape<DIM> vshape)
                                {
                                  sum += coefs(j) * Vec<DIM> (vshape);
                                }));
        vals.Row(i) = sum;
      }
  }


  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  EvaluateTrans (const IntegrationRule & ir, 
                 FlatMatrixFixWidth<DIM> vals,
                 FlatVector<double> coefs) const  
  {    
    coefs = 0;
    for (size_t i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM, AutoDiff<DIM>> adp = ir[i]; 

        Vec<DIM> val = vals.Row(i);
        static_cast<const FEL*> (this) -> 
          T_CalcShape (TIP<DIM,AutoDiff<DIM>>(adp),
                       SBLambda([coefs,val] (size_t j, THDiv2Shape<DIM> vshape)
                                {
                                  coefs(j) += InnerProduct (val, Vec<DIM> (vshape));
                                }));
      }
  }



  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  Evaluate (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
  {
    static Timer t("HDivFE::Evalaute(SIMD)");
    ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        // Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
        TIP<DIM,AutoDiffRec<DIM,SIMD<double>>> adp = GetTIP(mir[i]);                    
        Vec<DIM,SIMD<double>> sum(0.0);
        static_cast<const FEL*> (this) ->         
          T_CalcShape (adp, SBLambda ([coefs,&sum] (size_t j, THDiv2Shape<DIM,SIMD<double>> shape)
                                      {
                                            Vec<DIM,SIMD<double>> vshape = shape;
                                            sum += coefs(j) * vshape;
                                      }));
        for (size_t k = 0; k < DIM; k++)
          values(k,i) = sum(k);
      }
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  AddTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
            BareSliceVector<> coefs) const
  {
    static Timer t("HDivFE::AddTrans(SIMD)");
    ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        // Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
        TIP<DIM,AutoDiffRec<DIM,SIMD<double>>> adp = GetTIP(mir[i]);                    
        Vec<DIM, SIMD<double>> vali;
        for (int k = 0; k < DIM; k++)
          vali(k) = values(k,i);
        static_cast<const FEL*> (this) -> 
          T_CalcShape (adp, SBLambda ([vali,coefs] (size_t j, THDiv2Shape<DIM,SIMD<double>> shape)
                                      {
                                        Vec<DIM,SIMD<double>> vshape = shape;                                            
                                        SIMD<double> sum = 0.0;
                                        for (size_t k = 0; k < DIM; k++)
                                          sum += vali(k) * vshape(k);
                                        coefs(j) += HSum(sum);
                                      }));
      }
  }


  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  EvaluateDiv (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs,
               BareVector<SIMD<double>> values) const
  {
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        // Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
        TIP<DIM,AutoDiffRec<DIM,SIMD<double>>> adp = GetTIP(mir[i]);                            
        SIMD<double> sum(0.0);
        static_cast<const FEL*> (this) ->         
          T_CalcShape (adp, SBLambda ([&sum,coefs] (size_t j, THDiv2DivShape<DIM,SIMD<double>> divshape)
                                      {
                                        // SIMD<double> simdshape = divshape;
                                        SIMD<double> simdshape = divshape.Get();
                                        sum += coefs(j) * simdshape;
                                      }));
        values(i) = sum;
      }
  }
  
  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  AddDivTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareVector<SIMD<double>> values,
               BareSliceVector<> coefs) const
  {
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        // Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
        TIP<DIM,AutoDiffRec<DIM,SIMD<double>>> adp = GetTIP(mir[i]);                                    
        SIMD<double> vali = values(i);
        static_cast<const FEL*> (this) -> 
          T_CalcShape (adp, SBLambda ([coefs,vali] (size_t j, THDiv2DivShape<DIM,SIMD<double>> divshape)
                                      {
                                        SIMD<double> simdshape = divshape.Get();
                                        coefs(j) += HSum(simdshape*vali);
                                      }));
      }
  }

  
#endif
}
