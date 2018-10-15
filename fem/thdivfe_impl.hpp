#include <thdivfe.hpp>

namespace ngfem
{

  /*
  template <int DIM>
  INLINE auto GetTIPHDiv (const IntegrationPoint & ip)
  {
    TIP<DIM,AutoDiff<DIM>> tip = ip;
    return tip;
  }
  */
  
  template<int DIMS, int DIMR>
  INLINE auto GetTIPHDiv (const MappedIntegrationPoint<DIMS,DIMR> & mip)
  {
    return GetTIP (mip);
  }

  template<int DIMS, int DIMR>
  INLINE auto GetTIPHDiv (const SIMD<MappedIntegrationPoint<DIMS,DIMR>> & mip)
  {
    return GetTIP (mip);
  }

  template <int DIMR>
  INLINE auto GetTIPHDiv (const MappedIntegrationPoint<2,DIMR> & mip)
  {
    TIP<2,AutoDiffRec<DIMR>> adp;      
    Mat<DIMR,2> jac = mip.GetJacobian();
    jac *= 1/mip.GetJacobiDet();
    const auto &ip = mip.IP();
    adp.x.Value() = ip(0);
    adp.y.Value() = ip(1);
    for (int i = 0; i < DIMR; i++)
      {
        adp.x.DValue(i) = jac(i,1);
        adp.y.DValue(i) = -jac(i,0);
      }
    return adp;
  }
  
  template <int DIMR>
  INLINE auto GetTIPHDiv (const SIMD<MappedIntegrationPoint<2,DIMR>> & mip)
  {
    TIP<2,AutoDiffRec<DIMR,SIMD<double>>> adp;      
    Mat<DIMR,2,SIMD<double>> jac = mip.GetJacobian();
    jac *= 1/mip.GetJacobiDet();
    const auto &ip = mip.IP();
    adp.x.Value() = ip(0);
    adp.y.Value() = ip(1);
    for (int i = 0; i < DIMR; i++)
      {
        adp.x.DValue(i) = jac(i,1);
        adp.y.DValue(i) = -jac(i,0);
      }
    return adp;
  }


  
  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {    
    TIP<DIM,AutoDiff<DIM>> pt = ip;
    static_cast<const FEL*> (this) -> 
      T_CalcShape (pt, // GetTIPHDiv<DIM> (ip),
                   SBLambda( [shape] (size_t nr, THDiv2Shape<DIM> val)  LAMBDA_INLINE
                             {
                               FlatVec<DIM> (&shape(nr,0)) = Vec<DIM> (val);
                             }));
  }
  
  
  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  CalcDivShape (const IntegrationPoint & ip, 
		SliceVector<> divshape) const
  {  
    TIP<DIM,AutoDiff<DIM>> pt = ip;    
    static_cast<const FEL*> (this) -> 
      T_CalcShape (pt, // GetTIPHDiv<DIM>(ip),
                   SBLambda( [divshape] (size_t nr, THDiv2DivShape<DIM> val) LAMBDA_INLINE
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
    static_cast<const FEL*> (this) -> 
      T_CalcShape (GetTIPHDiv(mip),
                   SBLambda([shape] (size_t nr, auto s) LAMBDA_INLINE
                            {
                              FlatVec<DIM> (&shape(nr,0)) = HDiv2ShapeNew(s);
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
    static_cast<const FEL*> (this) ->                 
      T_CalcShape (GetTIPHDiv(mip), 
                   SBLambda ([shapes] (size_t j, auto s)
                             {
                               auto vshape = HDiv2ShapeNew (s);
                               for (size_t k = 0; k < vshape.Size(); k++)
                                 shapes(j*DIM+k, 0) = vshape(k);
                             }));
  }
  
  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET>::
  CalcMappedShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                   BareSliceMatrix<SIMD<double>> shapes) const
  {
    Iterate<4-DIM>
      ([this,&bmir,shapes](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         if (bmir.DimSpace() == DIMSPACE)
           {
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
             for (size_t i = 0; i < mir.Size(); i++)
               {
                 auto shapesi = shapes.Col(i);
                 static_cast<const FEL*> (this) ->                 
                   T_CalcShape (GetTIPHDiv(mir[i]),
                                SBLambda ([shapesi] (size_t j, auto s) 
                                          {
                                            auto vshape = HDiv2ShapeNew (s); 
                                            for (size_t k = 0; k < vshape.Size(); k++)
                                              shapesi(j*vshape.Size()+k) = vshape(k);
                                          }));
               }
           }
       });
  }
  
  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET>::
  CalcMappedDivShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                      BareSliceMatrix<SIMD<double>> divshapes) const
  {
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        static_cast<const FEL*> (this) ->                 
          T_CalcShape (GetTIP(mir[i]),
                       SBLambda ([divshapes,i] (size_t j, THDiv2DivShape<DIM,SIMD<double>> val)
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
    Iterate<4-DIM>
      ([this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         if (bmir.DimSpace() == DIMSPACE)
           {
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
             for (size_t i = 0; i < mir.Size(); i++)
               {
                 Vec<DIMSPACE,SIMD<double>> sum(0.0);
                 static_cast<const FEL*> (this) ->         
                   T_CalcShape (GetTIPHDiv(mir[i]),
                                SBLambda ([coefs,&sum] (size_t j, auto s)
                                          {
                                            sum += coefs(j) * HDiv2ShapeNew(s); 
                                          }));
                 for (size_t k = 0; k < DIMSPACE; k++)
                   values(k,i) = sum(k);
               }
           }
       });
    
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  AddTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
            BareSliceVector<> coefs) const
  {
    Iterate<4-DIM>
      ([this,&bmir,values,coefs](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         if (bmir.DimSpace() == DIMSPACE)
           {
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
             for (size_t i = 0; i < mir.Size(); i++)
               {
                 Vec<DIMSPACE, SIMD<double>> vali;
                 for (int k = 0; k < DIMSPACE; k++)
                   vali(k) = values(k,i);
                 static_cast<const FEL*> (this) -> 
                   T_CalcShape (GetTIPHDiv(mir[i]),
                                SBLambda ([vali,coefs] (size_t j, auto s)
                                          {
                                            auto vshape = HDiv2ShapeNew(s); 
                                            SIMD<double> sum = 0.0;
                                            for (size_t k = 0; k < vshape.Size(); k++)
                                              sum += vali(k) * vshape(k);
                                            coefs(j) += HSum(sum);
                                          }));
               }
           }
       });
  }


  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  EvaluateDiv (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs,
               BareVector<SIMD<double>> values) const
  {
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        SIMD<double> sum(0.0);
        static_cast<const FEL*> (this) ->         
          T_CalcShape (GetTIP(mir[i]),
                       SBLambda ([=,&sum] (size_t j, THDiv2DivShape<DIM,SIMD<double>> divshape)
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
        SIMD<double> vali = values(i);
        static_cast<const FEL*> (this) -> 
          T_CalcShape (GetTIP(mir[i]),
                       SBLambda ([coefs,vali] (size_t j, THDiv2DivShape<DIM,SIMD<double>> divshape)
                                 {
                                   SIMD<double> simdshape = divshape.Get();
                                   coefs(j) += HSum(simdshape*vali);
                                 }));
      }
  }

  
#endif
}
