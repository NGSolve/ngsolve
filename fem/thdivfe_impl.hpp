#include <thdivfe.hpp>
#include "hdivfe_utils.hpp"

namespace ngfem
{
  
  template<int DIM>
  INLINE auto GetTIPHDiv (const IntegrationPoint & ip);


  template<>
  INLINE auto GetTIPHDiv<2> (const IntegrationPoint & ip)
  {
    TIP<2,AutoDiff<2>> tip { ip, ip.FacetNr(), ip.VB() };
    tip.x.DValue(0) = 0;
    tip.x.DValue(1) = 1;
    tip.y.DValue(0) = -1;
    tip.y.DValue(1) = 0;
    return tip;
  }

  template<>
  INLINE auto GetTIPHDiv<3> (const IntegrationPoint & ip)
  {
    TIP<3,AutoDiff<3>> tip = ip;
    return tip;
  }  

  
  
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
    TIP<2,AutoDiff<DIMR>> adp(mip.IP().FacetNr(), mip.IP().VB());    
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
    // adp.facetnr = mip.IP().FacetNr();
    // adp.vb = mip.IP().VB();
    return adp;
  }
  
  template <int DIMR>
  INLINE auto GetTIPHDiv (const SIMD<MappedIntegrationPoint<2,DIMR>> & mip)
  {
    TIP<2,AutoDiff<DIMR,SIMD<double>>> adp(mip.IP().FacetNr(), mip.IP().VB());
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
    // adp.facetnr = mip.IP().FacetNr();
    // adp.vb = mip.IP().VB();
    return adp;
  }


  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  CalcShape (const IntegrationPoint & ip, 
	     BareSliceMatrix<> shape) const
  {
    /*
    static_cast<const FEL*> (this) -> 
      T_CalcShape (GetTIPGrad<DIM>(ip), 
                   SBLambda( [shape] (size_t nr, THDiv2Shape<DIM> val)  LAMBDA_INLINE
                             {
                               shape.Row(nr) = Vec<DIM> (val);
                             }));
    */

    static_cast<const FEL*> (this) -> 
      T_CalcShape (GetTIPHDiv<DIM>(ip), 
                   SBLambda( [shape] (size_t nr, auto val)  LAMBDA_INLINE
                   {
                     shape.Row(nr) = HDiv2ShapeNew(val);
                   }));
  }
  
  
  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  CalcDivShape (const IntegrationPoint & ip, 
		BareSliceVector<> divshape) const
  {  
    static_cast<const FEL*> (this) -> 
      T_CalcShape (GetTIPGrad<DIM>(ip), 
                   SBLambda( [divshape] (size_t nr, THDiv2DivShape<DIM> val) LAMBDA_INLINE
                             {
                               divshape(nr) = val;
                             }));
  }

#ifndef FASTCOMPILE

  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
                   BareSliceMatrix<> shape) const
  {
    Iterate<4-DIM>
      ([this,&bmip,shape](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         if (bmip.DimSpace() == DIMSPACE)
           {
             auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);
             static_cast<const FEL*> (this) ->
               T_CalcShape (GetTIPHDiv(mip),
                            SBLambda ([shape](size_t nr, auto s) 
                                      {
                                        shape.Row(nr) = HDiv2ShapeNew(s);
                                      }));
           }
       });
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET>::
  CalcMappedShape (const BaseMappedIntegrationRule & bmir, 
                   BareSliceMatrix<> shapes) const
  {
    Iterate<4-DIM>
      ([this,&bmir,shapes](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         if (bmir.DimSpace() == DIMSPACE)
           {
             auto & mir = static_cast<const MappedIntegrationRule<DIM,DIM+CODIM.value>&> (bmir);
             for (size_t i = 0; i < mir.Size(); i++)
               this->CalcMappedShape (mir[i], shapes.Cols(i*DIMSPACE,(i+1)*DIMSPACE));
           }
       });
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
                                            shapesi.Range(j*vshape.Size(), (j+1)*vshape.Size()) = vshape;
                                          }));
               }
           }
       });
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET>::
  CalcMappedNormalShape (const SIMD_BaseMappedIntegrationRule & bmir, 
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
                 auto nv = mir[i].GetNV();
                 auto shapesi = shapes.Col(i);
                 static_cast<const FEL*> (this) ->                 
                   T_CalcShape (GetTIPHDiv(mir[i]),
                                SBLambda ([shapesi, nv] (size_t j, auto s) 
                                          {
                                            auto vshape = HDiv2ShapeNew (s); 
                                            // shapesi.Range(j*vshape.Size(), (j+1)*vshape.Size()) = vshape;
                                            shapesi(j) = InnerProduct(nv, vshape);
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
        auto divshapesi = divshapes.Col(i);
        static_cast<const FEL*> (this) ->                 
          T_CalcShape (GetTIP(mir[i]),
                       SBLambda ([divshapesi] (size_t j, THDiv2DivShape<DIM,SIMD<double>> val)
                                 {
                                   divshapesi(j) = val;
                                 }));
      }
  }

  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, 
	    BareSliceMatrix<> vals) const  
  {
    /*
    for (size_t i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM, AutoDiff<DIM>> adp = ir[i]; 

        Vec<DIM> sum = 0;
        static_cast<const FEL*> (this) -> 
          T_CalcShape (TIP<DIM,AutoDiff<DIM>>(adp, ir[i].FacetNr(), ir[i].VB()),
                       SBLambda([coefs,&sum] (size_t j, THDiv2Shape<DIM> vshape)
                                {
                                  sum += coefs(j) * Vec<DIM> (vshape);
                                }));
        vals.Row(i) = sum;
      }
    */

    for (size_t i = 0; i < ir.GetNIP(); i++)
      {
        // Vec<DIM, AutoDiff<DIM>> adp = ir[i]; 

        Vec<DIM> sum = 0;
        static_cast<const FEL*> (this) -> 
          T_CalcShape (GetTIPHDiv<DIM>(ir[i]), 
                       SBLambda([coefs,&sum] (size_t j, auto s)
                                {
                                  sum += coefs(j) * HDiv2ShapeNew(s);
                                }));
        vals.Row(i) = sum;
      }
  }


  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  EvaluateTrans (const IntegrationRule & ir, 
                 BareSliceMatrix<> vals,
                 FlatVector<double> coefs) const  
  {
    /*
    coefs = 0;
    for (size_t i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM, AutoDiff<DIM>> adp = ir[i]; 

        Vec<DIM> val = vals.Row(i);
        static_cast<const FEL*> (this) -> 
          T_CalcShape (TIP<DIM,AutoDiff<DIM>>(adp, ir[i].FacetNr(), ir[i].VB()),
                       SBLambda([coefs,val] (size_t j, THDiv2Shape<DIM> vshape)
                                {
                                  coefs(j) += InnerProduct (val, Vec<DIM> (vshape));
                                }));
      }
    */

    coefs = 0;
    for (size_t i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM> val = vals.Row(i);
        static_cast<const FEL*> (this) -> 
          T_CalcShape (GetTIPHDiv<DIM>(ir[i]),
                       SBLambda([coefs,val] (size_t j, auto s)
                                {
                                  coefs(j) += InnerProduct (val, HDiv2ShapeNew(s));
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
                 values.Col(i).Range(DIMSPACE) = sum;
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
                 Vec<DIMSPACE, SIMD<double>> vali = values.Col(i);
                 // for (int k = 0; k < DIMSPACE; k++)
                 // vali(k) = values(k,i);
                 static_cast<const FEL*> (this) -> 
                   T_CalcShape (GetTIPHDiv(mir[i]),
                                SBLambda ([vali,coefs] (size_t j, auto s)
                                          {
                                            auto vshape = HDiv2ShapeNew(s); 
                                            coefs(j) += HSum(InnerProduct(vali,vshape));
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
                                   // SIMD<double> simdshape = divshape.Get();
                                   // sum += coefs(j) * simdshape;
                                   sum += coefs(j) * divshape.Get();
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
                                   coefs(j) += HSum(divshape.Get()*vali);
                                 }));
      }
  }

  


#endif
}
