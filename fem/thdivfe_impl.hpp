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
    /*
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
    */
    Iterate<4-DIM>
      ([this,&bmir,divshapes](auto CODIM)
      {
        constexpr int DIMSPACE = DIM+CODIM.value;
        if (bmir.DimSpace() == DIMSPACE)
          {
            // cout << "use new div, dim " << DIMSPACE << endl;
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
             for (size_t i = 0; i < mir.Size(); i++)
               {
                 auto divshapesi = divshapes.Col(i);
                 auto invJ = 1.0/mir[i].GetJacobiDet();
                 static_cast<const FEL*> (this) ->                 
                   T_CalcShape (GetTIPGrad<DIM>(mir[i].IP()),
                                SBLambda ([divshapesi,invJ] (size_t j, THDiv2DivShape<DIM,SIMD<double>> val)
                                          {
                                            divshapesi(j) = invJ*val.Get(); 
                                          }));
               }
           }
       });
   
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
    /*
    cout << "evaldiv, simd" << endl
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
    */

    /*
    Iterate<4-DIM>
      ([this,&bmir,coefs,values](auto CODIM)
      {
        constexpr int DIMSPACE = DIM+CODIM.value;
        if (bmir.DimSpace() == DIMSPACE)
          {
            auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
            for (size_t i = 0; i < mir.Size(); i++)
               {
                 SIMD<double> sum(0.0);                 
                 static_cast<const FEL*> (this) ->                 
                   T_CalcShape (GetTIPGrad<DIM>(mir[i].IP()),
                                SBLambda ([=,&sum] (size_t j, THDiv2DivShape<DIM,SIMD<double>> divshape)
                                          {
                                            // auto vshape = HDiv2DivShapeNew (s); 
                                            // shapesi.Range(j*vshape.Size(), (j+1)*vshape.Size()) = vshape;
                                            // divshapesi(j) = HDiv2DivShapeNew (s);
                                            sum += coefs(j) * divshape.Get();                                            
                                          }));
                 values(i) = 1/mir[i].GetJacobiDet() * sum;
               }
           }
       });
    */

    auto & mir = bmir;
    for (size_t i = 0; i < mir.Size(); i++)
      {
        SIMD<double> sum(0.0);                 
        static_cast<const FEL*> (this) ->                 
          T_CalcShape (GetTIPGrad<DIM>(mir[i].IP()),
                       SBLambda ([=,&sum] (size_t j, THDiv2DivShape<DIM,SIMD<double>> divshape)
                       {
                         sum += coefs(j) * divshape.Get();                                            
                       }));
        values(i) = 1/mir[i].GetJacobiDet() * sum;
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

  


  // reference point, outer AutoDiff differentiates the shape
  template <int DIM, typename SCAL, typename IP>
  INLINE auto GetTIPHDivNested (const IP & ip)
  {
    typedef AutoDiff<DIM,SCAL> TI;
    typedef AutoDiff<DIM,TI> T;
    auto coord = [&ip] (int m, Vec<DIM> d)
    {
      T c;
      c.Value() = TI(SCAL(ip(m)), m);
      for (int i = 0; i < DIM; i++)
        c.DValue(i) = TI(SCAL(d(i)));
      return c;
    };
    if constexpr (DIM == 2)
      return TIP<2,T> (coord(0, Vec<2>(0,1)), coord(1, Vec<2>(-1,0)), ip.FacetNr(), ip.VB());
    else
      return TIP<3,T> (coord(0, Vec<3>(1,0,0)), coord(1, Vec<3>(0,1,0)), coord(2, Vec<3>(0,0,1)),
                       ip.FacetNr(), ip.VB());
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> ::
  CalcDShape (const IntegrationPoint & ip, BareSliceMatrix<> dshape) const
  {
    if constexpr (DIM == 1)
      HDivFiniteElement<DIM>::CalcDShape (ip, dshape);
    else
      static_cast<const FEL*> (this) ->
        T_CalcShape (GetTIPHDivNested<DIM,double> (ip),
                     SBLambda ([dshape] (size_t nr, auto s)
                               {
                                 auto val = HDiv2ShapeNew(s);
                                 for (int j = 0; j < DIM; j++)
                                   for (int k = 0; k < DIM; k++)
                                     dshape(nr, j*DIM+k) = val(j).DValue(k);
                               }));
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> ::
  CalcDShape (const SIMD<IntegrationPoint> & ip, BareSliceVector<SIMD<double>> dshape) const
  {
    if constexpr (DIM == 1)
      HDivFiniteElement<DIM>::CalcDShape (ip, dshape);
    else
      static_cast<const FEL*> (this) ->
        T_CalcShape (GetTIPHDivNested<DIM,SIMD<double>> (ip),
                     SBLambda ([dshape] (size_t nr, auto s)
                               {
                                 auto val = HDiv2ShapeNew(s);
                                 size_t base = nr*(DIM*DIM+DIM);
                                 for (int j = 0; j < DIM; j++)
                                   {
                                     for (int k = 0; k < DIM; k++)
                                       dshape(base+j*DIM+k) = val(j).DValue(k);
                                     dshape(base+DIM*DIM+j) = val(j).Value();
                                   }
                               }));
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> ::
  EvaluateGradRef (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs,
                   BareSliceMatrix<SIMD<double>> values) const
  {
    if constexpr (DIM == 1)
      HDivFiniteElement<DIM>::EvaluateGradRef (ir, coefs, values);
    else
      for (size_t i = 0; i < ir.Size(); i++)
        {
          Vec<DIM,AutoDiff<DIM,SIMD<double>>> sum = AutoDiff<DIM,SIMD<double>> (SIMD<double>(0.0));
          static_cast<const FEL*> (this) ->
            T_CalcShape (GetTIPHDivNested<DIM,SIMD<double>> (ir[i]),
                         SBLambda ([coefs,&sum] (size_t j, auto s)
                                   {
                                     auto val = HDiv2ShapeNew(s);
                                     for (int k = 0; k < DIM; k++)
                                       sum(k) += coefs(j) * val(k);
                                   }));
          for (int j = 0; j < DIM; j++)
            {
              for (int k = 0; k < DIM; k++)
                values(j*DIM+k, i) = sum(j).DValue(k);
              values(DIM*DIM+j, i) = sum(j).Value();
            }
        }
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> ::
  AddGradRefTrans (const SIMD_IntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                   BareSliceVector<> coefs) const
  {
    if constexpr (DIM == 1)
      HDivFiniteElement<DIM>::AddGradRefTrans (ir, values, coefs);
    else
      for (size_t i = 0; i < ir.Size(); i++)
        {
          Vec<DIM*DIM+DIM,SIMD<double>> vali = values.Col(i);
          static_cast<const FEL*> (this) ->
            T_CalcShape (GetTIPHDivNested<DIM,SIMD<double>> (ir[i]),
                         SBLambda ([vali,coefs] (size_t nr, auto s)
                                   {
                                     auto val = HDiv2ShapeNew(s);
                                     SIMD<double> sum(0.0);
                                     for (int j = 0; j < DIM; j++)
                                       {
                                         for (int k = 0; k < DIM; k++)
                                           sum += vali(j*DIM+k) * val(j).DValue(k);
                                         sum += vali(DIM*DIM+j) * val(j).Value();
                                       }
                                     coefs(nr) += HSum(sum);
                                   }));
        }
  }


#endif
}
