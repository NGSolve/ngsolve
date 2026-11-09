#ifndef FILE_THCURLFE_IMPL
#define FILE_THCURLFE_IMPL


#include "thcurlfe.hpp"
#include "recursive_pol.hpp" // for SBLambda

namespace ngfem
{


  
  /*******************************************/
  /* T_HCurlHOFiniteElement                  */
  /*******************************************/

  
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcShape (const IntegrationPoint & ip, BareSliceMatrix<> shape) const
  {
    this->T_CalcShape (GetTIPGrad<DIM>(ip), 
                       SBLambda ([shape](size_t i, auto s)
                                 { shape.Row(i) = s.Value(); }));
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET, SHAPES,BASE> :: 
  CalcCurlShape (const IntegrationPoint & ip, BareSliceMatrix<> shape) const
  {
    this->T_CalcShape (GetTIPGrad<DIM>(ip), 
                       SBLambda ([shape](size_t i, auto s) 
                                 { shape.Row(i) = s.CurlValue(); }));
  } 

#ifndef FASTCOMPILE

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET, SHAPES, BASE> :: 
  CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
                   BareSliceMatrix<> shape) const
  {
    Switch<4-DIM>
      (bmip.DimSpace()-DIM,[this,&bmip,shape](auto CODIM)
       {
         auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);
         this->T_CalcShape (GetTIP(mip),
                            SBLambda ([shape](size_t i, auto s) 
                                      {
                                        shape.Row(i) = s.Value();
                                      }));
       });
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcMappedShape (const BaseMappedIntegrationRule & bmir, 
                   BareSliceMatrix<> shapes) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,shapes](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const MappedIntegrationRule<DIM,DIM+CODIM.value>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           this->CalcMappedShape (mir[i], shapes.Cols(i*DIMSPACE,(i+1)*DIMSPACE));
       });
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET, SHAPES, BASE> :: 
  CalcMappedShape (const SIMD<BaseMappedIntegrationPoint> & bmip,
                   BareSliceMatrix<SIMD<double>> shape) const
  {
    Switch<4-DIM>
      (bmip.DimSpace()-DIM,[this,&bmip,shape](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mip = static_cast<const SIMD<MappedIntegrationPoint<DIM,DIM+CODIM.value>>&> (bmip);
         this->T_CalcShape (GetTIP(mip),
                            SBLambda ([shape,DIMSPACE=DIMSPACE](size_t i, auto s) 
                                      {
                                        shape.Col(0).Range(i*DIMSPACE, (i+1)*DIMSPACE) = s.Value();
                                      }));
       });
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcMappedShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                   BareSliceMatrix<SIMD<double>> shapes) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,shapes](auto CODIM)
       {
         static constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             auto shapei = shapes.Col(i);
             this->T_CalcShape
               (GetTIP(mir[i]),
                SBLambda ([shapei] (size_t j, auto s)
                          {
                            shapei.Range(j*DIMSPACE, (j+1)*DIMSPACE) = s.Value();
                          }));
               }
       });
  }


  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcMappedCurlShape (const BaseMappedIntegrationPoint & bmip,
                       BareSliceMatrix<> curlshape) const
  {
    auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM>&> (bmip);

    Vec<DIM,AutoDiff<DIM> > adp = mip;
    TIP<DIM,AutoDiff<DIM>> tip(adp, bmip.IP().FacetNr(), bmip.IP().VB());
    this->T_CalcShape (GetTIP(mip), 
                       SBLambda ([&](size_t i, auto s) 
                                 { 
                                   curlshape.Row(i) = s.CurlValue();
                                 }));
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcMappedCurlShape (const BaseMappedIntegrationRule & mir, 
                       BareSliceMatrix<> curlshape) const
  {
    for (int i = 0; i < mir.Size(); i++)
      CalcMappedCurlShape (mir[i], 
                           curlshape.Cols(DIM_CURL_(DIM)*i, DIM_CURL_(DIM)*(i+1)));
  }    


  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcMappedCurlShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                       BareSliceMatrix<SIMD<double>> shapes) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,shapes](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         constexpr int DIM_CURL = DIM_CURL_(DIMSPACE);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             auto shapei = shapes.Col(i);
             this->T_CalcShape (GetTIP(mir[i]),
                                SBLambda ([shapei,DIM_CURL=DIM_CURL] (size_t j, auto s)
                                          {
                                            shapei.Range(j*DIM_CURL, (j+1)*DIM_CURL) = s.CurlValue();
                                          }));
           }
       });
  }

  

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  auto T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  EvaluateCurlShape (const IntegrationPoint & ip, 
                     BareSliceVector<double> x,
                     LocalHeap & lh) const -> Vec<DIM_CURL_(DIM)>
  {
    // Vec<DIM, AutoDiff<DIM> > adp = ip;
    // TIP<DIM,AutoDiff<DIM>> tip(adp);
    
    Vec<DIM_CURL_(DIM)> sum = 0.0;
    this->T_CalcShape (GetTIPGrad<DIM>(ip), 
                       SBLambda ([&sum, x](size_t i, auto s) 
                                 { sum += x(i) * s.CurlValue(); }));
    return sum;
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  EvaluateCurl (const IntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<> curl) const
  {
    LocalHeapMem<10000> lhdummy("evalcurl-heap");
    for (int i = 0; i < ir.Size(); i++)
      curl.Row(i) = EvaluateCurlShape (ir[i], coefs, lhdummy);
  }

  
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  Evaluate (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs,
            BareSliceMatrix<SIMD<double>> values) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIMSPACE,SIMD<double>> sum(0.0);
             this->T_CalcShape (GetTIP(mir[i]),
                                SBLambda ([&sum,coefs] (size_t j, auto shape)
                                          {
                                            sum += coefs(j) * shape.Value();
                                          }));
             values.Col(i).Range(DIMSPACE) = sum;
           }
       });
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  Evaluate (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<Complex> coefs,
            BareSliceMatrix<SIMD<Complex>> values) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIMSPACE,SIMD<Complex>> sum = SIMD<Complex>(0.0);
             this->T_CalcShape (GetTIP(mir[i]),
                                SBLambda ([&sum,coefs] (size_t j, auto shape)
                                          {
                                            sum += coefs(j) * shape.Value();
                                          }));
             values.Col(i).Range(DIMSPACE) = sum;
           }
       });
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  EvaluateCurl (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         constexpr int DIM_CURL = DIM_CURL_(DIMSPACE);
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIM_CURL,SIMD<double>> sum(0.0);            
             this->T_CalcShape (GetTIP(mir[i]),
                                SBLambda ([coefs,&sum] (size_t j, auto shape)
                                          {
                                            sum += coefs(j) * shape.CurlValue();
                                          }));
             values.Col(i).Range(DIM_CURL) = sum;
           }
       });
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  EvaluateCurl (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         constexpr int DIM_CURL = DIM_CURL_(DIMSPACE);
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIM_CURL,SIMD<Complex>> sum = SIMD<Complex>(0.0);
             this->T_CalcShape (GetTIP(mir[i]),
                                SBLambda ([coefs, &sum] (size_t j, auto shape)
                                          {
                                            sum += coefs(j) * shape.CurlValue();
                                          }));
             values.Col(i).Range(DIM_CURL) = sum;
           }
       });
  }

  
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  AddTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
            BareSliceVector<> coefs) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIMSPACE,SIMD<double>> vali = values.Col(i);
             this->T_CalcShape (GetTIP(mir[i]),
                                SBLambda ([vali,coefs] (size_t j, auto s)
                                          {
                                            coefs(j) += HSum(InnerProduct(s.Value(), vali));
                                          }));
           }
       });
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  AddTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<Complex>> values,
            BareSliceVector<Complex> coefs) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIMSPACE,SIMD<Complex>> vali = values.Col(i);
             this->T_CalcShape (GetTIP(mir[i]),
                                SBLambda ([vali,coefs] (size_t j, auto s)
                                          {
                                            coefs(j) += HSum(InnerProduct(s.Value(), vali));
                                          }));
           }
       });
  }

  
  
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  AddCurlTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
                BareSliceVector<> coefs) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         constexpr int DIM_CURL = DIM_CURL_(DIMSPACE);
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIM_CURL,SIMD<double>> vali = values.Col(i);
             this->T_CalcShape (GetTIP(mir[i]),
                                SBLambda ([vali,coefs] (size_t j, auto s)
                                          {
                                            coefs(j) += HSum(InnerProduct(s.CurlValue(), vali));
                                          }));
           }
       });
  }
  
  
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  AddCurlTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<Complex>> values,
                BareSliceVector<Complex> coefs) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         constexpr int DIM_CURL = DIM_CURL_(DIMSPACE);
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIM_CURL,SIMD<Complex>> vali = values.Col(i);
             this->T_CalcShape (GetTIP(mir[i]),
                                SBLambda ([vali,coefs] (size_t j, auto s)
                                          {
                                            coefs(j) += HSum(InnerProduct(s.CurlValue(), vali));
                                          }));
           }
       });
  }

  // reference point, outer AutoDiff differentiates the shape
  template <int DIM, typename SCAL, typename IP>
  INLINE auto GetTIPGradNested (const IP & ip)
  {
    typedef AutoDiff<DIM,SCAL> TI;
    typedef AutoDiff<DIM,TI> T;
    auto coord = [&ip] (int m)
    {
      T c;
      c.Value() = TI(SCAL(ip(m)), m);
      for (int i = 0; i < DIM; i++)
        c.DValue(i) = TI(SCAL(i == m ? 1.0 : 0.0));
      return c;
    };
    if constexpr (DIM == 1)
      return TIP<1,T> (coord(0), ip.FacetNr(), ip.VB());
    else if constexpr (DIM == 2)
      return TIP<2,T> (coord(0), coord(1), ip.FacetNr(), ip.VB());
    else
      return TIP<3,T> (coord(0), coord(1), coord(2), ip.FacetNr(), ip.VB());
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> ::
  CalcDShape (const IntegrationPoint & ip, BareSliceMatrix<> dshape) const
  {
    if constexpr (DIM == 0) return; else
    this->T_CalcShape (GetTIPGradNested<DIM,double> (ip),
                       SBLambda ([dshape] (size_t nr, auto s)
                                 {
                                   auto val = s.Value();
                                   for (int j = 0; j < DIM; j++)
                                     for (int k = 0; k < DIM; k++)
                                       dshape(nr, j*DIM+k) = val(j).DValue(k);
                                 }));
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> ::
  CalcDShape (const SIMD<IntegrationPoint> & ip, BareSliceVector<SIMD<double>> dshape) const
  {
    if constexpr (DIM == 0) return; else
    this->T_CalcShape (GetTIPGradNested<DIM,SIMD<double>> (ip),
                       SBLambda ([dshape] (size_t nr, auto s)
                                 {
                                   auto val = s.Value();
                                   size_t base = nr*(DIM*DIM+DIM);
                                   for (int j = 0; j < DIM; j++)
                                     {
                                       for (int k = 0; k < DIM; k++)
                                         dshape(base+j*DIM+k) = val(j).DValue(k);
                                       dshape(base+DIM*DIM+j) = val(j).Value();
                                     }
                                 }));
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> ::
  EvaluateGradRef (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs,
                   BareSliceMatrix<SIMD<double>> values) const
  {
    if constexpr (DIM == 0) return; else
    for (size_t i = 0; i < ir.Size(); i++)
      {
        Vec<DIM,AutoDiff<DIM,SIMD<double>>> sum = AutoDiff<DIM,SIMD<double>> (SIMD<double>(0.0));
        this->T_CalcShape (GetTIPGradNested<DIM,SIMD<double>> (ir[i]),
                           SBLambda ([coefs,&sum] (size_t j, auto s)
                                     {
                                       auto val = s.Value();
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

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> ::
  AddGradRefTrans (const SIMD_IntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                   BareSliceVector<> coefs) const
  {
    if constexpr (DIM == 0) return; else
    for (size_t i = 0; i < ir.Size(); i++)
      {
        Vec<DIM*DIM+DIM,SIMD<double>> vali = values.Col(i);
        this->T_CalcShape (GetTIPGradNested<DIM,SIMD<double>> (ir[i]),
                           SBLambda ([vali,coefs] (size_t nr, auto s)
                                     {
                                       auto val = s.Value();
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

#endif
