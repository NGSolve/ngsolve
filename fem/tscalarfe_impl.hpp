#ifdef __CUDA_ARCH__
#include "tscalarfe_impl_cuda.hpp"
#endif


#ifndef FILE_TSCALARFE_IMPL
#define FILE_TSCALARFE_IMPL


#include "tscalarfe.hpp"
#include "recursive_pol.hpp"
#include "shapefunction_utils.hpp"

namespace ngfem
{

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcShape (const IntegrationPoint & ip, BareSliceVector<> shape) const
  {
    T_CalcShape (GetTIP<DIM>(ip), shape);    
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcDShape (const IntegrationPoint & ip, 
              BareSliceMatrix<> dshape) const
  {
    T_CalcShape (GetTIPGrad<DIM> (ip),
                 SBLambda ([dshape] (int i, auto shape)
                 { dshape.Row(i) = ngfem::GetGradient(shape); }));
  }

#ifndef FASTCOMPILE

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcShape (const IntegrationRule & ir, BareSliceMatrix<> shape) const
  {
    for (int i = 0; i < ir.Size(); i++)
      T_CalcShape (GetTIP<DIM>(ir[i]), shape.Col(i));        
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcShape (const SIMD_IntegrationRule & ir, BareSliceMatrix<SIMD<double>> shapes) const
  {
    /*
    for (size_t i = 0; i < ir.Size(); i++)
      T_CalcShape (GetTIP<DIM>(ir[i]),
                   SBLambda([&](size_t j, SIMD<double> shape)
                            { shapes(j,i) = shape; } ));
    */
    for (size_t i = 0; i < ir.Size(); i++)
      T_CalcShape (GetTIP<DIM>(ir[i]), shapes.Col(i));
  }
  
  
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  double T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const IntegrationPoint & ip, BareSliceVector<double> x) const
  {
    double sum = 0;
    T_CalcShape (GetTIP<DIM>(ip),
                 SBLambda ([x,&sum](size_t i, double val) { sum += x(i)*val; } ));
    return sum;
  }  


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const IntegrationRule & ir, BareSliceVector<double> coefs, BareSliceVector<double> vals) const
  {
    for (size_t i = 0; i < ir.GetNIP(); i++)
      {
        double sum = 0;
        T_CalcShape (GetTIP<DIM>(ir[i]),
                     SBLambda ( [coefs,&sum](size_t j, double shape) { sum += coefs(j)*shape; } ));
        vals(i) = sum;
      }
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs, BareVector<SIMD<double>> values) const
  {
    FlatArray<SIMD<IntegrationPoint>> hir = ir;
    size_t i = 0;
    for ( ; i+2 <= hir.Size(); i+=2)
      {
        MultiSIMD<2,double> sum = 0;
        auto tip1 = GetTIP<DIM>(hir[i]);
        auto tip2 = GetTIP<DIM>(hir[i+1]);
        TIP<DIM,MultiSIMD<2,double>> tip(tip1,tip2);
        
        double * pcoefs = coefs.Data();
        size_t dist = coefs.Dist();
        T_CalcShape (tip, 
                     SBLambda ( [&](size_t j, MultiSIMD<2,double> shape)
                                {
                                  // sum += *pcoefs * shape;
                                  sum = FMA(MultiSIMD<2,double>(*pcoefs), shape, sum);
                                  pcoefs += dist; }
                                ));
        
        // std::tie(values(i), values(i+1)) = sum;
        values(i) = sum.Lo();
        values(i+1) = sum.Hi();
      }

    if (i < hir.Size())
      {
        SIMD<double> sum = 0;
        double * pcoefs = coefs.Data();
        size_t dist = coefs.Dist();
        T_CalcShape (GetTIP<DIM>(hir[i]),
                     SBLambda ( [&](int j, SIMD<double> shape)
                                { sum += (*pcoefs)*shape; pcoefs += dist; } ));
        values(i) = sum;
      }
  }
    

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const SIMD_IntegrationRule & ir,
            SliceMatrix<> coefs,
            BareSliceMatrix<SIMD<double>> values) const
  {
    FlatArray<SIMD<IntegrationPoint>> hir = ir;    
    size_t j = 0;
    for ( ; j+4 <= coefs.Width(); j+=4)
      {
        for (size_t i = 0; i < hir.Size(); i++)
          {
            SIMD<double> sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
            // double * pcoefs = &coefs(j);
            double * pcoefs = coefs.Addr(0,j);
            size_t dist = coefs.Dist();
            T_CalcShape (GetTIP<DIM>(hir[i]), 
                         SBLambda ( [&pcoefs, dist, &sum1, &sum2, &sum3, &sum4](int j, SIMD<double> shape)
                                    {
                                      sum1 += pcoefs[0]*shape;
                                      sum2 += pcoefs[1]*shape;
                                      sum3 += pcoefs[2]*shape;
                                      sum4 += pcoefs[3]*shape;
                                      pcoefs += dist;
                                    } ));
            values(j,i) = sum1;
            values(j+1,i) = sum2;
            values(j+2,i) = sum3;
            values(j+3,i) = sum4;
          }
      }
    switch (coefs.Width()&3)
      {
      case 0: break;
      case 1: Evaluate (ir, coefs.Col(j), values.Row(j)); break;
      case 2:
        {
          for (size_t i = 0; i < hir.Size(); i++)
            {
              SIMD<double> sum1 = 0, sum2 = 0;
              // double * pcoefs = &coefs(j);
              double * pcoefs = coefs.Addr(0,j);
              size_t dist = coefs.Dist();
              T_CalcShape (GetTIP<DIM>(hir[i]), 
                           SBLambda ( [&pcoefs, dist,&sum1, &sum2](int j, SIMD<double> shape)
                                      {
                                        sum1 += pcoefs[0]*shape;
                                        sum2 += pcoefs[1]*shape;
                                        pcoefs += dist;
                                      } ));
              values(j,i) = sum1;
              values(j+1,i) = sum2;
            }
          break;
        case 3:
          {
            for (size_t i = 0; i < hir.Size(); i++)
              {
                SIMD<double> sum1 = 0, sum2 = 0, sum3 = 0;
                // double * pcoefs = &coefs(j);
                double * pcoefs = coefs.Addr(0,j);
                size_t dist = coefs.Dist();
                T_CalcShape (GetTIP<DIM>(hir[i]), 
                             SBLambda ( [&pcoefs, dist, &sum1,&sum2,&sum3](int j, SIMD<double> shape)
                                        {
                                          sum1 += pcoefs[0]*shape;
                                          sum2 += pcoefs[1]*shape;
                                          sum3 += pcoefs[2]*shape;
                                          pcoefs += dist;
                                        } ));
                values(j,i) = sum1;
                values(j+1,i) = sum2;
                values(j+2,i) = sum3;
              }
            break;
          }
        }
      }
        /*
    for ( ; j < coefs.Width(); j++)
      Evaluate (ir, coefs.Col(j), values.Row(j));
        */
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const IntegrationRule & ir, SliceMatrix<> coefs, BareSliceMatrix<> values) const
  {
    for (size_t i = 0; i < ir.GetNIP(); i++)
      {
        auto hrow = values.Row(i).Range(coefs.Width());
        hrow = 0.0;
        T_CalcShape (GetTIP<DIM>(ir[i]), 
                     SBLambda ( [&](size_t j, double shape) 
                                { 
                                  hrow += shape * coefs.Row(j); 
                                }));
      }
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateTrans (const IntegrationRule & ir, BareSliceVector<> vals, BareSliceVector<double> coefs) const
  {
    coefs.Range(0,ndof) = 0.0;
    for (size_t i = 0; i < ir.GetNIP(); i++)
      {
        double vali = vals(i);
        T_CalcShape (GetTIP<DIM>(ir[i]), 
                     SBLambda ( [coefs, vali](size_t j, double shape) 
                                { coefs(j) += vali*shape; } ));
      }
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  AddTrans (const SIMD_IntegrationRule & ir, BareVector<SIMD<double>> values,
            BareSliceVector<> coefs) const
  {
    FlatArray<SIMD<IntegrationPoint>> hir = ir;
    /*
    for (int i = 0; i < hir.Size(); i++)
      {
        Vec<DIM,SIMD<double>> pt = hir[i];
        SIMD<double> val = values.Get(i);
        T_CalcShape (&pt(0), SBLambda ( [&](int j, SIMD<double> shape) { coefs(j) += HSum(val*shape); } ));
      }
    */

    size_t i = 0;
    for ( ; i+2 <= hir.Size(); i+=2)
      {
        TIP<DIM,SIMD<double>> tip1 = hir[i].TIp<DIM>();
        TIP<DIM,SIMD<double>> tip2 = hir[i+1].TIp<DIM>();
        TIP<DIM,MultiSIMD<2,double>> tip(tip1,tip2);

        MultiSIMD<2,double> val (values(i), values(i+1));

        double * pcoefs = coefs.Data();
        size_t dist = coefs.Dist();
        /*
        T_CalcShape (tip, 
                     SBLambda
                     ([&](int j, MultiSIMD<2,double> shape)
                      { *pcoefs += HSum(val*shape); pcoefs += dist; }
                      ));
        */
        T_CalcShape (tip, 
                     SBLambdaDuo
                     ([&](int j, MultiSIMD<2,double> shape)
                      {
                        *pcoefs += HSum(val*shape); pcoefs += dist;
                      },
                      [&](int j, MultiSIMD<2,double> shape, int j2, MultiSIMD<2,double> shape2)
                      {
                        auto v2 = HSum(val*shape, val*shape2);
                        *pcoefs += get<0>(v2); pcoefs += dist;
                        *pcoefs += get<1>(v2); pcoefs += dist;
                      }
                      ));        
      }

    for ( ; i < hir.Size(); i++)
      {
        TIP<DIM,SIMD<double>> tip = hir[i].TIp<DIM>();
        SIMD<double> val (values(i));

        double * pcoefs = coefs.Data();
        size_t dist = coefs.Dist();
        T_CalcShape (tip, 
                     SBLambdaDuo
                     ([&](int j, SIMD<double> shape)
                      {
                        *pcoefs += HSum(val*shape); pcoefs += dist;
                      },
                      [&](int j, SIMD<double> shape, int j2, SIMD<double> shape2)
                      {
                        auto v2 = HSum(val*shape, val*shape2);
                        *pcoefs += get<0>(v2); pcoefs += dist;
                        *pcoefs += get<1>(v2); pcoefs += dist;
                      }
                      ));
      }


    
    
    /*
    for (int i = 0; i < hir.Size(); i+=3)
      {
        Vec<DIM,SIMD<double>> pt1 = hir[i];
        Vec<DIM,SIMD<double>> pt2 = hir[(i+1 < hir.Size()) ? i+1 : i];
        Vec<DIM,SIMD<double>> pt3 = hir[(i+2 < hir.Size()) ? i+2 : i];
        
        Vec<DIM,MultiSIMD<3,double>> pt;
        for (int i = 0; i < DIM; i++)
          pt(i) = MultiSIMD<3,double> (pt1(i), pt2(i), pt3(i));          
        MultiSIMD<3,double> val (values.Get(i),
                                 i+1 < hir.Size() ? values.Get(i+1) : SIMD<double> (0.0),
                                 i+2 < hir.Size() ? values.Get(i+2) : SIMD<double> (0.0));

        // T_CalcShape (&pt(0), SBLambda ( [&](int j, MultiSIMD<3,double> shape) { coefs(j) += HSum(val*shape); } ));

        double * pcoefs = &coefs(0);
        size_t dist = coefs.Dist();
        T_CalcShape (TIP<DIM,MultiSIMD<3,double>> (pt),
                     SBLambda ( [&](int j, MultiSIMD<3,double> shape)
                                { *pcoefs += HSum(val*shape); pcoefs += dist; } ));
        
      }
    */
  }

  // #endif // FASTCOMPILE
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  AddDualTrans (const SIMD_IntegrationRule & ir, BareVector<SIMD<double>> values,
            BareSliceVector<> coefs) const
  {
    FlatArray<SIMD<IntegrationPoint>> hir = ir;
    for (int i = 0; i < hir.Size(); i++)
      {
        TIP<DIM,SIMD<double>> tip = hir[i].TIp<DIM>();
        SIMD<double> val = values(i);
        static_cast<const FEL*> (this)->        
          T_CalcDualShape (tip, SBLambda ( [&](int j, SIMD<double> shape) { coefs(j) += HSum(val*shape); } ));
      }
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  AddDualTrans (const IntegrationRule & ir, BareSliceVector<double> values,
            BareSliceVector<> coefs) const
  {
    FlatArray<IntegrationPoint> hir = ir;
    for (int i = 0; i < hir.Size(); i++)
      {
        TIP<DIM,double> tip = hir[i].TIp<DIM>();
        double val = values(i);
        static_cast<const FEL*> (this)->                
          T_CalcDualShape (tip, SBLambda ( [&](int j, double shape) { coefs(j) += val*shape; } ));
      }
  }


  // #ifndef FASTCOMPILE
  
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  AddTrans (const SIMD_IntegrationRule & ir,
            BareSliceMatrix<SIMD<double>> values,
            SliceMatrix<> coefs) const
  {
    FlatArray<SIMD<IntegrationPoint>> hir = ir;    
    size_t j = 0;
    for ( ; j+4 <= coefs.Width(); j+=4)
      {
        for (size_t i = 0; i < hir.Size(); i++)
          {
            TIP<DIM,SIMD<double>> pt = hir[i].TIp<DIM>();
            SIMD<double> val1 = values(j,i);
            SIMD<double> val2 = values(j+1,i);
            SIMD<double> val3 = values(j+2,i);
            SIMD<double> val4 = values(j+3,i);
            double * pcoefs = &coefs(j);
            size_t dist = coefs.Dist();
            T_CalcShape (pt, 
                         SBLambda ( [&](int j, SIMD<double> shape)
                                    {
                                      auto val = HSum(shape*val1, shape*val2, shape*val3, shape*val4);
                                      val += SIMD<double,4> (pcoefs);
                                      // _mm256_storeu_pd (pcoefs, val.Data());
                                      val.Store(pcoefs);
                                      pcoefs += dist;
                                    } ));
          }
      }
    switch (coefs.Width()&3)
      {
      case 0: break;
      case 1: AddTrans (ir, values.Row(j), coefs.Col(j)); break;
      case 2:
        {
          /*
          for (size_t i = 0; i < hir.Size(); i++)
            {
              TIP<DIM,SIMD<double>> pt = hir[i].TIp<DIM>();
              SIMD<double> val1 = values(j,i);
              SIMD<double> val2 = values(j+1,i);
              __m256i mask = _mm256_set_epi64x(0, 0, -1, -1);
              double * pcoefs = &coefs(j);
              size_t dist = coefs.Dist();
              T_CalcShape (pt, 
                           SBLambda ( [&](int j, SIMD<double> shape)
                                      {
                                        auto val = HSum(shape*val1, shape*val2, shape*val2, shape*val2);
                                        val += SIMD<double,4> (_mm256_maskload_pd (pcoefs, mask));
                                        _mm256_maskstore_pd (pcoefs, mask, val.Data());
                                        pcoefs += dist;
                                      } ));
            }
          */
          /*
          SIMD<mask64,4> mask(2);
          for (size_t i = 0; i < hir.Size(); i++)
            {
              TIP<DIM,SIMD<double>> pt = hir[i].TIp<DIM>();
              SIMD<double> val1 = values(j,i);
              SIMD<double> val2 = values(j+1,i);
              double * pcoefs = &coefs(j);
              size_t dist = coefs.Dist();
              T_CalcShape (pt, 
                           SBLambda ( [val1,val2,mask,&pcoefs,dist](int j, SIMD<double> shape)
                                      {
                                        auto val = HSum(shape*val1, shape*val2, shape*val2, shape*val2);
                                        val += SIMD<double,4> (pcoefs, mask);
                                        val.Store(pcoefs, mask);
                                        pcoefs += dist;
                                      } ));
            }
          */
          for (size_t i = 0; i < hir.Size(); i++)
            {
              TIP<DIM,SIMD<double>> pt = hir[i].TIp<DIM>();
              SIMD<double> val1 = values(j,i);
              SIMD<double> val2 = values(j+1,i);
              double * pcoefs = &coefs(j);
              size_t dist = coefs.Dist();
              T_CalcShape (pt, 
                           SBLambda ( [val1,val2,&pcoefs,dist](int j, SIMD<double> shape)
                                      {
                                        auto val = HSum(shape*val1, shape*val2);
                                        val += SIMD<double,2> (pcoefs);
                                        val.Store(pcoefs);
                                        pcoefs += dist;
                                      } ));
            }
          break;
        }
      case 3:
        {
          SIMD<mask64,4> mask(3);
          for (size_t i = 0; i < hir.Size(); i++)
            {
              TIP<DIM,SIMD<double>> pt = hir[i].TIp<DIM>();
              SIMD<double> val1 = values(j,i);
              SIMD<double> val2 = values(j+1,i);
              SIMD<double> val3 = values(j+2,i);
              double * pcoefs = &coefs(j);
              size_t dist = coefs.Dist();
              T_CalcShape (pt, 
                           SBLambda ( [val1,val2,val3,mask,dist,&pcoefs](int j, SIMD<double> shape)
                                      {
                                        auto val = HSum(shape*val1, shape*val2, shape*val3, shape*val3);
                                        val += SIMD<double,4> (pcoefs, mask); 
                                        val.Store(pcoefs, mask);
                                        pcoefs += dist;
                                      } ));
            }
          break;
        }
      }
        /*
    for ( ; j < coefs.Width(); j++)
      Evaluate (ir, coefs.Col(j), values.Row(j));
        */
  }
  

  
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  auto T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGrad (const IntegrationPoint & ip, BareSliceVector<double> coefs) const -> Vec<DIM>
  {
    Vec<DIM> sum = 0.0;
    T_CalcShape (GetTIPGrad<DIM>(ip), 
                 SBLambda ( [&](int i, auto val) 
                            { 
                              sum += coefs(i) * ngfem::GetGradient(val);
                            }));
    return sum;
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGrad (const IntegrationRule & ir, BareSliceVector<double> coefs, 
                BareSliceMatrix<> vals) const
  {
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM> sum = 0.0;
        T_CalcShape (GetTIPGrad<DIM>(ir[i]), 
                     SBLambda ([&sum, coefs] (size_t j, auto shape)
                               { sum += coefs(j) * ngfem::GetGradient(shape); }));
        vals.Row(i) = sum; 
      }
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGrad (const SIMD_BaseMappedIntegrationRule & bmir,
                BareSliceVector<> coefs,
                BareSliceMatrix<SIMD<double>> values) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM, [this,&bmir,coefs,values] (auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;         
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             double *pcoefs = &coefs(0);
             const size_t dist = coefs.Dist();
             
             Vec<DIMSPACE,SIMD<double>> sum(0.0);
             this->T_CalcShape (GetTIP(mir[i]), 
                                SBLambda ([&pcoefs,dist,&sum]
                                          (size_t j, auto shape)
                                          { 
                                            sum += *pcoefs * ngfem::GetGradient(shape);
                                            pcoefs += dist;
                                          }));
             values.Col(i).Range(DIMSPACE) = sum;
           }
       });
  }

  
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGrad (const SIMD_IntegrationRule & ir,
                BareSliceVector<> coefs,
                BareSliceMatrix<SIMD<double>> values) const
  {
    for (int i = 0; i < ir.Size(); i++)
      {
        Vec<DIM,SIMD<double>> sum(0.0);
        T_CalcShape (GetTIPGrad<DIM> (ir[i]),
                     SBLambda ([&sum, coefs] (size_t j, auto shape)
                               { sum += coefs(j) * ngfem::GetGradient(shape); }));
        values.Col(i).Range(DIM) = sum;
      }
  }

  
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGradTrans (const IntegrationRule & ir, 
                     BareSliceMatrix<> vals, BareSliceVector<double> coefs) const
  {
    coefs.Range(0,ndof) = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM> vali = vals.Row(i);
        T_CalcShape (GetTIPGrad<DIM>(ir[i]), 
                     SBLambda ([coefs, vali] (int j, auto shape)
                               { coefs(j) += InnerProduct (vali, ngfem::GetGradient(shape)); }));
      }
  }
  
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGradTrans (const IntegrationRule & ir, SliceMatrix<> values, SliceMatrix<> coefs) const
  {
    int nels = coefs.Width();
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        // Vec<DIM, AutoDiff<DIM>> adp = ir[i];  
        T_CalcShape (// TIP<DIM, AutoDiff<DIM>> (adp),
                     GetTIPGrad<DIM>(ir[i]),
                     SBLambda ([&] (int j, auto shape)
                               { 
                                 FlatMatrixFixWidth<DIM> mvals(nels, &values(i,0));
                                 coefs.Row(j) += mvals * ngfem::GetGradient(shape);
                               }));
      }
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  AddGradTrans (const SIMD_BaseMappedIntegrationRule & bmir,
                BareSliceMatrix<SIMD<double>> values,
                BareSliceVector<> coefs) const
  {
    if constexpr (DIM == 0) return;
    Iterate<4-DIM>
      ([&](auto CODIM)
       {
         constexpr auto DIMSPACE = DIM+CODIM.value;
         if (bmir.DimSpace() == DIMSPACE)
           {
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
             for (size_t i = 0; i < mir.Size(); i++)
               {
                 // Directional derivative
                 [[maybe_unused]]
                   Vec<DIM, SIMD<double>> jac_dir = mir[i].GetJacobianInverse() * values.Col(i);
                 
                 const auto &ip = mir[i].IP();
                 TIP<DIM,AutoDiff<1,SIMD<double>>>adp(ip.FacetNr(), ip.VB());
                 if constexpr(DIM>0)
                     adp.x = AutoDiff<1, SIMD<double>>( ip(0), jac_dir(0) );
                 if constexpr(DIM>1)
                     adp.y = AutoDiff<1, SIMD<double>>( ip(1), jac_dir(1) );
                 if constexpr(DIM>2)
                     adp.z = AutoDiff<1, SIMD<double>>( ip(2), jac_dir(2) );

                 double * pcoef = &coefs(0);
                 size_t dist = coefs.Dist();
                 this->T_CalcShape (adp,
                                    SBLambda ([dist,&pcoef] (size_t j, auto shape)
                                              {
                                                *pcoef += HSum(shape.DValue(0));
                                                pcoef += dist;
                                              }));
               }
           }
       });
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  AddGradTrans (const SIMD_BaseMappedIntegrationRule & bmir,
                BareSliceMatrix<SIMD<double>> values,
                SliceMatrix<> coefs) const
  {
    Iterate<4-DIM>
      ([&](auto CODIM)
       {
         constexpr auto DIMSPACE = DIM+CODIM.value;
         if (bmir.DimSpace() == DIMSPACE)
           {
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);

             size_t j = 0;
             for ( ; j+4 <= coefs.Width(); j+=4)
               {
                 for (size_t i = 0; i < mir.Size(); i++)
                   {
                     TIP<DIM,AutoDiff<DIMSPACE,SIMD<double>>>adp = GetTIP(mir[i]);
                     double * pcoef = &coefs(0,j);
                     size_t dist = coefs.Dist();
                     // Vec<4*DIMSPACE,SIMD<double>> vals = values.Col(i).Range(j*DIMSPACE, (j+4)*DIMSPACE);
                     Vec<DIMSPACE,SIMD<double>> vals1 = values.Col(i).Range(j*DIMSPACE, (j+1)*DIMSPACE);
                     Vec<DIMSPACE,SIMD<double>> vals2 = values.Col(i).Range((j+1)*DIMSPACE, (j+2)*DIMSPACE);
                     Vec<DIMSPACE,SIMD<double>> vals3 = values.Col(i).Range((j+2)*DIMSPACE, (j+3)*DIMSPACE);
                     Vec<DIMSPACE,SIMD<double>> vals4 = values.Col(i).Range((j+3)*DIMSPACE, (j+4)*DIMSPACE);

                     
                     this->T_CalcShape (adp,
                                        SBLambda ([=,&pcoef] (size_t j, auto shape)
                                                  {
                                                    auto grad = ngfem::GetGradient(shape);
                                                    SIMD<double> sum1 = InnerProduct(vals1, grad);
                                                    SIMD<double> sum2 = InnerProduct(vals2, grad);
                                                    SIMD<double> sum3 = InnerProduct(vals3, grad);
                                                    SIMD<double> sum4 = InnerProduct(vals4, grad);

                                                    SIMD<double,4> allsum = HSum(sum1, sum2, sum3, sum4);
                                                    allsum += SIMD<double,4> (pcoef);
                                                    allsum.Store(pcoef);
                                                    pcoef += dist;
                                                  }));
                   }
               }

             for ( ; j+1 <= coefs.Width(); j++)
               {
                 for (size_t i = 0; i < mir.Size(); i++)
                   {
                     // TIP<DIM,AutoDiff<DIMSPACE,SIMD<double>>>adp = GetTIP(mir[i]);
                     double * pcoef = &coefs(0,j);
                     size_t dist = coefs.Dist();
                     Vec<DIMSPACE,SIMD<double>> vals = values.Col(i).Range(j*DIMSPACE, (j+1)*DIMSPACE);
                     this->T_CalcShape (GetTIP(mir[i]),   // adp
                                        SBLambda ([=,&pcoef] (size_t j, auto shape)
                                                  {
                                                    *pcoef += HSum(InnerProduct(ngfem::GetGradient(shape), vals));
                                                    pcoef += dist;
                                                  }));
                   }
               }
           }
       });
  }

  
  /*
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcDShape (const IntegrationPoint & ip, 
	      const std::function<void(int,Vec<DIM>)> & callback) const
  {
    Vec<DIM, AutoDiff<DIM> > adp;
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);
      
    // DShapeAssign<DIM> ds(dshape); 
    // T_CalcShape (&adp(0), ds);


    T_CalcShape (&adp(0), SBLambda ([&] (int i, AutoDiff<DIM> shape)
                                    {
				      Vec<DIM> v;
				      shape.StoreGradient (&v(0));
				      callback (i,v);
				    }));
  }
  */



  /*
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcMappedDShape (const MappedIntegrationPoint<DIM,DIM> & mip, 
		    FlatMatrixFixWidth<DIM> dshape) const
  {
    Vec<DIM, AutoDiff<DIM> > adp;   
    for (int i = 0; i < DIM; i++)
      adp[i].Value() = mip.IP()(i);
      
    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
	adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);

    T_CalcShape (&adp(0), SBLambda ([&] (int i, AutoDiff<DIM> shape)
                                    { shape.StoreGradient (&dshape(i,0)) ; }));
  }
  */

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcMappedDShape (const BaseMappedIntegrationPoint & bmip, 
		    BareSliceMatrix<> dshape) const
  {
    Switch<4-DIM>
      (bmip.DimSpace()-DIM, [&bmip, dshape, this](auto CODIM)
       {
         constexpr int DIM_ = DIM;
         constexpr int DIMSPACE = int(DIM)+int(CODIM.value);
         static_assert(DIM<=DIMSPACE, "dim<=dimspace");
         
         auto & mip = static_cast<const MappedIntegrationPoint<DIM_,DIMSPACE> &> (bmip);
         auto dshapes = dshape.AddSize(ndof, DIMSPACE);
         
         this->T_CalcShape (GetTIP(mip),
                            SBLambda ([dshapes] (size_t i, auto shape)
                                      { dshapes.Row(i) = ngfem::GetGradient(shape); }));
       });

    /*
    if (bmip.DimSpace() == DIM)
      {
        auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM> &> (bmip);
        auto dshapes = dshape.AddSize(ndof, DIM);
        
        T_CalcShape (GetTIP(mip),
                     SBLambda ([dshapes] (int i, auto shape)
                               { dshapes.Row(i) = ngfem::GetGradient(shape); }));
      }
    else if (bmip.DimSpace() == DIM+1)
      {
        constexpr int DIM1 = DIM<3 ? DIM+1 : DIM;
        auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM1> &> (bmip);
        auto dshapes = dshape.AddSize(ndof, DIM1);
        
        T_CalcShape (GetTIP(mip),
                     SBLambda ([dshapes] (int i, auto shape)
                               {dshapes.Row(i) = ngfem::GetGradient(shape);}));
      }
    else
      {
        cout << "CalcMappedDShape called for bboundary (not implemented)" << endl;        
      }
    */
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcMappedDShape (const BaseMappedIntegrationRule & bmir, 
		    BareSliceMatrix<> dshape) const
  {
    /*
    // auto & mir = static_cast<const MappedIntegrationRule<DIM,DIM> &> (bmir);
    for (size_t i = 0; i < bmir.Size(); i++)
      T_ScalarFiniteElement::CalcMappedDShape (bmir[i], dshape.Cols(i*DIM,(i+1)*DIM));
    */

    Switch<4-DIM>
      (bmir.DimSpace()-DIM, [&bmir, dshape, this](auto CODIM)
       {
         constexpr int DIM_ = DIM;
         constexpr int DIMSPACE = int(DIM)+int(CODIM.value);
         auto & mir = static_cast<const MappedIntegrationRule<DIM_,DIMSPACE> &> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             auto dshapes = dshape.Cols(i*DIMSPACE, (i+1)*DIMSPACE).AddSize(ndof, DIMSPACE);
             this->T_CalcShape (GetTIP(mir[i]),
                                SBLambda ([dshapes] (size_t j, auto shape)
                                          { dshapes.Row(j) = ngfem::GetGradient(shape); }));
           }
       });
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcMappedDShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                    BareSliceMatrix<SIMD<double>> dshapes) const
  {
   if (bmir.DimSpace() == DIM)
      {
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            SIMD<double> * pdshapes = dshapes.Col(i).Data();
            size_t dist = dshapes.Dist();
            
            // TIP<DIM,AutoDiff<DIM,SIMD<double>>> adp = GetTIP(mir[i]);
            T_CalcShape (GetTIP(mir[i]), // adp,
                         SBLambda ([&] (size_t j, AutoDiff<DIM,SIMD<double>> shape)
                                   { 
                                     Iterate<DIM> ( [&] (size_t ii) {
                                         *pdshapes = shape.DValue(ii);
                                         pdshapes += dist;
                                       });
                                   }));
          }
      }
   else if (bmir.DimSpace() == DIM+1)
     {
       constexpr int DIM1 = DIM<3 ? DIM+1 : DIM;
       auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
       for (size_t i = 0; i < mir.Size(); i++)
         {
           SIMD<double> * pdshapes = dshapes.Col(i).Data();
           size_t dist = dshapes.Dist();
            
           // TIP<DIM,AutoDiff<DIM1,SIMD<double>>> adp = GetTIP(mir[i]);
           T_CalcShape (GetTIP(mir[i]), // adp,
                        SBLambda ([&] (size_t j, AutoDiff<DIM1,SIMD<double>> shape)
                                  {
                                    /*
                                    Iterate<DIM1> ( [&] (size_t ii) {
                                        *pdshapes = shape.DValue(ii);
                                        pdshapes += dist;
                                      });
                                    */
                                    for (size_t k = 0; k < DIM1; k++)
                                      {
                                        *pdshapes = shape.DValue(k);
                                        pdshapes += dist;
                                      }
                                  }));
         }
     }
   else
     {
       cout << "EvaluateGrad(simd) called for bboundary (not implemented)" << endl;        
     }
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>  
  void T_ScalarFiniteElement<FEL,ET,BASE> ::
  CalcDDShape (const IntegrationPoint & ip, 
               BareSliceMatrix<> ddshape) const
  {
    TIP<DIM, AutoDiff<DIM>> t1 = ip;
    TIP<DIM, AutoDiffDiff<DIM>> tip = t1;

    T_CalcShape (tip, 
                 SBLambda ([ddshape] (size_t i, auto shape)
                           {
                             auto row = ddshape.Row(i);
                             for (int d1 = 0; d1 < DIM; d1++)
                               for (int d2 = 0; d2 < DIM; d2++)
                                 row(d1*DIM+d2) = shape.DDValue(d1,d2);
                           }));
  }


  
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> ::   
  CalcMappedDDShape (const BaseMappedIntegrationPoint & bmip, 
                     BareSliceMatrix<> ddshape) const
  {
    auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM>&> (bmip);
    T_CalcShape (GetTIPHesse (mip),
                 SBLambda ([ddshape] (size_t i, auto shape)
                           {
                             auto row = ddshape.Row(i);
                             for (int d1 = 0; d1 < DIM; d1++)
                               for (int d2 = 0; d2 < DIM; d2++)
                                 row(d1*DIM+d2) = shape.DDValue(d1,d2);
                           }));
  }


#endif

  template <class FEL, ELEMENT_TYPE ET, class BASE>  
  bool T_ScalarFiniteElement<FEL,ET,BASE> :: GetDiagDualityMassInverse (FlatVector<> diag) const 
  {
    return static_cast<const FEL*>(this)->GetDiagDualityMassInverse2(diag);
  }
  
  
  
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcDualShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const
  {
    // static_cast<const FEL*>(this) -> CalcDualShape2 (mip, shape);    
    /*
    try
      {
        static_cast<const FEL*>(this) -> CalcDualShape2 (mip, shape);
      }
    catch (const Exception& e)
      {
        double imeas = 1.0/mip.GetMeasure();
        shape = 0.0;
        static_cast<const FEL*> (this)->        
          T_CalcDualShape (GetTIP<DIM>(mip.IP()), SBLambda ( [&](int j, double val) { shape(j) = imeas * val; }));
      }
    */
    double imeas = 1.0/mip.GetMeasure();
    shape.Range(ndof) = 0.0;
    static_cast<const FEL*> (this)->        
      T_CalcDualShape (GetTIP<DIM>(mip.IP()), SBLambda ( [&](int j, double val) { shape(j) = imeas * val; }));
  }
  
  


}



#endif
