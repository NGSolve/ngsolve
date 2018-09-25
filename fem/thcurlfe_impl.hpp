#ifndef FILE_THCURLFE_IMPL
#define FILE_THCURLFE_IMPL



namespace ngfem
{

  /*
  template <int DIM, typename SCAL = double>
  class HCurl_Shape : public Vec<DIM,SCAL>
  {
  public:
    template <typename T>
    HCurl_Shape (T shape) : Vec<DIM,SCAL>(shape.Value()) { ; }
  };

  template <int DIM, typename SCAL = double>
  class HCurl_CurlShape : public Vec<DIM_CURL_TRAIT<DIM>::DIM,SCAL>
  {
  public:
    template <typename T>
    HCurl_CurlShape (T shape) 
      : Vec<DIM_CURL_TRAIT<DIM>::DIM,SCAL> (shape.CurlValue()) { ; }
  };
  */


  
  /*******************************************/
  /* T_HCurlHOFiniteElement                  */
  /*******************************************/

  
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcShape (const IntegrationPoint & ip, SliceMatrix<> shape) const
  {    
    Vec<DIM, AutoDiff<DIM> > adp = ip; 
    T_CalcShape (&adp(0), SBLambda ([shape](size_t i, auto s) 
                                    { FlatVec<DIM> (&shape(i,0)) = s.Value(); }));
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET, SHAPES,BASE> :: 
  CalcCurlShape (const IntegrationPoint & ip, SliceMatrix<> shape) const
  {  
    Vec<DIM, AutoDiff<DIM> > adp = ip; 
    T_CalcShape (&adp(0), SBLambda ([shape](size_t i, auto s) 
                                    { FlatVec<DIM_CURL_(DIM)> (&shape(i,0)) = s.CurlValue(); }));
  } 

#ifndef FASTCOMPILE
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET, SHAPES, BASE> :: 
  CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
                   SliceMatrix<> shape) const
  {
    auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM>&> (bmip);
    Vec<DIM, AutoDiff<DIM> > adp = mip; 
    T_CalcShape (&adp(0), SBLambda ([shape](size_t i, auto s) 
				    { 
				      FlatVec<DIM> (&shape(i,0)) = s.Value(); 
				    }));
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcMappedShape (const MappedIntegrationRule<DIM,DIM> & mir, 
                   SliceMatrix<> shape) const
  {
    for (size_t i = 0; i < mir.Size(); i++)
      CalcMappedShape (mir[i], shape.Cols(i*DIM,(i+1)*DIM));
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
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
                 Vec<DIM, AutoDiff<DIMSPACE,SIMD<double>>> adp = mir[i];
                 auto shapei = shapes.Col(i);
                 T_CalcShape (&adp(0), SBLambda ([shapei,DIMSPACE] (size_t j, auto s)
                                                 {
                                                   auto shape = s.Value();
                                                   for (size_t k = 0; k < DIMSPACE; k++)
                                                     shapei(j*DIMSPACE+k) = shape(k);
                                                 }));
               }
             
           }
       });
    /*
    if ((DIM == 3) || (bmir.DimSpace() == DIM))
      {
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_Shape<DIM,SIMD<double>> shape)
                                            {
                                              for (size_t k = 0; k < DIM; k++)
                                                shapes(j*DIM+k, i) = shape(k);
                                            }));
          }
      }
    else
      {
        constexpr int DIM1 = DIM==3 ? DIM : DIM+1;
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_Shape<DIM1,SIMD<double>> shape)
                                            {
                                              for (size_t k = 0; k < DIM1; k++)
                                                shapes(j*DIM1+k,i) = shape(k);
                                            }));
          }
      }
    */
  }


  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcMappedCurlShape (const BaseMappedIntegrationPoint & bmip,
                       SliceMatrix<> curlshape) const
  {
    auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM>&> (bmip);    
    if (DIM == 2)
      {
        CalcCurlShape (mip.IP(), curlshape);
        curlshape /= mip.GetJacobiDet();        
      }
    else
      {
        Vec<DIM, AutoDiff<DIM> > adp = mip; 
        T_CalcShape (&adp(0), SBLambda ([&](size_t i, auto s) 
                                        { 
                                          FlatVec<DIM_CURL_(DIM)> (&curlshape(i,0)) = s.CurlValue(); 
                                        }));
      }
  }


  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcMappedCurlShape (const MappedIntegrationRule<DIM,DIM> & mir, 
                       SliceMatrix<> curlshape) const
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
    Iterate<4-DIM>
      ([this,&bmir,shapes](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         if (bmir.DimSpace() == DIMSPACE)
           {
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
             constexpr int DIM_CURL = DIM_CURL_(DIMSPACE);        
             for (size_t i = 0; i < mir.Size(); i++)
               {
                 Vec<DIM, AutoDiff<DIMSPACE,SIMD<double>>> adp = mir[i];
                 auto shapei = shapes.Col(i);
                 T_CalcShape (&adp(0), SBLambda ([shapei,DIM_CURL] (size_t j, auto s)
                                                 {
                                                   auto cshape = s.CurlValue();
                                                   for (size_t k = 0; k < DIM_CURL; k++)
                                                     shapei(j*DIM_CURL+k) = cshape(k);
                                                 }));
               }
           }
       });

    /*
    if ((DIM == 3) || (bmir.DimSpace() == DIM))
      {
        constexpr int DIM_CURL = DIM_CURL_(DIM);
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_CurlShape<DIM,SIMD<double>> cshape)
                                            {
                                              for (size_t k = 0; k < DIM_CURL; k++)
                                                shapes(j*DIM_CURL+k, i) = cshape(k);
                                            }));
          }
      }
    else
      {
        constexpr int DIM1 = DIM==3 ? DIM : DIM+1;
        constexpr int DIM_CURL = DIM_CURL_(DIM1);        
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_CurlShape<DIM1,SIMD<double>> cshape)
                                            {
                                              for (size_t k = 0; k < DIM_CURL; k++)
                                                shapes(j*DIM_CURL+k,i) = cshape(k);
                                            }));
          }
      }
    */
  }

  

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  auto T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  EvaluateCurlShape (const IntegrationPoint & ip, 
                     FlatVector<double> x,
                     LocalHeap & lh) const -> Vec<DIM_CURL_(DIM)>
  {
    Vec<DIM, AutoDiff<DIM> > adp = ip; 
    Vec<DIM_CURL_(DIM)> sum = 0.0;
    T_CalcShape (&adp(0), SBLambda ([&sum, x](size_t i, auto s) 
                                    { sum += x(i) * s.CurlValue(); }));
    return sum;
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  EvaluateCurl (const IntegrationRule & ir, FlatVector<> coefs, FlatMatrixFixWidth<DIM_CURL_(DIM)> curl) const
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
    Iterate<4-DIM>
      ([this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         if (bmir.DimSpace() == DIMSPACE)
           {
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
             for (size_t i = 0; i < mir.Size(); i++)
               {
                 Vec<DIM, AutoDiff<DIMSPACE,SIMD<double>>> adp = mir[i];
                 Vec<DIMSPACE,SIMD<double>> sum(0.0);
                 T_CalcShape (&adp(0), SBLambda ([&sum,coefs] (size_t j, auto shape)
                                                 {
                                                   sum += coefs(j) * shape.Value();
                                                 }));
                 for (size_t k = 0; k < DIMSPACE; k++)
                   values(k,i) = sum(k); 
               }
           }
       });

    /*
    if ((DIM == 3) || (bmir.DimSpace() == DIM))
      {
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
            Vec<DIM,SIMD<double>> sum(0.0);
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_Shape<DIM,SIMD<double>> shape)
                                            {
                                            // double coef = coefs(j);
                                            //  Iterate<DIM> ( [&] (auto ii) {
                                            //    sum(ii.value) += coef * shape(ii.value);
                                            //  });
                                              sum += coefs(j) * Vec<DIM,SIMD<double>> (shape);
                                            }));
            for (size_t k = 0; k < DIM; k++)
              values(k,i) = sum(k).Data();
          }
      }
    else if (bmir.DimSpace() == DIM+1)
      {
        constexpr int DIM1 = DIM<3 ? DIM+1 : DIM;
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            Vec<DIM1,SIMD<double>> sum(0.0);            
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_Shape<DIM1,SIMD<double>> shape)
                                            {
                                              double coef = coefs(j);
                                              for (int k = 0; k < DIM1; k++)
                                                sum(k) += coef * shape(k);
                                            }));
            for (size_t k = 0; k < DIM1; k++)
              values(k,i) = sum(k).Data();
          }
      }
    else
      {
        constexpr int DIM1 = DIM<2 ? DIM+2 : DIM;
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            Vec<DIM1,SIMD<double>> sum(0.0);            
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_Shape<DIM1,SIMD<double>> shape)
                                            {
                                              double coef = coefs(j);
                                              for (int k = 0; k < DIM1; k++)
                                                sum(k) += coef * shape(k);
                                            }));
            for (size_t k = 0; k < DIM1; k++)
              values(k,i) = sum(k).Data();
          }

      }
    */
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  Evaluate (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<Complex> coefs,
            BareSliceMatrix<SIMD<Complex>> values) const
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
                 Vec<DIM, AutoDiff<DIMSPACE,SIMD<double>>> adp = mir[i];
                 Vec<DIMSPACE,SIMD<Complex>> sum = SIMD<Complex>(0.0);
                 T_CalcShape (&adp(0), SBLambda ([&sum,coefs] (size_t j, auto shape)
                                                 {
                                                   sum += coefs(j) * shape.Value();
                                                 }));
                 for (size_t k = 0; k < DIMSPACE; k++)
                   values(k,i) = sum(k);
               }
           }
       });
    /*
    if ((DIM == 3) || (bmir.DimSpace() == DIM))
      {
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
            Vec<DIM,SIMD<Complex>> sum = SIMD<Complex>(0.0);
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_Shape<DIM,SIMD<double>> shape)
                                            {
                                              Iterate<DIM> ( [&] (auto ii) {
                                                  sum(ii.value) += shape(ii.value) * coefs(j);
                                                });
                                            }));
            for (size_t k = 0; k < DIM; k++)
              values(k,i) = sum(k);
          }
      }
    else if (bmir.DimSpace() == DIM+1)
      {
        constexpr int DIM1 = DIM<3 ? DIM+1 : DIM;
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            Vec<DIM1,SIMD<Complex>> sum = SIMD<Complex>(0.0);
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_Shape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<Complex> coef = coefs(j);
                                              for (int k = 0; k < DIM1; k++)
                                                sum(k) += shape(k) * coef;
                                            }));
            for (size_t k = 0; k < DIM1; k++)
              values(k,i) = sum(k);
          }
      }
    else
      {
        constexpr int DIM1 = DIM<2 ? DIM+2 : DIM;
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            Vec<DIM1,SIMD<Complex>> sum = SIMD<Complex>(0.0);
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_Shape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<Complex> coef = coefs(j);
                                              for (int k = 0; k < DIM1; k++)
                                                sum(k) += shape(k) * coef;
                                            }));
            for (size_t k = 0; k < DIM1; k++)
              values(k,i) = sum(k);
          }
      }
    */
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  EvaluateCurl (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
  {
    Iterate<4-DIM>
      ([this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         if (bmir.DimSpace() == DIMSPACE)
           {
             constexpr int DIM_CURL = DIM_CURL_(DIMSPACE);
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
             for (size_t i = 0; i < mir.Size(); i++)
               {
                 Vec<DIM, AutoDiff<DIMSPACE,SIMD<double>>> adp = mir[i];
                 Vec<DIM_CURL,SIMD<double>> sum(0.0);            
                 T_CalcShape (&adp(0), SBLambda ([coefs,&sum] (size_t j, auto shape)
                                                 {
                                                   sum += coefs(j) * shape.CurlValue();
                                                 }));
                 for (size_t k = 0; k < DIM_CURL; k++)
                   values(k,i) = sum(k).Data();
               }
           }
       });
    
    /*
    if ((DIM == 3) || (bmir.DimSpace() == DIM))
      {
        constexpr int DIM_CURL = DIM_CURL_(DIM);                
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
            Vec<DIM_CURL,SIMD<double>> sum(0.0);
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_CurlShape<DIM,SIMD<double>> shape)
                                            {
                                              double coef = coefs(j);
                                              for (int k = 0; k < DIM_CURL; k++)
                                                sum(k) += shape(k) * coef;
                                            }));
            for (size_t k = 0; k < DIM_CURL; k++)
              values(k,i) = sum(k).Data();
          }
      }
    else if (bmir.DimSpace() == DIM+1)
      {
        constexpr int DIM1 = DIM<3 ? DIM+1 : DIM;
        constexpr int DIM_CURL = DIM_CURL_(DIM1);                        
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            Vec<DIM_CURL,SIMD<double>> sum(0.0);            
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_CurlShape<DIM1,SIMD<double>> shape)
                                            {
                                              double coef = coefs(j);
                                              for (int k = 0; k < DIM_CURL; k++)
                                                sum(k) += coef * shape(k);
                                            }));
            for (size_t k = 0; k < DIM_CURL; k++)
              values(k,i) = sum(k).Data();
          }
      }
    else
      {
        constexpr int DIM1 = DIM<2 ? DIM+2 : DIM;
        constexpr int DIM_CURL = DIM_CURL_(DIM1);                        
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            Vec<DIM_CURL,SIMD<double>> sum(0.0);            
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_CurlShape<DIM1,SIMD<double>> shape)
                                            {
                                              double coef = coefs(j);
                                              for (int k = 0; k < DIM_CURL; k++)
                                                sum(k) += coef * shape(k);
                                            }));
            for (size_t k = 0; k < DIM_CURL; k++)
              values(k,i) = sum(k).Data();
          }

      }
    */
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  EvaluateCurl (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const
  {
    Iterate<4-DIM>
      ([this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         if (bmir.DimSpace() == DIMSPACE)
           {
             constexpr int DIM_CURL = DIM_CURL_(DIMSPACE);
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
             
             for (size_t i = 0; i < mir.Size(); i++)
               {
                 Vec<DIM, AutoDiff<DIMSPACE,SIMD<double>>> adp = mir[i];
                 Vec<DIM_CURL,SIMD<Complex>> sum = SIMD<Complex>(0.0);
                 T_CalcShape (&adp(0), SBLambda ([coefs, &sum] (size_t j, auto shape)
                                                 {
                                                   sum += coefs(j) * shape.CurlValue();
                                                 }));
                 for (size_t k = 0; k < DIM_CURL; k++)
                   values(k,i) = sum(k);
               }
           }
       });
  
    /*
    if ((DIM == 3) || (bmir.DimSpace() == DIM))
      {
        constexpr int DIM_CURL = DIM_CURL_(DIM);                
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
            Vec<DIM_CURL,SIMD<Complex>> sum = SIMD<Complex>(0.0);
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_CurlShape<DIM,SIMD<double>> shape)
                                            {
					      Iterate<DIM> ([&] (auto ii) {
						  sum(ii.value) += shape(ii.value) * coefs(j);
						});
                                            }));
            for (size_t k = 0; k < DIM_CURL; k++)
              values(k,i) = sum(k);
          }
      }
    else if (bmir.DimSpace() == DIM+1)
      {
        constexpr int DIM1 = DIM<3 ? DIM+1 : DIM;
        constexpr int DIM_CURL = DIM_CURL_(DIM1);                        
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            Vec<DIM_CURL,SIMD<Complex>> sum = SIMD<Complex>(0.0);            
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_CurlShape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<Complex> coef = coefs(j);
                                              for (int k = 0; k < DIM_CURL; k++)
                                                sum(k) += coef * shape(k);
                                            }));
            for (size_t k = 0; k < DIM_CURL; k++)
              values(k,i) = sum(k);
          }
      }
    else
      {
        constexpr int DIM1 = DIM<2 ? DIM+2 : DIM;
        constexpr int DIM_CURL = DIM_CURL_(DIM1);                        
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            Vec<DIM_CURL,SIMD<Complex>> sum = SIMD<Complex>(0.0);            
            T_CalcShape (&adp(0), SBLambda ([&] (int j, HCurl_CurlShape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<Complex> coef = coefs(j);
                                              for (int k = 0; k < DIM_CURL; k++)
                                                sum(k) += coef * shape(k);
                                            }));
            for (size_t k = 0; k < DIM_CURL; k++)
              values(k,i) = sum(k);
          }
	
      }
    */
  }

  
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  AddTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
            BareSliceVector<> coefs) const
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
                 Vec<DIM, AutoDiff<DIMSPACE,SIMD<double>>> adp = mir[i];
                 Vec<DIMSPACE,SIMD<double>> vali = values.Col(i);
                 T_CalcShape (&adp(0), SBLambda ([vali,coefs] (size_t j, auto s)
                                                 {
                                                   /*
                                                   auto shape = s.Value();
                                                   SIMD<double> sum = 0.0;
                                                   for (size_t k = 0; k < shape.Size(); k++)
                                                     sum += shape(k) * vali(k);
                                                   coefs(j) += HSum(sum);
                                                   */
                                                   coefs(j) += HSum(InnerProduct(s.Value(), vali));
                                                 }));
               }
           }
       });

    /*
    if ((DIM == 3) || (bmir.DimSpace() == DIM))
      {
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_Shape<DIM,SIMD<double>> shape)
                                            {
                                              SIMD<double> sum = 0.0;
                                              for (int k = 0; k < DIM; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }
      }
    else if(bmir.DimSpace() == DIM+1)
      {
        constexpr int DIM1 = DIM<3 ? DIM+1 : DIM;
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_Shape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<double> sum = 0.0;
                                              for (int k = 0; k < DIM1; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }
      }
    else
      {
	constexpr int DIM1 = DIM<2 ? DIM+2 : DIM;
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_Shape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<double> sum = 0.0;
                                              for (int k = 0; k < DIM1; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }
      }
    */
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  AddTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<Complex>> values,
            BareSliceVector<Complex> coefs) const
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
                 Vec<DIM, AutoDiff<DIMSPACE,SIMD<double>>> adp = mir[i];
                 Vec<DIMSPACE,SIMD<Complex>> vali = values.Col(i);
                 T_CalcShape (&adp(0), SBLambda ([vali,coefs] (size_t j, auto s)
                                                 {
                                                   /*
                                                   auto shape = s.Value();
                                                   SIMD<Complex> sum = 0.0;
                                                   for (size_t k = 0; k < shape.Size(); k++)
                                                     sum += shape(k) * vali(k);
                                                   coefs(j) += HSum(sum);
                                                   */
                                                   coefs(j) += HSum(InnerProduct(s.Value(), vali));
                                                 }));
               }
           }
       });

    /*
    if ((DIM == 3) || (bmir.DimSpace() == DIM))
      {
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_Shape<DIM,SIMD<double>> shape)
                                            {
                                              SIMD<Complex> sum = 0.0;
                                              for (int k = 0; k < DIM; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }
      }
    else if(bmir.DimSpace() == DIM+1)
      {
        constexpr int DIM1 = DIM<3 ? DIM+1 : DIM;
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_Shape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<Complex> sum = 0.0;
                                              for (int k = 0; k < DIM1; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }
      }
    else
      {
        constexpr int DIM1 = DIM<2 ? DIM+2 : DIM;
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_Shape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<Complex> sum = 0.0;
                                              for (int k = 0; k < DIM1; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }

      }
    */
  }



  
  
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  AddCurlTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
                BareSliceVector<> coefs) const
  {
    Iterate<4-DIM>
      ([this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         constexpr int DIM_CURL = DIM_CURL_(DIMSPACE);                                         
         if (bmir.DimSpace() == DIMSPACE)
           {
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
             for (size_t i = 0; i < mir.Size(); i++)
               {
                 Vec<DIM, AutoDiff<DIMSPACE,SIMD<double>>> adp = mir[i];
                 Vec<DIM_CURL,SIMD<double>> vali = values.Col(i);
                 T_CalcShape (&adp(0), SBLambda ([vali,coefs] (size_t j, auto s)
                                                 {
                                                   /*
                                                   auto cshape = s.CurlValue();
                                                   SIMD<double> sum = 0.0;
                                                   for (size_t k = 0; k < cshape.Size(); k++)
                                                     sum += cshape(k) * vali(k);
                                                   coefs(j) += HSum(sum);
                                                   */
                                                   coefs(j) += HSum(InnerProduct(s.CurlValue(), vali));
                                                 }));
               }
           }
       });


    
    /*
    // throw ExceptionNOSIMD ("thcurlfe - simd - add curl trans not implemented");        
    if ((DIM == 3) || (bmir.DimSpace() == DIM))
      {
        constexpr int DIM_CURL = DIM_CURL_(DIM);                        
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_CurlShape<DIM,SIMD<double>> shape)
                                            {
                                              SIMD<double> sum = 0.0;
                                              for (int k = 0; k < DIM_CURL; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }
      }
    else if (bmir.DimSpace() == DIM+1)
      {
        constexpr int DIM1 = DIM<3 ? DIM+1 : DIM;
        constexpr int DIM_CURL = DIM_CURL_(DIM1);                                
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_CurlShape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<double> sum = 0.0;
                                              for (int k = 0; k < DIM_CURL; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }
      }
    else
      {
        constexpr int DIM1 = DIM<2 ? DIM+2 : DIM;
        constexpr int DIM_CURL = DIM_CURL_(DIM1);                                
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_CurlShape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<double> sum = 0.0;
                                              for (int k = 0; k < DIM_CURL; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }

      }
    */
  }
  
  
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  AddCurlTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<Complex>> values,
                BareSliceVector<Complex> coefs) const
  {
    Iterate<4-DIM>
      ([this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         constexpr int DIM_CURL = DIM_CURL_(DIMSPACE);                                         
         if (bmir.DimSpace() == DIMSPACE)
           {
             auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
             for (size_t i = 0; i < mir.Size(); i++)
               {
                 Vec<DIM, AutoDiff<DIMSPACE,SIMD<double>>> adp = mir[i];
                 Vec<DIM_CURL,SIMD<Complex>> vali = values.Col(i);
                 T_CalcShape (&adp(0), SBLambda ([vali,coefs] (size_t j, auto s)
                                                 {
                                                   /*
                                                   auto cshape = s.CurlValue();
                                                   SIMD<Complex> sum = 0.0;
                                                   for (size_t k = 0; k < cshape.Size(); k++)
                                                     sum += cshape(k) * vali(k);
                                                   coefs(j) += HSum(sum);
                                                   */
                                                   coefs(j) += HSum(InnerProduct(s.CurlValue(), vali));
                                                 }));
               }
           }
       });

    /*
    if ((DIM == 3) || (bmir.DimSpace() == DIM))
      {
        constexpr int DIM_CURL = DIM_CURL_(DIM);                        
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_CurlShape<DIM,SIMD<double>> shape)
                                            {
                                              SIMD<Complex> sum = 0.0;
                                              for (int k = 0; k < DIM_CURL; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }
      }
    else if (bmir.DimSpace() == DIM+1)
      {
        constexpr int DIM1 = DIM<3 ? DIM+1 : DIM;
        constexpr int DIM_CURL = DIM_CURL_(DIM1);                                
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_CurlShape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<Complex> sum = 0.0;
                                              for (int k = 0; k < DIM_CURL; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }
      }
    else
      {
        constexpr int DIM1 = DIM<2 ? DIM+2 : DIM;
        constexpr int DIM_CURL = DIM_CURL_(DIM1);                                
        auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM1>&> (bmir);
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Vec<DIM, AutoDiff<DIM1,SIMD<double>>> adp = mir[i];
            T_CalcShape (&adp(0), SBLambda ([&] (size_t j, HCurl_CurlShape<DIM1,SIMD<double>> shape)
                                            {
                                              SIMD<Complex> sum = 0.0;
                                              for (int k = 0; k < DIM_CURL; k++)
                                                sum += shape(k) * values(k,i);
                                              coefs(j) += HSum(sum);
                                            }));
          }

      }
    */
  }

#endif
}


#endif
