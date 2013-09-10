#ifndef FILE_THCURLFE_IMPL
#define FILE_THCURLFE_IMPL



namespace ngfem
{

  
  /*******************************************/
  /* T_HCurlHOFiniteElement                  */
  /*******************************************/

  
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcShape (const IntegrationPoint & ip, SliceMatrix<> shape) const
  {    
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    T_CalcShape (adp, SBLambda ([&](int i, HCurl_Shape<DIM> s) 
                                { shape.Row(i) = s; }));
  }

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET, SHAPES,BASE> :: 
  CalcCurlShape (const IntegrationPoint & ip, SliceMatrix<> shape) const
  {  
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    T_CalcShape (adp, SBLambda ([&](int i, HCurl_CurlShape<DIM> s) 
                                { shape.Row(i) = s; }));
  } 

  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET, SHAPES, BASE> :: 
  CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                   SliceMatrix<> shape) const
  {
    AutoDiff<DIM> adp[DIM];

    for (int i = 0; i < DIM; i++)
      adp[i].Value() = mip.IP()(i);

    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
        adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);

    T_CalcShape (adp, SBLambda ([&](int i, HCurl_Shape<DIM> s) 
                                { shape.Row(i) = s; }));

  }

  /// compute curl of shape
  template <ELEMENT_TYPE ET, typename SHAPES, typename BASE>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES,BASE> :: 
  CalcMappedCurlShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                       SliceMatrix<> curlshape) const
  { 
    if (DIM == 2)
      {
        // not yet tested
        CalcCurlShape (mip.IP(), curlshape);
        curlshape /= mip.GetJacobiDet();        
        // MatrixFixWidth<DIM_CURL> hmat(curlshape.Height());
        // CalcCurlShape (mip.IP(), hmat);
        // curlshape = 1.0/mip.GetJacobiDet() * hmat;
      }
    else
      {
	AutoDiff<DIM> adp[DIM];
	
	for (int i = 0; i < DIM; i++)
	  adp[i].Value() = mip.IP()(i);
	
	for (int i = 0; i < DIM; i++)
	  for (int j = 0; j < DIM; j++)
	    adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);

        T_CalcShape (adp, SBLambda ([&](int i, HCurl_CurlShape<DIM> s) 
                                    { curlshape.Row(i) = s; }));
      }
  }
  
  /*
  template <ELEMENT_TYPE ET>
  Vec <DIM_CURL_TRAIT<ET_trait<ET>::DIM>::DIM>
  T_HCurlHighOrderFiniteElement<ET> :: 
  EvaluateCurlShape (const IntegrationPoint & ip, 
                     FlatVector<double> x,
                     LocalHeap & lh) const
  {
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    HCurlEvaluateCurl<DIM> ds(x); 
    static_cast<const HCurlHighOrderFE<ET>*> (this) -> T_CalcShape (adp, ds);
    return ds.Sum();
  }
  */

}


#endif
