/*********************************************************************/
/* File:   hcurlhofe.cpp                                             */
/* Author: Sabine Zaglmayr                                           */
/* Date:   20. Maerz 2003                                            */
/*                                                                   */
/* AutoDiff - revision: J. Schoeberl, March 2009                     */
/*********************************************************************/

#include <fem.hpp>    
#include <hcurlhofe_impl.hpp>


namespace ngfem
{

  /*******************************************/
  /* T_HCurlHOFiniteElement                  */
  /*******************************************/

  template <ELEMENT_TYPE ET, typename SHAPES>
  T_HCurlHighOrderFiniteElement<ET, SHAPES> :: 
  T_HCurlHighOrderFiniteElement (int aorder)
  {
    for (int i = 0; i < N_EDGE; i++) order_edge[i] = aorder;
    for (int i = 0; i < N_FACE; i++) order_face[i] = aorder;
    if (DIM == 3) order_cell = aorder;

    for(int i = 0; i < N_EDGE; i++) usegrad_edge[i] = 1;
    for(int i=0; i < N_FACE; i++) usegrad_face[i] = 1;
    if (DIM == 3) usegrad_cell = 1;

    for (int i = 0; i < N_VERTEX; i++) vnums[i] = i;
    // eltype = ET;
  }

  
  template <ELEMENT_TYPE ET, typename SHAPES>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES> :: 
  CalcShape (const IntegrationPoint & ip, FlatMatrixFixWidth<DIM> shape) const
  {    
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    /*
    HCurlShapeAssign<DIM> ds(shape); 
    static_cast<const SHAPES*> (this) -> T_CalcShape (adp, ds);
    */
    T_CalcShape (adp, SBLambda ([&](int i, HCurl_Shape<DIM> s) 
                                { shape.Row(i) = s; }));
  }

  template <ELEMENT_TYPE ET, typename SHAPES>
  void T_HCurlHighOrderFiniteElement<ET, SHAPES> :: 
  CalcCurlShape (const IntegrationPoint & ip, FlatMatrixFixWidth<DIM_CURL> shape) const
  {  
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);
    
    /*
    HCurlCurlShapeAssign<DIM> ds(shape); 
    static_cast<const SHAPES*> (this) -> T_CalcShape (adp, ds);
    */
    T_CalcShape (adp, SBLambda ([&](int i, HCurl_CurlShape<DIM> s) 
                                { shape.Row(i) = s; }));
  } 

  template <ELEMENT_TYPE ET, typename SHAPES>
  void T_HCurlHighOrderFiniteElement<ET, SHAPES> :: 
  CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                   FlatMatrixFixWidth<DIM> shape) const
  {
    AutoDiff<DIM> adp[DIM];

    for (int i = 0; i < DIM; i++)
      adp[i].Value() = mip.IP()(i);

    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
        adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);

    /*
    HCurlShapeAssign<DIM> ds(shape); 
    static_cast<const SHAPES*> (this) -> T_CalcShape (adp, ds);
    */
    T_CalcShape (adp, SBLambda ([&](int i, HCurl_Shape<DIM> s) 
                                { shape.Row(i) = s; }));

  }

  /// compute curl of shape
  template <ELEMENT_TYPE ET, typename SHAPES>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES> :: 
  CalcMappedCurlShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                       FlatMatrixFixWidth<DIM_CURL> curlshape) const
  { 
    if (DIM == 2)
      {
        // not yet tested
        CalcCurlShape (mip.IP(), curlshape);
        curlshape /= mip.GetJacobiDet();        
      }
    else
      {
	AutoDiff<DIM> adp[DIM];
	
	for (int i = 0; i < DIM; i++)
	  adp[i].Value() = mip.IP()(i);
	
	for (int i = 0; i < DIM; i++)
	  for (int j = 0; j < DIM; j++)
	    adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);

        /*
	HCurlCurlShapeAssign<DIM> ds(curlshape); 
	static_cast<const SHAPES*> (this) -> T_CalcShape (adp, ds);
        */
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

  template <ELEMENT_TYPE ET, typename SHAPES>
  void T_HCurlHighOrderFiniteElement<ET,SHAPES> :: 
  ComputeNDof()
  {
    ndof = N_EDGE;

    for (int i = 0; i < N_EDGE; i++)
      if(order_edge[i] > 0)
        ndof += usegrad_edge[i]*order_edge[i];
    
    for(int i = 0; i < N_FACE; i++)
      if (FaceType(i) == ET_TRIG)
        {
          if (order_face[i][0] > 1) 
            ndof += ((usegrad_face[i]+1)*order_face[i][0]+2)*(order_face[i][0]-1)/2 ;
        }
      else
        {
          if(order_face[i][0]>=0 && order_face[i][1]>=0)
            ndof +=  (usegrad_face[i]+1)*order_face[i][0]*order_face[i][1] 
              + order_face[i][0] + order_face[i][1]; 
        }

    switch (ET)
      {
      case ET_TET: 
        if(order_cell[0] > 2)
          ndof += ((usegrad_cell + 2) * order_cell[0] + 3) 
            * (order_cell[0]-2) * (order_cell[0]-1) / 6; 
        break;
      case ET_PRISM:
        if(order_cell[2] > 0 && order_cell[0] > 1)
          ndof += ((usegrad_cell+2)*order_cell[2] + 1) * order_cell[0]*(order_cell[0]-1)/2
            + (order_cell[0]-1)*order_cell[2]; 
        break;
      case ET_PYRAMID:
        {
          int pc = order_cell[0]; //SZ: no problem to do anisotropic, but for the moment 
          // is it worth getting crazy :-) 
          if(order_cell[0]>1)
            ndof += usegrad_cell*(pc-1)*pc*(2*pc-1)/6 + pc*(2*pc*pc+3*pc-2)/3; 
          break;
        }
      case ET_HEX:
        if(order_cell[0] >= 0 && order_cell[1]>= 0 && order_cell[2]>=0)
          ndof += (usegrad_cell + 2)* order_cell[0] * order_cell[1] * order_cell[2]
            + order_cell[1]*order_cell[2]  + order_cell[0]*(order_cell[1] + order_cell[2]);  
        break;
      default:
        ;
      }

    order = 0; // max(order_edges,order_face,order_cell);  
    for (int i = 0; i < N_EDGE; i++)
      order = max (order, order_edge[i]);

    for(int i=0; i < N_FACE; i++) 
      if (ET_trait<ET>::FaceType(i) == ET_TRIG)
        order = max (order, order_face[i][0]);
      else
        order = max (order, Max (order_face[i]));

    if (DIM == 3)
      order = max (order, Max(order_cell));

    // for integration order .. 
    if (ET == ET_PRISM || ET == ET_HEX || ET == ET_PYRAMID || ET == ET_QUAD)
      order++;
    else
      if (order==0) order++;
  }


  template class  HCurlHighOrderFiniteElement<1>;
  template class  HCurlHighOrderFiniteElement<2>;
  template class  HCurlHighOrderFiniteElement<3>; 

  template class HCurlHighOrderFE<ET_SEGM>;
  template class HCurlHighOrderFE<ET_TRIG>;
  template class HCurlHighOrderFE<ET_QUAD>;
  template class HCurlHighOrderFE<ET_TET>;
  template class HCurlHighOrderFE<ET_HEX>;
  template class HCurlHighOrderFE<ET_PRISM>;
  template class HCurlHighOrderFE<ET_PYRAMID>;
}
