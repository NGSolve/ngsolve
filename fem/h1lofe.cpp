/*********************************************************************/
/* File:   h1lofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

 
#include <fem.hpp>

namespace ngfem
{

  /*
  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  void T_ScalarFiniteElement2<FEL,ET,NDOF,ORDER> ::
  CalcShape (const IntegrationPoint & ip, 
	     FlatVector<> shape) const
  {
    double pt[DIM];
      for (int i = 0; i < DIM; i++) pt[i] = ip(i);
      Spec().T_CalcShape (pt, shape); 
    }

  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  void T_ScalarFiniteElement2<FEL,ET,NDOF,ORDER> ::
  CalcShapeStat (const IntegrationPoint & ip, 
                               FlatVector<> shape) 
    {
      double pt[DIM];
      for (int i = 0; i < DIM; i++) pt[i] = ip(i);
      FEL::T_CalcShape (pt, shape); 
    }
    
    
  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  void T_ScalarFiniteElement2<FEL,ET,NDOF,ORDER> ::
  CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<DIM> dshape) const
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);
      
      DShapeAssign<DIM> ds(dshape); 
      static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);
      
      // FEL::CalcDShapeStat (ip, dshape);
    }

  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  void T_ScalarFiniteElement2<FEL,ET,NDOF,ORDER> ::
  CalcDShapeStat (const IntegrationPoint & ip, 
				FlatMatrixFixWidth<DIM> dshape) 
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);
      
      DShapeAssign<DIM> ds(dshape); 
      FEL::T_CalcShape (adp, ds);
      
      // FEL::CalcDShapeStat (ip, dshape);
    }



  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  void T_ScalarFiniteElement2<FEL,ET,NDOF,ORDER> ::
  CalcMappedDShape (const SpecificIntegrationPoint<DIM,DIM> & sip, 
		    FlatMatrixFixWidth<DIM> dshape) const
  {
      AutoDiff<DIM> adp[DIM];
      
      for (int i = 0; i < DIM; i++)
        adp[i].Value() = sip.IP()(i);
      
      for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
          adp[i].DValue(j) = sip.GetJacobianInverse()(i,j);
      
      DShapeAssign<DIM> ds(dshape); 
      Spec().T_CalcShape (adp, ds);
    }
  */

  /*
  template class T_ScalarFiniteElement2<FE_Segm0,ET_SEGM,1,0>;
  template class T_ScalarFiniteElement2<FE_Segm1,ET_SEGM,2,1>;
  template class T_ScalarFiniteElement2<FE_Segm1,ET_SEGM,2,1>;
  template class T_ScalarFiniteElement2<FE_Segm1L2,ET_SEGM,2,1>;
  template class T_ScalarFiniteElement2<FE_Segm2,ET_SEGM,3,2>;
  template class T_ScalarFiniteElement2<FE_Segm2HB,ET_SEGM,3,2>;
  template class T_ScalarFiniteElement2<FE_Segm2L2,ET_SEGM,3,2>;
  template class T_ScalarFiniteElement2<FE_Trig0,ET_TRIG,1,0>;
  template class T_ScalarFiniteElement2<FE_Trig1,ET_TRIG,3,1>;
  template class T_ScalarFiniteElement2<FE_Trig2,ET_TRIG,6,2>;

  template class T_ScalarFiniteElement2<FE_Quad0,ET_QUAD,1,0>;
  template class T_ScalarFiniteElement2<FE_Quad1,ET_QUAD,4,1>;

  template class T_ScalarFiniteElement2<FE_Trig1,ET_TRIG,3,1>;
  */  
}
