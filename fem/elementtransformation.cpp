/*********************************************************************/
/* File:   finiteelement.cpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


/* 
   Finite Element Definitions
*/



#include <fem.hpp>


namespace ngfem
{
  template <int DIMS, int DIMR>
  FE_ElementTransformation<DIMS, DIMR> :: FE_ElementTransformation ()
    : pointmat(0,0,0), pointmat_ownmem(false), nvmat(0,0,0)
  { ; }
  

  template <int DIMS, int DIMR>
  FE_ElementTransformation<DIMS, DIMR> :: ~FE_ElementTransformation ()
  {
    if (pointmat_ownmem) delete [] &pointmat(0,0); 
  }
  
  
  ///
  template <int DIMS, int DIMR>
  void FE_ElementTransformation<DIMS, DIMR> :: 
  CalcJacobian (const IntegrationPoint & ip, FlatMatrix<> dxdxi) const
  {
    for (int i = 0; i < DIMR; i++)
      dxdxi.Row(i) = fel->EvaluateGrad (ip, pointmat.Row(i));
  }
  

  template <int DIMS, int DIMR>
  void FE_ElementTransformation<DIMS, DIMR> :: 
  CalcPoint (const IntegrationPoint & ip, 
	     FlatVector<> point) const
  {
    for (int i = 0; i < DIMR; i++)
      point(i) = fel->Evaluate (ip, pointmat.Row(i));
  }
  
  template <int DIMS, int DIMR>
  void FE_ElementTransformation<DIMS, DIMR> :: 
  CalcPointJacobian (const IntegrationPoint & ip,
		     FlatVector<> point, 
		     FlatMatrix<> dxdxi) const
  {
    CalcPoint (ip, point);
    CalcJacobian (ip, dxdxi);
  }


  template <int DIMS, int DIMR>
  void FE_ElementTransformation<DIMS, DIMR> :: 
  CalcMultiPointJacobian (const IntegrationRule & ir,
			  BaseMappedIntegrationRule & bmir) const
  {
    MappedIntegrationRule<DIMS,DIMR> & mir = 
      static_cast<MappedIntegrationRule<DIMS,DIMR> &>(bmir);
    
    Vector<> shapes(ir.Size());
    MatrixFixWidth<DIMS> grad(ir.Size());

    for (int j = 0; j < DIMR; j++)
      {
	fel->Evaluate (ir, pointmat.Row(j), shapes);
	fel->EvaluateGrad (ir, pointmat.Row(j), grad);
	
	for (int i = 0; i < ir.Size(); i++)
	  {
	    mir[i].Point()(j) = shapes(i);
	    mir[i].Jacobian().Row(j) = grad.Row(i);
	  }
      }

    for (int i = 0; i < ir.Size(); i++)
      mir[i].Compute();
  }
  
  
  template class FE_ElementTransformation<1,1>;
  template class FE_ElementTransformation<2,2>;
  template class FE_ElementTransformation<3,3>;
  template class FE_ElementTransformation<1,2>;
  template class FE_ElementTransformation<2,3>;
}
