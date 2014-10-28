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
    : /* pointmat(0,0,0), pointmat_ownmem(false), */ nvmat(0,0,0)
  { ; }
  

  template <int DIM> ScalarFiniteElement<DIM> * GetP1FE (ELEMENT_TYPE type);

  template <> ScalarFiniteElement<1> * GetP1FE (ELEMENT_TYPE type)
  {
    static FE_Segm1 segm;
    switch (type)
      {
      case ET_SEGM: return (&segm);
      default:
        throw ("FE_ElementTrafo, undefined 1D elementtype");
      }
  }
  template <> ScalarFiniteElement<2> * GetP1FE (ELEMENT_TYPE type)
  {
    static FE_Trig1 trig;
    static FE_Quad1 quad;
    switch (type)
      {
      case ET_TRIG: return (&trig);
      case ET_QUAD: return (&quad);
      default:
        throw ("FE_ElementTrafo, undefined 2D elementtype");
      }
  }
  template <> ScalarFiniteElement<3> * GetP1FE (ELEMENT_TYPE type)
  {
    static FE_Tet1 tet;
    static FE_Hex1 hex;
    static FE_Prism1 prism;
    static FE_Pyramid1 pyr;
    switch (type)
      {
      case ET_TET: return (&tet);
      case ET_HEX: return (&hex);
      case ET_PRISM: return (&prism);
      case ET_PYRAMID: return (&pyr);
      default:
        throw ("FE_ElementTrafo, undefined 3D elementtype");
      }
  }


  template <int DIMS, int DIMR>
  FE_ElementTransformation<DIMS, DIMR> :: FE_ElementTransformation (ELEMENT_TYPE type, FlatMatrix<> pmat)
    : pointmat(Trans(pmat))
  {
    fel = GetP1FE<DIMS> (type);
  }

  template <int DIMS, int DIMR>
  FE_ElementTransformation<DIMS, DIMR> :: ~FE_ElementTransformation ()
  {
    // if (pointmat_ownmem) delete [] &pointmat(0,0); 
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
