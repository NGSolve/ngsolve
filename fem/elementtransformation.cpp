/*********************************************************************/
/* File:   finiteelement.cpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


/* 
   Finite Element Definitions
*/



// #include <fem.hpp>
#include "elementtransformation.hpp"
#include "h1lofe.hpp"


namespace ngfem
{

  SIMD_BaseMappedIntegrationRule & ElementTransformation :: operator() (const SIMD_IntegrationRule & ir, Allocator & lh) const
  {
    throw ExceptionNOSIMD("ElementTransformation(SIMD_IR) not overloaded");
  }

  void ElementTransformation :: VCalcHesse (const SIMD<ngfem::IntegrationPoint> & ip, SIMD<double> * hesse) const
  {
    cout << "ElementTransformation::VCalcHesse not overloaded for " << typeid(*this).name() << endl;
  }
  
  void ElementTransformation :: CalcMultiPointJacobian (const SIMD_IntegrationRule & ir,
                                                        SIMD_BaseMappedIntegrationRule & mir) const
  {
    cout << "CalcMultiPointJacobian - SIMD not overloaded for class " << typeid(ir).name() << endl;
    throw Exception("CalcMultiPointJacobian (SIMD) not overloaded");
  }
  
  
  
  template <int DIMS, int DIMR>
  FE_ElementTransformation<DIMS, DIMR> :: FE_ElementTransformation ()
    : ElementTransformation(ET_POINT, VOL,-1,-1),  nvmat(0,0,0)
  { ; }

  template <int DIM> ScalarFiniteElement<DIM> * GetP1FE (ELEMENT_TYPE type);

  template <> ScalarFiniteElement<0> * GetP1FE (ELEMENT_TYPE type)
  {
    static FE_Point point;
    switch (type)
      {
      case ET_POINT: return (&point);
      default:
        throw ("FE_ElementTrafo, undefined 0D elementtype");
      }
  }
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
    static ScalarFE<ET_TRIG,1> trig;
    static ScalarFE<ET_QUAD,1> quad;
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
    static ScalarFE<ET_TET,1> tet;
    static FE_Hex1 hex;
    static FE_Prism1 prism;
    static ScalarFE<ET_HEXAMID,1> hexamid;
    static FE_Pyramid1 pyr;
    switch (type)
      {
      case ET_TET: return (&tet);
      case ET_HEX: return (&hex);
      case ET_PRISM: return (&prism);
      case ET_PYRAMID: return (&pyr);
      case ET_HEXAMID: return (&hexamid);        
      default:
        throw Exception("FE_ElementTrafo, undefined 3D elementtype");
      }
  }


  template <int DIMS, int DIMR>
  FE_ElementTransformation<DIMS, DIMR> :: FE_ElementTransformation (ELEMENT_TYPE type, SliceMatrix<> pmat)
    : ElementTransformation(type, VOL,-1,-1), pointmat(Trans(pmat))
  {
    fel = GetP1FE<DIMS> (type);
    // SetElementType (type);
  }

  template <int DIMS, int DIMR>
  FE_ElementTransformation<DIMS, DIMR> :: FE_ElementTransformation (ELEMENT_TYPE type)
    : FE_ElementTransformation (type, SliceMatrix<> (ElementTopology::GetNVertices(type), 
                                                     Dim(type), 3,
                                                     (double*)&ElementTopology::GetVertices(type)[0][0]))
  { ; }



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

  template <int DIMS, int DIMR>
  void FE_ElementTransformation<DIMS, DIMR> :: 
  CalcMultiPointJacobian (const SIMD_IntegrationRule & ir,
			  SIMD_BaseMappedIntegrationRule & bmir) const
  {
    SIMD_MappedIntegrationRule<DIMS,DIMR> & mir = 
      static_cast<SIMD_MappedIntegrationRule<DIMS,DIMR> &>(bmir);
    
    constexpr int SW = SIMD<IntegrationPoint>::Size();
    Vector<> shapes(ir.Size()*SW);
    MatrixFixWidth<DIMS> grad(shapes.Size());

    for (int j = 0; j < DIMR; j++)
      {
        // really slow for the moment ...
        for (int k = 0; k < ir.Size(); k++)
          {
            auto simd_ip = ir[k];
            for (int k2 = 0; k2 < SW; k2++)
              {
                shapes(k*simd_ip.Size()+k2) = fel->Evaluate(simd_ip[k2], pointmat.Row(j));
                grad.Row(k*simd_ip.Size()+k2) = fel->EvaluateGrad(simd_ip[k2], pointmat.Row(j));
              }
          }

	for (int i = 0; i < ir.Size(); i++)
	  {
	    mir[i].Point()(j) = &shapes(i*SW);
            for (int k = 0; k < DIMS; k++)
              mir[i].Jacobian()(j,k) = [&](int i2) { return grad(SW*i+i2, k); };
	  }
      }

    for (int i = 0; i < ir.Size(); i++)
      mir[i].Compute();
  }

  const ElementTransformation & GetFEElementTransformation (ELEMENT_TYPE et)
  {
    static FE_ElementTransformation<0,0> trafo_point(ET_POINT);
    static FE_ElementTransformation<1,1> trafo_segm(ET_SEGM);
    static FE_ElementTransformation<2,2> trafo_trig(ET_TRIG);
    static FE_ElementTransformation<2,2> trafo_quad(ET_QUAD);
    static FE_ElementTransformation<3,3> trafo_tet(ET_TET);
    static FE_ElementTransformation<3,3> trafo_prism(ET_PRISM);
    static FE_ElementTransformation<3,3> trafo_pyramid(ET_PYRAMID);
    static FE_ElementTransformation<3,3> trafo_hexamid(ET_HEXAMID);
    static FE_ElementTransformation<3,3> trafo_hex(ET_HEX);

    switch (et)
      {
      case ET_POINT: return trafo_point;
      case ET_SEGM: return trafo_segm;
      case ET_TRIG: return trafo_trig;
      case ET_QUAD: return trafo_quad;
      case ET_TET: return trafo_tet;
      case ET_PRISM: return trafo_prism;
      case ET_PYRAMID: return trafo_pyramid;
      case ET_HEXAMID: return trafo_hexamid;
      case ET_HEX: return trafo_hex;
      }
      
    throw Exception(string("no trafo for element type ") + ToString(et));
  }

  
  template class FE_ElementTransformation<1,1>;
  template class FE_ElementTransformation<2,2>;
  template class FE_ElementTransformation<3,3>;
  template class FE_ElementTransformation<1,2>;
  template class FE_ElementTransformation<2,3>;
  template class FE_ElementTransformation<1,3>;

  template class FE_ElementTransformation<0,0>;
  template class FE_ElementTransformation<0,1>;  
  template class FE_ElementTransformation<0,2>;
  template class FE_ElementTransformation<0,3>;
}
