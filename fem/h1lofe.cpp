/*********************************************************************/
/* File:   h1lofe.cpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

 
#include <fem.hpp>
#include <tscalarfe_impl.hpp>


namespace ngfem
{

  template class T_ScalarFiniteElementFO<ScalarDummyFE<ET_POINT>,ET_POINT,0,0>;
  template class T_ScalarFiniteElementFO<ScalarDummyFE<ET_SEGM>,ET_SEGM,0,0>;
  template class T_ScalarFiniteElementFO<ScalarDummyFE<ET_TRIG>,ET_TRIG,0,0>;
  template class T_ScalarFiniteElementFO<ScalarDummyFE<ET_QUAD>,ET_QUAD,0,0>;
  template class T_ScalarFiniteElementFO<ScalarDummyFE<ET_TET>,ET_TET,0,0>;
  template class T_ScalarFiniteElementFO<ScalarDummyFE<ET_PRISM>,ET_PRISM,0,0>;
  template class T_ScalarFiniteElementFO<ScalarDummyFE<ET_PYRAMID>,ET_PYRAMID,0,0>;
  template class T_ScalarFiniteElementFO<ScalarDummyFE<ET_HEX>,ET_HEX,0,0>;


  
  template class T_ScalarFiniteElementFO<FE_Point,ET_POINT,1,0>;
  template class T_ScalarFiniteElementFO<FE_Segm0,ET_SEGM,1,0>;
  template class T_ScalarFiniteElementFO<FE_Segm1,ET_SEGM,2,1>;
  template class T_ScalarFiniteElementFO<FE_Segm2,ET_SEGM,3,2>;

  template class T_ScalarFiniteElementFO<FE_Segm2HB,ET_SEGM,3,2>;
  template class T_ScalarFiniteElementFO<FE_Segm1L2,ET_SEGM,2,1>;

  template class T_ScalarFiniteElementFO<FE_Segm2L2,ET_SEGM,3,2>;

  template class T_ScalarFiniteElementFO<FE_NcSegm1,ET_SEGM,1,1>;
  template class T_ScalarFiniteElementFO<FE_Segm3Pot,ET_SEGM,4,3>;


  template class T_ScalarFiniteElementFO<FE_TSegmL2<0>,ET_SEGM,1,0>;
  template class T_ScalarFiniteElementFO<FE_TSegmL2<1>,ET_SEGM,2,1>;
  template class T_ScalarFiniteElementFO<FE_TSegmL2<2>,ET_SEGM,3,2>;
  template class T_ScalarFiniteElementFO<FE_TSegmL2<3>,ET_SEGM,4,3>;




  template class T_ScalarFiniteElement<FE_Trig0,ET_TRIG>;
  template class T_ScalarFiniteElement<FE_Trig1,ET_TRIG>;
  template class T_ScalarFiniteElement<FE_Trig2,ET_TRIG>;
  template class T_ScalarFiniteElement<FE_Trig2HB,ET_TRIG>;
  template class T_ScalarFiniteElementFO<FE_NcTrig1,ET_TRIG,3,1>;
  template class  T_ScalarFiniteElementFO<FE_Quad0,ET_QUAD,1,0>;
  template class  T_ScalarFiniteElementFO<FE_Quad1,ET_QUAD,4,1>;
  template class  T_ScalarFiniteElementFO<FE_Quad2,ET_QUAD,9,2>;
  template class  T_ScalarFiniteElementFO<FE_Quad2aniso,ET_QUAD,6,2>;

  template class  T_ScalarFiniteElementFO<FE_Tet0,ET_TET,1,0>;
  template class  T_ScalarFiniteElementFO<FE_Tet1,ET_TET,4,1>;
  template class  T_ScalarFiniteElementFO<FE_Tet2,ET_TET,10,2>;
  template class  T_ScalarFiniteElementFO<FE_Tet2HB,ET_TET,10,2>;
  template class  T_ScalarFiniteElementFO<FE_NcTet1,ET_TET,4,1>;


  template class  T_ScalarFiniteElementFO<FE_Prism0,ET_PRISM,1,0>;
  template class  T_ScalarFiniteElementFO<FE_Prism1,ET_PRISM,6,1>;
  template class  T_ScalarFiniteElementFO<FE_Prism2,ET_PRISM,18,2>;
  template class  T_ScalarFiniteElementFO<FE_Prism2aniso,ET_PRISM,12,2>;
  template class  T_ScalarFiniteElementFO<FE_Prism2HBaniso,ET_PRISM,12,2>;


  template class  T_ScalarFiniteElementFO<FE_Hex0,ET_HEX,1,0>;
  template class  T_ScalarFiniteElementFO<FE_Hex1,ET_HEX,8,1>;

  template class  T_ScalarFiniteElementFO<FE_Pyramid0,ET_PYRAMID,1,0>;
  template class  T_ScalarFiniteElementFO<FE_Pyramid1,ET_PYRAMID,5,1>;
}
