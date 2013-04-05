/*********************************************************************/
/* File:   h1hofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

 
#include <fem.hpp>

#include <h1hofe_impl.hpp>
#include <tscalarfe_impl.hpp>


namespace ngfem
{

  template <ELEMENT_TYPE ET>
  void T_H1HighOrderFiniteElement<ET> :: 
  ComputeNDof()
  {
    ndof = N_VERTEX;
    
    for (int i = 0; i < N_EDGE; i++)
      ndof += order_edge[i] - 1;
    
    for (int i = 0; i < N_FACE; i++)
      ndof += ::ngfem::PolBubbleDimension (FaceType(i), order_face[i]);

    if (DIM == 3)
      ndof += ::ngfem::PolBubbleDimension (ET, order_cell);

    order = 1;
    for (int i = 0; i < N_EDGE; i++) order = max(order, order_edge[i]);
    for (int i = 0; i < N_FACE; i++) order = max(order, Max (order_face[i])); 
    if (DIM == 3) order = max (order, Max (order_cell));
  }


  

  template class H1HighOrderFiniteElement<0>;
  template class H1HighOrderFiniteElement<1>;
  template class H1HighOrderFiniteElement<2>;
  template class H1HighOrderFiniteElement<3>;


  template NGS_DLL_HEADER class H1HighOrderFE<ET_POINT>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_SEGM>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_TRIG>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_QUAD>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_TET>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_PRISM>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_PYRAMID>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_HEX>;
 

  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_POINT>, ET_POINT>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_SEGM>, ET_SEGM>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_TRIG>, ET_TRIG>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_QUAD>, ET_QUAD>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_TET>, ET_TET>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_PRISM>, ET_PRISM>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_HEX>, ET_HEX>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_PYRAMID>, ET_PYRAMID>;


  template class T_H1HighOrderFiniteElement<ET_POINT>;
  template class T_H1HighOrderFiniteElement<ET_SEGM>;
  template class T_H1HighOrderFiniteElement<ET_TRIG>;
  template class T_H1HighOrderFiniteElement<ET_QUAD>;
  template class T_H1HighOrderFiniteElement<ET_TET>;
  template class T_H1HighOrderFiniteElement<ET_PRISM>;
  template class T_H1HighOrderFiniteElement<ET_PYRAMID>;
  template class T_H1HighOrderFiniteElement<ET_HEX>;
}
