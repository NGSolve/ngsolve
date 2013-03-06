/*********************************************************************/
/* File:   l2hofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

 
#include <fem.hpp>
#include <l2hofefo.hpp>


#include "tscalarfe.cpp"



namespace ngfem
{

  template <int ORDER>
  PrecomputedShapesContainer<PrecomputedScalShapes<2> > L2HighOrderFEFO<ET_TRIG, ORDER>::precomp;

  /*
  template class L2HighOrderFiniteElementFO<1>;
  template class L2HighOrderFiniteElementFO<2>;
  template class L2HighOrderFiniteElementFO<3>;
  */

  template class T_L2HighOrderFiniteElementFO<ET_TRIG,0>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,1>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,2>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,3>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,4>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,5>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,6>;


  template class L2HighOrderFEFO<ET_TRIG,0>;
  template class L2HighOrderFEFO<ET_TRIG,1>;
  template class L2HighOrderFEFO<ET_TRIG,2>;
  template class L2HighOrderFEFO<ET_TRIG,3>;
  template class L2HighOrderFEFO<ET_TRIG,4>;
  template class L2HighOrderFEFO<ET_TRIG,5>; 
  template class L2HighOrderFEFO<ET_TRIG,6>;

  template class T_ScalarFiniteElement2<L2HighOrderFEFO<ET_TRIG,0>, ET_TRIG>;
  template class T_ScalarFiniteElement2<L2HighOrderFEFO<ET_TRIG,1>, ET_TRIG>;
  template class T_ScalarFiniteElement2<L2HighOrderFEFO<ET_TRIG,2>, ET_TRIG>;
  template class T_ScalarFiniteElement2<L2HighOrderFEFO<ET_TRIG,3>, ET_TRIG>;
  template class T_ScalarFiniteElement2<L2HighOrderFEFO<ET_TRIG,4>, ET_TRIG>;
  template class T_ScalarFiniteElement2<L2HighOrderFEFO<ET_TRIG,5>, ET_TRIG>;
  template class T_ScalarFiniteElement2<L2HighOrderFEFO<ET_TRIG,6>, ET_TRIG>;
}
 
