/*********************************************************************/
/* File:   hdiv_equations.cpp                                        */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/

/*
   Finite Element Integrators
*/

#define FILE_HDIV_EQUATIONS_CPP 


#include <fem.hpp>
#include "hdiv_equations.hpp"
#include <diffop_impl.hpp>

namespace ngfem
{
  static RegisterBilinearFormIntegrator<MassHDivIntegrator<2>> init_mhd2("masshdiv", 2, 1);
  static RegisterBilinearFormIntegrator<MassHDivIntegrator<3>> init_mhd3("masshdiv", 3, 1);
  static RegisterBilinearFormIntegrator<DivDivHDivIntegrator<2>> init_ddhd2("divdivhdiv", 2, 1);
  static RegisterBilinearFormIntegrator<DivDivHDivIntegrator<3>> init_ddhd3("divdivhdiv", 3, 1);
  static RegisterBilinearFormIntegrator<RobinHDivIntegrator<2>> init_rhd2("robinhdiv", 2, 1);
  static RegisterBilinearFormIntegrator<RobinHDivIntegrator<3>> init_rhd3("robinhdiv", 3, 1);
  
  static RegisterLinearFormIntegrator<DivSourceHDivIntegrator<2>> init_dshd2("divsource", 2, 1);
  static RegisterLinearFormIntegrator<DivSourceHDivIntegrator<3>> init_dshd3("divsource", 3, 1);
  static RegisterLinearFormIntegrator<SourceHDivIntegrator<2>> init_shd2("sourcehdiv", 2, 2);
  static RegisterLinearFormIntegrator<SourceHDivIntegrator<3>> init_shd3("sourcehdiv", 3, 3);
  static RegisterLinearFormIntegrator<NeumannHDivIntegrator<2>> init_nhd2("neumannhdiv", 2, 1);
  static RegisterLinearFormIntegrator<NeumannHDivIntegrator<3>> init_nhd3("neumannhdiv", 3, 1);
}



