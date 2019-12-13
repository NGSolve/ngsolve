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
#include <diffop_impl.hpp>

namespace ngfem
{

 
  using namespace ngfem;

  /*
  template <int D> MassHDivIntegrator<D> ::
  MassHDivIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpIdHDiv<D>, DiagDMat<D>, HDivFiniteElement<D> > (DiagDMat<D> (coeff))
  { ; }

  template <int D> DivDivHDivIntegrator<D> ::
  DivDivHDivIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpDivHDiv<D>, DiagDMat<1>, HDivFiniteElement<D> > (DiagDMat<1> (coeff))
  { ; }
  
  template <int D> DivSourceHDivIntegrator<D> ::
  DivSourceHDivIntegrator (CoefficientFunction * coeff)
    : T_BIntegrator<DiffOpDivHDiv<D>, DVec<1>, HDivFiniteElement<D> > (DVec<1> (coeff))
  { ; }


  template <int D> 
  Integrator * BaseSourceHDivIntegrator<D> ::
  Create (Array<CoefficientFunction*> & coeffs)
  {
    if ((coeffs.Size()==1) && (coeffs[0]->Dimension() == D)){
      if(D == 2)
	return new SourceHDivIntegratorN<2> (coeffs[0]);
      else // if (D == 3)
	return new SourceHDivIntegratorN<3> (coeffs[0]);
    }
    if(D == 2)
      return new SourceHDivIntegrator<2> (coeffs[0], coeffs[1]);
    else // if (D == 3)
      return new SourceHDivIntegrator<3> (coeffs[0], coeffs[1], coeffs[2]);
  }

  
  SourceHDivIntegrator<3> ::
  SourceHDivIntegrator (CoefficientFunction * coeff1,
                        CoefficientFunction * coeff2,
                        CoefficientFunction * coeff3)
    : BaseSourceHDivIntegrator<3> (DVec<3> (coeff1, coeff2,coeff3))
  { ; }
  
  SourceHDivIntegrator<2> ::
  SourceHDivIntegrator (CoefficientFunction * coeff1,
                        CoefficientFunction * coeff2)
    : BaseSourceHDivIntegrator<2> (DVec<2> (coeff1, coeff2))
  { ; }
  */


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



