/*********************************************************************/
/* File:   hdiv_equations.cpp                                        */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/

/*
   Finite Element Integrators
*/



#include <fem.hpp>

namespace ngfem
{

 
  using namespace ngfem;

  
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


  template class MassHDivIntegrator<2>;
  template class DivDivHDivIntegrator<2>;
  // template class SourceHDivIntegrator<2>;
  template class SourceHDivIntegratorN<2>;
  template class DivSourceHDivIntegrator<2>;

  template class MassHDivIntegrator<3>;
  template class DivDivHDivIntegrator<3>;
  // template class SourceHDivIntegrator<3>;
  template class SourceHDivIntegratorN<3>;
  template class DivSourceHDivIntegrator<3>;




  namespace hdiv_equations_cpp
  {
    class Init
    {
    public:
      Init ();
    };


    
    Init::Init()
    {
      GetIntegrators().AddBFIntegrator ("masshdiv", 2, 1,
					MassHDivIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("masshdiv", 3, 1,
					MassHDivIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("divdivhdiv", 2, 1,
					DivDivHDivIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("divdivhdiv", 3, 1,
					DivDivHDivIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("robinhdiv", 2, 1,
					RobinHDivIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("robinhdiv", 3, 1,
					RobinHDivIntegrator<3>::Create);

      GetIntegrators().AddLFIntegrator ("divsource", 2, 1,
					DivSourceHDivIntegrator<2>::Create);
      GetIntegrators().AddLFIntegrator ("divsource", 3, 1,
					DivSourceHDivIntegrator<3>::Create);


      GetIntegrators().AddLFIntegrator ("sourcehdiv", 2, 2,
					SourceHDivIntegrator<2>::Create);
      GetIntegrators().AddLFIntegrator ("sourcehdiv", 3, 3,
					SourceHDivIntegrator<3>::Create);
      GetIntegrators().AddLFIntegrator ("neumannhdiv", 2, 1,
					NeumannHDivIntegrator<2>::Create);
      GetIntegrators().AddLFIntegrator ("neumannhdiv", 3, 1,
					NeumannHDivIntegrator<3>::Create);
    }

    Init init;
  }
}



