/*********************************************************************/
/* File:   maxwellintegrator.cpp                                     */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/


/* 
   Finite Element Integrators for Maxwell equations
*/

#define FILE_HCURL_EQUATIONS_CPP 
#pragma implementation "hcurl_equations.hpp"

#include <fem.hpp>
#include <diffop_impl.hpp>

int link_maxwellintegrator;

  
namespace ngfem
{
  /*
  template <int D, typename FEL>
  CurlCurlEdgeIntegrator<D,FEL> :: CurlCurlEdgeIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpCurlEdge<D>, DiagDMat<DIM_CURL_TRAIT<D>::DIM>, FEL>
  (DiagDMat<DIM_CURL_TRAIT<D>::DIM> (coeff))
  { ; }
  
  template <int D, typename FEL>
  CurlCurlEdgeIntegrator<D,FEL> :: CurlCurlEdgeIntegrator (Array<CoefficientFunction*> & coeffs)
    : T_BDBIntegrator<DiffOpCurlEdge<D>, DiagDMat<DIM_CURL_TRAIT<D>::DIM>, FEL> (coeffs)
  { ; }


  template <int D, typename FEL>
  MassEdgeIntegrator<D,FEL> :: MassEdgeIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpIdEdge<D>, DiagDMat<D>, FEL> (DiagDMat<D> (coeff))
  { ; }
  
  template <int D, typename FEL>
  MassEdgeIntegrator<D,FEL> :: MassEdgeIntegrator (Array<CoefficientFunction*> & coeffs)
    : T_BDBIntegrator<DiffOpIdEdge<D>, DiagDMat<D>, FEL> (coeffs)
  { ; }

  CurlCurlBoundaryEdgeIntegrator ::
  CurlCurlBoundaryEdgeIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpCurlBoundaryEdge<>, DiagDMat<1>, HCurlFiniteElement<2> > 
  (DiagDMat<1> (coeff))
  { ; }
  
  CurlCurlBoundaryEdgeIntegrator :: 
  CurlCurlBoundaryEdgeIntegrator (Array<CoefficientFunction*> & coeffs)
    : T_BDBIntegrator<DiffOpCurlBoundaryEdge<>, DiagDMat<1>, HCurlFiniteElement<2> > (coeffs)
  { ; }
  */
  

  // template <> 
  MassEdgeAnisotropicIntegrator<3, HCurlFiniteElement<3> > ::
  MassEdgeAnisotropicIntegrator (shared_ptr<CoefficientFunction> coeff00,
				 shared_ptr<CoefficientFunction> coeff10,
				 shared_ptr<CoefficientFunction> coeff11,
				 shared_ptr<CoefficientFunction> coeff20,
				 shared_ptr<CoefficientFunction> coeff21,
				 shared_ptr<CoefficientFunction> coeff22)
    : T_BDBIntegrator<DiffOpIdEdge<3>, SymDMat<3>, HCurlFiniteElement<3> >
  (SymDMat<3> (coeff00, coeff10, coeff11, coeff20, coeff21, coeff22))
  { ; }
  


  /*
  template <int D, typename FEL> SourceEdgeIntegrator<D,FEL> ::
  SourceEdgeIntegrator (CoefficientFunction * coeff)
    : T_BIntegrator<DiffOpIdEdge<D>, DVec<D>, FEL> (DVec<D> (coeff))
  { ; }

  template <int D, typename FEL> SourceEdgeIntegrator<D,FEL> ::
  SourceEdgeIntegrator (CoefficientFunction * coef1,
			CoefficientFunction * coef2)
    : T_BIntegrator<DiffOpIdEdge<D>, DVec<D>, FEL> (DVec<D> (coef1, coef2))
  { ; }

  template <int D, typename FEL> SourceEdgeIntegrator<D,FEL> ::
  SourceEdgeIntegrator (CoefficientFunction * coef1,
			CoefficientFunction * coef2,
			CoefficientFunction * coef3)
    : T_BIntegrator<DiffOpIdEdge<D>, DVec<D>, FEL> (DVec<D> (coef1, coef2, coef3))
  { ; }


  template <int D, typename FEL> SourceEdgeIntegrator<D,FEL> ::
  SourceEdgeIntegrator (Array<CoefficientFunction*> & coeffs)
    : T_BIntegrator<DiffOpIdEdge<D>, DVec<D>, FEL> (coeffs)
  { ; }

  */


  /*
  template class RobinEdgeIntegrator<2>;
  template class RobinEdgeIntegrator<3>;

  template class SourceEdgeIntegrator<2>;
  template class SourceEdgeIntegrator<3>;

  template class NeumannEdgeIntegrator<2>;
  template class NeumannEdgeIntegrator<3>;
  */

  namespace maxwellint {
    
    static RegisterBilinearFormIntegrator<CurlCurlEdgeIntegrator<2> > initcce2 ("curlcurledge", 2, 1);
    static RegisterBilinearFormIntegrator<CurlCurlEdgeIntegrator<3> > initcce3 ("curlcurledge", 3, 1);

    static RegisterBilinearFormIntegrator<MassEdgeIntegrator<2> > initmasse2 ("massedge", 2, 1);
    static RegisterBilinearFormIntegrator<MassEdgeIntegrator<3> > initmasse3 ("massedge", 3, 1);

    static RegisterBilinearFormIntegrator<RobinEdgeIntegrator<2> > initrobin2 ("robinedge", 2, 1);
    static RegisterBilinearFormIntegrator<RobinEdgeIntegrator<3> > initrobin3 ("robinedge", 3, 1);
    static RegisterBilinearFormIntegrator<CurlCurlBoundaryEdgeIntegrator> initccb ("curlcurlboundaryedge", 3, 1);

    static RegisterBilinearFormIntegrator<CurlCurlEdgeOrthoIntegrator<3>> initorthocc("orthocurlcurledge", 3, 3);
    static RegisterBilinearFormIntegrator<MassEdgeOrthoIntegrator<2>> initmoe2("orthomassedge", 2, 2);
    static RegisterBilinearFormIntegrator<MassEdgeOrthoIntegrator<3>> initmoe3("orthomassedge", 3, 3);



    static RegisterLinearFormIntegrator<SourceEdgeIntegrator<2> > initse2 ("sourceedge", 2, 2);
    static RegisterLinearFormIntegrator<SourceEdgeIntegrator<3> > initse3 ("sourceedge", 3, 3);

    static RegisterLinearFormIntegrator<NeumannEdgeIntegrator<2> > initneue2 ("neumannedge", 2, 2);
    static RegisterLinearFormIntegrator<NeumannEdgeIntegrator<3> > initneue3 ("neumannedge", 3, 3);

    static RegisterLinearFormIntegrator<CurlEdgeIntegrator<2> > initcurle2 ("curledge", 2, 1);
    static RegisterLinearFormIntegrator<CurlEdgeIntegrator<3> > initcurle3 ("curledge", 3, 3);



    static RegisterLinearFormIntegrator<TangentialSourceEdgeIntegrator<2> > initts2 ("tangentialsourceedge", 2, 1);
    static RegisterLinearFormIntegrator<TangentialSourceEdgeIntegrator<3> > initts3 ("tangentialsourceedge", 3, 1);
    static RegisterLinearFormIntegrator<CurlBoundaryEdgeIntegrator<> > initcbe3 ("curlboundaryedge", 3, 1);


    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetIntegrators().AddBFIntegrator ("massedgeanisotropic", 3, 6,
					MassEdgeAnisotropicIntegrator<3>::Create);
    }

    Init init;

  }
}

