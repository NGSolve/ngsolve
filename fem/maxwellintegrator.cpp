/*********************************************************************/
/* File:   maxwellintegrator.cpp                                     */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/


/* 
   Finite Element Integrators for Maxwell equations
*/
 

#include <fem.hpp>
int link_maxwellintegrator;

  
namespace ngfem
{
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




  template class MassEdgeIntegrator<2>;
  template class MassEdgeIntegrator<3>;

  template class CurlCurlEdgeIntegrator<2>;
  template class CurlCurlEdgeIntegrator<3>;

  template class SourceEdgeIntegrator<2>;
  template class SourceEdgeIntegrator<3>;


  namespace maxwellint {
    
    static RegisterBilinearFormIntegrator<CurlCurlEdgeIntegrator<2> > initcce2 ("curlcurledge", 2, 1);
    static RegisterBilinearFormIntegrator<CurlCurlEdgeIntegrator<3> > initcce3 ("curlcurledge", 3, 1);

    static RegisterBilinearFormIntegrator<MassEdgeIntegrator<2> > initmasse2 ("massedge", 2, 1);
    static RegisterBilinearFormIntegrator<MassEdgeIntegrator<3> > initmasse3 ("massedge", 3, 1);

    static RegisterBilinearFormIntegrator<RobinEdgeIntegrator<2> > initrobin2 ("robinedge", 2, 1);
    static RegisterBilinearFormIntegrator<RobinEdgeIntegrator<3> > initrobin3 ("robinedge", 3, 1);


    static RegisterLinearFormIntegrator<SourceEdgeIntegrator<2> > initse2 ("sourceedge", 2, 2);
    static RegisterLinearFormIntegrator<SourceEdgeIntegrator<3> > initse3 ("sourceedge", 3, 3);


    class Init
    { 
    public: 
      Init ();

    };
    
    Init::Init()
    {
      GetIntegrators().AddBFIntegrator ("curlcurlboundaryedge", 3, 1,
					CurlCurlBoundaryEdgeIntegrator::Create);
      GetIntegrators().AddBFIntegrator ("orthocurlcurledge", 3, 3,
					CurlCurlEdgeOrthoIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("orthomassedge", 2, 2,
					MassEdgeOrthoIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("orthomassedge", 3, 3,
					MassEdgeOrthoIntegrator<3>::Create);
      /*
      GetIntegrators().AddBFIntegrator ("robinedge", 3, 1,
					RobinEdgeIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("robinedge", 2, 1,
					RobinEdgeIntegrator<2>::Create);      
      */
      GetIntegrators().AddBFIntegrator ("massedgeanisotropic", 3, 6,
					MassEdgeAnisotropicIntegrator<3>::Create);

      GetIntegrators().AddLFIntegrator ("neumannedge", 3, 3,
					NeumannEdgeIntegrator<3>::Create);
      GetIntegrators().AddLFIntegrator ("neumannedge", 2, 2,
					NeumannEdgeIntegrator<2>::Create);
      GetIntegrators().AddLFIntegrator ("curlboundaryedge", 3, 1,
					CurlBoundaryEdgeIntegrator<>::Create);
      GetIntegrators().AddLFIntegrator ("curledge", 2, 1,
					CurlEdgeIntegrator<2>::Create); 
      GetIntegrators().AddLFIntegrator ("curledge", 3, 3,
					CurlEdgeIntegrator<3>::Create);
      GetIntegrators().AddLFIntegrator ("tangentialsourceedge", 3, 1,
					TangentialSourceEdgeIntegrator<3>::Create);
      GetIntegrators().AddLFIntegrator ("tangentialsourceedge", 2, 1,
					TangentialSourceEdgeIntegrator<2>::Create);
    }

    Init init;
  }
}
