/*********************************************************************/
/* File:   maxwellintegrator.cpp                                     */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/


/* 
   Finite Element Integrators for Maxwell equations
*/

#define FILE_HCURL_EQUATIONS_CPP 

// #include <fem.hpp>
#include "hcurl_equations.hpp"
#include "diffop_impl.hpp"

int link_maxwellintegrator;

  
namespace ngfem
{
 
 
  ///
  template <int D, typename FEL = HCurlFiniteElement<D> >
  class MassEdgeAnisotropicIntegrator 
    : public T_BDBIntegrator<DiffOpIdEdge<D>, SymDMat<D>, FEL>
  { 
  };



  template <> 
  class MassEdgeAnisotropicIntegrator<3, HCurlFiniteElement<3> >
    : public T_BDBIntegrator<DiffOpIdEdge<3>, SymDMat<3>, HCurlFiniteElement<3> >
  {
  public:
    ///
    MassEdgeAnisotropicIntegrator (shared_ptr<CoefficientFunction> coeff00,
				   shared_ptr<CoefficientFunction> coeff10,
				   shared_ptr<CoefficientFunction> coeff11,
				   shared_ptr<CoefficientFunction> coeff20,
				   shared_ptr<CoefficientFunction> coeff21,
				   shared_ptr<CoefficientFunction> coeff22);
    /*
      : T_BDBIntegrator<DiffOpIdEdge<3>, SymDMat<3>, HCurlFiniteElement<3> >
    (SymDMat<3> (coeff00, coeff10, coeff11, coeff20, coeff21, coeff22))
    { ; }
    */

    static shared_ptr<BilinearFormIntegrator> Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return make_shared<MassEdgeAnisotropicIntegrator> (coeffs[0], coeffs[1], coeffs[2],
                                                         coeffs[3], coeffs[4], coeffs[5]);
    }
  
    ///
    virtual string Name () const 
    { return "MassEdgeAnisotropic"; }
  };


  
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

