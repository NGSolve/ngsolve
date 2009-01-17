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

  using namespace ngfem;


  namespace maxwellint {

    class Init
    { 
    public: 
      Init ();

    };
    
    Init::Init()
    {
      GetIntegrators().AddBFIntegrator ("curlcurledge", 3, 1,
					CurlCurlEdgeIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("curlcurledge", 2, 1,
					CurlCurlEdgeIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("curlcurlboundaryedge", 3, 1,
					CurlCurlBoundaryEdgeIntegrator::Create);
      GetIntegrators().AddBFIntegrator ("orthocurlcurledge", 3, 3,
					CurlCurlEdgeOrthoIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("massedge", 3, 1,
					MassEdgeIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("massedge", 2, 1,
					MassEdgeIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("orthomassedge", 2, 2,
					MassEdgeOrthoIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("orthomassedge", 3, 3,
					MassEdgeOrthoIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("robinedge", 3, 1,
					RobinEdgeIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("robinedge", 2, 1,
					RobinEdgeIntegrator<2>::Create);      
      GetIntegrators().AddBFIntegrator ("massedgeanisotropic", 3, 6,
					MassEdgeAnisotropicIntegrator<3>::Create);


      GetIntegrators().AddLFIntegrator ("sourceedge", 3, 3,
					SourceEdgeIntegrator<3>::Create);
      GetIntegrators().AddLFIntegrator ("sourceedge", 2, 2,
					SourceEdgeIntegrator<2>::Create);
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
