/*********************************************************************/
/* File:   bdbequations.cpp                                          */
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
  


  template <int D>  MassIntegrator<D> :: 
  MassIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpId<D>, DiagDMat<1>, NodalFiniteElement<D> > (DiagDMat<1> (coeff))
  { ; }

  template <int D, typename FEL> LaplaceIntegrator<D,FEL> ::
  LaplaceIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpGradient<D>, DiagDMat<D>, FEL> (DiagDMat<D> (coeff))
  { ; }

  template <int D, typename FEL> LaplaceBoundaryIntegrator<D,FEL> ::
  LaplaceBoundaryIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpGradientBoundary<D>, DiagDMat<D>, FEL> (DiagDMat<D> (coeff))
  { ; }


  template <int D, typename FEL> RotSymLaplaceIntegrator<D,FEL> ::
  RotSymLaplaceIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpGradient<D>, RotSymLaplaceDMat<D>, FEL> (RotSymLaplaceDMat<D> (coeff))
  { ; }



  template <int D> RobinIntegrator<D> ::
  RobinIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpIdBoundary<D>, DiagDMat<1>, NodalFiniteElement<D-1> > (DiagDMat<1> (coeff))
  { ; }





  template <int D, typename FEL> SourceIntegrator<D,FEL> ::
  SourceIntegrator (CoefficientFunction * coeff)
    : T_BIntegrator<DiffOpId<D>, DVec<1>, FEL> (DVec<1> (coeff))
  { ; }

  template <int D, typename FEL> NeumannIntegrator<D,FEL> ::
  NeumannIntegrator (CoefficientFunction * coeff)
    : T_BIntegrator<DiffOpIdBoundary<D>, DVec<1>, FEL> (DVec<1> (coeff))
  { ; }
  



  template class MassIntegrator<1>;
  template class MassIntegrator<2>;
  template class MassIntegrator<3>;

  template class LaplaceIntegrator<1>;
  template class LaplaceIntegrator<2>;
  template class LaplaceIntegrator<3>;

  template class LaplaceBoundaryIntegrator<2>;
  template class LaplaceBoundaryIntegrator<3>;


  template class RobinIntegrator<2>;
  template class RobinIntegrator<3>;


  template class SourceIntegrator<1>;
  template class SourceIntegrator<2>;
  template class SourceIntegrator<3>;

  template class NeumannIntegrator<2>;
  template class NeumannIntegrator<3>;






  // standard integratos:
  
  namespace bdbequations_cpp
  {
    class Init
    { 
    public:  
      Init ();
    };        
    
    Init::Init()
    {
      GetIntegrators().AddBFIntegrator ("laplace", 2, 1,
					LaplaceIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("laplace", 3, 1,
					LaplaceIntegrator<3>::Create);
      /*
      GetIntegrators().AddBFIntegrator ("rotsymlaplace", 2, 1,
					RotSymLaplaceIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("ortholaplace", 2, 2,
					OrthoLaplaceIntegrator<2>::Create);
      */
      GetIntegrators().AddBFIntegrator ("ortholaplace", 2, 2,
					OrthoLaplaceIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("ortholaplace", 3, 3,
					OrthoLaplaceIntegrator<3>::Create);
      

      GetIntegrators().AddBFIntegrator ("mass", 2, 1,
					MassIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("mass", 3, 1,
					MassIntegrator<3>::Create);
      /*
      GetIntegrators().AddBFIntegrator ("divdiv", 2, 1,
					DivDivIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("curlcurl", 2, 1,
					CurlCurlIntegrator<>::Create);
      GetIntegrators().AddBFIntegrator ("curlcurl", 3, 1,
					CurlCurl3dIntegrator<>::Create);
      */

      GetIntegrators().AddBFIntegrator ("robin", 2, 1,
					RobinIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("robin", 3, 1,
					RobinIntegrator<3>::Create);

      GetIntegrators().AddBFIntegrator ("laplaceboundary", 2, 1,
					LaplaceBoundaryIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("laplaceboundary", 3, 1,
					LaplaceBoundaryIntegrator<3>::Create);
      /*
      GetIntegrators().AddBFIntegrator ("normalrobin", 2, 1,
					NormalRobinIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("normalrobin", 3, 1,
					NormalRobinIntegrator<3>::Create);
      */      
      GetIntegrators().AddBFIntegrator ("elasticity", 2, 2,
					ElasticityIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("elasticity", 3, 2,
					ElasticityIntegrator<3>::Create);
      
      GetIntegrators().AddBFIntegrator ("orthoelasticity", 2, 9,
					OrthotropicElasticityIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("orthoelasticity", 3, 9,
					OrthotropicElasticityIntegrator<3>::Create);
      
      GetIntegrators().AddBFIntegrator ("orthocylelasticity", 2, 10,
					OrthotropicCylElasticityIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("orthocylelasticity", 3, 10,
					OrthotropicCylElasticityIntegrator<3>::Create);
      
      /*
      GetIntegrators().AddBFIntegrator ("dirichletpenalty", 3, 1,
					DirichletPenaltyIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("dirichletpenalty", 2, 1,
					DirichletPenaltyIntegrator<2>::Create);
      */


      GetIntegrators().AddLFIntegrator ("source", 2, 1,
					SourceIntegrator<2>::Create);
      GetIntegrators().AddLFIntegrator ("source", 3, 1,
					SourceIntegrator<3>::Create);
      GetIntegrators().AddLFIntegrator("gradsource", 3, 3, 
				       GradSourceIntegrator<3>::Create); 

      GetIntegrators().AddLFIntegrator ("neumann", 2, 1,
					NeumannIntegrator<2>::Create);
      GetIntegrators().AddLFIntegrator ("neumann", 3, 1,
					NeumannIntegrator<3>::Create);

      GetIntegrators().AddLFIntegrator ("normalneumann", 2, 1,
					NormalNeumannIntegrator<2>::Create);
      GetIntegrators().AddLFIntegrator ("normalneumann", 3, 1,
					NormalNeumannIntegrator<3>::Create);
    }
    
    Init init;
    
  }
 
}
