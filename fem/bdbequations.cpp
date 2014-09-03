/*********************************************************************/
/* File:   bdbequations.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/
  
/*  
   Finite Element Integrators 
*/

#define FILE_BDBEQUATIONS_CPP


#include <fem.hpp>
#include <diffop_impl.hpp>



  
namespace ngfem
{ 

  /*
  template <int D, typename FEL> LaplaceBoundaryIntegrator<D,FEL> ::
  LaplaceBoundaryIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpGradientBoundary<D>, DiagDMat<D>, FEL> (DiagDMat<D> (coeff))
  { ; }
  */

  /*
  template <int D, typename FEL> RotSymLaplaceIntegrator<D,FEL> ::
  RotSymLaplaceIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpGradient<D>, RotSymLaplaceDMat<D>, FEL> (RotSymLaplaceDMat<D> (coeff))
  { ; }
  */

  // standard integratos:
  
  static RegisterBilinearFormIntegrator<LaplaceIntegrator<1> > initlap1 ("laplace", 1, 1);
  static RegisterBilinearFormIntegrator<LaplaceIntegrator<2> > initlap2 ("laplace", 2, 1);
  static RegisterBilinearFormIntegrator<LaplaceIntegrator<3> > initlap3 ("laplace", 3, 1);

  static RegisterBilinearFormIntegrator<MassIntegrator<1> > initmass1 ("mass", 1, 1);
  static RegisterBilinearFormIntegrator<MassIntegrator<2> > initmass2 ("mass", 2, 1);
  static RegisterBilinearFormIntegrator<MassIntegrator<3> > initmass3 ("mass", 3, 1);

  static RegisterBilinearFormIntegrator<RobinIntegrator<1> > initrobin1 ("robin", 1, 1);
  static RegisterBilinearFormIntegrator<RobinIntegrator<2> > initrobin2 ("robin", 2, 1);
  static RegisterBilinearFormIntegrator<RobinIntegrator<3> > initrobin3 ("robin", 3, 1);


  static RegisterBilinearFormIntegrator<LaplaceBoundaryIntegrator<2> > initlb2 ("laplaceboundary", 2, 1);
  static RegisterBilinearFormIntegrator<LaplaceBoundaryIntegrator<3> > initlb3 ("laplaceboundary", 3, 1);
  /*
      GetIntegrators().AddBFIntegrator ("laplaceboundary", 2, 1,
					LaplaceBoundaryIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("laplaceboundary", 3, 1,
					LaplaceBoundaryIntegrator<3>::Create);
  */

  static RegisterLinearFormIntegrator<SourceIntegrator<1> > initsource1 ("source", 1, 1);
  static RegisterLinearFormIntegrator<SourceIntegrator<2> > initsource2 ("source", 2, 1);
  static RegisterLinearFormIntegrator<SourceIntegrator<3> > initsource3 ("source", 3, 1);

  static RegisterLinearFormIntegrator<NeumannIntegrator<1> > initneumann1 ("neumann", 1, 1);
  static RegisterLinearFormIntegrator<NeumannIntegrator<2> > initneumann2 ("neumann", 2, 1);
  static RegisterLinearFormIntegrator<NeumannIntegrator<3> > initneumann3 ("neumann", 3, 1);
  

  static RegisterBilinearFormIntegrator<ElasticityIntegrator<2> > initelast2 ("elasticity", 2, 2);
  static RegisterBilinearFormIntegrator<ElasticityIntegrator<3> > initelast3 ("elasticity", 3, 2);

  static RegisterBilinearFormIntegrator<RotSymLaplaceIntegrator<2>> initrs2 ("rotsymlaplace", 2, 1);
  static RegisterBilinearFormIntegrator<RotSymLaplaceIntegrator<3>> initrs3 ("rotsymlaplace", 3, 1);

  static RegisterBilinearFormIntegrator<OrthoLaplaceIntegrator<2>> initolap2 ("ortholaplace", 2, 2);
  static RegisterBilinearFormIntegrator<OrthoLaplaceIntegrator<3>> initolap3 ("ortholaplace", 3, 3);




  static RegisterBilinearFormIntegrator<OrthotropicElasticityIntegrator<3>>  init_oelast3("orthoelasticity", 3, 9);
  static RegisterBilinearFormIntegrator<OrthotropicCylElasticityIntegrator<3>> init_coelast3  ("orthocylelasticity", 3, 10);

      
  static RegisterLinearFormIntegrator<GradSourceIntegrator<2>> init_gradsource2 ("gradsource", 2, 2);
  static RegisterLinearFormIntegrator<GradSourceIntegrator<3>> init_gradsource3 ("gradsource", 3, 3);

  static RegisterLinearFormIntegrator<NormalNeumannIntegrator<2>> init_normneu2("normalneumann", 2, 1);
  static RegisterLinearFormIntegrator<NormalNeumannIntegrator<3>> init_normneu3("normalneumann", 3, 1);
 
}

