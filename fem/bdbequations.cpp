/*********************************************************************/
/* File:   bdbequations.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/
  
/*  
   Finite Element Integrators 
*/

#define FILE_BDBEQUATIONS_CPP


// #include <fem.hpp>
#include "finiteelement.hpp"
#include "bdbequations.hpp"
#include "elasticity_equations.hpp"
#include "diffop_impl.hpp"

  
namespace ngfem
{ 





  template <int D, typename FEL>
  void DiffOpHesseBoundary<D,FEL> ::
  GenerateMatrixSIMDIR (const FiniteElement & bfel,
                        const SIMD_BaseMappedIntegrationRule & bmir,
                        BareSliceMatrix<SIMD<double>> mat)
  {
    auto & fel = static_cast<const FEL&>(bfel);
    size_t nd = fel.GetNDof();

    STACK_ARRAY(SIMD<double>, mem, nd*D*D);
    FlatMatrix<SIMD<double>> ddshape(nd, D*D, mem);

    for (size_t i = 0; i < bmir.Size(); i++)
      {
        fel.CalcMappedDDShape(bmir[i], ddshape);
        mat.Col(i) = ddshape.AsVector();
      }

    /*
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (bmir);
      
    STACK_ARRAY(SIMD<double>, mem1, 6*D*nd_u);
    FlatMatrix<SIMD<double>> shape_u(nd_u*D, 4, &mem1[0]);

    auto shape_ul = shape_u.Col(0);
    auto shape_ur = shape_u.Col(1);
    auto shape_ull = shape_u.Col(2);
    auto shape_urr = shape_u.Col(3);
    
    // FlatMatrix<SIMD<double>> dshape_u_ref(nd_u*(D-1), 1, &mem1[4*D*nd_u]);
    FlatMatrix<SIMD<double>> dshape_u(nd_u*D, 1, &mem1[5*D*nd_u]);

    LocalHeapMem<10000> lh("diffophesse-lh");

    auto & ir = mir.IR();
    for (size_t i = 0; i < mir.Size(); i++)
      {
        const SIMD<IntegrationPoint> & ip = ir[i];
        const ElementTransformation & eltrans = mir[i].GetTransformation();

        for (int j = 0; j < D-1; j++)   // d / dxj
          {
            HeapReset hr(lh);
            SIMD<IntegrationPoint> ipts[4];
            ipts[0] = ip;
            ipts[0](j) -= eps();
            ipts[1] = ip;
            ipts[1](j) += eps();              
            ipts[2] = ip;
            ipts[2](j) -= 2*eps();
            ipts[3] = ip;
            ipts[3](j) += 2*eps();

            SIMD_IntegrationRule ir(4, ipts);
            SIMD_MappedIntegrationRule<D-1,D> mirl(ir, eltrans, lh);

            fel.CalcMappedDShape (mirl, shape_u);
            
            dshape_u.Col(0) = (1.0/(12.0*eps())) * (8.0*shape_ur-8.0*shape_ul-shape_urr+shape_ull);
            for (size_t l = 0; l < D; l++)
              for (size_t k = 0; k < nd_u; k++)
                mat(k*D*D+j*D+l, i) = dshape_u(k*D+l, 0);
          }
          
          for (size_t j = 0; j < D; j++)
            for (size_t k = 0; k < nd_u; k++)
              {
                Vec<D-1,SIMD<double>> dshape_u_ref;
                Vec<D,SIMD<double>> dshape_u;
                
                for (size_t l = 0; l < D-1; l++)
                  dshape_u_ref(l) = mat(k*D*D+l*D+j, i);
                
                dshape_u = Trans(mir[i].GetJacobianInverse()) * dshape_u_ref;
                
                for (size_t l = 0; l < D; l++)
                  mat(k*D*D+l*D+j, i) = dshape_u(l);
              }
        }
    */
  }

  
  template <int D, typename FEL>
  void DiffOpHesseBoundary<D,FEL> ::
  ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
               BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
  {
    auto & fel_u = static_cast<const FEL&>(fel);
    size_t nd = fel_u.GetNDof();

    STACK_ARRAY(SIMD<double>, mem, nd*D*D);
    FlatMatrix<SIMD<double>> ddshape(nd, D*D, mem);

    y.Rows(D*D).Cols(bmir.Size()) = 0.0;
    for (size_t i = 0; i < bmir.Size(); i++)
      {
        fel_u.CalcMappedDDShape(bmir[i], ddshape);
        for (size_t j = 0; j < nd; j++)
          y.Col(i) += x(j)*ddshape.Row(j);
      }
    /*
      int size = (bmir.Size()+1)*500*SIMD<double>::Size();
      STACK_ARRAY(char, data, size);
      LocalHeap lh(data, size);

      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (bmir);
      auto & ir = mir.IR();
      const ElementTransformation & trafo = mir.GetTransformation();
      auto & fel_u = static_cast<const FEL&>(fel);
      FlatMatrix<SIMD<double>> hxl(D, mir.IR().Size(), lh);
      FlatMatrix<SIMD<double>> hxr(D, mir.IR().Size(), lh);
      FlatMatrix<SIMD<double>> hxll(D, mir.IR().Size(), lh);
      FlatMatrix<SIMD<double>> hxrr(D, mir.IR().Size(), lh);
      FlatMatrix<SIMD<double>> hx(D, mir.IR().Size(), lh);

      
      for (int k = 0; k < mir.Size(); k++)
        for (int m = 0; m < D*D; m++)
          y(m, k) = SIMD<double> (0.0);
      
      for (int j = 0; j < D-1; j++)
        {
          // hx = (F^-1 * x).Row(j)
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irl(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irl.Size(); k++)
              {
                irl[k] = ir[k];
                irl[k](j) -= eps();
              }
            SIMD_MappedIntegrationRule<D-1,D> mirl(irl, trafo, lh);
            fel_u.EvaluateGrad (mirl, x, hxl);
          }
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irr(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irr.Size(); k++)
              {
                irr[k] = ir[k];              
                irr[k](j) += eps();
              }
            SIMD_MappedIntegrationRule<D-1,D> mirr(irr, trafo, lh);
            fel_u.EvaluateGrad (mirr, x, hxr);
          }
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irll(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irll.Size(); k++)
              {
                irll[k] = ir[k];
                irll[k](j) -= 2*eps();
              }
            SIMD_MappedIntegrationRule<D-1,D> mirll(irll, trafo, lh);
            fel_u.EvaluateGrad (mirll, x, hxll);
          }
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irrr(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irrr.Size(); k++)
              {
                irrr[k] = ir[k];              
                irrr[k](j) += 2*eps();
              }
            SIMD_MappedIntegrationRule<D-1,D> mirrr(irrr, trafo, lh);
            fel_u.EvaluateGrad (mirrr, x, hxrr);
          }
          // hx = 1.0/(2*eps()) * (hxr-hxl);
          // dshape_u_ref = (1.0/(12.0*eps)) * (8.0*shape_ur-8.0*shape_ul-shape_urr+shape_ull);
          hx = 1.0/(12*eps()) * (8*hxr-8*hxl-hxrr+hxll);
          for (int k = 0; k < mir.Size(); k++)
            {
              auto jacinv = mir[k].GetJacobianInverse();
              for (int l = 0; l < D; l++)
                {
                  for (int m = 0; m < D; m++)
                    y(m*D+l, k) += jacinv(j,m) * hx(l, k);
                }
            }
        }
    */
  }


  template <int D, typename FEL>
  void DiffOpHesseBoundary<D,FEL> ::
  AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                  BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
  {
    auto & fel_u = static_cast<const FEL&>(fel);
    size_t nd = fel_u.GetNDof();

    STACK_ARRAY(SIMD<double>, mem, nd*D*D);
    FlatMatrix<SIMD<double>> ddshape(nd, D*D, mem);

    for (size_t i = 0; i < bmir.Size(); i++)
      {
        fel_u.CalcMappedDDShape(bmir[i], ddshape);
        for (size_t j = 0; j < nd; j++)
          y(j) += HSum(InnerProduct(ddshape.Row(j), x.Col(i)));
      }


    /*
    
      size_t size = (bmir.Size()+1)*500*SIMD<double>::Size();
      STACK_ARRAY(char, data, size);
      LocalHeap lh(data, size);

      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (bmir);
      auto & ir = mir.IR();
      const ElementTransformation & trafo = mir.GetTransformation();
      auto & fel_u = static_cast<const FEL&>(fel);

      FlatMatrix<SIMD<double>> hx1(D, mir.Size(), lh);
      FlatMatrix<SIMD<double>> hx2(D, mir.Size(), lh);

      for (size_t j = 0; j < D-1; j++)
        {
          // hx = (F^-1 * x).Row(j)
          for (size_t k = 0; k < mir.Size(); k++)
            {
              auto jacinv = mir[k].GetJacobianInverse();
              for (int l = 0; l < D; l++)
                {
                  SIMD<double> sum = 0;
                  for (int m = 0; m < D; m++)
                    sum += jacinv(j,m) * x(m*D+l, k);
                  // hx.Get(l,k) = (-(0.5/eps()) * sum).Data();
                  hx1(l,k) = (-(8/(12*eps())) * sum).Data();
                  hx2(l,k) = ( (1/(12*eps())) * sum).Data();
                }
            }
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irl(mir.IR().GetNIP(), lh);
            for (size_t k = 0; k < irl.Size(); k++)
              {
                irl[k] = ir[k];
                irl[k](j) -= eps();
              }
            SIMD_MappedIntegrationRule<D-1,D> mirl(irl, trafo, lh);
            fel_u.AddGradTrans (mirl, hx1, y);
            irl.NothingToDelete();
          }
          {
            HeapReset hr(lh);
            hx1 *= -1;
            SIMD_IntegrationRule irr(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irr.Size(); k++)
              {
                irr[k] = ir[k];              
                irr[k](j) += eps();
              }
            SIMD_MappedIntegrationRule<D-1,D> mirr(irr, trafo, lh);
            fel_u.AddGradTrans (mirr, hx1, y);
          }
          {
            HeapReset hr(lh);
            SIMD_IntegrationRule irl(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irl.Size(); k++)
              {
                irl[k] = ir[k];
                irl[k](j) -= 2*eps();
              }
            SIMD_MappedIntegrationRule<D-1,D> mirl(irl, trafo, lh);
            fel_u.AddGradTrans (mirl, hx2, y);
          }
          {
            HeapReset hr(lh);
            hx2 *= -1;
            SIMD_IntegrationRule irr(mir.IR().GetNIP(), lh);
            for (int k = 0; k < irr.Size(); k++)
              {
                irr[k] = ir[k];              
                irr[k](j) += 2*eps();
              }
            SIMD_MappedIntegrationRule<D-1,D> mirr(irr, trafo, lh);
            fel_u.AddGradTrans (mirr, hx2, y);
          }
        }
    */
  }

  template class DiffOpHesseBoundary<3,ScalarFiniteElement<2>>;
  template class DiffOpHesseBoundary<2,ScalarFiniteElement<1>>;
  
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

