/*********************************************************************/
/* File:   hybridDG.cpp                                              */
/* Author: RWTH, H. Egger, J. Schoeberl                              */
/* Date:   10. Feb. 2008                                             */
/*********************************************************************/
  
/*  
   Finite Element Integrators 
*/
  
#include <fem.hpp>
  
namespace ngfem
{
  using namespace ngfem;

  template <int D>
  class HDG_LaplaceIntegrator : public BilinearFormIntegrator
  {
  protected:
    CoefficientFunction *coef_lam;
  public:
    HDG_LaplaceIntegrator (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      coef_lam  = coeffs[0];
    }

    virtual ~HDG_LaplaceIntegrator () { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new HDG_LaplaceIntegrator (coeffs);
    }
    
    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void AssembleElementMatrix (const FiniteElement & fel,
                                        const ElementTransformation & eltrans, 
                                        FlatMatrix<double> & elmat,
                                        LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("HDG laplace");
      static int timer2 = NgProfiler::CreateTimer ("HDG laplace boundary");

      NgProfiler::RegionTimer reg (timer);


      const CompoundFiniteElement & cfel = 
        dynamic_cast<const CompoundFiniteElement&> (fel);
    
      const L2HighOrderFiniteElement<D> & fel_l2 = 
        dynamic_cast<const L2HighOrderFiniteElement<D>&> (cfel[0]);
          const FacetVolumeFiniteElement<D> & fel_facet = 
        dynamic_cast<const FacetVolumeFiniteElement<D> &> (cfel[1]);
    
  
      ELEMENT_TYPE eltype = cfel.ElementType();
      
      int nd_l2 = fel_l2.GetNDof();
      int nd_facet = fel_facet.GetNDof();
      int nd = nd_l2 + nd_facet;

      int base_l2 = 0;
      int base_facet = base_l2+nd_l2;
      
      elmat.AssignMemory(nd, nd, lh);
      elmat = 0.0;

        
      FlatVector<> mat_l2(nd_l2, lh);
      FlatVector<> mat_dudn(nd_l2, lh);
      FlatVector<> mat_facet(nd_facet, lh);
      
      FlatMatrixFixHeight<2> bmat(nd, lh);
      FlatMatrixFixHeight<2> dbmat(nd, lh);
      Mat<2> dmat;

      FlatMatrix<> mat_gradgrad (nd_l2, lh);
      FlatMatrixFixWidth<D> dshape(nd_l2, lh);
      FlatMatrixFixWidth<D> fac_dshape(nd_l2, lh);

      double lam = coef_lam->EvaluateConst();
          
      const IntegrationRule & ir_vol =
        SelectIntegrationRule (eltype, 2*fel_l2.Order());
      
      mat_gradgrad = 0.0;

      for (int l = 0; l < ir_vol.GetNIP(); l++)
        {
          HeapReset hr(lh);
          const SpecificIntegrationPoint<D,D> sip(ir_vol[l], eltrans, lh);
          
          fel_l2.CalcMappedDShape (sip, dshape);
          fac_dshape = (lam * sip.GetJacobiDet() * ir_vol[l].Weight()) * dshape;

          mat_gradgrad += fac_dshape * Trans (dshape);
        }

      elmat.HRange(base_l2, base_l2+nd_l2).VRange(base_l2,base_l2+nd_l2) 
        += mat_gradgrad;
      

      // The facet contribution
      NgProfiler::RegionTimer reg2 (timer2);     


      int nfacet = ElementTopology::GetNFacets(eltype);
      
      Facet2ElementTrafo transform(eltype); 
      const NORMAL * normals = ElementTopology::GetNormals(eltype);

      for (int k = 0; k < nfacet; k++)
        {
          HeapReset hr(lh);
          ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);

          Vec<D> normal_ref;
          for (int i=0; i<D; i++)
            normal_ref(i) = normals[k][i];

          const IntegrationRule & ir_facet =
            SelectIntegrationRule (etfacet, 2*max (0,fel_l2.Order()));
      

	  Array<int> facetdofs, fdofs;
	  for (int i = 0; i < nd_l2; i++)
	    facetdofs.Append(i);

	  fel_facet.GetFacetDofNrs(k, fdofs);
	  for (int i = 0; i < fdofs.Size(); i++)
	    facetdofs.Append(base_facet+fdofs[i]);

	  FlatMatrixFixHeight<2> comp_bmat(facetdofs.Size(), lh);
	  FlatMatrixFixHeight<2> comp_dbmat(facetdofs.Size(), lh);
	  FlatMatrix<> comp_elmat(facetdofs.Size(), facetdofs.Size(), lh);
	  
	  comp_elmat = 0;
          bmat = 0.0;

          for (int l = 0; l < ir_facet.GetNIP(); l++)
            {
              IntegrationPoint ip = transform(k, ir_facet[l]);
              SpecificIntegrationPoint<D,D> sip (ip, eltrans, lh);
	    
              Mat<D> jac = sip.GetJacobian();
              Mat<D> inv_jac = sip.GetJacobianInverse();
              double det = sip.GetJacobiDet();


              Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
              double len = L2Norm (normal);
              normal /= len;

              fel_facet.CalcFacetShape(k, ir_facet[l], mat_facet);
              fel_l2.CalcShape(sip.IP(), mat_l2);

              Vec<D> invjac_normal = inv_jac * normal;
              mat_dudn = fel_l2.GetDShape (sip.IP(), lh) * invjac_normal;
              
              bmat.Row(0).Range (base_l2   , base_l2   +nd_l2 )   = mat_dudn;
              bmat.Row(1).Range (base_l2   , base_l2   +nd_l2   ) = mat_l2;
              bmat.Row(1).Range (base_facet, base_facet+nd_facet) = -mat_facet;

	      for (int i = 0; i < facetdofs.Size(); i++)
		comp_bmat.Col(i) = bmat.Col(facetdofs[i]);

              dmat(0,0) = 0;
              dmat(1,0) = dmat(0,1) = -1;
              dmat(1,1) = 10 * sqr (fel_l2.Order()) * (len/det);


              dmat *= len * ir_facet[l].Weight();
              comp_dbmat = dmat * comp_bmat;
              comp_elmat += Trans (comp_bmat) * comp_dbmat;

              NgProfiler::AddFlops (timer2, nd*nd*4);
            }
	  
	  for (int i = 0; i < facetdofs.Size(); i++)
	    for (int j = 0; j < facetdofs.Size(); j++)
	      elmat(facetdofs[i], facetdofs[j]) += comp_elmat(i,j);
        }
    }
  };
  

  
  namespace inithdg
  {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetIntegrators().AddBFIntegrator ("HDG_laplace", 2, 1,
                                        HDG_LaplaceIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("HDG_laplace", 3, 1,
                                        HDG_LaplaceIntegrator<3>::Create);
    }
    
    Init init;
  }

}

