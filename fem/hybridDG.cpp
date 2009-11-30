/*********************************************************************/
/* File:   hybridDG.cpp                                              */
/* Author: H. Egger, J. Schoeberl, RWTH                              */
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
    double alpha;   // interior penalyty
    CoefficientFunction *coef_lam;
  public:
    HDG_LaplaceIntegrator (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      coef_lam  = coeffs[0];
      alpha = coeffs[1] -> EvaluateConst();
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
      const ScalarFiniteElement<D> & fel_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>&> (cfel[0]);
      const FacetVolumeFiniteElement<D> & fel_facet = 
        dynamic_cast<const FacetVolumeFiniteElement<D> &> (cfel[1]);
      
  
      ELEMENT_TYPE eltype = cfel.ElementType();
      
      int nd_l2 = fel_l2.GetNDof();
      int nd_facet = fel_facet.GetNDof();
      int nd = cfel.GetNDof();  

      int base_l2 = 0;
      int base_facet = base_l2+nd_l2;
      
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

          
      const IntegrationRule & ir_vol =
        SelectIntegrationRule (eltype, 2*fel_l2.Order());
      
      mat_gradgrad = 0.0;

      for (int l = 0; l < ir_vol.GetNIP(); l++)
        {
          HeapReset hr(lh);
          const SpecificIntegrationPoint<D,D> sip(ir_vol[l], eltrans, lh);
          double lam = coef_lam->Evaluate(sip);
          
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
              double lam = coef_lam->Evaluate(sip);
              
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
              dmat(1,1) = alpha * sqr (fel_l2.Order()) * (len/det);


              dmat *= lam * len * ir_facet[l].Weight();
              comp_dbmat = dmat * comp_bmat;
              comp_elmat += Trans (comp_bmat) * comp_dbmat;

              NgProfiler::AddFlops (timer2, nd*nd*4);
            }

	  for (int i = 0; i < facetdofs.Size(); i++)
	    for (int j = 0; j < facetdofs.Size(); j++)
	      elmat(facetdofs[i], facetdofs[j]) += comp_elmat(i,j);
        }


      // the vertex glue
      if (const_cast<CompoundFiniteElement&>(cfel).GetNComponents() == 3)
	{
	  const ScalarFiniteElement<D> & fel_h1 = 
	    dynamic_cast<const ScalarFiniteElement<D> &> (cfel[2]);

	  int nd_h1 = fel_h1.GetNDof();
	  int base_h1 = nd_l2 + nd_facet;
	  FlatVector<> vec_h1(nd_h1, lh);
	  FlatVector<> b_vec(nd, lh);
	 
	  if (D == 2)
	    {
	      
	      const POINT3D * verts = ElementTopology::GetVertices (eltype);
	      int nv = ElementTopology::GetNVertices(eltype);

	      double scale = elmat(0,0);
	      
	      for (int i = 0; i < nv; i++)
		{
		  IntegrationPoint ip;
		  for (int j = 0; j < D; j++)
		    ip(j) = verts[i][j];
		  
		  fel_l2.CalcShape(ip, mat_l2);
		  fel_h1.CalcShape(ip, vec_h1);
		  
		  b_vec = 0.0;
		  b_vec.Range(0, nd_l2) = mat_l2;
		  b_vec.Range(base_h1, nd) = -vec_h1;
		  
		  elmat += scale * (b_vec * Trans(b_vec));
		}
	    }
	  else
	    {
	      const POINT3D * verts = ElementTopology::GetVertices (eltype);
	      // int nv = ElementTopology::GetNVertices(eltype);
	      const EDGE * edges = ElementTopology::GetEdges (eltype);
	      int ned = ElementTopology::GetNEdges(eltype);

	      const IntegrationRule & ir1d = SelectIntegrationRule (ET_SEGM, 2 * fel_l2.Order());


	      double scale = 10 * elmat(0,0);
	      // *testout << "scale = " << scale << endl;
	      for (int i = 0; i < ned; i++)
		{
		  Vec<3> p1, p2;
		  for (int j = 0; j < 3; j++)
		    { 
		      p1(j) = verts[edges[i][0]][j];
		      p2(j) = verts[edges[i][1]][j];
		    }
		  
		  // *testout << "p1 = " << p1 << ", p2 = " << p2 << endl;
		  for (int k = 0; k < ir1d.Size(); k++)
		    {
		      Vec<3> p = p1 + ir1d[k](0) * (p2-p1);
		      IntegrationPoint ip (p(0), p(1), p(2), ir1d[k].Weight());
		      SpecificIntegrationPoint<D,D> sip (ip, eltrans, lh);
		      
		      fel_l2.CalcShape(ip, mat_l2);
		      fel_h1.CalcShape(ip, vec_h1);
		      
		      // *testout << "p = " << p << ", vec_h1 = " << vec_h1 << endl;
		      b_vec = 0.0;
		      b_vec.Range(0, nd_l2) = mat_l2;
		      b_vec.Range(base_h1, nd) = -vec_h1;
		      
		      // *testout << "bvec = " << b_vec << endl;

		      elmat += scale * ir1d[k].Weight() * (b_vec * Trans(b_vec));
		    }
		}
	    }
	}
    }
  };




  template <int D>
  class HDG_ConvectionIntegrator : public BilinearFormIntegrator
  {
  protected:
    CoefficientFunction * coef_conv[D];
  public:
    HDG_ConvectionIntegrator (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      for (int j = 0; j < D; j++)
        coef_conv[j] = coeffs[j];
    }

    virtual ~HDG_ConvectionIntegrator () { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new HDG_ConvectionIntegrator (coeffs);
    }
    
    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void AssembleElementMatrix (const FiniteElement & fel,
                                        const ElementTransformation & eltrans, 
                                        FlatMatrix<double> & elmat,
                                        LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("HDG convection");
      static int timer2 = NgProfiler::CreateTimer ("HDG convection boundary");

      NgProfiler::RegionTimer reg (timer);

      const CompoundFiniteElement & cfel = 
        dynamic_cast<const CompoundFiniteElement&> (fel);
    
      const ScalarFiniteElement<D> & fel_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>&> (cfel[0]);
          const FacetVolumeFiniteElement<D> & fel_facet = 
        dynamic_cast<const FacetVolumeFiniteElement<D> &> (cfel[1]);
    
  
      ELEMENT_TYPE eltype = cfel.ElementType();
      
      int nd_l2 = fel_l2.GetNDof();
      int nd_facet = fel_facet.GetNDof();
      int nd = nd_l2 + nd_facet;

      int base_l2 = 0;
      int base_facet = base_l2+nd_l2;
      
      elmat = 0.0;

      FlatVector<> shape(nd_l2, lh);
      FlatVector<> conv_dshape(nd_l2, lh);
      FlatVector<> mat_dudn(nd_l2, lh);
      FlatVector<> shape_facet(nd_facet, lh);
      
      FlatMatrixFixHeight<2> bmat(nd, lh);
      FlatMatrixFixHeight<2> dbmat(nd, lh);
      Mat<2> dmat;

      FlatMatrix<> mat_l2 (nd_l2, lh);
      FlatMatrixFixWidth<D> dshape(nd_l2, lh);
      FlatMatrixFixWidth<D> fac_dshape(nd_l2, lh);

          
      const IntegrationRule & ir_vol =
        SelectIntegrationRule (eltype, 2*fel_l2.Order());
      
      mat_l2 = 0.0;

      for (int l = 0; l < ir_vol.GetNIP(); l++)
        {
          HeapReset hr(lh);
          const SpecificIntegrationPoint<D,D> sip(ir_vol[l], eltrans, lh);
          Vec<D> conv;
          for (int j = 0; j < D; j++)
            conv(j) = coef_conv[j]->Evaluate(sip);
          
          fel_l2.CalcShape (sip.IP(), shape);
          fel_l2.CalcMappedDShape (sip, dshape);
          
          conv_dshape = dshape * conv;

          conv_dshape *= sip.GetJacobiDet() * ir_vol[l].Weight();
          
          // mat_l2 = shape * Trans (conv_dshape);
          mat_l2 -= conv_dshape * Trans (shape);
        }
      
      elmat.HRange(base_l2, base_l2+nd_l2).VRange(base_l2,base_l2+nd_l2) 
        = mat_l2;
      // = Trans (mat_l2);

      
      // The facet contribution
      int nfacet = ElementTopology::GetNFacets(eltype);
      
      Facet2ElementTrafo transform(eltype); 
      const NORMAL * normals = ElementTopology::GetNormals(eltype);

      NgProfiler::RegionTimer reg2 (timer2);     

      for (int k = 0; k < nfacet; k++)
        {
          HeapReset hr(lh);
          ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);

          Vec<D> normal_ref;
          for (int i=0; i<D; i++)
            normal_ref(i) = normals[k][i];

          const IntegrationRule & ir_facet =
            SelectIntegrationRule (etfacet, fel_l2.Order()+fel_facet.Order());
      

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

              Vec<D> conv;
              for (int j = 0; j < D; j++)
                conv(j) = coef_conv[j]->Evaluate(sip);

              Mat<D> jac = sip.GetJacobian();
              Mat<D> inv_jac = sip.GetJacobianInverse();
              double det = sip.GetJacobiDet();

              Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
              double len = L2Norm (normal);
              normal /= len;

              double bn = InnerProduct (conv, normal);
              bool inflow = (bn < 0);

              fel_facet.CalcFacetShape(k, ir_facet[l], shape_facet);
              fel_l2.CalcShape(sip.IP(), shape);
              
              bmat.Row(0).Range (base_l2   , base_l2   +nd_l2   ) = shape;
              bmat.Row(1).Range (base_facet, base_facet+nd_facet) = shape_facet;

	      for (int i = 0; i < facetdofs.Size(); i++)
		comp_bmat.Col(i) = bmat.Col(facetdofs[i]);


              dmat = 0.0;

              if (inflow)
                dmat(0,1) = bn;
              else
                dmat(0,0) = bn;

              if (!inflow)
                {
                  dmat(1,0) = -bn; 
                  dmat(1,1) = bn;
                }



              dmat *= len * ir_facet[l].Weight();
              comp_dbmat = dmat * comp_bmat;
              comp_elmat += Trans (comp_bmat) * comp_dbmat;
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
      GetIntegrators().AddBFIntegrator ("HDG_laplace", 2, 2,
                                        HDG_LaplaceIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("HDG_laplace", 3, 2,
                                        HDG_LaplaceIntegrator<3>::Create);

      GetIntegrators().AddBFIntegrator ("HDG_convection", 2, 2,
                                        HDG_ConvectionIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("HDG_convection", 3, 3,
                                        HDG_ConvectionIntegrator<3>::Create);

    }
    
    Init init;
  }

}

