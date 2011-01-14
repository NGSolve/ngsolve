/*********************************************************************/
/* File:   hybridDG.cpp                                              */
/* Author: H. Egger, J. Schoeberl, RWTH                              */
/* Date:   10. Feb. 2008                                             */
/*********************************************************************/
  
/*  
    Finite Element Integrators 
*/
  
#include <fem.hpp>
#include <facetfe.hpp>  

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

    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void CalcElementMatrix (const FiniteElement & fel,
				    const ElementTransformation & eltrans, 
				    FlatMatrix<double> & elmat,
				    LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("HDG laplace");
      static int timer1 = NgProfiler::CreateTimer ("HDG laplace volume");
      static int timer1a = NgProfiler::CreateTimer ("HDG laplace volume, lapack");
      static int timer2 = NgProfiler::CreateTimer ("HDG laplace boundary");
      static int timer3 = NgProfiler::CreateTimer ("HDG laplace edge glue");

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
      
      {
	NgProfiler::RegionTimer reg (timer1);     
	HeapReset hr(lh);

	FlatMatrix<> mat_gradgrad (nd_l2, lh);
	FlatMatrixFixWidth<D> dshape(nd_l2, lh);

	const IntegrationRule & ir_vol =
	  SelectIntegrationRule (eltype, 2*fel_l2.Order());
	MappedIntegrationRule<D,D> mir_vol(ir_vol, eltrans, lh);

	
	FlatMatrix<> bmats(ir_vol.GetNIP()*D, nd_l2, lh);
	FlatMatrix<> dbmats(ir_vol.GetNIP()*D, nd_l2, lh);

	for (int l = 0; l < ir_vol.GetNIP(); l++)
	  {
	    const SpecificIntegrationPoint<D,D> & sip = mir_vol[l];
	    double lam = coef_lam->Evaluate(sip);
	    
	    fel_l2.CalcMappedDShape (sip, dshape);

	    bmats.Rows(l*D, (l+1)*D) = Trans(dshape);
	    dbmats.Rows(l*D, (l+1)*D) = (lam * sip.GetWeight()) * Trans(dshape);
	      // (lam * sip.GetMeasure() * ir_vol[l].Weight()) * Trans(dshape);
	  }

	NgProfiler::RegionTimer reg1a (timer1a);     

	// mat_gradgrad = Trans (dbmats) * bmats;
        LapackMultAtB (dbmats, bmats, mat_gradgrad);

	elmat.Cols(base_l2, base_l2+nd_l2).Rows(base_l2,base_l2+nd_l2) 
	  += mat_gradgrad;
      }


      // The facet contribution
      {
	NgProfiler::RegionTimer reg2 (timer2);     


	int nfacet = ElementTopology::GetNFacets(eltype);
      
	Facet2ElementTrafo transform(eltype); 
	const NORMAL * normals = ElementTopology::GetNormals(eltype);

	FlatMatrixFixHeight<2> bmat(nd, lh);
	FlatMatrixFixHeight<2> dbmat(nd, lh);
	Mat<2> dmat;

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

	    FlatMatrix<> comp_bmats(2*ir_facet.GetNIP(), facetdofs.Size(), lh);
	    FlatMatrix<> comp_dbmats(2*ir_facet.GetNIP(), facetdofs.Size(), lh);

	    FlatMatrix<> comp_elmat(facetdofs.Size(), facetdofs.Size(), lh);
	  
	    comp_elmat = 0;
	    bmat = 0.0;

	    fel_facet.SelectFacet (k);

	    for (int l = 0; l < ir_facet.GetNIP(); l++)
	      {
		IntegrationPoint ip = transform(k, ir_facet[l]);
		SpecificIntegrationPoint<D,D> sip (ip, eltrans, lh);
		double lam = coef_lam->Evaluate(sip);
              
		Mat<D> jac = sip.GetJacobian();
		Mat<D> inv_jac = sip.GetJacobianInverse();
		double det = sip.GetMeasure();


		Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
		double len = L2Norm (normal);
		normal /= len;

		fel_facet.CalcFacetShape(k, ir_facet[l], mat_facet);
// 		fel_facet.CalcShape(ip, mat_facet);
		fel_l2.CalcShape(ip, mat_l2);

		Vec<D> invjac_normal = inv_jac * normal;
		mat_dudn = fel_l2.GetDShape (sip.IP(), lh) * invjac_normal;
              
		bmat.Row(0).Range (base_l2   , base_l2   +nd_l2 )   = mat_dudn;
		bmat.Row(1).Range (base_l2   , base_l2   +nd_l2   ) = mat_l2;
		bmat.Row(1).Range (base_facet, base_facet+nd_facet) = -mat_facet;

		for (int i = 0; i < facetdofs.Size(); i++)
		  comp_bmat.Col(i) = bmat.Col(facetdofs[i]);

		dmat(0,0) = 0;
		dmat(1,0) = dmat(0,1) = -1;
		// dmat(1,1) = alpha * sqr (fel_l2.Order()) * (len/det);
		dmat(1,1) = alpha * ((fel_l2.Order()+1)*(fel_l2.Order()+D)/D * len) *(1.0/det);


		dmat *= lam * len * ir_facet[l].Weight();

		comp_bmats.Rows (2*l, 2*l+2) = comp_bmat;
		comp_dbmats.Rows (2*l, 2*l+2) = dmat * comp_bmat;

		NgProfiler::AddFlops (timer2, nd*nd*4);
	      }
	    // comp_elmat = Trans (comp_bmats) * comp_dbmats;
	    LapackMultAtB (comp_bmats, comp_dbmats, comp_elmat);

	    for (int i = 0; i < facetdofs.Size(); i++)
	      for (int j = 0; j < facetdofs.Size(); j++)
		elmat(facetdofs[i], facetdofs[j]) += comp_elmat(i,j);
	  }
      }

      // the vertex glue
      if (const_cast<CompoundFiniteElement&>(cfel).GetNComponents() == 3)
	{
	  NgProfiler::RegionTimer reg3 (timer3);     
	  if (D == 2)
	    {
	      const ScalarFiniteElement<D> & fel_h1 = 
		dynamic_cast<const ScalarFiniteElement<D> &> (cfel[2]);

	      int nd_h1 = fel_h1.GetNDof();
	      int base_h1 = nd_l2 + nd_facet;
	      FlatVector<> vec_h1(nd_h1, lh);
	      FlatVector<> b_vec(nd, lh);
	      
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
	      HeapReset hr(lh);

	      const EdgeVolumeFiniteElement<D> & fel_edge = 
		dynamic_cast<const EdgeVolumeFiniteElement<D> &> (cfel[2]);

	      int nd_edge = fel_edge.GetNDof();
	      int base_edge = nd_l2 + nd_facet;
	      FlatVector<> vec_edge(nd_edge, lh);
	      FlatVector<> b_vec(nd, lh);
	      FlatVector<> db_vec(nd, lh);

	      const POINT3D * verts = ElementTopology::GetVertices (eltype);
	      const EDGE * edges = ElementTopology::GetEdges (eltype);
	      int ned = ElementTopology::GetNEdges(eltype);

	      const IntegrationRule & ir1d = SelectIntegrationRule (ET_SEGM, 2 * fel_l2.Order());

	      FlatMatrix<> bmats(ir1d.Size(), nd, lh);
	      FlatMatrix<> dbmats(ir1d.Size(), nd, lh);

	      for (int i = 0; i < ned; i++)
		{
		  Vec<3> p1, p2;
		  for (int j = 0; j < 3; j++)
		    { 
		      p1(j) = verts[edges[i][0]][j];
		      p2(j) = verts[edges[i][1]][j];
		    }
		  
		  for (int k = 0; k < ir1d.Size(); k++)
		    {
		      Vec<3> p = p1 + ir1d[k](0) * (p2-p1);
		      IntegrationPoint ip (p(0), p(1), p(2), ir1d[k].Weight());
		      SpecificIntegrationPoint<D,D> sip (ip, eltrans, lh);
		      double lam = coef_lam->Evaluate(sip);

		      
		      Vec<3> tau = sip.GetJacobian() * (p2-p1);
		      double h = L2Norm(tau);

		      fel_l2.CalcShape(ip, mat_l2);
		      fel_edge.CalcEdgeShape(i, ip, vec_edge);
		      
		      b_vec = 0.0;
		      b_vec.Range(0, nd_l2) = mat_l2;
		      b_vec.Range(base_edge, nd) = -vec_edge;
		      
		      db_vec = (alpha * lam * h * ir1d[k].Weight()) * b_vec;
		      // elmat += db_vec * Trans(b_vec);
		      bmats.Row(k) = b_vec;
		      dbmats.Row(k) = db_vec;
		    }
		  LapackMultAddAtB (bmats, dbmats, 1, elmat);
		  // elmat += Trans (bmats) * dbmats;
		}
	    }
	}
    }


    virtual void ApplyElementMatrix (const FiniteElement & fel,
				     const ElementTransformation & eltrans, 
				     const FlatVector<double> & elx, 
				     FlatVector<double> & ely,
				     void * precomputed,
				     LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("HDG apply laplace");
      static int timer1 = NgProfiler::CreateTimer ("HDG apply laplace volume");
      static int timer2 = NgProfiler::CreateTimer ("HDG apply laplace boundary");
      static int timer3 = NgProfiler::CreateTimer ("HDG apply laplace edge glue");

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


      ely = 0.0;

      {
	NgProfiler::RegionTimer reg (timer1);     
	HeapReset hr(lh);


	IntegrationRuleTP<D> ir_vol(eltrans, 2*fel_l2.Order(), false, lh);
	// const IntegrationRule & ir_vol = SelectIntegrationRule (eltype, 2*fel_l2.Order());


	MappedIntegrationRule<D,D> mir_vol(ir_vol, eltrans, lh);

	FlatMatrixFixWidth<D> grad(ir_vol.GetNIP(), lh);
	
	fel_l2.EvaluateGrad (ir_vol, elx.Range(base_l2, base_l2+nd_l2), grad);
	

	for (int l = 0; l < ir_vol.GetNIP(); l++)
	  {
	    const SpecificIntegrationPoint<D,D> & sip = mir_vol[l];
	    double lam = coef_lam->Evaluate(sip);
	    
	    Vec<D> gi = grad.Row(l);
	    Vec<D> hv1 = Trans (sip.GetJacobianInverse()) * gi;
	    Vec<D> hv2 = sip.GetJacobianInverse() * hv1;
	    gi = (lam * sip.GetJacobiDet() * ir_vol[l].Weight()) * hv2;

	    grad.Row(l) = gi;
	  }

	fel_l2.EvaluateGradTrans (ir_vol, grad, ely.Range(base_l2, base_l2+nd_l2));	
      }

      
      // The facet contribution
      {
	NgProfiler::RegionTimer reg2 (timer2);     


	int nfacet = ElementTopology::GetNFacets(eltype);

	int sort[D+1];
	eltrans.GetSort (FlatArray<int> (D+1,&sort[0]));
      
	Facet2ElementTrafo transform(eltype);
	const NORMAL * normals = ElementTopology::GetNormals(eltype);

	Mat<2> dmat;

	for (int k = 0; k < nfacet; k++)

	  {
	    HeapReset hr(lh);

	    fel_facet.SelectFacet (k);

	    Vec<D> normal_ref;
	    for (int i=0; i<D; i++)
	      normal_ref(i) = normals[k][i];

	    /*
	    ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);

	    const IntegrationRule & ir_facet =
	      SelectIntegrationRule (etfacet, 2*fel_l2.Order());

	    IntegrationRule ir_vol(ir_facet.Size());
	    

	    for (int l = 0; l < ir_facet.GetNIP(); l++)
	      {
		IntegrationPoint ip = transform(k, ir_facet[l]);
		ir_vol.Append (ip);
	      }
	    */

	    int sort[4];
	    FlatArray<int> fa_sort(D+1, &sort[0]);
	    eltrans.GetSort (fa_sort);
	    IntegrationRuleTP<D> ir_vol (eltype, fa_sort, NODE_TYPE(D-1), k, 2*fel_l2.Order(), lh);


	    MappedIntegrationRule<D,D> mir(ir_vol, eltrans, lh);

	    FlatVector<> shapes_l2(ir_vol.Size(), lh);
	    FlatVector<> shapes_facet(ir_vol.Size(), lh);
	    FlatMatrixFixWidth<D> grad_l2(ir_vol.Size(), lh);


	    fel_l2.Evaluate (ir_vol, elx.Range(base_l2, base_l2+nd_l2), shapes_l2);
	    fel_l2.EvaluateGrad (ir_vol, elx.Range(base_l2, base_l2+nd_l2), grad_l2);
	    fel_facet.Evaluate (ir_vol, elx.Range(base_facet, base_facet+nd_facet), shapes_facet);



	    for (int l = 0; l < ir_vol.GetNIP(); l++)
	      {
		const SpecificIntegrationPoint<D,D> & sip = mir[l];
		double lam = coef_lam->Evaluate(sip);
              
		Mat<D> jac = sip.GetJacobian();
		Mat<D> inv_jac = sip.GetJacobianInverse();
		double det = sip.GetJacobiDet();


		Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
		double len = L2Norm (normal);
		normal /= len;

		Vec<D> invjac_normal = inv_jac * normal;
              
		dmat(0,0) = 0;
		dmat(1,0) = dmat(0,1) = -1;
		dmat(1,1) = alpha * ((fel_l2.Order()+1)*(fel_l2.Order()+D)/D * len) *(1.0/det);

		dmat *= lam * len * ir_vol[l].Weight();

		Vec<2> hv1, hv2;
		hv1(0) = InnerProduct (grad_l2.Row(l), invjac_normal);
		hv1(1) = shapes_l2(l) - shapes_facet(l);

		hv2 = dmat * hv1;

		shapes_l2(l) = hv2(1);
		shapes_facet(l) = -hv2(1);
		grad_l2.Row(l) = hv2(0) * invjac_normal;
	      }


	    FlatVector<> hely(nd_l2, lh);
	    fel_l2.EvaluateTrans (ir_vol, shapes_l2, hely);
	    ely.Range (base_l2, base_l2+nd_l2) += hely;

	    fel_l2.EvaluateGradTrans (ir_vol, grad_l2, hely);
	    ely.Range (base_l2, base_l2+nd_l2) += hely;

	    FlatVector<> hely_facet(nd_facet, lh);
	    fel_facet.EvaluateTrans (ir_vol, shapes_facet, hely_facet);
	    ely.Range (base_facet, base_facet+nd_facet) += hely_facet;
	  }
      }


      // the vertex glue
      if (const_cast<CompoundFiniteElement&>(cfel).GetNComponents() == 3)
	{
	  NgProfiler::RegionTimer reg3 (timer3);     
	  if (D == 2)
	    {
	      const ScalarFiniteElement<D> & fel_h1 = 
		dynamic_cast<const ScalarFiniteElement<D> &> (cfel[2]);

	      int nd_h1 = fel_h1.GetNDof();
	      int base_h1 = nd_l2 + nd_facet;
	      FlatVector<> vec_h1(nd_h1, lh);
	      FlatVector<> b_vec(nd, lh);
	      FlatVector<> mat_l2(nd_l2, lh);
	      

	      const POINT3D * verts = ElementTopology::GetVertices (eltype);
	      int nv = ElementTopology::GetNVertices(eltype);

	      double scale = alpha; // elmat(0,0);
	      FlatMatrix<> elmat (nd, lh);
	      elmat = 0.0;

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
	      ely += elmat * elx;
	    }
	  else
	    {
	      HeapReset hr(lh);

	      const EdgeVolumeFiniteElement<D> & fel_edge = 
		dynamic_cast<const EdgeVolumeFiniteElement<D> &> (cfel[2]);

	      int nd_edge = fel_edge.GetNDof();
	      int base_edge = nd_l2 + nd_facet;


	      const POINT3D * verts = ElementTopology::GetVertices (eltype);
	      const EDGE * edges = ElementTopology::GetEdges (eltype);
	      int ned = ElementTopology::GetNEdges(eltype);


	      for (int i = 0; i < ned; i++)
		{
		  Vec<3> p1, p2;
		  for (int j = 0; j < 3; j++)
		    { 
		      p1(j) = verts[edges[i][0]][j];
		      p2(j) = verts[edges[i][1]][j];
		    }


		  /*
		    const IntegrationRule & ir1d = SelectIntegrationRule (ET_SEGM, 2 * fel_l2.Order());
		  IntegrationRule ir_vol(ir1d.Size());

		  for (int l = 0; l < ir1d.GetNIP(); l++)
		    {
		      Vec<3> p = p1 + ir1d[l](0) * (p2-p1);
		      IntegrationPoint ip (p(0), p(1), p(2), ir1d[l].Weight());
		      ir_vol.Append (ip);
		    }
		  */

		  int sort[4];
		  FlatArray<int> fa_sort(D+1, &sort[0]);
		  eltrans.GetSort (fa_sort);
		  IntegrationRuleTP<D> ir_vol (eltype, fa_sort, NODE_TYPE(1), i, 2*fel_l2.Order(), lh);



		  MappedIntegrationRule<D,D> mir(ir_vol, eltrans, lh);
			    
		  FlatVector<> shapes_l2(ir_vol.Size(), lh);
		  FlatVector<> shapes_edge(ir_vol.Size(), lh);



		  fel_edge.SelectEdge(i);

		  fel_l2.Evaluate (ir_vol, elx.Range(base_l2, base_l2+nd_l2), shapes_l2);
		  fel_edge.Evaluate (ir_vol, elx.Range(base_edge, base_edge+nd_edge), shapes_edge);

		  shapes_l2 -= shapes_edge;

		  for (int k = 0; k < ir_vol.Size(); k++)
		    {
		      const SpecificIntegrationPoint<D,D> & sip = mir[k];

		      double lam = coef_lam->Evaluate(sip);
		      
		      Vec<3> tau = sip.GetJacobian() * (p2-p1);
		      double h = L2Norm(tau);

		      shapes_l2(k) *=  alpha * lam * h * ir_vol[k].Weight();
		    }


		  FlatVector<> hv_l2(nd_l2, lh);
		  FlatVector<> hv_edge(nd_edge, lh);
		  fel_l2.EvaluateTrans (ir_vol, shapes_l2, hv_l2);
		  fel_edge.EvaluateTrans (ir_vol, shapes_l2, hv_edge);

		  ely.Range (0, nd_l2) += hv_l2;
		  ely.Range (base_edge, base_edge+nd_edge) -= hv_edge;
		}
	    }
	}
    }
  };






  template <int D>
  class HDG_ConvectionIntegrator : public BilinearFormIntegrator
  {
  protected:
    Array<CoefficientFunction *> coef_conv;
  public:
    HDG_ConvectionIntegrator (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      coef_conv.SetSize(coeffs.Size());
      for (int j = 0; j < coeffs.Size(); j++)
        coef_conv[j] = coeffs[j];
    }

    virtual ~HDG_ConvectionIntegrator () { ; }

    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void CalcElementMatrix (const FiniteElement & fel,
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

          if (coef_conv.Size()>1)
	    for (int j = 0; j < D; j++)
	      conv(j) = coef_conv[j]->Evaluate(sip);
	  else
	    coef_conv[0]->Evaluate(sip,conv);

          fel_l2.CalcShape (sip.IP(), shape);
          fel_l2.CalcMappedDShape (sip, dshape);
          
          conv_dshape = dshape * conv;

          conv_dshape *= sip.GetJacobiDet() * ir_vol[l].Weight();
          
          // mat_l2 = shape * Trans (conv_dshape);
          mat_l2 -= conv_dshape * Trans (shape);
        }
      
      elmat.Cols(base_l2, base_l2+nd_l2).Rows(base_l2,base_l2+nd_l2) 
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
	      if (coef_conv.Size()>1)
		for (int j = 0; j < D; j++)
		  conv(j) = coef_conv[j]->Evaluate(sip);
	      else
		coef_conv[0]->Evaluate(sip,conv);


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
    //TODO: Vertex glue
  };



  static RegisterBilinearFormIntegrator<HDG_LaplaceIntegrator<2> > initlap2 ("HDG_laplace", 2, 2);
  static RegisterBilinearFormIntegrator<HDG_LaplaceIntegrator<3> > initlap3 ("HDG_laplace", 3, 2);

  static RegisterBilinearFormIntegrator<HDG_ConvectionIntegrator<2> > initconv21 ("HDG_convection", 2, 1);
  static RegisterBilinearFormIntegrator<HDG_ConvectionIntegrator<2> > initconv22 ("HDG_convection", 2, 2);

  static RegisterBilinearFormIntegrator<HDG_ConvectionIntegrator<3> > initconv31 ("HDG_convection", 2, 1);
  static RegisterBilinearFormIntegrator<HDG_ConvectionIntegrator<3> > initconv33 ("HDG_convection", 3, 3);
}

