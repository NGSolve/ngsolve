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
      static Timer timer ("HDG laplace");
      static Timer timer1 ("HDG laplace volume");
      static Timer timer1a ("HDG laplace volume, lapack");
      static Timer timer2 ("HDG laplace boundary");

      RegionTimer reg (timer);

      const CompoundFiniteElement & cfel = 
        dynamic_cast<const CompoundFiniteElement&> (fel);
      const ScalarFiniteElement<D> & fel_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>&> (cfel[0]);
      const FacetVolumeFiniteElement<D> & fel_facet = 
        dynamic_cast<const FacetVolumeFiniteElement<D> &> (cfel[1]);
      
  
      ELEMENT_TYPE eltype = cfel.ElementType();
      
      IntRange l2_dofs = cfel.GetRange (0);
      IntRange facet_dofs = cfel.GetRange (1);

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

	IntegrationRule ir_vol(eltype, 2*fel_l2.Order());
	MappedIntegrationRule<D,D> mir_vol(ir_vol, eltrans, lh);
	
	FlatMatrix<> bmats(ir_vol.GetNIP()*D, nd_l2, lh);
	FlatMatrix<> dbmats(ir_vol.GetNIP()*D, nd_l2, lh);

	for (int l = 0; l < ir_vol.GetNIP(); l++)
	  {
	    const MappedIntegrationPoint<D,D> & mip = mir_vol[l];
	    double lam = coef_lam->Evaluate(mip);
	    
	    fel_l2.CalcMappedDShape (mip, dshape);

	    bmats.Rows(l*D, (l+1)*D) = Trans(dshape);
	    dbmats.Rows(l*D, (l+1)*D) = (lam * mip.GetWeight()) * Trans(dshape);
	  }

	NgProfiler::RegionTimer reg1a (timer1a);     

	// mat_gradgrad = Trans (dbmats) * bmats;
        LapackMultAtB (dbmats, bmats, mat_gradgrad);

	elmat.Cols(l2_dofs).Rows(l2_dofs) += mat_gradgrad;
      }


      // The facet contribution
      {
	NgProfiler::RegionTimer reg2 (timer2);     

	int nfacet = ElementTopology::GetNFacets(eltype);
      
	Facet2ElementTrafo transform(eltype); 
	FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);

	FlatMatrixFixHeight<2> bmat(nd, lh);
	FlatMatrixFixHeight<2> dbmat(nd, lh);
	Mat<2> dmat;

	for (int k = 0; k < nfacet; k++)
	  {
	    HeapReset hr(lh);
	    ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);

	    Vec<D> normal_ref = normals[k];

	    IntegrationRule ir_facet(etfacet, 2*fel_l2.Order());
	    IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);

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
	    
	    MappedIntegrationRule<D,D> mir(ir_facet_vol, eltrans, lh);

	    for (int l = 0; l < ir_facet.GetNIP(); l++)
	      {
		MappedIntegrationPoint<D,D> & mip = mir[l];

		double lam = coef_lam->Evaluate(mip);
              
		Mat<D> inv_jac = mip.GetJacobianInverse();
		double det = mip.GetMeasure();


		Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
		double len = L2Norm (normal);
		normal /= len;

		fel_facet.CalcFacetShape(k, ir_facet[l], mat_facet);
		// fel_facet.CalcShape(ip, mat_facet);
		fel_l2.CalcShape(ir_facet_vol[l], mat_l2);

		Vec<D> invjac_normal = inv_jac * normal;
		mat_dudn = fel_l2.GetDShape (mip.IP(), lh) * invjac_normal;
		
		bmat.Row(0).Range (l2_dofs)    = mat_dudn;
		bmat.Row(1).Range (l2_dofs)    = mat_l2;
		bmat.Row(1).Range (facet_dofs) = -mat_facet;

		for (int i = 0; i < facetdofs.Size(); i++)
		  comp_bmat.Col(i) = bmat.Col(facetdofs[i]);

		dmat(0,0) = 0;
		dmat(1,0) = dmat(0,1) = -1;
		// dmat(1,1) = alpha * sqr (fel_l2.Order()) * (len/det);
		dmat(1,1) = alpha * ((fel_l2.Order()+1)*(fel_l2.Order()+D)/D) * (len/det);

		dmat *= lam * len * ir_facet[l].Weight();

		comp_bmats.Rows (2*l, 2*l+2) = comp_bmat;
		comp_dbmats.Rows (2*l, 2*l+2) = dmat * comp_bmat;

		NgProfiler::AddFlops (timer2, nd*nd*4);
	      }
	    // comp_elmat = Trans (comp_bmats) * comp_dbmats;
	    LapackMultAtB (comp_bmats, comp_dbmats, comp_elmat);

	    elmat.Rows(facetdofs).Cols(facetdofs) += comp_elmat;
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
      // int nd = cfel.GetNDof();  

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
	    const MappedIntegrationPoint<D,D> & mip = mir_vol[l];
	    double lam = coef_lam->Evaluate(mip);
	    
	    Vec<D> gi = grad.Row(l);
	    Vec<D> hv1 = Trans (mip.GetJacobianInverse()) * gi;
	    Vec<D> hv2 = mip.GetJacobianInverse() * hv1;
	    gi = (lam * mip.GetJacobiDet() * ir_vol[l].Weight()) * hv2;

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
	// const NORMAL * normals = ElementTopology::GetNormals(eltype);
	FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);

	Mat<2> dmat;

	for (int k = 0; k < nfacet; k++)

	  {
	    HeapReset hr(lh);

	    fel_facet.SelectFacet (k);

	    Vec<D> normal_ref = normals(k);

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
		const MappedIntegrationPoint<D,D> & mip = mir[l];
		double lam = coef_lam->Evaluate(mip);
              
		Mat<D> inv_jac = mip.GetJacobianInverse();
		double det = mip.GetJacobiDet();

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
    }



    virtual void
    CalcFlux (const FiniteElement & fel,
	      const FiniteElement & felflux,
	      const ElementTransformation & eltrans,
	      const FlatVector<> & elu, 
	      FlatVector<> & elflux,
	      bool applyd,
	      LocalHeap & lh) const
    {
      static Timer timer1 ("hdg - calcflux");
      static Timer timer_solve ("hdg - calcflux solve");
      static Timer timer_el ("hdg - calcflux el");
      static Timer timer_facet ("hdg - calcflux facet");

      RegionTimer rt (timer1);


      const IntegrationRule & ir = 
	SelectIntegrationRule(fel.ElementType(), max(fel.Order(),felflux.Order())+felflux.Order());

      MappedIntegrationRule<D,D> mir (ir, eltrans, lh);

      const ScalarFiniteElement<D> & fel_l2 = 
	dynamic_cast<const ScalarFiniteElement<D>&> 
	(dynamic_cast<const CompoundFiniteElement&> (fel)[0]);

      const FacetVolumeFiniteElement<D> & fel_facet = 
        dynamic_cast<const FacetVolumeFiniteElement<D> &>
	(dynamic_cast<const CompoundFiniteElement&> (fel)[1]);


      const HDivFiniteElement<D> & hdfelflux = 
	dynamic_cast<const HDivFiniteElement<D>&> (felflux);


      int ndflux = felflux.GetNDof();
      int ndu = fel.GetNDof();
      int ndut = fel_l2.GetNDof();


      FlatMatrix<double> elmat(ndflux, lh);
      FlatMatrix<double> helmat(ndflux, lh);
      FlatVector<> elflux1 (ndflux, lh);
      FlatMatrixFixWidth<D> shape (ndflux, lh);
      FlatMatrixFixWidth<D> dshape (ndut, lh);
      FlatVector<> shapen (ndflux, lh);

      elmat = 0;
      elflux1 = 0;

      timer_el.Start();

      FlatMatrixFixWidth<D> gradm(ir.GetNIP(), lh);
      fel_l2.EvaluateGrad (ir, elu.Range(0, ndut), gradm);

      for (int i = 0; i < ir.GetNIP(); i++)
	{
	  Vec<D> grad = Trans (mir[i].GetJacobianInverse()) * gradm.Row(i);
	  grad *= coef_lam->Evaluate(mir[i]) * mir[i].GetWeight();
	  
	  hdfelflux.CalcMappedShape (mir[i], shape);
	  elmat += mir[i].GetWeight() * shape * Trans(shape);
	  elflux1 += shape * grad;
	}

      timer_el.Stop();

      timer_facet.Start();

      ELEMENT_TYPE eltype = fel_l2.ElementType();
      
      int nfacet = ElementTopology::GetNFacets(eltype);

      int sort[D+1];
      eltrans.GetSort (FlatArray<int> (D+1,&sort[0]));
      
      Facet2ElementTrafo transform(eltype);
      // const NORMAL * normals = ElementTopology::GetNormals(eltype);
      FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);
      Mat<2> dmat;
      
      for (int k = 0; k < nfacet; k++)
	{
	  HeapReset hr(lh);

	  fel_facet.SelectFacet (k);
	  
	  Vec<D> normal_ref = normals[k];
	  ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
	  
	  const IntegrationRule & ir_facet =
	    SelectIntegrationRule (etfacet, 2*fel_l2.Order());
	  
	  IntegrationRule ir_vol;
	  
	  for (int l = 0; l < ir_facet.GetNIP(); l++)
	    {
	      IntegrationPoint ip = transform(k, ir_facet[l]);
	      ip.SetWeight (ir_facet[l].Weight());
	      ir_vol.Append (ip);
	    }

	  MappedIntegrationRule<D,D> mir(ir_vol, eltrans, lh);


	  FlatVector<> shapes_l2(ir_vol.Size(), lh);
	  FlatVector<> shapes_facet(ir_vol.Size(), lh);
	  FlatMatrixFixWidth<D> grad_l2(ir_vol.Size(), lh);
	  

	  fel_l2.Evaluate (ir_vol, elu.Range(0, ndut), shapes_l2);
	  fel_l2.EvaluateGrad (ir_vol, elu.Range(0, ndut), grad_l2);
	  fel_facet.Evaluate (ir_vol, elu.Range(ndut, ndu), shapes_facet);
	  
	  helmat = 0.0;
	  for (int l = 0; l < ir_facet.GetNIP(); l++)
	    {
	      const MappedIntegrationPoint<D,D> & mip = mir[l];
	      double lam = coef_lam->Evaluate(mip);

              
	      Mat<D> inv_jac = mip.GetJacobianInverse();
	      double det = mip.GetJacobiDet();

	      
	      Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
	      double len = L2Norm (normal);
	      normal /= len;
	      
	      Vec<D> invjac_normal = inv_jac * normal;
		
	      dmat(0,0) = 0;
	      dmat(1,0) = dmat(0,1) = -1;
	      dmat(1,1) = alpha * ((fel_l2.Order()+1)*(fel_l2.Order()+D)/D * len) *(1.0/det);
	      // dmat *= lam * len * ir_vol[l].Weight();

	      Vec<2> hv1, hv2;
	      hv1(0) = InnerProduct (grad_l2.Row(l), invjac_normal);
	      hv1(1) = shapes_l2(l) - shapes_facet(l);
	      dmat *= lam;
	      
	      hv2 = dmat * hv1;
	      
	      hdfelflux.CalcMappedShape (mip, shape);
	      
	      shapen = shape * normal;

	      elmat += 1e6*mip.GetWeight() * shapen * Trans(shapen);
	      elflux1 -= 1e6*mip.GetWeight() * hv2(1) * shapen;
	    }
	}

      timer_facet.Stop();
      timer_solve.Start();
      FlatCholeskyFactors<double> invelmat(elmat, lh);
      invelmat.Mult (elflux1, elflux);
      timer_solve.Stop();
    }

  };





  template <int D>
  class HDG_LaplaceIntegrator2 : public BilinearFormIntegrator
  {
  protected:
    CoefficientFunction * coef_lam;
  public:
    HDG_LaplaceIntegrator2 (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      coef_lam  = coeffs[0];
    }

    virtual ~HDG_LaplaceIntegrator2 () { ; }

    virtual bool BoundaryForm () const { return 0; }



    virtual void CalcElementMatrix (const FiniteElement & fel,
				    const ElementTransformation & eltrans, 
				    FlatMatrix<double> & elmat,
				    LocalHeap & lh) const
    {
      double alpha = 0.1; // 0.01;

      static Timer timer ("HDG laplace");
      static Timer timer1 ("HDG laplace volume");
      static Timer timer1a ("HDG laplace volume, lapack");
      static Timer timer2 ("HDG laplace boundary");
      static Timer timer2a ("HDG laplace boundary - mult");
      static Timer timer3 ("HDG laplace inv/mult");

      RegionTimer reg (timer);

      const CompoundFiniteElement & cfel = 
        dynamic_cast<const CompoundFiniteElement&> (fel);
      const ScalarFiniteElement<D> & fel_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>&> (cfel[0]);
      const FacetVolumeFiniteElement<D> & fel_facet = 
        dynamic_cast<const FacetVolumeFiniteElement<D> &> (cfel[1]);
      
  
      ELEMENT_TYPE eltype = cfel.ElementType();

      IntRange l2_dofs = cfel.GetRange (0);
      IntRange facet_dofs = cfel.GetRange (1);
      
      int nd_l2 = fel_l2.GetNDof();
      int nd_facet = fel_facet.GetNDof();
      int nd = cfel.GetNDof();  

      // int base_l2 = 0;
      // int base_facet = base_l2+nd_l2;
      
      elmat = 0.0;

      FlatVector<> mat_l2(nd_l2, lh);
      FlatVector<> mat_dudn(nd_l2, lh);
      FlatVector<> mat_facet(nd_facet, lh);

      FlatMatrix<> mat_gradgrad (nd_l2, lh);
      FlatMatrix<> mat_gradgradinv (nd_l2, lh);
      FlatMatrix<> mat_robin (nd_l2, lh);

      {
	NgProfiler::RegionTimer reg (timer1);     
	HeapReset hr(lh);

	FlatMatrixFixWidth<D> dshape(nd_l2, lh);

	const IntegrationRule & ir_vol =
	  SelectIntegrationRule (eltype, 2*fel_l2.Order());
	MappedIntegrationRule<D,D> mir_vol(ir_vol, eltrans, lh);
	
	FlatMatrix<> bmats(ir_vol.GetNIP()*D, nd_l2, lh);
	FlatMatrix<> dbmats(ir_vol.GetNIP()*D, nd_l2, lh);

	for (int l = 0; l < ir_vol.GetNIP(); l++)
	  {
	    const MappedIntegrationPoint<D,D> & mip = mir_vol[l];
	    double lam = coef_lam->Evaluate(mip);
	    
	    fel_l2.CalcMappedDShape (mip, dshape);

	    bmats.Rows(l*D, (l+1)*D) = Trans(dshape);
	    dbmats.Rows(l*D, (l+1)*D) = (lam * mip.GetWeight()) * Trans(dshape);
	  }

	NgProfiler::RegionTimer reg1a (timer1a);     

	// mat_gradgrad = Trans (dbmats) * bmats;
        LapackMultAtB (dbmats, bmats, mat_gradgrad);
	elmat.Cols(l2_dofs).Rows(l2_dofs) = mat_gradgrad;
      }


      // The facet contribution
      {

	int nfacet = ElementTopology::GetNFacets(eltype);
      
	Facet2ElementTrafo transform(eltype); 
	// FlatVector<Vec<3> > normals (nfacet, ElementTopology::GetNormals(eltype));
	FlatVector<Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);

	FlatMatrix<> mat_mixed(nd_l2, nd, lh);
	FlatMatrix<> mat_mixedT(nd, nd_l2, lh);
	FlatVector<> jump(nd, lh);
	mat_mixed = 0.0;
	mat_robin = 0.0;
	
	{
	  NgProfiler::RegionTimer reg2 (timer2);     
	  for (int k = 0; k < nfacet; k++)
	    {
	      HeapReset hr(lh);
	      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);

	      const IntegrationRule & ir_facet = 
		SelectIntegrationRule (etfacet, 2*max (fel_l2.Order(), fel_facet.Order()));

	      fel_facet.SelectFacet (k);

	      FlatMatrix<> jumps(nd, ir_facet.GetNIP(), lh);
	      FlatMatrix<> fac_dudns(nd_l2, ir_facet.GetNIP(), lh);
	      FlatMatrix<> fac_jumps(nd, ir_facet.GetNIP(), lh);

	      for (int l = 0; l < ir_facet.GetNIP(); l++)
		{
		  IntegrationPoint ip = transform(k, ir_facet[l]);
		  MappedIntegrationPoint<D,D> mip (ip, eltrans);
		  double lam = coef_lam->Evaluate(mip);
              
		  Mat<D> inv_jac = mip.GetJacobianInverse();
		  double det = mip.GetMeasure();

		  Vec<D> normal = det * Trans (inv_jac) * normals(k);
		  double len = L2Norm (normal);
		  normal /= len;

		  fel_l2.CalcShape(ip, jump.Range (l2_dofs) );
		  fel_facet.CalcShape(ip, jump.Range (facet_dofs));
		  jump.Range (facet_dofs) *= -1;

		  Vec<D> invjac_normal = inv_jac * normal;
		  mat_dudn = fel_l2.GetDShape (mip.IP(), lh) * invjac_normal;

		  jumps.Col(l) = jump;
		  fac_jumps.Col(l) = alpha*lam * len/det * fel_l2.Order() * ir_facet[l].Weight() * jump;
		  fac_dudns.Col(l) = lam * len * ir_facet[l].Weight() * mat_dudn;
		}

	      NgProfiler::RegionTimer reg2a (timer2a);     
	      NgProfiler::AddFlops (timer2a, nd*nd*ir_facet.GetNIP());
	      NgProfiler::AddFlops (timer2a, nd_l2*nd*ir_facet.GetNIP());
	      NgProfiler::AddFlops (timer2a, nd_l2*nd_l2*ir_facet.GetNIP());
	      /*
		elmat += Symmetric (fac_jumps * Trans (jumps));
		mat_mixed += fac_dudns * Trans (jumps);
		mat_robin += Symmetric (fac_jumps.Rows(l2_dofs) * Trans (jumps.Rows(l2_dofs)));
	      */
	      LapackMultAddABt (fac_jumps, jumps, 1, elmat);
	      LapackMultAddABt (fac_dudns, jumps, 1, mat_mixed);
	      LapackMultAddABt (fac_jumps.Rows(l2_dofs), jumps.Rows(l2_dofs), 1, mat_robin);
	    }
	}
	
	elmat.Rows(l2_dofs) -= mat_mixed;
	elmat.Cols(l2_dofs) -= Trans(mat_mixed);


	double eps = 1e-6; // alpha;
	mat_gradgrad += eps * mat_robin;

	RegionTimer reg3 (timer3);

	mat_mixedT = Trans (mat_mixed);
	LapackInverse (mat_gradgrad);

	FlatMatrix<> hmT(nd, nd_l2, lh);
	/*
	  hmT = mat_mixedT * Trans (mat_gradgrad);
	  elmat += (1+alpha) * mat_mixedT * Trans (hmT);
	*/
	LapackMultABt (mat_mixedT, mat_gradgrad, hmT);
	LapackMultAddABt (mat_mixedT, hmT, 1+alpha, elmat);
      }
    }
  };









  


  template <int D>
  class HDG_LaplaceIntegrator3 : public BilinearFormIntegrator
  {
  protected:
    CoefficientFunction * coef_lam;
  public:
    HDG_LaplaceIntegrator3 (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      coef_lam  = coeffs[0];
    }

    virtual ~HDG_LaplaceIntegrator3 () { ; }

    virtual bool BoundaryForm () const { return 0; }



    virtual void CalcElementMatrix (const FiniteElement & fel,
				    const ElementTransformation & eltrans, 
				    FlatMatrix<double> & elmat,
				    LocalHeap & lh) const
    {
      // double alpha = 0.1; // 0.01;

      static Timer timer ("HDG laplace");
      static Timer timer1 ("HDG laplace volume");
      static Timer timer1a ("HDG laplace volume, lapack");
      static Timer timer2 ("HDG laplace boundary");
      static Timer timer2a ("HDG laplace boundary - mult");
      static Timer timer3 ("HDG laplace inv/mult");

      RegionTimer reg (timer);

      const CompoundFiniteElement & cfel = 
        dynamic_cast<const CompoundFiniteElement&> (fel);
      const ScalarFiniteElement<D> & fel_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>&> (cfel[0]);
      const FacetVolumeFiniteElement<D> & fel_facet = 
        dynamic_cast<const FacetVolumeFiniteElement<D> &> (cfel[1]);
      
  
      ELEMENT_TYPE eltype = cfel.ElementType();

      IntRange l2_dofs = cfel.GetRange (0);
      IntRange facet_dofs = cfel.GetRange (1);
      
      int nd_l2 = fel_l2.GetNDof();
      int nd_facet = fel_facet.GetNDof();
      int nd = cfel.GetNDof();  

      // int base_l2 = 0;
      // int base_facet = base_l2+nd_l2;
      
      elmat = 0.0;

      FlatVector<> mat_l2(nd_l2, lh);
      FlatVector<> mat_dudn(nd_l2, lh);
      FlatVector<> mat_facet(nd_facet, lh);

      FlatMatrix<> mat_gradgrad (nd_l2, lh);
      FlatMatrix<> mat_gradgradinv (nd_l2, lh);
      FlatMatrix<> mat_robin (nd_l2, lh);
      
      FlatMatrix<> mat_mass (nd_l2, lh);
      FlatMatrix<> mat_inv_mass (nd_l2, lh);
      mat_mass = 0.0;

      {
	NgProfiler::RegionTimer reg (timer1);     
	HeapReset hr(lh);

	FlatMatrixFixWidth<D> dshape(nd_l2, lh);

	const IntegrationRule & ir_vol =
	  SelectIntegrationRule (eltype, 2*fel_l2.Order());
	MappedIntegrationRule<D,D> mir_vol(ir_vol, eltrans, lh);
	
	FlatMatrix<> bmats(ir_vol.GetNIP()*D, nd_l2, lh);
	FlatMatrix<> dbmats(ir_vol.GetNIP()*D, nd_l2, lh);

	for (int l = 0; l < ir_vol.GetNIP(); l++)
	  {
	    const MappedIntegrationPoint<D,D> & mip = mir_vol[l];
	    double lam = coef_lam->Evaluate(mip);
	    
	    fel_l2.CalcMappedDShape (mip, dshape);

	    bmats.Rows(l*D, (l+1)*D) = Trans(dshape);
	    dbmats.Rows(l*D, (l+1)*D) = (lam * mip.GetWeight()) * Trans(dshape);

	    fel_l2.CalcShape (ir_vol[l], mat_l2);
	    mat_mass += (lam*mip.GetWeight()) * mat_l2 * Trans(mat_l2);
	  }

	NgProfiler::RegionTimer reg1a (timer1a);     

	// mat_gradgrad = Trans (dbmats) * bmats;
        LapackMultAtB (dbmats, bmats, mat_gradgrad);
	elmat.Cols(l2_dofs).Rows(l2_dofs) = mat_gradgrad;

	mat_inv_mass = mat_mass;
	LapackInverse (mat_inv_mass);
      }


      // The facet contribution
      {

	int nfacet = ElementTopology::GetNFacets(eltype);
      
	Facet2ElementTrafo transform(eltype); 
	// FlatVector<Vec<3> > normals (nfacet, ElementTopology::GetNormals(eltype));
	FlatVector<Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);

	FlatMatrix<> mat_mixed(nd_l2, nd, lh);
	FlatMatrix<> mat_mixed2(nd_l2, nd, lh);
	FlatMatrix<> mat_mixed2h(nd_l2, nd, lh);
	FlatMatrix<> mat_mixedT(nd, nd_l2, lh);
	FlatVector<> jump(nd, lh);
	mat_mixed = 0.0;
	mat_robin = 0.0;
	
	{
	  NgProfiler::RegionTimer reg2 (timer2);     
	  for (int k = 0; k < nfacet; k++)
	    {
	      HeapReset hr(lh);
	      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);

	      const IntegrationRule & ir_facet = 
		SelectIntegrationRule (etfacet, 2*max (fel_l2.Order(), fel_facet.Order()));

	      fel_facet.SelectFacet (k);

	      FlatMatrix<> jumps(nd, ir_facet.GetNIP(), lh);
	      FlatMatrix<> fac_dudns(nd_l2, ir_facet.GetNIP(), lh);
	      FlatMatrix<> fac_jumps(nd, ir_facet.GetNIP(), lh);

	      for (int l = 0; l < ir_facet.GetNIP(); l++)
		{
		  IntegrationPoint ip = transform(k, ir_facet[l]);
		  MappedIntegrationPoint<D,D> mip (ip, eltrans);
		  double lam = coef_lam->Evaluate(mip);
              
		  Mat<D> inv_jac = mip.GetJacobianInverse();
		  double det = mip.GetMeasure();

		  Vec<D> normal = det * Trans (inv_jac) * normals(k);
		  double len = L2Norm (normal);
		  normal /= len;

		  fel_l2.CalcShape(ip, jump.Range (l2_dofs) );
		  fel_facet.CalcShape(ip, jump.Range (facet_dofs));
		  jump.Range (facet_dofs) *= -1;

		  Vec<D> invjac_normal = inv_jac * normal;
		  mat_dudn = fel_l2.GetDShape (mip.IP(), lh) * invjac_normal;

		  jumps.Col(l) = jump;
		  fac_jumps.Col(l) = lam * len * ir_facet[l].Weight() * jump;
		  fac_dudns.Col(l) = lam * len * ir_facet[l].Weight() * mat_dudn;
		}

	      NgProfiler::RegionTimer reg2a (timer2a);     
	      NgProfiler::AddFlops (timer2a, nd*nd*ir_facet.GetNIP());
	      NgProfiler::AddFlops (timer2a, nd_l2*nd*ir_facet.GetNIP());
	      NgProfiler::AddFlops (timer2a, nd_l2*nd_l2*ir_facet.GetNIP());
	      /*
		elmat += Symmetric (fac_jumps * Trans (jumps));
		mat_mixed += fac_dudns * Trans (jumps);
		mat_robin += Symmetric (fac_jumps.Rows(l2_dofs) * Trans (jumps.Rows(l2_dofs)));
	      */
	      double h = pow (mat_mass(0,0), 1.0/D); 
	      LapackMultAddABt (fac_jumps, jumps, fel_l2.Order()/h, elmat);
	      LapackMultAddABt (fac_dudns, jumps, 1, mat_mixed);

	      LapackMultABt (fac_jumps.Rows(l2_dofs), jumps, mat_mixed2);
	      LapackMultAB (mat_inv_mass, mat_mixed2, mat_mixed2h);
	      LapackMultAddAtB (mat_mixed2, mat_mixed2h, 5, elmat);
	      // mat_mixed2 = fac_jumps.Rows(l2_dofs) * Trans (jumps);
	      // mat_mixed2h = mat_inv_mass * mat_mixed2;
	      // elmat += 2 * Trans (mat_mixed2) * mat_mixed2h;
	    }
	}
	
	elmat.Rows(l2_dofs) -= mat_mixed;
	elmat.Cols(l2_dofs) -= Trans(mat_mixed);

	/*
	double eps = 1e-6; // alpha;
	mat_gradgrad += eps * mat_robin;

	RegionTimer reg3 (timer3);

	mat_mixedT = Trans (mat_mixed);
	LapackInverse (mat_gradgrad);

	FlatMatrix<> hmT(nd, nd_l2, lh);
	LapackMultABt (mat_mixedT, mat_gradgrad, hmT);
	LapackMultAddABt (mat_mixedT, hmT, 1+alpha, elmat);
	*/
      }
    }
  };




  template <int D>
  class HDG_LaplaceIntegrator4 : public BilinearFormIntegrator
  {
  protected:
    CoefficientFunction * coef_lam;
  public:
    HDG_LaplaceIntegrator4 (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      coef_lam  = coeffs[0];
    }

    virtual ~HDG_LaplaceIntegrator4  () { ; }

    virtual bool BoundaryForm () const { return 0; }



    virtual void CalcElementMatrix (const FiniteElement & fel,
				    const ElementTransformation & eltrans, 
				    FlatMatrix<double> & elmat,
				    LocalHeap & lh) const
    {
      // double alpha = 0.1; // 0.01;

      static Timer timer ("HDG laplace");
      static Timer timer1 ("HDG laplace volume");
      static Timer timer1a ("HDG laplace volume, lapack");
      static Timer timer2 ("HDG laplace boundary");
      static Timer timer2a ("HDG laplace boundary - mult");
      static Timer timer3 ("HDG laplace inv/mult");

      RegionTimer reg (timer);

      const CompoundFiniteElement & cfel = 
        dynamic_cast<const CompoundFiniteElement&> (fel);
      const ScalarFiniteElement<D> & fel_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>&> (cfel[0]);
      const FacetVolumeFiniteElement<D> & fel_facet = 
        dynamic_cast<const FacetVolumeFiniteElement<D> &> (cfel[1]);
      
  
      ELEMENT_TYPE eltype = cfel.ElementType();

      IntRange l2_dofs = cfel.GetRange (0);
      IntRange facet_dofs = cfel.GetRange (1);
      
      int nd_l2 = fel_l2.GetNDof();
      int nd_facet = fel_facet.GetNDof();
      int nd = cfel.GetNDof();  

      // int base_l2 = 0;
      // int base_facet = base_l2+nd_l2;
      
      elmat = 0.0;

      FlatVector<> mat_l2(nd_l2, lh);
      FlatVector<> mat_dudn(nd_l2, lh);
      FlatVector<> mat_facet(nd_facet, lh);

      FlatMatrix<> mat_gradgrad (nd_l2, lh);
      FlatMatrix<> mat_gradgradinv (nd_l2, lh);
      FlatMatrix<> mat_robin (nd_l2, lh);
      
      FlatMatrix<> mat_mass (nd_l2, lh);
      FlatMatrix<> mat_inv_mass (nd_l2, lh);
      mat_mass = 0.0;

      {
	NgProfiler::RegionTimer reg (timer1);     
	HeapReset hr(lh);

	FlatMatrixFixWidth<D> dshape(nd_l2, lh);

	const IntegrationRule & ir_vol =
	  SelectIntegrationRule (eltype, 2*fel_l2.Order());
	MappedIntegrationRule<D,D> mir_vol(ir_vol, eltrans, lh);
	
	FlatMatrix<> bmats(ir_vol.GetNIP()*D, nd_l2, lh);
	FlatMatrix<> dbmats(ir_vol.GetNIP()*D, nd_l2, lh);

	for (int l = 0; l < ir_vol.GetNIP(); l++)
	  {
	    const MappedIntegrationPoint<D,D> & mip = mir_vol[l];
	    double lam = coef_lam->Evaluate(mip);
	    
	    fel_l2.CalcMappedDShape (mip, dshape);

	    bmats.Rows(l*D, (l+1)*D) = Trans(dshape);
	    dbmats.Rows(l*D, (l+1)*D) = (lam * mip.GetWeight()) * Trans(dshape);

	    fel_l2.CalcShape (ir_vol[l], mat_l2);
	    mat_mass += (lam*mip.GetWeight()) * mat_l2 * Trans(mat_l2);
	  }

	NgProfiler::RegionTimer reg1a (timer1a);     

	// mat_gradgrad = Trans (dbmats) * bmats;
        LapackMultAtB (dbmats, bmats, mat_gradgrad);
	elmat.Cols(l2_dofs).Rows(l2_dofs) = mat_gradgrad;

	mat_inv_mass = mat_mass;
	LapackInverse (mat_inv_mass);
      }


      // The facet contribution
      {

	int nfacet = ElementTopology::GetNFacets(eltype);
      
	Facet2ElementTrafo transform(eltype); 
	// FlatVector<Vec<3> > normals (nfacet, ElementTopology::GetNormals(eltype));
	FlatVector<Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);

	FlatMatrix<> mat_mixed(nd_l2, nd, lh);
	FlatMatrix<> mat_mixed2(nd_l2, nd, lh);
	FlatMatrix<> mat_mixed2h(nd_l2, nd, lh);
	FlatMatrix<> mat_mixedT(nd, nd_l2, lh);
	FlatVector<> jump(nd, lh);
	mat_mixed = 0.0;
	mat_robin = 0.0;
	
	{
	  NgProfiler::RegionTimer reg2 (timer2);     

	  for (int dir = 0; dir < D; dir++)
	    {
	      mat_mixed2 = 0.0;

	      for (int k = 0; k < nfacet; k++)
		{
		  HeapReset hr(lh);
		  ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
		  
		  const IntegrationRule & ir_facet = 
		    SelectIntegrationRule (etfacet, 2*max (fel_l2.Order(), fel_facet.Order()));
		  
		  fel_facet.SelectFacet (k);
		  
		  FlatMatrix<> jumps(nd, ir_facet.GetNIP(), lh);
		  FlatMatrix<> fac_dudns(nd_l2, ir_facet.GetNIP(), lh);
		  FlatMatrix<> fac_jumps(nd, ir_facet.GetNIP(), lh);
		  Vec<D> normal;

		  for (int l = 0; l < ir_facet.GetNIP(); l++)
		    {
		      IntegrationPoint ip = transform(k, ir_facet[l]);
		      MappedIntegrationPoint<D,D> mip (ip, eltrans);
		      double lam = coef_lam->Evaluate(mip);
		      
		      Mat<D> inv_jac = mip.GetJacobianInverse();
		      double det = mip.GetMeasure();
		      
		      normal = det * Trans (inv_jac) * normals(k);
		      double len = L2Norm (normal);
		      normal /= len;
		      
		      fel_l2.CalcShape(ip, jump.Range (l2_dofs) );
		      fel_facet.CalcShape(ip, jump.Range (facet_dofs)); 
		      jump.Range (facet_dofs) *= -1;
		      
		      Vec<D> invjac_normal = inv_jac * normal;
		      mat_dudn = fel_l2.GetDShape (mip.IP(), lh) * invjac_normal;
		      
		      jumps.Col(l) = jump;
		      fac_jumps.Col(l) = lam * len * ir_facet[l].Weight() * jump;
		      fac_dudns.Col(l) = lam * len * ir_facet[l].Weight() * mat_dudn;
		    }
		  
		  /*
		    elmat += Symmetric (fac_jumps * Trans (jumps));
		    mat_mixed += fac_dudns * Trans (jumps);
		    mat_robin += Symmetric (fac_jumps.Rows(l2_dofs) * Trans (jumps.Rows(l2_dofs)));
		  */
		  double h = pow (mat_mass(0,0), 1.0/D); 
		  if (dir == 0)
		    {
		      LapackMultAddABt (fac_jumps, jumps, 0.1*fel_l2.Order()/h, elmat);
		      LapackMultAddABt (fac_dudns, jumps, 1, mat_mixed);
		    }
		  LapackMultAddABt (fac_jumps.Rows(l2_dofs), jumps, normal(dir), mat_mixed2);
		}
	      
	      LapackMultAB (mat_inv_mass, mat_mixed2, mat_mixed2h);
	      LapackMultAddAtB (mat_mixed2, mat_mixed2h, 1.5, elmat);
	    }
	}
	
	elmat.Rows(l2_dofs) -= mat_mixed;
	elmat.Cols(l2_dofs) -= Trans(mat_mixed);
	  
	/*
	double eps = 1e-6; // alpha;
	mat_gradgrad += eps * mat_robin;

	RegionTimer reg3 (timer3);

	mat_mixedT = Trans (mat_mixed);
	LapackInverse (mat_gradgrad);

	FlatMatrix<> hmT(nd, nd_l2, lh);
	LapackMultABt (mat_mixedT, mat_gradgrad, hmT);
	LapackMultAddABt (mat_mixedT, hmT, 1+alpha, elmat);
	*/
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

      IntRange l2_dofs = cfel.GetRange (0);
      IntRange facet_dofs = cfel.GetRange (1);
      
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

          
      IntegrationRule ir_vol(eltype, 2*fel_l2.Order());
      
      mat_l2 = 0.0;

      for (int l = 0; l < ir_vol.GetNIP(); l++)
        {
          HeapReset hr(lh);
          const MappedIntegrationPoint<D,D> mip(ir_vol[l], eltrans);
          Vec<D> conv;

          if (coef_conv.Size()>1)
	    for (int j = 0; j < D; j++)
	      conv(j) = coef_conv[j]->Evaluate(mip);
	  else
	    coef_conv[0]->Evaluate(mip, conv);

          fel_l2.CalcShape (mip.IP(), shape);
          fel_l2.CalcMappedDShape (mip, dshape);
          
          conv_dshape = dshape * conv;

          conv_dshape *= mip.GetJacobiDet() * ir_vol[l].Weight();
          
          mat_l2 -= conv_dshape * Trans (shape);
        }
      
      elmat.Cols(l2_dofs).Rows(l2_dofs) = mat_l2; 


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
              MappedIntegrationPoint<D,D> mip (ip, eltrans);
	      
              Vec<D> conv;
	      if (coef_conv.Size()>1)
		for (int j = 0; j < D; j++)
		  conv(j) = coef_conv[j]->Evaluate(mip);
	      else
		coef_conv[0]->Evaluate(mip,conv);


              // Mat<D> jac = mip.GetJacobian();
              Mat<D> inv_jac = mip.GetJacobianInverse();
              double det = mip.GetJacobiDet();

              Vec<D> normal = det * Trans (inv_jac) * normal_ref; 
      
              double len = L2Norm (normal);
              normal /= len;

              double bn = InnerProduct (conv, normal);
              bool inflow = (bn < 0);

              fel_facet.CalcFacetShape(k, ir_facet[l], shape_facet);
              fel_l2.CalcShape(mip.IP(), shape);
              
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

	  elmat.Rows(facetdofs).Cols(facetdofs) += comp_elmat;
        }
    }
  };



  static RegisterBilinearFormIntegrator<HDG_LaplaceIntegrator<2> > initlap2 ("HDG_laplace", 2, 2);
  static RegisterBilinearFormIntegrator<HDG_LaplaceIntegrator<3> > initlap3 ("HDG_laplace", 3, 2);

  static RegisterBilinearFormIntegrator<HDG_LaplaceIntegrator3<2> > initlap2a ("HDG_laplace2", 2, 1);
  static RegisterBilinearFormIntegrator<HDG_LaplaceIntegrator3<3> > initlap3a ("HDG_laplace2", 3, 1);

  static RegisterBilinearFormIntegrator<HDG_LaplaceIntegrator4<2> > initlap24 ("HDG_laplace4", 2, 1);
  static RegisterBilinearFormIntegrator<HDG_LaplaceIntegrator4<3> > initlap34 ("HDG_laplace4", 3, 1);


  static RegisterBilinearFormIntegrator<HDG_ConvectionIntegrator<2> > initconv21 ("HDG_convection", 2, 1);
  static RegisterBilinearFormIntegrator<HDG_ConvectionIntegrator<2> > initconv22 ("HDG_convection", 2, 2);

  static RegisterBilinearFormIntegrator<HDG_ConvectionIntegrator<3> > initconv31 ("HDG_convection", 3, 1);
  static RegisterBilinearFormIntegrator<HDG_ConvectionIntegrator<3> > initconv33 ("HDG_convection", 3, 3);
}

