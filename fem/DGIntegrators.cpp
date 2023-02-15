/*********************************************************************/
/* File:   DGIntegrators.cpp                                         */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   17. Nov. 2009                                             */
/* see also: http://sourceforge.net/apps/mediawiki/ngsolve/index.php?title=Discontinuous_Galerkin_Methods_from_the_PDE_file */
/*********************************************************************/


#include <fem.hpp>

namespace ngfem
{
  /** 
      DG for scalar Laplace
      //TODO: Convection-Integrators: volume term integrator in B1DB2-Form ? etc..)
      //TODO: find out why bilinearform -symmetric is not working correctly
      //TODO: Gridfunctioncoefficientfunction for conv-bilinearform
  */
  
  using namespace ngfem;

  namespace DG_FORMULATIONS{
    enum DGTYPE{
      IP,
      NIPG,
      BO
    };
  }
  
  template <int D, DG_FORMULATIONS::DGTYPE dgtype>
  class DGInnerFacet_LaplaceIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    double alpha;   // interior penalyty
    shared_ptr<CoefficientFunction> coef_lam;
  public:
    DGInnerFacet_LaplaceIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
    // : FacetBilinearFormIntegrator(coeffs)
    { 
      coef_lam  = coeffs[0];
      if (dgtype!=DG_FORMULATIONS::BO)
	alpha = coeffs[1] -> EvaluateConst();
    }

    virtual ~DGInnerFacet_LaplaceIntegrator () { ; }

    virtual VorB VB() const
    { return VOL; }
    virtual xbool IsSymmetric () const { return true; }
    
    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new DGInnerFacet_LaplaceIntegrator (coeffs);
    }
    
    virtual void CalcElementMatrix (const FiniteElement & fel,
                                        const ElementTransformation & eltrans, 
                                        FlatMatrix<double> elmat,
                                        LocalHeap & lh) const
    {
      throw Exception("DGInnerFacet_LaplaceIntegrator::CalcElementMatrix - not implemented!");
    }
    
    virtual void CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                                  const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                                  const FiniteElement & volumefel2, int LocalFacetNr2,
                                  const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                                  FlatMatrix<double> elmat,
                                  LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("DGInnerFacet_LaplaceIntegrator");

      if (LocalFacetNr2==-1) throw Exception("DGFacetLaplaceIntegrator: LocalFacetNr2==1");

      NgProfiler::RegionTimer reg (timer);
      const ScalarFiniteElement<D> * fel1_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel1);
      ELEMENT_TYPE eltype1 = volumefel1.ElementType();
      int nd1 = fel1_l2->GetNDof();

      const ScalarFiniteElement<D> * fel2_l2 = NULL;
      ELEMENT_TYPE eltype2 = eltype1;

      fel2_l2 = dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel2);
      eltype2 = volumefel2.ElementType();
      int nd2 = fel2_l2->GetNDof();
      int maxorder = max2(fel1_l2->Order(),fel2_l2->Order());
      
      elmat = 0.0;

      FlatVector<> mat1_shape(nd1, lh);
      FlatVector<> mat1_dudn(nd1, lh);
      FlatVector<> mat2_shape(nd2, lh);
      FlatVector<> mat2_dudn(nd2, lh);
      
      FlatMatrixFixHeight<2> bmat(nd1+nd2, lh);
      FlatMatrixFixHeight<2> dbmat(nd1+nd2, lh);
      Mat<2> dmat;

      FlatMatrixFixWidth<D> dshape(nd1, lh);
      FlatMatrixFixWidth<D> fac_dshape(nd1, lh);
      
      Facet2ElementTrafo transform1(eltype1,ElVertices1); 
      Facet2ElementTrafo transform2(eltype2,ElVertices2); 

      const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);
      const NORMAL * normals2 = ElementTopology::GetNormals(eltype2);

      HeapReset hr(lh);
      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

      Vec<D> normal_ref1, normal_ref2;
      for (int i=0; i<D; i++){
	normal_ref1(i) = normals1[LocalFacetNr1][i];
	normal_ref2(i) = normals2[LocalFacetNr2][i];
      }
      const IntegrationRule & ir_facet =
	SelectIntegrationRule (etfacet, 2*maxorder);
	
      if (maxorder==0) maxorder=1;
   
      bmat = 0.0;
      for (int l = 0; l < ir_facet.GetNIP(); l++)
	{
	  IntegrationPoint ip1 = transform1(LocalFacetNr1, ir_facet[l]);
	  
	  MappedIntegrationPoint<D,D> sip1 (ip1, eltrans1);
	  double lam = coef_lam->Evaluate(sip1);

	  // Mat<D> jac1 = sip1.GetJacobian();
	  Mat<D> inv_jac1 = sip1.GetJacobianInverse();
	  double det1 = sip1.GetJacobiDet();

	  Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;       
	  double len1 = L2Norm (normal1);
	  normal1 /= len1;

	  fel1_l2->CalcShape(sip1.IP(), mat1_shape);
	  Vec<D> invjac_normal1 = inv_jac1 * normal1;
	  mat1_dudn = fel1_l2->GetDShape (sip1.IP(), lh) * invjac_normal1;
	  
	  IntegrationPoint ip2 = (LocalFacetNr2!=-1) ? transform2(LocalFacetNr2, ir_facet[l]) : ip1;
	  MappedIntegrationPoint<D,D> sip2 (ip2, eltrans2);
	  // double lam2 = coef_lam->Evaluate(sip2);
	  // Mat<D> jac2 = sip2.GetJacobian();
	  Mat<D> inv_jac2 = sip2.GetJacobianInverse();
	  double det2 = sip2.GetJacobiDet();
	  
	  Vec<D> normal2 = det2 * Trans (inv_jac2) * normal_ref2;       
	  double len2 = L2Norm (normal2); 
	  if(abs(len1-len2)>1e-6){
	    std::cout << "len :\t" << len1 << "\t=?=\t" << len2 << std::endl;
	    throw Exception ("DGInnerFacet_LaplaceIntegrator: len1!=len2");
	  }
	  normal2 /= len2;
	  Vec<D> invjac_normal2;;
	  fel2_l2->CalcShape(sip2.IP(), mat2_shape);
	  invjac_normal2 = inv_jac2 * normal2;
	  mat2_dudn = fel2_l2->GetDShape (sip2.IP(), lh) * invjac_normal2;
	  
	  bmat.Row(0).Range (0   , nd1)   = 0.5 * mat1_dudn;	    
	  bmat.Row(0).Range (nd1   , nd1+nd2)   = -0.5 * mat2_dudn;
	  bmat.Row(1).Range (0   , nd1)   = mat1_shape;
	  bmat.Row(1).Range (nd1   , nd1+nd2)   = -mat2_shape;

	  dmat(0,0) = 0;
	  dmat(1,0) = -1;
	  switch (dgtype){
	    case DG_FORMULATIONS::BO:
	      dmat(0,1) = 1;
	      dmat(1,1) = 0;
	      break;
	    case DG_FORMULATIONS::NIPG:
	      dmat(0,1) = 1; 
// 	      dmat(1,1) = alpha * sqr (maxorder) * (len1/det1);
	      dmat(1,1) = alpha * ((maxorder+1.0)*(maxorder+D)/D * len1) *(1.0/det1);
	      break;
	    case DG_FORMULATIONS::IP:
	    default:
	      dmat(0,1) = -1; 
// 	      dmat(1,1) = alpha * sqr (maxorder) * (len1/det1);
	      dmat(1,1) = alpha * ((maxorder+1.0)*(maxorder+D)/D * len1) *(1.0/det1);
	      break;	      
	  }
	  dmat *= lam * len1 * ir_facet[l].Weight();
	  dbmat = dmat * bmat;
	  elmat += Trans (bmat) * dbmat;
	}
	if (LocalFacetNr2==-1) elmat=0.0;
      }
  };



  template <int D>
  class ConvectionIntegrator : public BilinearFormIntegrator
  {
  protected:
    // shared_ptr<CoefficientFunction>  coef_conv[D];
    DVec<D> coef_conv;
  public:
    ConvectionIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : BilinearFormIntegrator(), coef_conv (coeffs)
    { 
      /*
      for (int j = 0; j < D; j++)
        coef_conv[j] = coeffs[j];
      */
    }

    virtual ~ConvectionIntegrator () { ; }

    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new ConvectionIntegrator (coeffs);
    }

    virtual VorB VB() const
    { return VOL; }
    virtual xbool IsSymmetric () const { return false; }

    virtual void CalcElementMatrix (const FiniteElement & fel,
                                        const ElementTransformation & eltrans, 
                                        FlatMatrix<double> elmat,
                                        LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("ConvectionIntegrator");

      NgProfiler::RegionTimer reg (timer);

      const ScalarFiniteElement<D> & fel_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>&> (fel);
      ELEMENT_TYPE eltype = fel.ElementType();
      
      int nd = fel_l2.GetNDof();

      FlatVector<> shape(nd, lh);
      FlatMatrixFixWidth<D> dshape(nd, lh);
      FlatVector<> conv_dshape(nd, lh);
      
      FlatMatrixFixHeight<2> bmat(nd, lh);
      FlatMatrixFixHeight<2> dbmat(nd, lh);
      // Mat<2> dmat;
          
      const IntegrationRule & ir_vol =
        SelectIntegrationRule (eltype, 2*fel_l2.Order());
      
      elmat = 0.0;

      for (int l = 0; l < ir_vol.GetNIP(); l++)
      {
	HeapReset hr(lh);
	const MappedIntegrationPoint<D,D> sip(ir_vol[l], eltrans);
	Vec<D> conv;
        coef_conv.GenerateVector (fel_l2, sip, conv, lh);
        /*
	for (int j = 0; j < D; j++)
	  conv(j) = coef_conv[j]->Evaluate(sip);
        */

	fel_l2.CalcShape (sip.IP(), shape);
	fel_l2.CalcMappedDShape (sip, dshape);
	
	conv_dshape = dshape * conv;

	conv_dshape *= sip.GetJacobiDet() * ir_vol[l].Weight();
	
        elmat  -= conv_dshape * Trans (shape);
      }
    }//end of CalcElementMatrix
    
  };//end of class 

  
  
  
  template <int D>
  class DGInnerFacet_ConvectionIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    Array<shared_ptr<CoefficientFunction> > coef_b;
  public:
    DGInnerFacet_ConvectionIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
    // : FacetBilinearFormIntegrator(coeffs)
    { 
      coef_b.SetSize(D);
      for (int i=0; i<D; i++)
	coef_b[i]  = coeffs[i];
    }

    virtual ~DGInnerFacet_ConvectionIntegrator () { ; }
    
    virtual VorB VB() const
    { return VOL; }
    virtual xbool IsSymmetric () const { return false; }
    
    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new DGInnerFacet_ConvectionIntegrator (coeffs);
    }
    
    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> elmat,
                                    LocalHeap & lh) const
    {
      throw Exception("DGInnerFacet_ConvectionIntegrator::CalcElementMatrix - not implemented!");
    }
    
    virtual void CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                         FlatMatrix<double> elmat,
                         LocalHeap & lh) const 
    {
      static int timer = NgProfiler::CreateTimer ("DGInnerFacet_ConvectionIntegrator");
      if (LocalFacetNr2==-1) throw Exception ("DGInnerFacet_ConvectionIntegrator: LocalFacetNr2==-1");
      NgProfiler::RegionTimer reg (timer);
      const ScalarFiniteElement<D> * fel1_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel1);
      ELEMENT_TYPE eltype1 = volumefel1.ElementType();
      int nd1 = fel1_l2->GetNDof();

      const ScalarFiniteElement<D> * fel2_l2 =  
	dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel2);
      ELEMENT_TYPE eltype2 = volumefel2.ElementType();
      eltype2 = volumefel2.ElementType();
      int nd2 = fel2_l2->GetNDof();
      int maxorder = max2(fel1_l2->Order(),fel2_l2->Order());
      
      elmat = 0.0;

      FlatVector<> mat1_shape(nd1, lh);
      FlatVector<> mat2_shape(nd2, lh);
      
      FlatVector<> b1mat(nd1+nd2, lh); //B(v)
      FlatVector<> b2mat(nd1+nd2, lh); //B(u)
//       double dmat; //bn
      FlatMatrixFixWidth<D> dshape(nd1, lh);
      FlatMatrixFixWidth<D> fac_dshape(nd1, lh);
      Facet2ElementTrafo transform1(eltype1,ElVertices1); 
      const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);
      Facet2ElementTrafo transform2(eltype2,ElVertices2); 
      const NORMAL * normals2 = ElementTopology::GetNormals(eltype2);

      HeapReset hr(lh);
      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

      Vec<D> normal_ref1, normal_ref2;
      for (int i=0; i<D; i++){
	normal_ref1(i) = normals1[LocalFacetNr1][i];
	normal_ref2(i) = normals2[LocalFacetNr2][i];
      }
      const IntegrationRule & ir_facet =
	SelectIntegrationRule (etfacet, 2*maxorder);
      if (maxorder==0) maxorder=1;
   
      b1mat = 0.0;
      b2mat = 0.0;
      for (int l = 0; l < ir_facet.GetNIP(); l++)
	{
	  IntegrationPoint ip1 = transform1(LocalFacetNr1, ir_facet[l]);
	  
	  MappedIntegrationPoint<D,D> sip1 (ip1, eltrans1);

	  // Mat<D> jac1 = sip1.GetJacobian();
	  Mat<D> inv_jac1 = sip1.GetJacobianInverse();
	  double det1 = sip1.GetJacobiDet();

	  Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;       
	  double len1 = L2Norm (normal1);
	  normal1 /= len1;

	  IntegrationPoint ip2 = transform2(LocalFacetNr2, ir_facet[l]);
	  MappedIntegrationPoint<D,D> sip2 (ip2, eltrans2);

	  // Mat<D> jac2 = sip2.GetJacobian();
	  Mat<D> inv_jac2 = sip2.GetJacobianInverse();
	  double det2 = sip2.GetJacobiDet();
	  
	  Vec<D> normal2 = det2 * Trans (inv_jac2) * normal_ref2;       
	  double len2 = L2Norm (normal2); 
	  
	  if(abs(len1-len2)>1e-6){
	    std::cout << "len :\t" << len1 << "\t=?=\t" << len2 << std::endl;
	    throw Exception("DGInnerFacet_ConvectionIntegrator: len1 != len2");
	  }
	  
	  normal2 /= len2;
	  fel1_l2->CalcShape(sip1.IP(), mat1_shape);
	  fel2_l2->CalcShape(sip2.IP(), mat2_shape);
	  
	  double bn = 0;//, bn1 = 0, bn2 = 0;
	  for (int i=0; i<D;i++){
	    bn += 0.5 * (coef_b[i]->Evaluate(sip1) + coef_b[i]->Evaluate(sip2) )  * normal1(i);
	  }
	  b1mat = 0.0; 
	  b1mat.Range (0   , nd1)   = mat1_shape;
	  b1mat.Range (nd1   , nd1+nd2)   = -mat2_shape; //sign for outer normal (bn2 = - bn1)
	  b2mat = 0.0;
	  if(bn>0) //flow from 1 to 2 => 1 is upwind direction
	    b2mat.Range (0   , nd1) = mat1_shape;
	  else
	    b2mat.Range (nd1   , nd1+nd2) = mat2_shape;
	  // dmat = bn;
	  b1mat *= bn * len1 * ir_facet[l].Weight();
	  elmat += b1mat * Trans (b2mat);
	}
      }
  };

// TODO: This doesn't have to be a FacetLinearFormIntegrator if the normal is available anyway!
  template <int D>
  class DGBoundaryFacet_ConvectionIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    double alpha;   // interior penalyty
    Array<shared_ptr<CoefficientFunction> > coef_b;
  public:
    DGBoundaryFacet_ConvectionIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
    // : FacetBilinearFormIntegrator(coeffs)
    { 
      coef_b.SetSize(D);
      for (int i=0; i<D; i++)
	coef_b[i]  = coeffs[i];
    }

    virtual ~DGBoundaryFacet_ConvectionIntegrator () { ; }

    virtual VorB VB() const
    { return BND; }
    virtual xbool IsSymmetric () const { return false; }
  
    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new DGBoundaryFacet_ConvectionIntegrator (coeffs);
    }
    
    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> elmat,
                                    LocalHeap & lh) const
    {
      throw Exception("DGBoundaryFacet_ConvectionIntegrator::CalcElementMatrix - not implemented!");
    }
    
    virtual void CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                                  const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                                  const ElementTransformation & seltrans, FlatArray<int> & SElVertices,
                                  FlatMatrix<double> elmat,
                                  LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("DGBoundaryFacet_ConvectionIntegrator");
      NgProfiler::RegionTimer reg (timer);
      const ScalarFiniteElement<D> & fel = 
        dynamic_cast<const ScalarFiniteElement<D>&> (volumefel);
      ELEMENT_TYPE eltype = volumefel.ElementType();
      int nd = fel.GetNDof();

      elmat = 0.0;

      FlatVector<> mat_shape(nd, lh);
      
      FlatVector<> b1mat(nd, lh); //B(v)
      FlatVector<> b2mat(nd, lh); //B(u)
      // double dmat; //bn
      Facet2ElementTrafo transform(eltype,ElVertices); 
      const NORMAL * normals = ElementTopology::GetNormals(eltype);

      HeapReset hr(lh);
      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, LocalFacetNr);

      Vec<D> normal_ref;
      for (int i=0; i<D; i++){
	normal_ref(i) = normals[LocalFacetNr][i];
      }
      const IntegrationRule & ir_facet =
	SelectIntegrationRule (etfacet, 2*fel.Order());

   

      for (int l = 0; l < ir_facet.GetNIP(); l++)
	{
	  IntegrationPoint ip = transform(LocalFacetNr, ir_facet[l]);
	  
	  MappedIntegrationPoint<D,D> sip (ip, eltrans);

	  // Mat<D> jac = sip.GetJacobian();
	  Mat<D> inv_jac = sip.GetJacobianInverse();
	  double det = sip.GetJacobiDet();

	  Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
	  double len = L2Norm (normal);
	  normal /= len;
	  
	  double bn = 0;
	  for (int i=0; i<D;i++){
	    bn += coef_b[i]->Evaluate(sip) * normal(i);
	  }
	  if (bn<0) continue; //on this integration point, there is only inflow
	  fel.CalcShape(sip.IP(), mat_shape);
	  b1mat = b2mat = 0.0;	    
	  b1mat = mat_shape;
	  b2mat = 0.0;
	  b2mat = mat_shape;
	  b1mat *= bn * len * ir_facet[l].Weight();
	  elmat += b1mat * Trans (b2mat);
	}
      }
  };
  
  
  
  template <int D, DG_FORMULATIONS::DGTYPE dgtype>
  class DGBoundaryFacet_LaplaceIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    double alpha;   // interior penalyty
    shared_ptr<CoefficientFunction> coef_lam;
  public:
    DGBoundaryFacet_LaplaceIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
    // : FacetBilinearFormIntegrator(coeffs)
    { 
      coef_lam  = coeffs[0];
      if (dgtype!=DG_FORMULATIONS::BO)
	alpha = coeffs[1] -> EvaluateConst();
    }

    virtual VorB VB() const
    { return BND; }
    virtual xbool IsSymmetric () const { return true; }

    virtual ~DGBoundaryFacet_LaplaceIntegrator () { ; }

    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new DGBoundaryFacet_LaplaceIntegrator (coeffs);
    }
    
    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> elmat,
                                    LocalHeap & lh) const
    {
      throw Exception("DGBoundaryFacet_LaplaceIntegrator::CalcElementMatrix - not implemented!");
    }
    
    virtual void CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                                  const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                                  const ElementTransformation & seltrans, FlatArray<int> & SElVertices,
                                  FlatMatrix<double> elmat,
                                  LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("DGBoundaryFacet_LaplaceIntegrator boundary");
      NgProfiler::RegionTimer reg (timer);
      const ScalarFiniteElement<D> * fel1_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel);
      ELEMENT_TYPE eltype1 = volumefel.ElementType();
      int nd1 = fel1_l2->GetNDof();

      // const ScalarFiniteElement<D> * fel2_l2 = NULL;
      // ELEMENT_TYPE eltype2 = eltype1;
      int nd2 = 0;

      int maxorder = fel1_l2->Order();

      
      elmat = 0.0;

      FlatVector<> mat1_shape(nd1, lh);
      FlatVector<> mat1_dudn(nd1, lh);
      
      FlatMatrixFixHeight<2> bmat(nd1+nd2, lh);
      FlatMatrixFixHeight<2> dbmat(nd1+nd2, lh);
      Mat<2> dmat;

      FlatMatrixFixWidth<D> dshape(nd1, lh);
      FlatMatrixFixWidth<D> fac_dshape(nd1, lh);
      Facet2ElementTrafo transform1(eltype1,ElVertices); 
      const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);
      
      HeapReset hr(lh);
      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

      Vec<D> normal_ref1, normal_ref2;
      for (int i=0; i<D; i++){
	normal_ref1(i) = normals1[LocalFacetNr][i];
      }
      const IntegrationRule & ir_facet =
	SelectIntegrationRule (etfacet, 2*maxorder);
      if (maxorder==0) maxorder=1;

      bmat = 0.0;
      for (int l = 0; l < ir_facet.GetNIP(); l++)
	{
	  IntegrationPoint ip1 = transform1(LocalFacetNr, ir_facet[l]);
	  
	  MappedIntegrationPoint<D,D> sip1 (ip1, eltrans);
	  double lam = coef_lam->Evaluate(sip1);

	  MappedIntegrationPoint<D-1,D> sips (ir_facet[l], seltrans);
	  
	  // Mat<D> jac1 = sip1.GetJacobian();
	  Mat<D> inv_jac1 = sip1.GetJacobianInverse();
	  double det1 = sip1.GetJacobiDet();

	  Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;       
	  double len1 = L2Norm (normal1);
	  normal1 /= len1;

	  fel1_l2->CalcShape(sip1.IP(), mat1_shape);
	  Vec<D> invjac_normal1 = inv_jac1 * normal1;
	  mat1_dudn = fel1_l2->GetDShape (sip1.IP(), lh) * invjac_normal1;
	  
	  bmat.Row(0).Range (0   , nd1)   = mat1_dudn;
	  bmat.Row(1).Range (0   , nd1)   = mat1_shape;
	  
	  dmat(0,0) = 0;
	  dmat(1,0) = -1;
	  
	  switch (dgtype){
	    case DG_FORMULATIONS::BO:
	      dmat(0,1) = 1;
	      dmat(1,1) = 0;
	      break;
	    case DG_FORMULATIONS::NIPG:
	      dmat(0,1) = 1; 
// 	      dmat(1,1) = alpha * sqr (maxorder) * (len1/det1);
	      dmat(1,1) = alpha * ((maxorder+1.0)*(maxorder+D)/D * len1) *(1.0/det1);
	      break;
	    case DG_FORMULATIONS::IP:
	    default:
	      dmat(0,1) = -1; 
// 	      dmat(1,1) = alpha * sqr (maxorder) * (len1/det1);
	      dmat(1,1) = alpha * ((maxorder+1.0)*(maxorder+D)/D * len1) *(1.0/det1);
	      break;	      
	  }
	  dmat *= lam * len1 * ir_facet[l].Weight();
	  dbmat = dmat * bmat;
	  elmat += Trans (bmat) * dbmat;
	}
      }
  };
  
  
  template <int D, DG_FORMULATIONS::DGTYPE dgtype>
  class DGFacet_DirichletBoundaryIntegrator : public FacetLinearFormIntegrator
  {
  protected:
    double alpha;   // interior penalyty
    shared_ptr<CoefficientFunction> coef_lam;
    shared_ptr<CoefficientFunction> coef_dir;
  public:
    DGFacet_DirichletBoundaryIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : FacetLinearFormIntegrator( /* coeffs */)
    { 
      coef_lam  = coeffs[0];
      coef_dir  = coeffs[1];
      alpha = coeffs[2] -> EvaluateConst();
    }


    virtual VorB VB() const
    { return BND; }
    
    virtual ~DGFacet_DirichletBoundaryIntegrator () { ; }

    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new DGFacet_DirichletBoundaryIntegrator (coeffs);
    }
    
    virtual void CalcElementVector (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatVector<double> elvec,
                                    LocalHeap & lh) const
    {
      throw Exception("DGFacet_DirichletBoundaryIntegrator::CalcElementVector - not implemented!");
    }
    
    virtual void CalcFacetVector (const FiniteElement & volumefel, int LocalFacetNr,
			 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			 const ElementTransformation & seltrans,
                         FlatVector<double> elvec, LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("DGFacet_DirichletBoundaryIntegrator");
      NgProfiler::RegionTimer reg (timer);
      const ScalarFiniteElement<D> * fel1_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel);
      ELEMENT_TYPE eltype1 = volumefel.ElementType();
      int nd1 = fel1_l2->GetNDof();
      elvec = 0.0;

      FlatVector<> mat1_shape(nd1, lh);
      FlatVector<> mat1_dudn(nd1, lh);
      
      Facet2ElementTrafo transform1(eltype1);//,ElVertices); at domain boundaries don't change orientation: orientation should still coincide with the orientation of the surface element transformation

      const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);

      HeapReset hr(lh);
      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

      Vec<D> normal_ref1;
      Vec<D> dvec(2);
      dvec=0.0;
      for (int i=0; i<D; i++){
	normal_ref1(i) = normals1[LocalFacetNr][i];
      }
      int maxorder = fel1_l2->Order();
      const IntegrationRule & ir_facet =
	SelectIntegrationRule (etfacet, 2*maxorder);
      if (maxorder==0) maxorder=1;

      for (int l = 0; l < ir_facet.GetNIP(); l++)
	{
	  IntegrationPoint ip1 = transform1(LocalFacetNr, ir_facet[l]);
	  
	  MappedIntegrationPoint<D,D> sip1 (ip1, eltrans);
	  double lam = coef_lam->Evaluate(sip1);

	  MappedIntegrationPoint<D-1,D> sips (ir_facet[l], seltrans);
	  double dir = coef_dir->Evaluate(sips);

	  // Mat<D> jac1 = sip1.GetJacobian();
	  Mat<D> inv_jac1 = sip1.GetJacobianInverse();
	  double det1 = sip1.GetJacobiDet();

	  Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;       
	  double len1 = L2Norm (normal1);
	  normal1 /= len1;

	  fel1_l2->CalcShape(sip1.IP(), mat1_shape);

	  Vec<D> invjac_normal1 = inv_jac1 * normal1;

	  mat1_dudn = fel1_l2->GetDShape (sip1.IP(), lh) * invjac_normal1;
	  
	  
	  switch (dgtype){
	    case DG_FORMULATIONS::BO:
	      dvec(0) = 1; 
	      dvec(1) = 0;
	      break;
	    case DG_FORMULATIONS::NIPG:
	      dvec(0) = 1; 
// 	      dvec(1) = alpha * sqr (maxorder) * (len1/det1);
	      dvec(1) = alpha * ((maxorder+1.0)*(maxorder+D)/D * len1) *(1.0/det1);	      
	      break;
	    case DG_FORMULATIONS::IP:
	    default:
	      dvec(0) = -1;
// 	      dvec(1) = alpha * sqr (maxorder) * (len1/det1);
	      dvec(1) = alpha * ((maxorder+1.0)*(maxorder+D)/D * len1) *(1.0/det1);	      
	      break;	      
	  }	  
	  
	  double fac = len1*ir_facet[l].Weight()*lam*dir;
	  elvec += (fac * dvec(0))*mat1_dudn;
	  elvec += (fac * dvec(1) ) * mat1_shape;
	}
      }
  };

 
    template <int D>
  class DGFacet_NeumannBoundaryIntegrator : public FacetLinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> coef_lam;
  public:
    DGFacet_NeumannBoundaryIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : FacetLinearFormIntegrator( /* coeffs */ )
    { 
      coef_lam  = coeffs[0];
    }

    virtual VorB VB() const
    { return BND; }

    virtual ~DGFacet_NeumannBoundaryIntegrator () { ; }

    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new DGFacet_NeumannBoundaryIntegrator (coeffs);
    }
    
    virtual void CalcElementVector (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatVector<double> elvec,
                                    LocalHeap & lh) const
    {
      throw Exception("DGFacet_NeumannBoundaryIntegrator::CalcElementVector - not implemented!");
    }
    
    virtual void CalcFacetVector (const FiniteElement & volumefel, int LocalFacetNr,
			 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			 const ElementTransformation & seltrans,
                         FlatVector<double> elvec, LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("DGFacet_NeumannBoundaryIntegrator");

      NgProfiler::RegionTimer reg (timer);
      const ScalarFiniteElement<D> * fel1_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel);
      ELEMENT_TYPE eltype1 = volumefel.ElementType();
      int nd1 = fel1_l2->GetNDof();
      elvec = 0.0;

      FlatVector<> mat1_shape(nd1, lh);
      FlatVector<> mat1_dudn(nd1, lh);
      
      Facet2ElementTrafo transform1(eltype1);//,ElVertices); at domain boundaries don't change orientation: orientation should still coincide with the orientation of the surface element transformation

      const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);

      HeapReset hr(lh);
      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

      Vec<D> normal_ref1;
      for (int i=0; i<D; i++){
	normal_ref1(i) = normals1[LocalFacetNr][i];
      }
      int maxorder = fel1_l2->Order();
      const IntegrationRule & ir_facet =
	SelectIntegrationRule (etfacet, 2*maxorder);
      if (maxorder==0) maxorder=1;

      for (int l = 0; l < ir_facet.GetNIP(); l++)
	{
	  IntegrationPoint ip1 = transform1(LocalFacetNr, ir_facet[l]);
	  
	  MappedIntegrationPoint<D,D> sip1 (ip1, eltrans);
	  double lam = coef_lam->Evaluate(sip1);

	  MappedIntegrationPoint<D-1,D> sips (ir_facet[l], seltrans);
	  
	  // Mat<D> jac1 = sip1.GetJacobian();
	  Mat<D> inv_jac1 = sip1.GetJacobianInverse();
	  double det1 = sip1.GetJacobiDet();

	  Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;       
	  double len1 = L2Norm (normal1);
	  normal1 /= len1;

	  fel1_l2->CalcShape(sip1.IP(), mat1_shape);

	  // Vec<D> invjac_normal1 = inv_jac1 * normal1;

	  double fac = len1*ir_facet[l].Weight()*lam;
	  elvec += fac *  mat1_shape;
	}
      }
  };

// TODO: This doesn't have to be a FacetLinearFormIntegrator if the normal is available anyway!
  template <int D>
  class DGFacet_ConvectionDirichletBoundaryIntegrator : public FacetLinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> coef_rob;
    Array<shared_ptr<CoefficientFunction> > coef_b;
  public:
    DGFacet_ConvectionDirichletBoundaryIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : FacetLinearFormIntegrator( /* coeffs */)
    { 
      coef_b.SetSize(D);
      for (int i=0; i<D; i++)
	coef_b[i]  = coeffs[i];
      coef_rob = coeffs[D];
    }

    virtual VorB VB() const
    { return BND; }
    
    virtual ~DGFacet_ConvectionDirichletBoundaryIntegrator () { ; }

    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new DGFacet_ConvectionDirichletBoundaryIntegrator (coeffs);
    }
    
    virtual void CalcElementVector (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatVector<double> elvec,
                                    LocalHeap & lh) const
    {
      throw Exception("DGFacet_ConvectionDirichletBoundaryIntegrator::CalcElementVector - not implemented!");
    }
    
    virtual void CalcFacetVector (const FiniteElement & volumefel, int LocalFacetNr,
			 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			 const ElementTransformation & seltrans,
                         FlatVector<double> elvec, LocalHeap & lh) const
    {
       static int timer = NgProfiler::CreateTimer ("DGFacet_ConvectionDirichletBoundaryIntegrator");
      NgProfiler::RegionTimer reg (timer);
      const ScalarFiniteElement<D> & fel = 
        dynamic_cast<const ScalarFiniteElement<D>&> (volumefel);
      ELEMENT_TYPE eltype = volumefel.ElementType();
      int nd = fel.GetNDof();

      elvec = 0.0;

      FlatVector<> mat_shape(nd, lh);
      
      FlatVector<> b1mat(nd, lh); //B(v)
      FlatVector<> b2mat(nd, lh); //B(u)
      // double dmat; //bn
//       Facet2ElementTrafo transform(eltype,ElVertices); 
      Facet2ElementTrafo transform(eltype);//,ElVertices); at domain boundaries don't change orientation: orientation should still coincide with the orientation of the surface element transformation

      const NORMAL * normals = ElementTopology::GetNormals(eltype);

      HeapReset hr(lh);
      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, LocalFacetNr);

      Vec<D> normal_ref;
      for (int i=0; i<D; i++){
	normal_ref(i) = normals[LocalFacetNr][i];
      }
      const IntegrationRule & ir_facet =
	SelectIntegrationRule (etfacet, 2*fel.Order());

   

      for (int l = 0; l < ir_facet.GetNIP(); l++)
	{
	  IntegrationPoint ip = transform(LocalFacetNr, ir_facet[l]);
	  
	  MappedIntegrationPoint<D,D> sip (ip, eltrans);
	  MappedIntegrationPoint<D-1,D> sipsurf (ir_facet[l], seltrans);

	  // Mat<D> jac = sip.GetJacobian();
	  Mat<D> inv_jac = sip.GetJacobianInverse();
	  double det = sip.GetJacobiDet();

	  Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
	  double len = L2Norm (normal);
	  normal /= len;
	  
	  double bn = 0;
	  for (int i=0; i<D;i++){
	    bn += coef_b[i]->Evaluate(sip) * normal(i);
	  }
	  double val = coef_rob->Evaluate(sipsurf);
	  if (bn>0) continue; //on this integration point, there is only inflow
	  fel.CalcShape(sip.IP(), mat_shape);
	  
	  b1mat = mat_shape;
	  
	  elvec -= val * bn * len * ir_facet[l].Weight() * b1mat;
	}
      }
  };


  template class DGFacet_NeumannBoundaryIntegrator<2>;
  template class DGFacet_NeumannBoundaryIntegrator<3>;

//convection:  
  static RegisterBilinearFormIntegrator<ConvectionIntegrator<2> > initconv_bf_inner_vol_2d ("convection", 2, 2);
  static RegisterBilinearFormIntegrator<ConvectionIntegrator<3> > initconv_bf_inner_vol_3d ("convection", 3, 3);
  static RegisterBilinearFormIntegrator<DGInnerFacet_ConvectionIntegrator<2> > initconv_bf_inner_fac_2d ("DG_innfac_convection", 2, 2);
  static RegisterBilinearFormIntegrator<DGInnerFacet_ConvectionIntegrator<3> > initconv_bf_inner_fac_3d ("DG_innfac_convection", 3, 3);
  static RegisterBilinearFormIntegrator<DGBoundaryFacet_ConvectionIntegrator<2> > initconv_bf_bound_fac_2d ("DG_bndfac_convection", 2, 2);
  static RegisterBilinearFormIntegrator<DGBoundaryFacet_ConvectionIntegrator<3> > initconv_bf_bound_fac_3d ("DG_bndfac_convection", 3, 3);
  static RegisterLinearFormIntegrator<DGFacet_ConvectionDirichletBoundaryIntegrator<2> > initconv_lf_bound_fac_2d ("DG_bndfac_convdir", 2, 3);
  static RegisterLinearFormIntegrator<DGFacet_ConvectionDirichletBoundaryIntegrator<3> > initconv_lf_bound_fac_3d ("DG_bndfac_convdir", 3, 4);

//laplace
  static RegisterBilinearFormIntegrator<DGInnerFacet_LaplaceIntegrator<2,DG_FORMULATIONS::IP> > initlap_bf_inner_fac_2d ("DGIP_innfac_laplace", 2, 2);
  static RegisterBilinearFormIntegrator<DGInnerFacet_LaplaceIntegrator<3,DG_FORMULATIONS::IP> > initlap_bf_inner_fac_3d ("DGIP_innfac_laplace", 3, 2);
  static RegisterBilinearFormIntegrator<DGBoundaryFacet_LaplaceIntegrator<2,DG_FORMULATIONS::IP> > initlap_bf_outer_fac_2d ("DGIP_bndfac_laplace", 2, 2);
  static RegisterBilinearFormIntegrator<DGBoundaryFacet_LaplaceIntegrator<3,DG_FORMULATIONS::IP> > initlap_bf_outer_fac_3d ("DGIP_bndfac_laplace", 3, 2);
  static RegisterLinearFormIntegrator<DGFacet_DirichletBoundaryIntegrator<2,DG_FORMULATIONS::IP> > initlap_lf_outer_fac_2d ("DGIP_bndfac_dir", 2, 3);
  static RegisterLinearFormIntegrator<DGFacet_DirichletBoundaryIntegrator<3,DG_FORMULATIONS::IP> > initlap_lf_outer_fac_3d ("DGIP_bndfac_dir", 3, 3);
  static RegisterLinearFormIntegrator<DGFacet_NeumannBoundaryIntegrator<2> > initneu_lf_outer_fac_2d ("DGIP_bndfac_neumann", 2, 2);
  static RegisterLinearFormIntegrator<DGFacet_NeumannBoundaryIntegrator<3> > initneu_lf_outer_fac_3d ("DGIP_bndfac_neumann", 3, 2);

// as nitsche  
//those integrators may also be useful/used without DG being involved and thereby the designation "nitsche" is more appropriate here:
  static RegisterBilinearFormIntegrator<DGBoundaryFacet_LaplaceIntegrator<2,DG_FORMULATIONS::IP> > initlap_bf_nitsche_2d ("nitsche", 2, 2);
  static RegisterBilinearFormIntegrator<DGBoundaryFacet_LaplaceIntegrator<3,DG_FORMULATIONS::IP> > initlap_bf_nitsche_3d ("nitsche", 3, 2);
  static RegisterLinearFormIntegrator<DGFacet_DirichletBoundaryIntegrator<2,DG_FORMULATIONS::IP> > initlap_lf_nitsche_2d ("nitsche", 2, 3);
  static RegisterLinearFormIntegrator<DGFacet_DirichletBoundaryIntegrator<3,DG_FORMULATIONS::IP> > initlap_lf_nitsche_3d ("nitsche", 3, 3);
  
/*    // nonsymmteric Interior Penalty method: consistent and stable
  static RegisterBilinearFormIntegrator<DGInnerFacet_LaplaceIntegrator<2,DG_FORMULATIONS::NIPG> > initlap_bf_inner_fac_2d_nipg ("DGNIPG_innerfac_laplace", 2, 2);
  static RegisterBilinearFormIntegrator<DGInnerFacet_LaplaceIntegrator<3,DG_FORMULATIONS::NIPG> > initlap_bf_inner_fac_3d_nipg ("DGNIPG_innerfac_laplace", 3, 2);
  static RegisterBilinearFormIntegrator<DGBoundaryFacet_LaplaceIntegrator<2,DG_FORMULATIONS::NIPG> > initlap_bf_outer_fac_2d_nipg ("DGNIPG_boundfac_laplace", 2, 2);
  static RegisterBilinearFormIntegrator<DGBoundaryFacet_LaplaceIntegrator<3,DG_FORMULATIONS::NIPG> > initlap_bf_outer_fac_3d_nipg ("DGNIPG_boundfac_laplace", 3, 2);
  static RegisterBilinearFormIntegrator<DGFacet_DirichletBoundaryIntegrator<2,DG_FORMULATIONS::NIPG> > initlap_lf_outer_fac_2d_nipg ("DGNIPG_boundfac_dir", 2, 3);
  static RegisterBilinearFormIntegrator<DGFacet_DirichletBoundaryIntegrator<3,DG_FORMULATIONS::NIPG> > initlap_lf_outer_fac_3d_nipg ("DGNIPG_boundfac_dir", 3, 3);

    // Baumann-Oden method: (only) consistent
  static RegisterBilinearFormIntegrator<DGInnerFacet_LaplaceIntegrator<2,DG_FORMULATIONS::BO> > initlap_bf_inner_fac_2d_bo ("DGBO_innerfac_laplace", 2, 2);
  static RegisterBilinearFormIntegrator<DGInnerFacet_LaplaceIntegrator<3,DG_FORMULATIONS::BO> > initlap_bf_inner_fac_3d_bo ("DGBO_innerfac_laplace", 3, 2);
  static RegisterBilinearFormIntegrator<DGBoundaryFacet_LaplaceIntegrator<2,DG_FORMULATIONS::BO> > initlap_bf_outer_fac_2d_bo ("DGBO_boundfac_laplace", 2, 1);
  static RegisterBilinearFormIntegrator<DGBoundaryFacet_LaplaceIntegrator<3,DG_FORMULATIONS::BO> > initlap_bf_outer_fac_3d_bo ("DGBO_boundfac_laplace", 3, 1);
  static RegisterBilinearFormIntegrator<DGFacet_DirichletBoundaryIntegrator<2,DG_FORMULATIONS::BO> > initlap_lf_outer_fac_2d_bo ("DGBO_boundfac_dir", 2, 3);
  static RegisterBilinearFormIntegrator<DGFacet_DirichletBoundaryIntegrator<3,DG_FORMULATIONS::BO> > initlap_lf_outer_fac_3d_bo ("DGBO_boundfac_dir", 3, 3);
*/					 

}; //end of namespace ngfem
