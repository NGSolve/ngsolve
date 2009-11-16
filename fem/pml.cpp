/*********************************************************************/
/* File:   pml.cpp                                                   */
/* Author: Joachim Schoeberl                                         */
/* Date:   Nov 3, 2003                                               */
/*********************************************************************/

/*
  Billinear form integrator over complex domain
*/


// **************************************************************** 
// Parameters are set in SetPMLParameters() (via access on pde), which is 
// called by the pml bdb integrators. 
// 
// 
// Define constants in pde:
// pml_alpha: factor for damping   
// For type of PML_layer (circular,square,rectangular) 
// pml_r = : PML-layer around circle with radius pml_r and center 0 
// pml_x = : PML-layer around square [-pml_x,pml_x]^d
// pml_xmin=, pml_ymin=, pml_zmin=, pml_xmax=, ... 
//     PML_layer around [pml_xmin,pml_xmax]x[pml_ymin,pml_ymax]x[pml_zmin,pml_zmax]
//
// *****************************************************************   


#include <fem.hpp>
#include "pml.hpp"


extern ngstd::SymbolTable<double> & GetConstantTable ();



namespace ngfem
{

  Complex alpha(0,1);
  double pml_r = 1;
  double pml_x = 1;

  double pml_xmin[3] = { -1, -1, -1};
  double pml_xmax[3] = {  1,  1,  1};

  Vec<3> pml_center (0, 0, 0);

  /*
    rect_pml = 0 .... circular pml with radius pml_r and center pml_center
    rect_pml = 1 .... square pml on square (-pml_x, pml_x)^d
    rect_pml = 2 .... rectangular pml on (pml_xmin, pml_xmax) x (pml_ymin, pml_ymax) x (pml_zmin, pml_zmax)
   */
  int rect_pml = 0;
  


  bool apply_deriv_alpha = 0;



  
  template <>
  SpecificIntegrationPoint<2,2,Complex> :: 
  SpecificIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,
			    LocalHeap & lh)
    : DimSpecificIntegrationPoint<2,Complex> (aip, aeltrans)
  {
    Mat<2,2> hdxdxi;
    Vec<2> hpoint;

    eltrans.CalcPointJacobian (ip, hpoint, hdxdxi, lh);

    
    // rect_pml = 3;

    switch  (rect_pml)
      {
      case 0:  // circle
	{
	  double abs_x = L2Norm (hpoint);
	  if (abs_x <= pml_r)  //  || eltrans.GetElementIndex() == 0)
	    {
	      point = hpoint;
	      dxdxi = hdxdxi;
	    }
	  else
	    {
	      Complex g = 1.+alpha*(1.0-pml_r/abs_x);
	      point = g * hpoint;
	      // SZ: sollte da nicht abs_x * abs_x anstelle  abs_x*abs_x * abs_x stehen? 
              // JS: das hat schon so gestimmt
	      Mat<2,2,Complex> trans =
		g * Id<2>() + (pml_r*alpha/(abs_x*abs_x*abs_x)) * (hpoint * Trans(hpoint));
	      dxdxi = trans * hdxdxi;
	    }
	  break;
	}
      case 1:
	{
	  Vec<2> dabs_dpoint;
	  double abs_x;
	  if (fabs (hpoint(0)) > fabs(hpoint(1)))
	    {
	      abs_x = fabs (hpoint(0));
	      dabs_dpoint(0) = (hpoint(0) > 0) ? 1 : -1;
	      dabs_dpoint(1) = 0;
	    }
	  else
	    {
	      abs_x = fabs (hpoint(1));
	      dabs_dpoint(0) = 0;
	      dabs_dpoint(1) = (hpoint(1) > 0) ? 1 : -1;
	    }
	  
	  if (abs_x <= pml_x)   //  || hpoint(1) < 0)
	    {
	      point = hpoint;
	      dxdxi = hdxdxi;
	    }
	  else
	    {
	      Complex g = 1.+alpha*(1.0-pml_x/abs_x);
	      point = g * hpoint;
	      Mat<2,2,Complex> trans =
		g * Id<2>() + (pml_x*alpha/(abs_x*abs_x)) * (hpoint * Trans(dabs_dpoint));
	      dxdxi = trans * hdxdxi;
	    }
	  break;
	}

      case 2:
	{
	  point = hpoint;
	  dxdxi = hdxdxi;

	  for (int i = 0; i < 2; i++)
	    {
	      if (hpoint(i) > pml_xmax[i])
		{
		  point(i) += alpha * (point(i) - pml_xmax[i]);

		  Mat<2,2,Complex> trans = Id<2>();
		  trans(i,i) += alpha;
		  Mat<2,2,Complex> hm = dxdxi;
		  dxdxi = trans * hm;
		}

	      else if (hpoint(i) < pml_xmin[i])
		{
		  point(i) -= alpha * (pml_xmin[i] - point(i));
		  
		  Mat<2,2,Complex> trans = Id<2>();
		  trans(i,i) += alpha;
		  Mat<2,2,Complex> hm = dxdxi;
		  dxdxi = trans * hm;
		}
	    }
	  break;
	}


      case 3:
	{
	  // double r2 = 1.2;
	  double abs_x = L2Norm (hpoint);
	  if (abs_x <= pml_r)
	    {
	      point = hpoint;
	      dxdxi = hdxdxi;
	    }
	  else
	    {
	      Complex g = 1.+alpha*(1.0-pml_r/abs_x);
	      point = g * hpoint;
	      Mat<2,2,Complex> trans =
		g * Id<2>() + (pml_r*alpha/(abs_x*abs_x*abs_x)) * (hpoint * Trans(hpoint));
	      dxdxi = trans * hdxdxi;
	    }
	  break;
	}


      default:
        break;
      }
	
    det = Det (dxdxi);
    dxidx = Inv (dxdxi);
  }

  
  template <>
  SpecificIntegrationPoint<3,3,Complex> :: 
  SpecificIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,
			    LocalHeap & lh)
    : DimSpecificIntegrationPoint<3,Complex> (aip, aeltrans)
  {
    Mat<3,3> hdxdxi;
    Vec<3> hpoint, hvec;

    eltrans.CalcPointJacobian (ip, hpoint, hdxdxi, lh);
    // eltrans.CalcJacobian (ip, hdxdxi, lh);
    // eltrans.CalcPoint (ip, hpoint, lh);

    switch (rect_pml)
      {
      case 0:
	{
	  hvec = hpoint - pml_center;
	  double abs_x = L2Norm (hvec);
	  if (abs_x <= pml_r)
	    {
	      point = hpoint;
	      dxdxi = hdxdxi;
	    }
	  else
	    {
	      Complex g = 1.+alpha*(1-pml_r/abs_x);
	      point = pml_center + g * hvec;
	      Mat<3,3,Complex> trans =
		g * Id<3>() + (pml_r*alpha/(abs_x*abs_x*abs_x)) * (hvec * Trans(hvec));
	      dxdxi = trans * hdxdxi;
	      
	      // dxdxi = g * Id<3>() + (alpha/(abs_x*abs_x*abs_x)) * (hpoint * Trans(hpoint));
	    }
	  break;
	}
      case 1:
	{
	  double abs_x = fabs (hpoint(0));
	  if (abs_x <= 3.0)
	    {
	      point = hpoint;
	      dxdxi = hdxdxi;
	    }
	  else
	    {
	      Complex g = 1.+ alpha * (1.0-3.0/abs_x);  
	      point = g * hpoint;
	    dxdxi = g * Id<3>();
	    dxdxi(0,0) += 3.0*alpha/(abs_x*abs_x*abs_x) * hpoint(0);
	    }
	  break;
	}
      case 2:
	{
	  point = hpoint;
	  dxdxi = hdxdxi;

	  for (int i = 0; i < 3; i++)
	    {
	      if (hpoint(i) > pml_xmax[i])
		{
		  point(i) += alpha * (point(i) - pml_xmax[i]);

		  Mat<3,3,Complex> trans = Id<3>();
		  trans(i,i) += alpha;
		  Mat<3,3,Complex> hm = dxdxi;
		  dxdxi = trans * hm;
		}

	      else if (hpoint(i) < pml_xmin[i])
		{
		  point(i) -= alpha * (pml_xmin[i] - point(i));
		  
		  Mat<3,3,Complex> trans = Id<3>();
		  trans(i,i) += alpha;
		  Mat<3,3,Complex> hm = dxdxi;
		  dxdxi = trans * hm;
		}
	    }
	  break;
	}

      default: 
	break;
      }
    
    det = Det (dxdxi);
    dxidx = Inv (dxdxi);
  }


  // SZ moved to intrule.cpp 
  template class SpecificIntegrationPoint<2,2,Complex>;
  template class SpecificIntegrationPoint<3,3,Complex>;







  
  template <>
  SpecificIntegrationPoint<2,2,AutoDiff<1,Complex> > :: 
  SpecificIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,
			    LocalHeap & lh)
    : DimSpecificIntegrationPoint<2,AutoDiff<1,Complex> > (aip, aeltrans)
  {
    Mat<2,2> hdxdxi;
    Vec<2> hpoint;

    eltrans.CalcPointJacobian (ip, hpoint, hdxdxi, lh);



    AutoDiff<1,Complex> ad_alpha(alpha, 0);

    switch  (rect_pml)
      {
      case 0:  // circle
	{
	  double abs_x = L2Norm (hpoint);
	  if (abs_x <= pml_r)
	    {
	      point = hpoint;
	      dxdxi = hdxdxi;
	    }
	  else
	    {
	      AutoDiff<1,Complex> g = 1.+ad_alpha*(1.0-pml_r/abs_x);
	      point = g * hpoint;
	      Mat<2,2,AutoDiff<1,Complex> > trans =
		g * Id<2>() + (pml_r*ad_alpha/(abs_x*abs_x*abs_x)) * (hpoint * Trans(hpoint));
	      dxdxi = trans * hdxdxi;
	    }
	  break;
	}
      case 1:
	{
	  Vec<2> dabs_dpoint;
	  double abs_x;
	  if (fabs (hpoint(0)) > fabs(hpoint(1)))
	    {
	      abs_x = fabs (hpoint(0));
	      dabs_dpoint(0) = (hpoint(0) > 0) ? 1 : -1;
	      dabs_dpoint(1) = 0;
	    }
	  else
	    {
	      abs_x = fabs (hpoint(1));
	      dabs_dpoint(0) = 0;
	      dabs_dpoint(1) = (hpoint(1) > 0) ? 1 : -1;
	    }
	  
	  if (abs_x <= pml_x)   //  || hpoint(1) < 0)
	    {
	      point = hpoint;
	      dxdxi = hdxdxi;
	    }
	  else
	    {
	      AutoDiff<1,Complex> g = 1.+ad_alpha*(1.0-pml_x/abs_x);
	      point = g * hpoint;
	      Mat<2,2,AutoDiff<1,Complex> > trans =
		g * Id<2>() + (pml_r*ad_alpha/(abs_x*abs_x)) * (hpoint * Trans(dabs_dpoint));
	      dxdxi = trans * hdxdxi;
	    }
	  break;
	}

      case 2:
	{
	  point = hpoint;
	  dxdxi = hdxdxi;

	  for (int i = 0; i < 2; i++)
	    {
	      if (hpoint(i) > pml_xmax[i])
		{
		  
		  point(i) += ad_alpha * (point(i) - pml_xmax[i]);

		  Mat<2,2,AutoDiff<1,Complex> > trans = Id<2>();
		  trans(i,i) += ad_alpha;
		  Mat<2,2,AutoDiff<1,Complex> > hm = dxdxi;
		  dxdxi = trans * hm;
		}

	      else if (hpoint(i) < pml_xmin[i])
		{
		  point(i) -= ad_alpha * (pml_xmin[i] - point(i));
		  
		  
		  Mat<2,2,AutoDiff<1,Complex> > trans = Id<2>();
		  trans(i,i) += ad_alpha;
		  Mat<2,2,AutoDiff<1,Complex> > hm = dxdxi;
		  dxdxi = trans * hm;
		}
	    }
	  break;
	}




      default:
        break;
      }
	
    det = Det (dxdxi);
    dxidx = Inv (dxdxi);
  }



  template <>
  SpecificIntegrationPoint<3,3,AutoDiff<1,Complex> > :: 
  SpecificIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,
			    LocalHeap & lh)
    : DimSpecificIntegrationPoint<3,AutoDiff<1,Complex> > (aip, aeltrans)
  {
    cout << "AD not implemented for 3D" << endl;
  };


  template class SpecificIntegrationPoint<2,2,AutoDiff<1,Complex> >;
  template class SpecificIntegrationPoint<3,3,AutoDiff<1,Complex> >;




  void SetPMLParameters()
  {
    if (GetConstantTable().Used ("pml_r"))
      pml_r = GetConstantTable()["pml_r"];

    if (GetConstantTable().Used ("pml_cx"))
      pml_center(0) = GetConstantTable()["pml_cx"];
    if (GetConstantTable().Used ("pml_cy"))
      pml_center(1) = GetConstantTable()["pml_cy"];
    if (GetConstantTable().Used ("pml_cz"))
      pml_center(2) = GetConstantTable()["pml_cz"];

    
    if (GetConstantTable().Used ("pml_x"))
      {
	pml_x = GetConstantTable()["pml_x"];
	rect_pml = 1;
      }


    if (GetConstantTable().Used ("pml_xmin"))
      {
	pml_xmin[0] = GetConstantTable()["pml_xmin"];
	rect_pml = 2;
      }
    if (GetConstantTable().Used ("pml_xmax"))
      {
	pml_xmax[0] = GetConstantTable()["pml_xmax"];
	rect_pml = 2;
      }

    if (GetConstantTable().Used ("pml_ymin"))
      {
	pml_xmin[1] = GetConstantTable()["pml_ymin"];
	rect_pml = 2;
      }
    if (GetConstantTable().Used ("pml_ymax"))
      {
	pml_xmax[1] = GetConstantTable()["pml_ymax"];
	rect_pml = 2;
      }

    if (GetConstantTable().Used ("pml_zmin"))
      {
	pml_xmin[2] = GetConstantTable()["pml_zmin"];
	rect_pml = 2;
      }
    if (GetConstantTable().Used ("pml_zmax"))
      {
	pml_xmax[2] = GetConstantTable()["pml_zmax"];
	rect_pml = 2;
      }


    switch (rect_pml)
      {
      case 0: cout << "circular pml of radius " << pml_r << endl; break;
      case 1: cout << "square pml on +/- " << pml_x << endl; break;
      case 2: cout << "rectangular pml on " 
		   << "(" << pml_xmin[0] << "," << pml_xmax[0] 
		   << ") x (" << pml_xmin[1] << "," << pml_xmax[1] 
		   << ") x (" << pml_xmin[2] << "," << pml_xmax[2] << ")"
		   << endl; 
	break;
      }

    if (GetConstantTable().Used ("pml_alpha"))
      alpha = Complex (0, GetConstantTable()["pml_alpha"]);
  }






/// 
template <int D, typename FEL = ScalarFiniteElement<D> >
class PML_LaplaceIntegrator 
  : public PML_BDBIntegrator<DiffOpGradient<D>, DiagDMat<D>, FEL>
{
public:
  ///
  PML_LaplaceIntegrator (CoefficientFunction * coeff)
    : PML_BDBIntegrator<DiffOpGradient<D>, DiagDMat<D>, FEL> (DiagDMat<D> (coeff))
  {
    ;
  }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new PML_LaplaceIntegrator (coeffs[0]);
  }

  ///
  virtual string Name () const { return "PML_Laplace"; }
};




/// 
template <int D, typename FEL = ScalarFiniteElement<D> >
class PML_ElasticityIntegrator 
  : public PML_BDBIntegrator<DiffOpStrain<D>, ElasticityDMat<D>, FEL>
{
public:
  ///
  PML_ElasticityIntegrator (CoefficientFunction * coefe,
                            CoefficientFunction * coefnu)
    : PML_BDBIntegrator<DiffOpStrain<D>, ElasticityDMat<D>, FEL> (ElasticityDMat<D> (coefe, coefnu))
  {
    ;
  }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new PML_ElasticityIntegrator (coeffs[0], coeffs[1]);
  }

  ///
  virtual string Name () const { return "PML_Elasticity"; }
};





/// 
template <int D, typename FEL = ScalarFiniteElement<D> >
class PML_MassIntegrator 
  : public PML_BDBIntegrator<DiffOpId<D>, DiagDMat<1>, FEL>
{
public:
  ///
  PML_MassIntegrator (CoefficientFunction * coeff)
    : PML_BDBIntegrator<DiffOpId<D>, DiagDMat<1>, FEL> (DiagDMat<1> (coeff))
  {
    ;
  }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new PML_MassIntegrator (coeffs[0]);
  }

  ///
  virtual string Name () const { return "PML_Mass"; }
};



/*
///
class PML_ElasticityIntegrator
  : public PML_BDBIntegrator<DiffOpStrain<2>, PlaneStressDMat, ScalarFiniteElement>
{
public:
  ///
  PML_ElasticityIntegrator (CoefficientFunction * coefe,
			    CoefficientFunction * coefnu)
    : PML_BDBIntegrator<DiffOpStrain<2>, PlaneStressDMat, ScalarFiniteElement> 
  (PlaneStressDMat (coefe, coefnu))
  { ; }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new PML_ElasticityIntegrator (coeffs[0], coeffs[1]);
  }



  ///
  virtual string Name () const { return "Elasticity"; }
};
*/



template <int D>
class PML_CurlCurlEdgeIntegrator
  : public PML_BDBIntegrator<DiffOpCurlEdge<D>, DiagDMat<DIM_CURL_TRAIT<D>::DIM>, HCurlFiniteElement<D> >
{
public:
  ///
  PML_CurlCurlEdgeIntegrator (CoefficientFunction * coef)
    : PML_BDBIntegrator<DiffOpCurlEdge<D>, DiagDMat<DIM_CURL_TRAIT<D>::DIM>, HCurlFiniteElement<D> > 
  (DiagDMat<DIM_CURL_TRAIT<D>::DIM> (coef))
  { ; }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new PML_CurlCurlEdgeIntegrator (coeffs[0]);
  }

  virtual string Name () const { return "PML_CurlCurlEdge"; }
};


template <int D>
class PML_MassEdgeIntegrator
  : public PML_BDBIntegrator<DiffOpIdEdge<D>, DiagDMat<D>, HCurlFiniteElement<D> >
{
public:
  PML_MassEdgeIntegrator (CoefficientFunction * coef)
    : PML_BDBIntegrator<DiffOpIdEdge<D>, DiagDMat<D>, HCurlFiniteElement<D> > 
  (DiagDMat<D> (coef))
  { ; }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new PML_MassEdgeIntegrator (coeffs[0]);
  }

  ///
  virtual string Name () const { return "PML_Massedge"; }
};

















///
template <int D>
class PML_DivDivHDivIntegrator
  : public PML_BDBIntegrator<DiffOpDivHDiv<D>, DiagDMat<1>, HDivFiniteElement<D> >
{
public:
  ///
  PML_DivDivHDivIntegrator (CoefficientFunction * coef)
    : PML_BDBIntegrator<DiffOpDivHDiv<D>, DiagDMat<1>, HDivFiniteElement<D> > 
  (DiagDMat<1> (coef))
  { ; }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new PML_DivDivHDivIntegrator (coeffs[0]);
  }

  ///
  virtual string Name () const { return "PML_DivDivHDiv"; }
};


///
template<int D>
class PML_MassHDivIntegrator
  : public PML_BDBIntegrator<DiffOpIdHDiv<D>, DiagDMat<D>, HDivFiniteElement<D> >
{
public:
  ///
  PML_MassHDivIntegrator (CoefficientFunction * coef)
    : PML_BDBIntegrator<DiffOpIdHDiv<D>, DiagDMat<D>, HDivFiniteElement<D> > 
  (DiagDMat<D> (coef))
  { ; }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new PML_MassHDivIntegrator (coeffs[0]);
  }

  ///
  virtual string Name () const { return "PML_MassHDiv"; }
};





///
template<int D>
class PML_RobinHDivIntegrator
  : public PML_BDBIntegrator<DiffOpIdHDivBoundary<D>, DiagDMat<1>, HDivNormalFiniteElement<D-1> >
{
public:
  ///
  PML_RobinHDivIntegrator (CoefficientFunction * coef)
    : PML_BDBIntegrator<DiffOpIdHDivBoundary<D>, DiagDMat<1>, HDivNormalFiniteElement<D-1> > 
  (DiagDMat<1> (coef))
  { ; }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new PML_RobinHDivIntegrator (coeffs[0]);
  }

  ///
  virtual string Name () const { return "PML_RobinHDiv"; }
};















#ifdef VERYSPECIAL




class PML_ReducedMassIntegrator : public BilinearFormIntegrator
{
  DiagDMat<1> dmatop;
  CoefficientFunction * coef;
public:
  
  enum { DIM_SPACE   = 2 };
  enum { DIM_ELEMENT = 2 };
  enum { DIM_DMAT    = 1 };
  enum { DIM         = 1 };
  

  PML_ReducedMassIntegrator (CoefficientFunction * acoef)
    : BilinearFormIntegrator(), dmatop(acoef), coef(acoef)
  { ; }

  ///
  virtual ~PML_ReducedMassIntegrator ()
  { ; }

  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new PML_ReducedMassIntegrator (coeffs[0]);
  }

  virtual bool BoundaryForm () const
  { return 0; }

  ///
  virtual void
  AssembleElementMatrix (const FiniteElement & bfel, 
			 const ElementTransformation & eltrans, 
			 FlatMatrix<double> & elmat,
			 LocalHeap & locheap) const
  {
    throw Exception ("PML cannot generae real matrices");
  }

  ///
  virtual void
  AssembleElementMatrix (const FiniteElement & bfel, 
			 const ElementTransformation & eltrans, 
			 FlatMatrix<Complex> & elmat,
			 LocalHeap & locheap) const
  {

    FE_Trig2 trig;
    ScalarFiniteElement<2> & felq = trig;
    const ScalarFiniteElement<2> & fel = dynamic_cast<const ScalarFiniteElement<2> &> (bfel);
    int ndof = fel.GetNDof();
    int ndofq = felq.GetNDof();

    elmat.AssignMemory (ndof*DIM, ndof*DIM, locheap);
    elmat = 0;
	
    FlatMatrix<Complex> matb(ndofq, ndof, locheap);
    FlatMatrix<Complex> matc(ndofq, ndofq, locheap);
    FlatMatrix<Complex> matcinv(ndofq, ndofq, locheap);


    const IntegrationRule & ir = 
      GetIntegrationRules().SelectIntegrationRule (fel.ElementType(), 6);
    
    matb = 0;
    matc = 0;
    void * heapp = locheap.GetPointer();
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,Complex> 
	  sip(ir[i], eltrans, locheap);
	SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,double> 
	  sip_real(ir[i], eltrans, locheap);
	
	
	const FlatVector<> shape = fel.GetShape (sip.IP(), locheap);
	const FlatVector<> shapeq = felq.GetShape (sip.IP(), locheap);
	

	double val = Evaluate (*coef, sip_real);
	Complex fac = sip.GetJacobiDet() * sip.IP().Weight();
	
	matb += fac * (shapeq * Trans (shape));
	matc += (fac/val) * (shapeq * Trans (shapeq));
	locheap.CleanUp (heapp);
      }

    CalcInverse (matc, matcinv);
    elmat = Trans (matb) * matcinv * matb;

    /*
    try
      {
	const ScalarFiniteElement & fel = dynamic_cast<const ScalarFiniteElement&> (bfel);
	int ndof = fel.GetNDof();

	elmat.AssignMemory (ndof*DIM, ndof*DIM, locheap);
	elmat = 0;
	
	FlatMatrixFixHeight<DIM_DMAT, Complex> bmat (ndof * DIM, locheap);
	FlatMatrixFixHeight<DIM_DMAT, Complex> dbmat (ndof * DIM, locheap);

	Mat<DIM_DMAT,DIM_DMAT, Complex> dmat;

	const IntegrationRule & ir = 
	  GetIntegrationRules().SelectIntegrationRule (fel.ElementType(), 2);

	void * heapp = locheap.GetPointer();
	for (int i = 0; i < ir.GetNIP(); i++)
	  {
	    SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,Complex> 
	      sip(ir[i], eltrans, locheap);
	    SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,double> 
	      sip_real(ir[i], eltrans, locheap);


	    const FlatVector<> shape = fel.GetShape (sip.IP(), locheap);
	    for (int j = 0; j < shape.Height(); j++)
	      bmat(0, j) = shape(j);
	    
	    dmatop.GenerateMatrix (fel, sip_real, dmat, locheap);
	    Complex fac = sip.GetJacobiDet() * sip.IP().Weight();

	    dbmat = dmat * bmat;
	    elmat += fac * (Trans (dbmat) * bmat);

	    locheap.CleanUp (heapp);
	  } 
      }

    catch (Exception & e)
      {
	e.Append ("in AssembleElementMatrix, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin AssembleElementMatrix, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }
    */
  }


  ///
  virtual int GetDimension () const
  { return 1; }

  ///
  virtual int Lumping () const
    { return 0; }
  ///
  virtual string Name () const { return "PML_ReducedMass integrator"; }
};












class PML_InterpolationMassIntegrator : public BilinearFormIntegrator
{
  DiagDMat<1> dmatop;
  CoefficientFunction * coef;
public:
  
  enum { DIM_SPACE   = 2 };
  enum { DIM_ELEMENT = 2 };
  enum { DIM_DMAT    = 1 };
  enum { DIM         = 1 };
  

  PML_InterpolationMassIntegrator (CoefficientFunction * acoef)
    : BilinearFormIntegrator(), dmatop(acoef), coef(acoef)
  { ; }

  ///
  virtual ~PML_InterpolationMassIntegrator ()
  { ; }

  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new PML_InterpolationMassIntegrator (coeffs[0]);
  }

  virtual bool BoundaryForm () const
  { return 0; }

  ///
  virtual void
  AssembleElementMatrix (const FiniteElement & bfel, 
			 const ElementTransformation & eltrans, 
			 FlatMatrix<double> & elmat,
			 LocalHeap & locheap) const
  {
    throw Exception ("PML cannot generae real matrices");
  }

  ///
  virtual void
  AssembleElementMatrix (const FiniteElement & bfel, 
			 const ElementTransformation & eltrans, 
			 FlatMatrix<Complex> & elmat,
			 LocalHeap & locheap) const
  {

    FE_Trig1 trig;
    ScalarFiniteElement<2> & felq = trig;
    const ScalarFiniteElement<2> & fel = dynamic_cast<const ScalarFiniteElement<2> &> (bfel);
    int ndof = fel.GetNDof();
    int ndofq = felq.GetNDof();

    elmat.AssignMemory (ndof*DIM, ndof*DIM, locheap);
    elmat = 0;

    FlatMatrix<Complex> matm(ndofq, ndofq, locheap);

    FlatMatrix<Complex> matb(ndofq, ndof, locheap);
    FlatMatrix<Complex> matb2(ndofq, ndof, locheap);
    FlatMatrix<Complex> matc(ndofq, ndofq, locheap);
    FlatMatrix<Complex> matcinv(ndofq, ndofq, locheap);

    const IntegrationRule & ir = 
      GetIntegrationRules().SelectIntegrationRule (fel.ElementType(), 6);
    const IntegrationRule & irseg = 
      GetIntegrationRules().SelectIntegrationRule (ET_SEGM, 0);
    
    matm = 0;
    matb = 0;
    matc = 0;

    FlatVector<> shapee (3, locheap);


    void * heapp = locheap.GetPointer();
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,Complex> 
	  sip(ir[i], eltrans, locheap);
	SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,double> 
	  sip_real(ir[i], eltrans, locheap);
	
	const FlatVector<> shape = fel.GetShape (sip.IP(), locheap);
	const FlatVector<> shapeq = felq.GetShape (sip.IP(), locheap);

	double val = Evaluate (*coef, sip_real);
	Complex fac = sip.GetJacobiDet() * sip.IP().Weight() * val;
	
	matm += fac * (shapeq * Trans (shapeq));
	locheap.CleanUp (heapp);
      }


    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    const POINT3D * points = ElementTopology::GetVertices (ET_TRIG);
    for (int en = 0; en < 3; en++)
      for (int i = 0; i < irseg.GetNIP(); i++)
	{
	  int v1 = edges[en][0];
	  int v2 = edges[en][1];
	  double lam = irseg[i](0);
	  IntegrationPoint ip(lam*points[v1][0]+(1-lam)*points[v2][0],
			      lam*points[v1][1]+(1-lam)*points[v2][1], 0, 0);
	  
	  SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,Complex> 
	    sip(ip, eltrans, locheap);
	  SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,double> 
	    sip_real(ip, eltrans, locheap);
	  
	  const FlatVector<> shape = fel.GetShape (sip.IP(), locheap);
	  const FlatVector<> shapeq = felq.GetShape (sip.IP(), locheap);

	  shapee = 0.0;
	  shapee(en) = 1;

	  Complex fac = irseg[i].Weight();

	  matb += fac * (shapee * Trans (shape));
	  matc += fac * (shapee * Trans (shapeq));
	  locheap.CleanUp (heapp);
	}


    CalcInverse (matc, matcinv);
    matb2 = matcinv * matb;
    elmat = Trans (matb2) * matm * matb2;

    (*testout) << "matb = " << endl << matb << endl
	       << "matc = " << endl << matc << endl
	       << "matb2 = " << endl << matb2 << endl
	       << "matm = " << endl << matm << endl;
  }


  ///
  virtual int GetDimension () const
  { return 1; }

  ///
  virtual int Lumping () const
    { return 0; }
  ///
  virtual string Name () const { return "PML_InterpolationMass integrator"; }
};






#endif













  namespace pml_cpp
 {
    class Init
    { 
    public: 
      Init ();
    };        
    

    Init::Init()
    {
      GetIntegrators().AddBFIntegrator ("PML_laplace", 2, 1,
					PML_LaplaceIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("PML_mass", 2, 1,
					PML_MassIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("PML_elasticity", 3, 2,
      					PML_ElasticityIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("PML_elasticity", 2, 2,
      					PML_ElasticityIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("PML_laplace", 3, 1,
					PML_LaplaceIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("PML_mass", 3, 1,
					PML_MassIntegrator<3>::Create);
      /*
      GetIntegrators().AddBFIntegrator ("PML_reducedmass", 2, 1,
					PML_ReducedMassIntegrator::Create);
      GetIntegrators().AddBFIntegrator ("PML_interpolationmass", 2, 1,
					PML_InterpolationMassIntegrator::Create);
      */
      GetIntegrators().AddBFIntegrator ("PML_curlcurledge", 3, 1,
					PML_CurlCurlEdgeIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("PML_curlcurledge", 2, 1,
					PML_CurlCurlEdgeIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("PML_massedge", 3, 1,
					PML_MassEdgeIntegrator<3>::Create);
      GetIntegrators().AddBFIntegrator ("PML_massedge", 2, 1,
					PML_MassEdgeIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("PML_divdivhdiv", 2, 1,
					PML_DivDivHDivIntegrator<2>::Create);
      GetIntegrators().AddBFIntegrator ("PML_masshdiv", 2, 1,
					PML_MassHDivIntegrator<2>::Create);
      //      GetIntegrators().AddBFIntegrator ("PML_robinhdiv", 2, 1,
      //					PML_RobinHDivIntegrator<2>::Create);
    }
    
    Init init;
    
  }
 }

