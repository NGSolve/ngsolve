#ifndef FILE_ELASTICITY_EQUATIONS
#define FILE_ELASTICITY_EQUATIONS

/*********************************************************************/
/* File:   elasticity_equations.hpp                                  */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/

namespace ngfem
{


  /*
    
  Elasticity integrators:

  */




  
  /* ********************  Elasticity ************************** */



  /// Elasticity operator $(e_{11},e_{22},2 e_{12})$
  template <int D, typename FEL = ScalarFiniteElement<D> > 
  class DiffOpStrain : public DiffOp<DiffOpStrain<D, FEL> >
  {
  };

  template <typename FEL>
  class DiffOpStrain<2, FEL> : public DiffOp<DiffOpStrain<2, FEL> >
  {
  public:
    enum { DIM = 2 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 3 };
    enum { DIFFORDER = 1 };

    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      typedef typename remove_reference_t<MAT>::TSCAL TSCAL;
      int nd = fel.GetNDof();

      FlatMatrixFixHeight<2, TSCAL> grad (nd, lh);
      FlatMatrixFixWidth<2> dshape(nd, lh);
      static_cast<const FEL&>(fel).CalcDShape(mip.IP(), dshape);
      grad = Trans (mip.GetJacobianInverse ()) * Trans (dshape);
      /*
      grad = Trans (mip.GetJacobianInverse ()) * 
	Trans (static_cast<const FEL&>(fel).GetDShape(mip.IP(), lh));
      */
      mat.AddSize(3, fel.GetNDof())  = TSCAL (0);
      for (int i = 0; i < nd; i++)
	{
	  mat(0, DIM*i  ) = grad(0, i);
	  mat(1, DIM*i+1) = grad(1, i);
	  mat(2, DIM*i  ) = grad(1, i);
	  mat(2, DIM*i+1) = grad(0, i);
	}
    }
  };




  template <typename FEL>
  class DiffOpStrain<3, FEL> : public DiffOp<DiffOpStrain<3, FEL> >
  {
  public:
    enum { DIM = 3 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 3 };
    enum { DIM_DMAT = 6 };
    enum { DIFFORDER = 1 };

    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      typedef typename remove_reference_t<MAT>::TSCAL TSCAL;

      HeapReset hr(lh);
      /*
      mat = TSCAL(0);
      static_cast<const FEL &>(fel).
	CalcDShape (mip.IP(),
		    [&](int i, Vec<3> gradref)
		    {
		      Vec<3,TSCAL> grad = 
			Trans (mip.GetJacobianInverse ()) * gradref;

		      mat(0, 3*i  ) = grad(0);
		      mat(1, 3*i+1) = grad(1);
		      mat(2, 3*i+2) = grad(2);
		      
		      mat(3, 3*i  ) = grad(1);
		      mat(3, 3*i+1) = grad(0);
		      
		      mat(4, 3*i  ) = grad(2);
		      mat(4, 3*i+2) = grad(0);
			    
		      mat(5, 3*i+1) = grad(2);
		      mat(5, 3*i+2) = grad(1);
		    });
      */      
      int nd = fel.GetNDof();
      FlatMatrixFixHeight<3,TSCAL> grad (nd, lh);
      FlatMatrixFixWidth<3> dshape(nd, lh);
      static_cast<const FEL&>(fel).CalcDShape(mip.IP(), dshape);
      grad = Trans (mip.GetJacobianInverse ()) * Trans (dshape);
      
      /*
      grad =  Trans (mip.GetJacobianInverse ()) * 
	Trans (static_cast<const FEL &>(fel).GetDShape(mip.IP(),lh));
      */
      mat.AddSize(6, nd) = TSCAL (0);
      for (int i = 0; i < nd; i++)
	{
	  mat(0, DIM*i  ) = grad(0, i);
	  mat(1, DIM*i+1) = grad(1, i);
	  mat(2, DIM*i+2) = grad(2, i);

	  mat(3, DIM*i  ) = grad(1, i);
	  mat(3, DIM*i+1) = grad(0, i);

	  mat(4, DIM*i  ) = grad(2, i);
	  mat(4, DIM*i+2) = grad(0, i);

	  mat(5, DIM*i+1) = grad(2, i);
	  mat(5, DIM*i+2) = grad(1, i);
	}
      
    }


#ifdef SERVICE
    template <typename AFEL, typename MAT>
    static void GenerateMatrix (const AFEL & fel, 
				const MappedIntegrationPoint<3,3> & mip,
				MAT & mat, LocalHeap & lh)
    {
      typedef typename MAT::TSCAL TSCAL;

      int nd = fel.GetNDof();
      HeapReset hr(lh);

      mat = TSCAL(0);
      static_cast<const FEL &>(fel).
	CalcMappedDShape (mip,
			  [&](int i, Vec<3> grad)
			  {
			    mat(0, 3*i  ) = grad(0);
			    mat(1, 3*i+1) = grad(1);
			    mat(2, 3*i+2) = grad(2);
			    
			    mat(3, 3*i  ) = grad(1);
			    mat(3, 3*i+1) = grad(0);
			    
			    mat(4, 3*i  ) = grad(2);
			    mat(4, 3*i+2) = grad(0);
			    
			    mat(5, 3*i+1) = grad(2);
			    mat(5, 3*i+2) = grad(1);
			  });

      /*
      FlatMatrixFixWidth<3> grad (nd, lh);
      static_cast<const FEL &>(fel).CalcMappedDShape (mip, grad);

      mat = TSCAL (0);
      for (int i = 0; i < nd; i++)
	{
	  mat(0, DIM*i  ) = grad(i, 0);
	  mat(1, DIM*i+1) = grad(i, 1);
	  mat(2, DIM*i+2) = grad(i, 2);

	  mat(3, DIM*i  ) = grad(i, 1);
	  mat(3, DIM*i+1) = grad(i, 0);

	  mat(4, DIM*i  ) = grad(i, 2);
	  mat(4, DIM*i+2) = grad(i, 0);

	  mat(5, DIM*i+1) = grad(i, 2);
	  mat(5, DIM*i+2) = grad(i, 1);
	}
      */
    }
#endif
  };



  /// 2D plane strain, and 3D
  template <int DIM>
  class ElasticityDMat : public DMatOp<ElasticityDMat<DIM>,DIM*(DIM+1)/2>
  {
  public:
    shared_ptr<CoefficientFunction> coefe;
    shared_ptr<CoefficientFunction> coefnu;
  public:
    enum { DIM_DMAT = (DIM * (DIM+1)) / 2 };  

    ElasticityDMat (shared_ptr<CoefficientFunction> acoefe,
		    shared_ptr<CoefficientFunction> acoefnu) 
      : coefe(acoefe), coefnu(acoefnu) { ; }

    template <typename SCAL>
    static Mat<DIM_DMAT,DIM_DMAT,SCAL> GetMatrixType(SCAL s) { return SCAL(0); }


    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      mat = 0;
      double nu = Evaluate (*coefnu, mip);
      double e = Evaluate (*coefe, mip);
      int i;
      for (i = 0; i < DIM; i++)
	{
	  mat(i,i) = 1-nu;
	  for (int j = 0; j < i; j++)
	    mat(i,j) = mat(j,i) = nu;
	}
      for (i = DIM; i < (DIM*(DIM+1)/2); i++)
	mat(i,i) = 0.5 * (1-2*nu);

      mat *= (e / ((1 + nu) * (1 - 2 * nu)));
    }  
  };


  ///
  template <int DIM>
  class OrthotropicElasticityDMat : public DMatOp<OrthotropicElasticityDMat<DIM>,
                                                  DIM*(DIM+1)/2>
  {
  public:
    CoefficientFunction * coefE1; // Young's moduli
    CoefficientFunction * coefE2;
    CoefficientFunction * coefE3;
    CoefficientFunction * coefnu12; // Poisson's ratios (nu21/E2 = nu12/E1, nu31/E3 = nu13/E1, nu32/E3 = nu23/E2)
    CoefficientFunction * coefnu13;
    CoefficientFunction * coefnu23;
    CoefficientFunction * coefG12; // shear moduil
    CoefficientFunction * coefG13;
    CoefficientFunction * coefG23;
  public:
    enum { DIM_DMAT = (DIM * (DIM+1)) / 2 };  

    OrthotropicElasticityDMat (const Array<shared_ptr<CoefficientFunction>> & coefs)
      {
        cerr << "OrthotropicElasticityDMat currently not available" << endl;
      }
    /*
    OrthotropicElasticityDMat (CoefficientFunction * acoefE1,
			       CoefficientFunction * acoefE2,
			       CoefficientFunction * acoefE3,
			       CoefficientFunction * acoefnu12,
			       CoefficientFunction * acoefnu13,
			       CoefficientFunction * acoefnu23,
			       CoefficientFunction * acoefG12,
			       CoefficientFunction * acoefG13,
			       CoefficientFunction * acoefG23) 
      : coefE1(acoefE1), coefE2(acoefE2), coefE3(acoefE3),
	coefnu12(acoefnu12), coefnu13(acoefnu13), coefnu23(acoefnu23),
	coefG12(acoefG12), coefG13(acoefG13), coefG23(acoefG23) { ; }
    */

    template <typename SCAL>
    static Mat<DIM_DMAT,DIM_DMAT,SCAL> GetMatrixType(SCAL s) { return SCAL(0); }

    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      mat = 0;
      const double E1 = Evaluate (*coefE1, mip);
      const double E2 = Evaluate (*coefE2, mip);
      const double E3 = Evaluate (*coefE3, mip);

      if(E1 < 1.e-5 || E2 < 1.e-5 || E3 < 1.e-5) return;

      const double nu12 = Evaluate (*coefnu12, mip);
      const double nu21 = nu12*(E2/E1);
      const double nu13 = Evaluate (*coefnu13, mip);
      const double nu31 = nu13*(E3/E1);
      const double nu23 = Evaluate (*coefnu23, mip);
      const double nu32 = nu23*(E3/E2);

      if(nu12 < 0 || nu12 > 0.5 || nu21 < 0 || nu21 > 0.5 || nu13 < 0 || nu13 > 0.5 || nu31 < 0 || nu31 > 0.5 || nu23 < 0 || nu23 > 0.5 || nu32 < 0 || nu32 > 0.5)
	{
	  cerr << "WARNING: Bad choice for elasticity constants: " << endl
	       << "E1 " << E1 << " E2 " << E2 << " E3 " << E3 << endl
	       << "nu12 " << nu12 << " nu21 " << nu21 << " nu13 " << nu13 << " nu31 " << nu31 << " nu23 " << nu23 << " nu32 " << nu32 <<endl;
	}

      const double denom = 1. - nu13*nu32*nu21 - nu12*nu23*nu31 - nu12*nu21 - nu13*nu31 - nu23*nu32;

      mat(0,0) = E1*(1.-nu23*nu32)/denom; 
      mat(1,0) = mat(0,1) = E2*(nu12+nu13*nu32)/denom; mat(1,1) = E2*(1.-nu13*nu31)/denom;
      mat(2,0) = mat(0,2) = E3*(nu13+nu12*nu23)/denom; mat(2,1) = mat(1,2) = E3*(nu23+nu13*nu21)/denom; mat(2,2) = E3*(1.-nu12*nu21)/denom;

      mat(3,3) = Evaluate (*coefG12, mip);
      mat(4,4) = Evaluate (*coefG13, mip);
      mat(5,5) = Evaluate (*coefG23, mip);
    }  
  };

  /// Orthotropic Elasticity DMat with Cylindrical Coordinates
  template <int DIM>
  class OrthotropicCylElasticityDMat : public DMatOp<OrthotropicElasticityDMat<DIM>,
                                                     DIM*(DIM+1)/2>
  {
  public:
    CoefficientFunction * coefE1; // Young's moduli
    CoefficientFunction * coefE2;
    CoefficientFunction * coefE3;
    CoefficientFunction * coefnu12; // Poisson's ratios (nu21/E2 = nu12/E1, nu31/E3 = nu13/E1, nu32/E3 = nu23/E2)
    CoefficientFunction * coefnu13;
    CoefficientFunction * coefnu23;
    CoefficientFunction * coefG12; // shear moduil
    CoefficientFunction * coefG13;
    CoefficientFunction * coefG23;
    CoefficientFunction * coefUseCyl; // if 1 ... use cylindrical coordinates, if 0 ... standard ortot.
  public:
    enum { DIM_DMAT = (DIM * (DIM+1)) / 2 };  
    
    OrthotropicCylElasticityDMat (const Array<shared_ptr<CoefficientFunction>> & coefs)
      {
        cerr << "OrthotropicCylElasticityDMat currently not available" << endl;
      }


    OrthotropicCylElasticityDMat (CoefficientFunction * acoefE1,
				  CoefficientFunction * acoefE2,
				  CoefficientFunction * acoefE3,
				  CoefficientFunction * acoefnu12,
				  CoefficientFunction * acoefnu13,
				  CoefficientFunction * acoefnu23,
				  CoefficientFunction * acoefG12,
				  CoefficientFunction * acoefG13,
				  CoefficientFunction * acoefG23,
				  CoefficientFunction * acoefUseCyl) 
      : coefE1(acoefE1), coefE2(acoefE2), coefE3(acoefE3),
	coefnu12(acoefnu12), coefnu13(acoefnu13), coefnu23(acoefnu23),
	coefG12(acoefG12), coefG13(acoefG13), coefG23(acoefG23), coefUseCyl(acoefUseCyl) { ; }


    // template <typename SCAL>
    // static Mat<DIM_DMAT,DIM_DMAT,SCAL> GetMatrixType(SCAL s) { return SCAL(0); }


    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
    
      double E1 = Evaluate (*coefE1, mip);
      double E2 = Evaluate (*coefE2, mip);
      double E3 = Evaluate (*coefE3, mip);

      if(E1 < 1.e-5 || E2 < 1.e-5 || E3 < 1.e-5) return;

      double nu12 = Evaluate (*coefnu12, mip);
      double nu21 = nu12*(E2/E1);
      double nu13 = Evaluate (*coefnu13, mip);
      double nu31 = nu13*(E3/E1);
      double nu23 = Evaluate (*coefnu23, mip);
      double nu32 = nu23*(E3/E2);

      const double useCyl = Evaluate (*coefUseCyl, mip);

      double G12 = Evaluate (*coefG12, mip);
      double G13 = Evaluate (*coefG13, mip);
      double G23 = Evaluate (*coefG23, mip);



      double n1 = mip.GetPoint()(0);
      double n2 = mip.GetPoint()(1);
      const double l = sqrt(n1*n1+n2*n2);

      n1 /= l; n2 /= l;

      if(nu12 < 0 || nu12 > 0.5 || nu21 < 0 || nu21 > 0.5 || nu13 < 0 || nu13 > 0.5 || nu31 < 0 || nu31 > 0.5 || nu23 < 0 || nu23 > 0.5 || nu32 < 0 || nu32 > 0.5)
	{
	  cerr << "WARNING: Bad choice for elasticity constants: " << endl
	       << "E1 " << E1 << " E2 " << E2 << " E3 " << E3 << endl
	       << "nu12 " << nu12 << " nu21 " << nu21 << " nu13 " << nu13 << " nu31 " << nu31 << " nu23 " << nu23 << " nu32 " << nu32 <<endl;
	}

      const double denom = 1. - nu13*nu32*nu21 - nu12*nu23*nu31 - nu12*nu21 - nu13*nu31 - nu23*nu32;


    

      MAT aux(mat),transf(mat);

      aux = 0;

      aux(0,0) = E1*(1.-nu23*nu32)/denom; 
      aux(1,0) = aux(0,1) = E2*(nu12+nu13*nu32)/denom; aux(1,1) = E2*(1.-nu13*nu31)/denom;
      aux(2,0) = aux(0,2) = E3*(nu13+nu12*nu23)/denom; aux(2,1) = aux(1,2) = E3*(nu23+nu13*nu21)/denom; aux(2,2) = E3*(1.-nu12*nu21)/denom;

      aux(3,3) = G12;
      aux(4,4) = G13;
      aux(5,5) = G23;

      if(fabs(useCyl) > 0.5)
	{
	  transf = 0;

	  transf(0,0) = transf(1,1) = n1*n1; transf(0,1) = transf(1,0) = n2*n2; transf(2,2) = 1.;
	  transf(0,3) = 2.*n1*n2; transf(1,3) = -2.*n1*n2;
	  transf(3,0) = -n1*n2; transf(3,1) = n1*n2;
	  transf(3,3) = n1*n1-n2*n2; transf(4,4) = transf(5,5) = n1; transf(4,5) = n2; transf(5,4) = -n2;
	
	  mat = Trans(transf)*aux*transf;
	}
      else
	{
	  mat = aux;
	}
    }  
  };



  ///
  class PlaneStressDMat : public DMatOp<PlaneStressDMat,3>
  {
    CoefficientFunction * coefe;
    CoefficientFunction * coefnu;
  public:
    enum { DIM_DMAT = 3 };
  
    PlaneStressDMat (CoefficientFunction * acoefe,
		     CoefficientFunction * acoefnu) 
      : coefe(acoefe), coefnu(acoefnu) { ; }
  
    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      mat = 0;
      double nu = Evaluate (*coefnu, mip);
      double e = Evaluate (*coefe, mip);

      mat(0,0) = mat(1,1) = 1;
      mat(0,1) = mat(1,0) = nu;
      mat(2,2) = (1-nu) / 2;

      mat *= (e / (1 - nu * nu));
    }  
  };

  ///
  template <int D>
  class ElasticityIntegrator 
    : public T_BDBIntegrator<DiffOpStrain<D>, ElasticityDMat<D>, ScalarFiniteElement<D> >
  {
    typedef T_BDBIntegrator<DiffOpStrain<D>, ElasticityDMat<D>, ScalarFiniteElement<D> > BASE;
  public:
    ElasticityIntegrator (shared_ptr<CoefficientFunction> coefe,
			  shared_ptr<CoefficientFunction> coefnu)
      : BASE(ElasticityDMat<D> (coefe, coefnu))
    { ; }

    ElasticityIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : BASE(ElasticityDMat<D> (coeffs[0], coeffs[1]))
    { ; }

    /*
    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new ElasticityIntegrator (coeffs[0], coeffs[1]);
    }
    */

    ///
    virtual string Name () const { return "Elasticity"; }
  };


  /*
  // for plane stress
  ///
  template <>
  class ElasticityIntegrator <2>
  : public T_BDBIntegrator<DiffOpStrain<2>, PlaneStressDMat, ScalarFiniteElement<D> >
  {
  public:
  ///
  ElasticityIntegrator (CoefficientFunction * coefe,
  CoefficientFunction * coefnu)
  : T_BDBIntegrator<DiffOpStrain<2>, PlaneStressDMat, ScalarFiniteElement<D> > 
  (PlaneStressDMat (coefe, coefnu))
  { ; }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
  return new ElasticityIntegrator (coeffs[0], coeffs[1]);
  }

  ///
  virtual string Name () const { return "Elasticity"; }
  };
  */


  ///
  template <int D>
  class OrthotropicElasticityIntegrator 
    : public T_BDBIntegrator<DiffOpStrain<D>, OrthotropicElasticityDMat<D>, ScalarFiniteElement<D> >
  {
    typedef T_BDBIntegrator<DiffOpStrain<D>, OrthotropicElasticityDMat<D>, ScalarFiniteElement<D> > BASE;
  public:
    using T_BDBIntegrator<DiffOpStrain<D>, OrthotropicElasticityDMat<D>, ScalarFiniteElement<D> >::T_BDBIntegrator;
    ///
    /*
    OrthotropicElasticityIntegrator (CoefficientFunction * coefE1,
				     CoefficientFunction * coefE2,
				     CoefficientFunction * coefE3,
				     CoefficientFunction * coefnu12,
				     CoefficientFunction * coefnu13,
				     CoefficientFunction * coefnu23,
				     CoefficientFunction * coefG12,
				     CoefficientFunction * coefG13,
				     CoefficientFunction * coefG23)
      : T_BDBIntegrator<DiffOpStrain<D>, OrthotropicElasticityDMat<D>, ScalarFiniteElement<D> > 
    (OrthotropicElasticityDMat<D> (coefE1, coefE2, coefE3, coefnu12, coefnu13, coefnu23, coefG12, coefG13, coefG23))
    { ; }
  
    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new OrthotropicElasticityIntegrator (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4], coeffs[5], coeffs[6], coeffs[7], coeffs[8]);
    }
    */


    ///
    virtual string Name () const { return "OrthotropicElasticity"; }
  };


  ///
  template <int D>
  class OrthotropicCylElasticityIntegrator 
    : public T_BDBIntegrator<DiffOpStrain<D>, OrthotropicCylElasticityDMat<D>, ScalarFiniteElement<D> >
  {
    typedef T_BDBIntegrator<DiffOpStrain<D>, OrthotropicCylElasticityDMat<D>, ScalarFiniteElement<D> > BASE;
  public:
    using T_BDBIntegrator<DiffOpStrain<D>, OrthotropicCylElasticityDMat<D>, ScalarFiniteElement<D> >::T_BDBIntegrator;
    ///
    /*
    OrthotropicCylElasticityIntegrator (CoefficientFunction * coefE1,
					CoefficientFunction * coefE2,
					CoefficientFunction * coefE3,
					CoefficientFunction * coefnu12,
					CoefficientFunction * coefnu13,
					CoefficientFunction * coefnu23,
					CoefficientFunction * coefG12,
					CoefficientFunction * coefG13,
					CoefficientFunction * coefG23,
					CoefficientFunction * coefUseCyl)
      : T_BDBIntegrator<DiffOpStrain<D>, OrthotropicCylElasticityDMat<D>, ScalarFiniteElement<D> > 
    (OrthotropicCylElasticityDMat<D> (coefE1, coefE2, coefE3, coefnu12, coefnu13, coefnu23, coefG12, coefG13, coefG23, coefUseCyl))
    { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new OrthotropicCylElasticityIntegrator (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4], coeffs[5], coeffs[6], coeffs[7], coeffs[8], coeffs[9]);
    }
    */


    ///
    virtual string Name () const { return "OrthotropicCylElasticity"; }
  };

}




#endif
