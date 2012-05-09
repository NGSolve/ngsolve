#ifndef FILE_HCURL_EQUATIONS
#define FILE_HCURL_EQUATIONS

/*********************************************************************/
/* File:   hcurl_equations.hpp                                       */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/

namespace ngfem
{


  /*
    
  Maxwell integrators:


  Finite Element Integrators for H(curl)

  Mapping with covariant transformation

  Requires H(curl) finite elements
  */






  /// Identity operator, covariant transformation
  template <int D, typename FEL = HCurlFiniteElement<D> >
  class DiffOpIdEdge : public DiffOp<DiffOpIdEdge<D, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 0 };

    template <typename FEL1, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL1 & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      mat = Trans (mip.GetJacobianInverse ()) * 
	Trans (static_cast<const FEL&> (fel).GetShape(mip.IP(), lh));
    }

    template <typename FEL1>
    static void GenerateMatrix (const FEL1 & fel, 
				const MappedIntegrationPoint<D,D> & mip,
				FlatMatrixFixHeight<D> & mat, LocalHeap & lh)
    {
      static_cast<const FEL&> (fel).CalcMappedShape (mip, mat.Trans());
    }


    template <typename FEL1, typename MIP, class TVX, class TVY>
    static void Apply (const FEL1 & fel, const MIP & mip,
		       const TVX & x, TVY & y,
		       LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCAL;

      Vec<D,TSCAL> hx;
      hx = Trans (static_cast<const FEL&> (fel).GetShape (mip.IP(), lh)) * x;
      y = Trans (mip.GetJacobianInverse()) * hx;
    }

    template <typename FEL1, typename MIP, class TVX, class TVY>
    static void ApplyTrans (const FEL1 & fel, const MIP & mip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCAL;

      Vec<D,TSCAL> hx;
      hx = mip.GetJacobianInverse() * x;
      y = static_cast<const FEL&> (fel).GetShape (mip.IP(),lh) * hx;
    }
  };








  /// Operator $curl$, Piola-transformation
  template <int D>
  class DiffOpCurlEdge : public DiffOp<DiffOpCurlEdge<D> >
  {
  };


  template <> class DiffOpCurlEdge<2> : public DiffOp<DiffOpCurlEdge<2> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      mat = 1.0/mip.GetJacobiDet() * 
	Trans (fel.GetCurlShape(mip.IP(), lh));
    }


    template <typename FEL, typename MIP, class TVX, class TVY>
    static void Apply (const FEL & fel, const MIP & mip,
		       const TVX & x, TVY & y,
		       LocalHeap & lh) 
    {
      y = (1.0/mip.GetJacobiDet()) * 
	(Trans (fel.GetCurlShape(mip.IP(), lh)) * x);
    }
  };

  template <> class DiffOpCurlEdge<3> : public DiffOp<DiffOpCurlEdge<3> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 3 };
    enum { DIM_DMAT = 3 };
    enum { DIFFORDER = 1 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      mat = (1.0/mip.GetJacobiDet())
	* (mip.GetJacobian() * Trans (fel.GetCurlShape(mip.IP(), lh)));
    }

    template <typename FEL>
    static void GenerateMatrix (const FEL & fel, 
				const MappedIntegrationPoint<3,3> & mip,
				FlatMatrixFixHeight<3> & mat, LocalHeap & lh)
    {
      fel.CalcMappedCurlShape (mip, mat.Trans());
    }



    template <typename FEL, typename MIP, class TVX, class TVY>
    static void Apply (const FEL & fel, const MIP & mip,
		       const TVX & x, TVY & y,
		       LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCAL;

      Vec<3,TSCAL> hx;
      hx = fel.EvaluateCurlShape (mip.IP(), x, lh);
      y = (1.0/mip.GetJacobiDet()) * (mip.GetJacobian() * hx);
    }


    template <typename FEL, typename MIP, class TVX, class TVY>
    static void ApplyTrans (const FEL & fel, const MIP & mip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCAL;

      Vec<3,TSCAL> hx;
      hx = (1.0/mip.GetJacobiDet()) * (Trans (mip.GetJacobian()) * x);
      y = fel.GetCurlShape(mip.IP(), lh) * hx;
    }
  };

















  // \int_{C} v.\tau
  template <int D>
  class DiffOpTangentialComponentEdge : public DiffOp<DiffOpTangentialComponentEdge<D> >
  {
  public:
    enum { DIM = D };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      Vec<D> tv = mip.GetTV();
      Vec<D> tv_JI = mip.GetJacobianInverse () * tv;
   
      mat = Trans ( fel.GetShape(mip.IP(), lh) * tv_JI );
    
    }

  };



  /// Identity on boundary
  template <int D, typename FEL = HCurlFiniteElement<D-1> >
  class DiffOpIdBoundaryEdge : public DiffOp<DiffOpIdBoundaryEdge<D,FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 0 };

    template <typename FEL1, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL1 & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      mat = Trans (mip.GetJacobianInverse ()) * 
	Trans (static_cast<const FEL&> (fel).GetShape(mip.IP(),lh));
    }

    template <typename FEL1, typename MIP, class TVX, class TVY>
    static void Apply (const FEL1 & fel, const MIP & mip,
		       const TVX & x, TVY & y,
		       LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCAL;

      Vec<D-1,TSCAL> hx;
      hx = Trans (static_cast<const FEL&> (fel).GetShape (mip.IP(),lh)) * x;
      y = Trans (mip.GetJacobianInverse()) * hx;
    }

    template <typename FEL1, typename MIP, class TVX, class TVY>
    static void ApplyTrans (const FEL1 & fel, const MIP & mip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCAL;

      Vec<DIM_ELEMENT,TSCAL> hx;
      hx = mip.GetJacobianInverse() * x;
      y = static_cast<const FEL&> (fel).GetShape (mip.IP(),lh) * hx;

      /*
      FlatMatrixFixWidth<DIM_ELEMENT> mshape (y.Height(), &hv(0)); 
      FlatMatrix<> mshape2 (y.Height(), DIM_ELEMENT, &hv(0)); 
      y = mshape2 * hx; 
      */
    }
  };



  /// Curl on boundary
  class DiffOpCurlBoundaryEdge : public DiffOp<DiffOpCurlBoundaryEdge>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      mat = 1.0/mip.GetJacobiDet() * Trans (fel.GetCurlShape(mip.IP(),lh));
    }


    template <typename FEL, typename MIP, class TVX, class TVY>
    static void Apply (const FEL & fel, const MIP & mip,
		       const TVX & x, TVY & y,
		       LocalHeap & lh) 
    {
      y = (1.0/mip.GetJacobiDet()) * (Trans (fel.GetCurlShape(mip.IP(),lh)) * x);
    }

    template <typename FEL, typename MIP, class TVX, class TVY>
    static void ApplyTrans (const FEL & fel, const MIP & mip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCAL;
      y = fel.GetCurlShape(mip.IP(),lh) * ((1.0/mip.GetJacobiDet()) * x);
    
    }
  };









  // bilinearform integrators





  /// 
  template <int D, typename FEL = HCurlFiniteElement<D> >
  class CurlCurlEdgeIntegrator 
    : public T_BDBIntegrator<DiffOpCurlEdge<D>, DiagDMat<DIM_CURL_TRAIT<D>::DIM>, FEL>
  {
  public:
    CurlCurlEdgeIntegrator (CoefficientFunction * coeff);
    CurlCurlEdgeIntegrator (Array<CoefficientFunction*> & coeffs);
  
    /*
    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new CurlCurlEdgeIntegrator (coeffs[0]);
    }
    */
    virtual string Name () const { return "CurlCurlEdge"; }
  };



  /// 
  class CurlCurlBoundaryEdgeIntegrator 
    : public T_BDBIntegrator<DiffOpCurlBoundaryEdge, DiagDMat<1>, HCurlFiniteElement<2> >
  {
  public:
    ///
    CurlCurlBoundaryEdgeIntegrator (CoefficientFunction * coeff)
      : T_BDBIntegrator<DiffOpCurlBoundaryEdge, DiagDMat<1>, HCurlFiniteElement<2> > 
    (DiagDMat<1> (coeff))
    { ; }
  
    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new CurlCurlBoundaryEdgeIntegrator (coeffs[0]);
    }

    ///
    virtual bool BoundaryForm () const { return 1; }
    ///
    virtual string Name () const { return "CurlCurlBoundaryEdge"; }
  };

  /// 
  template <int D, typename FEL = HCurlFiniteElement<D> >
  class CurlCurlEdgeOrthoIntegrator 
    : public T_BDBIntegrator<DiffOpCurlEdge<D>, OrthoDMat<DIM_CURL_TRAIT<D>::DIM>, FEL>
  {
  public:
    ///
    CurlCurlEdgeOrthoIntegrator (CoefficientFunction * coeff1,
				 CoefficientFunction * coeff2,
				 CoefficientFunction * coeff3)
      : T_BDBIntegrator<DiffOpCurlEdge<D>, OrthoDMat<DIM_CURL_TRAIT<D>::DIM>, FEL> 
    (OrthoDMat<DIM_CURL_TRAIT<D>::DIM> (coeff1, coeff2, coeff3))
    { ; }
  
    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new CurlCurlEdgeOrthoIntegrator (coeffs[0], coeffs[1], coeffs[2]);
    }

    ///
    virtual string Name () const { return "CurlCurlEdgeOrtho"; }
  };




  ///
  template <int D, typename FEL = HCurlFiniteElement<D> >
  class MassEdgeIntegrator 
    : public T_BDBIntegrator<DiffOpIdEdge<D>, DiagDMat<D>, FEL>
  {
  public:
    ///
    MassEdgeIntegrator (CoefficientFunction * coeff);
    MassEdgeIntegrator (Array<CoefficientFunction*> & coeffs);
    /*
      : T_BDBIntegrator<DiffOpIdEdge<D>, DiagDMat<D>, FEL> (DiagDMat<D> (coeff))
    { ; }
    */


    /*
    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new MassEdgeIntegrator (coeffs[0]);
    }
    */
    ///
    virtual string Name () const { return "MassEdge"; }
  };


  ///
  template <int D, typename FEL = HCurlFiniteElement<D> >
  class MassEdgeOrthoIntegrator 
    : public T_BDBIntegrator<DiffOpIdEdge<D>, OrthoDMat<D>, FEL>
  {
  public:
    ///
    MassEdgeOrthoIntegrator (CoefficientFunction * coeff1,
			     CoefficientFunction * coeff2)
      : T_BDBIntegrator<DiffOpIdEdge<D>, OrthoDMat<D>, FEL> (OrthoDMat<D> (coeff1, coeff2))
    { ; }

    MassEdgeOrthoIntegrator (CoefficientFunction * coeff1,
			     CoefficientFunction * coeff2,
			     CoefficientFunction * coeff3)
      : T_BDBIntegrator<DiffOpIdEdge<D>, OrthoDMat<D>, FEL> (OrthoDMat<D> (coeff1, coeff2, coeff3))
    { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      if (D == 2)
	return new MassEdgeOrthoIntegrator (coeffs[0], coeffs[1]);
      else
	return new MassEdgeOrthoIntegrator (coeffs[0], coeffs[1], coeffs[2]);
    }
  
    ///
    virtual string Name () const { return "MassEdgeOrtho"; }
  };


  ///
  template <int D, typename FEL = HCurlFiniteElement<D> >
  class MassEdgeAnisotropicIntegrator 
    : public T_BDBIntegrator<DiffOpIdEdge<D>, SymDMat<D>, FEL>
  { 
  };



  template <> 
  class MassEdgeAnisotropicIntegrator<3, HCurlFiniteElement<3> >
    : public T_BDBIntegrator<DiffOpIdEdge<3>, SymDMat<3>, HCurlFiniteElement<3> >
  {
  public:
    ///
    MassEdgeAnisotropicIntegrator (CoefficientFunction * coeff00,
				   CoefficientFunction * coeff10,
				   CoefficientFunction * coeff11,
				   CoefficientFunction * coeff20,
				   CoefficientFunction * coeff21,
				   CoefficientFunction * coeff22)
      : T_BDBIntegrator<DiffOpIdEdge<3>, SymDMat<3>, HCurlFiniteElement<3> >
    (SymDMat<3> (coeff00, coeff10, coeff11, coeff20, coeff21, coeff22))
    { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new MassEdgeAnisotropicIntegrator (coeffs[0], coeffs[1], coeffs[2],
						coeffs[3], coeffs[4], coeffs[5]);
    }
  
    ///
    virtual string Name () const 
    { return "MassEdgeAnisotropic"; }
  };





  ///
  template <int D, typename FEL = HCurlFiniteElement<D-1> >
  class RobinEdgeIntegrator 
    : public T_BDBIntegrator<DiffOpIdBoundaryEdge<D,FEL>, DiagDMat<D>, FEL>
  {
  public:
    ///
    RobinEdgeIntegrator (CoefficientFunction * coeff)
      : T_BDBIntegrator<DiffOpIdBoundaryEdge<D,FEL>, DiagDMat<D>, FEL> (DiagDMat<D> (coeff))
    { ; }
    RobinEdgeIntegrator (Array<CoefficientFunction*> & coeffs)
      : T_BDBIntegrator<DiffOpIdBoundaryEdge<D,FEL>, DiagDMat<D>, FEL> (DiagDMat<D> (coeffs[0]))
    { ; }

    ///
    virtual bool BoundaryForm () const { return 1; }

    ///
    virtual string Name () const { return "RobinEdge"; }
  };







  // Linearform integrators 

  template <int D, typename FEL = HCurlFiniteElement<D> >
  class NGS_DLL_HEADER SourceEdgeIntegrator
    : public T_BIntegrator<DiffOpIdEdge<D>, DVec<D>, FEL>
  {
  public:
    SourceEdgeIntegrator (Array<CoefficientFunction*> & coeffs);

    SourceEdgeIntegrator (CoefficientFunction * coeff1);

    SourceEdgeIntegrator (CoefficientFunction * coeff1,
			  CoefficientFunction * coeff2);

    SourceEdgeIntegrator (CoefficientFunction * coeff1,
			  CoefficientFunction * coeff2,
			  CoefficientFunction * coeff3);

    virtual string Name () const { return "SourceEdge"; }
  };








  ///
  template <int D, typename FEL = HCurlFiniteElement<D> >
  class TangentialSourceEdgeIntegrator 
    : public T_BIntegrator<DiffOpIdEdge<D>, TVec<D>, FEL>
  {
  public:
    ///
    TangentialSourceEdgeIntegrator (CoefficientFunction * coeff)
      : T_BIntegrator<DiffOpIdEdge<D>, TVec<D>, FEL> 
    (TVec<D> (coeff))
    { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new TangentialSourceEdgeIntegrator<D,FEL> (coeffs[0]);
    }
  
    ///
    virtual string Name () const { return "TangentialSourceEdge"; }
  };
  


  ///
  template <int D, typename FEL = HCurlFiniteElement<D-1> >
  class NeumannEdgeIntegrator
    : public T_BIntegrator<DiffOpIdBoundaryEdge<D,FEL>, DVec<D>, FEL>
  {
  public:
    ///
    NeumannEdgeIntegrator (CoefficientFunction * coeff1,
			   CoefficientFunction * coeff2,
			   CoefficientFunction * coeff3)
      : T_BIntegrator<DiffOpIdBoundaryEdge<D,FEL>,DVec<D>, FEL> 
    (DVec<D> (coeff1, coeff2, coeff3))
    { ; }
    NeumannEdgeIntegrator (CoefficientFunction * coeff1,
			   CoefficientFunction * coeff2)
      : T_BIntegrator<DiffOpIdBoundaryEdge<D,FEL>,DVec<D>, FEL> 
    (DVec<D> (coeff1, coeff2))
    { ; }

    NeumannEdgeIntegrator (Array<CoefficientFunction*> & coeffs)
      : T_BIntegrator<DiffOpIdBoundaryEdge<D,FEL>,DVec<D>, FEL> 
    (DVec<D> (coeffs))
    { ; }
      


    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      if (D == 3)
	return new NeumannEdgeIntegrator<3> (coeffs[0], coeffs[1], coeffs[2]);
      else
	return new NeumannEdgeIntegrator<2> (coeffs[0], coeffs[1]);
    }
  
    ///
    virtual bool BoundaryForm () const { return 1; }
    ///
    virtual string Name () const { return "NeumannEdge"; }
  };






  ///
  template <int D, typename FEL = HCurlFiniteElement<D> >
  class CurlEdgeIntegrator 
    : public T_BIntegrator<DiffOpCurlEdge<D>, DVec<DIM_CURL_TRAIT<D>::DIM>, FEL>
  {
  public:
    ///
    CurlEdgeIntegrator (CoefficientFunction * coeff1)
      : T_BIntegrator<DiffOpCurlEdge<D>, DVec<DIM_CURL_TRAIT<D>::DIM>, FEL> 
    (DVec<DIM_CURL_TRAIT<D>::DIM> (coeff1))
    { ; }

    CurlEdgeIntegrator (CoefficientFunction * coeffx,
			CoefficientFunction * coeffy,
			CoefficientFunction * coeffz)
      : T_BIntegrator<DiffOpCurlEdge<D>, DVec<DIM_CURL_TRAIT<D>::DIM>, FEL> 
    (DVec<DIM_CURL_TRAIT<D>::DIM> (coeffx, coeffy, coeffz))
    { ; }



    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      if (D == 2)
	return new CurlEdgeIntegrator<2,FEL> (coeffs[0]);
      else
	return new CurlEdgeIntegrator<3,FEL> (coeffs[0], coeffs[1], coeffs[2]);
    }
  
    ///
    virtual bool BoundaryForm () const { return 0; }
    ///
    virtual string Name () const { return "CurlEdge"; }
  };




  ///
  template <typename FEL = HCurlFiniteElement<2> >
  class CurlBoundaryEdgeIntegrator 
    : public T_BIntegrator<DiffOpCurlBoundaryEdge, DVec<1>, FEL>
  {
  public:
    ///
    CurlBoundaryEdgeIntegrator (CoefficientFunction * coeff1)
      : T_BIntegrator<DiffOpCurlBoundaryEdge, DVec<1>, FEL> 
    (DVec<1> (coeff1))
    { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new CurlBoundaryEdgeIntegrator (coeffs[0]);
    }
  
    ///
    virtual bool BoundaryForm () const { return 1; }
    ///
    virtual string Name () const { return "CurlBoundaryEdge"; }
  };






}

#endif
