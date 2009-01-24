#ifndef FILE_HDIV_EQUATIONS
#define FILE_HDIV_EQUATIONS

/*********************************************************************/
/* File:   hdiv_equations.hpp                                        */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/

/*
   Finite Element Integrators
*/



//Almedin



  // Id and div for H(div) elements, Piola transformation



  /// Identity operator, Piolas transformation
  template <int D>
  class DiffOpIdHDiv : public DiffOp<DiffOpIdHDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 0 };
    
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const SIP & sip,
				MAT & mat, LocalHeap & lh)
    {
      mat = (1.0/sip.GetJacobiDet()) * (sip.GetJacobian() * Trans (fel.GetShape(sip.IP(), lh)));
    }
    
    
    ///
    template <typename FEL, typename SIP, class TVX, class TVY>
    static void Apply (const FEL & fel, const SIP & sip,
		       const TVX & x, TVY & y,
		       LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCAL;
      
      Vec<D,TSCAL> hv = Trans (fel.GetShape(sip.IP(), lh)) * x;
      hv *= (1.0/sip.GetJacobiDet());
      y = sip.GetJacobian() * hv;
    }

  ///
  template <typename FEL, typename SIP, class TVX, class TVY>
  static void ApplyTrans (const FEL & fel, const SIP & sip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh) 
  {
    typedef typename TVX::TSCAL TSCAL;

    Vec<D,TSCAL> hv = Trans (sip.GetJacobian()) * x;
    hv *= (1.0/sip.GetJacobiDet());
    y = fel.GetShape(sip.IP(),lh) * hv;
  }
  };




  /// Operator $div$
  template <int D>
  class DiffOpDivHDiv : public DiffOp<DiffOpDivHDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };

    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const SIP & sip,
				MAT & mat, LocalHeap & lh)
    {
    
      mat = 1.0/sip.GetJacobiDet() *
	Trans (fel.GetDivShape(sip.IP(), lh));
    }


    template <typename FEL, typename SIP>
    static void GenerateMatrix (const FEL & fel, const SIP & sip,
				FlatVector<double> & mat, LocalHeap & lh)
    {
    
      mat = 1.0/sip.GetJacobiDet() *
        (fel.GetDivShape(sip.IP(), lh));
    }

    ///
    template <typename FEL, typename SIP, class TVX, class TVY>
    static void Apply (const FEL & fel, const SIP & sip,
		       const TVX & x, TVY & y,
		       LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCAL;
      
      Vec<DIM,TSCAL> hv = Trans (fel.GetDivShape(sip.IP(), lh)) * x;
      y = (1.0/sip.GetJacobiDet()) * hv;
    }

  ///
  template <typename FEL, typename SIP, class TVX, class TVY>
  static void ApplyTrans (const FEL & fel, const SIP & sip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh) 
  {
    typedef typename TVX::TSCAL TSCAL;
    Vec<DIM,TSCAL> hv = x;
    hv *= (1.0/sip.GetJacobiDet());
    y = fel.GetDivShape(sip.IP(),lh) * hv;
  }






  };


/// Identity on boundary
template <int D>
class DiffOpIdHDivBoundary : public DiffOp<DiffOpIdHDivBoundary<D> >
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = 1 };
  enum { DIFFORDER = 0 };

  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
     mat =  (1.0/sip.GetJacobiDet())*Trans(fel.GetShape (sip.IP(), lh));

  }

  template <typename FEL, typename SIP, class TVX, class TVY>
  static void Apply (const FEL & fel, const SIP & sip,
		     const TVX & x, TVY & y,
		     LocalHeap & lh)
  {

    y = (1.0/sip.GetJacobiDet())*(Trans (fel.GetShape (sip.IP(), lh)) * x);

  }

  template <typename FEL, typename SIP, class TVX, class TVY>
  static void ApplyTrans (const FEL & fel, const SIP & sip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh)
  {
    y = fel.GetShape (sip.IP(), lh)*((1.0/sip.GetJacobiDet())* x);
   
  }
};



  ///
  template <int D, typename FEL = HDivFiniteElement<D> >
  class MassHDivIntegrator
    : public T_BDBIntegrator<DiffOpIdHDiv<D>, DiagDMat<D>, FEL>
  {
  public:
    ///
    MassHDivIntegrator (CoefficientFunction * coeff)
      : T_BDBIntegrator<DiffOpIdHDiv<D>, DiagDMat<D>, FEL> (DiagDMat<D> (coeff))
    { (*testout)<<"MassHDivIntegrator"<<endl; }
    
    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new MassHDivIntegrator (coeffs[0]);
    }

    ///
    virtual string Name () const { return "MassHDiv"; }
  };

  

  ///
  template <int D, typename FEL = HDivFiniteElement<D> >
  class DivDivHDivIntegrator
    : public T_BDBIntegrator<DiffOpDivHDiv<D>, DiagDMat<1>, FEL>
  {
  public:
    ///
    DivDivHDivIntegrator (CoefficientFunction * coeff)
      : T_BDBIntegrator<DiffOpDivHDiv<D>, DiagDMat<1>, FEL> (DiagDMat<1> (coeff))
    { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new DivDivHDivIntegrator (coeffs[0]);
    }
    
    ///
    virtual string Name () const { return "DivDivHDiv"; }
  };




  ///
  template <int D, typename FEL = HDivFiniteElement<D> >
  class DivSourceHDivIntegrator 
    : public T_BIntegrator<DiffOpDivHDiv<D>, DVec<1>, FEL>
  {
  public:
    ///
    DivSourceHDivIntegrator (CoefficientFunction * coeff)
      : T_BIntegrator<DiffOpDivHDiv<D>, DVec<1>, FEL> (DVec<1> (coeff))
    { ; }
    
    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new DivSourceHDivIntegrator (coeffs[0]);
    }

    ///
    virtual string Name () const { return "DivSourceHDiv"; }
  };





  ///
  template <int D, typename FEL = HDivFiniteElement<D> >
  class SourceHDivIntegrator 
    : public T_BIntegrator<DiffOpIdHDiv<D>, DVec<D>, FEL>
  {
  public:
    ///
    SourceHDivIntegrator (CoefficientFunction * coeff1,
			  CoefficientFunction * coeff2,
			  CoefficientFunction * coeff3)
      : T_BIntegrator<DiffOpIdHDiv<D>, DVec<D>, FEL> (DVec<D> (coeff1, coeff2,coeff3))
    { ; }
    SourceHDivIntegrator (CoefficientFunction * coeff1,
			  CoefficientFunction * coeff2)
      : T_BIntegrator<DiffOpIdHDiv<D>, DVec<D>, FEL> (DVec<D> (coeff1, coeff2))
    { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      if(D == 2)
	return new SourceHDivIntegrator (coeffs[0], coeffs[1]);
      else if (D == 3)
	return new SourceHDivIntegrator (coeffs[0], coeffs[1], coeffs[2]);
    }
    
    ///
    virtual string Name () const { return "SourceHDiv"; }
  };

//***************************************Ali***************************************************************
///
template <int D, typename FEL = HDivNormalFiniteElement<D-1> >
class NeumannHDivIntegrator
  : public T_BIntegrator<DiffOpIdHDivBoundary<D>, DVec<1>, FEL>
{
public:
  ///
  NeumannHDivIntegrator (CoefficientFunction * coeff)
    : T_BIntegrator<DiffOpIdHDivBoundary<D>, DVec<1>, FEL> (DVec<1> (coeff))
  { ; }

  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new NeumannHDivIntegrator (coeffs[0]);
  }

  ///
  virtual bool BoundaryForm () const { return 1; }
  ///
  virtual string Name () const { return "NeumannHDiv"; }
};


/// integrator for $\int_\Gamma \sigma_n \tau_n \, ds$
template <int D>
class RobinHDivIntegrator
  : public T_BDBIntegrator<DiffOpIdHDivBoundary<D>, DiagDMat<1>, HDivNormalFiniteElement<D-1> >
{
public:
  RobinHDivIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpIdHDivBoundary<D>, DiagDMat<1>, HDivNormalFiniteElement<D-1> > (DiagDMat<1> (coeff))
  { ; }


  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new RobinHDivIntegrator (coeffs[0]);

  }

  virtual bool BoundaryForm () const { return 1; }
  virtual string Name () const { return "RobinHDiv"; }
};
#endif



