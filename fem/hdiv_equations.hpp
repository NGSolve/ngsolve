#ifndef FILE_HDIV_EQUATIONS
#define FILE_HDIV_EQUATIONS

/*********************************************************************/
/* File:   hdiv_equations.hpp                                        */
/* Author: Joachim Schoeberl, Almedin Becirovic                      */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/

namespace ngfem
{


/*
  Finite Element Integrators for H(div)

  Mapping with Piola transformation

  Requires H(div) finite elements
*/



/// Identity operator, Piola transformation
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


  template <typename FEL, class MIR>
  static void ApplyIR (const FEL & fel, const MIR & mir,
		       const FlatVector<double> & x, FlatMatrix<double> & y,
		       LocalHeap & lh)
  {
    fel.Evaluate (mir.IR(), x, FlatMatrixFixWidth<D> (y.Height(), &y(0,0)));
    for (int i = 0; i < mir.Size(); i++)
      {
	Vec<D> hy = (1.0/mir[i].GetJacobiDet()) * mir[i].GetJacobian() * y.Row(i);
	y.Row(i) = hy;
      }
  }



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




/// divergence Operator
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
    mat = 1.0/sip.GetJacobiDet() * (fel.GetDivShape(sip.IP(), lh));
  }


  template <typename FEL, typename SIP, class TVX, class TVY>
  static void Apply (const FEL & fel, const SIP & sip,
                     const TVX & x, TVY & y,
                     LocalHeap & lh) 
  {
    typedef typename TVX::TSCAL TSCAL;
      
    Vec<DIM,TSCAL> hv = Trans (fel.GetDivShape(sip.IP(), lh)) * x;
    y = (1.0/sip.GetJacobiDet()) * hv;
  }

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




/// Identity for boundary-normal elements
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


/// Identity for boundary-normal elements, gives q_n n
template <int D>
class DiffOpIdVecHDivBoundary : public DiffOp<DiffOpIdVecHDivBoundary<D> >
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 0 };

  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
    mat =  (1.0/sip.GetJacobiDet())*Trans(fel.GetShape (sip.IP(), lh)) * Trans (sip.GetNV());;
  }

  template <typename FEL, typename SIP, class TVX, class TVY>
  static void Apply (const FEL & fel, const SIP & sip,
		     const TVX & x, TVY & y,
		     LocalHeap & lh)
  {
    y = ( (1.0/sip.GetJacobiDet())*(InnerProduct (fel.GetShape (sip.IP(), lh), x) )) * sip.GetNV();
  }

  template <typename FEL, typename SIP, class TVX, class TVY>
  static void ApplyTrans (const FEL & fel, const SIP & sip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh)
  {
    y = ((1.0/sip.GetJacobiDet())* InnerProduct (x, sip.GetNV()) ) * fel.GetShape (sip.IP(), lh);
  }
};




/// Integrator for term of zero-th order
template <int D>
class MassHDivIntegrator
  : public T_BDBIntegrator<DiffOpIdHDiv<D>, DiagDMat<D>, HDivFiniteElement<D> >
{
public:
  MassHDivIntegrator (CoefficientFunction * coeff);

  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new MassHDivIntegrator (coeffs[0]);
  }

  virtual string Name () const { return "MassHDiv"; }
};

  

/// Integrator for div u div v
template <int D>
class DivDivHDivIntegrator
  : public T_BDBIntegrator<DiffOpDivHDiv<D>, DiagDMat<1>, HDivFiniteElement<D> >
{
public:
  DivDivHDivIntegrator (CoefficientFunction * coeff);

  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new DivDivHDivIntegrator (coeffs[0]);
  }
    
  virtual string Name () const { return "DivDivHDiv"; }
};




/// source term integrator for \ff div v
template <int D>
class DivSourceHDivIntegrator 
  : public T_BIntegrator<DiffOpDivHDiv<D>, DVec<1>, HDivFiniteElement<D> >
{
public:
  DivSourceHDivIntegrator (CoefficientFunction * coeff);

  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new DivSourceHDivIntegrator (coeffs[0]);
  }

  virtual string Name () const { return "DivSourceHDiv"; }
};

/*
/// source term for H(div)
template <int D>
class SourceHDivIntegrator 
  : public T_BIntegrator<DiffOpIdHDiv<D>, DVec<D>, HDivFiniteElement<D> >
{
public:
  SourceHDivIntegrator (CoefficientFunction * coeff1,
                        CoefficientFunction * coeff2,
                        CoefficientFunction * coeff3);

  SourceHDivIntegrator (CoefficientFunction * coeff1,
                        CoefficientFunction * coeff2);

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
*/


  template <int D> class SourceHDivIntegrator;


template <int D>
class BaseSourceHDivIntegrator 
  : public T_BIntegrator<DiffOpIdHDiv<D>, DVec<D>, HDivFiniteElement<D> >
{
public:
  BaseSourceHDivIntegrator(const DVec<D> & dvec)
    : T_BIntegrator<DiffOpIdHDiv<D>, DVec<D>, HDivFiniteElement<D> > (dvec) { ; }

  static Integrator * Create (Array<CoefficientFunction*> & coeffs);
  /*
  {
    if(D == 2)
      return new SourceHDivIntegrator<2> (coeffs[0], coeffs[1]);
    else if (D == 3)
      return new SourceHDivIntegrator<3> (coeffs[0], coeffs[1], coeffs[2]);
  }
  */
  ///
  virtual string Name () const { return "SourceHDiv"; }
};



template <>
class SourceHDivIntegrator<2>
  : public BaseSourceHDivIntegrator<2> 
{
public:
  SourceHDivIntegrator (CoefficientFunction * coeff1,
                        CoefficientFunction * coeff2);
};

template <>
class SourceHDivIntegrator<3>
  : public BaseSourceHDivIntegrator<3> 
{
public:
  SourceHDivIntegrator (CoefficientFunction * coeff1,
                        CoefficientFunction * coeff2,
                        CoefficientFunction * coeff3);
};









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


/// integrator for \f$\int_\Gamma \sigma_n \tau_n \, ds\f$
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


}


#endif



