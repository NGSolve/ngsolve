#ifndef FILE_BDBEQUATIONS
#define FILE_BDBEQUATIONS

/*********************************************************************/
/* File:   bdbequations.hpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngfem
{

/* 
   realizations of bdb integrators for many equations.
   The differential operators provide the B-matrix,
   the DMatOps provide the coefficient tensors
*/


/// Gradient operator of dimension D
template <int D>
class DiffOpGradient : public DiffOp<DiffOpGradient<D> >
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 1 };


  ///
  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
    mat = Trans (sip.GetJacobianInverse ()) * 
      Trans (fel.GetDShape(sip.IP(),lh));
  }

  template <typename FEL>
  static void GenerateMatrix (const FEL & fel, 
                              const SpecificIntegrationPoint<D,D> & sip,
			      FlatMatrixFixHeight<D> & mat, LocalHeap & lh)
  {
    fel.CalcMappedDShape (sip, mat.Trans());
  }

  ///
  template <typename FEL, typename SIP, class TVX, class TVY>
  static void Apply (const FEL & fel, const SIP & sip,
		     const TVX & x, TVY & y,
		     LocalHeap & lh) 
  {
    typedef typename TVX::TSCAL TSCAL;

    Vec<D,TSCAL> hv = Trans (fel.GetDShape(sip.IP(), lh)) * x;
    y = Trans (sip.GetJacobianInverse()) * hv;
  }

  template <typename FEL, class MIR>
  static void ApplyIR (const FEL & fel, const MIR & mir,
		       const FlatVector<double> & x, FlatMatrix<double> & y,
		       LocalHeap & lh)
  {
    FlatMatrixFixWidth<D> grad(mir.Size(), &y(0));
    fel.EvaluateGrad (mir.IR(), x, grad);
    for  (int i = 0; i < mir.Size(); i++)
      {
	Vec<D> hv = grad.Row(i);
	grad.Row(i) = Trans (mir[i].GetJacobianInverse()) * hv;
      }
  }


  ///
  template <typename FEL, typename SIP, class TVX, class TVY>
  static void ApplyTrans (const FEL & fel, const SIP & sip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh) 
  {
    typedef typename TVX::TSCAL TSCAL;

    Vec<D,TSCAL> hv = sip.GetJacobianInverse() * x;
    y = fel.GetDShape(sip.IP(),lh) * hv;
  }


  ///
  template <typename SIP, class TVX>
  static void Transform (const SIP & sip, TVX & x)
  {
    Vec<D> hx = Trans (sip.GetJacobianInverse()) * x;
    x = hx; 
  }

  ///
  template <typename SIP, class TVX>
  static void TransformT (const SIP & sip, TVX & x)
  {
    Vec<D> hx = sip.GetJacobianInverse() * x;
    x = hx; 
  }
  
  ///
  template <typename FEL, class TVD, class TVY>
  static void ApplyGrid (const FEL & fel, 
			 const IntegrationRuleTP<D> & ir,
			 const TVY & hv, 
			 TVD & dvecs, 
			 LocalHeap & lh)
  {
    FlatMatrix<double> dmat(dvecs.Size(), D, const_cast<double*> (&dvecs[0](0)));
    fel.EvaluateDShapeGrid (ir, hv, dmat, lh);
  }

  ///
  template <typename FEL, class TVD, class TVY>
  static void ApplyTransGrid (const FEL & fel, 
			      const IntegrationRuleTP<D> & ir,
			      const TVD & dvecs, 
			      TVY & hv, LocalHeap & lh)
  {
    FlatMatrix<double> dmat(dvecs.Size(), D, const_cast<double*> (&dvecs[0](0)));
    fel.EvaluateDShapeGridTrans (ir, dmat, hv, lh);
  }

};




/// Boundary Gradient operator of dimension D
template <int D>
class DiffOpGradientBoundary : public DiffOp<DiffOpGradientBoundary<D> >
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 1 };

  ///
  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
    mat = Trans (sip.GetJacobianInverse ()) * 
      Trans (fel.GetDShape(sip.IP(),lh));
  }
};





/// Gradient operator in r-z coordinates
template <int D>
class DiffOpGradientRotSym : 
  public DiffOp<DiffOpGradientRotSym<D> >
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 1 };

  ///
  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
    typedef typename MAT::TSCAL TSCAL;

    mat =  Trans (sip.GetJacobianInverse ()) * 
      Trans (fel.GetDShape(sip.IP(),lh));

    int i;
    double cx = sip.GetPoint()(0);
    if (cx == 0) cx = 1e-10;
    for (int i = 0; i < mat.Width(); i++)
      mat(0,i) += fel.GetShape(sip.IP(), lh)(i) / cx;

    // do the rot
    for (int i = 0; i < mat.Width(); i++)
      {
	TSCAL hv = mat(0,i);
	mat(0,i) = mat(1,i);
	mat(1,i) = -hv;
      }
  }
};



/// Identity
template <int D>
class DiffOpId : public DiffOp<DiffOpId<D> >
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D };
  enum { DIM_DMAT = 1 };
  enum { DIFFORDER = 0 };

  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
    mat.Row(0) = fel.GetShape(sip.IP(), lh);
  }

  static void GenerateMatrix (const ScalarFiniteElement<D> & fel, 
			      const SpecificIntegrationPoint<D,D> & sip,
			      FlatMatrixFixHeight<1> & mat, LocalHeap & lh)
  {
    fel.CalcShape (sip.IP(), FlatVector<> (fel.GetNDof(), &mat(0,0)));
  }


  template <typename FEL, typename SIP, class TVX, class TVY>
  static void Apply (const FEL & fel, const SIP & sip,
		     const TVX & x, TVY & y,
		     LocalHeap & lh) 
  {
    y = Trans (fel.GetShape (sip.IP(), lh)) * x;
  }

  static void Apply (const ScalarFiniteElement<D> & fel, const SpecificIntegrationPoint<D,D> & sip,
		     const FlatVector<double> & x, FlatVector<double> & y,
		     LocalHeap & lh) 
  {
    y(0) = fel.Evaluate(sip.IP(), x);
  }

  template <typename FEL, class MIR>
  static void ApplyIR (const FEL & fel, const MIR & mir,
		       const FlatVector<double> & x, FlatMatrix<double> & y,
		       LocalHeap & lh)
  {
    fel.Evaluate (mir.IR(), x, FlatVector<> (mir.Size(), &y(0,0)));
  }



  template <typename FEL, typename SIP, class TVX, class TVY>
  static void ApplyTrans (const FEL & fel, const SIP & sip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh) 
  {
    y = fel.GetShape (sip.IP(), lh) * x;
  }



  template <typename SIP, class TVX>
  static void Transform (const SIP & sip, TVX & x)
  {
    // do nothing
    ; 
  }

  template <typename SIP, class TVX>
  static void TransformT (const SIP & sip, TVX & x)
  {
    // do nothing
    ; 
  }


  template <typename FEL, class TVD, class TVY>
  static void ApplyGrid (const FEL & fel, 
			 const IntegrationRuleTP<D> & ir,
			 const TVY & hv, 
			 TVD & dvecs, 
			 LocalHeap & lh)
  {
    FlatVector<double> dvec(dvecs.Size(), const_cast<double*> (&dvecs[0](0)));
    fel.EvaluateShapeGrid (ir, hv, dvec, lh);
  }


  template <typename FEL, class TVD, class TVY>
  static void ApplyTransGrid (const FEL & fel, 
			      const IntegrationRuleTP<D> & ir,
			      const TVD & dvecs, 
			      TVY & hv, LocalHeap & lh)
  {
    FlatVector<double> dvec(dvecs.Size(), const_cast<double*> (&dvecs[0](0)));
    fel.EvaluateShapeGridTrans (ir, dvec, hv, lh);
  }
};



/// Identity
template <int D, int SYSDIM>
class DiffOpIdSys : public DiffOp<DiffOpIdSys<D,SYSDIM> >
{
public:
  enum { DIM = SYSDIM };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D };
  enum { DIM_DMAT = SYSDIM };
  enum { DIFFORDER = 0 };

  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
    const FlatVector<> shape = fel.GetShape (sip.IP(), lh);
  

    typedef typename MAT::TSCAL TSCAL;
    mat = TSCAL(0.); 
    for (int j = 0; j < shape.Height(); j++)
      for (int i = 0; i < SYSDIM; i++)
	mat(i, j*SYSDIM+i) = shape(j);
  }
};





/// Identity on boundary
template <int D>
class DiffOpIdBoundary : public DiffOp<DiffOpIdBoundary<D> >
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
    mat.Row(0) = fel.GetShape(sip.IP(), lh);
  }

  template <typename FEL, typename SIP, class TVX, class TVY>
  static void Apply (const FEL & fel, const SIP & sip,
		     const TVX & x, TVY & y,
		     LocalHeap & lh) 
  {
    y = Trans (fel.GetShape (sip.IP(), lh)) * x;
  }

  static void Apply (const ScalarFiniteElement<D-1> & fel, const SpecificIntegrationPoint<D-1,D> & sip,
		     const FlatVector<double> & x, FlatVector<double> & y,
		     LocalHeap & lh) 
  {
    y(0) = fel.Evaluate(sip.IP(), x);
  }


  template <typename FEL, typename SIP, class TVX, class TVY>
  static void ApplyTrans (const FEL & fel, const SIP & sip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh) 
  {
    y = fel.GetShape (sip.IP(), lh) * x;
  }



  template <typename SIP, class TVX>
  static void Transform (const SIP & sip, TVX & x)
  {
    // do nothing
    ; 
  }

  template <typename SIP, class TVX>
  static void TransformT (const SIP & sip, TVX & x)
  {
    // do nothing
    ; 
  }


  template <typename FEL, class TVD, class TVY>
  static void ApplyGrid (const FEL & fel, 
			 const IntegrationRuleTP<D-1> & ir,
			 const TVY & hv, 
			 TVD & dvecs, 
			 LocalHeap & lh)
  {
    void * heapp = lh.GetPointer();
    FlatVector<double> dvec(dvecs.Size(), const_cast<double*> (&dvecs[0](0)));
    fel.EvaluateShapeGrid (ir, hv, dvec, lh);
  }


  template <typename FEL, class TVD, class TVY>
  static void ApplyTransGrid (const FEL & fel, 
			      const IntegrationRuleTP<D-1> & ir,
			      const TVD & dvecs, 
			      TVY & hv, LocalHeap & lh)
  {
    void * heapp = lh.GetPointer();
    FlatVector<double> dvec(dvecs.Size(), const_cast<double*> (&dvecs[0](0)));
    fel.EvaluateShapeGridTrans (ir, dvec, hv, lh);
  }

};



/// Identity
template <int D, int SYSDIM>
class DiffOpIdBoundarySys : public DiffOp<DiffOpIdBoundarySys<D,SYSDIM> >
{
public:
  enum { DIM = SYSDIM };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = SYSDIM };
  enum { DIFFORDER = 0 };

  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
    mat = 0.; 
    const FlatVector<> shape = fel.GetShape (sip.IP(), lh);
    for (int j = 0; j < shape.Height(); j++)
      for (int i = 0; i < SYSDIM; i++)
	mat(i, j*SYSDIM+i) = shape(j);
  }
};























/*
/// diagonal tensor, all values are the same
template <int DIM, typename SCAL = double>
class DiagDMat : public DMatOp<DiagDMat<DIM,SCAL> >
{
  CoefficientFunction * coef;
public:
  typedef SCAL TSCAL;
  enum { DIM_DMAT = DIM };
  DiagDMat (CoefficientFunction * acoef) : coef(acoef) { ; }

  DiagDMat (Array<CoefficientFunction*> & acoefs) : coef(acoefs[0]) { ; }

  template <typename FEL, typename SIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const SIP & sip,
		       MAT & mat, LocalHeap & lh) const
  {
    typedef typename MAT::TSCAL TRESULT;
    mat = TRESULT(0);
    TSCAL val = coef -> T_Evaluate<TSCAL> (sip);
    for (int i = 0; i < DIM; i++)
      mat(i, i) = ConvertTo<TRESULT> (val);
  }  

  template <typename FEL, class VECX, class VECY>
  void Apply (const FEL & fel, const BaseSpecificIntegrationPoint & sip,
	      const VECX & x, VECY & y, LocalHeap & lh) const
  {
    typedef typename VECY::TSCAL TRESULT;
    TSCAL val = coef -> T_Evaluate<TSCAL> (sip);
    for (int i = 0; i < DIM; i++)
      y(i) = ConvertTo<TRESULT> (val * x(i));
  }

  template <typename FEL, class VECX>
  void Apply1 (const FEL & fel, const BaseSpecificIntegrationPoint & sip,
	      const VECX & x, LocalHeap & lh) const
  {
    TSCAL val = coef -> T_Evaluate<TSCAL> (sip);
    x *= ConvertTo<typename VECX::TSCAL> (val);
  }

  template <typename FEL, typename MIR, typename TVX>
  void ApplyIR (const FEL & fel, const MIR & mir,
		TVX & x, LocalHeap & lh) const
  {
    FlatMatrix<TSCAL> values(mir.Size(), 1, lh);
    coef -> Evaluate (mir, values);
    for (int i = 0; i < mir.Size(); i++)
      x.Row(i) *=  ConvertTo<typename TVX::TSCAL>  (values(i, 0));
  }
};
*/

/// diagonal tensor, all values are the same
template <int DIM>
class DiagDMat : public DMatOp<DiagDMat<DIM> >
{
  CoefficientFunction * coef;
public:
  // typedef SCAL TSCAL;
  enum { DIM_DMAT = DIM };
  DiagDMat (CoefficientFunction * acoef) : coef(acoef) { ; }

  DiagDMat (Array<CoefficientFunction*> & acoefs) : coef(acoefs[0]) { ; }

  template <typename FEL, typename SIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const SIP & sip,
		       MAT & mat, LocalHeap & lh) const
  {
    typedef typename MAT::TSCAL TRESULT;
    mat = TRESULT(0);
    TRESULT val = coef -> T_Evaluate<TRESULT> (sip);
    for (int i = 0; i < DIM; i++)
      mat(i, i) = val;
  }  

  template <typename FEL, class VECX, class VECY>
  void Apply (const FEL & fel, const BaseSpecificIntegrationPoint & sip,
	      const VECX & x, VECY & y, LocalHeap & lh) const
  {
    typedef typename VECY::TSCAL TRESULT;
    TRESULT val = coef -> T_Evaluate<TRESULT> (sip);
    for (int i = 0; i < DIM; i++)
      y(i) = val * x(i);
  }

  template <typename FEL, class VECX>
  void Apply1 (const FEL & fel, const BaseSpecificIntegrationPoint & sip,
	       const VECX & x, LocalHeap & lh) const
  {
    typedef typename VECX::TSCAL TSCAL;
    x *= coef -> T_Evaluate<TSCAL> (sip);
  }

  template <typename FEL, typename MIR, typename TVX>
  void ApplyIR (const FEL & fel, const MIR & mir,
		TVX & x, LocalHeap & lh) const
  {
    typedef typename TVX::TSCAL TSCAL;
    FlatMatrix<TSCAL> values(mir.Size(), 1, lh);
    coef -> Evaluate (mir, values);
    for (int i = 0; i < mir.Size(); i++)
      x.Row(i) *=  values(i, 0);
  }
};








/// orthotropic tensor. 
template <int N> 
class OrthoDMat
{
};

template <> class OrthoDMat<1> : public DMatOp<OrthoDMat<1> >
{
  CoefficientFunction * coef;
public:
  enum { DIM_DMAT = 1 };
  OrthoDMat (CoefficientFunction * acoef) : coef(acoef) { ; }

  template <typename FEL, typename SIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const SIP & sip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat(0,0) = Evaluate (*coef, sip);
  }  

  template <typename FEL, typename SIP, class VECX, class VECY>
  void Apply (const FEL & fel, const SIP & sip,
	      const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef, sip) * x(0);
  }

  template <typename FEL, typename SIP, class VECX, class VECY>
  void ApplyTrans (const FEL & fel, const SIP & sip,
		   const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef, sip) * x(0);
  }

};


template <> class OrthoDMat<2>: public DMatOp<OrthoDMat<2> >
{
  CoefficientFunction * coef1;
  CoefficientFunction * coef2;
public:
  enum { DIM_DMAT = 2 };
  OrthoDMat (CoefficientFunction * acoef1,
	CoefficientFunction * acoef2)
    : coef1(acoef1), coef2(acoef2) { ; }
  OrthoDMat (CoefficientFunction * acoef1,
	CoefficientFunction * acoef2,
	CoefficientFunction * acoef3)
    : coef1(acoef1), coef2(acoef2) { ; }
  
  template <typename FEL, typename SIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const SIP & sip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    mat(0,0) = Evaluate (*coef1, sip);
    mat(1,1) = Evaluate (*coef2, sip);
  }  

  template <typename FEL, typename SIP>
  void GetEigensystem (const FEL & fel, const SIP & sip, 
		  Array<double> & eigenvalues,
		  Array<double> & eigenvectors,
		  LocalHeap & lh) const
  {
    eigenvalues[0] = Evaluate (*coef1, sip);
    eigenvalues[1] = Evaluate (*coef2, sip);
    eigenvectors[0] = eigenvectors[3] = 1.;
    eigenvectors[1] = eigenvectors[2] = 0.;
  }

  template <typename FEL, typename SIP, class VECX, class VECY>
  void Apply (const FEL & fel, const SIP & sip,
	      const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef1, sip) * x(0);
    y(1) = Evaluate (*coef2, sip) * x(1);
  }

  template <typename FEL, typename SIP, class VECX, class VECY>
  void ApplyTrans (const FEL & fel, const SIP & sip,
		   const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef1, sip) * x(0);
    y(1) = Evaluate (*coef2, sip) * x(1);
  }
};

template <> class OrthoDMat<3> : public DMatOp<OrthoDMat<3> >
{
  CoefficientFunction * coef1;
  CoefficientFunction * coef2;
  CoefficientFunction * coef3;
public:
  enum { DIM_DMAT = 3 };
  OrthoDMat (CoefficientFunction * acoef1,
	     CoefficientFunction * acoef2)
    : coef1(acoef1), coef2(acoef2), coef3(acoef2) { ; }
  OrthoDMat (CoefficientFunction * acoef1,
	     CoefficientFunction * acoef2,
	     CoefficientFunction * acoef3)
    : coef1(acoef1), coef2(acoef2), coef3(acoef3) { ; }
  
  template <typename FEL, typename SIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const SIP & sip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    mat(0,0) = Evaluate (*coef1, sip);
    mat(1,1) = Evaluate (*coef2, sip);
    mat(2,2) = Evaluate (*coef3, sip);
  }  

  
  template <typename FEL, typename SIP>
  void GetEigensystem (const FEL & fel, const SIP & sip, 
		  Array<double> & eigenvalues,
		  Array<double> & eigenvectors,
		  LocalHeap & lh) const
  {
    
    eigenvalues[0] = Evaluate(*coef1,sip);
    eigenvalues[1] = Evaluate(*coef2,sip);
    eigenvalues[2] = Evaluate(*coef3,sip);

    eigenvectors = 0.;
    eigenvectors[0] = eigenvectors[4] = eigenvectors[8] = 1.;
  }


  template <typename FEL, typename SIP, class VECX, class VECY>
  void Apply (const FEL & fel, const SIP & sip,
	      const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef1, sip) * x(0);
    y(1) = Evaluate (*coef2, sip) * x(1);
    y(2) = Evaluate (*coef3, sip) * x(2);
  }

  template <typename FEL, typename SIP, class VECX, class VECY>
  void ApplyTrans (const FEL & fel, const SIP & sip,
		   const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef1, sip) * x(0);
    y(1) = Evaluate (*coef2, sip) * x(1);
    y(2) = Evaluate (*coef3, sip) * x(2);
  }

  void SetCoefficientFunctions( CoefficientFunction * acoef1,
				CoefficientFunction * acoef2,
				CoefficientFunction * acoef3 )
  {
    // NOTE: alte coefficient-functions werden nicht geloescht!
    coef1 = acoef1;
    coef2 = acoef2;
    coef3 = acoef3;
  }
};









/// full symmetric tensor
template <int N> 
class SymDMat : public DMatOp<SymDMat<N> >
{
};

template <> class SymDMat<1> : public DMatOp<SymDMat<1> >
{
  CoefficientFunction * coef;
public:
  enum { DIM_DMAT = 1 };

  SymDMat (CoefficientFunction * acoef) : coef(acoef) { ; }

  template <typename FEL, typename SIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const SIP & sip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat(0,0) = Evaluate (*coef, sip);
  }  
};


template <> class SymDMat<2> : public DMatOp<SymDMat<2> >
{
  CoefficientFunction * coef00;
  CoefficientFunction * coef01;
  CoefficientFunction * coef11;
public:
  enum { DIM_DMAT = 2 };

  SymDMat (CoefficientFunction * acoef00,
	   CoefficientFunction * acoef01,
	   CoefficientFunction * acoef11)
    : coef00(acoef00), coef01(acoef01), coef11(acoef11) { ; }
  
  template <typename FEL, typename SIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const SIP & sip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    mat(0,0) = Evaluate (*coef00, sip);
    mat(0,1) = mat(1,0) = Evaluate (*coef01, sip);
    mat(1,1) = Evaluate (*coef11, sip);
  }  
};

template <> class SymDMat<3> : public DMatOp<SymDMat<3> >
{
  CoefficientFunction * coef00;
  CoefficientFunction * coef10;
  CoefficientFunction * coef11;
  CoefficientFunction * coef20;
  CoefficientFunction * coef21;
  CoefficientFunction * coef22;
public:
  enum { DIM_DMAT = 3 };

  SymDMat (CoefficientFunction * acoef00,
	   CoefficientFunction * acoef10,
	   CoefficientFunction * acoef11,
	   CoefficientFunction * acoef20,
	   CoefficientFunction * acoef21,
	   CoefficientFunction * acoef22)
    : coef00(acoef00), coef10(acoef10), coef11(acoef11),
      coef20(acoef20), coef21(acoef21), coef22(acoef22) { ; }
  
  template <typename FEL, typename SIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const SIP & sip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    mat(0,0) = Evaluate (*coef00, sip);
    mat(1,0) = mat(0,1) = Evaluate (*coef10, sip);
    mat(1,1) = Evaluate (*coef11, sip);
    mat(2,0) = mat(0,2) = Evaluate (*coef20, sip);
    mat(2,1) = mat(1,2) = Evaluate (*coef21, sip);
    mat(2,2) = Evaluate (*coef22, sip);
  }  
};









///
template <int DIM>
class NormalDMat : public DMatOp<NormalDMat<DIM> >
{
  CoefficientFunction * coef;
public:
  enum { DIM_DMAT = DIM };
  NormalDMat (CoefficientFunction * acoef) : coef(acoef) { ; }

  template <typename FEL, typename SIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const SIP & sip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    double val = Evaluate (*coef, sip);
    Vec<DIM> nv = sip.GetNV();
    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
	mat(i, j) = val * nv(i) * nv(j);
  }  
};











template <int N, typename TSCAL> 
class DVecBase 
{
 protected:
  CoefficientFunction * coefs[N];
 public:
  template <typename FEL, typename SIP, typename VEC>
  void GenerateVector (const FEL & fel, const SIP & sip,
		       VEC & vec, LocalHeap & lh) const
  {
    for (int i = 0; i < N; i++)
      vec(i) = coefs[i] -> Evaluate (sip);
  }  
};

template <int N>
class DVecBase<N, Complex>
{
 protected:
  CoefficientFunction * coefs[N];
 public:
  template <typename FEL, typename SIP, typename VEC>
  void GenerateVector (const FEL & fel, const SIP & sip,
		       VEC & vec, LocalHeap & lh) const
  {
    typedef typename VEC::TSCAL TRESULT;
    
    for (int i = 0; i < N; i++)
      vec(i) = ConvertTo<TRESULT> (coefs[i] -> EvaluateComplex (sip));
  }  
};

// template <typename T = double> 
// class DVecVec
// {
//  protected:
//   CoefficientFunction * coefs;
//  public:
//   typedef T TSCAL; 
//   
//   DVecVec (CoefficientFunction * acoef)
//   { 
//     coefs = acoef;
//   }
//   
//   template <typename FEL, typename SIP, typename VEC>
//   void GenerateVector (const FEL & fel, const SIP & sip,
// 		       VEC & vec, LocalHeap & lh) const
//   {
//       coefs -> Evaluate (sip,vec);
//   }  
// };


template <int N, typename T = double> 
class NGS_DLL_HEADER DVec  : public DVecBase<N,T>
{
  using DVecBase<N,T>::coefs;
  typedef T TSCAL;

public:
  DVec (Array<CoefficientFunction*> & acoeffs)
  {
    for (int i = 0; i < N; i++)
      coefs[i] = acoeffs[i];
  }
  
  DVec (CoefficientFunction * acoef)
  { 
    coefs[0] = acoef;
  }

  DVec (CoefficientFunction * acoef1,
	CoefficientFunction * acoef2)
  { 
    coefs[0] = acoef1;
    coefs[1] = acoef2;
  }

  DVec (CoefficientFunction * acoef1,
	CoefficientFunction * acoef2,
	CoefficientFunction * acoef3)
  { 
    coefs[0] = acoef1;
    coefs[1] = acoef2;
    coefs[2] = acoef3;
  }
};



  /*
template <typename T> 
class DVec<1, T> : public DVecBase<1,T>
{
public:
  using DVecBase<1,T>::coefs;
  typedef T TSCAL;

  DVec (CoefficientFunction * acoef)
  { 
    coefs[0] = acoef;
  }
};

template <typename T> 
class DVec<2, T> : public DVecBase<2,T>
{
public:
  typedef T TSCAL;
  using DVecBase<2,T>::coefs;
  DVec (CoefficientFunction * acoef1,
	CoefficientFunction * acoef2)
  { 
    coefs[0] = acoef1;
    coefs[1] = acoef2;
  }
};

template <typename T> 
class DVec<3, T> : public DVecBase<3,T>
{
public:
  typedef T TSCAL;
  using DVecBase<3,T>::coefs;

  DVec (CoefficientFunction * acoef1,
	CoefficientFunction * acoef2,
	CoefficientFunction * acoef3)
  { 
    coefs[0] = acoef1;
    coefs[1] = acoef2;
    coefs[2] = acoef3;
  }
};
  */




///
		   /*
template <int N, typename TSCAL = double> 
class DVec { };

template <typename TSCAL> class DVec<1, TSCAL>
{
  CoefficientFunction * coef;
public:
  typedef double TSCAL;

  DVec (CoefficientFunction * acoef) : coef(acoef) { ; }

  template <typename FEL, typename SIP, typename VEC>
  void GenerateVector (const FEL & fel, const SIP & sip,
		       VEC & vec, LocalHeap & lh) const
  {
    vec(0) = coef -> Evaluate (sip);
  }  
};

template <> class DVec<1, Complex>
{
  CoefficientFunction * coef;
public:
  typedef Complex TSCAL;

  DVec (CoefficientFunction * acoef) : coef(acoef) { ; }

  template <typename FEL, typename SIP, typename VEC>
  void GenerateVector (const FEL & fel, const SIP & sip,
		       VEC & vec, LocalHeap & lh) const
  {
    vec(0) = coef -> EvaluateComplex (sip);
  }  
};


template <> class DVec<2, double>
{
  CoefficientFunction * coef1;
  CoefficientFunction * coef2;
public:
  typedef double TSCAL;

  DVec (CoefficientFunction * acoef1,
	CoefficientFunction * acoef2)
    : coef1(acoef1), coef2(acoef2) { ; }

  template <typename FEL, typename SIP, typename VEC>
  void GenerateVector (const FEL & fel, const SIP & sip,
		       VEC & vec, LocalHeap & lh) const
  {
    vec(0) = Evaluate (*coef1, sip);
    vec(1) = Evaluate (*coef2, sip);
  }  
};

template <> class DVec<3, double>
{
  CoefficientFunction * coef1;
  CoefficientFunction * coef2;
  CoefficientFunction * coef3;
public:
  typedef double TSCAL;

  DVec (CoefficientFunction * acoef1,
	CoefficientFunction * acoef2,
	CoefficientFunction * acoef3)
    : coef1(acoef1), coef2(acoef2), coef3(acoef3) { ; }
  
  template <typename FEL, typename SIP, typename VEC>
  void GenerateVector (const FEL & fel, const SIP & sip,
		       VEC & vec, LocalHeap & lh) const
  {
    vec(0) = Evaluate (*coef1, sip);
    vec(1) = Evaluate (*coef2, sip);
    vec(2) = Evaluate (*coef3, sip);
  }  
};
*/




template <> class DVec<6>
{
  CoefficientFunction * coef1;
  CoefficientFunction * coef2;
  CoefficientFunction * coef3;
  CoefficientFunction * coef4;
  CoefficientFunction * coef5;
  CoefficientFunction * coef6;
public:
  DVec (CoefficientFunction * acoef1,
	CoefficientFunction * acoef2,
	CoefficientFunction * acoef3,
	CoefficientFunction * acoef4,
	CoefficientFunction * acoef5,
	CoefficientFunction * acoef6)
    : coef1(acoef1), coef2(acoef2), coef3(acoef3),
      coef4(acoef4), coef5(acoef5), coef6(acoef6) { ; }
  
  template <typename FEL, typename SIP, typename VEC>
  void GenerateVector (const FEL & fel, const SIP & sip,
		       VEC & vec, LocalHeap & lh) const
  {
    vec(0) = Evaluate (*coef1, sip);
    vec(1) = Evaluate (*coef2, sip);
    vec(2) = Evaluate (*coef3, sip);
    vec(3) = Evaluate (*coef4, sip);
    vec(4) = Evaluate (*coef5, sip);
    vec(5) = Evaluate (*coef6, sip);
  }  
};





template <int N, typename T = double>  
class DVecN
{
  CoefficientFunction * coef;
public:
  typedef T TSCAL;
  DVecN (CoefficientFunction * acoef)
    : coef(acoef) { ; }
  
  template <typename FEL, typename SIP, typename VEC>
  void GenerateVector (const FEL & fel, const SIP & sip,
		       VEC & vec, LocalHeap & lh) const
  {
    Vec<N> hv;
    coef -> Evaluate (sip, hv);
    for (int i = 0; i < N; i++)
      vec(i) = hv(i);
  }  
};



template <int N>
class TVec
{
  CoefficientFunction * coef;

public:
  typedef double TSCAL;
  TVec (CoefficientFunction * acoef) : coef(acoef) {;}

  

  template <typename FEL, typename SIP, typename VEC>
  void GenerateVector (const FEL & fel, const SIP & sip,
		       VEC & vec, LocalHeap & lh) const
  {
    vec = 0;

    typedef typename VEC::TSCAL TSCAL;
    
    TSCAL length = 0.;
    for(int i=0; i<N; i++)
      {
	//vec(i) = sip.GetJacobian()(i,0);
	vec(i) = sip.GetTV()(i);
	length += vec(i)*vec(i);
      }
    //(*testout) << "point " << sip.GetPoint() << " tv " << vec;
    vec *= Evaluate (*coef, sip)/sqrt(length);
    //(*testout) << " retval " << vec << endl;
  }

};











/// DMat for rot.-sym. Laplace operator
template <int DIM>
class RotSymLaplaceDMat : public DMatOp<RotSymLaplaceDMat<DIM> >
{
  CoefficientFunction * coef;
public:
  enum { DIM_DMAT = DIM };

  RotSymLaplaceDMat (CoefficientFunction * acoef) : coef(acoef) { ; }

  template <typename FEL, typename SIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const SIP & sip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    const double r = sip.GetPoint()(0);
    double val = r*Evaluate (*coef, sip);
    for (int i = 0; i < DIM; i++)
      mat(i, i) = val;
  }  
  

  template <typename FEL, typename SIP, class VECX, class VECY>
  void Apply (const FEL & fel, const SIP & sip,
	      const VECX & x, VECY & y, LocalHeap & lh) const
  {
    const double r = sip.GetPoint()(0);
    double val = r*Evaluate (*coef, sip);
    y = val * x;
  }
};






/// Identity on boundary
template <int D>
class DiffOpNormal : public DiffOp<DiffOpNormal<D> >
{
public:
  enum { DIM = D };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = 1 };
  enum { DIFFORDER = 0 };

  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
    FlatVector<> shape = fel.GetShape (sip.IP(), lh);
    Vec<D> nv = sip.GetNV();
    //Vec<D> p = sip.GetPoint();
    for (int j = 0; j < shape.Size(); j++)
      for (int i = 0; i < D; i++)
	mat(0, j*D+i) =  shape(j) * nv(i);

    //(*testout) << "sip = " << p << ", nv = " << nv << endl;
    //p(0) = 0.0;
    //p /= L2Norm(p);
    //nv /= L2Norm(nv);
    //(*testout) << "normalized, sip = " << p << ", nv = " << nv << endl;
    //(*testout) << "mat = " << mat << endl;
  }
};





/// integrator for \f$\int_\Gamma u_n v_n \, ds\f$
template <int D>
class NormalRobinIntegrator 
  : public T_BDBIntegrator<DiffOpNormal<D>, DiagDMat<1>, ScalarFiniteElement<D-1> >
{
public:
  NormalRobinIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpNormal<D>, DiagDMat<1>, ScalarFiniteElement<D-1> > (DiagDMat<1> (coeff))
  { ; }


  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new NormalRobinIntegrator (coeffs[0]);
  }

  virtual bool BoundaryForm () const { return 1; }
  virtual string Name () const { return "NormalRobin"; }
};






// ********************************* Scalar integrators: ********************

/// Integrator for grad u grad v
template <int D, typename FEL = ScalarFiniteElement<D> >
class NGS_DLL_HEADER LaplaceIntegrator 
  : public T_BDBIntegrator<DiffOpGradient<D>, DiagDMat<D>, FEL>
{
public:
  ///
  LaplaceIntegrator (CoefficientFunction * coeff);
  LaplaceIntegrator (Array<CoefficientFunction*> & coeffs);

  /*
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new LaplaceIntegrator (coeffs[0]);
  }
  */
  virtual string Name () const { return "Laplace"; }
};



/// 
template <int D, typename FEL = ScalarFiniteElement<D-1> >
class LaplaceBoundaryIntegrator 
  : public T_BDBIntegrator<DiffOpGradientBoundary<D>, DiagDMat<D>, FEL>
{
public:
  ///
  LaplaceBoundaryIntegrator (CoefficientFunction * coeff);
  /*
    : T_BDBIntegrator<DiffOpGradientBoundary<D>, DiagDMat<D>, FEL> (DiagDMat<D> (coeff))
  {
    ;
  }
  */

  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new LaplaceBoundaryIntegrator (coeffs[0]);
  }

  ///
  virtual string Name () const { return "Laplace-Boundary"; }
};



///
template <int D, typename FEL = ScalarFiniteElement<D> >
class RotSymLaplaceIntegrator 
  : public T_BDBIntegrator<DiffOpGradient<D>, RotSymLaplaceDMat<D>, FEL>
{
public:
  ///
  RotSymLaplaceIntegrator (CoefficientFunction * coeff);
  /*
    : T_BDBIntegrator<DiffOpGradient<D>, RotSymLaplaceDMat<D>, FEL> (RotSymLaplaceDMat<D> (coeff))
  { ; }
  */
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new RotSymLaplaceIntegrator (coeffs[0]);
  }
  
  ///
  virtual string Name () const { return "RotSymLaplace"; }
};


///
template <int D, typename FEL = ScalarFiniteElement<D> >
class OrthoLaplaceIntegrator 
  : public T_BDBIntegrator<DiffOpGradient<D>, OrthoDMat<D>, FEL>
{
public:
  ///
  OrthoLaplaceIntegrator (CoefficientFunction * coeff1, CoefficientFunction * coeff2)
    : T_BDBIntegrator<DiffOpGradient<D>, OrthoDMat<D>, FEL> (OrthoDMat<D> (coeff1, coeff2))
  { ; }
  OrthoLaplaceIntegrator (CoefficientFunction * coeff1, CoefficientFunction * coeff2, CoefficientFunction * coeff3)
    : T_BDBIntegrator<DiffOpGradient<D>, OrthoDMat<D>, FEL> (OrthoDMat<D> (coeff1, coeff2, coeff3))
  { ; }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    if(coeffs.Size() == 2)
      return new OrthoLaplaceIntegrator (coeffs[0], coeffs[1]);
    else
      return new OrthoLaplaceIntegrator (coeffs[0], coeffs[1], coeffs[2]);
  }


  ///
  virtual string Name () const { return "OrthoLaplace"; }
};



///
template <int D>
class MassIntegrator 
  : public T_BDBIntegrator<DiffOpId<D>, DiagDMat<1>, ScalarFiniteElement<D> >
{
public:
  ///
  NGS_DLL_HEADER MassIntegrator (CoefficientFunction * coeff);
  NGS_DLL_HEADER MassIntegrator (Array<CoefficientFunction*> & coeffs);
    
  ///
  static Integrator * Create (Array<CoefficientFunction*> & coeffs);
  ///
  virtual string Name () const { return "Mass"; }
};


template <int D>
class ComplexMassIntegrator 
// : public T_BDBIntegrator<DiffOpId<D>, DiagDMat<1, Complex>, ScalarFiniteElement<D> >
  : public T_BDBIntegrator<DiffOpId<D>, DiagDMat<1>, ScalarFiniteElement<D> >
{
public:
  ///
  ComplexMassIntegrator (CoefficientFunction * coeff);
  ///
  virtual string Name () const { return "ComplexMass"; }
};





/// integrator for \f$\int_\Gamma u v \, ds\f$
template <int D>
class NGS_DLL_HEADER RobinIntegrator 
  : public T_BDBIntegrator<DiffOpIdBoundary<D>, DiagDMat<1>, ScalarFiniteElement<D-1> >
{
public:
  ///
  RobinIntegrator (CoefficientFunction * coeff);
  ///
  static Integrator * Create (Array<CoefficientFunction*> & coeffs);
  // {
  // return new RobinIntegrator (coeffs[0]);
  // }
  ///
  virtual bool BoundaryForm () const { return 1; }
  ///
  virtual string Name () const { return "Robin"; }
};

template <int D>
class ComplexRobinIntegrator 
// : public T_BDBIntegrator<DiffOpIdBoundary<D>, DiagDMat<1, Complex>, ScalarFiniteElement<D-1> >
  : public T_BDBIntegrator<DiffOpIdBoundary<D>, DiagDMat<1>, ScalarFiniteElement<D-1> >
{
public:
  ///
  ComplexRobinIntegrator (CoefficientFunction * coeff);
  ///
  virtual string Name () const { return "ComplexRobin"; }
};





/*
template <int D>
class NormalRobinIntegrator 
  : public T_BDBIntegrator<DiffOpIdBoundary<D,D>, NormalDMat<D>, ScalarFiniteElement<D> >
{
public:
  NormalRobinIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpIdBoundary<D,D>, NormalDMat<D>, ScalarFiniteElement<D> > (NormalDMat<D> (coeff))
  { ; }

  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new NormalRobinIntegrator (coeffs[0]);
  }

  virtual bool BoundaryForm () const { return 1; }
  virtual string Name () const { return "NormalRobin"; }
};
*/


/// Elasticity operator $(e_{11},e_{22},2 e_{12})$
template <int D> 
class DiffOpDiv : public DiffOp<DiffOpDiv<D> >
{
public:
  enum { DIM = D };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D };
  enum { DIM_DMAT = 1 };
  enum { DIFFORDER = 1 };

  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
    int nd = fel.GetNDof();

    FlatMatrix<> grad (D, nd, lh);
    grad = Trans (sip.GetJacobianInverse ()) * 
      Trans (fel.GetDShape(sip.IP(), lh));
    
    mat = 0;
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < DIM; j++)
	mat(0, DIM*i+j) = grad(j, i);
  }
};


template <int D, typename FEL = ScalarFiniteElement<D> >
class DivDivIntegrator 
  : public T_BDBIntegrator<DiffOpDiv<D>, DiagDMat<1>, FEL>
{
public:
  ///
  DivDivIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpDiv<D>, DiagDMat<1>, FEL> (DiagDMat<1> (coeff))
  {
    ;
  }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new DivDivIntegrator (coeffs[0]);
  }

  ///
  virtual string Name () const { return "DivDiv"; }
};












///
class DiffOpCurl : public DiffOp<DiffOpCurl>
{
public:
  enum { DIM = 2 };
  enum { DIM_SPACE = 2 };
  enum { DIM_ELEMENT = 2 };
  enum { DIM_DMAT = 1 };
  enum { DIFFORDER = 1 };

  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
    int nd = fel.GetNDof();

    FlatMatrix<> grad (2, nd, lh);
    grad = Trans (sip.GetJacobianInverse ()) * 
      Trans (fel.GetDShape(sip.IP(), lh));
    
    mat = 0;
    for (int i = 0; i < nd; i++)
      {
	mat(0, DIM*i  ) = grad(1, i);
	mat(0, DIM*i+1) = -grad(0, i);
      }
  }
};


template <typename FEL = ScalarFiniteElement<2> >
class CurlCurlIntegrator 
  : public T_BDBIntegrator<DiffOpCurl, DiagDMat<1>, FEL>
{
public:
  ///
  CurlCurlIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpCurl, DiagDMat<1>, FEL> (DiagDMat<1> (coeff))
  {
    ;
  }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new CurlCurlIntegrator (coeffs[0]);
  }

  ///
  virtual string Name () const { return "CurlCurl"; }
};




///
class DiffOpCurl3d : public DiffOp<DiffOpCurl3d>
{
public:
  enum { DIM = 3 };
  enum { DIM_SPACE = 3 };
  enum { DIM_ELEMENT = 3 };
  enum { DIM_DMAT = 3 };
  enum { DIFFORDER = 1 };

  template <typename FEL, typename SIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const SIP & sip,
			      MAT & mat, LocalHeap & lh)
  {
    int nd = fel.GetNDof();

    FlatMatrix<> grad (3, nd, lh);
    grad = Trans (sip.GetJacobianInverse ()) * 
      Trans (fel.GetDShape(sip.IP(), lh));
    
    mat = 0;
    for (int i = 0; i < nd; i++)
      {
	mat(0, DIM*i+2) =  grad(1, i);
	mat(0, DIM*i+1) = -grad(2, i);
	mat(1, DIM*i+0) =  grad(2, i);
	mat(1, DIM*i+2) = -grad(0, i);
	mat(2, DIM*i+1) =  grad(0, i);
	mat(2, DIM*i+0) = -grad(1, i);
      }
  }
};


template <typename FEL = ScalarFiniteElement<3> >
class CurlCurl3dIntegrator 
  : public T_BDBIntegrator<DiffOpCurl3d, DiagDMat<3>, FEL>
{
public:
  ///
  CurlCurl3dIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpCurl3d, DiagDMat<3>, FEL> (DiagDMat<3> (coeff))
  {
    ;
  }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new CurlCurl3dIntegrator (coeffs[0]);
  }

  ///
  virtual string Name () const { return "CurlCurl3d"; }
};

























/* ************************** Linearform Integrators ************************* */

/// integrator for \f$\int_\Omega f v \f$
template <int D, typename FEL = ScalarFiniteElement<D>  >
class NGS_DLL_HEADER SourceIntegrator 
  : public T_BIntegrator<DiffOpId<D>, DVec<1>, FEL>
{
public:
  SourceIntegrator (CoefficientFunction * coeff);
  SourceIntegrator (Array<CoefficientFunction*> & coeffs);
  // static Integrator * Create (Array<CoefficientFunction*> & coeffs);
  virtual string Name () const { return "Source"; }
};

template <int D, typename FEL = ScalarFiniteElement<D>  >
class ComplexSourceIntegrator 
  : public T_BIntegrator<DiffOpId<D>, DVec<1, Complex>, FEL>
{
public:
  ComplexSourceIntegrator (CoefficientFunction * coeff);
  virtual string Name () const { return "ComplexSource"; }
};



///
template <int D, typename FEL = ScalarFiniteElement<D-1> >
class NGS_DLL_HEADER NeumannIntegrator 
  : public T_BIntegrator<DiffOpIdBoundary<D>, DVec<1>, FEL>
{
public:
  ///
  NeumannIntegrator (CoefficientFunction * coeff);

  static Integrator * Create (Array<CoefficientFunction*> & coeffs);
  ///  
  virtual bool BoundaryForm () const { return 1; }
  ///
  virtual string Name () const { return "Neumann"; }
};

template <int D, typename FEL = ScalarFiniteElement<D-1>  >
class ComplexNeumannIntegrator 
  : public T_BIntegrator<DiffOpIdBoundary<D>, DVec<1, Complex>, FEL>
{
public:
  ComplexNeumannIntegrator (CoefficientFunction * coeff);
  virtual bool BoundaryForm () const { return 1; }
  virtual string Name () const { return "ComplexNeumann"; }
};





/// integrator for \f$\int_\Gamma v_n \, ds\f$
template <int D, typename FEL = ScalarFiniteElement<D-1> >
class NormalNeumannIntegrator 
  : public T_BIntegrator<DiffOpNormal<D>, DVec<1>, FEL>
{
public:
  ///
  NormalNeumannIntegrator (CoefficientFunction * coeff)
    : T_BIntegrator<DiffOpNormal<D>, DVec<1>, FEL> (DVec<1> (coeff))
  { ; }

  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new NormalNeumannIntegrator (coeffs[0]);
  }

  ///  
  virtual bool BoundaryForm () const { return 1; }
  ///
  virtual string Name () const { return "NormalNeumann"; }
};






///
template <int D, typename FEL = ScalarFiniteElement<D>  >
class GradSourceIntegrator 
  : public T_BIntegrator<DiffOpGradient<D>, DVec<D>, FEL>
{
public:
  ///
  GradSourceIntegrator (CoefficientFunction * coeff1, CoefficientFunction *coeff2, CoefficientFunction *coeff3)
    : T_BIntegrator<DiffOpGradient<D>, DVec<D>, FEL> (DVec<D> (coeff1, coeff2, coeff3))
  { ; }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new GradSourceIntegrator (coeffs[0],coeffs[1],coeffs[2]);
  }

  ///
  virtual string Name () const { return "GradSource"; }
};










#ifdef NOTAVAILABLE
///
class DivIntegrator : public B1DB2Integrator<>
{
  ///
  CoefficientFunction * coeff;
  ///
  int dim;

public:
  ///
  DivIntegrator (int adim, CoefficientFunction * acoeff);
  ///
  ~DivIntegrator ();

  ///
  virtual int GetDimension1 () const  
  { return dim; }
  ///
  virtual int GetDimension2 () const  
  { return 1; }

  ///
  virtual int GetDimensionD1 () const 
  { return 1; }
  ///
  virtual int GetDimensionD2 () const 
  { return 1; }

  ///
  virtual int DiffOrder1 () const { return 1; }
  ///
  virtual int DiffOrder2 () const { return 0; }


  ///
  virtual void GenerateB1Matrix (const FiniteElement & fel,
				 const SpecificIntegrationPoint<> & ip,
				 ngbla::FlatMatrix<> & bmat,
				 LocalHeap & lh) const;
  ///
  virtual void GenerateB2Matrix (const FiniteElement & fel,
				 const SpecificIntegrationPoint<> & ip,
				 ngbla::FlatMatrix<> & bmat,
				 LocalHeap & lh) const;
  ///
  virtual void GenerateDMatrix (const FiniteElement & fel,
				const SpecificIntegrationPoint<> & ip,
				ngbla::FlatMatrix<> & dmat,
				LocalHeap & lh) const;


  ///
  virtual string Name () const { return "Div"; }
};

#endif


}

#endif
