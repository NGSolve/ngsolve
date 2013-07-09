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
template <int D, typename FEL = ScalarFiniteElement<D> >
class DiffOpGradient : public DiffOp<DiffOpGradient<D, FEL> >
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 1 };
  
  static string Name() { return "grad"; }

  static const FEL & Cast (const FiniteElement & fel) 
  { return static_cast<const FEL&> (fel); }

  ///
  template <typename MIP, typename MAT>
  static void GenerateMatrix (const FiniteElement & fel, 
			      const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    mat = Trans (mip.GetJacobianInverse ()) * 
      Trans (Cast(fel).GetDShape(mip.IP(),lh));
  }

  static void GenerateMatrix (const FiniteElement & fel, 
                              const MappedIntegrationPoint<D,D> & mip,
			      FlatMatrixFixHeight<D> & mat, LocalHeap & lh)
  {
    Cast(fel).CalcMappedDShape (mip, mat.Trans());
  }

  ///
  template <typename MIP, class TVX, class TVY>
  static void Apply (const FiniteElement & fel, const MIP & mip,
		     const TVX & x, TVY & y,
		     LocalHeap & lh) 
  {
    typedef typename TVX::TSCAL TSCAL;

    Vec<D,TSCAL> hv = Trans (Cast(fel).GetDShape(mip.IP(), lh)) * x;
    y = Trans (mip.GetJacobianInverse()) * hv;
  }

  using DiffOp<DiffOpGradient<D, FEL> >::ApplyIR;
  
  template <class MIR>
  static void ApplyIR (const FiniteElement & fel, const MIR & mir,
		       const FlatVector<double> x, FlatMatrix<double> y,
		       LocalHeap & lh)
  {
    FlatMatrixFixWidth<D> grad(mir.Size(), &y(0));
    Cast(fel).EvaluateGrad (mir.IR(), x, grad);
    for  (int i = 0; i < mir.Size(); i++)
      {
	Vec<D> hv = grad.Row(i);
	grad.Row(i) = Trans (mir[i].GetJacobianInverse()) * hv;
      }
  }


  ///
  template <typename MIP, class TVX, class TVY>
  static void ApplyTrans (const FiniteElement & fel, const MIP & mip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh) 
  {
    typedef typename TVX::TSCAL TSCAL;

    Vec<D,TSCAL> hv = mip.GetJacobianInverse() * x;
    y = Cast(fel).GetDShape(mip.IP(),lh) * hv;
  }
};




/// Boundary Gradient operator of dimension D
template <int D, typename FEL = ScalarFiniteElement<D-1> >
class DiffOpGradientBoundary : public DiffOp<DiffOpGradientBoundary<D, FEL> >
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 1 };

  ///
  template <typename AFEL, typename MIP, typename MAT>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    mat = Trans (mip.GetJacobianInverse ()) * 
      Trans (static_cast<const FEL&>(fel).GetDShape(mip.IP(),lh));
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
  template <typename FEL, typename MIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    typedef typename MAT::TSCAL TSCAL;

    mat =  Trans (mip.GetJacobianInverse ()) * 
      Trans (fel.GetDShape(mip.IP(),lh));

    double cx = mip.GetPoint()(0);
    if (cx == 0) cx = 1e-10;
    for (int i = 0; i < mat.Width(); i++)
      mat(0,i) += fel.GetShape(mip.IP(), lh)(i) / cx;

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
template <int D, typename FEL = ScalarFiniteElement<D> >
class DiffOpId : public DiffOp<DiffOpId<D, FEL> >
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D };
  enum { DIM_DMAT = 1 };
  enum { DIFFORDER = 0 };

  static const FEL & Cast (const FiniteElement & fel) 
  { return static_cast<const FEL&> (fel); }

  template <typename MIP, typename MAT>
  static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    mat.Row(0) = Cast(fel).GetShape(mip.IP(), lh);
  }

  static void GenerateMatrix (const FiniteElement & fel, 
			      const MappedIntegrationPoint<D,D> & mip,
			      FlatMatrixFixHeight<1> & mat, LocalHeap & lh)
  {
    Cast(fel).CalcShape (mip.IP(), FlatVector<> (fel.GetNDof(), &mat(0,0)));
  }


  template <typename MIP, class TVX, class TVY>
  static void Apply (const FiniteElement & fel, const MIP & mip,
		     const TVX & x, TVY & y,
		     LocalHeap & lh) 
  {
    y = Trans (Cast(fel).GetShape (mip.IP(), lh)) * x;
  }

  static void Apply (const FiniteElement & fel, const MappedIntegrationPoint<D,D> & mip,
		     const FlatVector<double> & x, FlatVector<double> & y,
		     LocalHeap & lh) 
  {
    y(0) = Cast(fel).Evaluate(mip.IP(), x);
  }


  using DiffOp<DiffOpId<D, FEL> >::ApplyIR;

  template <class MIR>
  static void ApplyIR (const FiniteElement & fel, const MIR & mir,
		       const FlatVector<double> x, FlatMatrix<double> y,
		       LocalHeap & lh)
  {
    Cast(fel).Evaluate (mir.IR(), x, FlatVector<> (mir.Size(), &y(0,0)));
  }



  template <typename MIP, class TVX, class TVY>
  static void ApplyTrans (const FiniteElement & fel, const MIP & mip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh) 
  {
    y = Cast(fel).GetShape (mip.IP(), lh) * x;
  }


  // using DiffOp<DiffOpId<D, FEL> >::ApplyTransIR;
  template <class MIR>
  static void ApplyTransIR (const FiniteElement & fel, 
                            const MIR & mir,
                            FlatMatrix<double> x, FlatVector<double> y,
                            LocalHeap & lh)
  {
    Cast(fel).EvaluateTrans (mir.IR(), FlatVector<> (mir.Size(), &x(0,0)), y);
  }

  template <class MIR>
  static void ApplyTransIR (const FiniteElement & fel, 
                            const MIR & mir,
                            FlatMatrix<Complex> x, FlatVector<Complex> y,
                            LocalHeap & lh)
  {
    DiffOp<DiffOpId<D, FEL> > :: ApplyTransIR (fel, mir, x, y, lh);    
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

  template <typename FEL, typename MIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    const FlatVector<> shape = fel.GetShape (mip.IP(), lh);
  

    typedef typename MAT::TSCAL TSCAL;
    mat = TSCAL(0.); 
    for (int j = 0; j < shape.Height(); j++)
      for (int i = 0; i < SYSDIM; i++)
	mat(i, j*SYSDIM+i) = shape(j);
  }
};





/// Identity on boundary
template <int D, typename FEL = ScalarFiniteElement<D-1> >
class DiffOpIdBoundary : public DiffOp<DiffOpIdBoundary<D, FEL> >
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = 1 };
  enum { DIFFORDER = 0 };

  template <typename AFEL, typename MIP, typename MAT>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    mat.Row(0) = static_cast<const FEL&>(fel).GetShape(mip.IP(), lh);
  }

  template <typename AFEL, typename MIP, class TVX, class TVY>
  static void Apply (const AFEL & fel, const MIP & mip,
		     const TVX & x, TVY & y,
		     LocalHeap & lh) 
  {
    y = Trans (static_cast<const FEL&>(fel).GetShape (mip.IP(), lh)) * x;
  }

  static void Apply (const ScalarFiniteElement<D-1> & fel, const MappedIntegrationPoint<D-1,D> & mip,
		     const FlatVector<double> & x, FlatVector<double> & y,
		     LocalHeap & lh) 
  {
    y(0) = static_cast<const FEL&>(fel).Evaluate(mip.IP(), x);
  }


  template <typename AFEL, typename MIP, class TVX, class TVY>
  static void ApplyTrans (const AFEL & fel, const MIP & mip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh) 
  {
    y = static_cast<const FEL&>(fel).GetShape (mip.IP(), lh) * x;
  }




  // using DiffOp<DiffOpIdBoundary<D, FEL> >::ApplyTransIR;
  template <class MIR>
  static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
                            FlatMatrix<double> x, FlatVector<double> y,
                            LocalHeap & lh)
  {
    // static Timer t("applytransir - bnd");
    // RegionTimer reg(t);

    static_cast<const FEL&>(fel).
      EvaluateTrans (mir.IR(), FlatVector<> (mir.Size(), &x(0,0)), y);
  }

  template <class MIR>
  static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
                            FlatMatrix<Complex> x, FlatVector<Complex> y,
                            LocalHeap & lh)
  { 
    DiffOp<DiffOpIdBoundary<D, FEL> > :: ApplyTransIR (fel, mir, x, y, lh);    
    // static_cast<const FEL&>(fel).
    // EvaluateTrans (mir.IR(), FlatVector<> (mir.Size(), &x(0,0)), y);
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

  template <typename FEL, typename MIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    mat = 0.; 
    const FlatVector<> shape = fel.GetShape (mip.IP(), lh);
    for (int j = 0; j < shape.Height(); j++)
      for (int i = 0; i < SYSDIM; i++)
	mat(i, j*SYSDIM+i) = shape(j);
  }
};
























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

  template <typename FEL, typename MIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const MIP & mip,
		       MAT & mat, LocalHeap & lh) const
  {
    typedef typename MAT::TSCAL TRESULT;
    TRESULT val = coef -> T_Evaluate<TRESULT> (mip);
    /*
    mat = TRESULT(0);
    for (int i = 0; i < DIM; i++)
      mat(i, i) = val;
    */
    mat = val * Id<DIM>();
  }  

  template <typename FEL, typename MIR, typename MAT>
  void GenerateMatrixIR (const FEL & fel, const MIR & mir,
			 const FlatArray<MAT> & mats, LocalHeap & lh) const
  {
    typedef typename MAT::TSCAL TRESULT;
    FlatMatrix<TRESULT> vals(mir.IR().GetNIP(), 1, lh);
    coef -> Evaluate (mir, vals);
    
    for (int j = 0; j < mir.IR().GetNIP(); j++)
      {
	mats[j] = vals(j,0) * Id<DIM>();
	/*
	mats[j] = TRESULT(0.0);
	for (int i = 0; i < DIM; i++)
	  mats[j](i, i) = vals(j,0);
	*/
      }
  }  

  template <typename FEL, class VECX, class VECY>
  void Apply (const FEL & fel, const BaseMappedIntegrationPoint & mip,
	      const VECX & x, VECY & y, LocalHeap & lh) const
  {
    typedef typename VECY::TSCAL TRESULT;
    TRESULT val = coef -> T_Evaluate<TRESULT> (mip);
    for (int i = 0; i < DIM; i++)
      y(i) = val * x(i);
  }

  template <typename FEL, class VECX>
  void Apply1 (const FEL & fel, const BaseMappedIntegrationPoint & mip,
	       const VECX & x, LocalHeap & lh) const
  {
    typedef typename VECX::TSCAL TSCAL;
    x *= coef -> T_Evaluate<TSCAL> (mip);
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

  template <typename FEL, typename MIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const MIP & mip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat(0,0) = Evaluate (*coef, mip);
  }  

  template <typename FEL, typename MIP, class VECX, class VECY>
  void Apply (const FEL & fel, const MIP & mip,
	      const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef, mip) * x(0);
  }

  template <typename FEL, typename MIP, class VECX, class VECY>
  void ApplyTrans (const FEL & fel, const MIP & mip,
		   const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef, mip) * x(0);
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
  
  template <typename FEL, typename MIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const MIP & mip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    mat(0,0) = Evaluate (*coef1, mip);
    mat(1,1) = Evaluate (*coef2, mip);
  }  

  template <typename FEL, typename MIP>
  void GetEigensystem (const FEL & fel, const MIP & mip, 
		  Array<double> & eigenvalues,
		  Array<double> & eigenvectors,
		  LocalHeap & lh) const
  {
    eigenvalues[0] = Evaluate (*coef1, mip);
    eigenvalues[1] = Evaluate (*coef2, mip);
    eigenvectors[0] = eigenvectors[3] = 1.;
    eigenvectors[1] = eigenvectors[2] = 0.;
  }

  template <typename FEL, typename MIP, class VECX, class VECY>
  void Apply (const FEL & fel, const MIP & mip,
	      const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef1, mip) * x(0);
    y(1) = Evaluate (*coef2, mip) * x(1);
  }

  template <typename FEL, typename MIP, class VECX, class VECY>
  void ApplyTrans (const FEL & fel, const MIP & mip,
		   const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef1, mip) * x(0);
    y(1) = Evaluate (*coef2, mip) * x(1);
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
  
  template <typename FEL, typename MIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const MIP & mip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    mat(0,0) = Evaluate (*coef1, mip);
    mat(1,1) = Evaluate (*coef2, mip);
    mat(2,2) = Evaluate (*coef3, mip);
  }  

  
  template <typename FEL, typename MIP>
  void GetEigensystem (const FEL & fel, const MIP & mip, 
		  Array<double> & eigenvalues,
		  Array<double> & eigenvectors,
		  LocalHeap & lh) const
  {
    
    eigenvalues[0] = Evaluate(*coef1,mip);
    eigenvalues[1] = Evaluate(*coef2,mip);
    eigenvalues[2] = Evaluate(*coef3,mip);

    eigenvectors = 0.;
    eigenvectors[0] = eigenvectors[4] = eigenvectors[8] = 1.;
  }


  template <typename FEL, typename MIP, class VECX, class VECY>
  void Apply (const FEL & fel, const MIP & mip,
	      const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef1, mip) * x(0);
    y(1) = Evaluate (*coef2, mip) * x(1);
    y(2) = Evaluate (*coef3, mip) * x(2);
  }

  template <typename FEL, typename MIP, class VECX, class VECY>
  void ApplyTrans (const FEL & fel, const MIP & mip,
		   const VECX & x, VECY & y, LocalHeap & lh) const
  {
    y(0) = Evaluate (*coef1, mip) * x(0);
    y(1) = Evaluate (*coef2, mip) * x(1);
    y(2) = Evaluate (*coef3, mip) * x(2);
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

  template <typename FEL, typename MIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const MIP & mip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat(0,0) = Evaluate (*coef, mip);
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
  
  template <typename FEL, typename MIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const MIP & mip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    mat(0,0) = Evaluate (*coef00, mip);
    mat(0,1) = mat(1,0) = Evaluate (*coef01, mip);
    mat(1,1) = Evaluate (*coef11, mip);
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
  
  template <typename FEL, typename MIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const MIP & mip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    mat(0,0) = Evaluate (*coef00, mip);
    mat(1,0) = mat(0,1) = Evaluate (*coef10, mip);
    mat(1,1) = Evaluate (*coef11, mip);
    mat(2,0) = mat(0,2) = Evaluate (*coef20, mip);
    mat(2,1) = mat(1,2) = Evaluate (*coef21, mip);
    mat(2,2) = Evaluate (*coef22, mip);
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

  template <typename FEL, typename MIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const MIP & mip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    double val = Evaluate (*coef, mip);
    Vec<DIM> nv = mip.GetNV();
    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
	mat(i, j) = val * nv(i) * nv(j);
  }  
};






template <int N>
class NGS_DLL_HEADER DVec  // : public DVecBase<N,T>
{
  CoefficientFunction * coefs[N];
  bool vectorial;
public:
  DVec (Array<CoefficientFunction*> & acoeffs)
  {
    vectorial = (N > 1) && (N == acoeffs[0]->Dimension());

    if (vectorial)
      coefs[0] = acoeffs[0];
    else
      for (int i = 0; i < N; i++)
	coefs[i] = acoeffs[i];
  }
  
  DVec (CoefficientFunction * acoef)
  { 
    vectorial = (N > 1) && (N == acoef->Dimension());
    coefs[0] = acoef;
  }

  DVec (CoefficientFunction * acoef1,
	CoefficientFunction * acoef2,
	CoefficientFunction * acoef3 = NULL,
	CoefficientFunction * acoef4 = NULL,
	CoefficientFunction * acoef5 = NULL,
	CoefficientFunction * acoef6 = NULL)
  { 
    vectorial = false;

    coefs[0] = acoef1;
    coefs[1] = acoef2;
    if (N >= 3) coefs[2] = acoef3;
    if (N >= 4) coefs[3] = acoef4;
    if (N >= 5) coefs[4] = acoef5;
    if (N >= 6) coefs[5] = acoef6;
  }
    

  template <typename FEL, typename MIP, typename VEC>
  void GenerateVector (const FEL & fel, const MIP & mip,
		       VEC & vec, LocalHeap & lh) const
  {
    typedef typename VEC::TSCAL TSCAL;

    if (vectorial)
      {
	coefs[0] -> Evaluate (mip, vec);
      }
    else
      for (int i = 0; i < N; i++)
	{
	  CoefficientFunction * hp = coefs[i];
	  vec(i) = hp -> T_Evaluate<TSCAL> (mip);

	  // gcc 4.6 complains, why ???? (JS)
	  // vec(i) = coefs[i] -> T_Evaluate<TSCAL> (mip);  
	}
  }  



  template <typename FEL, typename MIR, typename VEC>
  void GenerateVectorIR (const FEL & fel, const MIR & mir,
			 VEC & vecs, LocalHeap & lh) const
  {
    typedef typename VEC::TSCAL TSCAL;

    if (vectorial || N == 1)
      {
	coefs[0] -> Evaluate (mir, vecs);
      }
    else
      for (int j = 0; j < mir.Size(); j++)
	for (int i = 0; i < N; i++)
	  {
	    CoefficientFunction * hp = coefs[i];
	    vecs(j,i) = hp -> T_Evaluate<TSCAL> (mir[j]);
	    
	    // gcc 4.6 complains, why ???? (JS)
	    // vec(i) = (coefs[i]) -> T_Evaluate<TSCAL> (mip);  
	  }
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
  
  template <typename FEL, typename MIP, typename VEC>
  void GenerateVector (const FEL & fel, const MIP & mip,
		       VEC & vec, LocalHeap & lh) const
  {
    Vec<N> hv;
    coef -> Evaluate (mip, hv);
    for (int i = 0; i < N; i++)
      vec(i) = hv(i);
  }  

  template <typename FEL, typename MIR, typename VEC>
  void GenerateVectorIR (const FEL & fel, const MIR & mir,
			 VEC & vecs, LocalHeap & lh) const
  {
    for (int i = 0; i < mir.Size(); i++)
      GenerateVector (fel, mir[i], vecs.Row(i), lh);
  }
};



template <int N>
class TVec
{
  CoefficientFunction * coef;

public:
  typedef double TSCAL;
  TVec (CoefficientFunction * acoef) : coef(acoef) {;}

  

  template <typename FEL, typename MIP, typename VEC>
  void GenerateVector (const FEL & fel, const MIP & mip,
		       VEC & vec, LocalHeap & lh) const
  {
    vec = 0.0;

    typedef typename VEC::TSCAL TSCAL;
    
    TSCAL length = 0.;
    for(int i=0; i<N; i++)
      {
	//vec(i) = mip.GetJacobian()(i,0);
	vec(i) = mip.GetTV()(i);
	length += vec(i)*vec(i);
      }
    //(*testout) << "point " << mip.GetPoint() << " tv " << vec;
    vec *= Evaluate (*coef, mip)/sqrt(length);
    //(*testout) << " retval " << vec << endl;
  }

  template <typename FEL, typename MIR, typename VEC>
  void GenerateVectorIR (const FEL & fel, const MIR & mir,
			 VEC & vecs, LocalHeap & lh) const
  {
    for (int i = 0; i < mir.Size(); i++)
      GenerateVector (fel, mir[i], vecs.Row(i), lh);
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

  template <typename FEL, typename MIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const MIP & mip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;
    const double r = mip.GetPoint()(0);
    double val = r*Evaluate (*coef, mip);
    for (int i = 0; i < DIM; i++)
      mat(i, i) = val;
  }  
  

  template <typename FEL, typename MIP, class VECX, class VECY>
  void Apply (const FEL & fel, const MIP & mip,
	      const VECX & x, VECY & y, LocalHeap & lh) const
  {
    const double r = mip.GetPoint()(0);
    double val = r*Evaluate (*coef, mip);
    y = val * x;
  }
};






/// Identity on boundary
template <int D, typename FEL = ScalarFiniteElement<D> >
class DiffOpNormal : public DiffOp<DiffOpNormal<D, FEL> >
{
public:
  enum { DIM = D };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = 1 };
  enum { DIFFORDER = 0 };

  template <typename MIP, typename MAT>
  static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    FlatVector<> shape = static_cast<const FEL&> (fel).GetShape (mip.IP(), lh);
    Vec<D> nv = mip.GetNV();
    //Vec<D> p = mip.GetPoint();
    for (int j = 0; j < shape.Size(); j++)
      for (int i = 0; i < D; i++)
	mat(0, j*D+i) =  shape(j) * nv(i);

    //(*testout) << "mip = " << p << ", nv = " << nv << endl;
    //p(0) = 0.0;
    //p /= L2Norm(p);
    //nv /= L2Norm(nv);
    //(*testout) << "normalized, mip = " << p << ", nv = " << nv << endl;
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
  virtual ~LaplaceIntegrator ();
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
template <int D, typename FEL = ScalarFiniteElement<D> >
class NGS_DLL_HEADER MassIntegrator 
  : public T_BDBIntegrator<DiffOpId<D>, DiagDMat<1>, FEL >
{
public:
  ///
  MassIntegrator (CoefficientFunction * coeff);
  ///
  MassIntegrator (Array<CoefficientFunction*> & coeffs);
  ///
  virtual ~MassIntegrator ();
  ///
  virtual string Name () const { return "Mass"; }
};






/// integrator for \f$\int_\Gamma u v \, ds\f$
  template <int D, typename FEL = ScalarFiniteElement<D-1> >
class NGS_DLL_HEADER RobinIntegrator 
  : public T_BDBIntegrator<DiffOpIdBoundary<D>, DiagDMat<1>, FEL>
{
  typedef T_BDBIntegrator<DiffOpIdBoundary<D>, DiagDMat<1>, FEL> BASE;
public:
  ///
  RobinIntegrator (CoefficientFunction * coeff);
  ///
  RobinIntegrator (Array<CoefficientFunction*> & coeffs);
  ///
  virtual ~RobinIntegrator ();
  ///
  virtual bool BoundaryForm () const { return 1; }
  ///
  virtual string Name () const { return "Robin"; }
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


/// 
template <int D> 
class DiffOpDiv : public DiffOp<DiffOpDiv<D> >
{
public:
  enum { DIM = D };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D };
  enum { DIM_DMAT = 1 };
  enum { DIFFORDER = 1 };

  template <typename FEL, typename MIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    int nd = fel.GetNDof();

    FlatMatrix<> grad (D, nd, lh);
    grad = Trans (mip.GetJacobianInverse ()) * 
      Trans (fel.GetDShape(mip.IP(), lh));
    
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

  template <typename FEL, typename MIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    int nd = fel.GetNDof();

    FlatMatrix<> grad (2, nd, lh);
    grad = Trans (mip.GetJacobianInverse ()) * 
      Trans (fel.GetDShape(mip.IP(), lh));
    
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

  template <typename FEL, typename MIP, typename MAT>
  static void GenerateMatrix (const FEL & fel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    int nd = fel.GetNDof();

    FlatMatrix<> grad (3, nd, lh);
    grad = Trans (mip.GetJacobianInverse ()) * 
      Trans (fel.GetDShape(mip.IP(), lh));
    
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
  virtual ~SourceIntegrator ();
  virtual string Name () const { return "Source"; }
};


///
template <int D, typename FEL = ScalarFiniteElement<D-1> >
class NGS_DLL_HEADER NeumannIntegrator 
  : public T_BIntegrator<DiffOpIdBoundary<D>, DVec<1>, FEL>
{
public:
  ///
  NeumannIntegrator (CoefficientFunction * coeff);
  NeumannIntegrator (Array<CoefficientFunction*> & coeffs);
  virtual ~NeumannIntegrator ();
  ///  
  virtual bool BoundaryForm () const { return true; }
  ///
  virtual string Name () const { return "Neumann"; }
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
  GradSourceIntegrator (CoefficientFunction * coeff1, 
			CoefficientFunction * coeff2, 
			CoefficientFunction * coeff3)
    : T_BIntegrator<DiffOpGradient<D>, DVec<D>, FEL> (DVec<D> (coeff1, coeff2, coeff3))
  { ; }
  
  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new GradSourceIntegrator (coeffs[0],coeffs[1],coeffs[2]);
  }

  ///
  virtual string Name () const { return "GradSource"; }
};









}

#endif
