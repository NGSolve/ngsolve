#ifndef FILE_BDBINTEGRATOR
#define FILE_BDBINTEGRATOR

/*********************************************************************/
/* File:   bdbintegrator.hpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


namespace ngfem
{






/**
   Coefficient tensor.
   Base-class for template-polymorphismus.
   Provides application and transpose-application
*/
  template <class DMO, int DIM_DMAT>
class DMatOp
{
public:
  // typedef double TSCAL;
  /// is coefficient tensor symmetric ?
  enum { SYMMETRIC = 1 };
  
  template <typename SCAL>
  static Mat<DIM_DMAT,DIM_DMAT,SCAL> GetMatrixType(SCAL s) { return SCAL(0); }

  template <typename FEL, typename MIR, typename MAT>
  void GenerateMatrixIR (const FEL & fel, const MIR & mir, 
			 const FlatArray<MAT> & dmats, 
			 LocalHeap & lh) const
  {
    for (int i = 0; i < mir.IR().GetNIP(); i++)
      static_cast<const DMO*>(this) -> GenerateMatrix (fel, mir[i], dmats[i], lh);
  }


  /// generate linearized matrix in linearization point vec
  template <typename FEL, typename MIP, typename VEC, typename MAT>
  void GenerateLinearizedMatrix (const FEL & fel, const MIP & mip, VEC & vec,
				 MAT & mat, LocalHeap & lh) const
  {
    static_cast<const DMO*>(this) -> GenerateMatrix (fel, mip, mat, lh);
  }

  /// apply coefficient matrix.
  template <typename FEL, typename MIP, class TVX, class TVY>
  void Apply (const FEL & fel, const MIP & mip,
	      const TVX & x, TVY & y,
	      LocalHeap & lh) const
  {
    Mat<DMO::DIM_DMAT, DMO::DIM_DMAT, double> mat;
    static_cast<const DMO*>(this) -> GenerateMatrix (fel, mip, mat, lh);
    y = mat * x;
  }

  template <typename FEL, typename MIP, class TVX>
  void Apply1 (const FEL & fel, const MIP & mip,
	       TVX & x, LocalHeap & lh) const
  {
    Vec<DMO::DIM_DMAT, typename TVX::TSCAL> y;
    static_cast<const DMO*>(this) -> Apply (fel, mip, x, y, lh);
    x = y;
  }

  template <typename FEL, typename MIR, typename TVX>
  void ApplyIR (const FEL & fel, const MIR & mir,
		TVX & x, LocalHeap & lh) const
  {
    for (int i = 0; i < mir.Size(); i++)
      static_cast<const DMO*>(this) -> Apply1 (fel, mir[i], x.Row(i), lh);
  }


  template <typename FEL, typename MIP, class TVX, class TVY>
  void ApplyInv (const FEL & fel, const MIP & mip,
		 const TVX & x, TVY & y,
		 LocalHeap & lh) const
  {
    Mat<DMO::DIM_DMAT, DMO::DIM_DMAT, double> mat;
    Mat<DMO::DIM_DMAT, DMO::DIM_DMAT, double> inv;

    static_cast<const DMO*>(this) -> GenerateMatrix (fel, mip, mat, lh);
    CalcInverse (mat, inv);
    y = inv * x;
  }

  /// apply transpose coefficient tensor
  template <typename FEL, typename MIP, class TVX, class TVY>
  void ApplyTrans (const FEL & fel, const MIP & mip,
		   const TVX & x, TVY & y,
		   LocalHeap & lh) const
  {
    Mat<DMO::DIM_DMAT, DMO::DIM_DMAT, double> mat;
    static_cast<const DMO*>(this) -> GenerateMatrix (fel, mip, mat, lh);
    y = Trans (mat) * x;
  }

  /// computes energy 
  template <typename FEL, typename MIP, class TVX>
  double Energy (const FEL & fel, const MIP & mip,
		 const TVX & x, LocalHeap & lh) const  
  {
    TVX y;
    static_cast<const DMO*>(this) -> Apply (fel, mip, x, y, lh);
    return 0.5 * InnerProduct (x,y);
  }
};






#ifdef WIN32
#define __restrict__ __restrict
#endif

template <int M> NGS_DLL_HEADER
void FastMat (int n, Complex * ba, Complex *  pb, Complex * pc);

template <int M> NGS_DLL_HEADER
void FastMat (int n, double * __restrict__ ba, double *  __restrict__ pb, double * __restrict__ pc);
  





 
/**
   Element assembling.
   Assembling for bilinear-forms of type $\int (B v) : D (B u) dx$. \\
   Template argument DiffOp provides differential operator, i.e. B matrix,
   (e.g. gradient, strain operator, curl,...) \\
   DmatOp provides d-matrix (e.g. diagonal, anisotropic, plane stress, ...) \\
   FEL is element type to assemble matrix for (ScalarFiniteElement, 
   HCurlFiniteElement, FE_Trig1, ...)
 */
template <class DIFFOP, class DMATOP, class FEL = FiniteElement>
class T_BDBIntegrator : public virtual BilinearFormIntegrator
{
protected:
  DMATOP dmatop;
public:
  // typedef double TSCAL;  // necessary ?

  enum { DIM_SPACE   = DIFFOP::DIM_SPACE };
  enum { DIM_ELEMENT = DIFFOP::DIM_ELEMENT };
  enum { DIM_DMAT    = DIFFOP::DIM_DMAT };
  enum { DIM         = DIFFOP::DIM };
  // typedef typename DMATOP::TSCAL TSCAL;

  ///
  T_BDBIntegrator  (Array<CoefficientFunction*> & coeffs)
  : dmatop(coeffs)
  { ; }

  ///
  T_BDBIntegrator (const DMATOP & admat)
    : dmatop(admat)
  { ; }


  ///
  virtual ~T_BDBIntegrator ()
  { ; }

  ///
  virtual bool BoundaryForm () const
  { return int(DIM_SPACE) > int(DIM_ELEMENT); }

  virtual int DimElement () const
  { return DIM_ELEMENT; }

  virtual int DimSpace () const
  { return DIM_SPACE; }

  virtual int DimFlux () const
  { return DIM_DMAT; }

  DMATOP & DMat () { return dmatop; }
  const DMATOP & DMat () const { return dmatop; }

  virtual void CheckElement (const FiniteElement & el) const
  {
    if (!dynamic_cast<const FEL*> (&el) )
      throw Exception (string ("Element does not match integrator\n") +
                       string ("element type is ") + typeid(el).name() +
                       string (" expected type is ") + typeid(FEL).name() +
                       string (" integrator is ") + Name());
  }


  virtual void ApplyDMat (const FiniteElement & bfel,
			  const BaseMappedIntegrationPoint & bmip,
			  const FlatVector<double> & elx, 
			  FlatVector<double> & eldx,
			  LocalHeap & lh) const
  {
    dmatop.Apply(static_cast<const FEL&> (bfel),
		 static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> &>(bmip),
		 elx, eldx ,lh);
  }

  virtual void ApplyDMat (const FiniteElement & bfel,
			  const BaseMappedIntegrationPoint & bmip,
			  const FlatVector<Complex> & elx, 
			  FlatVector<Complex> & eldx,
			  LocalHeap & lh) const
  {
    dmatop.Apply(static_cast<const FEL&> (bfel),
		 static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> &>(bmip),
		 elx,eldx,lh);
  }

  virtual void ApplyDMat (const FiniteElement & bfel,
			  const BaseMappedIntegrationRule & bmir,
			  const FlatMatrix<double> & elx, 
			  FlatMatrix<double> & eldx,
			  LocalHeap & lh) const
  {
    const FEL & fel = static_cast<const FEL&> (bfel);
    const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);

    for (int i = 0; i < mir.Size(); i++)
      dmatop.Apply (fel, mir[i], elx.Row(i), eldx.Row(i), lh);
  }


  virtual void ApplyDMatInv (const FiniteElement & bfel,
			     const BaseMappedIntegrationRule & bmir,
			     const FlatMatrix<double> & elx, 
			     FlatMatrix<double> & eldx,
			     LocalHeap & lh) const
  {
    const FEL & fel = static_cast<const FEL&> (bfel);
    const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);

    for (int i = 0; i < mir.Size(); i++)
      dmatop.ApplyInv (fel, mir[i], elx.Row(i), eldx.Row(i), lh);
  }






  virtual void ApplyDMat (const FiniteElement & bfel,
			  const BaseMappedIntegrationRule & bmir,
			  const FlatMatrix<Complex> & elx, 
			  FlatMatrix<Complex> & eldx,
			  LocalHeap & lh) const
  {
    const FEL & fel = static_cast<const FEL&> (bfel);
    const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);

    for (int i = 0; i < mir.Size(); i++)
      dmatop.Apply (fel, mir[i], elx.Row(i), eldx.Row(i), lh);
  }



  virtual void
  CalcElementMatrix (const FiniteElement & bfel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    T_CalcElementMatrix<double> (bfel, eltrans, elmat, lh);
  }

  virtual void
  CalcElementMatrix (const FiniteElement & bfel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<Complex> & elmat,
		     LocalHeap & lh) const
  {
    T_CalcElementMatrix<Complex> (bfel, eltrans, elmat, lh);
  }


#ifdef TEXT_BOOK_VERSION

  template <typename TSCAL>
  void T_CalcElementMatrix (const FiniteElement & bfel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<TSCAL> & elmat,
			    LocalHeap & lh) const
  {

    try
      {
	const FEL & fel = static_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();
        
	elmat = 0;

	FlatMatrixFixHeight<DIM_DMAT, TSCAL> bmat (ndof * DIM, lh);
	FlatMatrixFixHeight<DIM_DMAT, TSCAL> dbmat (ndof * DIM, lh);
	Mat<DIM_DMAT,DIM_DMAT,TSCAL> dmat;

	const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());

	for (int i = 0 ; i < ir.GetNIP(); i++)
	  {
	    HeapReset hr(lh);

	    MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> 
	      mip(ir[i], eltrans, lh);

	    dmatop.GenerateMatrix (fel, mip, dmat, lh);

	    DIFFOP::GenerateMatrix (fel, mip, bmat, lh);
	    double fac = mip.GetMeasure() * mip.IP().Weight();
	    
	    dbmat = fac * (dmat * bmat);
	    elmat += Trans (bmat) * dbmat; 
	  } 
      }

    catch (Exception & e)
      {
	e.Append ("in CalcElementMatrix - textbook, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin CalcElementMatrix - textbook, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }
  }
#endif




#ifdef __SSE3__
#define BLOCK_VERSION
#endif

#ifdef BLOCK_VERSION

  template <typename TSCAL>
  void T_CalcElementMatrix (const FiniteElement & fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<TSCAL> & elmat,
			    LocalHeap & lh) const
  {
    // static int timer = NgProfiler::CreateTimer (string ("Elementmatrix, ") + Name());
    // NgProfiler::RegionTimer reg (timer);

    try
      {
	int ndof = fel.GetNDof();
	elmat = 0;

        enum { BLOCK = 2 * (12 / DIM_DMAT + 1) };
	
        HeapReset hr1(lh);
	
        FlatMatrixFixHeight<DIM_DMAT, double> bmatr (ndof * DIM, lh);
        FlatMatrixFixHeight<DIM_DMAT, TSCAL> bmat1 (ndof * DIM, lh);
	FlatMatrixFixHeight<DIM_DMAT, TSCAL> dbmat1 (ndof * DIM, lh);

        FlatMatrixFixHeight<2*DIM_DMAT, TSCAL> bmat2 (ndof * DIM, lh);
	FlatMatrixFixHeight<2*DIM_DMAT, TSCAL> dbmat2 (ndof * DIM, lh);
	
        FlatMatrixFixHeight<DIM_DMAT*BLOCK, TSCAL> bbmat (ndof * DIM, lh);
        FlatMatrixFixHeight<DIM_DMAT*BLOCK, TSCAL> bdbmat (ndof * DIM, lh);

        typedef decltype (DMATOP::GetMatrixType(TSCAL(0))) TDMAT;
	
        const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());
	MappedIntegrationRule<DIM_ELEMENT, DIM_SPACE> mir(ir, eltrans, lh);

        FlatArray<TDMAT> dmats(ir.GetNIP(), lh);
	dmatop.GenerateMatrixIR (fel, mir, dmats, lh);

        int i = 0;
        for (int i1 = 0; i1 < ir.GetNIP() / BLOCK; i1++)
          {
            for (int i2 = 0; i2 < BLOCK; i++, i2++)
              {
                HeapReset hr(lh);

                DIFFOP::GenerateMatrix (fel, mir[i], bmatr, lh);
		TDMAT dmat = mir[i].GetWeight() * dmats[i];
                
		bbmat.Rows(i2*DIM_DMAT, (i2+1)*DIM_DMAT) = bmatr;
                bdbmat.Rows(i2*DIM_DMAT, (i2+1)*DIM_DMAT) = dmat * bmatr; 
              }
            
            if (DMATOP::SYMMETRIC)
              FastMat<DIM_DMAT*BLOCK> (elmat.Height(), &bdbmat(0,0), &bbmat(0,0), &elmat(0,0));
            else
              elmat += Trans (bbmat) * bdbmat; 
          } 

        for ( ; i < ir.GetNIP()-1; i+=2)
          {
            HeapReset hr(lh);

            DIFFOP::GenerateMatrix (fel, mir[i], bmatr, lh);
            TDMAT dmat = mir[i].GetWeight() * dmats[i];

            bmat2.Rows(0,DIM_DMAT) = bmatr;
            dbmat2.Rows(0,DIM_DMAT) = dmat * bmatr;

            DIFFOP::GenerateMatrix (fel, mir[i+1], bmatr, lh);
            dmat = mir[i+1].GetWeight() * dmats[i+1];

            bmat2.Rows(DIM_DMAT,2*DIM_DMAT) = bmatr;
            dbmat2.Rows(DIM_DMAT,2*DIM_DMAT) = dmat * bmatr;
            
            if (DMATOP::SYMMETRIC)
              FastMat<2*DIM_DMAT> (elmat.Height(), & dbmat2(0,0), &bmat2(0,0), &elmat(0,0));
            else
              elmat += Trans (bmat2) * dbmat2;
          } 

        for ( ; i < ir.GetNIP(); i++)
          {
            HeapReset hr(lh);

            DIFFOP::GenerateMatrix (fel, mir[i], bmatr, lh);
            TDMAT dmat = mir[i].GetWeight() * dmats[i];

            bmat1 = bmatr;
            dbmat1 = dmat * bmatr;
            
            if (DMATOP::SYMMETRIC)
              // elmat += Symmetric (Trans (dbmat) * bmat);
              FastMat<DIM_DMAT> (elmat.Height(), &dbmat1(0,0), &bmat1(0,0), &elmat(0,0));
            else
              elmat += Trans (bmat1) * dbmat1;
          } 
        
        if (DMATOP::SYMMETRIC)
          {
            for (int i = 0; i < elmat.Height(); i++)
              for (int j = 0; j < i; j++)
                elmat(j,i) = elmat(i,j);
          }
      }

    catch (Exception & e)
      {
	e.Append ("in CalcElementMatrix - blockversion, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin CalcElementMatrix - blockversion, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }
  }



#else // blockversion
  // one matrix matrix multiplication
  ///
  template <typename TSCAL>
  void T_CalcElementMatrix (const FiniteElement & fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<TSCAL> & elmat,
			    LocalHeap & lh) const
  {
    static Timer timer (string ("Elementmatrix, ") + Name(), 2);
    static Timer timer2 (string ("Elementmatrix, ") + Name() + ", Lapack", 3);
    RegionTimer reg (timer);

    try
      {
	// const FEL & fel = static_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();

	HeapReset hr(lh);

	const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());
	MappedIntegrationRule<DIM_ELEMENT, DIM_SPACE> mir(ir, eltrans, lh);

	FlatMatrixFixHeight<DIM_DMAT, TSCAL> bmat (ndof * DIM, lh);
	Mat<DIM_DMAT,DIM_DMAT,TSCAL> dmat;

	FlatMatrix<TSCAL> bbmat (ndof * DIM, DIM_DMAT*ir.GetNIP(), lh);
	FlatMatrix<TSCAL> bdbmat (ndof * DIM, DIM_DMAT*ir.GetNIP(), lh);

	for (int i = 0; i < ir.GetNIP(); i++)
	  {
	    HeapReset hr(lh);
	    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip = mir[i];

	    DIFFOP::GenerateMatrix (fel, mip, bmat, lh);
	    dmatop.GenerateMatrix (fel, mip, dmat, lh);
	    dmat *= mip.GetWeight();

	    bbmat.Cols(i*DIM_DMAT, (i+1)*DIM_DMAT) = Trans (bmat);
	    bdbmat.Cols(i*DIM_DMAT, (i+1)*DIM_DMAT) = Trans (dmat * bmat);
	  }

	RegionTimer reg2 (timer2);

	if (ndof < 20)
	  {
	    if (DMATOP::SYMMETRIC)
	      elmat = Symmetric ( bdbmat * Trans (bbmat));
	    else
	      elmat = bbmat * Trans (bdbmat); 
	  }
	else
	  elmat = bbmat * Trans(bdbmat) | Lapack;

	timer.AddFlops (long(elmat.Height())*long(elmat.Width())*bbmat.Width());
      } 
    

    catch (Exception & e)
      {
	e.Append (string ("in CalcElementMatrix - lapack, type = ") + typeid(*this).name() + "\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append (string ("\nin CalcElementMatrix - lapack, type = ") + typeid(*this).name() + "\n");
	throw e2;
      }
  }

#endif





  virtual void
  CalcElementMatrixDiag (const FiniteElement & bfel,
			     const ElementTransformation & eltrans, 
			     FlatVector<double> & diag,
			     LocalHeap & lh) const
  {
    try
      {
	const FEL & fel = dynamic_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();

	diag.AssignMemory (ndof*DIM, lh);
	diag = 0.0;

	FlatMatrixFixHeight<DIM_DMAT, double> bmat (ndof * DIM, lh);
	Mat<DIM_DMAT,DIM_DMAT> dmat;

	const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());

	void * heapp = lh.GetPointer();
	for (int i = 0; i < ir.GetNIP(); i++)
	  {
	    MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> mip(ir[i], eltrans);

	    DIFFOP::GenerateMatrix (fel, mip, bmat, lh);
	    dmatop.GenerateMatrix (fel, mip, dmat, lh);

	    double fac =  
	      fabs (mip.GetJacobiDet()) * mip.IP().Weight();

	    for (int j = 0; j < diag.Size(); j++)
	      {
		Vec<DIM_DMAT> hv = dmat * bmat.Col(j);
		diag(j) += fac * InnerProduct (bmat.Col(j), hv);
	      }

	    lh.CleanUp (heapp);
	  } 
      }

    catch (Exception & e)
      {
	e.Append ("in CalcElementMatrixDiag, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin CalcElementMatrixDiag, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }
  }


  /*
  ///
  virtual void 
  ApplyElementMatrix (const FiniteElement & bfel, 
		      const ElementTransformation & eltrans, 
		      const FlatVector<double> & elx, 
		      FlatVector<double> & ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    static Timer timer (string ("BDBIntegrator::ApplyElementMatrix, ") + Name());
    RegionTimer reg (timer);


    const FEL & fel = static_cast<const FEL&> (bfel);
    int ndof = fel.GetNDof ();
    
    ely = 0;
    
    Vec<DIM_DMAT,double> hv1;
    Vec<DIM_DMAT,double> hv2;

    FlatVector<double> hely (ndof*DIM, lh);

    const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());

    MappedIntegrationRule<DIM_ELEMENT, DIM_SPACE> mir(ir, eltrans, lh);

    for (int i = 0; i < ir.GetNIP(); i++)
      {
	HeapReset hr (lh);
	const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip = mir[i];

	DIFFOP::Apply (fel, mip, elx, hv1, lh);
	dmatop.Apply (fel, mip, hv1, hv2, lh);
	DIFFOP::ApplyTrans (fel, mip, hv2, hely, lh);

	double fac = fabs (mip.GetJacobiDet()) * ir[i].Weight();
	ely += fac * hely;
      }     
  }


  ///
  virtual void 
  ApplyElementMatrix (const FiniteElement & bfel, 
		      const ElementTransformation & eltrans, 
		      const FlatVector<Complex> & elx, 
		      FlatVector<Complex> & ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    static Timer timer (string ("BDBIntegrator::ApplyElementMatrix<Complex>, ") + Name(), 2);
    RegionTimer reg (timer);

    const FEL & fel = static_cast<const FEL&> (bfel);
    int ndof = fel.GetNDof ();
    
    ely = 0;
    
    Vec<DIM_DMAT,Complex> hv1;
    Vec<DIM_DMAT,Complex> hv2;
    
    FlatVector<Complex> hely (ndof*DIM, lh);
    const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());
    MappedIntegrationRule<DIM_ELEMENT, DIM_SPACE> mir(ir, eltrans, lh);

    for (int i = 0; i < ir.GetNIP(); i++)
      {
	HeapReset hr (lh);
	const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip = mir[i];

	DIFFOP::Apply (fel, mip, elx, hv1, lh);
	dmatop.Apply (fel, mip, hv1, hv2, lh);
	DIFFOP::ApplyTrans (fel, mip, hv2, hely, lh);

	double fac = fabs (mip.GetJacobiDet()) * mip.IP().Weight();
	ely += fac * hely;
      }     
  }
  */





  virtual void 
  ApplyElementMatrix (const FiniteElement & bfel, 
		      const ElementTransformation & eltrans, 
		      const FlatVector<double> & elx, 
		      FlatVector<double> & ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    T_ApplyElementMatrix<double> (bfel, eltrans, elx, ely, precomputed, lh);
  }

  virtual void 
  ApplyElementMatrix (const FiniteElement & bfel, 
		      const ElementTransformation & eltrans, 
		      const FlatVector<Complex> & elx, 
		      FlatVector<Complex> & ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    T_ApplyElementMatrix<Complex> (bfel, eltrans, elx, ely, precomputed, lh);
  }


  template <typename TSCAL>
  void T_ApplyElementMatrix (const FiniteElement & fel, 
			     const ElementTransformation & eltrans, 
			     const FlatVector<TSCAL> & elx, 
			     FlatVector<TSCAL> & ely,
			     void * precomputed,
			     LocalHeap & lh) const
  {
    static Timer timer (string ("BDBIntegrator::ApplyElementMatrix<") 
			+ typeid(TSCAL).name() + ">" + Name(), 2);
    RegionTimer reg (timer);

    // const FEL & fel = static_cast<const FEL&> (bfel);
    // int ndof = fel.GetNDof ();
    
    /*
    ely = 0;
    
    Vec<DIM_DMAT,TSCAL> hv1;
    Vec<DIM_DMAT,TSCAL> hv2;
    
    FlatVector<TSCAL> hely (ndof*DIM, lh);
    const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());
    MappedIntegrationRule<DIM_ELEMENT, DIM_SPACE> mir(ir, eltrans, lh);

    for (int i = 0; i < ir.GetNIP(); i++)
      {
	HeapReset hr (lh);
	const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip = mir[i];

	DIFFOP::Apply (fel, mip, elx, hv1, lh);
	dmatop.Apply (fel, mip, hv1, hv2, lh);
	DIFFOP::ApplyTrans (fel, mip, hv2, hely, lh);

	ely += mip.GetWeight() * hely;
      }     
    */

    const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());
    MappedIntegrationRule<DIM_ELEMENT, DIM_SPACE> mir(ir, eltrans, lh);
    
    FlatMatrixFixWidth<DIM_DMAT, TSCAL> hv1(ir.GetNIP(), lh);
    FlatMatrixFixWidth<DIM_DMAT, TSCAL> hv2(ir.GetNIP(), lh);

    DIFFOP::ApplyIR (fel, mir, elx, hv1, lh);

    for (int i = 0; i < ir.GetNIP(); i++)
      {
	dmatop.Apply (fel, mir[i], hv1.Row(i), hv2.Row(i), lh);       
        hv2.Row(i) *= mir[i].GetWeight();
      }

    DIFFOP::ApplyTransIR (fel, mir, hv2, ely, lh);    
  }







  ///
  virtual void 
  ApplyMixedElementMatrix (const FiniteElement & bfel1, 
			   const FiniteElement & bfel2, 
			   const ElementTransformation & eltrans, 
			   const FlatVector<double> & elx, 
			   FlatVector<double> & ely,
			   LocalHeap & lh) const
  {
    const FEL & fel1 = dynamic_cast<const FEL&> (bfel1);
    const FEL & fel2 = dynamic_cast<const FEL&> (bfel2);
    // int ndof1 = fel1.GetNDof ();
    int ndof2 = fel2.GetNDof ();
    
    ely = 0;
    
    Vec<DIM_DMAT,double> hv1;
    Vec<DIM_DMAT,double> hv2;
    
    FlatVector<double> hely (ndof2*DIM, lh);

    const IntegrationRule & ir = GetIntegrationRule (fel2,eltrans.HigherIntegrationOrderSet());

    void * heapp = lh.GetPointer();
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	const IntegrationPoint & ip = ir[i];
	
	MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> mip (ir[i], eltrans);

	DIFFOP::Apply (fel1, mip, elx, hv1, lh);
	dmatop.Apply (fel1, mip, hv1, hv2, lh);
	DIFFOP::ApplyTrans (fel2, mip, hv2, hely, lh);

	double fac = fabs (mip.GetJacobiDet()) * ip.Weight();
	ely += fac * hely;
	lh.CleanUp (heapp);
      }     
  }



  IntegrationRule GetIntegrationRule (const FiniteElement & fel, 
				      const bool use_higher_integration_order = false) const
  {
    int order = 2 * fel.Order();

    ELEMENT_TYPE et = fel.ElementType();

    if (et == ET_TET || et == ET_TRIG || et == ET_SEGM)
      order -= 2 * DIFFOP::DIFFORDER;

    if (common_integration_order >= 0)
      order = common_integration_order;

    if (integration_order >= 0)
      order = integration_order;

    if(use_higher_integration_order && higher_integration_order > order)
      {
	//(*testout) << "changing order " << order << " to " << higher_integration_order << endl;
	order = higher_integration_order;
      }
      
    // return SelectIntegrationRule (et, order);
    return IntegrationRule (et, order);
  }










  virtual void
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationPoint & bmip,
	    const FlatVector<double> & elx, 
	    FlatVector<double> & flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    // const FEL & fel = static_cast<const FEL&> (bfel);
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);


    DIFFOP::Apply (fel, mip, elx, flux, lh);
    if (applyd)
      dmatop.Apply1 (fel, mip, flux, lh);
  }

  virtual void
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationRule & bmir,
	    const FlatVector<double> & elx, 
	    FlatMatrix<double> & flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    // const FEL & fel = static_cast<const FEL&> (bfel);
    const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);


    DIFFOP::ApplyIR (fel, mir, elx, flux, lh);
    if (applyd)
      dmatop.ApplyIR (fel, mir, flux, lh);
  }

  virtual void
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationPoint & bmip,
	    const FlatVector<Complex> & elx, 
	    FlatVector<Complex> & flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    // const FEL & fel = static_cast<const FEL&> (bfel);
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);

    if (applyd)
      {
	Vec<DIM_DMAT,Complex> hv1;
	DIFFOP::Apply (fel, mip, elx, hv1, lh);
	dmatop.Apply (fel, mip, hv1, flux, lh);
      }
    else
      {
	DIFFOP::Apply (fel, mip, elx, flux, lh);
      }
  }
  


  


  virtual void
  CalcFluxMulti (const FiniteElement & bfel,
		 const BaseMappedIntegrationPoint & bmip,		
		 int m,
		 const FlatVector<double> & elx, 
		 FlatVector<double> & flux,
		 bool applyd,
		 LocalHeap & lh) const
  {
    const FEL & fel = static_cast<const FEL&> (bfel);
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip = 
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);

    int ndof = fel.GetNDof();
    FlatMatrixFixHeight<DIM_DMAT> bmat (ndof * DIM, lh);

    DIFFOP::GenerateMatrix (fel, mip, bmat, lh);

    if (applyd)
      {
	Vec<DIM_DMAT> hv1;
	Mat<DIM_DMAT,DIM_DMAT> dmat;
	dmatop.GenerateMatrix (fel, mip, dmat, lh);

	for (int i = 0; i < m; i++)
	  {
	    SliceVector<double> slice_x (ndof*DIM, m, &const_cast<double&> (elx(i)));
	    SliceVector<double> slice_flux (DIM_DMAT, m, &flux(i));
	    hv1 = bmat * slice_x;
	    slice_flux = dmat * hv1;
	  }
      }
    else
      {
	for (int i = 0; i < m; i++)
	  {
	    // SliceVector<double> slice_x (ndof*DIM, m, &const_cast<double&> (elx(i)));
            SliceVector<double> slice_x (ndof*DIM, m, &elx(i));
	    SliceVector<double> slice_flux (DIM_DMAT, m, &flux(i));
	    slice_flux = bmat * slice_x;
	  }
      }
  }





  virtual void
  ApplyBTrans (const FiniteElement & fel,
	       const BaseMappedIntegrationPoint & bmip,
	       const FlatVector<double> & elx, 
	       FlatVector<double> & ely,
	       LocalHeap & lh) const
  {
    // const FEL & fel = static_cast<const FEL&> (bfel);
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);

    DIFFOP::ApplyTrans (fel, mip, elx, ely, lh);
  }
  

  virtual void
  ApplyBTrans (const FiniteElement & fel,
	       const BaseMappedIntegrationPoint & bmip,
	       const FlatVector<Complex> & elx, 
	       FlatVector<Complex> & ely,
	       LocalHeap & lh) const
  {
    // const FEL & fel = static_cast<const FEL&> (bfel);
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);

    DIFFOP::ApplyTrans (fel, mip, elx, ely, lh);
  }


  virtual void
  ApplyBTrans (const FiniteElement & fel,
	       const BaseMappedIntegrationRule & bmir,
	       const FlatMatrix<double> & elx, 
	       FlatVector<double> & ely,
	       LocalHeap & lh) const
  {
    // const FEL & fel = static_cast<const FEL&> (bfel);
    const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);

    DIFFOP::ApplyTransIR (fel, mir, elx, ely, lh);
  }

  ///
  virtual int GetDimension () const { return DIM; }
};
























template <class DIFFOP, class DMATOP, class FEL = FiniteElement>
class T_NonlinearBDBIntegrator : public T_BDBIntegrator<DIFFOP, DMATOP, FEL>
{
protected:
  enum { DIM_SPACE   = DIFFOP::DIM_SPACE };
  enum { DIM_ELEMENT = DIFFOP::DIM_ELEMENT };
  enum { DIM_DMAT    = DIFFOP::DIM_DMAT };
  enum { DIM         = DIFFOP::DIM };

  using T_BDBIntegrator<DIFFOP,DMATOP,FEL>::GetIntegrationRule;

public:
  ///
  T_NonlinearBDBIntegrator (const DMATOP & admat)
    : T_BDBIntegrator<DIFFOP,DMATOP,FEL> (admat)
  { ; }

  ///
  virtual ~T_NonlinearBDBIntegrator ()
  { ; }



  virtual void
  CalcLinearizedElementMatrix (const FiniteElement & fel, 
				   const ElementTransformation & eltrans, 
				   FlatVector<double> & elveclin,
				   FlatMatrix<double> & elmat,
				   LocalHeap & lh) const 
  {
    static Timer maintimer (string ("NonlinearBDB, CalcLinearized, ") + this->Name(), 2);
    static Timer bdbtimer ("NonlinearBDB, bdb product", 3);
    RegionTimer reg(maintimer);

    try
      {
	// const FEL & fel = dynamic_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();

	elmat = 0;
	
	FlatMatrixFixHeight<DIM_DMAT, double> bmat (ndof * DIM, lh);
	FlatMatrixFixHeight<DIM_DMAT, double> dbmat (ndof * DIM, lh);
	Vec<DIM_DMAT,double> hvlin;

	Mat<DIM_DMAT,DIM_DMAT> dmat;

	const IntegrationRule & ir = GetIntegrationRule (fel);

	void * heapp = lh.GetPointer();
	for (int i = 0; i < ir.GetNIP(); i++)
	  {
	    MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> mip(ir[i], eltrans);

	    DIFFOP::Apply (fel, mip, elveclin, hvlin, lh);

	    DIFFOP::GenerateMatrix (fel, mip, bmat, lh);

	    this->dmatop . GenerateLinearizedMatrix (fel, mip, hvlin, dmat, lh);

	    double fac =  
	      fabs (mip.GetJacobiDet()) * mip.IP().Weight();

	    {
	      NgProfiler::RegionTimer reg(bdbtimer);
	      dbmat = fac * (dmat * bmat);
	      elmat += Trans (bmat) * dbmat; 
	    }
	    
	    lh.CleanUp (heapp);
	  } 
      }

    catch (Exception & e)
      {
	e.Append ("in CalcLinearizedElementMatrix, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin CalcLinearizedElementMatrix, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }    
  }


  
  virtual void
  CalcLinearizedElementMatrix (const FiniteElement & bfel, 
			       const ElementTransformation & eltrans, 
			       FlatVector<Complex> & elveclin,
			       FlatMatrix<Complex> & elmat,
			       LocalHeap & lh) const 
  {
    static Timer maintimer (string ("NonlinearBDB, CalcLinearized<Complex>, ") + this->Name(), 2);
    static Timer bdbtimer ("NonlinearBDB, bdb product", 3);
    RegionTimer reg(maintimer);

    try
      {
	const FEL & fel = dynamic_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();

	elmat = 0;
	
	FlatMatrixFixHeight<DIM_DMAT, Complex> bmat (ndof * DIM, lh);
	FlatMatrixFixHeight<DIM_DMAT, Complex> dbmat (ndof * DIM, lh);
	Vec<DIM_DMAT,Complex> hvlin;

	Mat<DIM_DMAT,DIM_DMAT, Complex> dmat;

	const IntegrationRule & ir = GetIntegrationRule (fel);
	MappedIntegrationRule<DIM_ELEMENT, DIM_SPACE> mir(ir, eltrans, lh);


	FlatMatrix<Complex> bbmat (ndof * DIM, DIM_DMAT*ir.GetNIP(), lh);
	FlatMatrix<Complex> bdbmat (ndof * DIM, DIM_DMAT*ir.GetNIP(), lh);


	for (int i = 0; i < ir.GetNIP(); i++)
	  {
	    HeapReset hr(lh);

	    MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip = mir[i];

	    DIFFOP::Apply (fel, mip, elveclin, hvlin, lh);

	    DIFFOP::GenerateMatrix (fel, mip, bmat, lh);
	    this->dmatop . GenerateLinearizedMatrix (fel, mip, hvlin, dmat, lh);

	    double fac = mip.GetWeight();
	    /*
	    dbmat = dmat * bmat;
	    elmat += fac * (Trans (bmat) * dbmat); 
	    */

	    dmat *= fac;
	    bbmat.Cols(i*DIM_DMAT, (i+1)*DIM_DMAT) = Trans (bmat);
	    bdbmat.Cols(i*DIM_DMAT, (i+1)*DIM_DMAT) = Trans (dmat * bmat);
	  } 


	RegionTimer reg(bdbtimer);

	if (ndof < 20)
	  {
	    if (DMATOP::SYMMETRIC)
	      elmat = Symmetric ( bdbmat * Trans (bbmat));
	    else
	      elmat = bbmat * Trans (bdbmat); 
	  }
	else
	  // LapackMultABt (bbmat, bdbmat, elmat);
	  LapackMult (bbmat, Trans(bdbmat), elmat);


      }

    catch (Exception & e)
      {
	e.Append ("in CalcLinearizedElementMatrix, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin CalcLinearizedElementMatrix, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }    
  }







  ///
  virtual void 
  ApplyLinearizedElementMatrix (const FiniteElement & fel, 
				const ElementTransformation & eltrans, 
				const FlatVector<double> & ellin, 
				const FlatVector<double> & elx, 
				FlatVector<double> & ely,
				LocalHeap & lh) const
  {
    // const FEL & fel = dynamic_cast<const FEL&> (bfel);
    int ndof = fel.GetNDof ();
    
    ely = 0;
    
    Vec<DIM_DMAT,double> hvlin;
    Vec<DIM_DMAT,double> hvx;
    Vec<DIM_DMAT,double> hvy;
    
    FlatVector<double> hely (ndof*DIM, lh);

    const IntegrationRule & ir = GetIntegrationRule (fel);

    for (int i = 0; i < ir.GetNIP(); i++)
      {
	const IntegrationPoint & ip = ir[i];
	
	MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> mip (ir[i], eltrans);

	DIFFOP::Apply (fel, mip, ellin, hvlin, lh);
	DIFFOP::Apply (fel, mip, elx, hvx, lh);
	this->dmatop.ApplyLinearized (fel, mip, hvlin, hvx, hvy, lh);
	DIFFOP::ApplyTrans (fel, mip, hvy, hely, lh);

	double fac = fabs (mip.GetJacobiDet()) * ip.Weight();
	ely += fac * hely;
      }     
  }



 ///
  virtual void 
  ApplyLinearizedElementMatrix (const FiniteElement & bfel, 
				const ElementTransformation & eltrans, 
				const FlatVector<Complex> & ellin, 
				const FlatVector<Complex> & elx, 
				FlatVector<Complex> & ely,
				LocalHeap & lh) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    int ndof = fel.GetNDof ();
    
    ely = 0;
    
    Vec<DIM_DMAT,Complex> hvlin;
    Vec<DIM_DMAT,Complex> hvx;
    Vec<DIM_DMAT,Complex> hvy;
    
    FlatVector<Complex> hely (ndof*DIM, lh);
    
    /*
    int order = IntegrationOrder (fel);
    const IntegrationRule & ir = 
      GetIntegrationRules().SelectIntegrationRule (fel.ElementType(), order);
    */
    const IntegrationRule & ir = GetIntegrationRule (fel);



    for (int i = 0; i < ir.GetNIP(); i++)
      {
	MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> mip (ir[i], eltrans);

	DIFFOP::Apply (fel, mip, ellin, hvlin, lh);
	DIFFOP::Apply (fel, mip, elx, hvx, lh);
	this->dmatop.ApplyLinearized (fel, mip, hvlin, hvx, hvy, lh);
	DIFFOP::ApplyTrans (fel, mip, hvy, hely, lh);
	
	double fac = fabs (mip.GetJacobiDet()) * mip.IP().Weight();
	ely += fac * hely;
      }     
  }



  virtual double Energy (const FiniteElement & bfel, 
			 const ElementTransformation & eltrans, 
			 const FlatVector<double> & elx, 
			 LocalHeap & lh) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    // int ndof = fel.GetNDof ();
    
    Vec<DIM_DMAT,double> hvx;
    const IntegrationRule & ir = GetIntegrationRule (fel);

    double energy = 0;

    for (int i = 0; i < ir.GetNIP(); i++)
      {
	const IntegrationPoint & ip = ir[i];
	
	MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> mip (ir[i], eltrans);
	DIFFOP::Apply (fel, mip, elx, hvx, lh);

	double fac = fabs (mip.GetJacobiDet()) * ip.Weight();
	energy += fac * this->dmatop.Energy (fel, mip, hvx, lh);
      }     

    return energy;
  }

};


























/**
   Element vector assembling.
   Assembling for linear-forms of type $\int D (B u) dx$.
 */
template <class DIFFOP, class DVecOp, class FEL = FiniteElement>
class T_BIntegrator : public LinearFormIntegrator
{
protected:
  DVecOp dvecop;

public:
  enum { DIM_SPACE = DIFFOP::DIM_SPACE };
  enum { DIM_ELEMENT = DIFFOP::DIM_ELEMENT };
  enum { DIM_DMAT = DIFFOP::DIM_DMAT };
  enum { DIM = DIFFOP::DIM };
  // typedef typename DVecOp::TSCAL TSCAL;

  ///
  T_BIntegrator (const DVecOp & advec)
    : dvecop(advec)
  { ; }

  ///
  virtual ~T_BIntegrator ()
  { ; }

  ///
  virtual void CheckElement (const FiniteElement & el) const
  {
    if (!dynamic_cast<const FEL*> (&el) )
      throw Exception (string ("Element does not match integrator\n") +
                       string ("element type is ") + typeid(el).name() +
                       string (" expected type is ") + typeid(FEL).name() +
                       string ("integrator is ") + Name());
  }

  ///
  virtual bool BoundaryForm () const
  { return int(DIM_SPACE) > int(DIM_ELEMENT); }

  virtual int DimElement () const
  { return DIM_ELEMENT; }

  virtual int DimSpace () const
  { return DIM_SPACE; }



  virtual void
  CalcElementVector (const FiniteElement & bfel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> & elvec,
		     LocalHeap & lh) const
  {
    T_CalcElementVector<double> (bfel, eltrans, elvec, lh);
  }

  virtual void
  CalcElementVector (const FiniteElement & bfel,
		     const ElementTransformation & eltrans, 
		     FlatVector<Complex> & elvec,
		     LocalHeap & lh) const
  {
    T_CalcElementVector<Complex> (bfel, eltrans, elvec, lh);
  }



  ///

  template <typename TSCAL>
  void T_CalcElementVector (const FiniteElement & fel,
			    const ElementTransformation & eltrans, 
			    FlatVector<TSCAL> & elvec,
			    LocalHeap & lh) const
  {
    try
      {
	IntegrationRule ir(fel.ElementType(), IntegrationOrder(fel));
	MappedIntegrationRule<DIM_ELEMENT, DIM_SPACE> mir(ir, eltrans, lh);
	FlatMatrixFixWidth<DIM_DMAT, TSCAL> dvecs(ir.GetNIP(), lh);
	dvecop.GenerateVectorIR (fel, mir, dvecs, lh);
        for (int i = 0; i < ir.GetNIP(); i++)
          dvecs.Row(i) *= mir[i].GetWeight();

        DIFFOP::ApplyTransIR (fel, mir, dvecs, elvec, lh);
      }
    
    catch (Exception & e)
      {
        e.Append ("in CalcElementVector, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
        Exception e2(e.what());
	e2.Append ("\nin CalcElementVector, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }
  }
    




  virtual void
  CalcElementVectorIndependent (const FiniteElement & gfel, 
				    const BaseMappedIntegrationPoint & s_mip,
				    const BaseMappedIntegrationPoint & g_mip,
				    FlatVector<double> & elvec,
				    LocalHeap & lh,
				    const bool curveint = false) const
  {
    T_CalcElementVectorIndependent (gfel, s_mip, g_mip, elvec, lh, curveint);
  }

  virtual void
  CalcElementVectorIndependent (const FiniteElement & gfel, 
				    const BaseMappedIntegrationPoint & s_mip,
				    const BaseMappedIntegrationPoint & g_mip,
				    FlatVector<Complex> & elvec,
				    LocalHeap & lh,
				    const bool curveint = false) const
  {
    T_CalcElementVectorIndependent (gfel, s_mip, g_mip, elvec, lh, curveint);
  }


  template <typename TSCAL>
  void T_CalcElementVectorIndependent (const FiniteElement & gfel, 
					   const BaseMappedIntegrationPoint & s_mip,
					   const BaseMappedIntegrationPoint & g_mip,
					   FlatVector<TSCAL> & elvec,
					   LocalHeap & lh,
					   const bool curveint = false) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (gfel);
    int ndof = fel.GetNDof();
    
    elvec.AssignMemory (ndof * DIM, lh);
    //elvec = 0;

    Vec<DIM_DMAT, TSCAL> dvec;
	
    const MappedIntegrationPoint< DIM_SPACE, DIM_SPACE > & d_g_mip
      (static_cast<const MappedIntegrationPoint< DIM_SPACE, DIM_SPACE > &>(g_mip));

    if(curveint)
      {
	const MappedIntegrationPoint< 1, DIM_SPACE > & d_s_mip
	  (static_cast<const MappedIntegrationPoint< 1, DIM_SPACE > &>(s_mip));

	dvecop.GenerateVector (fel, d_s_mip, dvec, lh);
      }
    else
      {
	enum { HDIM = (DIM_SPACE > 1) ? DIM_SPACE-1 : 1 };
	
	const MappedIntegrationPoint< HDIM, DIM_SPACE > & d_s_mip
	  (static_cast<const MappedIntegrationPoint< HDIM, DIM_SPACE > &>(s_mip));
	
	dvecop.GenerateVector (fel, d_s_mip, dvec, lh);
      }


    DIFFOP::ApplyTrans (fel, d_g_mip, dvec, elvec, lh);
  
    //(*testout) << "dvec " << dvec << " elvec " << elvec << endl;

  }  




  int IntegrationOrder (const FiniteElement & fel) const
  {
    // int order = fel.Order()+2;    // low order case
    int order = 2*fel.Order()+1;  // high order case
    
    ELEMENT_TYPE et = fel.ElementType();

    if (et == ET_TET || et == ET_TRIG || et == ET_SEGM)
      order -= DIFFOP::DIFFORDER;

    /*    
    if (common_integration_order >= 0)
      order = common_integration_order;
    */
    if (integration_order >= 0)
      order = integration_order;

    return order;
  }

  
  ///
  virtual int GetDimension () const { return DIM; }

  ///
  virtual string Name () const { return "BDB integrator"; }
};



}


#endif
