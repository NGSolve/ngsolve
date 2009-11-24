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
template <class DMO>
class DMatOp
{
public:
  typedef double TSCAL;
  /// is coefficient tensor symmetric ?
  enum { SYMMETRIC = 1 };

  /// generate linearized matrix in linearization point vec
  template <typename FEL, typename SIP, typename VEC, typename MAT>
  void GenerateLinearizedMatrix (const FEL & fel, const SIP & sip, VEC & vec,
				 MAT & mat, LocalHeap & lh) const
  {
    static_cast<const DMO*>(this) -> GenerateMatrix (fel, sip, mat, lh);
  }


  /// apply coefficient matrix.
  template <typename FEL, typename SIP, class TVX, class TVY>
  void Apply (const FEL & fel, const SIP & sip,
	      const TVX & x, TVY & y,
	      LocalHeap & lh) const
  {
    Mat<DMO::DIM_DMAT, DMO::DIM_DMAT, double> mat;
    static_cast<const DMO*>(this) -> GenerateMatrix (fel, sip, mat, lh);
    y = mat * x;
  }

  /// apply transpose coefficient tensor
  template <typename FEL, typename SIP, class TVX, class TVY>
  void ApplyTrans (const FEL & fel, const SIP & sip,
		   const TVX & x, TVY & y,
		   LocalHeap & lh) const
  {
    Mat<DMO::DIM_DMAT, DMO::DIM_DMAT, double> mat;
    static_cast<const DMO*>(this) -> GenerateMatrix (fel, sip, mat, lh);
    y = Trans (mat) * x;
  }

  /// computes energy 
  template <typename FEL, typename SIP, class TVX>
  double Energy (const FEL & fel, const SIP & sip,
		 const TVX & x, LocalHeap & lh) const  
  {
    TVX y;
    static_cast<const DMO*>(this) -> Apply (fel, sip, x, y, lh);
    return 0.5 * InnerProduct (x,y);
  }
};






#ifdef WIN32
#define __restrict__ __restrict
#endif

template <int M>
void FastMat (int n, Complex * ba, Complex *  pb, Complex * pc);

template <int M>
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
  typedef typename DMATOP::TSCAL TSCAL;

  /*
  T_BDBIntegrator  (Array<CoefficientFunction*> & coeffs)
  : dmatop(coeffs)
  { ; }
  */

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
			  const BaseSpecificIntegrationPoint & bsip,
			  const FlatVector<double> & elx, 
			  FlatVector<double> & eldx,
			  LocalHeap & lh) const
  {
    eldx.AssignMemory (DIM_DMAT, lh);
    dmatop.Apply(dynamic_cast<const FEL&> (bfel),
		 static_cast<const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> &>(bsip),
		 elx, eldx ,lh);
  }

  virtual void ApplyDMat (const FiniteElement & bfel,
			  const BaseSpecificIntegrationPoint & bsip,
			  const FlatVector<Complex> & elx, 
			  FlatVector<Complex> & eldx,
			  LocalHeap & lh) const
  {
    eldx.AssignMemory (DIM_DMAT, lh);
    dmatop.Apply(dynamic_cast<const FEL&> (bfel),
		 static_cast<const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> &>(bsip),
		 elx,eldx,lh);
  }



#ifdef TEXT_BOOK_VERSION

  virtual void
  AssembleElementMatrix (const FiniteElement & bfel,
			 const ElementTransformation & eltrans, 
			 FlatMatrix<double> & elmat,
			 LocalHeap & locheap) const
  {

    try
      {
	const FEL & fel = static_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();
        
	elmat = 0;

	FlatMatrixFixHeight<DIM_DMAT, double> bmat (ndof * DIM, locheap);
	FlatMatrixFixHeight<DIM_DMAT, double> dbmat (ndof * DIM, locheap);
	Mat<DIM_DMAT,DIM_DMAT> dmat;

	const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());

	for (int i = 0 ; i < ir.GetNIP(); i++)
	  {
	    HeapReset hr(locheap);

	    SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> 
	      sip(ir[i], eltrans, locheap);

	    dmatop.GenerateMatrix (fel, sip, dmat, locheap);

	    DIFFOP::GenerateMatrix (fel, sip, bmat, locheap);
	    double fac =  
	      fabs (sip.GetJacobiDet()) * sip.IP().Weight();
	    
	    dbmat = fac * (dmat * bmat);
	    elmat += Trans (bmat) * dbmat; 
	  } 
      }

    catch (Exception & e)
      {
	e.Append ("in AssembleElementMatrix - textbook, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin AssembleElementMatrix - textbook, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }
  }
#endif






#define BLOCK_VERSION
#ifdef BLOCK_VERSION


  /// this is my preferred one !


  virtual void
  AssembleElementMatrix (const FiniteElement & bfel,
			 const ElementTransformation & eltrans, 
			 FlatMatrix<double> & elmat,
			 LocalHeap & lh) const
  {
    // static int timer = NgProfiler::CreateTimer (string ("Elementmatrix, ") + Name());
    // NgProfiler::RegionTimer reg (timer);

    try
      {
        const FEL & fel = *static_cast<const FEL*> (&bfel);
	int ndof = fel.GetNDof();
        
	elmat = 0;
	

        enum { BLOCK = 2 * (12 / DIM_DMAT + 1) };
	
        HeapReset hr1(lh);
	
        FlatMatrixFixHeight<DIM_DMAT, double> bmat (ndof * DIM, lh);
        FlatMatrixFixHeight<DIM_DMAT, double> dbmat (ndof * DIM, lh);
	
        FlatMatrixFixHeight<DIM_DMAT*BLOCK, double> bbmat (ndof * DIM, lh);
        FlatMatrixFixHeight<DIM_DMAT*BLOCK, double> bdbmat (ndof * DIM, lh);
        Mat<DIM_DMAT,DIM_DMAT> dmat;
	
        const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());

	/*
        FlatArray<Vec<DIM_SPACE> > pts(ir.GetNIP(), lh);
        FlatArray<Mat<DIM_SPACE, DIM_ELEMENT> > dxdxi(ir.GetNIP(), lh);

        eltrans.CalcMultiPointJacobian (ir, pts, dxdxi, lh);
	*/
	MappedIntegrationRule<DIM_ELEMENT, DIM_SPACE> mir(ir, eltrans, lh);
	
        int i = 0;
        for (int i1 = 0; i1 < ir.GetNIP() / BLOCK; i1++)
          {
            for (int i2 = 0; i2 < BLOCK; i++, i2++)
              {
                HeapReset hr(lh);
		/*
                SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> 
                  sip(ir[i], eltrans, pts[i], dxdxi[i]);
		*/

                const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & sip = mir[i];


                DIFFOP::GenerateMatrix (fel, sip, bmat, lh);
                dmatop.GenerateMatrix (fel, sip, dmat, lh);
                double fac =  
                  fabs (sip.GetJacobiDet()) * sip.IP().Weight();
                
                dbmat = fac * (dmat * bmat);
		
                for (int l = 0; l < ndof*DIM; l++)
                  for (int k = 0; k < DIM_DMAT; k++)
                    {
                      bbmat(i2*DIM_DMAT+k, l) = bmat(k,l);
                      bdbmat(i2*DIM_DMAT+k, l) = dbmat(k,l);
                    }
              }
            
            if (DMATOP::SYMMETRIC)
              // elmat += Symmetric (Trans (bdbmat) * bbmat);
              FastMat<DIM_DMAT*BLOCK> (elmat.Height(), &bdbmat(0,0), &bbmat(0,0), &elmat(0,0));
            else
              elmat += Trans (bbmat) * bdbmat; 
          } 
        
        for ( ; i < ir.GetNIP(); i++)
          {
            HeapReset hr(lh);
            
	    const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & sip = mir[i];
	    /*
            SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE>
	      sip(ir[i], eltrans, pts[i], dxdxi[i]);
	    */

            DIFFOP::GenerateMatrix (fel, sip, bmat, lh);
            dmatop.GenerateMatrix (fel, sip, dmat, lh);
            double fac =  
              fabs (sip.GetJacobiDet()) * sip.IP().Weight();
            
            dbmat = fac * (dmat * bmat);
            
            
            if (DMATOP::SYMMETRIC)
              // elmat += Symmetric (Trans (dbmat) * bmat);
              FastMat<DIM_DMAT> (elmat.Height(), &dbmat(0,0), &bmat(0,0), &elmat(0,0));
            else
              elmat += Trans (bmat) * dbmat;
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
	e.Append ("in AssembleElementMatrix - blockversion, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin AssembleElementMatrix - blockversion, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }
  }



#else // blockversion
  // one matrix matrix multiplication
  ///
  virtual void
  AssembleElementMatrix (const FiniteElement & bfel,
			 const ElementTransformation & eltrans, 
			 FlatMatrix<double> & elmat,
			 LocalHeap & locheap) const
  {
    static int timer = NgProfiler::CreateTimer (string ("Elementmatrix, ") + Name());
    NgProfiler::RegionTimer reg (timer);

    try
      {
	// cout << "start assemble, free = " << locheap.Available() << endl;

	const FEL & fel = dynamic_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();

	const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());

	// elmat.AssignMemory (ndof*DIM, ndof*DIM, locheap);
	// elmat = 0;
	
	FlatMatrixFixHeight<DIM_DMAT, double> bmat (ndof * DIM, locheap);
	FlatMatrixFixHeight<DIM_DMAT, double> dbmat (ndof * DIM, locheap);
	Mat<DIM_DMAT,DIM_DMAT> dmat;

	FlatMatrix<double> bbmat (ndof * DIM, DIM_DMAT*ir.GetNIP(), locheap);
	FlatMatrix<double> bdbmat (ndof * DIM, DIM_DMAT*ir.GetNIP(), locheap);

	void * heapp = locheap.GetPointer();
	for (int i = 0; i < ir.GetNIP(); i++)
	  {
	    // cout << "assemble, ip loop, free = " << locheap.Available() << endl;

	    SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> 
	      sip(ir[i], eltrans, locheap);
            
	    DIFFOP::GenerateMatrix (fel, sip, bmat, locheap);
	    dmatop.GenerateMatrix (fel, sip, dmat, locheap);
	    double fac =  
	      fabs (sip.GetJacobiDet()) * sip.IP().Weight();

	    dmat *= fac;
	    dbmat = dmat * bmat;
	    // dbmat = fac * (dmat * bmat);

	    /*
	    for (int l = 0; l < ndof*DIM; l++)
	      for (int k = 0; k < DIM_DMAT; k++)
		{
		  bbmat(l, i*DIM_DMAT+k) = bmat(k,l);
		  bdbmat(l, i*DIM_DMAT+k) = dbmat(k,l);
		}
	    */
	    for (int k = 0; k < DIM_DMAT; k++)
	      bbmat.Col(i*DIM_DMAT+k) = bmat.Row(k);

	    for (int k = 0; k < DIM_DMAT; k++)
	      bdbmat.Col(i*DIM_DMAT+k) = dbmat.Row(k);


	    locheap.CleanUp (heapp);
	  }


	/*
	if (DMATOP::SYMMETRIC)
	  elmat += Symmetric ( bdbmat * Trans (bbmat));
	else
	  elmat += bbmat * Trans (bdbmat); 
	*/

	LapackMultABt (bbmat, bdbmat, elmat);
	NgProfiler::AddFlops (timer, long(elmat.Height())*long(elmat.Width())*bbmat.Width());

	/*
	cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasTrans,
		     ndof*DIM, ndof*DIM, bbmat.Width(), 1.0,
		     &bdbmat(0,0), bbmat.Width(), 
		     &bbmat(0,0), bbmat.Width(),
		     1.0, &elmat(0,0), ndof*DIM);
	*/
      } 
    

    catch (Exception & e)
      {
	e.Append ("in AssembleElementMatrix - lapack, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin AssembleElementMatrix - lapack, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }
  }

#endif





  virtual void
  AssembleElementMatrixDiag (const FiniteElement & bfel,
			     const ElementTransformation & eltrans, 
			     FlatVector<double> & diag,
			     LocalHeap & locheap) const
  {
    try
      {
	const FEL & fel = dynamic_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();

	diag.AssignMemory (ndof*DIM, locheap);
	diag = 0.0;

	FlatMatrixFixHeight<DIM_DMAT, double> bmat (ndof * DIM, locheap);
	Mat<DIM_DMAT,DIM_DMAT> dmat;

	const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());

	void * heapp = locheap.GetPointer();
	for (int i = 0; i < ir.GetNIP(); i++)
	  {
	    SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> 
	      sip(ir[i], eltrans, locheap);

	    DIFFOP::GenerateMatrix (fel, sip, bmat, locheap);
	    dmatop.GenerateMatrix (fel, sip, dmat, locheap);

	    double fac =  
	      fabs (sip.GetJacobiDet()) * sip.IP().Weight();

	    for (int j = 0; j < diag.Size(); j++)
	      {
		Vec<DIM_DMAT> hv = dmat * bmat.Col(j);
		diag(j) += fac * InnerProduct (bmat.Col(j), hv);
	      }

	    locheap.CleanUp (heapp);
	  } 
      }

    catch (Exception & e)
      {
	e.Append ("in AssembleElementMatrixDiag, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin AssembleElementMatrixDiag, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }
  }


  ///
  virtual void 
  ApplyElementMatrix (const FiniteElement & bfel, 
		      const ElementTransformation & eltrans, 
		      const FlatVector<double> & elx, 
		      FlatVector<double> & ely,
		      void * precomputed,
		      LocalHeap & lh) const
  {
    static int maintimer = NgProfiler::CreateTimer ("BDBIntegrator::ApplyElementMatrix");
    NgProfiler::RegionTimer reg(maintimer);

    const FEL & fel = static_cast<const FEL&> (bfel);
    int ndof = fel.GetNDof ();
    
    ely = 0;
    
    Vec<DIM_DMAT,double> hv1;
    Vec<DIM_DMAT,double> hv2;

    FlatVector<double> hely (ndof*DIM, lh);

    const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());


    FlatArray<Vec<DIM_SPACE> > pts(ir.GetNIP(), lh);
    FlatArray<Mat<DIM_SPACE, DIM_ELEMENT> > dxdxi(ir.GetNIP(), lh);
    
    eltrans.CalcMultiPointJacobian (ir, pts, dxdxi, lh);


    for (int i = 0; i < ir.GetNIP(); i++)
      {
	HeapReset hr (lh);

	// const IntegrationPoint & ip = ir[i];
	
	// SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip (ir[i], eltrans, lh);
	SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> 
	  sip(ir[i], eltrans, pts[i], dxdxi[i]);

	DIFFOP::Apply (fel, sip, elx, hv1, lh);
	dmatop.Apply (fel, sip, hv1, hv2, lh);
	DIFFOP::ApplyTrans (fel, sip, hv2, hely, lh);

	double fac = fabs (sip.GetJacobiDet()) * ir[i].Weight();
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
		      LocalHeap & locheap) const
  {
    const FEL & fel = static_cast<const FEL&> (bfel);
    int ndof = fel.GetNDof ();
    
    ely = 0;
    
    Vec<DIM_DMAT,Complex> hv1;
    Vec<DIM_DMAT,Complex> hv2;
    
    FlatVector<Complex> hely (ndof*DIM, locheap);
    const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());


    for (int i = 0; i < ir.GetNIP(); i++)
      {
	SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip (ir[i], eltrans, locheap);

	DIFFOP::Apply (fel, sip, elx, hv1, locheap);
	dmatop.Apply (fel, sip, hv1, hv2, locheap);
	DIFFOP::ApplyTrans (fel, sip, hv2, hely, locheap);

	double fac = fabs (sip.GetJacobiDet()) * sip.IP().Weight();
	ely += fac * hely;
      }     
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
	
	SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip (ir[i], eltrans, lh);

	DIFFOP::Apply (fel1, sip, elx, hv1, lh);
	dmatop.Apply (fel1, sip, hv1, hv2, lh);
	DIFFOP::ApplyTrans (fel2, sip, hv2, hely, lh);

	double fac = fabs (sip.GetJacobiDet()) * ip.Weight();
	ely += fac * hely;
	lh.CleanUp (heapp);
      }     
  }



  virtual const IntegrationRule & GetIntegrationRule (const FiniteElement & fel, const bool use_higher_integration_order = false) const
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
      
    return SelectIntegrationRule (et, order);
  }

  


  virtual void
  CalcFlux (const FiniteElement & bfel,
	    const ElementTransformation & eltrans,
	    const IntegrationPoint & ip,
	    const FlatVector<double> & elx, 
	    FlatVector<double> & flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip (ip, eltrans, lh);
    CalcFlux (bfel, sip, elx, flux, applyd, lh);
  }


  virtual void
  CalcFlux (const FiniteElement & bfel,
	    const ElementTransformation & eltrans,
	    const IntegrationPoint & ip,
	    const FlatVector<Complex> & elx, 
	    FlatVector<Complex> & flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip (ip, eltrans, lh);
    CalcFlux (bfel, sip, elx, flux, applyd, lh);
  }


  virtual void
  CalcFlux (const FiniteElement & bfel,
	    const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & sip,
	    const FlatVector<double> & elx, 
	    FlatVector<double> & flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    
    // flux.AssignMemory (DIM_DMAT, lh);
    if (applyd)
      {
	Vec<DIM_DMAT,double> hv1;
	DIFFOP::Apply (fel, sip, elx, hv1, lh);
	dmatop.Apply (fel, sip, hv1, flux, lh);
      }
    else
      {
	DIFFOP::Apply (fel, sip, elx, flux, lh);
      }
  }


  virtual void
  CalcFluxMulti (const FiniteElement & bfel,
		 const BaseSpecificIntegrationPoint & bsip,		
		 int m,
		 const FlatVector<double> & elx, 
		 FlatVector<double> & flux,
		 bool applyd,
		 LocalHeap & lh) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & sip = 
      static_cast<const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bsip);

    // flux.AssignMemory (DIM_DMAT * m, lh);

    int ndof = fel.GetNDof();
    FlatMatrixFixHeight<DIM_DMAT> bmat (ndof * DIM, lh);

    DIFFOP::GenerateMatrix (fel, sip, bmat, lh);

    if (applyd)
      {
	Vec<DIM_DMAT> hv1;
	Mat<DIM_DMAT,DIM_DMAT> dmat;
	dmatop.GenerateMatrix (fel, sip, dmat, lh);

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
  CalcFlux (const FiniteElement & bfel,
	    const BaseSpecificIntegrationPoint & bsip,
	    const FlatVector<double> & elx, 
	    FlatVector<double> & flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & sip =
      static_cast<const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bsip);

    
    // flux.AssignMemory (DIM_DMAT, lh);
    if (applyd)
      {
	Vec<DIM_DMAT,double> hv1;
	DIFFOP::Apply (fel, sip, elx, hv1, lh);
	dmatop.Apply (fel, sip, hv1, flux, lh);
      }
    else
      {
	DIFFOP::Apply (fel, sip, elx, flux, lh);
      }
  }
  

  virtual void
  CalcFlux (const FiniteElement & bfel,
	    const BaseMappedIntegrationRule & bmir,
	    const FlatVector<double> & elx, 
	    FlatMatrix<double> & flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);


    if (applyd)
      {
	Vec<DIM_DMAT,double> hv1, hv2;
	for (int i = 0; i < mir.Size(); i++)
	  {
	    DIFFOP::Apply (fel, mir[i], elx, hv1, lh);
	    dmatop.Apply (fel, mir[i], hv1, hv2, lh);
	    flux.Row(i) = hv2;
	  }
      }
    else
      {
	DIFFOP::ApplyIR (fel, mir, elx, flux, lh);
	/*
	Vec<DIM_DMAT,double> hv1;
	for (int i = 0; i < mir.Size(); i++)
	  {
	    DIFFOP::Apply (fel, mir[i], elx, hv1, lh);
	    flux.Row(i) = hv1;
	  }
	*/
      }
  }
  






  virtual void
  CalcFlux (const FiniteElement & bfel,
	    const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & sip,
	    const FlatVector<Complex> & elx, 
	    FlatVector<Complex> & flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    // flux.AssignMemory (DIM_DMAT, lh);
    if (applyd)
      {
	Vec<DIM_DMAT,Complex> hv1;
	DIFFOP::Apply (fel, sip, elx, hv1, lh);
	dmatop.Apply (fel, sip, hv1, flux, lh);
      }
    else
      {
	DIFFOP::Apply (fel, sip, elx, flux, lh);
      }
  }
  

  virtual void
  CalcFlux (const FiniteElement & bfel,
	    const BaseSpecificIntegrationPoint & bsip,
	    const FlatVector<Complex> & elx, 
	    FlatVector<Complex> & flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip =
      static_cast<const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bsip);
    // flux.AssignMemory (DIM_DMAT, lh);
    if (applyd)
      {
	Vec<DIM_DMAT,Complex> hv1;
	DIFFOP::Apply (fel, sip, elx, hv1, lh);
	dmatop.Apply (fel, sip, hv1, flux, lh);
      }
    else
      {
	DIFFOP::Apply (fel, sip, elx, flux, lh);
      }
  }
  



  virtual void
  ApplyBTrans (const FiniteElement & bfel,
	       const BaseSpecificIntegrationPoint & bsip,
	       const FlatVector<double> & elx, 
	       FlatVector<double> & ely,
	       LocalHeap & lh) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip =
      static_cast<const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bsip);
    DIFFOP::ApplyTrans (fel, sip, elx, ely, lh);
    // ely *= sip.IP().Weight() * fabs (sip.GetJacobiDet());
  }
  

  virtual void
  ApplyBTrans (const FiniteElement & bfel,
	       const BaseSpecificIntegrationPoint & bsip,
	       const FlatVector<Complex> & elx, 
	       FlatVector<Complex> & ely,
	       LocalHeap & lh) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip =
      static_cast<const SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bsip);
    DIFFOP::ApplyTrans (fel, sip, elx, ely, lh);
    // ely *= sip.IP().Weight() * fabs (sip.GetJacobiDet());
  }
  


  ///
  virtual int GetDimension () const
  { return DIM; }

  ///
  virtual int Lumping () const
    { return 0; }
  ///
  /*
  virtual string Name () const 
  {
    return "BDB integrator"; 
  }
  */
};





























template <class DIFFOP, class DMATOP, class FEL = FiniteElement>
class T_NonlinearBDBIntegrator : public T_BDBIntegrator<DIFFOP, DMATOP, FEL>
{
protected:
  enum { DIM_SPACE   = DIFFOP::DIM_SPACE };
  enum { DIM_ELEMENT = DIFFOP::DIM_ELEMENT };
  enum { DIM_DMAT    = DIFFOP::DIM_DMAT };
  enum { DIM         = DIFFOP::DIM };

public:
  ///
  T_NonlinearBDBIntegrator (const DMATOP & admat)
    : T_BDBIntegrator<DIFFOP,DMATOP,FEL> (admat)
  { ; }

  ///
  virtual ~T_NonlinearBDBIntegrator ()
  { ; }



  virtual void
  AssembleLinearizedElementMatrix (const FiniteElement & bfel, 
				   const ElementTransformation & eltrans, 
				   FlatVector<double> & elveclin,
				   FlatMatrix<double> & elmat,
				   LocalHeap & locheap) const 
  {
    static int maintimer = NgProfiler::CreateTimer ("NonlinearBDB, Assemblelinearized");
    static int bdbtimer = NgProfiler::CreateTimer ("NonlinearBDB, bdb product");
    NgProfiler::RegionTimer reg(maintimer);

    try
      {
	const FEL & fel = dynamic_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();

	elmat = 0;
	
	FlatMatrixFixHeight<DIM_DMAT, double> bmat (ndof * DIM, locheap);
	FlatMatrixFixHeight<DIM_DMAT, double> dbmat (ndof * DIM, locheap);
	Vec<DIM_DMAT,double> hvlin;

	Mat<DIM_DMAT,DIM_DMAT> dmat;

	const IntegrationRule & ir = GetIntegrationRule (fel);

	void * heapp = locheap.GetPointer();
	for (int i = 0; i < ir.GetNIP(); i++)
	  {
	    SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> 
	      sip(ir[i], eltrans, locheap);

	    DIFFOP::Apply (fel, sip, elveclin, hvlin, locheap);

	    DIFFOP::GenerateMatrix (fel, sip, bmat, locheap);

	    this->dmatop . GenerateLinearizedMatrix (fel, sip, hvlin, dmat, locheap);

	    double fac =  
	      fabs (sip.GetJacobiDet()) * sip.IP().Weight();

	    {
	      NgProfiler::RegionTimer reg(bdbtimer);
	      dbmat = fac * (dmat * bmat);
	      elmat += Trans (bmat) * dbmat; 
	    }

	    locheap.CleanUp (heapp);
	  } 
      }

    catch (Exception & e)
      {
	e.Append ("in AssembleLinearizedElementMatrix, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin AssembleLinearizedElementMatrix, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }    
  }



  virtual void
  AssembleLinearizedElementMatrix (const FiniteElement & bfel, 
				   const ElementTransformation & eltrans, 
				   FlatVector<Complex> & elveclin,
				   FlatMatrix<Complex> & elmat,
				   LocalHeap & locheap) const 
  {
    try
      {
	const FEL & fel = dynamic_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();

	elmat = 0;
	
	FlatMatrixFixHeight<DIM_DMAT, Complex> bmat (ndof * DIM, locheap);
	FlatMatrixFixHeight<DIM_DMAT, Complex> dbmat (ndof * DIM, locheap);
	Vec<DIM_DMAT,Complex> hvlin;

	Mat<DIM_DMAT,DIM_DMAT, Complex> dmat;

	const IntegrationRule & ir = GetIntegrationRule (fel);

	void * heapp = locheap.GetPointer();
	for (int i = 0; i < ir.GetNIP(); i++)
	  {
	    SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> 
	      sip(ir[i], eltrans, locheap);

	    DIFFOP::Apply (fel, sip, elveclin, hvlin, locheap);

	    DIFFOP::GenerateMatrix (fel, sip, bmat, locheap);
	    this->dmatop . GenerateLinearizedMatrix (fel, sip, hvlin, dmat, locheap);

	    double fac =  
	      fabs (sip.GetJacobiDet()) * sip.IP().Weight();

	    dbmat = dmat * bmat;
	    elmat += fac * (Trans (bmat) * dbmat); 
	    
	    locheap.CleanUp (heapp);
	  } 
      }

    catch (Exception & e)
      {
	e.Append ("in AssembleLinearizedElementMatrix, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
	Exception e2(e.what());
	e2.Append ("\nin AssembleLinearizedElementMatrix, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }    
  }







  ///
  virtual void 
  ApplyLinearizedElementMatrix (const FiniteElement & bfel, 
				const ElementTransformation & eltrans, 
				const FlatVector<double> & ellin, 
				const FlatVector<double> & elx, 
				FlatVector<double> & ely,
				LocalHeap & locheap) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    int ndof = fel.GetNDof ();
    
    ely = 0;
    
    Vec<DIM_DMAT,double> hvlin;
    Vec<DIM_DMAT,double> hvx;
    Vec<DIM_DMAT,double> hvy;
    
    FlatVector<double> hely (ndof*DIM, locheap);

    /*    
    int order = IntegrationOrder (fel);
    const IntegrationRule & ir = 
      GetIntegrationRules().SelectIntegrationRule (fel.ElementType(), order);
    */
    const IntegrationRule & ir = GetIntegrationRule (fel);


    for (int i = 0; i < ir.GetNIP(); i++)
      {
	const IntegrationPoint & ip = ir[i];
	
	SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip (ir[i], eltrans, locheap);

	DIFFOP::Apply (fel, sip, ellin, hvlin, locheap);
	DIFFOP::Apply (fel, sip, elx, hvx, locheap);
	this->dmatop.ApplyLinearized (fel, sip, hvlin, hvx, hvy, locheap);
	DIFFOP::ApplyTrans (fel, sip, hvy, hely, locheap);

	double fac = fabs (sip.GetJacobiDet()) * ip.Weight();
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
				LocalHeap & locheap) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    int ndof = fel.GetNDof ();
    
    ely = 0;
    
    Vec<DIM_DMAT,Complex> hvlin;
    Vec<DIM_DMAT,Complex> hvx;
    Vec<DIM_DMAT,Complex> hvy;
    
    FlatVector<Complex> hely (ndof*DIM, locheap);
    
    /*
    int order = IntegrationOrder (fel);
    const IntegrationRule & ir = 
      GetIntegrationRules().SelectIntegrationRule (fel.ElementType(), order);
    */
    const IntegrationRule & ir = GetIntegrationRule (fel);



    for (int i = 0; i < ir.GetNIP(); i++)
      {
	SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip (ir[i], eltrans, locheap);

	DIFFOP::Apply (fel, sip, ellin, hvlin, locheap);
	DIFFOP::Apply (fel, sip, elx, hvx, locheap);
	this->dmatop.ApplyLinearized (fel, sip, hvlin, hvx, hvy, locheap);
	DIFFOP::ApplyTrans (fel, sip, hvy, hely, locheap);
	
	double fac = fabs (sip.GetJacobiDet()) * sip.IP().Weight();
	ely += fac * hely;
      }     
  }



  virtual double Energy (const FiniteElement & bfel, 
			 const ElementTransformation & eltrans, 
			 const FlatVector<double> & elx, 
			 LocalHeap & locheap) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    int ndof = fel.GetNDof ();
    
    Vec<DIM_DMAT,double> hvx;
    const IntegrationRule & ir = GetIntegrationRule (fel);

    double energy = 0;

    for (int i = 0; i < ir.GetNIP(); i++)
      {
	const IntegrationPoint & ip = ir[i];
	
	SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip (ir[i], eltrans, locheap);
	DIFFOP::Apply (fel, sip, elx, hvx, locheap);

	double fac = fabs (sip.GetJacobiDet()) * ip.Weight();
	energy += fac * this->dmatop.Energy (fel, sip, hvx, locheap);
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
  typedef typename DVecOp::TSCAL TSCAL;

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



  ///
  virtual void
  AssembleElementVector (const FiniteElement & bfel,
			 const ElementTransformation & eltrans, 
			 FlatVector<TSCAL> & elvec,
			 LocalHeap & lh) const
  {
    try
      {
	const FEL & fel = static_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();
	
	elvec = 0;

	FlatVector<TSCAL> hv(ndof * DIM, lh);
	Vec<DIM_DMAT, TSCAL> dvec;
	
	int order = IntegrationOrder (fel);
	const IntegrationRule & ir = SelectIntegrationRule (fel.ElementType(), order);
	
        FlatArray<Vec<DIM_SPACE> > pts(ir.GetNIP(), lh);
        FlatArray<Mat<DIM_SPACE, DIM_ELEMENT> > dxdxi(ir.GetNIP(), lh);


        eltrans.CalcMultiPointJacobian (ir, pts, dxdxi, lh);

	for (int i = 0; i < ir.GetNIP(); i++)
	  {
            HeapReset hr(lh);
            SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE> sip(ir[i], eltrans, pts[i], dxdxi[i]); // , lh);

	    dvecop.GenerateVector (fel, sip, dvec, lh);
	    DIFFOP::ApplyTrans (fel, sip, dvec, hv, lh);

	    double fac = fabs (sip.GetJacobiDet()) * sip.IP().Weight();
	    elvec += fac * hv;
	  }
      }
    
    catch (Exception & e)
      {
        e.Append ("in AssembleElementVector, type = ");
	e.Append (typeid(*this).name());
	e.Append ("\n");
	throw;
      }
    catch (exception & e)
      {
        Exception e2(e.what());
	e2.Append ("\nin AssembleElementVector, type = ");
	e2.Append (typeid(*this).name());
	e2.Append ("\n");
	throw e2;
      }
  }
    



  virtual void
  AssembleElementVectorIndependent (const FiniteElement & gfel, 
				    const BaseSpecificIntegrationPoint & s_sip,
				    const BaseSpecificIntegrationPoint & g_sip,
				    FlatVector<TSCAL> & elvec,
				    LocalHeap & locheap,
				    const bool curveint = false) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (gfel);
    int ndof = fel.GetNDof();
    
    elvec.AssignMemory (ndof * DIM, locheap);
    //elvec = 0;

    Vec<DIM_DMAT, TSCAL> dvec;
	
    const SpecificIntegrationPoint< DIM_SPACE, DIM_SPACE > & d_g_sip
      (static_cast<const SpecificIntegrationPoint< DIM_SPACE, DIM_SPACE > &>(g_sip));

    if(curveint)
      {
	const SpecificIntegrationPoint< 1, DIM_SPACE > & d_s_sip
	  (static_cast<const SpecificIntegrationPoint< 1, DIM_SPACE > &>(s_sip));

	dvecop.GenerateVector (fel, d_s_sip, dvec, locheap);
      }
    else
      {
	enum { HDIM = (DIM_SPACE > 1) ? DIM_SPACE-1 : 1 };
	
	const SpecificIntegrationPoint< HDIM, DIM_SPACE > & d_s_sip
	  (static_cast<const SpecificIntegrationPoint< HDIM, DIM_SPACE > &>(s_sip));
	
	dvecop.GenerateVector (fel, d_s_sip, dvec, locheap);
      }


    DIFFOP::ApplyTrans (fel, d_g_sip, dvec, elvec, locheap);
  
    //(*testout) << "dvec " << dvec << " elvec " << elvec << endl;

  }  




  int IntegrationOrder (const FEL & fel) const
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
  virtual int GetDimension () const
  { return DIM; }

  ///
  virtual int Lumping () const
    { return 0; }

  ///
  virtual string Name () const { return "BDB integrator"; }
};



}


#endif
