#ifndef FILE_PML
#define FILE_PML

/*********************************************************************/
/* File:   pml.hpp                                                   */
/* Author: Joachim Schoeberl                                         */
/* Date:   1. Jul. 2004                                              */
/*********************************************************************/

namespace ngfem
{

  extern Complex alpha;
  extern double pml_r;
  extern double pml_x;

  extern double pml_xmin[3];
  extern double pml_xmax[3]; 

  /*
    rect_pml = 0 .... circular pml with radius pml_r
    rect_pml = 1 .... square pml on square (-pml_x, pml_x)^d
    rect_pml = 2 .... rectangular pml on (pml_xmin, pml_xmax) x (pml_ymin, pml_ymax) x (pml_zmin, pml_zmax)
  */
  extern int rect_pml;
  
  extern bool apply_deriv_alpha;



  


  
  template <>
  SpecificIntegrationPoint<2,2,Complex> :: 
  SpecificIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,
			    LocalHeap & lh);

  
  template <>
  SpecificIntegrationPoint<2,2,AutoDiff<1,Complex> > :: 
  SpecificIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,
			    LocalHeap & lh);



  extern void SetPMLParameters();


  template <class DIFFOP, class DMATOP, class FEL = FiniteElement>
  class PML_BDBIntegrator : public T_BDBIntegrator<DIFFOP,DMATOP,FEL>
  {
  public:
  
    enum { DIM_SPACE   = DIFFOP::DIM_SPACE };
    enum { DIM_ELEMENT = DIFFOP::DIM_ELEMENT };
    enum { DIM_DMAT    = DIFFOP::DIM_DMAT };
    enum { DIM         = DIFFOP::DIM };
  
    ///
    PML_BDBIntegrator (const DMATOP & admat)
      : T_BDBIntegrator<DIFFOP,DMATOP,FEL> (admat)
    { 
      SetPMLParameters();

    }

    ///
    virtual ~PML_BDBIntegrator ()
    { ; }


    ///
    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
			   const ElementTransformation & eltrans, 
			   FlatMatrix<double> & elmat,
			   LocalHeap & locheap) const
    {
      throw Exception ("PML cannot generate real matrices");
    }

    ///
    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
			   const ElementTransformation & eltrans, 
			   FlatMatrix<Complex> & elmat,
			   LocalHeap & locheap) const
    {
      try
	{
	  const FEL & fel = static_cast<const FEL&> (bfel);
	  int ndof = fel.GetNDof();

	  elmat = 0;
	
	  FlatMatrixFixHeight<DIM_DMAT, Complex> bmat (ndof * DIM, locheap);
	  FlatMatrixFixHeight<DIM_DMAT, Complex> dbmat (ndof * DIM, locheap);

	  Mat<DIM_DMAT,DIM_DMAT, Complex> dmat;

	  const IntegrationRule & ir = GetIntegrationRule (fel);

	  for (int i = 0; i < ir.GetNIP(); i++)
	    {
              HeapReset hr (locheap);

	      SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,Complex> 
		sip(ir[i], eltrans, locheap);
	      SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,double> 
		sip_real(ir[i], eltrans, locheap);

	      DIFFOP::GenerateMatrix (fel, sip, bmat, locheap);
	      this->dmatop.GenerateMatrix (fel, sip_real, dmat, locheap);

	      Complex fac = sip.GetJacobiDet() * sip.IP().Weight();

              dmat *= fac;
	      dbmat = dmat * bmat;
	      // elmat += Trans (bmat) * dbmat;
              if (DMATOP::SYMMETRIC)
                FastMat<DIM_DMAT> (elmat.Height(), &dbmat(0,0), &bmat(0,0), &elmat(0,0));
              else
                elmat += Trans (bmat) * dbmat;
	    } 
	}

      catch (Exception & e)
	{
	  e.Append ("in CalcElementMatrix, type = ");
	  e.Append (typeid(*this).name());
	  e.Append ("\n");
	  throw;
	}
      catch (exception & e)
	{
	  Exception e2(e.what());
	  e2.Append ("\nin CalcElementMatrix, type = ");
	  e2.Append (typeid(*this).name());
	  e2.Append ("\n");
	  throw e2;
	}
    }



    virtual void 
    ApplyElementMatrix (const FiniteElement & bfel, 
			const ElementTransformation & eltrans, 
			const FlatVector<double> & elx, 
			FlatVector<double> & ely,
			void * precomputed,
			LocalHeap & locheap) const
    { ; }


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


      if (!apply_deriv_alpha)
	{
	  Vec<DIM_DMAT,Complex> hv1;
	  Vec<DIM_DMAT,Complex> hv2;
	
	  FlatVector<Complex> hely (ndof*DIM, locheap);
	  const IntegrationRule & ir = GetIntegrationRule (fel);
	
	
	  for (int i = 0; i < ir.GetNIP(); i++)
	    {
	      SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,Complex> 
		sip(ir[i], eltrans, locheap);

	      SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,double> 
		sip_real (ir[i], eltrans, locheap);
	    
	      DIFFOP::Apply (fel, sip, elx, hv1, locheap);
	      this -> dmatop.Apply (fel, sip_real, hv1, hv2, locheap);
	      DIFFOP::ApplyTrans (fel, sip, hv2, hely, locheap);
	    
	      Complex fac = sip.GetJacobiDet() * sip.IP().Weight();
	      ely += fac * hely;
	    }     
	}
      else
	{
	  Vec<DIM_DMAT, AutoDiff<1,Complex> > hv1;
	  Vec<DIM_DMAT, AutoDiff<1,Complex> > hv2;
	
	  FlatVector<AutoDiff<1,Complex> > hely (ndof*DIM, locheap);
	  FlatVector<AutoDiff<1,Complex> > helx (ndof*DIM, locheap);
	  const IntegrationRule & ir = GetIntegrationRule (fel);
	
	  for (int j = 0; j < helx.Size(); j++)
	    helx(j) = elx(j);
	
	  for (int i = 0; i < ir.GetNIP(); i++)
	    {
	      SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,AutoDiff<1,Complex> >
		sip(ir[i], eltrans, locheap);

	      SpecificIntegrationPoint<DIM_ELEMENT,DIM_SPACE,double> 
		sip_real (ir[i], eltrans, locheap);
	    
	      DIFFOP::Apply (fel, sip, helx, hv1, locheap);
	      this -> dmatop.Apply (fel, sip_real, hv1, hv2, locheap);
	      DIFFOP::ApplyTrans (fel, sip, hv2, hely, locheap);
	    
	      AutoDiff<1,Complex> fac = sip.GetJacobiDet() * sip.IP().Weight();

	      for (int j = 0; j < hely.Size(); j++)
		ely(j) += (fac * hely(j)).DValue(0);
	    }     
	}
    }


    ///
    virtual int GetDimension () const { return DIM; }
    ///
    virtual int Lumping () const { return 0; }
    ///
    virtual string Name () const { return "PML-BDB integrator"; }
  };

}

#endif
