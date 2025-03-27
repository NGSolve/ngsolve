#ifndef FILE_BDBINTEGRATOR
#define FILE_BDBINTEGRATOR

/*********************************************************************/
/* File:   bdbintegrator.hpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

// #include "fastmat.hpp"

#include "integrator.hpp"
#include "diffop.hpp"


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
	      const TVX & x, TVY && y,
	      LocalHeap & lh) const
  {
    Mat<DMO::DIM_DMAT, DMO::DIM_DMAT, double> mat;
    static_cast<const DMO*>(this) -> GenerateMatrix (fel, mip, mat, lh);
    y = mat * x;
  }


  template <typename FEL, typename MIP, class TVX, class TVY>
  void ApplyLinearized (const FEL & fel, const MIP & mip,
                        const TVX & lin, const TVX & x, TVY & y,
                        LocalHeap & lh) const
  {
    Mat<DMO::DIM_DMAT, DMO::DIM_DMAT, typename TVX::TSCAL> mat;
    static_cast<const DMO*>(this) -> GenerateLinearizedMatrix (fel, mip, lin, mat, lh);
    y = mat * x;
  }

  template <typename FEL, typename MIP, class TVX>
  void Apply1 (const FEL & fel, const MIP & mip,
	       TVX && x, LocalHeap & lh) const
  {
    Vec<DMO::DIM_DMAT, typename remove_reference<TVX>::type::TSCAL> y;
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
		 const TVX & x, TVY && y,
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


  /*
  template <int H, int DIST, typename T1, typename T2, typename T3>
  void FastMat (FlatMatrixFixHeight<H,T1,DIST> a,
                FlatMatrixFixHeight<H,T2,DIST> b,
                FlatMatrix<T3> c)
  {
    FastMat<H> (a.Width(), DIST, a.Data(), b.Data(), c.Data());
  }
  */
  

  template <class DMATOP> // , int DIM_ELEMENT, int DIM_SPACE>
  class T_BDBIntegrator_DMat : public BilinearFormIntegrator
  {
  protected:
    DMATOP dmatop;
    DifferentialOperator * diffop = NULL;
    enum { DIM_DMAT    = DMATOP::DIM_DMAT };
    
  public:
  
    T_BDBIntegrator_DMat  (const Array<shared_ptr<CoefficientFunction>> & coeffs);
    // : dmatop(coeffs) { ; }

    
    /*
    T_BDBIntegrator_DMat  (shared_ptr<CoefficientFunction> & c1)
      : dmatop(c1) { ; }
    */

    /*
    template <typename ... TMORE>
    T_BDBIntegrator_DMat  (shared_ptr<CoefficientFunction> c1, TMORE ... more_coefs);
    // : dmatop(c1, more_coefs ...) { ; }
    */

    /*
    T_BDBIntegrator_DMat  (const CoefficientFunction * coef)
      : dmatop(shared_ptr<CoefficientFunction> (const_cast<CoefficientFunction*>(coef), NOOP_Deleter))
    { ; }
    */
    ///
    T_BDBIntegrator_DMat (const DMATOP & admat);
    // : dmatop(admat) { ; }

    virtual ~T_BDBIntegrator_DMat () { delete diffop; }
    
    virtual xbool IsSymmetric () const override
    { return DMATOP::SYMMETRIC; }
    
    virtual int DimFlux () const override
    { return DMATOP::DIM_DMAT; }
    
    DMATOP & DMat () { return dmatop; }
    const DMATOP & DMat () const { return dmatop; }

    
    int GetIntegrationOrder (const FiniteElement & fel, 
                             const bool use_higher_integration_order = false) const
    {
      int order = 2 * fel.Order();
      
      ELEMENT_TYPE et = fel.ElementType();
      
      if (et == ET_TET || et == ET_TRIG || et == ET_SEGM)
        order -= 2 * diffop->DiffOrder();

      if (common_integration_order >= 0)
        order = common_integration_order;
      
      if (integration_order >= 0)
      order = integration_order;
      
      if(use_higher_integration_order && higher_integration_order > order)
        order = higher_integration_order;
      
      return order;
    }
    

    IntegrationRule GetIntegrationRule (const FiniteElement & fel, 
                                        const bool use_higher_integration_order = false) const
    {
      // return std::move(IntegrationRule (fel.ElementType(), GetIntegrationOrder(fel, use_higher_integration_order)));
      // return IntegrationRule (fel.ElementType(), GetIntegrationOrder(fel, use_higher_integration_order));
      return { fel.ElementType(), GetIntegrationOrder(fel, use_higher_integration_order) };
    }


    
    virtual void ApplyDMat (const FiniteElement & bfel,
                            const BaseMappedIntegrationPoint & bmip,
                            FlatVector<double> elx, 
                            FlatVector<double> eldx,
                            LocalHeap & lh) const override
    {
      dmatop.Apply(bfel, bmip,
                   // static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> &>(bmip),
                   elx, eldx ,lh);
    }

    virtual void ApplyDMat (const FiniteElement & bfel,
                            const BaseMappedIntegrationPoint & bmip,
                            FlatVector<Complex> elx, 
                            FlatVector<Complex> eldx,
                            LocalHeap & lh) const override
    {
      dmatop.Apply(bfel, bmip,
                   // static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> &>(bmip),
                   elx,eldx,lh);
    }
    
    virtual void ApplyDMat (const FiniteElement & bfel,
                            const BaseMappedIntegrationRule & bmir,
                            FlatMatrix<double> elx, 
                            FlatMatrix<double> eldx,
                            LocalHeap & lh) const override
    {
      // const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      // static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);
      
      for (int i = 0; i < bmir.Size(); i++)
        dmatop.Apply (bfel, bmir[i], elx.Row(i), eldx.Row(i), lh);
    }
    
    virtual void ApplyDMatInv (const FiniteElement & bfel,
                               const BaseMappedIntegrationRule & bmir,
                               FlatMatrix<double> elx, 
                               FlatMatrix<double> eldx,
                               LocalHeap & lh) const override
    {
      // const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      // static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);
      
      for (int i = 0; i < bmir.Size(); i++)
        dmatop.ApplyInv (bfel, bmir[i], elx.Row(i), eldx.Row(i), lh);
    }

    
    virtual void ApplyDMat (const FiniteElement & bfel,
                            const BaseMappedIntegrationRule & bmir,
                            FlatMatrix<Complex> elx, 
                            FlatMatrix<Complex> eldx,
                            LocalHeap & lh) const override
    {
      // const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      // static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);
      
      for (int i = 0; i < bmir.Size(); i++)
        dmatop.Apply (bfel, bmir[i], elx.Row(i), eldx.Row(i), lh);
    }
    




    virtual void
    CalcFlux (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & bmip,
              BareSliceVector<double> elx, 
              FlatVector<double> flux,
              bool applyd,
              LocalHeap & lh) const override
    {
      diffop -> Apply (fel, bmip, elx, flux, lh);

      FlatVec<DMATOP::DIM_DMAT,double> hflux(&flux(0));
      if (applyd)
        dmatop.Apply1 (fel, bmip, hflux, lh);
    }

    virtual void
    CalcFlux (const FiniteElement & fel,
              const BaseMappedIntegrationRule & bmir,
              BareSliceVector<double> elx, 
              BareSliceMatrix<double> flux,
              bool applyd,
              LocalHeap & lh) const override
    {
      diffop->Apply (fel, bmir, elx, flux, lh);
      
      FlatMatrixFixWidth<DMATOP::DIM_DMAT,double> hflux(bmir.Size(), &flux(0,0));
      if (applyd)
        dmatop.ApplyIR (fel, bmir, hflux, lh);
    }

    virtual void
    CalcFlux (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & bmip,
              BareSliceVector<Complex> elx, 
              FlatVector<Complex> flux,
              bool applyd,
              LocalHeap & lh) const override
    {
      diffop->Apply (fel, bmip, elx, flux, lh);

      FlatVec<DMATOP::DIM_DMAT,Complex> hflux(&flux(0));
      if (applyd)
        dmatop.Apply1 (fel, bmip, hflux, lh);
    }

    virtual void
    CalcFlux (const FiniteElement & fel,
              const BaseMappedIntegrationRule & bmir,
              BareSliceVector<Complex> elx, 
              BareSliceMatrix<Complex> flux,
              bool applyd,
              LocalHeap & lh) const override
    {
      diffop->Apply (fel, bmir, elx, flux, lh);
      
      FlatMatrixFixWidth<DMATOP::DIM_DMAT,Complex> hflux(bmir.Size(), &flux(0,0));
      if (applyd)
        dmatop.ApplyIR (fel, bmir, hflux, lh);
    }




  virtual void
  CalcFluxMulti (const FiniteElement & fel,
		 const BaseMappedIntegrationPoint & mip,		
		 int m,
                 FlatVector<double> elx, 
		 FlatVector<double> flux,
		 bool applyd,
		 LocalHeap & lh) const override
  {
    int ndof = fel.GetNDof();
    int dimension = this->GetDimension();
    FlatMatrixFixHeight<DIM_DMAT> bmat (ndof * dimension, lh);

    diffop->CalcMatrix (fel, mip, bmat, lh);

    if (applyd)
      {
	Vec<DIM_DMAT> hv1;
	Mat<DIM_DMAT,DIM_DMAT> dmat;
	dmatop.GenerateMatrix (fel, mip, dmat, lh);

	for (int i = 0; i < m; i++)
	  {
	    SliceVector<double> slice_x (ndof*dimension, m, &const_cast<double&> (elx(i)));
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
            SliceVector<double> slice_x (ndof*dimension, m, &elx(i));
	    SliceVector<double> slice_flux (DIM_DMAT, m, &flux(i));
	    slice_flux = bmat * slice_x;
	  }
      }
  }



    virtual void
    ApplyBTrans (const FiniteElement & fel,
                 const BaseMappedIntegrationPoint & bmip,
                 FlatVector<double> elx, 
                 FlatVector<double> ely,
                 LocalHeap & lh) const override
    {
      diffop->ApplyTrans (fel, bmip, elx, ely, lh);
    }
  
    
    virtual void
    ApplyBTrans (const FiniteElement & fel,
                 const BaseMappedIntegrationPoint & bmip,
                 FlatVector<Complex> elx, 
                 FlatVector<Complex> ely,
                 LocalHeap & lh) const override
    {
      diffop->ApplyTrans (fel, bmip, elx, ely, lh);
    }


    virtual void
    ApplyBTrans (const FiniteElement & fel,
                 const BaseMappedIntegrationRule & bmir,
                 FlatMatrix<double> elx, 
                 FlatVector<double> ely,
                 LocalHeap & lh) const override
    {
      // const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      // static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);
      diffop->ApplyTrans (fel, bmir, elx, ely, lh);
    }

    
    virtual void 
    ApplyElementMatrix (const FiniteElement & bfel, 
                        const ElementTransformation & eltrans, 
                        const FlatVector<double> elx, 
                        FlatVector<double> ely,
                        void * precomputed,
                        LocalHeap & lh) const override
    {
      T_ApplyElementMatrix<double> (bfel, eltrans, elx, ely, precomputed, lh);
    }
    
    virtual void 
    ApplyElementMatrix (const FiniteElement & bfel, 
                        const ElementTransformation & eltrans, 
                        FlatVector<Complex> elx, 
                        FlatVector<Complex> ely,
                        void * precomputed,
                        LocalHeap & lh) const override
    {
      T_ApplyElementMatrix<Complex> (bfel, eltrans, elx, ely, precomputed, lh);
    }
    
    
    template <typename TSCAL>
    void T_ApplyElementMatrix (const FiniteElement & fel, 
                               const ElementTransformation & eltrans, 
                               FlatVector<TSCAL> elx, 
                               FlatVector<TSCAL> ely,
                               void * precomputed,
                               LocalHeap & lh) const
    {
      const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());
      BaseMappedIntegrationRule & mir = eltrans(ir, lh); 
      
      FlatMatrixFixWidth<DMATOP::DIM_DMAT, TSCAL> hv1(ir.GetNIP(), lh);
      diffop->Apply (fel, mir, elx, hv1, lh);
      dmatop.ApplyIR (fel, mir, hv1, lh);
      for (int i = 0; i < mir.Size(); i++)
        hv1.Row(i) *= mir[i].GetWeight();
      diffop->ApplyTrans (fel, mir, hv1, ely, lh);    
    }
    

    ///
    virtual void 
    ApplyMixedElementMatrix (const FiniteElement & bfel1, 
                             const FiniteElement & bfel2, 
                             const ElementTransformation & eltrans, 
                             FlatVector<double> elx, 
                             FlatVector<double> ely,
                             LocalHeap & lh) const
    {
      HeapReset hr1 (lh);

      ely = 0;
      
      Vec<DIM_DMAT,double> hv1;
      Vec<DIM_DMAT,double> hv2;
      
      FlatVector<double> hely (ely.Size(), lh);
      
      const IntegrationRule & ir = GetIntegrationRule (bfel2,eltrans.HigherIntegrationOrderSet());
      
      for (int i = 0; i < ir.GetNIP(); i++)
        {
          HeapReset hr (lh);
          BaseMappedIntegrationPoint & mip = eltrans(ir[i], lh);
          // MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> mip (ir[i], eltrans);
          
          diffop->Apply (bfel1, mip, elx, hv1, lh);
          dmatop.Apply (bfel1, mip, hv1, hv2, lh);
          diffop->ApplyTrans (bfel2, mip, hv2, hely, lh);
          
          ely += mip.GetWeight() * hely;
        }     
    }




    
    virtual void
    CalcElementMatrixDiag (const FiniteElement & fel,
                           const ElementTransformation & eltrans, 
                           FlatVector<double> diag,
                           LocalHeap & lh) const override
    {
      try
        {
          // diag.AssignMemory (ndof*DIM, lh);
          diag = 0.0;
          
          FlatMatrixFixHeight<DIM_DMAT, double> bmat (diag.Size(), lh);
          Mat<DIM_DMAT,DIM_DMAT> dmat;
          
          const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());

          for (int i = 0; i < ir.GetNIP(); i++)
            {
              HeapReset hr(lh);
              BaseMappedIntegrationPoint & mip = eltrans(ir[i], lh);
              
              diffop->CalcMatrix (fel, mip, bmat, lh);
              dmatop.GenerateMatrix (fel, mip, dmat, lh);
              
              double fac =  mip.GetWeight();
              for (int j = 0; j < diag.Size(); j++)
                {
                  Vec<DIM_DMAT> hv = dmat * bmat.Col(j);
                  diag(j) += fac * InnerProduct (bmat.Col(j), hv);
                }
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

    

  };



  /*
  template <class DMATOP, int DIM_ELEMENT, int DIM_SPACE>
  T_BDBIntegrator_DMat<DMATOP,DIM_ELEMENT, DIM_SPACE>:: 
  T_BDBIntegrator_DMat  (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : dmatop(coeffs) { ; }
    
  template <class DMATOP, int DIM_ELEMENT, int DIM_SPACE> template <typename ... TMORE>
  T_BDBIntegrator_DMat<DMATOP,DIM_ELEMENT, DIM_SPACE>:: 
  T_BDBIntegrator_DMat  (shared_ptr<CoefficientFunction> c1, TMORE ... more_coefs)
    : dmatop(c1, more_coefs ...) { ; }
    
  
  template <class DMATOP, int DIM_ELEMENT, int DIM_SPACE>
  T_BDBIntegrator_DMat<DMATOP,DIM_ELEMENT, DIM_SPACE>:: 
  T_BDBIntegrator_DMat (const DMATOP & admat) 
    : dmatop(admat) { ; }
  */

  template <class DMATOP>
  T_BDBIntegrator_DMat<DMATOP>:: 
  T_BDBIntegrator_DMat  (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : dmatop(coeffs) { ; }
    
  /*
  template <class DMATOP> template <typename ... TMORE>
  T_BDBIntegrator_DMat<DMATOP>:: 
  T_BDBIntegrator_DMat  (shared_ptr<CoefficientFunction> c1, TMORE ... more_coefs)
    : dmatop(c1, more_coefs ...) { ; }
  */
  
  template <class DMATOP>
  T_BDBIntegrator_DMat<DMATOP>:: 
  T_BDBIntegrator_DMat (const DMATOP & admat) 
    : dmatop(admat) { ; }




 
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
class T_BDBIntegrator : public T_BDBIntegrator_DMat<DMATOP /* , DIFFOP::DIM_ELEMENT, DIFFOP::DIM_SPACE */ > 
{
protected:
  typedef T_BDBIntegrator_DMat<DMATOP /* , DIFFOP::DIM_ELEMENT, DIFFOP::DIM_SPACE */ > BASE;

  using BASE::diffop;

public:

  enum { DIM_SPACE   = DIFFOP::DIM_SPACE };
  enum { DIM_ELEMENT = DIFFOP::DIM_ELEMENT };
  enum { DIM_DMAT    = DIFFOP::DIM_DMAT };
  enum { DIM         = DIFFOP::DIM };

  using BASE::Name;
  using BASE::integration_order;
  using BASE::higher_integration_order;
  using BASE::common_integration_order;
  using BASE::GetIntegrationOrder;
  using BASE::GetIntegrationRule;

  using BASE::dmatop;
  
  /// inherited constructors
  // using BASE::T_BDBIntegrator_DMat;


  T_BDBIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : BASE(coeffs)
  {
    diffop = new T_DifferentialOperator<DIFFOP>; 
  }

  T_BDBIntegrator  (const shared_ptr<CoefficientFunction> & c1)
    : BASE(c1) 
  { 
    diffop = new T_DifferentialOperator<DIFFOP>; 
  }

  /*
  template <typename ... TMORE>
  T_BDBIntegrator  (shared_ptr<CoefficientFunction> c1, TMORE ... more_coefs)
    : BASE(c1, more_coefs ...) 
  { 
    diffop = new T_DifferentialOperator<DIFFOP>; 
  }
  */

  /*
  T_BDBIntegrator  (const CoefficientFunction * coef)
    : BASE (coef) 
  { 
    diffop = new T_DifferentialOperator<DIFFOP>; 
  } 
  */

  T_BDBIntegrator (const DMATOP & admat)
    : BASE(admat) 
  { 
    diffop = new T_DifferentialOperator<DIFFOP>; 
  }


  
  /*
  template <typename ... ARGS>
  T_BDBIntegrator (ARGS ... args)
    : BASE(args...)
  {
    diffop = make_shared<T_DifferentialOperator<DIFFOP>>(); 
  }
  */

  ///
  virtual ~T_BDBIntegrator () { ; }

  ///
  virtual int GetDimension () const { return DIM; }
  
  ///
  //virtual bool BoundaryForm () const
  //{ return int(DIM_SPACE) > int(DIM_ELEMENT); }
  virtual VorB VB() const
  { return VorB(int(DIM_SPACE)-int(DIM_ELEMENT)); }
    virtual int DimElement () const
  { return DIM_ELEMENT; }
  
  virtual int DimSpace () const
  { return DIM_SPACE; }

  virtual void CheckElement (const FiniteElement & el) const
  {
    if (!dynamic_cast<const FEL*> (&el) )
      {
        string err("Element does not match integrator\n");
        err += "element type is ";
        err += typeid(el).name();
        err += " expected type is ";
        err += typeid(FEL).name();
        err += " integrator is ";
        err += Name();
        throw Exception (err);
        /*
        throw Exception (string ("Element does not match integrator\n") +
                         string ("element type is ") + typeid(el).name() +
                         string (" expected type is ") + typeid(FEL).name() +
                         string (" integrator is ") + Name());
        */
      }
  }

  virtual void
  CalcElementMatrix (const FiniteElement & bfel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    T_CalcElementMatrix<double> (bfel, eltrans, elmat, lh);
  }

  virtual void
  CalcElementMatrix (const FiniteElement & bfel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<Complex> elmat,
		     LocalHeap & lh) const
  {
    T_CalcElementMatrix<Complex> (bfel, eltrans, elmat, lh);
  }


#ifdef TEXT_BOOK_VERSION

  template <typename TSCAL>
  void T_CalcElementMatrix (const FiniteElement & fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<TSCAL> elmat,
			    LocalHeap & lh) const
  {
    try
      {
	// const FEL & fel = static_cast<const FEL&> (bfel);
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
  // #define BLOCK_VERSION
#endif

  // #ifdef __MIC__
  // #define BLOCK_VERSION
  // #endif



#ifdef BLOCK_VERSION

  template <typename TSCAL>
  void T_CalcElementMatrix (const FiniteElement & fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<TSCAL> elmat,
			    LocalHeap & lh) const
  {
    // static Timer timer(string ("Elementmatrix, ") + Name()); RegionTimer reg (timer);
    try
      {
	int ndof = fel.GetNDof();
	elmat = 0;

        enum { BLOCK = 4 * (6 / DIM_DMAT + 1) };
	
        HeapReset hr1(lh);

#ifdef __MIC__
	enum { ROUNDUP = (DIM_DMAT*BLOCK+7) & (-8) };
#else
	enum { ROUNDUP = (DIM_DMAT*BLOCK+3) & (-4) };
#endif

        FlatMatrixFixHeight<DIM_DMAT*BLOCK, double, ROUNDUP> bbmat (ndof * DIM, lh);
        FlatMatrixFixHeight<DIM_DMAT*BLOCK, TSCAL, ROUNDUP> bdbmat (ndof * DIM, lh);

        typedef decltype (DMATOP::GetMatrixType(TSCAL(0))) TDMAT;
	
        IntegrationRule ir(fel.ElementType(), 
                           GetIntegrationOrder(fel,eltrans.HigherIntegrationOrderSet()));

	MappedIntegrationRule<DIM_ELEMENT, DIM_SPACE> mir(ir, eltrans, lh);


        FlatArray<TDMAT> dmats(ir.GetNIP(), lh);
	dmatop.GenerateMatrixIR (fel, mir, dmats, lh);

        int i = 0;
        for (int i1 = 0; i1 < ir.GetNIP() / BLOCK; i1++)
          {
            DIFFOP::GenerateMatrixIR (fel, mir.Range(i,i+BLOCK), bbmat, lh);

            for (int i2 = 0; i2 < BLOCK; i2++)
              {
		IntRange rows (i2*DIM_DMAT, (i2+1)*DIM_DMAT);
		TDMAT dmat = mir[i+i2].GetWeight() * dmats[i+i2];
		bdbmat.Rows(rows) = dmat * bbmat.Rows(rows);
              }

            i += BLOCK;
	    
            if (DMATOP::SYMMETRIC)
              FastMat (bdbmat, bbmat, elmat);
            else
              elmat += Trans (bbmat.Rows(0,DIM_DMAT*BLOCK)) * bdbmat.Rows(0,DIM_DMAT*BLOCK); 
          } 

        int rest = ir.GetNIP()-i;
        if (rest > 0)
          {
            DIFFOP::GenerateMatrixIR (fel, mir.Range(i,ir.GetNIP()), bbmat, lh);
            for (int i2 = 0; i2 < rest; i2++)
              {
		IntRange rows (i2*DIM_DMAT, (i2+1)*DIM_DMAT);
		TDMAT dmat = mir[i+i2].GetWeight() * dmats[i+i2];
		bdbmat.Rows(rows) = dmat * bbmat.Rows(rows);
              }

            if (DMATOP::SYMMETRIC)
              {
                /*
                int j = 0;
                for ( ; j <= rest-4; j += 4)
                  FastMat (FlatMatrixFixHeight<4*DIM_DMAT,TSCAL, ROUNDUP> (ndof*DIM, &bdbmat(j*DIM_DMAT)), 
                           FlatMatrixFixHeight<4*DIM_DMAT,double, ROUNDUP> (ndof*DIM, &bbmat(j*DIM_DMAT)), 
                           elmat);
                if (j <= rest-2)
                  {
                    FastMat (FlatMatrixFixHeight<2*DIM_DMAT,TSCAL, ROUNDUP> (ndof*DIM, &bdbmat(j*DIM_DMAT)), 
                             FlatMatrixFixHeight<2*DIM_DMAT,double, ROUNDUP> (ndof*DIM, &bbmat(j*DIM_DMAT)), 
                             elmat);
                    j += 2;
                  }
                if (j <= rest-1)
                  FastMat (FlatMatrixFixHeight<DIM_DMAT,TSCAL, ROUNDUP> (ndof*DIM, &bdbmat(j*DIM_DMAT)), 
                           FlatMatrixFixHeight<DIM_DMAT,double, ROUNDUP> (ndof*DIM, &bbmat(j*DIM_DMAT)), 
                           elmat);
                */


                int j = 0;
                int rd = rest*DIM_DMAT;

                for ( ; j <= rd-8; j += 8)
                  /*
                  FastMat (FlatMatrixFixHeight<8,TSCAL, ROUNDUP> (ndof*DIM, &bdbmat(j)), 
                           FlatMatrixFixHeight<8,double, ROUNDUP> (ndof*DIM, &bbmat(j)), 
                           elmat);
                  */
                  FastMat(bdbmat.template Rows<8>(j), bbmat.template Rows<8>(j), elmat);

                switch (rd - j)
                  {
                  case 1:
                    /*
                    FastMat (FlatMatrixFixHeight<1,TSCAL, ROUNDUP> (ndof*DIM, &bdbmat(j)), 
                             FlatMatrixFixHeight<1,double, ROUNDUP> (ndof*DIM, &bbmat(j)), 
                             elmat);
                    */
                    FastMat(bdbmat.template Rows<1>(j), bbmat.template Rows<1>(j), elmat);
                    break;
                  case 2:
                    /*
                    FastMat (FlatMatrixFixHeight<2,TSCAL, ROUNDUP> (ndof*DIM, &bdbmat(j)), 
                             FlatMatrixFixHeight<2,double, ROUNDUP> (ndof*DIM, &bbmat(j)), 
                             elmat);
                    */
                    FastMat(bdbmat.template Rows<2>(j), bbmat.template Rows<2>(j), elmat);
                    break;
                  case 3:
                    /*
                    FastMat (FlatMatrixFixHeight<3,TSCAL, ROUNDUP> (ndof*DIM, &bdbmat(j)), 
                             FlatMatrixFixHeight<3,double, ROUNDUP> (ndof*DIM, &bbmat(j)), 
                             elmat);
                    */
                    FastMat (bdbmat.template Rows<3>(j), bbmat.template Rows<3>(j), elmat);
                    break;
                  case 4:
                    /*
                    FastMat (FlatMatrixFixHeight<4,TSCAL, ROUNDUP> (ndof*DIM, &bdbmat(j)), 
                             FlatMatrixFixHeight<4,double, ROUNDUP> (ndof*DIM, &bbmat(j)), 
                             elmat);
                    */
                    FastMat(bdbmat.template Rows<4>(j), bbmat.template Rows<4>(j), elmat);
                    break;
                  case 5:
                    /*
                    FastMat (FlatMatrixFixHeight<5,TSCAL, ROUNDUP> (ndof*DIM, &bdbmat(j)), 
                             FlatMatrixFixHeight<5,double, ROUNDUP> (ndof*DIM, &bbmat(j)), 
                             elmat);
                    */
                    FastMat(bdbmat.template Rows<5>(j), bbmat.template Rows<5>(j), elmat);
                    break;
                  case 6:
                    /*
                    FastMat (FlatMatrixFixHeight<6,TSCAL, ROUNDUP> (ndof*DIM, &bdbmat(j)), 
                             FlatMatrixFixHeight<6,double, ROUNDUP> (ndof*DIM, &bbmat(j)), 
                             elmat);
                    */
                    FastMat(bdbmat.template Rows<6>(j), bbmat.template Rows<6>(j), elmat);
                    break;
                  case 7:
                    /*
                    FastMat (FlatMatrixFixHeight<7,TSCAL, ROUNDUP> (ndof*DIM, &bdbmat(j)), 
                             FlatMatrixFixHeight<7,double, ROUNDUP> (ndof*DIM, &bbmat(j)), 
                             elmat);
                    */
                    FastMat(bdbmat.template Rows<7>(j), bbmat.template Rows<7>(j), elmat);
                    break;
                  default:
                    ;
                  }
              }
            else
              elmat += Trans (bbmat.Rows(0,rest*DIM_DMAT)) * bdbmat.Rows(0,rest*DIM_DMAT);
          }

        
        if (DMATOP::SYMMETRIC)
          {
            for (int i = 0; i < elmat.Height(); i++)
              for (int j = 0; j < i; j++)
                elmat(j,i) = elmat(i,j);
          }

        ir.NothingToDelete();
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
			    FlatMatrix<TSCAL> elmat,
			    LocalHeap & lh) const
  {
    // static Timer timer (string ("Elementmatrix, ") + Name(), NoTracing);
    // static Timer timer2 (string ("Elementmatrix, ") + Name() + ", Lapack", NoTracing, NoTiming);
    // RegionTimer reg (timer);

    // try
      {
	// const FEL & fel = static_cast<const FEL&> (bfel);
	int ndof = fel.GetNDof();

	HeapReset hr(lh);

	const IntegrationRule & ir = GetIntegrationRule (fel,eltrans.HigherIntegrationOrderSet());
	MappedIntegrationRule<DIM_ELEMENT, DIM_SPACE> mir(ir, eltrans, lh);

	FlatMatrixFixHeight<DIM_DMAT> bmat (ndof * DIM, lh);
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

	// RegionTimer reg2 (timer2);

	if (ndof < 20)
	  {
	    if (DMATOP::SYMMETRIC)
	      elmat = Symmetric ( bdbmat * Trans (bbmat));
	    else
	      elmat = bbmat * Trans (bdbmat); 
	  }
	else
	  elmat = bbmat * Trans(bdbmat) | Lapack;

	// timer.AddFlops (long(elmat.Height())*long(elmat.Width())*bbmat.Width());
      } 
    
      /*
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
      */
  }

#endif



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
  using T_BDBIntegrator<DIFFOP,DMATOP,FEL>::diffop;
  
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
				   FlatVector<double> elveclin,
				   FlatMatrix<double> elmat,
				   LocalHeap & lh) const 
  {
    // static Timer maintimer (string ("NonlinearBDB, CalcLinearized, ") + this->Name(), NoTracing);
    // static Timer bdbtimer ("NonlinearBDB, bdb product", NoTracing, NoTiming);
    // RegionTimer reg(maintimer);

    try
      {
	int ndof = fel.GetNDof();

	elmat = 0;
	
	FlatMatrixFixHeight<DIM_DMAT, double> bmat (ndof * DIM, lh);
	FlatMatrixFixHeight<DIM_DMAT, double> dbmat (ndof * DIM, lh);
	Vec<DIM_DMAT,double> hvlin;

	Mat<DIM_DMAT,DIM_DMAT> dmat;

	const IntegrationRule & ir = GetIntegrationRule (fel);

	for (int i = 0; i < ir.GetNIP(); i++)
	  {
            HeapReset hr(lh);

	    MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> mip(ir[i], eltrans);

	    DIFFOP::Apply (fel, mip, elveclin, hvlin, lh);
	    DIFFOP::GenerateMatrix (fel, mip, bmat, lh);

	    this->dmatop . GenerateLinearizedMatrix (fel, mip, hvlin, dmat, lh);

	    double fac = mip.GetWeight();

	    {
	      // NgProfiler::RegionTimer reg(bdbtimer);
	      dbmat = fac * (dmat * bmat);
	      elmat += Trans (bmat) * dbmat; 
	    }
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
			       FlatVector<Complex> elveclin,
			       FlatMatrix<Complex> elmat,
			       LocalHeap & lh) const 
  {
    // static Timer maintimer (string ("NonlinearBDB, CalcLinearized<Complex>, ") + this->Name(), NoTracing);
    // static Timer bdbtimer ("NonlinearBDB, bdb product", NoTracing, NoTiming);
    // RegionTimer reg(maintimer);

    try
      {
	const FEL & fel = static_cast<const FEL&> (bfel);
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

	    diffop->Apply (fel, mip, elveclin, hvlin, lh);
	    DIFFOP::GenerateMatrix (fel, mip, bmat, lh);
	    this->dmatop . GenerateLinearizedMatrix (fel, mip, hvlin, dmat, lh);

	    dmat *= mip.GetWeight();
	    bbmat.Cols(i*DIM_DMAT, (i+1)*DIM_DMAT) = Trans (bmat);
	    bdbmat.Cols(i*DIM_DMAT, (i+1)*DIM_DMAT) = Trans (dmat * bmat);
	  } 


	// RegionTimer reg(bdbtimer);

	if (ndof < 20)
	  {
	    if (DMATOP::SYMMETRIC)
	      elmat = Symmetric ( bdbmat * Trans (bbmat));
	    else
	      elmat = bbmat * Trans (bdbmat); 
	  }
	else
	  // LapackMultABt (bbmat, bdbmat, elmat);
	  // LapackMult (bbmat, Trans(bdbmat), elmat);
          elmat = bbmat * Trans (bdbmat) | Lapack;
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
                                FlatVector<double> ellin, 
                                FlatVector<double> elx, 
				FlatVector<double> ely,
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

	diffop->Apply (fel, mip, ellin, hvlin, lh);
	diffop->Apply (fel, mip, elx, hvx, lh);
	this->dmatop.ApplyLinearized (fel, mip, hvlin, hvx, hvy, lh);
	diffop->ApplyTrans (fel, mip, hvy, hely, lh);

	double fac = fabs (mip.GetJacobiDet()) * ip.Weight();
	ely += fac * hely;
      }     
  }



 ///
  virtual void 
  ApplyLinearizedElementMatrix (const FiniteElement & bfel, 
				const ElementTransformation & eltrans, 
                                FlatVector<Complex> ellin, 
                                FlatVector<Complex> elx, 
				FlatVector<Complex> ely,
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
                         FlatVector<double> elx, 
			 LocalHeap & lh) const
  {
    const FEL & fel = dynamic_cast<const FEL&> (bfel);
    
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
  // DifferentialOperator * diffop = new T_DifferentialOperator<DIFFOP>;
  unique_ptr<DifferentialOperator> diffop = make_unique<T_DifferentialOperator<DIFFOP>>();  
public:
  enum { DIM_SPACE = DIFFOP::DIM_SPACE };
  enum { DIM_ELEMENT = DIFFOP::DIM_ELEMENT };
  enum { DIM_DMAT = DIFFOP::DIM_DMAT };
  enum { DIM = DIFFOP::DIM };
  // typedef typename DVecOp::TSCAL TSCAL;

  ///
  T_BIntegrator  (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : dvecop(coeffs)
  { ; }

  T_BIntegrator  (shared_ptr<CoefficientFunction> coef)
    : dvecop(coef)
  { ; }

  T_BIntegrator  (const CoefficientFunction * coef)
    : dvecop(shared_ptr<CoefficientFunction> (const_cast<CoefficientFunction*>(coef), NOOP_Deleter))
  { ; }

  T_BIntegrator (const DVecOp & advec)
    : dvecop(advec)
  { ; }

  ///
  virtual ~T_BIntegrator ()
  { ; }
  ///
  void CheckElement (const FiniteElement & el) const override
  {
    if (!dynamic_cast<const FEL*> (&el) )
      {
        string err("Element does not match integrator\n");
        err += "element type is ";
        err += typeid(el).name();
        err += " expected type is ";
        err += typeid(FEL).name();
        err += " integrator is ";
        err += Name();
        throw Exception (err);
        /*
        throw Exception (string ("Element does not match integrator\n") +
                       string ("element type is ") + typeid(el).name() +
                       string (" expected type is ") + typeid(FEL).name() +
                       string ("integrator is ") + Name());
        */
      }

  }
  ///
  //virtual bool BoundaryForm () const
  //{ return int(DIM_SPACE) > int(DIM_ELEMENT); }
  VorB VB() const override
  { return VorB(int(DIM_SPACE)-int(DIM_ELEMENT)); }
  int DimElement () const override
  { return DIM_ELEMENT; }

  int DimSpace () const override
  { return DIM_SPACE; }



  void
  CalcElementVector (const FiniteElement & bfel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> elvec,
		     LocalHeap & lh) const override
  {
    T_CalcElementVector<double> (bfel, eltrans, elvec, lh);
  }

  void
  CalcElementVector (const FiniteElement & bfel,
		     const ElementTransformation & eltrans, 
		     FlatVector<Complex> elvec,
		     LocalHeap & lh) const override
  {
    T_CalcElementVector<Complex> (bfel, eltrans, elvec, lh);
  }


  ///
  template <typename TSCAL>
  void T_CalcElementVector (const FiniteElement & fel,
			    const ElementTransformation & eltrans, 
			    FlatVector<TSCAL> elvec,
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

        // DIFFOP::ApplyTransIR (fel, mir, dvecs, elvec, lh);
        diffop->ApplyTrans (fel, mir, dvecs, elvec, lh);
      }
    catch (Exception & e)
      {
        e.Append (string ("in CalcElementVector <")+typeid(TSCAL).name()+
                  ">, type = " + typeid(*this).name() + "\n");
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
				    const bool curveint = false) const override
  {
    T_CalcElementVectorIndependent (gfel, s_mip, g_mip, elvec, lh, curveint);
  }

  virtual void
  CalcElementVectorIndependent (const FiniteElement & gfel, 
				    const BaseMappedIntegrationPoint & s_mip,
				    const BaseMappedIntegrationPoint & g_mip,
				    FlatVector<Complex> & elvec,
				    LocalHeap & lh,
				    const bool curveint = false) const override
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

    diffop->ApplyTrans (fel, d_g_mip, dvec, elvec, lh);
  
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
  int GetDimension () const override { return DIM; }

  ///
  string Name () const override { return "B integrator"; }
};



}


#endif
