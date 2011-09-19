#ifndef FILE_DIFFOP
#define FILE_DIFFOP

/*********************************************************************/
/* File:   diffop.hpp                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   24. Nov. 2009                                             */
/*********************************************************************/

namespace ngfem
{



  /**
     Differential Operator.
     Base-class for template-polymorphismus.
     Provides application and transpose-application.
     Operations can be applied for one integration point, or for the whole integration rule at once.
  */
  template<class DOP>
  class NGS_DLL_HEADER DiffOp
  {
  public:
    // enum { DIM_ELEMENT = TDOP::DIM_ELEMENT };
    // enum { DIM_SPACE = TDOP::DIM_SPACE };
  
    /**
       Computes the B-matrix. 
       The height is DIM_DMAT, the width is fel.GetNDof(). 
       FEL is the FiniteElement type specified in the BDB-Integrator 
       mip is the mapped integration point containing the Jacobi-Matrix 
       MAT is the resulting matrix (usually a FixedHeightMatrix)
    */
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      ;
    }

    /**
       Applies the B-matrix.
       Computes matrix-vector product with the B-matrix
    */
    template <typename FEL, typename MIP, class TVX, class TVY>
    static void Apply (const FEL & fel, const MIP & mip,
		       const TVX & x, TVY & y,
		       LocalHeap & lh)
    {
      typedef typename TVY::TSCAL TSCAL;

      HeapReset hr(lh);

      FlatMatrixFixHeight<DOP::DIM_DMAT, TSCAL> mat(DOP::DIM*fel.GetNDof(), lh);
      DOP::GenerateMatrix (fel, mip, mat, lh);
      y = mat * x;
    }

    /// Computes B-matrix times element vector in many points
    template <typename FEL, class MIR, class TVX, class TVY>
    static void ApplyIR (const FEL & fel, const MIR & mir,
			 const TVX & x, TVY & y,
			 LocalHeap & lh)
    {
      for (int i = 0; i < mir.Size(); i++)
	Apply (fel, mir[i], x, y.Row(i), lh);
    }


    /// Computes Transpose (B-matrix) times point value
    template <typename FEL, typename MIP, class TVX, class TVY>
    static void ApplyTrans (const FEL & fel, const MIP & mip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh) 
    {
      typedef typename TVY::TSCAL TSCAL;

      HeapReset hr(lh);

      FlatMatrixFixHeight<DOP::DIM_DMAT, TSCAL> mat(DOP::DIM*fel.GetNDof(), lh);
      DOP::GenerateMatrix (fel, mip, mat, lh);
      y = Trans (mat) * x;
    }


    /// Computes Transpose (B-matrix) times point value
    template <typename FEL, typename MIR, class TVX, class TVY>
    static void ApplyTransIR (const FEL & fel, const MIR & mir,
			      const TVX & x, TVY & y,
			      LocalHeap & lh) 
    {
      typedef typename TVY::TSCAL TSCAL;

      HeapReset hr(lh);
      FlatVector<TSCAL> hy(y.Size(), lh);

      y = 0.0;
      for (int i = 0; i < mir.Size(); i++)
	{
	  ApplyTrans (fel, mir[i], x.Row(i), hy, lh);
	  y += hy;
	}
    }





    /// old style ???
    template <typename MIP, class TVX>
    static void Transform (const MIP & mip, TVX & x)
    {
      ;
    }

    /// old style ???
    template <typename MIP, class TVX>
    static void TransformT (const MIP & mip, TVX & x)
    {
      ;
    }

    /// old style
    template <typename FEL, class TVD, class TVY, int D>
    static void ApplyGrid (const FEL & fel, 
			   const IntegrationRuleTP<D> & ir,
			   const TVY & hv, 
			   TVD & dvecs, 
			   LocalHeap & lh)
    {
      ;
    }

    /// old style
    template <typename FEL, class TVD, class TVY, int D>
    static void ApplyTransGrid (const FEL & fel, 
				const IntegrationRuleTP<D> & ir,
				const TVD & dvecs, 
				TVY & hv, LocalHeap & lh)
    {
      ;
    }

  };







  /**
     Differential Operator.
     Base-class for run-time polymorphismus.
     Provides application and transpose-application
  */
  class DifferentialOperator
  {
  public:
    /// destructor
    virtual ~DifferentialOperator ();
    /// dimension of range
    virtual int Dim() const = 0;   
    /// does it live on the boundary ?
    virtual bool Boundary() const = 0;
    /// calculates the matrix
    virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		FlatMatrix<double> mat, 
		LocalHeap & lh) const = 0;

    virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   FlatVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const
    {
      FlatMatrix<> mat(Dim(), fel.GetNDof(), lh);
      CalcMatrix (fel, mip, mat, lh);
      flux = mat * x;
    }

    virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   FlatVector<Complex> x, 
	   FlatVector<Complex> flux,
	   LocalHeap & lh) const
    {
      FlatMatrix<> mat(Dim(), fel.GetNDof(), lh);
      CalcMatrix (fel, mip, mat, lh);
      flux = mat * x;
    }

    virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationRule & mir,
	   FlatVector<double> x, 
	   FlatMatrix<double> flux,
	   LocalHeap & lh) const
    {
      for (int i = 0; i < mir.Size(); i++)
	Apply (fel, mir[i], x, flux.Row(i), lh);
    }

    virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationRule & mir,
	   FlatVector<Complex> x, 
	   FlatMatrix<Complex> flux,
	   LocalHeap & lh) const
    {
      for (int i = 0; i < mir.Size(); i++)
	Apply (fel, mir[i], x, flux.Row(i), lh);
    }




    virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		FlatVector<double> flux,
		FlatVector<double> x, 
		LocalHeap & lh) const 
    {
      FlatMatrix<> mat(Dim(), fel.GetNDof(), lh);
      CalcMatrix (fel, mip, mat, lh);
      flux = mat * x;
    }


    virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		FlatVector<Complex> flux,
		FlatVector<Complex> x, 
		LocalHeap & lh) const 
    {
      FlatMatrix<> mat(Dim(), fel.GetNDof(), lh);
      CalcMatrix (fel, mip, mat, lh);
      flux = mat * x;
    }


    virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		FlatMatrix<double> flux,
		FlatVector<double> x, 
		LocalHeap & lh) const 
    {
      FlatVector<double> hx(x.Size(), lh);
      x = 0.0;
      for (int i = 0; i < mir.Size(); i++)
	{
	  ApplyTrans (fel, mir[i], flux.Row(i), hx, lh);
	  x += hx;
	}
    }

    virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		FlatMatrix<Complex> flux,
		FlatVector<Complex> x, 
		LocalHeap & lh) const 
    {
      FlatVector<Complex> hx(x.Size(), lh);
      x = 0.0;
      for (int i = 0; i < mir.Size(); i++)
	{
	  ApplyTrans (fel, mir[i], flux.Row(i), hx, lh);
	  x += hx;
	}
    }
  };





  /**
     Connect compile-time polymorph DiffOp to run-time polymorph DifferentialOperator.
   */
  template <typename DIFFOP, typename FEL>
  class T_DifferentialOperator : public DifferentialOperator
  {
  protected:
    enum { DIM_SPACE   = DIFFOP::DIM_SPACE };
    enum { DIM_ELEMENT = DIFFOP::DIM_ELEMENT };
    enum { DIM_DMAT    = DIFFOP::DIM_DMAT };
    enum { DIM         = DIFFOP::DIM };

  public:
    virtual int Dim() const { return DIFFOP::DIM_DMAT; }
    virtual bool Boundary() const { return int(DIM_SPACE) > int(DIM_ELEMENT); }

    virtual void
    CalcMatrix (const FiniteElement & bfel,
		const BaseMappedIntegrationPoint & bmip,
		FlatMatrix<double> mat, 
		LocalHeap & lh) const
    {
      const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
	static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);

      const FEL & fel = static_cast<const FEL&> (bfel);

      DIFFOP::GenerateMatrix (fel, mip, mat, lh);
    }

    virtual void
    Apply (const FiniteElement & bfel,
	   const BaseMappedIntegrationPoint & bmip,
	   FlatVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const
    {
      const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
	static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);
      const FEL & fel = static_cast<const FEL&> (bfel);
      DIFFOP::Apply (fel, mip, x, flux, lh);
    }


    virtual void
    Apply (const FiniteElement & bfel,
	   const BaseMappedIntegrationRule & bmir,
	   FlatVector<double> x, 
	   FlatMatrix<double> flux,
	   LocalHeap & lh) const
    {
      const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
	static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);
      const FEL & fel = static_cast<const FEL&> (bfel);
      DIFFOP::ApplyIR (fel, mir, x, flux, lh);
    }




    virtual void
    ApplyTrans (const FiniteElement & bfel,
		const BaseMappedIntegrationPoint & bmip,
		FlatVector<double> flux,
		FlatVector<double> x, 
		LocalHeap & lh) const 
    {
      const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
	static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);
      const FEL & fel = static_cast<const FEL&> (bfel);
      DIFFOP::ApplyTrans (fel, mip, flux, x, lh);
    }    
  };




}


#endif
