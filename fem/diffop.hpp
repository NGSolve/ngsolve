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
     Provides application and transpose-application
  */
  template<class DOP>
  class DiffOp
  {
  public:
    // enum { DIM_ELEMENT = TDOP::DIM_ELEMENT };
    // enum { DIM_SPACE = TDOP::DIM_SPACE };
  
    /**
       Computes the B-matrix. 
       The height is DIM_DMAT, the width is fel.GetNDof(). 
       FEL is the FiniteElement type specified in the BDB-Integrator 
       sip is the mapped integration point containing the Jacobi-Matrix 
       MAT is the resulting matrix (usually a FixedHeightMatrix)
    */
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const SIP & sip,
				MAT & mat, LocalHeap & lh)
    {
      ;
    }

    /// Computes B-matrix times element vector
    template <typename FEL, typename SIP, class TVX, class TVY>
    static void Apply (const FEL & fel, const SIP & sip,
		       const TVX & x, TVY & y,
		       LocalHeap & lh)
    {
      typedef typename TVY::TSCAL TSCAL;

      HeapReset hr(lh);

      FlatMatrixFixHeight<DOP::DIM_DMAT, TSCAL> mat(DOP::DIM*fel.GetNDof(), lh);
      DOP::GenerateMatrix (fel, sip, mat, lh);
      y = mat * x;
    }

    /// Computes B-matrix times element vector in many points
    template <typename FEL, class MIR, class TVX, class TVY>
    static void ApplyIR (const FEL & fel, const MIR & mir,
			 const TVX & x, TVY & y,
			 LocalHeap & lh)
    {
      typedef typename TVY::TSCAL TSCAL;

      for (int i = 0; i < mir.Size(); i++)
	{
	  HeapReset hr(lh);
	  Apply (fel, mir[i], x, y.Row(i), lh);
	}
    }




    /// Computes Transpose (B-matrix) times point value
    template <typename FEL, typename SIP, class TVX, class TVY>
    static void ApplyTrans (const FEL & fel, const SIP & sip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh) 
    {
      typedef typename TVY::TSCAL TSCAL;

      HeapReset hr(lh);

      FlatMatrixFixHeight<DOP::DIM_DMAT, TSCAL> mat(DOP::DIM*fel.GetNDof(), lh);
      DOP::GenerateMatrix (fel, sip, mat, lh);
      y = Trans (mat) * x;
    }



    /// 
    template <typename SIP, class TVX>
    static void Transform (const SIP & sip, TVX & x)
    {
      ;
    }

    template <typename SIP, class TVX>
    static void TransformT (const SIP & sip, TVX & x)
    {
      ;
    }

    template <typename FEL, class TVD, class TVY, int D>
    static void ApplyGrid (const FEL & fel, 
			   const IntegrationRuleTP<D> & ir,
			   const TVY & hv, 
			   TVD & dvecs, 
			   LocalHeap & lh)
    {
      ;
    }

    template <typename FEL, class TVD, class TVY, int D>
    static void ApplyTransGrid (const FEL & fel, 
				const IntegrationRuleTP<D> & ir,
				const TVD & dvecs, 
				TVY & hv, LocalHeap & lh)
    {
      ;
    }

  };







  class DifferentialOperator
  {
  public:
    virtual int Dim() const = 0;    
    virtual bool Boundary() const = 0;

    virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseSpecificIntegrationPoint & sip,
		FlatMatrix<double> mat, 
		LocalHeap & lh) const = 0;

    virtual void
    Apply (const FiniteElement & fel,
	   const BaseSpecificIntegrationPoint & sip,
	   FlatVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const
    {
      FlatMatrix<> mat(Dim(), fel.GetNDof(), lh);
      CalcMatrix (fel, sip, mat, lh);
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
    ApplyTrans (const FiniteElement & fel,
		const BaseSpecificIntegrationPoint & sip,
		FlatVector<double> flux,
		FlatVector<double> x, 
		LocalHeap & lh) const 
    {
      FlatMatrix<> mat(Dim(), fel.GetNDof(), lh);
      CalcMatrix (fel, sip, mat, lh);
      flux = mat * x;
    }
  };






  template <typename DIFFOP>
  class T_DifferentialOperator : public DifferentialOperator
  {
  public:
    virtual int Dim() const { return DIFFOP::DIM_DMAT; }
    virtual bool Boundary() const { return DIFFOP::DIM_SPACE > DIFFOP::DIM_ELEMENT; }

    virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseSpecificIntegrationPoint & sip,
		FlatMatrix<double> mat, 
		LocalHeap & lh) const
    {
      DIFFOP::GenerateMatrix (fel, sip, mat, lh);
    }

    virtual void
    Apply (const FiniteElement & fel,
	   const BaseSpecificIntegrationPoint & sip,
	   FlatVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const
    {
      DIFFOP::Apply (fel, sip, x, flux, lh);
    }


    virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationRule & mir,
	   FlatVector<double> x, 
	   FlatMatrix<double> flux,
	   LocalHeap & lh) const
    {
      DIFFOP::ApplyIR (fel, mir, x, flux, lh);
    }




    virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseSpecificIntegrationPoint & sip,
		FlatVector<double> flux,
		FlatVector<double> x, 
		LocalHeap & lh) const 
    {
      DIFFOP::ApplyTrans (fel, sip, flux, x, lh);
    }    
  };




}


#endif
