#ifndef FILE_SCALARFE
#define FILE_SCALARFE

/*********************************************************************/
/* File:   scalarfe.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

#include "finiteelement.hpp"
#include "fe_interfaces.hpp"

namespace ngfem
{
  
  class BaseScalarFiniteElement : public FiniteElement 
  {
  public:
    using FiniteElement::FiniteElement;

    /// the name
    NGS_DLL_HEADER
    virtual string ClassName() const override;

    /// compute shape
    HD NGS_DLL_HEADER 
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const = 0;

    HD NGS_DLL_HEADER 
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<Complex> shape) const;
    
    /// compute dshape, matrix: ndof x spacedim
    HD NGS_DLL_HEADER 
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     BareSliceMatrix<> dshape) const = 0;
    

    /**
       returns shape functions in point ip.
    */
    INLINE FlatVector<> GetShape (const IntegrationPoint & ip, 
				  LocalHeap & lh) const
    {
      FlatVector<> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

    /// compute shape, row is shape nr, col is ip nr
    HD NGS_DLL_HEADER 
    virtual void CalcShape (const IntegrationRule & ir, 
                            BareSliceMatrix<> shape) const;
  
    /// compute shape, row is shape nr, col is ip nr
    HD NGS_DLL_HEADER 
    virtual void CalcShape (const SIMD_IntegrationRule & ir, 
                            BareSliceMatrix<SIMD<double>> shape) const;

    HD NGS_DLL_HEADER 
    virtual void CalcMappedDShape (const BaseMappedIntegrationPoint & mip, 
                                   BareSliceMatrix<> dshape) const = 0;


    HD NGS_DLL_HEADER 
    virtual void CalcMappedDShape (const BaseMappedIntegrationRule & mir, 
                                   BareSliceMatrix<> dshapes) const = 0;
    
    // rows dim*ndof, cols .. nip
    // rows:  phi0/dx, phi0/dy, phi0/dz, phi1/dx ... 
    HD NGS_DLL_HEADER 
    virtual void CalcMappedDShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                   BareSliceMatrix<SIMD<double>> dshapes) const;


    /**
       returns second derivatives in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    virtual const FlatMatrix<> GetDDShape (const IntegrationPoint & ip, LocalHeap & lh) const = 0;

    /// compute dshape, matrix: ndof x (spacedim spacedim)
    NGS_DLL_HEADER virtual void CalcDDShape (const IntegrationPoint & ip, 
                                             BareSliceMatrix<> ddshape) const = 0;
    
    /// compute dshape, matrix: ndof x (spacedim spacedim)
    NGS_DLL_HEADER virtual void CalcMappedDDShape (const BaseMappedIntegrationPoint & mip, 
                                                   BareSliceMatrix<> ddshape) const = 0;

    /// compute ddshape, matrix: ndof x (spacedim spacedim)
    NGS_DLL_HEADER virtual void CalcMappedDDShape (const SIMD<BaseMappedIntegrationPoint> & mip, 
                                                   BareSliceMatrix<SIMD<double>> ddshape) const = 0;


    
    /**
       Evaluates function in integration point ip.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual double Evaluate (const IntegrationPoint & ip, BareSliceVector<> x) const;
    HD NGS_DLL_HEADER virtual Complex Evaluate (const IntegrationPoint & ip, BareSliceVector<Complex> x) const;    


    /**
       Evaluate function in points of integrationrule ir.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, BareSliceVector<> coefs, BareSliceVector<> values) const;
    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs, BareVector<SIMD<double>> values) const;
    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_IntegrationRule & ir, SliceMatrix<> coefs, BareSliceMatrix<SIMD<double>> values) const;
    NGS_DLL_HEADER void Evaluate (const SIMD_IntegrationRule & ir, BareSliceVector<Complex> coefs, BareVector<SIMD<Complex>> values) const;
    /**
       Each column a vector ...
     */
    HD NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, SliceMatrix<> coefs, BareSliceMatrix<> values) const;
    
    /**
       Evaluate function in points of integrationrule ir, transpose operation.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual void EvaluateTrans (const IntegrationRule & ir, BareSliceVector<> values, BareSliceVector<> coefs) const;
    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_IntegrationRule & ir, BareVector<SIMD<double>> values, BareSliceVector<> coefs) const;
    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_IntegrationRule & ir, BareSliceMatrix<SIMD<double>> values, SliceMatrix<> coefs) const;
    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_IntegrationRule & ir, BareVector<SIMD<Complex>> values, BareSliceVector<Complex> coefs) const;

    /**
       Evaluate gradient in points of integrationrule ir.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual void EvaluateGrad (const IntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<> values) const = 0;
    
    
    HD NGS_DLL_HEADER virtual void EvaluateGrad (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const;
    HD NGS_DLL_HEADER virtual void EvaluateGrad (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const;
    // needed for ALE-trafo
    HD NGS_DLL_HEADER virtual void EvaluateGrad (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const;

    /**
       Evaluate gradient in points of integrationrule ir, transpose operation.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, BareSliceMatrix<> values, BareSliceVector<> coefs) const = 0;
    HD NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, SliceMatrix<> values, SliceMatrix<> coefs) const = 0;
    HD NGS_DLL_HEADER virtual void AddGradTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                                 BareSliceVector<> coefs) const;
    /// input du1/dx du1/dy du1/dz du2/dx ...
    /// output: ndof x components
    HD NGS_DLL_HEADER virtual void AddGradTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                                 SliceMatrix<> coefs) const;


    NGS_DLL_HEADER virtual void CalcDualShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const;
    NGS_DLL_HEADER virtual void AddDualTrans (const IntegrationRule & ir, BareSliceVector<double> values, BareSliceVector<> coefs) const;
    NGS_DLL_HEADER virtual void AddDualTrans (const SIMD_IntegrationRule & ir, BareVector<SIMD<double>> values, BareSliceVector<> coefs) const;
    
    HD NGS_DLL_HEADER virtual void GetDiagMassMatrix (FlatVector<> mass) const;
    NGS_DLL_HEADER virtual bool GetDiagDualityMassInverse (FlatVector<> diag) const { return false; }
    NGS_DLL_HEADER virtual bool DualityMassDiagonal () const { return false; }
  };

  /**
     Scalar finite element.
     Provides shape functions and derivatives.
  */
  template <int D>
  class ScalarFiniteElement : public BaseScalarFiniteElement
  {
  public:
    using BaseScalarFiniteElement::BaseScalarFiniteElement;

    HD int Dim () const final { return D; } 
    
    /**
       returns derivatives in point ip.
    */
    INLINE const FlatMatrixFixWidth<D> 
    GetDShape (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      FlatMatrixFixWidth<D> dshape(ndof, lh);
      CalcDShape (ip, dshape);
      return dshape;
    }

    using BaseScalarFiniteElement::CalcMappedDShape;    

    /// compute dshape, matrix: ndof x spacedim
    HD NGS_DLL_HEADER 
    virtual void CalcMappedDShape (const BaseMappedIntegrationPoint & mip, 
                                   BareSliceMatrix<> dshape) const override;


    HD NGS_DLL_HEADER 
    virtual void CalcMappedDShape (const BaseMappedIntegrationRule & mir, 
                                   BareSliceMatrix<> dshapes) const override;


    /**
       returns second derivatives in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatMatrix<> GetDDShape (const IntegrationPoint & ip, LocalHeap & lh) const override
    {
      FlatMatrix<> ddshape(ndof, D*D, lh);
      CalcDDShape (ip, ddshape);
      return ddshape;
    }

    /// compute dshape, matrix: ndof x (spacedim spacedim)
    NGS_DLL_HEADER virtual void CalcDDShape (const IntegrationPoint & ip, 
                                             BareSliceMatrix<> ddshape) const override;
    
    /// compute ddshape, matrix: ndof x (spacedim spacedim)
    NGS_DLL_HEADER virtual void CalcMappedDDShape (const BaseMappedIntegrationPoint & mip, 
                                                   BareSliceMatrix<> ddshape) const override;

    /// compute ddshape, matrix: ndof x (spacedim spacedim)
    NGS_DLL_HEADER virtual void CalcMappedDDShape (const SIMD<BaseMappedIntegrationPoint> & mip, 
                                                   BareSliceMatrix<SIMD<double>> ddshape) const override;


    /**
       Evaluates gradient in integration point ip.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual Vec<D> EvaluateGrad (const IntegrationPoint & ip, BareSliceVector<> x) const;

    using BaseScalarFiniteElement::EvaluateGrad;

    /**
       Evaluate gradient in points of integrationrule ir.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER void EvaluateGrad (const IntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<> values) const override;
    
    /**
       Evaluate gradient in points of integrationrule ir, transpose operation.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER void EvaluateGradTrans (const IntegrationRule & ir, BareSliceMatrix<> values, BareSliceVector<> coefs) const override;

    HD NGS_DLL_HEADER void EvaluateGradTrans (const IntegrationRule & ir, SliceMatrix<> values, SliceMatrix<> coefs) const override;

    NGS_DLL_HEADER virtual void Interpolate (const ElementTransformation & trafo, 
                                             const class CoefficientFunction & func, SliceMatrix<> coefs,
                                             LocalHeap & lh) const override;
    
  public:
    NGS_DLL_HEADER virtual std::list<std::tuple<std::string,double>> Timing () const override;
    NGS_DLL_HEADER virtual bool SolveDuality (SliceVector<> rhs, SliceVector<> u, LocalHeap & lhr) const override;
  };




















  template<ELEMENT_TYPE ET>
  class DGFiniteElement : public ScalarFiniteElement<ET_trait<ET>::DIM>,
                          public VertexOrientedFE<ET>
  {
  protected:
    // int vnums[1<<D];  
    enum { D = ET_trait<ET>::DIM };
    
    using ScalarFiniteElement<D>::ndof;
    using ScalarFiniteElement<D>::order;
    using VertexOrientedFE<ET>::vnums;
    
  public:
    /// global vertex numbers define ordering of vertices
    template <typename TA>
    void SetVertexNumbers (const TA & avnums)
    { 
      for (int i = 0; i < avnums.Size(); i++) vnums[i] = avnums[i]; 
    }
    DGFiniteElement * SetVertexNumbers (FlatArray<int> vnums) override
    { VertexOrientedFE<ET>::SetVertexNumbers(vnums); return this; }
    

    /// assign vertex number
    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }
    NGS_DLL_HEADER virtual void SetOrder (IVec<D> p) = 0;
    NGS_DLL_HEADER virtual void ComputeNDof() = 0;


    NGS_DLL_HEADER virtual void PrecomputeTrace () = 0; 
    NGS_DLL_HEADER virtual void PrecomputeGrad () = 0;

    NGS_DLL_HEADER void CalcTraceMatrix (int facet, FlatMatrix<> trace) const;
    NGS_DLL_HEADER void CalcGradientMatrix (FlatMatrix<> gmat) const;

    HD NGS_DLL_HEADER virtual void GetDiagMassMatrix (FlatVector<> mass) const override;

    NGS_DLL_HEADER virtual void GetGradient (FlatVector<> coefs, FlatMatrixFixWidth<D> grad) const;
    NGS_DLL_HEADER virtual void GetGradientTrans (FlatMatrixFixWidth<D> grad, FlatVector<> coefs) const;

    NGS_DLL_HEADER virtual void GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const;
    NGS_DLL_HEADER virtual void GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const;
  };
  



  extern template class ScalarFiniteElement<0>;
  extern template class ScalarFiniteElement<1>;
  extern template class ScalarFiniteElement<2>;
  extern template class ScalarFiniteElement<3>;

  extern template class DGFiniteElement<ET_POINT>;
  extern template class DGFiniteElement<ET_SEGM>;
  extern template class DGFiniteElement<ET_TRIG>;
  extern template class DGFiniteElement<ET_QUAD>;
  extern template class DGFiniteElement<ET_TET>;
  extern template class DGFiniteElement<ET_PRISM>;
  extern template class DGFiniteElement<ET_PYRAMID>;
  extern template class DGFiniteElement<ET_HEXAMID>;
  extern template class DGFiniteElement<ET_HEX>;
}

#endif
