#ifndef FILE_SCALARFE
#define FILE_SCALARFE

/*********************************************************************/
/* File:   scalarfe.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

#include "avector.hpp"

namespace ngfem
{
  
  class BaseScalarFiniteElement : public FiniteElement 
  {
  public:
    // using FiniteElement::FiniteElement;

    INLINE BaseScalarFiniteElement () { ; } 
    INLINE BaseScalarFiniteElement (int andof, int aorder)
      : FiniteElement (andof, aorder) { ; }


    /// compute shape
    HD NGS_DLL_HEADER 
    virtual void CalcShape (const IntegrationPoint & ip, 
                            SliceVector<> shape) const = 0;

    /// compute dshape, matrix: ndof x spacedim
    HD NGS_DLL_HEADER 
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     SliceMatrix<> dshape) const = 0;
  };

  /**
     Scalar finite element.
     Provides shape functions and derivaties.
  */
  template <int D>
  class ScalarFiniteElement : public BaseScalarFiniteElement
  {
  public:
    using BaseScalarFiniteElement::BaseScalarFiniteElement;

    /// the name
    NGS_DLL_HEADER virtual string ClassName() const;

    HD NGS_DLL_HEADER virtual int Dim () const { return D; }

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

    using  BaseScalarFiniteElement::CalcShape;
    using  BaseScalarFiniteElement::CalcDShape;

    /// compute shape, row is shape nr, col is ip nr
    HD NGS_DLL_HEADER 
    virtual void CalcShape (const IntegrationRule & ir, 
                            SliceMatrix<> shape) const;
  
    
/*    
    NGS_DLL_HEADER 
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     const std::function<void(int,Vec<D>)> & callback) const;
  */  


    /// compute dshape, matrix: ndof x spacedim
    HD NGS_DLL_HEADER 
    virtual void CalcMappedDShape (const MappedIntegrationPoint<D,D> & mip, 
                                   SliceMatrix<> dshape) const;


    HD NGS_DLL_HEADER 
    virtual void CalcMappedDShape (const MappedIntegrationRule<D,D> & mir, 
                                   SliceMatrix<> dshapes) const;



    /*
    template <typename ANY_MIP, typename T>
    INLINE void CalcMappedDShape (const ANY_MIP & mip, MatExpr<T> & mat) const
    {
      cout << "calc dshape from any matrix, type = " << typeid(T).name() << endl;
      CalcDShape (mip.IP(), 
		  [&](int i, Vec<3> gradref)
		  { mat.Row(i) = Trans(mip.GetJacobianInverse()) * gradref; });
    }
    */



    /**
       returns second derivatives in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatMatrix<> GetDDShape (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      FlatMatrix<> ddshape(ndof, D*D, lh);
      CalcDDShape (ip, ddshape);
      return ddshape;
    }

    /// compute dshape, matrix: ndof x (spacedim spacedim)
    NGS_DLL_HEADER virtual void CalcDDShape (const IntegrationPoint & ip, 
                                             FlatMatrix<> ddshape) const;
    
    /// compute dshape, matrix: ndof x (spacedim spacedim)
    NGS_DLL_HEADER virtual void CalcMappedDDShape (const MappedIntegrationPoint<D,D> & mip, 
                                                   FlatMatrix<> ddshape) const;


    /**
       Evaluates function in integration point ip.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual double Evaluate (const IntegrationPoint & ip, SliceVector<> x) const;

    /**
       Evaluates gradient in integration point ip.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual Vec<D> EvaluateGrad (const IntegrationPoint & ip, SliceVector<> x) const;

    
    /**
       Evaluate function in points of integrationrule ir.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, SliceVector<> coefs, FlatVector<> values) const;
    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_IntegrationRule & ir, SliceVector<> coefs, AFlatVector<double> values) const;
    /**
       Each column a vector ...
     */
    HD NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, SliceMatrix<> coefs, SliceMatrix<> values) const;

    /**
       Evaluate gradient in points of integrationrule ir.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual void EvaluateGrad (const IntegrationRule & ir, SliceVector<> coefs, FlatMatrixFixWidth<D> values) const;
    

    /**
       Evaluate function in points of integrationrule ir, transpose operation.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual void EvaluateTrans (const IntegrationRule & ir, FlatVector<> values, SliceVector<> coefs) const;
    HD NGS_DLL_HEADER virtual void EvaluateTrans (const SIMD_IntegrationRule & ir, AFlatVector<double> values, SliceVector<> coefs) const;

    /**
       Evaluate gradient in points of integrationrule ir, transpose operation.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<D> values, SliceVector<> coefs) const;

    HD NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, SliceMatrix<> values, SliceMatrix<> coefs) const;


    HD NGS_DLL_HEADER virtual void GetPolOrders (FlatArray<PolOrder<D> > orders) const;
  };




















  template<int D>
  class DGFiniteElement : public ScalarFiniteElement<D>
  {
  protected:
    int vnums[1<<D];  

    using ScalarFiniteElement<D>::ndof;
    using ScalarFiniteElement<D>::order;

  public:
    /// global vertex numbers define ordering of vertices
    template <typename TA>
    void SetVertexNumbers (const TA & avnums)
    { 
      for (int i = 0; i < avnums.Size(); i++) vnums[i] = avnums[i]; 
    }

    /// assign vertex number
    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }
    NGS_DLL_HEADER virtual void SetOrder (INT<D> p) = 0;
    NGS_DLL_HEADER virtual void ComputeNDof() = 0;


    NGS_DLL_HEADER virtual void PrecomputeTrace () = 0; 
    NGS_DLL_HEADER virtual void PrecomputeGrad () = 0;

    NGS_DLL_HEADER void CalcTraceMatrix (int facet, FlatMatrix<> trace) const;
    NGS_DLL_HEADER void CalcGradientMatrix (FlatMatrix<> gmat) const;

    HD NGS_DLL_HEADER virtual void GetDiagMassMatrix (FlatVector<> mass) const;

    NGS_DLL_HEADER virtual void GetGradient (FlatVector<> coefs, FlatMatrixFixWidth<D> grad) const;
    NGS_DLL_HEADER virtual void GetGradientTrans (FlatMatrixFixWidth<D> grad, FlatVector<> coefs) const;

    NGS_DLL_HEADER virtual void GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const;
    NGS_DLL_HEADER virtual void GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const;
  };
  



#ifdef FILE_SCALARFE_CPP
#define SCALARFE_EXTERN
#else
#define SCALARFE_EXTERN extern

  SCALARFE_EXTERN template class ScalarFiniteElement<0>;
  SCALARFE_EXTERN template class ScalarFiniteElement<1>;
  SCALARFE_EXTERN template class ScalarFiniteElement<2>;
  SCALARFE_EXTERN template class ScalarFiniteElement<3>;

  SCALARFE_EXTERN template class DGFiniteElement<0>;
  SCALARFE_EXTERN template class DGFiniteElement<1>;
  SCALARFE_EXTERN template class DGFiniteElement<2>;
  SCALARFE_EXTERN template class DGFiniteElement<3>;

#endif
}

#endif
