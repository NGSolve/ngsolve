#ifndef FILE_SCALARFE
#define FILE_SCALARFE

/*********************************************************************/
/* File:   scalarfe.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngfem
{


  /**
     Scalar finite element.
     Provides shape functions and derivaties.
  */
  template <int D>
  class ScalarFiniteElement : public FiniteElement
  {
  public:
    /// empty constructor
    NGS_DLL_HEADER ScalarFiniteElement () { ; } 
    /// provides type, number of dofs, maximal order of shapes
    NGS_DLL_HEADER ScalarFiniteElement (ELEMENT_TYPE aeltype, 
			 int andof = 0, int aorder = 0)
      : FiniteElement (aeltype, andof, aorder) 
    { ; }

    // destructor
    NGS_DLL_HEADER virtual ~ScalarFiniteElement () = 0;
    
    /// the name
    virtual string ClassName() const;

    /**
       returns shape functions in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    FlatVector<> GetShape (const IntegrationPoint & ip, 
                           LocalHeap & lh) const
    {
      FlatVector<> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

    /**
       returns derivatives in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatMatrixFixWidth<D> 
    GetDShape (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      FlatMatrixFixWidth<D> dshape(ndof, lh);
      CalcDShape (ip, dshape);
      return dshape;
    }


    /// compute shape
    NGS_DLL_HEADER virtual void CalcShape (const IntegrationPoint & ip, 
                                           FlatVector<> shape) const = 0;
  
    /// compute dshape, matrix: ndof x spacedim
    NGS_DLL_HEADER virtual void CalcDShape (const IntegrationPoint & ip, 
                                            FlatMatrixFixWidth<D> dshape) const;

    /// compute dshape, matrix: ndof x spacedim
    NGS_DLL_HEADER 
    virtual void CalcMappedDShape (const MappedIntegrationPoint<D,D> & mip, 
                                   FlatMatrixFixWidth<D> dshape) const;
    

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



    /**
       Evaluates function in integration point ip.
       Vector x provides coefficient vector.
     */
    NGS_DLL_HEADER virtual double Evaluate (const IntegrationPoint & ip, FlatVector<> x) const;

    /**
       Evaluates gradient in integration point ip.
       Vector x provides coefficient vector.
     */
    NGS_DLL_HEADER virtual Vec<D> EvaluateGrad (const IntegrationPoint & ip, FlatVector<> x) const;

    
    /**
       Evaluate function in points of integrationrule ir.
       Vector x provides coefficient vector.
     */
    NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, FlatVector<> coefs, FlatVector<> values) const;

    /**
       Evaluate gradient in points of integrationrule ir.
       Vector x provides coefficient vector.
     */
    NGS_DLL_HEADER virtual void EvaluateGrad (const IntegrationRule & ir, FlatVector<> coefs, FlatMatrixFixWidth<D> values) const;


    /**
       Evaluate function in points of integrationrule ir, transpose operation.
       Vector x provides coefficient vector.
     */
    NGS_DLL_HEADER virtual void EvaluateTrans (const IntegrationRule & ir, FlatVector<> values, FlatVector<> coefs) const;


    /**
       Evaluate gradient in points of integrationrule ir, transpose operation.
       Vector x provides coefficient vector.
     */
    NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<D> values, FlatVector<> coefs) const;

    NGS_DLL_HEADER virtual void GetPolOrders (FlatArray<PolOrder<D> > orders) const;
  };




















  template<int D>
  class DGFiniteElement : public ScalarFiniteElement<D>
  {
  protected:
    int vnums[1<<D];  

    using ScalarFiniteElement<D>::ndof;
    using ScalarFiniteElement<D>::order;
    // using ScalarFiniteElement<D>::eltype;

  public:
    /// global vertex numbers define ordering of vertices
    template <typename TA>
    void SetVertexNumbers (const TA & avnums)
    { for (int i = 0; i < avnums.Size(); i++) vnums[i] = avnums[i]; }

    /// assign vertex number
    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }
    NGS_DLL_HEADER virtual void SetOrder (INT<D> p) = 0;
    NGS_DLL_HEADER virtual void ComputeNDof() = 0;


    NGS_DLL_HEADER virtual void PrecomputeTrace () = 0; 
    NGS_DLL_HEADER virtual void PrecomputeGrad () = 0;

    NGS_DLL_HEADER void CalcTraceMatrix (int facet, FlatMatrix<> trace) const;
    NGS_DLL_HEADER void CalcGradientMatrix (FlatMatrix<> gmat) const;

    NGS_DLL_HEADER virtual void GetDiagMassMatrix (FlatVector<> mass) const;

    virtual void GetGradient (FlatVector<> coefs, FlatMatrixFixWidth<D> grad) const;
    virtual void GetGradientTrans (FlatMatrixFixWidth<D> grad, FlatVector<> coefs) const;

    virtual void GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const;
    virtual void GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const;
  };
  






}

#endif
