#ifndef FILE_L2HOFE
#define FILE_L2HOFE

/*********************************************************************/
/* File:   l2hofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#include "tscalarfe.hpp"
#include "precomp.hpp"


namespace ngfem
{



  template <ELEMENT_TYPE ET> class L2HighOrderFE_Shape;

  template <int DIM>
  class PrecomputedScalShapes
  {
  public:
    Matrix<> shapes;
    Matrix<> dshapes;
    
    PrecomputedScalShapes (int nip, int ndof)
      : shapes(nip, ndof), dshapes (DIM*nip, ndof)
    { ; }

  };


  template <ELEMENT_TYPE ET, 
            template <ELEMENT_TYPE ET> class SHAPES = L2HighOrderFE_Shape,
            class BASE = DGFiniteElement<ET_trait<ET>::DIM> > 
  class NGS_DLL_HEADER L2HighOrderFE :  
    public T_ScalarFiniteElement2<SHAPES<ET>, ET, BASE >,
    public ET_trait<ET>
  { 
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using ScalarFiniteElement<DIM>::eltype;
    using DGFiniteElement<DIM>::vnums;


    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::PolDimension;

    INT<DIM> order_inner; 

    typedef PrecomputedShapesContainer<PrecomputedScalShapes<DIM> > TPRECOMP;
    static TPRECOMP precomp;

    typedef HashTable<INT<2>, Matrix<>*> TPRECOMP_TRACE;
    static TPRECOMP_TRACE precomp_trace;

    typedef HashTable<INT<2>, Matrix<>*> TPRECOMP_GRAD;
    static TPRECOMP_GRAD precomp_grad;

  public:
    L2HighOrderFE () 
    { 
      eltype = ET; 
    }

    L2HighOrderFE (int aorder)
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++) vnums[i] = i;
      eltype = ET;

      order = aorder;
      order_inner = aorder;
      ndof = PolDimension (aorder);
    }

    /// global vertex numbers define ordering of vertices
    template <typename TA>
    void SetVertexNumbers (const TA & avnums)
    { for (int i = 0; i < N_VERTEX; i++) vnums[i] = avnums[i]; }

    /// different orders in differnt directions
    virtual void SetOrder (INT<DIM> p)  { order_inner = p; }

    virtual void ComputeNDof()
    {
      ndof = PolDimension (order_inner);
      order = 0;
      for (int i = 0; i < DIM; i++)
        order = max(order, order_inner[i]);
    }

    virtual void PrecomputeTrace ();
    virtual void PrecomputeGrad ();
    virtual void PrecomputeShapes (const IntegrationRule & ir);

    virtual void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const;
    virtual void EvaluateTrans (const IntegrationRule & ir, FlatVector<> values, FlatVector<> coefs) const;

    virtual void EvaluateGrad (const IntegrationRule & ir, FlatVector<> coefs, FlatMatrixFixWidth<DIM> values) const;
    virtual void EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> values, FlatVector<> coefs) const;

    virtual void GetGradient (FlatVector<> coefs, FlatMatrixFixWidth<DIM> grad) const;
    virtual void GetGradientTrans (FlatMatrixFixWidth<DIM> grad, FlatVector<> coefs) const;

    virtual void GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const;
    virtual void GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const;

    virtual void GetDiagMassMatrix (FlatVector<> mass) const;
  };





}

#endif
