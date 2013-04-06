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

  /**
     Base class for L2 - high order finite elements, i.e., a
     discontinuous approximation
  */
  template<int D>
  class L2HighOrderFiniteElement : virtual public ScalarFiniteElement<D>
  {
  protected:

    enum { DIM = D };
    int vnums[8];  
    INT<DIM> order_inner; 

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using ScalarFiniteElement<DIM>::eltype;


  public:
    /// global vertex numbers define ordering of vertices
    template <typename TA>
    void SetVertexNumbers (const TA & avnums)
    { for (int i = 0; i < avnums.Size(); i++) vnums[i] = avnums[i]; }

    /// assign vertex number
    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }

    /// different orders in differnt directions
    void SetOrder (INT<DIM> p) { order_inner = p; }

    /// calculate number of dofs
    virtual void ComputeNDof () = 0; 
  
    virtual void PrecomputeTrace () = 0; 
    virtual void PrecomputeGrad () = 0;

    void CalcTraceMatrix (int facet, FlatMatrix<> trace) const;
    void CalcGradientMatrix (FlatMatrix<> gmat) const;

    virtual void GetDiagMassMatrix (FlatVector<> mass) const;

    virtual void GetGradient (FlatVector<> coefs, FlatMatrixFixWidth<DIM> grad) const
    {
      Matrix<> gmat(DIM*grad.Height(), coefs.Size());
      CalcGradientMatrix (gmat);
      FlatVector<> vgrad(gmat.Height(), &grad(0,0));
      vgrad = gmat * coefs;
    }

    virtual void GetGradientTrans (FlatMatrixFixWidth<DIM> grad, FlatVector<> coefs) const 
    {
      Matrix<> gmat(DIM*grad.Height(), coefs.Size());
      CalcGradientMatrix (gmat);
      FlatVector<> vgrad(gmat.Height(), &grad(0,0));
      coefs = Trans (gmat) * vgrad;
    }

    virtual void GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const
    {
      Matrix<> trace(fcoefs.Size(), coefs.Size());
      CalcTraceMatrix(facet, trace);
      fcoefs = trace * coefs;
    }

    virtual void GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const
    {
      Matrix<> trace(fcoefs.Size(), coefs.Size());
      CalcTraceMatrix(facet, trace);
      coefs = Trans (trace) * fcoefs;
    }
  };



  /**
     Template family of L2 - high order finite elements.
     The template argument is the element shape
  */
  // template <ELEMENT_TYPE ET> class L2HighOrderFE;


  /**
     Barton-Nackman base class for L2 - high order finite elements
  */
  template <ELEMENT_TYPE ET>
  class T_L2HighOrderFiniteElement : 
    public L2HighOrderFiniteElement<ET_trait<ET>::DIM>,
    public ET_trait<ET> 
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using ScalarFiniteElement<DIM>::eltype;

    using L2HighOrderFiniteElement<DIM>::vnums;
    using L2HighOrderFiniteElement<DIM>::order_inner;

    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::N_EDGE;
    using ET_trait<ET>::N_FACE;
    using ET_trait<ET>::FaceType;
    using ET_trait<ET>::GetEdgeSort;
    using ET_trait<ET>::GetFaceSort;
    using ET_trait<ET>::PolDimension;


  public:

    T_L2HighOrderFiniteElement () 
    { eltype = ET; }

    T_L2HighOrderFiniteElement (int aorder) 
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++) vnums[i] = i;
      eltype = ET;

      order = aorder;
      order_inner = aorder;
      ndof = PolDimension (aorder);
    }

    virtual void ComputeNDof();
  };






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


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES = L2HighOrderFE_Shape> 
  class NGS_DLL_HEADER L2HighOrderFE :  public T_L2HighOrderFiniteElement<ET>,
					public T_ScalarFiniteElement2<SHAPES<ET>, ET >
  { 
  protected:
    using T_L2HighOrderFiniteElement<ET>::ndof;
    using T_L2HighOrderFiniteElement<ET>::order;
    using T_L2HighOrderFiniteElement<ET>::vnums;

    enum { DIM = ET_trait<ET>::DIM };
    typedef PrecomputedShapesContainer<PrecomputedScalShapes<DIM> > TPRECOMP;
    static TPRECOMP precomp;

    typedef HashTable<INT<2>, Matrix<>*> TPRECOMP_TRACE;
    static TPRECOMP_TRACE precomp_trace;

    typedef HashTable<INT<2>, Matrix<>*> TPRECOMP_GRAD;
    static TPRECOMP_GRAD precomp_grad;

  public:
    L2HighOrderFE () { ; }

    L2HighOrderFE (int aorder)
      : T_L2HighOrderFiniteElement<ET> (aorder) { ; }


    virtual void PrecomputeTrace ();
    virtual void PrecomputeGrad ();
    virtual void PrecomputeShapes (const IntegrationRule & ir);

    virtual void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const;

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
