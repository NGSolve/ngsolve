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
	    class SHAPES = L2HighOrderFE_Shape<ET>,
	    class BASE = T_ScalarFiniteElement<SHAPES, ET, DGFiniteElement<ET_trait<ET>::DIM> > >
	    
  class L2HighOrderFE : public BASE, public ET_trait<ET>
  { 
  protected:
    typedef BASE T_IMPL;
    typedef SHAPES T_SHAPES;
    
    using ET_trait<ET>::DIM;
    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::PolDimension;

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using DGFiniteElement<DIM>::vnums;


    INT<DIM> order_inner; 

    typedef PrecomputedShapesContainer<PrecomputedScalShapes<DIM> > TPRECOMP;
    static TPRECOMP precomp;

    typedef HashTable<INT<2>, Matrix<>*> TPRECOMP_TRACE;
    static TPRECOMP_TRACE precomp_trace;

    typedef HashTable<INT<2>, Matrix<>*> TPRECOMP_GRAD;
    static TPRECOMP_GRAD precomp_grad;

  public:
    INLINE L2HighOrderFE () { ; }
    INLINE L2HighOrderFE (int aorder)
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++) vnums[i] = i;
      order = aorder;
      order_inner = aorder;
      ndof = PolDimension (aorder);
    }

    virtual NGS_DLL_HEADER ~L2HighOrderFE () { ; }
    
    /// global vertex numbers define ordering of vertices
    template <typename TA>
    INLINE void SetVertexNumbers (const TA & avnums)
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

    NGS_DLL_HEADER virtual void PrecomputeTrace ();
    NGS_DLL_HEADER virtual void PrecomputeGrad ();
    NGS_DLL_HEADER virtual void PrecomputeShapes (const IntegrationRule & ir);

    NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const;
    NGS_DLL_HEADER virtual void EvaluateTrans (const IntegrationRule & ir, FlatVector<> values, FlatVector<> coefs) const;

    NGS_DLL_HEADER virtual void EvaluateGrad (const IntegrationRule & ir, FlatVector<> coefs, FlatMatrixFixWidth<DIM> values) const;
    NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> values, FlatVector<> coefs) const;

    NGS_DLL_HEADER virtual void GetGradient (FlatVector<> coefs, FlatMatrixFixWidth<DIM> grad) const;
    NGS_DLL_HEADER virtual void GetGradientTrans (FlatMatrixFixWidth<DIM> grad, FlatVector<> coefs) const;

    NGS_DLL_HEADER virtual void GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const;
    NGS_DLL_HEADER virtual void GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const;

    NGS_DLL_HEADER virtual void GetDiagMassMatrix (FlatVector<> mass) const;
  };

}


  
#ifdef FILE_L2HOFE_CPP
#define L2HOFE_EXTERN

#include <tscalarfe_impl.hpp>
#include <l2hofe_impl.hpp>

namespace ngfem
{
  template <> inline void L2HighOrderFE<ET_POINT> ::
  GetDiagMassMatrix(FlatVector<> mass) const
  {
    mass(0) = 1;
  }
  
  template <> inline void L2HighOrderFE<ET_SEGM> ::
  GetDiagMassMatrix(FlatVector<> mass) const
  {
    for (int ix = 0; ix <= order; ix++)
      mass(ix) = 1.0 / (2 * ix + 1);
  }
  
  template <> inline void L2HighOrderFE<ET_TRIG> ::
  GetDiagMassMatrix(FlatVector<> mass) const
  {
    for (int ix = 0, ii = 0; ix <= order; ix++)
      for (int iy = 0; iy <= order - ix; iy++, ii++)
	mass(ii) = 1.0 / ((2 * ix + 1) * (2 * ix + 2 * iy + 2));
  }
}

#else
#define L2HOFE_EXTERN extern
#endif




namespace ngfem
{
  L2HOFE_EXTERN template class L2HighOrderFE<ET_POINT>;
  L2HOFE_EXTERN template class L2HighOrderFE<ET_SEGM>;
  L2HOFE_EXTERN template class L2HighOrderFE<ET_TRIG>;
  L2HOFE_EXTERN template class L2HighOrderFE<ET_QUAD>;

  L2HOFE_EXTERN template class L2HighOrderFE<ET_TET>;
  L2HOFE_EXTERN template class L2HighOrderFE<ET_PRISM>;
  L2HOFE_EXTERN template class L2HighOrderFE<ET_PYRAMID>;
  L2HOFE_EXTERN template class L2HighOrderFE<ET_HEX>;

  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_POINT>, ET_POINT, DGFiniteElement<0> >;
  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_SEGM>, ET_SEGM, DGFiniteElement<1> >;
  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_TRIG>, ET_TRIG, DGFiniteElement<2> >;
  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_QUAD>, ET_QUAD, DGFiniteElement<2> >;

  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_TET>, ET_TET, DGFiniteElement<3> >;
  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_PRISM>, ET_PRISM, DGFiniteElement<3> >;
  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_PYRAMID>, ET_PYRAMID, DGFiniteElement<3> >;
  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_HEX>, ET_HEX, DGFiniteElement<3> >;

}

#endif
