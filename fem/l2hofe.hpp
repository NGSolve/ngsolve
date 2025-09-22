#ifndef FILE_L2HOFE
#define FILE_L2HOFE

/*********************************************************************/
/* File:   l2hofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#include "tscalarfe.hpp"
#include "precomp.hpp"

#include "recursive_pol.hpp"

namespace ngfem
{

  template <ELEMENT_TYPE ET>
  ScalarFiniteElement<Dim(ET)> * CreateL2HighOrderFE (int order, FlatArray<int> vnums, Allocator & lh);


  template <ELEMENT_TYPE ET> class L2HighOrderFE_Shape;

  template <int DIM>
  class PrecomputedScalShapes
  {
  public:
    Matrix<> shapes;
    Matrix<> dshapes;
    
    INLINE PrecomputedScalShapes (int nip, int ndof)
      : shapes(nip, ndof), dshapes (DIM*nip, ndof)
    { ; }

  };


  template <ELEMENT_TYPE ET, 
	    class SHAPES = L2HighOrderFE_Shape<ET>,
	    class BASE = T_ScalarFiniteElement<SHAPES, ET, DGFiniteElement<ET> > >
	    
  class L2HighOrderFE : public BASE, public ET_trait<ET>
  { 
  protected:
    typedef BASE T_IMPL;
    typedef SHAPES T_SHAPES;

    static constexpr int DIM = ngfem::Dim(ET);
    // enum { DIM = ET_trait<ET>::DIM };
    enum { N_VERTEX = ET_trait<ET>::N_VERTEX };

    using ET_trait<ET>::PolDimension;
    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using DGFiniteElement<ET>::vnums;


    IVec<DIM> order_inner; 

#ifndef __CUDA_ARCH__
    typedef PrecomputedShapesContainer<PrecomputedScalShapes<DIM> > TPRECOMP;
    static TPRECOMP precomp;

    typedef HashTable<IVec<2>, Matrix<>*> TPRECOMP_TRACE;
    static TPRECOMP_TRACE precomp_trace;

    typedef HashTable<IVec<2>, Matrix<>*> TPRECOMP_GRAD;
    static TPRECOMP_GRAD precomp_grad;
#endif

  public:
    using ET_trait<ET>::ElementType;
    using BaseScalarFiniteElement::CalcShape;
    
    INLINE L2HighOrderFE () { ; }
    INLINE L2HighOrderFE (int aorder)
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++) vnums[i] = i;
      order = aorder;
      order_inner = aorder;
      ndof = PolDimension (aorder);
    }

    // virtual NGS_DLL_HEADER ~L2HighOrderFE () { ; }
    
    /// global vertex numbers define ordering of vertices
    template <typename TA>
    INLINE void SetVertexNumbers (const TA & avnums)
    { for (int i = 0; i < N_VERTEX; i++) vnums[i] = avnums[i]; }

    /// different orders in different directions
    virtual void SetOrder (IVec<DIM> p) override { order_inner = p; }

    virtual void ComputeNDof() override
    {
      ndof = PolDimension (order_inner);
      order = 0;
      for (int i = 0; i < DIM; i++)
        order = max2(order, order_inner[i]);
    }

    HD virtual tuple<int,int,int,int> GetNDofVEFC () const override
    {
      switch (DIM)
        {
        case 0: return { 1, 0, 0, 0 };
        case 1: return { 0, ndof, 0, 0 };
        case 2: return { 0, 0, ndof, 0 };
        case 3: return { 0, 0, 0, ndof };
        }
      return { 0, 0, 0, 0 }; // for the compiler
    }

    NGS_DLL_HEADER virtual void PrecomputeTrace () override;
    NGS_DLL_HEADER virtual void PrecomputeGrad () override;
    NGS_DLL_HEADER virtual void PrecomputeShapes (const IntegrationRule & ir) override;

    using BASE::Evaluate;
    HD NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, BareSliceVector<double> coefs, FlatVector<double> vals) const;
    HD NGS_DLL_HEADER virtual void EvaluateTrans (const IntegrationRule & ir, BareSliceVector<> values, BareSliceVector<> coefs) const override;

    using BASE::EvaluateGrad;    
    HD NGS_DLL_HEADER virtual void EvaluateGrad (const IntegrationRule & ir, BareSliceVector<> coefs, FlatMatrixFixWidth<DIM> values) const;

    using BASE::EvaluateGradTrans;
    HD NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, BareSliceMatrix<> values, BareSliceVector<> coefs) const override;

    NGS_DLL_HEADER virtual void GetGradient (FlatVector<> coefs, FlatMatrixFixWidth<DIM> grad) const override;
    NGS_DLL_HEADER virtual void GetGradientTrans (FlatMatrixFixWidth<DIM> grad, FlatVector<> coefs) const override;

    NGS_DLL_HEADER virtual void GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const override;
    NGS_DLL_HEADER virtual void GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const override;

    HD NGS_DLL_HEADER virtual void GetDiagMassMatrix (FlatVector<> mass) const override;
    NGS_DLL_HEADER virtual bool DualityMassDiagonal () const override { return true; }
    NGS_DLL_HEADER virtual bool GetDiagDualityMassInverse (FlatVector<> diag) const override
    {
      GetDiagMassMatrix(diag);
      for (auto & d : diag) d = 1.0/d;
      return true;
    }
  };

}


  
#ifdef FILE_L2HOFE_CPP
#define L2HOFE_EXTERN

#include <tscalarfe_impl.hpp>
#include <l2hofe_impl.hpp>

#else
#define L2HOFE_EXTERN extern
#endif



// #ifdef FILE_L2HOFE_TRIG_CPP
// #include <tscalarfe_impl.hpp>
// #include <l2hofe_impl.hpp>
// #endif

namespace ngfem
{
  L2HOFE_EXTERN template class L2HighOrderFE<ET_POINT>;
  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_POINT>, ET_POINT, DGFiniteElement<ET_POINT> >;
  
  extern template class L2HighOrderFE<ET_SEGM>;
  extern template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_SEGM>, ET_SEGM, DGFiniteElement<ET_SEGM> >;
  
  extern template class L2HighOrderFE<ET_TRIG>;
  extern template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_TRIG>, ET_TRIG, DGFiniteElement<ET_TRIG> >;
  
  L2HOFE_EXTERN template class L2HighOrderFE<ET_QUAD>;
  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_QUAD>, ET_QUAD, DGFiniteElement<ET_QUAD> >;
  
  extern template class L2HighOrderFE<ET_TET>;
  extern template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_TET>, ET_TET, DGFiniteElement<ET_TET> >;
  
  L2HOFE_EXTERN template class L2HighOrderFE<ET_PRISM>;
  L2HOFE_EXTERN template class L2HighOrderFE<ET_PYRAMID>;
  L2HOFE_EXTERN template class L2HighOrderFE<ET_HEXAMID>;
  L2HOFE_EXTERN template class L2HighOrderFE<ET_HEX>;

  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_PRISM>, ET_PRISM, DGFiniteElement<ET_PRISM> >;
  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_PYRAMID>, ET_PYRAMID, DGFiniteElement<ET_PYRAMID> >;
  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_HEXAMID>, ET_HEXAMID, DGFiniteElement<ET_HEXAMID> >;
  L2HOFE_EXTERN template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_HEX>, ET_HEX, DGFiniteElement<ET_HEX> >;
}

#endif
