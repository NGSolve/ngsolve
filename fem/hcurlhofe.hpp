#ifndef FILE_HCURLHOFE_
#define FILE_HCURLHOFE_  

/*********************************************************************/
/* File:   hcurlhofe.hpp                                             */
/* Author: Sabine Zaglmayr                                           */
/* Date:   20. Maerz 2003                                            */
/*                                                                   */
/* AutoCurl - revision: J. Schoeberl, March 2009                     */
/*********************************************************************/


#include "thcurlfe.hpp"
#include "recursive_pol.hpp"

namespace ngfem
{



  /** 
      HCurlHighOrderFE of shape ET.
      The template specialization provides the shape functions.
  */
  template <ELEMENT_TYPE ET> class HCurlHighOrderFE_Shape;


  /**
     High order finite elements for H(curl).  
     These are the actual finite element classes to be used.  
     Shape functions are provided by the shape template
   */

  template <ELEMENT_TYPE ET, 
            template <ELEMENT_TYPE ET2> class TSHAPES = HCurlHighOrderFE_Shape,
            typename BASE = T_HCurlHighOrderFiniteElement<ET, TSHAPES<ET>> >

  class HCurlHighOrderFE : public BASE, public ET_trait<ET>, public VertexOrientedFE<ET>
  {
  protected:
    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::N_EDGE;
    using ET_trait<ET>::N_FACE;
    using ET_trait<ET>::N_CELL;
    using ET_trait<ET>::FaceType;
    using ET_trait<ET>::DIM;

    typedef short TORDER;

    using VertexOrientedFE<ET>::vnums;

    IVec<N_EDGE, TORDER> order_edge;
    IVec<N_FACE, IVec<2, TORDER>> order_face;
    IVec<3, TORDER> order_cell;
    
    //bool usegrad_edge[N_EDGE]; 
    IVec<N_EDGE, bool> usegrad_edge;
    IVec<N_FACE, bool> usegrad_face;
    bool usegrad_cell;
    bool type1;   

    using BASE::ndof;
    using BASE::order;

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

    // Typedefs should match with h1hofe, otherwise HCurlHighOrderFESpace::CreateGradient() fails
    // typedef LegendrePolynomial EdgeOrthoPol;
    typedef IntLegNoBubble EdgeOrthoPol;  // Integrated Legendre divided by bubble
    // typedef ChebyPolynomial EdgeOrthoPol;
    
    // typedef ChebyPolynomial QuadOrthoPol;
    typedef IntLegNoBubble QuadOrthoPol;

    // typedef DubinerBasis TrigOrthoPolGrad;
    typedef DubinerBasisOrthoBub TrigOrthoPolGrad;  // for bi-orthogonal dual shapes (Jan 2021)
    

  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    /* INLINE */ HCurlHighOrderFE () { ; } 

    /* INLINE */ HCurlHighOrderFE (int aorder) 
    {
      order_edge = aorder;
      order_face = aorder;
      type1 = false;
      
      if (DIM == 3) order_cell = aorder;
      
      for(int i = 0; i < N_EDGE; i++) usegrad_edge[i] = 1;
      for(int i = 0; i < N_FACE; i++) usegrad_face[i] = 1;
      if (DIM == 3) usegrad_cell = 1;
      
      for (int i = 0; i < N_VERTEX; i++) vnums[i] = i;
      
      ComputeNDof();
    }

    INLINE void SetOrderEdge (int nr, TORDER order) { order_edge[nr] = order; }
    INLINE void SetOrderFace (int nr, IVec<2,TORDER> order) { order_face[nr] = order; }

    INLINE void SetUseGradEdge(int nr, bool uge) { usegrad_edge[nr] = uge; }
    INLINE void SetUseGradFace(int nr, bool ugf) { usegrad_face[nr] = ugf; }


    INLINE void SetOrderCell (IVec<3> oi) { order_cell = oi; }

    /// set isotropic or anisotropic face orders
    template <typename TA>
    INLINE void SetOrderFace (const TA & of)
    { for (int i = 0; i < N_FACE; i++) order_face[i] = of[i]; }

    /// set edge orders
    template <typename TA>
    INLINE void SetOrderEdge (const TA & oe)
    { for (int i = 0; i < N_EDGE; i++) order_edge[i] = oe[i]; }

    /// use edge-gradients
    template <typename TA>
    INLINE void SetUseGradEdge (const TA & uge)
    { for (int i = 0; i < N_EDGE; i++) usegrad_edge[i] = uge[i]; }

    /// use face-gradients
    template <typename TA>
    INLINE void SetUseGradFace (const TA & ugf)
    { for (int i = 0; i < N_FACE; i++) usegrad_face[i] = ugf[i]; }

    INLINE void SetUseGradCell (bool ugc) 
    { usegrad_cell = ugc; }

    INLINE void SetType1 (bool t1)
    { type1 = t1; }
    
    void ComputeNDof();

    virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const override;
    virtual void CalcDualShape (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> shape) const override;
    virtual void EvaluateDual (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const override;
    virtual void AddDualTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
                               BareSliceVector<double> coefs) const override;
  };


  

  extern template class HCurlHighOrderFE<ET_SEGM>;
  extern template class HCurlHighOrderFE<ET_TRIG>;
  extern template class HCurlHighOrderFE<ET_QUAD>;
  extern template class HCurlHighOrderFE<ET_TET>;
  extern template class HCurlHighOrderFE<ET_PRISM>;
  extern template class HCurlHighOrderFE<ET_PYRAMID>;
  extern template class HCurlHighOrderFE<ET_HEX>;

  extern template class 
  T_HCurlHighOrderFiniteElement<ET_SEGM, HCurlHighOrderFE_Shape<ET_SEGM>>;
  extern template class 
  T_HCurlHighOrderFiniteElement<ET_TRIG, HCurlHighOrderFE_Shape<ET_TRIG>>;
  extern template class 
  T_HCurlHighOrderFiniteElement<ET_QUAD, HCurlHighOrderFE_Shape<ET_QUAD>>;

  extern template class 
  T_HCurlHighOrderFiniteElement<ET_TET, HCurlHighOrderFE_Shape<ET_TET>>;
  extern template class 
  T_HCurlHighOrderFiniteElement<ET_PRISM, HCurlHighOrderFE_Shape<ET_PRISM>>;
  extern template class 
  T_HCurlHighOrderFiniteElement<ET_PYRAMID, HCurlHighOrderFE_Shape<ET_PYRAMID>>;
  extern template class 
  T_HCurlHighOrderFiniteElement<ET_HEX, HCurlHighOrderFE_Shape<ET_HEX>>;
  
}


#endif

