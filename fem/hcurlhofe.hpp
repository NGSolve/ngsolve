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
#include "dualinverse.hpp"

namespace ngfem
{

  /** 
      HCurlHighOrderFE of shape ET.
      The template specialization provides the shape functions.
  */
  template <ELEMENT_TYPE ET> class HCurlHighOrderFE_Shape;

  /*
    HCurl dofs per node, nominal order p on every node.
    Gradients come from H1(p+1-type1), the rest (curl part) is the
    same for type1 and default, and its curl lies in HDiv(p-1).
    Within each node the gradients come first.
   */
  namespace hcurl_ndof
  {
    inline int TrigGrad (int p, bool t1) { int g = p-t1; return g > 1 ? g*(g-1)/2 : 0; }
    inline int TrigCurl (int p) { return p > 1 ? (p+2)*(p-1)/2 : 0; }
    inline int QuadGrad (int p0, int p1, bool t1) { return max(p0-t1,0) * max(p1-t1,0); }
    inline int QuadCurl (int p0, int p1) { return (p0 > 0 && p1 > 0) ? p0*p1-1 : 0; }
    inline int TetGrad (int p, bool t1) { int g = p-t1; return g > 2 ? (g-2)*(g-1)*g/6 : 0; }
    inline int TetCurl (int p) { return p > 2 ? (p-2)*(p-1)*(2*p+3)/6 : 0; }
    inline int PrismGrad (int p, int pz, bool t1)
    { int g = p-t1, gz = pz-t1; return (g > 1 && gz > 0) ? g*(g-1)/2*gz : 0; }
    inline int PrismCurl (int p, int pz)
    { return (p > 1 && pz > 0) ? (pz-1)*(p-1)*p + (p-1)*(p-2)/2 : 0; }
    inline int HexGrad (int p0, int p1, int p2, bool t1)
    { return max(p0-t1,0) * max(p1-t1,0) * max(p2-t1,0); }
    inline int HexCurl (int p0, int p1, int p2)
    {
      if (p0 < 1 || p1 < 1 || p2 < 1) return 0;
      int q0 = p0-1, q1 = p1-1, q2 = p2-1;
      return 2*q0*q1*q2 + q0*q1 + q0*q2 + q1*q2;
    }
    inline int PyramidGrad (int p, bool t1) { int g = p-t1; return g > 1 ? (g-1)*g*(2*g-1)/6 : 0; }
    inline int PyramidCurl (int p) { return p > 1 ? p*(2*p*p+3*p-2)/3 : 0; }

    template <typename TP>
    inline int Face (ELEMENT_TYPE ft, TP p, bool grad, bool t1)
    {
      if (ft == ET_TRIG) return grad*TrigGrad(p[0], t1) + TrigCurl(p[0]);
      return grad*QuadGrad(p[0], p[1], t1) + QuadCurl(p[0], p[1]);
    }
    template <typename TP>
    inline int FaceGrad (ELEMENT_TYPE ft, TP p, bool grad, bool t1)
    {
      if (!grad) return 0;
      return (ft == ET_TRIG) ? TrigGrad(p[0], t1) : QuadGrad(p[0], p[1], t1);
    }
    // cell dofs, (p0,p1,p2) = order_cell; trig/quad: the 2D cell
    template <typename TP>
    inline int Cell (ELEMENT_TYPE et, TP p, bool grad, bool t1)
    {
      switch (et)
        {
        case ET_TRIG: return grad*TrigGrad(p[0], t1) + TrigCurl(p[0]);
        case ET_QUAD: return grad*QuadGrad(p[0], p[1], t1) + QuadCurl(p[0], p[1]);
        case ET_TET: return grad*TetGrad(p[0], t1) + TetCurl(p[0]);
        case ET_PRISM: return grad*PrismGrad(p[0], p[2], t1) + PrismCurl(p[0], p[2]);
        case ET_HEX: return grad*HexGrad(p[0], p[1], p[2], t1) + HexCurl(p[0], p[1], p[2]);
        case ET_PYRAMID: return grad*PyramidGrad(p[0], t1) + PyramidCurl(p[0]);
        default: return 0;
        }
    }
    template <typename TP>
    inline int CellGrad (ELEMENT_TYPE et, TP p, bool grad, bool t1)
    {
      if (!grad) return 0;
      switch (et)
        {
        case ET_TRIG: return TrigGrad(p[0], t1);
        case ET_QUAD: return QuadGrad(p[0], p[1], t1);
        case ET_TET: return TetGrad(p[0], t1);
        case ET_PRISM: return PrismGrad(p[0], p[2], t1);
        case ET_HEX: return HexGrad(p[0], p[1], p[2], t1);
        case ET_PYRAMID: return PyramidGrad(p[0], t1);
        default: return 0;
        }
    }
  }


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
    bool biorth = false;  // trig cell: biorthogonal basis (Haubold-Beuchler-Schoeberl 2024)

    // lazily filled D^-1 for full-order elements, owned by the space:
    // a table with one slot per VertexOrderClass
    std::atomic<double*> * dual_inverse = nullptr;

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
    HCurlHighOrderFE * SetVertexNumbers (FlatArray<int> vnums) override
    { VertexOrientedFE<ELEMENT_TYPE(ET)>::SetVertexNumbers(vnums); return this; }
    /* INLINE */ HCurlHighOrderFE () { ; } 

    /* INLINE */ HCurlHighOrderFE (int aorder) 
    {
      order_edge = aorder;
      order_face = IVec<2, TORDER>(aorder);
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
    INLINE void SetBiorth (bool b)
    { biorth = b; }

    void ComputeNDof();

    // same counts as ComputeNDof
    virtual tuple<int,int,int,int> GetNDofVEFC () const override
    {
      int ne = 0, nf = 0;
      for (int i = 0; i < N_EDGE; i++)
        ne += 1 + ((order_edge[i] > 0) ? usegrad_edge[i]*order_edge[i] : 0);
      for (int i = 0; i < N_FACE; i++)
        nf += hcurl_ndof::Face (FaceType(i), order_face[i], usegrad_face[i], type1);
      return { 0, ne, nf, ndof-ne-nf };
    }

    virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const override;
    virtual void CalcDualShape (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> shape) const override;
    virtual void EvaluateDual (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const override;
    virtual void AddDualTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
                               BareSliceVector<double> coefs) const override;

    void SetDualInverse (std::atomic<double*> * table) { dual_inverse = table; }

    // the class is taken from the vertex numbers the element carries now, since a
    // wrapper space may renumber them after GetFE (PeriodicFESpace does)
    std::atomic<double*> * DualInverseSlot () const
    {
      int cls = this->VertexOrderClass();
      return (dual_inverse && cls >= 0) ? dual_inverse + cls : nullptr;
    }

    using BASE::SolveDuality;
    bool SolveDuality (SliceMatrix<> x, LocalHeap & lh) const override
    {
      auto slot = DualInverseSlot();
      if (!slot) return false;
      FlatMatrix<> inv(ndof, ndof, GetDualInverse<ET,DIM> (*this, *slot, lh));
      HeapReset hr(lh);
      FlatMatrix<> tmp(x.Height(), x.Width(), lh);
      tmp = x;
      if (x.Width() == 1)  // gemv instead of n x 1 gemm
        x.Col(0) = inv * tmp.Col(0);
      else
        x = inv * tmp;
      return true;
    }

    bool SolveDualityTrans (SliceMatrix<> x, LocalHeap & lh) const override
    {
      auto slot = DualInverseSlot();
      if (!slot) return false;
      FlatMatrix<> inv(ndof, ndof, GetDualInverse<ET,DIM> (*this, *slot, lh));
      HeapReset hr(lh);
      FlatMatrix<> tmp(x.Height(), x.Width(), lh);
      tmp = x;
      if (x.Width() == 1)
        x.Col(0) = Trans(inv) * tmp.Col(0);
      else
        x = Trans(inv) * tmp;
      return true;
    }
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

