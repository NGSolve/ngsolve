#ifndef FILE_H1HOFE
#define FILE_H1HOFE

/*********************************************************************/
/* File:   h1hofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

#include "tscalarfe.hpp"

namespace ngfem
{

  
  /// default shape function engine for high order h1-elements
  template <ELEMENT_TYPE ET> class H1HighOrderFE_Shape;


  /**
     High order finite elements for H1.  These are the actual finite
     element classes to be used.  
     shape functions are provided by the shape template
   */

  template <ELEMENT_TYPE ET, 
            class SHAPES = H1HighOrderFE_Shape<ET>,
            class BASE = T_ScalarFiniteElement< SHAPES, ET> >
  

  class H1HighOrderFE : public BASE, public ET_trait<ET>, public VertexOrientedFE<ET>
  {
  protected:

    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;

    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::N_EDGE;
    using ET_trait<ET>::N_FACE;
    using ET_trait<ET>::N_CELL;
    using ET_trait<ET>::FaceType;
    using ET_trait<ET>::GetEdgeSort;
    using ET_trait<ET>::GetFaceSort;
    using ET_trait<ET>::PolDimension;
    using ET_trait<ET>::PolBubbleDimension;

    typedef unsigned char TORDER;

    /// order of edge shapes
    std::array<TORDER, N_EDGE> order_edge; 

    /// order of face shapes
    std::array<IVec<2,TORDER>, N_FACE> order_face; 

    /// order of internal shapes (3d only)
    std::array<IVec<3,TORDER>, N_CELL > order_cell;
    
    bool nodalp2 = false;

    // lazily filled D^-1 for full-order elements, owned by the space
    std::atomic<double*> * dual_inverse = nullptr;

  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    H1HighOrderFE * SetVertexNumbers (FlatArray<int> vnums) override
    { VertexOrientedFE<ET>::SetVertexNumbers(vnums); return this; }

    // using ET_trait<ET>::ElementType;
    using BASE::ElementType;
    
    INLINE void SetNodalp2() { nodalp2 = true; }
    
    /// minimal constructor, orders will be set later
    INLINE H1HighOrderFE () { ; } 

    /// builds a functional element of order aorder.
    INLINE H1HighOrderFE (int aorder)
    { 
      ndof = PolDimension (aorder);
      for (int i = 0; i < N_VERTEX; i++) this->SetVertexNumber(i,i);
      for (int i = 0; i < N_EDGE; i++) order_edge[i] = aorder;
      for (int i = 0; i < N_FACE; i++) order_face[i] = aorder;   
      if (DIM == 3) order_cell[0] = aorder; 
      
      order = aorder;
    }

    // virtual NGS_DLL_HEADER ~H1HighOrderFE () { ; }
    /// set edge orders
    template <typename TA>
    void SetOrderEdge (const TA & oe)
    { for (int i = 0; i < N_EDGE; i++) order_edge[i] = oe[i]; }

    /// set edge order for edge nr
    void SetOrderEdge (int nr, int order) { order_edge[nr] = order; }

    /// set isotropic or anisotropic face orders
    template <typename TA>
    void SetOrderFace (const TA & of)
    { for (int i = 0; i < N_FACE; i++) order_face[i] = of[i]; }

    /// set anisotropic face order for face nr
    void SetOrderFace (int nr, IVec<2> order) { order_face[nr] = order; }

    /// set anisotropic cell order
    void SetOrderCell (IVec<3> oi)  { order_cell[0] = oi; }

    /// compute the element space dimension
    void ComputeNDof()
    {
      ndof = N_VERTEX;
      for (int i = 0; i < N_EDGE; i++)
        ndof += order_edge[i] - 1;
      for (int i = 0; i < N_FACE; i++)
        ndof += ::ngfem::PolBubbleDimension (FaceType(i), order_face[i]);
      if (DIM == 3)
        ndof += PolBubbleDimension (order_cell[0]);

      TORDER ho = 1;
      for (int i = 0; i < N_EDGE; i++) ho = max2(ho, order_edge[i]);
      for (int i = 0; i < N_FACE; i++) ho = max2(ho, Max (order_face[i])); 
      if (DIM == 3) order = max2 (ho, Max (order_cell[0]));
      order = ho;
    }

    virtual tuple<int,int,int,int> GetNDofVEFC () const override
    {
      int nv = N_VERTEX;
      int ne = 0, nf = 0, nc = 0;
      
      for (int i = 0; i < N_EDGE; i++)
        ne += order_edge[i] - 1;
      
      for (int i = 0; i < N_FACE; i++)
        nf += ::ngfem::PolBubbleDimension (FaceType(i), order_face[i]);
      
      if (DIM == 3)
        nc += PolBubbleDimension (order_cell[0]);
      return { nv, ne, nf, nc };
    }
      
    
    virtual bool DualityMassDiagonal () const override
    {
      return (ET == ET_SEGM) || (ET == ET_TRIG) || (ET == ET_QUAD)
        || (ET == ET_HEX) || (ET == ET_TET) || (ET == ET_PYRAMID) || (ET == ET_PRISM) || (ET == ET_POINT);
    }

    void SetDualInverse (std::atomic<double*> * slot) { dual_inverse = slot; }

    using BASE::SolveDuality;
    bool SolveDuality (SliceMatrix<> x, LocalHeap & lh) const override
    {
      if (!dual_inverse)
        return BASE::SolveDuality (x, lh);

      // D^-1 is block lower triangular in V,E,F,C order, with diagonal diagonal blocks
      FlatMatrix<> inv(ndof, ndof, GetDualInverse(lh));
      auto [nv, ne, nf, nc] = GetNDofVEFC();
      int st[] = { 0, nv, nv+ne, nv+ne+nf, int(ndof) };
      if (x.Width() == 1)  // gemv instead of n x 1 gemm
        {
          auto xv = x.Col(0);
          for (int b = 3; b >= 0; b--)
            {
              int first = st[b], last = st[b+1];
              for (int i = first; i < last; i++)
                xv(i) *= inv(i,i);
              if (first > 0 && last > first)
                xv.Range(first, last) += inv.Rows(first, last).Cols(0, first) * xv.Range(0, first);
            }
          return true;
        }

      for (int b = 3; b >= 0; b--)
        {
          int first = st[b], last = st[b+1];
          for (int i = first; i < last; i++)
            x.Row(i) *= inv(i,i);
          if (first > 0 && last > first)
            x.Rows(first, last) += inv.Rows(first, last).Cols(0, first) * x.Rows(0, first);
        }
      return true;
    }

    bool SolveDualityTrans (SliceMatrix<> x, LocalHeap & lh) const override
    {
      if (!dual_inverse) return false;
      FlatMatrix<> inv(ndof, ndof, GetDualInverse(lh));
      auto [nv, ne, nf, nc] = GetNDofVEFC();
      int st[] = { 0, nv, nv+ne, nv+ne+nf, int(ndof) };
      if (x.Width() == 1)
        {
          auto xv = x.Col(0);
          for (int b = 0; b < 4; b++)
            {
              int first = st[b], last = st[b+1];
              for (int i = first; i < last; i++)
                xv(i) *= inv(i,i);
              if (last < ndof && last > first)
                xv.Range(first, last) += Trans(inv.Rows(last, ndof).Cols(first, last)) * xv.Range(last, ndof);
            }
          return true;
        }

      for (int b = 0; b < 4; b++)
        {
          int first = st[b], last = st[b+1];
          for (int i = first; i < last; i++)
            x.Row(i) *= inv(i,i);
          if (last < ndof && last > first)
            x.Rows(first, last) += Trans(inv.Rows(last, ndof).Cols(first, last)) * x.Rows(last, ndof);
        }
      return true;
    }

  private:
    double * GetDualInverse (LocalHeap & lh) const
    {
      double * pinv = dual_inverse->load (std::memory_order_acquire);
      if (!pinv)
        {
          // columns of D^-1 from the matrix-free solve
          pinv = new double[ndof*ndof];
          FlatMatrix<> inv(ndof, ndof, pinv);
          Vector<> ej(ndof);
          for (int j = 0; j < ndof; j++)
            {
              HeapReset hr(lh);
              ej = 0.0; ej(j) = 1.0;
              this->SolveDuality1 (ej, inv.Col(j), lh);
            }
          double * expected = nullptr;
          if (!dual_inverse->compare_exchange_strong (expected, pinv, std::memory_order_acq_rel))
            {
              delete [] pinv;
              pinv = expected;
            }
        }
      return pinv;
    }
  };

}  



#ifdef FILE_H1HOFE_CPP

#define H1HOFE_EXTERN
#include <h1hofe_impl.hpp>
#include <tscalarfe_impl.hpp>

#else

#define H1HOFE_EXTERN extern

#endif

namespace ngfem
{
  H1HOFE_EXTERN template class H1HighOrderFE<ET_POINT>;
  extern template class H1HighOrderFE<ET_SEGM>;
  extern template class H1HighOrderFE<ET_TRIG>;
  extern template class H1HighOrderFE<ET_QUAD>;

  extern template class H1HighOrderFE<ET_TET>;
  extern template class H1HighOrderFE<ET_PRISM>;
  extern template class H1HighOrderFE<ET_PYRAMID>;
  H1HOFE_EXTERN template class H1HighOrderFE<ET_HEXAMID>;
  extern template class H1HighOrderFE<ET_HEX>;

  H1HOFE_EXTERN template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_POINT>, ET_POINT>;
  extern template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_SEGM>, ET_SEGM>;
  extern template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_TRIG>, ET_TRIG>;
  extern template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_QUAD>, ET_QUAD>;

  extern template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_TET>, ET_TET>;
  extern template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_PRISM>, ET_PRISM>;
  extern template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_PYRAMID>, ET_PYRAMID>;
  H1HOFE_EXTERN template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_HEXAMID>, ET_HEXAMID>;
  extern template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_HEX>, ET_HEX>;
}


#endif
