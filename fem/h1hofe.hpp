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

  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;    
    using ET_trait<ET>::ElementType;
    
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

    HD virtual tuple<int,int,int,int> GetNDofVEFC () const override
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
        || (ET == ET_HEX) || (ET == ET_TET) || (ET == ET_POINT);
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
