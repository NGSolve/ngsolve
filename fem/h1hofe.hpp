#ifndef FILE_H1HOFE
#define FILE_H1HOFE

/*********************************************************************/
/* File:   h1hofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


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
            template <ELEMENT_TYPE ET2> class TSHAPES = H1HighOrderFE_Shape,
            class BASE = T_ScalarFiniteElement< TSHAPES<ET>, ET> >
            
  class H1HighOrderFE : public BASE, 
                        public ET_trait<ET>
  {
  protected:

    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    // using ScalarFiniteElement<DIM>::eltype;


    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::N_EDGE;
    using ET_trait<ET>::N_FACE;
    using ET_trait<ET>::N_CELL;
    using ET_trait<ET>::FaceType;
    using ET_trait<ET>::GetEdgeSort;
    using ET_trait<ET>::GetFaceSort;
    using ET_trait<ET>::PolDimension;
    using ET_trait<ET>::PolBubbleDimension;

    /// global vertex numbers used of edge/face orientation
    Vec<N_VERTEX, int> vnums;

    /// order of edge shapes
    Vec<N_EDGE, int> order_edge; 

    /// order of face shapes
    Vec<N_FACE, INT<2> > order_face; 

    /// order of internal shapes (3d only)
    Vec<N_CELL, INT<3> > order_cell;

  public:
    /// minimal constructor, orders will be set later
    H1HighOrderFE () { ; }

    /// builds a functional element of order aorder.
    H1HighOrderFE (int aorder)
    { 
      ndof = PolDimension (aorder);
      
      for (int i = 0; i < N_VERTEX; i++) vnums[i] = i;
      for (int i = 0; i < N_EDGE; i++) order_edge[i] = aorder;
      for (int i = 0; i < N_FACE; i++) order_face[i] = aorder;   
      if (DIM == 3) order_cell[0] = aorder; 
      
      order = aorder;
    }

    /// assignes vertex numbers
    template <typename TA> 
    void SetVertexNumbers (const TA & avnums)
    { for (int i = 0; i < N_VERTEX; i++) vnums[i] = avnums[i]; }

    /// assign vertex number
    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }

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
    void SetOrderFace (int nr, INT<2> order) { order_face[nr] = order; }


    /// set anisotropic cell order
    void SetOrderCell (INT<3> oi)  { order_cell[0] = oi; }

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
      
      order = 1;
      for (int i = 0; i < N_EDGE; i++) order = max(order, order_edge[i]);
      for (int i = 0; i < N_FACE; i++) order = max(order, Max (order_face[i])); 
      if (DIM == 3) order = max (order, Max (order_cell[0]));
    }


  };



}


#endif
