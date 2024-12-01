#ifndef FILE_NODALHOFE
#define FILE_NODALHOFE

#include "tscalarfe.hpp"

namespace ngfem
{



  template <ELEMENT_TYPE ET>
  class NodalHOFE : public  T_ScalarFiniteElement<NodalHOFE<ET>, ET>,
                    public ET_trait<ET>, public VertexOrientedFE<ET>
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

    
  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;    
    NodalHOFE * SetVertexNumbers (FlatArray<int> vnums) override
    { VertexOrientedFE<ELEMENT_TYPE(ET)>::SetVertexNumbers(vnums); return this; }   // cast for msvc ?
    using ET_trait<ET>::ElementType;

    /// builds a functional element of order aorder.
    INLINE NodalHOFE (int aorder)
    { 
      ndof = PolDimension (aorder);
      order = aorder;
    }

    
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (TIP<DIM,Tx> ip, TFA & shape) const;
    
    /*
    virtual tuple<int,int,int,int> GetNDofVEFC () const override
    {
      int nv = N_VERTEX;
      int ne = 0, nf = 0, nc = 0;
      
      for (int i = 0; i < N_EDGE; i++)
        ne += int(order) - 1;
      
      for (int i = 0; i < N_FACE; i++)
        nf += ::ngfem::PolBubbleDimension (FaceType(i), INT<2>(order));
      
      if (DIM == 3)
        nc += PolBubbleDimension (INT<3>(order));
      return { nv, ne, nf, nc };
    }
    */
  };

}  



#ifdef FILE_NODALHOFE_CPP

#define NODALHOFE_EXTERN
#include <nodalhofe_impl.hpp>
#include <tscalarfe_impl.hpp>

#else

#define NODALHOFE_EXTERN extern

#endif

namespace ngfem
{
  NODALHOFE_EXTERN template class NodalHOFE<ET_SEGM>;
  NODALHOFE_EXTERN template class NodalHOFE<ET_TRIG>;
  NODALHOFE_EXTERN template class NodalHOFE<ET_TET>;

  NODALHOFE_EXTERN template class T_ScalarFiniteElement<NodalHOFE<ET_SEGM>, ET_SEGM>;
  NODALHOFE_EXTERN template class T_ScalarFiniteElement<NodalHOFE<ET_TRIG>, ET_TRIG>;
  NODALHOFE_EXTERN template class T_ScalarFiniteElement<NodalHOFE<ET_TET>, ET_TET>;
}

#endif
