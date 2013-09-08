#ifndef FILE_FACETHOFE
#define FILE_FACETHOFE

/*********************************************************************/
/* File:   facethofe.hpp                                             */
/* Author: A. Sinwel, H. Egger, J. Schoeberl                         */
/* Date:   2008                                                      */
/*********************************************************************/

// #include "tscalarfe.cpp"

namespace ngfem
{
  
  template <ELEMENT_TYPE ET> class FacetFE;


  template <ELEMENT_TYPE ET> 
  class FacetFiniteElement_Family :
    public FacetVolumeFiniteElement<ET_trait<ET>::DIM>,
    public ET_trait<ET> 
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };
    using ET_trait<ET>::N_FACET;

    using FacetVolumeFiniteElement<ET_trait<ET>::DIM>::ndof;
    using FacetVolumeFiniteElement<ET_trait<ET>::DIM>::first_facet_dof;
    using FacetVolumeFiniteElement<ET_trait<ET>::DIM>::facet_order;

  public:
    FacetFiniteElement_Family () { ; }

    virtual ELEMENT_TYPE ElementType() const { return ET; }

    virtual void ComputeNDof() 
    {
      ndof = 0;
      for (int i = 0; i < N_FACET; i++)
	{
	  first_facet_dof[i] = ndof;
	  int fo = facet_order[i];
	  switch (ElementTopology::GetFacetType (ET, i))
	    {
	    case ET_SEGM: ndof += fo+1; break;
	    case ET_TRIG: ndof += ( (fo+1) * (fo+2) ) / 2; break;
	    case ET_QUAD: ndof += sqr (fo+1); break;
	    default: ;
	    }
	}
      first_facet_dof[N_FACET] = ndof;
    }
    

    virtual void CalcFacetShapeVolIP(int fnr, const IntegrationPoint & ip, 
				     SliceVector<> shape) const
    {
      double pt[DIM];
      for (int i = 0; i < DIM; i++) pt[i] = ip(i);
      static_cast<const FacetFE<ET>*>(this)->T_CalcShapeFNr (fnr, pt, shape); 
    }
  };

  template <ELEMENT_TYPE ET>
  class FacetFE : public FacetFiniteElement_Family<ET>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShapeFNr (int fnr, Tx x[ET_trait<ET>::DIM], TFA & shape) const;
  };


#ifdef FILE_FACETHOFE_CPP
#else
  extern template class FacetFE<ET_SEGM>;
  extern template class FacetFE<ET_TRIG>;
  extern template class FacetFE<ET_QUAD>;
  extern template class FacetFE<ET_TET>;
  extern template class FacetFE<ET_HEX>;
  extern template class FacetFE<ET_PRISM>;
  extern template class FacetFE<ET_PYRAMID>;

  extern template class FacetFiniteElement_Family<ET_SEGM>;
  extern template class FacetFiniteElement_Family<ET_TRIG>;
  extern template class FacetFiniteElement_Family<ET_QUAD>;
  extern template class FacetFiniteElement_Family<ET_TET>;
  extern template class FacetFiniteElement_Family<ET_HEX>;
  extern template class FacetFiniteElement_Family<ET_PRISM>;
  extern template class FacetFiniteElement_Family<ET_PYRAMID>;

#endif

}

#endif
