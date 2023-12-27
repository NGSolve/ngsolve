#ifndef FILE_FACETHOFE
#define FILE_FACETHOFE

/*********************************************************************/
/* File:   facethofe.hpp                                             */
/* Author: A. Sinwel, H. Egger, J. Schoeberl                         */
/* Date:   2008                                                      */
/*********************************************************************/


#include "facetfe.hpp"


namespace ngfem
{
  
  template <ELEMENT_TYPE ET>
  class FacetFE : 
    public FacetVolumeFiniteElement<ET_trait<ET>::DIM>,
    public ET_trait<ET>,
    public VertexOrientedFE<ET>
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };
    using ET_trait<ET>::N_FACET;

    using FacetVolumeFiniteElement<ET_trait<ET>::DIM>::order;
    using FacetVolumeFiniteElement<ET_trait<ET>::DIM>::ndof;
    using FacetVolumeFiniteElement<ET_trait<ET>::DIM>::first_facet_dof;
    using FacetVolumeFiniteElement<ET_trait<ET>::DIM>::facet_order;
    using VertexOrientedFE<ET>::vnums;
    using VertexOrientedFE<ET>::GetVertexOrientedEdge;
    using VertexOrientedFE<ET>::GetVertexOrientedFace;
    bool nodal;
  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
  public:
    FacetFE (bool anodal = false) : nodal(anodal) { ; }

    virtual ELEMENT_TYPE ElementType() const override { return ET; }
    
    using FacetVolumeFiniteElement<ET_trait<ET>::DIM>::SetOrder;
    template <typename TA>
    void SetOrder (const TA & of)
    {
      for (int i = 0; i < N_FACET; i++) facet_order[i] = of[i];

      order = facet_order[0];      
      for (int i = 1; i < N_FACET; i++)
	order = max2(order, of[i]);
    }

    
    virtual void ComputeNDof() override
    {
      ndof = 0;
      for (int i = 0; i < N_FACET; i++)
	{
	  first_facet_dof[i] = ndof;
	  int fo = facet_order[i];
	  switch (ElementTopology::GetFacetType (ET, i))
	    {
	    case ET_POINT: ndof += 1; break;
	    case ET_SEGM: ndof += fo+1; break;
	    case ET_TRIG: ndof += ( (fo+1) * (fo+2) ) / 2; break;
	    case ET_QUAD: ndof += sqr (fo+1); break;
	    default: ;
	    }
	}      
      first_facet_dof[N_FACET] = ndof;
      order = facet_order[0];      
      for (int i = 1; i < N_FACET; i++)
	order = max2(order, facet_order[i]);
    }
    

    virtual void CalcFacetShapeVolIP (int fnr, const IntegrationPoint & ip,
                                      BareSliceVector<> shape) const override;

    virtual void CalcFacetShapeVolIR (int fnr, const SIMD_IntegrationRule & ir, 
                                      BareSliceMatrix<SIMD<double>> shape) const override;
    
    virtual void EvaluateFacetVolIp (int fnr, const SIMD_IntegrationRule & ir,
                                     BareSliceVector<> coefs, BareVector<SIMD<double>> values) const override;
    
    virtual void AddTransFacetVolIp (int fnr, const SIMD_IntegrationRule & ir,
                                     BareVector<SIMD<double>> values, BareSliceVector<> coefs) const override;
    
    virtual void CalcFacetDShapeVolIP (int fnr, const IntegrationPoint & ip, 
                                       BareSliceMatrix<> shape) const override;

  private:
    template<typename Tx, typename TFA>  
    //void T_CalcShapeFNr (int fnr, Tx x[ET_trait<ET>::DIM], TFA & shape) const;
    void T_CalcShapeFNr (int fnr, TIP<ET_trait<ET>::DIM,Tx> ip, TFA && shape) const;
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
#endif

}

#endif
