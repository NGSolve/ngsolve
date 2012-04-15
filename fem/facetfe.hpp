#ifndef FILE_FACETFE
#define FILE_FACETFE

/*********************************************************************/
/* File:   facetfe.hpp                                               */
/* Author: A. Sinwel, H. Egger, J. Schoeberl                         */
/* Date:   2008                                                      */
/*********************************************************************/

namespace ngfem
{

  /*
   * Facet Finite Elements
   */ 

  template <int D>
  class NGS_DLL_HEADER FacetVolumeFiniteElement;
  
  
  template <int D>
  class FacetFEFacet : virtual public ScalarFiniteElement<D>
  {
    int fnr;
    const FacetVolumeFiniteElement<D> & fe;
  public:
    FacetFEFacet (int afnr,
		  const FacetVolumeFiniteElement<D> & afe,
		  ELEMENT_TYPE aeltype, int andof, int aorder)
      : ScalarFiniteElement<D> (aeltype, andof, aorder), fnr(afnr), fe(afe) { ; }

    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const
    {
      fe.CalcFacetShapeVolIP(fnr, ip, shape);
    }
  };

  

  template <int D>
  class NGS_DLL_HEADER FacetVolumeFiniteElement : public FiniteElement
  {
  protected:
    int vnums[8];
    int facet_order[6]; 
    int first_facet_dof[7];

    // bool highest_order_dc;

    using FiniteElement::eltype;
    using FiniteElement::order;

  public:
    FacetVolumeFiniteElement (ELEMENT_TYPE aeltype)
      : FiniteElement (D, aeltype,-1,-1)
    {
      // highest_order_dc = false;
      for (int i=0; i<8; i++) 
	vnums[i] = i;
      for (int i = 0; i < 6; i++) 
	facet_order[i] = -1;
      for (int i=0; i < 7; i++) 
	first_facet_dof[i] = 0;
    }


    void SetVertexNumbers (FlatArray<int> & avnums)
    {
      for (int i = 0; i < avnums.Size(); i++)
	vnums[i] = avnums[i];
    }

    // void SetHighestOrderDC(bool set){highest_order_dc=set;}

    void SetOrder (int ao)  
    {
      order = ao;
      for (int i = 0; i < 6; i++)
	facet_order[i] = ao;
    }
    
    void SetOrder (FlatArray<int> & ao)
    {
      for (int i=0; i<ao.Size(); i++)
	facet_order[i] = ao[i];
      
      order = facet_order[0];        // integration order
      for (int i = 1; i < ao.Size(); i++)
	order = max(order, ao[i]);
    }
    
    
    FacetFEFacet<D> Facet (int fnr) const { return FacetFEFacet<D> (fnr, *this, eltype, ndof, facet_order[fnr]); }
    virtual void CalcFacetShapeVolIP (int fnr, const IntegrationPoint & ip, FlatVector<> shape) const = 0;


    /*
      // please use Facet(k).CalcShape (ip, shape)   instead
      // ip .. vol point

    void CalcFacetShape(int fnr, const IntegrationPoint & ip, FlatVector<> shape) const;
    void EvaluateFacet (int fnr, const IntegrationRule & ir, FlatVector<> coefs, FlatVector<> values) const;
    void EvaluateFacetTrans (int fnr, const IntegrationRule & ir, FlatVector<> values, FlatVector<> coefs) const;
    void SelectFacet (int afn) const { facetnr = afn; }
    */


    IntRange GetFacetDofs(int fnr) const
    {
      return IntRange (first_facet_dof[fnr], first_facet_dof[fnr+1]);
    }

    void GetFacetDofNrs(int fnr, Array<int>& fdnums) const
    {
      fdnums = GetFacetDofs(fnr);
    }

    int GetFacetNDof(int afnr) const { return first_facet_dof[afnr+1] - first_facet_dof[afnr]; };
    int GetFirstFacetDof(int afnr) const { return first_facet_dof[afnr]; } 
  
    virtual void ComputeNDof () = 0;
  };

}



#endif
