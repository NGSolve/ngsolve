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
  class FacetVolumeFiniteElement : public FiniteElement
  {
  protected:
    int vnums[8];
    int facet_order[6]; 
    int first_facet_dof[7];

    ScalarFiniteElement<D-1> * facets[6];
  public:
    FacetVolumeFiniteElement (ELEMENT_TYPE aeltype);

    void SetVertexNumbers (FlatArray<int> & avnums);
    void SetOrder (int ao);
    void SetOrder (FlatArray<int> & ao);
    int GetFacetOrder (int j) const { return facet_order[j]; }
    int GetVertexNumber (int j) const { return vnums[j]; }
   
    void CalcFacetShape(int fnr, const IntegrationPoint & ip, FlatVector<> shape) const;
    void EvaluateFacet (int fnr, const IntegrationRule & ir, FlatVector<> coefs, FlatVector<> values) const;
    void EvaluateFacetTrans (int fnr, const IntegrationRule & ir, FlatVector<> values, FlatVector<> coefs) const;
    
    


    void GetFacetDofNrs(int afnr, Array<int>& fdnums) const; 
    int GetFacetNDof(int afnr) const { return first_facet_dof[afnr+1] - first_facet_dof[afnr]; };
    int GetFirstFacetDof(int afnr) const { return first_facet_dof[afnr]; } 
  
    const ScalarFiniteElement<D-1> & GetFacetFE(int fnr) const { return *facets[fnr]; }

    virtual void ComputeNDof () = 0;
  };
}



#endif
