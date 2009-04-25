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

    L2HighOrderFiniteElement<D-1> * facets[6];
  public:
    FacetVolumeFiniteElement (ELEMENT_TYPE aeltype);

    void SetVertexNumbers (FlatArray<int> & avnums);
    void SetOrder (int ao);
    void SetOrder (FlatArray<int> & ao);
    int GetFacetOrder (int j) const { return facet_order[j]; }
    int GetVertexNumber (int j) const { return vnums[j]; }
   
    void CalcFacetShape(int fnr, const IntegrationPoint & ip, FlatVector<> shape) const;
  
    void GetFacetDofNrs(int afnr, Array<int>& fdnums) const; 
    int GetFacetNDof(int afnr) const { return first_facet_dof[afnr+1] - first_facet_dof[afnr]; };
    int GetFirstFacetDof(int afnr) const { return first_facet_dof[afnr]; } 
  
    const ScalarFiniteElement<D-1> & GetFacetFE(int fnr) const { return *facets[fnr]; }

    virtual void ComputeNDof () = 0;
  };


  //------------------------------------------------------------
  class FacetVolumeTrig : public FacetVolumeFiniteElement<2>
  {
  protected:
    L2HighOrderFE<ET_SEGM> facets2[3];
  public:
    FacetVolumeTrig(); 
    virtual void ComputeNDof();
  };


  // --------------------------------------------------------
  class FacetVolumeQuad : public FacetVolumeFiniteElement<2>
  {
  protected:
    L2HighOrderFE<ET_SEGM> facets2[4];
  public:
    FacetVolumeQuad();
    virtual void ComputeNDof();
  };

  // --------------------------------------------------------
  class FacetVolumeTet : public FacetVolumeFiniteElement<3>
  {
  protected:
    L2HighOrderFE<ET_TRIG> facets2[4];
  public:
    FacetVolumeTet(); 
    virtual void ComputeNDof();
  };

  // --------------------------------------------------------
  class FacetVolumeHex : public FacetVolumeFiniteElement<3>
  {
  protected:
    L2HighOrderFE<ET_QUAD> facetsq[6];
  public:
    FacetVolumeHex();
    virtual void ComputeNDof();
  };

  // --------------------------------------------------------
  class FacetVolumePrism : public FacetVolumeFiniteElement<3>
  {
  protected:
    L2HighOrderFE<ET_TRIG> facetst[2];
    L2HighOrderFE<ET_QUAD> facetsq[3];

  public:
    FacetVolumePrism();
    virtual void ComputeNDof();
  };


  // --------------------------------------------------------
  class FacetVolumePyramid : public FacetVolumeFiniteElement<3>
  {
  protected:
    int tnr, qnr; // active facets
    L2HighOrderFE<ET_TRIG> trig;
    L2HighOrderFE<ET_QUAD> quad;
  public:
    FacetVolumePyramid() : FacetVolumeFiniteElement<3> (ET_PYRAMID) {  qnr=tnr=-1; };
  
    // virtual void CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const; // maybe convenient for shape tester??
    virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const; 
    virtual void SetFacet(int afnr) const;
  
    virtual void ComputeNDof();
    //     virtual const FiniteElement & GetFacetFE(int fnr, LocalHeap& lh) const;
  };

}

#endif
