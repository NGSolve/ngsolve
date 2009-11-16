#ifndef FILE_FACETHOFE
#define FILE_FACETHOFE

/*********************************************************************/
/* File:   facethofe.hpp                                             */
/* Author: A. Sinwel, H. Egger, J. Schoeberl                         */
/* Date:   2008                                                      */
/*********************************************************************/


#include "l2hofe.hpp"

namespace ngfem
{

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
