#ifndef VECTOR_FACET_FE_HPP__
#define VECTOR_FACET_FE_HPP__

/*********************************************************************/
/* File:   vectorfacetfe.hpp                                         */
/* Author: A. Sinwel, (J. Schoeberl)                                 */
/* Date:   2008                                                      */
/*********************************************************************/

#include <fem.hpp>

namespace ngfem 
{
  /*
    facet element with tangential facet components.
    to be changed similar to scalar facetfe
  */
  
  template <ELEMENT_TYPE ET>
  class VectorFacetFacetFE : public HCurlFiniteElement<ET_trait<ET>::DIM>, public VertexOrientedFE<ET>
  {
  protected:
    using VertexOrientedFE<ET>::GetVertexOrientedEdge;
    using VertexOrientedFE<ET>::GetVertexOrientedFace;
    using VertexOrientedFE<ET>::vnums;

    INT<2> order_inner;
    // using HCurlFiniteElement<D>::eltype;
    using HCurlFiniteElement<ET_trait<ET>::DIM>::order;
 
  public:
    using VertexOrientedFE<ET>::SetVertexNumber;
    using VertexOrientedFE<ET>::SetVertexNumbers;

    VectorFacetFacetFE (int aorder)
    {
      order = aorder;
      order_inner = INT<2>(aorder,aorder);
      ComputeNDof();
    }

    VectorFacetFacetFE ()
      : order_inner (0)
    {
      // for(int i=0; i<VertexOrientedFE<ET>::N_VERTEX; i++)
      //   SetVertexNumber(i,-1);
    }

    HD virtual ELEMENT_TYPE ElementType() const { return ET; }

    INLINE void SetOrder (int aorder)
    {
      order = aorder;
      order_inner = aorder;
      ComputeNDof();
    }
  
    INLINE void SetOrder (INT<2> oi)
    {
      order = max2 (oi[0], oi[1]);
      order_inner = oi;
      ComputeNDof();
    }

    virtual void ComputeNDof ();

    virtual void CalcShape(const IntegrationPoint & ip,
         		    SliceMatrix<> shape) const;

  };



  template <int D>
  class VectorFacetVolumeFiniteElement : public HCurlFiniteElement<D>
  {
  protected:
    int vnums[8];
    INT<2> facet_order[6]; 
    int first_facet_dof[7];
    bool highest_order_dc;

    ELEMENT_TYPE eltype;
    // using HCurlFiniteElement<D>::eltype;
    using HCurlFiniteElement<D>::order;
  
  public:
    VectorFacetVolumeFiniteElement () // : nv(0), nf(0)
    { highest_order_dc=false; }

    VectorFacetVolumeFiniteElement (ELEMENT_TYPE aeltype);
    HD virtual ELEMENT_TYPE ElementType() const { return eltype; }

    void SetHighestOrderDC(bool set){highest_order_dc=set;}
    void SetVertexNumbers (FlatArray<int> & avnums);

    void SetOrder(int ao);

    void SetOrder(FlatArray<int> & ao);
    void SetOrder(FlatArray<INT<2> > & ao);

    INT<2> GetFacetOrder(int j) const 
    { return facet_order[j]; }
    int GetVertexNumber(int j) const 
    { return vnums[j]; }
  

    //   virtual void TransformFacetToVolumeShape ( int fanr, FlatMatrix<> shape1d, 
    // 					     FlatMatrix<> shape ) const;

    virtual void CalcShape (const IntegrationPoint & ip, SliceMatrix<> shape) const;
    virtual void CalcShape (const IntegrationPoint & ip, int facet, SliceMatrix<> shape) const = 0;
    
    virtual int GetNExtraShapes( int facet) const {return 0;}
    virtual void CalcExtraShape (const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<D> xshape) const {xshape = 0.0;}

    virtual void GetFacetDofNrs(int afnr, Array<int>& fdnums) const; 

    virtual int GetFacetNDof(int afnr) const 
    { return first_facet_dof[afnr+1] - first_facet_dof[afnr]; };

    virtual int GetFirstFacetDof(int afnr) const { 
      return first_facet_dof[afnr];}; 
  
    virtual void ComputeNDof () = 0;

    /// degrees of freedom sitting inside the element, used for static condensation
    virtual void GetInternalDofs (Array<int> & idofs) const;
    
    // utility
    // virtual int GetNF() const { return nf; };
    // virtual int GetNV() const { return nv; };
    // virtual void GetVertexNumbers(Array<int>&) const;
    // virtual void GetFacetOrders(Array<INT<2> >&) const;
  };



  class VectorFacetVolumeTrig : public VectorFacetVolumeFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    using VertexOrientedFE<ET_TRIG>::vnums;
  public:
    VectorFacetVolumeTrig() : VectorFacetVolumeFiniteElement<2>(ET_TRIG) { ; }
    virtual void ComputeNDof();
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
   
    virtual void CalcShape ( const IntegrationPoint & ip, int facet, SliceMatrix<> shape) const;
    virtual int GetNExtraShapes( int facet) const {return 1;};
    virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<2> xshape) const;
  };



  class VectorFacetVolumeQuad : public VectorFacetVolumeFiniteElement<2>
  {
  public:
    VectorFacetVolumeQuad() : VectorFacetVolumeFiniteElement<2> (ET_QUAD) { ; }
    virtual void ComputeNDof();
   
    virtual void CalcShape (const IntegrationPoint & ip, int facet, SliceMatrix<> shape) const;
    virtual int GetNExtraShapes( int facet) const {return 1;};
    virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<2> xshape) const;
    
  };


  class VectorFacetVolumeTet : public VectorFacetVolumeFiniteElement<3>, public VertexOrientedFE<ET_TET>
  {
    using VertexOrientedFE<ET_TET>::vnums;
  public:
    VectorFacetVolumeTet() : VectorFacetVolumeFiniteElement<3> (ET_TET) { ; }
    virtual void ComputeNDof();
    using VertexOrientedFE<ET_TET>::SetVertexNumbers;

    virtual void CalcShape (const IntegrationPoint & ip, int facet, SliceMatrix<> shape) const;
    virtual int GetNExtraShapes( int facet) const {return 2*(facet_order[facet][0]+2);};
    virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> xshape) const;    
  };



  class VectorFacetVolumeHex : public VectorFacetVolumeFiniteElement<3>
  {
  public:
    VectorFacetVolumeHex() : VectorFacetVolumeFiniteElement<3>(ET_HEX) { ; };
    virtual void ComputeNDof();
    HD virtual ELEMENT_TYPE ElementType() const { return ET_HEX; }   
    virtual void CalcShape ( const IntegrationPoint & ip, int facet, SliceMatrix<> shape ) const;
//     virtual int GetNExtraShapes( int facet) const {return 2*(2*facet_order[facet][0]+3);};
//     virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> xshape) const;    
  };



  class VectorFacetVolumePrism : public VectorFacetVolumeFiniteElement<3>
  {
  public:
    VectorFacetVolumePrism() : VectorFacetVolumeFiniteElement<3> (ET_PRISM) { ; };
    virtual void ComputeNDof();
   
    virtual void CalcShape ( const IntegrationPoint & ip, int facet, SliceMatrix<> shape) const;
    virtual int GetNExtraShapes( int facet) const;
    virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> xshape) const;    
  };


  class VectorFacetVolumePyramid : public VectorFacetVolumeFiniteElement<3>
  {
  public:
    VectorFacetVolumePyramid() : VectorFacetVolumeFiniteElement<3> (ET_PYRAMID) { ; };
    virtual void ComputeNDof();
   
    virtual void CalcShape (const IntegrationPoint & ip, int facet, SliceMatrix<> shape ) const;
//     virtual int GetNExtraShapes( int facet) const;
//     virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> xshape) const;    
    
  };

}



#endif


