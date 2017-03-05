#ifndef VECTOR_FACET_FE_HPP__
#define VECTOR_FACET_FE_HPP__

/*********************************************************************/
/* File:   vectorfacetfe.hpp                                         */
/* Author: A. Sinwel, (J. Schoeberl)                                 */
/* Date:   2008                                                      */
/*********************************************************************/


namespace ngfem 
{
  /*
    facet element with tangential facet components.
    to be changed similar to scalar facetfe
  */
  
  template <int D>
  //   class NGS_DLL_HEADER VectorFacetFacetFiniteElement : public FiniteElement
  class VectorFacetFacetFiniteElement : public HCurlFiniteElement<D>
  {
  protected:
    int vnums[8];

    INT<2> order_inner;
    // using HCurlFiniteElement<D>::eltype;
    using HCurlFiniteElement<D>::order;
 
  public:

    VectorFacetFacetFiniteElement ()
      : order_inner (0)
    {
      for ( int i = 0; i < 8; i++ )
	vnums[i] = -1;
    }

    VectorFacetFacetFiniteElement (int dim, ELEMENT_TYPE aeltype)
      : HCurlFiniteElement<D> (-1, -1)
                                  // : FiniteElement (aeltype, -1, -1 )
    {
      for (int i=0; i<8; i++) vnums[i] = -1; 
      order_inner = -1;
    }

    INLINE void SetVertexNumbers (FlatArray<int> & avnums)
    {
      for ( int i = 0; i < avnums.Size(); i++ ) vnums[i] = avnums[i];
    }

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

    virtual void ComputeNDof () = 0; 

    virtual void CalcShape (const IntegrationPoint & ip, SliceMatrix<> shape) const = 0;

    virtual const FlatMatrixFixWidth<D> GetShape (const IntegrationPoint & ip, 
                                                  LocalHeap & lh) const = 0;
  };



  /**
     High order 1D finite element
  */
  class VectorFacetFacetSegm : public VectorFacetFacetFiniteElement<1>, public VertexOrientedFE<ET_SEGM>
  {
    using VertexOrientedFE<ET_SEGM>::vnums;
  public:
    VectorFacetFacetSegm (int aorder=0);

    virtual void ComputeNDof();
    using VertexOrientedFE<ET_SEGM>::SetVertexNumbers;

    HD virtual ELEMENT_TYPE ElementType() const { return ET_SEGM; }

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
			    SliceMatrix<> shape) const;
			    
    virtual const FlatMatrixFixWidth<1> GetShape (const IntegrationPoint & ip, 
						  LocalHeap & lh) const
    {
      FlatMatrixFixWidth<1> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

  };


  /**
     High order triangular finite element
  */
  class VectorFacetFacetTrig : public VectorFacetFacetFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    using VertexOrientedFE<ET_TRIG>::vnums;
  public:
    VectorFacetFacetTrig (int aorder=0);
    virtual void ComputeNDof();
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;

    HD virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }
    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
			    SliceMatrix<> shape) const;

    virtual const FlatMatrixFixWidth<2> GetShape (const IntegrationPoint & ip, 
						  LocalHeap & lh) const
    {
      FlatMatrixFixWidth<2> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }
  };


  /**
     High order quadrilateral finite element
  */
  class VectorFacetFacetQuad : public VectorFacetFacetFiniteElement<2>
  {
  public:
    VectorFacetFacetQuad (int aorder=0);
    virtual void ComputeNDof();
    HD virtual ELEMENT_TYPE ElementType() const { return ET_QUAD; }
    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
			    SliceMatrix<> shape) const;

    virtual const FlatMatrixFixWidth<2> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const
    {
      FlatMatrixFixWidth<2> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

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


