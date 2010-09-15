#ifndef VECTOR_FACET_FE_HPP__
#define VECTOR_FACET_FE_HPP__

/*********************************************************************/
/* File:   vectorfacetfe.hpp                                         */
/* Author: A. Sinwel, (J. Schoeberl)                                 */
/* Date:   2008                                                      */
/*********************************************************************/


namespace ngfem 
{

  // #include <fem.hpp>

  class NGS_DLL_HEADER VectorFacetFacetFiniteElement : public FiniteElement
  {
  protected:
    int vnums[8];

    INT<2> order_inner;

  public:

    VectorFacetFacetFiniteElement ():
      FiniteElement(),
      order_inner( INT<2>(0,0) )
    {
      for ( int i = 0; i < 8; i++ )
	vnums[i] = -1;
    }

    VectorFacetFacetFiniteElement (int dim, ELEMENT_TYPE aeltype);

    void SetVertexNumbers (FlatArray<int> & avnums);

    void SetOrder (int o);
  
    void SetOrder (INT<2> oi);

    virtual void ComputeNDof () = 0; 

    virtual void CalcShape ( const IntegrationPoint & ip, FlatMatrix<> shape ) const = 0;

    virtual const FlatMatrix<> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const = 0;
    // const int * GetVNums() const { return vnums; }
 
  };



  /**
     High order 1D finite element
  */
  class NGS_DLL_HEADER VectorFacetFacetSegm : public VectorFacetFacetFiniteElement
  {
  public:
    VectorFacetFacetSegm (int aorder=0);

    virtual void ComputeNDof();

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrix<> shape) const;
			    
    virtual const FlatMatrix<> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const
    {
      FlatMatrix<> shape(ndof, 1, lh);
      CalcShape (ip, shape);
      return shape;
    }

  };


  /**
     High order triangular finite element
  */
  class VectorFacetFacetTrig : public VectorFacetFacetFiniteElement
  {
  public:
    VectorFacetFacetTrig (int aorder=0);
    virtual void ComputeNDof();

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrix<> shape) const;

    virtual const FlatMatrix<> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const
    {
      FlatMatrix<> shape(ndof, 2, lh);
      CalcShape (ip, shape);
      return shape;
    }
  };


  /**
     High order quadrilateral finite element
  */
  class VectorFacetFacetQuad : public VectorFacetFacetFiniteElement
  {
  public:
    VectorFacetFacetQuad (int aorder=0);
    virtual void ComputeNDof();

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrix<> shape) const;

    virtual const FlatMatrix<> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const
    {
      FlatMatrix<> shape(ndof, 2, lh);
      CalcShape (ip, shape);
      return shape;
    }

  };


  template <int D>
  class NGS_DLL_HEADER VectorFacetVolumeFiniteElement : public HCurlFiniteElement<D>
  {
  protected:
    int vnums[8];
    INT<2> facet_order[6]; 
    int first_facet_dof[7];
    // int nv; // num of vertices
    // int nf; // num of facets
    bool highest_order_dc;
    using HCurlFiniteElement<D>::eltype;
    using HCurlFiniteElement<D>::order;
  
  public:
    VectorFacetVolumeFiniteElement () // : nv(0), nf(0)
    { highest_order_dc=false; }

    VectorFacetVolumeFiniteElement (ELEMENT_TYPE aeltype);

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

    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrixFixWidth<D> shape) const = 0;
    virtual void CalcShape (const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<D> shape) const = 0;
    
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



  class VectorFacetVolumeTrig : public VectorFacetVolumeFiniteElement<2>
  {
  public:
    VectorFacetVolumeTrig() : VectorFacetVolumeFiniteElement<2>(ET_TRIG) { ; }
    virtual void ComputeNDof();
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrixFixWidth<2> shape) const; 
    virtual void CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<2> shape) const;
    virtual int GetNExtraShapes( int facet) const {return 1;};
    virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<2> xshape) const;
  };



  class VectorFacetVolumeQuad : public VectorFacetVolumeFiniteElement<2>
  {
  public:
    VectorFacetVolumeQuad() : VectorFacetVolumeFiniteElement<2> (ET_QUAD) { ; }
    virtual void ComputeNDof();
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrixFixWidth<2> shape) const; 
    virtual void CalcShape (const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<2> shape) const;
    virtual int GetNExtraShapes( int facet) const {return 1;};
    virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<2> xshape) const;
    
  };


  class VectorFacetVolumeTet : public VectorFacetVolumeFiniteElement<3>
  {
  public:
    VectorFacetVolumeTet() : VectorFacetVolumeFiniteElement<3> (ET_TET) { ; }
    virtual void ComputeNDof();
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrixFixWidth<3> shape) const; 
    virtual void CalcShape (const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> shape) const;
    virtual int GetNExtraShapes( int facet) const {return 2*(facet_order[facet][0]+2);};
    virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> xshape) const;    
  };



  class VectorFacetVolumeHex : public VectorFacetVolumeFiniteElement<3>
  {
  public:
    VectorFacetVolumeHex() : VectorFacetVolumeFiniteElement<3>(ET_HEX) { ; };
    virtual void ComputeNDof();
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrixFixWidth<3> shape) const; 
    virtual void CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> shape ) const;
//     virtual int GetNExtraShapes( int facet) const {return 2*(2*facet_order[facet][0]+3);};
//     virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> xshape) const;    
  };



  class VectorFacetVolumePrism : public VectorFacetVolumeFiniteElement<3>
  {
  public:
    VectorFacetVolumePrism() : VectorFacetVolumeFiniteElement<3> (ET_PRISM) { ; };
    virtual void ComputeNDof();
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrixFixWidth<3> shape) const; 
    virtual void CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> shape) const;
//     virtual int GetNExtraShapes( int facet) const;
//     virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> xshape) const;    
  };


  class VectorFacetVolumePyramid : public VectorFacetVolumeFiniteElement<3>
  {
  public:
    VectorFacetVolumePyramid() : VectorFacetVolumeFiniteElement<3> (ET_PYRAMID) { ; };
    virtual void ComputeNDof();
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrixFixWidth<3> shape) const; 
    virtual void CalcShape (const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> shape ) const;
//     virtual int GetNExtraShapes( int facet) const;
//     virtual void CalcExtraShape ( const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<3> xshape) const;    
    
  };

}



#endif


