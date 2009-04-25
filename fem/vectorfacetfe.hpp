#ifndef VECTOR_FACET_FE_HPP__
#define VECTOR_FACET_FE_HPP__

namespace ngfem 
{

  // #include <fem.hpp>

  class VectorFacetFacetFiniteElement : public FiniteElement
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
  class VectorFacetFacetSegm : public VectorFacetFacetFiniteElement
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



  class VectorFacetVolumeFiniteElement : public FiniteElement  
  {
  protected:
    int vnums[8];
    INT<2> facet_order[6]; 
    int first_facet_dof[7];
    int nv; // num of vertices
    int nf; // num of facets
  
  public:
    VectorFacetVolumeFiniteElement () :
      FiniteElement (),
      nv(0),
      nf(0)
    { ; }

    VectorFacetVolumeFiniteElement (int adim, ELEMENT_TYPE aeltype);

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

    virtual void CalcShape(const IntegrationPoint & ip, FlatMatrix<> shape) const = 0;
      
    //   virtual void SetFacet(int afnr) const = 0; 

    //   virtual void CalcFacetShape(int fnr, const IntegrationPoint & ip, 
    // 			      FlatMatrix<> shape) const = 0;

    virtual void GetFacetDofNrs(int afnr, Array<int>& fdnums) const; 

    virtual int GetFacetNDof(int afnr) const 
    { return first_facet_dof[afnr+1] - first_facet_dof[afnr]; };

    virtual int GetFirstFacetDof(int afnr) const { 
      return first_facet_dof[afnr];}; 
  
    virtual void ComputeNDof () = 0;
  
    // utility
    virtual int GetNF() const { return nf; };

    virtual int GetNV() const { return nv; };

    virtual void GetVertexNumbers(Array<int>&) const;

    virtual void GetFacetOrders(Array<INT<2> >&) const;
  
    virtual const FlatMatrix<> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const = 0;

    virtual void CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrix<> shape ) const = 0;

  
  };



  class VectorFacetVolumeTrig : public VectorFacetVolumeFiniteElement
  {
  protected:
    //     int fnr; // active facet
    //     VectorFacetFacetSegm facet;
  public:
    VectorFacetVolumeTrig() : VectorFacetVolumeFiniteElement(2, ET_TRIG) { ; };
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const; 
    //     virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, 
    // 				 FlatMatrix<> shape) const; 
    //     virtual void SetFacet(int afnr) const;
    
    virtual void ComputeNDof();

    //   virtual void TransformFacetToVolumeShape ( int fanr, FlatMatrix<> shape1d, 
    // 					     FlatMatrix<> shape ) const;
    virtual const FlatMatrix<> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const
    {
      FlatMatrix<> shape(ndof, 2, lh);
      CalcShape (ip, shape);
      return shape;
    }

    virtual void CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrix<> shape ) const;

  };

  // --------------------------------------------------------
  class VectorFacetVolumeQuad : public VectorFacetVolumeFiniteElement
  {
  protected:
    //     int fnr; // active facet
    //     VectorFacetFacetSegm facet;
  public:
    VectorFacetVolumeQuad() : VectorFacetVolumeFiniteElement(2, ET_QUAD) { ; };
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const; 
    //     virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, 
    // 				 FlatMatrix<> shape) const; 
    //     virtual void SetFacet(int afnr) const;
   
    virtual void ComputeNDof();

    //   virtual void TransformFacetToVolumeShape ( int fanr, FlatMatrix<> shape1d, 
    // 					     FlatMatrix<> shape ) const;

    virtual const FlatMatrix<> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const
    {
      FlatMatrix<> shape(ndof, 2, lh);
      CalcShape (ip, shape);
      return shape;
    }

    virtual void CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrix<> shape ) const;
  };

  // --------------------------------------------------------
  class VectorFacetVolumeTet : public VectorFacetVolumeFiniteElement
  {
  protected:
    //     int fnr; // active facet
    //     VectorFacetFacetTrig facet;
  public:
    VectorFacetVolumeTet() : VectorFacetVolumeFiniteElement(3, ET_TET) { ; };
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const; 

    //     virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, 
    // 				 FlatMatrix<> shape) const; 
    //     virtual void SetFacet(int afnr) const;
   
    virtual void ComputeNDof();
    //     virtual const FiniteElement & GetFacetFE(int fnr, LocalHeap& lh) const;

    virtual const FlatMatrix<> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const
    {
      FlatMatrix<> shape(ndof, 3, lh);
      CalcShape (ip, shape);
      return shape;
    }
    virtual void CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrix<> shape ) const;
  };

  // --------------------------------------------------------
  class VectorFacetVolumeHex : public VectorFacetVolumeFiniteElement
  {
  protected:
    //     int fnr; // active facet
    //     VectorFacetFacetQuad facet;
  public:
    VectorFacetVolumeHex() : VectorFacetVolumeFiniteElement(3, ET_HEX) { ; };
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const; 

    //     virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, 
    // 				 FlatMatrix<> shape) const; 
    //     virtual void SetFacet(int afnr) const;
   
    virtual void ComputeNDof();
    //     virtual const FiniteElement & GetFacetFE(int fnr, LocalHeap& lh) const;

    virtual const FlatMatrix<> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const
    {
      FlatMatrix<> shape(ndof, 3, lh);
      CalcShape (ip, shape);
      return shape;
    }
    virtual void CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrix<> shape ) const;
  };

  // --------------------------------------------------------
  class VectorFacetVolumePrism : public VectorFacetVolumeFiniteElement
  {
  protected:
    //     int tnr, qnr; // active facets
    //     VectorFacetFacetTrig trig;
    //     VectorFacetFacetQuad quad;
  public:
    VectorFacetVolumePrism() : VectorFacetVolumeFiniteElement(3, ET_PRISM) 
    { ; };
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const; 

    //     virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, 
    // 				 FlatMatrix<> shape) const; 
    //     virtual void SetFacet(int afnr) const;
   
    virtual void ComputeNDof();
    //     virtual const FiniteElement & GetFacetFE(int fnr, LocalHeap& lh) const;

    virtual const FlatMatrix<> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const
    {
      FlatMatrix<> shape(ndof, 3, lh);
      CalcShape (ip, shape);
      return shape;
    }
    virtual void CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrix<> shape ) const;
  };

  // --------------------------------------------------------
  class VectorFacetVolumePyramid : public VectorFacetVolumeFiniteElement
  {
  protected:
    //     int tnr, qnr; // active facets
    //     VectorFacetFacetTrig trig;
    //     VectorFacetFacetQuad quad;
  public:
    VectorFacetVolumePyramid() : VectorFacetVolumeFiniteElement(3, ET_PYRAMID) 
    { ; };
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const; 

    //     virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, 
    // 				 FlatMatrix<> shape) const; 
    //     virtual void SetFacet(int afnr) const;
   
    virtual void ComputeNDof();

    virtual const FlatMatrix<> GetShape (const IntegrationPoint & ip, 
					 LocalHeap & lh) const
    {
      FlatMatrix<> shape(ndof, 3, lh);
      CalcShape (ip, shape);
      return shape;
    }
    virtual void CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrix<> shape ) const;
  };

}



#endif


