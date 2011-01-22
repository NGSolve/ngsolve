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
  class NGS_DLL_HEADER FacetVolumeFiniteElement : virtual public ScalarFiniteElement<D>
  {
  protected:
    int vnums[8];
    int facet_order[6]; 
    int first_facet_dof[7];

    // bool highest_order_dc;
    using ScalarFiniteElement<D>::eltype;
    using ScalarFiniteElement<D>::order;

    ScalarFiniteElement<D-1> * facets[6];

    mutable int facetnr;

  public:
    FacetVolumeFiniteElement (ELEMENT_TYPE aeltype);

    void SetVertexNumbers (FlatArray<int> & avnums);
    // void SetHighestOrderDC(bool set){highest_order_dc=set;}

    void SetOrder (int ao);
    void SetOrder (FlatArray<int> & ao);
    int GetFacetOrder (int j) const { return facet_order[j]; }
    int GetVertexNumber (int j) const { return vnums[j]; }
   
    void CalcFacetShape(int fnr, const IntegrationPoint & ip, FlatVector<> shape) const;
    void EvaluateFacet (int fnr, const IntegrationRule & ir, FlatVector<> coefs, FlatVector<> values) const;
    void EvaluateFacetTrans (int fnr, const IntegrationRule & ir, FlatVector<> values, FlatVector<> coefs) const;
    
    
    void SelectFacet (int afn) const { facetnr = afn; }


    void GetFacetDofNrs(int afnr, Array<int>& fdnums) const; 
    int GetFacetNDof(int afnr) const { return first_facet_dof[afnr+1] - first_facet_dof[afnr]; };
    int GetFirstFacetDof(int afnr) const { return first_facet_dof[afnr]; } 
  
    const ScalarFiniteElement<D-1> & GetFacetFE(int fnr) const { return *facets[fnr]; }

    virtual void ComputeNDof () = 0;

    // virtual void GetInternalDofs (Array<int> & idofs) const; 

  };




  
  

  template <int DIM>
  class NGS_DLL_HEADER EdgeVolumeFiniteElement : virtual public ScalarFiniteElement<DIM>
  {
  protected:
    int vnums[8];
    mutable int edgenr;

    using ScalarFiniteElement<DIM> :: order;
    using ScalarFiniteElement<DIM> :: eltype;
    using ScalarFiniteElement<DIM> :: ndof;
  public:
    EdgeVolumeFiniteElement (ELEMENT_TYPE aeltype, int aorder);
    void SetVertexNumbers (FlatArray<int> & avnums)
    {
      for (int i = 0; i < avnums.Size(); i++)
	vnums[i] = avnums[i];
    }
    void CalcEdgeShape(int enr, const IntegrationPoint & ip, FlatVector<> shape) const;

    void SelectEdge (int nr) const { edgenr = nr; }
  };

  class NGS_DLL_HEADER EdgeVolumeTet : public EdgeVolumeFiniteElement<3>,
				       public T_ScalarFiniteElement2<EdgeVolumeTet, ET_TET>
  {
  public:
    EdgeVolumeTet (int aorder)
      : EdgeVolumeFiniteElement<3>(ET_TET, aorder) { ; }


    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const
    {
      Tx lam[4] = { hx[0], hx[1], hx[2], 1-hx[0]-hx[1]-hx[2] };

      INT<2> e = ET_trait<ET_TET> :: GetEdgeSort (edgenr, vnums);

      SetZero (shape, 0, ndof);

      int base = (order+1) * edgenr;
      LegendrePolynomial (order, 
			  lam[e[1]] - lam[e[0]],
			  shape.Addr(base));
    }
  };


  class NGS_DLL_HEADER EdgeVolumeTrig : public EdgeVolumeFiniteElement<2>,
					public T_ScalarFiniteElement2<EdgeVolumeTrig, ET_TRIG>
  {
  public:
    EdgeVolumeTrig (int aorder)
      : EdgeVolumeFiniteElement<2>(ET_TRIG, aorder) { ; }


    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[2], TFA & shape) const
    {
      Tx lam[4] = { hx[0], hx[1], 1-hx[0]-hx[1] };

      INT<2> e = ET_trait<ET_TRIG> :: GetEdgeSort (edgenr, vnums);
      
      SetZero (shape, 0, ndof);

      int base = (order+1) * edgenr;
      LegendrePolynomial (order, 
			  lam[e[1]] - lam[e[0]],
			  shape.Addr(base));
    }
  };


}



#endif
