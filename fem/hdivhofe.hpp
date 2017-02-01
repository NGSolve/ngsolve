#ifndef FILE_HDIVHOFE_
#define FILE_HDIVHOFE_ 

/*********************************************************************/
/* File:   hdivhofe.hpp                                              */
/* Author: A. Becirovic, S. Zaglmayr, J. Schoeberl                   */
/* Date:   15. Feb. 2003                                             */
/*********************************************************************/


#include "thdivfe.hpp"

namespace ngfem
{
  








  template <int D>
  class HDivHighOrderNormalFiniteElement : public HDivNormalFiniteElement<D>
  {
  protected:
    int vnums[4];
    INT<2> order_inner;

  public:
    ///
    HDivHighOrderNormalFiniteElement ();

    void SetVertexNumbers (FlatArray<int> & avnums);

    void SetOrderInner (int oi);
    void SetOrderInner (INT<2> oi);

    virtual void ComputeNDof () = 0;

    int EdgeOrientation (int enr) const
    {
      const EDGE * edges = ElementTopology::GetEdges (this->ElementType());
      return (vnums[edges[enr][1]] > vnums[edges[enr][0]]) ? 1 : -1;
    }
  };


  template <class T_ORTHOPOL = TrigExtensionMonomial>
  class HDivHighOrderNormalSegm : public HDivHighOrderNormalFiniteElement<1>
  {
  public:

    HDivHighOrderNormalSegm (int aorder);
    virtual void ComputeNDof();
    virtual ELEMENT_TYPE ElementType() const { return ET_SEGM; }

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip,
                            FlatVector<> shape) const;

  };

  template <class T_ORTHOPOL = TrigExtensionMonomial>
  class HDivHighOrderNormalTrig : public HDivHighOrderNormalFiniteElement<2>
  {
  public:

    HDivHighOrderNormalTrig (int aorder);
    virtual void ComputeNDof();
    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }
    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip,
                            FlatVector<> shape) const;
  };

  template <class T_ORTHOPOL = TrigExtensionMonomial>
  class HDivHighOrderNormalQuad : public HDivHighOrderNormalFiniteElement<2>
  {
  public:

    HDivHighOrderNormalQuad (int aorder);
    virtual void ComputeNDof();
    virtual ELEMENT_TYPE ElementType() const { return ET_QUAD; }
    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip,
                            FlatVector<> shape) const;
  };




  template <ELEMENT_TYPE ET> class HDivHighOrderFE_Shape;

  template <ELEMENT_TYPE ET> 
  class NGS_DLL_HEADER HDivHighOrderFE : 
    public T_HDivFiniteElement< HDivHighOrderFE_Shape<ET>, ET > , public ET_trait<ET>, public VertexOrientedFE<ET>
  {
  protected:
    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::N_FACET;
    using ET_trait<ET>::DIM;

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;  

    using HDivFiniteElement<DIM>::ndof;
    using HDivFiniteElement<DIM>::order;

    using VertexOrientedFE<ET>::vnums;
    

    INT<DIM> order_inner;
    INT<N_FACET,INT<DIM-1>> order_facet;  

    bool ho_div_free;
    bool only_ho_div;

  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    /// minimal constructor, orders will be set later
    HDivHighOrderFE () 
      : ho_div_free(false), only_ho_div(false)
    { ; }
  
    /// builds a functional element of order aorder.
    HDivHighOrderFE (int aorder)
      : ho_div_free(false), only_ho_div(false)
    { 
      for (int i = 0; i < N_VERTEX; i++) vnums[i] = i;

      order_inner = aorder;
      order_facet = aorder;

      ComputeNDof();
    }


    void SetOrderInner (INT<DIM> oi)
    { 
      order_inner = oi; 
    }

    template <typename TA>
    void SetOrderFacet (const TA & oe)
    { 
      for (int i = 0; i < N_FACET; i++) 
        order_facet[i] = oe[i]; 
    }

    void SetHODivFree (bool aho_div_free) 
    { 
      ho_div_free = aho_div_free; 
      only_ho_div = only_ho_div && !ho_div_free;
    };  

    void SetOnlyHODiv (bool aonly_ho_div) 
    { 
      only_ho_div = aonly_ho_div; 
      ho_div_free = ho_div_free && !only_ho_div;
    };  

    virtual void ComputeNDof();
    virtual ELEMENT_TYPE ElementType() const { return ET; }
    virtual void GetFacetDofs(int i, Array<int> & dnums) const;

    /// calc normal components of facet shapes, ip has facet-nr
    virtual void CalcNormalShape (const IntegrationPoint & ip, 
                                  SliceVector<> nshape) const;

  };


  
  // still to be changed ....

#ifdef HDIVHEX
  template<> 
  class HDivHighOrderFE<ET_HEX> : 
    public HDivHighOrderFiniteElement<3>
  {
  public:
    HDivHighOrderFE () { ; }
    HDivHighOrderFE (int aorder);



    virtual void ComputeNDof();
    virtual ELEMENT_TYPE ElementType() const { return ET_HEX; }

    // virtual void GetInternalDofs (Array<int> & idofs) const;
  

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip,
                            SliceMatrix<> shape) const;

    /// compute Div of shape
    virtual void CalcDivShape (const IntegrationPoint & ip,
                               SliceVector<> shape) const;
    /// compute Div numerical diff
    //void CalcNumDivShape( const IntegrationPoint & ip,
    //			FlatVector<> divshape) const;
    virtual void GetFacetDofs(int i, Array<int> & dnums) const; 

  };
#endif

}



#ifdef FILE_HDIVHOFE_CPP

#define HDIVHOFE_EXTERN
#include <thdivfe_impl.hpp>
#include <hdivhofe_impl.hpp>
#include <hdivhofefo.hpp>

#else

#define HDIVHOFE_EXTERN extern

#endif


namespace ngfem
{
  HDIVHOFE_EXTERN template class HDivHighOrderFE<ET_TRIG>;
  HDIVHOFE_EXTERN template class HDivHighOrderFE<ET_QUAD>;
  HDIVHOFE_EXTERN template class HDivHighOrderFE<ET_TET>;
  HDIVHOFE_EXTERN template class HDivHighOrderFE<ET_PRISM>;

  HDIVHOFE_EXTERN template class T_HDivFiniteElement<HDivHighOrderFE_Shape<ET_TRIG>, ET_TRIG>;
  HDIVHOFE_EXTERN template class T_HDivFiniteElement<HDivHighOrderFE_Shape<ET_QUAD>, ET_QUAD>;
  HDIVHOFE_EXTERN template class T_HDivFiniteElement<HDivHighOrderFE_Shape<ET_TET>, ET_TET>;
  HDIVHOFE_EXTERN template class T_HDivFiniteElement<HDivHighOrderFE_Shape<ET_PRISM>, ET_PRISM>;
}

#endif



