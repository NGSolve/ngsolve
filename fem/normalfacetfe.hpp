#ifndef NORMAL_FACET_FE_HPP__
#define NORMAL_FACET_FE_HPP__

/*********************************************************************/
/* File:   normalfacetfe.hpp                                         */
/* Author: J. Schoeberl                                              */
/* Date:   2019                                                      */
/*********************************************************************/


// #include <fem.hpp>
#include "thdivfe.hpp"
#include "recursive_pol.hpp"
#include <cassert>

namespace ngfem 
{
  /*
    facet element with normal facet components.
  */
  
  template <ELEMENT_TYPE ET>
  class NormalFacetFacetFE :
    public HDivNormalFiniteElement<ET_trait<ET>::DIM>,
    public VertexOrientedFE<ET>,
    public ET_trait<ET>
  {
  protected:
    IVec<2> order_inner;
    using VertexOrientedFE<ET>::vnums;
    using ET_trait<ET>::DIM;
    using HDivNormalFiniteElement<ET_trait<ET>::DIM>::order;
 
  public:
    using VertexOrientedFE<ET>::SetVertexNumber;
    using VertexOrientedFE<ET>::SetVertexNumbers;

    NormalFacetFacetFE (int aorder) : HDivNormalFiniteElement<ET_trait<ET>::DIM>(aorder+1,aorder)
    {
      order = aorder;
      order_inner = IVec<2>(aorder,aorder);
      ComputeNDof();
    }

    NormalFacetFacetFE () : HDivNormalFiniteElement<ET_trait<ET>::DIM>(0,0) { ; }

    HD virtual ELEMENT_TYPE ElementType() const override { return ELEMENT_TYPE(ET); }

    INLINE void SetOrder (int aorder)
    {
      order = aorder;
      order_inner = IVec<2>(aorder,aorder);
      ComputeNDof();
    }
  
    INLINE void SetOrder (IVec<2> oi)
    {
      order = max2 (oi[0], oi[1]);
      order_inner = oi;
      ComputeNDof();
    }

    virtual void ComputeNDof ();

    virtual void CalcShape(const IntegrationPoint & ip,
         		    FlatVector<> shape) const override;

    template<typename Tx, typename TFA>  
    void T_CalcShape (TIP<DIM,Tx> tip, TFA & shape) const;
  };


  

  template <ELEMENT_TYPE ET> class NormalFacetVolumeFE_Shape;
  
  template <ELEMENT_TYPE ET>
  class NormalFacetVolumeFE : public T_HDivFiniteElement<NormalFacetVolumeFE_Shape<ET>, ET>,
                              public ET_trait<ET>, public VertexOrientedFE<ET>
  {
  protected:
    using ET_T = ET_trait<ET>;
    IVec<2> facet_order[ET_T::N_FACET];
    int first_facet_dof[ET_T::N_FACET+1];
    bool highest_order_dc;
    using HDivFiniteElement<ET_trait<ET>::DIM>::order;
    using VertexOrientedFE<ET>::vnums;
    enum { DIM = ET_trait<ET>::DIM };
    typedef T_HDivFiniteElement<NormalFacetVolumeFE_Shape<ET>, ET> TBASE;
    
  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    
    NormalFacetVolumeFE () { highest_order_dc=false; }
    
    HD virtual ELEMENT_TYPE ElementType() const override { return ELEMENT_TYPE(ET); }

    void SetHighestOrderDC(bool set) { highest_order_dc=set; }

    void SetOrder(int ao)
    {
      order = ao;
      for ( int i = 0; i < ET_T::N_FACET; i++ )
        facet_order[i] = IVec<2> (ao, ao);
      ComputeNDof();
    }

    void SetOrder(FlatArray<int> & ao)
    {
      order = 0;
      assert(ao.Size()==ET_T::N_FACET);
      for ( int i = 0; i < ET_T::N_FACET; i++ )
        {
          order = max2 ( order, ao[i] );
          facet_order[i] = IVec<2> (ao[i], ao[i]);
        }
      ComputeNDof();
    }

    void SetOrder(FlatArray<IVec<2> > & ao)
    {
      order = 0;
      assert(ao.Size()==ET_T::N_FACET);
      for ( int i = 0; i < ET_T::N_FACET; i++ )
        {
          order = max3 ( order, ao[i][0], ao[i][1] );
          facet_order[i] = ao[i];
        }
      ComputeNDof();
    }

    IVec<2> GetFacetOrder(int j) const { return facet_order[j]; }
    int GetVertexNumber(int j) const { return vnums[j]; }

    /*
    virtual void CalcShape (const IntegrationPoint & ip, SliceMatrix<> shape) const override
    {
      int fnr = ip.FacetNr();
      if (fnr >= 0)
        {
          CalcShape (ip, fnr, shape);
          return;
        }
      throw Exception("NormalFacetVolumeFiniteElement<D>::CalcShape in global coordinates disabled");
    }

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                  BareSliceMatrix<SIMD<double>> shapes) const override;
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs,
                           BareSliceMatrix<SIMD<double>> values) const override;
    
    virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                           BareSliceVector<> coefs) const override;
    */
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[DIM], int fnr, TFA & shape) const;

    template<typename MIP, typename TFA>  
    void CalcDualShape2 (MIP mip, int fnr, TFA & shape) const;

    using TBASE::CalcShape;
    virtual void CalcShape (const IntegrationPoint & ip, int facet, SliceMatrix<> shape) const;

    using TBASE::CalcDualShape;
    virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const override;

    virtual int GetNExtraShapes( int facet) const {return 0;}
    virtual void CalcExtraShape (const IntegrationPoint & ip, int facet, FlatMatrixFixWidth<ET_T::DIM> xshape) const {xshape = 0.0;}

    virtual void GetFacetDofNrs(int afnr, Array<int>& fdnums) const
    {
      int first = first_facet_dof[afnr];
      int next = first_facet_dof[afnr+1];
      assert(next > first);
      fdnums.SetSize(next-first);
      for (auto i : Range(next-first))
        fdnums[i] = first+i;
    }

    virtual int GetFacetNDof(int afnr) const
    { return first_facet_dof[afnr+1] - first_facet_dof[afnr]; };
    virtual int GetFirstFacetDof(int afnr) const {
      return first_facet_dof[afnr];};

    virtual void ComputeNDof ();

    /// degrees of freedom sitting inside the element, used for static condensation
    virtual void GetInternalDofs (Array<int> & idofs) const
    {
      idofs.SetSize(0);
      if (highest_order_dc)
        {
          if (ET_T::DIM==2)
            {
              for (int i=0; i < ET_T::N_FACET; i++)
                idofs.Append(first_facet_dof[i+1]-1);
            }
          else
            {
              throw Exception ("NormalFacetFE with hodc not ready in 3D");
              for (int i=0; i < ET_T::N_FACET; i++)
                {
                  int pos = first_facet_dof[i]-2;
                  int fac = 4 - ElementTopology::GetNVertices(ElementTopology::GetFaceType(ET,i));
                  for (int k = 0; k <= facet_order[i][0]; k++){
                    pos += 2*(facet_order[i][0]+1-fac*k);
                    idofs.Append(pos);
                    idofs.Append(pos+1);
                  }
                }
            }//end if Dimension
        }
    }
  };

}



#endif


