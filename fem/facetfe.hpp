#ifndef FILE_FACETFE
#define FILE_FACETFE

/*********************************************************************/
/* File:   facetfe.hpp                                               */
/* Author: A. Sinwel, H. Egger, J. Schoeberl                         */
/* Date:   2008                                                      */
/*********************************************************************/


#include "scalarfe.hpp"

namespace ngfem
{

  /*
   * Facet Finite Elements
   */ 

  template <int D>
  class NGS_DLL_HEADER FacetVolumeFiniteElement;
  

  template <int D>
  class FacetFEFacet : public ScalarFiniteElement<D>
  {
    int fnr;
    const FacetVolumeFiniteElement<D> & fe;
  public:
    FacetFEFacet (int afnr,
		  const FacetVolumeFiniteElement<D> & afe,
		  int andof, int aorder)
      : ScalarFiniteElement<D> (andof, aorder), fnr(afnr), fe(afe) 
    { 
      ; // cout << "created facetfefacet" << endl;
    }

    HD virtual ELEMENT_TYPE ElementType() const override { return fe.ElementType(); }

    using ScalarFiniteElement<D>::CalcShape;
    HD virtual void CalcShape (const IntegrationPoint & ip, 
                               BareSliceVector<> shape) const override
    {
      fe.CalcFacetShapeVolIP(fnr, ip, shape);
    }
    
    HD virtual void CalcShape (const SIMD_IntegrationRule & ir, 
                               BareSliceMatrix<SIMD<double>> shape) const override
    {
      fe.CalcFacetShapeVolIR(fnr, ir, shape);
    }
    
    HD virtual void CalcDShape (const IntegrationPoint & ip, 
                                BareSliceMatrix<> dshape) const override
    {
      fe.CalcFacetDShapeVolIP(fnr, ip, dshape);
      //throw Exception ("facetfe - calcdshape not olverloaded");
    }

    using ScalarFiniteElement<D>::Evaluate;
    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs, BareVector<SIMD<double>> values) const override
    {
      fe.EvaluateFacetVolIp (fnr, ir, coefs, values);
    }
    
    using ScalarFiniteElement<D>::AddTrans;    
    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_IntegrationRule & ir, BareVector<SIMD<double>> values, BareSliceVector<> coefs) const override
    {
      fe.AddTransFacetVolIp (fnr, ir, values, coefs);
    }
    
  };

  
  

  template <int D>
  class FacetVolumeFiniteElement : public FiniteElement
  {
  protected:
    // int vnums[8];
    int facet_order[6]; 
    int first_facet_dof[7];

    using FiniteElement::ndof;
    using FiniteElement::order;

  public:

    /*
    void SetVertexNumbers (FlatArray<int> & avnums)
    {
      for (int i = 0; i < avnums.Size(); i++)
	vnums[i] = avnums[i];
    }
    */

    /*
    template <typename T>
    void SetVertexNumbers (const BaseArrayObject<T> & avnums)
    {
      for (int i = 0; i < avnums.Size(); i++)
	vnums[i] = avnums[i];
    }
    */
    
    void SetOrder (int ao)  
    {
      order = ao;
      for (int i = 0; i < 6; i++)
	facet_order[i] = ao;
    }

    void SetOrder (int nr, int o) { facet_order[nr] = o; order = max2(order,o); }
    
    void SetOrder (FlatArray<int> & ao)
    {
      for (int i=0; i<ao.Size(); i++)
	facet_order[i] = ao[i];
      
      order = facet_order[0];        // integration order
      for (int i = 1; i < ao.Size(); i++)
	order = max2(order, ao[i]);
    }

    FacetFEFacet<D> Facet (int fnr) const 
    {
      return FacetFEFacet<D> (fnr, *this, 
			      GetFacetDofs(fnr).Size(), facet_order[fnr]); 
    }


    virtual void CalcFacetShapeVolIP (int fnr, const IntegrationPoint & ip, 
				      BareSliceVector<> shape) const = 0;
    virtual void CalcFacetShapeVolIR (int fnr, const SIMD_IntegrationRule & ir, 
                                      BareSliceMatrix<SIMD<double>> shape) const = 0;

    virtual void EvaluateFacetVolIp (int fnr, const SIMD_IntegrationRule & ir, BareSliceVector<> coefs, BareVector<SIMD<double>> values) const = 0;
    virtual void AddTransFacetVolIp (int fnr, const SIMD_IntegrationRule & ir, BareVector<SIMD<double>> values, BareSliceVector<> coefs) const = 0;
    
    virtual void CalcFacetDShapeVolIP (int fnr, const IntegrationPoint & ip, 
                                       BareSliceMatrix<> shape) const = 0;

    IntRange GetFacetDofs(int fnr) const
    {
      return IntRange (first_facet_dof[fnr], first_facet_dof[fnr+1]);
    }

    virtual string ClassName() const { return "FacetVolumeFiniteElement"; }

    virtual void ComputeNDof () = 0;
  };




#ifdef FILE_FACETHOFE_CPP
#define FACETHOFE_EXTERN
#else
#define FACETHOFE_EXTERN extern
#endif

  FACETHOFE_EXTERN template class FacetVolumeFiniteElement<1>;
  FACETHOFE_EXTERN template class FacetVolumeFiniteElement<2>;
  FACETHOFE_EXTERN template class FacetVolumeFiniteElement<3>;

  FACETHOFE_EXTERN template class FacetFEFacet<1>;
  FACETHOFE_EXTERN template class FacetFEFacet<2>;
  FACETHOFE_EXTERN template class FacetFEFacet<3>;

}



#endif
