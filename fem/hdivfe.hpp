#ifndef FILE_HDIVFE
#define FILE_HDIVFE

/*********************************************************************/
/* File:   hdivfe.hpp                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   5. Jul. 2001                                              */
/*********************************************************************/

namespace ngfem
{

  /**
     Finite Elements for H(div)
     Raviart-Thomas, BDM, BDFM
  */


  class NGS_DLL_HEADER  BaseHDivFiniteElement : public FiniteElement
  {
  public:
    INLINE BaseHDivFiniteElement () { ; } 
    INLINE BaseHDivFiniteElement (int andof, int aorder)
      : FiniteElement (andof, aorder) { ; }

    virtual void CalcShape (const IntegrationPoint & ip, 
			    SliceMatrix<> shape) const = 0;

    virtual void CalcDivShape (const IntegrationPoint & ip,
			       SliceVector<> divshape) const = 0;
  };

  template <int D>
  class NGS_DLL_HEADER HDivFiniteElement : public BaseHDivFiniteElement
  {
  public:
    enum { DIM = D };

  public:
    ///
    INLINE HDivFiniteElement (int andof, int aorder)
      : BaseHDivFiniteElement (andof, aorder) { ; } 

    ///
    INLINE HDivFiniteElement () { ; }

    ///
    HD virtual ~HDivFiniteElement () { ; }

    /// 
    virtual string ClassName() const;

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip,
			    SliceMatrix<> shape) const = 0;

    /// compute div of shape
    virtual void CalcDivShape (const IntegrationPoint & ip,
			       SliceVector<> divshape) const;

    /// calc normal components of facet shapes, ip has facet-nr
    virtual void CalcNormalShape (const IntegrationPoint & ip, 
                                  SliceVector<> nshape) const;

    /// compute shape
    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
                                  SliceMatrix<> shape) const;

    virtual void CalcMappedShape (const SIMD<MappedIntegrationPoint<DIM,DIM>> & mip,
				  BareSliceMatrix<SIMD<double>> shape) const;

    virtual void CalcMappedShape (const BaseMappedIntegrationRule & bmir, SliceMatrix<> shapes) const;

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                  BareSliceMatrix<SIMD<double>> shapes) const;

    virtual void CalcMappedNormalShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                        BareSliceMatrix<SIMD<double>> shapes) const;

    /// compute div of shape
    virtual void CalcMappedDivShape (const MappedIntegrationPoint<DIM,DIM> & sip,
				     SliceVector<> divshape) const;

    virtual void CalcMappedDivShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                     BareSliceMatrix<SIMD<double>> divshapes) const;


    INLINE const FlatMatrixFixWidth<DIM> GetShape (const IntegrationPoint & ip,
                                                   LocalHeap & lh) const
    {
      FlatMatrixFixWidth<DIM> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

    INLINE const FlatVector<> GetDivShape (const IntegrationPoint & ip,
                                           LocalHeap & lh) const
    {
      FlatVector<> divshape(ndof, lh);
      CalcDivShape (ip, divshape);
      return divshape;
    }

    
    virtual void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, 
			   FlatMatrixFixWidth<D> vals) const;
    /*
    {
      MatrixFixWidth<D> shape(ndof);
      for (int i = 0; i < ir.GetNIP(); i++)
	{
	  CalcShape (ir[i], shape);
	  vals.Row(i) = Trans(shape) * coefs;
	}
    }
    */
    virtual void EvaluateTrans (const IntegrationRule & ir, 
                                FlatMatrixFixWidth<D> vals,
                                FlatVector<double> coefs) const;
    /*
    {
      MatrixFixWidth<D> shape(ndof);
      coefs = 0;
      for (int i = 0; i < ir.GetNIP(); i++)
	{
	  CalcShape (ir[i], shape);
          coefs += shape * vals.Row(i);
	}
    }
    */
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs,
                           BareSliceMatrix<SIMD<double>> values) const;
    virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                           BareSliceVector<> coefs) const;
    virtual void EvaluateDiv (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs,
                              BareVector<SIMD<double>> values) const;
    virtual void AddDivTrans (const SIMD_BaseMappedIntegrationRule & ir, BareVector<SIMD<double>> values,
                              BareSliceVector<> coefs) const;
    
    virtual void GetFacetDofs(int i, Array<int> & dnums) const;
    // { cout  << " GetFacetDofs for nothing " << endl; dnums.SetSize(0);}; 

    virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, SliceMatrix<> shape) const;

  protected:

    /// compute basis, will be orthogonalized
    virtual void CalcShape1 (const IntegrationPoint & ip,
			     FlatMatrixFixWidth<DIM> shape) const { ; }

    ///
    virtual void CalcShape2 (const IntegrationPoint & ip,
			     FlatMatrixFixWidth<DIM> shape) const { ; }

    ///
    void ComputeFaceMoments (int fnr, ScalarFiniteElement<DIM-1> & testfe,
			     FlatMatrix<> & moments,
			     int order, int shape = 1) const;


    virtual std::list<std::tuple<std::string,double>> Timing () const;
  };




  /**
    HDivNormalFiniteElement
  */

  template <int D>
  class NGS_DLL_HEADER HDivNormalFiniteElement : public FiniteElement
  {
  public:
    enum { DIM = D };

  public:
    ///
    HDivNormalFiniteElement (int andof, int aorder)
      : FiniteElement (andof, aorder){;}

    ///
    HD virtual ~HDivNormalFiniteElement () { ; }

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip,
			    FlatVector<> shape) const = 0;

    ///
    const FlatVector<> GetShape (const IntegrationPoint & ip,
				 LocalHeap & lh) const
    {
      FlatVector<> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
                                  SliceMatrix<> shape) const
    {
      throw Exception(string("CalcMappedShape not implemented for H(div) normal element ")+typeid(*this).name());
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           BareSliceVector<> coefs,
                           BareSliceMatrix<SIMD<double>> values) const
    {
      throw ExceptionNOSIMD ("HDivNormalFE::Evaluate (simd) not overloaded");
    }
    
    // values * normal
    virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir,
                           BareSliceMatrix<SIMD<double>> values,
                           BareSliceVector<> coefs) const
    {
      throw ExceptionNOSIMD ("HDivNormalFE::AddTrans (simd) not overloaded");
    }
    
  };



}

#endif



