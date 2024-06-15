#ifndef FILE_HCURLFE
#define FILE_HCURLFE

/*********************************************************************/
/* File:   hcurlfe.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   16. Apr. 2000                                             */
/*********************************************************************/

#include "finiteelement.hpp"


namespace ngfem
{


  /*
    HCurl Finite Element Definitions
  */


  template <int D> class HDivFiniteElement;

  
  constexpr int DIM_CURL_ (int D) { return (D*(D-1))/2; }

  /**
     H(Curl) finite element of dimension D
  */

  
  class BaseHCurlFiniteElement : public FiniteElement
  {
  public:
    using FiniteElement :: FiniteElement;
    
    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceMatrix<> shape) const = 0;

    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & mip,
				  BareSliceMatrix<> shape) const = 0;

    virtual void CalcMappedShape (const BaseMappedIntegrationRule & bmir, BareSliceMatrix<> shapes) const = 0;

    virtual void CalcMappedShape (const SIMD<BaseMappedIntegrationPoint> & bmip,
				  BareSliceMatrix<SIMD<double>> shape) const = 0;
    
    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                  BareSliceMatrix<SIMD<double>> shapes) const = 0;
    
    /// compute curl of shape, default: numerical diff
    virtual void CalcCurlShape (const IntegrationPoint & ip, 
				BareSliceMatrix<> curlshape) const = 0;
    
    virtual void CalcMappedCurlShape (const BaseMappedIntegrationPoint & mip,
				      BareSliceMatrix<> curlshape) const = 0;

    virtual void CalcMappedCurlShape (const BaseMappedIntegrationRule & mir, 
                                      BareSliceMatrix<> curlshape) const = 0;

    virtual void CalcMappedCurlShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                      BareSliceMatrix<SIMD<double>> curlshapes) const = 0;



    NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                                          BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const;
    NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                                          BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const;
    NGS_DLL_HEADER virtual void EvaluateCurl (const SIMD_BaseMappedIntegrationRule & ir,
                                              BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const;
    
    NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir,
                                          BareSliceMatrix<SIMD<double>> values,
                                          BareSliceVector<> coefs) const;
    NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir,
                                          BareSliceMatrix<SIMD<Complex>> values,
                                          BareSliceVector<Complex> coefs) const;
    NGS_DLL_HEADER virtual void AddCurlTrans (const SIMD_BaseMappedIntegrationRule & ir,
                                              BareSliceMatrix<SIMD<double>> values,
                                              BareSliceVector<> coefs) const;
    NGS_DLL_HEADER virtual void AddCurlTrans (const SIMD_BaseMappedIntegrationRule & ir,
                                              BareSliceMatrix<SIMD<Complex>> values,
                                              BareSliceVector<Complex> coefs) const;
  };

  
  template <int D>
  class HCurlFiniteElement : public BaseHCurlFiniteElement
  {

  public:
    enum { DIM = D };
    enum { DIM_CURL = DIM_CURL_(D) };


  public:
    using BaseHCurlFiniteElement::BaseHCurlFiniteElement;

    HD virtual ~HCurlFiniteElement () { ; }

    virtual string NGS_DLL_HEADER ClassName() const override;

  
    /// compute curl of shape, default: numerical diff
    void NGS_DLL_HEADER CalcCurlShape (const IntegrationPoint & ip,
                        BareSliceMatrix<> curlshape) const override;
    
    /// compute shape
    void NGS_DLL_HEADER CalcMappedShape (const BaseMappedIntegrationPoint & mip, BareSliceMatrix<> shape) const override;


    void NGS_DLL_HEADER CalcMappedShape (const BaseMappedIntegrationRule & bmir, BareSliceMatrix<> shapes) const override;

    void NGS_DLL_HEADER CalcMappedShape (const SIMD<BaseMappedIntegrationPoint> & bmip,
                          BareSliceMatrix<SIMD<double>> shape) const override;
    
    void NGS_DLL_HEADER CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir,
                          BareSliceMatrix<SIMD<double>> shapes) const override;
    
    /// compute curl of shape
    void NGS_DLL_HEADER CalcMappedCurlShape (const BaseMappedIntegrationPoint & mip,
                              BareSliceMatrix<> curlshape) const override;
    
    void NGS_DLL_HEADER CalcMappedCurlShape (const BaseMappedIntegrationRule & mir,
                              BareSliceMatrix<> curlshape) const override;
    
    void NGS_DLL_HEADER CalcMappedCurlShape (const SIMD_BaseMappedIntegrationRule & mir,
                              BareSliceMatrix<SIMD<double>> curlshapes) const override;

    ///
    const FlatMatrixFixWidth<DIM> GetShape (const IntegrationPoint & ip, 
					    LocalHeap & lh) const
    {
      FlatMatrixFixWidth<DIM> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }


    virtual Vec<D>
    EvaluateShape (const IntegrationPoint & ip, 
                   BareSliceVector<double> x, LocalHeap & lh) const
    {
      HeapReset hr(lh);
      return Trans (GetShape(ip, lh)) * x;
    }  

    using BaseHCurlFiniteElement::Evaluate;
    NGS_DLL_HEADER virtual void 
    Evaluate (const IntegrationRule & ir, BareSliceVector<> coefs, SliceMatrix<> values) const;

    NGS_DLL_HEADER virtual void 
    Evaluate (const MappedIntegrationRule<D,D> & mir, BareSliceVector<> coefs, SliceMatrix<> values) const;
    
    const FlatMatrixFixWidth<DIM_CURL_(D)> GetCurlShape (const IntegrationPoint & ip, 
                                                         LocalHeap & lh) const
    {
      FlatMatrixFixWidth<DIM_CURL_(D)> curlshape(ndof, lh);
      CalcCurlShape (ip, curlshape);
      return curlshape;
    }  
    
    template <typename TVX>
    Vec<DIM_CURL_(D), typename TVX::TSCAL> 
    EvaluateCurlShape (const IntegrationPoint & ip, 
		       const TVX & x, LocalHeap & lh) const
    {
      HeapReset hr(lh);
      return Trans (GetCurlShape(ip, lh)) * x;
    } 

    virtual Vec<DIM_CURL_(D)>
    EvaluateCurlShape (const IntegrationPoint & ip, 
		       BareSliceVector<double> x, LocalHeap & lh) const
    {
      HeapReset hr(lh);
      return Trans (GetCurlShape(ip, lh)) * x;
    }  

    using BaseHCurlFiniteElement::EvaluateCurl;    
    
    NGS_DLL_HEADER virtual void 
    EvaluateCurl (const IntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<> curl) const;

    NGS_DLL_HEADER virtual void 
    EvaluateMappedCurl (const MappedIntegrationRule<D,D> & mir, 
                        BareSliceVector<> coefs, FlatMatrixFixWidth<DIM_CURL_(D)> curl) const;


    
    NGS_DLL_HEADER virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const;

    NGS_DLL_HEADER virtual void CalcDualShape (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> shape) const;

    NGS_DLL_HEADER virtual void EvaluateDual (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
    { throw ExceptionNOSIMD("HCurlFE - simd evaldual not overloaded"); }      

    NGS_DLL_HEADER virtual void AddDualTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
                                              BareSliceVector<double> coefs) const
    { throw ExceptionNOSIMD("HCurlFE - simd adddualtrans not overloaded"); }
    
  protected:
    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<D> shape) const
    { ; }
    ///
    virtual void CalcShape2 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<D> shape) const
    { ; }
  
    virtual void CalcShape3 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<D> shape) const
    { ; }
  
    virtual void CalcShape4 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<D> shape) const
    { ; }
    
    ///
    void ComputeEdgeMoments (int enr, ScalarFiniteElement<1> & testfe,
			     FlatMatrix<> moments, int order, int shape = 1) const;
    ///
    void ComputeFaceMoments (int fnr, HDivFiniteElement<2> & testfe,
			     FlatMatrix<> moments, int order, int shape = 1) const;
    ///
    void ComputeVolMoments (HDivFiniteElement<3> & testfe,
			    FlatMatrix<> moments, int order, int shape = 1) const;

    
    NGS_DLL_HEADER virtual std::list<std::tuple<std::string,double>> Timing () const override;
  };







  


  template <int D>
  extern void ComputeGradientMatrix (const ScalarFiniteElement<D> & h1fe,
				     const HCurlFiniteElement<D> & hcurlfe,
				     FlatMatrix<> gradient);


  
  
  extern template class HCurlFiniteElement<1>;
  extern template class HCurlFiniteElement<2>;
  extern template class HCurlFiniteElement<3>;
}



#endif
