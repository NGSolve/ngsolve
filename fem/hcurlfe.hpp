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
    INLINE BaseHCurlFiniteElement () { ; } 
    INLINE BaseHCurlFiniteElement (int andof, int aorder)
      : FiniteElement (andof, aorder) { ; }

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
    ///
    INLINE HCurlFiniteElement () { ; }

    /// 
    INLINE HCurlFiniteElement (int andof, int aorder)
      : BaseHCurlFiniteElement (andof, aorder) { ; } 
  
    HD virtual ~HCurlFiniteElement () { ; }

    virtual string ClassName() const override;

  
    /// compute curl of shape, default: numerical diff
    void CalcCurlShape (const IntegrationPoint & ip, 
                        BareSliceMatrix<> curlshape) const override;
    
    /// compute shape
    void CalcMappedShape (const BaseMappedIntegrationPoint & mip, BareSliceMatrix<> shape) const override;


    void CalcMappedShape (const BaseMappedIntegrationRule & bmir, BareSliceMatrix<> shapes) const override;

    void CalcMappedShape (const SIMD<BaseMappedIntegrationPoint> & bmip,
                          BareSliceMatrix<SIMD<double>> shape) const override;
    
    void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                          BareSliceMatrix<SIMD<double>> shapes) const override;
    
    /// compute curl of shape
    void CalcMappedCurlShape (const BaseMappedIntegrationPoint & mip,
                              BareSliceMatrix<> curlshape) const override;
    
    void CalcMappedCurlShape (const BaseMappedIntegrationRule & mir, 
                              BareSliceMatrix<> curlshape) const override;
    
    void CalcMappedCurlShape (const SIMD_BaseMappedIntegrationRule & mir, 
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













  // hv.DValue() = (grad u) x (grad v) 
  template <typename SCAL>
  INLINE AutoDiff<3,SCAL> Cross (const AutoDiff<3,SCAL> & u,
                                 const AutoDiff<3,SCAL> & v)
  {
    SCAL hv[3];
    hv[0] = u.DValue(1)*v.DValue(2)-u.DValue(2)*v.DValue(1);
    hv[1] = u.DValue(2)*v.DValue(0)-u.DValue(0)*v.DValue(2);
    hv[2] = u.DValue(0)*v.DValue(1)-u.DValue(1)*v.DValue(0);
    return AutoDiff<3,SCAL> (0, hv);
  }

  template <typename SCAL>
  INLINE AutoDiff<1,SCAL> Cross (const AutoDiff<2,SCAL> & u,
                                 const AutoDiff<2,SCAL> & v)
  {
    AutoDiff<1,SCAL> hv;
    hv.Value() = 0.0;
    hv.DValue(0) = u.DValue(0)*v.DValue(1)-u.DValue(1)*v.DValue(0);
    return hv;
  }

  template <typename SCAL>
  INLINE AutoDiff<0,SCAL> Cross (const AutoDiff<1,SCAL> & u,
                                 const AutoDiff<1,SCAL> & v)
  {
    AutoDiff<0,SCAL> hv;
    hv.Value() = 0.0;
    return hv;
  }




  template <int DIM, typename SCAL> 
  class Du
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM,SCAL> u;

    Du (const AutoDiff<DIM,SCAL> au) : u(au) { }

    Vec<DIM,SCAL> Value () const
    {
      return GetGradient(u);
    }

    Vec<DIM_CURL,SCAL> CurlValue () const
    {
      return Vec<DIM_CURL,SCAL> (0.0);
    }
  };



  template <int DIM, typename SCAL>
  class uDv
  {
  public:
    const AutoDiff<DIM,SCAL> u, v;

    uDv (AutoDiff<DIM,SCAL> au, AutoDiff<DIM,SCAL> av)
      : u(au), v(av) { ; }

    auto Value () const
    {
      return u.Value() * GetGradient(v);
    }

    auto CurlValue () const
    {
      return Cross (GetGradient(u), GetGradient(v));
    }
  };



  template <int DIM, typename SCAL>
  class uDv_minus_vDu
  {
  public:
    const AutoDiff<DIM, SCAL> u, v;
    
    uDv_minus_vDu (const AutoDiff<DIM,SCAL> au, 
                   const AutoDiff<DIM,SCAL> av)
      : u(au), v(av) { }

    auto Value () const
    {
      return u.Value()*GetGradient(v)-v.Value()*GetGradient(u);
    }

    auto CurlValue () const    
    {
      return 2 * Cross (GetGradient(u), GetGradient(v));
    }
  };




  template <int DIM, typename SCAL>
  class wuDv_minus_wvDu
  {
  public:
    const AutoDiff<DIM,SCAL> u, v, w;

    wuDv_minus_wvDu (const AutoDiff<DIM,SCAL> au, 
                     const AutoDiff<DIM,SCAL> av,
                     const AutoDiff<DIM,SCAL> aw)
      : u(au), v(av), w(aw) { ; }

    auto Value () const
    {
      return w.Value()*u.Value()*GetGradient(v) - w.Value()*v.Value()*GetGradient(u);
    }

    auto CurlValue () const
    {
      return Cross(GetGradient(u*w),GetGradient(v)) - Cross(GetGradient(v*w), GetGradient(u));
    }
  };



  template <int D>
  extern void ComputeGradientMatrix (const ScalarFiniteElement<D> & h1fe,
				     const HCurlFiniteElement<D> & hcurlfe,
				     FlatMatrix<> gradient);


  

#ifdef FILE_HCURLFE_CPP
#define HCURLFE_EXTERN
#else
#define HCURLFE_EXTERN extern
  
  HCURLFE_EXTERN template class HCurlFiniteElement<1>;
  HCURLFE_EXTERN template class HCurlFiniteElement<2>;
  HCURLFE_EXTERN template class HCurlFiniteElement<3>;

#endif

}




#endif
