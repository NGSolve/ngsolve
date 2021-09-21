#ifndef FILE_HCURLFE
#define FILE_HCURLFE

/*********************************************************************/
/* File:   hcurlfe.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   16. Apr. 2000                                             */
/*********************************************************************/



namespace ngfem
{


  /*
    HCurl Finite Element Definitions
  */


  template <int D> class HDivFiniteElement;

  /*
  template <int D>
  class DIM_CURL_TRAIT
  {
  public:
    enum { DIM = (D*(D-1))/2 };
  };
  */
  
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
			    SliceMatrix<> shape) const = 0;

    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & mip,
				  SliceMatrix<> shape) const = 0;

    virtual void CalcMappedCurlShape (const BaseMappedIntegrationPoint & mip,
				      SliceMatrix<> curlshape) const = 0;
  };

  
  template <int D>
  class HCurlFiniteElement : public BaseHCurlFiniteElement
  {

  public:
    enum { DIM = D };
    // enum { DIM_CURL = DIM_CURL_TRAIT<D>::DIM };
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
    virtual void CalcCurlShape (const IntegrationPoint & ip, 
				SliceMatrix<> curlshape) const;

    /// compute shape
    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & mip,
				  SliceMatrix<> shape) const override;


    virtual void CalcMappedShape (const BaseMappedIntegrationRule & bmir, SliceMatrix<> shapes) const;

    virtual void CalcMappedShape (const SIMD<BaseMappedIntegrationPoint> & bmip,
				  BareSliceMatrix<SIMD<double>> shape) const;
    
    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                  BareSliceMatrix<SIMD<double>> shapes) const;
    
    /// compute curl of shape
    virtual void CalcMappedCurlShape (const BaseMappedIntegrationPoint & mip,
				      SliceMatrix<> curlshape) const override;

    virtual void CalcMappedCurlShape (const MappedIntegrationRule<DIM,DIM> & mir, 
                                      SliceMatrix<> curlshape) const;

    virtual void CalcMappedCurlShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                      BareSliceMatrix<SIMD<double>> curlshapes) const;

    ///
    const FlatMatrixFixWidth<DIM> GetShape (const IntegrationPoint & ip, 
					    LocalHeap & lh) const
    {
      FlatMatrixFixWidth<DIM> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

    ///
    /*
    const FlatMatrixFixWidth<DIM_CURL_TRAIT<D>::DIM> GetCurlShape (const IntegrationPoint & ip, 
                                                         LocalHeap & lh) const
    {
      FlatMatrixFixWidth<DIM_CURL_TRAIT<D>::DIM> curlshape(ndof, lh);
      CalcCurlShape (ip, curlshape);
      return curlshape;
    }  
    */
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

    // virtual Vec<DIM_CURL_TRAIT<D>::DIM>
    virtual Vec<DIM_CURL_(D)>
    EvaluateCurlShape (const IntegrationPoint & ip, 
		       BareSliceVector<double> x, LocalHeap & lh) const
    {
      HeapReset hr(lh);
      return Trans (GetCurlShape(ip, lh)) * x;
    }  

    NGS_DLL_HEADER virtual void 
    EvaluateCurl (const IntegrationRule & ir, BareSliceVector<> coefs, FlatMatrixFixWidth<DIM_CURL_(D)> curl) const;

    NGS_DLL_HEADER virtual void 
    EvaluateMappedCurl (const MappedIntegrationRule<D,D> & mir, 
                        BareSliceVector<> coefs, FlatMatrixFixWidth<DIM_CURL_(D)> curl) const;


    NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
    { throw ExceptionNOSIMD(string("HCurlFE - simd eval not overloaded, eltype = ")+typeid(*this).name()); }
    NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const
    { throw ExceptionNOSIMD(string("HCurlFE - simd<complex> eval not overloaded")+typeid(*this).name()); }
    NGS_DLL_HEADER virtual void EvaluateCurl (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
    { throw ExceptionNOSIMD(string("HCurlFE - simd evalcurl not overloaded")+typeid(*this).name()); }      

    NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                             BareSliceVector<> coefs) const
    { throw ExceptionNOSIMD(string("HCurlFE - simd addtrans not overloaded")+typeid(*this).name()); }
    NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values,
                                             BareSliceVector<Complex> coefs) const
    { throw ExceptionNOSIMD(string("HCurlFE - simd addtrans complex not overloaded")+typeid(*this).name()); }
    NGS_DLL_HEADER virtual void AddCurlTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                                 BareSliceVector<> coefs) const
    { throw ExceptionNOSIMD(string("HCurlFE - simd addcurltrans not overloaded")+typeid(*this).name()); }
    NGS_DLL_HEADER virtual void AddCurlTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values,
                                              BareSliceVector<Complex> coefs) const
    { throw ExceptionNOSIMD(string("HCurlFE - simd addcurltrans complex not overloaded")+typeid(*this).name()); }      

    
    NGS_DLL_HEADER virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, SliceMatrix<> shape) const;

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




  template <int DIM, typename SCAL> //  = double>
  class Du
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM,SCAL> u;

    Du (const AutoDiff<DIM,SCAL> au) : u(au) { }
    // Du (const AutoDiff<DIM,SCAL> au) : u(au) { }

    Vec<DIM,SCAL> Value () const
    {
      return GetGradient(u);
    }

    Vec<DIM_CURL,SCAL> CurlValue () const
    {
      return Vec<DIM_CURL,SCAL> (0.0);
    }
  };



  template <int DIM, typename SCAL = double>
  class uDv
  {
    // enum { DIM_CURL = (DIM * (DIM-1))/2 };    
  public:
    const AutoDiff<DIM,SCAL> u, v;

    uDv (AutoDiff<DIM,SCAL> au, AutoDiff<DIM,SCAL> av)
      : u(au), v(av) { ; }

    // uDv (AutoDiff<DIM,SCAL> au, AutoDiff<DIM,SCAL> av)
    // : u(au), v(av) { ; }
    
    // Vec<DIM,SCAL> Value () const
    auto Value () const
    {
      return u.Value() * GetGradient(v);
    }

    // Vec<DIM_CURL,SCAL> CurlValue () const
    auto CurlValue () const
    {
      return Cross (GetGradient(u), GetGradient(v));
    }
  };



  template <int DIM, typename SCAL = double>
  class uDv_minus_vDu
  {
    // enum { DIM_CURL = (DIM * (DIM-1))/2 };
  public:
    const AutoDiff<DIM, SCAL> u, v;
    
    uDv_minus_vDu (const AutoDiff<DIM,SCAL> au, 
                   const AutoDiff<DIM,SCAL> av)
      : u(au), v(av) { }

    /*
    uDv_minus_vDu (const AutoDiff<DIM,SCAL> au, 
                   const AutoDiff<DIM,SCAL> av)
      : u(au), v(av) { }
    */
    
    // Vec<DIM,SCAL> Value () const
    auto Value () const
    {
      return u.Value()*GetGradient(v)-v.Value()*GetGradient(u);
    }

    // Vec<DIM_CURL,SCAL> CurlValue () const
    auto CurlValue () const    
    {
      return 2 * Cross (GetGradient(u), GetGradient(v));
    }
  };




  template <int DIM, typename SCAL = double>
  class wuDv_minus_wvDu
  {
    // enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM,SCAL> u, v, w;

    wuDv_minus_wvDu (const AutoDiff<DIM,SCAL> au, 
                     const AutoDiff<DIM,SCAL> av,
                     const AutoDiff<DIM,SCAL> aw)
      : u(au), v(av), w(aw) { ; }
    /*
    wuDv_minus_wvDu (const AutoDiff<DIM,SCAL> au, 
                     const AutoDiff<DIM,SCAL> av,
                     const AutoDiff<DIM,SCAL> aw)
      : u(au), v(av), w(aw) { ; }
    */
    
    // Vec<DIM,SCAL> Value () const
    auto Value () const
    {
      /*
      Vec<DIM,SCAL> val;
      for (int j = 0; j < DIM; j++)
	val(j) = w.Value() * (u.Value() * v.DValue(j) - v.Value() * u.DValue(j));
      return val;
      */
      return w.Value()*u.Value()*GetGradient(v) - w.Value()*v.Value()*GetGradient(u);
    }

    // Vec<DIM_CURL,SCAL> CurlValue () const
    auto CurlValue () const
    {
      /*
      AutoDiff<DIM_CURL,SCAL> hd = Cross (u*w, v) + Cross(u, v*w);
      Vec<DIM_CURL,SCAL> val;
      for (int i = 0; i < DIM_CURL; i++) 
        val(i) = hd.DValue(i);
      return val;
      */
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
