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


  template <int D>
  class DIM_CURL_TRAIT
  {
  public:
    enum { DIM = (D*(D-1))/2 };
  };

  constexpr int DIM_CURL_ (int D) { return (D*(D-1))/2; }

  /**
     H(Curl) finite element of dimension D
  */
  template <int D>
  class HCurlFiniteElement : public FiniteElement
  {

  public:
    enum { DIM = D };
    enum { DIM_CURL = DIM_CURL_TRAIT<D>::DIM };
    // enum { DIM_CURL = DIM_CURL_(D) };


  public:
    ///
    INLINE HCurlFiniteElement () { ; }

    /// 
    INLINE HCurlFiniteElement (int andof, int aorder)
      : FiniteElement (andof, aorder) { ; } 
  
    HD virtual ~HCurlFiniteElement () { ; }

    virtual string ClassName() const;

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
			    SliceMatrix<> shape) const = 0;
  
    /// compute curl of shape, default: numerical diff
    virtual void CalcCurlShape (const IntegrationPoint & ip, 
				SliceMatrix<> curlshape) const;

    /// compute shape
    virtual void CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
				  SliceMatrix<> shape) const;

    virtual void CalcMappedShape (const MappedIntegrationRule<DIM,DIM> & mir, 
                                  SliceMatrix<> shape) const;

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                  BareSliceMatrix<SIMD<double>> shapes) const;
    
    /// compute curl of shape
    virtual void CalcMappedCurlShape (const MappedIntegrationPoint<DIM,DIM> & mip,
				      SliceMatrix<> curlshape) const;

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
    const FlatMatrixFixWidth<DIM_CURL_TRAIT<D>::DIM> GetCurlShape (const IntegrationPoint & ip, 
                                                         LocalHeap & lh) const
    {
      FlatMatrixFixWidth<DIM_CURL_TRAIT<D>::DIM> curlshape(ndof, lh);
      CalcCurlShape (ip, curlshape);
      return curlshape;
    }  
    
    template <typename TVX>
    Vec<DIM_CURL_TRAIT<D>::DIM, typename TVX::TSCAL> 
    EvaluateCurlShape (const IntegrationPoint & ip, 
		       const TVX & x, LocalHeap & lh) const
    {
      HeapReset hr(lh);
      return Trans (GetCurlShape(ip, lh)) * x;
    } 

    virtual Vec<DIM_CURL_TRAIT<D>::DIM> 
    EvaluateCurlShape (const IntegrationPoint & ip, 
		       FlatVector<double> x, LocalHeap & lh) const
    {
      HeapReset hr(lh);
      return Trans (GetCurlShape(ip, lh)) * x;
    }  

    NGS_DLL_HEADER virtual void 
    EvaluateCurl (const IntegrationRule & ir, FlatVector<> coefs, FlatMatrixFixWidth<DIM_CURL_TRAIT<D>::DIM> curl) const;

    NGS_DLL_HEADER virtual void 
    EvaluateMappedCurl (const MappedIntegrationRule<D,D> & mir, 
                        FlatVector<> coefs, FlatMatrixFixWidth<DIM_CURL_TRAIT<D>::DIM> curl) const;


    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
    { throw ExceptionNOSIMD("HCurlFE - simd eval not overloaded"); }
    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const
    { throw ExceptionNOSIMD("HCurlFE - simd<complex> eval not overloaded"); }
    HD NGS_DLL_HEADER virtual void EvaluateCurl (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
    { throw ExceptionNOSIMD("HCurlFE - simd evalcurl not overloaded"); }      

    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                             BareSliceVector<> coefs) const
    { throw ExceptionNOSIMD("HCurlFE - simd addtrans not overloaded"); }
    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values,
                                             BareSliceVector<Complex> coefs) const
    { throw ExceptionNOSIMD("HCurlFE - simd addtrans complex not overloaded"); }
    HD NGS_DLL_HEADER virtual void AddCurlTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                                 BareSliceVector<> coefs) const
    { throw ExceptionNOSIMD("HCurlFE - simd addcurltrans not overloaded"); }

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






  template <int DIM, typename SCAL = double>
  class Class_Du
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM,SCAL> u;

    Class_Du (const AutoDiff<DIM,SCAL> au)
      : u(au) { ; }

    Vec<DIM,SCAL> Value () const
    {
      Vec<DIM,SCAL> val;
      for (int j = 0; j < DIM; j++)
	val(j) = u.DValue(j);
      return val;
    }

    Vec<DIM_CURL,SCAL> CurlValue () const
    {
      return Vec<DIM_CURL,SCAL> (0.0);
    }
  };

  template <int DIM, typename SCAL>
  INLINE Class_Du<DIM,SCAL> Du (AutoDiff<DIM,SCAL> u)
  { return Class_Du<DIM,SCAL> (u); }




  template <int DIM, typename SCAL = double>
  class Class_uDv
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM,SCAL> u, v;
    
    Class_uDv (const AutoDiff<DIM,SCAL> au, 
               const AutoDiff<DIM,SCAL> av)
      : u(au), v(av) { ; }
    
    Vec<DIM,SCAL> Value () const
    {
      Vec<DIM,SCAL> val;
      for (int j = 0; j < DIM; j++)
	val(j) = u.Value() * v.DValue(j);
      return val;
    }

    Vec<DIM_CURL,SCAL> CurlValue () const
    {
      AutoDiff<DIM_CURL,SCAL> hd = Cross (u, v);
      Vec<DIM_CURL,SCAL> val;
      for (int i = 0; i < DIM_CURL; i++) 
        val(i) = hd.DValue(i);
      return val;
    }
  };

  template <int DIM, typename SCAL>
  INLINE Class_uDv<DIM,SCAL> uDv (AutoDiff<DIM,SCAL> u, AutoDiff<DIM,SCAL> v)
  { return Class_uDv<DIM,SCAL> (u,v); }


  template <int DIM, typename SCAL = double>
  class Class_uDv_minus_vDu
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM, SCAL> u, v;

    Class_uDv_minus_vDu (const AutoDiff<DIM,SCAL> au, 
                         const AutoDiff<DIM,SCAL> av)
      : u(au), v(av) { ; }

    Vec<DIM,SCAL> Value () const
    {
      Vec<DIM,SCAL> val;
      for (int j = 0; j < DIM; j++)
	val(j) = u.Value() * v.DValue(j) - v.Value() * u.DValue(j);
      return val;
    }

    Vec<DIM_CURL,SCAL> CurlValue () const
    {
      AutoDiff<DIM_CURL,SCAL> hd = Cross (u, v);
      Vec<DIM_CURL,SCAL> val;
      for (int i = 0; i < DIM_CURL; i++) 
        val(i) = 2 * hd.DValue(i);
      return val;
    }
  };

  template <int DIM, typename SCAL>
  INLINE Class_uDv_minus_vDu<DIM,SCAL> 
  uDv_minus_vDu (AutoDiff<DIM,SCAL> u, AutoDiff<DIM,SCAL> v)
  { return Class_uDv_minus_vDu<DIM,SCAL> (u,v); }





  template <int DIM, typename SCAL = double>
  class Class_wuDv_minus_wvDu
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM,SCAL> u, v, w;

    Class_wuDv_minus_wvDu (const AutoDiff<DIM,SCAL> au, 
                           const AutoDiff<DIM,SCAL> av,
                           const AutoDiff<DIM,SCAL> aw)
      : u(au), v(av), w(aw) { ; }
    
    Vec<DIM,SCAL> Value () const
    {
      Vec<DIM,SCAL> val;
      for (int j = 0; j < DIM; j++)
	val(j) = w.Value() * (u.Value() * v.DValue(j) - v.Value() * u.DValue(j));
      return val;
    }

    Vec<DIM_CURL,SCAL> CurlValue () const
    {
      AutoDiff<DIM_CURL,SCAL> hd = Cross (u*w, v) + Cross(u, v*w);
      Vec<DIM_CURL,SCAL> val;
      for (int i = 0; i < DIM_CURL; i++) 
        val(i) = hd.DValue(i);
      return val;
    }
  };


  template <int DIM, typename SCAL>
  INLINE Class_wuDv_minus_wvDu<DIM,SCAL> 
  wuDv_minus_wvDu (AutoDiff<DIM,SCAL> u, AutoDiff<DIM,SCAL> v, AutoDiff<DIM,SCAL> w)
  { return Class_wuDv_minus_wvDu<DIM,SCAL> (u,v,w); }







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
