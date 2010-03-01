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


  template <> class DIM_CURL_TRAIT<1>
  {
  public:
    enum { DIM = 1 };  // should be 0; changed to make gcc34 feel ok
  };


  /**
     H(Curl) finite element of dimension D
  */
  template <int D>
  class NGS_DLL_HEADER HCurlFiniteElement : public FiniteElement
  {

  public:
    enum { DIM = D };
    enum { DIM_CURL = DIM_CURL_TRAIT<D>::DIM };


  public:
    ///
    HCurlFiniteElement () { ; }

    /// 
    HCurlFiniteElement (ELEMENT_TYPE aeltype, int andof, int aorder)
      : FiniteElement (DIM, aeltype, andof, aorder) { ; } 
  
    virtual ~HCurlFiniteElement () { ; }

    virtual string ClassName(void) const;

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<DIM> shape) const = 0;
  
    /// compute curl of shape, default: numerical diff
    virtual void CalcCurlShape (const IntegrationPoint & ip, 
				FlatMatrixFixWidth<DIM_CURL> curlshape) const;

    /// compute shape
    virtual void CalcMappedShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
				  FlatMatrixFixWidth<DIM> shape) const;

    /// compute curl of shape
    virtual void CalcMappedCurlShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
				      FlatMatrixFixWidth<DIM_CURL> curlshape) const;


    ///
    const FlatMatrixFixWidth<DIM> GetShape (const IntegrationPoint & ip, 
					    LocalHeap & lh) const
    {
      FlatMatrixFixWidth<DIM> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

    ///
    const FlatMatrixFixWidth<DIM_CURL> GetCurlShape (const IntegrationPoint & ip, 
						     LocalHeap & lh) const
    {
      FlatMatrixFixWidth<DIM_CURL> curlshape(ndof, lh);
      CalcCurlShape (ip, curlshape);
      return curlshape;
    }  
    
    template <typename TVX>
    Vec<DIM_CURL, typename TVX::TSCAL> 
    EvaluateCurlShape (const IntegrationPoint & ip, 
		       const TVX & x, LocalHeap & lh) const
    {
      return Trans (GetCurlShape(ip, lh)) * x;
    } 

    virtual Vec<DIM_CURL> 
    EvaluateCurlShape (const IntegrationPoint & ip, 
		       FlatVector<double> x, LocalHeap & lh) const
    {
      return Trans (GetCurlShape(ip, lh)) * x;
    }  


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








  template <int DIM>
  class Du
  {
  public:
    const AutoDiff<DIM> & u;
    Du (const AutoDiff<DIM> & au)
      : u(au) { ; }
  };

  template <int DIM>
  class uDv
  {
  public:
    const AutoDiff<DIM> & u, v;
    uDv (const AutoDiff<DIM> & au, 
         const AutoDiff<DIM> & av)
      : u(au), v(av) { ; }
  };


  template <int DIM>
  class uDv_minus_vDu
  {
  public:
    const AutoDiff<DIM> & u, v;
    uDv_minus_vDu (const AutoDiff<DIM> & au, 
                   const AutoDiff<DIM> & av)
      : u(au), v(av) { ; }
  };

  template <int DIM>
  class wuDv_minus_wvDu
  {
  public:
    const AutoDiff<DIM> & u, v, w;
    wuDv_minus_wvDu (const AutoDiff<DIM> & au, 
                     const AutoDiff<DIM> & av,
                     const AutoDiff<DIM> & aw)
      : u(au), v(av), w(aw) { ; }
  };




  template <int DIM>
  class HCurlShapeElement
  {
    double * data;
  public:
    HCurlShapeElement (double * adata) : data(adata) { ; }

    void operator= (const Du<DIM> & uv) 
    { for (int i = 0; i < DIM; i++) 
        data[i] = uv.u.DValue(i); }

    void operator= (const uDv<DIM> & uv) 
    { for (int i = 0; i < DIM; i++) 
        data[i] = uv.u.Value() * uv.v.DValue(i); }

    void operator= (const uDv_minus_vDu<DIM> & uv) 
    { for (int i = 0; i < DIM; i++) 
        data[i] = 
          uv.u.Value() * uv.v.DValue(i)
          -uv.v.Value() * uv.u.DValue(i); }

    void operator= (const wuDv_minus_wvDu<DIM> & uv) 
    { for (int i = 0; i < DIM; i++) 
        data[i] = 
          uv.w.Value() * uv.u.Value() * uv.v.DValue(i)
          -uv.w.Value() * uv.v.Value() * uv.u.DValue(i); }
  };

  
  // hv.DValue() = (grad u) x (grad v) 
  inline AutoDiff<3> Cross (const AutoDiff<3> & u,
			    const AutoDiff<3> & v)
  {
    AutoDiff<3> hv;
    hv.Value() = 0.0;
    hv.DValue(0) = u.DValue(1)*v.DValue(2)-u.DValue(2)*v.DValue(1);
    hv.DValue(1) = u.DValue(2)*v.DValue(0)-u.DValue(0)*v.DValue(2);
    hv.DValue(2) = u.DValue(0)*v.DValue(1)-u.DValue(1)*v.DValue(0);
    return hv;
  }

  inline AutoDiff<1> Cross (const AutoDiff<2> & u,
			    const AutoDiff<2> & v)
  {
    AutoDiff<1> hv;
    hv.Value() = 0.0;
    hv.DValue(0) = u.DValue(0)*v.DValue(1)-u.DValue(1)*v.DValue(0);
    return hv;
  }



  template <int DIM>
  class HCurlCurlShapeElement
  {
    double * data;
    enum { DIM_CURL = (DIM * (DIM-1))/2 };
  public:
    HCurlCurlShapeElement (double * adata) : data(adata) { ; }

    void operator= (const Du<DIM> & uv) 
    { for (int i = 0; i < DIM; i++) 
        data[i] = 0; }

    void operator= (const uDv<DIM> & uv) 
    {
      AutoDiff<DIM_CURL> hd = Cross (uv.u, uv.v);
      for (int i = 0; i < DIM_CURL; i++) 
        data[i] = hd.DValue(i);
    }

    void operator= (const uDv_minus_vDu<DIM> & uv) 
    {
      AutoDiff<DIM_CURL> hd = Cross (uv.u, uv.v);
      for (int i = 0; i < DIM_CURL; i++) 
        data[i] = 2*hd.DValue(i);
    }

    void operator= (const wuDv_minus_wvDu<DIM> & uv) 
    {
      AutoDiff<DIM_CURL> hd = Cross (uv.u*uv.w, uv.v) + Cross(uv.u, uv.v*uv.w);
      for (int i = 0; i < DIM_CURL; i++) 
        data[i] = hd.DValue(i);
    }
  };

  template <int DIM>
  class HCurlEvaluateCurlElement
  {
    const double * coefs;
    enum { DIM_CURL = (DIM * (DIM-1))/2 };
    Vec<DIM_CURL> & sum;
  public:
    HCurlEvaluateCurlElement (const double * acoefs, Vec<DIM_CURL> & asum)
      : coefs(acoefs), sum(asum) { ; }

    void operator= (const Du<DIM> & uv) 
    { ; }

    void operator= (const uDv<DIM> & uv) 
    {
      AutoDiff<DIM_CURL> hd = Cross (uv.u, uv.v);
      for (int i = 0; i < DIM_CURL; i++) 
        sum[i] += *coefs * hd.DValue(i);
    }

    void operator= (const uDv_minus_vDu<DIM> & uv) 
    {
      AutoDiff<DIM_CURL> hd = Cross (uv.u, uv.v);
      for (int i = 0; i < DIM_CURL; i++) 
        sum[i] += 2 * *coefs * hd.DValue(i);
    }

    void operator= (const wuDv_minus_wvDu<DIM> & uv) 
    {
      AutoDiff<DIM_CURL> hd = Cross (uv.u*uv.w, uv.v) + Cross(uv.u, uv.v*uv.w);
      for (int i = 0; i < DIM_CURL; i++) 
        sum[i] += *coefs * hd.DValue(i);
    }
  };



  template <int DIM>
  class HCurlShapeAssign
  {
    double * dshape;
  public:
    HCurlShapeAssign (FlatMatrixFixWidth<DIM> mat)
    { dshape = &mat(0,0); }

    HCurlShapeElement<DIM> operator[] (int i) const
    { return HCurlShapeElement<DIM> (dshape + i*DIM); }
  };


  template <int DIM>
  class HCurlCurlShapeAssign
  {
    double * dshape;
    enum { DIM_CURL = (DIM * (DIM-1))/2 };
  public:
    HCurlCurlShapeAssign (FlatMatrixFixWidth<DIM_CURL> mat)
    { dshape = &mat(0,0); }

    HCurlCurlShapeElement<DIM> operator[] (int i) const
    { return HCurlCurlShapeElement<DIM> (dshape + i*DIM_CURL); }
  };

  template <int DIM>
  class HCurlEvaluateCurl
  {
    const double * coefs;
    enum { DIM_CURL = (DIM * (DIM-1))/2 };
    Vec<DIM_CURL> sum;
  public:
    HCurlEvaluateCurl (FlatVector<> acoefs)
    { coefs = &acoefs(0); sum = 0.0; }

    HCurlEvaluateCurlElement<DIM> operator[] (int i) 
    { return HCurlEvaluateCurlElement<DIM> (coefs+i, sum); }

    Vec<DIM_CURL> Sum() { return sum; }
  };







  /**
     Base-element for template polymorphism.
     Barton and Nackman Trick
  */
  
  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  class T_HCurlFiniteElement : public HCurlFiniteElement<ET_trait<ET>::DIM>
  {

  public:
    
  protected:
    enum { DIM = ET_trait<ET>::DIM };
    using HCurlFiniteElement<DIM>::DIM_CURL;

    T_HCurlFiniteElement ()
      : HCurlFiniteElement<DIM> (ET, NDOF, ORDER) { ; }

    virtual ~T_HCurlFiniteElement() { ; }

  public:

    virtual void CalcShape (const IntegrationPoint & ip, 
                            FlatMatrixFixWidth<DIM> shape) const
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);
      
      HCurlShapeAssign<DIM> ds(shape); 
      FEL::T_CalcShape (adp, ds);
    }

    virtual void
    CalcMappedShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
                     FlatMatrixFixWidth<DIM> shape) const
    {
      AutoDiff<DIM> adp[DIM];
      
      for (int i = 0; i < DIM; i++)
        adp[i].Value() = sip.IP()(i);
      
      for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
          adp[i].DValue(j) = sip.GetJacobianInverse()(i,j);
      
      HCurlShapeAssign<DIM> ds(shape); 
      FEL::T_CalcShape (adp, ds);
    }



    virtual void CalcCurlShape (const IntegrationPoint & ip, 
                                FlatMatrixFixWidth<DIM_CURL> curlshape) const
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);

      HCurlCurlShapeAssign<DIM> ds(curlshape); 
      FEL::T_CalcShape (adp, ds);
    }

    virtual void
    CalcMappedCurlShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
                         FlatMatrixFixWidth<DIM_CURL> curlshape) const
    {
      AutoDiff<DIM> adp[DIM];

      for (int i = 0; i < DIM; i++)
        adp[i].Value() = sip.IP()(i);
      
      for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
          adp[i].DValue(j) = sip.GetJacobianInverse()(i,j);

      HCurlCurlShapeAssign<DIM> ds(curlshape); 
      FEL::T_CalcShape (adp, ds);
    }

    virtual Vec <DIM_CURL>
    EvaluateCurlShape (const IntegrationPoint & ip, 
                       FlatVector<double> x,
                       LocalHeap & lh) const
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);
      
      HCurlEvaluateCurl<DIM> ds(x); 
      FEL::T_CalcShape (adp, ds);
      return ds.Sum();
    }

  };








  template <int D>
  extern void ComputeGradientMatrix (const ScalarFiniteElement<D> & h1fe,
				     const HCurlFiniteElement<D> & hcurlfe,
				     FlatMatrix<> gradient);






}




#endif
