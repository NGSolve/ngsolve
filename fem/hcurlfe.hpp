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
    INLINE HCurlFiniteElement () { ; }

    /// 
    INLINE HCurlFiniteElement (int andof, int aorder)
      : FiniteElement (andof, aorder) { ; } 
  
    virtual ~HCurlFiniteElement () { ; }

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
    
    /// compute curl of shape
    virtual void CalcMappedCurlShape (const MappedIntegrationPoint<DIM,DIM> & mip,
				      SliceMatrix<> curlshape) const;


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

  inline AutoDiff<0> Cross (const AutoDiff<1> & u,
			    const AutoDiff<1> & v)
  {
    AutoDiff<0> hv;
    hv.Value() = 0.0;
    return hv;
  }






  template <int DIM>
  class Du
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM> & u;

    Du (const AutoDiff<DIM> & au)
      : u(au) { ; }

    Vec<DIM> Value () const
    {
      Vec<DIM> val;
      for (int j = 0; j < DIM; j++)
	val(j) = u.DValue(j);
      return val;
    }

    Vec<DIM_CURL> CurlValue () const
    {
      return Vec<DIM> (0.0);
    }
  };




  template <int DIM>
  class uDv
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM> & u, v;

    uDv (const AutoDiff<DIM> & au, 
         const AutoDiff<DIM> & av)
      : u(au), v(av) { ; }

    Vec<DIM> Value () const
    {
      Vec<DIM> val;
      for (int j = 0; j < DIM; j++)
	val(j) = u.Value() * v.DValue(j);
      return val;
    }

    Vec<DIM_CURL> CurlValue () const
    {
      AutoDiff<DIM_CURL> hd = Cross (u, v);
      Vec<DIM_CURL> val;
      for (int i = 0; i < DIM_CURL; i++) 
        val(i) = hd.DValue(i);
      return val;
    }
  };



  template <int DIM>
  class uDv_minus_vDu
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM> & u, v;

    uDv_minus_vDu (const AutoDiff<DIM> & au, 
                   const AutoDiff<DIM> & av)
      : u(au), v(av) { ; }

    Vec<DIM> Value () const
    {
      Vec<DIM> val;
      for (int j = 0; j < DIM; j++)
	val(j) = u.Value() * v.DValue(j) - v.Value() * u.DValue(j);
      return val;
    }

    Vec<DIM_CURL> CurlValue () const
    {
      AutoDiff<DIM_CURL> hd = Cross (u, v);
      Vec<DIM_CURL> val;
      for (int i = 0; i < DIM_CURL; i++) 
        val(i) = 2 * hd.DValue(i);
      return val;
    }
  };


  template <int DIM>
  class wuDv_minus_wvDu
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM> & u, v, w;

    wuDv_minus_wvDu (const AutoDiff<DIM> & au, 
                     const AutoDiff<DIM> & av,
                     const AutoDiff<DIM> & aw)
      : u(au), v(av), w(aw) { ; }

    Vec<DIM> Value () const
    {
      Vec<DIM> val;
      for (int j = 0; j < DIM; j++)
	val(j) = w.Value() * (u.Value() * v.DValue(j) - v.Value() * u.DValue(j));
      return val;
    }

    Vec<DIM_CURL> CurlValue () const
    {
      AutoDiff<DIM_CURL> hd = Cross (u*w, v) + Cross(u, v*w);
      Vec<DIM_CURL> val;
      for (int i = 0; i < DIM_CURL; i++) 
        val(i) = hd.DValue(i);
      return val;
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
#endif
  
  HCURLFE_EXTERN template class HCurlFiniteElement<1>;
  HCURLFE_EXTERN template class HCurlFiniteElement<2>;
  HCURLFE_EXTERN template class HCurlFiniteElement<3>;

}




#endif
