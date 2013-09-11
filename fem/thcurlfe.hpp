#ifndef FILE_THCURLFE
#define FILE_THCURLFE

/*********************************************************************/
/* File:   thcurlfe.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   5. Sep. 2013                                              */
/*********************************************************************/

namespace ngfem
{


  
  /**
     HCurlHighOrderFE of shape ET.
     provides access functions, shape funcitons are provided by CalcShape template
  */
  template <ELEMENT_TYPE ET, typename SHAPES,
            typename BASE = HCurlFiniteElement<ET_trait<ET>::DIM>>
  class T_HCurlHighOrderFiniteElement : public BASE 
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };
  
    using BASE::DIM_CURL;
    using BASE::ndof;
    using BASE::order;

    // using BASE::xx_discontinuous;

  public:

    NGS_DLL_HEADER T_HCurlHighOrderFiniteElement () { ; }
    virtual ELEMENT_TYPE ElementType() const { return ET; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[2], TFA & shape) const
    { 
      static_cast<const SHAPES*> (this) -> T_CalcShape (hx, shape);
    }

    virtual void CalcShape (const IntegrationPoint & ip, 
                            SliceMatrix<> shape) const;

    virtual void CalcCurlShape (const IntegrationPoint & ip, 
                                SliceMatrix<> curlshape) const;

    virtual void CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                                  SliceMatrix<> shape) const;

    virtual void CalcMappedCurlShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                                      SliceMatrix<> curlshape) const;

    /*
      virtual Vec <DIM_CURL_TRAIT<ET_trait<ET>::DIM>::DIM>
      EvaluateCurlShape (const IntegrationPoint & ip, 
      FlatVector<double> x,
      LocalHeap & lh) const;
    */
  };








  template <int DIM>
  class HCurlShapeElement
  {
    FlatVec<DIM> data;
  public:
    HCurlShapeElement (FlatVector<> adata) : data(adata.Data()) { ; }

    template <typename F1>
    void operator= (F1 form1)  { data = form1.Value(); }
  };

  template <int DIM>
  class HCurlCurlShapeElement
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };
    FlatVec<DIM_CURL> data;
    // FlatVector<> data;
  public:
    // HCurlCurlShapeElement (FlatVec<DIM_CURL> adata) : data(adata) { ; }
    HCurlCurlShapeElement (FlatVector<> adata) : data(adata.Data()) { ; }

    template <typename F1>
    void operator= (F1 form1)  { data = form1.CurlValue(); }
  };


  template <int DIM>
  class HCurlEvaluateCurlElement
  {
    const double & coef;
    enum { DIM_CURL = (DIM * (DIM-1))/2 };
    Vec<DIM_CURL> & sum;
  public:
    HCurlEvaluateCurlElement (const double & acoef, Vec<DIM_CURL> & asum)
      : coef(acoef), sum(asum) { ; }

    template <typename F1>
    void operator= (F1 form1)  { sum += coef * form1.CurlValue(); }
  };



  template <int DIM>
  class HCurlShapeAssign
  {
    SliceMatrix<> shape; 
  public:
    HCurlShapeAssign (SliceMatrix<> mat) : shape(mat) { ; }

    HCurlShapeElement<DIM> operator[] (int i) const
    { return HCurlShapeElement<DIM> (shape.Row(i)); }
  };

  template <int DIM>
  class HCurlCurlShapeAssign
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };
    SliceMatrix<> cshape; 
  public:
    HCurlCurlShapeAssign (SliceMatrix<> mat) : cshape(mat) { ; }

    HCurlCurlShapeElement<DIM> operator[] (int i) const 
    { return HCurlCurlShapeElement<DIM> (cshape.Row(i)); }
  };

  template <int DIM>
  class HCurlEvaluateCurl
  {
    FlatVector<> coefs;
    enum { DIM_CURL = (DIM * (DIM-1))/2 };
    Vec<DIM_CURL> sum;
  public:
    HCurlEvaluateCurl (FlatVector<> acoefs) : coefs(acoefs), sum(0.0) { ; }

    HCurlEvaluateCurlElement<DIM> operator[] (int i) 
    { return HCurlEvaluateCurlElement<DIM> (coefs(i), sum); }

    Vec<DIM_CURL> Sum() { return sum; }
  };




  template <int DIM>
  class HCurl_Shape : public Vec<DIM>
  {
  public:
    template <typename T>
    HCurl_Shape (T shape) : Vec<DIM>(shape.Value()) { ; }
    // operator Vec<DIM> () const { return vec; }
  };

  template <int DIM>
  class HCurl_CurlShape : public Vec<DIM_CURL_TRAIT<DIM>::DIM>
  {
  public:
    template <typename T>
    HCurl_CurlShape (T shape) 
      : Vec<DIM_CURL_TRAIT<DIM>::DIM> (shape.CurlValue()) { ; }
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
      : HCurlFiniteElement<DIM> (NDOF, ORDER) { ; }

    virtual ~T_HCurlFiniteElement() { ; }

  public:
    virtual ELEMENT_TYPE ElementType() const { return ET; }
    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            SliceMatrix<> shape) const
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);
      
      HCurlShapeAssign<DIM> ds(shape); 
      FEL::T_CalcShape (adp, ds);
    }

    virtual void
    CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                     SliceMatrix<> shape) const
    {
      AutoDiff<DIM> adp[DIM];
      
      for (int i = 0; i < DIM; i++)
        adp[i].Value() = mip.IP()(i);
      
      for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
          adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);
      
      HCurlShapeAssign<DIM> ds(shape); 
      FEL::T_CalcShape (adp, ds);
    }


    virtual void CalcCurlShape (const IntegrationPoint & ip, 
                                SliceMatrix<> curlshape) const
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);

      HCurlCurlShapeAssign<DIM> ds(curlshape); 
      FEL::T_CalcShape (adp, ds);
    }

    virtual void
    CalcMappedCurlShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                         SliceMatrix<> curlshape) const
    {
      AutoDiff<DIM> adp[DIM];

      for (int i = 0; i < DIM; i++)
        adp[i].Value() = mip.IP()(i);
      
      for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
          adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);

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




}



#endif
