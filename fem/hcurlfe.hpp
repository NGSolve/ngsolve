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
  class HCurlFiniteElement : public FiniteElement
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



  /* **************************** Segm Elements *************** */


  ///
  class FE_NedelecSegm1 : public HCurlFiniteElement<1>
  {
  public:

    enum { NDOF = 1 };

  private:
    ///
    // static Array<IPData> ipdata;
    // bool ipdatadestructed;

  public:
    ///
    FE_NedelecSegm1();
    ///
    virtual ~FE_NedelecSegm1();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<1> shape) const;
  };



  ///
  class FE_NedelecSegm2 : public HCurlFiniteElement<1>
  {
  public:

    enum { NDOF = 2 };

  private:
    ///
    // static Array<IPData> ipdata;
    // bool ipdatadestructed;

  public:
    ///
    FE_NedelecSegm2();
    ///
    virtual ~FE_NedelecSegm2();

    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<1> shape) const;
  };



  ///
  class FE_NedelecSegm3 : public HCurlFiniteElement<1>
  {
    ///
// static Array<IPData> ipdata;
// bool ipdatadestructed;

  public:
    ///
    FE_NedelecSegm3();
    ///
    virtual ~FE_NedelecSegm3();

    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<1> shape) const;
  };






  /* *********************** Quad elements ******************* */

  class FE_NedelecQuad1 : public T_HCurlFiniteElement<FE_NedelecQuad1,ET_QUAD,4,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      Tx x = hx[0], y = hx[1];
      
      AutoDiff<2> lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
      AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
      
      const EDGE * edges = ElementTopology::GetEdges (ET_QUAD);
      for (int i = 0; i < 4; i++)
        {
          int es = edges[i][0], ee = edges[i][1];
          
          AutoDiff<2> xi  = sigma[ee]-sigma[es];
          AutoDiff<2> lam_e = lami[ee]+lami[es];  
          
          shape[i] = uDv<2> (0.5 * lam_e, xi); 
        }
    }
  };

  /*
  /// Gradients of Q1
  class FE_NedelecQuad1 : public HCurlFiniteElement<2>
  {
  public:
    enum { NDOF = 4 };

  protected:
    ///
    // static Array<IPData> ipdata;
    // bool ipdatadestructed;
  public:
    ///
    FE_NedelecQuad1();
    ///
    virtual ~FE_NedelecQuad1();

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<2> shape) const;

  };
  */


  /*
    template <int ORDER, int ZORDER>
    class FE_TNedelecQuadTraits
    {
    public:
    enum { NDOF = ORDER * (ZORDER+1) + (ORDER+1) * ZORDER };
    enum { NEDGEDOF = 2 * (ORDER + ZORDER) - 4 };
    };
  */

  template <int ORDER, int ZORDER>
  class FE_TNedelecQuad : public HCurlFiniteElement<2>
  {
  public:
    enum { NDOF = ORDER * (ZORDER+1) + (ORDER+1) * ZORDER };
    enum { NEDGEDOF = 2 * (ORDER + ZORDER) - 4 };

  protected:
    ///
    // static Array<IPData> ipdata;
    // bool ipdatadestructed;
    ///
    // static Mat<FE_TNedelecQuad<ORDER,ZORDER>::NDOF> trans;

    // static Mat<FE_TNedelecQuadTraits<ORDER,ZORDER>::NDOF> trans;
    // static Mat<FE_TNedelecQuadTraits<ORDER,ZORDER>::NEDGEDOF> trans2;

    static Matrix<> trans;
    static Matrix<> trans2;

    FE_NedelecQuad1 quad1;

  public:
    enum { MAXORDER = (ORDER > ZORDER) ? ORDER : ZORDER };

    ///
    FE_TNedelecQuad();
    ///
    virtual ~FE_TNedelecQuad();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<2> shape) const;
    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<2> shape) const;
    ///
    virtual void CalcShape2 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<2> shape) const;
    ///
    void Orthogonalize();  
  };




  /* ******************** triangular elements *********************** */

  class FE_NedelecTrig1 : public T_HCurlFiniteElement<FE_NedelecTrig1,ET_TRIG,3,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      Tx x = hx[0], y = hx[1];
      Tx lami[3] = { x, y, 1-x-y };
      
      const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
      for (int i = 0; i < 3; i++)
        shape[i] = uDv_minus_vDu<2> (lami[edges[i][0]], lami[edges[i][1]]);
    }
  };

  class FE_NedelecTrig2 : public T_HCurlFiniteElement<FE_NedelecTrig2,ET_TRIG,6,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      Tx x = hx[0], y = hx[1];
      Tx lami[3] = { x, y, 1-x-y };
      
      const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
      for (int i = 0; i < 3; i++)
        {
          shape[i] = uDv_minus_vDu<2> (lami[edges[i][0]], lami[edges[i][1]]);
          shape[i+3] = Du<2> (lami[edges[i][0]]*lami[edges[i][1]]);
        }
    }
  };

  class FE_NedelecTrig3 : public T_HCurlFiniteElement<FE_NedelecTrig3,ET_TRIG,12,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      Tx x = hx[0], y = hx[1];
      Tx lami[3] = { x, y, 1-x-y };
      
      const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
      for (int i = 0; i < 3; i++)
        {
          Tx lam1 = lami[edges[i][0]];
          Tx lam2 = lami[edges[i][1]];

          shape[i] = uDv_minus_vDu<2> (lam1, lam2);
          shape[i+3] = Du<2> (lam1*lam2);
          shape[i+6] = Du<2> (lam1*lam2*(lam1-lam2));
        }

      const FACE * faces = ElementTopology::GetFaces (ET_TRIG); 
      for (int k = 0; k < 3; k++)
        {
          int k1 = (k+1)%3, k2 = (k+2)%3;
          shape[9+k] = uDv_minus_vDu<2> (lami[faces[0][k]],
                                         lami[faces[0][k1]]*lami[faces[0][k2]]);
        }

    }
  };



  /*
  /// Lowest order Nedelec
  class FE_NedelecTrig1 : public HCurlFiniteElement<2>
  {
    ///
    // static Array<IPData> ipdata;
    // bool ipdatadestructed;

  public:

    ///
    FE_NedelecTrig1();
    ///
    virtual ~FE_NedelecTrig1();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<2> shape) const;
  };



  /// Nedelec type 2, order 1, gradients of P2
  class FE_NedelecTrig2 : public HCurlFiniteElement<2>
  {
  public:
    enum { NDOF = 6 };

  private:
    ///
    // static Array<IPData> ipdata;
    // bool ipdatadestructed;
    ///
    static Mat<NDOF> trans;

  public:
    ///
    FE_NedelecTrig2();
    ///
    virtual ~FE_NedelecTrig2();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<2> shape) const;

    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<2> shape) const;

    ///
    void Orthogonalize();
  };


  /// Nedelec type 2, order 2, gradients of P3
  class FE_NedelecTrig3 : public HCurlFiniteElement<2>
  {
  public:
    enum { NDOF = 12 };
    enum { NEDGEDOF = 6 };
    ///
    // static Array<IPData> ipdata;
    // bool ipdatadestructed;
    ///
    static Mat<NDOF> trans;
    ///
    static Mat<NEDGEDOF> trans2;
    ///
    FE_NedelecTrig2 trig1;
  public:
    ///
    FE_NedelecTrig3();
    ///
    virtual ~FE_NedelecTrig3();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<2> shape) const;

    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<2> shape) const;

    ///
    virtual void CalcShape2 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<2> shape) const;
  
    ///
    void Orthogonalize();
  };
  */










  /* *********************** Tetrahedral elements ********************** */
  
  class FE_NedelecTet1 : public T_HCurlFiniteElement<FE_NedelecTet1,ET_TET,6,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[3], TFA & shape) 
    {
      // Tx x = hx[0], y = hx[1], z = hx[2];
      // Tx lami[4] = { x, y, z, 1-x-y-z };
      Tx lami[4] = { x[0], x[1], x[2], 1-x[0]-x[1]-x[2] };      

      const EDGE * edges = ElementTopology::GetEdges (ET_TET);
      for (int i = 0; i < 6; i++)
        shape[i] = uDv_minus_vDu<3> (lami[edges[i][0]], lami[edges[i][1]]);
    }
  };

  class FE_NedelecTet2 : public T_HCurlFiniteElement<FE_NedelecTet2,ET_TET,12,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[3], TFA & shape) 
    {
      // Tx x = hx[0], y = hx[1], z = hx[2];
      // Tx lami[4] = { x, y, z, 1-x-y-z };
      Tx lami[4] = { x[0], x[1], x[2], 1-x[0]-x[1]-x[2] };      
      
      const EDGE * edges = ElementTopology::GetEdges (ET_TET);
      for (int i = 0; i < 6; i++)
        {
          shape[i] = uDv_minus_vDu<3> (lami[edges[i][0]], lami[edges[i][1]]);
          shape[i+6] = Du<3> (lami[edges[i][0]]*lami[edges[i][1]]);
        }
    }
  };


  class FE_NedelecTet3 : public T_HCurlFiniteElement<FE_NedelecTet3,ET_TET,30,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[3], TFA & shape)
    {
      Tx lami[4] = { x[0], x[1], x[2], 1-x[0]-x[1]-x[2] };      
      
      const EDGE * edges = ElementTopology::GetEdges (ET_TET);
      for (int i = 0; i < 6; i++)
        {
          Tx lam1 = lami[edges[i][0]];
          Tx lam2 = lami[edges[i][1]];
          shape[i] = uDv_minus_vDu<3> (lam1, lam2);
          shape[i+6] = Du<3> (lam1*lam2);
          shape[i+12] = Du<3> (lam1*lam2*(lam1-lam2));
        }

      const FACE * faces = ElementTopology::GetFaces (ET_TET); 
      for (int i = 0; i < 4; i++)
        for (int k = 0; k < 3; k++)
          {
            int k1 = (k+1)%3, k2 = (k+2)%3;
            shape[18+3*i+k] = uDv_minus_vDu<3> (lami[faces[i][k]],
                                                lami[faces[i][k1]]*lami[faces[i][k2]]);
          }
    }
  };


  /*
    class FE_NedelecTet1o : public HCurlFiniteElement<3>
    {
    public:
    enum { NDOF = 6 };

    private:
    ///
      // static Array<IPData> ipdata;
      // bool ipdatadestructed;

    public:

    ///
    FE_NedelecTet1o();
    ///
    virtual ~FE_NedelecTet1o();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<3> shape) const;

    virtual void CalcCurlShape (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<3> curlshape) const;
    };


  ///
  class FE_NedelecTet2o : public HCurlFiniteElement<3>
  {
  public:

    enum { NDOF = 12 };

  private:
    ///
    // static Array<IPData> ipdata;
    // bool ipdatadestructed;
    ///
    static Mat<NDOF> trans;

  public:
    ///
    FE_NedelecTet2o();
    ///
    virtual ~FE_NedelecTet2o();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
                            FlatMatrixFixWidth<3> shape) const;


    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
                             FlatMatrixFixWidth<3> shape) const;

    ///
    void Orthogonalize();
  };


  /// 2nd order Nedelec element of class II
  class FE_NedelecTet3 : public HCurlFiniteElement<3>
  {
  public:
    enum { NDOF = 30 };
    enum { NEDGEDOF = 12 };
    enum { NFACEDOF = 12 };

  protected:
    ///
    //     static Array<IPData> ipdata;
    // bool ipdatadestructed;
    ///
    static Mat<NDOF> trans;
    ///
    static Mat<NEDGEDOF> trans2;
    ///
    static Mat<NFACEDOF> trans3;

    FE_NedelecTet1 tet1;
  public:

    ///
    FE_NedelecTet3();
    ///
    virtual ~FE_NedelecTet3();


    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<3> shape) const;

    virtual void CalcCurlShape (const IntegrationPoint & ip, 
                                FlatMatrixFixWidth<3> curlshape) const;
  
    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
                             FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcShape2 (const IntegrationPoint & ip, 
                             FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcShape3 (const IntegrationPoint & ip, 
                             FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcCurlShape3 (const IntegrationPoint & ip, 
                                 FlatMatrixFixWidth<3> shape) const;

    ///
    void Orthogonalize();
  };
  */


  /// 2nd order Nedelec element of class II, without gradient fields
  class FE_NedelecTet3NoGrad : public HCurlFiniteElement<3>
  {
  public:
    enum { NDOF = 18 };
    enum { NFACEDOF = 12 };

  protected:
    static Mat<NFACEDOF> trans3;

    FE_NedelecTet1 tet1;
  public:
    FE_NedelecTet3NoGrad();
    virtual ~FE_NedelecTet3NoGrad();


    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<3> shape) const;

    virtual void CalcCurlShape (const IntegrationPoint & ip, 
                                FlatMatrixFixWidth<3> curlshape) const;

    virtual void CalcShape3 (const IntegrationPoint & ip, 
                             FlatMatrixFixWidth<3> shape) const;

    virtual void CalcCurlShape3 (const IntegrationPoint & ip, 
                                 FlatMatrixFixWidth<3> shape) const;

    void Orthogonalize();
  };



  /* *********************** Hex elements ************************ */ 



  /// 
  class FE_NedelecHex1 : public HCurlFiniteElement<3> 
  {
    /// 
// static Array<IPData> ipdata; 
// bool ipdatadestructed;
  
  public: 
    ///
    FE_NedelecHex1(); 
    ///
    virtual ~FE_NedelecHex1(); 
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
                            FlatMatrixFixWidth<3> shape) const; 
  }; 


  /* *********************** Prism elements ********************** */

  class FE_NedelecPrism1 : public T_HCurlFiniteElement<FE_NedelecPrism1,ET_PRISM,9,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      Tx x = hx[0], y = hx[1], z = hx[2];

      Tx lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
      Tx muz[6]  = { 1-z, 1-z, 1-z, z, z, z };
       
      const EDGE * edges = ElementTopology::GetEdges (ET_PRISM);
  
      // horizontal edge shapes
      for (int i = 0; i < 6; i++)
        {
          int es = edges[i][0], ee = edges[i][1]; 
          shape[i] = wuDv_minus_wvDu<3> (lami[es], lami[ee], muz[ee]);
        }
      
      //Vertical Edge Shapes
      for (int i = 6; i < 9; i++)
        {
          int es = edges[i][0], ee = edges[i][1];
          shape[i] = wuDv_minus_wvDu<3> (muz[es], muz[ee], lami[ee]);
        }
    }
  };
  
  /*
  ///
  class FE_NedelecPrism1 : public HCurlFiniteElement<3>
  {
    ///
    // static Array<IPData> ipdata;
    // bool ipdatadestructed;

  public:

    ///
    FE_NedelecPrism1();
    ///
    virtual ~FE_NedelecPrism1();

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<3> shape) const;
  };
  */

  /// \f$ \nabla Q (2,ZORDER) \f$
  template <int ZORDER>
  class FE_TNedelecPrism2 : public HCurlFiniteElement<3>
  {
    ///
// static Array<IPData> ipdata;
// bool ipdatadestructed;
    ///
    static Matrix<> trans;
    ///
    static Matrix<> trans2;
    ///
    static Matrix<> trans3;

    FE_NedelecPrism1 prism1;

  public:
    enum { NDOF = 6 * (ZORDER+1) + 6 * ZORDER };
    enum { NEDGEDOF = 6 + 3 * (ZORDER-1) };
    enum { NFACEDOF = 9 * ZORDER - 6} ;
    enum { MAXORDER = (2 > ZORDER) ? 2 : ZORDER };

    ///
    FE_TNedelecPrism2();
    ///
    virtual ~FE_TNedelecPrism2();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcShape2 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcShape3 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    ///
    void Orthogonalize();

  };



  /// potential space for Nedelec IIb
  class FE_Trig3Pot : public ScalarFiniteElement<2>
  {
    ///
// static IPDataArray ipdata;
  public:
  ///
    FE_Trig3Pot();
    ///
    virtual ~FE_Trig3Pot();

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
  }; 




  /// \f$ \nabla Q (3,ZORDER) \f$
  template <int ZORDER>
  class FE_TNedelecPrism3 : public HCurlFiniteElement<3>
  {
    ///
// static Array<IPData> ipdata;
// bool ipdatadestructed;
    ///
    static Matrix<> trans;
    ///
    static Matrix<> trans2;
    ///
    static Matrix<> trans_quad;
    ///
    static Matrix<> trans_trig;

    FE_NedelecPrism1 prism1;
    FE_NedelecTrig3 trig3;
    static class FE_Trig2 h1trig2;
    FE_Trig3Pot h1trig3;
    FE_TSegmL2<ZORDER> segm;
  public:
    enum { NDOF = 12 * (ZORDER+1) + 10 * ZORDER };
    enum { NEDGEDOF = 12 + 3 * (ZORDER-1) };
    enum { NQUADFACEDOF = 3 * (5*ZORDER-3) };
    enum { NTRIGFACEDOF = 6 };
    enum { MAXORDER = (3 > ZORDER) ? 3 : ZORDER };
    enum { NINNERDOF = 3 * (ZORDER-1) + ZORDER };

    ///
    FE_TNedelecPrism3();
    ///
    virtual ~FE_TNedelecPrism3();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcShape2 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    /// quad face dofs
    virtual void CalcShape3 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    /// trig face dofs
    virtual void CalcShape4 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcInner (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void GetInternalDofs (Array<int> & idofs) const;

    ///
    void Orthogonalize();
  };





  /// \f$ \nabla Q (3,ZORDER) \f$
  template <int ZORDER>
  class FE_TNedelecPrism3NoGrad : public HCurlFiniteElement<3>
  {
    ///
// static Array<IPData> ipdata;
// bool ipdatadestructed;
    ///
    static Matrix<> trans_quad;
    ///
    static Matrix<> trans_trig;

    FE_NedelecPrism1 prism1;
    FE_NedelecTrig3 trig3;
    static class FE_Trig2 h1trig2;
    FE_Trig3Pot h1trig3;
    FE_TSegmL2<ZORDER> segm;
  public:
    //  enum { NDOF = 12 * (ZORDER+1) + 10 * ZORDER };
    //  enum { NEDGEDOF = 12 + 3 * (ZORDER-1) };
    // 12 z + 12 + 10 z - 12 - 3z + 3 = 19 z + 3
    enum { NDOF = 19 * ZORDER + 3 };
    enum { NQUADFACEDOF = 3 * (5*ZORDER-3) };
    enum { NTRIGFACEDOF = 6 };
    enum { MAXORDER = (3 > ZORDER) ? 3 : ZORDER };
    // enum { NINNERDOF = 3 * (ZORDER-1) + ZORDER };
    enum { NINNERDOF = 3 * (ZORDER-1) + 1 };

    ///
    FE_TNedelecPrism3NoGrad();
    ///
    virtual ~FE_TNedelecPrism3NoGrad();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcShape2 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    /// quad face dofs
    virtual void CalcShape3 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    /// trig face dofs + inner dofs
    virtual void CalcShape4 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void CalcInner (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void GetInternalDofs (Array<int> & idofs) const;

    ///
    void Orthogonalize();
  };














  ///
  class FE_NedelecPyramid1 : public HCurlFiniteElement<3>
  {
    ///
// static Array<IPData> ipdata;
// bool ipdatadestructed;
    ///
    static Matrix<> trans;
  public:
    ///
    FE_NedelecPyramid1();
    ///
    virtual ~FE_NedelecPyramid1();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<3> shape) const;

    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    ///
    void Orthogonalize();

  };



  ///
  class FE_NedelecPyramid2 : public HCurlFiniteElement<3>
  {
  public:
    enum { NDOF = 20 };
    enum { NEDGEDOF = 8 };

  private:
    ///
    // static Array<IPData> ipdata;
    // bool ipdatadestructed;
    ///
    static Matrix<> trans;
    static Matrix<> trans2;
    static Matrix<> trans3;
  
    ///
    // static class FE_Quad1 quad1;
    // static class FE_Quad2 quad2;
    typedef FE_Quad1 quad1;
    typedef FE_Quad2 quad2;

    FE_NedelecPyramid1 pyramid1;
  public:
    ///
    FE_NedelecPyramid2();
    ///
    virtual ~FE_NedelecPyramid2();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<3> shape) const;

    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    virtual void CalcShape2 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;


    ///
    void Orthogonalize();

  };


  /// quad of order 3
  class FE_Quad3 : public ScalarFiniteElement<2>
  {
    // static IPDataArray ipdata;

  public:
    FE_Quad3();
    virtual ~FE_Quad3();
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<2> dshape) const;
  }; 



  ///
  class FE_NedelecPyramid3 : public HCurlFiniteElement<3>
  {
  public:
    enum { NDOF = 57 };
    enum { NEDGEDOF = 16 };
    enum { NFACEDOF = 24 };
    enum { NINNERDOF = 9 };
  private:
    ///
    // static Array<IPData> ipdata;
    // bool ipdatadestructed;
    ///
    static Mat<NDOF> trans;
    static Mat<NEDGEDOF> trans2;
    static Mat<NFACEDOF> trans3;
  
    ///
    // static class FE_Quad1 quad1;
    // static class FE_Quad2 quad2;
    typedef FE_Quad1 quad1;
    typedef FE_Quad2 quad2;

    FE_Quad3 quad3;
    FE_NedelecPyramid1 pyramid1;
  public:
    ///
    FE_NedelecPyramid3();
    ///
    virtual ~FE_NedelecPyramid3();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatMatrixFixWidth<3> shape) const;

    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    virtual void CalcShape2 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    virtual void CalcShape3 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const;

    ///
    virtual void GetInternalDofs (Array<int> & idofs) const;
    ///
    void Orthogonalize();
  };





  /*
 ///
 class FE_NedelecPyramid3NoGrad : public HCurlFiniteElement<3>
 {
 public:
 //  enum { NDOF = 57 };
 // enum { NEDGEDOF = 16 };
 enum { NDOF = 41 };
 enum { NFACEDOF = 24 };
 enum { NINNERDOF = 9 };
 private:
 ///
 static Array<IPData> ipdata;
 ///
 static Mat<NDOF> trans;
 // static Mat<NEDGEDOF> trans2;
 static Mat<NFACEDOF> trans3;
  
 ///
 FE_Quad1 quad1;
 FE_Quad2 quad2;
 FE_Quad3 quad3;
 FE_NedelecPyramid1 pyramid1;
 public:
 ///
 FE_NedelecPyramid3();
 ///
 virtual ~FE_NedelecPyramid3();
 ///
 virtual void CalcShape (const IntegrationPoint & ip, 
 FlatMatrixFixWidth<3> shape) const;

 virtual void CalcShape1 (const IntegrationPoint & ip, 
 FlatMatrixFixWidth<3> shape) const;

 virtual void CalcShape2 (const IntegrationPoint & ip, 
 FlatMatrixFixWidth<3> shape) const;

 virtual void CalcShape3 (const IntegrationPoint & ip, 
 FlatMatrixFixWidth<3> shape) const;

 ///
 virtual void GetInternalDofs (Array<int> & idofs) const;
 ///
 void Orthogonalize();
 };
  */


}




#endif
