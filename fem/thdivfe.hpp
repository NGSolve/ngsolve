#ifndef FILE_THDIVFE
#define FILE_THDIVFE

/*********************************************************************/
/* File:   thdivfe.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   5. Jul. 2001                                              */
/*********************************************************************/

namespace ngfem
{


  template <int D>
  inline double Dot (const AutoDiff<D> & u, const AutoDiff<D> & v)
  {
    double sum = 0;
    for (int i = 0; i < D; i++)
      sum += u.DValue(i) * v.DValue(i);
    return sum;
  }


  // rotated gradient
  template <int DIM> class DuRot;

  template <> class DuRot<2>
  {

  public:
    const AutoDiff<2> & u;

    DuRot (const AutoDiff<2> & au)
      : u(au) { ; }

    Vec<2> Value () const
    {
      Vec<2> val;
      val(0) = u.DValue(1);
      val(1) = -u.DValue(0);
      return val;
    }

    /*
    Vec<DIM_CURL> CurlValue () const
    {
      return Vec<DIM> (0.0);
    }
    */
  };





  template <int DIM>
  class uDvDw_Cyclic
  {
  public:
  public:
    const AutoDiff<DIM> & u, v, w;
    uDvDw_Cyclic (const AutoDiff<DIM> & au, 
                  const AutoDiff<DIM> & av,
                  const AutoDiff<DIM> & aw)
      : u(au), v(av), w(aw) { ; }
  };

  template <int DIM>
  class Du_Cross_Dv
  {
  public:
  public:
    const AutoDiff<DIM> & u, v;
    Du_Cross_Dv (const AutoDiff<DIM> & au, 
                 const AutoDiff<DIM> & av)
      : u(au), v(av) { ; }
  };

  template <int DIM>
  class wDu_Cross_Dv
  {
  public:
  public:
    const AutoDiff<DIM> & u, v, w;
    wDu_Cross_Dv (const AutoDiff<DIM> & au, 
                  const AutoDiff<DIM> & av,
                  const AutoDiff<DIM> & aw)
      : u(au), v(av), w(aw) { ; }
  };


  template <int DIM>
  class uDvDw_minus_DuvDw
  {
  public:
  public:
    const AutoDiff<DIM> & u, v, w;
    uDvDw_minus_DuvDw (const AutoDiff<DIM> & au, 
                       const AutoDiff<DIM> & av,
                       const AutoDiff<DIM> & aw)
      : u(au), v(av), w(aw) { ; }
  };

  template <int DIM>
  class curl_uDvw_minus_Duvw
  {
  public:
  public:
    const AutoDiff<DIM> & u, v, w;
    curl_uDvw_minus_Duvw (const AutoDiff<DIM> & au, 
                          const AutoDiff<DIM> & av,
                          const AutoDiff<DIM> & aw)
      : u(au), v(av), w(aw) { ; }
  };



  template <int DIM> class THDiv2Shape
  {
  public:
    INLINE operator Vec<DIM> () { return 0.0; }
  };


  template <> class THDiv2Shape<2>
  {
    Vec<2> data;
  public:
    INLINE THDiv2Shape (Du<2> uv)
    {
      data = Vec<2> (uv.u.DValue(1), -uv.u.DValue(0));
    }
    
    INLINE THDiv2Shape (uDv<2> uv)
    {
      data = Vec<2> (-uv.u.Value()*uv.v.DValue(1), 
                     uv.u.Value()*uv.v.DValue(0));
    }

    INLINE THDiv2Shape (const uDv_minus_vDu<2> & uv) 
    { 
      data(0) = -uv.u.Value() * uv.v.DValue(1) + uv.u.DValue(1) * uv.v.Value();
      data(1) =  uv.u.Value() * uv.v.DValue(0) - uv.u.DValue(0) * uv.v.Value();
    }

    INLINE THDiv2Shape (const wuDv_minus_wvDu<2> & uv) 
    { 
      data[0] = -uv.u.Value() * uv.v.DValue(1) + uv.u.DValue(1) * uv.v.Value();
      data[1] =  uv.u.Value() * uv.v.DValue(0) - uv.u.DValue(0) * uv.v.Value();
      data[0] *= uv.w.Value();
      data[1] *= uv.w.Value();
    }
    
    INLINE operator Vec<2> () { return data; }
  };


  template <> class THDiv2Shape<3>
  {
    Vec<3> data;
  public:

    INLINE THDiv2Shape (const uDvDw_Cyclic<3> & uvw) 
    { 
      AutoDiff<3> hv =
        uvw.u.Value() * Cross (uvw.v, uvw.w) +
        uvw.v.Value() * Cross (uvw.w, uvw.u) +
        uvw.w.Value() * Cross (uvw.u, uvw.v);

      for (int i = 0; i < 3; i++)
        data[i] = hv.DValue(i);
    }

    INLINE THDiv2Shape (const Du_Cross_Dv<3> & uv) 
    { 
      AutoDiff<3> hv = Cross (uv.u, uv.v);
      for (int i = 0; i < 3; i++)
        data[i] = hv.DValue(i);
    }

    INLINE THDiv2Shape (const wDu_Cross_Dv<3> & uvw) 
    { 
      AutoDiff<3> hv = Cross (uvw.u, uvw.v);
      for (int i = 0; i < 3; i++)
        data[i] = uvw.w.Value() * hv.DValue(i);
    }


    INLINE THDiv2Shape (const uDvDw_minus_DuvDw<3> & uvw) 
    { 
      AutoDiff<3> hv =
        uvw.u.Value() * Cross (uvw.v, uvw.w) +
        uvw.v.Value() * Cross (uvw.w, uvw.u);

      for (int i = 0; i < 3; i++)
        data[i] = hv.DValue(i);
    }

    INLINE THDiv2Shape (const curl_uDvw_minus_Duvw<3> & uvw) 
    { 
      AutoDiff<3> hv = Cross (uvw.u*uvw.w, uvw.v) - Cross (uvw.v*uvw.w, uvw.u);
      for (int i = 0; i < 3; i++)
        data[i] = hv.DValue(i);
    }


    INLINE operator Vec<3> () { return data; }
  };




  // 2D 
  template <int DIM>
  class HDivShapeElement
  {
    double * data;
  public:
    HDivShapeElement (double * adata) : data(adata) { ; }

    void operator= (THDiv2Shape<DIM> hd2vec)
    {
      Vec<DIM> v = hd2vec;
      for (int j = 0; j < DIM; j++)
        data[j] = v(j);
    }


    /*
    void operator= (const Du<DIM> & uv) 
    {
      data[0] =  uv.u.DValue(1);
      data[1] = -uv.u.DValue(0);
    }

    void operator= (const uDv<DIM> & uv) 
    { 
      data[0] = -uv.u.Value() * uv.v.DValue(1);
      data[1] =  uv.u.Value() * uv.v.DValue(0);
    }

    void operator= (const uDv_minus_vDu<DIM> & uv) 
    { 
      data[0] = -uv.u.Value() * uv.v.DValue(1) + uv.u.DValue(1) * uv.v.Value();
      data[1] =  uv.u.Value() * uv.v.DValue(0) - uv.u.DValue(0) * uv.v.Value();
    }
    
    void operator= (const wuDv_minus_wvDu<DIM> & uv) 
    { 
      data[0] = -uv.u.Value() * uv.v.DValue(1) + uv.u.DValue(1) * uv.v.Value();
      data[1] =  uv.u.Value() * uv.v.DValue(0) - uv.u.DValue(0) * uv.v.Value();
      data[0] *= uv.w.Value();
      data[1] *= uv.w.Value();
    }

    
    void operator= (const uDvDw_Cyclic<DIM> & uvw) 
    { 
      AutoDiff<3> hv =
        uvw.u.Value() * Cross (uvw.v, uvw.w) +
        uvw.v.Value() * Cross (uvw.w, uvw.u) +
        uvw.w.Value() * Cross (uvw.u, uvw.v);

      for (int i = 0; i < 3; i++)
        data[i] = hv.DValue(i);
    }

    void operator= (const Du_Cross_Dv<DIM> & uv) 
    { 
      AutoDiff<3> hv = Cross (uv.u, uv.v);
      for (int i = 0; i < 3; i++)
        data[i] = hv.DValue(i);
    }

    void operator= (const wDu_Cross_Dv<DIM> & uvw) 
    { 
      AutoDiff<3> hv = Cross (uvw.u, uvw.v);
      for (int i = 0; i < 3; i++)
        data[i] = uvw.w.Value() * hv.DValue(i);
    }


    void operator= (const uDvDw_minus_DuvDw<DIM> & uvw) 
    { 
      AutoDiff<3> hv =
        uvw.u.Value() * Cross (uvw.v, uvw.w) +
        uvw.v.Value() * Cross (uvw.w, uvw.u);

      for (int i = 0; i < 3; i++)
        data[i] = hv.DValue(i);
    }

    void operator= (const curl_uDvw_minus_Duvw<DIM> & uvw) 
    { 
      AutoDiff<3> hv = Cross (uvw.u*uvw.w, uvw.v) - Cross (uvw.v*uvw.w, uvw.u);
      for (int i = 0; i < 3; i++)
        data[i] = hv.DValue(i);
    }
    */
  };










  template <int DIM>
  class HDivEvaluateShapeElement
  {
    const double * coefs;
    Vec<DIM> & sum;
  public:
    HDivEvaluateShapeElement (const double * acoefs, Vec<DIM> & asum)
      : coefs(acoefs), sum(asum) { ; }


    void operator= (THDiv2Shape<DIM> hd2vec)
    {
      sum += *coefs * Vec<DIM> (hd2vec);
    }

    /*
    void operator= (const Du<DIM> & uv) 
    {
      sum(0) += *coefs * uv.u.DValue(1);
      sum(1) -= *coefs * uv.u.DValue(0);
    }

    void operator= (const uDv<DIM> & uv) 
    { 
      sum(0) -= (*coefs) * uv.u.Value() * uv.v.DValue(1);
      sum(1) += (*coefs) * uv.u.Value() * uv.v.DValue(0);
    }

    void operator= (const uDv_minus_vDu<DIM> & uv) 
    { 
      sum(0) += (*coefs) * (-uv.u.Value() * uv.v.DValue(1) + uv.u.DValue(1) * uv.v.Value());
      sum(1) += (*coefs) * ( uv.u.Value() * uv.v.DValue(0) - uv.u.DValue(0) * uv.v.Value());
    }

    void operator= (const wuDv_minus_wvDu<DIM> & uv) 
    { 
      double fac = *coefs * uv.w.Value();
      sum(0) += fac * (-uv.u.Value() * uv.v.DValue(1) + uv.u.DValue(1) * uv.v.Value());
      sum(1) += fac * ( uv.u.Value() * uv.v.DValue(0) - uv.u.DValue(0) * uv.v.Value());
    }


    void operator= (const uDvDw_Cyclic<DIM> & uvw) 
    { 
      AutoDiff<3> hv =
        uvw.u.Value() * Cross (uvw.v, uvw.w) +
        uvw.v.Value() * Cross (uvw.w, uvw.u) +
        uvw.w.Value() * Cross (uvw.u, uvw.v);
      
      for (int i = 0; i < 3; i++)
        sum(i) += (*coefs) * hv.DValue(i);
    }

    void operator= (const Du_Cross_Dv<DIM> & uv) 
    { 
      AutoDiff<3> hv = Cross (uv.u, uv.v);
      for (int i = 0; i < 3; i++)
        sum(i) += (*coefs) * hv.DValue(i);
    }

    void operator= (const wDu_Cross_Dv<DIM> & uvw) 
    { 
      AutoDiff<3> hv = Cross (uvw.u, uvw.v);
      for (int i = 0; i < 3; i++)
        sum(i) += (*coefs) * uvw.w.Value() * hv.DValue(i);

    }

    void operator= (const uDvDw_minus_DuvDw<DIM> & uvw) 
    { 
      AutoDiff<3> hv =
        uvw.u.Value() * Cross (uvw.v, uvw.w) +
        uvw.v.Value() * Cross (uvw.w, uvw.u);

      for (int i = 0; i < 3; i++)
        sum(i) += (*coefs) * hv.DValue(i);
    }

    void operator= (const curl_uDvw_minus_Duvw<DIM> & uvw) 
    { 
      AutoDiff<3> hv = Cross (uvw.u*uvw.w, uvw.v) - Cross (uvw.v*uvw.w, uvw.u);
      for (int i = 0; i < 3; i++)
        sum(i) += (*coefs) * hv.DValue(i);

    }
    */
  };








  template <int DIM>
  class HDivDivShapeElement
  {
    double * data;
  public:
    HDivDivShapeElement (double * adata) : data(adata) { ; }

    void operator= (const Du<DIM> & uv) 
    { 
      data[0] = 0;
    }

    void operator= (const uDv<DIM> & uv) 
    {
      AutoDiff<1> hd = Cross (uv.u, uv.v);
      data[0] = -hd.DValue(0);
    }

    void operator= (const uDv_minus_vDu<DIM> & uv) 
    {
      data[0] = -2*uv.u.DValue(0) * uv.v.DValue(1) 
        + 2*uv.u.DValue(1) * uv.v.DValue(0);
    }

    void operator= (const wuDv_minus_wvDu<DIM> & uv) 
    {
      AutoDiff<1> hd = Cross (uv.u*uv.w, uv.v) + Cross(uv.u, uv.v*uv.w);
      data[0] = -hd.DValue(0);
    }


    void operator= (const uDvDw_Cyclic<DIM> & uvw) 
    { 
      data[0] = 
        Dot (uvw.u, Cross (uvw.v, uvw.w)) +
        Dot (uvw.v, Cross (uvw.w, uvw.u)) +
        Dot (uvw.w, Cross (uvw.u, uvw.v));
    }


    void operator= (const Du_Cross_Dv<DIM> & uv) 
    { 
      data[0] = 0;
    }

    void operator= (const wDu_Cross_Dv<DIM> & uvw) 
    { 
      data[0] = Dot (uvw.w, Cross (uvw.u, uvw.v));
    }

    void operator= (const uDvDw_minus_DuvDw<DIM> & uv) 
    { 
      data[0] = 
        Dot (uv.u, Cross (uv.v, uv.w)) +
        Dot (uv.v, Cross (uv.w, uv.u));
    }

    void operator= (const curl_uDvw_minus_Duvw<DIM> & uvw) 
    { 
      data[0] = 0;
    }
      
  };


  template <int DIM>
  class HDivShapeAssign
  {
    double * dshape;
  public:
    HDivShapeAssign (FlatMatrixFixWidth<DIM> mat)
    { dshape = &mat(0,0); }
    
    HDivShapeElement<DIM> operator[] (int i) const
    { return HDivShapeElement<DIM> (dshape + i*DIM); }
  };


  template <int DIM>
  class HDivDivShapeAssign
  {
    SliceVector<> dshape;
  public:
    HDivDivShapeAssign (SliceVector<>  mat)
      : dshape(mat) { ; }

    HDivDivShapeElement<DIM> operator[] (int i) const
    { return HDivDivShapeElement<DIM> (&dshape(i)); }
  };

  template <int DIM>
  class HDivEvaluateShape
  {
    const double * coefs;
    Vec<DIM> sum;
  public:
    HDivEvaluateShape (FlatVector<> acoefs)
    { coefs = &acoefs(0); sum = 0.0; }

    HDivEvaluateShapeElement<DIM> operator[] (int i) 
    { return HDivEvaluateShapeElement<DIM> (coefs+i, sum); }

    Vec<DIM> Sum() { return sum; }
  };






  







  template <class FEL, ELEMENT_TYPE ET>
  class T_HDivFiniteElement 
    : virtual public HDivFiniteElement<ET_trait<ET>::DIM>
    
  {
    enum { DIM = ET_trait<ET>::DIM };

  public:

    virtual void CalcShape (const IntegrationPoint & ip, 
			    SliceMatrix<> shape) const;
    
    virtual void CalcDivShape (const IntegrationPoint & ip, 
			       SliceVector<> divshape) const;
    
    virtual void CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
				  SliceMatrix<> shape) const;
    
    virtual void Evaluate (const IntegrationRule & ir, 
			   FlatVector<double> coefs, 
			   FlatMatrixFixWidth<DIM> vals) const;
  };

}


#endif
