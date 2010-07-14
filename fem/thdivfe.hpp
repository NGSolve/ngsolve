#ifndef FILE_THDIVFE
#define FILE_THDIVFE

/*********************************************************************/
/* File:   thdivfe.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   5. Jul. 2001                                              */
/*********************************************************************/

namespace ngfem
{






  /*
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
  */

  template <int D>
  inline double Dot (const AutoDiff<D> & u, const AutoDiff<D> & v)
  {
    double sum = 0;
    for (int i = 0; i < D; i++)
      sum += u.DValue(i) * v.DValue(i);
    return sum;
  }



  /*
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
  */

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








  // 2D 
  template <int DIM>
  class HDivShapeElement
  {
    double * data;
  public:
    HDivShapeElement (double * adata) : data(adata) { ; }
    
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
  };










  template <int DIM>
  class HDivEvaluateShapeElement
  {
    const double * coefs;
    Vec<DIM> & sum;
  public:
    HDivEvaluateShapeElement (const double * acoefs, Vec<DIM> & asum)
      : coefs(acoefs), sum(asum) { ; }

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
      // cout << "not implemented 63601" << endl;
      
      AutoDiff<3> hv =
        uvw.u.Value() * Cross (uvw.v, uvw.w) +
        uvw.v.Value() * Cross (uvw.w, uvw.u) +
        uvw.w.Value() * Cross (uvw.u, uvw.v);
      
      for (int i = 0; i < 3; i++)
        sum(i) += (*coefs) * hv.DValue(i);
    }

    void operator= (const Du_Cross_Dv<DIM> & uv) 
    { 
      // cout << "not implemented 63602" << endl;
      
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
    double * dshape;
  public:
    HDivDivShapeAssign (FlatVector<>  mat)
    { dshape = &mat(0); }

    HDivDivShapeElement<DIM> operator[] (int i) const
    { return HDivDivShapeElement<DIM> (dshape + i); }
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
			    FlatMatrixFixWidth<DIM> shape) const
    {    
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
	adp[i] = AutoDiff<DIM> (ip(i), i);

      HDivShapeAssign<DIM> ds(shape); 
      static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);
    }


    
    virtual void CalcDivShape (const IntegrationPoint & ip, 
			       FlatVector<> divshape) const
    {  
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
	adp[i] = AutoDiff<DIM> (ip(i), i);

      HDivDivShapeAssign<DIM> ds(divshape); 
      static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);
    }

    
    virtual void CalcMappedShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
				  FlatMatrixFixWidth<DIM> shape) const
    {   
      AutoDiff<DIM> adp[DIM];

      for (int i = 0; i < DIM; i++)
	adp[i].Value() = sip.IP()(i);

      for (int i = 0; i < DIM; i++)
	for (int j = 0; j < DIM; j++)
	  adp[i].DValue(j) = sip.GetJacobianInverse()(i,j);

      HDivShapeAssign<DIM> ds(shape); 
      static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);

      /*
      MatrixFixWidth<DIM> hshape(shape.Height());
      HDivFiniteElement<DIM>::CalcMappedShape (sip, hshape);
      *testout << "mapped = " << endl << shape << endl;
      *testout << "base = " << endl << hshape << endl;
      *testout << "diff = " << endl << shape-hshape << endl;
      */
    }

      
    
    virtual void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, 
			   FlatMatrixFixWidth<DIM> vals) const
    {    
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < ir.GetNIP(); i++)
	{
	  for (int j = 0; j < DIM; j++)
	    adp[j] = AutoDiff<DIM> (ir[i](j), j);
	
	  HDivEvaluateShape<DIM> ds(coefs);
	  static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);
	  vals.Row(i) = ds.Sum(); 
	}
    }

    
  };

}


#endif
