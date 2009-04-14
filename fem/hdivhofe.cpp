
#include <fem.hpp>
#define DEBUG

namespace ngfem
{  
  using namespace ngfem;

  
  // hv.DValue() = (grad u) x (grad v) 
  AutoDiff<3> GradCrossGrad (const AutoDiff<3> & u,
                             const AutoDiff<3> & v)
  {
    AutoDiff<3> hv;
    hv.Value() = 0.0;
    hv.DValue(0) = u.DValue(1)*v.DValue(2)-u.DValue(2)*v.DValue(1);
    hv.DValue(1) = u.DValue(2)*v.DValue(0)-u.DValue(0)*v.DValue(2);
    hv.DValue(2) = u.DValue(0)*v.DValue(1)-u.DValue(1)*v.DValue(0);
    return hv;
  }
  
  Vec<3> GradCrossVec (const AutoDiff<3> & u,
                       const Vec<3> & v)
  {
    Vec<3> hv;
    hv(0) = u.DValue(1)*v(2)-u.DValue(2)*v(1);
    hv(1) = u.DValue(2)*v(0)-u.DValue(0)*v(2);
    hv(2) = u.DValue(0)*v(1)-u.DValue(1)*v(0);
    return hv;
  }
  
  
  
  // 2d (gradu) x (grad v)
  double  GradCrossGrad2 (const AutoDiff<2> & u, const AutoDiff<2> & v) 
  {
    return(u.DValue(0)*v.DValue(1) - u.DValue(2)*v.DValue(1)); 
  }  
  
 
   
  //------------------------------------------------------------------------
  // HDivHighOrderFiniteElement
  //------------------------------------------------------------------------
  template <int D>
  HDivHighOrderFiniteElement<D> ::
  HDivHighOrderFiniteElement (ELEMENT_TYPE aeltype)
    : HDivFiniteElement<D> (aeltype, -1, -1)
  {
    for (int i = 0; i < 8; i++)
      vnums[i] = i;

    ho_div_free = 0;
  }

  template <int D>
  void HDivHighOrderFiniteElement<D>::
  SetVertexNumbers (FlatArray<int> & avnums)
  {
    for (int i = 0; i < avnums.Size(); i++)
      vnums[i] = avnums[i];
  }

  template <int D>
  void HDivHighOrderFiniteElement<D>::
  SetOrderEdge (FlatArray<int> & oe)
  {
    for (int i = 0; i < oe.Size(); i++)
      order_edge[i] = oe[i];
    ComputeNDof();
  }

  template <int D>
  void HDivHighOrderFiniteElement<D>::
  SetOrderFace (FlatArray<int> & of)
  {
    for (int i = 0; i < of.Size(); i++)
      order_face[i] = INT<2>(of[i],of[i]);
    ComputeNDof();
   
  }
  
  template <int D>
  void HDivHighOrderFiniteElement<D>::
  SetOrderFace (FlatArray<INT<2> > & of)
  {
    for (int i = 0; i < of.Size(); i++)
      order_face[i] = of[i];
    ComputeNDof();
  }

  template <int D>
  void HDivHighOrderFiniteElement<D>::
  SetOrderInner (int oi)
  {
    order_inner = INT<3>(oi,oi,oi);
    ComputeNDof();

  }
  
  template <int D>
  void HDivHighOrderFiniteElement<D>::
  SetOrderInner (INT<3> oi)
  {
    order_inner = oi;
    ComputeNDof();
  }


  //------------------------------------------------------------------------
  // HDivHighOrderNormalFiniteElement
  //------------------------------------------------------------------------
  template <int D>
  HDivHighOrderNormalFiniteElement<D> ::
  HDivHighOrderNormalFiniteElement (ELEMENT_TYPE aeltype)
    : HDivNormalFiniteElement<D> (aeltype, -1, -1)
  {
    for (int i = 0; i < 4; i++)
      vnums[i] = i;
  }

  template <int D>
  void HDivHighOrderNormalFiniteElement<D>::
  SetVertexNumbers (FlatArray<int> & avnums)
  {

    for (int i = 0; i < avnums.Size(); i++)
      vnums[i] = avnums[i];
  }

  template <int D>
  void HDivHighOrderNormalFiniteElement<D>::
  SetOrderInner (int oi)
  {
    order_inner = INT<2>(oi,oi);
  }
  
  template <int D>
  void HDivHighOrderNormalFiniteElement<D>::
  SetOrderInner (INT<2> oi)
  {
    order_inner = oi;
  }

  //------------------------------------------------------------------------
  // HDivHighOrderNormalSegm
  //------------------------------------------------------------------------
  template <class T_ORTHOPOL>
  HDivHighOrderNormalSegm<T_ORTHOPOL> :: HDivHighOrderNormalSegm (int aorder)
    : HDivHighOrderNormalFiniteElement<1>(ET_SEGM)
  {
    order_inner = INT<2>(aorder,aorder);
    ComputeNDof();
  }

  template <class T_ORTHOPOL>
  void HDivHighOrderNormalSegm<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = order_inner[0] + 1;
    order = order_inner[0];
  }

  template <class T_ORTHOPOL>
  void HDivHighOrderNormalSegm<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip,
							 FlatVector<> shape) const
  {
    double x = ip(0);

    int p = order_inner[0];

    int fac = 1;
    if (vnums[0] > vnums[1])
      {
	fac *= -1;
	x=1-x;
      }


    ArrayMem<double, 10> rec_pol(p+1), leg(2);

    MatrixFixWidth<2> DExt(p);
    T_ORTHOPOL::CalcTrigExtDeriv(p+1, 1-2*x, 0, DExt);

    //IntegratedLegendrePolynomial (order_inner, 1-2*x, rec_pol);
    // LegendrePolynomial (order_inner, 2*x-1, rec_pol);
    //T_ORTHOPOL::CalcDeriv(order_inner+1, 1-2*x, rec_pol);
    shape(0) = -1.*fac;
    for(int j=1; j<=p;j++)
      shape(j) = 2.* fac * DExt(j-1,0);

  }
  template class HDivHighOrderNormalSegm<IntegratedLegendreMonomialExt>;
  template class HDivHighOrderNormalSegm<TrigExtensionMonomial>;
  //template class HDivHighOrderNormalSegm<TrigExtensionOptimal>;
  //template class HDivHighOrderNormalSegm<TrigExtensionMin>;

  //------------------------------------------------------------------------
  // HDivHighOrderNormalQuad
  //------------------------------------------------------------------------
  template <class T_ORTHOPOL>
  HDivHighOrderNormalQuad<T_ORTHOPOL> :: HDivHighOrderNormalQuad (int aorder)
    : HDivHighOrderNormalFiniteElement<2>(ET_QUAD)
  {
    order_inner = INT<2>(aorder,aorder);
    ComputeNDof();
  }

  template <class T_ORTHOPOL>
  void HDivHighOrderNormalQuad<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = (order_inner[0] < 0) ? 0 : (1 + order_inner[0]*order_inner[1] + order_inner[0] + order_inner[1]);
    order = max(order_inner[0],order_inner[1]);
    order++; // order used for numerical integration
  }

  template <class T_ORTHOPOL>
  void HDivHighOrderNormalQuad<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip,
							 FlatVector<> shape) const
  {
    int i, j, k, l, m, ii;
    AutoDiff<2> x (ip(0), 0);
    AutoDiff<2> y (ip(1), 1);

    AutoDiff<2> lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};
    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};

    shape = 0.0;

    ii = 1;

    INT<2> p = order_inner;
    int pp = max(p[0],p[1]); 
    
    ArrayMem<AutoDiff<2>,20> pol_xi(p[0]+1), pol_eta(p[1]+1);
    AutoDiff<2> polprod;

    int fmax = 0;
    for (int j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax])
	fmax = j;

    int f1 = (fmax+3)%4;
    int f2 = (fmax+1)%4;

    int fac = 1;
    if(vnums[f2] > vnums[f1])
      {
	swap(f1,f2); // fmax > f1 > f2;
	fac *= -1;
      }

    AutoDiff<2> xi  = sigma[fmax]-sigma[f1];
    AutoDiff<2> eta = sigma[fmax]-sigma[f2];

    shape(0) = fac;

    T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
    T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);

    // Typ 1
    for (k = 0; k < p[0]; k++)
      {
	for (l = 0; l < p[1]; l++, ii++)
	  {
            shape(ii) = 2.*(pol_eta[l].DValue(0)*pol_xi[k].DValue(1)-pol_eta[l].DValue(1)*pol_xi[k].DValue(0));
	  }
      }

    //Typ 2

    for (k = 0; k < p[0]; k++)
      shape(ii++) = -eta.DValue(0)*pol_xi[k].DValue(1) + eta.DValue(1)*pol_xi[k].DValue(0); 
    for (k = 0; k < p[1]; k++)
      shape(ii++)   = -xi.DValue(0)*pol_eta[k].DValue(1) + xi.DValue(1)*pol_eta[k].DValue(0);

    return;
  }

  template class HDivHighOrderNormalQuad<IntegratedLegendreMonomialExt>;
  template class HDivHighOrderNormalQuad<TrigExtensionMonomial>;
  // template class HDivHighOrderNormalQuad<TrigExtensionOptimal>;
  template class HDivHighOrderNormalQuad<TrigExtensionMin>;

  //------------------------------------------------------------------------
  // HDivHighOrderNormalTrig
  //------------------------------------------------------------------------

  template <class T_ORTHOPOL>
  HDivHighOrderNormalTrig<T_ORTHOPOL> :: HDivHighOrderNormalTrig (int aorder)
    : HDivHighOrderNormalFiniteElement<2>(ET_TRIG)
  {
    order_inner = INT<2>(aorder,aorder);
    ComputeNDof();
  }

  template <class T_ORTHOPOL>
  void HDivHighOrderNormalTrig<T_ORTHOPOL> :: ComputeNDof()
  {

    ndof = 1 + (order_inner[0]*order_inner[0]+3*order_inner[0])/2;
    order = order_inner[0];
    order++;
  }

  template <class T_ORTHOPOL>
  void HDivHighOrderNormalTrig<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip,
							 FlatVector<> shape) const
  {
    double x = ip(0);
    double y = ip(1);

    double lami[3];

    lami[0] = x;
    lami[1] = y;
    lami[2] = 1-x-y;

    Mat<3,2> dlami(0);
    dlami(0,0) = 1.;
    dlami(1,1) = 1.;
    dlami(2,0) = -1.;
    dlami(2,1) = -1.;

    int ii, i, j, k, l, is, ie, iop;

    int p = order_inner[0];
    ii = 1;

    int fav[3];
    for(i=0;i<3;i++) fav[i] = i;


    //Sort vertices  first edge op minimal vertex
    int fswap = 1;
    if(vnums[fav[0]] > vnums[fav[1]]) { swap(fav[0],fav[1]); fswap *= -1; }
    if(vnums[fav[1]] > vnums[fav[2]]) { swap(fav[1],fav[2]); fswap *= -1; }
    if(vnums[fav[0]] > vnums[fav[1]]) { swap(fav[0],fav[1]); fswap *= -1; }

    is = fav[0]; ie = fav[1]; iop = fav[2];

    AutoDiff<2> ls = lami[is];
    AutoDiff<2> le = lami[ie];
    AutoDiff<2> lo = lami[iop];

    //AutoDiff<3> lsle = lami[is]*lami[ie];
    for (j = 0; j < 2; j++)
      {
        ls.DValue(j) = dlami(is,j);
        le.DValue(j) = dlami(ie,j);
        lo.DValue(j) = dlami(iop,j);
      }

    Vec<2> nedelec;
    for (j = 0; j < 2; j++)
      nedelec(j) = ls.Value()*le.DValue(j) - le.Value()*ls.DValue(j);

    AutoDiff<2> lsle=ls*le;

    // RT_0-normal low order shapes
    shape(0) = fswap;

    Vec<2> grad1, grad2;
    ArrayMem<AutoDiff<2>, 100> ad_rec_pol1(100);
    ArrayMem<AutoDiff<2>, 100> ad_rec_pol2(100);

    ScaledLegendrePolynomial(p-1, le-ls, 1-lo, ad_rec_pol1);
    LegendrePolynomial(p-1, 2*lo-1, ad_rec_pol2);

    for (k = 0; k <= p-1; k++)
      {
	ad_rec_pol1[k] *= lsle;
	ad_rec_pol2[k] *= lo;
      }
    // Typ 1
    for (k = 0; k <= p-1; k++)
      {
	for (l = 0; l <= p-1-k; l++, ii++)
	  {
	    for (j = 0; j < 2; j++)
	      {
		grad1(j) = ad_rec_pol1[k].DValue(j);
		grad2(j) = ad_rec_pol2[l].DValue(j);
	      }
	    shape(ii) = 2. * (grad1(1)*grad2(0) - grad1(0)*grad2(1));

	  }
      }

    // Typ 2
    double curlned;
    curlned = 2.* (ls.DValue(0)*le.DValue(1) - ls.DValue(1)*le.DValue(0));
    for (k = 0; k <= p-1; k++, ii++)
      {
        for (j = 0; j < 2; j++)
          grad2(j) = ad_rec_pol2[k].DValue(j);
        shape(ii) = (grad2(0)*nedelec(1) - grad2(1)*nedelec(0)) + ad_rec_pol2[k].Value()*curlned;
      }       
  }
  template class HDivHighOrderNormalTrig<IntegratedLegendreMonomialExt>;
  template class HDivHighOrderNormalTrig<TrigExtensionMonomial>;
  // template class HDivHighOrderNormalTrig<TrigExtensionOptimal>;
  template class HDivHighOrderNormalTrig<TrigExtensionMin>;











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


  template <int D>
  inline double Dot (const AutoDiff<D> & u, const AutoDiff<D> & v)
  {
    double sum = 0;
    for (int i = 0; i < D; i++)
      sum += u.DValue(i) * v.DValue(i);
    return sum;
  }




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
      
      // for (int i = 0; i < DIM; i++) 
      // data[i] = uv.u.Value() * uv.v.DValue(i); 
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

  /*
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
  */


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


  /*
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

  */





  template <ELEMENT_TYPE ET>
  void T_HDivHighOrderFiniteElement<ET> :: 
  CalcShape (const IntegrationPoint & ip, FlatMatrixFixWidth<DIM> shape) const
  {    
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    HDivShapeAssign<DIM> ds(shape); 
    static_cast<const HDivHighOrderFE<ET>*> (this) -> T_CalcShape (adp, ds);
  }


  template <ELEMENT_TYPE ET>
  void T_HDivHighOrderFiniteElement<ET> :: 
  CalcDivShape (const IntegrationPoint & ip, FlatVector<> shape) const
  {  
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    HDivDivShapeAssign<DIM> ds(shape); 
    static_cast<const HDivHighOrderFE<ET>*> (this) -> T_CalcShape (adp, ds);
  }



















  //------------------------------------------------------------------------
  // HDivHighOrderTrig
  //------------------------------------------------------------------------

  HDivHighOrderFE<ET_TRIG> :: HDivHighOrderFE (int aorder)
  // : T_HDivHighOrderFiniteElement<ET_TRIG> ()>(ET_TRIG)
  {
    order_inner = INT<3>(aorder,0,0);
    for (int i = 0; i < 3; i++)
      order_edge[i] = aorder;

    ComputeNDof();
  }

  void HDivHighOrderFE<ET_TRIG> :: ComputeNDof()
  {
    // ndof_edge =0;
    ndof = 3; // Thomas-Raviart 
    // if oder_edge < 1 --> Thomas Raviart
    int i;

    for (i = 0; i < 3; i++)
      {  
	ndof += order_edge[i];
	// ndof_edge += order_edge[i];
      }
    
    if (order_inner[0] > 1)
      { 
        if (ho_div_free)
          ndof += order_inner[0]*(order_inner[0]-1)/2;
        else
          ndof += order_inner[0]*order_inner[0]-1;
	//ndof += order_inner_curl*(order_inner_curl-1)/2;
	//ndof += order_inner*(order_inner-1)/2 + order_inner-1;
      }

    order = 0; // max(order_edges_normal,order_inner);
    for (i = 0; i < 3; i++)
      {
	if (order_edge[i] > order)
	  order = order_edge[i];
      }

    if (order_inner[0] > order) 
      order = order_inner[0];

    order++;
  }


  void HDivHighOrderFE<ET_TRIG> ::
  GetInternalDofs (Array<int> & idofs) const
  {
    if (discontinuous)
      {
        idofs.SetSize(0);
        for (int i=0; i < ndof; i++)
          idofs.Append(i);
        return ;
      }
    else
      {
        idofs.SetSize (0);

        int base = 3;
        for (int i = 0; i < 3; i++)
          {
            base += order_edge[i];
          }


        if(order_inner[0] > 1)
          {
            int p = order_inner[0];
            int ni = p*p -1 ;
            for (int i = 0; i < ni; i++)
              idofs.Append (base+i);
          }
        //(*testout) << "idofs = " << idofs << endl;
      }
  }
  

  void HDivHighOrderFE<ET_TRIG> :: GetFacetDofs(int fa, Array<int> & dnums) const 
  {
    if (fa > 2 ) 
      {
        cout << " Warning HDIVHighOrderTrigSZ::GetFacetDofNrs() index out of range" << endl; 
      
        dnums.SetSize(0); 
        return; 
      } 
  
    dnums.SetSize(0); 
    dnums.Append(fa); 
    
    int ii = 3; // Thomas-Raviart 
    for (int i = 0; i < fa; i++)
      ii += order_edge[i];

    for(int i = 0; i < order_edge[fa]; i++)  
      dnums.Append(ii+i); 
  }                  



  template<typename Tx, typename TFA>  
  void  HDivHighOrderFE<ET_TRIG> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1];
    Tx lami[3] = { x, y, 1-x-y };

    ArrayMem<AutoDiff<2>,10> adpol1(order),adpol2(order);	
	
    int ii = 3; 
    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    for (int i = 0; i < 3; i++)
      {
	int es = edges[i][0], ee = edges[i][1];
	if (vnums[es] > vnums[ee])  swap (es, ee);

	//Nedelec low order edge shape function 
        shape[i] = uDv_minus_vDu<2> (lami[es], lami[ee]);

	int p = order_edge[i]; 
	//HO-Edge shapes (Gradient Fields)   
	if(p > 0) //  && usegrad_edge[i]) 
	  { 
	    AutoDiff<2> xi = lami[ee] - lami[es]; 
	    AutoDiff<2> eta = 1 - lami[ee] - lami[es]; 
	    T_ORTHOPOL::CalcTrigExt(p+1, xi, eta, adpol1); 
	   
	    for(int j = 0; j < p; j++) 
              shape[ii++] = Du<2> (adpol1[j]);
	  }
      }   

    //Inner shapes (Face) 
    int p = order_inner[0];      
    if(p > 1) 
      {
	int fav[3] = { 0, 1, 2 }; 
	//Sort vertices ... v(f0) < v(f1) < v(f2) 
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	  

	AutoDiff<2> xi  = lami[fav[2]]-lami[fav[1]];
	AutoDiff<2> eta = lami[fav[0]]; 

        TrigShapesInnerLegendre::CalcSplitted(p+1, xi, eta, adpol1,adpol2);
	
	// gradients:
        for (int j = 0; j < p-1; j++)
          for (int k = 0; k < p-1-j; k++, ii++)
            shape[ii] = Du<2> (adpol1[j] * adpol2[k]);
        
        
        if (!ho_div_free)
          {
            // other combination
            for (int j = 0; j < p-1; j++)
              for (int k = 0; k < p-1-j; k++, ii++)
                shape[ii] = uDv_minus_vDu<2> (adpol2[k], adpol1[j]);
            
            // rec_pol * Nedelec0 
            for (int j = 0; j < p-1; j++, ii++)
              shape[ii] = wuDv_minus_wvDu<2> (lami[fav[1]], lami[fav[2]], adpol2[j]);
          }
      }
  }



#ifdef HDIV_V2  
  void HDivHighOrderFE<ET_TRIG> :: CalcShape (const IntegrationPoint & ip,
                                              FlatMatrixFixWidth<2> shape) const
  {
    T_HDivHighOrderFiniteElement<ET_TRIG>::CalcShape (ip, shape);
    return;

    double x = ip(0);
    double y = ip(1);

    double lami[3];

    lami[0] = x;
    lami[1] = y;
    lami[2] = 1-x-y;


    Mat<3,2> dlami(0);
    dlami(0,0) = 1.;
    dlami(1,1) = 1.;
    dlami(2,0) = -1.;
    dlami(2,1) = -1.;

    // Curl(lami_i) 
    Mat<3,2> clami(0); 
    clami(0,1) = -1.;
    clami(1,0) = 1.;
    clami(2,1) = 1.;
    clami(2,0) = -1.;

    shape = 0.0;

    int p, p_curl;

    int i, j, k, l;
    int ii = 3;


    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);

    int is, ie, io;
    for (i = 0; i < 3; i++)
      {
	p = order_edge[i];
	is = edges[i][0];
	ie = edges[i][1];
	io = 3-is-ie;

	if (vnums[is] > vnums[ie])
	  swap (is, ie);

	//RT low order edge shape function
	for(j=0; j<2; j++)
	  shape(i,j) = lami[ie]*clami(is,j) - lami[is]*clami(ie,j);
	
	//Edge normal shapes
	if(p>0)
	  {
	    MatrixFixWidth<2> DExt(p+2);

	    T_ORTHOPOL::CalcTrigExtDeriv(p+1, lami[ie]-lami[is], lami[io], DExt);

	    for(j=0; j<= p-1;j++)
	      {
		shape(ii,0) = (clami(ie,0) - clami(is,0)) * DExt(j,0) + clami(io,0) * DExt(j,1);
		shape(ii++,1) = (clami(ie,1) - clami(is,1)) * DExt(j,0) + clami(io,1) * DExt(j,1);

	      }
	  }

      }

    if(order_inner[0] < 2) return;

    p = order_inner[0];
    p_curl = order_inner[0];

    //Inner shapes
    if(p>1)
      {

	ArrayMem<double,10> rec_pol(order-1), drec_polx(order-1),  drec_polt(order-1);
	ArrayMem<double,10> rec_pol2(order-1), drec_pol2x(order-1), drec_pol2t(order-1);
	ArrayMem<double,30> mem(3*order*sizeof(double));
	ArrayMem<double,30> mem2(3*order*sizeof(double));
	FlatMatrixFixWidth<3> curl_recpol(order, &mem[0]);
	FlatMatrixFixWidth<3> curl_recpol2(order, &mem2[0]);

	int fav[3];
	for(i=0;i<3;i++) fav[i] = i;

	//Sort vertices  first edge op minimal vertex
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);


	double ls = lami[fav[0]];
	double le = lami[fav[1]];
        double lo = lami[fav[2]];

	ScaledLegendrePolynomialandDiff(p_curl-2, le-ls, 1-lo,
					rec_pol, drec_polx, drec_polt);
	LegendrePolynomialandDiff(p_curl-2, 2*lo-1,rec_pol2, drec_pol2x);
	// \curl (ls * le * rec_pol), and \curl (lo rec_pol2)
	for (j = 0; j <= p_curl-2; j++)
	  {
	    for (l = 0; l < 2; l++)
	      {
		curl_recpol(j,l) = ls * le * (drec_polx[j] * (clami(fav[1],l) - clami(fav[0],l)) -
					      drec_polt[j] * (clami(fav[2],l))) +
		  (ls * clami(fav[1],l) + le * clami(fav[0],l)) * rec_pol[j];

		curl_recpol2(j,l) = lo * (drec_pol2x[j] * 2*clami(fav[2],l)) +
		  clami(fav[2],l) * rec_pol2[j];
	      }
	  }

	// curls:
	for (j = 0; j <= p_curl-2; j++)
	  for (k = 0; k <= p_curl-2-j; k++, ii++)
	    for (l = 0; l < 2; l++)
	      shape(ii, l) = ls*le*rec_pol[j]*curl_recpol2(k,l) + curl_recpol(j,l)*lo*rec_pol2[k];

        if (!ho_div_free)
          {
            // rotations of curls
            for (j = 0; j <= p-2; j++)
              for (k = 0; k <= p-2-j; k++, ii++)
                for (l = 0; l < 2; l++)
                  shape(ii, l) = ls*le*rec_pol[j]*curl_recpol2(k,l) - curl_recpol(j,l) * lo * rec_pol2[k];
            
            // rec_pol2 * RT_0
            for (j = 0; j <= p-2; j++, ii++)
              for (l = 0; l < 2; l++)
                shape(ii,l) = lo*rec_pol2[j] * (ls*clami(fav[1],l)-le*clami(fav[0],l));
          }
      }
  }

  void HDivHighOrderFE<ET_TRIG> :: CalcDivShape (const IntegrationPoint & ip,
                                                 FlatVector<> shape) const
  {
    T_HDivHighOrderFiniteElement<ET_TRIG>::CalcDivShape (ip, shape);
    return;
    cout << "divshape,new = " << endl << shape << endl;


    double lami[3];
    lami[0] = ip(0);
    lami[1] = ip(1);
    lami[2] = 1-ip(0)-ip(1);

    Mat<3,2> dlami(0.);
    dlami(0,0) = 1.;
    dlami(1,1) = 1.;
    dlami(2,0) = -1.;
    dlami(2,1) = -1.;



    Mat<3,2> clami(0.);
    clami(0,1) = -1.;
    clami(1,0) = 1.;
    clami(2,1) = 1.;
    clami(2,0) = -1.;

    shape = 0.0;

    int p, p_curl;

    int i, j, k, l;
    int ii = 3; 



    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    for (i = 0; i < 3; i++)
      {
        p = order_edge[i];

	int is = edges[i][0]; 
	int ie = edges[i][1];
	
	if (vnums[is] > vnums[ie])
	  swap (is, ie);

	//RT low order edge shape function  
	//shape(i) = 2*(dlami(ie,0)*clami(is,0) +  dlami(ie,1)*clami(ie,1));
	for(j=0; j<= 1; j++)
	  shape(i) += dlami(ie,j)*clami(is,j) - dlami(is,j)*clami(ie,j);
	

	//Edge normal shapes (curls) 
	if(p>0)
	  {
	    for(j=0; j<= p-1;j++)
	      ii++; 
	  }
      }
    
    if(order_inner[0] < 2) return;
    
    p = order_inner[0];
    p_curl = order_inner[0];
        
    //Inner shapes (Face) 
    if(p>1) 
      {
	
	ArrayMem<double,10> rec_pol(order-1), drec_polx(order-1),  drec_polt(order-1);
	ArrayMem<double,10> rec_pol2(order-1), drec_pol2x(order-1), drec_pol2t(order-1);  
	ArrayMem<double,30> mem(3*order*sizeof(double));
	ArrayMem<double,30> mem2(3*order*sizeof(double));
	FlatMatrixFixWidth<2> grad_recpol(order, &mem[0]);
	FlatMatrixFixWidth<2> grad_recpol2(order, &mem2[0]);

	int fav[3];
	for(i=0;i<3;i++) fav[i] = i;

	//Sort vertices  first edge op minimal vertex
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);

	
	double ls = lami[fav[0]]; 
	double le = lami[fav[1]]; 
       	double lo = lami[fav[2]];  
	
	ScaledLegendrePolynomialandDiff(p_curl-2, le-ls, 1-lo,
					rec_pol, drec_polx, drec_polt);  	    
	LegendrePolynomialandDiff(p_curl-2, 2*lo-1,rec_pol2, drec_pol2x);
	// \curl (ls * le * rec_pol), and \curl (lo rec_pol2)
	for (j = 0; j <= p_curl-2; j++)
	  {
	    for (l = 0; l < 2; l++)
	      {
		grad_recpol(j,l) =
		  ls * le * (drec_polx[j] * (dlami(fav[1],l) - dlami(fav[0],l)) -
			     drec_polt[j] * (dlami(fav[2],l))) +
		  (ls * dlami(fav[1],l) + le * dlami(fav[0],l)) * rec_pol[j];
		grad_recpol2(j,l) =
		  lo * (drec_pol2x[j] * 2*dlami(fav[2],l)) + 
		  dlami(fav[2],l) * rec_pol2[j];
	      }
	  }



	// curls:
	for (j = 0; j <= p_curl-2; j++)
	  for (k = 0; k <= p_curl-2-j; k++)
	    ii++;
	     
        if (!ho_div_free)
          {
            // rotations of curls  
            for (j = 0; j <= p-2; j++)
              for (k = 0; k <= p-2-j; k++, ii++)
                for (l = 0; l < 2; l++)
                  shape(ii) = 2 * (grad_recpol(j,0) * grad_recpol2(k,1)
                                   - grad_recpol(j,1) * grad_recpol2(k,0));
            // shape (ii) = grad u . curl v = 2 (u_x v_y - u_y v_x) ; 
            
            // div(rec_pol * RT_0)  =
            double divrt0 = 2*(dlami(fav[0],0)*clami(fav[1],0) +  dlami(fav[0],1)*clami(fav[1],1));
            Vec<2> rt0;
            for (l = 0; l<2; l++)
              rt0(l) = lami[fav[0]]*clami(fav[1],l)-lami[fav[1]]*clami(fav[0],l);
            for (j = 0; j <= p-2; j++, ii++)
              shape(ii) = divrt0 * lo*rec_pol2[j] + rt0(0)*grad_recpol2(j,0) +
                rt0(1)*grad_recpol2(j,1); 
          }
      }

    cout << "shape,orig = " << endl << shape << endl;
  }
#endif

  /*
    template class HDivHighOrderTrig<IntegratedLegendreMonomialExt>;
    template class HDivHighOrderTrig<TrigExtensionMonomial>;
    //  template class HDivHighOrderTrig<TrigExtensionOptimal>;
    template class HDivHighOrderTrig<TrigExtensionMin>;
  */

  //------------------------------------------------------------------------
  // HDivHighOrderQuad
  //------------------------------------------------------------------------




  HDivHighOrderFE<ET_QUAD> :: HDivHighOrderFE (int aorder)
  // : HDivHighOrderFiniteElement<2>(ET_QUAD)
  {
    order_inner = INT<3>(aorder,aorder,0);
    for (int i = 0; i < 4; i++)
      order_edge[i] = aorder;
    ComputeNDof();
  }

  void HDivHighOrderFE<ET_QUAD> :: ComputeNDof()
  {
    int i;
    ndof = 4;

    for (i = 0; i < 4; i++)
      ndof += order_edge[i];

    INT<2> p = INT<2>(order_inner[0],order_inner[1]);
    int ni = ho_div_free ? 
      p[0]*p[1] : 2*p[0]*p[1] + p[0] + p[1];

    ndof += ni; // 2*order_inner[0]*order_inner[1]+order_inner[0]+order_inner[1];

    order = 0;
    for (i = 0; i < 4; i++)
      if (order_edge[i] > order)
	order = order_edge[i];
    if (order_inner[0] > order)
      order = order_inner[0];
    order++;
  }


  void HDivHighOrderFE<ET_QUAD> ::
  GetInternalDofs (Array<int> & idofs) const
  {
    if (discontinuous)
      {
        idofs.SetSize(0);
        for (int i=0; i<ndof; i++)
          idofs.Append(i);
        return ;
      }
    else
      {
        idofs.SetSize (0);

        int base = 4;
        for (int i = 0; i < 4; i++)
          {
            base += order_edge[i];
          }

        INT<2> p = INT<2>(order_inner[0],order_inner[1]);
        int ni = ho_div_free ? 
          p[0]*p[1] : 2*p[0]*p[1] + p[0] + p[1];

        for (int i = 0; i < ni; i++)
          idofs.Append (base+i);

        //(*testout) << "idofs = " << idofs << endl;
      }
  }
  
  void HDivHighOrderFE<ET_QUAD> :: GetFacetDofs(int fa, Array<int> & dnums) const
  {
    if (fa > 3 ) 
      {
        cout << " Warning HDIVHighOrderQuadSZ::GetFacetDofNrs() index out of range" << endl; 
      
        dnums.SetSize(0); 
        return; 
      } 

   
    dnums.SetSize(0); 
    dnums.Append(fa);  
   
    int ii = 4; // Thomas-Raviart 
    for (int i = 0; i < fa; i++)
      ii += order_edge[i];

    for(int i = 0; i < order_edge[fa]; i++)  
      dnums.Append(ii+i); 
   

  } 
  



  template<typename Tx, typename TFA>  
  void  HDivHighOrderFE<ET_QUAD> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1];

    AutoDiff<2> lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  

    int ii = 4;
    ArrayMem<AutoDiff<2>, 10> pol_xi(order+2), pol_eta(order+2);

    // edges
    const EDGE * edges = ElementTopology::GetEdges (ET_QUAD);
    for (int i = 0; i < 4; i++)
      {
	int p = order_edge[i]; 
	int es = edges[i][0], ee = edges[i][1];
	if (vnums[es] > vnums[ee]) swap (es, ee);

	AutoDiff<2> xi  = sigma[ee]-sigma[es];
	AutoDiff<2> lam_e = lami[ee]+lami[es];  // attention in [0,1]

	// Nedelec0-shapes
        shape[i] = uDv<2> (0.5 * lam_e, xi); 

	// High Order edges ... Gradient fields 
	// if(usegrad_edge[i])
        {
          T_ORTHOPOL::Calc (p+1, xi, pol_xi);  
          for (int j = 0; j < p; j++)
            shape[ii++] = Du<2> (pol_xi[j] * lam_e);
        }
      }
     
    // INT<2> p = order_face[0]; // (order_cell[0],order_cell[1]);
    INT<2> p (order_inner[0], order_inner[1]); // (order_cell[0],order_cell[1]);
    int fmax = 0; 
    for (int j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax])
	fmax = j;
    
    int f1 = (fmax+3)%4; 
    int f2 = (fmax+1)%4; 
    if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f2 > f1; 

    AutoDiff<2> xi = sigma[fmax]-sigma[f1];  // in [-1,1]
    AutoDiff<2> eta = sigma[fmax]-sigma[f2]; // in [-1,1]
    
    T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
    T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);
    
    //Gradient fields 
    // if(usegrad_face[0])
    for (int k = 0; k < p[0]; k++)
      for (int j= 0; j < p[1]; j++)
        shape[ii++] = Du<2> (pol_xi[k]*pol_eta[j]);
    
    if (!ho_div_free)
      {
        //Rotation of Gradient fields 
        for (int k = 0; k < p[0]; k++)
          for (int j= 0; j < p[1]; j++)
            shape[ii++] = uDv_minus_vDu<2> (pol_eta[j], pol_xi[k]);
        
        //Missing ones 
        for(int j = 0; j< p[0]; j++)
          shape[ii++] = uDv<2> (0.5*pol_xi[j], eta);
        
        for(int j = 0; j < p[1]; j++)
          shape[ii++] = uDv<2> (0.5*pol_eta[j], xi); 
      }
  }




  void HDivHighOrderFE<ET_QUAD> :: CalcShape (const IntegrationPoint & ip,
                                              FlatMatrixFixWidth<2> shape) const
  {
    T_HDivHighOrderFiniteElement<ET_QUAD>::CalcShape (ip, shape);
    return;
    cout << "shape,new = " << endl << shape << endl;


    AutoDiff<2> x (ip(0), 0);
    AutoDiff<2> y (ip(1), 1);

    AutoDiff<2> lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};
    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};

    AutoDiff<2> ext[4] = {1-y, y, 1-x, x};

    Mat<4,2> can(0);
    can(0,1) = -1.;
    can(1,1) = 1.;
    can(2,0) = -1.;
    can(3,0) = 1.;

    shape = 0.0;
    int ii = 4;

    ArrayMem<AutoDiff<2>,20> rec_pol(order+2);

    const EDGE * edges = ElementTopology::GetEdges (ET_QUAD);

    int is, ie;
    for (int i = 0; i < 4; i++)
      {
        int p = order_edge[i];
        
        is = edges[i][0];
        ie = edges[i][1];
        
        int fac=1;
        if (vnums[is] > vnums[ie])
          { swap (is, ie); fac *= -1;}

        AutoDiff<2> prod;
        AutoDiff<2> xi  = sigma[ie]-sigma[is]; // in [-1,1]
        AutoDiff<2> eta = lami[ie]+lami[is];  // attention in [0,1]

        T_ORTHOPOL::Calc(p+1, xi,rec_pol);

        //RT low order edge shape function
        for(int k=0; k<2; k++)
	  shape(i,k) = -1.*fac*ext[i].Value()*can(i,k);

	//Edge normal shapes
        for (int j=0; j < p; j++, ii++)
	  {
	    prod = eta*rec_pol[j];
	    shape(ii,0) = prod.DValue(1);
	    shape(ii,1) = -prod.DValue(0);
	  }

	/*
	  AutoDiff<2> ss = sigma[is];
	  AutoDiff<2> se = sigma[ie];
	  AutoDiff<2> xi  = sigma[ie]-sigma[is]; // in [-1,1]
	  AutoDiff<2> eta = lami[ie]+lami[is];  // attention in [0,1]
	  LegendrePolynomial(p+1, se-ss, rec_pol);
	  //RT low order edge shape function
	  for(int k=0; k<2; k++)
	  shape(i,k) = -1.*fac*ext[i].Value()*can(i,k);
	  // Edge normal shapes
	  for (int j=1; j<=p; j++, ii++)
	  for (int k=0; k<2; k++)
	  shape(ii,k) = 2.*fac*rec_pol[j].Value()*ext[i].Value()*can(i,k);
	*/
      }


    // Inner shapes
    INT<2> p = INT<2>(order_inner[0],order_inner[1]);

    ArrayMem<AutoDiff<2>,20> polxi(p[0]+1), poleta(p[1]+1);
    AutoDiff<2> polprod;

    int fmax = 0;
    for (int j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax])
	fmax = j;

    int f1 = (fmax+3)%4;
    int f2 = (fmax+1)%4;


    if(vnums[f2] > vnums[f1])
      swap(f1,f2);  // fmax > f1 > f2;


    AutoDiff<2> xi  = sigma[fmax]-sigma[f1];
    AutoDiff<2> eta = sigma[fmax]-sigma[f2];

    //LegendrePolynomial(p+1, xi, polxi);
    //LegendrePolynomial(p+1, eta, poleta);

    T_ORTHOPOL::Calc(p[0]+1, xi,polxi);
    T_ORTHOPOL::Calc(p[1]+1,eta,poleta);


    for (int i = 0; i < p[0]; i++)
      {
	for (int j = 0; j < p[1]; j++)
	  {
	    polprod = polxi[i]*poleta[j];
	    shape(ii,0) = polprod.DValue(1);
	    shape(ii++,1) = -polprod.DValue(0);
	  }
      }
    for (int i = 0; i < p[0]; i++)
      {
	for (int j = 0; j < p[1]; j++)
	  {
	    shape(ii,0) = (polxi[i].DValue(1)*poleta[j].Value() - polxi[i].Value()*poleta[j].DValue(1));
	    shape(ii++,1) = (-polxi[i].DValue(0)*poleta[j].Value() + polxi[i].Value()*poleta[j].DValue(0));
	  }
      }
    for (int i= 0; i < p[0]; i++, ii++)
      {
        shape(ii,0) = (eta.DValue(1)*polxi[i].Value());
        shape(ii,1) = (-eta.DValue(0)*polxi[i].Value());
      }
    for (int i= 0; i < p[1]; i++, ii++)
      {
        shape(ii,0) = (xi.DValue(1)*poleta[i].Value());
        shape(ii,1) = (-xi.DValue(0)*poleta[i].Value());
      }


    cout << "shape,old = " << endl << shape << endl;
    return;
  }


  void HDivHighOrderFE<ET_QUAD> :: CalcDivShape (const IntegrationPoint & ip,
                                                 FlatVector<> shape ) const
  {
    T_HDivHighOrderFiniteElement<ET_QUAD>::CalcDivShape (ip, shape);
    return;


    AutoDiff<2> x (ip(0), 0);
    AutoDiff<2> y (ip(1), 1);


    AutoDiff<2> lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};
    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};
    AutoDiff<2> ext[4] = {(1-y), y, (1-x), x};
    int ind[4] = {1, 1, 0, 0};
    int can[4] = {-1, 1, -1, 1};
    shape = 0.0;


    int ii = 4;

    ArrayMem<AutoDiff<2>,10> rec_pol(order+1);
    AutoDiff<2> pol_ext;

    const EDGE * edges = ElementTopology::GetEdges (ET_QUAD);

    int is, ie;
    for (int i = 0; i < 4; i++)
      {
	int p = order_edge[i];
	is = edges[i][0];
	ie = edges[i][1];

        int fac=1;
	if (vnums[is] > vnums[ie])
	  {  swap (is, ie); fac*=-1;
	  }

        AutoDiff<2> ss = sigma[is];
	AutoDiff<2> se = sigma[ie];

        LegendrePolynomial(p+1, se-ss, rec_pol);

        //RT low order edge shape function

	shape(i) = -1.*fac*ext[i].DValue(ind[i])*can[i];

	//Edge normal shapes
	/*for (int j=1; j<=p; j++)
	  {
	  pol_ext = rec_pol[j]*ext[i];
	  shape(ii++) = 2.*fac*pol_ext.DValue(ind[i])*can[i];
	  //shape(ii++) = rec_pol[j].DValue(ind[i]);

	  }*/
	ii += p;

      }
    // Inner shapes
    INT<2> p = INT<2> (order_inner[0],order_inner[1]);

    ArrayMem<AutoDiff<2>,20> polxi(p[0]+1), poleta(p[1]+1);
    AutoDiff<2> polprod;

    int fmax = 0;
    for (int j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax])
	fmax = j;

    int f1 = (fmax+3)%4;
    int f2 = (fmax+1)%4;

    if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f1 > f2;

    AutoDiff<2> xi  = sigma[fmax]-sigma[f1];
    AutoDiff<2> eta = sigma[fmax]-sigma[f2];


    //LegendrePolynomial(p+1, xi, polxi);
    //LegendrePolynomial(p+1, eta, poleta);
    T_ORTHOPOL::Calc(p[0]+1, xi,polxi);
    T_ORTHOPOL::Calc(p[1]+1,eta,poleta);

    ii += p[0]*p[1];

    for (int i = 0; i < p[0]; i++)
      {
	for (int j = 0; j < p[1]; j++)
	  {
	    shape(ii++) = (2*(polxi[i].DValue(1)*poleta[j].DValue(0) - polxi[i].DValue(0)*poleta[j].DValue(1)));
          }
      }
    for (int i= 0; i < p[0]; i++, ii++)
      shape(ii) = (eta.DValue(1)*polxi[i].DValue(0)-eta.DValue(0)*polxi[i].DValue(1));
    for (int i= 0; i < p[1]; i++, ii++)
      shape(ii) = (xi.DValue(1)*poleta[i].DValue(0)-xi.DValue(0)*poleta[i].DValue(1));
    return;

  }


#ifdef V2  
  void HDivHighOrderFE<ET_QUAD> :: CalcNumDivShape (const IntegrationPoint & ip,
                                                    FlatVector<> divshape) const
  {
    double x=ip(0);
    double y=ip(1);
    LocalHeap lh(10000000);
    FlatMatrixFixWidth<2> cshape(ndof, lh);
    FlatMatrixFixWidth<2> cshape1(ndof, lh);
    FlatMatrixFixWidth<2> cshape_(ndof, lh);
    FlatMatrixFixWidth<2> cshape1_(ndof, lh);
    FlatMatrixFixWidth<2> cshape2(ndof, lh);
    FlatMatrixFixWidth<2> cshape3(ndof, lh);
    FlatMatrixFixWidth<2> cshape2_(ndof, lh);
    FlatMatrixFixWidth<2> cshape3_(ndof, lh);

    double h = 0.001;
    IntegrationPoint ipp(x+h,y,0,0);
    IntegrationPoint ipp1(x-h,y,0,0);
    IntegrationPoint ipp_(x+2*h,y,0,0);
    IntegrationPoint ipp1_(x-2*h,y,0,0);
    IntegrationPoint ipp2(x,y+h,0,0);
    IntegrationPoint ipp3(x,y-h,0,0);
    IntegrationPoint ipp2_(x,y+2*h,0,0);
    IntegrationPoint ipp3_(x,y-2*h,0,0);

    CalcShape(ipp,cshape);
    CalcShape(ipp1,cshape1);
    CalcShape(ipp_,cshape_);
    CalcShape(ipp1_,cshape1_);
    CalcShape(ipp2,cshape2);
    CalcShape(ipp3,cshape3);
    CalcShape(ipp2_,cshape2_);
    CalcShape(ipp3_,cshape3_);

    for (int i = 0; i< ndof; i++)
      divshape(i) = 2/(3*h)*(cshape(i,0)-cshape1(i,0))-1/(12*h)*(cshape_(i,0)-cshape1_(i,0)) + 2/(3*h)*(cshape2(i,1)-cshape3(i,1))-1/(12*h)*(cshape2_(i,1)-cshape3_(i,1));

    return;
  }
#endif


  /*
    template class HDivHighOrderQuad<IntegratedLegendreMonomialExt>;
    template class HDivHighOrderQuad<TrigExtensionMonomial>;
    // template class HDivHighOrderQuad<TrigExtensionOptimal>;
    template class HDivHighOrderQuad<TrigExtensionMin>;
  */


  //------------------------------------------------------------------------
  // HDivHighOrderTet
  //------------------------------------------------------------------------


  HDivHighOrderFE<ET_TET> :: HDivHighOrderFE (int aorder)
  {
    order_inner = INT<3>(aorder,aorder,aorder);
    for (int i = 0; i < 4; i++)
      order_face[i] = INT<2>(aorder,aorder);

    ComputeNDof();
  }

 

  void HDivHighOrderFE<ET_TET> :: ComputeNDof()
  {
    ndof = 4;

    for (int i = 0; i < 4; i++)
      {
        if (order_face[i][0] > 0)
	  { 
	    int p = order_face[i][0];
	    ndof += (p*p+3*p)/2;
	  }
      }
    
    int p = order_inner[0];
    int pc = order_inner[0]; // should be order_inner_curl!!!  
    
    if(pc > 1) 
      ndof += pc*(pc+1)*(pc-1)/3 + pc*(pc-1)/2;
    if(p > 1 && !ho_div_free) 
      ndof += p*(p+1)*(p-1)/6 + p*(p-1)/2 + p-1;


    order = 0; // max(order_face_normal,order_inner);
    for (int i = 0; i < 4; i++)
      {
	if (order_face[i][0] > order)
	  order = order_face[i][0];
      }

    if (order_inner[0] > order)
      order = order_inner[0];

    // order++;
  }



  // template <class T_ORTHOPOL>
  void HDivHighOrderFE<ET_TET> ::
  GetInternalDofs (Array<int> & idofs) const
  {
    if (discontinuous)
      {
        idofs.SetSize(0);
        for (int i=0; i<ndof; i++)
          idofs.Append(i);
        return ;
      }
    else
      {
        idofs.SetSize (0);

        int base = 4;
        for (int i = 0; i < 4; i++)
          {
            int p = order_face[i][0];
            base += (p*p+3*p)/2;
          }


        if(order_inner[0] >= 2)
          {
            int p = order_inner[0];
            int pc = order_inner[0];

            //int ni = 6*(p-1) + 4*(p-1)*(p-2) + (p-1)*(p-2)*(p-3)/2;
            // int ni = p*(p+1)*(p-1)/2 + p*(p-1)+(p-1);

            int ni;
            if(pc > 1) 
              ni += pc*(pc+1)*(pc-1)/3 + pc*(pc-1)/2;
            if(p > 1 && !ho_div_free) 
              ni += p*(p+1)*(p-1)/6 + p*(p-1)/2 + p-1;
            
            for (int i = 0; i < ni; i++)
              idofs.Append (base+i);
          }
        //(*testout) << "idofs = " << idofs << endl;
      }
  }
  
  

  void HDivHighOrderFE<ET_TET> :: GetFacetDofs(int fa, Array<int> & dnums) const 
  {
    if (fa >= 4 ) 
      {
        cout << " Warning HDIVHighOrderTet::GetFacetDofNrs() index out of range" << endl; 
      
        dnums.SetSize(0); 
        return; 
      } 
  
    dnums.SetSize(0); 
    dnums.Append(fa);  // lowest-order
    
    int ii = 4; // Thomas-Raviart 
     
    for (int i = 0; i < fa; i++)
      {     
        int p = order_face[i][0];
        ii+= (p*p+3*p)/2; 
      }
    
    int p = order_face[fa][0];
    int nf = (p*p+3*p)/2;
    
    for(int i = 0; i < nf; i++)  
      dnums.Append(ii+i); 
  }                  





  template<typename Tx, typename TFA>  
  void  HDivHighOrderFE<ET_TET> :: T_CalcShape (Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];
    Tx lami[4] = { x, y, z, 1-x-y-z };

    ArrayMem<Tx,10> adpol1(order), adpol2(order), adpol3(order);
	
    int ii = 4; 
    const FACE * faces = ElementTopology::GetFaces (ET_TET);
    for (int i = 0; i < 4; i++)
      {
	int p = order_face[i][0];

        int fav[3];
        for(int j = 0; j < 3; j++) fav[j]=faces[i][j];
        
        //Sort vertices  first edge op minimal vertex
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        int fop = 6 - fav[0] - fav[1] - fav[2];
        
	// RT lowest order
        shape[i] = uDvDw_Cyclic<3> (lami[fav[0]], lami[fav[1]], lami[fav[2]]);

        Tx xi = lami[fav[1]]-lami[fav[0]];
        Tx eta = lami[fav[2]];
        Tx zeta = lami[fop];  
      
        T_FACESHAPES::CalcSplitted (p+2, xi, eta, zeta, adpol1, adpol2); 

        // Compability WITH TRIG!! 
        for (int k = 0; k < adpol1.Size(); k++)
          adpol1[k] *= 0.5; 
          
        // Curl (Type 2) 2*grad v x grad u
        for (int j = 0; j <= p-1; j++) 
          for (int k = 0; k <= p-1-j; k++)
            shape[ii++] = Du_Cross_Dv<3> (adpol2[k], adpol1[j]);

        // Curl (Type 3) //curl( * v) = nabla v x ned + curl(ned)*v
        for (int j = 0; j <= p-1; j++)
          shape[ii++] = curl_uDvw_minus_Duvw<3> (lami[fav[0]], lami[fav[1]], adpol2[j]);
      }

    
    // cell-based shapes 
    int p = order_inner[0];
    int pc = order_inner[0]; // should be order_inner_curl  
    int pp = max(p,pc); 
    if ( pp >= 2 )
      {
        T_INNERSHAPES::CalcSplitted(pp+2, lami[0]-lami[3], lami[1], lami[2], adpol1, adpol2, adpol3 );
      
        // Curl-Fields 
        for (int i = 0; i <= pc-2; i++)
          for (int j = 0; j <= pc-2-i; j++)
            for (int k = 0; k <= pc-2-i-j; k++)
              {
                // grad v  x  grad (uw)
                shape[ii++] = Du_Cross_Dv<3> (adpol2[j], adpol1[i]*adpol3[k]);
      
                // grad w  x  grad (uv)
                shape[ii++] = Du_Cross_Dv<3> (adpol3[k], adpol1[i]*adpol2[j]);
              }     


        // Type 1 : Curl(T3)
        // ned = lami[0] * nabla(lami[3]) - lami[3] * nabla(lami[0]) 
        for (int j= 0; j <= pc-2; j++)
          for (int k = 0; k <= pc-2-j; k++)
            shape[ii++] = curl_uDvw_minus_Duvw<3> (lami[0], lami[3], adpol2[j]*adpol3[k]);


        if (!ho_div_free)
          { 
            // Type 2:  
            // (grad u  x  grad v) w 
            for (int i = 0; i <= p-2; i++)
              for (int j = 0; j <= p-2-i; j++)
                for (int k = 0; k <= p-2-i-j; k++)
                  shape[ii++] = wDu_Cross_Dv<3> (adpol1[i], adpol2[j], adpol3[k]);

            // (ned0 x grad v) w    
            for (int j = 0; j <= p-2; j++)
              for (int k= 0; k <= p-2-j; k++)
                shape[ii++] = wDu_Cross_Dv<3> (lami[0], adpol2[j], lami[3]*adpol3[k]);
            
            // Type 3: 
            // (ned0 x e_z) v = (N_y, -N_x,0)^T * v ) 
            for (int j=0; j<=p-2; j++) 
              shape[ii++] = wDu_Cross_Dv<3> (lami[0], z, lami[3]*adpol2[j]);
          }
      }
  }








#define oldv
#ifdef oldv 


  void HDivHighOrderFE<ET_TET> :: CalcShape (const IntegrationPoint & ip,
                                             FlatMatrixFixWidth<3> shape) const
  {
    shape = 0.0;
    T_HDivHighOrderFiniteElement<ET_TET>::CalcShape (ip, shape);
    return;
    // *testout << "newshape = " << endl << shape << endl;


    double x = ip(0);
    double y = ip(1);
    double z = ip(2);


    int i, j, ii, k,l, m, n;

    int print=0;

    double lami[4];
    for (j=0;j<3;j++) lami[j] = ip(j);
    lami[3] = 1-ip(0)-ip(1)-ip(2);

    int index[3][3];

    int node[3];

    Mat<4,3> dlami;
    dlami = 0.0;
    dlami(0,0) = 1.;
    dlami(1,1) = 1.;
    dlami(2,2) = 1.;
    dlami(3,0) = -1.;
    dlami(3,1) = -1.;
    dlami(3,2) = -1.;

    shape = 0.0;

    const EDGE * edges = ElementTopology::GetEdges (ET_TET);
    const FACE * faces = ElementTopology::GetFaces (ET_TET);
    const POINT3D * vertices = ElementTopology::GetVertices (ET_TET);


    Vec<3> tang, tang1, normal, vp;

    int is, ie, iop, iop1, iop2;

    int p=(order+1)*(order+2)*(order+3)/2;
    

    ArrayMem<double,10> rec_pol(p), rec_pol_1(p), rec_pol_2(p);

    //********************************FACE SHAPES********************************************************
	// RT_0
	// Typ 1
	// Typ 2
    
	ii = 4;
    for (i = 0; i < 4; i++)
      {
	p = order_face[i][0];
	
	if (p >= 0)
	  {
            int fav[3];
	    for(j=0; j<3; j++) fav[j]=faces[i][j];

	    //Sort vertices  first edge op minimal vertex
	    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
	    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
	    int fop = 6 - fav[0] - fav[1] - fav[2];

            is = fav[0]; ie = fav[1]; iop = fav[2];

	    AutoDiff<3> ls = lami[is];
	    AutoDiff<3> le = lami[ie];
	    AutoDiff<3> lo = lami[iop];
	    AutoDiff<3> lz = lami[fop];
           

	    for (j = 0; j < 3; j++)
	      {
		ls.DValue(j) = dlami(is,j);
		le.DValue(j) = dlami(ie,j);
		lo.DValue(j) = dlami(iop,j);
		lz.DValue(j) = dlami(fop,j);

	      }
	    AutoDiff<3> lsle=ls*le;

            Vec<3> vp, vp1, vp2, vp3, nedelec;
	    Vec<3> dlamis, dlamie, dlamio;

	    for (n = 0; n < 3; n++)
	      {
		dlamis(n) = dlami(is,n);
		dlamie(n) = dlami(ie,n);
		dlamio(n) = dlami(iop,n);
		nedelec(n) = ls.Value()*dlamie(n) - le.Value()*dlamis(n);
	      }
	    vp1 = Cross (dlamis, dlamie);
	    vp2 = Cross (dlamio, dlamis);
	    vp3 = Cross (dlamie, dlamio);



	    // RT_0 low order shapes
	    for (j = 0; j < 3; j++)
	      shape(i,j) = ls.Value()*vp3(j)+le.Value()*vp2(j)+lo.Value()*vp1(j);


            Vec<3> grad1, grad2;
	    ArrayMem<AutoDiff<3>, 10> ad_rec_pol1(order+1);
	    ArrayMem<AutoDiff<3>, 10> ad_rec_pol2(order+1);

            // *testout << " p " << p << endl; 
            
          
	    ScaledLegendrePolynomial(p-1, le-ls, ls+le, ad_rec_pol1);
	    ScaledLegendrePolynomial(p-1, lo-le-ls, lo+ls+le, ad_rec_pol2);

            for (k = 0; k <= p-1; k++)
	      { 
		ad_rec_pol1[k] *= lsle;
		ad_rec_pol2[k] *= lo;
	      }
           
            // Typ 1
            for (k = 0; k <= p-1; k++)
	      {
		for (l = 0; l <= p-1-k; l++, ii++)
		  {

		    for (j = 0; j < 3; j++)
		      {
			grad1(j) = ad_rec_pol1[k].DValue(j);
			grad2(j) = ad_rec_pol2[l].DValue(j);
		      }
		    vp = Cross(grad2,grad1);
		    for (m = 0; m < 3; m++)
		      shape(ii,m) = 2.*vp(m);
		  }
	      }

	    // Typ 2
	    for (k = 0; k <= p-1; k++, ii++)
	      {
		for (j = 0; j < 3; j++)
                  grad2(j) = ad_rec_pol2[k].DValue(j);
		vp = Cross(nedelec,grad2);
		for (m = 0; m < 3; m++)
	          shape(ii,m) = (-1.*vp(m) + 2 * vp1(m) * ad_rec_pol2[k].Value());
	      }

	  }
      }
    


    // if(order_inner[0] < 2) return;
    // **********************************INNER*********************************************

    
      
    // **************curl inner shapes****************************************************
    p = order_inner[0];
    Vec<3> grad1, grad2, grad3, vp1, vp2, vp3;
    ArrayMem<AutoDiff<3>, 10> ad_poli(order-1);
    ArrayMem<AutoDiff<3>, 10> ad_polj(order-1);
    ArrayMem<AutoDiff<3>, 10> ad_polk(order-1);
    AutoDiff<3> pol_prod;
    AutoDiff<3> x1 (ip(0),0);
    AutoDiff<3> y1 (ip(1),1);
    AutoDiff<3> z1 (ip(2),2);
    
    //AutoDiff<3> lz = lami[fop];
    T_INNERSHAPES::CalcSplitted(p+2, x1-(1-x1-y1-z1), y1, z1,ad_poli, ad_polj, ad_polk );
    //   ScaledLegendrePolynomial(p-1, le-ls, ls+le, ad_rec_pol1);
    //   ScaledLegendrePolynomial(p-1, lo-le-ls, lo+ls+le, ad_rec_pol2);
    for (i = 0; i <= p-2; i++)
      {
	
	for (j = 0; j <= p-2-i; j++)
	  {
	    
	    for (k = 0; k <= p-2-i-j; k++, ii+=3)
	      {
		
		for (l = 0; l < 3; l++)
		  {
		    grad1(l) = ad_poli[i].DValue(l);
		    grad2(l) = ad_polj[j].DValue(l);
		    grad3(l) = ad_polk[k].DValue(l);
		  }
		
		vp1 = Cross(grad2,grad1);
		vp2 = Cross(grad2,grad3);
		vp3 = Cross(grad3,grad1);

		
		for (m = 0; m < 3; m++)
		  {
		    shape(ii,m)   = 2.* (ad_polk[k].Value()*vp1(m) + ad_poli[i].Value()*vp2(m));
		    shape(ii+1,m) = 2.* (ad_polj[j].Value()*vp3(m) - ad_poli[i].Value()*vp2(m));
		    shape(ii+2,m) = ad_polk[k].Value()*vp1(m);
		  }
              
             
            
	      }
	  }
      }
    AutoDiff<3> le = 1-x1-y1-z1;
    
    Vec<3> ned=0.;;
    ned(0) = - x1.Value() - le.Value(); 
    ned(1) =  - x1.Value(); 
    ned(2)= - x1.Value();
   
    Vec<3> curlned=0.;
    curlned(1) = (-x1.DValue(0)*le.DValue(2));
    curlned(2) = (x1.DValue(0)*le.DValue(1));
    // Typ 4  curl(Nedelec*ad_polk*ad_polj)
    for (i = 0; i <= p-2; i++)
      {
	for (j = 0; j <= p-2-i; j++, ii+=2)
	  {
	    pol_prod = ad_polj[j]*ad_polk[i];
	    for (k = 0; k < 3; k++)
	      {
		grad1(k) = pol_prod.DValue(k);
		grad2(k) = ad_polj[i].DValue(k);
	      }
	    vp1 = Cross(grad1,ned);
	    vp2 = Cross(ned,grad2);
	    
         
          
	    for (m = 0; m < 3; m++)
	      {
		shape(ii,m) = 2.*(curlned(m)*pol_prod.Value()) + vp1(m);
		shape(ii+1,m) = vp2(m)*ad_polk[j].Value();
	      }
                       
	  }
      }
    for (i = 0; i <= p-2; i++, ii++)
      {
	for (j = 0; j < 3; j++)
	  grad1(j) = z1.DValue(j);
        vp1 = Cross(ned,grad1);   
	for (m = 0; m < 3; m++)
	  {
	    shape(ii,m) = vp1(m)* ad_polj[i].Value();
          
	    
	  }
      }

    *testout << "oldshape = " << endl << shape << endl;
   
    // (*testout)<<"shape tet="<<shape<<endl;
    // cout << "shape, old = " << shape << endl;     
  }


  void HDivHighOrderFE<ET_TET> :: CalcDivShape (const IntegrationPoint & ip,
                                                FlatVector<> divshape) const
  {
    divshape = 0.0;
    T_HDivHighOrderFiniteElement<ET_TET>::CalcDivShape (ip, divshape);
    return;


    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    int i, j, ii, k, l, m, n;
    int is, ie, iop, iop1, iop2;




    double lami[4];
    for (j=0;j<3;j++) lami[j] = ip(j);
    lami[3] = 1-ip(0)-ip(1)-ip(2);

    Mat<4,3> dlami;
    dlami = 0.0;
    dlami(0,0) = 1.;
    dlami(1,1) = 1.;
    dlami(2,2) = 1.;
    dlami(3,0) = -1.;
    dlami(3,1) = -1.;
    dlami(3,2) = -1.;

    divshape = 0.0;

    const EDGE * edges = ElementTopology::GetEdges (ET_TET);
    const FACE * faces = ElementTopology::GetFaces (ET_TET);
    const POINT3D * vertices = ElementTopology::GetVertices (ET_TET);

    Vec<3> tang, tang1, normal, vp;

    int p=(order+1)*(order+2)*(order+3)/2;




    ArrayMem<double,10> rec_pol(p), rec_pol_1(p), rec_pol_2(p);
    ArrayMem<double,10> drec_pol(p), drec_pol_1(p), drec_pol_2(p);
    ArrayMem<AutoDiff<3>, 100> ad_rec_pol(100);


    //*******************************DIV FACE-SHAPES*********************************************



	// RT_0
	// Typ 1
	// Typ 2

	ii = 4;
    for (i = 0; i < 4; i++)
      {
	p = order_face[i][0];

	if (p >= 0)
	  {
            int fav[3];
	    for(j=0; j<3; j++) fav[j]=faces[i][j];

	    //Sort vertices  first edge op minimal vertex
	    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
	    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
	    int fop = 6 - fav[0] - fav[1] - fav[2];

            is = fav[0]; ie = fav[1]; iop = fav[2];

	    AutoDiff<3> ls = lami[is];
	    AutoDiff<3> le = lami[ie];
	    AutoDiff<3> lo = lami[iop];
	    AutoDiff<3> lz = lami[fop];


	    for (j = 0; j < 3; j++)
	      {
		ls.DValue(j) = dlami(is,j);
		le.DValue(j) = dlami(ie,j);
		lo.DValue(j) = dlami(iop,j);
		lz.DValue(j) = dlami(fop,j);

	      }
	    AutoDiff<3> lsle=ls*le;

            Vec<3> vp, vp1, vp2, vp3, nedelec;
	    Vec<3> dlamis, dlamie, dlamio;
	    for (n = 0; n < 3; n++)
	      {
		dlamis(n) = dlami(is,n);
		dlamie(n) = dlami(ie,n);
		dlamio(n) = dlami(iop,n);
		nedelec(n) = ls.Value()*dlamie(n) - le.Value()*dlamis(n);
	      }
	    vp1 = Cross (dlamis, dlamie);
	    vp2 = Cross (dlamio, dlamis);
	    vp3 = Cross (dlamie, dlamio);



	    // Divergenz RT_0 low order divshapes
	    for (j = 0; j < 3; j++)
	      divshape(i) += dlamis(j)*vp3(j) + dlamie(j)*vp2(j) + dlamio(j)*vp1(j);

	    // Divergenz Typ 1 = 0

	    // Divergenz Typ 2 = 0


	  }
      }
    ii += 4*(p*p+3*p)/2;


    //**********************************DIV INNER-SHAPES*********************************************
	if(order_inner[0] < 2) return;

    
    //***************DIV curl shapes*****************************************************
        p = order_inner[0];
    Vec<3> grad1, grad2, grad3, vp1, vp2, vp3;
    ArrayMem<AutoDiff<3>, 100> ad_poli(100);
    ArrayMem<AutoDiff<3>, 100> ad_polj(100);
    ArrayMem<AutoDiff<3>, 100> ad_polk(100);
    AutoDiff<3> pol_prod;
    AutoDiff<3> x1 (ip(0),0);
    AutoDiff<3> y1 (ip(1),1);
    AutoDiff<3> z1 (ip(2),2);

    T_INNERSHAPES::CalcSplitted(p+2, x1-(1-x1-y1-z1), y1, z1,ad_poli, ad_polj, ad_polk );
    //   ScaledLegendrePolynomial(p-1, le-ls, ls+le, ad_rec_pol1);
    //   ScaledLegendrePolynomial(p-1, lo-le-ls, lo+ls+le, ad_rec_pol2);

    for (i = 0; i <= p-2; i++)
      {

	for (j = 0; j <= p-2-i; j++)
	  {

	    for (k = 0; k <= p-2-i-j; k++, ii+=3)
	      {

		for (l = 0; l < 3; l++)
		  {
		    grad1(l) = ad_poli[i].DValue(l);
		    grad2(l) = ad_polj[j].DValue(l);
		    grad3(l) = ad_polk[k].DValue(l);
		  }
		vp1 = Cross(grad2,grad1);
		//  Divergenz Typ 1 = 0
		//  Divergenz Typ 2 = 0
		//  Divergenz Typ 3 div((adpolj.D X adpoli.D)adpolk) = - adpolj.D*(adpolk.D X adpoli.D)

            
            
		divshape(ii+2) = InnerProduct(grad3,vp1);
           
        
	      }
	  }
      }

  

    AutoDiff<3> le = 1-x1-y1-z1;

    Vec<3> ned=0.;;
    ned(0) = - x1.Value() - le.Value(); ned(1) =  - x1.Value(); ned(2)= - x1.Value();
    Vec<3> curlned=0.;
    curlned(1) = 2; // (-x1.DValue(0)*le.DValue(2));
    curlned(2) = -2; // (x1.DValue(0)*le.DValue(1));
    // Divergenz Typ 4  div(curl(Nedelec*ad_polk*ad_polj)) = 0
    // Divergenz Typ 5 div((Nedelec X ad_polj)ad_polk) = ad_polj.D*ad_polk*curlned + ad_polk.D*(ned X ad_polj.D)
    for (i = 0; i <= p-2; i++)
      {
	for (j = 0; j <= p-2-i; j++, ii+=2)
	  {

	    for (k = 0; k < 3; k++)
	      {
		grad1(k) = ad_polj[i].DValue(k);
		grad2(k) = ad_polk[j].DValue(k);

	      }
	    vp1 = Cross(ned,grad1);
	    double h1 = InnerProduct(grad2,vp1);
	    double h2 = InnerProduct(grad1,curlned);

	    divshape(ii+1) = h2*ad_polk[j].Value() + h1;
         
          
	  }
      }
    // Divergenz Typ 6 div((Nedelec X grad(lam4))*ad_polj) = grad(lam4)*ad_polj*curlned + ad_polj.D*(ned X grad(lam4))
    for (i = 0; i <= p-2; i++, ii++)
      {
	for (j = 0; j < 3; j++)
	  {
	    grad1(j) = z1.DValue(j);
	    grad2(j) = ad_polj[i].DValue(j);
	  }
	vp1 = Cross(ned,grad1);

	double h1 = InnerProduct(grad2,vp1);
	double h2 = InnerProduct(grad1,curlned);

	divshape(ii) = ad_polj[i].Value()*h2 + h1;
     
     
      }

    //(*testout)<<"divshape analytisch tet="<<divshape<<endl;
	    
    //Vector<> divshape_num(divshape.Size());
    //HDivFiniteElement<3> :: CalcDivShape (ip, divshape_num);
	    
    //(*testout) << "divshape numerisch = " << endl << divshape_num << endl;

    // (*testout) << "difference = " << endl << divshape-divshape_num << endl;

  }

#else 
  

  void HDivHighOrderFE<ET_TET> :: CalcShape (const IntegrationPoint & ip,
                                             FlatMatrixFixWidth<3> shape) const
  {
    /*
    shape = 0.0;
    T_HDivHighOrderFiniteElement<ET_TET>::CalcShape (ip, shape);
    cout << "shape xx, new = " << shape << endl;
    */

    AutoDiff<3> x(ip(0),0);
    AutoDiff<3> y(ip(1),1);
    AutoDiff<3> z(ip(2),2);

    AutoDiff<3> lami[4] = {x, y, z, 1-x-y-z}; 
    
    int i, j, ii, k,l, m, n;
    int print=0;
    shape = 0.0; 
  
        

        
        
    const FACE * faces = ElementTopology::GetFaces (ET_TET);
  
    ArrayMem<AutoDiff<3>,10> adpol1(order+2),adpol2(order+2),adpol3(order+2); 
     
    ii =4; 
    for (i = 0; i < 4; i++)
      {
        int p = order_face[i][0]; 

        int fav[3] =  { faces[i][0], faces[i][1], faces[i][2] };
      
        //Sort vertices  first edge op minimal vertex
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        int fop = 6 - fav[0] - fav[1] - fav[2];
      
        AutoDiff<3> cr0 = GradCrossGrad(lami[fav[1]],lami[fav[2]]); 
        AutoDiff<3> cr1 = GradCrossGrad(lami[fav[2]],lami[fav[0]]);
        AutoDiff<3> cr2 = GradCrossGrad(lami[fav[0]],lami[fav[1]]);
      
        AutoDiff<3> rt = cr0 * lami[fav[0]].Value() + cr1 *lami[fav[1]].Value() + cr2 *lami[fav[2]].Value();  
      
        //RT_0 low order
        for(l=0;l<3;l++) 
          shape(i,l) = rt.DValue(l); 
      
        // Curl-Shapes   
        AutoDiff<3> xi = lami[fav[1]]-lami[fav[0]];
        AutoDiff<3> eta = lami[fav[2]];
        AutoDiff<3> zeta = lami[fop];  // in [0,1] // lam_F
      
        T_FACESHAPES::CalcSplitted (p+2, xi, eta, zeta, adpol1, adpol2); 

     
        // Compability WITH TRIG!! 
        for (k=0;k<adpol1.Size();k++)
          adpol1[k]*=0.25; 
      
        // Curl (Type 2) 2*grad v x grad u
        for (j = 0; j <= p-1; j++) 
          for (k = 0; k <= p-1-j; k++, ii++)
            { 
              AutoDiff<3> hv = GradCrossGrad(adpol2[k],adpol1[j]); 
              for(l=0;l<3;l++)
                shape(ii,l) = 2*hv.DValue(l); 
            }
      
        // Curl (Type 3) //curl( * v) = nabla v x ned + curl(ned)*v
        Vec<3> ned; 
        for(l=0;l<3;l++) 
          ned(l) = lami[fav[0]].Value()*lami[fav[1]].DValue(l) - lami[fav[1]].Value()*lami[fav[0]].DValue(l);
        AutoDiff<3> curlned =  
          2*GradCrossGrad(lami[fav[0]],lami[fav[1]]);
      
        for (j = 0; j <= p-1; j++, ii++)
          {
            Vec<3> hv1 = GradCrossVec(adpol2[j],ned); 
          
            for(l=0;l<3;l++)
              shape(ii,l) = curlned.DValue(l) * adpol2[j].Value()+ hv1(l);  
          }
      }
    
      
    // cell-based shapes 
    int p = order_inner[0];
    int pc = order_inner[0]; // should be order_inner_curl  
    int pp = max(p,pc); 
    if ( p >= 2 || pc >= 2 ) 
      {
        T_INNERSHAPES::CalcSplitted(pp+2, lami[0]-lami[3], lami[1], lami[2],adpol1, adpol2, adpol3 );
      
        // Curl-Fields 
        for (i = 0; i <= pc-2; i++)
          for (j = 0; j <= pc-2-i; j++)
            for (k = 0; k <= pc-2-i-j; k++, ii+=2)
              {
                // 2 grad v  x  grad (uw)
                AutoDiff<3> hv = 2*GradCrossGrad(adpol2[j],adpol1[i]*adpol3[k]);
      
                for (l = 0; l < 3; l++)
                  shape(ii,l) = hv.DValue(l);
                        
                // 2 grad w  x  grad (uv)
                hv = 2*GradCrossGrad(adpol3[k], adpol1[i] * adpol2[j]);
      
                for (l = 0; l < 3; l++)
                  shape(ii+1,l) = hv.DValue(l); 
     
          
              }     

        // Type 1 : Curl(T3)
        // ned = lami[0] * nabla(lami[3]) - lami[3] * nabla(lami[0]) 
        double curlned[3] = {0,2,-2}; 
        Vec<3> curlnedv; curlnedv(0) = 0; curlnedv(1)=2; curlnedv(2)=-2; // = {0,2,-2}; 
        Vec<3> nedv; 
     
        for(l=0;l<3;l++) 
          nedv(l) = lami[0].Value()*lami[3].DValue(l) - lami[3].Value()*lami[0].DValue(l);
      
        for (j= 0; j <= pc-2; j++)
          for (k = 0; k <= pc-2-j; k++, ii++)
            {
              // Curl(Ned*vj*wk) = vj*wk*Curl(ned) + nabla(vj*wk)xNed0
              AutoDiff<3> pjk = adpol2[j] * adpol3[k]; 
              Vec<3> vv = GradCrossVec(pjk,nedv); 
              vv +=pjk.Value()*curlnedv; 
             
              for(l=0; l<3; l++) 
                shape(ii,l) = vv(l); 
            }
       
        // Type 2:  
        // (grad u  x  grad v) w 
        for (i = 0; i <= p-2; i++)
          for (j = 0; j <= p-2-i; j++)
            {
              AutoDiff<3> hv = GradCrossGrad(adpol1[i],adpol2[j]);
              for (k = 0; k <= p-2-i-j; k++,ii++)
                for (l = 0; l < 3; l++)
                  shape(ii,l) = hv.DValue(l) * adpol3[k].Value(); 
            }    
          
        // (ned0 x grad v) w    
        for (j = 0; j <= p-2; j++)
          {
            Vec<3> vv = GradCrossVec(adpol2[j],nedv);
            for (k= 0; k <= p-2-j; k++, ii++)
              for(l=0; l<3; l++) 
                shape(ii,l) = -adpol3[k].Value()*vv(l);
          } 
      
          
        // Type 3: 
        // (ned0 x e_z) v = (N_y, -N_x,0)^T * v ) 
        for (j=0; j<=p-2; j++,ii++) 
          {        
            Vec<3> hv = GradCrossVec(z,nedv);
            for(int l=0;l<3;l++) shape(ii,l) = -hv(l)*adpol2[j].Value();
          }
      
      }
      
    // cout << " shape, old " << shape << endl; 
        
    return; 

  }
      
  
  void HDivHighOrderFE<ET_TET> :: CalcDivShape (const IntegrationPoint & ip,
                                                FlatVector<> shape) const             
  {
    /*
      LocalHeap lh(100000); 
      FlatVector<> shape2(shape.Size(),lh); 
      HDivHighOrderFiniteElement<3> :: CalcDivShape(ip,shape2); 
    */
    AutoDiff<3> x(ip(0),0);
    AutoDiff<3> y(ip(1),1);
    AutoDiff<3> z(ip(2),2);
    AutoDiff<3> lami[4] = {x, y, z, 1-x-y-z}; 
    
            
            
    shape = 0.0; 
  
    const FACE * faces = ElementTopology::GetFaces (ET_TET);
    const POINT3D * vertices = ElementTopology::GetVertices (ET_TET);
    
    ArrayMem< AutoDiff<3> , 10> adpol1(order+2), adpol2(order+2), adpol3(order+2); 
     
    int ii = 4;
    for (int i = 0; i < 4; i++)
      {
        int p = order_face[i][0]; 
        ii += (p*p+3*p)/2; 

        int fav[3] =  { faces[i][0], faces[i][1], faces[i][2] };
      
        //Sort vertices  first edge op minimal vertex
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        int fop = 6 - fav[0] - fav[1] - fav[2];
      
   
        AutoDiff<3> cr0 = GradCrossGrad(lami[fav[1]],lami[fav[2]]); 
        AutoDiff<3> cr1 = GradCrossGrad(lami[fav[2]],lami[fav[0]]);
        AutoDiff<3> cr2 = GradCrossGrad(lami[fav[0]],lami[fav[1]]);
      
    
   
        for(int l=0;l<3;l++)
          {
            shape(i) += cr0.DValue(l)*lami[fav[0]].DValue(l) + cr1.DValue(l)*lami[fav[1]].DValue(l) +cr2.DValue(l)*lami[fav[2]].DValue(l) ; 
          }
      
    
      }
    // cell-based shapes 
    int p = order_inner[0];
    int pc = order_inner[0]; // should be order_inner_curl  
    int pp = max(p,pc); 
    if (pc >= 2 ) 
      ii += pc*(pc+1)*(pc-1)/3 + pc*(pc-1)/2; 
      
    if (p >=2) 
      { 
        
        
        T_INNERSHAPES::CalcSplitted(pp+2, x-(1-x-y-z), y, z,adpol1, adpol2, adpol3 );
      
        // Type 2:  
        // (grad u  x  grad (v)) . grad w 
        for (int i = 0; i <= p-2; i++)
          for (int j = 0; j <= p-2-i; j++)
            for (int k = 0; k <= p-2-i-j; k++,ii++)             
              { 
       
          
                AutoDiff<3> hv = GradCrossGrad(adpol1[i],adpol2[j]);
                for (int l = 0; l < 3; l++)
                  shape(ii) += hv.DValue(l) * adpol3[k].DValue(l); 
  
              }
      
            
        Vec<3> curlnedv; curlnedv(0) = 0; curlnedv(1)=2; curlnedv(2)=-2; // = {0,2,-2}; 
        Vec<3> nedv; 
     
        for(int l=0;l<3;l++) 
          nedv(l) = lami[0].Value()*lami[3].DValue(l) - lami[3].Value()*lami[0].DValue(l);
        AutoDiff<3> ned[3]; 
        for(int l=0;l<3;l++) 
          ned[l] = lami[0]*lami[3].DValue(l) - lami[3]*lami[0].DValue(l);
      
      
        // div((ned0 x grad v) w) = grad v * curlned *w + gradw *(ned0xgradv)    
        for (int j = 0; j <= p-2; j++)
          {
        
        
            Vec<3> vv = GradCrossVec(adpol2[j],nedv);
            double sp = 0; 
            for(int l=0;l<3;l++) sp+=curlnedv(l)*adpol2[j].DValue(l);
            for (int k= 0; k <= p-2-j; k++, ii++)
              { 
         
          
                shape(ii) = adpol3[k].Value()*sp; 
                for(int l=0; l<3; l++) 
                  shape(ii) -=  adpol3[k].DValue(l)*vv(l);
         
              }
          } 
      
        // div((ned0 x grad z) v) = grad z * curlned *v + gradv *(ned0xgradz)    
     
        Vec<3> vv = GradCrossVec(z,nedv);
        double sp = 0; 
        for(int l=0;l<3;l++) sp+=curlnedv(l)*z.DValue(l);
        for (int k= 0; k <= p-2; k++, ii++)
          { 
            shape(ii) = adpol2[k].Value()*sp; 
            for(int l=0; l<3; l++) 
              shape(ii) -=  adpol2[k].DValue(l)*vv(l);
          }
      
      
      }
      
    //*testout << " divshape new " << shape << endl; 
    //*testout << " shape2 " << shape2 << endl;
    //*testout << " diff " << shape-shape2 << endl;   
  }
          
#endif

  //------------------------------------------------------------------------
  // HDivHighOrderPrism
  //------------------------------------------------------------------------

  template <class T_ORTHOPOL>
  HDivHighOrderPrism<T_ORTHOPOL> :: HDivHighOrderPrism (int aorder)
    : HDivHighOrderFiniteElement<3>(ET_PRISM)
  {
    order_inner = INT<3> (aorder,aorder,aorder);

    for (int i = 0; i < 5; i++)
      order_face[i] = INT<2> (aorder,aorder);

    ComputeNDof();
  }

  template <class T_ORTHOPOL>
  void HDivHighOrderPrism<T_ORTHOPOL> :: ComputeNDof()
  {
    int i;
    ndof = 5; // RT_0


    // trig_faces
    for (i=0; i < 2; i++)
      ndof += (order_face[i][0]*order_face[i][0]+3*order_face[i][0])/2;

    // quad_faces
    for (i=2; i < 5; i++)
      ndof +=  order_face[i][0]*order_face[i][1] + order_face[i][0] + order_face[i][1];

    // SZ: ATTENTION PRISM up to now only using for order_inner[0] !!  
    if (order_inner[0]>0 )
      {
	// inner_dof horizontal
	ndof += (order_inner[0]+1)*(3*(order_inner[0]-1)+(order_inner[0]-2)*(order_inner[0]-1));

	// inner dof vertical
	ndof += (order_inner[0]-1)*(order_inner[0]+1)*(order_inner[0]+2)/2;
      }

    //(*testout)<<"Prismen ndof = "<<ndof<<endl;
    order = 0; // max(order_edges_tang,order_inner);

    for(i=0; i<2; i++)
      {
        int pp = order_face[i][0];  
        if (pp > order)
          order = pp;
      }
      
    for(i=3; i<5; i++)
      {
        int pp = max(order_face[i][0],order_face[i][1]);
        if (pp > order)
          order = pp;
      }

    if (order_inner[0] > order)
      order = order_inner[0];
   
    order++;
  }

  
  
  
  //SZ : Attention PRISMA has still isotropic inner_order
  template <class T_ORTHOPOL>
  void HDivHighOrderPrism<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip,
						    FlatMatrixFixWidth<3> shape) const
  {
    AutoDiff<2> x (ip(0), 0);
    AutoDiff<2> y (ip(1), 1);
    AutoDiff<1> z (ip(2), 0);
    AutoDiff<2> lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    AutoDiff<1> muz[6]  = { 1-z, 1-z, 1-z, z, z, z };

    Mat<6,2> clami(0);
    clami(0,1) = -1.; clami(3,1) = -1.;
    clami(1,0) = 1.; clami(4,0) = 1.;
    clami(2,1) = 1.; clami(5,1) = 1.;
    clami(2,0) = -1.; clami(5,0) = -1.;

    int i, j, k, l, m, ii;

    int p;

    shape = 0.0;


    const FACE * faces = ElementTopology::GetFaces (ET_PRISM);
    const EDGE * edges = ElementTopology::GetEdges (ET_PRISM);


    ArrayMem<AutoDiff<2>,20> adpolxy1(order+1), adpolxy2(order+1), adpolx(order+1), adpoly(order+1);
    ArrayMem<AutoDiff<1>,20> adpolz(order+1);

    ii = 5;

    // trig face shapes
    // RT_0
    // shape(0,2) = -(1-z.Value());
    // shape(1,2) = z.Value();
    //(*testout)<<"shape trig RT_0="<<shape<<endl<<endl;
    for (i = 0; i < 2; i++)
      {
	p = order_face[i][0];
        int fac = 1;
        int fav[3] = {faces[i][0], faces[i][1], faces[i][2]};

	if(vnums[fav[0]] > vnums[fav[1]]) { swap(fav[0],fav[1]); fac *= -1;}
	if(vnums[fav[1]] > vnums[fav[2]]) { swap(fav[1],fav[2]); fac *= -1;}
	if(vnums[fav[0]] > vnums[fav[1]]) { swap(fav[0],fav[1]); fac *= -1;}

	if ( i == 0)
	  shape(i,2) = -1.*fac*muz[fav[0]].Value();
	else
	  shape(i,2) = fac*muz[fav[0]].Value();


        int is = fav[0]; int ie = fav[1]; int iop = fav[2];
        AutoDiff<2> ls = lami[is];
        AutoDiff<2>  le = lami[ie];
        AutoDiff<2>  lo = lami[iop];
	AutoDiff<2> xi = lami[ie]-lami[is];
	AutoDiff<2> eta = lami[iop]; // 1-lami[f2]-lami[f1];



	//T_TRIGFACESHAPES::CalcSplitted(p+2,xi,eta,adpolxy1,adpolxy2);

	ScaledLegendrePolynomial(p-1, xi, 1-eta, adpolxy1);
	LegendrePolynomial(p-1, 2*eta-1, adpolxy2);

	for (k = 0; k <= p-1; k++)
	  {
	    adpolxy1[k] *= ls*le;
	    adpolxy2[k] *= lo;
	  }



        Vec<3> gradmu, grad1, grad2, grad3, vp1, vp2, vp3;
	gradmu = 0; grad1 = 0; grad2=0; grad3 = 0; vp1 = 0; vp2 = 0; vp3 = 0;



	// Typ 1
	for (j = 0; j <= p-1; j++)
	  {
	    for (k = 0; k <= p-1-j; k++, ii++)
	      {

		gradmu(2) = muz[iop].DValue(0);

		for (m = 0; m < 2; m++)
		  {
		    grad1(m) = adpolxy1[j].DValue(m);
		    grad2(m) = adpolxy2[k].DValue(m) * muz[iop].Value();
		  }
		grad2(2) = adpolxy2[k].Value() * muz[iop].DValue(0);


		vp3 = Cross(grad2,grad1);

		for (l = 0; l < 3; l++)
		  shape(ii, l) = 2*vp3(l);

	      }
	  }

	Vec<3> ned=0.;
	ned(0) = ls.Value()*le.DValue(0) - le.Value()*ls.DValue(0);
	ned(1) = ls.Value()*le.DValue(1) - le.Value()*ls.DValue(1);

        Vec<3> curlned=0.;
	curlned(2) = 2.0*(lami[is].DValue(0)*lami[ie].DValue(1) -
			  lami[ie].DValue(0)*lami[is].DValue(1));

	Vec<3> grad_pol2_mu = 0;
        Vec<3> vp = 0;

	// Typ 2
	//  curl(Ned0*adpolxy2[j]*muz) = grad(adpolxy2[j]*muz) X Ned0 + curlNed0(3)*adpolxy2[j]*muz
	for (j = 0; j <= p-1; j++,ii++)
	  {
	    for (k = 0; k < 2; k++)
	      grad_pol2_mu(k) = adpolxy2[j].DValue(k)*muz[iop].Value();
	    grad_pol2_mu(2) = adpolxy2[j].Value()*muz[iop].DValue(0);
	    vp = Cross(grad_pol2_mu,ned);
	    for (l = 0; l < 3; l++)
	      shape(ii,l) = vp(l) + curlned(l)* adpolxy2[j].Value()*muz[iop].Value();
	    
	  }

      }



    // quad faces
    for (i = 2; i < 5; i++)
      {
	INT<2> p = order_face[i];
	int fmax = 0;

	for (j = 1; j < 4; j++)
	  {
	    if (vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;
	  }
	int fz = 3-fmax;
	int ftrig = fmax^1;

	int fac = 1;
	int f = faces[i][fmax];
	int f1 = faces[i][ftrig];
	int f2 = faces[i][fz];

	AutoDiff<2> xi = lami[f]-lami[f1];
	AutoDiff<2> eta = 1-lami[f]-lami[f1];
	AutoDiff<1> zeta = muz[f]-muz[f2];


        int pp = max(p[0],p[1]); 
	T_ORTHOPOL::CalcTrigExt(pp+1,xi,eta,adpolxy1);
	T_ORTHOPOL::Calc(pp+1,zeta,adpolz);



	int hf1 = faces[i][(fmax+3)%4];
	int hf2 = faces[i][(fmax+1)%4];


        if (vnums[hf1] > vnums[hf2])
	  {
	    fac *= -1;
	  }

	//RT low order shape function
	for(j=0; j<2; j++)
	  {
	    shape(i,j) = fac*(lami[faces[i][1]].Value()*clami(faces[i][0],j) -
			      lami[faces[i][0]].Value()*clami(faces[i][1],j));

	  }



	// curl(nabla(polxy)*polz - polxy*nabla(polz)) = 2 * (grad(polz) X grad(polxy))
	if (vnums[f1] > vnums[f2])
	  {
	    for (k = 0; k <= p[0]-1; k++)
	      {
		for (j = 0; j <= p[1]-1; j++, ii++)
		  {

		    shape(ii,0) = -2.*adpolxy1[k].DValue(1)*adpolz[j].DValue(0);
		    shape(ii,1) =  2.*adpolxy1[k].DValue(0)*adpolz[j].DValue(0);

		  }
	      }
	  }
	else
	  {
	    for (j = 0; j <= p[0]-1; j++)
	      {
		for (k = 0; k <= p[1]-1; k++, ii++)
		  {

		    shape(ii,0) =  2.*adpolxy1[k].DValue(1)*adpolz[j].DValue(0);
		    shape(ii,1) = -2.*adpolxy1[k].DValue(0)*adpolz[j].DValue(0);

		  }
	      }
	  }


	// curl ((ned0trig)*adpolz) = 2*(grad(adpolz) X ned0trig) + 2*curlned0trig*adpolz,
	// curl((ned0_quad)* adpolxy1) = grad(adpolxy1) X ned0quad     (ned0quad = grad(zeta))
	Vec<2> ned0trig = 0.;
	for (j = 0; j < 2; j++)
	  ned0trig(j) =  lami[f1].Value()*lami[f].DValue(j)  -
	    lami[f1].DValue(j)*lami[f].Value();
        double curlned0trig =
	  -2*(lami[f].DValue(0)*lami[f1].DValue(1)  -
	      lami[f].DValue(1)*lami[f1].DValue(0));

	
        if (vnums[f1] > vnums[f2])
          {
            for (j= 0; j <= p[0]-1; j++, ii++)
              {
                shape(ii,0) =  adpolxy1[j].DValue(1)*zeta.DValue(0);
                shape(ii,1) = -adpolxy1[j].DValue(0)*zeta.DValue(0);
              }
            for(j=0; j<= p[1]-1;j++,ii++)
              {
                shape(ii,0) = -2*adpolz[j].DValue(0)*ned0trig(1);
                shape(ii,1) = 2*adpolz[j].DValue(0)*ned0trig(0);
                shape(ii,2) = 2*curlned0trig * adpolz[j].Value();
              }
          }  
        else
          {
          
            for(j=0; j<= p[0]-1;j++,ii++)
              {
                shape(ii,0) = -2*adpolz[j].DValue(0)*ned0trig(1);
                shape(ii,1) = 2*adpolz[j].DValue(0)*ned0trig(0);
                shape(ii,2) = 2*curlned0trig * adpolz[j].Value();
              }
             
            for (j= 0; j <= p[1]-1; j++, ii++)
              {
                shape(ii,0) =  adpolxy1[j].DValue(1)*zeta.DValue(0);
                shape(ii,1) = -adpolxy1[j].DValue(0)*zeta.DValue(0);
              }
          }   

      }

    p = order_inner[0];

    Mat<3,2> cdlami(0);
    cdlami(0,1) = -1.;
    cdlami(1,0) = 1.;
    cdlami(2,1) = 1.;
    cdlami(2,0) = -1.;

    if(p>=2)
      {
	//ArrayMem<double,10> rec_pol(order-1), drec_polx(order-1),  drec_polt(order-1);
	//ArrayMem<double,10> rec_pol2(order-1), drec_pol2x(order-1), drec_pol2t(order-1);

	// HACK, performed by MW, authorized by SZ
	ArrayMem<double,10> rec_pol(max2(order-1,p+1)), drec_polx(max2(order-1,p+1)),  drec_polt(max2(order-1,p+1));
	ArrayMem<double,10> rec_pol2(max2(order-1,p+1)), drec_pol2x(max2(order-1,p+1)), drec_pol2t(max2(order-1,p+1));
	ArrayMem<double,30> mem(3*order*sizeof(double));
	ArrayMem<double,30> mem2(3*order*sizeof(double));
	FlatMatrixFixWidth<3> curl_recpol(order, &mem[0]);
	FlatMatrixFixWidth<3> curl_recpol2(order, &mem2[0]);


	int fav[3] = {0, 1, 2};



	double ls = x.Value();
	double le = y.Value();
	double lo = 1-x.Value()-y.Value();




	ScaledLegendrePolynomialandDiff(p, le-ls, 1-lo,
					rec_pol, drec_polx, drec_polt);
	LegendrePolynomialandDiff(p, 2*lo-1,rec_pol2, drec_pol2x);
	LegendrePolynomial(p+1,2*z-1,adpolz);
	// \curl (ls * le * rec_pol), and \curl (lo rec_pol2)
	for (j = 0; j <= p-1; j++)
	  {
	    for (l = 0; l < 2; l++)
	      {
		curl_recpol(j,l) = ls * le * (drec_polx[j] * (cdlami(fav[1],l) - cdlami(fav[0],l)) -
					      drec_polt[j] * (cdlami(fav[2],l))) +
		  (ls * cdlami(fav[1],l) + le * cdlami(fav[0],l)) * rec_pol[j];

		curl_recpol2(j,l) = lo * (drec_pol2x[j] * 2*cdlami(fav[2],l)) +
		  cdlami(fav[2],l) * rec_pol2[j];
	      }
	  }

	// curls:
	for (i = 0; i <= p; i++)
	  for (j = 0; j <= p-2; j++)
	    for (k = 0; k <= p-2-j; k++, ii++)
	      for (l = 0; l < 2; l++)
		shape(ii, l) = (ls*le*rec_pol[j]*curl_recpol2(k,l) + curl_recpol(j,l)*lo*rec_pol2[k])*adpolz[i].Value();

	// rotations of curls
	for (i = 0; i <= p; i++)
	  for (j = 0; j <= p-2; j++)
	    for (k = 0; k <= p-2-j; k++, ii++)
	      for (l = 0; l < 2; l++)
		shape(ii, l) = (ls*le*rec_pol[j]*curl_recpol2(k,l) - curl_recpol(j,l) * lo * rec_pol2[k])*adpolz[i].Value();

	// rec_pol2 * RT_0
	for (i = 0; i <= p; i++)
	  for (j = 0; j <= p-2; j++, ii++)
	    for (l = 0; l < 2; l++)
	      shape(ii,l) = (lo*rec_pol2[j] * (ls*clami(fav[1],l)-le*clami(fav[0],l)))*adpolz[i].Value();




      }


    if (p>=2)
      {
	T_ORTHOPOL::Calc (p, 2*z-1, adpolz);

	ScaledLegendrePolynomial (p+1, x, 1-y, adpolx);
	LegendrePolynomial (p+1, 2*y-1, adpoly);
	// const. lowest-order * IntLegendre + 2 X (linear * IntLegendre)
	for( i = 0; i <= p-2; i++)
	  for ( j = 0; j <= p; j++)
	    for (k = 0; k <= p-j; k++, ii++)
	      shape(ii,2) = adpolx[j].Value()*adpoly[k].Value()*adpolz[i].Value();
      }


    //(*testout)<<" Prismen ii = "<<ii<<endl;
    return;


  }
  /*
    template <class T_ORTHOPOL>
    void HDivHighOrderPrism<T_ORTHOPOL> :: CalcDivShape (const IntegrationPoint & ip,
    FlatVector<> divshape) const
    {

    AutoDiff<2> x (ip(0), 0);
    AutoDiff<2> y (ip(1), 1);
    AutoDiff<1> z (ip(2), 0);
    AutoDiff<2> lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    AutoDiff<1> muz[6]  = { 1-z, 1-z, 1-z, z, z, z };

    Mat<6,2> clami(0);
    clami(0,1) = -1.; clami(3,1) = -1.;
    clami(1,0) = 1.; clami(4,0) = 1.;
    clami(2,1) = 1.; clami(5,1) = 1.;
    clami(2,0) = -1.; clami(5,0) = -1.;

    int i, j, k, l, m, ii, p;

    divshape = 0.0;


    const FACE * faces = ElementTopology::GetFaces (ET_PRISM);
    const EDGE * edges = ElementTopology::GetEdges (ET_PRISM);


    ArrayMem<AutoDiff<2>,20> adpolxy1(order+1), adpolxy2(order+1), adpolx(order+1), adpoly(order+1);
    ArrayMem<AutoDiff<1>,20> adpolz(order+1);

    ii = 5;

    // trig faces
    for (i = 0; i < 2; i++)
    {
    p = order_face[i];
    int fac = 1;
    int fav[3] = {faces[i][0], faces[i][1], faces[i][2]};

    if(vnums[fav[0]] > vnums[fav[1]]) { swap(fav[0],fav[1]); fac *= -1;}
    if(vnums[fav[1]] > vnums[fav[2]]) { swap(fav[1],fav[2]); fac *= -1;}
    if(vnums[fav[0]] > vnums[fav[1]]) { swap(fav[0],fav[1]); fac *= -1;}

    if ( i == 0)
    divshape(i) = -1.*fac*muz[fav[0]].DValue(0);
    else
    divshape(i) = fac*muz[fav[0]].DValue(0);

    ii += (p*p +3*p)/2;

    // div(Typ 1) = 0
    // div ( 2*(grad(adpolxy2[j]*muz) X adpolxy1[i] ) = 0
    // div (Typ 2) = 0
    //  div(curl(Ned0*adpolxy2[j]*muz)) = div(grad(adpolxy2[j]*muz) X Ned0 + curlNed0(3)*adpolxy2[j]*muz) = 0
    }



    // quad faces
    for (i = 2; i < 5; i++)
    {
    p = order_face[i];
    int fmax = 0;
    for (j = 1; j < 4; j++)
    {
    if (vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;
    }
    int hf1 = faces[i][(fmax+3)%4];
    int hf2 = faces[i][(fmax+1)%4];

    int fac = 1;
    if (vnums[hf1] > vnums[hf2])
    {
    fac *= -1;
    }
    //div (RT low order shape function)
    for(j=0; j<2; j++)
    {
    divshape(i) = fac*(lami[faces[i][1]].DValue(j)*clami(faces[i][0],j) -
    lami[faces[i][0]].DValue(j)*clami(faces[i][1],j));
    }
    // div (Typ 1) = 0
    // div (curl(nabla(polxy)*polz - polxy*nabla(polz))) = div (2 * (grad(polz) X grad(polxy)) ) = 0

    // div (Typ 2 ) and (Typ 3) = 0
    // curl ((ned0trig)*adpolz) = 2*(grad(adpolz) X ned0trig) + 2*curlned0trig*adpolz,
    // curl((ned0_quad)* adpolxy1) = grad(adpolxy1) X ned0quad     (ned0quad = grad(zeta))

    ii += p*p + 2*p;
    }

    p = order_inner;

    Mat<3,2> cdlami(0);
    cdlami(0,1) = -1.;
    cdlami(1,0) = 1.;
    cdlami(2,1) = 1.;
    cdlami(2,0) = -1.;

    if(p>=2)
    {
    ArrayMem<double,10> rec_pol(order-1), drec_polx(order-1),  drec_polt(order-1);
    ArrayMem<double,10> rec_pol2(order-1), drec_pol2x(order-1), drec_pol2t(order-1);
    ArrayMem<double,30> mem(3*order*sizeof(double));
    ArrayMem<double,30> mem2(3*order*sizeof(double));
    FlatMatrixFixWidth<3> curl_recpol(order, &mem[0]);
    FlatMatrixFixWidth<3> curl_recpol2(order, &mem2[0]);

    int fav[3] = {0, 1, 2};


    double ls = x.Value();
    double le = y.Value();
    double lo = 1-x.Value()-y.Value();




    ScaledLegendrePolynomialandDiff(p, le-ls, 1-lo,
    rec_pol, drec_polx, drec_polt);
    LegendrePolynomialandDiff(p, 2*lo-1,rec_pol2, drec_pol2x);
    LegendrePolynomial(p+1,2*z-1,adpolz);
    // \curl (ls * le * rec_pol), and \curl (lo rec_pol2)
    for (j = 0; j <= p-1; j++)
    {
    for (l = 0; l < 2; l++)
    {
    curl_recpol(j,l) = ls * le * (drec_polx[j] * (cdlami(fav[1],l) - cdlami(fav[0],l)) -
    drec_polt[j] * (cdlami(fav[2],l))) +
    (ls * cdlami(fav[1],l) + le * cdlami(fav[0],l)) * rec_pol[j];

    curl_recpol2(j,l) = lo * (drec_pol2x[j] * 2*cdlami(fav[2],l)) +
    cdlami(fav[2],l) * rec_pol2[j];
    }
    }

    // curls:
    for (i = 0; i <= p; i++)
    for (j = 0; j <= p-2; j++)
    for (k = 0; k <= p-2-j; k++, ii++)
    for (l = 0; l < 2; l++)
    shape(ii, l) = (ls*le*rec_pol[j]*curl_recpol2(k,l) + curl_recpol(j,l)*lo*rec_pol2[k])*adpolz[i].Value();

    // rotations of curls
    for (i = 0; i <= p; i++)
    for (j = 0; j <= p-2; j++)
    for (k = 0; k <= p-2-j; k++, ii++)
    for (l = 0; l < 2; l++)
    shape(ii, l) = (ls*le*rec_pol[j]*curl_recpol2(k,l) - curl_recpol(j,l) * lo * rec_pol2[k])*adpolz[i].Value();

    // rec_pol2 * RT_0
    for (i = 0; i <= p; i++)
    for (j = 0; j <= p-2; j++, ii++)
    for (l = 0; l < 2; l++)
    shape(ii,l) = (lo*rec_pol2[j] * (ls*clami(fav[1],l)-le*clami(fav[0],l)))*adpolz[i].Value();

    }


    if (p>=2)
    {
    T_ORTHOPOL::Calc (p, 2*z-1, adpolz);

    ScaledLegendrePolynomial (p+1, x, 1-y, adpolx);
    LegendrePolynomial (p+1, 2*y-1, adpoly);
    // const. lowest-order * IntLegendre + 2 X (linear * IntLegendre)
    for( i = 0; i <= p-2; i++)
    for ( j = 0; j <= p; j++)
    for (k = 0; k <= p-j; k++, ii++)
    shape(ii,2) = adpolx[j].Value()*adpoly[k].Value()*adpolz[i].Value();
    }
    //cout<<"shape.h ="<<shape.Height()<< " ii = "<<ii<<endl;


    /*  for( i = 0; i <= p-1; i++, ii++)
    {
    shape(ii,2) = 1.* adpolz[i].Value();
    }

    if (p >= 1)
    for( i = 0; i <= p-1; i++, ii+=2)
    {
    shape(ii,2) = x.Value() * adpolz[i].Value();
    shape(ii+1,2) = y.Value() * adpolz[i].Value();
    }

    // edge dofs H1 * IntLegendre
    if (p >=2)
    for (i = 0; i < 3; i++)
    {
    int es = edges[i][0];
    int ee = edges[i][1];
    if (vnums[es] > vnums[ee]) swap (es, ee);
    T_ORTHOPOL::CalcTrigExt (p-1, lami[ee]-lami[es], 1-lami[es]-lami[ee], adpolxy1);
    for (j = 0; j <= p-1; j++)
    for (k = 0; k <= p-2; k++, ii++)
    shape(ii,2) = adpolxy1[k].Value()*adpolz[j].Value();
    }
    if (p>=3)
    {
    T_TRIGFACESHAPES::Calc (p, x-y, 1-x-y, adpolxy1);
    for (j = 0; j <= p-1; j++)
    for (k = 0; k <= p-3; k++, ii++)
    shape(ii,2) = adpolxy1[k].Value()*adpolz[j].Value();
    }*/

  //(*testout)<<"shape inner="<<shape<<endl<<endl;

  /*   return;


  }*/

  template <class T_ORTHOPOL>
  void HDivHighOrderPrism<T_ORTHOPOL> ::
  GetInternalDofs (Array<int> & idofs) const
  {
    if (discontinuous)
      {
        idofs.SetSize(0);
        for (int i=0; i<ndof; i++)
          idofs.Append(i);
        return ;
      }
    else
      {
        idofs.SetSize (0);

        if(order_inner[0] >= 2) // else no inner dofs
          {
            int base = 5; // low order
      
            // trig faces
            for (int i = 0; i < 2; i++)
              {
                int p = order_face[i][0];
                base += (p*p+3*p)/2;  // see ComputeNDof
              }

            // quad faces
            for (int i=2; i<5; i++)
              {
                INT<2> p = order_face[i];
                base += p[0]*p[1]+p[0]+p[1];  // see ComputeNDof
              }

            for (int i=base; i<ndof; i++)
              idofs.Append(i);
          }
        //(*testout) << "idofs = " << idofs << endl;
      }
  }

  
  template <class T_ORTHOPOL>
  void HDivHighOrderPrism<T_ORTHOPOL> :: GetFacetDofs(int fa, Array<int> & dnums) const 
  {
    if (fa >= 5 ) 
      {
        cout << " Warning HDIVHighOrderPrism::GetFacetDofNrs() index out of range" << endl; 
      
        dnums.SetSize(0); 
        return; 
      } 
  
    dnums.SetSize(0); 
    dnums.Append(fa);  // lowest-order
    
    int base = 5; // low order
      
    int nf; 
    // trig faces
    for (int i = 0; i < 2 && i <fa; i++)
      {
        int p = order_face[i][0];
        base += (p*p+3*p)/2;  // see ComputeNDof
      }

   
    if(fa<2) 
      {
        int p = order_face[fa][0];
        nf = (p*p+3*p)/2;
      }
    else
      {
        // quad faces
        for (int i=2; i<fa; i++)
          {
            INT<2> p = order_face[i];
            base += p[0]*p[1]+p[0]+p[1];  // see ComputeNDof
          }
        INT<2> p = order_face[fa][0];
        nf = p[0]*p[1]+p[0]+p[1];
      }
      
    for (int i=0; i<nf; i++)
      dnums.Append(i+base);
  }   
  

  //------------------------------------------------------------------------
  // HDivHighOrderHex
  //------------------------------------------------------------------------

  template <class T_ORTHOPOL>
  HDivHighOrderHex<T_ORTHOPOL> :: HDivHighOrderHex (int aorder)
    : HDivHighOrderFiniteElement<3>(ET_HEX)
  {
    int i;
    // nf = 6;
    order_inner = aorder;
    for (i = 0; i < 6; i++)
      order_face[i] = aorder;
    ComputeNDof();
  }

  template <class T_ORTHOPOL>
  void HDivHighOrderHex<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 6; // RT_0



    //ndof = 0;
    int i;
   
    for (i = 0; i < 6; i++)
      {
	//if (order_face[i][0] > 0)
        {
          INT<2> p = order_face[i];
          // ndof_face += p[0]*p[1]+p[0]+p[1];
          ndof += p[0]*p[1]+p[0]+p[1];
        }
      }
    INT<3> p = order_inner;
    int ndof_inner = 3*p[0]*p[1]*p[2] + 2*p[0]*p[1] + 2 *p[1]*p[2] + 2 * p[0]*p[2] + p[0] + p[1] + p[2]; 
    //3*p*(p+1)*(p+1);
    ndof += ndof_inner;

    order = 0; // max(order_face_normal,order_inner);
    for (i = 0; i < 6; i++)
      {
        int pp = max(order_face[i][0],order_face[i][1]);
	if (pp > order)
	  order = pp;

      }
    int pp = max(order_inner[0],max(order_inner[1],order_inner[2]));
    order=max(pp,order); 

    if (order == 0) order = 1;
    order++; // integration order
  }


  template <class T_ORTHOPOL> 
  void HDivHighOrderHex<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip,
						  FlatMatrixFixWidth<3> shape) const
  {

    int i, j, k, l, m;
    AutoDiff<3> x (ip(0),0);
    AutoDiff<3> y (ip(1),1);
    AutoDiff<3> z (ip(2),2);

    AutoDiff<3> lami[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
			 (1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z};
    AutoDiff<3> sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
			  (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z};

    AutoDiff<3> ext[6] = {1-z, z, 1-y, x, y, 1-x};

    Mat<6,3> can(0);
    can(0,2) = -1.;
    can(1,2) = 1.;
    can(2,1) = -1.;
    can(3,0) = 1.;
    can(4,1) = 1.;
    can(5,0) = -1.;

    const FACE * faces = ElementTopology::GetFaces (ET_HEX);


    shape = 0.0;

    int ii = 6;

    ArrayMem<AutoDiff<3>, 20> pol_xi(order+2),pol_eta(order+2),pol_zeta(order+2);

    // (*testout)<<"x="<<x<<" y="<<y<<" z="<<z<<endl<<endl;;
    //Faces
    for (i = 0; i<6; i++)
      {
	INT<2> p = order_face[i];

	AutoDiff<3> lam_f = 0;
	for (j = 0; j < 4; j++)
	  lam_f += lami[faces[i][j]];

	int fmax = 0;
	int fac=1;
	for (j = 1; j < 4; j++)
	  {
	    if (vnums[faces[i][j]] > vnums[faces[i][fmax]])
	      {  fmax = j; //fac = -1;
	      }
	  }

	int f1 = faces[i][(fmax+3)%4];
	int f2 = faces[i][(fmax+1)%4];
	fmax = faces[i][fmax];


	if(vnums[f2] > vnums[f1])
	  {
	    swap(f1,f2);  // fmax > f1 > f2
	    fac *= -1;
	  }

	AutoDiff<3> xi = sigma[fmax]-sigma[f1];
	AutoDiff<3> eta = sigma[fmax]-sigma[f2];

	T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
	T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);

	Vec<3> grad1, grad2, grad3, grad4, vp, vp1, vp2;
	// RT_0
	for(k = 0; k < 3; k++)
	  shape(i,k) = fac*ext[i].Value()*can(i,k);
	//(*testout)<<"shape="<<shape<<endl<<endl;
	// Typ 1
	for (k = 0; k < p[0]; k++)
	  {
	    for (l = 0; l < p[1]; l++, ii++)
	      {
		for (j = 0; j < 3; j++)
		  {
		    grad1(j) = pol_xi[k].DValue(j);
		    grad2(j) = pol_eta[l].DValue(j);
		    grad3(j) = lam_f.DValue(j);
		    grad4(j) = pol_xi[k].DValue(j)*pol_eta[l].Value()-pol_xi[k].Value()*pol_eta[l].DValue(j);
		  }

		vp = Cross(grad2,grad1);
		vp1 = Cross(grad3,grad4);

		for (m = 0; m < 3; m++)
		  shape(ii,m) = 2.*lam_f.Value()*vp(m) + vp1(m);
	      }
	  }


	//Typ 2
	for (k = 0; k < p[0]; k++)
	  pol_xi[k]  *= lam_f;
        for (k = 0; k < p[1]; k++)
	  pol_eta[k] *= lam_f;
	  

	for (j = 0; j < 3; j++)
	  {

	    grad3(j) = xi.DValue(j);
	    grad4(j) = eta.DValue(j);
	  }

	for (k = 0; k < p[0]; k++, ii++)
	  {
	    for (j = 0; j < 3; j++)
	      {
		grad1(j) = pol_xi[k].DValue(j);
		
	      }
	   
	    vp2=Cross(grad1,grad4);
	    for (m = 0; m < 3; m++)
	      {
		
		shape(ii,m) = vp2(m);
	      }
	  }
        
        for (k = 0; k < p[1]; k++, ii++)
          {
            for (j = 0; j < 3; j++)
              {
            
                grad2(j) = pol_eta[k].DValue(j);
              }
            vp1=Cross(grad2,grad3);
          
            for (m = 0; m < 3; m++)
              {
                shape(ii,m)   = vp1(m);
            
              }
          }

	//(*testout)<<"shape="<<shape<<endl<<endl;
      }

    //Inner
    INT<3> p = order_inner;

    ArrayMem<AutoDiff<3>, 20> leg_x(p[0]+2), leg_y(p[1]+2), leg_z(p[2]+2);
    LegendrePolynomial (p[0]+1, 2*x-1, leg_x);
    LegendrePolynomial (p[1]+1, 2*y-1, leg_y);
    LegendrePolynomial (p[2]+1, 2*z-1, leg_z);
    T_ORTHOPOL::Calc(p[0]+1,2*x-1,pol_xi);
    T_ORTHOPOL::Calc(p[1]+1,2*y-1,pol_eta);
    T_ORTHOPOL::Calc(p[2]+1,2*z-1,pol_zeta);

    for (i = 0; i < p[0]; i++)
      for(j = 0; j <= p[1]; j++)
	for(k = 0; k <= p[2]; k++, ii++)
	  shape(ii,0)   = pol_xi[i].Value()*leg_y[j].Value()*leg_z[k].Value();
    
    for (i = 0; i <= p[0]; i++)
      for(j = 0; j < p[1]; j++)
	for(k = 0; k <= p[2]; k++, ii++)
	  shape(ii,1) = leg_x[i].Value()*pol_eta[j].Value()*leg_z[k].Value();

    for (i = 0; i <= p[0]; i++)
      for(j = 0; j <= p[1]; j++)
	for(k = 0; k < p[2]; k++, ii++)
	  shape(ii,2) = leg_x[i].Value()*leg_y[j].Value()*pol_zeta[k].Value();
    //(*testout)<<"Hex ii ="<<ii<<endl;



    //    if (order_face[4] == 3 && order_face[5] == 2)
    //      (*testout) << "shape = " << endl << shape << endl;
    return;

  }

  template <class T_ORTHOPOL>
  void HDivHighOrderHex<T_ORTHOPOL> :: CalcDivShape (const IntegrationPoint & ip,
						     FlatVector<> divshape) const
  {
    int i, j, k, l, m;
    AutoDiff<3> x (ip(0),0);
    AutoDiff<3> y (ip(1),1);
    AutoDiff<3> z (ip(2),2);

    AutoDiff<3> ext[6] = {1-z, z, 1-y, x, y, 1-x};

    int ind[6] = {2, 2, 1, 0, 1, 0};
    int can[6] = {-1, 1, -1, 1, 1, -1};

    const FACE * faces = ElementTopology::GetFaces (ET_HEX);

    divshape = 0.0;

    int ii = 6;

    ArrayMem<AutoDiff<3>, 20> pol_xi(order+2),pol_eta(order+2),pol_zeta(order+2);

    //Faces
    for (i = 0; i<6; i++)
      {
	INT<2> p = order_face[i];
	int fmax = 0;
	for (j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] > vnums[faces[i][fmax]])
	    fmax = j;
	int f1 = faces[i][(fmax+3)%4]; 
	int f2 = faces[i][(fmax+1)%4];
	fmax = faces[i][fmax];

	int fac=1;
	if(vnums[f2] > vnums[f1])
	  {
	    swap(f1,f2);  // fmax > f1 > f2
	    fac *= -1;
	  }

	// Divergenz RT_0
	divshape(i) = fac*ext[i].DValue(ind[i])*can[i];

	// Divergenz Typ 1 = 0

	// Divergenz Typ 2 = 0

	ii += (p[0]*p[1]+p[0]+p[1]);

      }


    //Inner
    INT<3> p = order_inner;
    ArrayMem<AutoDiff<3>, 20> leg_x(p[0]+2), leg_y(p[1]+2), leg_z(p[2]+2);
    AutoDiff<3> te_1=0.; AutoDiff<3> te_2=0.; AutoDiff<3> te_3=0.;
    LegendrePolynomial (p[0]+1, 2*x-1, leg_x);
    LegendrePolynomial (p[1]+1, 2*y-1, leg_y);
    LegendrePolynomial (p[2]+1, 2*z-1, leg_z);
    T_ORTHOPOL::Calc(p[0]+1,2*x-1,pol_xi);
    T_ORTHOPOL::Calc(p[1]+1,2*y-1,pol_eta);
    T_ORTHOPOL::Calc(p[2]+1,2*z-1,pol_zeta);


    for (i = 0; i < p[0]; i++)
      for(j = 0; j <= p[1]; j++)
	for(k = 0; k <= p[2]; k++, ii++)
	  {
	    te_1 = pol_xi[i]*leg_y[j]*leg_z[k];
	    divshape(ii) = te_1.DValue(0);
	  }
    for (i = 0; i <= p[0]; i++)
      for(j = 0; j < p[1]; j++)
	for(k = 0; k <= p[2]; k++, ii++)
	  {
	    te_2 = leg_x[i]*pol_eta[j]*leg_z[k];
	    divshape(ii) = te_2.DValue(1);
	  }

    for (i = 0; i <= p[0]; i++)
      for(j = 0; j <= p[1]; j++)
	for(k = 0; k < p[2]; k++, ii++)
	  {
	    te_3 = leg_x[i]*leg_y[j]*pol_zeta[k];
	    divshape(ii) = te_3.DValue(2);
	  }

    return;

  }


  template <class T_ORTHOPOL>
  void HDivHighOrderHex<T_ORTHOPOL> ::
  GetInternalDofs (Array<int> & idofs) const
  {
    if (discontinuous)
      {
        idofs.SetSize(0);
        for (int i=0; i<ndof; i++)
          idofs.Append(i);
        return ;
      }
    else 
      {
        idofs.SetSize (0);

        //if(order_inner >= 2) // else no inner dofs
        {
          int base = 6; // low order
      
          // quad faces
          for (int i=0; i<6; i++)
            {
              INT<2> p = order_face[i];
              base += p[0]*p[1]+p[0]+p[1];  // see ComputeNDof
            }

          for (int i=base; i<ndof; i++)
            idofs.Append(i);
        }
        //(*testout) << "idofs = " << idofs << endl;
      }
  }

  template <class T_ORTHOPOL>
  void HDivHighOrderHex<T_ORTHOPOL> :: GetFacetDofs(int fa, Array<int> & dnums) const 
  {
    if (fa >= 6 ) 
      {
        cout << " Warning HDIVHighOrderHex::GetFacetDofNrs() index out of range" << endl; 
      
        dnums.SetSize(0); 
        return; 
      } 
  
    dnums.SetSize(0); 
    dnums.Append(fa);  // lowest-order
    
    int base = 6; // low order
        
    // quad faces
    for (int i=0; i<fa; i++)
      {
        INT<2> p = order_face[i];
        base += p[0]*p[1]+p[0]+p[1];  // see ComputeNDof
      }
    INT<2> p = order_face[fa][0];
    int nf = p[0]*p[1]+p[0]+p[1];
         
    for (int i=0; i<nf; i++)
      dnums.Append(base+i);
  }   
  

  template class HDivHighOrderHex<IntegratedLegendreMonomialExt>;
  // template class HDivHighOrderTet<IntegratedLegendreMonomialExt>;
  template class HDivHighOrderPrism<IntegratedLegendreMonomialExt>;

  /*
    template class HDivHighOrderHex<TrigExtensionMonomial>;
    template class HDivHighOrderTet<TrigExtensionMonomial>;
    template class HDivHighOrderPrism<TrigExtensionMonomial>;
  */

  /*
    template class HDivHighOrderHex<TrigExtensionMin>;
    template class HDivHighOrderTet<TrigExtensionMin>;
    template class HDivHighOrderPrism<TrigExtensionMin>;
  */

  //  template class HDivHighOrderHex<TrigExtensionOptimal>;
  //  template class HDivHighOrderTet<TrigExtensionOptimal>;
  //  template class HDivHighOrderPrism<TrigExtensionOptimal>;


  template class  HDivHighOrderFiniteElement<2>;
  template class  HDivHighOrderFiniteElement<3>;
  template class  HDivHighOrderNormalFiniteElement<1>;
  template class  HDivHighOrderNormalFiniteElement<2>;

}




