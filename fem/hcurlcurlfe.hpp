#ifndef FILE_HCURLCURLFE
#define FILE_HCURLCURLFE

/*********************************************************************/
/* File:   hcurlcurlfe.hpp                                           */
/* Author: Michael Neunteufel                                        */
/* Date:   June 2018                                                 */
/*********************************************************************/


namespace ngfem
{

  template <int DIM>
  class HCurlCurlFiniteElement : public FiniteElement
  {
  public:
    using FiniteElement::FiniteElement;
    using FiniteElement::ndof;
    using FiniteElement::order;

    // old style
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const = 0;


    // new implementation
    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM> & mip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedShape_Vector (const MappedIntegrationPoint<DIM,DIM> & mip,
      BareSliceMatrix<double> shape) const = 0;

    //virtual void CalcMappedShape_Matrix (const SIMD_BaseMappedIntegrationRule & mir, 
    //                                     BareSliceMatrix<SIMD<double>> shapes) const = 0;
    
    //virtual void Evaluate_Matrix (const SIMD_BaseMappedIntegrationRule & ir,
    //                              BareSliceVector<> coefs,
    //                              BareSliceMatrix<SIMD<double>> values) const = 0;

    //virtual void AddTrans_Matrix (const SIMD_BaseMappedIntegrationRule & ir,
    //                              BareSliceMatrix<SIMD<double>> values,
    //                              BareSliceVector<> coefs) const = 0;
  };

  template <int D,typename VEC,typename MAT>
  void VecToSymMat(const VEC & vec,MAT & mat)
  {
    switch(D)
    {
    case 2:
      mat(0) = vec(0);
      mat(3) = vec(1);
      mat(1) = mat(2) = vec(2);
      break;
    case 3:
      mat(0) = vec(0);
      mat(4) = vec(1);
      mat(8) = vec(2);
      mat(1) = mat(3) = vec(5);
      mat(2) = mat(6) = vec(4);
      mat(5) = mat(7) = vec(3);
      break;
    }
  }

  template <typename T>
  Vec<6> SymDyadProd(AutoDiff<3,T> a, AutoDiff<3,T> b)
  {
    return Vec<6>(2*a.DValue(0)*b.DValue(0),2*a.DValue(1)*b.DValue(1),2*a.DValue(2)*b.DValue(2), a.DValue(1)*b.DValue(2)+a.DValue(2)*b.DValue(1), a.DValue(0)*b.DValue(2)+a.DValue(2)*b.DValue(0),a.DValue(1)*b.DValue(0)+a.DValue(0)*b.DValue(1));
  }

  template <ELEMENT_TYPE ET> class HCurlCurlFE;

  
  template <ELEMENT_TYPE ET>
  class T_HCurlCurlFE : public HCurlCurlFiniteElement<ET_trait<ET>::DIM>,
    public VertexOrientedFE<ET>
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };
    enum { DIM_STRESS = (DIM*(DIM+1))/2 };
    
    using VertexOrientedFE<ET>::vnums;
    using HCurlCurlFiniteElement<ET_trait<ET>::DIM>::ndof;
    using HCurlCurlFiniteElement<ET_trait<ET>::DIM>::order;
    

    INT<1> order_edge[ET_trait<ET>::N_EDGE];
    INT<DIM-1> order_facet[ET_trait<ET>::N_FACET];
    INT<DIM> order_inner;

    
  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    
    T_HCurlCurlFE (int aorder)
    {
      order = aorder;
      for (auto & of : order_facet) of = aorder;
      order_inner = aorder;

    }
    
    virtual ELEMENT_TYPE ElementType() const override { return ET; }
    const HCurlCurlFE<ET> * Cast() const { return static_cast<const HCurlCurlFE<ET>*> (this); } 
    
    INLINE void SetOrderFacet (int nr, INT<DIM-1,int> order) { order_facet[nr] = order; }
    INLINE void SetOrderEdge (int nr, INT<1,int> order) { order_edge[nr] = order; }
    INLINE void SetOrderInner (INT<DIM,int> order) { order_inner = order; }

    virtual void ComputeNDof()
    {
      cout << "Error, T_HCurlCurlFE<ET>:: ComputeNDof not available, only for ET == TRIG" << endl;
    }

    // old style
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const override
    {
      Vec<DIM, AutoDiffDiff<DIM>> adp;
      for ( int i=0; i<DIM; i++)
        adp(i) = AutoDiff<DIM>(ip(i),i);
      Cast() -> T_CalcShape (TIP<DIM, AutoDiffDiff<DIM>> (adp), SBLambda([&] (int nr, auto val)
                                          {
                                            shape.Row(nr).AddSize(DIM_STRESS) = val;
                                          }));
    }


    // new style
    virtual void CalcMappedShape_Vector (const MappedIntegrationPoint<DIM,DIM> & mip,
                            BareSliceMatrix<double> shape) const override
    {
      Vec<DIM, AutoDiff<DIM>> adp = mip;
   
      Cast() -> T_CalcShape (TIP<DIM, AutoDiffDiff<DIM>> (adp), SBLambda([shape] (int nr, auto val)
                                          {
                                            shape.Row(nr).AddSize(DIM_STRESS) = val;
                                          }));
    }


    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM> & mip,
                            BareSliceMatrix<double> shape) const override
    {
      Vec<DIM, AutoDiff<DIM>> adp = mip;
      Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM>> (adp),SBLambda([shape](int nr,auto val)
      {
        VecToSymMat<DIM> (val, shape.Row(nr));
      }));
    }
   
  };



  
#ifdef FILE_HCURLCURLFE_CPP
#define HCURLCURLFE_EXTERN
#else
#define HCURLCURLFE_EXTERN extern
#endif
  
  HCURLCURLFE_EXTERN template class HCurlCurlFiniteElement<2>;
  HCURLCURLFE_EXTERN template class HCurlCurlFiniteElement<3>;
  
  
  template <> class HCurlCurlFE<ET_TRIG> : public T_HCurlCurlFE<ET_TRIG> 
  {
    
  public:
    using T_HCurlCurlFE<ET_TRIG> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<3; i++)
      {
        ndof += order_facet[i][0]+1;
        order = max2(order, order_facet[i][0]);
      }
      int ninner = 3*order_inner[0]*(order_inner[0]+1)/2 ;
      order = max2(order, order_inner[0]);
      
      ndof += ninner;


    }
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      auto x = ip.x, y = ip.y;
      Tx ddlami[3] ={ x, y, 1-x-y };
      
      int ii = 0;
      const EDGE * edges = ElementTopology::GetEdges(ET_TRIG);

      for (int i = 0; i < 3; i++)
        {
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es,ee);

 
          Tx ls = ddlami[es], le = ddlami[ee];
          Vec<3> symdyadic = Vec<3>(2*ls.DValue(0)*le.DValue(0),2*ls.DValue(1)*le.DValue(1), ls.DValue(1)*le.DValue(0)+ls.DValue(0)*le.DValue(1));

          LegendrePolynomial::EvalScaled(order, ls-le,ls+le, SBLambda([&] (size_t nr, Tx val)
                            {
                              auto value = val.Value()*symdyadic;
                              shape[ii] = value;
                              ii++;
                            }));
        }



      if (order > 0)
        {
          for (int i = 0; i < 3; i++)
            {
              int es = edges[i][0], ee = edges[i][1];
              if (vnums[es] > vnums[ee]) swap (es,ee);
              
              Tx ls = ddlami[es], le = ddlami[ee], lt = ddlami[3-ee-es];
              Vec<3> symdyadic = lt.Value()*Vec<3>(2*ls.DValue(0)*le.DValue(0),2*ls.DValue(1)*le.DValue(1), ls.DValue(1)*le.DValue(0)+ls.DValue(0)*le.DValue(1));
              
              DubinerBasis3::Eval(order-1, ls,le,
                                  SBLambda([&] (size_t nr, Tx val)
                                           {
                                             auto value = val.Value()*symdyadic;
                                             shape[ii] = value;
                                             ii++;
                                           }));
            }
        }
      
    };
  };
  
  /*template <> class HCurlCurlFE<ET_QUAD> : public T_HCurlCurlFE<ET_QUAD> 
  {
    
  public:
    using T_HCurlCurlFE<ET_QUAD> :: T_HCurlCurlFE;

    enum {incsg = -1};
    enum {incsugv = -1};

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<4; i++)
      {
        ndof += order_facet[i][0]+1;
        order = max2(order, order_facet[i][0]);
      }
      int ninner = (order_inner[0]+1+incsg)*(order_inner[0]+1+incsg) + 
        (order_inner[0]+2)*(order_inner[0]) *2 +
        2*(order_inner[0]+1+incsugv) +1;
      order = max2(order, order_inner[0]);
      order += 5;
      ndof += ninner;

    }
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      auto x = ip.x, y = ip.y;
      Tx lx[4] ={1-x, x, x, 1-x};
      Tx ly[4] = {1-y, 1-y, y, y};
      
      int ii = 0;

      const EDGE * edges = ElementTopology::GetEdges(ET_QUAD);

      ArrayMem<Tx,20> u(order+2), v(order+2);
      
      for (int i = 0; i < 4; i++)
        {
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es,ee);
          
          Tx xi = lx[ee]+ly[ee]-lx[es]-ly[es];
          Tx eta = lx[es]*ly[es]+lx[ee]*ly[ee];

	  IntegratedLegendreMonomialExt::Calc(order_facet[i][0]+2,xi,u);

          
          for (int l = 0; l <= order_facet[i][0]; l++)
            shape[ii++] = SigmaGrad (eta*u[l]);
        }


      int oi=order_inner[0];


      IntegratedLegendreMonomialExt::Calc(oi+3,lx[0]-lx[1],u);
      IntegratedLegendreMonomialExt::Calc(oi+3,ly[0]-ly[2],v);
      
      
      for(int i = 0; i <= oi+incsg; i++)
      {
        for(int j = 0; j <= oi+incsg; j++)
        {
          shape[ii++] = SigmaGrad(u[i]*v[j]);
        }
      }
      for(int i = 0; i <= oi+1; i++)
      {
        for(int j = 0; j <= oi-1; j++)
        {
          shape[ii++] = vSigmaGradu(u[i],v[j]);
          shape[ii++] = vSigmaGradu(v[i],u[j]);
        }
      }

      shape[ii++] = Sigma_u_Gradv(lx[0], ly[0]);

      for(int i = 0; i <= oi+incsugv; i++)
      {
        shape[ii++] = Sigma_u_Gradv(u[i], ly[0]);
        shape[ii++] = Sigma_u_Gradv(v[i], lx[0]); //
      }
    };
  };

  // ***************** S_zz(uvw) ****************************** *
  // write uvw into zz component
  template <int D> class T_S_zz;
  template <> class T_S_zz<3>
  {
    AutoDiff<2> u, v;
    AutoDiff<1> w;
  public:
    T_S_zz ( AutoDiff<2> au, AutoDiff<2> av, AutoDiff<1> aw) : u(au), v(av), w(aw) { ; }
    Vec<6> Shape() 
    { 
      Vec<6> sigma(0.);
      sigma[2] = u.Value()*v.Value()*w.Value();
      return sigma;
    }

    Vec<3> DivShape()
    {
      return Vec<3> (0., 0., u.Value()*v.Value()*w.DValue(0));
    }

  };
  
  template <int D>
  auto S_zz (AutoDiff<D> au, AutoDiff<D> av, AutoDiff<1> aw)
  { return T_S_zz<D+1>(au, av, aw); }
    
  // ***************** S_xz ****************************** *
  template <int D, typename T> class T_S_xz;
  template <typename T> class T_S_xz<3,T>
  {
    AutoDiff<2,T> uv;
    AutoDiff<1,T> w;

    int comp;
  public:
    T_S_xz ( int acomp, AutoDiff<2,T> auv, AutoDiff<1,T> aw) : uv(auv), w(aw), comp(acomp) { ; }
    Vec<6,T> Shape() 
    { 
      Vec<6,T> sigma;
      sigma = 0.;
      if (comp==0)
        sigma[4] = uv.Value()*w.Value();
      else
        sigma[3] = uv.Value()*w.Value();
      return sigma;
    }

    Vec<3,T> DivShape()
    {
      if (comp == 0)
        return Vec<3,T> (uv.Value()*w.DValue(0), 0, uv.DValue(0)*w.Value() );
      else
        return Vec<3,T> (0, uv.Value()*w.DValue(0), uv.DValue(1)*w.Value() );
    }

  };
  
  template <int D, typename T>
  auto S_xz (int comp, AutoDiff<D,T> auv, AutoDiff<1,T> aw)
  { return T_S_xz<D+1,T>(comp,auv, aw); }



  template <typename T>
  class T_Prism_wSigmaGradu
  {
    AutoDiffDiff<2,T> u;
    AutoDiff<1,T> w;
  public:
    T_Prism_wSigmaGradu ( AutoDiffDiff<2,T> au, AutoDiff<1,T> aw) : u(au), w(aw) { ; }
    Vec<6,T> Shape() 
    { 
      Vec<3,T> sigma2d = T_SigmaGrad<2,T>(u).Shape();
      Vec<6,T> sigma(0.);
      sigma[0] = w.Value()*sigma2d[0];
      sigma[1] = w.Value()*sigma2d[1];
      sigma[5] = w.Value()*sigma2d[2];
      return sigma;
    }

    Vec<3,T> DivShape()
    {
      return Vec<3,T> (0., 0., 0);
    }

  };

  template <typename T>
  auto Prism_wSigmaGradu ( AutoDiffDiff<2,T> au, AutoDiff<1,T> aw)
  { return T_Prism_wSigmaGradu<T>(au, aw); }


  template <typename T>
  class T_Prism_wType2
  {
    AutoDiffDiff<2,T> u, v;
    AutoDiff<1,T> w;
  public:
    T_Prism_wType2 ( AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiff<1,T> aw) : u(au), v(av),  w(aw) { ; }
    Vec<6,T> Shape() 
    { 
      Vec<3,T> sigma2d = T_Type2<2,T>(u,v).Shape();
      Vec<6,T> sigma(0.);
      sigma[0] = w.Value()*sigma2d[0];
      sigma[1] = w.Value()*sigma2d[1];
      sigma[5] = w.Value()*sigma2d[2];
      return sigma;
    }

    Vec<3,T> DivShape()
    {
      Vec<2> divsigma2d = w.Value()*T_Type2<2,T>(u,v).DivShape();
      return Vec<3,T> (divsigma2d[0], divsigma2d[1], 0);
    }

  };

  template <typename T>
  auto Prism_wType2 (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiff<1,T> aw)
  { return T_Prism_wType2<T>(au, av, aw); }


  template <typename T>
  class T_Prism_wType3
  {
    AutoDiffDiff<2,T> u, v;
    AutoDiff<1,T> w;
  public:
    T_Prism_wType3 ( AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiff<1,T> aw) : u(au), v(av),  w(aw) { ; }
    Vec<6> Shape() 
    { 
      Vec<3,T> sigma2d = T_Type3<2,T>(u,v).Shape();
      Vec<6,T> sigma(0.);
      sigma[0] = w.Value()*sigma2d[0];
      sigma[1] = w.Value()*sigma2d[1];
      sigma[5] = w.Value()*sigma2d[2];
      return sigma;
    }

    Vec<3,T> DivShape()
    {
      Vec<2,T> divsigma2d = w.Value()*T_Type3<2,T>(u,v).DivShape();
      return Vec<3,T> (divsigma2d[0], divsigma2d[1], 0);
    }
  };

  template <typename T>
  auto Prism_wType3 (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiff<1,T> aw)
  { return T_Prism_wType3<T>(au, av, aw); }


  template <typename T>
  class T_Prism_wType4
  {
    AutoDiffDiff<2,T> u, v, w;
    AutoDiff<1,T> wz;
  public:
    T_Prism_wType4 ( AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiffDiff<2,T> aw, AutoDiff<1,T> awz) : u(au), v(av),  w(aw), wz(awz) { ; }
    Vec<6,T> Shape() 
    { 
      Vec<3,T> sigma2d = wz.Value()*T_Sigma_Duv_minus_uDv_w<2,T>(u,v,w).Shape();
      Vec<6,T> sigma(0.);
      sigma[0] = sigma2d[0];
      sigma[1] = sigma2d[1];
      sigma[5] = sigma2d[2];
      return sigma;
    }

    Vec<3,T> DivShape()
    {
      Vec<2,T> divsigma2d = wz.Value()*T_Sigma_Duv_minus_uDv_w<2,T>(u,v,w).DivShape();
      return Vec<3,T> (divsigma2d[0], divsigma2d[1], 0);
    }

  };

  template <typename T>
  auto Prism_wType4 (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, AutoDiffDiff<2,T> aw, AutoDiff<1,T> awz)
  { return T_Prism_wType4<T>(au, av, aw, awz); }

  

  class Prism_SymRotRot_Dl2xDl1_vw
  {
    AutoDiff<2> l1,l2,v;
    AutoDiff<1> wz;
  public:
    Prism_SymRotRot_Dl2xDl1_vw ( AutoDiff<2> lam1, AutoDiff<2> lam2, AutoDiff<2> av, AutoDiff<1> awz) : l1(lam1), l2(lam2), v(av), wz(awz) { ; }
    Vec<6> Shape() 
    { 
      Vec<3> sigma2d = wz.Value()*T_SymRotRot_Dl2xDl1_v(l1,l2,v).Shape();
      Vec<6> sigma(0.);
      sigma[0] = sigma2d[0];
      sigma[1] = sigma2d[1];
      sigma[5] = sigma2d[2];
      return sigma;
    }

    Vec<3> DivShape()
    {
      Vec<2> divsigma2d = wz.Value()*T_SymRotRot_Dl2xDl1_v(l1,l2,v).DivShape();
      return Vec<3> (divsigma2d[0], divsigma2d[1], 0);
    }

  };

  template <typename T>
  class Prism_Dl1xDl3_symtensor_Dl2xDl4_u
  {
    AutoDiff<3,T> l1,l2,l3, l4;
    AutoDiff<3,T> u;
  public:
    Prism_Dl1xDl3_symtensor_Dl2xDl4_u ( AutoDiff<3,T> lam1, AutoDiff<3,T> lam2, AutoDiff<3,T> alz1, AutoDiff<3,T> alz2, AutoDiff<3,T> av) 
      : l1(lam1), l2(lam2), l3(alz1), l4(alz2), u(av) { ; }
    Vec<6,T> Shape() 
    { 
      auto rotlam1 = Cross(l1, l3);
      auto rotlam2 = Cross(l2, l4);

      Vec<6,T> sigma(0.);
      sigma[0] = u.Value()*rotlam1.DValue(0)*rotlam2.DValue(0);
      sigma[1] = u.Value()*rotlam1.DValue(1)*rotlam2.DValue(1);
      sigma[2] = u.Value()*rotlam1.DValue(2)*rotlam2.DValue(2);
      sigma[3] = 0.5*u.Value()*(rotlam1.DValue(2)*rotlam2.DValue(1) + rotlam2.DValue(2)*rotlam1.DValue(1));
      sigma[4] = 0.5*u.Value()*(rotlam1.DValue(2)*rotlam2.DValue(0) + rotlam2.DValue(2)*rotlam1.DValue(0));
      sigma[5] = 0.5*u.Value()*(rotlam1.DValue(0)*rotlam2.DValue(1) + rotlam2.DValue(0)*rotlam1.DValue(1));
      return sigma;
    }

    INLINE AutoDiff<3,T> Cross (const AutoDiff<3,T> & x,
                                const AutoDiff<3,T> & y)
    {
      T hv[3];
      hv[0] = x.DValue(1)*y.DValue(2)-x.DValue(2)*y.DValue(1);
      hv[1] = x.DValue(2)*y.DValue(0)-x.DValue(0)*y.DValue(2);
      hv[2] = x.DValue(0)*y.DValue(1)-x.DValue(1)*y.DValue(0);
      return AutoDiff<3,T> (0,hv);
    }
    Vec<3,T> DivShape()
    {
      auto lam1 = Cross(l1, l3);
      auto lam2 = Cross(l2, l4);
      return Vec<3,T> (u.DValue(0)*lam1.DValue(0)*lam2.DValue(0) + 
        0.5*(u.DValue(1)*(lam1.DValue(0)*lam2.DValue(1)+lam2.DValue(0)*lam1.DValue(1)) + u.DValue(2)*(lam1.DValue(0)*lam2.DValue(2)+lam2.DValue(0)*lam1.DValue(2))),
        u.DValue(1)*lam1.DValue(1)*lam2.DValue(1) + 
        0.5*(u.DValue(0)*(lam1.DValue(0)*lam2.DValue(1)+lam2.DValue(0)*lam1.DValue(1)) + u.DValue(2)*(lam1.DValue(1)*lam2.DValue(2)+lam2.DValue(1)*lam1.DValue(2))),
        u.DValue(2)*lam1.DValue(2)*lam2.DValue(2) + 
        0.5*(u.DValue(0)*(lam1.DValue(0)*lam2.DValue(2)+lam2.DValue(0)*lam1.DValue(2)) + u.DValue(1)*(lam1.DValue(1)*lam2.DValue(2)+lam2.DValue(1)*lam1.DValue(2)))
        );
    }

  };


  template <> class HCurlCurlFE<ET_PRISM> : public T_HCurlCurlFE<ET_PRISM> 
  {
  public:
    // order k+1 for certain components, for inner and boundary shapes
    // analysis from TDNNS paper for case xx1=0, zz1=xx2=zz2=1 for inner and boundary shapes
    // however, works also when boundary order is not increased.. check
    enum { incrorder_xx1 = 0};
    enum { incrorder_zz1 = 1};
    enum { incrorder_xx2 = 1};
    enum { incrorder_zz2 = 1};
    enum { incrorder_xx1_bd = 0};
    enum { incrorder_zz1_bd = 0};
    enum { incrorder_xx2_bd = 0};
    enum { incrorder_zz2_bd = 0};
    using T_HCurlCurlFE<ET_PRISM> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<2; i++)
      {
        ndof += (order_facet[i][0]+1+incrorder_zz1_bd)*(order_facet[i][0]+2+incrorder_zz1_bd)/2;
        order = max2(order, order_facet[i][0]+incrorder_zz1_bd);
      }
      for (int i=2; i<5; i++)
      {
        ndof += (order_facet[i][0]+1+incrorder_xx1_bd)*(order_facet[i][1]+1+incrorder_xx2_bd);
        order = max2(order, order_facet[i][0]+incrorder_xx2_bd);
      }
      int oi0 = order_inner[0];
      int oi2 = order_inner[2];
      int ninner = 3*((oi0+1+incrorder_xx1)*(oi0+incrorder_xx1))/2 *(oi2+1+incrorder_xx2) 
        + (oi0+1)*(oi0+2)*(oi2+1) 
        + (oi0+1+incrorder_zz1)*(oi0+2+incrorder_zz1)*(oi2-1+incrorder_zz2)/2;
      ndof += ninner; 

      order = max3(order, oi0+incrorder_zz1, oi2+incrorder_zz2);

    }

    // works only with old-style Transformation
    // does not work with CalcMappedShape
   template <typename Tx, typename TFA> 
    void T_CalcShape_Complex (TIP<3,Tx> ip, TFA & shape) const
    {
      AutoDiffDiff<2> x(ip.x.Value(),0);
      AutoDiffDiff<2> y(ip.y.Value(),1);
      AutoDiff<2> xd(ip.x.Value(),0);
      AutoDiff<2> yd(ip.y.Value(),1);
      AutoDiff<1> z(ip.z.Value(), 0);
      AutoDiffDiff<2> lami[6] ={ x,y,1-x-y,x,y,1-x-y };
      AutoDiff<2> lamid[6] ={ xd,yd,1-xd-yd,xd,yd,1-xd-yd };
      AutoDiff<1> lamiz[6] ={ 1-z,1-z,1-z,z,z,z };

      int ii = 0;
      
      int maxorder_facet =
        max2(order_facet[0][0],max2(order_facet[1][0],order_facet[2][0]));

      const FACE * faces = ElementTopology::GetFaces(ET_PRISM);

      ArrayMem<AutoDiffDiff<2>,20> ha(maxorder_facet+2);
      ArrayMem<AutoDiffDiff<2>,20> u(order+2), v(order+3);
      ArrayMem<AutoDiff<2>,20> leg_u(order+2), leg_v(order+3);
      ArrayMem<AutoDiff<1>,20> leg_w(order+2);

      
      // Trig faces, (p+1)(p+2)/2
      for (int fa=0; fa<2; fa++)
      {
        int fav[3] ={faces[fa][0],faces[fa][1],faces[fa][2]};

        //Sort vertices  first edge op minimal vertex
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        
        leg_u.SetSize(order_facet[fa][0]+incrorder_zz1_bd+1);
        leg_v.SetSize(order_facet[fa][0]+incrorder_zz1_bd+1);
        ScaledLegendrePolynomial(order_facet[fa][0]+incrorder_zz1_bd,lamid[fav[0]]-lamid[fav[1]],1-lamid[fav[0]]-lamid[fav[1]],leg_u);
        LegendrePolynomial::Eval(order_facet[fa][0]+incrorder_zz1_bd,2 * lamid[fav[2]] - 1,leg_v);

        for(int j = 0; j <= order_facet[fa][0]+incrorder_zz1_bd; j++)
          for(int k = 0; k <= order_facet[fa][0]+incrorder_zz1_bd-j; k++)
            shape[ii++] = S_zz(leg_u[j],leg_v[k],lamiz[fav[0]]);
      }
      // quad faces -- use face bubbles of trig multiplied by leg_w
      // (px+1)(pz+1)
      for(int fa = 2; fa < 5; fa++)
      {
        int fmax = 0;
        for(int j = 1; j < 4; j++)
          if(vnums[faces[fa][j]] > vnums[faces[fa][fmax]]) fmax = j;

        int fz,ftrig;

        fz = 3 - fmax;
        ftrig = fmax^1;

        fmax = faces[fa][fmax];
        fz = faces[fa][fz];
        ftrig = faces[fa][ftrig];


        // int orderz = order_facet[fa][1];

        bool rotate = false;
        if(vnums[fz] > vnums[ftrig]) rotate = true;
        leg_w.SetSize(order_facet[fa][1]+incrorder_xx2_bd+1);
        ha.SetSize(order_facet[fa][0]+incrorder_xx1_bd+1);
        LegendrePolynomial::Eval(order_facet[fa][1]+incrorder_xx2_bd,lamiz[fmax]*2-1,leg_w);


        // edge functions are all div-free!
        IntegratedLegendreMonomialExt::CalcTrigExt(order_facet[fa][0]+incrorder_xx1_bd+2,
          lami[fmax]-lami[ftrig],1-lami[fmax]-lami[ftrig],ha);

        if(rotate)
          for(int k = 0; k <= order_facet[fa][1]+incrorder_xx2_bd; k++)
            for(int l = 0; l <= order_facet[fa][0]+incrorder_xx1_bd; l++)
            {
              shape[ii++] = Prism_wSigmaGradu(ha[l],leg_w[k]);
            }

        else
          for(int l = 0; l <= order_facet[fa][0]+incrorder_xx1_bd; l++)
            for(int k = 0; k <= order_facet[fa][1]+incrorder_xx2_bd; k++)
            {
              shape[ii++] = Prism_wSigmaGradu(ha[l],leg_w[k]);
            }


      }



      int oi = order_inner[0];
      leg_u.SetSize(oi+incrorder_zz1+1);
      leg_v.SetSize(oi+incrorder_zz1+1);
      leg_w.SetSize(oi+incrorder_xx2+1);
      u.SetSize(oi-1+incrorder_xx1+1);
      v.SetSize(oi-1+incrorder_xx1+1);

      ScaledLegendrePolynomial(oi+incrorder_zz1, lamid[0]-lamid[1], 1-lamid[0]-lamid[1], leg_u);
      LegendrePolynomial::Eval(oi+incrorder_zz1, 2*lamid[2]-1, leg_v);
      LegendrePolynomial::Eval(oi+incrorder_xx2, 2*lamiz[0]-1, leg_w);

      // ------------------------------------
      // based on elasticity-complex-based triangle shapes
      IntegratedLegendreMonomialExt::CalcTrigExt(oi-1+incrorder_xx1+2,lami[0]-lami[1],1-lami[0]-lami[1],u);
      LegendrePolynomial::EvalMult(oi-1+incrorder_xx1,2*lami[2]-1, lami[2], v);
      for(int k=0; k<=oi+incrorder_xx2; k++)
      {
        for(int i = 0; i <= oi-1+incrorder_xx1; i++)
        {
          for(int j = 0; j+i <= oi-1+incrorder_xx1; j++)
          {
            shape[ii++] = Prism_wSigmaGradu(u[i]*v[j],leg_w[k]);
            shape[ii++] = Prism_wType2(u[i],v[j],leg_w[k]);
          }
        }
        for(int i = 0; i <= oi-1+incrorder_xx1; i++)
        {
          for(int j = 0; j+i <= oi-1+incrorder_xx1; j++)
          {
            if(j > 0)
            {
              shape[ii++] = Prism_wType3(u[i],v[j],leg_w[k]);
            }
          }
        }
        for (int i = 0; i < oi+incrorder_xx1; i++)
        {
          shape[ii++] = Prism_wType4 (lami[0], -lami[1], v[i],leg_w[k]);
        }

      }

      // S_xz
      for (int i=0; i<=oi; i++)
      {
        for (int j=0; j+i<=oi; j++)
        {
          AutoDiff<2> uv = leg_u[i]*leg_v[j];
          for (int k=0; k<=oi; k++)
          {
            shape[ii++] = S_xz(0,uv, leg_w[k]);
            shape[ii++] = S_xz(1,uv, leg_w[k]);
          }
        }
      }

      // S_zz
      for(int k=0; k<=oi-2+incrorder_zz2; k++)
      {
        AutoDiff<1> bubw = leg_w[k]*lamiz[0]*(1-lamiz[0]);
        for(int i=0; i<=oi+incrorder_zz1; i++)
        {
          for(int j=0; j<=oi+incrorder_zz1-i; j++)
          {
            shape[ii++] = S_zz(leg_u[i],leg_v[j],bubw);
          }
        }
      }


      };




    // alternative to T_CalcShape, with "simpler" shape functions,
    // that are described in anisotropic paper
    // works with CalcMappedShape etc. routines, also for curved elements
    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y, z = ip.z;
      typedef decltype(x.Value()+x.Value()) T;
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> zz(z.Value(), &z.DValue(0));
      AutoDiff<3,T> lx[6] ={ xx, yy, 1-xx-yy, xx, yy, 1-xx-yy };
      AutoDiff<3,T> lz[6] ={ 1-zz,1-zz,1-zz,zz,zz,zz };
      int ii = 0;
      
      // int maxorder_facet =
      // max2(order_facet[0][0],max2(order_facet[1][0],order_facet[2][0]));

      const FACE * faces = ElementTopology::GetFaces(ET_PRISM);

      ArrayMem<AutoDiff<3,T>,20> leg_u(order+2), leg_v(order+3);
      ArrayMem<AutoDiff<3,T>,20> leg_w(order+2);

      
      // Trig faces, (p+1)(p+2)/2
      for (int fa=0; fa<2; fa++)
      {
        int fav[3] ={faces[fa][0],faces[fa][1],faces[fa][2]};

        //Sort vertices  first edge op minimal vertex
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
        if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        
        leg_u.SetSize(order_facet[fa][0]+incrorder_zz1_bd+1);
        leg_v.SetSize(order_facet[fa][0]+incrorder_zz1_bd+1);
        ScaledLegendrePolynomial(order_facet[fa][0]+incrorder_zz1_bd,lx[fav[0]]-lx[fav[1]],lx[fav[0]]+lx[fav[1]],leg_u);
        LegendrePolynomial::Eval(order_facet[fa][0]+incrorder_zz1_bd,2 * lx[fav[2]] - 1,leg_v);

        for(int j = 0; j <= order_facet[fa][0]+incrorder_zz1_bd; j++)
          for(int k = 0; k <= order_facet[fa][0]+incrorder_zz1_bd-j; k++)
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[fav[0]], lx[fav[1]], lx[fav[2]], lx[fav[2]], leg_u[j]*leg_v[k]*lz[fav[0]]);
      }
      // quad faces -- use face bubbles of trig multiplied by leg_w
      // (px+1)(pz+1)
      for(int fa = 2; fa < 5; fa++)
      {
        int fmax = 0;
        for(int j = 1; j < 4; j++)
          if(vnums[faces[fa][j]] > vnums[faces[fa][fmax]]) fmax = j;

        int fz,ftrig;

        fz = 3 - fmax;
        ftrig = fmax^1;

        fmax = faces[fa][fmax];
        fz = faces[fa][fz];
        ftrig = faces[fa][ftrig];


        // int orderz = order_facet[fa][1];

        bool rotate = false;
        if(vnums[fz] > vnums[ftrig]) rotate = true;
        leg_w.SetSize(order_facet[fa][1]+incrorder_xx2_bd+1);
        LegendrePolynomial::Eval(order_facet[fa][1]+incrorder_xx2_bd,lz[fmax]*2-1,leg_w);


        ScaledLegendrePolynomial(order_facet[fa][0]+incrorder_xx1_bd, lx[fmax]-lx[ftrig], lx[fmax]+lx[ftrig], leg_u);      

        if(rotate)
          for(int k = 0; k <= order_facet[fa][1]+incrorder_xx2_bd; k++)
            for(int l = 0; l <= order_facet[fa][0]+incrorder_xx1_bd; l++)
            {
              shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lz[fmax], lz[fz], lx[fmax], lx[ftrig], leg_u[l]* leg_w[k]);
            }

        else
          for(int l = 0; l <= order_facet[fa][0]+incrorder_xx1_bd; l++)
            for(int k = 0; k <= order_facet[fa][1]+incrorder_xx2_bd; k++)
            {
              shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[fmax], lx[ftrig], lz[fmax], lz[fz], leg_u[l]* leg_w[k]);
            }


      }



      int oi = order_inner[0];
      leg_u.SetSize(oi+incrorder_zz1+1);
      leg_v.SetSize(oi+incrorder_zz1+1);
      leg_w.SetSize(oi+incrorder_xx2+1);

      ScaledLegendrePolynomial(oi+incrorder_zz1, lx[0]-lx[1], lx[0]+lx[1], leg_u);
      LegendrePolynomial::Eval(oi+incrorder_zz1, 2*lx[2]-1, leg_v);
      LegendrePolynomial::Eval(oi+incrorder_xx2, 2*lz[0]-1, leg_w);

      for(int k=0; k<=oi+incrorder_xx2; k++)
      {
        for(int i = 0; i <= oi-1+incrorder_xx1; i++)
        {
          for(int j = 0; j+i <= oi-1+incrorder_xx1; j++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[0], lx[1], lz[0], lz[0], lx[2]*leg_u[i]*leg_v[j]* leg_w[k]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[2], lx[0], lz[0], lz[0], lx[1]*leg_u[i]*leg_v[j]* leg_w[k]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[1], lx[2], lz[0], lz[0], lx[0]*leg_u[i]*leg_v[j]* leg_w[k]);
          }
        }
      }


      // S_xz
      for (int i=0; i<=oi; i++)
      {
        for (int j=0; j+i<=oi; j++)
        {
          for (int k=0; k<=oi; k++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lz[0], lx[1], lx[0], lx[0], leg_u[i]*leg_v[j]*leg_w[k]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lz[0], lx[0], lx[1], lx[1], leg_u[i]*leg_v[j]*leg_w[k]);
          }
        }
      }

      // S_zz
      for(int k=0; k<=oi-2+incrorder_zz2; k++)
      {
        AutoDiff<3,T> bubw = leg_w[k]*lz[0]*(1-lz[0]);
        for(int i=0; i<=oi+incrorder_zz1; i++)
        {
          for(int j=0; j<=oi+incrorder_zz1-i; j++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[0], lx[2], lx[1], lx[1], leg_u[i]*leg_v[j]*bubw);
          }
        }
      }


    };

  };

  */

  
  template <> class HCurlCurlFE<ET_TET> : public T_HCurlCurlFE<ET_TET> 
  {
  public:
    using T_HCurlCurlFE<ET_TET> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;

      for (int i=0; i<6; i++)
      {
        ndof += order_edge[i][0]+1;
        order = max2(order, order_edge[i][0]);
      }
      
      for (int i=0; i<4; i++)
      {
        ndof += 3*(order_facet[i][0])*(order_facet[i][0]+1)/2;
        order = max2(order, order_facet[i][0]);
      }
     
      int p = order_inner[0];
      int ninner = p > 1 ? 6*(p+1)*(p)*(p-1)/6 : 0;
      ndof += ninner; 

      order = max2(order, p);
    }


    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y, z = ip.z;
      typedef decltype(x.Value()+x.Value()) T;      
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> zz(z.Value(), &z.DValue(0));
      AutoDiff<3,T> lam[4] = {xx, yy, zz, 1-xx-yy-zz};
      int ii = 0;

      for (int i = 0; i < 6; i++)
        {
          INT<2> e = ET_trait<ET_TET>::GetEdgeSort (i, vnums);
          auto ls = lam[e[0]], le = lam[e[1]];

          Vec<6> symdyadic = SymDyadProd(ls,le);

          LegendrePolynomial::EvalScaled(order, ls-le,ls+le, SBLambda([&] (size_t nr, auto val)
                            {
                              auto value = val.Value()*symdyadic;
                              shape[ii] = value;
                              ii++;
                            }));
        }

      if (order > 0)
        {
          for(int fa = 0; fa < 4; fa++)
            {
              INT<4> f = ET_trait<ET_TET>::GetFaceSort(fa, vnums);
              auto ls = lam[f[0]], le = lam[f[1]], lt = lam[f[2]];
              
              Vec<6> symdyadic1 = lt.Value()*SymDyadProd(ls,le);
              Vec<6> symdyadic2 = ls.Value()*SymDyadProd(lt,le);
              Vec<6> symdyadic3 = le.Value()*SymDyadProd(ls,lt);
              
              DubinerBasis3::Eval(order-1, ls,le,
                                  SBLambda([&] (size_t nr, auto val)
                                           {
                                             auto value = val.Value()*symdyadic1;
                                             shape[ii] = value;
                                             ii++;
                                             value = val.Value()*symdyadic2;
                                             shape[ii] = value;
                                             ii++;
                                             value = val.Value()*symdyadic3;
                                             shape[ii] = value;
                                             ii++;
                                           }));
            }
        }

      if (order > 1)
        {
          auto li = lam[0], lj = lam[1], lk = lam[2], ll = lam[3];

          Vec<6> symdyadic1 = li.Value()*lj.Value()*SymDyadProd(lk,ll);
          Vec<6> symdyadic2 = lj.Value()*lk.Value()*SymDyadProd(ll,li);
          Vec<6> symdyadic3 = lk.Value()*ll.Value()*SymDyadProd(li,ll);
          Vec<6> symdyadic4 = ll.Value()*li.Value()*SymDyadProd(lj,lk);
          Vec<6> symdyadic5 = li.Value()*lk.Value()*SymDyadProd(lj,ll);
          Vec<6> symdyadic6 = lj.Value()*ll.Value()*SymDyadProd(li,lk);


          LegendrePolynomial leg;
          JacobiPolynomialAlpha jac1(1);    
          leg.EvalScaled1Assign 
            (order-2, lam[2].Value()-lam[3].Value(), lam[2].Value()+lam[3].Value(),
             SBLambda ([&](size_t k, auto polz) LAMBDA_INLINE
                       {
                         // JacobiPolynomialAlpha jac(2*k+1);
                         JacobiPolynomialAlpha jac2(2*k+2);
                         
                         jac1.EvalScaledMult1Assign
                           (order-2-k, lam[1].Value()-lam[2].Value()-lam[3].Value(), 1-lam[0].Value(), polz, 
                            SBLambda ([&] (size_t j, double polsy) LAMBDA_INLINE
                                      {
                                        // JacobiPolynomialAlpha jac(2*(j+k)+2);
                                        jac2.EvalMult(order-2 - k - j, 2 * lam[0].Value() - 1, polsy, 
                                                      SBLambda([&](size_t j, auto val) LAMBDA_INLINE
                                                               {
                                                                 auto value = val*symdyadic1;
                                                                 shape[ii] = value;
                                                                 ii++;
                                                                 value = val*symdyadic2;
                                                                 shape[ii] = value;
                                                                 ii++;
                                                                 value = val*symdyadic3;
                                                                 shape[ii] = value;
                                                                 ii++;
                                                                 value = val*symdyadic4;
                                                                 shape[ii] = value;
                                                                 ii++;
                                                                 value = val*symdyadic5;
                                                                 shape[ii] = value;
                                                                 ii++;
                                                                 value = val*symdyadic6;
                                                                 shape[ii] = value;
                                                               }));
                                        jac2.IncAlpha2();
                                      }));
                         jac1.IncAlpha2();
                       }));
        }
      

    };

  };
  

  /*


  template <> class HCurlCurlFE<ET_HEX> : public T_HCurlCurlFE<ET_HEX> 
  {
  public:
    using T_HCurlCurlFE<ET_HEX> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<6; i++)
      {
        ndof += (order_facet[i][0]+1)*(order_facet[i][0]+1);
        order = max2(order, order_facet[i][0]);
      }
      int p = order_inner[0];
      int ninner = 3*p*(p+2)*(p+2) + 3*(p+2)*(p+1)*(p+1);
      ndof += ninner; 

      order = max2(order, p);

    }


    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y, z = ip.z;
      typedef decltype(x.Value()+x.Value()) T;            
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> zz(z.Value(), &z.DValue(0));
      AutoDiff<3,T> lx[2] ={ 1-xx, xx};
      AutoDiff<3,T> ly[2] ={ 1-yy, yy};
      AutoDiff<3,T> lz[2] ={ 1-zz, zz};
      AutoDiff<3,T> sigma[8] = {1-xx + 1-yy + 1-zz,
        xx + 1-yy + 1-zz,
        xx + yy + 1-zz,
        1-xx + yy + 1-zz,
        1-xx + 1-yy + zz,
        xx + 1-yy + zz,
        xx + yy + zz,
        1-xx + yy + zz};
      int ii = 0;
      
      // int maxorder_facet =
      //     max2(order_facet[0][0],max2(order_facet[1][0],order_facet[2][0]));

      const FACE * faces = ElementTopology::GetFaces(ET_HEX);

      ArrayMem<AutoDiff<3,T>,20> leg_u(order+2), leg_v(order+3);
      ArrayMem<AutoDiff<3,T>,20> leg_w(order+2);
      AutoDiff<3,T> lam_face;
      
      for(int fa = 0; fa < 6; fa++)
      {
        int fmax = 0;
        lam_face = -1 + 0.25*sigma[faces[fa][0]];
        for(int j = 1; j < 4; j++)
        {
          if(vnums[faces[fa][j]] > vnums[faces[fa][fmax]]) fmax = j;
          lam_face += sigma[faces[fa][j]]*0.25;
        }

        int fz,ftrig;

        fz = 3 - fmax;
        ftrig = fmax^1;

        fmax = faces[fa][fmax];
        fz = faces[fa][fz];
        ftrig = faces[fa][ftrig];


        // int orderz = order_facet[fa][1];

        if(vnums[fz] < vnums[ftrig]) swap(fz, ftrig);
        int p = order_facet[fa][0];
        LegendrePolynomial::Eval(p, sigma[fmax] - sigma[ftrig],leg_u);
        LegendrePolynomial::Eval(p, sigma[fmax] - sigma[fz],leg_v);

        for(int k = 0; k <= p; k++)
          for(int l = 0; l <= p; l++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(0.5*(sigma[fmax] - sigma[fz]),-0.5*(sigma[fmax] - sigma[fz]),
              0.5*(sigma[fmax] - sigma[ftrig]),-0.5*(sigma[fmax] - sigma[ftrig]),leg_u[l]* leg_v[k]*lam_face);
          }



      }



      int oi = order_inner[0];
      leg_u.SetSize(oi+1);
      leg_v.SetSize(oi+1);
      leg_w.SetSize(oi+1);

      LegendrePolynomial::Eval(oi+1,sigma[0] - sigma[1],leg_u);
      LegendrePolynomial::Eval(oi+1,sigma[0] - sigma[3],leg_v);
      LegendrePolynomial::Eval(oi+1,sigma[0] - sigma[4],leg_w);

      for(int k=0; k<=oi-1; k++)
      {
        for(int i = 0; i <= oi+1; i++)
        {
          for(int j = 0; j <= oi+1; j++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[0], lx[1], lz[0], lz[1], ly[0]*ly[1]*leg_u[i]*leg_v[k]* leg_w[j]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[0], lx[1], ly[0], ly[1], lz[0]*lz[1]*leg_u[i]*leg_v[j]* leg_w[k]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(ly[0], ly[1], lz[0], lz[1], lx[0]*lx[1]*leg_u[k]*leg_v[j]* leg_w[i]);
          }
        }
      }


      // S_xz
      for (int i=0; i<=oi; i++)
      {
        for (int j=0; j<=oi; j++)
        {
          for (int k=0; k<=oi+1; k++)
          {
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lx[0], lx[0], ly[0], lz[0], leg_u[k]*leg_v[j]*leg_w[i]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(ly[0], ly[0], lx[0], lz[0], leg_u[i]*leg_v[k]*leg_w[j]);
            shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(lz[0], lz[0], lx[0], ly[0], leg_u[i]*leg_v[j]*leg_w[k]);
          }
        }
      }


    };

  };

  */


  ////////////////////// SURFACE ////////////////////////////
    template <int DIM>
  class HCurlCurlSurfaceFiniteElement : public FiniteElement
  {
  public:
    using FiniteElement::FiniteElement;
    using FiniteElement::ndof;
    using FiniteElement::order;

    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM+1> & mip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedShape_Vector (const MappedIntegrationPoint<DIM,DIM+1> & mip,
      BareSliceMatrix<double> shape) const = 0;

  };


  template <ELEMENT_TYPE ET> class HCurlCurlSurfaceFE;

  
  template <ELEMENT_TYPE ET>
  class T_HCurlCurlSurfaceFE : public HCurlCurlSurfaceFiniteElement<ET_trait<ET>::DIM>,
    public VertexOrientedFE<ET>
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };
    enum { DIM_STRESS = ((DIM+2)*(DIM+1))/2 };
    
    using VertexOrientedFE<ET>::vnums;
    using HCurlCurlSurfaceFiniteElement<ET_trait<ET>::DIM>::ndof;
    using HCurlCurlSurfaceFiniteElement<ET_trait<ET>::DIM>::order;

    INT<DIM> order_inner;


  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    
    T_HCurlCurlSurfaceFE (int aorder)
    {
      order = aorder;
      order_inner = aorder;
    }
    
    virtual ELEMENT_TYPE ElementType() const { return ET; }
    const HCurlCurlSurfaceFE<ET> * Cast() const { return static_cast<const HCurlCurlSurfaceFE<ET>*> (this); } 
    
    INLINE void SetOrderInner (INT<DIM,int> order) { order_inner = order; }

    virtual void ComputeNDof()
    {
      cout << "Error, T_HCurlCurlSurfaceFE<ET>:: ComputeNDof not available for base class" << endl;
    }

    virtual void CalcMappedShape_Vector (const MappedIntegrationPoint<DIM,DIM+1> & mip,
                            BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiff<DIM+1>> adp = mip;
      Vec<DIM, AutoDiffDiff<DIM+1>> addp;
      for (int i=0; i<DIM+1; i++)
      {
        addp[i] = adp[i].Value();
        addp[i].LoadGradient(&adp[i].DValue(0));
      }
      Cast() -> T_CalcShape (TIP<DIM, AutoDiffDiff<DIM+1>> (addp), SBLambda([shape] (int nr, auto val)
                                          {
                                            shape.Row(nr).AddSize(DIM_STRESS) = val;
                                          }));
    }


    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM+1> & mip,
                            BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiff<DIM+1>> adp = mip;
      Vec<DIM, AutoDiffDiff<DIM+1>> addp;
      for (int i=0; i<DIM+1; i++)
      {
        addp[i] = adp[i].Value();
        addp[i].LoadGradient(&adp[i].DValue(0));
      }
      Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM+1>> (addp),SBLambda([shape](int nr, auto val)//Capture
      {
        Vec<DIM_STRESS> vecshape = val;
        BareVector<double> matshape = shape.Row(nr);
        VecToSymMat<DIM+1> (vecshape, matshape);
      }));
    }


  };

  template <> class HCurlCurlSurfaceFE<ET_SEGM> : public T_HCurlCurlSurfaceFE<ET_SEGM> 
  {
    
  public:
    using T_HCurlCurlSurfaceFE<ET_SEGM> :: T_HCurlCurlSurfaceFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      ndof += order_inner[0]+1;
      order = max2(order,order_inner[0]);

    }
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<1,Tx> ip, TFA & shape) const
    {
      auto x = ip.x;
      AutoDiffDiff<2> ddlami[2] ={ x, 1-x };
      
      int ii = 0;

      int es = 0,ee = 1;
      if(vnums[es] > vnums[ee]) swap (es,ee);

      AutoDiffDiff<2> ls = ddlami[es],le = ddlami[ee];
   
      Vec<3> symdyadic = Vec<3>(2*ls.DValue(0)*le.DValue(0),2*ls.DValue(1)*le.DValue(1), ls.DValue(1)*le.DValue(0)+ls.DValue(0)*le.DValue(1));
          
      LegendrePolynomial::EvalScaled(order, ls.Value()-le.Value(),ls.Value()+le.Value(),
                                     SBLambda([shape, &ii,symdyadic] (size_t nr, auto val)
                                              {
                                                auto value = val*symdyadic;
                                                shape[ii] = value;
                                                ii++;
                                              }));
      
    };
  };


  template <> class HCurlCurlSurfaceFE<ET_TRIG> : public T_HCurlCurlSurfaceFE<ET_TRIG> 
  {
    
  public:
    using T_HCurlCurlSurfaceFE<ET_TRIG> :: T_HCurlCurlSurfaceFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      ndof += 3*(order_inner[0]+1)*(order_inner[0]+2)/2;
      order = max2(order, order_inner[0]);
    }


    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y;
      typedef decltype(x.Value()+x.Value()) T;                  
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> lx[6] ={ xx, yy, 1-xx-yy};
      int ii = 0;

      /*const EDGE * edges = ElementTopology::GetEdges(ET_TRIG);
      for (int i = 0; i < 3; i++)
        {
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es,ee);
          auto ls = lx[es], le = lx[ee];
          Vec<6> symdyadic =  SymDyadProd(ls,le);
          LegendrePolynomial::EvalScaled(order, ls-le,ls+le,
                                         SBLambda([&] (size_t nr, auto val)
                                                  {
                                                    shape[ii++] = val.Value()*symdyadic;
                                                  }));
                                                  }*/
      //const FACE * faces = ElementTopology::GetFaces(ET_TET);
      //int fav[3] ={0,1,2};
      //Sort vertices  first edge op minimal vertex
      //if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
      //if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
      //if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
        

      for( int i = 0; i < 3; i++)
        {
          INT<2> e = ET_trait<ET_TRIG>::GetEdgeSort (i, vnums);
          auto ls = lx[e[0]], le = lx[e[1]];
          Vec<6> symdyadic =  SymDyadProd(ls,le);
              
          LegendrePolynomial::EvalScaled(order, ls-le,ls+le,
                                         SBLambda([&] (size_t nr, auto val)
                                                  {
                                                    auto value = val.Value()*symdyadic;
                                                    shape[ii] = value;
                                                    ii++;
                                                  }));
        }


      if (order > 0)
        {
          INT<4> f =  ET_trait<ET_TRIG>::GetFaceSort(0, vnums);
          
          auto ls = lx[f[0]], le = lx[f[1]], lt = lx[f[2]];
          Vec<6> symdyadic1 = lt.Value()*SymDyadProd(ls,le);
          Vec<6> symdyadic2 = ls.Value()*SymDyadProd(lt,le);
          Vec<6> symdyadic3 = le.Value()*SymDyadProd(ls,lt);
          
          DubinerBasis3::Eval(order-1, ls,le,
                              SBLambda([&] (size_t nr, auto val)
                                       {
                                         auto value = val.Value()*symdyadic1;
                                         shape[ii] = value;
                                         ii++;
                                         value = val.Value()*symdyadic2;
                                         shape[ii] = value;
                                         ii++;
                                         value = val.Value()*symdyadic3;
                                         shape[ii] = value;
                                         ii++;
                                       }));
        }
     
    }
  };


  /*template <> class HCurlCurlSurfaceFE<ET_QUAD> : public T_HCurlCurlSurfaceFE<ET_QUAD> 
  {
    
  public:
    using T_HCurlCurlSurfaceFE<ET_QUAD> :: T_HCurlCurlSurfaceFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      ndof += (order_inner[0]+1+HCurlCurlFE<ET_PRISM>::incrorder_xx1_bd)*(order_inner[1]+1+HCurlCurlFE<ET_PRISM>::incrorder_xx2_bd);
      order = max2(order, order_inner[0]+HCurlCurlFE<ET_PRISM>::incrorder_xx1_bd);
      order = max2(order, order_inner[1]+HCurlCurlFE<ET_PRISM>::incrorder_xx2_bd);
    }


    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      AutoDiffDiff<3> x = ip.x, z = ip.y;
      typedef decltype(x.Value()+x.Value()) T;                  
      AutoDiff<3> xx(x.Value(), &x.DValue(0));
      AutoDiff<3> zz(z.Value(), &z.DValue(0));
      AutoDiff<3> sigma[4] = {1-xx+1-zz, xx+1-zz, xx+zz, 1-xx+zz};
      int ii = 0;
      

      ArrayMem<AutoDiff<3>,20> leg_u(order_inner[0]+2);
      ArrayMem<AutoDiff<3>,20> leg_w(order_inner[1]+2);

      
      int fmax = 0;
      for(int j = 1; j < 4; j++)
        if(vnums[j] > vnums[fmax]) fmax = j;

      int f1, f2;
      f1 = (fmax+1)%4;
      f2 = (fmax+3)%4;


      if(vnums[f1] > vnums[f2])
      {
        swap(f1,f2);
      }

      LegendrePolynomial::Eval(order_inner[0],sigma[fmax] - sigma[f1],leg_u);
      LegendrePolynomial::Eval(order_inner[0],sigma[fmax] - sigma[f2],leg_w);

      for(int k = 0; k <= order_inner[0]+HCurlCurlFE<ET_PRISM>::incrorder_xx2_bd; k++)
        for(int l = 0; l <= order_inner[0]+HCurlCurlFE<ET_PRISM>::incrorder_xx1_bd; l++)
        {
          shape[ii++] = Prism_Dl1xDl3_symtensor_Dl2xDl4_u<T>(0.5*(sigma[fmax]-sigma[f2]),-0.5*(sigma[fmax]-sigma[f2]),
            0.5*(sigma[fmax]-sigma[f1]),-0.5*(sigma[fmax]-sigma[f1]),leg_u[l]* leg_w[k]);
        }
                

    }
    };*/


  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_TRIG>;
  //HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_QUAD>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_TET>;
  //HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_PRISM>;
  //HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_HEX>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlSurfaceFE<ET_SEGM>;
  //HCURLCURLFE_EXTERN template class T_HCurlCurlSurfaceFE<ET_QUAD>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlSurfaceFE<ET_TRIG>;

}


#endif
  
