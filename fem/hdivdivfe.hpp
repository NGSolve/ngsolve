#ifndef FILE_HDIVDIVFE
#define FILE_HDIVDIVFE

/*********************************************************************/
/* File:   hdivdivfe.hpp                                             */
/* Author: Astrid Pechstein, Joachim Schoeberl                       */
/* Date:   orig 2006, redesign Dec 2016                              */
/*********************************************************************/


namespace ngfem
{

  
  class HDivDivFiniteElement : public FiniteElement
  {
  public:
    using FiniteElement::FiniteElement;
    
    // sigma_xx, sigma_yy, sigma_xy
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const = 0;

    virtual void CalcDivShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<double> divshape) const = 0;

    /*
    virtual void CalcDivDivShape (const IntegrationPoint & ip, 
                                  FlatVector<> ddshape) const;
    */
  };
  


  template <ELEMENT_TYPE ET> class HDivDivFE;

  
  template <ELEMENT_TYPE ET>
  class T_HDivDivFE : public HDivDivFiniteElement
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };
    enum { DIM_STRESS = (DIM*(DIM+1))/2 };
    
    enum { N_VERTEX = ET_trait<ET>::N_VERTEX };
    enum { N_FACET   = ET_trait<ET>::N_FACET };    
    
    size_t vnums[N_VERTEX];
    INT<DIM-1> order_facet[N_FACET];
    INT<DIM> order_inner;

    // additional div-div free bubbles
    bool plus;

  public:

    T_HDivDivFE (int aorder, bool _plus = false)
      : plus(_plus)
    {
      order = aorder;
      for (auto & of : order_facet) of = aorder;
      order_inner = aorder;
      ndof = DIM*(DIM+1)/2 * ET_trait<ET>::PolDimension(aorder);

      if (plus)
        {
          order++;
          ndof += 2*aorder;
        }
    }
    
    virtual ELEMENT_TYPE ElementType() const { return ET; }
    const HDivDivFE<ET> * Cast() const { return static_cast<const HDivDivFE<ET>*> (this); } 
    
    template <typename TA> 
    void SetVertexNumbers (const TA & avnums)
    { for (int i = 0; i < N_VERTEX; i++) vnums[i] = avnums[i]; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const
    {
      AutoDiffDiff<2> hx[] = { AutoDiffDiff<2>(ip(0), 0), AutoDiffDiff<2>(ip(1), 1) } ;

      Cast() -> T_CalcShape (hx, SBLambda([&] (int nr, auto val)
                                          {
                                            shape.Row(nr).AddSize(DIM_STRESS) = val.Shape();
                                          }));
    }


    virtual void CalcDivShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<double> shape) const
    {
      AutoDiffDiff<2> hx[] = { AutoDiffDiff<2>(ip(0), 0), AutoDiffDiff<2>(ip(1), 1) } ;
      
      Cast() -> T_CalcShape (hx, SBLambda([&] (int nr, auto val)
                                          {
                                            shape.Row(nr).AddSize(DIM) = val.DivShape();
                                          }));
    }

  };


  
  // ***************** SigmaGrad ****************************** */
  // sigma (nabla u)
  
  template <int D> class T_SigmaGrad;
  template <> class T_SigmaGrad<2>
  {
    AutoDiffDiff<2> u;
  public:
    T_SigmaGrad  (AutoDiffDiff<2> au) : u(au) { ; }
    Vec<3> Shape() { return Vec<3> (u.DDValue(1,1), u.DDValue(0,0), -u.DDValue(1,0)); }
    Vec<2> DivShape() { return Vec<2> (0.0, 0.0); }
  };
  
  template <int D>
  auto SigmaGrad (AutoDiffDiff<D> au) { return T_SigmaGrad<D>(au); }
  
  
  // ***************** Sigma_u_Gradv ****************************** */
  // sigma (u nabla v)
  
  template <int D> class T_Sigma_u_Gradv;
  template <> class T_Sigma_u_Gradv<2>
  {
    AutoDiffDiff<2> u, v;
  public:
    T_Sigma_u_Gradv  (AutoDiffDiff<2> au, AutoDiffDiff<2> av) : u(au), v(av) { ; }
    Vec<3> Shape() { return Vec<3> (u.Value()*v.DDValue(1,1) + u.DValue(1)*v.DValue(1),
                                    u.Value()*v.DDValue(0,0) + u.DValue(0)*v.DValue(0),
                                    -u.Value()*v.DDValue(1,0) - 0.5 * (u.DValue(0)*v.DValue(1)+u.DValue(1)*v.DValue(0))); }
    Vec<2> DivShape()
    {
      double uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      double ux = u.DValue(0), uy = u.DValue(1);
      double vxx = v.DDValue(0,0), vyy = v.DDValue(1,1), vxy = v.DDValue(0,1);
      double vx = v.DValue(0), vy = v.DValue(1);
      
      return -0.5 * Vec<2> (uyy*vx - uxy*vy + uy*vxy - ux*vyy,
                            -uxy*vx + uxx*vy - uy*vxx + ux*vxy);
    }
  };
  
  template <int D>
  auto Sigma_u_Gradv (AutoDiffDiff<D> au, AutoDiffDiff<D> av) { return T_Sigma_u_Gradv<D>(au, av); }
  
  
  // ***************** Type2 ****************************** */
  // ????
  
  template <int D> class T_Type2;
  template <> class T_Type2<2>
  {
    AutoDiffDiff<2> u,v;
  public:
    T_Type2  (AutoDiffDiff<2> au, AutoDiffDiff<2> av) : u(au), v(av) { ; }
    Vec<3> Shape() { return Vec<3> (u.DDValue(1,1)*v.Value() - 2*u.DValue(1)*v.DValue(1) + u.Value()*v.DDValue(1,1),
                                    u.DDValue(0,0)*v.Value() - 2*u.DValue(0)*v.DValue(0) + u.Value()*v.DDValue(0,0),
                                    -(u.DDValue(0,1)*v.Value() - u.DValue(0)*v.DValue(1) -
                                      u.DValue(1)*v.DValue(0) + u.Value()*v.DDValue(1,0))); }

    Vec<2> DivShape()
    {
      double uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      double ux = u.DValue(0), uy = u.DValue(1);
      double vxx = v.DDValue(0,0), vyy = v.DDValue(1,1), vxy = v.DDValue(0,1);
      double vx = v.DValue(0), vy = v.DValue(1);
      
      return Vec<2> (2*uyy*vx + 2*ux*vyy - 2*uxy*vy - 2*uy*vxy, 2*uxx*vy + 2*uy*vxx - 2*uxy*vx - 2*ux*vxy);
    }
  };
  
  template <int D>
  auto Type2 (AutoDiffDiff<D> au, AutoDiffDiff<D> av) { return T_Type2<D>(au, av); }
  
  // ***************** Type3 ****************************** */
  // ????
  
  template <int D> class T_Type3;
  template <> class T_Type3<2>
  {
    AutoDiffDiff<2> u,v;
  public:
    T_Type3  (AutoDiffDiff<2> au, AutoDiffDiff<2> av) : u(au), v(av) { ; }
    Vec<3> Shape() { return Vec<3> (u.DDValue(1,1)*v.Value() - u.Value()*v.DDValue(1,1),
                                    u.DDValue(0,0)*v.Value() - u.Value()*v.DDValue(0,0),
                                    -(u.DDValue(0,1)*v.Value() - u.Value()*v.DDValue(1,0))); }
    Vec<2> DivShape()
    {
      double uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      double ux = u.DValue(0), uy = u.DValue(1);
      double vxx = v.DDValue(0,0), vyy = v.DDValue(1,1), vxy = v.DDValue(0,1);
      double vx = v.DValue(0), vy = v.DValue(1);

      return Vec<2> (uyy*vx - uxy*vy - ux*vyy + uy*vxy, uxx*vy - uxy*vx - uy*vxx + ux*vxy);
    }
  };
  
  template <int D>
  auto Type3 (AutoDiffDiff<D> au, AutoDiffDiff<D> av) { return T_Type3<D>(au, av); }

  
  // ***************** Sigma ((Duv - vDu) w) ****************************** */
  // where u, v are linear ????
  template <int D> class T_Sigma_Duv_minus_uDv_w;
  template <> class T_Sigma_Duv_minus_uDv_w<2>
  {
    AutoDiffDiff<2> u,v,w;
  public:
    T_Sigma_Duv_minus_uDv_w  (AutoDiffDiff<2> au, AutoDiffDiff<2> av, AutoDiffDiff<2> aw) : u(au), v(av), w(aw) { ; }
    Vec<3> Shape() { return Vec<3> ( u.DValue(1)*w.Value() + u.Value()*w.DValue(1),
                                     v.DValue(0)*w.Value() + v.Value()*w.DValue(0),
                                     - (u.DValue(0)*w.Value() + u.Value()*w.DValue(0) +
                                        v.DValue(1)*w.Value() + v.Value()*w.DValue(1)) * 0.5
                                     ); }
    Vec<2> DivShape()
    {
      return Vec<2> ( 0.5*(u.DDValue(0,1)*w.Value() + u.DValue(1)*w.DValue(0) + u.DValue(0)*w.DValue(1) + u.Value()*w.DDValue(0,1)
                           - v.DDValue(1,1)*w.Value() - 2*v.DValue(1)*w.DValue(1) - v.Value()*w.DDValue(1,1)),
                      0.5*(v.DDValue(0,1)*w.Value() + v.DValue(1)*w.DValue(0) + v.DValue(0)*w.DValue(1) + v.Value()*w.DDValue(0,1)
                           - u.DDValue(0,0)*w.Value() - 2*u.DValue(0)*w.DValue(0) - u.Value()*w.DDValue(0,0)));
    }
  };
  
  template <int D>
  auto Sigma_Duv_minus_uDv_w (AutoDiffDiff<D> au, AutoDiffDiff<D> av, AutoDiffDiff<D> aw)
  { return T_Sigma_Duv_minus_uDv_w<D>(au, av, aw); }
  
  


  
  template <> class HDivDivFE<ET_TRIG> : public T_HDivDivFE<ET_TRIG> 
  {
    typedef T_HDivDivFE<ET_TRIG> BASE;
    
  public:
    using T_HDivDivFE<ET_TRIG> :: T_HDivDivFE;

    template <typename TFA> 
    void T_CalcShape (AutoDiffDiff<2> hx[2], TFA & shape) const
    {
      auto x = hx[0], y = hx[1];
      AutoDiffDiff<2> ddlami[3] ={ x, y, 1-x-y };
      
      int ii = 0;
      
      int maxorder_facet =
        max2(order_facet[0][0],max2(order_facet[1][0],order_facet[2][0]));

      const EDGE * edges = ElementTopology::GetEdges(ET_TRIG);

      ArrayMem<AutoDiffDiff<2>,20> ha(maxorder_facet+1);
      ArrayMem<AutoDiffDiff<2>,20> u(order_inner[0]+1), v(order_inner[0]+1);
      
      for (int i = 0; i < 3; i++)
        {
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es,ee);
          
          AutoDiffDiff<2> ls = ddlami[es], le = ddlami[ee];
          
          // edge functions are all div-free!
          IntegratedLegendreMonomialExt::CalcTrigExt(maxorder_facet+2,
                                                     le-ls, 1-le-ls, ha);
          
          for (int l = 0; l <= order_facet[i][0]; l++)
            shape[ii++] = SigmaGrad (ha[l]);
        }
      
      int es = 0; int ee = 1; int et = 2;
      AutoDiffDiff<2> ls = ddlami[es];
      AutoDiffDiff<2> le = ddlami[ee];
      AutoDiffDiff<2> lt = ddlami[et];
      
      IntegratedLegendreMonomialExt::CalcTrigExt(order_inner[0]+2,le-ls,1-le-ls,u);
      LegendrePolynomial::EvalMult(order_inner[0], 2*lt-1, lt, v);
      
      int oi=order_inner[0];
      
      for (int i = 0; i <= oi-1; i++)
        for (int j = 0; j+i <= oi-1; j++)
          {
            shape[ii++] = SigmaGrad(u[i]*v[j]);
            shape[ii++] = Type2(u[i], v[j]);
            if (j > 0)
              shape[ii++] = Type3(u[i], v[j]);
          }

      AutoDiffDiff<2> phie0 = le, phie1 = -ls;
      for (int i = 0; i < oi; i++)
        shape[ii++] = Sigma_Duv_minus_uDv_w (phie1, phie0, v[i]);
      
      // element bubbles for Sigma+ space
      if (plus)
        for (int i = 0; i <= oi-1; i++)
          {
            AutoDiffDiff<2> bubble = u[i]*v[oi-1-i];
            shape[ii++] = Sigma_u_Gradv(bubble, x);
            shape[ii++] = Sigma_u_Gradv(bubble, y);
          }
    };
  };
  

}


#endif
  
