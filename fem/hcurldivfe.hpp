#ifndef FILE_HCURLDIVFE
#define FILE_HCURLDIVFE

/*********************************************************************/
/* File:   hcurldivfe.hpp                                            */
/* Author: Philip Lederer                                            */
/* Date:   2017/2018                                                 */
/*********************************************************************/


namespace ngfem
{

  template <int DIM>
  class HCurlDivFiniteElement : public FiniteElement
  {
  public:
    using FiniteElement::FiniteElement;
    using FiniteElement::ndof;
    using FiniteElement::order;

    // old style
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const = 0;

    virtual void CalcDivShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<double> divshape) const = 0;

    // new implementation
    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM> & mip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedShape_Vector (const MappedIntegrationPoint<DIM,DIM> & mip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedDivShape (const MappedIntegrationPoint<DIM,DIM> & mip,
      BareSliceMatrix<double> shape) const = 0;

  };

  template <int D,typename VEC,typename MAT>
  void VecToMat(const VEC & vec,MAT & mat)
  {
    switch(D)
    {
    case 2:
      mat(0) = vec(0);
      mat(1) = vec(1);
      mat(2) = vec(2);
      mat(3) = vec(3);
      break;
    case 3:
      mat(0) = vec(0);
      mat(1) = vec(1);
      mat(2) = vec(2);
      mat(3) = vec(3);
      mat(4) = vec(4);
      mat(5) = vec(5);
      mat(6) = vec(6);
      mat(7) = vec(7);
      mat(8) = vec(8);
      break;
    }

  }

  template <ELEMENT_TYPE ET> class HCurlDivFE;

  
  template <ELEMENT_TYPE ET>
  class T_HCurlDivFE : public HCurlDivFiniteElement<ET_trait<ET>::DIM>,
    public VertexOrientedFE<ET>
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };
    enum { DIM_STRESS = DIM*DIM };
    
    using VertexOrientedFE<ET>::vnums;
    using HCurlDivFiniteElement<ET_trait<ET>::DIM>::ndof;
    using HCurlDivFiniteElement<ET_trait<ET>::DIM>::order;

    INT<DIM-1> order_facet[ET_trait<ET>::N_FACET];
    INT<DIM> order_inner;

    // additional curl-div free bubbles
    bool plus;

  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    
    T_HCurlDivFE (int aorder, bool _plus = false)
      : plus(_plus)
    {
      order = aorder;
      for (auto & of : order_facet) of = aorder;
      order_inner = aorder;
    }
    
    virtual ELEMENT_TYPE ElementType() const { return ET; }
    const HCurlDivFE<ET> * Cast() const { return static_cast<const HCurlDivFE<ET>*> (this); } 
    
    INLINE void SetOrderFacet (int nr, INT<DIM-1,int> order) { order_facet[nr] = order; }
    INLINE void SetOrderInner (INT<DIM,int> order) { order_inner = order; }

    virtual void ComputeNDof()
    {
      cout << "Error, T_HCurlDivFE<ET>:: ComputeNDof not available, only for ET == TRIG" << endl;
    }

    // old style
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiffDiff<DIM>> adp;
      for ( int i=0; i<DIM; i++)
      {
        adp(i) = AutoDiffDiff<DIM>(ip(i),i);
      }

      Cast() -> T_CalcShape (TIP<DIM, AutoDiffDiff<DIM>> (adp), SBLambda([&] (int nr, auto val)
                                          {
                                            shape.Row(nr).AddSize(DIM_STRESS) = val.Shape();
                                          }));
    }

    virtual void CalcDivShape (const IntegrationPoint & ip,
                               BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiffDiff<DIM>> adp;
      for ( int i=0; i<DIM; i++)
      {
        adp[i] = AutoDiffDiff<DIM>(ip(i),i);
      }
      
      Cast() -> T_CalcShape (TIP<DIM, AutoDiffDiff<DIM>> (adp), SBLambda([&] (int nr, auto val)
                                          {
                                            shape.Row(nr).AddSize(DIM) = val.DivShape();
                                          }));
    }

    // new style
    virtual void CalcMappedShape_Vector (const MappedIntegrationPoint<DIM,DIM> & mip,
                            BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiff<DIM>> adp = mip;
      Vec<DIM, AutoDiffDiff<DIM>> addp;
      for (int i=0; i<DIM; i++)
      {
        addp[i] = adp[i].Value();
        addp[i].LoadGradient(&adp[i].DValue(0));
      }
      Cast() -> T_CalcShape (TIP<DIM, AutoDiffDiff<DIM>> (addp), SBLambda([&] (int nr, auto val)
                                          {
                                            shape.Row(nr).AddSize(DIM_STRESS) = val.Shape();
                                          }));
    }


    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM> & mip,
                            BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiff<DIM>> adp = mip;
      Vec<DIM, AutoDiffDiff<DIM>> addp;
      for (int i=0; i<DIM; i++)
      {
        addp[i] = adp[i].Value();
        addp[i].LoadGradient(&adp[i].DValue(0));
      }
      Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM>> (addp),SBLambda([&](int nr,auto val)
      {
        Vec<DIM_STRESS> vecshape = val.Shape();
        BareVector<double> matshape = shape.Row(nr);
        VecToMat<DIM> (vecshape, matshape);
      }));
    }


    virtual void CalcMappedDivShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                            BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiff<DIM>> adp = mip;
      Vec<DIM, AutoDiffDiff<DIM>> addp;
      for (int i=0; i<DIM; i++)
      {
        addp[i] = adp[i].Value();
        addp[i].LoadGradient(&adp[i].DValue(0));
      }

      if(!mip.GetTransformation().IsCurvedElement()) // non-curved element
      {
        Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM>> (addp),SBLambda([&](int nr,auto val)
        {
          shape.Row(nr).AddSize(DIM) = val.DivShape();
        }));
      }
      else // curved element
      {
	throw Exception ("not implemented yet!");
        /*Mat<DIM> jac = mip.GetJacobian();
        Mat<DIM> inv_jac = mip.GetJacobianInverse();
        Mat<DIM> hesse[3],finvT_h_tilde_finv[3];
        mip.CalcHesse (hesse[0],hesse[1],hesse[2]);

        Mat<DIM,DIM,AutoDiff<DIM> > f_tilde;
        for(int i = 0; i < DIM; i++)
        {
          for(int j = 0; j < DIM; j++)
          {
            f_tilde(i,j).Value() = jac(i,j);
            for(int k = 0; k < DIM; k++)
              f_tilde(i,j).DValue(k) = hesse[i](j,k);
          }
        }

        AutoDiff<DIM> ad_det = Det (f_tilde);
        AutoDiff<DIM> iad_det = 1.0 / ad_det;
        f_tilde *= iad_det;

        for(int i=0; i<DIM; i++)
        {
          finvT_h_tilde_finv[i] = 0;
          for(int alpha=0; alpha<DIM; alpha++)
            for(int beta=0; beta<DIM; beta++)
              for(int gamma=0; gamma<DIM; gamma++)
                for(int delta=0; delta<DIM; delta++)
                  finvT_h_tilde_finv[i](alpha,beta) += inv_jac(gamma,alpha)*f_tilde(i,gamma).DValue(delta)*inv_jac(delta,beta);
        }

        Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM>> (addp),SBLambda([&](int nr,auto val)
                                  {
                                    shape.Row(nr).AddSize(DIM) = val.DivShape();
                                    BareVector<double> divshape = shape.Row(nr);
                                    Vec<DIM_STRESS> vecshape = val.Shape();
                                    Vec<DIM*DIM> matshape;
                                    VecToSymMat<DIM> (vecshape, matshape);

                                    for(int k=0; k<DIM; k++)
                                    {
                                      for(int j=0; j<DIM*DIM; j++)
                                      {
                                        divshape(k) += mip.GetJacobiDet() * finvT_h_tilde_finv[k](j) * matshape(j);
                                      }
                                    }
                                    
                                  }));
	*/
      }
    }

  };
  
  /* calculates [(grad l1) o-times (rot-grad l2) ] * intlegendre */
  class T_Dl2xRotDl1_v
  {
    AutoDiffDiff<2> l1,l2,v;
  public:
    T_Dl2xRotDl1_v  (AutoDiffDiff<2> lam1, AutoDiffDiff<2> lam2, AutoDiffDiff<2> av) : l1(lam1), l2(lam2), v(av) { ; }

    //Shape returns (sig_xx, sig_xy, sig_yx, sig_yy)
    Vec<4> Shape() {
      return Vec<4> (-v.Value()*(l1.DValue(0)*l2.DValue(1)),
		     v.Value()*(l1.DValue(0)*l2.DValue(0)),
		     -v.Value()*(l1.DValue(1)*l2.DValue(1)),
		     v.Value()*(l1.DValue(1)*l2.DValue(0))
		     );
    }

    Vec<2> DivShape()
    {
      double vx = v.DValue(0), vy = v.DValue(1);

      double lam1x = l1.DValue(0);
      double lam1y = l1.DValue(1);
      double lam1xx = l1.DDValue(0,0);
      double lam1xy = l1.DDValue(0,1);
      double lam1yx = l1.DDValue(1,0);
      double lam1yy = l1.DDValue(1,1);
      
      double lam2x = l2.DValue(0);
      double lam2y = l2.DValue(1);
      double lam2xx = l2.DDValue(0,0);
      double lam2xy = l2.DDValue(0,1);
      double lam2yx = l2.DDValue(1,0);
      double lam2yy = l2.DDValue(1,1);
      
      return Vec<2> (
		     v.Value() * ( - (lam1xx*lam2y + lam1x*lam2yx) + (lam1yx*lam2x + lam1x*lam2xy)) - vx *lam1x*lam2y + vy*lam1x*lam2x, 
		     v.Value() * ( - (lam1xy*lam2y + lam1y*lam2yx) + (lam1yy*lam2x + lam1y*lam2xy)) - vx *lam1y*lam2y + vy*lam1y*lam2x
        ); 
    }

  };

  
  /* calculates [(grad l1) o-times (rot-grad l2) - (grad l2) o-times (rot-grad l1)  ] * intlegendre */
  class T_Dl2xRotDl1_minus_Dl1xRotDl2_v
  {
    AutoDiffDiff<2> l1,l2,v;
  public:
    T_Dl2xRotDl1_minus_Dl1xRotDl2_v  (AutoDiffDiff<2> lam1, AutoDiffDiff<2> lam2, AutoDiffDiff<2> av) : l1(lam1), l2(lam2), v(av) { ; }

    //Shape returns (sig_xx, sig_xy, sig_yx, sig_yy)
    Vec<4> Shape() {
      return Vec<4> (-v.Value()*(l1.DValue(0)*l2.DValue(1) - l2.DValue(0)*l1.DValue(1) ) ,
		     v.Value()*(l1.DValue(0)*l2.DValue(0) - l2.DValue(0)*l1.DValue(0)),
		     -v.Value()*(l1.DValue(1)*l2.DValue(1) -  l2.DValue(1)*l1.DValue(1)),
		     v.Value()*(l1.DValue(1)*l2.DValue(0) - l2.DValue(1)*l1.DValue(0))
		     );
    }

    Vec<2> DivShape()
    {
      double vx = v.DValue(0), vy = v.DValue(1);

      double lam1x = l1.DValue(0);
      double lam1y = l1.DValue(1);
      double lam1xx = l1.DDValue(0,0);
      double lam1xy = l1.DDValue(0,1);
      double lam1yx = l1.DDValue(1,0);
      double lam1yy = l1.DDValue(1,1);
      
      double lam2x = l2.DValue(0);
      double lam2y = l2.DValue(1);
      double lam2xx = l2.DDValue(0,0);
      double lam2xy = l2.DDValue(0,1);
      double lam2yx = l2.DDValue(1,0);
      double lam2yy = l2.DDValue(1,1);
      
      return Vec<2> (
	v.Value() * ( - (lam1xx*lam2y + lam1x*lam2yx - lam2xx*lam1y - lam2x*lam1yx)
		      + (lam1yx*lam2x + lam1x*lam2xy - lam2yx*lam1x - lam2x*lam1xy)) - vx * (lam1x*lam2y - lam2x*lam1y) + vy*(lam1x*lam2x - lam2x*lam1x), 
      	     v.Value() * ( - (lam1xy*lam2y + lam1y*lam2yx - lam2xy*lam1y - lam2y*lam1yx)
			   + (lam1yy*lam2x + lam1y*lam2xy - lam2yy*lam1x - lam2y*lam1xy)) - vx *(lam1y*lam2y - lam2y*lam1y) + vy*(lam1y*lam2x - lam2y*lam1x)
       ); 
    }

  };

  

  template <> class HCurlDivFE<ET_TRIG> : public T_HCurlDivFE<ET_TRIG> 
  {
    
  public:
    using T_HCurlDivFE<ET_TRIG> :: T_HCurlDivFE;

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
      if (plus)
      { 
        order ++;
        ninner += 2*order_inner[0]; 
      }
      ndof += ninner+1;

    }
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {
      auto x = ip.x, y = ip.y;
      AutoDiffDiff<2> ddlami[3] ={ x, y, 1-x-y };
      
      int ii = 0;
      
      int maxorder_facet =
        max2(order_facet[0][0],max2(order_facet[1][0],order_facet[2][0]));

      const EDGE * edges = ElementTopology::GetEdges(ET_TRIG);

      ArrayMem<AutoDiffDiff<2>,20> ha(maxorder_facet+1);
      ArrayMem<AutoDiffDiff<2>,20> u(order_inner[0]+2), v(order_inner[0]+2);
      
      for (int i = 0; i < 3; i++)
        {
          int es = edges[i][0], ee = edges[i][1];
	  //cout<<"Edge i = "<<i<<endl;
	  
          if (vnums[es] > vnums[ee]) swap (es,ee);

	  //cout<<"es = "<< es<<endl;
	  //cout<<"ee = "<< ee<<endl;
	  
          AutoDiffDiff<2> ls = ddlami[es], le = ddlami[ee];

	  // Ask joachim about this magic !!!
          // edge functions are all div-free!
          //IntegratedLegendreMonomialExt::CalcTrigExt(maxorder_facet,le-ls, 1-le-ls, ha);
	  
          ScaledLegendrePolynomial(maxorder_facet,le-ls, 1-le-ls,ha);

          
          for (int l = 0; l <= order_facet[i][0]; l++)
	    shape[ii++] =  T_Dl2xRotDl1_v(le, ls, ha[l]);
            
        }
      
      int es = 0; int ee = 1; int et = 2;
      AutoDiffDiff<2> ls = ddlami[es];
      AutoDiffDiff<2> le = ddlami[ee];
      AutoDiffDiff<2> lt = ddlami[et];
      
      //int oi=order_inner[0];
      //int oi_plus = oi; //plus ? oi+1 : oi;

      ScaledLegendrePolynomial(0,le-ls, 1-le-ls,ha);
      //IntegratedLegendreMonomialExt::CalcTrigExt(maxorder_facet,le-ls, 1-le-ls, ha);
      
      shape[ii++] =  T_Dl2xRotDl1_minus_Dl1xRotDl2_v(le, ls, ha[0]);
      
      
      //IntegratedLegendreMonomialExt::CalcTrigExt(oi_plus+3,le-ls,1-le-ls,u);
      //LegendrePolynomial::EvalMult(oi_plus+1, 2*lt-1, lt, v);

    };
  };
  
  ////////////////////// SURFACE ////////////////////////////
    template <int DIM>
  class HCurlDivSurfaceFiniteElement : public FiniteElement
  {
  public:
    using FiniteElement::FiniteElement;
    using FiniteElement::ndof;
    using FiniteElement::order;

    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM+1> & mip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedShape_Vector (const MappedIntegrationPoint<DIM,DIM+1> & mip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const = 0;

  };


  template <ELEMENT_TYPE ET> class HCurlDivSurfaceFE;

  
  template <ELEMENT_TYPE ET>
  class T_HCurlDivSurfaceFE : public HCurlDivSurfaceFiniteElement<ET_trait<ET>::DIM>,
    public VertexOrientedFE<ET>
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };
    enum { DIM_STRESS = (DIM+1)*(DIM+1) }; //check this!!!
    
    using VertexOrientedFE<ET>::vnums;
    using HCurlDivSurfaceFiniteElement<ET_trait<ET>::DIM>::ndof;
    using HCurlDivSurfaceFiniteElement<ET_trait<ET>::DIM>::order;

    INT<DIM> order_inner;


  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    
    T_HCurlDivSurfaceFE (int aorder)
    {
      order = aorder;
      order_inner = aorder;
    }
    
    virtual ELEMENT_TYPE ElementType() const { return ET; }
    const HCurlDivSurfaceFE<ET> * Cast() const { return static_cast<const HCurlDivSurfaceFE<ET>*> (this); } 
    
    INLINE void SetOrderInner (INT<DIM,int> order) { order_inner = order; }

    virtual void ComputeNDof()
    {
      cout << "Error, T_HCurlDivSurfaceFE<ET>:: ComputeNDof not available for base class" << endl;
    }

    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiffDiff<DIM+1>> adp;
      for ( int i=0; i<DIM+1; i++)
      {
        adp(i) = AutoDiffDiff<DIM+1>(ip(i),i);
      }

      Cast() -> T_CalcShape (TIP<DIM, AutoDiffDiff<DIM+1>> (adp), SBLambda([&] (int nr, auto val)
                                          {
                                            shape.Row(nr).AddSize(DIM_STRESS) = val.Shape();
                                          }));
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
      Cast() -> T_CalcShape (TIP<DIM, AutoDiffDiff<DIM+1>> (addp), SBLambda([&] (int nr, auto val)
                                          {
                                            shape.Row(nr).AddSize(DIM_STRESS) = val.Shape();
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
      Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM+1>> (addp),SBLambda([&](int nr,auto val)
      {
        Vec<DIM_STRESS> vecshape = val.Shape();
        BareVector<double> matshape = shape.Row(nr);
        VecToMat<DIM+1> (vecshape, matshape);
      }));
    }


  };

  template <> class HCurlDivSurfaceFE<ET_SEGM> : public T_HCurlDivSurfaceFE<ET_SEGM> 
  {
    
  public:
    using T_HCurlDivSurfaceFE<ET_SEGM> :: T_HCurlDivSurfaceFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      ndof += order_inner[0]+1;
      order = max2(order,order_inner[0]);

    }
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<1,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {
      auto x = ip.x;
      AutoDiffDiff<2> ddlami[2] ={ x, 1-x };
      
      int ii = 0;
      
      ArrayMem<AutoDiffDiff<2>,20> ha(order_inner[0]+2);
      
      int es = 0,ee = 1;
      if(vnums[es] > vnums[ee]) swap (es,ee);

      AutoDiffDiff<2> ls = ddlami[es],le = ddlami[ee];

      //IntegratedLegendreMonomialExt::Calc(order_inner[0]+2, le-ls, ha);
      ScaledLegendrePolynomial(order_inner[0],le-ls,1,ha);

      for(int l = 0; l <= order_inner[0]; l++)
	shape[ii++] =  T_Dl2xRotDl1_v(le, ls, ha[l]);      
    };
  };

}


#endif
  
