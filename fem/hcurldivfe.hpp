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

    virtual void CalcCurlShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<double> divshape) const = 0;
    
    virtual void CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
      BareSliceMatrix<double> shape) const = 0;    

    virtual void CalcMappedDivShape (const MappedIntegrationPoint<DIM,DIM> & mip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedCurlShape (const MappedIntegrationPoint<DIM,DIM> & mip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & ir,
                                         BareSliceMatrix<SIMD<double>> shapes) const = 0;

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                                  BareSliceVector<> coefs,
                                  BareSliceMatrix<SIMD<double>> values) const = 0;

    virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir,
                                  BareSliceMatrix<SIMD<double>> values,
                                  BareSliceVector<> coefs) const = 0;   
  };
  

  
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

    int order_facet[ET_trait<ET>::N_FACET];
    int order_inner;

    // additional curl-div free bubbles
    bool plus;
    bool withtrace;

  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    
    T_HCurlDivFE (int aorder, bool _plus = false, bool _withtrace = false)
      : plus(_plus), withtrace(_withtrace)
    {
      order = aorder;
      for (auto & of : order_facet) of = aorder;
      order_inner = aorder;
    }
    
    virtual ELEMENT_TYPE ElementType() const { return ET; }
    const HCurlDivFE<ET> * Cast() const { return static_cast<const HCurlDivFE<ET>*> (this); } 
    
    INLINE void SetOrderFacet (int nr, int order) { order_facet[nr] = order; }
    INLINE void SetOrderInner (int order) { order_inner = order; }

    virtual void ComputeNDof()
    {
      cout << "Error, T_HCurlDivFE<ET>:: ComputeNDof not available, only for ET == TRIG" << endl;
    }

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

    virtual void CalcCurlShape (const IntegrationPoint & ip,
                               BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiffDiff<DIM>> adp;
      for ( int i=0; i<DIM; i++)
      {
        adp[i] = AutoDiffDiff<DIM>(ip(i),i);
      }
      
      Cast() -> T_CalcShape (TIP<DIM, AutoDiffDiff<DIM>> (adp), SBLambda([&] (int nr, auto val)
                                          {
                                            shape.Row(nr).AddSize(DIM) = val.CurlShape();
                                          }));
    }
    


    virtual void CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
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
	shape.Row(nr).AddSize(DIM_STRESS) = val.Shape();
      }));
    }

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                                         BareSliceMatrix<SIMD<double>> shapes) const override
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
      for (size_t i = 0; i < mir.Size(); i++)
        {
          Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir[i];
          TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          
          this->Cast() -> T_CalcShape (addp,
                                       SBLambda ([i,shapes] (size_t j, auto val) 
                                                 {

                                                   Vec<DIM*DIM,SIMD<double>> vecshape = val.Shape();
                                                   for (size_t k = 0; k < DIM*DIM; k++)
                                                     shapes(j*(DIM*DIM)+k,i) = vecshape(k);
                                                 }));
        }
    }


    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & bmir,
                                  BareSliceVector<> coefs,
                                  BareSliceMatrix<SIMD<double>> values) const override
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
      for (size_t i = 0; i < bmir.Size(); i++)
        {
          double *pcoefs = &coefs(0);
          const size_t dist = coefs.Dist();
          
          Vec<DIM_STRESS,SIMD<double>> sum(0.0);
          Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir[i];
          TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          
          Cast() -> T_CalcShape (addp,
                                 SBLambda ([&sum,&pcoefs,dist] (size_t j, auto val)
                                           {
                                             sum += (*pcoefs)*val.Shape();
                                             pcoefs += dist;
                                           }));

	  for (size_t k = 0; k < DIM*DIM; k++)
	    values(k,i) = sum(k);
        }
    }

    virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & bmir,
                                  BareSliceMatrix<SIMD<double>> values,
                                  BareSliceVector<> coefs) const override
    {
       for (size_t i = 0; i < bmir.Size(); i++)
        {
          Mat<DIM,DIM,SIMD<double>> mat;
	  
	  auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);

	  for (size_t k = 0; k < DIM*DIM; k++)
	    mat(k) = values(k,i);
	  
	  Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir[i];
          TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          double *pcoefs = &coefs(0);
          const size_t dist = coefs.Dist();

          Cast() -> T_CalcShape (addp,
                                 SBLambda ([mat,&pcoefs,dist] (size_t j, auto val)
                                           {                                          
					     Vec<DIM*DIM,SIMD<double>> vecshape = val.Shape();
                                             
                                             SIMD<double> sum = 0.0;
                                             for (size_t k = 0; k < DIM*DIM; k++)
                                               sum += mat(k) * vecshape(k);
                                             
                                             *pcoefs += HSum(sum);
                                             pcoefs += dist;
                                           }));
        }
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
        Mat<DIM> jac = mip.GetJacobian();
        Mat<DIM> inv_jac = mip.GetJacobianInverse();        
	Mat<DIM> hesse_FinvT[3], F_HFinvT_Finv[3];
		
	double eps = 1e-4;
		
	Mat<DIM> jacrinv, jaclinv,jacrrinv, jacllinv;
	for (int dir = 0; dir < DIM; dir++)
	  {
	    IntegrationPoint ipr = mip.IP();
	    IntegrationPoint ipl = mip.IP();
	    IntegrationPoint iprr = mip.IP();
	    IntegrationPoint ipll = mip.IP();
	    
	    ipr(dir) += eps;
	    ipl(dir) -= eps;
	    iprr(dir) += 2*eps;
	    ipll(dir) -= 2*eps;
	    
	    MappedIntegrationPoint<DIM,DIM> mipr(ipr, mip.GetTransformation());
	    MappedIntegrationPoint<DIM,DIM> mipl(ipl, mip.GetTransformation());
	    MappedIntegrationPoint<DIM,DIM> miprr(iprr, mip.GetTransformation());
	    MappedIntegrationPoint<DIM,DIM> mipll(ipll, mip.GetTransformation());
	    
	    jacrinv = Trans(mipr.GetJacobianInverse());    
	    jaclinv = Trans(mipl.GetJacobianInverse());

	    jacrrinv = Trans(miprr.GetJacobianInverse());    
	    jacllinv = Trans(mipll.GetJacobianInverse());
	    
	    for (int j = 0; j < DIM; j++)
	      {
		hesse_FinvT[0](j,dir) = (8.0*jacrinv(0,j) - 8.0*jaclinv(0,j) - jacrrinv(0,j) + jacllinv(0,j) ) / (12.0*eps);
		hesse_FinvT[1](j,dir) = (8.0*jacrinv(1,j) - 8.0*jaclinv(1,j) - jacrrinv(1,j) + jacllinv(1,j) ) / (12.0*eps);
		hesse_FinvT[2](j,dir) = (8.0*jacrinv(2,j) - 8.0*jaclinv(2,j) - jacrrinv(2,j) + jacllinv(2,j) ) / (12.0*eps);

	      }
	  }
	
	for(int i=0; i<DIM; i++)
	  F_HFinvT_Finv[i] = jac * hesse_FinvT[i] * inv_jac;

        Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM>> (addp),SBLambda([&](int nr,auto val)
                                  {
                                    shape.Row(nr).AddSize(DIM) = val.DivShape();
                                    BareVector<double> divshape = shape.Row(nr);				    
                                    Vec<DIM*DIM> matshape = val.Shape();				    
				    
                                    for(int k=0; k<DIM; k++)
                                    {
                                      for(int j=0; j<DIM*DIM; j++)
                                      {
					divshape(k) += F_HFinvT_Finv[k](j) * matshape(j);
                                      }
                                    }
                                    
                                  }));
	
      }
    }

     virtual void CalcMappedCurlShape (const MappedIntegrationPoint<DIM,DIM> & mip,
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
          shape.Row(nr).AddSize(DIM) = val.CurlShape();
        }));	
      }
      else
	{	
        Mat<DIM> jac = mip.GetJacobian();
        Mat<DIM> inv_jac = mip.GetJacobianInverse();        
	Mat<DIM> curl_FT[2], F_curlFT_Finv[2];
	Vec<DIM> curl_Jinv;
		
	double eps = 1e-6;
	
	Mat<DIM> jacr, jacl;
	for (int dir = 0; dir < DIM; dir++)
	  {
	    IntegrationPoint ipr = mip.IP();
	    IntegrationPoint ipl = mip.IP();
    
	    ipr(dir) += eps;
	    ipl(dir) -= eps;    	    

	    mip.GetTransformation().CalcJacobian(ipr, jacr);
	    mip.GetTransformation().CalcJacobian(ipl, jacl);
	    
	    jacr = Trans(jacr);
	    jacl = Trans(jacl);

	    for (int j = 0; j < DIM; j++)
	     {	       
	       curl_FT[0](DIM-1-dir,j) = pow(-1.0,dir) * (jacr(j,0) - jacl(j,0)) / (2.0*eps);
	       curl_FT[1](DIM-1-dir,j) = pow(-1.0,dir)  * (jacr(j,1) - jacl(j,1)) / (2.0*eps);	      
	     }
	  }

	F_curlFT_Finv[0] = jac * curl_FT[0] * inv_jac;
	F_curlFT_Finv[1] = jac * curl_FT[1] * inv_jac;
	
	Mat<DIM> hesse[3];
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
	curl_Jinv(0) = -iad_det.DValue(1);
	curl_Jinv(1) = iad_det.DValue(0);
	
	Vec<DIM> curl_Jinv_FT;
	curl_Jinv_FT(0) = curl_Jinv(0) * Trans(jac)(0,0) + curl_Jinv(1) * Trans(jac)(1,0);
	curl_Jinv_FT(1) = curl_Jinv(0) * Trans(jac)(0,1) + curl_Jinv(1) * Trans(jac)(1,1);

        Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM>> (addp),SBLambda([&](int nr,auto val)
                                  {
                                    shape.Row(nr).AddSize(DIM) = val.CurlShape();
                                    BareVector<double> curlshape = shape.Row(nr);				    
                                    Vec<DIM*DIM> matshape = val.Shape();				    				    
                                    for(int k=0; k<DIM; k++)
                                    {
                                      for(int j=0; j<DIM*DIM; j++)
				      {
					curlshape(k) += 1.0/mip.GetJacobiDet() * F_curlFT_Finv[k](j) * matshape(j);
                                      }
				      for(int j=0; j<DIM; j++)
				      {
					curlshape(k) += curl_Jinv_FT(j) * matshape(k+j*DIM);
				      }
                                    }
                                    
                                  }));
	
      }
    }

  };

  /* ############### edge basis functions - div-free ############### */
  /* sigma(grad v) = Curl(grad v), where Curl is the 1D to 2D curl operator */
  template <int D, typename T>  class T_Sigma_gradv;
  template <typename T>  class T_Sigma_gradv<2,T>
  {
    AutoDiffDiff<2,T> v;
  public:
    T_Sigma_gradv  (AutoDiffDiff<2,T> av) : v(av){ ; }
    
    Vec<4,T> Shape() {
      return Vec<4,T> (-v.DDValue(0,1), v.DDValue(0,0),
		     -v.DDValue(1,1),v.DDValue(0,1)
		     );
    }

    Vec<2,T> DivShape()
    {      
      return Vec<2,T> (0.0,0.0);     
    }

    Vec<2,T> CurlShape()
    {      
      return Vec<2,T> (0.0,0.0);     
    }
  };

  template <int D, typename T>
  auto Sigma_gradv (AutoDiffDiff<D,T> av) { return T_Sigma_gradv<D,T>(av); }

  /* ############### div-free basis function WITH trace ############### */
  /* For basis functions including the trace */
  
  template <int D, typename T> class T_Sigma_gradu_v;
  template <typename T> class T_Sigma_gradu_v<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_Sigma_gradu_v  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av){ ; }

    Vec<4,T> Shape() {      
      return Vec<4,T> (-u.DDValue(1,0) * v.Value()  -  v.DValue(1)*u.DValue(0),
		     u.DDValue(0,0) * v.Value() + v.DValue(0)*u.DValue(0),
		     -u.DDValue(1,1) * v.Value() - v.DValue(1)*u.DValue(1),
		     u.DDValue(0,1) * v.Value() +  v.DValue(0)*u.DValue(1)
		     );
    }

    Vec<2,T> DivShape()
    {     
      return Vec<2,T> (0.0,0.0);
    }

    Vec<2,T> CurlShape()
    {
      T vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T uxx = u.DDValue(0,0), uxy = u.DDValue(0,1), uyy = u.DDValue(1,1);
      T vx = v.DValue(0), vy = v.DValue(1);

      return Vec<2,T> ( (vy*uxy - vx*uyy) +  (ux * vyy - uy * vxy), (-vy*uxx + vx*uxy) + (-ux*vxy + uy*vxx));     
    }    
  };

  template <int D, typename T>
  auto Sigma_gradu_v (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av ) { return T_Sigma_gradu_v<D,T>(au,av); }

  /* ############### Type 1.1 (for QUAD) - inner basis functions - not div-free ############### */
  /* grad(u) * Curl(v) */
  template <int D, typename T> class T_Gradu_Curlv;
  template <typename T>  class T_Gradu_Curlv<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_Gradu_Curlv  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av){ ; }

    Vec<4,T> Shape() {
      return Vec<4,T> (-  v.DValue(1)*u.DValue(0),
		      v.DValue(0)*u.DValue(0),
		     - v.DValue(1)*u.DValue(1),
		     v.DValue(0)*u.DValue(1)
		     );
    }

    Vec<2,T> DivShape()
    {
      T vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);
      T uxx = u.DDValue(0,0), uxy = u.DDValue(0,1), uyy = u.DDValue(1,1);      
      T ux = u.DValue(0), uy = u.DValue(1);
      T vx = v.DValue(0), vy = v.DValue(1);
      
      return Vec<2,T> (-vy*uxx + vx*uxy,-vy*uxy + vx * uyy);
    }

    Vec<2,T> CurlShape()
    {
      T vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T uxx = u.DDValue(0,0), uxy = u.DDValue(0,1), uyy = u.DDValue(1,1);
      T vx = v.DValue(0), vy = v.DValue(1);

      return Vec<2,T> (vyy*ux - vxy*uy,-vxy*ux+vxx*uy);     
    }    
  };

  template <int D, typename T>
  auto Gradu_Curlv (AutoDiffDiff<D,T> au,AutoDiffDiff<D,T> av) { return T_Gradu_Curlv<D, T>(au,av); }

  
  /* ############### Type 2 (QUAD) - inner basis functions - div-free ############### */
  /* u * sigma(grad v) = Curl(grad v), where Curl is the 1D to 2D curl operator */
  template <int D, typename T> class T_u_Sigma_gradv;
  template <typename T> class T_u_Sigma_gradv<2,T>
  {
    AutoDiffDiff<2,T> v,u;
  public:
    T_u_Sigma_gradv  (AutoDiffDiff<2,T> au,AutoDiffDiff<2,T> av) : u(au),v(av){ ; }
    
    Vec<4,T> Shape() {
      return Vec<4,T> (-u.Value() * v.DDValue(0,1), u.Value() *  v.DDValue(0,0),
		     -u.Value() * v.DDValue(1,1), u.Value() * v.DDValue(0,1)
		     );
    }

    Vec<2,T> DivShape()
    {      
      return Vec<2,T> (-u.DValue(0) *v.DDValue(0,1) + u.DValue(1) * v.DDValue(0,0) , -u.DValue(0) *v.DDValue(1,1) + u.DValue(1) * v.DDValue(0,1));     
    }

    Vec<2,T> CurlShape()
    {      
      return Vec<2,T> (u.DValue(1) *v.DDValue(0,1) - u.DValue(0) * v.DDValue(1,1) ,-u.DValue(1) *v.DDValue(0,0) + u.DValue(0) * v.DDValue(0,1) );     
    }
  };

  template <int D, typename T>
  auto u_Sigma_gradv (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_u_Sigma_gradv<D,T>(au,av); }


  /* ############### Type 2 - inner basis functions - NOT div-free ############### */
  /* Curl(grad(u)) * v - grad(u) * Curl(v) */
  template <int D, typename T>  class T_Curlgraduv_graducurlv;
  template <typename T>  class T_Curlgraduv_graducurlv<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_Curlgraduv_graducurlv  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av){ ; }

    Vec<4,T> Shape() {

      auto trace = (v.DValue(1)*u.DValue(0) -  v.DValue(0)*u.DValue(1)  )/2.0;
      
      return Vec<4,T> (-u.DDValue(1,0) * v.Value()  +  v.DValue(1)*u.DValue(0) - trace,
		     u.DDValue(0,0) * v.Value() - v.DValue(0)*u.DValue(0),
		     -u.DDValue(1,1) * v.Value() + v.DValue(1)*u.DValue(1),
		     u.DDValue(0,1) * v.Value() -  v.DValue(0)*u.DValue(1) - trace
		     );
    }

    Vec<2,T> DivShape()
    {
      T vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T uxx = u.DDValue(0,0), uxy = u.DDValue(0,1), uyy = u.DDValue(1,1);
      T vx = v.DValue(0), vy = v.DValue(1);
            
      return 0.5 * Vec<2,T> (3 * uxx * vy - 3 * uxy * vx - vxy*ux + vxx*uy, 3* uxy * vy - 3 * uyy * vx - vyy*ux + vxy*uy);

    }

    Vec<2,T> CurlShape()
    {
      T vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T uxx = u.DDValue(0,0), uxy = u.DDValue(0,1), uyy = u.DDValue(1,1);
      T vx = v.DValue(0), vy = v.DValue(1);
      
      return Vec<2,T> ( (vy*uxy - vx*uyy) - (ux * vyy - uy * vxy), (-vy*uxx + vx*uxy) - (-ux*vxy + uy*vxx));     
    }
    
  };

  template <int D, typename T>
  auto Curlgraduv_graducurlv (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_Curlgraduv_graducurlv<D,T>(au, av); }

  /* ############### Type 3 - inner basis functions - div-free ############### */
  /*  Curl( [grad(l1) l2 - l1 grad(l2)] * v ) */
  template <int D, typename T> class T_type4;
  template <typename T> class T_type4<2,T>
  {
    AutoDiffDiff<2,T> l1,l2,v;
  public:
    T_type4  (AutoDiffDiff<2,T> lam1, AutoDiffDiff<2,T> lam2, AutoDiffDiff<2,T> av) : l1(lam1), l2(lam2), v(av){ ; }

    Vec<4,T> Shape() {
      T lam1x = l1.DValue(0), lam1y = l1.DValue(1);
      T lam2x = l2.DValue(0), lam2y = l2.DValue(1);
      T vx = v.DValue(0), vy = v.DValue(1);                      

      //auto trace = ( (v.Value() * ( - lam1x * lam2y + lam2x * lam1y) - (lam1x*l2.Value() - lam2x*l1.Value()) * vy)
      //		     + ( v.Value() * ( lam1y * lam2x - lam2y * lam1x) + (lam1y*l2.Value() - lam2y*l1.Value()) * vx))/2.0;
      
      return Vec<4,T> (v.Value() * ( - lam1x * lam2y + lam2x * lam1y) - (lam1x*l2.Value() - lam2x*l1.Value()) * vy,
		      (lam1x*l2.Value() - lam2x*l1.Value()) * vx,
		     -(lam1y*l2.Value() - lam2y*l1.Value()) * vy,
		     v.Value() * ( lam1y * lam2x - lam2y * lam1x) + (lam1y*l2.Value() - lam2y*l1.Value()) * vx
		     ); 
    }

    Vec<2,T> DivShape()
    {     
      return Vec<2,T> (0.0,0.0);     
    }

    Vec<2,T> CurlShape()
    {
      T lam1x = l1.DValue(0), lam1y = l1.DValue(1); 
      T lam2x = l2.DValue(0), lam2y = l2.DValue(1); 
      T vx = v.DValue(0), vy = v.DValue(1), vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);

      return Vec<2,T> ( vyy*(lam1x*l2.Value() - lam2x*l1.Value()) - vxy * (lam1y*l2.Value() - lam2y*l1.Value()) - 3* vy*(-lam1x*lam2y+lam2x*lam1y),
      		      vxx*(lam1y*l2.Value() - lam2y*l1.Value()) - vxy * (lam1x*l2.Value() - lam2x*l1.Value()) + 3* vx*(-lam1x*lam2y+lam2x*lam1y)
      				);
    }    
  };

  template <int D, typename T>
  auto type4 (AutoDiffDiff<D,T> al1, AutoDiffDiff<D,T> al2, AutoDiffDiff<D,T> av) { return T_type4<D,T>(al1, al2, av); }
  
  /* ############### Special functions for curld-div bubbles ############### */
  /* is equal to type 1 - tr(type 1)  */
  /* this produces a divergence, but curldiv is equal to zero */
  template <int D, typename T> class T_curldivfreebubbles;
  template <typename T> class T_curldivfreebubbles<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_curldivfreebubbles  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av){ ; }

    Vec<4,T> Shape() {
      return Vec<4,T> (-u.DDValue(0,1) * v.Value() -  v.DValue(0)*u.DValue(1),
		     u.DDValue(0,0) * v.Value() + v.DValue(0)*u.DValue(0),
		     -u.DDValue(1,1) * v.Value() - v.DValue(1)*u.DValue(1),
		     u.DDValue(1,0) * v.Value() + v.DValue(1)*u.DValue(0)
		     );
    }

    Vec<2,T> DivShape()
    {
      T uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T vxx = v.DDValue(0,0), vyy = v.DDValue(1,1), vxy = v.DDValue(0,1);
      T vx = v.DValue(0), vy = v.DValue(1);
      
      return Vec<2,T> (-vxx*uy - vx*uxy + uxx*vy + ux*vxy,
		     -uyy*vx - vxy*uy + uxy*vy + vyy*ux);

    }

    Vec<2,T> CurlShape()
    {     
      throw Exception("not implemented for curldivfreebubbles");
    }
    
  };

  template <int D, typename T>
  auto curldivfreebubbles (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_curldivfreebubbles<D, T>(au, av); }

  /* Edge basis functions which are normal-tangential continuous */
  /* calculates [(grad l1) o-times (rot-grad l2) ] * scaledlegendre */
  /* DivShape assumes that phi_12 = [(grad l1) o-times (rot-grad l2) ] is constant!!! */
  /* This is the old implementation! */
  //
  //class T_Dl2xRotDl1_v
  //{
  //  AutoDiffDiff<2> l1,l2,v;
  //public:
  //  T_Dl2xRotDl1_v  (AutoDiffDiff<2> lam1, AutoDiffDiff<2> lam2, AutoDiffDiff<2> av) : l1(lam1), l2(lam2), v(av) { ; }
  //
  //  Vec<4> Shape() {
  //    return Vec<4> (-v.Value()*(l1.DValue(0)*l2.DValue(1)),
  //		     v.Value()*(l1.DValue(0)*l2.DValue(0)),
  //		     -v.Value()*(l1.DValue(1)*l2.DValue(1)),
  //		     v.Value()*(l1.DValue(1)*l2.DValue(0))
  //		     );
  //  }
  //
  //  Vec<2> DivShape()
  //  {
  //    double vx = v.DValue(0), vy = v.DValue(1);
  //
  //    double lam1x = l1.DValue(0);
  //    double lam1y = l1.DValue(1);      
  //    
  //    double lam2x = l2.DValue(0);
  //    double lam2y = l2.DValue(1);
  //    
  //    return Vec<2> (  - vx *lam1x*lam2y + vy*lam1x*lam2x, -vx *lam1y*lam2y + vy*lam1y*lam2x) ;      
  //  }
  //
  //  Vec<2> CurlShape()
  //  {     
  //    throw Exception("not implemented for T_Dl2xRotDl1_v");
  //  }
  //
  //};

  
  /* Face basis functions which are normal-tangential continuous */
  /* calculates [(grad l1) o-times (grad l2 x grad l3)] * legendre */
  /* DivShape assumes that phi_12 =  [(grad l1) o-times (grad l2 x grad l3)] is constant!!! */
  template <typename T>
  class T_Dl1_o_Dl2xDl3_v
  {
    AutoDiff<3,T> l1,l2,l3,v;
  public:
    T_Dl1_o_Dl2xDl3_v  (AutoDiff<3,T> lam1, AutoDiff<3,T> lam2, AutoDiff<3,T> lam3, AutoDiff<3,T> av) : l1(lam1), l2(lam2), l3(lam3), v(av) { ; }
    
    Vec<9,T> Shape() {

      T cross1 = l2.DValue(1)*l3.DValue(2) - l2.DValue(2)*l3.DValue(1);
      T cross2 = -(l2.DValue(0)*l3.DValue(2) - l2.DValue(2)*l3.DValue(0));
      T cross3 = l2.DValue(0)*l3.DValue(1) - l2.DValue(1)*l3.DValue(0);

      Vec<9,T> sigmaref;
      
      for (int i=0; i<3; i++)
	{
	  sigmaref(i*3)= v.Value() * l1.DValue(i) * cross1;
	  sigmaref(i*3+1)= v.Value() * l1.DValue(i) * cross2;
	  sigmaref(i*3+2)= v.Value() * l1.DValue(i) * cross3;
	}
      return sigmaref;  

    }

    Vec<3,T> DivShape()
    {
      T vx = v.DValue(0), vy = v.DValue(1), vz = v.DValue(2);

      T cross1 = l2.DValue(1)*l3.DValue(2) - l2.DValue(2)*l3.DValue(1);
      T cross2 = -(l2.DValue(0)*l3.DValue(2) - l2.DValue(2)*l3.DValue(0));
      T cross3 = l2.DValue(0)*l3.DValue(1) - l2.DValue(1)*l3.DValue(0);

      return Vec<3,T> (vx * l1.DValue(0) * cross1 + vy * l1.DValue(0) * cross2 + vz * l1.DValue(0) * cross3,
		     vx * l1.DValue(1) * cross1 + vy * l1.DValue(1) * cross2 + vz * l1.DValue(1) * cross3,
		     vx * l1.DValue(2) * cross1 + vy * l1.DValue(2) * cross2 + vz * l1.DValue(2) * cross3
		     );
    }

    Vec<3,T> CurlShape()
    {     
      throw Exception("not implemented for T_Dl1_o_Dl2xDl3_v");
    }

  };

  
  /* Identity = Inner bubble function (normatl-tangential component is zero) */
  /* calculates I * legendre */
  template <int D, typename T>  class T_Id_v;
  template <typename T>  class T_Id_v<3,T>
  {
     AutoDiff<3,T> v;
  public:
    T_Id_v  (AutoDiff<3,T> av) : v(av) { ; }
    
    Vec<9,T> Shape() {
      T zero = 0.0;
      Vec<9,T> Id_v= zero;

      for (int i=0; i<3; i++)
	Id_v(i*(4))= v.Value();      
      return Id_v;
    }

    Vec<3,T> DivShape()
    {
      Vec<3,T> div_Id_v=0;
      for (int i=0; i<3; i++)
	div_Id_v(i) = v.DValue(i);
      return div_Id_v;
    }

    Vec<3,T> CurlShape()
    {     
      Vec<3, T> curl_Id_v=0;
      curl_Id_v(0) = v.DValue(2) -v.DValue(1);
      curl_Id_v(1) = v.DValue(0) -v.DValue(2);
      curl_Id_v(2) = v.DValue(1) -v.DValue(0); 
          
      return curl_Id_v;
    }

  };

  template <int D, typename T>
  auto Id_v (AutoDiff<D,T> av) { return T_Id_v<D,T>(av); }
  
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
        ndof += order_facet[i]+1;
        order = max2(order, order_facet[i]);
      }      
      int ninner = 3 * ((order_inner +1) * (order_inner))/2; 
      order = max2(order, order_inner);
      if (plus)
      {
	//throw Exception(" please update this first - ComputeNdof in HCurlDiv<ET_TRIG>");
        order ++;	
        ninner += 2 *(order_inner+1); 
      }
      ndof += ninner;
      if (withtrace)
	ndof += (order_inner +1) * (order_inner+2)/2.0;
    }
    
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {
      auto x = ip.x, y = ip.y;
      Tx ddlami[3] ={ x, y, 1-x-y };
      
      int ii = 0;
      
      int oi=order_inner;      
      int maxorder_facet =
        max2(order_facet[0],max2(order_facet[1],order_facet[2]));      

      ArrayMem<Tx,20> ha(maxorder_facet);
      ArrayMem<Tx,20> v(oi), u(oi);

      for (int i = 0; i < 3; i++)
        {
	  INT<2> e = ET_trait<ET_TRIG>::GetEdgeSort (i, vnums);	  	  
          Tx ls = ddlami[e[0]], le = ddlami[e[1]];
	  	 
	  IntLegNoBubble::EvalMult (maxorder_facet, le-ls, 0.25*le*ls, ha);
	  
          for (int l = 0; l <= order_facet[i]; l++)	    
	    shape[ii++] = Sigma_gradv(ha[l]);	 
        }
      
      Tx ls = ddlami[0];
      Tx le = ddlami[1];
      Tx lt = ddlami[2];
                  
      IntLegNoBubble::EvalMult (oi, le-lt, 0.25*le*lt, u);
      LegendrePolynomial::EvalMult(oi, 2*ls-1, ls, v);
      
      for(int i = 0; i <= oi-1; i++)
      {
        for(int j = 0; j+i <= oi-1; j++)
        {	  
	  shape[ii++] = Curlgraduv_graducurlv(u[i],v[j]);	  	  
        }	
      }
     
      IntLegNoBubble::EvalMult (oi, le-ls, 0.25*le*ls, u);
      LegendrePolynomial::EvalMult(oi, 2*lt-1, lt, v);
      
      for(int i = 0; i <= oi-1; i++)
      {
        for(int j = 0; j+i <= oi-1; j++)
        {
	  shape[ii++] = Sigma_gradv(u[i]*v[j]); //divfree!
	  shape[ii++] = Curlgraduv_graducurlv(u[i],v[j]); 	 
        }	
      }
      
      if (withtrace)
	{
	  LegendrePolynomial::Eval(oi, 2*lt-1, v);
	  for (int i = 0; i <= oi; i++)
	    shape[ii++] = type4(le, ls, v[i]);

	  IntLegNoBubble::EvalMult (oi, le-lt, 0.25*le*lt, u);
	  LegendrePolynomial::EvalMult(oi, 2*ls-1, ls, v);
	  for(int i = 0; i <= oi-1; i++)
	      for(int j = 0; j+i <= oi-1; j++)       
		    shape[ii++] = Sigma_gradu_v(u[i],v[j]);	  		  	    
	}
    };
  };


    template <> class HCurlDivFE<ET_QUAD> : public T_HCurlDivFE<ET_QUAD> 
  {
    
  public:
    using T_HCurlDivFE<ET_QUAD> :: T_HCurlDivFE;

    virtual void ComputeNDof()
    {     
      order = 0;
      ndof = 0;
      for (int i=0; i<4; i++)
	{
	  ndof += order_facet[i]+1;
	  order = max2(order, order_facet[i]);
	}
      
      int ninner = 2*(order_inner+1)*(order_inner+1) + (order_inner+2)*(order_inner) *2;
      
      order = max2(order, order_inner);
      order += 4;
      ndof += ninner;     
     
    }
    
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {
      auto x = ip.x, y = ip.y;
      Tx lx[4] ={1-x, x, x, 1-x};
      Tx ly[4] = {1-y, 1-y, y, y};
      //Tx lam[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};
      Tx edgebubbles[4] = {(1-x)*x, x*(1-x), y*(1-y), (1-y)*y}; 
      
      int ii = 0;
      
      int oi=order_inner;
      
      int maxorder_facet =
        max2(order_facet[3],max2(order_facet[0],max2(order_facet[1],order_facet[2])));

      const EDGE * edges = ElementTopology::GetEdges(ET_QUAD);

      ArrayMem<Tx,20> ha(maxorder_facet);
      ArrayMem<Tx,20> v(oi), u(oi);
      for (int i = 0; i < 4; i++)
        {
	  INT<2> e = ET_trait<ET_QUAD>::GetEdgeSort (i, vnums);
	  
          Tx xi = lx[e[1]]+ly[e[1]]-lx[e[0]]-ly[e[0]];
          Tx eta = lx[e[0]]*ly[e[0]]+lx[e[1]]*ly[e[1]];
	  	 
	  IntLegNoBubble::EvalMult (maxorder_facet , xi, 0.25*edgebubbles[i], ha);	  
	  
          for (int l = 0; l <= order_facet[i]; l++)	    
	    shape[ii++] = Sigma_gradv(eta*ha[l]);	 
        }
           
      IntLegNoBubble::EvalMult (oi+2, lx[0]-lx[1], 0.25*lx[0]*lx[1], u);
      IntLegNoBubble::EvalMult (oi+2, ly[0]-ly[2], 0.25*ly[0]*ly[2], v);

      // constants in diagonal
      // off-diagonal constant functions are provided by edge functions
      shape[ii++] = Sigma_gradu_v(lx[0],ly[0]);      
      shape[ii++] = Sigma_gradu_v(ly[0],lx[0]);

      //provide mixed functions in the diagonal
      for(int i = 0; i <= oi-1; i++)
      {
        for(int j = 0; j <= oi-1; j++)
        {
          shape[ii++] = Sigma_gradv(u[i]*v[j]);
	  shape[ii++] = Sigma_gradu_v(u[i],v[j]);
	}
      }

      //are needed to compensate the terms in the off diagonal from the block before
      for(int i = 0; i <= oi+1; i++)
      {
        for(int j = 0; j <= oi-1; j++)
        {
          shape[ii++] = u_Sigma_gradv(u[j],v[i]);
          shape[ii++] = u_Sigma_gradv(v[j],u[i]);	  	  
        }
      }
       
      // lienar (and high order) parts in the diagonal
      for(int i = 0; i <= oi-1; i++)
       {
	 shape[ii++] = Sigma_gradu_v(ly[0],u[i]);
	 shape[ii++] = Sigma_gradu_v(lx[0],v[i]);

	 shape[ii++] = Gradu_Curlv(u[i],ly[0]);
	 shape[ii++] = Gradu_Curlv(v[i],lx[0]);	 
       }
      
    };
  };
  
 
     template <> class HCurlDivFE<ET_TET> : public T_HCurlDivFE<ET_TET> 
  {
    
  public:
    using T_HCurlDivFE<ET_TET> :: T_HCurlDivFE;

    virtual void ComputeNDof()
    {     
      order = 0;
      ndof = 0;
      for (int i=0; i<4; i++)
      {
        ndof += (order_facet[i]+1)*(order_facet[i]+2);
        order = max2(order, order_facet[i]);
      }
      //first type + second type
      int ninner = (order_inner + 1)*(order_inner + 2)*(order_inner + 3)/6.0 + 8.0/6.0* ((order_inner +2) * (order_inner +1) * (order_inner));
      
      order = max2(order, order_inner);
      //if (plus)
      //{	
      //  order ++;
      //  ninner += 2*order_inner; 
      //}
      ndof += ninner;     
    }
    
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {
      auto x = ip.x, y = ip.y, z = ip.z ;
      //AutoDiffDiff<3> ddlami[4] ={ x, y, z, 1-x-y-z };
      //Tx ddlami[4] ={ x, y, z, 1-x-y-z };

      typedef decltype(x.Value()+x.Value()) T;      
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> zz(z.Value(), &z.DValue(0));
      AutoDiff<3,T> ddlami[4] = {xx, yy, zz, 1-xx-yy-zz};
      
      int ii = 0;
      
      int maxorder_facet =
        max2(order_facet[0],max2(order_facet[1],order_facet[2]));


       const FACE * faces = ElementTopology::GetFaces(ET_TET);

       
       ArrayMem<AutoDiff<3,T>,20> ha((maxorder_facet+1)*(maxorder_facet+2)/2.0); 
       ArrayMem<AutoDiff<3,T>,20> v((order_inner + 1)*(order_inner + 2)*(order_inner + 3)/6.0);
       ArrayMem<AutoDiff<3,T>,20> dubb(order_inner*(order_inner+1)*(order_inner+2)/6.0);
       
      /* Edge based basis functions for tangential-normal continuity */
      for(int fa = 0; fa < 4; fa++)
        {
	  int fav[3] = {faces[fa][0], faces[fa][1], faces[fa][2]};
	  
	  int p = order_facet[fa];
	  //Sort vertices  first edge op minimal vertex
	  if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0], fav[1]);
	  if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1], fav[2]);
	  if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0], fav[1]);


          AutoDiff<3,T> ls = ddlami[fav[0]], le = ddlami[fav[1]], lt = ddlami[fav[2]];
	  
	  DubinerBasis3::Eval (maxorder_facet, ls, le, ha);

          for (int l = 0; l < (order_facet[fa]+1)*(order_facet[fa]+2)/2.0; l++)
	    {	      
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(le,ls,lt,ha[l]);	      
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(ls,lt,le,ha[l]);
	    }
        }

      int oi=order_inner;
      
      int es = 0; int ee = 1; int et = 2; int eo = 3 ;
      AutoDiff<3,T> ls = ddlami[es]; 
      AutoDiff<3,T> le = ddlami[ee];
      AutoDiff<3,T> lt = ddlami[et];
      AutoDiff<3,T> lo = ddlami[eo];

      LegendrePolynomial leg;
      JacobiPolynomialAlpha jac1(1);

      //############ type 1 ############
      
      leg.EvalScaled1Assign 
	(oi, lt-lo, lt+lo,
	 SBLambda ([&](size_t k, AutoDiff<3,T> polz) LAMBDA_INLINE
		   {
		     JacobiPolynomialAlpha jac2(2*k+2);
 
		     jac1.EvalScaledMult1Assign
		       (oi-k, le-lt-lo, 1-ls, polz, 
			SBLambda ([&] (size_t j, AutoDiff<3,T> polsy) LAMBDA_INLINE
				  {				    
				    jac2.EvalMult(oi - k - j, 2 * ls - 1, polsy, 
						  SBLambda([&](size_t j, AutoDiff<3,T> val) LAMBDA_INLINE
							   {
							     shape[ii++] =  Id_v(val);	      							     					     							   }));
				    jac2.IncAlpha2();
				  }));
		     jac1.IncAlpha2();
		   }));

      
      //############ type 2 ############
            
      leg.EvalScaled1Assign 
	(oi-1, lt-lo, lt+lo,
	 SBLambda ([&](size_t k, AutoDiff<3,T> polz) LAMBDA_INLINE
		   {
		     JacobiPolynomialAlpha jac2(2*k+2);
 
		     jac1.EvalScaledMult1Assign
		       (oi-1-k, le-lt-lo, 1-ls, polz, 
			SBLambda ([&] (size_t j, AutoDiff<3,T> polsy) LAMBDA_INLINE
				  {				    
				    jac2.EvalMult(oi-1 - k - j, 2 * ls - 1, polsy, 
						  SBLambda([&](size_t j, AutoDiff<3,T> val) LAMBDA_INLINE
							   {
							     shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(le,ls,lt,lo*val);	      
							     shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(ls,lt,le,lo*val);	  
							     shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(le,ls,lo,lt*val);	      
							     shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(ls,lo,le,lt*val);	  
							     shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(le,lo,lt,ls*val);	      
							     shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(lo,lt,le,ls*val);	  
							     shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(lo,ls,lt,le*val);	      
							     shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(lt,ls,lo,le*val);							     
							   }));
				    jac2.IncAlpha2();
				  }));
		     jac1.IncAlpha2();
		   }));
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

    virtual void CalcMappedShape (const MappedIntegrationPoint<DIM,DIM+1> & mip,
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

    int order_inner;


  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    
    T_HCurlDivSurfaceFE (int aorder)
    {
      order = aorder;
      order_inner = aorder;
    }
    
    virtual ELEMENT_TYPE ElementType() const { return ET; }
    const HCurlDivSurfaceFE<ET> * Cast() const { return static_cast<const HCurlDivSurfaceFE<ET>*> (this); } 
    
    INLINE void SetOrderInner (int order) { order_inner = order; }

    virtual void ComputeNDof()
    {
      cout << "Error, T_HCurlDivSurfaceFE<ET>:: ComputeNDof not available for base class" << endl;
    }

    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const
    {
      Vec<DIM, AutoDiffDiff<DIM+1>> adp;
      for ( int i=0; i<DIM; i++)
      {
        adp(i) = AutoDiffDiff<DIM+1>(ip(i),i);
      }

      Cast() -> T_CalcShape (TIP<DIM, AutoDiffDiff<DIM+1>> (adp), SBLambda([&] (int nr, auto val)
                                          {
                                            //shape.Row(nr).AddSize(DIM_STRESS) = val.Shape();
					    shape.Row(nr).AddSize(DIM_STRESS) = val;
                                          }));
    }

    
    virtual void CalcMappedShape (const MappedIntegrationPoint<DIM,DIM+1> & mip,
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
	if (DIM==1)
	  shape.Row(nr).AddSize(DIM_STRESS) = val;
	else
	  shape.Row(nr).AddSize(1) = val;
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
      ndof += order_inner+1;
      order = max2(order,order_inner);

    }
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<1,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {      
      auto x = ip.x;
      AutoDiffDiff<2> ddlami[2] ={ x, 1-x };
      
      int ii = 0;
      
      ArrayMem<AutoDiffDiff<2>,20> ha(order_inner+1);
      
      int es = 0,ee = 1;
      if(vnums[es] > vnums[ee]) swap (es,ee);

      AutoDiffDiff<2> ls = ddlami[es],le = ddlami[ee];
      
      IntLegNoBubble::EvalMult (order_inner, le-ls, 0.25*le*ls, ha);
      
      for(int l = 0; l <= order_inner; l++)	
	shape[ii++] = Sigma_gradv(ha[l]).Shape();
    };
  };


  
  /* Face basis functions which are normal-tangential continuous */
  /* calculates [(grad l1) o-times (grad l2 x grad l3)] * legendre */
  /* (grad l2 x grad l3) is a scalar!! (cross product on surface */
  
  class T_Dl1_o_Dl2xDl3_v_surf
  {
    AutoDiffDiff<3> l1,l2,l3,v;
  public:
    T_Dl1_o_Dl2xDl3_v_surf  (AutoDiffDiff<3> lam1, AutoDiffDiff<3> lam2, AutoDiffDiff<3> lam3, AutoDiffDiff<3> av) : l1(lam1), l2(lam2), l3(lam3), v(av) { ; }
    
    Vec<2> Shape() {
      double cross = l2.DValue(0)*l3.DValue(1) - l2.DValue(1)*l3.DValue(0);      
      return Vec<2> (v.Value()*l1.DValue(0) * cross,  v.Value()*l1.DValue(1) * cross);
    }

    Vec<2> DivShape()
    {
      throw Exception("not available on surface");
    }

  };
  
  template <> class HCurlDivSurfaceFE<ET_TRIG> : public T_HCurlDivSurfaceFE<ET_TRIG> 
  {
    
  public:
    using T_HCurlDivSurfaceFE<ET_TRIG> :: T_HCurlDivSurfaceFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      ndof += (order_inner+1) * (order_inner+2);
      order = max2(order,order_inner);
    }
    
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {            
      auto x = ip.x, y= ip.y;
      AutoDiffDiff<3> ddlami[3] ={ x, y, 1-x-y };
      
      int ii = 0;

      ArrayMem<AutoDiffDiff<3>,20> ha((order_inner+1)*(order_inner+2)/2.0);
      
      int es = 0, ee = 1, et = 2;
      if(vnums[es] > vnums[ee]) swap(es, ee);
      if(vnums[ee] > vnums[et]) swap(ee, et);
      if(vnums[es] > vnums[et]) swap(es,et);
            
      AutoDiffDiff<3> ls = ddlami[es],le = ddlami[ee], lt = ddlami[et];
      
      DubinerBasis3::Eval (order_inner, ls, le, ha);

      for (int l = 0; l < (order_inner+1)*(order_inner+2)/2.0; l++)
	    {
	      shape[ii++] = T_Dl1_o_Dl2xDl3_v_surf(le,ls,lt,ha[l]).Shape();
	      shape[ii++] = T_Dl1_o_Dl2xDl3_v_surf(ls,lt,le,ha[l]).Shape();
	    }     
    };
      
  };

}


#endif
  
