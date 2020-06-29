#ifndef FILE_HCURLDIVFE
#define FILE_HCURLDIVFE

/*********************************************************************/
/* File:   hcurldivfe.hpp                                            */
/* Author: Philip Lederer                                            */
/* Date:   2017/2018                                                 */
/*********************************************************************/

//include recursive_pol.hpp;
#include "recursive_pol_tet.hpp"
#include "hcurldivfe_impl.hpp"
//#include "l2orth3d.hpp"

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

    virtual void CalcMappedShape (const SIMD<MappedIntegrationPoint<DIM,DIM>> & mip,
                                         BareSliceMatrix<SIMD<double>> shapes) const = 0;

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & ir,
				      BareSliceMatrix<SIMD<double>> shapes) const = 0;
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                                  BareSliceVector<> coefs,
                                  BareSliceMatrix<SIMD<double>> values) const = 0;

    virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir,
                                  BareSliceMatrix<SIMD<double>> values,
                                  BareSliceVector<> coefs) const = 0;

    virtual void CalcMappedDivShape (const SIMD_BaseMappedIntegrationRule & bmir, 
				     BareSliceMatrix<SIMD<double>> divshapes) const = 0;    
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
    int order_trace;
    //bool curlbubbles;
    bool GGbubbles;
    //L2orthTet l2orth;
    
  public:          
    T_HCurlDivFE (int aorder, bool aGGbubbles = false)
    {
      order = aorder;
      for (auto & of : order_facet) of = aorder;
      order_inner = aorder;
      order_trace = -1;
      GGbubbles = aGGbubbles;
    }   
    
    using VertexOrientedFE<ET>::SetVertexNumbers;

    /*template <typename TA>
    void SetVertexNumbers(const TA & avnums)
    {
      //FlatArray<int> aavnums;
      for (int i=0; i<avnums.Size(); i++)
	{
	  //aavnums[i] = avnums[i];
	  vnums[i] = avnums[i];
	}

      l2orth.SetVertexNumbers(avnums);
      }*/
    
    
    virtual ELEMENT_TYPE ElementType() const override { return ET; }
    const HCurlDivFE<ET> * Cast() const { return static_cast<const HCurlDivFE<ET>*> (this); } 
    
    INLINE void SetOrderFacet (int nr, int order) { order_facet[nr] = order; }
    INLINE void SetOrderInner (int order) { order_inner = order; }
    INLINE void SetOrderTrace (int order) { order_trace = order; }

    virtual void ComputeNDof()
    {
      cout << "Error, T_HCurlDivFE<ET>:: ComputeNDof not available, only for ET == TRIG,TET,QUAD" << endl;
    }

    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<double> shape) const override
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
                               BareSliceMatrix<double> shape) const override
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
                               BareSliceMatrix<double> shape) const override
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
                            BareSliceMatrix<double> shape) const override
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

     virtual void CalcMappedShape (const SIMD<MappedIntegrationPoint<DIM,DIM>> & mip, 
                                         BareSliceMatrix<SIMD<double>> shapes) const override
    {
      Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mip;
      TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          
      this->Cast() -> T_CalcShape (addp,
				   SBLambda ([shapes] (size_t j, auto val) 
					     {

					       Vec<DIM*DIM,SIMD<double>> vecshape = val.Shape();
					       for (size_t k = 0; k < DIM*DIM; k++)
						 shapes(j*(DIM*DIM)+k,0) = vecshape(k);
					     }));
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
          Vec<DIM*DIM,SIMD<double>> mat;
	  
	  auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);

	  //for (size_t k = 0; k < DIM*DIM; k++)
	  //  mat(k) = values(k,i);
          mat = values.Col(i);
	  
	  Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir[i];
          TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp(adp);
          double *pcoefs = &coefs(0);
          const size_t dist = coefs.Dist();

          Cast() -> T_CalcShape (addp,
                                 SBLambda ([mat,&pcoefs,dist] (size_t j, auto val)
                                           {                                          
					     Vec<DIM*DIM,SIMD<double>> vecshape = val.Shape();
                                             *pcoefs += HSum(InnerProduct(mat,vecshape));
                                             pcoefs += dist;
                                           }));
        }
    }

    virtual void CalcMappedDivShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                            BareSliceMatrix<double> shape) const override
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

	// see calcMappedDivShape for SIMD for better description
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
                            BareSliceMatrix<double> shape) const override
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


    virtual void CalcMappedDivShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                      BareSliceMatrix<SIMD<double>> divshapes) const override
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);

      if(!mir.GetTransformation().IsCurvedElement()) // non-curved element
      {		
	for (size_t i = 0; i < mir.Size(); i++)
	{
	  /*
          auto jacinv = mir[i].GetJacobianInverse();
          auto d = mir[i].GetJacobiDet();
          
          Vec<DIM,SIMD<double>> vec;
          SIMD<double> mem[DIM*DIM_STRESS];
          FlatMatrix<SIMD<double>> trans(DIM,DIM,&mem[0]);
          trans = 1/d*Trans(jacinv);
	  */
	  
          //Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir.ip()[i];
	  Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir[i];
	  Vec<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp;
	  for (int j=0; j<DIM; j++)
	    {
	      addp[j] = adp[j].Value();
	      addp[j].LoadGradient(&adp[j].DValue(0));
	    }
	  Cast() -> T_CalcShape
            (TIP<DIM,AutoDiffDiff<DIM,SIMD<double>>>(addp),
             //SBLambda([divshapes,i,trans](int j, auto val)
	     SBLambda([divshapes,i](int j, auto val)
                      {
			//divshapes.Rows(j*DIM,(j+1)*DIM).Col(i).AddSize(DIM) = trans * val.DivShape();
                        divshapes.Rows(j*DIM,(j+1)*DIM).Col(i).Range(0,DIM) = val.DivShape();
                      }));
	}
      }
      else
	{	  
	  //throw ExceptionNOSIMD(string("HCurlDiv - CalcMappedDivShape SIMD only for noncurved elements"));
	  for (size_t i = 0; i < mir.Size(); i++)
          {
	    auto mip = mir[i];
	    Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir[i];
	    Mat<DIM,DIM,SIMD<double>> jac = mip.GetJacobian();
	    Mat<DIM,DIM,SIMD<double>> inv_jac = mip.GetJacobianInverse();        
	    Mat<DIM,DIM,SIMD<double>> F_HFinvT_Finv[3];

	    Vec<DIM, Mat<DIM,DIM,SIMD<double>>> hesse;
            mir.GetTransformation().CalcHesse (mir.IR()[i], hesse);

	    Vec<DIM, Mat<DIM,DIM,SIMD<double>>> hesseinvT;	    
	    Vec<DIM, Mat<DIM,DIM,SIMD<double>>> dd_of_F_xi;
	    Vec<DIM, Mat<DIM,DIM,SIMD<double>>> dd_of_FinvT_xi;

	    //revert ordering of CalcHesse
	    //dd_of_J_xi contains der derivative of the Jacobian F with respect to x_i
	    for (int l = 0; l < DIM; l++)
	      {
		for (int j = 0; j < DIM; j++)
		  for (int k = 0; k < DIM; k++)
		    dd_of_F_xi[l](j,k)=hesse[j](l,k);
	      }

	    //use the formula d_xi(F^-T) = - F^-T * Trans(d_xi(F)) * F^-T 
	    for (int l = 0; l < DIM; l++)
	      dd_of_FinvT_xi[l] = -Trans(inv_jac) * Trans(dd_of_F_xi[l]) * Trans(inv_jac);

	    //reorder such that hesseinvT contains derivatives in all directions x_i of one component
	    for (int l = 0; l < DIM; l++)
	      {
		for (int j = 0; j < DIM; j++)
		  for (int k = 0; k < DIM; k++)
		    hesseinvT[l](j,k) = dd_of_FinvT_xi[j](l,k);
	      }	    	    

	    //This comes from the formula of div(sigma) on curved elements
	    for(int j=0; j<DIM; j++)	    
	      F_HFinvT_Finv[j] = jac * Trans(hesseinvT[j]) *  inv_jac;
	    

	    //Vec<DIM,AutoDiff<DIM,SIMD<double>>> adp = mir[i];
	    Vec<DIM,AutoDiffDiff<DIM,SIMD<double>>> addp;

	    for (int j=0; j<DIM; j++)
	    {
	      addp[j] = adp[j].Value();
	      addp[j].LoadGradient(&adp[j].DValue(0));
	    }
	    
	    Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM, SIMD<double>>> (addp),SBLambda([&](int nr,auto val)
									      {
										BareSliceVector<SIMD<double>> divshape = divshapes.Rows(nr*DIM,(nr+1)*DIM).Col(i);
										Vec<DIM,SIMD<double>> div1 = val.DivShape();										
										Vec<DIM*DIM, SIMD<double>> matshape = val.Shape();				    
				    
										for(int k=0; k<DIM; k++)
										  {
										    SIMD<double> sum = div1(k);
										    for(int j=0; j<DIM*DIM; j++)
										      sum += F_HFinvT_Finv[k](j) * matshape(j);
										    divshape(k) = sum;
										  }
									      }));
	    
	  }
	  
	}
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
        ndof += order_facet[i]+1;
        order = max2(order, order_facet[i]);
      }      
      int ninner = 3 * ((order_inner +1) * (order_inner))/2; 
      order = max2(order, order_inner);
      
      ndof += ninner;
      if (order_trace > -1)
	{
	  ndof += (order_trace +1) * (order_trace+2)/2.0;
	  order = max2(order, order_trace);
	}

      if(GGbubbles)
	{
	  ndof += order_inner+1;
	  order +=1;
	}
      
    }
    
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {
      auto x = ip.x, y = ip.y;
      Tx ddlami[3] ={ x, y, 1-x-y };
      
      int ii = 0;
      
      int oi=order_inner;
      int ot=order_trace; 
      int maxorder_facet =
        max2(order_facet[0],max2(order_facet[1],order_facet[2]));      

      ArrayMem<Tx,20> ha(maxorder_facet+1);
      ArrayMem<Tx,20> v(oi+1), u(oi+1);

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
      
      //if trace order is not equal -1
      if (ot>-1)
	{
	  LegendrePolynomial::Eval(ot, 2*lt-1, v);
	  for (int i = 0; i <= ot; i++)
	    shape[ii++] = type4(le, ls, v[i]);

	  IntLegNoBubble::EvalMult (ot, le-lt, 0.25*le*lt, u);
	  LegendrePolynomial::EvalMult(ot, 2*ls-1, ls, v);
	  for(int i = 0; i <= ot-1; i++)
	      for(int j = 0; j+i <= ot-1; j++)       
		    shape[ii++] = Sigma_gradu_v(u[i],v[j]);	  		  	    
	}
                  
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

      
      if(GGbubbles)
	{
	  /*
	  IntLegNoBubble::EvalMult (oi, le-ls, 0.25*le*ls, u);
	  LegendrePolynomial::EvalMult(oi, 2*lt-1, lt, v);

	  //in 2 dimensions u[i]*v[oi-1-i] should be a H1_0 bubble
      
	  for(int i = 0; i <= oi-1; i++)
	    {
	      shape[ii++] = CurlBubble2D(u[i]*v[oi-1-i]);
	    }
	  */
	  //Vector<AutoDiffDiff<2, T> >  l2shape( (oi+1)*(oi+2) / 2 );
	  
	  Vector<Tx>  l2shape( (oi+1)*(oi+2) / 2);
	  DubinerBasis::Eval ( oi , ddlami[0], ddlami[1], l2shape);

	  Vector<Tx>  S(oi+1);

	  S(0) = l2shape(oi);	    
	  for ( int i = 1; i<oi+1; i++)
	      S(i) = l2shape((i+1)*oi-(i-1)*(i)/2);
	  
	  Tx b = ddlami[0]*ddlami[1]*ddlami[2];
	  
	  for (int i=0; i<oi+1;i++)
	    {	     
	      shape[ii++] = GGbubble(S(i), b);
	    }
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
      order += 1; //lo facet functions are linear
      //int ninner = (order_inner+1)*(order_inner+1) + (order_inner+2)*(order_inner) *2;

      int ninner = (order_inner+1)*(order_inner+1);

      if (order_inner > 0)
	ninner += 2 * (order_inner+2)*(order_inner);
      else
	ninner += 2; 
      
      order = max2(order, order_inner);
      order += 2; // we include higher inner bubbles and  tensor space 
      ndof += ninner;

      if (order_trace > -1)
	{
	  ndof += (order_trace+1)*(order_trace+1);
	  order = max2(order, order_trace);
	}
    }
    
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {      
      auto x = ip.x, y = ip.y;
      Tx lx[4] ={1-x, x, x, 1-x};
      Tx ly[4] = {1-y, 1-y, y, y};
      //Tx lam[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};
      Tx edgebubbles[4] = {(1-x)*x, x*(1-x), y*(1-y), (1-y)*y};

      // typedef decltype(x.Value()+x.Value()) T;
     		  
      int ii = 0;
      
      int oi = order_inner;
      int ot = order_trace;
      
      int maxorder_facet  =
        max2(order_facet[3],max2(order_facet[0],max2(order_facet[1],order_facet[2])));
      
      ArrayMem<Tx,20> ha(maxorder_facet+1);
      ArrayMem<Tx,20> v(oi+3), u(oi+3);
      for (int i = 0; i < 4; i++)
        {
	  INT<2> e = ET_trait<ET_QUAD>::GetEdgeSort (i, vnums);
	  
          Tx xi = lx[e[1]]+ly[e[1]]-lx[e[0]]-ly[e[0]];
          Tx eta = lx[e[0]]*ly[e[0]]+lx[e[1]]*ly[e[1]];
	  	 
	  IntLegNoBubble::EvalMult (maxorder_facet , xi, 0.25*edgebubbles[i], ha);	  

	  // edge basis functions given by (on the unit square)
	  // (0, u(x) (1-y),0,0) u(x) \in P^k_x -> basis on lower edge
	  // (0, u(x) y,0,0) u(x) \in P^k_x -> basis on upper edge
	  // (0, v(y) (1-x),0,0) u(x) \in P^k_x -> basis on left edge
	  // (0, v(y) x,0,0) u(x) \in P^k_x -> basis on right edge
	  
          for (int l = 0; l <= order_facet[i]; l++)
	    shape[ii++] = u_Sigma_gradv(eta,ha[l]); //on QUADS edge basisfunctions are NOT div-free
	    //shape[ii++] = Sigma_gradv(eta*ha[l]); //div-free	old version but produces high order terms on the diagonal
        }

      // polynomials in x direction. First one is  quadratic
      IntLegNoBubble::EvalMult (oi+2, lx[0]-lx[1], 0.25*lx[0]*lx[1], u);
      // polynomials in y direction. First one is  quadratic
      IntLegNoBubble::EvalMult (oi+2, ly[0]-ly[2], 0.25*ly[0]*ly[2], v);

      //!!!!!!!!!!!!!!
      // produces a sigma space that is nt-continuous and trace free and includes polynomials up to order k
      // plus certain bubbles
      //!!!!!!!!!!!!!!
      
      //////////////////////////////////////////////////////////////////////////////////////
      // functions on the diagonal split into dev-free part and orthogonal part (if ot > -1)
      //
      // 1.) constants      
      // gives the constant function (-1,0,0,1) -> div-free      
      shape[ii++] = Sigma_gradv(lx[0]*ly[0]); 
      
      //
      // 2.) linear and high order 
      for(int i = 0; i <= oi-1; i++)
       {	 
	 // gives (-d_x u(x)/2,0,0, d_x u(x)/2) -> produces a divergence
	 shape[ii++] = Gradu_Curlv(u[i],ly[0], 1.0); // trace free!
	 // gives (d_y v(x)/2,0,0, -d_y v(x)/2) -> produces a divergence
	 shape[ii++] = Gradu_Curlv(v[i],lx[0], 1.0);

	 // block 1 and 2 produce already (P^k_x,0,0,-P^k_x), (P^k_y,0,0,-P^k_y)
	 // and (ot > -1):  (0,0,0,P^k_x), (0,0,0,P^k_y)
	 // 3.) high order mixed functions
	 // both parts (x or y) of the mixed polynomials on the diagonal are at least linear
	 for(int j = 0; j <= oi-1; j++)
	   {
	     shape[ii++] = Gradu_Curlv(u[i],v[j], 1.0);  
	   }
       }

      // gives the constant function (0,0,0,1) -> div-free
      if (ot>-1)
	{
	  //at least constant
	  shape[ii++] = Sigma_gradu_v(ly[0],lx[0]);
	  // linear and higher
	  for(int i = 0; i <= ot-1; i++)
	    {
	      // gives (0,0,0,- d_x u(x))
	      shape[ii++] = Sigma_gradu_v(ly[0],u[i]); //div-free
	      // gives (-d_y v(y), 0, 0,0)
	      shape[ii++] = Sigma_gradu_v(lx[0],v[i]); //div-free

	      // gives (-d_x u(x) d_y v(y), 0, 0, 0)
	      for(int j = 0; j <= ot-1; j++)
		shape[ii++] = Gradu_Curlv(u[i],v[j], 0.0);
	    }	 
	}
	
      //////////////////////////////////////////////////////////////////////////////////////

      // edges and shapes one the diagonal give:
      // edges: (0,P^k_x o-times P^1_y,0,0)
      //        (0,0,P^1_x o-times P^k_y,0)
      // hence we still need
      //       (0,P^{k}_x o-times P^{k}_y\P^{1}_y, 0,0)
      //       (0,0,P^{k}_y o-times P^{k}_x\P^{1}_y,0)
      // further we include bubbles such that
      //       (0,P^{k+1}_x o-times P^{k-1}_y * (1-y)*y, 0,0)
      //       (0,0,P^{k+1}_y o-times P^{k-1}_x * (1-x)*x,0)
      // is included. Needed for LBB-proof of MCS

      if (oi == 0)
	{
	  shape[ii++] = u_Sigma_gradv(v[0],u[0]);
	  shape[ii++] = u_Sigma_gradv(u[0],v[0]);
	}
      else
	{
	  for(int i = 0; i <= oi+1; i++)
	    {
	      for(int j = 0; j <= oi-1; j++)
		{	  
		  // gives (0 ,-d_xx u(x) v(y),0,0) with v(y) bubble (j=0) or higher (hence inner bubble)
		  shape[ii++] = u_Sigma_gradv(v[j],u[i]);
		  // gives  (0 ,0 , d_yy v(y) u(x), 0) with u(x) bubble or higher (hence inner bubble)	  
		  shape[ii++] = u_Sigma_gradv(u[j],v[i]);
		}
	    }
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
      int ninner = 8.0/6.0* ((order_inner +2) * (order_inner +1) * (order_inner));
      
      order = max2(order, order_inner);
      ndof += ninner;
      
      if (order_trace > -1)
	{
	  ndof += (order_trace +1) * (order_trace+2)* (order_trace+3)/6.0;
	  order = max2(order, order_trace);
	}

      if (GGbubbles)
	{
	  //GG bubbles of jay
	  ndof += 3*(order_inner+1)*(order_inner+2)/2;
	  order +=1;
	}
    }
    
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {
      auto x = ip.x, y = ip.y, z = ip.z ;     
      
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
	  
	  DubinerBasis::Eval (maxorder_facet, ls, le, ha);

          for (int l = 0; l < (p+1)*(p+2)/2.0; l++)
	    {	      
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(le,ls,lt,ha[l]);	      
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(ls,lt,le,ha[l]);
	    }
        }

      int oi=order_inner;
      int ot=order_trace; 
      
      int es = 0; int ee = 1; int et = 2; int eo = 3 ;
      AutoDiff<3,T> ls = ddlami[es]; 
      AutoDiff<3,T> le = ddlami[ee];
      AutoDiff<3,T> lt = ddlami[et];
      AutoDiff<3,T> lo = ddlami[eo];

      LegendrePolynomial leg;
      JacobiPolynomialAlpha jac1(1);


      //############ type 1 ############
      if (ot>-1)
	{
	  ArrayMem<AutoDiff<3,T>,20> dub_vals_inner((ot+1)*(ot+2)*(ot+3)/6.0);      
	  DubinerBasis3D::Eval(ot,ls,le,lt ,dub_vals_inner);
	  
	  for (int l = 0; l < (ot+1)*(ot+2)*(ot+3)/6.0; l++)
	     shape[ii++] =  Id_v(dub_vals_inner[l]); 	  	  
	}
      
      //############ type 2 ############
      int ndof_inner = 1.0/6.0* ((order_inner +2) * (order_inner +1) * (order_inner));
      
      ArrayMem<AutoDiff<3,T>,20> dub_vals(ndof_inner);
      
      DubinerBasis3D::Eval(oi-1,ls,le,lt ,dub_vals);
      
      for (int l = 0; l < ndof_inner; l++)
	    {
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(le,ls,lt,lo*dub_vals[l]);	      
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(ls,lt,le,lo*dub_vals[l]);	  
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(le,ls,lo,lt*dub_vals[l]);	      
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(ls,lo,le,lt*dub_vals[l]);	  
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(le,lo,lt,ls*dub_vals[l]);	      
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(lo,lt,le,ls*dub_vals[l]);	  
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(lo,ls,lt,le*dub_vals[l]);	      
	      shape[ii++] =  T_Dl1_o_Dl2xDl3_v<T>(lt,ls,lo,le*dub_vals[l]);
	    }
      
      if(GGbubbles)
	{
	  AutoDiffDiff<3, T> ax[4] = { x, y, z, 1-x-y-z};
	  
	  Mat<3,3,T> B;
	  for (int i = 0; i<3; i++)
	    for (int j = 0; j<3; j++)
	      B(i,j) = ax[0].Value()*ax[1].Value()*ax[2].Value()*ax[3].DValue(i)*ax[3].DValue(j) + ax[1].Value()*ax[2].Value()*ax[3].Value()*ax[0].DValue(i)*ax[0].DValue(j) +
		       ax[2].Value()*ax[3].Value()*ax[0].Value()*ax[1].DValue(i)*ax[1].DValue(j) + ax[3].Value()*ax[0].Value()*ax[1].Value()*ax[2].DValue(i)*ax[2].DValue(j);

	  AutoDiffDiff<3,T> curlB[3];      
	  for(int i = 0; i<3; i++)
	    {
	      curlB[i] = ax[3].DValue(i) * Cross(ax[0] * ax[1] * ax[2],ax[3]) + ax[0].DValue(i) * Cross(ax[1] * ax[2] * ax[3],ax[0]) +
                         ax[1].DValue(i) * Cross(ax[2] * ax[3] * ax[0],ax[1]) + ax[2].DValue(i) * Cross(ax[3] * ax[0] * ax[1],ax[2]);
	    }

	  Mat<3,3,T> S[3];

	  AutoDiffDiff<3, T> al[6] = { x, y, x, z, y, z};
	  
	  for (int i = 0; i<3; i++)	   
	    {
	      S[i](0,0) = 0; S[i](1,1) = 0; S[i](2,2) = 0;
	      S[i](0,1)= al[2*i].DValue(0) * al[2*i+1].DValue(1) - al[2*i].DValue(1) * al[2*i+1].DValue(0); 
	      S[i](0,2)= al[2*i].DValue(0) * al[2*i+1].DValue(2) - al[2*i].DValue(2) * al[2*i+1].DValue(0);
	      S[i](1,2)= al[2*i].DValue(1) * al[2*i+1].DValue(2) - al[2*i].DValue(2) * al[2*i+1].DValue(1);
	      
	      S[i](1,0) = -S[i](0,1); S[i](2,0) = -S[i](0,2); S[i](2,1) = -S[i](1,2);
	    }	  
	  ArrayMem<AutoDiffDiff<3,T>,20> highest_dub_vals_inner((oi+1)*(oi+2)/2);
	  
	  DubinerBasis3D::EvalHighestOrder(oi,ax[0],ax[1],ax[2] ,highest_dub_vals_inner);
	  
	  for (int l = 0; l < (oi+1)*(oi+2)/2; l++)
	    {
	      shape[ii++] = GGbubble_3D(highest_dub_vals_inner[l], S[0], B, curlB);
	      shape[ii++] = GGbubble_3D(highest_dub_vals_inner[l], S[1], B, curlB);
	      shape[ii++] = GGbubble_3D(highest_dub_vals_inner[l], S[2], B, curlB);
	    }	 
	}
     };
  };

  template <> class HCurlDivFE<ET_HEX> : public T_HCurlDivFE<ET_HEX> 
  {
    
  public:
    using T_HCurlDivFE<ET_HEX> :: T_HCurlDivFE;

    virtual void ComputeNDof()
    {      
      order = 0;
      ndof = 0;
      for (int i=0; i<6; i++)
      {
        ndof += 2*(order_facet[i]+1)*(order_facet[i]+1);
        order = max2(order, order_facet[i]+1);
      }
      int p = order_inner;
      int ninner = 2 * (p+1) * (p+1) * (p+1);

      if (p > 0)
	ninner += 6 * (p+2) * (p+1) * (p);
      else
	ninner += 6;
      
      
      order = max2(order, order_inner);
      ndof += ninner;
      
      if (order_trace > -1)
	{
	  ndof += (order_trace +1) * (order_trace+1)* (order_trace+1);
	  order = max2(order, order_trace);
	}
      order+=2;
      if (GGbubbles)
	{
	  throw Exception ("GGBubbles not implemented for Hcurldiv on HEXES");
	}
    }
    
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {
      auto x = ip.x, y = ip.y, z = ip.z ;     
      
      typedef decltype(x.Value()+x.Value()) T;      
      AutoDiff<3,T> xx(x.Value(), &x.DValue(0));
      AutoDiff<3,T> yy(y.Value(), &y.DValue(0));
      AutoDiff<3,T> zz(z.Value(), &z.DValue(0));      
      AutoDiff<3,T> lx[2] ={ 1-xx, xx};
      AutoDiff<3,T> ly[2] ={ 1-yy, yy};
      AutoDiff<3,T> lz[2] ={ 1-zz, zz};
      AutoDiff<3,T> lami[8]={(1-xx)*(1-yy)*(1-zz),xx*(1-yy)*(1-zz),xx*yy*(1-zz),(1-xx)*yy*(1-zz),
                  (1-xx)*(1-yy)*zz,xx*(1-yy)*zz,xx*yy*zz,(1-xx)*yy*zz};
      
      AutoDiff<3,T> sigma[8] = {1-xx + 1-yy + 1-zz, xx + 1-yy + 1-zz, xx + yy + 1-zz, 1-xx + yy + 1-zz,
				1-xx + 1-yy + zz, xx + 1-yy + zz, xx + yy + zz, 1-xx + yy + zz};
      
      int ii = 0;
      
      int maxorder_facet =
        max2(order_facet[0],max2(order_facet[1],order_facet[2]));

      const FACE * faces = ElementTopology::GetFaces(ET_HEX);
             
      ArrayMem<AutoDiff<3,T>,20> leg_u(order+2), leg_v(order+2), leg_w(order+2);

      for (int i = 0; i<6; i++)
        {
          int p = order_facet[i];
          
          AutoDiff<3,T> lam_f(0);
          for (int j = 0; j < 4; j++)
            lam_f += lami[faces[i][j]];
          
          INT<4> f = ET_trait<ET_HEX>::GetFaceSort (i, vnums);	  
          AutoDiff<3,T> xi  = sigma[f[0]] - sigma[f[1]]; 
          AutoDiff<3,T> eta = sigma[f[0]] - sigma[f[3]];
          //auto nv = GetGradient(lam_f);
	  
	  LegendrePolynomial (p+1, eta, leg_u);
          LegendrePolynomial (p+1, xi, leg_v);
          
          for (int j = 0; j <= p; j++)
            for (int k = 0; k <= p; k++)
	      {
		shape[ii++] = T_Dl1_o_Dl2xDl3_v<T>(xi,xi,eta,lam_f*leg_u[j]*leg_v[k]);
		shape[ii++] = T_Dl1_o_Dl2xDl3_v<T>(eta,xi,eta,lam_f*leg_u[j]*leg_v[k]);
		//shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(xi,lam_f,lam_f*leg_u[j]*leg_v[k]);  // not div-free
		//shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(eta,lam_f,lam_f*leg_u[j]*leg_v[k]); // not div-free
	      }	            
        }

      int oi=order_inner;
      int ot=order_trace;
      
      int p = max(oi,ot);
      
      AutoDiff<3,T> xi  = sigma[0] - sigma[1];
      AutoDiff<3,T> eta = sigma[0] - sigma[3];
      AutoDiff<3,T> nv = sigma[0] - sigma[4];
          
      LegendrePolynomial (p+2, xi,  leg_u);
      LegendrePolynomial (p+2, eta, leg_v);
      LegendrePolynomial (p+2, nv,  leg_w);

      
      //polynomials on the diagonal
      // 2 * (oi+1)**3
      // with trace + (ot+1)**3
      for (int i = 0; i <= oi; i++)
	for (int j = 0; j <= oi; j++)
	  for (int k = 0; k <= oi; k++)
                {
		  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(xi,xi,leg_u[i]*leg_v[j]*leg_w[k]);
		  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(eta,eta,leg_u[i]*leg_v[j]*leg_w[k]);
		  //if (ot>-1)
		  //  shape[ii++] = Id_v(leg_u[i]*leg_v[j]*leg_w[k]);
		}

      
      if (ot>-1)
	{
	  for (int i = 0; i <= ot; i++)
	    for (int j = 0; j <= ot; j++)
	      for (int k = 0; k <= ot; k++)
		shape[ii++] = Id_v(leg_u[i]*leg_v[j]*leg_w[k]);
	}
      
      if (oi == 0)
	{
	  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(xi,nv, leg_u[0]*leg_v[0]*leg_w[0] * lz[0] * lz[1]);
	  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(eta,nv,leg_u[0]*leg_v[0]*leg_w[0] * lz[0] * lz[1]);
	  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(xi,eta,leg_u[0]*leg_v[0]*leg_w[0] * ly[0] * ly[1]);
	  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(nv,eta,leg_u[0]*leg_v[0]*leg_w[0] * ly[0] * ly[1]);
	  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(eta,xi,leg_u[0]*leg_v[0]*leg_w[0] * lx[0] * lx[1]);
	  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(nv,xi, leg_u[0]*leg_v[0]*leg_w[0] * lx[0] * lx[1]);
	}
      else
	{       
	  //polynomials on the off diagonal that are nt-bubbles
	  for (int i = 0; i <= oi+1; i++)
	    for (int j = 0; j <= oi; j++)
	      for (int k = 0; k <= oi-1; k++)
                {
		  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(xi,nv, leg_u[i]*leg_v[j]*leg_w[k] * lz[0] * lz[1]);
		  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(eta,nv,leg_u[j]*leg_v[i]*leg_w[k] * lz[0] * lz[1]);

		  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(xi,eta,leg_u[i]*leg_v[k]*leg_w[j]  * ly[0] * ly[1]);
		  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(nv,eta,leg_u[j]*leg_v[k]*leg_w[i]  * ly[0] * ly[1]);

		  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(eta,xi,leg_u[k]*leg_v[i]*leg_w[j]  * lx[0] * lx[1]);
		  shape[ii++] = T_dev_Dl1_o_Dl2_v<T>(nv,xi, leg_u[k]*leg_v[j]*leg_w[i]  * lx[0] * lx[1]);
		  
		}
      
	}	  
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
					    if (DIM==1)
					      shape.Row(nr).AddSize(DIM_STRESS) = val;
					    else
					      shape.Row(nr).AddSize(DIM+1) = val;					    
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
	  throw ExceptionNOSIMD(string("HCurlDiv - CalcMappedShape on surface elements only on (surface)DIM==1"));
	//shape.Row(nr).AddSize(1) = val;
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
    AutoDiff<3> l1,l2,l3,v;
  public:
    T_Dl1_o_Dl2xDl3_v_surf  (AutoDiff<3> lam1, AutoDiff<3> lam2, AutoDiff<3> lam3, AutoDiff<3> av) : l1(lam1), l2(lam2), l3(lam3), v(av) { ; }
    
    Vec<2> Shape() {
      double cross = l2.DValue(0)*l3.DValue(1) - l2.DValue(1)*l3.DValue(0);      
      return Vec<2> (v.Value()*l1.DValue(0) * cross,  v.Value()*l1.DValue(1) * cross);
    }

    Vec<2> DivShape()
    {
      throw Exception("not available on surface");
    }

  };

  
  class T_dev_Dl1_o_Dl2_v_surf 
  {
    AutoDiff<3> l1,l2,v;
  public:
    T_dev_Dl1_o_Dl2_v_surf   (AutoDiff<3> lam1,AutoDiff<3> lam2,  AutoDiff<3> av) : l1(lam1), l2(lam2), v(av) { ; }
    
    Vec<2> Shape() {
      //double cross = l2.DValue(2);      
      return Vec<2> (v.Value()*l1.DValue(0) * l2.DValue(2),  v.Value()*l1.DValue(1) * l2.DValue(2));
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
      AutoDiff<3> xx(x.Value(), &x.DValue(0));
      AutoDiff<3> yy(y.Value(), &y.DValue(0));
      
      AutoDiff<3> ddlami[3] ={ xx, yy, 1-xx-yy };      
      
      int ii = 0;

      ArrayMem<AutoDiff<3>,20> ha((order_inner+1)*(order_inner+2)/2.0);
      
      int es = 0, ee = 1, et = 2;
      if(vnums[es] > vnums[ee]) swap(es, ee);
      if(vnums[ee] > vnums[et]) swap(ee, et);
      if(vnums[es] > vnums[et]) swap(es,et);
            
      AutoDiff<3> ls = ddlami[es],le = ddlami[ee], lt = ddlami[et];
      
      DubinerBasis::Eval (order_inner, ls, le, ha);

      for (int l = 0; l < (order_inner+1)*(order_inner+2)/2.0; l++)
	    {
	      shape[ii++] = T_Dl1_o_Dl2xDl3_v_surf(le,ls,lt,ha[l]).Shape();
	      shape[ii++] = T_Dl1_o_Dl2xDl3_v_surf(ls,lt,le,ha[l]).Shape();
	    }     
    };
      
  };
  
  template <> class HCurlDivSurfaceFE<ET_QUAD> : public T_HCurlDivSurfaceFE<ET_QUAD> 
  {
    
  public:
    using T_HCurlDivSurfaceFE<ET_QUAD> :: T_HCurlDivSurfaceFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      ndof += 2* (order_inner+1) * (order_inner+1);
      order = max2(order,order_inner);
    }
    
   template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip/*AutoDiffDiff<2> hx[2]*/, TFA & shape) const
    {            
      auto x = ip.x, z= ip.y;      
      //AutoDiffDiff<3> ddlami[3] ={ x, y, 1-x-y };      
      typedef decltype(x.Value()+x.Value()) T;   
      AutoDiff<3> xx(x.Value(), &x.DValue(0));
      AutoDiff<3> zz(z.Value(), &z.DValue(0));
      
      AutoDiff<3,T> lami[4]={(1-xx)*(1-zz),xx*(1-zz),xx*zz,(1-xx)*zz};
      
      AutoDiff<3> sigma[4] = {1-xx+1-zz, xx+1-zz, xx+zz, 1-xx+zz};

      ArrayMem<AutoDiff<3>,20> leg_u(order_inner+2);
      ArrayMem<AutoDiff<3>,20> leg_v(order_inner+2);
      
      int ii = 0;

      INT<4> f = ET_trait<ET_QUAD>::GetFaceSort (0, vnums);

      //AutoDiff<3,T> lam_f(0);
      //for (int j = 0; j < 4; j++)
      //lam_f += lami[j];

      AutoDiff<3,T> xi  = sigma[f[0]] - sigma[f[1]]; 
      AutoDiff<3,T> eta = sigma[f[0]] - sigma[f[3]];
	       	    
      LegendrePolynomial::Eval(order_inner+1,eta,leg_u);
      LegendrePolynomial::Eval(order_inner+1,xi,leg_v);

      for(int j = 0; j <= order_inner; j++)
        for(int k = 0; k <= order_inner; k++)
	  {	     
	    shape[ii++] = T_Dl1_o_Dl2xDl3_v_surf(xi,xi,eta,leg_u[j] * leg_v[k]).Shape();
	    shape[ii++] =  T_Dl1_o_Dl2xDl3_v_surf(eta,xi,eta,leg_u[j] * leg_v[k]).Shape();	    
	    //shape[ii++] =  T_dev_Dl1_o_Dl2_v_surf(xi,lam_f,leg_u[k] * leg_v[k]).Shape();
	    //shape[ii++] =  T_dev_Dl1_o_Dl2_v_surf(eta,lam_f,leg_u[k] * leg_v[k]).Shape();	    
	  }
    };
      
  };

}


#endif
  
