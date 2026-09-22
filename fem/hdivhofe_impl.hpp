#ifndef FILE_HDIVHOFE_IMPL
#define FILE_HDIVHOFE_IMPL

/*********************************************************************/
/* File:   hdivhofe_impl.hpp                                         */
/* Author: A. Becirovic, S. Zaglmayr, J. Schoeberl                   */
/* Date:   15. Feb. 2003                                             */
/*********************************************************************/

#include "recursive_pol_tet.hpp"

namespace ngfem
{


  template <>
  class HDivHighOrderFE_Shape<ET_TRIG> : public HDivHighOrderFE<ET_TRIG>
  {
  public:
    template<typename Tx, typename TFA>  
    // void T_CalcShape (Tx hx[2], TFA & shape) const;
    INLINE void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const;

    template <typename MIP, typename TFA>
    inline void CalcDualShape2 (const MIP & mip, TFA & shape) const;
  };

  template <>
  class HDivHighOrderFE_Shape<ET_QUAD> : public HDivHighOrderFE<ET_QUAD>
  {
  public:
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const;

    template <typename MIP, typename TFA>
    inline void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception(string("CalcDualShape missing for HighOrderHDiv element ")+ElementTopology::GetElementName(ET_QUAD));
    }
  };

  template<> 
  class HDivHighOrderFE_Shape<ET_TET> : public HDivHighOrderFE<ET_TET>
  {
    typedef TetShapesInnerLegendre T_INNERSHAPES;
    typedef TetShapesFaceLegendre T_FACESHAPES; 
  public:
    template<typename Tx, typename TFA>  
    inline void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const;
    template <typename MIP, typename TFA>
    inline void CalcDualShape2 (const MIP & mip, TFA & shape) const;
  };

  template<>
  class HDivHighOrderFE_Shape<ET_PRISM> : public HDivHighOrderFE<ET_PRISM>
  {
    typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const;

    template <typename MIP, typename TFA>
    inline void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception(string("CalcDualShape missing for HighOrderHDiv element ")+ElementTopology::GetElementName(ET_PRISM));
    }
  };

  template<>
  class HDivHighOrderFE_Shape<ET_PYRAMID> : public HDivHighOrderFE<ET_PYRAMID>
  {
  public:
    template<typename Tx, typename TFA>
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const;

    template <typename MIP, typename TFA>
    inline void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception(string("CalcDualShape missing for HighOrderHDiv element ")+ElementTopology::GetElementName(ET_PYRAMID));
    }
  };

  template<>
  class HDivHighOrderFE_Shape<ET_HEX> : public HDivHighOrderFE<ET_HEX>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const;

    template <typename MIP, typename TFA>
    inline void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception(string("CalcDualShape missing for HighOrderHDiv element ")+ElementTopology::GetElementName(ET_HEX));
    }
  };






  template<typename Tx, typename TFA>  
  INLINE void  HDivHighOrderFE_Shape<ET_TRIG> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {
    if (only_ho_div && (order_inner[0] <= 1)) return;

    Tx lam[3] = { ip.x, ip.y, 1-ip.x-ip.y };

    int ii = 3; 
    if (!only_ho_div)
      {
        for (int i = 0; i < 3; i++)
          {
            IVec<2> e = ET_trait<ET_TRIG>::GetEdgeSort (i, vnums);

            //Nedelec low order edge shape function 
            shape[i] = uDv_minus_vDu (lam[e[0]], lam[e[1]]);

            int p = order_facet[i][0]; 
            //HO-Edge shapes (Gradient Fields)   
            if(p > 0) //  && usegrad_edge[i]) 
              { 
                Tx xi = lam[e[1]] - lam[e[0]]; 
                
                // LegendrePolynomial::
                IntLegNoBubble::
                  EvalScaledMult (p-1, xi, lam[e[0]]+lam[e[1]], 
                                  lam[e[0]]*lam[e[1]], 
                                  SBLambda([&](int i, Tx v)
                                           {
                                             shape[ii++] = Du(v);
                                           }));
              }
          }   
      }
    else
      ii = 0;

    //Inner shapes (Face) 
    int p = order_inner[0];
    int pd = p;
    if (RT) pd++;
    
    if (pd > 1) 
      {
        IVec<4> fav = ET_trait<ET_TRIG>::GetFaceSort (0, vnums);

 	// rotated gradients:
	//if(!only_ho_div)
	if(p>1 && !only_ho_div)
          {
            DubinerBasis::EvalMult (p-2, lam[fav[0]], lam[fav[1]], 
                                    lam[fav[0]]*lam[fav[1]]*lam[fav[2]], 
                                    SBLambda ([&](int nr, Tx val)
                                              {
                                                shape[ii++] = Du(val);
                                              }));
          }

        if (!ho_div_free)
          {
            // Tx x = lam[fav[0]];
            // Tx y = lam[fav[1]];
            LegendrePolynomial leg;
            // IntLegNoBubble leg;
            /*
            leg.EvalScaledMult1Assign 
              (pd-2, y-(1-x-y), 1-x, y*(1-x-y),
               SBLambda ([&] (int i, Tx val1) LAMBDA_INLINE 
                         {
                           JacobiPolynomialAlpha jac(1+2*i);
                           jac.EvalMult1Assign 
                             (pd-2-i, 2*x-1, x, 
                              SBLambda([&](int j, Tx val2) 
                                       {
                                         shape[ii++] = uDv_minus_vDu (val1,val2);
                                       }));
                         }));
            
            // rec_pol * Nedelec0 
            leg.EvalMult 
              (pd-2, 2*x-1, x, 
               SBLambda([&] (int j, Tx val)
                        {
                          shape[ii++] = wuDv_minus_wvDu (lam[fav[1]], lam[fav[2]], val);
                        }));
            */
            // the correct RT !
            auto l0 = lam[0];
            auto l1 = lam[1];
            auto l2 = lam[2];
            
            DubinerBasis::EvalMult (pd-2, lam[0], lam[1], 
                                    lam[0], 
                                    SBLambda ([l1,l2,shape,&ii](int nr, Tx val)
                                               {
                                                 shape[ii++] = wuDv_minus_wvDu (l1, l2, val);
                                               }));

            LegendrePolynomial::EvalScaledMult 
              (pd-2, l2-l1, l1+l2, l2, 
               SBLambda([shape,l0,l1,&ii] (int j, Tx val)
                        {
                          shape[ii++] = wuDv_minus_wvDu (l1, l0, val);
                        }));
          }

      }
  }




template <typename MIP, typename TFA>
  inline void HDivHighOrderFE_Shape<ET_TRIG>::CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    auto & ip = mip.IP();
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
    T x = ip(0), y = ip(1);
    T lam[3] = { x, y, 1-x-y };
    Vec<2,T> pnts[3] = { { 1, 0 }, { 0, 1 } , { 0, 0 } };
    int facetnr = ip.FacetNr();

    int ii = 3;

    if (ip.VB() == BND)
      { // facet shapes
        for (int i = 0; i < 3; i++)
          {
            int p = order_facet[i][0];
            if (i == facetnr)
              {

                IVec<2> e = GetEdgeSort (i, vnums);
                T xi = lam[e[1]]-lam[e[0]];
                Vec<2,T> tauref = pnts[e[1]] - pnts[e[0]];
                
		Vec<2,T> nvref = Vec<2,T>(tauref[1],-tauref[0]);
                auto nv = L2Norm(nvref) / L2Norm(Trans(mip.GetJacobianInverse())*nvref) * Trans(mip.GetJacobianInverse())*nvref; // latest version, see TET for derivation
                
		LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              auto vshape = val * nv;
                              if (nr==0)
                                shape[i] = vshape;
                              else
                                shape[ii+nr-1] = vshape;
				}));
              }
            ii += p;
          }
      }
    else
      {
	for (int i = 0; i < 3; i++)
          ii += order_facet[i][0];
      }
    if (ip.VB() == VOL)
      {
	DubinerBasis::Eval(order_inner[0]-2+RT, x, y,
                           SBLambda([&] (size_t nr, auto val)
                                    {
                                      shape[ii++] = Trans(mip.GetJacobianInverse())*Vec<2,T> (val, 0);
				      if(RT)
					shape[ii++] = Trans(mip.GetJacobianInverse())*Vec<2,T> (0,val);
				      else
					shape[ii++] = Trans(mip.GetJacobianInverse())*Vec<2,T> (val*y, -val*x);
                                    }));

	if(!RT)
	  LegendrePolynomial::Eval(order_inner[0]-2,y,
				 SBLambda([&] (size_t nr, auto val)
					  {
					    shape[ii++] = Trans(mip.GetJacobianInverse())*Vec<2,T>(0,val);
					  }));

      }
  }





  template<typename Tx, typename TFA>  
  void  HDivHighOrderFE_Shape<ET_QUAD> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {

    if (only_ho_div && (order_inner[0]<=1 && order_inner[1]<=1)) return;
    Tx x = ip.x, y = ip.y;

    Tx lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  

    size_t ii = 4;
    STACK_ARRAY(Tx, mem, 2*order+4);
    Tx * pol_xi = &mem[0];
    Tx * pol_eta = &mem[order+2];
    
    if (!only_ho_div){
      // edges
      // const EDGE * edges = ElementTopology::GetEdges (ET_QUAD);
      for (int i = 0; i < 4; i++)
        {
          int p = order_facet[i][0]; 
          // int es = edges[i][0], ee = edges[i][1];
          // if (vnums[es] > vnums[ee]) swap (es, ee);
          IVec<2> e = GetVertexOrientedEdge (i);	  

          Tx xi  = sigma[e[1]]-sigma[e[0]];
          Tx lam_e = lami[e[1]]+lami[e[0]]; 

          // Nedelec0-shapes
          shape[i] = uDv (0.5 * lam_e, xi); 
          
          // High Order edges ... Gradient fields 
          // if(usegrad_edge[i])
          {
            IntLegNoBubble::
              EvalMult (p-1, xi, 0.25*(1-xi*xi)*lam_e,
                        SBLambda([&] (size_t j, Tx val)
                                 {
                                   shape[ii++] = Du(val);
                               }));
          }
        }
    }
    else
      ii = 0;
    
    // IVec<2> p = order_face[0]; // (order_cell[0],order_cell[1]);
    IVec<2> p (order_inner[0], order_inner[1]); // (order_cell[0],order_cell[1]);
    int fmax = 0; 
    for (int j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax])
	fmax = j;
    
    int f1 = (fmax+3)%4; 
    int f2 = (fmax+1)%4; 
    if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f2 > f1; 
    
    Tx xi = sigma[fmax]-sigma[f1];  // in [-1,1]
    Tx eta = sigma[fmax]-sigma[f2]; // in [-1,1]

    /*
    T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
    T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);
    */
    IntLegNoBubble::EvalMult (p[0]-1, xi, 0.25*(1-xi*xi), pol_xi);
    IntLegNoBubble::EvalMult (p[1]-1, eta, 0.25*(1-eta*eta), pol_eta);
    if (!only_ho_div)
    {    
      //Gradient fields 
      // if(usegrad_face[0])
      for (size_t k = 0; k < p[0]; k++)
        {
          auto val_k = pol_xi[k];
          for (size_t j= 0; j < p[1]; j++)
            shape[ii++] = Du(val_k*pol_eta[j]);
        }
    }

    if (!ho_div_free)
      {
        //Rotation of Gradient fields 
        for (int k = 0; k < p[0]; k++)
          for (int j= 0; j < p[1]; j++)
            shape[ii++] = uDv_minus_vDu (pol_eta[j], pol_xi[k]);

        //Missing ones 
        for(int j = 0; j< p[0]; j++)
          shape[ii++] = uDv (0.5*pol_xi[j], eta);
        
        for(int j = 0; j < p[1]; j++)
          shape[ii++] = uDv (0.5*pol_eta[j], xi); 
      }
  }














  
  /// compute shape
  template<typename Tx, typename TFA>  
  inline void HDivHighOrderFE_Shape<ET_TET> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    if (only_ho_div && order_inner[0]<=1) return;
    Tx x = ip.x, y = ip.y, z = ip.z;
    Tx lami[4] = { x, y, z, 1-x-y-z };

	
    size_t ii = 4; 
    if (!only_ho_div)
    {
      const FACE * faces = ElementTopology::GetFaces (ET_TET).Data();
      for (size_t i = 0; i < 4; i++)
        {
          int p = order_facet[i][0];

          int fav[3];
          for(int j = 0; j < 3; j++) fav[j]=faces[i][j];
          
          //Sort vertices  first edge op minimal vertex
          if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
          if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
          if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);
          int fop = 6 - fav[0] - fav[1] - fav[2];
          
          // RT lowest order
          shape[i] = uDvDw_Cyclic (lami[fav[0]], lami[fav[1]], lami[fav[2]]);

          Tx xi = lami[fav[1]]-lami[fav[0]];
          Tx sum = lami[fav[1]]+lami[fav[0]];
          Tx bub = lami[fav[1]]*lami[fav[0]];
          Tx eta = lami[fav[2]];
          Tx zeta = lami[fop];  

          // Typ 1
          IntLegNoBubble::
            EvalScaledMult (p-1, xi, sum, bub, 
                            SBLambda([&](int k, Tx polk) LAMBDA_INLINE
                                     {
                                       IntegratedJacobiPolynomialAlpha jac(2*k+3);
                                       jac.EvalScaledMult(p-k-1, eta-sum, 1-zeta, eta, 
                                                          SBLambda ([&](int i, Tx val) LAMBDA_INLINE
                                                                    {
                                                                      shape[ii++] = 
                                                                        Du_Cross_Dv (val, polk);
                                                                    }));
                                     }));
          
          IntegratedJacobiPolynomialAlpha jac(3);
          jac.EvalScaledMult(p-1, eta-sum, 1-zeta, eta, 
                             SBLambda ([&](int i, Tx val) LAMBDA_INLINE
                                       {
                                         shape[ii++] = 
                                           curl_uDvw_minus_Duvw (lami[fav[0]], lami[fav[1]], val);
                                       }));
        }
    }
    else
      ii = 0;
    // cell-based shapes 
    int p = order_inner[0];
    int pc = order_inner[0]; // should be order_inner_curl  
    // int pp = max2(p,pc);
    int pd = p;
    if (RT) pd++;
    if ( order >= 2 )
      {
        STACK_ARRAY(Tx, mem, 3*order);
        Tx * adpol1 = &mem[0];
        Tx * adpol2 = &mem[order];
        Tx * adpol3 = &mem[2*order];

        T_INNERSHAPES::CalcSplitted(order+2, lami[0]-lami[3], lami[1], lami[2], adpol1, adpol2, adpol3 );
      
        if (pc >= 2 && !only_ho_div){
          // Curl-Fields 
          for (int i = 0; i <= pc-2; i++)
            for (int j = 0; j <= pc-2-i; j++)
              for (int k = 0; k <= pc-2-i-j; k++)
                {
                  // grad v  x  grad (uw)
                  shape[ii++] = Du_Cross_Dv (adpol2[j], adpol1[i]*adpol3[k]);
        
                  // grad w  x  grad (uv)
                  shape[ii++] = Du_Cross_Dv (adpol3[k], adpol1[i]*adpol2[j]);
                }     


          // Type 1 : Curl(T3)
          // ned = lami[0] * nabla(lami[3]) - lami[3] * nabla(lami[0]) 
          for (int j= 0; j <= pc-2; j++)
            for (int k = 0; k <= pc-2-j; k++)
              shape[ii++] = curl_uDvw_minus_Duvw (lami[0], lami[3], adpol2[j]*adpol3[k]);
        }

        if (!ho_div_free)
          {
            // Type 2:  
            // (grad u  x  grad v) w
            /*
            for (int i = 0; i <= pd-2; i++)
              for (int j = 0; j <= pd-2-i; j++)
                for (int k = 0; k <= pd-2-i-j; k++)
                  shape[ii++] = wDu_Cross_Dv<3> (adpol1[i], adpol2[j], adpol3[k]);

            // (ned0 x grad v) w    
            for (int j = 0; j <= pd-2; j++)
              for (int k= 0; k <= pd-2-j; k++)
                shape[ii++] = wDu_Cross_Dv<3> (lami[0], adpol2[j], lami[3]*adpol3[k]);
            
            // Type 3: 
            // (ned0 x e_z) v = (N_y, -N_x,0)^T * v ) 
            for (int j=0; j<=pd-2; j++) 
              shape[ii++] = wDu_Cross_Dv<3> (lami[0], z, lami[3]*adpol2[j]);
            */



            DubinerBasis3D::EvalMult (pd-2, lami[0], lami[1], lami[2],
                                      lami[0], 
                                      SBLambda ([&](int nr, Tx val)
                                                {
                                                  shape[ii++] = z_times_uDvDw_Cyclic (lami[1], lami[2], lami[3], val);
                                                }));


            DubinerBasis::EvalScaledMult (pd-2, lami[1], lami[2], 1-lami[0], 
                                          lami[1],
                                          SBLambda ([&](int nr, Tx val)
                                                    {
                                                      shape[ii++] = z_times_uDvDw_Cyclic (lami[0], lami[2], lami[3], val);
                                                    }));
            
            LegendrePolynomial::EvalScaledMult 
              (pd-2, lami[2]-lami[3], lami[2]+lami[3], lami[2], 
               SBLambda([&] (int j, Tx val)
                        {
                          shape[ii++] = z_times_uDvDw_Cyclic (lami[0], lami[1], lami[3], val);                          
                        }));
          }
      }
  }

  



template <typename MIP, typename TFA>
  inline void HDivHighOrderFE_Shape<ET_TET> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;        
    auto & ip = mip.IP();
    T x = ip(0), y = ip(1), z = ip(2);
    T lam[4] = { x, y, z, 1-x-y-z };
    Vec<3> pnts[4] = { { 1, 0, 0 }, { 0, 1, 0 } , { 0, 0, 1 }, { 0, 0, 0 } };
    int facetnr = ip.FacetNr();
    int ii = 4;

    if (ip.VB() == BND)
      { // facet shapes
        for (int i = 0; i < 4; i++)
          {
            int p = order_facet[i][0];
	    if (i == facetnr)
              {
		IVec<4> fav = GetFaceSort (i, vnums);
                T xi = lam[fav[0]]-lam[fav[2]];
                T eta = lam[fav[1]]-lam[fav[2]];
                Vec<3,T> tauref1 = pnts[fav[1]] - pnts[fav[0]];
		Vec<3,T> tauref2 = pnts[fav[2]] - pnts[fav[0]];
		Vec<3,T> nvref   = Cross(tauref1,tauref2);
                


                // Rewrite new version: use often that Cof(F) = Det(F)*F^{-T}
                // Cross(F*tauref1, F*tauref2): Cross(F*tauref1,F*tauref2)*v = Det(F*tauref1,F*tauref2,v) = Det(F)*Det(tauref1,tauref2,F^{-1}v) = Det(F)F^{-T}Cross(tauref1,tauref2)*v = Cof(F)*nvref*v
                // -> Cross(F*tauref1, F*tauref2) = Cof(F)*nvref
                // -> nv = L2Norm(nvref) / L2Norm(Cof(F)*nvref) * Cof(F) * nvref;
                // 1/L2Norm(Cof(F)*nvref) * Cof(F)*nvref = Det(F)/Det(F) * 1/L2Norm(Cof(F)*nvref) * Cof(F)*nvref = 1/L2Norm(1/Det(F)*Cof(F)*nvref) * 1/Det(F)*Cof(F)*nvref
                //                                       = 1/L2Norm(F^{-T}*nvref) * F^{-T}*nvref
                // -> nv = L2Norm(nvref)/L2Norm(F^{-T}*nvref) * F^{-T}*nvref
                Vec<3,T> nv = L2Norm(nvref) / L2Norm(Trans(mip.GetJacobianInverse())*nvref) * Trans(mip.GetJacobianInverse())*nvref; // latest version 
                
		DubinerBasis::Eval(order_facet[i][0], xi, eta,
                                   SBLambda([&] (size_t nr, auto val)
                                            {
                                              Vec<3,T> vshape = val * nv;
                                              if (nr==0)
                                                shape[i] = vshape;
                                              else
                                                shape[ii+nr-1] = vshape;
					     }));
              }
            ii += (p+1)*(p+2)/2-1;
          }
      }
    else
      {
	for (int i = 0; i < 4; i++)
          ii += (order_facet[i][0]+1)*(order_facet[i][0]+2)/2-1;
      }

    if (ip.VB() == VOL && order > 1)
      {
	DubinerBasis3D::Eval (order-2, lam[0], lam[1], lam[2],
			      SBLambda([&](size_t j, T val) LAMBDA_INLINE
			      {
				if(RT)
				  {
				    shape[ii++] = Trans(mip.GetJacobianInverse())*Vec<3,T>(val,0,0);
				    shape[ii++] = Trans(mip.GetJacobianInverse())*Vec<3,T>(0,val,0);
				    shape[ii++] = Trans(mip.GetJacobianInverse())*Vec<3,T>(0,0,val);
				  }
				else
				  {
				    shape[ii++] = Trans(mip.GetJacobianInverse())*Cross(Vec<3,T>(val,0,0),Vec<3,T>(x,y,z));
				    shape[ii++] = Trans(mip.GetJacobianInverse())*Cross(Vec<3,T>(0,val,0),Vec<3,T>(x,y,z));
				    shape[ii++] = Trans(mip.GetJacobianInverse())*Vec<3,T>(val, 0, 0);
				  }
			      }));
	if (!RT)
	  {
	    DubinerBasis::Eval(order_inner[0]-2, x, y,
                           SBLambda([&] (size_t nr, auto val)
                                    {
                                      shape[ii++] = Trans(mip.GetJacobianInverse())*Cross(Vec<3,T>(0,0,val),Vec<3,T>(x,y,z));
                                    }));
        
	    DubinerBasis::Eval(order_inner[0]-2, y, z,
                           SBLambda([&] (size_t nr, auto val)
                                    {
                                      shape[ii++] = Trans(mip.GetJacobianInverse())*Vec<3,T>(0,val,0);
                                    }));
        
	    LegendrePolynomial::Eval(order_inner[0]-2,z,
				 SBLambda([&] (size_t nr, auto val)
					  {
					    shape[ii++] = Trans(mip.GetJacobianInverse())*Vec<3,T>(0,0,val);
					  }));
	  }
      }
  }







  template<typename Tx, typename TFA>  
  void HDivHighOrderFE_Shape<ET_PRISM> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x, y = ip.y, z = ip.z;

    Tx lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    Tx muz[6]  = { 1-z, 1-z, 1-z, z, z, z };
       
    const FACE * faces = ElementTopology::GetFaces (ET_PRISM).Data(); 

    ArrayMem<Tx,20> adpolxy1(order+4),adpolxy2(order+4); 
    ArrayMem<Tx,20> adpolz(order+4);   
    ArrayMem<Tx,10> adpol1(order), adpol2(order), adpol3(order);
    
    // trig faces

    int ii = 5;
    for (int i = 0; i < 2; i++)
      {
	int p = order_facet[i][0];
	int fav[3] = {faces[i][0], faces[i][1], faces[i][2]};

	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	  	
	
	shape[i] = wDu_Cross_Dv (lami[fav[0]], lami[fav[1]], muz[fav[0]]);

	Tx xi = lami[fav[1]]-lami[fav[0]];
	Tx eta = lami[fav[2]];

        Tx sum = lami[fav[1]]+lami[fav[0]];
        Tx bub = lami[fav[1]]*lami[fav[0]];
        // AutoDiff<3> zeta = lami[fop];  


        /*
	T_TRIGFACESHAPES::CalcSplitted(p+2,xi,eta,adpol1,adpol2); 
	for (int k = 0; k < adpol1.Size(); k++)
          adpol1[k] *= 0.5; 
	
        // Curl (Type 2) 2*grad v x grad u
        for (int j = 0; j <= p-1; j++) 
          for (int k = 0; k <= p-1-j; k++)
	    shape[ii++] = Du_Cross_Dv<3> (adpol2[k]*muz[fav[2]], adpol1[j]);

        // Curl (Type 3) //curl( * v) = nabla v x ned + curl(ned)*v
        for (int j = 0; j <= p-1; j++)
	  shape[ii++] = curl_uDvw_minus_Duvw<3> (lami[fav[0]], lami[fav[1]], adpol2[j]*muz[fav[2]]);
        */

        IntLegNoBubble::EvalScaledMult (p-1, xi, sum, bub, adpol1); 
        
        // Typ 1
        for (int k = 0; k <= p-1; k++)
          {
            IntegratedJacobiPolynomialAlpha jac(2*k+3);
            // jac.EvalScaledMult(p-1, eta-sum, 1-zeta, eta, adpol2);
            jac.EvalMult(p-1, eta-sum, eta, adpol2);
            for (int l = 0; l <= p-1-k; l++)
              shape[ii++] = Du_Cross_Dv (adpol2[l]*muz[fav[2]], adpol1[k]);
          }
        
        IntegratedJacobiPolynomialAlpha jac(3);
        // jac.EvalScaledMult(p-1, eta-sum, 1-zeta, eta, adpol2);
        jac.EvalMult(p-1, eta-sum, eta, adpol2);
        for (int j = 0; j <= p-1; j++)
          shape[ii++] = curl_uDvw_minus_Duvw (lami[fav[0]], lami[fav[1]], adpol2[j]*muz[fav[2]]);
      }    
    
    // quad faces
    for (int i = 2; i < 5; i++)
      {
	IVec<2> p = order_facet[i];
	 
	int fmax = 0;
	for (int j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] < vnums[faces[i][fmax]]) fmax = j;

	int fz = 3-fmax; 
	int ftrig = fmax^1; 
	int f = faces[i][fmax];
	int f1 = faces[i][ftrig];
	int f2 = faces[i][fz];

	Tx xi = lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
	Tx eta = 1-lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
	Tx zeta = muz[faces[i][fmax]]-muz[faces[i][fz]]; 
	
	int pp = int(max2(p[0],p[1]))+1;
        /*
	T_ORTHOPOL::CalcTrigExt(pp,xi,eta,adpolxy1); 
	T_ORTHOPOL::Calc(pp,zeta,adpolz); 
        */

	Tx bub = lami[faces[i][fmax]]*lami[faces[i][ftrig]]; 
        IntLegNoBubble::EvalScaledMult (pp-1, xi, 1-eta, 4*bub, adpolxy1); 
        IntLegNoBubble::EvalMult (pp-1, zeta, 1-zeta*zeta, adpolz); 

	double fac = (vnums[faces[i][fz]] > vnums[faces[i][ftrig]]) ? 1 : -1;

	shape[i] = uDvDw_minus_DuvDw (lami[faces[i][fmax]],
                                      lami[faces[i][ftrig]], -0.5*fac*zeta);
					    

	if (vnums[f1] < vnums[f2])
	  {
	    for (int k = 0; k <= p[0]-1; k++)
	      for (int j = 0; j <= p[1]-1; j++, ii++)
		shape[ii] = Du_Cross_Dv (adpolxy1[k], -2*adpolz[j]);
	  }
	else
	  {
	    for (int j = 0; j <= p[1]-1; j++)
	      for (int k = 0; k <= p[0]-1; k++, ii++)
		shape[ii] = Du_Cross_Dv (adpolxy1[k],  2*adpolz[j]);
	  }
	  

        if (vnums[f1] < vnums[f2])
          {
            for (int j= 0; j <= p[0]-1; j++, ii++)
	      shape[ii] = Du_Cross_Dv (adpolxy1[j], zeta);
            for(int j=0; j<= p[1]-1;j++,ii++)
	      shape[ii] = curl_uDvw_minus_Duvw (lami[f1], lami[f], 2*adpolz[j]);
          }  
        else
          {
            for(int j = 0; j <= p[0]-1; j++,ii++)
	      shape[ii] = curl_uDvw_minus_Duvw (lami[f1], lami[f], 2*adpolz[j]);
            for (int j= 0; j <= p[1]-1; j++, ii++)
	      shape[ii] = Du_Cross_Dv (adpolxy1[j], zeta);
          }   
      }    

    int p = order_inner[0];
    int pz = order_inner[2];
    if(p >= 1 && pz >= 1)
      {
	T_TRIGFACESHAPES::CalcSplitted(p+2,x-y,1-x-y,adpolxy1,adpolxy2);
	T_ORTHOPOL::Calc(pz+2,2*z-1,adpolz); 
	

	for(int i=0;i<=p-1;i++)
	  for(int j=0;j<=p-1-i;j++)
	    for(int k=0;k<=pz-1;k++)
	      shape[ii++] = Du_Cross_Dv (adpolxy1[i],adpolxy2[j]*adpolz[k]);

	for(int i=0;i<=p-1;i++)
	  for(int j=0;j<=p-1-i;j++)
	    for(int k=0;k<=pz-1;k++)
	      shape[ii++] = curl_uDvw_minus_Duvw (adpolxy1[i],adpolxy2[j],adpolz[k]);

	for(int j=0;j<=p-1;j++) 
	  for (int k=0;k<=pz-1;k++) 
            shape[ii++] = curl_uDvw_minus_Duvw (x,y, adpolxy2[j]*adpolz[k]);

    	for(int i = 0; i <= p-1; i++) 
	  for(int j = 0; j <= p-1-i; j++) 
            shape[ii++] = curl_uDvw_minus_Duvw (z,1-z, adpolxy1[i]*adpolxy2[j]);

        if (!ho_div_free)
          {  // not yet verified
            ScaledLegendrePolynomial (p, x-y, x+y, adpolxy1);
            LegendrePolynomial (p, 1-2*x, adpolxy2);
            LegendrePolynomial (pz, 1-2*z, adpolz);

            /*
            for (int i = 0; i <= p; i++)
              for (int j = 0; j <= p-i; j++)
                for (int k = 0; k <= pz; k++)
                  if (i+j+k > 0)
                    shape[ii++] = wDu_Cross_Dv<3> ((x-y)*adpolxy1[i], x*adpolxy2[j], z*(1-z)*adpolz[k]);
            */

            for (int i = 0; i <= p; i++)
              for (int j = 0; j <= p-i; j++)
                for (int k = 0; k < pz; k++)
                  shape[ii++] = wDu_Cross_Dv ((x-y)*adpolxy1[i], x*adpolxy2[j], z*(1-z)*adpolz[k]);

            for (int i = 0; i < p; i++)
              for (int j = 0; j < p-i; j++)
                shape[ii++] = wDu_Cross_Dv (z, x*y*adpolxy1[i], (1-x-y)*adpolxy2[j]);

            for (int i = 0; i < p; i++)
              shape[ii++] = wDu_Cross_Dv (z, x, y*(1-x-y)*adpolxy1[i]);


            /*
            for (int i = 0; i <= p-1; i++)
              for (int k = 0; k <= pz; k++)
                shape[ii++] = wDu_Cross_Dv<3> (x*y*adpolxy1[i], z*adpolz[k],  1-x-y);
            */
          }
      }

    if (ii != ndof) cout << "hdiv-prism: dofs missing, ndof = " << ndof << ", ii = " << ii << endl;
  }




  /*
    Pyramid, rational shape functions (J. Schoeberl, 2026)

    Space: Nigam-Phillips H(div) pyramid space R^(2)_{p+1}  (M2AN 46 (2012)), in the
    collapsed-hex coordinates xt = x/(1-z), yt = y/(1-z), lam = 1-z:
        V_p = curl( F1 grad(xt) ) + curl( F2 grad(yt) ) + H grad(xt) x grad(yt),
        F1 = span{ xt^a yt^b lam^c : a<=c-2, b<=c-1, c<=p+2 },  F2 symmetric,
        H  = span{ xt^a yt^b lam^c : max(a,b)<=c-3, c<=p+3 },
    enlarged by 2p div-free interior functions curl( z x^(j+1) (1-x-z) grad(yt) ), curl( z y^(j+1) (1-y-z) grad(xt) ),
    which are needed so that the triangular face functions can mirror the tet face functions
    for every vertex sorting.
    Properties: contains P_p^3, normal traces P_p on trig faces / Q_p on the quad face (matching tet, hex),
    div maps onto the L2 pyramid space (order p-1, or p for RT).
  */
  template<typename Tx, typename TFA>  
  void HDivHighOrderFE_Shape<ET_PYRAMID> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x, y = ip.y, z = ip.z;
    z.Value() = z.Value()*(1-1e-12);

    Tx lam = 1-z;
    Tx xt = x/lam;
    Tx yt = y/lam;

    Tx sigma2d[4]  = { (1-xt)+(1-yt), xt+(1-yt), xt+yt, (1-xt)+yt };
    Tx lambda2d[4] = { (1-xt)*(1-yt), xt*(1-yt), xt*yt, (1-xt)*yt };

    Tx lampow[20];   // lam^k
    lampow[0] = Tx(1.0);
    for (int k = 1; k < 20; k++) lampow[k] = lampow[k-1]*lam;

    int p = order_inner[0];
    int pd = RT ? p : p-1;     // order of divergence space
    
    size_t ii = 5;
    if (only_ho_div) ii = 0;

    ArrayMem<Tx,20> Lxi(order+3), Leta(order+3), Pxi(order+3), Peta(order+3);
    IntegratedLegendrePolynomial::Eval (order+2, 2*xt-1, Lxi);
    IntegratedLegendrePolynomial::Eval (order+2, 2*yt-1, Leta);
    LegendrePolynomial::Eval (order+2, 2*xt-1, Pxi);
    LegendrePolynomial::Eval (order+2, 2*yt-1, Peta);
    
    if (!only_ho_div)
      {
        const FACE * faces = ElementTopology::GetFaces (ET_PYRAMID).Data();

        // trig faces
        for (int i = 0; i < 4; i++)
          {
            int pf = order_facet[i][0];
            int f0 = faces[i][0], f1 = faces[i][1];
            
            Tx L = lambda2d[f0] + lambda2d[f1];    // 1 on the face, 0 on the opposite trig face
            Tx bary[5];
            bary[f0] = (sigma2d[f0]-L)*lam;
            bary[f1] = (sigma2d[f1]-L)*lam;
            bary[4] = z;
            
            // sort vertices, first edge opposite minimal vertex
            int fav[3] = { f0, f1, 4 };
            double sgn = 1;
            if(vnums[fav[0]] > vnums[fav[1]]) { swap(fav[0],fav[1]); sgn = -sgn; }
            if(vnums[fav[1]] > vnums[fav[2]]) { swap(fav[1],fav[2]); sgn = -sgn; }
            if(vnums[fav[0]] > vnums[fav[1]]) { swap(fav[0],fav[1]); sgn = -sgn; }

            Tx ba = bary[fav[0]], bb = bary[fav[1]], bc = bary[fav[2]];

            // RT_0:   sgn * ( L n_i - z/2 (xt, yt, -1) ),  n_i outward direction of the base edge 
            // in hex coordinates: 
            Tx ha(0.0), hb(0.0);
            switch (i)
              {
              case 0: hb = -L*lam; break;
              case 1: ha =  L*lam; break;
              case 2: hb =  L*lam; break;
              case 3: ha = -L*lam; break;
              }
            shape[i] = aDvDw_bDwDu_cDuDv (sgn*ha, sgn*hb, sgn*0.5*z*lam*lam, xt, yt, z);

            if (pf < 1) continue;

            // same as tet, with linear extension of the face barycentrics, times L
            Tx xi = bb-ba;
            Tx sum = bb+ba;
            Tx bub = bb*ba;
            Tx eta = bc;

            // Typ 1
            IntLegNoBubble::
              EvalScaledMult (pf-1, xi, sum, bub, 
                              SBLambda([&](int k, Tx polk) LAMBDA_INLINE
                                       {
                                         IntegratedJacobiPolynomialAlpha jac(2*k+3);
                                         jac.EvalMult(pf-k-1, eta-sum, eta, 
                                                      SBLambda ([&](int j, Tx val) LAMBDA_INLINE
                                                                {
                                                                  shape[ii++] = Du_Cross_Dv (val, L*polk);
                                                                }));
                                       }));
            // Typ 2
            IntegratedJacobiPolynomialAlpha jac(3);
            jac.EvalMult(pf-1, eta-sum, eta, 
                         SBLambda ([&](int j, Tx val) LAMBDA_INLINE
                                   {
                                     shape[ii++] = curl_uDvw_minus_Duvw (ba, bb, L*val);
                                   }));
          }

        // quad face, same as hex, with lam_f replaced by powers of lam
        {
          IVec<2> pf = order_facet[4];
          IVec<4> f = GetFaceSort (4, vnums);	  
          Tx xi  = sigma2d[f[0]]-sigma2d[f[1]];
          Tx eta = sigma2d[f[0]]-sigma2d[f[3]];

          shape[4] = wDu_Cross_Dv (eta, xi, -0.25*lampow[3]);

          ArrayMem<Tx, 20> L_xi(order+2), L_eta(order+2);
          IntegratedLegendrePolynomial::Eval (pf[0]+1, xi, L_xi);
          IntegratedLegendrePolynomial::Eval (pf[1]+1, eta, L_eta);

          for (int k = 0; k < pf[0]; k++)
            for (int l = 0; l < pf[1]; l++, ii++)
              shape[ii] = curl_uDvw_minus_Duvw (L_xi[k+2], L_eta[l+2], -lampow[max2(k,l)+3]);
          
          for (int k = 0; k < pf[0]; k++)
            shape[ii++] = Du_Cross_Dv (L_xi[k+2]*lampow[k+3], eta);
          
          for (int k = 0; k < pf[1]; k++)
            shape[ii++] = Du_Cross_Dv (L_eta[k+2]*lampow[k+3], xi);
        }

        // inner, div-free
        if (p >= 1)
          {
            // curl ( lam^c grad (Lxi Leta) )
            for (int r = 0; r < p; r++)
              for (int s = 0; s < p; s++)
                shape[ii++] = Du_Cross_Dv (lampow[max2(r,s)+3], Lxi[r+2]*Leta[s+2]);

            // curl ( z lam^c P_A(xt) L_{b+2}(yt) grad xt ),  and symmetric
            for (int c = 2; c <= p+1; c++)
              for (int b = 0; b <= c-2; b++)
                for (int a = 0; a <= c-2; a++)
                  if (b <= c-3 || a == 0)
                    {
                      shape[ii++] = Du_Cross_Dv (z*lampow[c]*Pxi[a]*Leta[b+2], xt);
                      shape[ii++] = Du_Cross_Dv (z*lampow[c]*Peta[a]*Lxi[b+2], yt);
                    }
          }
      }

    // inner, non-div-free: divergence spans L2 pyramid space of order pd
    if (!ho_div_free && pd >= 1)
      {
        // base trace P_a(xt) P_b(yt), b >= 1, compensated by curl( P_a L_{b+1} lam^(M+2) grad xt)
        for (int a = 0; a <= pd; a++)
          for (int b = 1; b <= pd; b++)
            {
              int M = max2(a,b);
              shape[ii++] = aDvDw_bDwDu_cDuDv (Tx(0.0), (M+2)*Pxi[a]*Leta[b+1]*lampow[M+1],
                                               2*z*lampow[M+2]*Pxi[a]*Peta[b], xt, yt, z);
            }
        // base trace P_a(xt), a >= 1
        for (int a = 1; a <= pd; a++)
          shape[ii++] = aDvDw_bDwDu_cDuDv ((a+2)*Lxi[a+1]*lampow[a+1], Tx(0.0),
                                           2*z*lampow[a+2]*Pxi[a], xt, yt, z);
        // zero base trace
        for (int a = 0; a <= pd; a++)
          for (int b = 0; b <= pd; b++)
            {
              int M = max2(a,b);
              for (int m = 0; m < pd-M; m++)
                shape[ii++] = wDu_Cross_Dv (xt, yt, Pxi[a]*Peta[b]*z*lampow[M+3+m]);
            }
      }
    
    if (ii != ndof) cout << "hdiv-pyramid: dofs missing, ndof = " << ndof << ", ii = " << ii << endl;
  }

  template<typename Tx, typename TFA>  
  void HDivHighOrderFE_Shape<ET_HEX> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x, y = ip.y, z = ip.z;

    Tx lami[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
			 (1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z};
    Tx sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
			  (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z};

    ArrayMem<Tx, 20> pol_xi(order+2),pol_eta(order+2),pol_zeta(order+2);

    int ii = 6;

    //Faces
    const FACE * faces = ElementTopology::GetFaces (ET_HEX).Data();
    for (int i = 0; i < 6; i++)
      {
	IVec<2> p = order_facet[i];

	Tx lam_f(0);
	for (int j = 0; j < 4; j++)
	  lam_f += lami[faces[i][j]];

        IVec<4> f = GetFaceSort (i, vnums);	  
        Tx xi  = sigma[f[0]]-sigma[f[1]];
        Tx eta = sigma[f[0]]-sigma[f[3]];

        shape[i] = wDu_Cross_Dv(eta, xi, -0.25*lam_f);

        ArrayMem<Tx, 20> L_xi(order+2),L_eta(order+2);

        IntegratedLegendrePolynomial::Eval(p[0]+1,xi,L_xi);
        IntegratedLegendrePolynomial::Eval(p[1]+1,eta,L_eta);

          // if (p[0] >= 1)
          //    IntLegNoBubble::EvalMult(p[0]-1, xi, 1-xi*xi, pol_xi);
          // if (p[1] >= 1)
          //    IntLegNoBubble::EvalMult(p[1]-1,eta, 1-eta*eta, pol_eta);

        for (int k = 0; k < p[0]; k++)
          for (int l = 0; l < p[1]; l++, ii++)
            //shape[ii] = wDu_Cross_Dv<3> (pol_eta[l], pol_xi[k], 2*lam_f);
            shape[ii] = curl_uDvw_minus_Duvw(L_xi[k+2],L_eta[l+2],-lam_f); //divfree
        
        for (int k = 0; k < p[0]; k++)
          //shape[ii++] = wDu_Cross_Dv<3> (-eta, pol_xi[k], lam_f);
          shape[ii++] = Du_Cross_Dv(L_xi[k+2]*lam_f,eta); //divfree

        for (int k = 0; k < p[1]; k++)
          //shape[ii++] = wDu_Cross_Dv<3> (-xi, pol_eta[k], lam_f);
          shape[ii++] = Du_Cross_Dv(L_eta[k+2]*lam_f,xi); //divfree
      }

    auto p = order_inner;

    LegendrePolynomial::Eval(p[0], 2*x-1, pol_xi);
    LegendrePolynomial::Eval(p[1], 2*y-1, pol_eta);
    LegendrePolynomial::Eval(p[2], 2*z-1, pol_zeta);
   
    // ----------------------- hodivfree-SZ-basis -----------------------
    Tx xi  = 2*x-1;
    Tx eta = 2*y-1;
    Tx zeta = 2*z-1;

    ArrayMem<Tx, 20> L_xi(order+2),L_eta(order+2),L_zeta(order+2);

    IntegratedLegendrePolynomial::Eval(p[0]+1,xi,L_xi);
    IntegratedLegendrePolynomial::Eval(p[1]+1,eta,L_eta);
    IntegratedLegendrePolynomial::Eval(p[2]+1,zeta,L_zeta); 

    // Type 1 (divergence-free):
    for (int i = 0; i < p[0]; i++)
      for (int j = 0; j < p[1]; j++)
        for (int k = 0; k < p[2]; k++)
          shape[ii++] = wDu_Cross_Dv (L_xi[i+2]*L_zeta[k+2], y, -4*pol_eta[j+1]);

    for (int i = 0; i < p[0]; i++)
      for (int j = 0; j < p[1]; j++)
        for (int k = 0; k < p[2]; k++)
          shape[ii++] = wDu_Cross_Dv (L_eta[j+2]*L_zeta[k+2], x, 4*pol_xi[i+1]);

    for (int j = 0; j < p[1]; j++)
      for (int k = 0; k < p[2]; k++)
        shape[ii++] = Du_Cross_Dv (L_eta[j+2]*L_zeta[k+2], 2*xi);

    for (int i = 0; i < p[0]; i++)
      for (int k = 0; k < p[2]; k++)
        shape[ii++] = Du_Cross_Dv (L_xi[i+2]*L_zeta[k+2], 2*eta);

    for (int i = 0; i < p[0]; i++)
      for (int j = 0; j < p[1]; j++)
        shape[ii++] = Du_Cross_Dv (L_xi[i+2]*L_eta[j+2], 2*zeta);
    

    if (!ho_div_free)
    // ----------------------- remaining-SZ-basis -----------------------
    {
      // Type 2:
      for (int i = 0; i < p[0]; i++)
        for (int j = 0; j < p[1]; j++)
          for (int k = 0; k < p[2]; k++)
            shape[ii++] = uDvDw_minus_DuvDw (L_xi[i+2], L_eta[j+2], L_zeta[k+2]);                                 

      for (int j = 0; j < p[1]; j++)
        for (int k = 0; k < p[2]; k++)
          shape[ii++] = uDvDw_minus_DuvDw (L_eta[j+2], L_zeta[k+2], x);     

      for (int i = 0; i < p[0]; i++)
        for (int k = 0; k < p[2]; k++)
          shape[ii++] = uDvDw_minus_DuvDw (L_zeta[k+2], L_xi[i+2], y);          

      for (int i = 0; i < p[0]; i++)
        for (int j = 0; j < p[1]; j++)
          shape[ii++] = uDvDw_minus_DuvDw (L_xi[i+2], L_eta[j+2], z);  

      // Type 3:
      for (int i = 0; i < p[0]; i++)
        shape[ii++] = wDu_Cross_Dv (y, z, L_xi[i+2]);

      for (int j = 0; j < p[1]; j++)
        shape[ii++] = wDu_Cross_Dv (x, -z, L_eta[j+2]); 
        
      for (int k = 0; k < p[2]; k++)
        shape[ii++] = wDu_Cross_Dv (x, y, L_zeta[k+2]);        
    }

      /*// ----------------------- full-non-SZ-basis -----------------------

      for (int i = 0; i <= p[0]; i++)
        for (int j = 0; j <= p[1]; j++)
          for (int k = 0; k < p[2]; k++)
            shape[ii++] = wDu_Cross_Dv (x, y, pol_xi[i]*pol_eta[j]*pol_zeta[k]*z*(1-z));
      for (int i = 0; i <= p[0]; i++)
        for (int j = 0; j < p[1]; j++)
          for (int k = 0; k <= p[2]; k++)
            shape[ii++] = wDu_Cross_Dv (x, z, pol_xi[i]*pol_eta[j]*pol_zeta[k]*y*(1-y));
      for (int i = 0; i < p[0]; i++)
        for (int j = 0; j <= p[1]; j++)
          for (int k = 0; k <= p[2]; k++)
            shape[ii++] = wDu_Cross_Dv (y, z, pol_xi[i]*pol_eta[j]*pol_zeta[k]*x*(1-x));

      */

  }

    
  template<>
  void HDivHighOrderFE<ET_TRIG> :: 
  CalcNormalShape (const IntegrationPoint & ip, 
                   BareSliceVector<> nshape) const
  {
    // Vector<> nshape1(nshape.Size());
    // HDivFiniteElement<2>::CalcNormalShape (ip, nshape1);

    int fnr = ip.FacetNr();
    double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };

    IVec<2> e0 = ET_trait<ET_TRIG>::GetEdge (fnr);
    double fac = vnums[e0[0]] > vnums[e0[1]] ? 1 : -1;
    IVec<2> e = ET_trait<ET_TRIG>::GetEdgeSort (fnr, vnums);
    AutoDiff<1> xi (lam[e[1]]-lam[e[0]], 0);
    
    ArrayMem<AutoDiff<1>,10> adpol1(order);
    
    nshape[0] = fac*xi.DValue(0);

    int p = order_inner[0]; 
    IntLegNoBubble:: EvalMult (p-1, xi, 0.25*(1-xi*xi), adpol1);
    
    for(int j = 0; j < p; j++) 	      
      nshape[j+1] = -2*fac*adpol1[j].DValue(0);
    /*
    cout << "nshape = " << nshape << endl;
    cout << "nshape1 = " << nshape1 << endl;
    cout << "************************************ diff = " << L2Norm (nshape-nshape1) << endl;
    */
  }

  
  template<>
  void HDivHighOrderFE<ET_QUAD> :: 
  CalcNormalShape (const IntegrationPoint & ip, 
                   BareSliceVector<> nshape) const
  {
    // Vector<> nshape1(nshape.Size());
    // HDivFiniteElement<2>::CalcNormalShape (ip, nshape1);

    int fnr = ip.FacetNr();
    double lam[] = { (1-ip(0))*(1-ip(1)), ip(0)*(1-ip(1)), ip(0)*ip(1),(1-ip(0))*ip(1) };

    IVec<2> e0 = ET_trait<ET_QUAD>::GetEdge (fnr);
    double fac = vnums[e0[0]] > vnums[e0[1]] ? 1 : -1;
    IVec<2> e = ET_trait<ET_QUAD>::GetEdgeSort (fnr, vnums);
    AutoDiff<1> xi (lam[e[1]]-lam[e[0]], 0);
    
    ArrayMem<AutoDiff<1>,10> adpol1(order);
    
    nshape[0] = fac*xi.DValue(0);

    int p = order_inner[0]; 
    IntLegNoBubble:: EvalMult (p-1, xi, 0.25*(1-xi*xi), adpol1);
    
    for(int j = 0; j < p; j++) 	      
      nshape[j+1] = -2*fac*adpol1[j].DValue(0);
    /*
    cout << "nshape = " << nshape << endl;
    cout << "nshape1 = " << nshape1 << endl;
    cout << "************************************ diff = " << L2Norm (nshape-nshape1) << endl;
    */
  }



  template<>
  void HDivHighOrderFE<ET_TET> :: 
  CalcNormalShape (const IntegrationPoint & ip, 
                   BareSliceVector<> nshape) const
  {
    // Vector<> nshape1(nshape.Size());
    // HDivFiniteElement<3>::CalcNormalShape (ip, nshape1);

    int fnr = ip.FacetNr();
    double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
    
    IVec<4> face = ET_trait<ET_TET>::GetFace(fnr);
    
    IntegrationPoint ip2d(lam[face[0]], lam[face[1]], lam[face[2]]);
    ArrayMem<int,3> vnumsf(3);
    for (int j = 0; j < 3; j++) vnumsf[j] = vnums[face[j]];
    
    HDivHighOrderNormalTrig<> trig(order_facet[fnr][0]);
    trig.SetVertexNumbers (vnumsf);

    VectorMem<20> tmp(ndof); // nshape.Size());
    trig.CalcShape (ip2d, tmp);
    nshape = -tmp;
    // cout << "nshape1 = " << endl << nshape1 << endl;
    // cout << "nshape = " << endl << nshape << endl;
    // cout << "******************************************** diff = " << L2Norm (nshape-nshape1) << endl;
  }


  template <ELEMENT_TYPE ET>
  void HDivHighOrderFE<ET> ::
  CalcNormalShape(const IntegrationPoint & ip,
                  BareSliceVector<> nshape) const
  {
    cout << "HDivHOFE, calcnormalshape not overloaded" << endl;
  }
  

  template <ELEMENT_TYPE ET>
  void HDivHighOrderFE<ET> ::
  CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const
  {
    shape.AddSize(this->ndof, bmip.DimSpace()) = 0.0;
    Switch<4-DIM>
      (bmip.DimSpace()-DIM,[this,&bmip,shape](auto CODIM)
       {
         auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);

         static_cast<const HDivHighOrderFE_Shape<ET>*> (this)
           -> CalcDualShape2 (mip, SBLambda([shape] (size_t i, auto val)
                                            {
                                              shape.Row(i) = val;
                                            }));
       });
  }



}



#endif
