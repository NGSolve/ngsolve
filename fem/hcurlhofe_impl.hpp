#ifndef FILE_HCURLHOFE_IMPL
#define FILE_HCURLHOFE_IMPL

/*********************************************************************/
/* File:   hcurlhofe.hpp                                             */
/* Author: Sabine Zaglmayr, Joachim Schoeber                         */
/* Date:   20. Maerz 2003                                            */
/*                                                                   */
/* AutoCurl - revision: J. Schoeberl, March 2009                     */
/*********************************************************************/

#include "recursive_pol.hpp"
// #include "thdivfe.hpp"
#include "hcurlhofe.hpp"
#include "hcurlfe_utils.hpp"

namespace ngfem
{

  // declaration of the shapes ...

  

  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET2> class TSHAPES, typename BASE>
  void HCurlHighOrderFE<ET,TSHAPES,BASE> :: ComputeNDof()
  {
    ndof = N_EDGE;

    for (int i = 0; i < N_EDGE; i++)
      if(order_edge[i] > 0)
        ndof += usegrad_edge[i]*order_edge[i];
    
    for(int i = 0; i < N_FACE; i++)
      if (FaceType(i) == ET_TRIG)
        {
          if (order_face[i][0] > 1)
            {
              int p = order_face[i][0];
              int pg = p - (type1 ? 1 : 0);
              ndof += usegrad_face[i]*pg*(pg-1)/2;
              ndof += (p+2)*(p-1)/2;
                // ndof += ((usegrad_face[i]+1)*order_face[i][0]+2)*(order_face[i][0]-1)/2;
            }
        }
      else
        {
          if(order_face[i][0]>=0 && order_face[i][1]>=0)
            ndof +=  (usegrad_face[i]+1)*order_face[i][0]*order_face[i][1] 
              + order_face[i][0] + order_face[i][1]; 
        }

    switch (ET)
      {
      case ET_TET: 
        if(order_cell[0] > 2) {
	  if (type1)	    
	    ndof += usegrad_cell*(order_cell[0]-3)*(order_cell[0]-2)*(order_cell[0]-1)/6  + (order_cell[0]-2)*(order_cell[0]-1)*(2*order_cell[0]+3)/6 ;
	  else 
	    ndof += ((usegrad_cell + 2) * order_cell[0] + 3) 
	      * (order_cell[0]-2) * (order_cell[0]-1) / 6;
	}
        break;
      case ET_PRISM:
        if(order_cell[2] > 0 && order_cell[0] > 1)
          ndof += ((usegrad_cell+2)*order_cell[2] + 1) * order_cell[0]*(order_cell[0]-1)/2
            + (order_cell[0]-1)*order_cell[2]; 
        break;
      case ET_PYRAMID:
        {
          int pc = order_cell[0]; //SZ: no problem to do anisotropic, but for the moment 
          // is it worth getting crazy :-) 
          if(order_cell[0]>1)
            ndof += usegrad_cell*(pc-1)*pc*(2*pc-1)/6 + pc*(2*pc*pc+3*pc-2)/3; 
          break;
        }
      case ET_HEX:
        if(order_cell[0] >= 0 && order_cell[1]>= 0 && order_cell[2]>=0)
          ndof += (usegrad_cell + 2)* order_cell[0] * order_cell[1] * order_cell[2]
            + order_cell[1]*order_cell[2]  + order_cell[0]*(order_cell[1] + order_cell[2]);  
        break;
      default:
        ;
      }
    
    TORDER horder = 0; 
    for (int i = 0; i < N_EDGE; i++)
      horder = max2 (horder, order_edge[i]);

    for(int i=0; i < N_FACE; i++) 
      if (ET_trait<ET>::FaceType(i) == ET_TRIG)
        horder = max2 (horder, order_face[i][0]);
      else
        horder = max2 (horder, Max (order_face[i]));

    if (DIM == 3)
      horder = max2 (horder, Max(order_cell));

    // for integration order .. 
    if (ET == ET_PRISM || ET == ET_HEX || ET == ET_PYRAMID || ET == ET_QUAD)
      horder++;
    else
      if (horder==0) horder++;
    order = horder;
  }



  
  template <ELEMENT_TYPE ET> 
  class HCurlHighOrderFE_Shape :  public HCurlHighOrderFE<ET>
  {
    using ET_trait<ET>::DIM;    
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (TIP<DIM,Tx> ip, TFA & shape) const;

    template <typename MIP, typename TFA>
    inline void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception(string("CalcDualShape missing for HighOrderHCurl element ")+ElementTopology::GetElementName(ET));
    }
  };



  
  //------------------------------------------------------------------------
  // HCurlHighOrderSegm
  //------------------------------------------------------------------------

  
  template<> template<typename Tx, typename TFA>  
  void HCurlHighOrderFE_Shape<ET_SEGM> :: T_CalcShape (TIP<1,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x;
    Tx lam[2] = { x, 1-x };

    IVec<2> e = GetEdgeSort (0, vnums);	  
    
    //Nedelec low order edge shape function 
    shape[0] = uDv_minus_vDu (lam[e[0]], lam[e[1]]);

    int p = order_edge[0]; //order_cell[0]; 
    //HO-Edge shapes (Gradient Fields)   
    if(p > 0 && usegrad_cell)
      { 
        // LegendrePolynomial::
        size_t ii = 1;
        EdgeOrthoPol::
          EvalScaledMult (p-1, 
                          lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
                          lam[e[0]]*lam[e[1]], 
                          SBLambda ([&](int nr, Tx val)
                                    {
                                      shape[ii++] = Du (val);
                                    }));
      }
  }
  


  
  //------------------------------------------------------------------------
  // HCurlHighOrderTrig
  //------------------------------------------------------------------------
  
  template<> template<typename Tx, typename TFA>  
  INLINE void HCurlHighOrderFE_Shape<ET_TRIG> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x, y = ip.y;
    Tx lam[3] = { x, y, 1-x-y };

	
    int ii = 3; 
    for (int i = 0; i < 3; i++)
      {
        IVec<2> e = GetEdgeSort (i, vnums);	  

	//Nedelec low order edge shape function 
        shape[i] = uDv_minus_vDu (lam[e[0]], lam[e[1]]);

	int p = order_edge[i]; 
	//HO-Edge shapes (Gradient Fields)   
	if(p > 0 && usegrad_edge[i]) 
	  { 
	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalScaledMult (order_edge[i]-1, 
			      lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
			      lam[e[0]]*lam[e[1]], 
                              // adpol1
                              SBLambda ([&](int nr, Tx val)
                                        {
                                          shape[ii++] = Du (val);
                                        }));
	  }
      }   
     
    //Inner shapes (Face) 
    int p = order_face[0][0];      
    if(p > 1) 
      {
	IVec<4> fav = GetFaceSort (0, vnums);

	// Tx xi  = lam[fav[2]]-lam[fav[1]];
	// Tx eta = lam[fav[0]]; 
	int pg = p - 2 - (type1 ? 1 : 0);

	// gradients:	
	if (usegrad_face[0] && pg >=0)
            
	  // DubinerBasis
          TrigOrthoPolGrad::EvalMult
	    (pg, lam[fav[0]], lam[fav[1]], 
	     lam[fav[0]]*lam[fav[1]]*lam[fav[2]], 
	     SBLambda
	     ([&](int nr, Tx val)
	      {
		shape[ii++] = Du (val);
	      }));

	// now ready: compatibility with tets/prisms ! non-gradient inner shapes, same for type 1 and 2
        if (true) // if (type1)
          {
            DubinerBasis::EvalMult
              (p-2, lam[fav[0]], lam[fav[1]], 
               lam[fav[0]], 
               SBLambda
               ([&](int nr, Tx val)
                {
                  shape[ii++] = wuDv_minus_wvDu (lam[fav[1]], lam[fav[2]], val);
                }));
            
            LegendrePolynomial::EvalMult 
              (p-2, lam[fav[2]]-lam[fav[1]], lam[fav[2]], 
               SBLambda
               ([&] (int j, Tx val)
                {
                  shape[ii++] = wuDv_minus_wvDu (lam[fav[1]], lam[fav[0]], val);
                }));
          }
        else
          {
            Tx xi  = lam[fav[2]]-lam[fav[1]];
            Tx eta = lam[fav[0]];
            
            ArrayMem<Tx,20> adpol1(order),adpol2(order);	
            TrigShapesInnerLegendre::CalcSplitted(p+1, xi, eta, adpol1,adpol2);
	
            // other combination
            for (int j = 0; j < p-1; j++)
              for (int k = 0; k < p-1-j; k++, ii++)
                shape[ii] = uDv_minus_vDu (adpol2[k], adpol1[j]);     
            
            // rec_pol * Nedelec0 
            for (int j = 0; j < p-1; j++, ii++)
              shape[ii] = wuDv_minus_wvDu (lam[fav[1]], lam[fav[2]], adpol2[j]);
          }
      }
  }


  //------------------------------------------------------------------------
  // HCurlHighOrderQuad
  //------------------------------------------------------------------------
  

  template<> template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE_Shape<ET_QUAD> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x, y = ip.y;
    Tx hx[2] = { x, y };
    Tx lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  

    int ii = 4;
    ArrayMem<Tx, 10> pol_xi(order+2), pol_eta(order+2);

    for (int i = 0; i < 4; i++)
      {
	// int p = order_edge[i]; 
        IVec<2> e = GetEdgeSort (i, vnums);	  
	
	Tx xi  = sigma[e[1]]-sigma[e[0]];
	Tx lam_e = lami[e[0]]+lami[e[1]];  
        Tx bub = 0.25 * lam_e * (1 - xi*xi);

	// Nedelec0-shapes
	shape[i] = uDv (0.5 * lam_e, xi); 

	// High Order edges ... Gradient fields 
	if(usegrad_edge[i])
	  {
            /*
	    LegendrePolynomial::
	      EvalMult (order_edge[i]-1, 
			xi, bub, pol_xi);

	    for (int j = 0; j < p; j++)
              shape[ii++] = Du<2> (pol_xi[j]);
            */
	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalMult (order_edge[i]-1, 
			xi, bub, SBLambda ([&](int i, Tx val)
                                           {
                                             shape[ii++] = Du (val);
                                           }));

	  }
      }



     
    IVec<2> p = order_face[0]; // (order_cell[0],order_cell[1]);


    if (usegrad_face[0] && p[0] >= 1 && p[1] >= 1)
      {
        Vec<2,Tx> xi = ET_trait<ET_QUAD>::XiFace(0, hx, vnums);
	Tx bub = 1.0/16 * (1-xi(0)*xi(0))*(1-xi(1)*xi(1));
        
        QuadOrthoPol::EvalMult1Assign(p[0]-1, xi(0), bub,
              SBLambda ([&](int i, Tx val) LAMBDA_INLINE 
                    {  
                      QuadOrthoPol::EvalMult (p[1]-1, xi(1), val, 
                                                    SBLambda([&](int i2, Tx v2)
                                                             {
                                                               shape[ii++] = Du (v2);
                                                             }));
                    }));
      }





    int fmax = 0; 
    for (int j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax])
	fmax = j;
    
    int f1 = (fmax+3)%4; 
    int f2 = (fmax+1)%4; 
    if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f2 > f1; 

    Tx xi = sigma[fmax]-sigma[f1];  // in [-1,1]
    Tx eta = sigma[fmax]-sigma[f2]; // in [-1,1]
    
    T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
    T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);

    /*
    //Gradient fields 
    if(usegrad_face[0])
      for (int k = 0; k < p[0]; k++)
	for (int j= 0; j < p[1]; j++)
          shape[ii++] = Du<2> (pol_xi[k]*pol_eta[j]);
    */

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


  
  //------------------------------------------------------------------------
  //        Tetrahedron
  //------------------------------------------------------------------------
 

  template<> template<typename Tx, typename TFA>  
  void HCurlHighOrderFE_Shape<ET_TET> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x, y = ip.y, z = ip.z;    
    Tx lam[4] = { x, y, z, 1-x-y-z };

    ArrayMem<Tx,20> adpol1(order+2),adpol2(order+2),adpol3(order+2); 
    int ii = 6; 

    for (int i = 0; i < N_EDGE; i++)
      { 
	int p = order_edge[i]; 
        IVec<2> e = GetEdgeSort (i, vnums);	  
	
	//Nedelec low order edge shape function 
        shape[i] = uDv_minus_vDu (lam[e[0]], lam[e[1]]);

	//HO-Edge shape functions (Gradient Fields) 	
	if (p > 0 && usegrad_edge[i]) 
	  {     
	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalScaledMult (p-1, 
			      lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
			      lam[e[0]]*lam[e[1]], 
                              SBLambda ([&](int i, Tx val)
                                        {
                                          shape[ii++] = Du (val);
                                        }));
	  }
      }

    // face shape functions
    for(int i = 0; i < N_FACE; i++) 
      if (order_face[i][0] >= 2)
        {
          IVec<4> fav = GetFaceSort (i, vnums);
          
          int vop = 6 - fav[0] - fav[1] - fav[2];  	
          int p = order_face[i][0];
          int pg = p - (type1 ? 1 : 0);
	  
          Tx xi = lam[fav[2]]-lam[fav[1]];
          Tx eta = lam[fav[0]]; // lo 
          Tx zeta = lam[vop];   // lz 
          
          TetShapesFaceLegendre::CalcSplitted (p+1, xi, eta, zeta, adpol1, adpol2); 
          
          // gradients 
          if (usegrad_face[i] && pg >= 2)
            
	    // DubinerBasis
            TrigOrthoPolGrad::EvalScaledMult
	      (pg-2, lam[fav[0]], lam[fav[1]],
	       1-lam[vop],
	       lam[fav[0]]*lam[fav[1]]*lam[fav[2]], 
	       SBLambda
	       ([&](int nr, Tx val)
		{
		  shape[ii++] = Du (val);
		}));
            

	  // non-gradient face shapes
          // if (type1) {
	  if (true) {

	    DubinerBasis::EvalMult
	      (p-2, lam[fav[0]], lam[fav[1]],
	       lam[fav[0]],
	       SBLambda
	       ([&](int nr, Tx val)
		{
		  shape[ii++] =	wuDv_minus_wvDu(lam[fav[1]], lam[fav[2]], val);
		}));
	    
            LegendrePolynomial::EvalMult
	      (p-2, lam[fav[2]]-lam[fav[1]], lam[fav[2]], 
	       SBLambda
	       ([&] (int j, Tx val)
		{
		  shape[ii++] = wuDv_minus_wvDu(lam[fav[1]], lam[fav[0]], val);
		}));
	  }
	  else {

	    // other combination
	    for (int j = 0; j <= p-2; j++)
	      for (int k = 0; k <= p-2-j; k++, ii++)
		shape[ii] = uDv_minus_vDu (adpol2[k], adpol1[j]);
          
	    // type 3
	    for (int j = 0; j <= p-2; j++, ii++)
	      shape[ii] = wuDv_minus_wvDu (lam[fav[1]], lam[fav[2]], adpol2[j]);
	  }
        }
    
    
    int p = order_cell[0];
    int pg = p - (type1 ? 1 : 0);

    // gradient inner shapes
	
    if (usegrad_cell && pg >= 3)

      DubinerBasis3DOrthoBub::EvalMult
	(pg-3, lam[0], lam[1], lam[2], lam[0]*lam[1]*lam[2]*lam[3], 
	 SBLambda
	 ([&](int nr, Tx val)
	  {
	    shape[ii++] = Du(val);
	  }));

    // non-gradient inner shapes, same for type 1 and type 2 Nedelec:

    if (p >= 3)  {

      Tx lam21 = lam[2]*lam[1], lam30 = lam[3]*lam[0];

      DubinerBasis3D::Eval
	(p-3, lam[0], lam[1], lam[2],
	 SBLambda
	 ([&](int nr, Tx val)
	  {
	    shape[ii++] = wuDv_minus_wvDu(lam[3], lam[0], lam21*val);
	    shape[ii++] = wuDv_minus_wvDu(lam[2], lam[1], lam30*val);
	  }));


      DubinerBasis::EvalMult
	(p-3, lam[0], lam[1], lam[0]*lam[1], 
	 SBLambda
	 ([&](int nr, Tx val)
	  {
	    shape[ii++] = wuDv_minus_wvDu(lam[3], lam[2], val);
	  }));
    }
  }


		        
  //------------------------------------------------------------------------
  //                   Prism
  //------------------------------------------------------------------------

  template<> template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE_Shape<ET_PRISM> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;

    // Tx x = hx[0], y = hx[1], z = hx[2];
    Tx x = ip.x, y = ip.y, z = ip.z;
    
    Tx lam[6] = { x, y, 1-x-y, x, y, 1-x-y };
    Tx muz[6]  = { 1-z, 1-z, 1-z, z, z, z };

    Tx sigma[6];
    for (int i = 0; i < 6; i++) sigma[i] = lam[i] + muz[i];

    ArrayMem<Tx,20> adpolxy1(order+3),adpolxy2(order+3); 
    ArrayMem<Tx,20> adpolz(order+3);   

    int ii = 9;
    
    // horizontal edge shapes
    for (int i = 0; i < 6; i++)
      {
	int p = order_edge[i]; 
        IVec<2> e = GetEdgeSort (i, vnums);	  
	
	//Nedelec0
        shape[i] = wuDv_minus_wvDu (lam[e[0]], lam[e[1]], muz[e[1]]);
   
	//high order \nabla (P_edge(x,y) * muz)
	if (p > 0 && usegrad_edge[i])
	  {
	    /*
              T_ORTHOPOL::CalcTrigExt(p+1, lam[e[1]]-lam[e[0]],
              1-lam[e[0]]-lam[e[1]],adpolxy1);
              for(int j = 0; j <= p-1; j++)
              shape[ii++] = Du<3> (adpolxy1[j] * muz[e[1]]);
	    */
            Tx xi = lam[e[1]]-lam[e[0]]; 
            Tx eta = lam[e[0]]+lam[e[1]]; 
            Tx bub = lam[e[0]]*lam[e[1]]*muz[e[1]];

	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalScaledMult (p-1, xi, eta, bub, adpolxy1);
	    for(int j = 0; j <= p-1; j++)
              shape[ii++] = Du (adpolxy1[j]);
	  }
      }

    //Vertical Edge Shapes
    for (int i = 6; i < 9; i++)
      {
	int p = order_edge[i]; 
        IVec<2> e = GetEdgeSort (i, vnums);	  

        shape[i] = wuDv_minus_wvDu (muz[e[0]], muz[e[1]], lam[e[1]]);
	
	//high order edges:  \nabla (T_ORTHOPOL^{p+1}(2z-1) * lam(x,y))
	if(p > 0 && usegrad_edge[i])
	  {
	    // T_ORTHOPOL::Calc (p+1, muz[e[1]]-muz[e[0]], adpolz);
	    // for (int j = 0; j < p; j++)
	    //   shape[ii++] = Du<3> (adpolz[j] * lam[e[1]]);

	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalMult (p-1, 
			muz[e[1]]-muz[e[0]], 
			muz[e[0]]*muz[e[1]]*lam[e[1]], adpolz);
	    
	    for (int j = 0; j < p; j++)
              shape[ii++] = Du (adpolz[j]);
	  }
      }


    const FACE * faces = ElementTopology::GetFaces (ET_PRISM); 

    // trig face shapes
    for (int i = 0; i < 2; i++)
      {
	int p = order_face[i][0];
	if (p < 2) continue;

	IVec<4> fav = GetFaceSort (i, vnums);

        {
          // gradients 
          if (usegrad_face[i])
            {
              TrigOrthoPolGrad::EvalMult (p-2, lam[fav[0]], lam[fav[1]], 
                                       lam[fav[0]]*lam[fav[1]]*lam[fav[2]]*muz[fav[2]], 
                                       SBLambda ([&](int nr, Tx val)
                                                 {
                                                   shape[ii++] = Du (val);
                                                 }));
            }
        }

        // if (type1) {        
        if (true) {
          
          DubinerBasis::EvalMult
            (p-2, lam[fav[0]], lam[fav[1]],
             lam[fav[0]]*muz[fav[2]],
             SBLambda
             ([&](int nr, Tx val)
             {
               shape[ii++] =  wuDv_minus_wvDu(lam[fav[1]], lam[fav[2]], val);
             }));
	    
          LegendrePolynomial::EvalMult
            (p-2, lam[fav[2]]-lam[fav[1]], lam[fav[2]]*muz[fav[2]], 
             SBLambda
             ([&] (int j, Tx val)
             {
               shape[ii++] = wuDv_minus_wvDu(lam[fav[1]], lam[fav[0]], val);
             }));
        }

        else

          {
            
            Tx xi = lam[fav[2]]-lam[fav[1]];
            Tx eta = lam[fav[0]]; // 1-lam[f2]-lam[f1];
            
            T_TRIGFACESHAPES::CalcSplitted(p+1,xi,eta,adpolxy1,adpolxy2); 
            /*
              if(usegrad_face[i])
              // gradient-fields =>  \nabla( adpolxy1*adpolxy2*muz )
              for (int j = 0; j <= p-2; j++)
              for (int k = 0; k <= p-2-j; k++)
              shape[ii++] = Du<3> (adpolxy1[j]*adpolxy2[k] * muz[fav[2]]);
        */
            
            
            // rotations of grad-fields => grad(uj)*vk*w -  uj*grad(vk)*w 
            for (int j = 0; j <= p-2; j++)
              for (int k = 0; k <= p-2-j; k++)
                shape[ii++] = wuDv_minus_wvDu (adpolxy2[k], adpolxy1[j], muz[fav[2]]);
            
            //  Ned0*adpolxy2[j]*muz 
            for (int j = 0; j <= p-2; j++,ii++)
              shape[ii] = wuDv_minus_wvDu (lam[fav[1]], lam[fav[2]], adpolxy2[j]*muz[fav[2]]);
          }
      }
    

    // quad faces
    for (int i = 2; i < 5; i++)
      {
	IVec<2> p = order_face[i];
        IVec<4> f = GetFaceSort (i, vnums);	  

        {
          Tx xi  = sigma[f[0]] - sigma[f[1]]; 
          Tx eta = sigma[f[0]] - sigma[f[3]];
        
          Tx scalexi(1.0), scaleeta(1.0);
          if (f[0] / 3 == f[1] / 3)  
            scalexi = lam[f[0]]+lam[f[1]];  // xi is horizontal
          else
            scaleeta = lam[f[0]]+lam[f[3]];
          
          Tx bub = (1.0/16)*(scaleeta*scaleeta-eta*eta)*(scalexi*scalexi-xi*xi);
          QuadOrthoPol::EvalScaled     (p[0]-1, xi, scalexi, adpolxy1);
          QuadOrthoPol::EvalScaledMult (p[1]-1, eta, scaleeta, bub, adpolz);

          
          if(usegrad_face[i])
            {
              // Gradientfields nabla(polxy*polz) 
              for (int k = 0; k <= p[0]-1; k++)
                for (int j = 0; j <= p[1]-1; j++)
                  shape[ii++] = Du (adpolxy1[k] * adpolz[j]);
            }
        }
        

	int fmax = 0;
	for (int j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;

	int fz = 3-fmax; 
	int ftrig = fmax^1; 
	Tx xi = lam[faces[i][fmax]]-lam[faces[i][ftrig]]; 
	Tx eta = 1-lam[faces[i][fmax]]-lam[faces[i][ftrig]]; 
	Tx zeta = muz[faces[i][fmax]]-muz[faces[i][fz]]; 
	
	int pp = int(max2(p[0],p[1]))+1;
	T_ORTHOPOL::CalcTrigExt(pp,xi,eta,adpolxy1); 
	T_ORTHOPOL::Calc(pp,zeta,adpolz); 
        
        
        /*
	if(usegrad_face[i])
	  {
	    // Gradientfields nabla(polxy*polz) 
	    if (vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) 
	      for (int k = 0; k <= p[0]-1; k++)
		for (int j = 0; j <= p[1]-1; j++)
                  shape[ii++] = Du<3> (adpolxy1[k] * adpolz[j]);
	    else
	      for (int j = 0; j <= p[0]-1; j++)
		for (int k = 0; k <= p[1]-1; k++)
                  shape[ii++] = Du<3> (adpolxy1[k] * adpolz[j]);
	  }
        */

	// Rotations of GradFields => nabla(polxy)*polz - polxy*nabla(polz)
	if (vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) 
	  for (int k = 0; k <= p[0]-1; k++)
	    for (int j = 0; j <= p[1]-1; j++)
              shape[ii++] = uDv_minus_vDu (adpolz[j], adpolxy1[k]);
	else
	  for (int j = 0; j <= p[0]-1; j++)
	    for (int k = 0; k <= p[1]-1; k++)
              shape[ii++] = uDv_minus_vDu (adpolxy1[k], adpolz[j]);
	
	// Type 3 
	// (ned0_trig)*polz, (ned0_quad)* polxy 

	if(vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) // p = (p_trig,p_z) 
	  {
	    for(int j=0;j<=p[0]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu (muz[faces[i][fz]], muz[faces[i][fmax]], adpolxy1[j]);
	    for(int j=0;j<=p[1]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu (lam[faces[i][ftrig]], lam[faces[i][fmax]], adpolz[j]);
	  }
	else 
	  {
	    for(int j=0;j<=p[0]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu (lam[faces[i][ftrig]], lam[faces[i][fmax]], adpolz[j]);
	    for(int j=0;j<=p[1]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu (muz[faces[i][fz]], muz[faces[i][fmax]], adpolxy1[j]);
	  }
      }
    
    if(order_cell[0] > 1 && order_cell[2] > 0) 
      {
        IVec<3> p = order_cell[0];
        if (usegrad_cell && p[0] > 1 && p[2] > 0)
          {
            // gradientfields
            int nf = (p[0]-0)*(p[0]-1)/2;
            ArrayMem<Tx,20> pol_trig(nf);
            
            DubinerBasis::EvalMult (p[0]-2, x, y, x*y*(1-x-y),pol_trig);
            LegendrePolynomial::EvalMult (p[2]-1, 2*z-1, z*(1-z), adpolz);
            
            for (int i = 0; i < nf; i++)
              for (int k = 0; k <= p[2]-1; k++)
                shape[ii++] = Du (pol_trig[i] * adpolz[k]);
          }


	T_TRIGFACESHAPES::CalcSplitted(order_cell[0]+1,x-y,1-x-y,adpolxy1,adpolxy2);
	T_ORTHOPOL::Calc(order_cell[2]+1,2*z-1,adpolz); 
	
        /*
	// gradientfields
	if(usegrad_cell)
	  for(int i=0;i<=order_cell[0]-2;i++)
	    for(int j=0;j<=order_cell[0]-2-i;j++)
	      for(int k=0;k<=order_cell[2]-1;k++)
                shape[ii++] = Du<3> (adpolxy1[i]*adpolxy2[j]*adpolz[k]);
        */

	// Rotations of gradientfields
	for(int i=0;i<=order_cell[0]-2;i++)
	  for(int j=0;j<=order_cell[0]-2-i;j++)
	    for(int k=0;k<=order_cell[2]-1;k++)
	      {
                shape[ii++] = wuDv_minus_wvDu (adpolxy1[i],adpolxy2[j],adpolz[k]);
                shape[ii++] = uDv_minus_vDu (adpolxy1[i],adpolxy2[j]*adpolz[k]);
	      }

	// Type 3 
	// ned0(trig) * polxy2[j]*polz 
	// z.DValue(0) * polxy1[i] * polxy2[j] 
	// double ned_trig[2] = {y.Value(),-x.Value()};  
	for(int j=0;j<=order_cell[0]-2;j++) 
	  for (int k=0;k<=order_cell[2]-1;k++) 
            shape[ii++] = wuDv_minus_wvDu (x,y, adpolxy2[j]*adpolz[k]);

    	for(int i = 0; i <= order_cell[0]-2; i++) 
	  for(int j = 0; j <= order_cell[0]-2-i; j++) 
            shape[ii++] = wuDv_minus_wvDu (z,1-z, adpolxy1[i]*adpolxy2[j]);
      }
  }



  //------------------------------------------------------------------------
  // HCurlHighOrderHex
  //------------------------------------------------------------------------


  template<> template<typename Tx, typename TFA>  
  void HCurlHighOrderFE_Shape<ET_HEX> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    // Tx x = hx[0], y = hx[1], z = hx[2];
    Tx x = ip.x, y = ip.y, z = ip.z;

    Tx lami[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
                (1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
    Tx sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
                 (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 
    
    int ii = 12; 
    ArrayMem<Tx, 20> pol_xi(order+2),pol_eta(order+2),pol_zeta(order+2);
   
    // edges
    for (int i = 0; i < 12; i++)
      {
	int p = order_edge[i]; 
        IVec<2> e = GetEdgeSort (i, vnums);	  
	
	Tx xi  = sigma[e[1]]-sigma[e[0]];
	Tx lam_e = lami[e[0]]+lami[e[1]];  
        Tx bub = 0.25 * lam_e * (1 - xi*xi);

	// Nedelec0-shapes
	shape[i] = uDv (0.5 * lam_e, xi); 

	// High Order edges ... Gradient fields 
	if(p > 0 && usegrad_edge[i])
	  {
	    //LegendrePolynomial::
            EdgeOrthoPol::
	      EvalMult (p-1, 
			xi, bub, pol_xi);

	    for (int j = 0; j < p; j++)
              shape[ii++] = Du (pol_xi[j]);
	  }
      }
    
    //Faces 
    const FACE * faces = ElementTopology::GetFaces (ET_HEX);
    for (int i = 0; i<6; i++)
      {
	IVec<2> p = order_face[i];

	Tx lam_f(0);
	for (int j = 0; j < 4; j++)
	  lam_f += lami[faces[i][j]];

        {
          IVec<4> f = GetFaceSort (i, vnums);	  
          Tx xi  = sigma[f[0]] - sigma[f[1]]; 
          Tx eta = sigma[f[0]] - sigma[f[3]];
        
          Tx bub = lam_f*(1.0/16)*(1.0-eta*eta)*(1.0-xi*xi);
          QuadOrthoPol::Eval     (p[0]-1, xi, pol_xi);
          QuadOrthoPol::EvalMult (p[1]-1, eta, bub, pol_eta);
          
          if(usegrad_face[i])
            {
              // Gradientfields nabla(polxy*polz) 
              for (int k = 0; k <= p[0]-1; k++)
                for (int j = 0; j <= p[1]-1; j++)
                  shape[ii++] = Du (pol_xi[k] * pol_eta[j]);
            }
        }

	
	
	int qmax = 0;
	for (int j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] > vnums[faces[i][qmax]])
	    qmax = j;
	
	int q1 = (qmax+3)%4; 
	int q2 = (qmax+1)%4; 

	if(vnums[faces[i][q2]] > vnums[faces[i][q1]])
	  swap(q1,q2);  // fmax > f1 > f2

	int fmax = faces[i][qmax]; 
	int f1 = faces[i][q1]; 
	int f2 = faces[i][q2]; 
	      
	Tx xi = sigma[fmax]-sigma[f1]; 
	Tx eta = sigma[fmax]-sigma[f2]; 
    
	T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
	T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);

        /*
	//Gradient fields 
	if(usegrad_face[i])
	  for (int k = 0; k < p[0]; k++)
	    for (int j= 0; j < p[1]; j++)
              shape[ii++] = Du<3> (lam_f * pol_xi[k] * pol_eta[j]);
        */

	//Rotation of Gradient fields 
	for (int k = 0; k < p[0]; k++)
	  for (int j= 0; j < p[1]; j++)
            shape[ii++] = uDv_minus_vDu (pol_eta[j], lam_f * pol_xi[k]);

	// Missing ones 
	for(int j = 0; j < p[0];j++) 
          shape[ii++] = wuDv_minus_wvDu (Tx(0.5), eta, pol_xi[j]*lam_f); 

	for(int j = 0; j < p[1];j++) 
          shape[ii++] = wuDv_minus_wvDu (Tx(0.5), xi, pol_eta[j]*lam_f); 
      }


    
    {
    IVec<3> p = order_cell[0];
    if(usegrad_cell)
      if (p[0] >= 1 && p[1] >= 1 && p[2] >= 1)
        {
          QuadOrthoPol::EvalMult (p[0]-1, 2*x-1, x*(1-x), pol_xi);
          QuadOrthoPol::EvalMult (p[1]-1, 2*y-1, y*(1-y), pol_eta);
          QuadOrthoPol::EvalMult (p[2]-1, 2*z-1, z*(1-z), pol_zeta);
          
          for (int i = 0; i < p[0]; i++)
            for (int j = 0; j < p[1]; j++)
              {
                Tx pxy = pol_xi[i] * pol_eta[j];
                for (int k = 0; k < p[2]; k++)
                  shape[ii++] = Du (pxy * pol_zeta[k]);
              }
        }
    }

    
    // Element-based shapes
    T_ORTHOPOL::Calc(order_cell[0]+1,2*x-1,pol_xi);
    T_ORTHOPOL::Calc(order_cell[1]+1,2*y-1,pol_eta);
    T_ORTHOPOL::Calc(order_cell[2]+1,2*z-1,pol_zeta); 
    
    /*
    //Gradient fields
    if(usegrad_cell)
      for (int i=0; i<order_cell[0]; i++)
	for(int j=0; j<order_cell[1]; j++) 
	  for(int k=0; k<order_cell[2]; k++)
            shape[ii++] = Du<3> (pol_xi[i] * pol_eta[j] * pol_zeta[k]);
    */

    //Rotations of gradient fields
    for (int i=0; i<order_cell[0]; i++)
      for(int j=0; j<order_cell[1]; j++) 
	for(int k=0; k<order_cell[2]; k++)
	  {
            shape[ii++] = uDv_minus_vDu (pol_xi[i] * pol_eta[j], pol_zeta[k]);
            shape[ii++] = uDv_minus_vDu (pol_xi[i], pol_eta[j] * pol_zeta[k]);
	  } 
    
    for(int i = 0; i < order_cell[0]; i++) 
      for(int j = 0; j < order_cell[1]; j++)
        shape[ii++] = wuDv_minus_wvDu (z,1-z,pol_xi[i] * pol_eta[j]);

    for(int i = 0; i < order_cell[0]; i++) 
      for(int k = 0; k < order_cell[2]; k++)
        shape[ii++] = wuDv_minus_wvDu (y,1-y,pol_xi[i] * pol_zeta[k]);

    for(int j = 0; j < order_cell[1]; j++)
      for(int k = 0; k < order_cell[2]; k++)
        shape[ii++] = wuDv_minus_wvDu (x,1-x,pol_eta[j] * pol_zeta[k]);
  }



  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_HEX> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;        
    auto & ip = mip.IP();
    T x = ip(0), y = ip(1), z = ip(2);
    // T lam[4] = { x, y, z, 1-x-y-z };
    // Vec<3> pnts[4] = { { 1, 0, 0 }, { 0, 1, 0 } , { 0, 0, 1 }, { 0, 0, 0 } };
    // T lam[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
    // (1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
    T sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
                (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z};
    Vec<3> pnts[8] = 
      { 
	{ 0, 0, 0 },
	{ 1, 0, 0 },
	{ 1, 1, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 1 },
	{ 1, 0, 1 },
	{ 1, 1, 1 },
	{ 0, 1, 1 }
      };
    
    
    int facetnr = ip.FacetNr();
    int ii = 12;

    if (ip.VB() == BBND)
      { // edge shapes
        for (int i = 0; i < 12; i++)
          {
            int p = order_edge[i];
            if (i == facetnr)
              {
                IVec<2> e = GetEdgeSort (i, vnums);
                T xi = sigma[e[1]]-sigma[e[0]];
                Vec<3> tauref = pnts[e[1]] - pnts[e[0]];
                Vec<3,T> tau = mip.GetJacobian()*tauref;
                tau /= mip.GetMeasure();
                LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              Vec<3,T> vshape = val * tau;
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
        throw Exception ("H(curl)-hex: dual shapes supported only on edges");
        for (int i = 0; i < 12; i++)
          ii += order_edge[i];
      }

    /*
      // just copied from tet, needs adaption ...
    if (ip.VB() == BND)
      {
	//AutoDiff<3,T> xa(ip(0), 0), ya(ip(1),1), za(ip(2),2);
	//AutoDiff<3,T> lami[4] = { xa, ya, za, (T)(1.0) };

	for (int f = 0; f < 4; f++)
	  {
	    int p = order_face[f][0];
	     if (f == facetnr)
	       {
		 IVec<4> fav = GetFaceSort (facetnr, vnums);
		 //AutoDiff<3,T> adxi = lami[fav[0]]-lami[fav[2]];
		 //AutoDiff<3,T> adeta = lami[fav[1]]-lami[fav[2]];
		 Vec<3> adxi = pnts[fav[0]] - pnts[fav[2]];
		 Vec<3> adeta = pnts[fav[1]] - pnts[fav[2]];
		 T xi = lam[fav[0]];
		 T eta = lam[fav[1]];
		 
		 Matrix<> F(3,2);
		 F.Cols(0,1) = adxi;//Vec<3,T>(adxi.DValue(0),adxi.DValue(1),adxi.DValue(2));
		 F.Cols(1,2) = adeta;//Vec<3,T>(adeta.DValue(0),adeta.DValue(1),adeta.DValue(2));
		 
		 Matrix<> Ftmp(2,2);
		 Ftmp = Trans(F)*F;
		 auto det = sqrt(Ftmp(0,0)*Ftmp(1,1)-Ftmp(1,0)*Ftmp(0,1));
		 
		 DubinerBasis::Eval(order-2, xi, eta,
                                    SBLambda([&] (size_t nr, auto val)
                                             {
                                               shape[ii++] = 1/(det*mip.GetMeasure())*mip.GetJacobian()*(F*Vec<2,T> (val, 0));
                                               shape[ii++] = 1/(det*mip.GetMeasure())*mip.GetJacobian()*(F*Vec<2,T> (val*xi, val*eta));		
                                             }));
		 LegendrePolynomial::Eval(order-2,xi,
					  SBLambda([&] (size_t nr, auto val)
						   {
						     shape[ii++] = 1/(det*mip.GetMeasure())*mip.GetJacobian()*(F*Vec<2,T>(0, val));
						   }));
	       }
	     else
	      ii += (p+1)*(p-1);
	  }
      }
    else
      {
        for (int i = 0; i < 4; i++)
	  {
	    int p = order_face[i][0];
	    ii += (p+1)*(p-1);
	  }
      }
    if (ip.VB() == VOL)
      {
	// auto xphys = mip.GetPoint()(0);
        // auto yphys = mip.GetPoint()(1);
	// auto zphys = mip.GetPoint()(2);
	
	LegendrePolynomial leg;
	JacobiPolynomialAlpha jac1(1);    
	leg.EvalScaled1Assign 
	  (order-3, lam[2]-lam[3], lam[2]+lam[3],
	   SBLambda ([&](size_t k, T polz) LAMBDA_INLINE
		     {
		       // JacobiPolynomialAlpha jac(2*k+1);
		       JacobiPolynomialAlpha jac2(2*k+2);
		       
		       jac1.EvalScaledMult1Assign
			 (order-3-k, lam[1]-lam[2]-lam[3], 1-lam[0], polz, 
			  SBLambda ([&] (size_t j, T polsy) LAMBDA_INLINE
				    {
				      // JacobiPolynomialAlpha jac(2*(j+k)+2);
				      jac2.EvalMult(order-3 - k - j, 2 * lam[0] - 1, polsy, 
						    SBLambda([&](size_t j, T val) LAMBDA_INLINE
							     {
							       shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*Vec<3,T>(val*x, val*y, val*z);
							       shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*Vec<3,T>(val, 0, 0);
							       shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*Vec<3,T>(0, val, 0);
							     }));
				      jac2.IncAlpha2();
				    }));
		       jac1.IncAlpha2();
		     }));


        DubinerBasis::Eval(order-3, x, y,
                           SBLambda([&] (size_t nr, auto val)
                                    {
                                      shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*Vec<3,T> (0, 0, val);
                                    }));
      }
    */
  }








  

  //------------------------------------------------------------------------
  //            Pyramid
  //------------------------------------------------------------------------


  template<> template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE_Shape<ET_PYRAMID> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;  

    // Tx x = hx[0], y = hx[1], z = hx[2];
    Tx x = ip.x, y = ip.y, z = ip.z;
    
    //if(z.Value()==1.) z.Value() -=1.e-8; 
    z.Value() = z.Value()*(1-1e-12);

    Tx xt = x/(1-z); 
    Tx yt = y/(1-z); 
    Tx sigma[5] = {(1-xt)+(1-yt)+(1-z),xt+(1-yt)+(1-z), xt + yt + (1-z), 
			    (1-xt)+yt+(1-z),z}; 

    Tx lami[5] = {(1-xt)*(1-yt)*(1-z),xt*(1-yt)*(1-z), xt * yt * (1-z), 
			   (1-xt)*yt*(1-z),z}; 

    Tx lambda[5] = {(1-xt)*(1-yt),xt*(1-yt), xt * yt, 
                             (1-xt)*yt,z}; 
        
       
    ArrayMem<Tx, 20> pol_xi(order+2), pol_eta(order+2), pol_zeta(order+2);
    
    int ii =8; 
 
    // horizontal edges incl. Nedelec 0
    for (int i = 0; i < 4; i++)
      {
        int p = order_edge[i];
        IVec<2> e = GetEdgeSort (i, vnums);	  
	
	Tx xi  = sigma[e[1]] - sigma[e[0]];   
	Tx lam_t = lambda[e[1]] + lambda[e[0]]; 
        
        shape[i] = uDv (0.5 * (1-z)*(1-z)*lam_t, xi);

	if(p > 0 && usegrad_edge[i])
	  {
            Tx bub = 0.25*(1-xi*xi)*(1-z)*(1-z)*lam_t;
	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalScaledMult (p-1,
			      xi*(1-z), 1-z, bub,
                              SBLambda ([&](int i, Tx val)
                                        {
                                          shape[ii++] = Du(val);
                                        }));
	  }
      }
    
    // vertical edges incl. Nedelec 0  
    for(int i = 4; i < 8; i++)
      {
        int p = order_edge[i];
        IVec<2> e = GetEdgeSort (i, vnums);	  

        shape[i] = uDv_minus_vDu (lami[e[0]], lami[e[1]]);

	if (p > 0 && usegrad_edge[i])
	  {
            Tx xi = lami[e[1]]-lami[e[0]]; 
            Tx lam_e = lami[e[0]]+lami[e[1]];
            Tx bub = 0.25 * (lam_e*lam_e-xi*xi);

	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalScaledMult (p-1,
                              xi, lam_e, bub, 
                              SBLambda ([&](int i, Tx val)
                                        {
                                          shape[ii++] = Du(val);
                                        }));
	  }
      }

    const FACE * faces = ElementTopology::GetFaces (ET_PYRAMID); 

    // trig face dofs
    for (int i = 0; i < 4; i++)
      if (order_face[i][0] >= 2)
	{
	  int p = order_face[i][0];
	  Tx lam_face = lambda[faces[i][0]] + lambda[faces[i][1]];  
	  Tx bary[3] = 
	    {(sigma[faces[i][0]]-(1-z)-lam_face)*(1-z), 
	     (sigma[faces[i][1]]-(1-z)-lam_face)*(1-z), z}; 
			     
	  int fav[3] = {0, 1, 2};
	  if(vnums[faces[i][fav[0]]] > vnums[faces[i][fav[1]]]) swap(fav[0],fav[1]); 
	  if(vnums[faces[i][fav[1]]] > vnums[faces[i][fav[2]]]) swap(fav[1],fav[2]);
	  if(vnums[faces[i][fav[0]]] > vnums[faces[i][fav[1]]]) swap(fav[0],fav[1]); 	

	  if(usegrad_face[i])
            {
              Tx bub = lam_face * bary[fav[0]]*bary[fav[1]]*bary[fav[2]];
              TrigOrthoPolGrad::
                EvalMult (p-2, bary[fav[0]], bary[fav[1]], bub, 
                          SBLambda ([&](int nr, Tx val)
                                    {
                                      shape[ii++] = Du (val);
                                    }));
            }
          
          /*
	  // phi = pol_xi * pol_eta * lam_face; 
	  // Type 1: Gradient Functions 
	  if(usegrad_face[i])
	    for(int j=0;j<= p-2; j++)
	      for(int k=0;k<=p-2-j; k++)
                shape[ii++] = Du<3> (pol_xi[j] * pol_eta[k]);
          */

          // if (type1) {          
          if (true) {
            
            DubinerBasis::EvalMult
              (p-2, bary[fav[0]], bary[fav[1]],
               bary[fav[0]]*lam_face,
               SBLambda
               ([&](int nr, Tx val)
               {
                 shape[ii++] =  wuDv_minus_wvDu(bary[fav[1]], bary[fav[2]], val);
               }));
	    
            LegendrePolynomial::EvalMult
              (p-2, bary[fav[2]]-bary[fav[1]], bary[fav[2]]*lam_face, 
               SBLambda
               ([&] (int j, Tx val)
               {
                 shape[ii++] = wuDv_minus_wvDu(bary[fav[1]], bary[fav[0]], val);
             }));
          }

          else {
            T_TRIGFACESHAPES::CalcSplitted(p+1, bary[fav[2]]-bary[fav[1]], 
                                           bary[fav[0]],pol_xi,pol_eta);
            
            for(int j=0;j<=p-2;j++) pol_eta[j] *= lam_face;
            
            // Type 2:  
            for(int j=0;j<= p-2; j++)
              for(int k=0;k<=p-2-j; k++)
                shape[ii++] = uDv_minus_vDu (pol_eta[k], pol_xi[j]);
            
            // Type 3: Nedelec-based ones (Ned_0*v_j)
            for(int j=0;j<=p-2;j++)
              shape[ii++] = wuDv_minus_wvDu (bary[fav[1]], bary[fav[2]], pol_eta[j]);
          }
	}


    // quad face 
    if (order_face[4][0] >= 1)
      {
	int px = order_face[4][0];
	int py = order_face[4][0]; // SZ-Attentione 
	int p = max2(px, py);

	Tx fac(1.0);
	for (int k = 1; k <= p+1; k++) fac *= (1-z);

	IVec<4> f = GetFaceSort (4, vnums);	  
	Tx xi  = sigma[f[0]] - sigma[f[1]]; 
	Tx eta = sigma[f[0]] - sigma[f[3]];

        if (usegrad_face[4])
          {
            // Type 1: Gradient-fields 
            
            QuadOrthoPol::
              EvalMult (px-1, xi, fac*0.25*(1-xi*xi), pol_xi);
            QuadOrthoPol::
              EvalMult (py-1, eta, 0.25*(1-eta*eta), pol_eta);

            for (int k = 0; k <= px-1; k++) 
              for (int j = 0; j <= py-1; j++, ii++) 
                shape[ii] = Du (pol_xi[k] * pol_eta[j]);
          }
        
	int fmax = 0;
	for (int l=1; l<4; l++) 
	  if (vnums[l] > vnums[fmax]) fmax = l;  

	int f1 = (fmax+3)%4;
	int f2 = (fmax+1)%4; 
	if(vnums[f1]>vnums[f2]) swap(f1,f2);  // fmax > f2 > f1 

 	xi  = sigma[fmax] - sigma[f2]; 
	eta = sigma[fmax] - sigma[f1];

	T_ORTHOPOL::Calc (px+1, xi, pol_xi);	
	T_ORTHOPOL::Calc (py+1, eta, pol_eta);
	
	for(int k = 0; k < py; k++) pol_eta[k] *= fac;

	// Type 2: 
	for (int k = 0; k < px; k++) 
	  for (int j = 0; j < py; j++) 
            shape[ii++] = uDv_minus_vDu (pol_eta[j], pol_xi[k]);

	// Type 3:
	for (int k = 0; k < px; k++)
          shape[ii++] = uDv (0.5*pol_xi[k]*fac, eta);

	for (int k = 0; k < py; k++)
          shape[ii++] = uDv (0.5*pol_eta[k] /* *fac */, xi); 
      }

    if (order_cell[0] >= 2)
      {
	int pp = order_cell[0];
	// According H^1 terms: 
	// u_i = L_i+2(2xt-1)
	// v_j = L_j+2(2yt-1) 
	// w_k = z * (1-z)^(k+2)  with 0 <= i,j <= k, 0<= k <= p-2  
	
	LegendrePolynomial::EvalMult (pp-1, 2*xt-1, xt*(1-xt), pol_xi);
        LegendrePolynomial::EvalMult (pp-1, 2*yt-1, yt*(1-yt), pol_eta);
	// T_ORTHOPOL::Calc (pp+3, 2*xt-1, pol_xi);
	// T_ORTHOPOL::Calc (pp+3, 2*yt-1, pol_eta);		
		
	pol_zeta[0] = z*(1-z)*(1-z);
	for (int k=1;k<=pp-2;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];

	if(usegrad_cell)
	  {
	    for(int k=0;k<= pp-2;k++)
              {
                for(int i=0;i<=k;i++)
                  for(int j=0;j<=k;j++)
                    shape[ii++] = Du (pol_xi[i]*pol_eta[j]*pol_zeta[k]);
              }
	  }
	
	// Type 2a: l.i. combinations of grad-terms   
	// shape = u_i \nabla(v_j) w_k 
	// shape = u_i v_j \nabla(w_k) 
	for(int k=0;k<= pp-2;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
              shape[ii] = uDv (pol_xi[i]*pol_zeta[k], pol_eta[j]);

	// Type 2b: shape = v_j w_k \nabla (xt) 
	//          shape = u_i w_k \nabla (yt)
	for(int k = 0;k<= pp-2;k++)
	  for(int j=0;j<=k;j++) 
            shape[ii++] = uDv (pol_eta[j]*pol_zeta[k], xt);
	
	for(int  k = 0;k<= pp-2;k++)
	  for (int i=0;i<=k;i++)
            shape[ii++] = uDv (pol_xi[i]*pol_zeta[k], yt);
	
	// 3rd component spans xi^i eta^j zeta^(k-1), i,j <= k
	// pol_zeta starts linear in zeta 
	// pol_xi and pol_eta quadratic in xi resp. eta 
	pol_zeta[0] = (1-z);
	for (int k=1;k<=pp;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];
     
	for(int k=0;k<= pp-1;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
              shape[ii] = uDv (pol_eta[j] * pol_xi[i] * pol_zeta[k], z);
      }
  }




  // dual shapes

  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_TRIG> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    // shape = 0;
    auto & ip = mip.IP();
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
    T x = ip(0), y = ip(1);
    T lam[3] = { x, y, 1-x-y };
    Vec<2,T> pnts[3] = { { 1, 0 }, { 0, 1 } , { 0, 0 } };
    int facetnr = ip.FacetNr();

    if (ip.VB() == BND)
      { // facet shapes
        int ii = 3;
        for (int i = 0; i < 3; i++)
          {
            int p = order_edge[i];
            if (i == facetnr)
              {
                IVec<2> e = GetEdgeSort (i, vnums);
                T xi = lam[e[1]]-lam[e[0]];
                Vec<2,T> tauref = pnts[e[1]] - pnts[e[0]];
                auto tau = mip.GetJacobian()*tauref;
                tau /= mip.GetMeasure();

                LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              auto vshape = val * tau;
                              if (nr==0)
                                shape[i] = vshape;
                              else
                                shape[ii+nr-1] = vshape;
                            }));
              }
            ii += p;
          }
      }
    if (ip.VB() == VOL)
      { // inner shapes
        int ii = 3;
        for (int i = 0; i < 3; i++)
          ii += order_edge[i];


        auto trafo = 1/mip.GetMeasure()*mip.GetJacobian();
        DubinerBasis::Eval(order-2, x, y,
                           SBLambda([&] (size_t nr, auto val)
                                    {
                                      shape[ii++] = trafo * Vec<2,T> (val, 0);
                                      if (type1)
                                        shape[ii++] = trafo * Vec<2,T> (0, val);
                                      else
                                        shape[ii++] = trafo * Vec<2,T> (val*x, val*y);
                                    }));
        if (!type1)
          LegendrePolynomial::Eval(order-2,x,
                                   SBLambda([&] (size_t nr, auto val)
                                            {
                                              shape[ii++] = trafo * Vec<2,T>(0,val);
                                            }));
      }
  }

  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_QUAD> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    auto & ip = mip.IP();
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
    T x = ip(0), y = ip(1);
    Vec<2,T> pnts[4] = { { 0, 0 },
	{ 1, 0 },
	{ 1, 1 },
	{ 0, 1 } };
    T sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    int facetnr = ip.FacetNr();

    if (ip.VB() == BND)
      { // facet shapes
        int ii = 4;
        for (int i = 0; i < 4; i++)
          {
            int p = order_edge[i];
            if (i == facetnr)
              {
                IVec<2> e = GetEdgeSort (i, vnums);
                T xi = sigma[e[1]]-sigma[e[0]];
                Vec<2,T> tauref = pnts[e[1]] - pnts[e[0]];
                auto tau = mip.GetJacobian()*tauref;
                tau /= mip.GetMeasure();

                LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              auto vshape = val * tau;
                              if (nr==0)
                                shape[i] = vshape;
                              else
                                shape[ii+nr-1] = vshape;
                            }));
              }
            ii += p;
          }
      }
    if (ip.VB() == VOL)
      { // inner shapes
        int ii = 4;
        for (int i = 0; i < 4; i++)
          ii += order_edge[i];

        //do not sort face!
        //IVec<4> f = GetFaceSort (0, vnums);  
        //T xi = sigma[f[0]]-sigma[f[1]]; 
        //T eta = sigma[f[0]]-sigma[f[3]]; 
        T xi = sigma[0]-sigma[1]; 
        T eta = sigma[0]-sigma[3]; 
        ArrayMem<T, 20> polx(order+2), poly(order+2);
        LegendrePolynomial (order, xi, polx);
        LegendrePolynomial (order, eta, poly);

        int p = order_face[0][0];
        
        for (int i = 0; i < p+1; i++)
          for (int j = 0; j < p; j++)
            {
              shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*Vec<2,T>(polx[i] * poly[j],0);
              shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*Vec<2,T>(0,polx[j] * poly[i]);
            }	
      }
  }

  

  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_TET> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    // shape = 0;
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;        
    auto & ip = mip.IP();
    T x = ip(0), y = ip(1), z = ip(2);
    T lam[4] = { x, y, z, 1-x-y-z };
    Vec<3> pnts[4] = { { 1, 0, 0 }, { 0, 1, 0 } , { 0, 0, 1 }, { 0, 0, 0 } };
    int facetnr = ip.FacetNr();
    int ii = 6;

    if (ip.VB() == BBND)
      { // edge shapes
        for (int i = 0; i < 6; i++)
          {
            int p = order_edge[i] * usegrad_edge[i];
            if (i == facetnr)
              {
                IVec<2> e = GetEdgeSort (i, vnums);
                T xi = lam[e[1]]-lam[e[0]];
                Vec<3> tauref = pnts[e[1]] - pnts[e[0]];
                Vec<3,T> tau = mip.GetJacobian()*tauref;
                tau /= mip.GetMeasure();
                LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              Vec<3,T> vshape = val * tau;
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
        for (int i = 0; i < 6; i++)
          ii += order_edge[i] * usegrad_edge[i];
      }
    if (ip.VB() == BND)
      {
	//AutoDiff<3,T> xa(ip(0), 0), ya(ip(1),1), za(ip(2),2);
	//AutoDiff<3,T> lami[4] = { xa, ya, za, (T)(1.0) };

	for (int f = 0; f < 4; f++)
	  {
	    int p = order_face[f][0];
	     if (f == facetnr)
	       {
		 IVec<4> fav = GetFaceSort (facetnr, vnums);
		 //AutoDiff<3,T> adxi = lami[fav[0]]-lami[fav[2]];
		 //AutoDiff<3,T> adeta = lami[fav[1]]-lami[fav[2]];
		 Vec<3> adxi = pnts[fav[0]] - pnts[fav[2]];
		 Vec<3> adeta = pnts[fav[1]] - pnts[fav[2]];
		 T xi = lam[fav[0]];
		 T eta = lam[fav[1]];
		 
		 Matrix<> F(3,2);
		 F.Col(0) = adxi;//Vec<3,T>(adxi.DValue(0),adxi.DValue(1),adxi.DValue(2));
		 F.Col(1) = adeta;//Vec<3,T>(adeta.DValue(0),adeta.DValue(1),adeta.DValue(2));
		 
		 Matrix<> Ftmp(2,2);
		 Ftmp = Trans(F)*F;
		 auto det = sqrt(Ftmp(0,0)*Ftmp(1,1)-Ftmp(1,0)*Ftmp(0,1));
		 
		 DubinerBasis::Eval(order-2, xi, eta,
                                    SBLambda([&] (size_t nr, auto val)
                                             {						 
                                               shape[ii++] = 1/(det*mip.GetMeasure())*mip.GetJacobian()*(F*Vec<2,T> (val, 0));
					       if (type1)
						 shape[ii++] = 1/(det*mip.GetMeasure())*mip.GetJacobian()*(F*Vec<2,T>(0, val));
					       else
						 shape[ii++] = 1/(det*mip.GetMeasure())*mip.GetJacobian()*(F*Vec<2,T> (val*xi, val*eta));		
                                             }));
		 if (!type1)
		   LegendrePolynomial::Eval(order-2,xi,
					  SBLambda([&] (size_t nr, auto val)
						   {
						     shape[ii++] = 1/(det*mip.GetMeasure())*mip.GetJacobian()*(F*Vec<2,T>(0, val));
						   }));
	       }
	     else
	       ii += (p+!type1)*(p-1);
	     
	  }
      }
    else
      {
        for (int i = 0; i < 4; i++)
	  {
	    int p = order_face[i][0];
	    ii += (p+!type1)*(p-1);
	  }
      }
    if (ip.VB() == VOL)
      {
	// auto xphys = mip.GetPoint()(0);
        // auto yphys = mip.GetPoint()(1);
	// auto zphys = mip.GetPoint()(2);

	auto trafo = 1/mip.GetMeasure()*mip.GetJacobian();

	DubinerBasis3D::Eval (order-3, lam[0], lam[1], lam[2],
			      SBLambda([&](size_t j, T val) LAMBDA_INLINE
	 						     {
							       if (type1)
								 shape[ii++] = trafo*Vec<3,T> (0, 0, val);
							       else
								 shape[ii++] = trafo*Vec<3,T>(val*x, val*y, val*z);
	 						       shape[ii++] = trafo*Vec<3,T>(val, 0, 0);
	 						       shape[ii++] = trafo*Vec<3,T>(0, val, 0);
	 						     }));
	
	if (!type1)
	  DubinerBasis::Eval(order-3, x, y,
                           SBLambda([&] (size_t nr, auto val)
                                    {
                                      shape[ii++] = trafo*Vec<3,T> (0, 0, val);
                                    }));
      }
  }

  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_PRISM> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    // shape = 0;
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;
    auto & ip = mip.IP();
    T x = ip(0), y = ip(1), z = ip(2);
    T lam[6] = { x, y, 1-x-y, x, y, 1-x-y };
    T muz[6]  = { 1-z, 1-z, 1-z, z, z, z };

    T sigma[6];
    for (int i = 0; i < 6; i++) sigma[i] = lam[i] + muz[i];

    Vec<3> pnts[6] = { { 1, 0, 0 }, { 0, 1, 0 } , { 0, 0, 0 },
                       { 1, 0, 1 }, { 0, 1, 1 }, { 0, 0, 1 }};
    int facetnr = ip.FacetNr();
    int ii = 9;

    if (ip.VB() == BBND)
      {
        for (int i = 0; i < 9; i++)
          {
            if(order_edge[i] > 1) throw Exception("Dual shapes for prisms for order > 1 not implemented!");
            int p = order_edge[i] * usegrad_edge[i];
            if (i == facetnr)
              {
                IVec<2> e = GetEdgeSort (i, vnums);
                T xi = sigma[e[1]] - sigma[e[0]];
                Vec<3> tauref = pnts[e[1]] - pnts[e[0]];
                Vec<3,T> tau = mip.GetJacobian()*tauref;
                tau /= mip.GetMeasure();
                LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              Vec<3,T> vshape = val * tau;
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
        for (int i = 0; i < 9; i++)
          ii += order_edge[i] * usegrad_edge[i];
      }
    if (ip.VB() == BND)
      {
        // TODO
      }
    else
      {
        for (int i = 0; i < 4; i++)
	  {
	    int p = order_face[i][0];
	    ii += (p+1)*(p-1);
	  }
      }
    if (ip.VB() == VOL)
      {
        // TODO
      }
  }
  
  template <ELEMENT_TYPE ET, 
            template <ELEMENT_TYPE ET2> class TSHAPES, 
            typename BASE>
  void HCurlHighOrderFE<ET,TSHAPES,BASE> ::
  CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const
  {
    shape.AddSize(ndof, bmip.DimSpace()) = 0.0;
    Switch<4-DIM>
      (bmip.DimSpace()-DIM,[this,&bmip,shape](auto CODIM)
       {
         auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);
         static_cast<const HCurlHighOrderFE_Shape<ET>*> (this)
           -> CalcDualShape2 (mip, SBLambda([shape] (size_t i, auto val)
                                            {
                                              shape.Row(i) = val;
                                            }));
       });
  }
  

  template <ELEMENT_TYPE ET, 
            template <ELEMENT_TYPE ET2> class TSHAPES, 
            typename BASE>
  void HCurlHighOrderFE<ET,TSHAPES,BASE> ::
  CalcDualShape (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> shapes) const
  {
 
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,shapes](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         shapes.AddSize(ndof*DIMSPACE, mir.Size()) = 0.0;
         for (size_t i = 0; i < mir.Size(); i++)
           static_cast<const HCurlHighOrderFE_Shape<ET>*> (this)
             -> CalcDualShape2 (mir[i], SBLambda([shapes,i,DIMSPACE] (size_t j, auto val)
                                                 {
                                                   for (size_t k = 0; k < DIMSPACE; k++)
                                                     shapes(j*DIMSPACE+k, i) = val(k);
                                                 }));
       });
  }

  template <ELEMENT_TYPE ET, 
            template <ELEMENT_TYPE ET2> class TSHAPES, 
            typename BASE>
  void HCurlHighOrderFE<ET,TSHAPES,BASE> ::
  EvaluateDual (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIMSPACE,SIMD<double>> sum (SIMD<double>(0.0));
             static_cast<const HCurlHighOrderFE_Shape<ET>*> (this)
               -> CalcDualShape2 (mir[i], SBLambda([&sum, coefs] (size_t j, auto val)
                                                   {
                                                     sum += coefs(j) * val;
                                                   }));
             for (size_t k = 0; k < DIMSPACE; k++)
               values(k, i) = sum(k);
           }
       });
  }

  template <ELEMENT_TYPE ET, 
            template <ELEMENT_TYPE ET2> class TSHAPES, 
            typename BASE>
  void HCurlHighOrderFE<ET,TSHAPES,BASE> ::
  AddDualTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
                BareSliceVector<double> coefs) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIMSPACE,SIMD<double>> value = values.Col(i);
             // for (size_t k = 0; k < DIMSPACE; k++)
             // value(k) = values(k, i);
             
             static_cast<const HCurlHighOrderFE_Shape<ET>*> (this)
               -> CalcDualShape2 (mir[i], SBLambda([value, coefs] (size_t j, auto val)
                                                   {
                                                     coefs(j) += HSum(InnerProduct(value,val));
                                                   }));
           }
       });
  }
}

#endif
