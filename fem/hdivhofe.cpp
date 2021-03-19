#define FILE_HDIVHOFE_CPP

#include <fem.hpp>
#include <hdivhofe.hpp>
/*
#include <thdivfe_impl.hpp>
#include <hdivhofe_impl.hpp>
#include <hdivhofefo.hpp>
*/

namespace ngfem
{  

   
  //------------------------------------------------------------------------
  // HDivHighOrderFiniteElement
  //------------------------------------------------------------------------

  /*
  template <int D>
  HDivHighOrderFiniteElement<D> ::
  HDivHighOrderFiniteElement (ELEMENT_TYPE aeltype)
    : HDivFiniteElement<D> (-1, -1)
  {
    for (int i = 0; i < 8; i++)
      vnums[i] = i;

    ho_div_free = 0;
    only_ho_div = 0;
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
  SetOrderInner (INT<D> oi)
  {
    order_inner = oi;
    ComputeNDof();
  }

  template <int D>
  void HDivHighOrderFiniteElement<D> :: Print (ostream & ost) const
  {
    ELEMENT_TYPE et = this->ElementType();
    ost << "HDivHighOrderFiniteElement<" << ElementTopology::GetElementName(et) << ">:" << endl;
    if (D == 2)
      for (int j = 0; j < ElementTopology::GetNEdges(et); j++)
	ost << "order_edge[" << j << "] = " << order_edge[j] << endl;
    else
      for (int j = 0; j < ElementTopology::GetNFaces(et); j++)
	ost << "order_face[" << j << "] = " << order_face[j] << endl;
  }
  */


  //------------------------------------------------------------------------
  // HDivHighOrderNormalFiniteElement
  //------------------------------------------------------------------------
  /*
  template <int D>
  HDivHighOrderNormalFiniteElement<D> ::
  HDivHighOrderNormalFiniteElement ()
    : HDivNormalFiniteElement<D> (-1, -1)
  {
    // for (int i = 0; i < 4; i++)
    // vnums[i] = i;
    ;
  }
  */
  /*
  template <int D>
  void HDivHighOrderNormalFiniteElement<D>::
  SetVertexNumbers (FlatArray<int> & avnums)
  {

    for (int i = 0; i < avnums.Size(); i++)
      vnums[i] = avnums[i];
  }
  */

  /*
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
  */
  
  //------------------------------------------------------------------------
  // HDivHighOrderNormalSegm
  //------------------------------------------------------------------------
  template <class T_ORTHOPOL>
  HDivHighOrderNormalSegm<T_ORTHOPOL> :: HDivHighOrderNormalSegm (int aorder)
  // : HDivHighOrderNormalFiniteElement<1>()
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

  /*
  template <class T_ORTHOPOL>
  void HDivHighOrderNormalSegm<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip,
							 FlatVector<> shape) const
  {
    AutoDiff<1> x (ip(0), 0);
    AutoDiff<1> lam[2] = { x, 1-x };

    ArrayMem<AutoDiff<1>,10> adpol1(order);
	
    INT<2> e = ET_trait<ET_SEGM>::GetEdgeSort (0, vnums);	  
    
    shape[0] = -lam[e[0]].DValue(0);

    int p = order_inner[0]; 
    // LegendrePolynomial::
    IntLegNoBubble::
      EvalMult (p-1, 
                lam[e[1]]-lam[e[0]], lam[e[0]]*lam[e[1]], adpol1);
    
    for(int j = 0; j < p; j++) 	      
      shape[j+1] = -adpol1[j].DValue(0);
      }
  */



  template class HDivHighOrderNormalSegm<IntegratedLegendreMonomialExt>;
  template class HDivHighOrderNormalSegm<TrigExtensionMonomial>;
  //template class HDivHighOrderNormalSegm<TrigExtensionOptimal>;
  //template class HDivHighOrderNormalSegm<TrigExtensionMin>;

  //------------------------------------------------------------------------
  // HDivHighOrderNormalQuad
  //------------------------------------------------------------------------
  template <class T_ORTHOPOL>
  HDivHighOrderNormalQuad<T_ORTHOPOL> :: HDivHighOrderNormalQuad (int aorder)
  // : HDivHighOrderNormalFiniteElement<2>()
  {
    order_inner = INT<2>(aorder,aorder);
    ComputeNDof();
  }

  template <class T_ORTHOPOL>
  void HDivHighOrderNormalQuad<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = (order_inner[0] < 0) ? 0 : (1 + order_inner[0]*order_inner[1] + order_inner[0] + order_inner[1]);
    order = max2(order_inner[0],order_inner[1]);
    order++; // order used for numerical integration
  }

#ifdef OLD
  template <class T_ORTHOPOL>
  void HDivHighOrderNormalQuad<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip,
							 FlatVector<> shape) const
  {
    AutoDiff<2> x (ip(0), 0);
    AutoDiff<2> y (ip(1), 1);

    // AutoDiff<2> lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};
    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};

    shape = 0.0;

    int ii = 1;

    INT<2> p = order_inner;
    // int pp = max2(p[0],p[1]); 
    
    ArrayMem<AutoDiff<2>,20> pol_xi(p[0]+1), pol_eta(p[1]+1);

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

    /*
    T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
    T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);
    */
    IntLegNoBubble::EvalMult (p[0]-1, xi, 1-xi*xi, pol_xi);
    IntLegNoBubble::EvalMult (p[1]-1, eta, 1-eta*eta, pol_eta);

    // Typ 1
    for (int k = 0; k < p[0]; k++)
      for (int l = 0; l < p[1]; l++, ii++)
	shape(ii) = 2.*(pol_eta[l].DValue(0)*pol_xi[k].DValue(1)-pol_eta[l].DValue(1)*pol_xi[k].DValue(0));

    //Typ 2
    for (int k = 0; k < p[0]; k++)
      shape(ii++) = -eta.DValue(0)*pol_xi[k].DValue(1) + eta.DValue(1)*pol_xi[k].DValue(0); 
    for (int k = 0; k < p[1]; k++)
      shape(ii++)   = -xi.DValue(0)*pol_eta[k].DValue(1) + xi.DValue(1)*pol_eta[k].DValue(0);
  }
#endif 

  template class HDivHighOrderNormalQuad<IntegratedLegendreMonomialExt>;
  template class HDivHighOrderNormalQuad<TrigExtensionMonomial>;
  // template class HDivHighOrderNormalQuad<TrigExtensionOptimal>;
  template class HDivHighOrderNormalQuad<TrigExtensionMin>;

  //------------------------------------------------------------------------
  // HDivHighOrderNormalTrig
  //------------------------------------------------------------------------

  template <class T_ORTHOPOL>
  HDivHighOrderNormalTrig<T_ORTHOPOL> :: HDivHighOrderNormalTrig (int aorder)
  // : HDivHighOrderNormalFiniteElement<2>()
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

#ifdef OLD
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

    int ii, is, ie, iop;

    int p = order_inner[0];
    ii = 1;

    int fav[3];
    for(int i=0;i<3;i++) fav[i] = i;

    //Sort vertices first edge op minimal vertex
    int fswap = 1;
    if(vnums[fav[0]] > vnums[fav[1]]) { swap(fav[0],fav[1]); fswap *= -1; }
    if(vnums[fav[1]] > vnums[fav[2]]) { swap(fav[1],fav[2]); fswap *= -1; }
    if(vnums[fav[0]] > vnums[fav[1]]) { swap(fav[0],fav[1]); fswap *= -1; }

    is = fav[0]; ie = fav[1]; iop = fav[2];

    AutoDiff<2> ls = lami[is];
    AutoDiff<2> le = lami[ie];
    AutoDiff<2> lo = lami[iop];

    //AutoDiff<3> lsle = lami[is]*lami[ie];
    for (int j = 0; j < 2; j++)
      {
        ls.DValue(j) = dlami(is,j);
        le.DValue(j) = dlami(ie,j);
        lo.DValue(j) = dlami(iop,j);
      }

    Vec<2> nedelec;
    for (int j = 0; j < 2; j++)
      nedelec(j) = ls.Value()*le.DValue(j) - le.Value()*ls.DValue(j);

    // RT_0-normal low order shapes
    shape(0) = fswap;


    /*
    AutoDiff<2> lsle=ls*le;
    Vec<2> grad1, grad2;
    ArrayMem<AutoDiff<2>, 100> ad_rec_pol1(p);
    ArrayMem<AutoDiff<2>, 100> ad_rec_pol2(p);

    ScaledLegendrePolynomial(p-1, le-ls, 1-lo, ad_rec_pol1);
    // IntLegNoBubble::EvalScaled (p-1, le-ls, 1-lo, ad_rec_pol1); 
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
    */

    ArrayMem<AutoDiff<2>, 20> adpol1(p);
    ArrayMem<AutoDiff<2>, 20> adpol2(p);

    IntLegNoBubble::EvalScaledMult (p-1, le-ls, le+ls, ls*le, adpol1); 
    /*
    // Typ 1
    for (int k = 0; k <= p-1; k++)
      {
        IntegratedJacobiPolynomialAlpha jac(2*k+3);
        jac.EvalMult(p-1-k, 2*lo-1, lo, adpol2);
	for (int l = 0; l <= p-1-k; l++, ii++)
          shape(ii) = Cross (adpol2[l], adpol1[k]).DValue(0);
      }
    */

    // Vector<> hshape(shape.Size());
    // ArrayMem<double, 20> adpol2b(p);
    // hshape = 0.0;
    // int iii = 1;
    for (int k = 0; k <= p-1; k++)
      {
        JacobiPolynomialAlpha jac(2*k+3);
	double factor = Cross (lo, adpol1[k]).DValue(0);
        jac.EvalMult(p-1-k, 2*lo.Value()-1, factor, shape+ii);
	ii += p-k;
      }
    // cout << "shape = " << endl << shape << endl << " hshape = " << hshape << endl;

    IntegratedJacobiPolynomialAlpha jac(3);
    jac.EvalMult(p-1, 2*lo-1, lo, adpol2);

    // Typ 2
    double curlned;
    curlned = 2.* (ls.DValue(0)*le.DValue(1) - ls.DValue(1)*le.DValue(0));
    for (int k = 0; k <= p-1; k++, ii++)
      shape(ii) = adpol2[k].DValue(0)*nedelec(1) - adpol2[k].DValue(1)*nedelec(0) 
        + adpol2[k].Value()*curlned;
  }
#endif

  template class HDivHighOrderNormalTrig<IntegratedLegendreMonomialExt>;
  template class HDivHighOrderNormalTrig<TrigExtensionMonomial>;
  // template class HDivHighOrderNormalTrig<TrigExtensionOptimal>;
  template class HDivHighOrderNormalTrig<TrigExtensionMin>;










  template <ELEMENT_TYPE ET>
  void HDivHighOrderFE<ET> :: 
  ComputeNDof()
  {
    if (DIM == 2)
      {
        if (only_ho_div)
	  {
	    if (ET == ET_TRIG)
              {
                ndof = order_inner[0]*(order_inner[0]+1)/2 - 1;
                order = order_inner[0];
              }
	    else
              {
                ndof = order_inner[0]*order_inner[1] + order_inner[0] + order_inner[1];
                order = max2(order_inner[0], order_inner[1])+1;
              }
	    return;
	  }
        else
	  {
	    ndof = ET_trait<ET>::N_EDGE;
	    
	    for(int i = 0; i < ET_trait<ET>::N_EDGE; i++)
	      ndof += order_facet[i][0];
	    
	    if (ET == ET_TRIG)
	      {
		if (order_inner[0] > 1)
		  { 
		    if (ho_div_free)
		      ndof += order_inner[0]*(order_inner[0]-1)/2;
		    else
		      ndof += order_inner[0]*order_inner[0]-1;
		  }
                if (RT && order_inner[0] > 0)
		  ndof += order_inner[0] + 1;
	      }
	    else
	      {  // quad
		INT<2> p(order_inner[0], order_inner[1]);
		
		int ni = ho_div_free
		  ? p[0]*p[1] 
		  : 2*p[0]*p[1] + p[0] + p[1];
		
		ndof += ni; 
	      }
	  }

        order = 0; 
        for (int i = 0; i < ET_trait<ET>::N_EDGE; i++)
          if (order_facet[i][0] > order)
            order = order_facet[i][0];
        
        for (int j = 0; j < 2; j++)
          if (order_inner[j] > order) 
            order = order_inner[0];

        if (ET != ET_TRIG) order++;
        if(RT) order ++;
      }
    else
      {

        int p = order_inner[0];
        int pc = order_inner[0]; // should be order_inner_curl!!!  
        int pz = p; // order_inner[2];

        if (only_ho_div){
          switch (ET)
            {   
            case ET_TRIG: case ET_QUAD: // for the compiler
              break;
            case ET_TET: 
                ndof = p*(p+1)*(p-1)/6 + p*(p-1)/2 + p-1;
              break;
            case ET_PRISM:
              if (order_inner[0]>0 )
                ndof = (p+1)*(p+2)*(pz+1)/2 - 1;
              /*
                  // inner_dof horizontal
                  ndof += (order_inner[0]+1)*(3*(order_inner[0]-1)+(order_inner[0]-2)*(order_inner[0]-1));
                  // inner dof vertical
                  ndof += (order_inner[0]-1)*(order_inner[0]+1)*(order_inner[0]+2)/2;
              */
              break;
            case ET_HEX:
              ndof = 3*(p+1)*(p+1)*p;
              break; 
            }
        }
        else
        {
          ndof = ET_trait<ET>::N_FACE;

          for(int i = 0; i < ET_trait<ET>::N_FACE; i++)
            {
              INT<2> p = order_facet[i];
              if (ET_trait<ET>::FaceType(i) == ET_TRIG)
                ndof += (p[0]*p[0]+3*p[0])/2;
              else
                ndof +=  p[0]*p[1] + p[0] + p[1];
            }

        // cout << "ndof, bound = " << ndof << endl;
          switch (ET)
            {   
            case ET_TRIG: case ET_QUAD: // for the compiler
              break;

            case ET_TET: 
              if(pc > 1) 
                ndof += pc*(pc+1)*(pc-1)/3 + pc*(pc-1)/2;
              if(p > 1 && !ho_div_free) 
                ndof += p*(p+1)*(p-1)/6 + p*(p-1)/2 + p-1;
	      if(RT && p >= 1)
		ndof += (p+1)*(p+2)/2;
              break;

            case ET_PRISM:
              // SZ: ATTENTION PRISM up to now only using for order_inner[0] !!  
              if (order_inner[0]>0 )
                {
              ndof += (p+2)*p*(pz+1) + (p+1)*(p+2)*pz/2;
              if (ho_div_free)
                ndof -= (p+1)*(p+2)*(pz+1)/2 - 1;
              /*
                  // inner_dof horizontal
                  ndof += (order_inner[0]+1)*(3*(order_inner[0]-1)+(order_inner[0]-2)*(order_inner[0]-1));
                  // inner dof vertical
                  ndof += (order_inner[0]-1)*(order_inner[0]+1)*(order_inner[0]+2)/2;
              */
                }
              break;
            case ET_HEX:
              ndof += 3*(p+1)*(p+1)*p;       
              if (ho_div_free)     
                ndof -= p*p*p+3*p*p+3*p;
              break; 
            }
        }
	// cout << "ndof, tot = " << ndof << endl;

        order = 0; 
        for (int i = 0; i < ET_trait<ET>::N_FACE; i++)
          {
            int p = max2(order_facet[i][0], order_facet[i][1]);
            if (p > order) order = p;
          }

        int pi = max3(order_inner[0], order_inner[1], order_inner[2]);
        if (pi > order) order = pi;

	if (ET != ET_TET) order++;
	
	if(RT) order ++;
      }
  }


  template <ELEMENT_TYPE ET>
  tuple<int,int,int,int> HDivHighOrderFE<ET> :: 
  GetNDofVEFC() const
  {
    int nv=0, ne = 0, nf = 0, nc = 0;
    if (DIM == 2)
      {
        // if (only_ho_div)
	//   {
	//     if (ET == ET_TRIG)
        //       {
        //         ndof = order_inner[0]*(order_inner[0]+1)/2 - 1;
        //         order = order_inner[0];
        //       }
	//     else
        //       {
        //         ndof = order_inner[0]*order_inner[1] + order_inner[0] + order_inner[1];
        //         order = max2(order_inner[0], order_inner[1])+1;
        //       }
	//     return;
	//   }
        // else
	  {
	    ne = ET_trait<ET>::N_EDGE;
	    
	    for(int i = 0; i < ET_trait<ET>::N_EDGE; i++)
	      ne += order_facet[i][0];
	    
	    if (ET == ET_TRIG)
	      {
		if (order_inner[0] > 1)
		  { 
		    if (ho_div_free)
		      nf += order_inner[0]*(order_inner[0]-1)/2;
		    else
		      nf += order_inner[0]*order_inner[0]-1;
		  }
                if (RT && order_inner[0] > 0)
		  nf += order_inner[0] + 1;
	      }
	    else
	      {  // quad
		INT<2> p(order_inner[0], order_inner[1]);
		
		int ni = ho_div_free
		  ? p[0]*p[1] 
		  : 2*p[0]*p[1] + p[0] + p[1];
		
		nf += ni; 
	      }
	  }
      }
    else
      {

        int p = order_inner[0];
        int pc = order_inner[0]; // should be order_inner_curl!!!  
        int pz = p; // order_inner[2];

        // if (only_ho_div){
        //   switch (ET)
        //     {   
        //     case ET_TRIG: case ET_QUAD: // for the compiler
        //       break;
        //     case ET_TET: 
        //         ndof = p*(p+1)*(p-1)/6 + p*(p-1)/2 + p-1;
        //       break;
        //     case ET_PRISM:
        //       if (order_inner[0]>0 )
        //         ndof = (p+1)*(p+2)*(pz+1)/2 - 1;
        //       /*
        //           // inner_dof horizontal
        //           ndof += (order_inner[0]+1)*(3*(order_inner[0]-1)+(order_inner[0]-2)*(order_inner[0]-1));
        //           // inner dof vertical
        //           ndof += (order_inner[0]-1)*(order_inner[0]+1)*(order_inner[0]+2)/2;
        //       */
        //       break;
        //     case ET_HEX:
        //       ndof = 3*(p+1)*(p+1)*p;
        //       break; 
        //     }
        // }
        // else
        {
          nf = ET_trait<ET>::N_FACE;

          for(int i = 0; i < ET_trait<ET>::N_FACE; i++)
            {
              INT<2> p = order_facet[i];
              if (ET_trait<ET>::FaceType(i) == ET_TRIG)
                nf += (p[0]*p[0]+3*p[0])/2;
              else
                nf +=  p[0]*p[1] + p[0] + p[1];
            }

          switch (ET)
            {   
            case ET_TRIG: case ET_QUAD: // for the compiler
              break;

            case ET_TET: 
              if(pc > 1) 
                nc += pc*(pc+1)*(pc-1)/3 + pc*(pc-1)/2;
              if(p > 1 && !ho_div_free) 
                nc += p*(p+1)*(p-1)/6 + p*(p-1)/2 + p-1;
	      if(RT && p >= 1)
		nc += (p+1)*(p+2)/2;
              break;

            case ET_PRISM:
              // SZ: ATTENTION PRISM up to now only using for order_inner[0] !!  
              if (order_inner[0]>0 )
                {
                  nc += (p+2)*p*(pz+1) + (p+1)*(p+2)*pz/2;
                  if (ho_div_free)
                    nc -= (p+1)*(p+2)*(pz+1)/2 - 1;
                }
              break;
            case ET_HEX:
              nc += 3*(p+1)*(p+1)*p;       
              if (ho_div_free)     
                nc -= p*p*p+3*p*p+3*p;
              break; 
            }
        }
      }

    return { nv, ne, nf, nc };
  }


  template <ELEMENT_TYPE ET>
  void HDivHighOrderFE<ET> :: 
  GetFacetDofs(int i, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    dnums.Append (i);

    int base = ET_trait<ET>::N_FACET;
    for (int j = 0; j < ET_trait<ET>::N_FACET; j++)
      {
	int nf = 0;
	switch (ElementTopology::GetFacetType(ET,j))
	  {
	  case ET_SEGM: nf = order_facet[j][0]; break;
	    // case ET_TRIG: nf = (sqr (order_face[0])+3*order_face[0])/2; break;
	  case ET_TRIG: nf = (order_facet[j][0]+1)*(order_facet[j][0]+2)/2-1; break;
	  case ET_QUAD: nf = (order_facet[j][0]+1)*(order_facet[j][1]+1)-1; break;
	  default:
	    throw Exception("what kind of face is that ?");
	  }
	
	if (i == j)
	  {
	    dnums += IntRange (base, base+nf);
	    // cout << "i = " << i << ", dnums = " << dnums << endl;
	    return;
	  }
	base += nf;
      }
    throw Exception("illegal facet index");
  }


  //------------------------------------------------------------------------
  // HDivHighOrderTrig
  //------------------------------------------------------------------------

  
  /*
  HDivHighOrderFE<ET_TRIG> :: HDivHighOrderFE (int aorder)
  {
    order_inner = INT<3>(aorder,0,0);
    for (int i = 0; i < 3; i++)
      order_edge[i] = aorder;
    ComputeNDof();
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
  */


  //------------------------------------------------------------------------
  // HDivHighOrderQuad
  //------------------------------------------------------------------------


  /*

  HDivHighOrderFE<ET_QUAD> :: HDivHighOrderFE (int aorder)
  {
    order_inner = INT<3>(aorder,aorder,0);
    for (int i = 0; i < 4; i++)
      order_edge[i] = aorder;
    ComputeNDof();
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
  */



  //------------------------------------------------------------------------
  // HDivHighOrderTet
  //------------------------------------------------------------------------

  /*
  HDivHighOrderFE<ET_TET> :: HDivHighOrderFE (int aorder)
  {
    order_inner = INT<3>(aorder,aorder,aorder);
    for (int i = 0; i < 4; i++)
      order_face[i] = INT<2>(aorder,aorder);

    ComputeNDof();
  }
  */


  /*
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
  */

  /*
    // in the header 
  
  template<typename Tx, typename TFA>  
  void  HDivHighOrderFE<ET_TET> :: T_CalcShape (Tx hx[], TFA & shape) const
  {
    if (only_ho_div && order_inner[0]<=1) return;
    Tx x = hx[0], y = hx[1], z = hx[2];
    Tx lami[4] = { x, y, z, 1-x-y-z };

    ArrayMem<Tx,10> adpol1(order), adpol2(order), adpol3(order);
	
    int ii = 4; 
    if (!only_ho_div)
    {
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

          // Compatibility WITH TRIG!! 
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
    }
    else
      ii = 0;
    // cell-based shapes 
    int p = order_inner[0];
    int pc = order_inner[0]; // should be order_inner_curl  
    int pp = max2(p,pc); 
    if ( pp >= 2 )
      {
        T_INNERSHAPES::CalcSplitted(pp+2, lami[0]-lami[3], lami[1], lami[2], adpol1, adpol2, adpol3 );
      
        if (!only_ho_div){
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
        }

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
  */





  //------------------------------------------------------------------------
  // HDivHighOrderPrism
  //------------------------------------------------------------------------

 
  
  #ifdef ABC
  
  //SZ : Attention PRISMA has still isotropic inner_order
  void HDivHighOrderFE<ET_PRISM> :: CalcShape (const IntegrationPoint & ip,
                                               FlatMatrixFixWidth<3> shape) const
  {
    // shape = 0.0;
    T_HDivHighOrderFiniteElement<ET_PRISM>::CalcShape (ip, shape);
    return;
    // *testout << "t_shape = " << endl << shape << endl;

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
    // const EDGE * edges = ElementTopology::GetEdges (ET_PRISM);


    ArrayMem<AutoDiff<2>,20> adpolxy1(order+1), adpolxy2(order+1), adpolx(order+1), adpoly(order+1);
    ArrayMem<AutoDiff<1>,20> adpolz(order+2);

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


        int pp = max2(p[0],p[1]); 
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


    if(p>=1)
      {
	//ArrayMem<double,10> rec_pol(order-1), drec_polx(order-1),  drec_polt(order-1);
	//ArrayMem<double,10> rec_pol2(order-1), drec_pol2x(order-1), drec_pol2t(order-1);

	// HACK, performed by MW, authorized by SZ
	ArrayMem<double,10> rec_pol(max2(order,p+2)), drec_polx(max2(order,p+2)),  drec_polt(max2(order,p+2));
	ArrayMem<double,10> rec_pol2(max2(order,p+2)), drec_pol2x(max2(order,p+2)), drec_pol2t(max2(order,p+2));
	ArrayMem<double,30> mem(3*(order+1)*sizeof(double));
	ArrayMem<double,30> mem2(3*(order+1)*sizeof(double));
	FlatMatrixFixWidth<3> curl_recpol(order+1, &mem[0]);
	FlatMatrixFixWidth<3> curl_recpol2(order+1, &mem2[0]);

	int fav[3] = {0, 1, 2};


	double ls = x.Value();
	double le = y.Value();
	double lo = 1-x.Value()-y.Value();




	ScaledLegendrePolynomialandDiff(p+1, le-ls, 1-lo,
					rec_pol, drec_polx, drec_polt);
	LegendrePolynomialandDiff(p+1, 2*lo-1,rec_pol2, drec_pol2x);
	LegendrePolynomial(p+2,2*z-1,adpolz);
	// \curl (ls * le * rec_pol), and \curl (lo rec_pol2)
	for (j = 0; j <= p; j++)
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
	  for (j = 0; j <= p-1; j++)
	    for (k = 0; k <= p-1-j; k++, ii++)
	      for (l = 0; l < 2; l++)
		shape(ii, l) = (ls*le*rec_pol[j]*curl_recpol2(k,l) + curl_recpol(j,l)*lo*rec_pol2[k])*adpolz[i].Value();

	// rotations of curls
	for (i = 0; i <= p; i++)
	  for (j = 0; j <= p-1; j++)
	    for (k = 0; k <= p-1-j; k++, ii++)
	      for (l = 0; l < 2; l++)
		shape(ii, l) = (ls*le*rec_pol[j]*curl_recpol2(k,l) - curl_recpol(j,l) * lo * rec_pol2[k])*adpolz[i].Value();

	// rec_pol2 * RT_0
	for (i = 0; i <= p; i++)
	  for (j = 0; j <= p-1; j++, ii++)
	    for (l = 0; l < 2; l++)
	      shape(ii,l) = (lo*rec_pol2[j] * (ls*clami(fav[1],l)-le*clami(fav[0],l)))*adpolz[i].Value();
      }

    // *testout << "horiz, ii = " << ii << endl;

    if (p>=1)
      {
	T_ORTHOPOL::Calc (p+1, 2*z-1, adpolz);

	ScaledLegendrePolynomial (p+1, x, 1-y, adpolx);
	LegendrePolynomial (p+1, 2*y-1, adpoly);
	// const. lowest-order * IntLegendre + 2 X (linear * IntLegendre)
	for( i = 0; i <= p-1; i++)
	  for ( j = 0; j <= p; j++)
	    for (k = 0; k <= p-j; k++, ii++)
	      shape(ii,2) = adpolx[j].Value()*adpoly[k].Value()*adpolz[i].Value();
      }


    // *testout << "total, ii = " << ii << endl;
    // *testout << "orig shape = " << endl << shape << endl;

    T_HDivHighOrderFiniteElement<ET_PRISM>::CalcShape (ip, shape);
  }


  void HDivHighOrderFE<ET_PRISM> :: CalcDivShape (const IntegrationPoint & ip,
						       FlatVector<> divshape) const
  {
    // HDivFiniteElement<3>::CalcDivShape (ip, divshape);
    // *testout << "num divshape = " << endl << divshape << endl;

    T_HDivHighOrderFiniteElement<ET_PRISM>::CalcDivShape (ip, divshape);
    // *testout << "some exact divshape = " << endl << divshape << endl;
  }
#endif



#ifdef GONE_TO_THE_HEADER

  //in the header, due to undefined reference problems (ICC)

  template<typename Tx, typename TFA>  
  void  HDivHighOrderFE<ET_PRISM> :: T_CalcShape (Tx hx[], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];

    AutoDiff<3> lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    AutoDiff<3> muz[6]  = { 1-z, 1-z, 1-z, z, z, z };
       
    const FACE * faces = ElementTopology::GetFaces (ET_PRISM); 

    ArrayMem<AutoDiff<3>,20> adpolxy1(order+4),adpolxy2(order+4); 
    ArrayMem<AutoDiff<3>,20> adpolz(order+4);   


    ArrayMem<Tx,10> adpol1(order), adpol2(order), adpol3(order);
    
    // trig faces

    int ii = 5;
    for (int i = 0; i < 2; i++)
      {
	int p = order_face[i][0];
	int fav[3] = {faces[i][0], faces[i][1], faces[i][2]};

	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	  	
	
	shape[i] = wDu_Cross_Dv<3> (lami[fav[0]], lami[fav[1]], muz[fav[0]]);

	AutoDiff<3> xi = lami[fav[1]]-lami[fav[0]];
	AutoDiff<3> eta = lami[fav[2]];

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
      }    
    
    // quad faces
    for (int i = 2; i < 5; i++)
      {
	INT<2> p = order_face[i];
	 
	int fmax = 0;
	for (int j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;

	int fz = 3-fmax; 
	int ftrig = fmax^1; 
	int f = faces[i][fmax];
	int f1 = faces[i][ftrig];
	int f2 = faces[i][fz];

	AutoDiff<3> xi = lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
	AutoDiff<3> eta = 1-lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
	AutoDiff<3> zeta = muz[faces[i][fmax]]-muz[faces[i][fz]]; 
	
	int pp = int(max2(p[0],p[1]))+1;
	T_ORTHOPOL::CalcTrigExt(pp,xi,eta,adpolxy1); 
	T_ORTHOPOL::Calc(pp,zeta,adpolz); 

	double fac = (vnums[faces[i][fz]] > vnums[faces[i][ftrig]]) ? 1 : -1;

	shape[i] = uDvDw_minus_DuvDw<3> (lami[faces[i][fmax]],
					 lami[faces[i][ftrig]], -0.5*fac*zeta);
					    

	if (vnums[f1] > vnums[f2])
	  {
	    for (int k = 0; k <= p[0]-1; k++)
	      for (int j = 0; j <= p[1]-1; j++, ii++)
		shape[ii] = Du_Cross_Dv<3> (adpolxy1[k], -2*adpolz[j]);
	  }
	else
	  {
	    for (int j = 0; j <= p[1]-1; j++)
	      for (int k = 0; k <= p[0]-1; k++, ii++)
		shape[ii] = Du_Cross_Dv<3> (adpolxy1[k],  2*adpolz[j]);
	  }
	  

        if (vnums[f1] > vnums[f2])
          {
            for (int j= 0; j <= p[0]-1; j++, ii++)
	      shape[ii] = Du_Cross_Dv<3> (adpolxy1[j], zeta);
            for(int j=0; j<= p[1]-1;j++,ii++)
	      shape[ii] = curl_uDvw_minus_Duvw<3> (lami[f1], lami[f], 2*adpolz[j]);
          }  
        else
          {
            for(int j = 0; j <= p[0]-1; j++,ii++)
	      shape[ii] = curl_uDvw_minus_Duvw<3> (lami[f1], lami[f], 2*adpolz[j]);
            for (int j= 0; j <= p[1]-1; j++, ii++)
	      shape[ii] = Du_Cross_Dv<3> (adpolxy1[j], zeta);
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
	      shape[ii++] = Du_Cross_Dv<3> (adpolxy1[i],adpolxy2[j]*adpolz[k]);

	for(int i=0;i<=p-1;i++)
	  for(int j=0;j<=p-1-i;j++)
	    for(int k=0;k<=pz-1;k++)
	      shape[ii++] = curl_uDvw_minus_Duvw<3> (adpolxy1[i],adpolxy2[j],adpolz[k]);

	for(int j=0;j<=p-1;j++) 
	  for (int k=0;k<=pz-1;k++) 
            shape[ii++] = curl_uDvw_minus_Duvw<3> (x,y, adpolxy2[j]*adpolz[k]);

    	for(int i = 0; i <= p-1; i++) 
	  for(int j = 0; j <= p-1-i; j++) 
            shape[ii++] = curl_uDvw_minus_Duvw<3> (z,1-z, adpolxy1[i]*adpolxy2[j]);

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
                  shape[ii++] = wDu_Cross_Dv<3> ((x-y)*adpolxy1[i], x*adpolxy2[j], z*(1-z)*adpolz[k]);

            for (int i = 0; i < p; i++)
              for (int j = 0; j < p-i; j++)
                shape[ii++] = wDu_Cross_Dv<3> (z, x*y*adpolxy1[i], (1-x-y)*adpolxy2[j]);

            for (int i = 0; i < p; i++)
              shape[ii++] = wDu_Cross_Dv<3> (z, x, y*(1-x-y)*adpolxy1[i]);


            /*
            for (int i = 0; i <= p-1; i++)
              for (int k = 0; k <= pz; k++)
                shape[ii++] = wDu_Cross_Dv<3> (x*y*adpolxy1[i], z*adpolz[k],  1-x-y);
            */
          }
      }

    if (ii != ndof) cout << "hdiv-prism: dofs missing, ndof = " << ndof << ", ii = " << ii << endl;
  }
#endif
 


  /*
  void HDivHighOrderFE<ET_PRISM> :: GetFacetDofs(int fa, Array<int> & dnums) const 
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
        INT<2> p = order_face[fa];
        nf = p[0]*p[1]+p[0]+p[1];
      }
      
    for (int i=0; i<nf; i++)
      dnums.Append(i+base);
  }   
  */


  //------------------------------------------------------------------------
  // HDivHighOrderHex
  //------------------------------------------------------------------------


#ifdef HDIVHEX

  HDivHighOrderFE<ET_HEX> :: HDivHighOrderFE (int aorder)
    : HDivHighOrderFiniteElement<3>(ET_HEX)
  {
    order_inner = aorder;
    for (int i = 0; i < 6; i++)
      order_face[i] = aorder;
    ComputeNDof();
  }


  void HDivHighOrderFE<ET_HEX> :: ComputeNDof()
  {
    ndof = 6; // RT_0
   
    for (int i = 0; i < 6; i++)
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

    order = 0; // max2(order_face_normal,order_inner);
    for (int i = 0; i < 6; i++)
      {
        int pp = max2(order_face[i][0],order_face[i][1]);
	if (pp > order)
	  order = pp;

      }
    int pp = max2(order_inner[0],max2(order_inner[1],order_inner[2]));
    order=max2(pp,order); 

    if (order == 0) order = 1;
    order++; // integration order
  }


  void HDivHighOrderFE<ET_HEX> :: CalcShape (const IntegrationPoint & ip,
                                             SliceMatrix<> shape) const
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

  void HDivHighOrderFE<ET_HEX> :: CalcDivShape (const IntegrationPoint & ip,
                                                SliceVector<> divshape) const
  {
    int i, j, k; 
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

  
  void HDivHighOrderFE<ET_HEX> :: GetFacetDofs(int fa, Array<int> & dnums) const 
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
    INT<2> p = order_face[fa];
    int nf = p[0]*p[1]+p[0]+p[1];
         
    for (int i=0; i<nf; i++)
      dnums.Append(base+i);
  }   
#endif


  // template class  HDivHighOrderFiniteElement<2>;
  // template class  HDivHighOrderFiniteElement<3>;
  template class  HDivHighOrderNormalFiniteElement<1>;
  template class  HDivHighOrderNormalFiniteElement<2>;

  /*
  template class T_HDivHighOrderFiniteElement<ET_TRIG>;
  template class T_HDivHighOrderFiniteElement<ET_QUAD>;
  template class T_HDivHighOrderFiniteElement<ET_TET>;
  template class T_HDivHighOrderFiniteElement<ET_PRISM>;
  */




  template class HDivHighOrderFEFO<ET_TRIG,1>;
  template class HDivHighOrderFEFO<ET_TRIG,2>;
  template class HDivHighOrderFEFO<ET_TRIG,3>;
  template class HDivHighOrderFEFO<ET_TRIG,4>;
  template class HDivHighOrderFEFO<ET_TRIG,5>;
  template class HDivHighOrderFEFO<ET_TRIG,6>;

}




