#ifndef FILE_L2HOFE_IMPL
#define FILE_L2HOFE_IMPL

/*********************************************************************/
/* File:   l2hofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/



namespace ngfem
{





  template <>
  class L2HighOrderFE_Shape<ET_POINT> : public L2HighOrderFE<ET_POINT>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const;
  };

  template <>
  class L2HighOrderFE_Shape<ET_SEGM> : public L2HighOrderFE<ET_SEGM>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const;
  };

  template <> 
  class L2HighOrderFE_Shape<ET_TRIG> : public L2HighOrderFE<ET_TRIG>
  {
  public:
    template<typename Tx, typename TFA> 
    void T_CalcShape (Tx x[], TFA & shape) const;
  };
  
  template <> 
  class L2HighOrderFE_Shape<ET_QUAD> : public L2HighOrderFE<ET_QUAD>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[], TFA & shape) const;
  };

  template <> 
  class L2HighOrderFE_Shape<ET_TET> : public L2HighOrderFE<ET_TET>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const;
  };

  template <> 
  class L2HighOrderFE_Shape<ET_PRISM> : public L2HighOrderFE<ET_PRISM>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[], TFA & shape) const;
  };

  template <> 
  class L2HighOrderFE_Shape<ET_PYRAMID> : public L2HighOrderFE<ET_PYRAMID>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const;
  };

  template <> 
  class L2HighOrderFE_Shape<ET_HEX> : public L2HighOrderFE<ET_HEX>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const;
  };









  /* *********************** Point  **********************/
  

  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_POINT> ::
  T_CalcShape (Tx x[], TFA & shape) const
  {
    shape[0] = 1.0;
  }



  /* *********************** Segment  **********************/
  

  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_SEGM> ::
  T_CalcShape (Tx x[], TFA & shape) const
  {
    Tx lam[2] = { x[0], 1-x[0] };
    INT<2> e = GetEdgeSort (0, vnums);
    LegendrePolynomial (order, lam[e[1]]-lam[e[0]], shape);
  }



  /* *********************** Triangle  **********************/


  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_TRIG> ::
  T_CalcShape (Tx x[], TFA & shape) const
  {
    Tx lam[3] = { x[0], x[1], 1-x[0]-x[1] };
    INT<4> f = GetFaceSort (0, vnums);
    int p = order_inner[0];
    DubinerBasis::Eval (p, lam[f[0]], lam[f[1]], shape);
  }


  /* *********************** Quad  **********************/

  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_QUAD> ::
  T_CalcShape (Tx hx[], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1];
    
    Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    
    INT<4> f = GetFaceSort (0, vnums);  
    
    Tx xi = sigma[f[0]]-sigma[f[1]]; 
    Tx eta = sigma[f[0]]-sigma[f[3]]; 
    
    int n = max(order_inner[0],order_inner[1]);
    ArrayMem<Tx, 20> polx(n+1), poly(n+1);
      
    LegendrePolynomial (n, xi, polx);
    LegendrePolynomial (n, eta, poly);
      
    for (int i = 0, ii = 0; i <= order_inner[0]; i++)
      for (int j = 0; j <= order_inner[1]; j++)
        shape[ii++] = polx[i] * poly[j];
  }


  /* *********************** Tet  **********************/

  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_TET> ::
  T_CalcShape (Tx x[], TFA & shape) const
  {
    Tx lami[4] = { x[0], x[1], x[2], 1-x[0]-x[1]-x[2] };
    
    int sort[4];
    for (int i = 0; i < 4; i++) sort[i] = i;
    
    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
    if (vnums[sort[2]] > vnums[sort[3]]) Swap (sort[2], sort[3]);
    if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
    if (vnums[sort[1]] > vnums[sort[3]]) Swap (sort[1], sort[3]);
    if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);

    Tx lamis[4];
    for (int i = 0; i < 4; i++)
      lamis[i] = lami[sort[i]];

    ArrayMem<Tx, 20> memx(sqr(order+1));
    ArrayMem<Tx, 20> memy(sqr(order+1));

    FlatMatrix<Tx> polsx(order+1, &memx[0]);
    FlatMatrix<Tx> polsy(order+1, &memy[0]);
    VectorMem<10, Tx> polsz(order+1);
    
    for (int i = 0; i <= order; i++)
      JacobiPolynomial (order, 2*lamis[0]-1, 2*i+2, 0, polsx.Row(i));
    for (int i = 0; i <= order; i++)
      ScaledJacobiPolynomial (order, lamis[1]-lamis[2]-lamis[3], 1-lamis[0], 2*i+1, 0, polsy.Row(i));

    ScaledLegendrePolynomial (order, lamis[2]-lamis[3], lamis[2]+lamis[3], polsz);

    for (int i = 0, ii = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        for (int k = 0; k <= order-i-j; k++, ii++)
          shape[ii] = polsz(k) * polsy(k, j) * polsx(j+k, i);
  }





  /* *********************** Prism  **********************/




  template<typename Tx, typename TFA>  
  void  L2HighOrderFE_Shape<ET_PRISM> ::
  T_CalcShape (Tx hx[], TFA & shape) const
  {
    Tx lami[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };

    int sort[3];
    for (int i = 0; i < 3; i++) sort[i] = i;
    
    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
    if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);

    Tx lamis[3];
    for (int i = 0; i < 3; i++)
      lamis[i] = lami[sort[i]];

    Tx x = lamis[0];
    // Tx y = lamis[1];
    Tx z = hx[2];

    int p=order_inner[0];
    int q=order_inner[1];

    ArrayMem<Tx, 20> memx(sqr(p+1));
    FlatMatrix<Tx> polsx(p+1, &memx[0]);

    VectorMem<10, Tx> polsy(p+1);
    VectorMem<10, Tx> polsz(q+1);
    
    for (int i = 0; i <= p; i++)
      JacobiPolynomial (p, 2*x-1, 2*i+1, 0, polsx.Row(i));

    ScaledLegendrePolynomial (order, lamis[1]-lamis[2], lamis[1]+lamis[2], polsy);
    LegendrePolynomial (order, 2*z-1, polsz);

    int ii = 0;
    for (int k = 0; k <= q; k++)
      for (int i = 0; i <= p; i++)
        for (int j = 0; j <= p-i; j++)
          shape[ii++] = polsx(j,i) * polsy(j) * polsz(k);
  }




  /* *********************** Pyramid  **********************/


  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_PYRAMID> :: 
  T_CalcShape (Tx hx[], TFA & shape) const
  {
    Tx x = hx[0];
    Tx y = hx[1];
    Tx z = hx[2];

    // if (z == 1) z -= 1e-8;
    z *= (1-1e-8);
    Tx xt = 2 * (x / (1-z)) - 1;
    Tx yt = 2 * (y / (1-z)) - 1;

    VectorMem<10, Tx> polsx(order+1);
    VectorMem<10, Tx> polsy(order+1);
    
    ArrayMem<Tx, 20> memz(sqr(order+1));
    FlatMatrix<Tx> polsz(order+1, &memz[0]);

    Tx fac = 1.0;
    for (int i = 0; i <= order; i++)
      {
	JacobiPolynomial (order, 2*z-1, 2*i+2, 0, polsz.Row(i));
	polsz.Row(i) *= fac;
	fac *= (1-z);
      }

    LegendrePolynomial (order, xt, polsx);
    LegendrePolynomial (order, yt, polsy);
    
    int ii = 0;
    for (int iz = 0; iz <= order; iz++)
      for (int ix = 0; ix <= order-iz; ix++)
	for (int iy = 0; iy <= order-iz; iy++, ii++)
	  shape[ii] = polsx(ix) * polsy(iy) * polsz(max(ix,iy), iz);
  }


  /* *********************** Hex  **********************/


  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_HEX> :: T_CalcShape (Tx hx[], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];
    
    // no orientation necessary
    int p=order_inner[0];
    int q=order_inner[1];
    int r=order_inner[2];
    
    ArrayMem<Tx, 20> polx(p+1), poly(q+1), polz(r+1);

    LegendrePolynomial (p, 2*x-1, polx);
    LegendrePolynomial (q, 2*y-1, poly);
    LegendrePolynomial (r, 2*z-1, polz);
  
    for (int i = 0, ii = 0; i <= p; i++)
      for (int j = 0; j <= q; j++)
        for (int k = 0; k <= r; k++)
          shape[ii++] = polx[i] * poly[j] * polz[k];
  }



}

#endif
