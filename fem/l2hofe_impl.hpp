#ifndef FILE_L2HOFE_IMPL
#define FILE_L2HOFE_IMPL

/*********************************************************************/
/* File:   l2hofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/



namespace ngfem
{


#ifndef __CUDA_ARCH__
  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  typename L2HighOrderFE<ET,SHAPES,BASE>::TPRECOMP L2HighOrderFE<ET,SHAPES,BASE>::precomp;

  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  typename L2HighOrderFE<ET,SHAPES,BASE>::TPRECOMP_TRACE L2HighOrderFE<ET,SHAPES,BASE>::precomp_trace(320);

  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  typename L2HighOrderFE<ET,SHAPES,BASE>::TPRECOMP_GRAD L2HighOrderFE<ET,SHAPES,BASE>::precomp_grad(40);
#endif





  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  PrecomputeTrace ()
  {
#ifndef __CUDA_ARCH__
    for (int f = 0; f < ElementTopology::GetNFacets(ET); f++)
      {
        int classnr =  ET_trait<ET>::GetFacetClassNr (f, vnums);
        if (precomp_trace.Used (INT<2> (order, classnr)))
          continue;
        
        ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (ET, f);
        int nf;
        switch (etfacet)
          {
          case ET_POINT: nf = 1; break;
          case ET_SEGM: nf = order+1; break;
          case ET_TRIG: nf = (order+1)*(order+2)/2; break;
          case ET_QUAD: nf = sqr(order+1); break;
          default: nf = 0;
          }

        Matrix<> * trace = new Matrix<>(nf, ndof);
        DGFiniteElement<ET>::CalcTraceMatrix (f, *trace);
        precomp_trace.Set (INT<2> (order, classnr), trace);
      }
#endif
  }
  
  
  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  PrecomputeGrad ()
  {
#ifndef __CUDA_ARCH__
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    
    if (precomp_grad.Used (INT<2> (order, classnr)))
      return;

    Matrix<> * gmat = new Matrix<>(ndof*DIM, ndof);
    DGFiniteElement<ET>::CalcGradientMatrix (*gmat);
    precomp_grad.Set (INT<2> (order, classnr), gmat);
#endif
  }
    


  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  PrecomputeShapes (const IntegrationRule & ir) 
  {
#ifndef __CUDA_ARCH__
    int classnr =  ET_trait<ET>::GetClassNr (vnums);

    if (precomp.Get (classnr, order, ir.GetNIP())) return;
    
    PrecomputedScalShapes<DIM> * pre = new  PrecomputedScalShapes<DIM> (ir.GetNIP(), ndof);
    
    MatrixFixWidth<DIM> dshapes(ndof);
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	this->CalcShape (ir[i], pre->shapes.Row(i));
	this->CalcDShape (ir[i], dshapes);
	pre->dshapes.Rows (DIM*i, DIM*(i+1)) = Trans (dshapes);
      }
    
    precomp.Add (classnr, order, ir.GetNIP(), pre);
#endif
  }
  


  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  Evaluate (const IntegrationRule & ir, BareSliceVector<double> coefs, FlatVector<double> vals) const
  {
#ifndef __CUDA_ARCH__
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());
    if (pre)
      vals = pre->shapes * coefs;
    else
#endif
      BASE :: Evaluate (ir, coefs, vals);
  }

  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  EvaluateTrans (const IntegrationRule & ir, FlatVector<> values, BareSliceVector<> coefs) const
  {
#ifndef __CUDA_ARCH__
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());

    if (pre)
      coefs.Range(0,ndof) = Trans(pre->shapes)*values;
    else
#endif
      BASE :: EvaluateTrans (ir, values, coefs);
  }

  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  EvaluateGrad (const IntegrationRule & ir, BareSliceVector<> coefs, FlatMatrixFixWidth<DIM> values) const
    {
#ifndef __CUDA_ARCH__
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
	{
	  FlatVector<> vval(DIM*values.Height(), &values(0,0));
	  vval = pre->dshapes * coefs;
	}
      else
#endif
	BASE :: EvaluateGrad (ir, coefs, values);
    }

  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> values, BareSliceVector<> coefs) const
    {
#ifndef __CUDA_ARCH__
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
	coefs.Range(0,ndof) = Trans (pre->dshapes) * FlatVector<> (DIM*ndof, &values(0,0));  // values.Height !!!
      else
#endif
	BASE :: EvaluateGradTrans (ir, values, coefs);
    }


  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  GetGradient (FlatVector<> coefs, FlatMatrixFixWidth<DIM> grad) const
  {
#ifndef __CUDA_ARCH__
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    int bnr, pos;
    if (precomp_grad.Used (INT<2> (order, classnr), bnr, pos))
      {
        FlatMatrix<> gmat = *precomp_grad.Get (bnr, pos);
        FlatVector<> vgrad(grad.Height()*DIM, &grad(0,0));
        // vgrad = gmat * coefs;
        MultMatVec (gmat, coefs, vgrad);
      }
    else
#endif
      DGFiniteElement<ET>::GetGradient (coefs, grad);
  }
  
  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  GetGradientTrans (FlatMatrixFixWidth<DIM> grad, FlatVector<> coefs) const
  {
#ifndef __CUDA_ARCH__
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    int bnr, pos;
    if (precomp_grad.Used (INT<2> (order, classnr), bnr, pos))
      {
        FlatMatrix<> gmat = *precomp_grad.Get (bnr, pos);
        FlatVector<> vgrad(grad.Height()*DIM, &grad(0,0));
        // coefs = Trans(gmat) * vgrad;
        MultMatTransVec (gmat, vgrad, coefs);        
      }
    else
#endif
      DGFiniteElement<ET>::GetGradientTrans (grad, coefs);
  }









  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const
    {
#ifndef __CUDA_ARCH__
      int classnr =  ET_trait<ET>::GetFacetClassNr (facet, vnums);
      int bnr, pos;
      if (precomp_trace.Used (INT<2> (order, classnr), bnr, pos))
	{
	  FlatMatrix<> trace = *precomp_trace.Get (bnr, pos);
	  // fcoefs = trace * coefs;
          MultMatVec (trace, coefs, fcoefs);
	}
      else
#endif
	DGFiniteElement<ET>::GetTrace (facet, coefs, fcoefs);
    }

  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const
  {
#ifndef __CUDA_ARCH__
      int classnr =  ET_trait<ET>::GetFacetClassNr (facet, vnums);
      if (precomp_trace.Used (INT<2> (order, classnr)))
	{
	  FlatMatrix<> trace = *precomp_trace.Get (INT<2> (order, classnr));          
	  // coefs = Trans(*precomp_trace.Get (INT<2> (order, classnr))) * fcoefs;
          MultMatTransVec (trace, fcoefs, coefs);
	}
      else
#endif
	DGFiniteElement<ET>::GetTraceTrans (facet, fcoefs, coefs);
    }


  template <ELEMENT_TYPE ET, class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
    switch (ET)
      {
      case ET_POINT:
        mass(0) = 1;
        break;
        
      case ET_SEGM:
        for (int ix = 0; ix <= order; ix++)
          mass(ix) = 1.0 / (2 * ix + 1);
        break;
        
      case ET_TRIG:
        for (int ix = 0, ii = 0; ix <= order; ix++)
          for (int iy = 0; iy <= order - ix; iy++, ii++)
            mass(ii) = 1.0 / ((2 * ix + 1) * (2 * ix + 2 * iy + 2));
        break;

      case ET_QUAD:
        for (int ix = 0, ii = 0; ix <= order; ix++)
          for (int iy = 0; iy <= order; iy++, ii++)
            mass(ii) = 1.0 / ((2 * ix + 1) * (2 * iy + 1));
        break;

      case ET_TET:
        for (int ix = 0, ii = 0; ix <= order; ix++)
          for (int iy = 0; iy <= order - ix; iy++)
            for (int iz = 0; iz <= order - ix-iy; iz++, ii++)
              mass(ii) = 1.0 / ((2 * ix + 1) * (2 * ix + 2 * iy + 2) * (2 * ix + 2 * iy + 2 * iz + 3));
        break;
        
      case ET_HEX:
        for (int ix = 0, ii = 0; ix <= order; ix++)
          for (int iy = 0; iy <= order; iy++)
            for (int iz = 0; iz <= order; iz++, ii++)
              mass(ii) = 1.0 / ((2 * ix + 1) * (2 * iy + 1) * (2 * iz + 1));
        break;
        
      default:
        DGFiniteElement<ET>::GetDiagMassMatrix (mass);
      }
  }


  









  template <ELEMENT_TYPE ET> 
  class L2HighOrderFE_Shape : public L2HighOrderFE<ET>
  {
    using L2HighOrderFE<ET>::vnums;
    using L2HighOrderFE<ET>::order;
    using L2HighOrderFE<ET>::order_inner;
    using L2HighOrderFE<ET>::GetFaceSort;
    using L2HighOrderFE<ET>::GetEdgeSort;
  public:
    // template<typename Tx, typename TFA>  
    // INLINE void T_CalcShape (Tx hx[], TFA & shape) const;

    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (TIP<ET_trait<ET>::DIM,Tx> ip, TFA & shape) const;

    template<typename Tx, typename TFA>  
    INLINE void T_CalcDualShape (TIP<ET_trait<ET>::DIM,Tx> ip, TFA & shape) const
    {
      if (ip.vb == VOL)
        T_CalcShape (ip, shape);
    }
    
    /*
    {
      throw Exception (string("TIP not implemented, fe = ")+typeid(*this).name());
    }
    */
  };




  /* *********************** Point  **********************/
  

  template<> template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_POINT> ::
  T_CalcShape (TIP<0,Tx> ip, TFA & shape) const
  {
    shape[0] = Tx(1.0);
  }



  /* *********************** Segment  **********************/
  

  template<> template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_SEGM> ::
  T_CalcShape (TIP<1,Tx> ip, TFA & shape) const
  {
    Tx lam[2] = { ip.x, 1-ip.x };
    INT<2> e = GetEdgeSort (0, vnums);
    LegendrePolynomial (order, lam[e[1]]-lam[e[0]], shape);
  }



  /* *********************** Triangle  **********************/

  /*
  template<> template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_TRIG> ::
  T_CalcShape (Tx x[], TFA & shape) const
  {
#ifndef __CUDA_ARCH__
    Tx lam[3] = { x[0], x[1], 1-x[0]-x[1] };
    INT<4> f = GetFaceSort (0, vnums);
    size_t p = order_inner[0];
    DubinerBasis3::Eval (p, lam[f[0]], lam[f[1]], shape);
#else
    int p = order_inner[0];
    DubinerBasis3::Eval (p, x[0], x[1], shape);
#endif
  }
  */


  template<> template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_TRIG> ::
  T_CalcShape (TIP<DIM,Tx> ip, TFA & shape) const
  {
    Tx lam[3] = { ip.x, ip.y, 1-ip.x-ip.y };
    INT<4> f = GetFaceSort (0, vnums);
    size_t p = order_inner[0];
    DubinerBasis::Eval (p, lam[f[0]], lam[f[1]], shape);

    /*
    Tx lam[3] = { ip.x, ip.y, 1-ip.x-ip.y };
    size_t p = order_inner[0];
    DubinerBasis3::Eval (p, lam[loclam[0]], lam[loclam[1]], shape);
    */
    /*
    size_t p = order_inner[0];
    DubinerBasis3::Eval (p, ip.x, ip.y, shape);    
    */
  }



  
  /* *********************** Quad  **********************/

  template<> template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_QUAD> ::
  T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x, y = ip.y;
    Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    
    INT<4> f = GetFaceSort (0, vnums);  
    
    Tx xi = sigma[f[0]]-sigma[f[1]]; 
    Tx eta = sigma[f[0]]-sigma[f[3]]; 

    int p=order_inner[0];
    int q=order_inner[1];

    STACK_ARRAY(Tx, mem, p+q+2);
    Tx * polx = &mem[0];
    Tx * poly = &mem[p+1];
      
    LegendrePolynomial (p, xi, polx);
    LegendrePolynomial (q, eta, poly);

    for (size_t i = 0, ii = 0; i <= p; i++)
      for (size_t j = 0; j <= q; j++)
        shape[ii++] = polx[i] * poly[j];
  }


  /* *********************** Tet  **********************/

  template<> template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_TET> ::
  T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx lami[4] = { ip.x, ip.y, ip.z, 1-ip.x-ip.y-ip.z };

    unsigned char sort[4] = { 0, 1, 2, 3 };
    
    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
    if (vnums[sort[2]] > vnums[sort[3]]) Swap (sort[2], sort[3]);
    if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
    if (vnums[sort[1]] > vnums[sort[3]]) Swap (sort[1], sort[3]);
    if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
    
    Tx lamis[4];
    for (int i = 0; i < 4; i++)
      lamis[i] = lami[sort[i]];

    /*
    ArrayMem<Tx, 20> memx(sqr(order+1));
    ArrayMem<Tx, 20> memy(sqr(order+1));
    VectorMem<10, Tx> polsz(order+1);
    FlatMatrix<Tx> polsx(order+1, &memx[0]);
    FlatMatrix<Tx> polsy(order+1, &memy[0]);

    for (int i = 0; i <= order; i++)
      {
        JacobiPolynomialAlpha jac(2*i+2);
        jac.Eval (order, 2*lamis[0]-1, polsx.Row(i));
      }

    for (int i = 0; i <= order; i++)
      {
        JacobiPolynomialAlpha jac(2*i+1);
        jac.EvalScaled (order, lamis[1]-lamis[2]-lamis[3], 1-lamis[0], polsy.Row(i));
      }

    // ScaledLegendrePolynomial (order, lamis[2]-lamis[3], lamis[2]+lamis[3], polsz);

    for (int i = 0, ii = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        for (int k = 0; k <= order-i-j; k++, ii++)
          shape[ii] = polsz(k) * polsy(k, j) * polsx(j+k, i);
    */

    DubinerBasis3D::Eval (this->order, lamis[0], lamis[1], lamis[2], shape);

    /*
    size_t ii = 0;
    size_t order = this->order;
    LegendrePolynomial leg;
    JacobiPolynomialAlpha jac1(1);    
    leg.EvalScaled1Assign 
      (order, lamis[2]-lamis[3], lamis[2]+lamis[3],
       SBLambda ([&](size_t k, Tx polz) LAMBDA_INLINE
                 {
                   // JacobiPolynomialAlpha jac(2*k+1);
                   JacobiPolynomialAlpha jac2(2*k+2);
 
                   jac1.EvalScaledMult1Assign
                     (order-k, lamis[1]-lamis[2]-lamis[3], 1-lamis[0], polz, 
                      SBLambda ([&] (size_t j, Tx polsy) LAMBDA_INLINE
                                {
                                  // JacobiPolynomialAlpha jac(2*(j+k)+2);
                                  jac2.EvalMult(order - k - j, 2 * lamis[0] - 1, polsy, 
                                                SBLambda([&](size_t j, Tx val) LAMBDA_INLINE
                                                         {
                                                           shape[ii] = val; 
                                                           ii++;
                                                         }));
                                  jac2.IncAlpha2();
                                }));
                   jac1.IncAlpha2();
                 }));
    */
  }





  /* *********************** Prism  **********************/




  template<> template<typename Tx, typename TFA>  
  void  L2HighOrderFE_Shape<ET_PRISM> ::
  T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx lami[3] = { ip.x, ip.y, 1-ip.x-ip.y }; // hx[0], hx[1], 1-hx[0]-hx[1] };

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
    Tx z = ip.z; // hx[2];

    int p=order_inner[0];
    int q=order_inner[1];

    ArrayMem<Tx, 20> memx(sqr(p+1));
    FlatMatrix<Tx> polsx(p+1, &memx[0]);

    VectorMem<10, Tx> polsy(p+1);
    VectorMem<10, Tx> polsz(q+1);

    for (int i = 0; i <= p; i++)
      // JacobiPolynomial (p, 2*x-1, 2*i+1, 0, polsx.Row(i));
      JacobiPolynomialAlpha (p, 2*x-1, 2*i+1, polsx.Row(i));

    // ScaledLegendrePolynomial (order, lamis[1]-lamis[2], lamis[1]+lamis[2], polsy);
    LegendrePolynomial::EvalScaled (p, lamis[1]-lamis[2], lamis[1]+lamis[2], polsy);
    LegendrePolynomial (q, 2*z-1, polsz);

    int ii = 0;
    for (int k = 0; k <= q; k++)
      for (int i = 0; i <= p; i++)
        for (int j = 0; j <= p-i; j++)
          shape[ii++] = polsx(j,i) * polsy(j) * polsz(k);
  }




  /* *********************** Pyramid  **********************/


  template<> template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_PYRAMID> :: 
  T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x;
    Tx y = ip.y;
    Tx z = ip.z;

    // if (z == 1) z -= 1e-8;
    z *= (1-1e-8);
    Tx xt = 2 * (x / (1-z)) - 1;
    Tx yt = 2 * (y / (1-z)) - 1;

    VectorMem<10, Tx> polsx(order+1);
    VectorMem<10, Tx> polsy(order+1);
    
    ArrayMem<Tx, 20> memz(sqr(order+1));
    FlatMatrix<Tx> polsz(order+1, &memz[0]);

    Tx fac(1.0);
    for (int i = 0; i <= order; i++)
      {
	// JacobiPolynomial (order, 2*z-1, 2*i+2, 0, polsz.Row(i));
        JacobiPolynomialAlpha (order, 2*z-1, 2*i+2, polsz.Row(i));
	polsz.Row(i) *= fac;
	fac *= (1-z);
      }

    LegendrePolynomial (order, xt, polsx);
    LegendrePolynomial (order, yt, polsy);
    
    int ii = 0;
    for (int iz = 0; iz <= order; iz++)
      for (int ix = 0; ix <= order-iz; ix++)
	for (int iy = 0; iy <= order-iz; iy++, ii++)
	  shape[ii] = polsx(ix) * polsy(iy) * polsz(max2(ix,iy), iz);
  }


  /* *********************** Hex  **********************/


  template<> template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_HEX> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x, y = ip.y, z = ip.z; 
    
    // no orientation necessary
    int p=order_inner[0];
    int q=order_inner[1];
    int r=order_inner[2];
    
    STACK_ARRAY(Tx, mem, p+q+r+3);
    Tx * polx = &mem[0];
    Tx * poly = &mem[p+1];
    Tx * polz = &mem[p+q+2];

    LegendrePolynomial (p, 2*x-1, polx);
    LegendrePolynomial (q, 2*y-1, poly);
    LegendrePolynomial (r, 2*z-1, polz);
  
    for (int i = 0, ii = 0; i <= p; i++)
      for (int j = 0; j <= q; j++)
        {
          Tx hval = polx[i] * poly[j];
          for (int k = 0; k <= r; k++)
            shape[ii++] = hval * polz[k];
        }

  }







}

#endif
