/*********************************************************************/
/* File:   l2hofe.cpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/



#include <fem.hpp>


namespace ngfem
{
  using namespace ngfem;

  template <int D>
  L2HighOrderFiniteElement<D> ::   
  L2HighOrderFiniteElement (int dim, ELEMENT_TYPE aeltype)
    : NodalFiniteElement<D> (dim, aeltype, -1, -1) 
  { 
    for (int i = 0; i < 8; i++)
      vnums[i] = i;
    
  }
  
  template <int D>
  void L2HighOrderFiniteElement<D> :: 
  SetVertexNumbers (FlatArray<int> & avnums)
  {
    for (int i=0; i<avnums.Size(); i++)
      vnums[i] = avnums[i];
    for (int i=avnums.Size(); i<8; i++)
      vnums[i] = -1;
  }

  template <int D>
  void L2HighOrderFiniteElement<D> :: 
  SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    for (int i=0; i<avnums.Size(); i++)
      vnums[i] = avnums[i];
    for (int i=avnums.Size(); i<8; i++)
      vnums[i] = -1;
  }


  template <int D>
  void L2HighOrderFiniteElement<D>::
  SetOrder (int o)
  {
    order_inner = INT<3>(o,o,o); 
    ComputeNDof();
  }

  template <int D>
  void L2HighOrderFiniteElement<D>::
  SetOrderInner (INT<3> oi)
  {
    order_inner = oi;   
    ComputeNDof();
  }

  template <int D>
  void L2HighOrderFiniteElement<D>:: 
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize(0);
    for (int i = 0; i < NodalFiniteElement<D>::GetNDof(); i++)
      idofs.Append (i);
  }


  /* *********************** Segment  **********************/

  L2HighOrderSegm :: L2HighOrderSegm (int aorder)
    : L2HighOrderFiniteElement<1>(1, ET_SEGM)
  {
    order = aorder;
    order_inner = INT<3> (aorder, aorder, aorder);
    ComputeNDof();
  }


  void L2HighOrderSegm :: ComputeNDof()
  {
    order = order_inner[0];
    ndof = order+1;
  }

  void L2HighOrderSegm :: CalcShape (const IntegrationPoint & ip, 
				     FlatVector<> shape) const
  {
    double x = ip(0);
    // orient
    if ( vnums[0] > vnums[1])
      x = 1-x;
    
    LegendrePolynomial (order, 2*x-1, shape);
  }
  
  void L2HighOrderSegm :: CalcDShape (const IntegrationPoint & ip, 
				     FlatMatrix<> dshape) const
  {
    AutoDiff<1> x(ip(0),0);
    if ( vnums[0] > vnums[1])
      x = 1-x;
    
    ArrayMem<AutoDiff<1>, 20> pol(order_inner[0]+1); 
    LegendrePolynomial (order_inner[0], 2*x-1, pol);

    for(int i=0; i<order_inner[0]+1; i++) 
      dshape(i,0) = pol[i].DValue(0); 
  }

  /* *********************** Triangle  **********************/

  L2HighOrderTrig :: L2HighOrderTrig (int aorder)
    : L2HighOrderFiniteElement<2>(2, ET_TRIG)
  {
    order_inner = INT<3> (aorder,aorder,aorder); 
    ComputeNDof();
  }


  void L2HighOrderTrig :: ComputeNDof()
  {
    ndof = (order_inner[0]+1) * (order_inner[0]+2) / 2;
    order = order_inner[0]; 
  }

  /*
  void L2HighOrderTrig :: CalcShape (const IntegrationPoint & ip, 
				     FlatVector<> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    int ii = 0;

    int n = order_inner[0];

    ArrayMem<double, 20> polx(n+1), poly(n+1);

    ScaledLegendrePolynomial (n, 2*x+y-1, 1-y, polx);
    LegendrePolynomial (n, 2*y-1, poly);
    
    for (int i = 0; i <= order_inner[0]; i++)
      for (int j = 0; j <= order_inner[0]-i; j++)
	shape(ii++) = polx[i] * poly[j];
  }
#ifdef DSHAPE_L2
  void L2HighOrderTrig :: CalcDShape (const IntegrationPoint & ip, 
				     FlatMatrix<> dshape) const
  {
    AutoDiff<2> x(ip(0),0), y(ip(1),1);
    int ii = 0;
    int n = order_inner[0];

    ArrayMem<AutoDiff<2>, 20> polx(n+1), poly(n+1);

    ScaledLegendrePolynomial (n, 2*x+y-1, 1-y, polx);
    LegendrePolynomial (n, 2*y-1, poly);
    
    for (int i = 0; i <= order_inner[0]; i++)
      for (int j = 0; j <= order_inner[0]-i; j++, ii++)
	for(int l=0;l<2;l++) 
	  dshape(ii,l) = polx[i].DValue(l) * poly[j].Value() + polx[i].Value()*poly[j].DValue(l);
	  
  }
#endif
*/



  void L2HighOrderTrig :: CalcShape (const IntegrationPoint & ip, 
				     FlatVector<> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    
    double lami[3] = { x, y, 1-x-y };
    int fav[3] = { 0, 1, 2};
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
    x = lami[fav[0]]; y=lami[fav[1]];
    
    
     
    int ii = 0;

    int n = order_inner[0];

    ArrayMem<double, 20> polx(n+1), poly(n+1);

    ScaledLegendrePolynomial (n, 2*x+y-1, 1-y, polx);
    
    for (int i = 0; i <= order_inner[0]; i++)
      {
        JacobiPolynomial (n-i, 2*y-1, 2*i+1, 0, poly);
        for (int j = 0; j <= order_inner[0]-i; j++)
          shape(ii++) = polx[i] * poly[j];
      }
  }

  void L2HighOrderTrig :: CalcDShape (const IntegrationPoint & ip, 
				     FlatMatrix<> dshape) const
  {
    AutoDiff<2> x(ip(0),0), y(ip(1),1);
      
    AutoDiff<2> lami[3] = { x, y, 1-x-y };
    int fav[3] = { 0, 1, 2};
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
    x = lami[fav[0]]; y=lami[fav[1]];
    
    int ii = 0;
    int n = order_inner[0];

    ArrayMem<AutoDiff<2>, 20> polx(n+1), poly(n+1);

    ScaledLegendrePolynomial (n, 2*x+y-1, 1-y, polx);
    
    for (int i = 0; i <= order_inner[0]; i++)
    {
      JacobiPolynomial (n, 2*y-1, 2*i+1, 0, poly);
      for (int j = 0; j <= order_inner[0]-i; j++, ii++)
        for(int l=0;l<2;l++) 
          dshape(ii,l) = polx[i].DValue(l) * poly[j].Value() + polx[i].Value()*poly[j].DValue(l);
    }
  }





  /* *********************** Quadrilateral  **********************/

  L2HighOrderQuad :: L2HighOrderQuad (int aorder)
    : L2HighOrderFiniteElement<2>(2, ET_QUAD)
  {
    order_inner = INT<3>(aorder,aorder,aorder);
    ComputeNDof();
  }

  void L2HighOrderQuad :: ComputeNDof()
  {
    ndof = (order_inner[0]+1) * (order_inner[1]+1); 
    order = max(order_inner[0],order_inner[1]); 
  }

  void L2HighOrderQuad :: CalcShape (const IntegrationPoint & ip, 
				     FlatVector<> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    
    // orient: copied from h1
    double lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    double sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    int fmax = 0; 
    for (int j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax]) fmax = j;
    int f1 = (fmax+3)%4; 
    int f2 = (fmax+1)%4; 
    if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f1 > f2; 
    x = sigma[fmax]-sigma[f1]; 
    y = sigma[fmax]-sigma[f2]; 
    
    int ii = 0;
    
    int n = max(order_inner[0],order_inner[1]);
    ArrayMem<double, 20> polx(n+1), poly(n+1);

    LegendrePolynomial (n, x, polx);
    LegendrePolynomial (n, y, poly);
    
    for (int i = 0; i <= order_inner[0]; i++)
      for (int j = 0; j <= order_inner[1]; j++)
	shape(ii++) = polx[i] * poly[j];
  }
  
  void L2HighOrderQuad :: CalcDShape(const IntegrationPoint & ip, 
				     FlatMatrix<> dshape) const  
  {
    AutoDiff<2> x(ip(0), 0);
    AutoDiff<2> y(ip(1), 1);

      // orient: copied from h1
    AutoDiff<2> lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    int fmax = 0; 
    for (int j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax]) fmax = j;
    int f1 = (fmax+3)%4; 
    int f2 = (fmax+1)%4; 
    if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f1 > f2; 
    x = sigma[fmax]-sigma[f1]; 
    y = sigma[fmax]-sigma[f2]; 
     
    
    int ii=0; 
    int n = max(order_inner[0],order_inner[1]); 
    
    ArrayMem<AutoDiff<2>,20> polx(n+1), poly(n+1); 
    LegendrePolynomial (n, x, polx); 
    LegendrePolynomial (n, y, poly); 
    
    for (int i = 0; i <= order_inner[0]; i++)
      for (int j = 0; j <= order_inner[1]; j++)
    { 
      dshape(ii,0) = polx[i].DValue(0) * poly[j].Value() + polx[i].Value() * poly[j].DValue(0);
      dshape(ii++,1) = polx[i].DValue(1) * poly[j].Value() + polx[i].Value() * poly[j].DValue(1);
    } 
        
  }
  
  
  /* *********************** Tetrahedron  **********************/

  L2HighOrderTet :: L2HighOrderTet (int aorder) 
  : L2HighOrderFiniteElement<3>(3, ET_TET)
  {
    order = aorder;
    order_inner = INT<3>(aorder, aorder, aorder);
    ComputeNDof();
  }


  void L2HighOrderTet :: ComputeNDof()
  {
    ndof = ((order+1) * (order+2) * (order+3)) / 6; // P_k
    order = order_inner[0];
  }


  void L2HighOrderTet :: CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    
    // no orientation necessary
    int n=order;
   // double eps=0.0;
    ArrayMem<double, 20> polx(n+1), poly(n+1), polz(n+1);

//     LegendrePolynomial (n, 2*x-(1+eps), polx);
//     ScaledLegendrePolynomial (n, 2*y+x-(1+eps), (1+eps)-x, poly);
//     ScaledLegendrePolynomial (n, 2*z+x+y-(1+eps), (1+eps)-x-y, polz);
    
//     int ii = 0;
//     for (int i = 0; i <= order; i++)
//       for (int j = 0; j <= order-i; j++)
//         for (int k = 0; k <= order-i-j; k++, ii++)
//           for (int l=0; l<3; l++)
//	      shape(ii++) = polx[i] * poly[j] * polz[k];



// Polynomials orthogonal w.r.t L2 inner product
    ScaledLegendrePolynomial ( n, 2*x+y+z-1, 1-y-z, polx);
  
    int ii = 0;
    for (int i = 0; i <= order; i++)
      {
	ScaledJacobiPolynomial (n, 2*y+z-1, 1-z, 2*i+1, 0, poly);
	for (int j = 0; j <= order-i; j++)
	  {
	    JacobiPolynomial(n, 2*z-1, 2*i+2*j+2, 0, polz);
	    for (int k = 0; k <= order-i-j; k++)
	      shape(ii++) = polx[i] * poly[j] * polz[k];
	  }
      }
  }

  void L2HighOrderTet :: CalcDShape (const IntegrationPoint & ip, FlatMatrix<> dshape) const
  {
    AutoDiff<3> x(ip(0),0);
    AutoDiff<3> y(ip(1),1); 
    AutoDiff<3> z(ip(2),2);
    
      
   // double eps=1e-10;
    
    int n = order;
    ArrayMem<AutoDiff<3>, 20> polx(n+1), poly(n+1), polz(n+1);

//     LegendrePolynomial (n, 2*x-(1+eps), polx);
//     ScaledLegendrePolynomial (n, 2*y+x-(1+eps), (1+eps)-x, poly);
//     ScaledLegendrePolynomial (n, 2*z+x+y-(1+eps), (1+eps)-x-y, polz);
  
//     int ii = 0;
//     for (int i = 0; i <= order; i++)
//       for (int j = 0; j <= order-i; j++)
//         for (int k = 0; k <= order-i-j; k++, ii++)
//           for (int l=0; l<3; l++)
//             dshape(ii,l) = polx[i].DValue(l) * poly[j].Value()   * polz[k].Value() 
//                 +  polx[i].Value()   * poly[j].DValue(l) * polz[k].Value() 
//                 +  polx[i].Value()   * poly[j].Value()   * polz[k].DValue(l);

// Polynomials orthogonal w.r.t L2 inner product
    ScaledLegendrePolynomial ( n, 2*x+y+z-1, 1-y-z, polx);
  
    int ii = 0;
    for (int i = 0; i <= order; i++)
      {
	ScaledJacobiPolynomial (n, 2*y+z-1, 1-z, 2*i+1, 0, poly);
	for (int j = 0; j <= order-i; j++)
	  {
	    JacobiPolynomial(n, 2*z-1, 2*i+2*j+2, 0, polz);
	    for (int k = 0; k <= order-i-j; k++,ii++) 
          { 
            AutoDiff<3> prod = polx[i] * poly[j] * polz[k]; 
           
	      for (int l=0; l<3; l++)
		dshape(ii,l) = prod.DValue(l);
          }

	  }
      }
     
  }
  


  /* *********************** Prism  **********************/

  L2HighOrderPrism :: L2HighOrderPrism (int aorder) 
  : L2HighOrderFiniteElement<3>(3, ET_PRISM)
  {
    order = aorder;
    order_inner = INT<3>(aorder, aorder, aorder);
    ComputeNDof();
  }


  void L2HighOrderPrism :: ComputeNDof()
  {
    ndof = ((order_inner[0]+1) * (order_inner[0]+2) * (order_inner[2]+1)) / 2; // P_k x Q_k
    order = max(order_inner[0], order_inner[2]);
  }


  void L2HighOrderPrism :: CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    
    // no orientation necessary
    int p=order_inner[0];
    int q=order_inner[1];
    
    ArrayMem<double, 20> polx(p+1), poly(p+1), polz(q+1);

    LegendrePolynomial (p, 2*x-1, polx);
    ScaledLegendrePolynomial (p, 2*y+x-1, 1-x, poly);
    
    LegendrePolynomial (q, 2*z-1, polz);
  
    int ii = 0;
    for (int i = 0; i <= p; i++)
      for (int j = 0; j <= p-i; j++)
        for (int k = 0; k <= q; k++)
          shape(ii++) = polx[i] * poly[j] * polz[k];
    
  }

  void L2HighOrderPrism :: CalcDShape (const IntegrationPoint & ip, FlatMatrix<> dshape) const
  {
    AutoDiff<3> x(ip(0),0), y(ip(1),1), z(ip(2),2);
    
      
    int p = order_inner[0];
    int q = order_inner[2];
      
    ArrayMem<AutoDiff<3>, 20> polx(p+1), poly(p+1), polz(q+1);

    LegendrePolynomial (p, 2*x-1, polx);
    ScaledLegendrePolynomial (p, 2*y+x-1, 1-x, poly);
    LegendrePolynomial (q, 2*z-1, polz);
      
    int ii = 0;
    for (int i = 0; i <= p; i++)
      for (int j = 0; j <= p-i; j++)
        for (int k = 0; k <= q; k++, ii++)
          for (int l=0; l<3; l++)
            dshape(ii,l) = polx[i].DValue(l) * poly[j].Value()   * polz[k].Value() 
                +  polx[i].Value()   * poly[j].DValue(l) * polz[k].Value() 
                +  polx[i].Value()   * poly[j].Value()   * polz[k].DValue(l);
  }
  


  /* *********************** Hex  **********************/

  L2HighOrderHex :: L2HighOrderHex (int aorder) 
  : L2HighOrderFiniteElement<3>(3, ET_HEX)
  {
    order = aorder;
    order_inner = INT<3>(aorder, aorder, aorder);
    ComputeNDof();
  }


  void L2HighOrderHex :: ComputeNDof()
  {
    ndof = (order_inner[0]+1) * (order_inner[1]+1) * (order_inner[2]+1); 
    order = max(order_inner[0], max(order_inner[1], order_inner[2]) );
  }


  void L2HighOrderHex :: CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    
    // no orientation necessary
    int p=order_inner[0];
    int q=order_inner[1];
    int r=order_inner[2];
    
    ArrayMem<double, 20> polx(p+1), poly(q+1), polz(r+1);

    LegendrePolynomial (p, 2*x-1, polx);
    LegendrePolynomial (q, 2*y-1, poly);
    LegendrePolynomial (r, 2*z-1, polz);
  
    int ii = 0;
    for (int i = 0; i <= p; i++)
      for (int j = 0; j <= q; j++)
        for (int k = 0; k <= r; k++)
          shape(ii++) = polx[i] * poly[j] * polz[k];
  }

  void L2HighOrderHex :: CalcDShape (const IntegrationPoint & ip, FlatMatrix<> dshape) const
  {
    AutoDiff<3> x(ip(0),0), y(ip(1),1), z(ip(2),2);
    
      
    int p = order_inner[0];
    int q = order_inner[1];
    int r = order_inner[2];
      
    ArrayMem<AutoDiff<3>, 20> polx(p+1), poly(q+1), polz(r+1);

    LegendrePolynomial (p, 2*x-1, polx);
    LegendrePolynomial (q, 2*y-1, poly);
    LegendrePolynomial (r, 2*z-1, polz);
      
    int ii = 0;
    for (int i = 0; i <= p; i++)
      for (int j = 0; j <= q; j++)
        for (int k = 0; k <= r; k++, ii++)
          for (int l=0; l<3; l++)
            dshape(ii,l) = polx[i].DValue(l) * poly[j].Value() * polz[k].Value() 
                +  polx[i].Value()   * poly[j].DValue(l) * polz[k].Value() 
                +  polx[i].Value()   * poly[j].Value()   * polz[k].DValue(l);
  }







  template class L2HighOrderFiniteElement<1>;
  template class L2HighOrderFiniteElement<2>;
  template class L2HighOrderFiniteElement<3>;

} // namespace








