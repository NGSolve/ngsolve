/*********************************************************************/
/* File:   l2hofe.cpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/



#include <fem.hpp>


namespace ngfem
{
  using namespace ngfem;

  /*
  template <int D>
  L2HighOrderFiniteElement<D> ::   
  L2HighOrderFiniteElement ()  // (ELEMENT_TYPE aeltype)
    : ScalarFiniteElement<D> (ET_TRIG, -1, -1) 
  { 
    for (int i = 0; i < 8; i++)
      vnums[i] = i;
  }
  */

  template <int D>
  void L2HighOrderFiniteElement<D> :: 
  SetVertexNumbers (FlatArray<int> & avnums)
  {
    for (int i=0; i<avnums.Size(); i++)
      vnums[i] = avnums[i];
    for (int i=avnums.Size(); i<8; i++)
      vnums[i] = -1;
  }

  /*
  template <int D>
  void L2HighOrderFiniteElement<D> :: 
  SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    for (int i=0; i<avnums.Size(); i++)
      vnums[i] = avnums[i];
    for (int i=avnums.Size(); i<8; i++)
      vnums[i] = -1;
  }
  */

  template <int D>
  void L2HighOrderFiniteElement<D>::
  SetOrder (int o)
  {
    for (int j = 0; j < DIM; j++)
      order_inner[j] = o;
    ComputeNDof();
  }

  template <int D>
  void L2HighOrderFiniteElement<D>::
  SetOrder (INT<DIM> oi)
  {
    order_inner = oi;   
    ComputeNDof();
  }

  template <int D>
  void L2HighOrderFiniteElement<D>:: 
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize(0);
    for (int i = 0; i < ScalarFiniteElement<D>::GetNDof(); i++)
      idofs.Append (i);
  }







  template <ELEMENT_TYPE ET>
  void T_L2HighOrderFiniteElement<ET> :: 
  CalcShape (const IntegrationPoint & ip, 
             FlatVector<> shape) const
  {
    double pt[DIM];
    for (int i = 0; i < DIM; i++) pt[i] = ip(i);
    static_cast<const L2HighOrderFE<ET>*> (this) -> T_CalcShape (pt, shape); 
  }

  template <ELEMENT_TYPE ET>
  void T_L2HighOrderFiniteElement<ET> :: 
  CalcDShape (const IntegrationPoint & ip, 
              FlatMatrixFixWidth<DIM> dshape) const
  {
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    ArrayMem<AutoDiff<DIM>,40> sds(ndof);
    static_cast<const L2HighOrderFE<ET>*> (this) -> T_CalcShape (adp, sds);

    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < DIM; j++)
	dshape(i,j) = sds[i].DValue(j);
  }








  template <ELEMENT_TYPE ET>
  void T_L2HighOrderFiniteElement<ET> :: 
  ComputeNDof()
  {
    switch (ET)
      {
      case ET_SEGM:
	ndof = order+1;
	break;
      case ET_TRIG:
	ndof = (order_inner[0]+1) * (order_inner[0]+2) / 2;
	break;
      case ET_QUAD:
	ndof = (order_inner[0]+1) * (order_inner[1]+1); 
	break;
      case ET_TET: 
	ndof = ((order+1) * (order+2) * (order+3)) / 6; // P_k
        break;
      case ET_PRISM:
	ndof = ((order_inner[0]+1) * (order_inner[0]+2) * (order_inner[2]+1)) / 2; // P_k x Q_k
        break;
      case ET_PYRAMID:
	cout << "pyramid not implemented" << endl;
        break;
      case ET_HEX:
	ndof = (order_inner[0]+1) * (order_inner[1]+1) * (order_inner[2]+1); 
        break;
      default:
        ;
      }
    
    order = 0;
    for (int i = 0; i < DIM; i++)
      order = max(order, order_inner[i]);
  }




  /* *********************** Segment  **********************/


  L2HighOrderFE<ET_SEGM> :: L2HighOrderFE (int aorder)
  {
    order_inner = INT<1> (aorder);
    ComputeNDof();
  }


  template<typename Tx, typename TFA>  
  void L2HighOrderFE<ET_SEGM> :: T_CalcShape (Tx hx[1], TFA & shape) const
  {
    Tx x = hx[0];
    // orient
    if ( vnums[0] > vnums[1])
      x = 1-x;
    
    LegendrePolynomial (order, 2*x-1, shape);
  }
  

  /* *********************** Triangle  **********************/

  L2HighOrderFE<ET_TRIG> :: L2HighOrderFE (int aorder)
  {
    order_inner = INT<2> (aorder,aorder); 
    ComputeNDof();
  }


  template<typename Tx, typename TFA>  
  void L2HighOrderFE<ET_TRIG> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx x = hx[0];
    Tx y = hx[1];
    
    Tx lami[3] = { x, y, 1-x-y };
    int fav[3] = { 0, 1, 2};
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
    x = lami[fav[0]]; y=lami[fav[1]];
    
     
    int ii = 0;
    int n = order_inner[0];

    ArrayMem<Tx, 20> polx(n+1), poly(n+1);

    ScaledLegendrePolynomial (n, 2*x+y-1, 1-y, polx);
    
    for (int i = 0; i <= order_inner[0]; i++)
      {
        JacobiPolynomial (n-i, 2*y-1, 2*i+1, 0, poly);
        for (int j = 0; j <= order_inner[0]-i; j++)
          shape[ii++] = polx[i] * poly[j];
      }
  }


  /* *********************** Quadrilateral  **********************/

  L2HighOrderFE<ET_QUAD> :: L2HighOrderFE (int aorder)
  {
    order_inner = INT<2>(aorder,aorder);
    ComputeNDof();
  }


  template<typename Tx, typename TFA>  
  void L2HighOrderFE<ET_QUAD> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx x = hx[0];
    Tx y = hx[1];
    
    // orient: copied from h1
    Tx lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
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
    ArrayMem<Tx, 20> polx(n+1), poly(n+1);

    LegendrePolynomial (n, x, polx);
    LegendrePolynomial (n, y, poly);
    
    for (int i = 0; i <= order_inner[0]; i++)
      for (int j = 0; j <= order_inner[1]; j++)
	shape[ii++] = polx[i] * poly[j];
  }
  
  /* *********************** Tetrahedron  **********************/

  L2HighOrderFE<ET_TET> :: L2HighOrderFE (int aorder)
  {
    order_inner = INT<3>(aorder, aorder, aorder);
    ComputeNDof();
  }


  template<typename Tx, typename TFA>  
  void L2HighOrderFE<ET_TET> :: T_CalcShape (Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0];
    Tx y = hx[1];
    Tx z = hx[2];
    
    // no orientation necessary
    int n=order;
   // double eps=0.0;
    ArrayMem<Tx, 20> polx(n+1), poly(n+1), polz(n+1);

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
	      shape[ii++] = polx[i] * poly[j] * polz[k];
	  }
      }
  }


  /* *********************** Prism  **********************/

  L2HighOrderFE<ET_PRISM> :: L2HighOrderFE (int aorder)
  {
    order_inner = INT<3>(aorder, aorder, aorder);
    ComputeNDof();
  }



  template<typename Tx, typename TFA>  
  void L2HighOrderFE<ET_PRISM> :: T_CalcShape (Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0];
    Tx y = hx[1];
    Tx z = hx[2];
    
    // no orientation necessary
    int p=order_inner[0];
    int q=order_inner[1];
    
    ArrayMem<Tx, 20> polx(p+1), poly(p+1), polz(q+1);

    LegendrePolynomial (p, 2*x-1, polx);
    ScaledLegendrePolynomial (p, 2*y+x-1, 1-x, poly);
    
    LegendrePolynomial (q, 2*z-1, polz);
  
    int ii = 0;
    for (int i = 0; i <= p; i++)
      for (int j = 0; j <= p-i; j++)
        for (int k = 0; k <= q; k++)
          shape[ii++] = polx[i] * poly[j] * polz[k];
  }


  /* *********************** Hex  **********************/

  L2HighOrderFE<ET_HEX> :: L2HighOrderFE (int aorder)
  {
    order_inner = INT<3>(aorder, aorder, aorder);
    ComputeNDof();
  }



  template<typename Tx, typename TFA>  
  void L2HighOrderFE<ET_HEX> :: T_CalcShape (Tx hx[3], TFA & shape) const
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
  
    int ii = 0;
    for (int i = 0; i <= p; i++)
      for (int j = 0; j <= q; j++)
        for (int k = 0; k <= r; k++)
          shape[ii++] = polx[i] * poly[j] * polz[k];
  }

  template class L2HighOrderFiniteElement<1>;
  template class L2HighOrderFiniteElement<2>;
  template class L2HighOrderFiniteElement<3>;

} // namespace








