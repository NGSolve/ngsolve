/*********************************************************************/
/* File:   l2hofe.cpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#include <fem.hpp>
#include "l2hofe.hpp"


#include "tscalarfe.cpp"


namespace ngfem
{

  template <int D>
  void L2HighOrderFiniteElement<D>:: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
    IntegrationRule ir(eltype, 2*order);
    VectorMem<50> shape(ndof);
    mass = 0;
    for (int i = 0; i < ir.Size(); i++)
      {
        this -> CalcShape (ir[i], shape);
        for (int j = 0; j < ndof; j++)
          mass(j) += ir[i].Weight() * sqr (shape(j));
      }
  }


  template <int D>
  void L2HighOrderFiniteElement<D>:: 
  CalcTraceMatrix (int facet, FlatMatrix<> trace) const
  {
    ELEMENT_TYPE ftype = ElementTopology::GetFacetType (eltype, facet);
    Facet2ElementTrafo f2el(eltype, FlatArray<int> (8, const_cast<int*> (vnums)) );
    const IntegrationRule & ir = SelectIntegrationRule (ftype, 2*order);

    L2HighOrderFiniteElement<1> * facetfe1 = NULL;
    L2HighOrderFiniteElement<2> * facetfe2 = NULL;
    switch (ftype)
      {
      case ET_SEGM : facetfe1 = new L2HighOrderFE<ET_SEGM> (order); break;
      case ET_TRIG : facetfe2 = new L2HighOrderFE<ET_TRIG> (order); break;
      case ET_QUAD : facetfe2 = new L2HighOrderFE<ET_QUAD> (order); break;
      default:
	;
      }

    int ndof_facet = trace.Height();
    Vector<> shape(ndof);
    Vector<> fshape(ndof_facet);
    Vector<> norms(ndof_facet);

    trace = 0.0;
    norms = 0.0;
    for (int i = 0; i < ir.Size(); i++)
      {
	if (D == 2)
	  facetfe1 -> CalcShape (ir[i], fshape);
	else
	  facetfe2 -> CalcShape (ir[i], fshape);

	this -> CalcShape (f2el (facet, ir[i]), shape);

	trace += ir[i].Weight() * fshape * Trans (shape);
	for (int j = 0; j < norms.Size(); j++)
	  norms(j) += ir[i].Weight() * sqr (fshape(j));
      }

    for (int j = 0; j < fshape.Size(); j++)
      trace.Row(j) /= norms(j);

    delete facetfe1;
    delete facetfe2;
  }


  template <int D>
  void L2HighOrderFiniteElement<D>:: 
  CalcGradientMatrix (FlatMatrix<> gmat) const
  {
    IntegrationRule ir (eltype, 2*order);

    Vector<> shape(ndof);
    MatrixFixWidth<D> dshape(ndof);
    Vector<> norms(ndof);
    
    gmat = 0.0;
    norms = 0.0;
    for (int i = 0; i < ir.Size(); i++)
      {
	this -> CalcShape (ir[i], shape);
	this -> CalcDShape (ir[i], dshape);
        
        for (int j = 0; j < ndof; j++)
          for (int k = 0; k < ndof; k++)
            for (int l = 0; l < D; l++)
              gmat(k*D+l, j) += ir[i].Weight() * dshape(j,l) * shape(k);

	for (int j = 0; j < norms.Size(); j++)
	  norms(j) += ir[i].Weight() * sqr (shape(j));
      }
    for (int j = 0; j < ndof; j++)
      gmat.Rows(D*j, D*(j+1)) /= norms(j);
  }


  template <ELEMENT_TYPE ET>
  void T_L2HighOrderFiniteElement<ET> :: 
  ComputeNDof()
  {
    switch (ET)
      {
      case ET_SEGM:
	ndof = order_inner[0]+1;
	break;
      case ET_TRIG:
	ndof = (order_inner[0]+1) * (order_inner[0]+2) / 2;
	break;
      case ET_QUAD:
	ndof = (order_inner[0]+1) * (order_inner[1]+1); 
	break;
      case ET_TET: 
	ndof = ((order_inner[0]+1) * (order_inner[0]+2) * (order_inner[0]+3)) / 6; // P_k
        break;
      case ET_PRISM:
	ndof = ((order_inner[0]+1) * (order_inner[0]+2) * (order_inner[2]+1)) / 2; // P_k x Q_k
        break;
      case ET_PYRAMID:
        ndof = (order_inner[0]+2)*(order_inner[0]+1)*(2*order_inner[0]+3) / 6; 
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



  


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  typename L2HighOrderFE<ET,SHAPES>::TPRECOMP L2HighOrderFE<ET,SHAPES>::precomp;

  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  typename L2HighOrderFE<ET,SHAPES>::TPRECOMP_TRACE L2HighOrderFE<ET,SHAPES>::precomp_trace(320);

  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  typename L2HighOrderFE<ET,SHAPES>::TPRECOMP_GRAD L2HighOrderFE<ET,SHAPES>::precomp_grad(40);




  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  void L2HighOrderFE<ET,SHAPES> :: 
  PrecomputeTrace ()
  {
      for (int f = 0; f < ElementTopology::GetNFacets(ET); f++)
	{
	  int classnr =  ET_trait<ET>::GetFacetClassNr (f, vnums);

	  if (precomp_trace.Used (INT<2> (order, classnr)))
	    continue;

	  ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (ET, f);
	  int nf;
	  switch (etfacet)
	    {
	    case ET_SEGM: nf = order+1; break;
	    case ET_TRIG: nf = (order+1)*(order+2)/2; break;
	    case ET_QUAD: nf = sqr(order+1); break;
	    default: nf = 0;
	    }
	  

	  Matrix<> * trace = new Matrix<>(nf, ndof);
	  L2HighOrderFiniteElement<DIM>::CalcTraceMatrix (f, *trace);
	  precomp_trace.Set (INT<2> (order, classnr), trace);
	}
    }
    

  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  void L2HighOrderFE<ET,SHAPES> :: 
  PrecomputeGrad ()
  {
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    
    if (precomp_grad.Used (INT<2> (order, classnr)))
      return;

    Matrix<> * gmat = new Matrix<>(ndof*DIM, ndof);
    L2HighOrderFiniteElement<DIM>::CalcGradientMatrix (*gmat);
    precomp_grad.Set (INT<2> (order, classnr), gmat);
  }
    


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  void L2HighOrderFE<ET,SHAPES> :: 
  PrecomputeShapes (const IntegrationRule & ir) 
  {
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
  }
  


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  void L2HighOrderFE<ET,SHAPES> :: 
  Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const
  {
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    
    PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());
    if (pre)
      vals = pre->shapes * coefs;
    else
      T_ScalarFiniteElement2< SHAPES<ET>, ET > :: Evaluate (ir, coefs, vals);
  }



  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  void L2HighOrderFE<ET,SHAPES> :: 
  EvaluateGrad (const IntegrationRule & ir, FlatVector<> coefs, FlatMatrixFixWidth<DIM> values) const
    {
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
	{
	  FlatVector<> vval(DIM*values.Height(), &values(0,0));
	  vval = pre->dshapes * coefs;
	}
      else
	T_ScalarFiniteElement2< SHAPES<ET>, ET > :: EvaluateGrad (ir, coefs, values);
    }


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  void L2HighOrderFE<ET,SHAPES> :: 
  EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> values, FlatVector<> coefs) const
    {
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
	coefs = Trans (pre->dshapes) * FlatVector<> (DIM*ndof, &values(0,0));  // values.Height !!!
      else
	T_ScalarFiniteElement2< SHAPES<ET>, ET > :: EvaluateGradTrans (ir, values, coefs);
    }






  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  void L2HighOrderFE<ET,SHAPES> :: 
  GetGradient (FlatVector<> coefs, FlatMatrixFixWidth<DIM> grad) const
  {
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    int bnr, pos;
    if (precomp_grad.Used (INT<2> (order, classnr), bnr, pos))
      {
        FlatMatrix<> gmat = *precomp_grad.Get (bnr, pos);
        FlatVector<> vgrad(grad.Height()*DIM, &grad(0,0));
        vgrad = gmat * coefs;
      }
    else
      L2HighOrderFiniteElement<DIM>::GetGradient (coefs, grad);
  }
  
  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  void L2HighOrderFE<ET,SHAPES> :: 
  GetGradientTrans (FlatMatrixFixWidth<DIM> grad, FlatVector<> coefs) const
  {
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    int bnr, pos;
    if (precomp_grad.Used (INT<2> (order, classnr), bnr, pos))
      {
        FlatMatrix<> gmat = *precomp_grad.Get (bnr, pos);
        FlatVector<> vgrad(grad.Height()*DIM, &grad(0,0));
        coefs = Trans(gmat) * vgrad;
      }
    else
      L2HighOrderFiniteElement<DIM>::GetGradientTrans (grad, coefs);
  }









  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  void L2HighOrderFE<ET,SHAPES> :: 
  GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const
    {
      int classnr =  ET_trait<ET>::GetFacetClassNr (facet, vnums);
      int bnr, pos;
      if (precomp_trace.Used (INT<2> (order, classnr), bnr, pos))
	{
	  FlatMatrix<> trace = *precomp_trace.Get (bnr, pos);
	  fcoefs = trace * coefs;
	}
      else
	L2HighOrderFiniteElement<DIM>::GetTrace (facet, coefs, fcoefs);
    }

  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  void L2HighOrderFE<ET,SHAPES> :: 
  GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const
  {
      int classnr =  ET_trait<ET>::GetFacetClassNr (facet, vnums);
      if (precomp_trace.Used (INT<2> (order, classnr)))
	{
	  coefs = Trans(*precomp_trace.Get (INT<2> (order, classnr))) * fcoefs;
	}
      else
	L2HighOrderFiniteElement<DIM>::GetTraceTrans (facet, fcoefs, coefs);
    }


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES>
  void L2HighOrderFE<ET,SHAPES> :: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
    L2HighOrderFiniteElement<ET_trait<ET>::DIM>::GetDiagMassMatrix (mass);
  }

  template <> void L2HighOrderFE<ET_SEGM> :: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
    // L2HighOrderFiniteElement<1>::GetDiagMassMatrix (mass);

    for (int ix = 0; ix <= order; ix++)
      mass(ix) = 1.0 / (2*ix+1);
  }

  template <> void L2HighOrderFE<ET_TRIG> :: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
    // L2HighOrderFiniteElement<2>::GetDiagMassMatrix (mass);

    for (int ix = 0, ii = 0; ix <= order; ix++)
      for (int iy = 0; iy <= order-ix; iy++, ii++)
        mass(ii) = 1.0 / ( (2*iy+1) * (2*ix+2*iy+2));
  }










  /* *********************** Segment  **********************/

  /*
  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_SEGM> :: T_CalcShape (Tx x[], TFA & shape) const
  {
    Tx lam[2] = { x[0], 1-x[0] };

    INT<2> e = GetEdgeSort (0, vnums);

    LegendrePolynomial (order, lam[e[1]]-lam[e[0]], shape);
  }
  */


  /* *********************** Triangle  **********************/

  /*
  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_TRIG> :: T_CalcShape (Tx x[], TFA & shape) const
  {
    Tx lam[3] = { x[0], x[1], 1-x[0]-x[1] };
    
    INT<4> f = GetFaceSort (0, vnums);
    
    int p = order_inner[0];

    DubinerBasis::Eval (p, lam[f[0]], lam[f[1]], shape);
  }
  */

  /*
  */

  /* *********************** Quadrilateral  **********************/

  /*
  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_QUAD> :: T_CalcShape (Tx hx[], TFA & shape) const
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
b  }
  */


  /* *********************** Tetrahedron  **********************/

  /*
  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_TET> :: T_CalcShape (Tx x[], TFA & shape) const
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
  */

  /* *********************** Prism  **********************/


  /*
  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_PRISM> :: T_CalcShape (Tx hx[], TFA & shape) const
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
  */

  /* *********************** Pyramid  **********************/


  template<typename Tx, typename TFA>  
  void L2HighOrderFE_Shape<ET_PYRAMID> :: T_CalcShape (Tx hx[], TFA & shape) const
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



  template class L2HighOrderFiniteElement<1>;
  template class L2HighOrderFiniteElement<2>;
  template class L2HighOrderFiniteElement<3>;

  template NGS_DLL_HEADER class L2HighOrderFE<ET_SEGM>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_TRIG>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_QUAD>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_TET>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_PRISM>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_PYRAMID>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_HEX>;
 

  template class T_L2HighOrderFiniteElement<ET_SEGM>;
  template class T_L2HighOrderFiniteElement<ET_TRIG>;
  template class T_L2HighOrderFiniteElement<ET_QUAD>; 
  template class T_L2HighOrderFiniteElement<ET_TET>;
  template class T_L2HighOrderFiniteElement<ET_PRISM>;
  template class T_L2HighOrderFiniteElement<ET_PYRAMID>;
  template class T_L2HighOrderFiniteElement<ET_HEX>;



  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_SEGM>, ET_SEGM>;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_TRIG>, ET_TRIG>;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_QUAD>, ET_QUAD>;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_TET>, ET_TET>;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_PRISM>, ET_PRISM>;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_HEX>, ET_HEX>;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_PYRAMID>, ET_PYRAMID>;







  class testing
  { 
  public:
    testing()
    {
      cout << "test calcgrad" << endl;
      LocalHeap lh(1000000, "heaptest");
      L2HighOrderFE<ET_TRIG> fel(2);

      Vector<> mass(fel.GetNDof());
      fel.GetDiagMassMatrix (mass);
      cout << "mass = " << endl << mass << endl;

      IntegrationRule ir(ET_TRIG, 3);
      Vector<> coefs(fel.GetNDof());
      Matrix<> cgrad(fel.GetNDof(), 2);

      coefs = 0;
      coefs(2) = 1;
      
      fel.GetGradient (coefs, cgrad);
      cout << "grad = " << endl << cgrad << endl;
      
      Matrix<> vgrad(ir.GetNIP(), 2);
      fel.EvaluateGrad (ir, coefs, vgrad);
      
      cout << "grad(ir) = " << endl << vgrad << endl;

      FlatVector<> hv(ir.GetNIP(), lh);
      for (int j = 0; j < 2; j++)
        {
          fel.Evaluate (ir, cgrad.Col(j)|lh, hv);
          vgrad.Col(j) = hv;
        }
      
      cout << "grad(ir2) = " << endl << vgrad << endl;
    }
  };

  // static testing test;
  
} // namespace







