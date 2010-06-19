/*********************************************************************/
/* File:   h1hofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

 
#include <fem.hpp>
#include <h1hofe.hpp>

#include "tscalarfe.cpp"

namespace ngfem
{
  using namespace ngfem;

  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetOrderFace (FlatArray<int> & of)
  {
    for (int i = 0; i < of.Size(); i++)
      order_face[i] = INT<2>(of[i],of[i]);
  }

  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetOrderFace (FlatArray<INT<2> > & of)
  {
    for (int i = 0; i < of.Size(); i++)
      order_face[i] = of[i];
  }
  
  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetOrderEdge (FlatArray<int> & oe)
  {
    for (int i = 0; i < oe.Size(); i++)
      order_edge[i] = oe[i];
  }

  template <ELEMENT_TYPE ET>
  void T_H1HighOrderFiniteElement<ET> :: 
  GetDofs (Array<Dof> & dofs) const
  {
    dofs.SetSize (0);
    
    for (int i = 0; i < N_VERTEX; i++)
      dofs.Append (Dof (Node (NT_VERTEX, i), 0));

    for (int i = 0; i < N_EDGE; i++)
      for (int j = 0; j < order_edge[i]-1; j++)
	dofs.Append (Dof (Node (NT_EDGE, i), j));

    for (int i = 0; i < N_FACE; i++)
      {
        INT<2> p = order_face[i];

        int nf;
        if (FaceType (i) == ET_TRIG)
          nf = (p[0]-1)*(p[0]-2) / 2;
        else
          nf = (p[0]-1)*(p[1]-1);

	for (int j = 0; j < nf; j++)
	  dofs.Append (Dof (Node(NT_FACE, i), j));
      }

    int ni = 0;
    INT<3> p = order_cell;
    switch (ET)
      {
      case ET_TET:  if(p[0] > 3) ni = (p[0]-1)*(p[0]-2)*(p[0]-3)/6;   
	break;
      case ET_PRISM: if(p[0]>2 && p[2]>1) ni = (p[0]-1)*(p[0]-2)*(p[2]-1)/2;
	break;
      case ET_PYRAMID: if(p[0]>2) ni = (p[0]-1)*(p[0]-2)*(2*p[0]-3)/6;
	break;
        case ET_HEX: if(p[0]>1 && p[1] > 1 && p[2]>1) ni = (p[0]-1)*(p[1]-1)*(p[2]-1);
          break;
      default:
	;
      }
    
    for (int j = 0; j < ni; j++)
      dofs.Append (Dof (Node(NT_CELL, 0), j));
  }


  template <ELEMENT_TYPE ET>
  void T_H1HighOrderFiniteElement<ET> :: 
  ComputeNDof()
  {
    ndof = N_VERTEX;
    
    for (int i = 0; i < N_EDGE; i++)
      ndof += order_edge[i] - 1;
    
    for (int i = 0; i < N_FACE; i++)
      if (FaceType(i) == ET_TRIG)
        ndof += (order_face[i][0]-1)*(order_face[i][0]-2) / 2;
      else
        ndof += (order_face[i][0]-1)*(order_face[i][1]-1);
    
    switch (ET)
      {
      case ET_TET: 
        ndof += (order_cell[0]-1)*(order_cell[0]-2)*(order_cell[0]-3) / 6; 
        break;
      case ET_PRISM:
        ndof += (order_cell[0]-1)*(order_cell[0]-2)/2*(order_cell[2]-1);
        break;
      case ET_PYRAMID:
        ndof += (order_cell[0]-1)*(order_cell[0]-2)*(2*order_cell[0]-3) / 6; 
        break;
      case ET_HEX:
        ndof += (order_cell[0]-1)*(order_cell[1]-1)*(order_cell[2]-1);
        break;
      default:
        ;
      }
    
    order = 1;
    for (int i = 0; i < N_EDGE; i++)
      order = max(order, order_edge[i]);

    for (int i = 0; i < N_FACE; i++)
      order = max(order, max (order_face[i][0], order_face[i][1]));
    
    if (DIM == 3)
      {
	order = max(order, order_cell[0]);
	order = max(order, order_cell[1]);
	order = max(order, order_cell[2]);
      }
  }


  template <ELEMENT_TYPE ET>
  void T_H1HighOrderFiniteElement<ET> :: 
  GetInternalDofs (Array<int> & idofs) const
  {
    int ni = 0;
    switch (ET)
      {
      case ET_SEGM: 
        ni = order_edge[0]-1;
        break;
      case ET_TRIG: 
        ni = (order_face[0][0]-1)*(order_face[0][0]-2) / 2; 
        break;
      case ET_QUAD: 
        ni = (order_face[0][0]-1)*(order_face[0][1]-1);
        break;
      case ET_TET: 
        ni = (order_cell[0]-1)*(order_cell[0]-2)*(order_cell[0]-3) / 6;
        break;
      case ET_PRISM: 
        ni = (order_cell[0]-1)*(order_cell[0]-2)/2*(order_cell[2]-1);
        break;
      case ET_PYRAMID:
        ni = (order_cell[0]-1)*(order_cell[0]-2)*(2*order_cell[0]-3) / 6; 
        break;
      case ET_HEX: 
        ni = (order_cell[0]-1)*(order_cell[1]-1)*(order_cell[2]-1);
        break;
      }
    
    idofs.SetSize (0);
    for (int i = 0; i < ni; i++)
      idofs.Append (ndof-ni+i);
  }






  /* *********************** Segment  **********************/

  template<typename Tx, typename TFA>  
  void H1HighOrderFE<ET_SEGM> :: T_CalcShape (Tx hx[1], TFA & shape) const
  {
    Tx x = hx[0];
    Tx lam[2] = { x, 1-x };

    shape[0] = lam[0];
    shape[1] = lam[1];

    INT<2> e = GetEdgeSort (0, vnums);

    LegendrePolynomial leg;
    leg.EvalMult (order_edge[0]-2, 
		  lam[e[1]]-lam[e[0]], lam[e[0]]*lam[e[1]], shape.Addr(2));
  }

  /* *********************** Triangle  **********************/

  template<typename Tx, typename TFA>  
  void H1HighOrderFE<ET_TRIG> :: T_CalcShape (Tx x[2], TFA & shape) const
  {
    Tx lam[3] = { x[0], x[1], 1-x[0]-x[1] };

    for (int i = 0; i < 3; i++)
      shape[i] = lam[i];

    int ii = 3;
    
    // edge dofs
    for (int i = 0; i < N_EDGE; i++)
      if (order_edge[i] >= 2)
	{ 
          INT<2> e = GetEdgeSort (i, vnums);
	  LegendrePolynomial leg;
	  leg.EvalScaledMult (order_edge[i]-2, 
			      lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
			      lam[e[0]]*lam[e[1]], shape.Addr(ii));
	  ii += order_edge[i]-1;
	}

    int p = order_face[0][0];
    if (p >= 3)
      {
        INT<4> f = GetFaceSort (0, vnums);

	DubinerBasis dub;
	dub.EvalMult (p-3, lam[f[0]], lam[f[1]], lam[f[0]]*lam[f[1]]*lam[f[2]], shape.Addr(ii));
      }
  }


  /* *********************** Quadrilateral  **********************/

  template<typename Tx, typename TFA>  
  void H1HighOrderFE<ET_QUAD> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1];
    Tx lam[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    

    // vertex shapes
    for(int i=0; i < N_VERTEX; i++) shape[i] = lam[i]; 
    int ii = 4;

    ArrayMem<Tx,20> polxi(order+1), poleta(order+1);
     
    // edge dofs
    for (int i = 0; i < N_EDGE; i++)
      {
	int p = order_edge[i];
        INT<2> e = GetEdgeSort (i, vnums);	  
        
        Tx xi = sigma[e[1]]-sigma[e[0]]; 
        Tx lam_e = lam[e[0]]+lam[e[1]];
        Tx bub = 0.25 * lam_e * (1 - xi*xi);

	LegendrePolynomial leg;
	leg.EvalMult (p-2, xi, bub, shape.Addr(ii));
	ii += p-1;
      }    
    
    INT<2> p = order_face[0];
    if (p[0] >= 2 && p[1] >= 2)
      {
        INT<4> f = GetFaceSort (0, vnums);  // vnums[f[0]] > vnums[f[1]] > vnums[f[3]]

	Tx xi = sigma[f[0]]-sigma[f[1]]; 
	Tx eta = sigma[f[0]]-sigma[f[3]]; 
	
	T_ORTHOPOL::Calc(p[0], xi,polxi);
	T_ORTHOPOL::Calc(p[1],eta,poleta);
	
	for (int k = 0; k <= p[0]-2; k++)
	  for (int j = 0; j <= p[1]-2; j++)
	    shape[ii++] = polxi[k] * poleta[j];
      }
  }


  /* *********************** Tetrahedron  **********************/

  template<typename Tx, typename TFA>  
  void  H1HighOrderFE<ET_TET> :: T_CalcShape (Tx x[3], TFA & shape) const
  {
    Tx lam[4] = { x[0], x[1], x[2], 1-x[0]-x[1]-x[2] };

    ArrayMem<Tx, 20> polx(order+1), poly(order+1), polz(order+1); 

    // vertex shapes
    for (int i = 0; i < 4; i++)
      shape[i] = lam[i];
    int ii = 4; 

    // edge dofs
    for (int i = 0; i < N_EDGE; i++)
      if (order_edge[i] >= 2)
	{
          INT<2> e = GetEdgeSort (i, vnums);
	  LegendrePolynomial leg;
	  leg.EvalScaledMult (order_edge[i]-2, 
			      lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
			      lam[e[0]]*lam[e[1]], shape.Addr(ii));
	  ii += order_edge[i]-1;
	}

    // face dofs
    for (int i = 0; i < N_FACE; i++)
      if (order_face[i][0] >= 3)
	{
          INT<4> f = GetFaceSort (i, vnums);
	  int vop = 6 - f[0] - f[1] - f[2];  	
          
	  DubinerBasis dub;
	  int p = order_face[i][0];
	  dub.EvalScaledMult (p-3, lam[f[0]], lam[f[1]], 1-lam[vop], 
			      lam[f[0]]*lam[f[1]]*lam[f[2]], shape.Addr(ii));
	  ii += (p-2)*(p-1)/2;
	}

    if (order_cell[0] >= 4)
      ii += T_INNERSHAPES::Calc (order_cell[0], 
                                 lam[0]-lam[3], lam[1], lam[2], 
                                 shape.Addr(ii) );
  }




  /* *********************** Prism  **********************/

  template<typename Tx, typename TFA>  
  void  H1HighOrderFE<ET_PRISM> :: T_CalcShape (Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];
    Tx lam[6] = { x, y, 1-x-y, x, y, 1-x-y };
    Tx muz[6]  = { 1-z, 1-z, 1-z, z, z, z };
    
    // vertex shapes
    for (int i = 0; i < 6; i++)
      shape[i] = lam[i] * muz[i];

    int ii = 6;

    // horizontal edge dofs
    for (int i = 0; i < 6; i++)
      if (order_edge[i] >= 2)
	{
          INT<2> e = GetEdgeSort (i, vnums);

	  Tx xi = lam[e[1]]-lam[e[0]]; 
	  Tx eta = lam[e[0]]+lam[e[1]]; 

	  LegendrePolynomial leg;
	  leg.EvalScaledMult (order_edge[i]-2, xi, eta, 
			      lam[e[0]]*lam[e[1]]*muz[e[1]], shape.Addr(ii));
	  ii += order_edge[i]-1;
	}
    
    // vertical edges
    for (int i = 6; i < 9; i++)
      if (order_edge[i] >= 2)
	{
          INT<2> e = GetEdgeSort (i, vnums);
	  LegendrePolynomial leg;
	  leg.EvalMult (order_edge[i]-2, 
			muz[e[1]]-muz[e[0]], 
			muz[e[0]]*muz[e[1]]*lam[e[1]], shape.Addr(ii));
	  ii += order_edge[i]-1;
	}
    

    ArrayMem<Tx,20> polx(order+1), polz(order+1);

    const FACE * faces = ElementTopology::GetFaces (ET_PRISM);
    // trig face dofs
    for (int i = 0; i < 2; i++)
      if (order_face[i][0] >= 3)
	{
          INT<4> f = GetFaceSort (i, vnums);

	  DubinerBasis dub;
	  int p = order_face[i][0];
	  dub.EvalMult (p-3, lam[f[0]], lam[f[1]],
			lam[0]*lam[1]*lam[2]*muz[f[2]], shape.Addr(ii));
	  ii += (p-2)*(p-1)/2;
	}
   
    // quad face dofs
    for (int i = 2; i < 5; i++)
      if (order_face[i][0] >= 2 && order_face[i][1] >= 2)
	{
	  INT<2> p = order_face[i];
	  int pp = int(max2(p[0], p[1]))+1;
	 
	  int fmax = 0;
	  for(int j = 1; j < 4; j++) 
	    if(vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;
	  
	  int fz = 3-fmax; 
	  int ftrig = fmax^1; 
	  
	  Tx xi = lam[faces[i][fmax]] - lam[faces[i][ftrig]]; 
	  Tx eta = 1-lam[faces[i][fmax]]-lam[faces[i][ftrig]]; 
	  Tx zeta = muz[faces[i][fmax]]-muz[faces[i][fz]]; 
	  
	  T_ORTHOPOL::CalcTrigExt (pp, xi, eta, polx);  
	  T_ORTHOPOL::Calc (pp, zeta, polz);	
	  
	  // global x-direction is towards second-largest vertex 
	  if (vnums[faces[i][ftrig]] > vnums[faces[i][fz]])  
	    for (int k = 0; k < p[0]-1; k++)
	      for (int j = 0; j < p[1]-1; j++)
		shape[ii++] = polx[k] * polz[j];
	  else
	    for (int j = 0; j < p[0]-1; j++)
	      for (int k = 0; k < p[1]-1; k++)
		shape[ii++] = polx[k] * polz[j];
	}
    
    // volume dofs:
    if (order_cell[0] > 2 && order_cell[2] > 1)
      {
	ArrayMem<Tx,20> pol_trig((order_cell[0]-1)*(order_cell[0]-2)/2);
	int nf = T_TRIGSHAPES::Calc (order_cell[0], x-y, 1-x-y, pol_trig);

	T_ORTHOPOL:: Calc(order_cell[2], 2*z-1, polz);
	for (int i = 0; i < nf; i++)
	  for (int k = 0; k < order_cell[2]-1; k++)
	    shape[ii++] = pol_trig[i] * polz[k];
      }
  }





  /* *********************** Hex  **********************/


  template<typename Tx, typename TFA>  
  void  H1HighOrderFE<ET_HEX> :: T_CalcShape (Tx hx[3], TFA & shape) const
  { 
    Tx x = hx[0], y = hx[1], z = hx[2];

    Tx lam[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
		(1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
    Tx sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
		 (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 

    // vertex shapes
    for(int i=0; i<8; i++) shape[i] = lam[i]; 
    int ii = 8;

    ArrayMem<Tx,20> polx(order+1), poly(order+1), polz(order+1);
    
    // edge dofs
    for (int i = 0; i < N_EDGE; i++)
      if (order_edge[i] >= 2)
	{
	  int p = order_edge[i];
          INT<2> e = GetEdgeSort (i, vnums);	  

	  Tx xi = sigma[e[1]]-sigma[e[0]]; 
	  Tx lam_e = lam[e[0]]+lam[e[1]];
	  Tx bub = 0.25 * lam_e * (1 - xi*xi);
	  
	  LegendrePolynomial leg;
	  leg.EvalMult (p-2, xi, bub, shape.Addr(ii));
	  ii += p-1;

	/*
	  Tx xi = sigma[e[1]]-sigma[e[0]]; 
	  Tx lam_e = lam[e[0]]+lam[e[1]];

	  ii += T_ORTHOPOL::CalcMult (order_edge[i], 
                                      xi, lam_e, 
                                      shape.Addr(ii));
	*/
	}
     
    for (int i = 0; i < N_FACE; i++)
      if (order_face[i][0] >= 2 && order_face[i][1] >= 2)
	{
	  INT<2> p = order_face[i];
	 
          INT<4> f = ET_trait<ET_HEX>::GetFaceSort (i, vnums);	  

	  Tx lam_f = 0;
	  for (int j = 0; j < 4; j++)
	    lam_f += lam[f[j]];

	  Tx xi  = sigma[f[0]] - sigma[f[1]]; 
	  Tx eta = sigma[f[0]] - sigma[f[3]];
	  
	  T_ORTHOPOL::CalcMult (p[0], xi, lam_f, polx);	
	  T_ORTHOPOL::Calc (p[1], eta, poly);	
	  
	  for (int k = 0; k < p[0]-1; k++) 
            for (int j = 0; j < p[1]-1; j++) 
              shape[ii++]= polx[k] * poly[j];
	}
    
    // volume dofs:
    if (order_cell[0] >= 2 && order_cell[1] >= 2 && order_cell[2] >= 2)
      {
	T_ORTHOPOL::Calc (order_cell[0], 2*x-1, polx);
	T_ORTHOPOL::Calc (order_cell[1], 2*y-1, poly);
	T_ORTHOPOL::Calc (order_cell[2], 2*z-1, polz);
	for (int i = 0; i < order_cell[0]-1; i++)
	  for (int j = 0; j < order_cell[1]-1; j++)
	    {
	      Tx pxy = polx[i] * poly[j];
	      for (int k = 0; k < order_cell[2]-1; k++)
		shape[ii++] = pxy * polz[k];
	    }
      }
  }

  /* ******************************** Pyramid  ************************************ */

  template<typename Tx, typename TFA>  
  void  H1HighOrderFE<ET_PYRAMID> :: T_CalcShape (Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];

    if (z == 1.) z -= 1e-10;

    Tx xt = x / (1-z);
    Tx yt = y / (1-z);
    
    Tx sigma[4]  = { (1-xt)+(1-yt), xt+(1-yt), xt+yt, (1-xt)+yt };
    Tx lambda[4] = { (1-xt)*(1-yt), xt*(1-yt), xt*yt, (1-xt)*yt };
    Tx lambda3d[5];

    for (int i = 0; i < 4; i++)  
      lambda3d[i] = lambda[i] * (1-z);
    lambda3d[4] = z;


    for (int i = 0; i < 5; i++)  
      shape[i] = lambda3d[i];

    int ii = 5;


    //horizontal edge dofs 
    // const EDGE * edges = ElementTopology::GetEdges (ET_PYRAMID);
    for (int i = 0; i < 4; i++)
      if (order_edge[i] >= 2)
	{
	  /*
	  int es = edges[i][0], ee = edges[i][1];
	  if (vnums[es] > vnums[ee]) swap (es, ee);

	  Tx lam = sigma[ee]-sigma[es]; 
	  Tx lam_edge = lambda[es] + lambda[ee];

          ii += T_ORTHOPOL::CalcTrigExtMult (order_edge[i], lam*(1-z), z, lam_edge, 
                                             shape.Addr(ii));
	  */

	  int p = order_edge[i];
	  INT<2> e = GetEdgeSort (i, vnums);	  

	  Tx xi = sigma[e[1]]-sigma[e[0]]; 
	  Tx lam_e = lambda[e[0]]+lambda[e[1]];
	  Tx bub = 0.25 * lam_e * (1 - xi*xi)*(1-z)*(1-z);
	  
	  LegendrePolynomial leg;
	  leg.EvalScaledMult (p-2, xi*(1-z), 1-z, bub, shape.Addr(ii));
	  ii += p-1;
	}
    
    // vertical edges
    for (int i = 4; i < 8; i++) 
      if (order_edge[i] >= 2)
	{
	  /*
	  int es = edges[i][0], ee = edges[i][1];
	  if (vnums[es] > vnums[ee]) swap (es, ee);
	  
	  ii += T_ORTHOPOL::CalcTrigExt (order_edge[i], lambda3d[ee]-lambda3d[es],  
					 1-lambda3d[es]-lambda3d[ee], shape.Addr(ii));
	  */


	  int p = order_edge[i];
	  INT<2> e = GetEdgeSort (i, vnums);	  

	  Tx xi = lambda3d[e[1]]-lambda3d[e[0]]; 
	  Tx lam_e = lambda3d[e[0]]+lambda3d[e[1]];
	  Tx bub = 0.25 * (lam_e*lam_e-xi*xi);
	  
	  LegendrePolynomial leg;
	  leg.EvalScaledMult (p-2, xi, lam_e, bub, shape.Addr(ii));
	  ii += p-1;
	}


    ArrayMem<Tx,20> polx(order+1), poly(order+1), polz(order+1);
    const FACE * faces = ElementTopology::GetFaces (ET_PYRAMID);

    // trig face dofs
    for (int i = 0; i < 4; i++)
      if (order_face[i][0] >= 3)
	{
	  int p = order_face[i][0];
	  Tx lam_face = lambda[faces[i][0]] + lambda[faces[i][1]];  // vertices on quad    
	  Tx bary[3] = 
	    {(sigma[faces[i][0]]-lam_face)*(1-z), (sigma[faces[i][1]]-lam_face)*(1-z), z};  
	  int fav[3] = {0, 1, 2};
	  if(vnums[faces[i][fav[0]]] > vnums[faces[i][fav[1]]]) swap(fav[0],fav[1]); 
	  if(vnums[faces[i][fav[1]]] > vnums[faces[i][fav[2]]]) swap(fav[1],fav[2]);
	  if(vnums[faces[i][fav[0]]] > vnums[faces[i][fav[1]]]) swap(fav[0],fav[1]); 	

	  /*
	  int ndf = T_TRIGSHAPES::CalcMult
	    (p, bary[fav[2]]-bary[fav[1]], bary[fav[0]], lam_face, shape.Addr(ii));
	  
	  ii += ndf;
	  */
	  Tx bub = lam_face * bary[0]*bary[1]*bary[2];
	  
	  DubinerBasis dub;
	  dub.EvalMult (p-3, bary[fav[0]], bary[fav[1]], bub, shape.Addr(ii));
	  ii += (p-2)*(p-1)/2;


	}
    
    // quad face dof
    if (order_face[4][0] >= 2 && order_face[4][1] >= 2)
      {
	int px = order_face[4][0];
	int py = order_face[4][1];
	int p = max2(px, py);

	Tx fac = 1.0;
	for (int k = 1; k <= p; k++)
	  fac *= (1-z);

	int fmax = 0;
	for (int j=1; j<4; j++) 
	  if (vnums[j] > vnums[fmax]) fmax = j;  

	int f1 = (fmax+3)%4;
	int f2 = (fmax+1)%4; 
	if(vnums[f2]>vnums[f1]) swap(f1,f2);  // fmax > f1 > f2 
       
 	Tx xi  = sigma[fmax] - sigma[f1]; 
	Tx eta = sigma[fmax] - sigma[f2];

	T_ORTHOPOL::Calc (p, xi, polx);	
	T_ORTHOPOL::Calc (p, eta, poly);	

	for (int k = 0; k < p-1; k++) 
	  for (int j = 0; j < p-1; j++) 
	    shape[ii++] = polx[k] * poly[j] * fac; 
      }

    
    if (order_cell[0] >= 3)
      {
	T_ORTHOPOL::Calc (order_cell[0], 2*xt-1, poly);
	T_ORTHOPOL::Calc (order_cell[0], 2*yt-1, polx);		
	
	Tx pz = z*(1-z)*(1-z);
	
	for(int k = 0; k <= order_cell[0]-3; k++)
	  {
	    for(int i = 0; i <= k; i++)
	      { 
		Tx bubpik = pz * polx[i];
		for (int j = 0; j <= k; j++)
		  shape[ii++] = bubpik * poly[j];
	      }
	    pz *= 1-z;  
	  }
      }
  }
  

  template class H1HighOrderFiniteElement<1>;
  template class H1HighOrderFiniteElement<2>;
  template class H1HighOrderFiniteElement<3>;


  template class T_ScalarFiniteElement2<H1HighOrderFE<ET_SEGM>, ET_SEGM>;
  template class T_ScalarFiniteElement2<H1HighOrderFE<ET_TRIG>, ET_TRIG>;
  template class T_ScalarFiniteElement2<H1HighOrderFE<ET_QUAD>, ET_QUAD>;
  template class T_ScalarFiniteElement2<H1HighOrderFE<ET_TET>, ET_TET>;
  template class T_ScalarFiniteElement2<H1HighOrderFE<ET_PRISM>, ET_PRISM>;
  template class T_ScalarFiniteElement2<H1HighOrderFE<ET_HEX>, ET_HEX>;
  template class T_ScalarFiniteElement2<H1HighOrderFE<ET_PYRAMID>, ET_PYRAMID>;


  template class T_H1HighOrderFiniteElement<ET_SEGM>;
  template class T_H1HighOrderFiniteElement<ET_TRIG>;
  template class T_H1HighOrderFiniteElement<ET_QUAD>;
  template class T_H1HighOrderFiniteElement<ET_TET>;
  template class T_H1HighOrderFiniteElement<ET_PRISM>;
  template class T_H1HighOrderFiniteElement<ET_PYRAMID>;
  template class T_H1HighOrderFiniteElement<ET_HEX>;

  // template class H1HighOrderFE<ET_SEGM>;
  // template class H1HighOrderFE<ET_TET>;
}

