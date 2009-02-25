/*********************************************************************/
/* File:   h1hofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


 
#include <fem.hpp>
namespace ngfem
{
  using namespace ngfem;

  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetVertexNumbers (FlatArray<int> & avnums)
  {
    for (int i = 0; i < avnums.Size(); i++)
      vnums[i] = avnums[i];
  }
    
  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetOrderInner (int oi)
  {
    order_inner = INT<3> (oi,oi,oi); 
  }

  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetOrderInner (INT<3> oi)
  {
    order_inner = oi;
  }

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

  template<int D>
  void H1HighOrderFiniteElement<D> :: GetDofs (Array<Dof> & dofs) const
  {
    dofs.SetSize (0);
    
    for (int i = 0; i < ElementTopology::GetNNodes(eltype, NT_VERTEX); i++)
      dofs.Append (Dof (Node (NT_VERTEX, i), 0));

    for (int i = 0; i < ElementTopology::GetNNodes(eltype, NT_EDGE); i++)
      for (int j = 0; j < order_edge[i]-1; j++)
	dofs.Append (Dof (Node (NT_EDGE, i), j));

    for (int i = 0; i < ElementTopology::GetNNodes(eltype, NT_FACE); i++)
      {
        INT<2> p;
        if (D == 2)
          { p[0] = order_inner[0]; p[1] = order_inner[1]; }
        else
          p = order_face[i];
        
        int nf;
        if (D == 2 && eltype == ET_TRIG || 
            D == 3 && ElementTopology::GetFacetType (eltype, i) == ET_TRIG)
          nf = (p[0]-1)*(p[0]-2) / 2;
        else
          nf = (p[0]-1)*(p[1]-1);

	for (int j = 0; j < nf; j++)
	  dofs.Append (Dof (Node(NT_FACE, i), j));
      }

    int ni = 0;
    INT<3> p = order_inner;
    if (D == 3)
      switch (eltype)
        {
        case ET_TET:  if(p[0] > 3) ni = (p[0]-1)*(p[0]-2)*(p[0]-3)/6;   
          break;
        case ET_PRISM: if(p[0]>2 && p[2]>1) ni = (p[0]-1)*(p[0]-2)*(p[2]-1)/2;
          break;
        case ET_PYRAMID: if(p[0]>2) ni = (p[0]-1)*(p[0]-2)*(2*p[0]-3)/6;
          break;
        case ET_HEX: if(p[0]>1 && p[1] > 1 && p[2]>1) ni = (p[0]-1)*(p[1]-1)*(p[2]-1);
          break;
        }
    
    for (int j = 0; j < ni; j++)
      dofs.Append (Dof (Node(NT_CELL, 0), j));
  }


  template <ELEMENT_TYPE ET>
  void T_H1HighOrderFiniteElement<ET> :: 
  ComputeNDof()
  {
    ndof = ET_trait<ET>::N_VERTEX;
    
    for (int i = 0; i < ET_trait<ET>::N_EDGE; i++)
      ndof += order_edge[i] -1;
    
    for (int i = 0; i < ET_trait<ET>::N_FACE; i++)
      if (ET_trait<ET>::FaceType(i) == ET_TRIG)
        ndof += (order_face[i][0]-1)*(order_face[i][0]-2) / 2;
      else
        ndof += (order_face[i][0]-1)*(order_face[i][1]-1);
    
    switch (ET)
      {
      case ET_TET: 
        ndof += (order_inner[0]-1)*(order_inner[0]-2)*(order_inner[0]-3) / 6; 
        break;
      case ET_PRISM:
        ndof += (order_inner[0]-1)*(order_inner[0]-2)/2*(order_inner[2]-1);
        break;
      case ET_PYRAMID:
        ndof += (order_inner[0]-1)*(order_inner[0]-2)*(2*order_inner[0]-3) / 6; 
        break;
      case ET_HEX:
        ndof += (order_inner[0]-1)*(order_inner[1]-1)*(order_inner[2]-1);
        break;
      default:
        ;
      }
    
    order = 1;
    for (int i = 0; i < ET_trait<ET>::N_EDGE; i++)
      order = max(order, order_edge[i]);

    for (int i = 0; i < ET_trait<ET>::N_FACE; i++)
      order = max(order, max (order_face[i][0], order_face[i][1]));
    
    order = max(order, order_inner[0]);
    order = max(order, order_inner[1]);
    order = max(order, order_inner[2]);
  }


  template <ELEMENT_TYPE ET>
  void T_H1HighOrderFiniteElement<ET> :: 
  GetInternalDofs (Array<int> & idofs) const
  {
    int ni = 0;
    switch (ET)
      {
      case ET_TRIG: 
        ni = (order_face[0][0]-1)*(order_face[0][0]-2) / 2; 
        break;
      case ET_QUAD: 
        ni = (order_face[0][0]-1)*(order_face[0][1]-1);
        break;
      case ET_TET: 
        ni = (order_inner[0]-1)*(order_inner[0]-2)*(order_inner[0]-3) / 6;
        break;
      case ET_PRISM: 
        ni = (order_inner[0]-1)*(order_inner[0]-2)/2*(order_inner[2]-1);
        break;
      case ET_PYRAMID:
        ni = (order_inner[0]-1)*(order_inner[0]-2)*(2*order_inner[0]-3) / 6; 
        break;
      case ET_HEX: 
        ni = (order_inner[0]-1)*(order_inner[1]-1)*(order_inner[2]-1);
        break;
      }
    
    idofs.SetSize (0);
    for (int i = 0; i < ni; i++)
      idofs.Append (ndof-ni+i);
  }





  template <ELEMENT_TYPE ET>
  void T_H1HighOrderFiniteElement<ET> :: 
  CalcShape (const IntegrationPoint & ip, 
             FlatVector<> shape) const
  {
    double pt[DIM];
    for (int i = 0; i < DIM; i++) pt[i] = ip(i);
    static_cast<const H1HighOrderFE<ET>*> (this) -> T_CalcShape (pt, shape); 
  }

  template <ELEMENT_TYPE ET>
  void T_H1HighOrderFiniteElement<ET> :: 
  CalcDShape (const IntegrationPoint & ip, 
              FlatMatrix<> dshape) const
  {
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    ArrayMem<AutoDiff<DIM>,40> sds(ndof);
    static_cast<const H1HighOrderFE<ET>*> (this) -> T_CalcShape (adp, sds);

    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < DIM; j++)
	dshape(i,j) = sds[i].DValue(j);
  }









  /* *********************** Segment  **********************/



  H1HighOrderFE<ET_SEGM> :: H1HighOrderFE (int aorder)
    : T_H1HighOrderFiniteElement<ET_SEGM> (), ScalarFiniteElement<1>(ET_SEGM)
  {
    order_inner = aorder;
    ComputeNDof();
  }

  /*
  void H1HighOrderFE<ET_SEGM> :: ComputeNDof()
  {
    ndof = 2;
    ndof += order_inner[0]-1;
    order = order_inner[0];
  }
  */

  template<typename Tx, typename TFA>  
  void H1HighOrderFE<ET_SEGM> :: T_CalcShape (Tx hx[1], TFA & sds) const
  {
    Tx x = hx[0];
    Tx lami[2] = { x, 1-x };
 
    sds = 0.0;

    sds[0] = lami[0];
    sds[1] = lami[1];

    int ii = 2;
    int ee=1, es=0;
    if (vnums[es] > vnums[ee]) swap(es,ee);
    T_ORTHOPOL::Calc (order_inner[0], lami[ee]-lami[es], sds.Addr(ii));
  }



  /* *********************** Triangle  **********************/


  H1HighOrderFE<ET_TRIG> :: H1HighOrderFE(int aorder)
    : T_H1HighOrderFiniteElement<ET_TRIG> (), ScalarFiniteElement<2>(ET_TRIG)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    order_face[0] = INT<2> (aorder,aorder);
    for (int i = 0; i < 3; i++)
      order_edge[i] = aorder;

    ComputeNDof();
  }

  /*
  void H1HighOrderFE<ET_TRIG> :: ComputeNDof()
  {
    ndof = 3;

    for (int i = 0; i < 3; i++)
      ndof += order_edge[i] - 1;
    int oi = order_inner[0]; 
    ndof += (oi-1)*(oi-2) / 2;

    order = 1;
    for (int i = 0; i < 3; i++)
      if (order_edge[i] > order)
	order = order_edge[i];
    if (oi > order)
      order = oi;
  }
  */

  /*
  void H1HighOrderFE<ET_TRIG> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 3;
    for (int i = 0; i < 3; i++)
      base += order_edge[i] - 1;

    if(order_inner[0] > 2)
      {
	int ni = (order_inner[0]-1)*(order_inner[0]-2) / 2;
	for (int i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }
  */


  template<typename Tx, typename TFA>  
  void H1HighOrderFE<ET_TRIG> :: T_CalcShape (Tx hx[2], TFA & sds) const
  {
    Tx x = hx[0], y = hx[1];
    Tx lami[3] = { x, y, 1-x-y };

    sds = 0.0;

    for (int i = 0; i < 3; i++)
      sds[i] = lami[i];

    int ii = 3;
    
    // edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    for (int i = 0; i < 3; i++)
      if (order_edge[i] >= 2)
	{ 
	  int es = edges[i][0], ee = edges[i][1];
	  if (vnums[es] > vnums[ee]) swap (es, ee);

	  ii += T_ORTHOPOL::CalcTrigExt (order_edge[i], 
					 lami[ee]-lami[es], 1-lami[es]-lami[ee], 
                                         sds.Addr(ii));
	}
    
    int fav[3] = { 0, 1, 2 }; 
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	

    if (order_inner[0] >= 3)
      {
	ii += T_INNERSHAPES::Calc (order_inner[0], lami[fav[2]]-lami[fav[1]],
				   lami[fav[0]], sds.Addr(ii));
      }
  }


  /* *********************** Quadrilateral  **********************/

  H1HighOrderFE<ET_QUAD> :: H1HighOrderFE (int aorder)
  // : T_H1HighOrderFiniteElement<ET_QUAD> (), ScalarFiniteElement<2>(ET_QUAD)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    order_face[0] = INT<2> (aorder,aorder);
    for (int i = 0; i < 4; i++)
      order_edge[i] = aorder;
    ComputeNDof();
  }

  /*
  void H1HighOrderFE<ET_QUAD> :: ComputeNDof()
  {
    cout << "compute ndof, quad" << endl;
    order = 1;
    for (int i = 0; i < 4; i++)
      if (order_edge[i] > order)
	order = order_edge[i];

    for (int m = 0; m < 2; m++)
      if (order_inner[m] > order)
	order = order_inner[m];
    
    ndof = 4;
    for (int i = 0; i < 4; i++)
      if(order_edge[i]>1)
	ndof += order_edge[i] - 1;
    
    cout << "ndof, e = " << ndof << endl;
    if(order_inner[0] > 1 && order_inner[1] >1)
      ndof += (order_inner[0]-1)*(order_inner[1]-1);

    cout << "order_face = " <<order_face[0] << endl;
    cout << "order_inner = " <<order_inner << endl;

    cout << "ndof, f = " << ndof << endl;
  }

  void H1HighOrderFE<ET_QUAD> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 4;

    for (int i = 0; i < 4; i++)
      base += order_edge[i] - 1;

    if(order_inner[0] >= 2 && order_inner[1] >= 2)
      {
	int ni = (order_inner[0]-1)*(order_inner[1]-1);
	for (int i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }
  */

  template<typename Tx, typename TFA>  
  void H1HighOrderFE<ET_QUAD> :: T_CalcShape (Tx hx[2], TFA & sds) const
  {
    Tx x = hx[0], y = hx[1];
    Tx lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    
    // sds = 0.0;

    // vertex shapes
    for(int i=0; i < 4; i++) sds[i] = lami[i]; 
    int ii = 4;

    ArrayMem<Tx,20> polxi(order+1), poleta(order+1);
     
    // edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_QUAD);
    for (int i = 0; i < 4; i++)
      {
	int p = order_edge[i];
	int es = edges[i][0];
	int ee = edges[i][1];
	if (vnums[es] > vnums[ee]) swap (es, ee);

	Tx xi  = sigma[ee]-sigma[es]; // in [-1,1] 
	Tx eta = lami[ee]+lami[es];  // attention in [0,1]

	T_ORTHOPOL::Calc (p, xi, polxi);
	
	for (int j = 0; j <= p-2; j++)
	  sds[ii++] = eta * polxi[j];
      }    
    
    
    INT<2> p (order_face[0][0], order_face[0][1]);
    if (p[0] >= 2 && p[1] >= 2)
      {
	int fmax = 0; 
	for (int j = 1; j < 4; j++)
	  if (vnums[j] > vnums[fmax]) fmax = j;
	
	int f1 = (fmax+3)%4; 
	int f2 = (fmax+1)%4; 
	if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f1 > f2; 
	
	Tx xi = sigma[fmax]-sigma[f1]; 
	Tx eta = sigma[fmax]-sigma[f2]; 
	
	T_ORTHOPOL::Calc(p[0]+1, xi,polxi);
	T_ORTHOPOL::Calc(p[1]+1,eta,poleta);
	
	for (int k = 0; k <= p[0]-2; k++)
	  for (int j = 0; j <= p[1]-2; j++)
	    sds[ii++] = polxi[k] * poleta[j];
      }
  }



  /* *********************** Tetrahedron  **********************/


  H1HighOrderFE<ET_TET> :: H1HighOrderFE (int aorder)
    : T_H1HighOrderFiniteElement<ET_TET> (), ScalarFiniteElement<3>(ET_TET)
  {
    order_inner = INT<3> (aorder,aorder,aorder);

    for (int i = 0; i < 6; i++)
      order_edge[i] = aorder;
    for (int i = 0; i < 4; i++)
      order_face[i] = INT<2> (aorder,aorder);

    ComputeNDof();
  }

  /*
  void H1HighOrderFE<ET_TET> :: ComputeNDof()
  {
    ndof = 4;

    for (int i = 0; i < 6; i++)
      ndof += order_edge[i] - 1;
    for (int i = 0; i < 4; i++)
      ndof += (order_face[i][0]-1)*(order_face[i][0]-2) / 2;
    ndof += (order_inner[0]-1)*(order_inner[0]-2)*(order_inner[0]-3) / 6;

    order = 1;
    for (int i = 0; i < 6; i++)
      if (order_edge[i] > order)
	order = order_edge[i];
    for (int i = 0; i < 4; i++)
      if (order_face[i][0] > order)
	order = order_face[i][0];
    if (order_inner[0] > order)
      order = order_inner[0];
  }
  */
  /*
  void H1HighOrderFE<ET_TET> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 4;

    for (int i = 0; i < 6; i++)
      base += order_edge[i] - 1;
    for (int i = 0; i < 4; i++)
      base += (order_face[i][0]-1)*(order_face[i][0]-2) / 2;

    if(order_inner[0] > 2)
      {
	int ni = (order_inner[0]-1)*(order_inner[0]-2)*(order_inner[0]-3) / 6;
	for (int i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }
  */

  template<typename Tx, typename TFA>  
  void  H1HighOrderFE<ET_TET> :: T_CalcShape (Tx hx[3], TFA & sds) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];
    sds = 0.0;

    Tx lami[4] = { x, y, z, 1-x-y-z };

    ArrayMem<Tx, 20> polx(order+1), poly(order+1), polz(order+1); 

    // vertex shapes
    for (int i = 0; i < 4; i++)
      sds[i] = lami[i];
    int ii = 4; 

    // edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_TET);
    for (int i = 0; i < 6; i++)
      if (order_edge[i] >= 2)
	{
	  int es = edges[i][0], ee = edges[i][1];
	  if (vnums[es] > vnums[ee]) swap (es, ee);
	  
	  ii += T_ORTHOPOL::CalcTrigExt (order_edge[i], lami[ee]-lami[es], 
					 1-lami[es]-lami[ee], sds.Addr(ii) );

	}

    // face dofs
    const FACE * faces = ElementTopology::GetFaces (ET_TET);
    for (int i = 0; i < 4; i++)
      if (order_face[i][0] >= 3)
	{
	  int fav[3] = { faces[i][0], faces[i][1], faces[i][2] };
	  
	  if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	  if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	  if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
	  
	  int vop = 6 - fav[0] - fav[1] - fav[2];  	
	  
	  // Tx * hp = &sds[ii];
	  ii += T_FACESHAPES::Calc (order_face[i][0], 
				    lami[fav[2]]-lami[fav[1]],
				    lami[fav[0]], lami[vop],  sds.Addr(ii));
	}

    if (order_inner[0] >= 4)
      {
	// Tx * psds = &sds[ii];
	ii += T_INNERSHAPES::Calc (order_inner[0], x-(1-x-y-z), y, z, sds.Addr(ii) );
      }
  }




  /* *********************** Prism  **********************/


  H1HighOrderFE<ET_PRISM> :: H1HighOrderFE (int aorder)
    : T_H1HighOrderFiniteElement<ET_PRISM> (), ScalarFiniteElement<3>(ET_PRISM)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    for (int i = 0; i < 9; i++)
      order_edge[i] = aorder;
    for (int i = 0; i < 5; i++)
      order_face[i] = INT<2> (aorder,aorder);

    ComputeNDof();
  }

  /*
  void H1HighOrderFE<ET_PRISM> :: ComputeNDof()
  {
    ndof = 6;

    for (int i = 0; i < 9; i++)
      if(order_edge[i] > 1)
	ndof += order_edge[i] - 1;
    for (int i = 0; i < 2; i++)
      if(order_face[i][0] > 2)
	ndof += (order_face[i][0]-1)*(order_face[i][0]-2) / 2;
    for (int i = 2; i < 5; i++)
      if(order_face[i][0] > 1 && order_face[i][1]>1)
	ndof += (order_face[i][0]-1)*(order_face[i][1]-1);
    if(order_inner[0] > 2 && order_inner[2] >1) 
      ndof += (order_inner[0]-1)*(order_inner[0]-2)/2*(order_inner[2]-1);

    order = 1;
    for (int i = 0; i < 9; i++)
      if (order_edge[i] > order)
	order = order_edge[i];
    for (int i = 0; i < 5; i++)
      for (int k = 0; k < 2; k++)
	if (order_face[i][k] > order)
	  order = order_face[i][k];
    if (order_inner[0] > order)
      order = order_inner[0];
    if (order_inner[2] > order)
      order = order_inner[2];
  }
  

  void H1HighOrderFE<ET_PRISM> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 6;

    for (int i = 0; i < 9; i++)
      if(order_edge[i] > 1 ) 
	base += order_edge[i] - 1;

    for (int i = 0; i < 2; i++)
      {
	if(order_face[i][0] >2)
	  {
	    int ni = (order_face[i][0]-1)*(order_face[i][0]-2) / 2;
	    base += ni;
	  }
      }

    for (int i = 2; i < 5; i++)
      if(order_face[i][0]>1 && order_face[i][1]>1)
	base += (order_face[i][0]-1)*(order_face[i][1]-1);
    
    if(order_inner[0] > 2 && order_inner[2] > 1)
      {
	int ni = (order_inner[0]-1)*(order_inner[0]-2)/2*(order_inner[2]-1);
	for (int i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }
  */


  template<typename Tx, typename TFA>  
  void  H1HighOrderFE<ET_PRISM> :: T_CalcShape (Tx hx[3], TFA & sds) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];
    Tx lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    Tx muz[6]  = { 1-z, 1-z, 1-z, z, z, z };
    
    sds = 0.0;

    // vertex shapes
    for (int i = 0; i < 6; i++)
      sds[i] = lami[i] * muz[i];

    int ii = 6;

    // horizontal edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_PRISM);
    for (int i = 0; i < 6; i++)
      if (order_edge[i] >= 2)
	{
	  int p = order_edge[i];
	  int es = edges[i][0], ee = edges[i][1];
	  if (vnums[es] > vnums[ee]) swap (es, ee);

	  Tx xi = lami[ee]-lami[es]; 
	  Tx eta = lami[es]+lami[ee]; 
	  
	  Tx * hp = &sds[ii];
	  T_ORTHOPOL::CalcTrigExt (p, xi, 1-eta, hp);
	  
	  for (int k = 0; k < p-1; k++)
	    sds[ii++] *= muz[ee];  
	}
    
    // vertical edges
    for (int i = 6; i < 9; i++)
      if (order_edge[i] >= 2)
	{
	  int p = order_edge[i];
	  int es = edges[i][0], ee = edges[i][1];
	  if (vnums[es] > vnums[ee]) swap (es, ee);
	  
	  Tx * hp = &sds[ii];
	  T_ORTHOPOL::Calc (p, muz[ee]-muz[es], hp);
	  
	  for (int j = 0; j < p-1; j++)
	    sds[ii++] *= lami[ee];
	}
    

    ArrayMem<Tx,20> polx(order+1) /* , poly(order+1) */ , polz(order+1);

    const FACE * faces = ElementTopology::GetFaces (ET_PRISM);
    // trig face dofs
    for (int i = 0; i < 2; i++)
      if (order_face[i][0] >= 3)
	{
	  int fav[3] = { faces[i][0], faces[i][1], faces[i][2] };
	  if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	  if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	  if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
	  
	  Tx * pface = &sds[ii];
	  int ndf = T_TRIGFACESHAPES::Calc 
	    (order_face[i][0], lami[fav[2]]-lami[fav[1]], lami[fav[0]], pface);
	  
	  for (int j = 0; j < ndf; j++)
	    pface[j] *= muz[fav[2]];
	  ii += ndf;
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
	  
	  Tx xi = lami[faces[i][fmax]] - lami[faces[i][ftrig]]; 
	  Tx eta = 1-lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
	  Tx zeta = muz[faces[i][fmax]]-muz[faces[i][fz]]; 
	  
	  T_ORTHOPOL::CalcTrigExt (pp, xi, eta, polx);  
	  T_ORTHOPOL::Calc (pp, zeta, polz);	
	  
	  // global x-direction is towards second-largest vertex (JS)
	  if (vnums[faces[i][ftrig]] > vnums[faces[i][fz]])  
	    for (int k = 0; k < p[0]-1; k++)
	      for (int j = 0; j < p[1]-1; j++)
		sds[ii++] = polx[k] * polz[j];
	  else
	    for (int j = 0; j < p[0]-1; j++)
	      for (int k = 0; k < p[1]-1; k++)
		sds[ii++] = polx[k] * polz[j];
	}
    
    // volume dofs:
    if (order_inner[0] > 2 && order_inner[2] > 1)
      {
	ArrayMem<Tx,20> pol_trig((order_inner[0]-1)*(order_inner[0]-2)/2);
	int nf = T_TRIGFACESHAPES::Calc (order_inner[0], x-y, 1-x-y, pol_trig);

	T_ORTHOPOL:: Calc(order_inner[2], 2*z-1, polz);
	for (int i = 0; i < nf; i++)
	  for (int k = 0; k < order_inner[2]-1; k++)
	    sds[ii++] = pol_trig[i] * polz[k];
      }
  }





  /* *********************** Hex  **********************/

  H1HighOrderFE<ET_HEX> :: H1HighOrderFE (int aorder)
    : T_H1HighOrderFiniteElement<ET_HEX> (), ScalarFiniteElement<3>(ET_HEX)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    for (int i = 0; i < 12; i++)
      order_edge[i] = aorder;
    for (int i = 0; i < 6; i++)
      order_face[i] = INT<2> (aorder,aorder);
    ComputeNDof();
  }

  /*
  void H1HighOrderFE<ET_HEX> :: ComputeNDof()
  {
    ndof = 8;

    for (int i = 0; i < 12; i++)
      ndof += order_edge[i] - 1;
    for (int i = 0; i < 6; i++)
      ndof +=  (order_face[i][0]-1)*(order_face[i][1]-1);
    ndof += (order_inner[0]-1)*(order_inner[1]-1)*(order_inner[2]-1);

    order = 1;
    for (int i = 0; i < 12; i++)
      if (order_edge[i] > order)
	order = order_edge[i];
    for (int i = 0; i < 6; i++)
      for (int k = 0; k < 2; k++)
	if (order_face[i][k] > order)
	  order = order_face[i][k];
    for(int m=0; m<3; m++)
      if (order_inner[m] > order)
	order = order_inner[m];
  }

  void H1HighOrderFE<ET_HEX> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 8;

    for (int i = 0; i < 12; i++)
      base += order_edge[i] - 1;
    for (int i = 0; i < 6; i++)
      base += (order_face[i][0]-1)*(order_face[i][1]-1);

    if(order_inner[0] >= 2 && order_inner[1] >= 2 && order_inner[2] >= 2)
      {
	int ni = (order_inner[0]-1)*(order_inner[1]-1)*(order_inner[2]-1);
	for (int i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }
  */

  template<typename Tx, typename TFA>  
  void  H1HighOrderFE<ET_HEX> :: T_CalcShape (Tx hx[3], TFA & sds) const
  { 
    Tx x = hx[0], y = hx[1], z = hx[2];

    sds = 0.0;

    Tx lami[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
		(1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
    Tx sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
		 (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 

    // vertex shapes
    for(int i=0; i<8; i++) sds[i] = lami[i]; 
    int ii = 8;

    ArrayMem<Tx,20> polx(order+1), poly(order+1), polz(order+1);
    
    // edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_HEX);
    // const POINT3D * points = ElementTopology::GetVertices (ET_HEX);
    for (int i = 0; i < 12; i++)
      if (order_edge[i] >= 2)
	{
	  int p = order_edge[i];
	  int es = edges[i][0], ee = edges[i][1];
	  if (vnums[es] > vnums[ee]) swap (es, ee);
	  
	  Tx xi = sigma[ee]-sigma[es]; 
	  Tx lam_e = lami[es]+lami[ee];
	  T_ORTHOPOL::Calc (p, xi, polx);
	  
	  for (int j = 0; j < p-1; j++)
	    sds[ii++] = lam_e * polx[j];
	}
     
    const FACE * faces = ElementTopology::GetFaces (ET_HEX);
    for (int i = 0; i < 6; i++)
      if (order_face[i][0] >= 2 && order_face[i][1] >= 2)
	{
	  int px = order_face[i][0];
	  int py = order_face[i][1];
	  int p = max2(px, py);

	  Tx lam_f = 0;
	  for (int j = 0; j < 4; j++)
	    lam_f += lami[faces[i][j]];
	  
	  int fmax = 0;
	  for (int j=1; j<4; j++) 
	    if (vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;  
	  
	  int f1 = faces[i][(fmax+3)%4];
	  int f2 = faces[i][(fmax+1)%4]; 
	  fmax = faces[i][fmax]; 
	  
	  if(vnums[f2]>vnums[f1]) swap(f1,f2);  // fmax > f1 > f2 
	  
	  Tx xi  = sigma[fmax] - sigma[f1]; 
	  Tx eta = sigma[fmax] - sigma[f2];
	  
	  T_ORTHOPOL::Calc (p, xi, polx);	
	  T_ORTHOPOL::Calc (p, eta, poly);	
	  
	  for (int k = 0; k < px-1; k++) 
	    {
	      Tx pxl = polx[k] * lam_f;
	      for (int j = 0; j < py-1; j++) 
		sds[ii++]= pxl * poly[j];
	    }
	}
    
    // volume dofs:
    if (order_inner[0] >= 2 && order_inner[1] >= 2 && order_inner[2] >= 2)
      {
	T_ORTHOPOL::Calc (order_inner[0], 2*x-1, polx);
	T_ORTHOPOL::Calc (order_inner[1], 2*y-1, poly);
	T_ORTHOPOL::Calc (order_inner[2], 2*z-1, polz);
	for (int i = 0; i < order_inner[0]-1; i++)
	  for (int j = 0; j < order_inner[1]-1; j++)
	    {
	      Tx pxy = polx[i] * poly[j];
	      for (int k = 0; k < order_inner[2]-1; k++)
		sds[ii++] = pxy * polz[k];
	    }
      }
  }

  /* ******************************** Pyramid  ************************************ */

  H1HighOrderFE<ET_PYRAMID> :: H1HighOrderFE (int aorder)
    : T_H1HighOrderFiniteElement<ET_PYRAMID>(), ScalarFiniteElement<3>(ET_PYRAMID)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    for (int i = 0; i < 8; i++)
      order_edge[i] = aorder;
    for (int i = 0; i < 5; i++)
      order_face[i] = INT<2> (aorder,aorder);

    ComputeNDof();
  }

  void H1HighOrderFE<ET_PYRAMID> :: ComputeNDof()
  {
    ndof = 5;

    for (int i = 0; i < 8; i++)
      ndof += order_edge[i] - 1;
    for (int i = 0; i < 4; i++)
      ndof += (order_face[i][0]-1)*(order_face[i][0]-2) / 2;
    ndof += (order_face[4][0]-1)*(order_face[4][1]-1);
    ndof += (order_inner[0]-1)*(order_inner[0]-2)*(2*order_inner[0]-3) / 6; 


    order = 1;
    for (int i = 0; i < 8; i++)
      if (order_edge[i] > order)
	order = order_edge[i];
    for (int i = 0; i < 5; i++)
      for (int k = 0; k < 2; k++)
	if (order_face[i][k] > order)
	  order = order_face[i][k];
    if (order_inner[0] > order)
      order = order_inner[0];
  }


  void H1HighOrderFE<ET_PYRAMID> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 5;

    for (int i = 0; i < 8; i++)
      base += order_edge[i] - 1;
    for (int i = 0; i < 4; i++)
      base += (order_face[i][0]-1)*(order_face[i][0]-2) / 2;
    base += (order_face[4][0]-1)*(order_face[4][1]-1);

    if(order_inner[0] > 2)
      {
	int ni = (order_inner[0]-1)*(order_inner[0]-2)*(2*order_inner[0]-3) / 6;
	
	for (int i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }

  template<typename Tx, typename TFA>  
  void  H1HighOrderFE<ET_PYRAMID> :: T_CalcShape (Tx hx[3], TFA & sds) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];
    sds = 0.0;

    if (z == 1.) z -= 1e-10;

    Tx xt = x / (1-z);
    Tx yt = y / (1-z);
    
    // Tx mux[4] = { 1-xt, xt, xt, 1-xt };
    // Tx muy[4] = { 1-yt, 1-yt, yt, yt };

    Tx sigma[4]  = { (1-xt)+(1-yt), xt+(1-yt), xt+yt, (1-xt)+yt };
    Tx lambda[4] = { (1-xt)*(1-yt), xt*(1-yt), xt*yt, (1-xt)*yt };

    for (int i = 0; i < 4; i++)  
      sds[i] = lambda[i] * (1-z);
    sds[4] = z;

    int ii = 5;


    //horizontal edge dofs 
    const EDGE * edges = ElementTopology::GetEdges (ET_PYRAMID);
    for (int i = 0; i < 4; i++)
      if (order_edge[i] >= 2)
	{
	  int p = order_edge[i];

	  int es = edges[i][0], ee = edges[i][1];
	  if (vnums[es] > vnums[ee]) swap (es, ee);

	  Tx lam = sigma[ee]-sigma[es]; 
	  Tx lam_edge = lambda[es] + lambda[ee];
	  
	  Tx * pedge = &sds[ii];
	  T_ORTHOPOL::CalcTrigExt (p, lam*(1-z), z, pedge);

	  for (int j = 0; j <= p-2; j++)
	    sds[ii++] *= lam_edge;
	}
    
    // vertical edges
    for (int i = 4; i < 8; i++) 
      if (order_edge[i] >= 2)
	{
	  int es = edges[i][0], ee = edges[i][1];
	  if (vnums[es] > vnums[ee]) swap (es, ee);
	  
	  Tx * pedge = &sds[ii];
	  ii += T_ORTHOPOL::CalcTrigExt (order_edge[i], sds[ee]-sds[es],  
					 1-sds[es]-sds[ee], pedge);
	}


    ArrayMem<Tx,20> polx(order+1), poly(order+1), polz(order+1);
 
    const FACE * faces = ElementTopology::GetFaces (ET_PYRAMID);
    const POINT3D * points = ElementTopology::GetVertices (ET_PYRAMID);

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
	 
	  Tx * pface = &sds[ii];
	  
	  int ndf = T_TRIGSHAPES::Calc 
	    (p, bary[fav[2]]-bary[fav[1]], bary[fav[0]], pface);
	  
	  for (int j = ii; j < ii+ndf; j++)
	    sds[j] *= lam_face; 
	  ii += ndf;
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
	    sds[ii++]= polx[k] * poly[j] * fac; 
      }

    //#ifdef SABINE  
    if (order_inner[0] >= 3)
      {
	T_ORTHOPOL::Calc (order_inner[0], 2*xt-1, poly);
	T_ORTHOPOL::Calc (order_inner[0], 2*yt-1, polx);		
	
	Tx pz = z*(1-z)*(1-z);
	
	for(int k=0; k<= order_inner[0]-3; k++)
	  {
	    for(int i = 0; i <= k; i++)
	      { 
		Tx bubpik = pz * polx[i];
		for (int j = 0; j <= k; j++)
		  sds[ii++] = bubpik * poly[j];
	      }
	    pz *= 1-z;  
	  }
      }
    //#else

//     if (order_inner[0] >= 3)
//       {
// 	LegendrePolynomial (order_inner[0]-2, -1+2*xt, polx);
// 	LegendrePolynomial (order_inner[0]-2, -1+2*yt, poly); 
	
// 	Tx bub = xt*(1-xt)*yt*(1-yt)*z*(1-z)*(1-z);
	
// 	for(int i = 0; i <= order_inner[0]-3; i++)
// 	  {
// 	    for (int j = 0; j <= i; j++)
// 	      {
// 		 Tx bubpj = bub * polx[j];
// 		 for (int k = 0; k <= i; k++)
// 		  sds[ii++] = bubpj * poly[k];
// 	      }
// 	    bub *= 1-z;  
// 	  }
//       }
//#endif 
  }


  template class H1HighOrderFiniteElement<1>;
  template class H1HighOrderFiniteElement<2>;
  template class H1HighOrderFiniteElement<3>;

  // template class T_H1HighOrderFiniteElement<ET_SEGM>;
  // template class H1HighOrderFE<ET_SEGM>;
}

