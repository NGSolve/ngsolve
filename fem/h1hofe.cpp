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
  H1HighOrderFiniteElement<D> ::   
  H1HighOrderFiniteElement (ELEMENT_TYPE aeltype)
    : NodalFiniteElement<D> (D, aeltype, -1, -1) 
  { 
    for (int i = 0; i < 8; i++)
      vnums[i] = i;

    for (int i = 0; i < 8; i++)
      order_vertex[i] = NodalFiniteElement<D>::order;

    augmented = 0;
    NodalFiniteElement<D>::needs_update = true;
  }
   
  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    for (int i = 0; i < avnums.Size(); i++)
      vnums[i] = avnums[i];
  }
    
  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetOrderInner (int oi)
  {
    order_inner = INT<3> (oi,oi,oi); 
    NodalFiniteElement<D>::needs_update = true;
    // ComputeNDof();
  }

  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetOrderInner (INT<3> oi)
  {
    order_inner = oi;
    NodalFiniteElement<D>::needs_update = true;
    // ComputeNDof();
  }

  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetAugmented (int aa)
  {
    augmented = aa;
    NodalFiniteElement<D>::needs_update = true;
    // ComputeNDof();
  }

  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetOrderFace (FlatArray<int> & of)
  {
    for (int i = 0; i < of.Size(); i++)
      order_face[i] = INT<2>(of[i],of[i]);
    NodalFiniteElement<D>::needs_update = true;
    // ComputeNDof();
  }

  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetOrderFace (FlatArray<INT<2> > & of)
  {
    for (int i = 0; i < of.Size(); i++)
      order_face[i] = of[i];
    NodalFiniteElement<D>::needs_update = true;
    // ComputeNDof();
  }

  
  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetOrderEdge (FlatArray<int> & oe)
  {
    for (int i = 0; i < oe.Size(); i++)
      order_edge[i] = oe[i];
    NodalFiniteElement<D>::needs_update = true;
    // ComputeNDof();
  }

  template<int D>
  void H1HighOrderFiniteElement<D>::
  SetOrderVertex (FlatArray<int> & ov)
  {
    for (int i = 0; i < ov.Size(); i++)
      order_vertex[i] = ov[i];
    NodalFiniteElement<D>::needs_update = true;
    // ComputeNDof();
  }


  template<int D>
  void H1HighOrderFiniteElement<D> :: GetDofs (Array<Dof> & dofs) const
  {
    dofs.SetSize (0);
    
    for (int i = 0; i < ElementTopology::GetNNodes(NodalFiniteElement<D>::eltype, NT_VERTEX); i++)
      dofs.Append (Dof (Node (NT_VERTEX, i), 0));

    for (int i = 0; i < ElementTopology::GetNNodes(NodalFiniteElement<D>::eltype, NT_EDGE); i++)
      for (int j = 0; j < order_edge[i]-1; j++)
	dofs.Append (Dof (Node (NT_EDGE, i), j));

    for (int i = 0; i < ElementTopology::GetNNodes(NodalFiniteElement<D>::eltype, NT_FACE); i++)
      {
        INT<2> p;
        if (D == 2)
          { p[0] = order_inner[0]; p[1] = order_inner[1]; }
        else
          p = order_face[i];
        
        int nf;
        if (D == 2 && NodalFiniteElement<D>::eltype == ET_TRIG || 
            D == 3 && 
            ElementTopology::GetFacetType (NodalFiniteElement<D>::eltype, i) == ET_TRIG)
          nf = (p[0]-1)*(p[0]-2) / 2;
        else
          nf = (p[0]-1)*(p[1]-1);

	for (int j = 0; j < nf; j++)
	  dofs.Append (Dof (Node(NT_FACE, i), j));
      }

    int ni = 0;
    INT<3> p = order_inner;
    if (D == 3)
      switch (NodalFiniteElement<D>::eltype)
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


  /* *********************** Segment  **********************/


  template <class T_ORTHOPOL>
  H1HighOrderSegm<T_ORTHOPOL> :: H1HighOrderSegm (int aorder)
    : H1HighOrderFiniteElement<1>(ET_SEGM)
  {
    order_inner = aorder;
    ComputeNDof();
  }


  template <class T_ORTHOPOL>
  void H1HighOrderSegm<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 2;
    if (augmented == 1) ndof += 2;
    if (augmented == 2) ndof = 2 * order_inner[0];

    ndof += order_inner[0]-1;
    order = order_inner[0];

    needs_update = 0;
  }



  template <class T_ORTHOPOL> 
  void H1HighOrderSegm<T_ORTHOPOL> ::CalcShape( const IntegrationPoint & ip,
						FlatVector<> shape) const
  {
    double x = ip(0);
    double lami[2] = { x, 1-x };
 
    shape = 0.0;

    shape(0) = lami[0];
    shape(1) = lami[1];

    int ii = 2;
    
    // vertex shapes
    if (augmented == 1)
      for (int i = 0; i < 2; i++)
	shape(ii++) = T_VERTEXSHAPES::Calc (order, lami[i]);
    else
      if (augmented == 2)
	for (int i = 0; i < 2; i++)
	  {
	    ArrayMem<double,20> pol(order+1);
	    LegendrePolynomial (order-1, -1+2*lami[i], pol);	
	    for (int j = 1; j < order; j++)
	      shape(ii++) = pol[j] * lami[i];
	  }
    
    int ee=1, es=0;
    if (vnums[es] > vnums[ee]) swap(es,ee);
    double * pedge = &shape(ii);
    T_ORTHOPOL::Calc (order_inner[0], lami[ee]-lami[es], pedge);
  }
 
  template <class T_ORTHOPOL>
  void H1HighOrderSegm<T_ORTHOPOL> :: CalcDShape (const IntegrationPoint & ip, 
						  FlatMatrix<> dshape) const
  {
    AutoDiff<1> x(ip(0),0); 
    AutoDiff<1> lami[2] = {x,1-x}; 
    dshape(0,0) = lami[0].DValue(0);
    dshape(1,0) = lami[1].DValue(0);
    
    int ee=1,es=0;
    if (vnums[es] > vnums[ee]) swap(es,ee);
    
    int ii=2; 
    ArrayMem<AutoDiff<1>,10> pedge(order_inner[0]);  
    T_ORTHOPOL::Calc (order_inner[0], lami[ee]-lami[es], pedge);
    for(int i=0;i<order_inner[0]-1;i++)
      dshape(ii++,0)=pedge[i].DValue(0); 
  }



  /* *********************** Triangle  **********************/

  template <class T_ORTHOPOL>
  H1HighOrderTrig<T_ORTHOPOL> :: H1HighOrderTrig(int aorder)
    : H1HighOrderFiniteElement<2> (ET_TRIG)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    for (int i = 0; i < 3; i++)
      order_edge[i] = aorder;

    ComputeNDof();
  }


  template <class T_ORTHOPOL>
  void H1HighOrderTrig<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 3;
    if (augmented == 1) ndof += 3;
    if (augmented == 2) ndof = 3 * order_inner[0];

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

    needs_update = 0;
  }



  template <class T_ORTHOPOL>
  void H1HighOrderTrig<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 3;
    if (augmented == 1) base *= 2;
    if (augmented == 2) base *= order_inner[0];

    int i;
    for (i = 0; i < 3; i++)
      base += order_edge[i] - 1;

    if(order_inner[0] > 2)
      {
	int ni = (order_inner[0]-1)*(order_inner[0]-2) / 2;
	for (i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }



  template <class T_ORTHOPOL>
  void H1HighOrderTrig<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
						 FlatVector<> shape) const
  {
    // T_CalcShape<double, double, FlatVector<> > (ip(0), ip(1), shape); 
    T_CalcShape (ip(0), ip(1), shape); 
  }

  template <class T_ORTHOPOL>
  void H1HighOrderTrig<T_ORTHOPOL> :: CalcDShape (const IntegrationPoint & ip, 
						  FlatMatrix<> dshape) const
  {
    AutoDiff<2> x(ip(0), 0);
    AutoDiff<2> y(ip(1), 1);

    ArrayMem<AutoDiff<2>,40> sds(ndof);
    // T_CalcShape<AutoDiff<2>, AutoDiff<2>, FlatArray<AutoDiff<2> > > (x,y,sds);
    T_CalcShape (x,y,sds);

    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 2; j++)
	dshape(i,j) = sds[i].DValue(j);
  }



  template<typename T_ORTHOPOL> 
  template<typename Tx, typename Ty, typename TFA>  
  void  H1HighOrderTrig<T_ORTHOPOL> :: T_CalcShape (Tx x, Ty y, TFA & sds) const
  {
    if (needs_update) 
      throw Exception ("H1HighOrderTrig needs update");


    Tx lami[3] = { x, y, 1-x-y };

    sds = 0.0;

    for (int i = 0; i < 3; i++)
      sds[i] = lami[i];

    int ii = 3;

    // vertex shapes
    if (augmented == 1)
      for (int i = 0; i < 3; i++)
	sds[ii++] = T_VERTEXSHAPES::Calc (order, lami[i]);
    else
      if (augmented == 2)
	for (int i = 0; i < 3; i++)
	  {
	    ArrayMem<Tx,20> pol(order+1);
	    LegendrePolynomial (order-1, -1+2*lami[i], pol);	
	    for (int j = 1; j < order; j++)
	      sds[ii++] = pol[j] * lami[i];
	  }
    
    // edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    for (int i = 0; i < 3; i++)
      if (order_edge[i] >= 2)
	{ 
	  int es = edges[i][0], ee = edges[i][1];
	  if (vnums[es] > vnums[ee]) swap (es, ee);

          /*
	  Tx * hp = &sds[ii];
	  ii += T_ORTHOPOL::CalcTrigExt (order_edge[i], 
					 lami[ee]-lami[es], 1-lami[es]-lami[ee], hp);
          */
	  ii += T_ORTHOPOL::CalcTrigExt (order_edge[i], 
					 lami[ee]-lami[es], 1-lami[es]-lami[ee], 
                                         sds.Addr(ii));
          // Range(ii, ii+order_edge[i]-1));
	}
    
    int fav[3] = { 0, 1, 2 }; 
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	

    if (order_inner[0] >= 3)
      {
	// Tx * hp = &sds[ii];
	ii += T_INNERSHAPES::Calc (order_inner[0], lami[fav[2]]-lami[fav[1]],
				   lami[fav[0]], sds.Addr(ii));
        // Range (ii, ii+order_inner[0]-1));
      }
  }

  /* *********************** Quadrilateral  **********************/

  template <class T_ORTHOPOL>
  H1HighOrderQuad<T_ORTHOPOL> :: H1HighOrderQuad (int aorder)
    : H1HighOrderFiniteElement<2> (ET_QUAD)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    for (int i = 0; i < 4; i++)
      order_edge[i] = aorder;
    ComputeNDof();
  }


  template <class T_ORTHOPOL>
  void H1HighOrderQuad<T_ORTHOPOL> :: ComputeNDof()
  {
    order = 1;
    for (int i = 0; i < 4; i++)
      if (order_edge[i] > order)
	order = order_edge[i];

    for (int m = 0; m < 2; m++)
      if (order_inner[m] > order)
	order = order_inner[m];
    
    ndof = 4;
    if (augmented == 1) ndof *= 2;
    if (augmented == 2) ndof *= order;

    for (int i = 0; i < 4; i++)
      if(order_edge[i]>1)
	ndof += order_edge[i] - 1;
    
    if(order_inner[0] > 1 && order_inner[1] >1)
      ndof += (order_inner[0]-1)*(order_inner[1]-1);

    needs_update = 0;
  }

  template <class T_ORTHOPOL>
  void H1HighOrderQuad<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 4;
    if (augmented == 1) base *= 2;
    if (augmented == 2) base *= order_inner[0];

    int i;
    for (i = 0; i < 4; i++)
      base += order_edge[i] - 1;

    if(order_inner[0] >= 2 && order_inner[1] >= 2)
      {
	int ni = (order_inner[0]-1)*(order_inner[1]-1);
	for (i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }

  template <class T_ORTHOPOL>
  void H1HighOrderQuad<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
						 FlatVector<> shape) const
  {  
    T_CalcShape<double, double, FlatVector<> > (ip(0), ip(1), shape); 

    /*
      ArrayMem<AutoDiff<2>, 40> sds(ndof);
      CalcShapeDShape (ip, sds);
      for (int i = 0; i < ndof; i++)
      shape(i) = sds[i].Value();
    */
  }

  template <class T_ORTHOPOL>
  void H1HighOrderQuad<T_ORTHOPOL> :: CalcDShape (const IntegrationPoint & ip, 
						  FlatMatrix<> dshape) const
  {
    AutoDiff<2> x(ip(0), 0);
    AutoDiff<2> y(ip(1), 1);

    ArrayMem<AutoDiff<2>,40> sds(ndof);
    T_CalcShape<AutoDiff<2>, AutoDiff<2>, FlatArray<AutoDiff<2> > >
      (x,y,sds);
    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 2; j++)
	dshape(i, j) = sds[i].DValue(j);
    /*
      ArrayMem<AutoDiff<2>,25> sds(ndof);
      CalcShapeDShape (ip, sds);
      for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 2; j++)
      dshape(i, j) = sds[i].DValue(j);
    */
  }

  template<typename T_ORTHOPOL> 
  template<typename Tx, typename Ty, typename TFA>  
  void  H1HighOrderQuad<T_ORTHOPOL> :: T_CalcShape (Tx x, Ty y, TFA & sds) const
  {
    Tx lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    
    sds = 0.0;

    // vertex shapes
    for(int i=0; i < 4; i++) sds[i] = lami[i]; 
    int ii = 4;

    if (augmented == 1)
      for (int i = 0; i < 4; i++)
	sds[ii++] = T_VERTEXSHAPES::Calc (order, sds[i]);
    else
      if (augmented == 2)
	for (int i = 0; i < 4; i++)
	  {
	    ArrayMem<Tx,20> pol(order+1);
	    LegendrePolynomial (order-1, -1+2*sds[i], pol);	
	    for (int j = 1; j < order; j++)
	      sds[ii++] = pol[j] * sds[i];
	  }

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
    
    int inci; 
    if (order_inner[0] >= 2 && order_inner[1] >= 2)
      {
	inci = ii ; 
	INT<2> p (order_inner[0], order_inner[1]);

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

  template <class T_ORTHOPOL>
  H1HighOrderTet<T_ORTHOPOL> :: H1HighOrderTet (int aorder)
    : H1HighOrderFiniteElement<3> (ET_TET)
  {
    order_inner = INT<3> (aorder,aorder,aorder);

    for (int i = 0; i < 6; i++)
      order_edge[i] = aorder;
    for (int i = 0; i < 4; i++)
      order_face[i] = INT<2> (aorder,aorder);

    ComputeNDof();
  }


  template <class T_ORTHOPOL>
  void H1HighOrderTet<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 4;
    if (augmented == 1) ndof *= 2;
    if (augmented == 2) ndof *= order_inner[0];

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

    needs_update = 0;
  }


  template <class T_ORTHOPOL>
  void H1HighOrderTet<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 4;
    if (augmented == 1) base *= 2;
    if (augmented == 2) base *= order_inner[0];

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





  template <class T_ORTHOPOL>
  void H1HighOrderTet<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
						FlatVector<> shape) const
  {
    T_CalcShape<double, double, double, FlatVector<> > (ip(0), ip(1), ip(2), shape); 
  }
  
  template <class T_ORTHOPOL>
  void H1HighOrderTet<T_ORTHOPOL> :: CalcDShape (const IntegrationPoint & ip, 
						 FlatMatrix<> dshape) const
  {
    AutoDiff<3> x(ip(0), 0);
    AutoDiff<3> y(ip(1), 1);
    AutoDiff<3> z(ip(2), 2);

    ArrayMem<AutoDiff<3>,40> sds(ndof);
    T_CalcShape<AutoDiff<3>, AutoDiff<3>, AutoDiff<3>, FlatArray<AutoDiff<3> > >
      (x,y,z,sds);
    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 3; j++)
	dshape(i, j) = sds[i].DValue(j);
  }


  template <class T_ORTHOPOL>
  void H1HighOrderTet<T_ORTHOPOL> :: 
  CalcMappedDShape (const BaseSpecificIntegrationPoint & bsip, 
                    FlatMatrix<> dshape) const
  {
    const SpecificIntegrationPoint<3,3> & sip = 
      static_cast<const SpecificIntegrationPoint<3,3> &> (bsip);
    
    const IntegrationPoint & ip = sip.IP();
    AutoDiff<3> x(ip(0)), y(ip(1)), z(ip(2));
    for (int j = 0; j < 3; j++)
      {
        x.DValue(j) = sip.GetJacobianInverse()(0,j);
        y.DValue(j) = sip.GetJacobianInverse()(1,j);
        z.DValue(j) = sip.GetJacobianInverse()(2,j);
      }

    ArrayMem<AutoDiff<3>,40> sds(ndof);
    T_CalcShape<AutoDiff<3>, AutoDiff<3>, AutoDiff<3>, FlatArray<AutoDiff<3> > >
      (x,y,z,sds);
    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 3; j++)
	dshape(i, j) = sds[i].DValue(j);
  }
   

  template<typename T_ORTHOPOL> 
  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void  H1HighOrderTet<T_ORTHOPOL> :: T_CalcShape (Tx x, Ty y, Tz z, TFA & sds) const
  {
    sds = 0.0;

    Tx lami[4] = { x, y, z, 1-x-y-z };

    ArrayMem<Tx, 20> polx(order+1), poly(order+1), polz(order+1); 

    // vertex shapes
    for (int i = 0; i < 4; i++)
      sds[i] = lami[i];
    int ii = 4; 

    if (augmented == 1)
      for (int i = 0; i < 4; i++)
	sds[ii++] = T_VERTEXSHAPES::Calc (order, lami[i]);
    else
      if (augmented == 2)
	for (int i = 0; i < 4; i++)
	  {
	    ArrayMem<Tx,20> pol(order+1);
	    LegendrePolynomial (order-1, -1+2*lami[i], pol);	
	    for (int j = 1; j < order; j++)
	      sds[ii++] = pol[j] * lami[i];
	  }
 
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


  template <class T_ORTHOPOL>
  H1HighOrderPrism<T_ORTHOPOL> :: H1HighOrderPrism (int aorder)
    : H1HighOrderFiniteElement<3> (ET_PRISM)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    for (int i = 0; i < 9; i++)
      order_edge[i] = aorder;
    for (int i = 0; i < 5; i++)
      order_face[i] = INT<2> (aorder,aorder);

    plate = 0;

    ComputeNDof();
  }


  template <class T_ORTHOPOL>  
  void H1HighOrderPrism<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 6;
    if (augmented == 1) ndof *= 2;
    if (augmented == 2) ndof *= order_inner[0];

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

    needs_update = 0;
  }
  

  template <class T_ORTHOPOL>
  void H1HighOrderPrism<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 6;
    if (augmented == 1) base *= 2;
    if (augmented == 2) base *= order_inner[0];

    for (int i = 0; i < 9; i++)
      if(order_edge[i] > 1 ) 
	base += order_edge[i] - 1;

    for (int i = 0; i < 2; i++)
      {
	if(order_face[i][0] >2)
	  {
	    int ni = (order_face[i][0]-1)*(order_face[i][0]-2) / 2;
	    
	    if (plate)
	      for (int j = 0; j < ni; j++)
		idofs.Append (base+j);
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



  template <class T_ORTHOPOL>
  void H1HighOrderPrism<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
						  FlatVector<> shape) const
  {
    T_CalcShape<double, double, double, FlatVector<> > (ip(0), ip(1), ip(2), shape); 
  }

  template <class T_ORTHOPOL>
  void H1HighOrderPrism<T_ORTHOPOL> :: CalcDShape (const IntegrationPoint & ip, 
						   FlatMatrix<> dshape) const
  {
    AutoDiff<3> x(ip(0), 0);
    AutoDiff<3> y(ip(1), 1);
    AutoDiff<3> z(ip(2), 2);

    ArrayMem<AutoDiff<3>,40> sds(ndof);
    T_CalcShape<AutoDiff<3>, AutoDiff<3>, AutoDiff<3>, FlatArray<AutoDiff<3> > >
      (x,y,z,sds);
    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 3; j++)
	dshape(i, j) = sds[i].DValue(j);
  }


  template<typename T_ORTHOPOL> 
  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void  H1HighOrderPrism<T_ORTHOPOL> :: T_CalcShape (Tx x, Ty y, Tz z, TFA & sds) const
  {
    Tx lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    Tx muz[6]  = { 1-z, 1-z, 1-z, z, z, z };
    
    sds = 0.0;

    // vertex shapes
    for (int i = 0; i < 6; i++)
      sds[i] = lami[i] * muz[i];

    int ii = 6;

    if (augmented == 1)
      for (int i = 0; i < 6; i++)
	sds[ii++] = T_VERTEXSHAPES::Calc (order, sds[i]);
    else
      if (augmented == 2)
	for (int i = 0; i < 6; i++)
	  {
	    ArrayMem<Tx,20> pol(order+1);
	    LegendrePolynomial (order-1, -1+2*sds[i], pol);	
	    for (int j = 1; j < order; j++)
	      sds[ii++] = pol[j] * sds[i];
	  }

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

  template <class T_ORTHOPOL>
  H1HighOrderHex<T_ORTHOPOL> :: H1HighOrderHex (int aorder)
    : H1HighOrderFiniteElement<3> (ET_HEX)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    for (int i = 0; i < 12; i++)
      order_edge[i] = aorder;
    for (int i = 0; i < 6; i++)
      order_face[i] = INT<2> (aorder,aorder);
    ComputeNDof();
  }

  template <class T_ORTHOPOL>
  void H1HighOrderHex<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 8;

    if (augmented == 1) ndof *= 2;
    if (augmented == 2) ndof *= order_inner[0];

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

    needs_update = 0;
  }

  template <class T_ORTHOPOL>
  void H1HighOrderHex<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 8;
    if (augmented == 1) base *= 2;
    if (augmented == 2) base *= order_inner[0];

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



  
  template <class T_ORTHOPOL>
  void H1HighOrderHex<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
						FlatVector<> shape) const
  {
    T_CalcShape<double, double, double, FlatVector<> > (ip(0), ip(1), ip(2), shape); 
  }

  template <class T_ORTHOPOL>
  void H1HighOrderHex<T_ORTHOPOL> :: CalcDShape (const IntegrationPoint & ip, 
						 FlatMatrix<> dshape) const
  {
    AutoDiff<3> x(ip(0), 0);
    AutoDiff<3> y(ip(1), 1);
    AutoDiff<3> z(ip(2), 2);

    ArrayMem<AutoDiff<3>,40> sds(ndof);
    T_CalcShape<AutoDiff<3>, AutoDiff<3>, AutoDiff<3>, FlatArray<AutoDiff<3> > >
      (x,y,z,sds);
    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 3; j++)
	dshape(i, j) = sds[i].DValue(j);
  }




  template<typename T_ORTHOPOL> 
  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void  H1HighOrderHex<T_ORTHOPOL> :: T_CalcShape (Tx x, Ty y, Tz z, TFA & sds) const
  { 

    sds = 0.0;

    Tx lami[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
		(1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
    Tx sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
		 (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 

    // vertex shapes
    for(int i=0; i<8; i++) sds[i] = lami[i]; 
    int ii = 8;

    if (augmented == 1)
      for (int i = 0; i < 8; i++)
	sds[ii++] = T_VERTEXSHAPES::Calc (order, sds[i]);
    else
      if (augmented == 2)
	for (int i = 0; i < 8; i++)
	  {
	    ArrayMem<Tx,20> pol(order+1);
	    LegendrePolynomial (order-1, -1+2*sds[i], pol);	
	    for (int j = 1; j < order; j++)
	      sds[ii++] = pol[j] * sds[i];
	  }

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

  template <class T_ORTHOPOL>
  H1HighOrderPyramid<T_ORTHOPOL> :: H1HighOrderPyramid (int aorder)
    : H1HighOrderFiniteElement<3> (ET_PYRAMID)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    for (int i = 0; i < 8; i++)
      order_edge[i] = aorder;
    for (int i = 0; i < 5; i++)
      order_face[i] = INT<2> (aorder,aorder);

    ComputeNDof();
  }


  template <class T_ORTHOPOL>
  void H1HighOrderPyramid<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 5;

    if (augmented == 1) ndof *= 2;
    if (augmented == 2) ndof *= order_inner[0];

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

    needs_update = 0;
  }


  template <class T_ORTHOPOL>
  void H1HighOrderPyramid<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);

    int base = 5;
    if (augmented == 1) base *= 2;
    if (augmented == 2) base *= order_inner[0];

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

  template <class T_ORTHOPOL>
  void H1HighOrderPyramid<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
						    FlatVector<> shape) const
  {
    T_CalcShape<double, double, double, FlatVector<> > (ip(0), ip(1), ip(2), shape); 
    /*
      ArrayMem<AutoDiff<3>, 40> sds(ndof);
      CalcShapeDShape (ip, sds);
      for (int i = 0; i < ndof; i++)
      shape(i) = sds[i].Value();
    */
  }

  template <class T_ORTHOPOL>  
  void H1HighOrderPyramid<T_ORTHOPOL> :: CalcDShape (const IntegrationPoint & ip, 
						     FlatMatrix<> dshape) const
  {
    AutoDiff<3> x(ip(0), 0);
    AutoDiff<3> y(ip(1), 1);
    AutoDiff<3> z(ip(2), 2);

    ArrayMem<AutoDiff<3>,40> sds(ndof);
    T_CalcShape<AutoDiff<3>, AutoDiff<3>, AutoDiff<3>, FlatArray<AutoDiff<3> > >
      (x,y,z,sds);
    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 3; j++)
	dshape(i, j) = sds[i].DValue(j);
    /*
      ArrayMem<AutoDiff<3>,40> sds(ndof);
      CalcShapeDShape (ip, sds);
      for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 3; j++)
      dshape(i, j) = sds[i].DValue(j);
    */
  }



  template<typename T_ORTHOPOL> 
  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void  H1HighOrderPyramid<T_ORTHOPOL> :: T_CalcShape (Tx x, Ty y, Tz z, TFA & sds) const
  {
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

    if (augmented == 1)
      for (int i = 0; i < 5; i++)
	sds[ii++] = T_VERTEXSHAPES::Calc (order, sds[i]);
    else
      if (augmented == 2)
	for (int i = 0; i < 5; i++)
	  {
	    ArrayMem<Tx,20> pol(order+1);
	    LegendrePolynomial (order-1, -1+2*sds[i], pol);	
	    for (int j = 1; j < order; j++)
	      sds[ii++] = pol[j] * sds[i];
	  }



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

  template class H1HighOrderSegm<IntegratedLegendreMonomialExt>;
  template class H1HighOrderTrig<IntegratedLegendreMonomialExt>;
  template class H1HighOrderQuad<IntegratedLegendreMonomialExt>;
  template class H1HighOrderTet<IntegratedLegendreMonomialExt>;
  template class H1HighOrderPrism<IntegratedLegendreMonomialExt>;
  template class H1HighOrderHex<IntegratedLegendreMonomialExt>;
  template class H1HighOrderPyramid<IntegratedLegendreMonomialExt>;

  /*
  template class H1HighOrderSegm<TrigExtensionMin>;
  template class H1HighOrderTrig<TrigExtensionMin>;
  template class H1HighOrderQuad<TrigExtensionMin>;
  template class H1HighOrderTet<TrigExtensionMin>;
  template class H1HighOrderPrism<TrigExtensionMin>;
  template class H1HighOrderHex<TrigExtensionMin>;
  template class H1HighOrderPyramid<TrigExtensionMin>;
  */

  /*
  template class H1HighOrderSegm<TrigExtensionOptimal>;
  template class H1HighOrderTrig<TrigExtensionOptimal>;
  template class H1HighOrderQuad<TrigExtensionOptimal>;
  template class H1HighOrderTet<TrigExtensionOptimal>;
  template class H1HighOrderPrism<TrigExtensionOptimal>;
  template class H1HighOrderHex<TrigExtensionOptimal>;
  template class H1HighOrderPyramid<TrigExtensionOptimal>;
  */




  /*
  template class H1HighOrderSegm<TrigExtensionMonomial>;
  template class H1HighOrderSegm<TrigExtensionOptimal>;
  template class H1HighOrderSegm<TrigExtensionMin>;

  template class H1HighOrderTrig<TrigExtensionMonomial>;
  template class H1HighOrderTrig<TrigExtensionOptimal>;
  template class H1HighOrderTrig<TrigExtensionMin>;

  template class H1HighOrderQuad<TrigExtensionMonomial>;
  template class H1HighOrderQuad<TrigExtensionOptimal>;
  template class H1HighOrderQuad<TrigExtensionMin>;

  template class H1HighOrderTet<TrigExtensionMonomial>;
  template class H1HighOrderTet<TrigExtensionOptimal>;
  template class H1HighOrderTet<TrigExtensionMin>;

  template class H1HighOrderPrism<TrigExtensionMonomial>;
  template class H1HighOrderPrism<TrigExtensionOptimal>;
  template class H1HighOrderPrism<TrigExtensionMin>;

  template class H1HighOrderHex<TrigExtensionMonomial>;
  template class H1HighOrderHex<TrigExtensionOptimal>;
  template class H1HighOrderHex<TrigExtensionMin>;

  template class H1HighOrderPyramid<TrigExtensionMonomial>;
  template class H1HighOrderPyramid<TrigExtensionOptimal>;
  template class H1HighOrderPyramid<TrigExtensionMin>;
  */

}

