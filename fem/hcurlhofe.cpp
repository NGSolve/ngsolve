#include <fem.hpp>    

namespace ngfem
{
  using namespace ngfem;


  
  // hv.DValue() = (grad u) x (grad v) 
  inline AutoDiff<3> Cross (const AutoDiff<3> & u,
			    const AutoDiff<3> & v)
  {
    AutoDiff<3> hv;
    hv.Value() = 0.0;
    hv.DValue(0) = u.DValue(1)*v.DValue(2)-u.DValue(2)*v.DValue(1);
    hv.DValue(1) = u.DValue(2)*v.DValue(0)-u.DValue(0)*v.DValue(2);
    hv.DValue(2) = u.DValue(0)*v.DValue(1)-u.DValue(1)*v.DValue(0);
    return hv;
  }
  
  /* inline AutoDiff<2>  Cross (const AutoDiff<2> & u, const AutoDiff<2> & v) 
  {
    AutoDiff<2> hv; 
    hv.Value() = 0.0; 
    hv.DValue(0) =  u.DValue(0)*v.DValue(1)-u.DValue(1)*v.DValue(0);
    return hv; 
    } */ 


  //------------------------------------------------------------------------
  // HCurlHighOrderFiniteElement
  //------------------------------------------------------------------------

  template <int D>
  HCurlHighOrderFiniteElement<D> ::
  HCurlHighOrderFiniteElement (ELEMENT_TYPE aeltype)
    : HCurlFiniteElement<D> (aeltype, -1, -1) 
  { 
    for (int i = 0; i < 8; i++)
      vnums[i] = i;

    usegrad_cell = 1;
    for(int i=0; i<6; i++)
      usegrad_face[i] = 1;
    for(int i=0; i<12; i++)
      usegrad_edge[i] = 1;

  }

  template <int D>
  void HCurlHighOrderFiniteElement<D>::
  SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    for (int i = 0; i < avnums.Size(); i++)
      vnums[i] = avnums[i];
  }
    
  template <int D>
  void HCurlHighOrderFiniteElement<D>::
  SetOrderInner (int oi)
  {
    order_inner = INT<3> (oi,oi,oi); 
  }

  template <int D>
  void HCurlHighOrderFiniteElement<D>::
  SetOrderInner (INT<3> oi)
  {
    order_inner = oi; 
  }

  template <int D>
  void HCurlHighOrderFiniteElement<D>::
  SetOrderFace (FlatArray<int> & of)
  {
    for (int i = 0; i < of.Size(); i++)
      order_face[i] = INT<2> (of[i],of[i]);
  }

  template <int D>
  void HCurlHighOrderFiniteElement<D>::
  SetOrderFace (FlatArray<INT<2> > & of)
  {
    for (int i = 0; i < of.Size(); i++)
	order_face[i] = of[i];
  }

  template <int D>
  void HCurlHighOrderFiniteElement<D>::
  SetOrderEdge (FlatArray<int> & oen)
  {
    for (int i = 0; i < oen.Size(); i++)
      order_edge[i] = oen[i];
  }


  template <int D> 
  void HCurlHighOrderFiniteElement<D>:: 
  SetUsegradEdge (FlatArray<int> & uge)
  {
    for (int i=0; i<uge.Size(); i++) usegrad_edge[i]=uge[i]; 
  }
  
  template <int D> 
  void HCurlHighOrderFiniteElement<D>:: 
  SetUsegradFace (FlatArray<int> & ugf)
  {
    for (int i=0; i<ugf.Size(); i++) usegrad_face[i]=ugf[i]; 
  }
  
  template <int D> 
  void HCurlHighOrderFiniteElement<D>:: 
  SetUsegradCell (int ugc)
  {
    usegrad_cell=ugc; 
  }
  
  template <int D> 
  void HCurlHighOrderFiniteElement<D>:: 
  SetAugmented (int aa)
  {
    augmented = aa;
    ComputeNDof();
  }

  template <int D> 
  void HCurlHighOrderFiniteElement<D>:: 
  SetOrderVertex (FlatArray<int> & ov)
  {
    for (int i = 0; i < ov.Size(); i++)
      order_vertex[i] = ov[i];
    ComputeNDof();
  }


  template <int D> 
  void HCurlHighOrderFiniteElement<D>:: 
  PrintInfo() const
  {
    (*testout) << "order_inner " << order_inner << " order_face ";
    for(int i=0; i<6; i++)
      (*testout) << order_face[i] << " ";
    (*testout) << "order_edge ";
    for(int i=0; i<12; i++)
      (*testout) << order_edge[i] << " ";
    (*testout) << "usegrad_cell " << usegrad_cell << " usgrad_face ";
    for(int i=0; i<6; i++)
      (*testout) << usegrad_face[i] << " ";
    (*testout) << "usegrad_edge ";
    for(int i=0; i<12; i++)
      (*testout) << usegrad_edge[i] << " ";
  }

  //------------------------------------------------------------------------
  // HCurlHighOrderSegm
  //------------------------------------------------------------------------
  
  template <class T_ORTHOPOL> 
  HCurlHighOrderSegm<T_ORTHOPOL> :: HCurlHighOrderSegm (int aorder)
    : HCurlHighOrderFiniteElement<1>(ET_SEGM)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    usegrad_cell = 1; 
    augmented =0; 
    for(int i=0;i<2;i++) order_vertex[i] = aorder+1; 
    ComputeNDof(); 
  }

  template <class T_ORTHOPOL>
  void HCurlHighOrderSegm<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof =1; 
    if(augmented ==1) ndof+=2; 
    if(usegrad_cell)
      ndof += order_inner[0];    
    order = order_inner[0];
    order++; // integration order 
  }

  
  template <class T_ORTHOPOL>
  void HCurlHighOrderSegm<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
				     FlatMatrixFixWidth<1> shape) const
  {
    AutoDiff<1> x (ip(0),0); 
    AutoDiff<1> lami[2] = {1-x,x}; 
    
    int es = 0, ee =1;
    if (vnums[es] > vnums[ee]) swap(es,ee);  
    AutoDiff<1> xi = lami[ee] - lami[es]; 
    
    shape = 0.0; //!
    
    // Nedelec0-shapes
    shape(0,0) = 0.5*(xi.DValue(0)); 

    int ii = 1; 
    if(augmented==1)
      {
	ArrayMem<double,20> pol(order); 
	for(int i=0;i<2;i++)	
	  {
	    AutoDiff<1> v = T_VERTEXSHAPES::Calc(order_vertex[i],lami[i]);
	    shape(ii++,0)= v.DValue(0); 
	  }
      }
      
    if (order_inner[0]>=1)
      if(usegrad_cell)
	{ 
	  ArrayMem<AutoDiff<1>, 10> pol_xi(order_inner[0]+2);
	  T_ORTHOPOL::Calc (order_inner[0]+1, xi,pol_xi);  
	  
	  for (int j = 0; j < order_inner[0]; j++)
	    shape(ii++,0) = pol_xi[j].DValue(0); 
	}    
  }


   //------------------------------------------------------------------------
  // HCurlHighOrderTrig
  //------------------------------------------------------------------------
  
  template <class T_ORTHOPOL>
  HCurlHighOrderTrig<T_ORTHOPOL> :: HCurlHighOrderTrig (int aorder)
    : HCurlHighOrderFiniteElement<2>(ET_TRIG)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    for (int i = 0; i < 3; i++)
      order_edge[i] = aorder; 
    
    augmented=0; 
    for(int i=0;i<3;i++) order_vertex[i] = aorder+1; 
    ComputeNDof();
  }
  
  template <class T_ORTHOPOL>
  void HCurlHighOrderTrig<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 3; // Nedelec 
    if(augmented==1) ndof+=3; 

    int i; 
    for(i=0; i<3; i++)
      ndof += usegrad_edge[i]*order_edge[i]; 
 	
    if (order_inner[0] > 1)  
      ndof += ((usegrad_cell + 1) * order_inner[0] +2) * (order_inner[0]-1) /2; //      - (order_inner[0]-1);
	
    order = 1; // max(order_edges_tang,order_inner);  
    for (i = 0; i < 3; i++)
      {
	if (order_edge[i] > order)
	  order = order_edge[i];
      }
   
    if (order_inner[0] > order) 
      order = order_inner[0];
   
    if(order==0) order++; 
  }

  template <class T_ORTHOPOL>
  void HCurlHighOrderTrig<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {

    if ( discontinuous )
      {
	int ni = ndof;
	idofs.SetSize (ndof);
	for ( int i = 0; i < ni; i++ )
	  idofs[i] = i;
	return;
      }

    idofs.SetSize (0);
    int i, base = 3; // Nedelec 
    if(augmented == 1 ) base +=3; 

    for (i = 0; i < 3; i++)
      base += usegrad_edge[i]*order_edge[i];
    
    if(order_inner[0] > 1)
      {
	int ni = ((usegrad_cell + 1) * order_inner[0] +2) * (order_inner[0]-1) /2; //        - (order_inner[0]-1);
	
	for (i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }



  template <class T_ORTHOPOL>
  void HCurlHighOrderTrig<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
						    FlatMatrixFixWidth<2> shape) const
  {
    
    AutoDiff<2> x(ip(0),0); 
    AutoDiff<2> y(ip(1),1); 
    AutoDiff<2> lami[3] = {x,y,1-x-y};

    shape = 0.0; 
    
    int i, j, k, l;
    int ii = 3; 
    
    if(augmented==1)
      for(i=0;i<3;i++,ii++)
	{
	  AutoDiff<2> v= T_VERTEXSHAPES::Calc (order_vertex[i], lami[i]);
	  for(l=0;l<2;l++) shape(ii,l) = v.DValue(l); 
	}

    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    ArrayMem<double,30> mem(2*(order+2)); 
    ArrayMem<AutoDiff<2>,10> adpol1(order),adpol2(order);	
	
    for (i = 0; i < 3; i++)
      {
	int es = edges[i][0]; 
	int ee = edges[i][1];
	
	if (vnums[es] > vnums[ee])  swap (es, ee);

	//Nedelec low order edge shape function 
	for(l=0; l<2; l++) 
	  shape(i,l) = lami[ee].DValue(l)*lami[es].Value() - 
	    lami[es].DValue(l)*lami[ee].Value();  


	int p = order_edge[i]; 
	//HO-Edge shapes (Gradient Fields)   
	if(p>0 && usegrad_edge[i]) 
	  { 
	    AutoDiff<2> xi = lami[ee] - lami[es]; 
	    AutoDiff<2> eta = 1 - lami[ee] - lami[es]; 
	    T_ORTHOPOL::CalcTrigExt(p+1,xi, eta, adpol1); 
	   
	    for(j=0; j< p;j++,ii++) 
	      for(l=0;l<2;l++)
		shape(ii,l) = adpol1[j].DValue(l);
	  }
      }   
     
    int p = order_inner[0];      
    //Inner shapes (Face) 
    if(p>1) 
      {
	int fav[3] = {0,1,2}; 
	//Sort vertices ... v(f0) < v(f1) < v(f2) 
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	  

	AutoDiff<2> xi  = lami[fav[2]]-lami[fav[1]];
	AutoDiff<2> eta = lami[fav[0]]; 

	T_INNERSHAPES::CalcSplitted(p+1, xi, eta, adpol1,adpol2);
	
	// gradients:
	if(usegrad_cell)
	  for (j = 0; j < p-1; j++)
	    for (k = 0; k < p-1-j; k++, ii++)
	      for (l = 0; l < 2; l++)
		shape(ii, l) = adpol1[j].DValue(l) * adpol2[k].Value()
		  + adpol1[j].Value() * adpol2[k].DValue(l); 
	
	// other combination
	for (j = 0; j < p-1; j++)
	  for (k = 0; k < p-1-j; k++, ii++)
	    for (l = 0; l < 2; l++)
	       shape(ii, l) = adpol1[j].DValue(l) * adpol2[k].Value()
		 - adpol1[j].Value() * adpol2[k].DValue(l); 
		
	
	double ned0[2]; 
	for(j=0;j<2;j++) 
	  ned0[j] = lami[fav[1]].Value()*lami[fav[2]].DValue(j) - lami[fav[2]].Value()*lami[fav[1]].DValue(j) ;  

	// rec_pol * Nedelec0 
	for (j = 0; j < p-1; j++, ii++)
	  for (l = 0; l < 2; l++)
	    shape(ii,l) = adpol2[j].Value() * ned0[l]; 

	// (*testout) << "nedelec trig, shape = " << endl << shape << endl;
      }
  }
 


  template <class T_ORTHOPOL>
  void HCurlHighOrderTrig<T_ORTHOPOL> :: CalcCurlShape (const IntegrationPoint & ip, 
		       			FlatMatrixFixWidth<1> shape) const
  { 
    shape = 0.0; 
    
    AutoDiff<2> x(ip(0),0); 
    AutoDiff<2> y(ip(1),1); 
    AutoDiff<2> lami[3] = {x,y,1-x-y};
    
    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    
    int i, j, k, l, ii;
    ii=3; 
    if(augmented==1) ii = 6; 
      
    for (i = 0; i < 3; i++)
      {
	int es = edges[i][0]; 
	int ee = edges[i][1];
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	//Nedelec low order edge shape function 2*grad(lam_s) x grad(lam_e)
	shape(i,0) = 2*(lami[es].DValue(0)*lami[ee].DValue(1)
			 - lami[es].DValue(1)*lami[ee].DValue(0));  
	// Curl(Edge tangential shapes)  = 0.0; 
 	if(usegrad_edge[i]) ii += order_edge[i];
      }    
      
    int p = order_inner[0]; 
    
    // Inner shapes (Face) 
    
    int fav[3] ={0,1,2};    
    //Sort vertices  first edge op minimal vertex 
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	  
    
    AutoDiff<2> xi  = lami[fav[2]]-lami[fav[1]];
    AutoDiff<2> eta = lami[fav[0]]; // 1- lami[f2] - lami[f1]; 

    ArrayMem<AutoDiff<2>,10> adpol1(order),adpol2(order); 
    T_INNERSHAPES::CalcSplitted(p+1,xi,eta,adpol1,adpol2); 
    
    // curl(gradients) = 0.0
    if(usegrad_cell) ii += (p-1)*p/2; //  - (p-1);

    // curl(other combinations) = 2 grad(v) x grad(u) 
    for (j = 0; j < p-1; j++)
      for (k = 0; k < p-1-j; k++)
	shape(ii++, 0) = 2*(adpol2[k].DValue(0) * adpol1[j].DValue(1) 
			    - adpol2[k].DValue(1) * adpol1[j].DValue(0));

    // curl(adpol2 * Nedelec0)   
    double ned0[2] = { lami[fav[1]].Value()*lami[fav[2]].DValue(0) - 
	 	       lami[fav[2]].Value()*lami[fav[1]].DValue(0), 
		       lami[fav[1]].Value()*lami[fav[2]].DValue(1) - 
		       lami[fav[2]].Value()*lami[fav[1]].DValue(1) };
        
    double curlned0 = 2*(lami[fav[1]].DValue(0)*lami[fav[2]].DValue(1)
			 - lami[fav[1]].DValue(1)*lami[fav[2]].DValue(0)); 
 
    for (j = 0; j <= p-2; j++)  
      shape(ii++,0) = curlned0 * adpol2[j].Value() + adpol2[j].DValue(0)*ned0[1] - 
 	adpol2[j].DValue(1)*ned0[0]; 
  }


  //------------------------------------------------------------------------
  // HCurlHighOrderQuad
  //------------------------------------------------------------------------
  
  template <class T_ORTHOPOL>  
  HCurlHighOrderQuad<T_ORTHOPOL> :: HCurlHighOrderQuad (int aorder)
    : HCurlHighOrderFiniteElement<2>(ET_QUAD)
  {
    order_inner = INT<3> (aorder,aorder,aorder);
    for (int i = 0; i < 4; i++)
	order_edge[i] = aorder; 

    for(int i=0;i<4;i++) order_vertex[i] = aorder+1; 
    augmented = 0; 
    ComputeNDof();
  }
  
  template <class T_ORTHOPOL>
  void HCurlHighOrderQuad<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 4; // Nedelec 
    if(augmented==1) ndof +=4; 
    
    for(int i=0;i<4;i++)
      ndof += usegrad_edge[i]*order_edge[i];
          
    if (order_inner[0] >= 0 && order_inner[1] >= 0) 
    	ndof += (usegrad_cell+ 1) * order_inner[0] * order_inner[1] +
	  order_inner[0] + order_inner[1]; 
		  
    order = 0; // max(order_edges_tang,order_inner);  
    for (int i = 0; i < 4; i++)
      if (order_edge[i] > order) order = order_edge[i];
   
    for(int k=0; k<2; k++)
      if (order_inner[k] > order) // order_edge_normal = order_inner; 
	order = order_inner[k];

    
    order++; // for integration order
  }


  template <class T_ORTHOPOL>
  void HCurlHighOrderQuad<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    if ( discontinuous )
      {
	int ni = ndof;
	idofs.SetSize (ndof);
	for ( int i = 0; i < ni; i++ )
	  idofs[i] = i;
	return;
      }

    idofs.SetSize (0);

    int i, base = 4; // Nedelec 
    if(augmented == 1 ) base +=4; 

    for (i = 0; i < 4; i++)
      base += usegrad_edge[i]*order_edge[i];
    
    if(order_inner[0] >= 0 && order_inner[1] >= 0)
      {
	int ni = (usegrad_cell+ 1) * order_inner[0] * order_inner[1] +
	  order_inner[0] + order_inner[1];
	
	for (i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }



  template <class T_ORTHOPOL>
  void HCurlHighOrderQuad<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
				       FlatMatrixFixWidth<2> shape) const
  {
    int i, j, k, l;
    AutoDiff<2> x (ip(0),0);
    AutoDiff<2> y (ip(1),1);

    AutoDiff<2> lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
 
    const EDGE * edges = ElementTopology::GetEdges (ET_QUAD);

    shape = 0;
    int ii = 4;

    ArrayMem<AutoDiff<2>, 10> pol_xi(order+2);
    ArrayMem<AutoDiff<2>, 10> pol_eta(order+2);
    
    if(augmented==1)
      for(i=0;i<4;i++,ii++)
	{
	  AutoDiff<2> v= T_VERTEXSHAPES::Calc (order_vertex[i], lami[i]);
	  for(l=0;l<2;l++) shape(ii,l) = v.DValue(l); 
	}


    // edges
    for (i = 0; i < 4; i++)
      {
	int p = order_edge[i]; 
	int es = edges[i][0];
	int ee = edges[i][1];
	if (vnums[es] > vnums[ee]) swap (es, ee);

	AutoDiff<2> xi  = sigma[ee]-sigma[es];
	AutoDiff<2> lam_e = lami[ee]+lami[es];  // attention in [0,1]

	// Nedelec0-shapes
	for(l=0; l<2; l++) 
	  shape(i,l) = 0.5*xi.DValue(l)*lam_e.Value(); 
	
	// High Order edges ... Gradient fields 
	if(usegrad_edge[i])
	  {
	    T_ORTHOPOL::Calc (p+1, xi,pol_xi);  
	    for (j = 0; j < p; j++, ii++)
	      for(l=0; l<2; l++)
		shape(ii,l) =  lam_e.DValue(l)*pol_xi[j].Value() + 
		  lam_e.Value()*pol_xi[j].DValue(l);
	  }
      }
     
    INT<2> p(order_inner[0],order_inner[1]);
    int fmax = 0; 
    for (j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax])
	fmax = j;
    
    int f1 = (fmax+3)%4; 
    int f2 = (fmax+1)%4; 
    if(vnums[f1] > vnums[f2]) swap(f1,f2);  // fmax > f2 > f1; 

    AutoDiff<2> xi = sigma[fmax]-sigma[f2];  // in [-1,1]
    AutoDiff<2> eta = sigma[fmax]-sigma[f1]; // in [-1,1]
    
    T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
    T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);
    
    //Gradient fields 
    if(usegrad_cell)
      for (k = 0; k < p[0]; k++)
	for (j= 0; j < p[1]; j++, ii++)
	  for(l=0; l<2; l++)
	    shape(ii,l) = pol_xi[k].Value()*pol_eta[j].DValue(l) 
	      + pol_xi[k].DValue(l)*pol_eta[j].Value(); 
    
    //Rotation of Gradient fields 
    for (k = 0; k < p[0]; k++)
      for (j= 0; j < p[1]; j++, ii++)
	for(l=0; l<2; l++)
	  shape(ii,l) = pol_xi[k].DValue(l)*pol_eta[j].Value() 
	    - pol_xi[k].Value()*pol_eta[j].DValue(l); 

    //Missing ones 
    for(j=0; j< p[0]; j++, ii++)
      for(l=0;l<2;l++)
	shape(ii,l) = 0.5*eta.DValue(l)*pol_xi[j].Value();
    for(j=0; j< p[1]; j++, ii++)
      for(l=0;l<2;l++)
	shape(ii,l) = 0.5*xi.DValue(l)*pol_eta[j].Value(); 
    
    return;
  }
  
  template <class T_ORTHOPOL>
  void HCurlHighOrderQuad<T_ORTHOPOL> :: CalcCurlShape (const IntegrationPoint & ip,
 							FlatMatrixFixWidth<1> shape) const
  {
    int i, j, k, l;
    AutoDiff<2> x (ip(0),0);
    AutoDiff<2> y (ip(1),1);

    AutoDiff<2> lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
 
    const EDGE * edges = ElementTopology::GetEdges (ET_QUAD);

    shape = 0.0;
    int ii = 4;

    ArrayMem<AutoDiff<2>, 10> pol_xi(order+2);
    ArrayMem<AutoDiff<2>, 10> pol_eta(order+2);
    
    if(augmented==1) ii+=4; 
    
    // edges
    for (i = 0; i < 4; i++)
      {
	int p = order_edge[i]; 
	int es = edges[i][0];
	int ee = edges[i][1];
	if (vnums[es] > vnums[ee]) swap (es, ee);

	AutoDiff<2> xi  = sigma[ee]-sigma[es];
	AutoDiff<2> lam_e = lami[ee]+lami[es];  

	// Nedelec0-shapes
	for(l=0; l<2; l++) 
	  shape(i,0) = 0.5*(lam_e.DValue(0)*xi.DValue(1) - 
			    lam_e.DValue(1)*xi.DValue(0));

	// High Order edges ... Gradient fields 
	if(usegrad_edge[i]) ii+=p; 
      }
     
    INT<2> p(order_inner[0],order_inner[1]);
    int fmax = 0; 
    for (j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax])
	fmax = j;
    
    int f1 = (fmax+3)%4; 
    int f2 = (fmax+1)%4; 
    if(vnums[f1] > vnums[f2]) swap(f1,f2);  // fmax > f2 > f1; 

    AutoDiff<2> xi = sigma[fmax]-sigma[f2];  // in [-1,1]
    AutoDiff<2> eta = sigma[fmax]-sigma[f1]; // in [-1,1]
    
    T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
    T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);
    
    //Gradient fields 
    if(usegrad_cell) ii+=p[0]*p[1]; 
        
    //Rotation of Gradient fields
    for (k = 0; k < p[0]; k++)
      for (j= 0; j < p[1]; j++, ii++)
	shape(ii,0)=2*(pol_eta[j].DValue(0)*pol_xi[k].DValue(1)-pol_eta[j].DValue(1)*pol_xi[k].DValue(0));

    //Missing ones 
    for(j=0; j< p[0]; j++, ii++)
      shape(ii,0) = 0.5*(pol_xi[j].DValue(0)*eta.DValue(1)-pol_xi[j].DValue(1)*eta.DValue(0)); 
    
    for(j=0; j< p[1]; j++, ii++)
      shape(ii,0) = 0.5*(pol_eta[j].DValue(0)*xi.DValue(1)-pol_eta[j].DValue(1)*xi.DValue(0)); 

    return;
  }

  
  //------------------------------------------------------------------------
  // HCurlHighOrderTet
  //------------------------------------------------------------------------
 
  template <class T_ORTHOPOL>
  HCurlHighOrderTet<T_ORTHOPOL> :: HCurlHighOrderTet (int aorder)
    : HCurlHighOrderFiniteElement<3>(ET_TET)
  {
    nv = 4; 
    ned = 6; 
    nf = 4; 

    int i;
    order_inner = INT<3> (aorder,aorder,aorder);
    for (i = 0; i < ned; i++)
      order_edge[i] = aorder;
    
    for (i=0; i<nf; i++)  
      order_face[i] = INT<2> (aorder,aorder); 

    augmented = 0; 
    for(int i=0;i<4;i++) order_vertex[i] = aorder+1; 
    ComputeNDof();
  }
  
  template <class T_ORTHOPOL>
  void HCurlHighOrderTet<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 6; // Nedelec 
    if(augmented==1) ndof +=4; 

    int i;
    for (i = 0; i < 6; i++)
      if(order_edge[i] > 0)
	ndof += usegrad_edge[i]*order_edge[i];
    
    for(i=0; i<4; i++)
      if (order_face[i][0] > 1) 
	ndof += ((usegrad_face[i]+1)*order_face[i][0]+2)*(order_face[i][0]-1)/2 ;
    
    if(order_inner[0] > 2)
      ndof += ((usegrad_cell + 2) * order_inner[0] + 3) 
	* (order_inner[0]-2) * (order_inner[0]-1) / 6; 

    order = 0; // max(order_edges,order_face,order_inner);  
    for (i = 0; i < 6; i++)
      if (order_edge[i] > order)  order = order_edge[i];
    for(i=0; i<4; i++) 
      if (order_face[i][0] > order) 
	order = order_face[i][0]; 
    if (order_inner[0] > order) // order_edge_normal = order_inner; 
      order = order_inner[0];
    
    if(order==0) order++; // for integration order .. 
  }

  template <class T_ORTHOPOL>
  void HCurlHighOrderTet<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    if ( discontinuous )
      {
	int ni = ndof;
	idofs.SetSize (ndof);
	for ( int i = 0; i < ni; i++ )
	  idofs[i] = i;
	return;
      }

    idofs.SetSize (0);
    /*
    if (ndof == 12)
      {
	for (int i = 6; i < 12; i++)
	  idofs.Append (i);
      }

    if (ndof == 30)
      {
	for (int i = 6; i < 18; i++)
	  idofs.Append (i);
	for (int i = 0; i < 4; i++)
	  idofs.Append (18+3*i);
      }

    static int cnt = 0;
    cnt++;
    if (cnt < 10)
      cout << "special hack for Hcurlhighordertet, idofs" << endl;
    return;
    */

    int i, base = 6; // Nedelec 
    if(augmented == 1 ) base +=4; 

    for (i = 0; i < 6; i++)
      base += usegrad_edge[i]*order_edge[i];
    for(i=0; i<4; i++)
      if (order_face[i][0] > 1) 
	base += ((usegrad_face[i]+1)*order_face[i][0]+2)*(order_face[i][0]-1)/2 ;
    
    if(order_inner[0] > 2)
      {
	int ni = ((usegrad_cell + 2) * order_inner[0] + 3) 
	  * (order_inner[0]-2) * (order_inner[0]-1) / 6; 
	
	for (i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }
  
  template <class T_ORTHOPOL>
  void HCurlHighOrderTet<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
						   FlatMatrixFixWidth<3> shape) const
  {    
    AutoDiff<3> x(ip(0), 0);
    AutoDiff<3> y(ip(1), 1);
    AutoDiff<3> z(ip(2), 2);
    
    AutoDiff<3> lami[4] = { x, y, z, 1-x-y-z };

    int i, j, ii, k,l;
    
    shape = 0.0; 
       
    const EDGE * edges = ElementTopology::GetEdges (ET_TET);
    const FACE * faces = ElementTopology::GetFaces (ET_TET); 
    
    ArrayMem<AutoDiff<3>,20> polxi(order+1), poleta(order+1), polzeta(order+1);
    // ArrayMem<double,30> mem(2*order*sizeof(double)); 
    ArrayMem<AutoDiff<3>,10> adpol1(order+2),adpol2(order+2),adpol3(order+2); 
    ii = 6; 

    if(augmented==1)
      for(int i=0;i<4;i++,ii++)
	{
	  AutoDiff<3> v= T_VERTEXSHAPES::Calc (order_vertex[i], lami[i]);
	  for(l=0;l<3;l++) shape(ii,l) = v.DValue(l); 
	}

    for (i = 0; i < 6; i++)
      { 
	int p = order_edge[i]; 
	int es = edges[i][0]; 
	int ee = edges[i][1]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	//Nedelec low order edge shape function 
	for(l=0; l<3; l++) 
	  shape(i,l) = lami[ee].DValue(l)*lami[es].Value()  
	    - lami[ee].Value()*lami[es].DValue(l) ;
	
	//HO-Edge shape functions (Gradient Fields) 	
	if(p>0 && usegrad_edge[i]) 
	  {
	    AutoDiff<3> xi = lami[ee]-lami[es]; 
	    AutoDiff<3> eta = 1-lami[ee]-lami[es]; 
	    T_ORTHOPOL::CalcTrigExt(p+1,xi,eta, adpol1); 

	    for(j=0; j< p;j++,ii++) 	      
	      for(l=0; l<3; l++) 
		shape(ii,l) = adpol1[j].DValue(l); 
	  }
      }

    // face shape functions
    for(i=0; i<4; i++) 
      {
	int fav[3] = { faces[i][0], faces[i][1], faces[i][2] };

	//Sort vertices  first edge op minimal vertex 
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	  	

	int vop = 6 - fav[0] - fav[1] - fav[2];  	
	int p = order_face[i][0];

	if (p >= 2)
	  {
	    AutoDiff<3> xi = lami[fav[2]]-lami[fav[1]];
	    AutoDiff<3> eta = lami[fav[0]]; // lo 
	    AutoDiff<3> zeta = lami[vop];   // lz 

	    T_FACESHAPES::CalcSplitted (p+1, xi, eta, zeta, adpol1, adpol2); 

	    // gradients 
	    if (usegrad_face[i])
	      for (j = 0; j <= p-2; j++)
		for (k = 0; k <= p-2-j; k++, ii++)
		  for (l = 0; l < 3; l++)
		    shape(ii, l) =  adpol1[j].DValue(l)*adpol2[k].Value() 
		      + adpol1[j].Value()*adpol2[k].DValue(l); 

	    // other combination
	    for (j = 0; j <= p-2; j++)
	      for (k = 0; k <= p-2-j; k++, ii++)
		  for (l = 0; l < 3; l++)
		    shape(ii, l) =  adpol1[j].DValue(l)*adpol2[k].Value() 
		      - adpol1[j].Value()*adpol2[k].DValue(l); 

	    double ned[3]; 
	    for(l=0;l<3;l++) 
	      ned[l] = 	lami[fav[1]].Value()*lami[fav[2]].DValue(l) 
		- lami[fav[2]].Value()*lami[fav[1]].DValue(l);  
	
	    // rec_pol2 * Nedelec0
	    for (j = 0; j <= p-2; j++, ii++)
	      for (l = 0; l < 3; l++)
		shape(ii,l) = adpol2[j].Value()*ned[l]; 
	  }
      }
    
    int p = order_inner[0]; 
    
    T_INNERSHAPES::CalcSplitted(p+1, x-(1-x-y-z), y, z,adpol1, adpol2, adpol3 );
    
    //gradient fields 
    if(usegrad_cell)
      for (i = 0; i <= p-3; i++)
	for (j = 0; j <= p-3-i; j++)
	  for (k = 0; k <= p-3-i-j; k++, ii++)
	    for(l=0; l<3; l++) 
	      shape(ii,l) = adpol1[i].DValue(l) * adpol2[j].Value()* adpol3[k].Value() 
		+  adpol1[i].Value() * adpol2[j].DValue(l)* adpol3[k].Value() 
		+  adpol1[i].Value() * adpol2[j].Value()* adpol3[k].DValue(l) ;

    // other combinations
    for (i = 0; i <= p-3; i++)
      for (j = 0; j <= p-3-i; j++)
	for (k = 0; k <= p-3-i-j; k++, ii++, ii++)
	  for(l=0; l<3; l++) 
	    {
	      shape(ii,l) = adpol1[i].DValue(l) * adpol2[j].Value()* adpol3[k].Value() 
		- adpol1[i].Value() * adpol2[j].DValue(l)* adpol3[k].Value() 
		+  adpol1[i].Value() * adpol2[j].Value()* adpol3[k].DValue(l) ; 
	      shape(ii+1,l) = adpol1[i].DValue(l) * adpol2[j].Value()* adpol3[k].Value() 
		+  adpol1[i].Value() * adpol2[j].DValue(l)* adpol3[k].Value() 
		-  adpol1[i].Value() * adpol2[j].Value()* adpol3[k].DValue(l) ; 
	    } 
       
     // Nedelec-types  l4.DValue()*x.Value() - l4.Value()*x.DValue(); 
     
    AutoDiff<3> l4 = 1-x-y-z; 
    double ned[3] = {- x.Value() - l4.Value(), - x.Value(), - x.Value()};
        
     for (j= 0; j <= p-3; j++)
       for (k = 0; k <= p-3-j; k++, ii++)
	 for(l=0; l<3; l++) 
	   shape(ii,l) = adpol2[j].Value()*adpol3[k].Value()*ned[l]; 
  }



  #ifndef NUMCURLTET 
  template <class T_ORTHOPOL>
  void HCurlHighOrderTet<T_ORTHOPOL> :: CalcCurlShape (const IntegrationPoint & ip, 
					 	       FlatMatrixFixWidth<3> shape) const
  {  
    // MatrixFixWidth<3> shape2(shape.Height()); 
    // HCurlHighOrderFiniteElement<3> :: CalcCurlShape(ip,shape2); 
   
  
    AutoDiff<3> x (ip(0), 0); 
    AutoDiff<3> y (ip(1), 1); 
    AutoDiff<3> z (ip(2), 2); 
    
    AutoDiff<3> lami[4] = { x, y, z, 1.0-x-y-z };
    
    const int c[3][2] = {{1,2},{2,0},{0,1}}; 
      
    shape = 0.0; 
    
    const EDGE * edges = ElementTopology::GetEdges (ET_TET);
    const FACE * faces = ElementTopology::GetFaces (ET_TET); 
  
    int ii = 6; 
    if(augmented == 1) ii+=4; 

    for (int i = 0; i < 6; i++)
      {
	int es = edges[i][0]; 
	int ee = edges[i][1]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	//Nedelec low order edge shape function  
        AutoDiff<3> hv = Cross(lami[es],lami[ee]);
	for(int l=0; l<3; l++) 
          shape(i,l) = 2*hv.DValue(l); 
	
	// Tangential Edge Shapes are grad-fields ->  0.0
 	if(usegrad_edge[i]) ii += order_edge[i]; 
      }  
    
    ArrayMem<AutoDiff<3>,20 > adpol1(order+2),adpol2(order+2) ,adpol3(order+2); 

    // face shape functions
    for(int i=0; i<4; i++) // faces 
      {
	int fav[3];
	for(int j=0; j<3; j++) fav[j]=faces[i][j];

	//Sort vertices  first edge op minimal vertex  
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	  	

	int vop = 6 - fav[0] - fav[1] - fav[2];  	
	int p = order_face[i][0];

	if (p >= 2)
	  {
            AutoDiff<3> lam1 = lami[fav[1]];
            AutoDiff<3> lam2 = lami[fav[2]];

	    AutoDiff<3> xi = lam2-lam1;
	    AutoDiff<3> eta = lami[fav[0]]; // in [0,1]
	    AutoDiff<3> zeta = lami[vop];  // in [0,1] // lam_F
	     
	    T_FACESHAPES::CalcSplitted (p+1, xi, eta, zeta, adpol1, adpol2); 
	    // gradients  =>  0 
	    if(usegrad_face[i]) ii+=(p-1)*p/2; 
	    
	    // other combination => 2 * (nabla pol2) x (nabla pol1) 
	    for (int j = 0; j <= p-2; j++)
	      for (int k = 0; k <= p-2-j; k++, ii++)
		{ 
		  AutoDiff<3> hv = Cross(adpol2[k],adpol1[j]); 
		  for(int l=0;l<3;l++)
		    shape(ii,l) = 2*hv.DValue(l); 
		}
	    
            AutoDiff<3> adned;
	    for(int l=0;l<3;l++) 
	      adned.DValue(l) = lam1.Value()*lam2.DValue(l) - lam2.Value()*lam1.DValue(l);  
	    AutoDiff<3> curlned = 2*Cross(lam1,lam2);
 	    
	    for (int j = 0; j <= p-2; j++, ii++)
	      {
                AutoDiff<3> hv = Cross(adpol2[j],adned);
		for(int l=0;l<3;l++)
		  shape(ii,l) = curlned.DValue(l) * adpol2[j].Value()
                    + hv.DValue(l);
	      }
	  } 
      }
	
    // element-based shapes 
    int p = order_inner[0];
    if (p >= 3)
      {
	T_INNERSHAPES::CalcSplitted(p+1, x-(1-x-y-z), y, z,adpol1, adpol2, adpol3 );

	//gradient fields  => 0.0
	if(usegrad_cell) ii+= (p-2)*(p-1)*p/6; 

	// other  combinations => 
	for (int i = 0; i <= p-3; i++)
	  for (int j = 0; j <= p-3-i; j++)
	    for (int k = 0; k <= p-3-i-j; k++, ii+=2)
	      {
		// 2 grad v  x  grad (uw)
		AutoDiff<3> hv = 2*Cross(adpol2[j],adpol1[i]*adpol3[k]);
		for (int l = 0; l < 3; l++)
		  shape(ii,l) = hv.DValue(l); 
		
		// 2 grad w  x  grad (uv)
		hv = 2*Cross(adpol3[k], adpol1[i] * adpol2[j]);
		for (int l = 0; l < 3; l++)
		  shape(ii+1,l) = hv.DValue(l); 
	      }

	// Type 3: 
	
	// oldversion 
	// ned = lami[0] * nabla(lami[3]) - lami[3] * nabla(lami[4]) 
	double ned[3] = {- x.Value() - lami[3].Value(), - x.Value(), - x.Value()};
	// AutoDiff<3> curlned = 2*Cross(x,l4); 
	double curlned[3] = {0,2,-2}; 
          
 	for (int j= 0; j <= p-3; j++)
	  for (int k = 0; k <= p-3-j; k++, ii++)
	    for(int l=0; l<3; l++) 
	      shape(ii,l) = adpol2[j].Value()*adpol3[k].Value()*curlned[l]
		+ adpol3[k].Value()*(adpol2[j].DValue(c[l][0])*ned[c[l][1]]
				     -adpol2[j].DValue(c[l][1])*ned[c[l][0]]) 
		+ adpol2[j].Value()*(adpol3[k].DValue(c[l][0])*ned[c[l][1]] 
				     -adpol3[k].DValue(c[l][1])*ned[c[l][0]]); 

     	// new version (not tested yet) 
	/* 
	   double curlned[3] = {0,2,-2}; 
	   for (j= 0; j <= p-3; j++)
	   for (k = 0; k <= p-3-j; k++, ii++)
	   {
	   AutoDiff<3> pjk = adpol2[j] * adpol3[k]; 
	   AutoDiff<3> hv1 = Cross(lami[3],pjk); 
	   AutoDiff<3> hv2 = Cross(lami[0],pjk); 	   
	   for(l=0; l<3; l++) 
	   shape(ii,l) = pjk.Value()*curlned[l] + hv1.DValue(l)* lami[0].Value() 
	   - hv2.DValue(l) * lami[3].Value(); 
	   }   
	*/
      }
    
    // shape2 -=shape; 
    /* (*testout) << " curlshape code  " << endl << shape << endl ; 
    (*testout) << " curlshape num " << endl << shape2 << endl; 
    shape2 -= shape; 
    (*testout) << " curlshape diff " << endl << shape2 << endl; 
    */
    }		
  #endif







  
  template <class T_ORTHOPOL>
  Vec<3> HCurlHighOrderTet<T_ORTHOPOL> ::
  EvaluateCurlShape (const IntegrationPoint & ip, FlatVector<double> vec, LocalHeap & lh) const
  {
    // return Trans (GetCurlShape(ip, lh)) * vec;

    Vec<3> curl = 0.0;
  
    AutoDiff<3> x (ip(0), 0); 
    AutoDiff<3> y (ip(1), 1); 
    AutoDiff<3> z (ip(2), 2); 
    
    AutoDiff<3> lami[4] = { x, y, z, 1.0-x-y-z };
    
    // shape = 0.0; 
    
    const EDGE * edges = ElementTopology::GetEdges (ET_TET);
    const FACE * faces = ElementTopology::GetFaces (ET_TET); 
  
    int ii = 6; 
    if(augmented == 1) ii+=4; 

    for (int i = 0; i < 6; i++)
      {
	int es = edges[i][0]; 
	int ee = edges[i][1]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	//Nedelec low order edge shape function  
        AutoDiff<3> hv = Cross(lami[es],lami[ee]);
	for(int l=0; l<3; l++) 
          curl(l) += vec(i) * 2 * hv.DValue(l);
	
	// Tangential Edge Shapes are grad-fields ->  0.0
 	if(usegrad_edge[i]) ii += order_edge[i]; 
      }  
    
    ArrayMem<AutoDiff<3>,8> adpol1(order+2),adpol2(order+2) ,adpol3(order+2); 

    // face shape functions
    for(int i=0; i<4; i++) // faces 
      {
	int p = order_face[i][0];
	if (p >= 2)
	  {	
            int fav[3];
            for(int j=0; j<3; j++) 
              fav[j]=faces[i][j];
            
            //Sort vertices  first edge op minimal vertex  
            if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
            if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
            if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	  	
            
            int vop = 6 - fav[0] - fav[1] - fav[2];  	

            AutoDiff<3> lam1 = lami[fav[1]];
            AutoDiff<3> lam2 = lami[fav[2]];
	    AutoDiff<3> eta = lami[fav[0]]; // in [0,1]
	    AutoDiff<3> zeta = lami[vop];  // in [0,1] // lam_F
	     
	    T_FACESHAPES::CalcSplitted (p+1, lam2-lam1, eta, zeta, adpol1, adpol2); 

	    // gradients  =>  0 
	    if(usegrad_face[i]) ii+=(p-1)*p/2; 
	    
	    // other combination => 2 * (nabla pol2) x (nabla pol1) 
	    for (int j = 0; j <= p-2; j++)
              {
                AutoDiff<3> sum = 0.0;
                for (int k = 0; k <= p-2-j; k++, ii++)
                  sum += vec(ii) * adpol2[k];
                
                AutoDiff<3> hv = Cross(sum,adpol1[j]); 
                for(int l=0;l<3;l++)
                  curl(l) += 2*hv.DValue(l);
              }
	    
            AutoDiff<3> adned;
	    for(int l=0;l<3;l++) 
	      adned.DValue(l) = lam1.Value()*lam2.DValue(l) - lam2.Value()*lam1.DValue(l);  
	    AutoDiff<3> curlned = 2*Cross(lam1,lam2);

            AutoDiff<3> sum = 0.0;
            for (int j = 0; j <= p-2; j++, ii++)
              sum += vec(ii) * adpol2[j];

            AutoDiff<3> hv = Cross(sum,adned);
            for(int l=0;l<3;l++)
              curl(l) += curlned.DValue(l) * sum.Value() + hv.DValue(l);
	  } 
      }

    // element-based shapes 
    int p = order_inner[0];
    if (p >= 3)
      {
	T_INNERSHAPES::CalcSplitted(p+1, x-(1-x-y-z), y, z,adpol1, adpol2, adpol3 );

	//gradient fields  => 0.0
	if(usegrad_cell) ii+= (p-2)*(p-1)*p/6; 

	// other  combinations => 
	for (int i = 0; i <= p-3; i++)
	  for (int j = 0; j <= p-3-i; j++)
            {
              AutoDiff<3> sum1 = 0;
              AutoDiff<3> sum2 = 0;

              for (int k = 0; k <= p-3-i-j; k++, ii+=2)
                {
                  sum1 += vec(ii) * adpol3[k];
                  sum2 += vec(ii+1) * adpol3[k];
                }

              AutoDiff<3> hv = 2*Cross(adpol2[j],adpol1[i]*sum1);
              for (int l = 0; l < 3; l++)
                curl(l) += hv.DValue(l);
                  
              // 2 grad w  x  grad (uv)
              hv = 2*Cross(sum2, adpol1[i] * adpol2[j]);
              for (int l = 0; l < 3; l++)
                curl(l) += hv.DValue(l);
            }

	// Type 3: 
	
	// oldversion 
	// ned = lami[0] * nabla(lami[3]) - lami[3] * nabla(lami[4]) 
	// AutoDiff<3> curlned = 2*Cross(x,l4); 

	// double ned[3] = {- x.Value() - lami[3].Value(), - x.Value(), - x.Value()};
        AutoDiff<3> adned;
        adned.DValue(0) = -x.Value()-lami[3].Value();
        adned.DValue(1) = -x.Value();
        adned.DValue(2) = -x.Value();

        // for(int l=0; l<3; l++) adned.DValue(l) = ned[l];
	double curlned[3] = {0,2,-2}; 
          
 	for (int j= 0; j <= p-3; j++)
          {
            AutoDiff<3> sum = 0.0;
            for (int k = 0; k <= p-3-j; k++, ii++)
              sum += vec(ii) * adpol3[k];
            sum *= adpol2[j];
            AutoDiff<3> hv = Cross (sum, adned);

	    for(int l=0; l<3; l++) 
              curl(l) += sum.Value()*curlned[l] + hv.DValue(l);
          }
      }
    return curl;
  }




		         
  
  //------------------------------------------------------------------------
  // HCurlHighOrderPrism 
  //------------------------------------------------------------------------

  template <class T_ORTHOPOL>
  HCurlHighOrderPrism<T_ORTHOPOL> :: HCurlHighOrderPrism (int aorder)
    : HCurlHighOrderFiniteElement<3>(ET_PRISM)
  {
    nv = 6; 
    ned = 9; 
    nf = 5; 
    
    int i;
    order_inner = INT<3> (aorder,aorder,aorder);

    for (i = 0; i < ned; i++)
      order_edge[i] = aorder;
    
    for (i=0; i<nf; i++) 
      order_face[i] = INT<2> (aorder,aorder); 
    

    augmented = 0; 
    for(int i=0;i<6;i++) order_vertex[i] = aorder +1; 
    ComputeNDof();
  }

  template <class T_ORTHOPOL>
  void HCurlHighOrderPrism<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 9; // Nedelec 
    if (augmented == 1) ndof+=6; 
    // if oder_edge <= 1 --> Nedelec
    int i;
    
    for (i = 0; i < 9; i++)
      ndof += usegrad_edge[i]*order_edge[i];
    
    // trig_faces 
    for (i=0; i < 2; i++) 
      if (order_face[i][0] > 1)  
	ndof += ((usegrad_face[i]+1)*order_face[i][0] + 2)*(order_face[i][0]-1)/2; 

    // quad_faces
    for (i=2; i < 5; i++)
      if(order_face[i][0]>=0 && order_face[i][1]>=0)
	  ndof +=  (usegrad_face[i]+1)*order_face[i][0]*order_face[i][1] 
	    + order_face[i][0] + order_face[i][1]; 

    // inner_dof
    if(order_inner[2] > 0 && order_inner[0] > 1)
      ndof += ((usegrad_cell+2)*order_inner[2] + 1) * order_inner[0]*(order_inner[0]-1)/2
	+ (order_inner[0]-1)*order_inner[2]; 
    
    order = 0; // max(order_edges_tang,order_inner);  
    for (i = 0; i < 9; i++)
      if (order_edge[i] > order)
	  order = order_edge[i];
    
    for(i=0; i<5; i++) 
      for(int j=0;j<2;j++)
	if (order_face[i][j] > order) 
	  order = order_face[i][j]; 
    
    if (order_inner[0] > order) 
      order = order_inner[0];
    if (order_inner[2] > order) 
      order = order_inner[2];
   
    order++; 
  }

  template <class T_ORTHOPOL>
  void HCurlHighOrderPrism<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    if ( discontinuous )
      {
	int ni = ndof;
	idofs.SetSize (ndof);
	for ( int i = 0; i < ni; i++ )
	  idofs[i] = i;
	return;
      }

    idofs.SetSize (0);

    int i, base = 9; // Nedelec 
    if(augmented==1) base+=6; 

    for (i = 0; i < 9; i++)
      base += usegrad_edge[i]*order_edge[i];

    // trig_faces 
    for (i=0; i < 2; i++) 
      if (order_face[i][0] > 1)  
	base += ((usegrad_face[i]+1)*order_face[i][0] + 2)* (order_face[i][0]-1)/2; 

    // quad_faces
    for (i=2; i < 5; i++)
      if(order_face[i][0]>=0 &&  order_face[i][1]>=0)
	base +=  (usegrad_face[i]+1)*order_face[i][0]*order_face[i][1]
	  + order_face[i][0] + order_face[i][1];

    if(order_inner[0] > 1 && order_inner[2]> 0)
      {
	int ni = 
	  ((usegrad_cell+2)*order_inner[2] + 1) * order_inner[0]*(order_inner[0]-1)/2
	  + (order_inner[0]-1)*order_inner[2]; 
	for (i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }

 template <class T_ORTHOPOL>
  void HCurlHighOrderPrism<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
		  			     FlatMatrixFixWidth<3> shape) const
  {
    AutoDiff<2> x (ip(0), 0);
    AutoDiff<2> y (ip(1), 1);
    AutoDiff<1> z (ip(2), 0);
    AutoDiff<2> lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    AutoDiff<1> muz[6]  = { 1-z, 1-z, 1-z, z, z, z };
       
    int i, j, k,l, ii;
    
    shape = 0.0; 
       
    const EDGE * edges = ElementTopology::GetEdges (ET_PRISM);
    const FACE * faces = ElementTopology::GetFaces (ET_PRISM); 

    // order+1 auf order+3 erh�ht, test JS
    ArrayMem<AutoDiff<2>,20> adpolxy1(order+3),adpolxy2(order+3); 
    ArrayMem<AutoDiff<1>,20> adpolz(order+3);   

    ii = 9;
    
    if(augmented==1)
      {
	AutoDiff<3> xx(ip(0),0); 
	AutoDiff<3> yy(ip(1),0); 
	AutoDiff<3> zz(ip(2),0); 
	AutoDiff<3> llami[6] =  { xx*(1-zz), yy*(1-zz), (1-xx-yy)*(1-zz), xx*zz, yy*zz, (1-xx-yy)*zz };
	      
	for(int i=9;i<15;i++,ii++)
	  {
	    AutoDiff<3> v = T_VERTEXSHAPES::Calc (order_vertex[i], llami[i]);
	    for(l=0;l<3;l++) shape(ii,l) = v.DValue(l); 
	  }
      }
    
    // horizontal edge shapes
    for (i = 0; i < 6; i++)
      {
	int p = order_edge[i]; 
        int es = edges[i][0]; 
	int ee = edges[i][1]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	//Nedelec0
	for(l=0; l<2; l++) 
	  shape(i,l) = muz[ee].Value()*
	    (lami[es].Value()*lami[ee].DValue(l) - lami[ee].Value()*lami[es].DValue(l));
	     
	//high order \nabla (P_edge(x,y) * muz)
	if(usegrad_edge[i])
	  {
	    T_ORTHOPOL::CalcTrigExt(p+1, lami[ee]-lami[es],
				    1-lami[es]-lami[ee],adpolxy1);
	    
	    for(j=0;j<=p-1;j++)
	      {
		for(l=0;l<2;l++)
		  shape(ii,l) = adpolxy1[j].DValue(l)*muz[ee].Value(); 
		shape(ii++,2) = adpolxy1[j].Value()*muz[ee].DValue(0);
	      }
	  }
      }
    
    //Vertical Edge Shapes
    for (i = 6; i < 9; i++)
      {
	int p = order_edge[i]; 
        int es = edges[i][0]; 
	int ee = edges[i][1];
	if (vnums[es] > vnums[ee]) swap (es, ee);

	shape(i,2) = 0.5*(muz[ee].DValue(0)-muz[es].DValue(0))*lami[ee].Value(); 
	
	//high order edges:  \nabla (T_ORTHOPOL^{p+1}(2z-1) * lami(x,y))
	if(usegrad_edge[i])
	  {
	    T_ORTHOPOL::Calc (p+1, muz[ee]-muz[es], adpolz);
	    
	    for (j = 0; j < p; j++)
	      {
		for (l = 0; l < 2; l++)
		  shape(ii,l) = adpolz[j].Value() * lami[ee].DValue(l);
		shape(ii++,2) = adpolz[j].DValue(0) * lami[ee].Value();
	      }
	  }
      }

    // trig face shapes
    for (i = 0; i < 2; i++)
      {
	int p = order_face[i][0];
	if (p < 2) continue;

	int fav[3] = {faces[i][0], faces[i][1], faces[i][2]};

	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	  	
	
	AutoDiff<2> xi = lami[fav[2]]-lami[fav[1]];
	AutoDiff<2> eta = lami[fav[0]]; // 1-lami[f2]-lami[f1];
	
	T_TRIGFACESHAPES::CalcSplitted(p+1,xi,eta,adpolxy1,adpolxy2); 

	if(usegrad_face[i])
	  // gradient-fields =>  \nabla( adpolxy1*adpolxy2*muz )
	  for (j = 0; j <= p-2; j++)
	    for (k = 0; k <= p-2-j; k++,ii++)
	      {
		for (l = 0; l < 2; l++)
		  shape(ii, l) = (adpolxy1[j].DValue(l)*adpolxy2[k].Value() 
				  + adpolxy1[j].Value()*adpolxy2[k].DValue(l))
		    *muz[fav[2]].Value();
		
		shape(ii,2) = muz[fav[2]].DValue(0)*(adpolxy1[j].Value()
						     * adpolxy2[k].Value());
	      }

	
	// rotations of grad-fields => grad(uj)*vk*w -  uj*grad(vk)*w 

	for (j = 0; j <= p-2; j++)
	  for (k = 0; k <= p-2-j; k++)
	    {
	      for (l = 0; l < 2; l++)
		shape(ii, l) = (adpolxy1[j].DValue(l)*adpolxy2[k].Value() 
				- adpolxy1[j].Value()*adpolxy2[k].DValue(l))
		  *muz[fav[2]].Value();
	      
	    
	       shape(ii++,2) = muz[fav[2]].DValue(0)*(adpolxy1[j].Value()
	      		     *adpolxy2[k].Value());
	     
	    }
	
	double ned[2]= {lami[fav[1]].Value()*lami[fav[2]].DValue(0) 
			- lami[fav[2]].Value()*lami[fav[1]].DValue(0),
			lami[fav[1]].Value()*lami[fav[2]].DValue(1)
			-lami[fav[2]].Value()*lami[fav[1]].DValue(1)}; 
	
	//  Ned0*adpolxy2[j]*muz 
	for (j = 0; j <= p-2; j++,ii++)
	  for (l = 0; l < 2; l++)
	    shape(ii,l) = ned[l] * adpolxy2[j].Value()*muz[fav[2]].Value(); 
      }
    
    
    // quad faces
    for (i = 2; i < 5; i++)
      {
	INT<2> p = order_face[i];
	 
	int fmax = 0;
	for (j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;

	int fz = 3-fmax; 
	int ftrig = fmax^1; 
	AutoDiff<2> xi = lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
	AutoDiff<2> eta = 1-lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
	AutoDiff<1> zeta = muz[faces[i][fmax]]-muz[faces[i][fz]]; 
	
	int pp = int(max2(p[0],p[1]))+1;
	T_ORTHOPOL::CalcTrigExt(pp,xi,eta,adpolxy1); 
	T_ORTHOPOL::Calc(pp,zeta,adpolz); 

	if(usegrad_face[i])
	  {
	    // Gradientfields nabla(polxy*polz) 
	    if (vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) 
	      for (k = 0; k <= p[0]-1; k++)
		for (j = 0; j <= p[1]-1; j++)
		  { 
		    for (l=0; l<2;l++) 
		      shape(ii,l) = adpolxy1[k].DValue(l) * adpolz[j].Value();
		    shape(ii++,2) = adpolxy1[k].Value()*adpolz[j].DValue(0); 
		  } 
	    else
	      for (j = 0; j <= p[0]-1; j++)
		for (k = 0; k <= p[1]-1; k++)
		  {
		    for (l=0; l<2;l++) 
		      shape(ii,l) = adpolxy1[k].DValue(l) * adpolz[j].Value();
		    shape(ii++,2) = adpolxy1[k].Value()*adpolz[j].DValue(0); 
		  }
	  }
	  
	// Rotations of GradFields => nabla(polxy)*polz - polxy*nabla(polz)
	if (vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) 
	  for (k = 0; k <= p[0]-1; k++)
	    for (j = 0; j <= p[1]-1; j++)
	      { 
		for (l=0; l<2;l++) 
		      shape(ii,l) = adpolxy1[k].DValue(l) * adpolz[j].Value();
		shape(ii++,2) = -adpolxy1[k].Value()*adpolz[j].DValue(0); 
	      }
	else
	  for (j = 0; j <= p[0]-1; j++)
	    for (k = 0; k <= p[1]-1; k++)
	      {
		for (l=0; l<2;l++) 
		  shape(ii,l) = -adpolxy1[k].DValue(l) * adpolz[j].Value();
		shape(ii++,2) = adpolxy1[k].Value()*adpolz[j].DValue(0); 
	      }
	
	// Type 3 
	// (ned0_trig)*polz, (ned0_quad)* polxy 
	double ned0trig[2] = { 
	  lami[faces[i][fmax]].DValue(0)*lami[faces[i][ftrig]].Value() 
	  -lami[faces[i][fmax]].Value()*lami[faces[i][ftrig]].DValue(0),
	  lami[faces[i][fmax]].DValue(1)*lami[faces[i][ftrig]].Value() 
	  -lami[faces[i][fmax]].Value()*lami[faces[i][ftrig]].DValue(1)};
	  
	if(vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) // p = (p_trig,p_z) 
	  {
	    for(j=0;j<=p[0]-1;j++,ii++) 
	      shape(ii,2) = 0.5*zeta.DValue(0)*adpolxy1[j].Value();  
	    for(j=0;j<=p[1]-1;j++,ii++) 
	      for(l=0; l<2; l++)
		shape(ii,l) =  ned0trig[l]*adpolz[j].Value(); 
	  }
	else 
	  {
	    for(j=0;j<=p[0]-1;j++,ii++) 
	      for(l=0; l<2; l++)
		shape(ii,l) =  ned0trig[l]*adpolz[j].Value(); 
	    for(j=0;j<=p[1]-1;j++,ii++) 
	      shape(ii,2) = 0.5*zeta.DValue(0)*adpolxy1[j].Value();  
	  }
	
      }
    

    if(order_inner[0] > 1&& order_inner[2]>0) 
      {
	T_TRIGFACESHAPES::CalcSplitted(order_inner[0]+1,x-y,1-x-y,adpolxy1,adpolxy2);
	T_ORTHOPOL::Calc(order_inner[2]+1,2*z-1,adpolz); 
	
	// gradientfields
	if(usegrad_cell)
	  for(i=0;i<=order_inner[0]-2;i++)
	    for(j=0;j<=order_inner[0]-2-i;j++)
	      for(k=0;k<=order_inner[2]-1;k++)
		{
		  for(l=0;l<2;l++)
		    shape(ii,l) = (adpolxy1[i].DValue(l)*adpolxy2[j].Value() +
				   adpolxy1[i].Value()*adpolxy2[j].DValue(l))
		      *adpolz[k].Value(); 
		  shape(ii++,2) = adpolxy1[i].Value()*adpolxy2[j].Value()*
		    adpolz[k].DValue(0); 
		}
	 

	// Rotations of gradientfields
	for(i=0;i<=order_inner[0]-2;i++)
	  for(j=0;j<=order_inner[0]-2-i;j++)
	    for(k=0;k<=order_inner[2]-1;k++)
	      {
		for(l=0;l<2;l++)
		  shape(ii,l) = (adpolxy1[i].DValue(l)*adpolxy2[j].Value() -
				 adpolxy1[i].Value()*adpolxy2[j].DValue(l))
		    *adpolz[k].Value(); 
		shape(ii++,2) = adpolxy1[i].Value()*adpolxy2[j].Value()*adpolz[k].DValue(0); 
		
		for(l=0;l<2;l++)
		  shape(ii,l) = (adpolxy1[i].DValue(l)*adpolxy2[j].Value() +
				  adpolxy1[i].Value()*adpolxy2[j].DValue(l))
		    *adpolz[k].Value(); 
		shape(ii++,2) = -adpolxy1[i].Value()*adpolxy2[j].Value()
		*adpolz[k].DValue(0); 
	      }

	// Type 3 
	// ned0(trig) * polxy2[j]*polz 
	// z.DValue(0) * polxy1[i] * polxy2[j] 
	double ned_trig[2] = {y.Value(),-x.Value()};  
	for(j=0;j<=order_inner[0]-2;j++) 
	  for (k=0;k<=order_inner[2]-1;k++,ii++) 
	      for(l=0;l<2;l++) 
		shape(ii,l) = ned_trig[l]*adpolxy2[j].Value()*adpolz[k].Value(); 
	    
    	for(i=0;i<=order_inner[0]-2;i++) 
	  for(j=0;j<=order_inner[0]-2-i;j++) 
	      shape(ii++,2) = adpolxy1[i].Value()*adpolxy2[j].Value(); 
	
      }
    
    //    (*testout) << "shape = " << shape << endl;
  }

#ifndef NUMCURLPRISM
  template <class T_ORTHOPOL>
  void HCurlHighOrderPrism<T_ORTHOPOL> :: CalcCurlShape (const IntegrationPoint & ip, 
							 FlatMatrixFixWidth<3> shape) const


  {
    /*
    MatrixFixWidth<3> shape2(shape.Height()); 
    HCurlHighOrderFiniteElement<3> :: CalcCurlShape(ip,shape2); 
    */
	
    AutoDiff<3> x (ip(0), 0);
    AutoDiff<3> y (ip(1), 1);
    AutoDiff<3> z (ip(2), 2);
    AutoDiff<3> lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    AutoDiff<3> muz[6]  = { 1-z, 1-z, 1-z, z, z, z };
    
    int i, j, k,l, ii;
    
    shape = 0.0; 
        
    const EDGE * edges = ElementTopology::GetEdges (ET_PRISM);
    const FACE * faces = ElementTopology::GetFaces (ET_PRISM); 

    // order+1 auf order+3 erh�ht, test JS
    ArrayMem<AutoDiff<3>,20> adpolxy1(order+3),adpolxy2(order+3); 
    ArrayMem<AutoDiff<3>,20> adpolz(order+3);   
    
    ii = 9;
    if(augmented==1) ii+=6; 
    // horizontal edge shapes
    for (i = 0; i < 6; i++)
      {
	int p = order_edge[i]; 
        int es = edges[i][0]; 
	int ee = edges[i][1]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	AutoDiff<3> hd = Cross(muz[ee]*lami[es],lami[ee])-Cross(muz[ee]*lami[ee],
							    lami[es]); 
	
	//Nedelec0
	for(l=0; l<3; l++) 
	  shape(i,l) = hd.DValue(l); 
	
	//high order edges = gradfields --> curl = 0.0 
	if(usegrad_edge[i])
	  ii+=p; 
      }
     
    //Vertical Edge Shapes
    for (i = 6; i < 9; i++)
      {
	int p = order_edge[i]; 
        int es = edges[i][0]; 
	int ee = edges[i][1];
	if (vnums[es] > vnums[ee]) swap (es, ee);

	AutoDiff<3> hd = 0.5*Cross(lami[ee],muz[ee]-muz[es]); 
		
	for(l=0;l<3;l++)
	  shape(i,l) = hd.DValue(l); 
		
	//high order edges:  \nabla (T_ORTHOPOL^{p+1}(2z-1) * lami(x,y))
	if(usegrad_edge[i]) ii+=p; 
      }
    
    // trig face shapes
    for (i = 0; i < 2; i++)
      {
	int p = order_face[i][0];
	if (p < 2) continue;

	int fav[3] = {faces[i][0], faces[i][1], faces[i][2]};

	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	  	
	
	AutoDiff<3> xi = lami[fav[2]]-lami[fav[1]];
	AutoDiff<3> eta = lami[fav[0]]; // 1-lami[f2]-lami[f1];
	
	T_TRIGFACESHAPES::CalcSplitted(p+1,xi,eta,adpolxy1,adpolxy2); 

	// gradient-fields =>  curl = 0.0
	if(usegrad_face[i])
	  ii+=(p-1)*p/2; 


	// rotations of grad-fields => grad(p1)*p2 -  p1*grad(p2)
	// --> 2 * grad(p2) x grad(p1)
	for (j = 0; j <= p-2; j++)
	  for (k = 0; k <= p-2-j; k++)
	    {
	      AutoDiff<3> hd = 2*Cross(adpolxy2[k],adpolxy1[j]*muz[fav[2]]); 
	      for (l = 0; l < 3; l++)
		shape(ii,l) = hd.DValue(l); 
	      ii++;
	    }

	//  curl(Ned0*adpolxy2[j]*muz) 
	for (j = 0; j <= p-2; j++,ii++)
	  {
	    AutoDiff<3> hd = 
	      Cross(lami[fav[1]]*adpolxy2[j]*muz[fav[2]],lami[fav[2]])
	      - Cross(lami[fav[2]]*adpolxy2[j]*muz[fav[2]],lami[fav[1]]); 
	    
	    for (l = 0; l < 3; l++)
	      shape(ii,l) = hd.DValue(l); 
	  }
      }
    
    // quad faces
    for (i = 2; i < 5; i++)
      {
	INT<2> p = order_face[i]; 

	int fmax = 0;
	for (j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;

	int fz = 3-fmax; 
	int ftrig = fmax^1; 
	AutoDiff<3> xi = lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
	AutoDiff<3> eta = 1-lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
	AutoDiff<3> zeta = muz[faces[i][fmax]]-muz[faces[i][fz]]; 
	
	int pp = int(max2(p[0],p[1]))+1; 
	T_ORTHOPOL::CalcTrigExt(pp,xi,eta,adpolxy1); 
	T_ORTHOPOL::Calc(pp,zeta,adpolz); 
	
	// Gradientfields nabla(polxy*polz) 
	if(usegrad_face[i])
	  ii+=p[0]*p[1]; 
	  
	// Rotations of GradFields => nabla(polxy)*polz - polxy*nabla(polz)
	// --> 2* grad(polz) * grad(polxy1) 

	double sign=-1; 
	if (vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) //SZ org: <
	  sign = 1; 
	 
	for (k = 0; k <= p[0]-1; k++)
	  for (j = 0; j <= p[1]-1; j++,ii++)
	    {
	      AutoDiff<3> hd = sign*2*Cross(adpolz[j],adpolxy1[k]);
	      for (l=0; l<3;l++) 
		shape(ii,l) = hd.DValue(l); 
	    }

	if(vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) 
	  {
	    for(j=0;j<=p[0]-1;j++, ii++)
	      {  
		AutoDiff<3> hd = 0.5*Cross(adpolxy1[j],zeta);
		for(l=0; l<3; l++)
		  shape(ii,l) = hd.DValue(l);  
	      }
	    for(j=0;j<=p[1]-1;j++, ii++) 
	      {
		AutoDiff<3> hd = 
		  Cross(lami[faces[i][ftrig]]*adpolz[j],lami[faces[i][fmax]])
		  -Cross(lami[faces[i][fmax]]*adpolz[j],lami[faces[i][ftrig]]); 
		for(l=0; l<3; l++)
		  shape(ii,l) = hd.DValue(l);
	      } 
	  } 
	else 
	  {
	    for(j=0;j<=p[0]-1;j++, ii++)
	      {  
		AutoDiff<3> hd = 
		  Cross(lami[faces[i][ftrig]]*adpolz[j],lami[faces[i][fmax]])
		  -Cross(lami[faces[i][fmax]]*adpolz[j],lami[faces[i][ftrig]]); 
		for(l=0; l<3; l++)
		  shape(ii,l) = hd.DValue(l);
	      }
	    for(j=0;j<=p[1]-1;j++, ii++) 
	      {
		AutoDiff<3> hd = 0.5*Cross(adpolxy1[j],zeta);
		for(l=0; l<3; l++)
		  shape(ii,l) = hd.DValue(l);  
	      } 
	  } 
      }
 
   
    if(order_inner[0]>1 && order_inner[2]>0) 
      {
	T_TRIGFACESHAPES::CalcSplitted(order_inner[0]+1,x-y,1-x-y,adpolxy1,adpolxy2);
	T_ORTHOPOL::Calc(order_inner[2]+1,2*z-1,adpolz); 
	
	// gradientfields
	if(usegrad_cell)
	  ii+=order_inner[2]*(order_inner[0]-1)*order_inner[0]/2; 

	// Rotations of gradientfields
	for(i=0;i<=order_inner[0]-2;i++)
	  for(j=0;j<=order_inner[0]-2-i;j++)
	    for(k=0;k<=order_inner[2]-1;k++,ii++)
	      {
		AutoDiff<3> hd = 2*Cross(adpolxy2[j],adpolxy1[i]*adpolz[k]); 
		for(l=0;l<3;l++)
		  shape(ii,l) = hd.DValue(l); 
		ii++; 
		hd = 2*Cross(adpolz[k],adpolxy1[i]*adpolxy2[j]); 
		for(l=0;l<3;l++)
		  shape(ii,l) = hd.DValue(l); 
	      }
	
	// Type 3
	// ned0(trig) * polxy2[j]*polz 
	// z.DValue(0) * polxy1[i] * polxy2[j] 
	double ned0trig[3] = {y.Value(),-x.Value(),0.};  
	double curlned[3] = {0.,0.,-2.};

	for(j=0;j<=order_inner[0]-2;j++) 
	  for (k=0;k<=order_inner[2]-1;k++,ii++)
	    {
	      AutoDiff<3> hd = 
		Cross(adpolz[k]*adpolxy2[j]*y, x) 
		- Cross(adpolz[k]*adpolxy2[j]*x, y);  
		
	      for(l=0;l<3;l++) 
		shape(ii,l) = hd.DValue(l); 
	    }
    	for(i=0;i<=order_inner[0]-2;i++) 
	  for(j=0;j<=order_inner[0]-2-i;j++,ii++) 
	    {
	      AutoDiff<3> hd =	Cross(adpolxy1[i]*adpolxy2[j],z); 
	      
	      for(l=0;l<3;l++)
		shape(ii,l) = hd.DValue(l); 
	    } 
      }

       /*
       // shape2 -=shape; 
       (*testout) << " curlshape code  " << endl << shape << endl ; 
    (*testout) << " curlshape num " << endl << shape2 << endl; 
    shape2 -= shape; 
    (*testout) << " curlshape diff " << endl << shape2 << endl; 
       */
      } 
#endif 

  //------------------------------------------------------------------------
  // HCurlHighOrderHex
  //------------------------------------------------------------------------
  template <class T_ORTHOPOL>
  HCurlHighOrderHex<T_ORTHOPOL> :: HCurlHighOrderHex (int aorder)
    : HCurlHighOrderFiniteElement<3>(ET_HEX)
  {
    int i; 

    order_inner = INT<3> (aorder,aorder,aorder);
    for (i = 0; i < 12; i++)
      order_edge[i] = aorder; 
    for (i=0; i<6; i++) 
      order_face[i] = INT<2> (aorder,aorder); 
    augmented = 0; 
    for(int i=0;i<8;i++) order_vertex[i] = aorder+1; 
    ComputeNDof();
  }
  
  template <class T_ORTHOPOL> 
  void HCurlHighOrderHex<T_ORTHOPOL> :: ComputeNDof()
  { 
    ndof = 12; // Nedelec 
    if(augmented==1) ndof +=8; 
    int i;
    for (i = 0; i < 12; i++)
      if(order_edge[i]>0)
	ndof += usegrad_edge[i]*order_edge[i];
    
    for(i=0; i< 6; i++)   
      if(order_face[i][0]>=0 && order_face[i][1]>=0)
	ndof += (usegrad_face[i]+1)*order_face[i][0]*order_face[i][1]
	  + order_face[i][0] + order_face[i][1];
      
    if(order_inner[0] >= 0 && order_inner[1]>= 0 && order_inner[2]>=0)
      ndof += (usegrad_cell + 2)* order_inner[0] * order_inner[1] * order_inner[2]
	+ order_inner[1]*order_inner[2]  + order_inner[0]*(order_inner[1] + order_inner[2]);  
    
    order = 0; // max(order_edges_tang,order_faces,order_inner);  
    for (i = 0; i < 12; i++)
      if (order_edge[i] > order)
	order = order_edge[i];
    for (i = 0; i < 6; i++)
      for(int j= 0; j< 2; j++ )
	if (order_face[i][j] > order)
	  order = order_face[i][j];
    
    for(i=0;i<3;i++)
      if (order_inner[i] > order) 
	order = order_inner[i];
    
    
    order++; // integration order 
  
  }

  template <class T_ORTHOPOL>
  void HCurlHighOrderHex<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    if ( discontinuous )
      {
	int ni = ndof;
	idofs.SetSize (ndof);
	for ( int i = 0; i < ni; i++ )
	  idofs[i] = i;
	return;
      }

    idofs.SetSize (0);

    int i, base = 12; // Nedelec 
    if(augmented==1) base +=8; 
    for (i = 0; i < 12; i++)
      base += usegrad_edge[i]*order_edge[i];
    for(i=0; i<6; i++)
      if (order_face[i][0] >= 0 && order_face[i][1] >= 0) 
	base += (usegrad_face[i]+1)*order_face[i][0] * order_face[i][1] 
	  + order_face[i][0] + order_face[i][1]; 
    
    if(order_inner[0] > 0 && order_inner[1] > 0  && order_inner[2]>0)
      {
	int ni = (usegrad_cell + 2)* order_inner[0]*order_inner[1]*order_inner[2]
	  + order_inner[1]*order_inner[2] 
	     + order_inner[0]*(order_inner[1] + order_inner[2]);  
	
	for (i = 0; i < ni; i++)
	  idofs.Append (base+i);
      }
  }

  template <class T_ORTHOPOL> 
  void HCurlHighOrderHex<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
				       FlatMatrixFixWidth<3> shape) const
  {
    int i, j, k, l;
    AutoDiff<3> x (ip(0),0);
    AutoDiff<3> y (ip(1),1);
    AutoDiff<3> z (ip(2),2); 
    
    AutoDiff<3> lami[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
			 (1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
    AutoDiff<3> sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
			  (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 
     
    const EDGE * edges = ElementTopology::GetEdges (ET_HEX);

    shape = 0.0;
    int ii = 12; 

    ArrayMem<AutoDiff<3>, 20> pol_xi(order+2),pol_eta(order+2),pol_zeta(order+2);
   
    if(augmented==1)
      for(int i=12;i<20;i++,ii++)
	{
	  AutoDiff<3> v= T_VERTEXSHAPES::Calc (order_vertex[i], lami[i]);
	  for(l=0;l<3;l++) shape(ii,l) = v.DValue(l); 
	}

    // edges
    for (i = 0; i < 12; i++)
      {
	int es = edges[i][0];
	int ee = edges[i][1];
	int p = order_edge[i]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);

	AutoDiff<3> xi = sigma[ee]-sigma[es]; 
	AutoDiff<3> lam_e = lami[ee]+lami[es]; 
	
	// Nedelec0-shapes
	for(l=0; l<3; l++) 
	  shape(i,l) = 0.5*xi.DValue(l)*lam_e.Value(); 
	  	
	if(usegrad_edge[i])
	  {
	    // High Order edges ... Gradient fields 
	    T_ORTHOPOL::Calc (p+1, xi, pol_xi);
	    
	    for (j = 0; j < p; j++, ii++)
	      for(l=0; l<3; l++)
		shape(ii,l) = lam_e.DValue(l)*pol_xi[j].Value() 
		  + lam_e.Value()*pol_xi[j].DValue(l); 
	  }
      }
    
    const FACE * faces = ElementTopology::GetFaces (ET_HEX);
    //Faces 
    for (i = 0; i<6; i++)
      {
	INT<2> p = order_face[i];
	
	AutoDiff<3> lam_f = 0;
	for (j = 0; j < 4; j++)
	  lam_f += lami[faces[i][j]];
	
	int qmax = 0;
	for (j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] > vnums[faces[i][qmax]])
	    qmax = j;
	
	int q1 = (qmax+3)%4; 
	int q2 = (qmax+1)%4; 

	if(vnums[faces[i][q2]] > vnums[faces[i][q1]])
	  swap(q1,q2);  // fmax > f1 > f2

	// horiz = 1 if sigma_{fmax}-sigma_{f1} is horizontal coordinate 
	// for anisotropic meshes (vertical(esp. z) is the longer one) 
	double horiz=1.; 
	if( (qmax & 2) == (q2 & 2)) 
	  horiz = -1.; 
		
	int fmax = faces[i][qmax]; 
	int f1 = faces[i][q1]; 
	int f2 = faces[i][q2]; 
	      
	AutoDiff<3> xi = sigma[fmax]-sigma[f1]; 
	AutoDiff<3> eta = sigma[fmax]-sigma[f2]; 
    
	T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
	T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);
	
	//Gradient fields 
	if(usegrad_face[i])
	  for (k = 0; k < p[0]; k++)
	    for (j= 0; j < p[1]; j++, ii++)
	      for(l=0; l<3; l++)
		shape(ii,l) = lam_f.Value()*(pol_xi[k].DValue(l)*pol_eta[j].Value() + 
					     pol_xi[k].Value()*pol_eta[j].DValue(l))
		  + lam_f.DValue(l)*pol_xi[k].Value()*pol_eta[j].Value();
	
	//Rotation of Gradient fields 
	for (k = 0; k < p[0]; k++)
	  for (j= 0; j < p[1]; j++, ii++)
	    for(l=0; l<3; l++)
	      shape(ii,l) =  
		lam_f.Value()*(pol_xi[k].DValue(l)*pol_eta[j].Value()
			       -pol_xi[k].Value()*pol_eta[j].DValue(l)) 
		+ horiz*lam_f.DValue(l)*pol_xi[k].Value()*pol_eta[j].Value();

	// Missing ones 
	for(j=0;j<p[0];j++, ii++) 
	  for(l=0;l<3;l++)
	    shape(ii,l) = 0.5* eta.DValue(l)*pol_xi[j].Value()*lam_f.Value(); 
	for(j=0;j<p[1];j++, ii++) 
	  for(l=0;l<3;l++)
	    shape(ii,l) = 0.5* xi.DValue(l)*pol_eta[j].Value()*lam_f.Value();

      }
    

    
    // Element-based shapes
    T_ORTHOPOL::Calc(order_inner[0]+1,2*x-1,pol_xi);
    T_ORTHOPOL::Calc(order_inner[1]+1,2*y-1,pol_eta);
    T_ORTHOPOL::Calc(order_inner[2]+1,2*z-1,pol_zeta); 
    
    //Gradient fields
    if(usegrad_cell)
      for (i=0; i<order_inner[0]; i++)
	for(j=0; j<order_inner[1]; j++) 
	  for(k=0; k<order_inner[2]; k++,ii++)
	    {
	      shape(ii,0)= pol_xi[i].DValue(0)*pol_eta[j].Value()*pol_zeta[k].Value();
	      shape(ii,1)= pol_xi[i].Value()*pol_eta[j].DValue(1)*pol_zeta[k].Value();
	      shape(ii,2)= pol_xi[i].Value()*pol_eta[j].Value()*pol_zeta[k].DValue(2); 
	    }
    //Rotations of gradient fields
    for (i=0; i<order_inner[0]; i++)
      for(j=0; j<order_inner[1]; j++) 
	for(k=0; k<order_inner[2]; k++,ii+=2)
	  {
	    shape(ii,0) = pol_xi[i].DValue(0)*pol_eta[j].Value()*pol_zeta[k].Value(); 
	    shape(ii,1)= - pol_xi[i].Value()*pol_eta[j].DValue(1)*pol_zeta[k].Value();
	    shape(ii,2) =  pol_xi[i].Value()*pol_eta[j].Value()*pol_zeta[k].DValue(2);
	    
	    shape(ii+1,0) = pol_xi[i].DValue(0)*pol_eta[j].Value()*pol_zeta[k].Value();
	    shape(ii+1,1) = pol_xi[i].Value()*pol_eta[j].DValue(1)*pol_zeta[k].Value();
	    shape(ii+1,2) = - pol_xi[i].Value()*pol_eta[j].Value()*pol_zeta[k].DValue(2); 
	  } 
    
    for(i=0; i<order_inner[0]; i++) 
      for(j=0; j<order_inner[1]; j++)
	shape(ii++,2) = pol_xi[i].Value()*pol_eta[j].Value();       
   

    for(i=0; i<order_inner[0]; i++) 
      for(k=0; k<order_inner[2]; k++)
	shape(ii++,1) = pol_xi[i].Value()*pol_zeta[k].Value(); 

    for(j=0; j<order_inner[1]; j++)
      for(k=0; k<order_inner[2]; k++)
	shape(ii++,0)=pol_eta[j].Value()*pol_zeta[k].Value();   
       
    return;
  }
 
#ifndef NUMCURLHEX

 template <class T_ORTHOPOL> 
  void HCurlHighOrderHex<T_ORTHOPOL> :: CalcCurlShape (const IntegrationPoint & ip, 
				       FlatMatrixFixWidth<3> shape) const
  {
    // MatrixFixWidth<3> shape2(shape.Height()); 
    // HCurlHighOrderFiniteElement<3> :: CalcCurlShape(ip,shape2); 
    
    int i, j, k, l;
    AutoDiff<3> x (ip(0),0);
    AutoDiff<3> y (ip(1),1);
    AutoDiff<3> z (ip(2),2); 
    
    AutoDiff<3> lami[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
			 (1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
    AutoDiff<3> sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
			  (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 
     
    const EDGE * edges = ElementTopology::GetEdges (ET_HEX);

    shape = 0.0;
     int c[3][2] = {{1,2},{2,0},{0,1}}; 

    int ii = 12; 
    if(augmented==1) ii+=8; 

    ArrayMem<AutoDiff<3>, 20> pol_xi(order+2),pol_eta(order+2),pol_zeta(order+2);
    ArrayMem<AutoDiff<3>, 20> pol_lfxi(order+2),pol_lfeta(order+2);
   
     // edges
    for (i = 0; i < 12; i++)
      {
	int es = edges[i][0];
	int ee = edges[i][1];
	int p = order_edge[i]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);

	AutoDiff<3> xi = sigma[ee]-sigma[es]; 
	AutoDiff<3> lam_e = lami[ee]+lami[es]; 
	
	// Nedelec0-shapes
	AutoDiff<3> hv =  0.5 * Cross(lam_e,xi); 
	for(l=0; l<3; l++) 
	  shape(i,l) = hv.DValue(l); 
	// shape (i,l) = (lam_e.DValue(c[l][0])*xi.DValue(c[l][1]) - 
	// lam_e.DValue(c[l][1])*xi.DValue(c[l][0]));
	    
	  	
	if(usegrad_edge[i])
	    ii+=p; 
      }
    
    const FACE * faces = ElementTopology::GetFaces (ET_HEX);
    //Faces 
    for (i = 0; i<6; i++)
      {
	INT<2> p = order_face[i];
	AutoDiff<3> lam_f = 0;
;	for (j = 0; j < 4; j++)
	  lam_f += lami[faces[i][j]];
	
	int qmax = 0;
	for (j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] > vnums[faces[i][qmax]])
	    qmax = j;
	
	int q1 = (qmax+3)%4; 
	int q2 = (qmax+1)%4; 
       

	if(vnums[faces[i][q2]] > vnums[faces[i][q1]])
	  swap(q1,q2);  // fmax > f1 > f2

	// horiz = 1 if sigma_{fmax}-sigma_{f1} is horizontal coordinate 
	// for anisotropic meshes (vertical(esp. z) is the longer one) 
	double horiz=1.; 
	if( (qmax & 2) == (q2 & 2)) 
	  horiz = -1.; 
		
	int fmax = faces[i][qmax]; 
	int f1 = faces[i][q1]; 
	int f2 = faces[i][q2]; 
	      
	AutoDiff<3> xi = sigma[fmax]-sigma[f1]; 
	AutoDiff<3> eta = sigma[fmax]-sigma[f2]; 
    
	T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
	T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);
	
	//Gradient fields 
	if(usegrad_face[i])
	  ii+=p[0]*p[1]; 
	

	for (j = 0; j < p[0]; j++)
	  pol_lfxi[j] = lam_f*pol_xi[j];
	for (j = 0; j < p[1]; j++)
	  pol_lfeta[j] = lam_f*pol_eta[j];
	
	//Rotation of Gradient fields 
	if (horiz == 1)
	  {	
	    for (k = 0; k < p[0]; k++)
	      for (j= 0; j < p[1]; j++, ii++)
		{
		  AutoDiff<3> hd = 2.0 * Cross (pol_eta[j], pol_lfxi[k]);
		  for(l=0; l<3; l++) shape(ii,l) = hd.DValue(l);
		}
	  }
	else
	  {
	    for (k = 0; k < p[0]; k++)
	      for (j= 0; j < p[1]; j++, ii++)
		{
		  AutoDiff<3> hd = 2.0 * Cross (pol_lfeta[j], pol_xi[k]);
		  for(l=0; l<3; l++) shape(ii,l) = hd.DValue(l);
		}
	  }

	//Missing ones 
	for (j= 0; j < p[0]; j++, ii++)
	  {
	    AutoDiff<3> hd = 0.5 * Cross(pol_lfxi[j],eta);
	    for(l=0; l<3; l++)
	      shape(ii,l) = hd.DValue(l);  
	  }
	for(j=0;j<p[1];j++, ii++)
	  {
	    AutoDiff<3> hd = 0.5 * Cross (pol_lfeta[j], xi);
	    for(l=0; l<3; l++)
	      shape(ii,l) = hd.DValue(l);
	  }
      }
        
    // Element-based shapes
        
    T_ORTHOPOL::Calc(order_inner[0]+1,2*x-1,pol_xi);
    T_ORTHOPOL::Calc(order_inner[1]+1,2*y-1,pol_eta);
    T_ORTHOPOL::Calc(order_inner[2]+1,2*z-1,pol_zeta); 
    
    //Gradient fields
    if(usegrad_cell)
      ii+=order_inner[0]*order_inner[1]*order_inner[2]; 
    
    //Rotations of gradient fields
    for (i=0; i<order_inner[0]; i++)
      for(j=0; j<order_inner[1]; j++) 
	for(k=0; k<order_inner[2]; k++,ii++)
	  
	    {
	      AutoDiff<3> uw = pol_xi[i]*pol_zeta[k]; 
	      AutoDiff<3> hd = 2* Cross(pol_eta[j],uw); 
	     
	      for (l=0; l<3; l++) shape(ii,l)=hd.DValue(l); 

	      AutoDiff<3> uv = pol_xi[i]*pol_eta[j]; 
	      hd = 2* Cross(pol_zeta[k],uv); 
	      ii++;
	      for (l=0; l<3; l++) shape(ii,l)=hd.DValue(l); 
	      
	    } 
    
    //Missing ones
    for(i=0;i<order_inner[0];i++)
      for(j=0;j<order_inner[1];j++)
	{
	  AutoDiff<3> uv = pol_xi[i]*pol_eta[j]; 
	  AutoDiff<3> hd = Cross(uv,z); 
	  for (l=0; l<3; l++) shape(ii,l)=hd.DValue(l);
	  ii++;
	}
    for(i=0;i<order_inner[0];i++)
      for(k=0;k<order_inner[2];k++)
	{
	  AutoDiff<3> uw = pol_xi[i]*pol_zeta[k]; 
	  AutoDiff<3> hd = Cross(uw,y); 
	  for (l=0; l<3; l++) shape(ii,l)=hd.DValue(l); 
	  ii++; 
	}
    for(j=0;j<order_inner[1];j++)
      for(k=0;k<order_inner[2];k++)
	{
	  AutoDiff<3> vw = pol_eta[j]*pol_zeta[k]; 
	  AutoDiff<3> hd = Cross(vw,x); 
	  for (l=0; l<3; l++) shape(ii,l)=hd.DValue(l); 
	  ii++; 
	}
    
    /* (*testout) << " CURL SHAPE HEX " << endl; 
      (*testout) << " curlshape code  " << endl << shape << endl ;  
      (*testout) << " curlshape num " << endl << shape2 << endl; 
      shape2 -= shape;  
      (*testout) << " curlshape diff " << endl << shape2 << endl; 
    */
    return; 
  }
 #endif 
//------------------------------------------------------------------------
// HCurlHighOrderPyr 
//------------------------------------------------------------------------
  
  template <class T_ORTHOPOL>
  HCurlHighOrderPyr<T_ORTHOPOL> :: HCurlHighOrderPyr (int aorder)
    : HCurlHighOrderFiniteElement<3>(ET_PYRAMID)
  {
    for(int m=0;m<3; m++) order_inner[m] = aorder;
  
    for (int i = 0; i < 8; i++)
      order_edge[i] = aorder;
    
    for (int i=0; i< 5; i++) 
      order_face[i] = INT<2> (aorder,aorder); 
    
    order_inner =INT<3>(aorder,aorder,aorder); 
    ComputeNDof();
  }

  template <class T_ORTHOPOL>
  void HCurlHighOrderPyr<T_ORTHOPOL> :: ComputeNDof()
  {
    ndof = 8; // Nedelec 
    for (int i = 0; i < 8; i++)
      ndof += usegrad_edge[i]* order_edge[i];

    // trig_faces 
    for (int i=0; i < 4; i++) 
      if (order_face[i][0] > 1)  
	ndof += (usegrad_face[i]+1)*(order_face[i][0]-1)*order_face[i][0]/2 
	  + (order_face[i][0] - 1); 
    
    // quad_face
    if(order_face[4][0]>=0 && order_face[4][1]>=0)
      ndof +=  (usegrad_face[4]+1)*order_face[4][0]*order_face[4][1]  
		+ order_face[4][0] + order_face[4][1];
    
    int pc = order_inner[0]; //SZ: no problem to do anisotropic, but for the moment 
                            // is it worth getting crazy :-) 
    if(order_inner[0]>1)
      {
	ndof += usegrad_cell*(pc-1)*pc*(2*pc-1)/6 + pc*(2*pc*pc+3*pc-2)/3; 
      }
      // ndof+= usegrad_cell*(order_inner[0]-1)*order_inner[0]*(2*order_inner[0]-1)/6+ order_inner[0]*(2*order_inner[0]*order_inner[0]+3*order_inner[0]-2)/3;
    
    //ndof -= 0.5*order_inner[0]*(order_inner[0]+1); 

    order = 0; // max(order_edges_tang,order_inner);  
    for (int i = 0; i < 8; i++)
      if (order_edge[i] > order)
	order = order_edge[i];
    
    for(int i=0; i<5; i++) 
      if (order_face[i][0] > order) 
	order = order_face[i][0]; 
    
    if (order_inner[0] > order) // order_edge_normal = order_inner; 
      order = order_inner[0];
 
    order++;
  }

  template <class T_ORTHOPOL>
  void HCurlHighOrderPyr<T_ORTHOPOL> ::  
  GetInternalDofs (Array<int> & idofs) const
  {
    if ( discontinuous )
      {
	int ni = ndof;
	idofs.SetSize (ndof);
	for ( int i = 0; i < ni; i++ )
	  idofs[i] = i;
	return;
      }

    idofs.SetSize (0);

    int base = 8; // Nedelec 
    for (int i = 0; i < 8; i++)
      base += usegrad_edge[i]*order_edge[i];
    
    // trig_faces 
    for (int i=0; i < 4; i++) 
      if (order_face[i][0] > 1)  
	base += (usegrad_face[i]+1)*(order_face[i][0]-1)*order_face[i][0]/2 
	  + (order_face[i][0] - 1); 

    // quad_faces
    base +=  ((usegrad_face[4]+1)*order_face[4][0] + 2)*order_face[4][0];
   

    int pc = order_inner[0];
   
    if(pc>1)
      {
        int ni = usegrad_cell*(pc-1)*pc*(2*pc-1)/6 + pc*(2*pc*pc+3*pc-2)/3; 
	
	for(int i=0;i<ni;i++)
	  idofs.Append(base+i);
      }
  }

  template <class T_ORTHOPOL>
  void HCurlHighOrderPyr<T_ORTHOPOL> :: CalcShape (const IntegrationPoint & ip, 
						   FlatMatrixFixWidth<3> shape) const
  {
   
    AutoDiff<3> x (ip(0), 0);
    AutoDiff<3> y (ip(1), 1);
    AutoDiff<3> z (ip(2), 2);

    if(z.Value()==1.) z.Value() -=1.e-8; 


    AutoDiff<3> xt = x/(1-z); 
    AutoDiff<3> yt = y/(1-z); 
    AutoDiff<3> sigma[5] = {(1-xt)+(1-yt)+(1-z),xt+(1-yt)+(1-z), xt + yt + (1-z), 
			    (1-xt)+yt+(1-z),z}; 

    AutoDiff<3> lami[5] = {(1-xt)*(1-yt)*(1-z),xt*(1-yt)*(1-z), xt * yt * (1-z), 
			   (1-xt)*yt*(1-z),z}; 

    AutoDiff<3> sig[5] = {(1-xt)+(1-yt),xt+(1-yt), xt + yt, 
			    (1-xt)+yt,z}; 

    AutoDiff<3> lambda[5] = {(1-xt)*(1-yt),xt*(1-yt), xt * yt, 
			   (1-xt)*yt,z}; 
        
    
    shape = 0.0; 
   
       
    const EDGE * edges = ElementTopology::GetEdges (ET_PYRAMID);
    const FACE * faces = ElementTopology::GetFaces (ET_PYRAMID); 
    
    ArrayMem<AutoDiff<3>, 20> pol_xi(order+2), pol_eta(order+2), pol_zeta(order+2),pol2_zeta(order+2); 
    
    int ii =8; 
 
    // horizontal edges incl. Nedelec 0
    for (int i = 0; i < 4; i++)
      {
	int es = edges[i][0]; 
	int ee = edges[i][1]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	AutoDiff<3> xi  = sigma[ee] - sigma[es];   
	AutoDiff<3> lam_e = lami[ee] + lami[es];
	AutoDiff<3> lam_t = lambda[ee] + lambda[es]; // lam_t = lam_e/(1-z); 
	
	for(int l=0; l<3; l++) 
	  shape(i,l) = (lami[ee].DValue(l) * lami[es].Value() - lami[ee].Value()* lami[es].DValue(l))
	    /(lam_e.Value()+1e-8) * (1-z.Value());   
	    
	if(usegrad_edge[i])
	  {
	    T_ORTHOPOL::CalcTrigExt (order_edge[i]+1, xi*(1-z), z, pol_xi);
	    
	    for(int j=0;j<order_edge[i];j++, ii++)
	      for(int l=0;l<3;l++) 
		shape(ii,l) = pol_xi[j].Value()*lam_t.DValue(l) + 
		  pol_xi[j].DValue(l)*lam_t.Value(); 
	  }
      }
    
    // vertical edges incl. Nedelec 0  
    for(int i=4;i<8;i++)
      {
	int es = edges[i][0];
	int ee = edges[i][1]; 
	
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	for( int l=0;l<3;l++)
	  shape(i,l) =  lami[ee].DValue(l)*lami[es].Value() - lami[es].DValue(l)*lami[ee].Value(); 
	
	if (usegrad_edge[i])
	  {
	    int ne = T_ORTHOPOL::CalcTrigExt (order_edge[i]+1, lami[ee]-lami[es],  
					      1-lami[es]-lami[ee], pol_xi);
	    for(int j=0; j < order_edge[i]; j++,ii++)
	      for(int l=0;l<3;l++)
		shape(ii,l) = pol_xi[j].DValue(l);
	  }
      }
     
  // trig face dofs
  for (int i = 0; i < 4; i++)
      if (order_face[i][0] >= 2)
	{
	  int p = order_face[i][0];
	  AutoDiff<3> lam_face = lambda[faces[i][0]] + lambda[faces[i][1]];  
	  AutoDiff<3> bary[3] = 
	    {(sigma[faces[i][0]]-(1-z)-lam_face)*(1-z), 
	     (sigma[faces[i][1]]-(1-z)-lam_face)*(1-z), z}; 
			     
	  int fav[3] = {0, 1, 2};
	  if(vnums[faces[i][fav[0]]] > vnums[faces[i][fav[1]]]) swap(fav[0],fav[1]); 
	  if(vnums[faces[i][fav[1]]] > vnums[faces[i][fav[2]]]) swap(fav[1],fav[2]);
	  if(vnums[faces[i][fav[0]]] > vnums[faces[i][fav[1]]]) swap(fav[0],fav[1]); 	
	 
	  T_TRIGFACESHAPES::CalcSplitted(p+1, bary[fav[2]]-bary[fav[1]], 
					 bary[fav[0]],pol_xi,pol_eta);
	  
	  for(int j=0;j<=p-2;j++) pol_eta[j] *= lam_face;  
	  
	  // phi = pol_xi * pol_eta * lam_face; 
	  // Type 1: Gradient Functions 
	  if(usegrad_face[i])
	    for(int j=0;j<= p-2; j++)
	      for(int k=0;k<=p-2-j; k++, ii++)
		for(int l=0;l<3;l++)
		  shape(ii,l) = pol_xi[j].DValue(l)*pol_eta[k].Value() 
				  + pol_xi[j].Value()*pol_eta[k].DValue(l); 
	  // Type 2:  
	  for(int j=0;j<= p-2; j++)
	    for(int k=0;k<=p-2-j; k++, ii++)
	      for(int l=0;l<3;l++)
		shape(ii,l) = pol_xi[j].DValue(l)*pol_eta[k].Value() 
		  - pol_xi[j].Value()*pol_eta[k].DValue(l);
	 
	  // Type 3: Nedelec-based ones (Ned_0*v_j)
	  double ned[3]; 
	  for(int l=0;l<3;l++) 
	    ned[l] = - bary[fav[2]].Value()*bary[fav[1]].DValue(l) 
	      + bary[fav[1]].Value()*bary[fav[2]].DValue(l);
	  for(int j=0;j<=p-2;j++,ii++)
	    for(int l=0;l<3;l++)
	      shape(ii,l) = ned[l]*pol_eta[j].Value(); 
	}
 
    // quad face 
    if (order_face[4][0] >= 1)
      {
	int px = order_face[4][0];
	int py = order_face[4][0]; // SZ-Attentione 
	int p = max2(px, py);
	int pp = p+1;

	AutoDiff<3> fac = 1.0;
	for (int k = 1; k <= p; k++)
	  fac *= (1-z);
	
	int fmax = 0;
	for (int l=1; l<4; l++) 
	  if (vnums[l] > vnums[fmax]) fmax = l;  
    // 
	int f1 = (fmax+3)%4;
	int f2 = (fmax+1)%4; 
	if(vnums[f1]>vnums[f2]) swap(f1,f2);  // fmax > f2 > f1 
       
 	AutoDiff<3> xi  = sigma[fmax] - sigma[f2]; 
	AutoDiff<3> eta = sigma[fmax] - sigma[f1];
	
	pol_eta[pp-1] = 0; //!
	pol_xi[pp-1] = 0;  //!

	T_ORTHOPOL::Calc (pp+2, xi, pol_xi);	
	T_ORTHOPOL::Calc (pp+2, eta, pol_eta);
	
	for(int k=0;k<pp;k++) pol_eta[k] = fac*pol_eta[k]; 

	// Type 1: Gradient-fields 
	if (usegrad_face[4])
	  for (int k = 0; k <= px-1; k++) 
	    for (int j = 0; j <= py-1; j++, ii++) 
	      for (int l = 0; l < 3; l++)
		shape(ii,l) = pol_xi[k].DValue(l) * pol_eta[j].Value() 
		  + pol_xi[k].Value() * pol_eta[j].DValue(l); 

	// Type 2: 
	for (int k = 0; k < px; k++) 
	  for (int j = 0; j < py; j++, ii++) 
	    for (int l = 0; l < 3; l++)
	      shape(ii,l)= pol_xi[k].DValue(l) * pol_eta[j].Value() 
		- pol_xi[k].Value() * pol_eta[j].DValue(l);

	// Type 3:
	for (int k=0; k< px; k++,ii++)
	  for(int l=0;l<3;l++)
	    shape(ii,l) = 0.5*eta.DValue(l)*pol_xi[k].Value()*fac.Value();

	for (int k=0; k< py; k++,ii++)
	  for(int l=0;l<3;l++)
	    shape(ii,l) = 0.5*xi.DValue(l)*pol_eta[k].Value(); 
      }
 
#ifdef newinner 
  if (order_inner[0] >= 2)
      {
	int pp = order_inner[0];
	// According H^1 terms: 
	// u_i = L_i+2(2xt-1)
	// v_j = L_j+2(2yt-1) 
	// w_k = z * (1-z)^(k+2)  with 0 <= i,j <= k, 0<= k <= p-2  
	
	T_ORTHOPOL::Calc (pp+3, 2*xt-1, pol_xi);
	T_ORTHOPOL::Calc (pp+3, 2*yt-1, pol_eta);		
 		
	pol_zeta[0] = z*(1-z)*(1-z);
	for (int k=1;k<=pp-2;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];

	// This one includes constant and linear polynomial
	AutoDiff<3> zz = 2*z-1; 
	IntegratedLegendrePolynomial (pp, zz, pol2_zeta); 
	

	// --> Attention: Use Integrated Legendre also for H^1
	if(usegrad_cell)
	  {
	    for(int k=0;k<= pp-2;k++)
	      for(int i=0;i<=k;i++)
		for(int j=0;j<=k;j++,ii++)
		  { 
		    if(i+j+2>=k)
		      {
			for(int l=0;l<3;l++)
			  shape(ii,l) = pol_xi[i].DValue(l)* pol_eta[j].Value()* pol_zeta[k].Value() 
			    + pol_xi[i].Value()* pol_eta[j].DValue(l)* pol_zeta[k].Value() 
			    + pol_xi[i].Value()* pol_eta[j].Value()* pol_zeta[k].DValue(l);
		      }
		    else 
		      {
			AutoDiff<3> pz = pol_zeta[i+j+2].Value()*pol2_zeta[k-i-j-2].Value();  
			for(int l=0;l<3;l++)
			  shape(ii,l) = pol_xi[i].DValue(l)* pol_eta[j].Value()* pz.Value() 
			    + pol_xi[i].Value()* pol_eta[j].DValue(l)* pz.Value() 
			    + pol_xi[i].Value()* pol_eta[j].Value()* pz.DValue(l);
		      }
		  }
	  }
	// Type 2a: l.i. combinations of grad-terms   
	// shape = u_i \nabla(v_j) w_k 
	// shape = u_i v_j \nabla(w_k) 
	for(int k=0;k<= pp-2;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
	      {
		if(i+j+2>=k)
		  {
		    for(int l=0;l<3;l++)
		      shape(ii,l) = pol_eta[j].DValue(l)*pol_xi[i].Value()*pol_zeta[k].Value(); 
		  }
		else 
		  {
		    AutoDiff<3> pz = pol_zeta[i+j+2].Value()*pol2_zeta[k-i-j-2].Value();  
		    for(int l=0;l<3;l++)
		      shape(ii,l) = pol_eta[j].DValue(l)*pol_xi[i].Value()*pz.Value(); 
		  }
		    
		//	shape(ii++,2) = pol_xi[i].Value()*pol_eta[j].Value()* pol_zeta[k].DValue(2);
	      }


	// Type 2b: shape = v_j w_k \nabla (xt) 
	//          shape = u_i w_k \nabla (yt)
	for(int k = 0;k<= pp-2;k++)
	  for(int j=0;j<=k;j++,ii++) 
	    {
	      double vw; 
	      if(j>=k)
		vw = pol_eta[j].Value()*pol_zeta[k].Value();
	      else
		vw = pol_eta[j].Value()*pol_zeta[j].Value()*pol2_zeta[k-j].Value();
	
	      shape(ii,0)=xt.DValue(0)*vw;
	      shape(ii,2)=xt.DValue(2)*vw;
	    }
	
	for(int  k = 0;k<= pp-2;k++)
	  for (int i=0;i<=k;i++,ii++)
	    {
	      double uw; 
	      if(i>=k)
		uw = pol_xi[i].Value()*pol_zeta[k].Value();
	      else
		uw = pol_xi[i].Value()*pol_zeta[i].Value()*pol2_zeta[k-i].Value();
	      
	      shape(ii,1)=yt.DValue(1)* uw; 
	      shape(ii,2)=yt.DValue(2)* uw;
	    }
	

	// 3rd component spans xi^i eta^j zeta^(k-1), i,j <= k
	// pol_zeta starts linear in zeta 
	// pol_xi and pol_eta quadratic in xi resp. eta 

	pol_zeta[0] = (1-z);
	for (int k=1;k<=pp;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];
	
	for(int k=0;k<= pp-1;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
	      {
		double pz; 
		if(i+j+3>=k)
		  pz = pol_zeta[k].Value(); 
		else 
		  pz = pol_zeta[i+j+3].Value()*pol2_zeta[k-i-j-3].Value();
		
		shape(ii,2) = pol_eta[j].Value()*pol_xi[i].Value()*pz; 
	      }
      }
    	
#else 
    if (order_inner[0] >= 2)
      {
	int pp = order_inner[0];
	// According H^1 terms: 
	// u_i = L_i+2(2xt-1)
	// v_j = L_j+2(2yt-1) 
	// w_k = z * (1-z)^(k+2)  with 0 <= i,j <= k, 0<= k <= p-2  
	
	T_ORTHOPOL::Calc (pp+3, 2*xt-1, pol_xi);
	T_ORTHOPOL::Calc (pp+3, 2*yt-1, pol_eta);		
		
	pol_zeta[0] = z*(1-z)*(1-z);
	for (int k=1;k<=pp-2;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];

	// This one includes constant and linear polynomial
	AutoDiff<3> zz = 2*z-1; 
	IntegratedLegendrePolynomial (pp, zz, pol2_zeta); 
	

	// --> Attention: Use Integrated Legendre also for H^1
	if(usegrad_cell)
	  {
	    for(int k=0;k<= pp-2;k++)
	      for(int i=0;i<=k;i++)
		for(int j=0;j<=k;j++,ii++)
		  { 
		    for(int l=0;l<3;l++)
		      shape(ii,l) = pol_xi[i].DValue(l)* pol_eta[j].Value()* pol_zeta[k].Value() 
			+ pol_xi[i].Value()* pol_eta[j].DValue(l)* pol_zeta[k].Value() 
			+ pol_xi[i].Value()* pol_eta[j].Value()* pol_zeta[k].DValue(l);
		    
		  }
	  }
	
      
	
	// Type 2a: l.i. combinations of grad-terms   
	// shape = u_i \nabla(v_j) w_k 
	// shape = u_i v_j \nabla(w_k) 
	for(int k=0;k<= pp-2;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
	      {
		for(int l=0;l<3;l++)
		  shape(ii,l) = pol_eta[j].DValue(l)*pol_xi[i].Value()*pol_zeta[k].Value(); 
	      }
	

	// Type 2b: shape = v_j w_k \nabla (xt) 
	//          shape = u_i w_k \nabla (yt)
	for(int k = 0;k<= pp-2;k++)
	  for(int j=0;j<=k;j++,ii++) 
	    {
	      double vw = pol_eta[j].Value()*pol_zeta[k].Value();
	      shape(ii,0)=xt.DValue(0)*vw;
	      shape(ii,2)=xt.DValue(2)*vw;
	    }
	
	for(int  k = 0;k<= pp-2;k++)
	  for (int i=0;i<=k;i++)
	    {
	      double uw = pol_xi[i].Value()*pol_zeta[k].Value();
	      shape(ii,1)=yt.DValue(1)* uw; 
	      shape(ii++,2)=yt.DValue(2)* uw;
	    }
	
	
	// 3rd component spans xi^i eta^j zeta^(k-1), i,j <= k
	// pol_zeta starts linear in zeta 
	// pol_xi and pol_eta quadratic in xi resp. eta 
	pol_zeta[0] = (1-z);
	for (int k=1;k<=pp;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];
     
	for(int k=0;k<= pp-1;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
	      {
		double pz = pol_zeta[k].Value(); 
		shape(ii,2) = pol_eta[j].Value()*pol_xi[i].Value()*pz; 
	      }
	
      }
#endif 	  
      
          



    /*      *testout << " pyramid ii " <<  ii << endl; 
     *testout << " pyramid ndof " << ndof << endl; 
    
     *testout << " pyramid shapes x = " << x.Value() << "\t" << y.Value() << "\t" << z.Value() << endl;  
     
     *testout << shape << endl; */ 
     
  }

      
  

       
   
#ifndef NUMCURLPYR
  template <class T_ORTHOPOL>
  void HCurlHighOrderPyr<T_ORTHOPOL> :: CalcCurlShape (const IntegrationPoint & ip, 
						   FlatMatrixFixWidth<3> shape) const
  {
 
    //    MatrixFixWidth<3> shape2(shape.Height()); 
    // HCurlHighOrderFiniteElement<3> :: CalcCurlShape(ip,shape2); 
   
     

    AutoDiff<3> x (ip(0), 0);
    AutoDiff<3> y (ip(1), 1);
    AutoDiff<3> z (ip(2), 2);

    int i, j, k,l;
    
    AutoDiff<3> xt = x/(1-z); 
    AutoDiff<3> yt = y/(1-z); 
    AutoDiff<3> sigma[5] = {(1-xt)+(1-yt)+(1-z),xt+(1-yt)+(1-z), xt + yt + (1-z), 
			    (1-xt)+yt+(1-z),z}; 

    AutoDiff<3> lami[5] = {(1-xt)*(1-yt)*(1-z),xt*(1-yt)*(1-z), xt * yt * (1-z), 
			   (1-xt)*yt*(1-z),z}; 

    AutoDiff<3> sig[5] = {(1-xt)+(1-yt),xt+(1-yt), xt + yt, 
			    (1-xt)+yt,z}; 

    AutoDiff<3> lambda[5] = {(1-xt)*(1-yt),xt*(1-yt), xt * yt, 
			   (1-xt)*yt,z}; 

    shape = 0.0; 
    int ii = 8; 
       

    
    const EDGE * edges = ElementTopology::GetEdges (ET_PYRAMID);
    const FACE * faces = ElementTopology::GetFaces (ET_PYRAMID); 
    
    ArrayMem<AutoDiff<3>, 20> pol_xi(order+2), pol_eta(order+2), pol_zeta(order+2), pol2_zeta(order+4); 
        
    // horizontal edges incl. Nedelec 0
    for (i = 0; i < 4; i++)
      {
      	int es = edges[i][0]; 
	int ee = edges[i][1]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	AutoDiff<3> hd1 = Cross(lami[es]/(lambda[ee]+lambda[es]+1.e-8),lami[ee]); 
	AutoDiff<3> hd2 = Cross(lami[es],lami[ee]/(lambda[ee]+lambda[es]+1.e-8));
 
	for(l=0; l<3; l++) 
	  shape(i,l) = hd1.DValue(l) + hd2.DValue(l); 

	if(usegrad_edge[i]) ii += order_edge[i]; 
      }
    // vertical edges incl. Nedelec 0  
    for(i=4;i<8;i++)
      {
	int es = edges[i][0];
	int ee = edges[i][1]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	AutoDiff<3> hv = 2*Cross(lami[es],lami[ee]); 
	for(l=0;l<3;l++)
	  shape(i,l) = hv.DValue(l); 

	if (usegrad_edge[i]) ii+=order_edge[i]; 
      }
  
    // trig face dofs
    for (i = 0; i < 4; i++)
      if (order_face[i][0] >= 2)
	{
	  int p = order_face[i][0];
	  AutoDiff<3> lam_face = lambda[faces[i][0]] + lambda[faces[i][1]];  
	  AutoDiff<3> bary[3] = 
	    {(sigma[faces[i][0]]-(1-z)-lam_face)*(1-z), 
	     (sigma[faces[i][1]]-(1-z)-lam_face)*(1-z), z}; 
			     
	  int fav[3] = {0, 1, 2};
	  if(vnums[faces[i][fav[0]]] > vnums[faces[i][fav[1]]]) swap(fav[0],fav[1]); 
	  if(vnums[faces[i][fav[1]]] > vnums[faces[i][fav[2]]]) swap(fav[1],fav[2]);
	  if(vnums[faces[i][fav[0]]] > vnums[faces[i][fav[1]]]) swap(fav[0],fav[1]); 	
	 
	  T_TRIGFACESHAPES::CalcSplitted(p+1, bary[fav[2]]-bary[fav[1]], 
					 bary[fav[0]],pol_xi,pol_eta);
	  
	  for(j=0;j<=p-2;j++) pol_eta[j] *= lam_face;  
	  
	  // Type 1: Gradient Functions 
	  if(usegrad_face[i]) ii += (p-1)*p/2; 

	  // Type 2:  
	  for(j=0;j<= p-2; j++)
	    for(k=0;k<=p-2-j; k++, ii++) 
	      {
		AutoDiff<3> hv =2* Cross(pol_eta[k],pol_xi[j]);
		  for(l=0;l<3;l++)
		    shape(ii,l) =  hv.DValue(l); 
	      }
	  
	  // Type 3: Nedelec-based ones (Ned_0*v_j)
	  AutoDiff<3> curlned = 2*Cross(bary[fav[1]],bary[fav[2]]);
	  	  
	  for(j=0;j<=p-2;j++,ii++)
	    {
	      AutoDiff<3> hv1 = Cross(pol_eta[j],bary[fav[2]]); 
	      AutoDiff<3> hv2 = Cross(pol_eta[j],bary[fav[1]]);
	      
	      for(l=0;l<3;l++)
		shape(ii,l) = curlned.DValue(l)*pol_eta[j].Value() 
		  + bary[fav[1]].Value()*hv1.DValue(l) - bary[fav[2]].Value()*hv2.DValue(l); 
	    }
 	}
    
    // quad face 
    if (order_face[4][0] >= 1)
      {
	int px = order_face[4][0];
	int py = order_face[4][1];
	int p = max2(px, py);
	int pp = p+1;

	AutoDiff<3> fac = 1.0;
	for (k = 1; k <= p; k++)
	  fac *= (1-z);
	
	int fmax = 0;
	for (l=1; l<4; l++) 
	  if (vnums[l] > vnums[fmax]) fmax = l;  

	int f1 = (fmax+3)%4;
	int f2 = (fmax+1)%4; 
	if(vnums[f1]>vnums[f2]) swap(f1,f2);  // fmax > f2 > f1 
       
 	AutoDiff<3> xi  = sigma[fmax] - sigma[f2]; 
	AutoDiff<3> eta = sigma[fmax] - sigma[f1];

	T_ORTHOPOL::Calc (pp, xi, pol_xi);	
	T_ORTHOPOL::Calc (pp, eta, pol_eta);
	
	for(k=0;k<p;k++) pol_eta[k] = fac*pol_eta[k]; 

	// Type 1: Gradient-fields 
	if (usegrad_face[4]) ii+=px*py; 
	  
	// Type 2: 
	for (k = 0; k < px; k++) 
	  for (j = 0; j < py; j++, ii++) 
	    {
	      AutoDiff<3> hv1 = 2*Cross(pol_eta[j],pol_xi[k]); 
	      for (l = 0; l < 3; l++)
		shape(ii,l)= hv1.DValue(l); 
	    }

	// Type 3:
	for (k=0; k< px; k++,ii++)
	  {
	    AutoDiff<3> hv = 0.5*Cross(pol_xi[k]*fac, eta); 
	    for(l=0;l<3;l++)
	      shape(ii,l) = hv.DValue(l);
	  }
	
	for (k=0; k< py; k++,ii++)
	  {
	    AutoDiff<3> hv = 0.5*Cross(pol_eta[k], xi); 
	    for(l=0;l<3;l++)
	      shape(ii,l) = hv.DValue(l); 
	  }
      }

#ifdef newinner 
  if (order_inner[0] >= 2)
      {
	int pp = order_inner[0];
	// According H^1 terms: 
	// u_i = L_i+2(2xt-1)
	// v_j = L_j+2(2yt-1) 
	// w_k = z * (1-z)^(k+2)  with 0 <= i,j <= k, 0<= k <= p-2  
	
	T_ORTHOPOL::Calc (pp+3, 2*xt-1, pol_xi);
	T_ORTHOPOL::Calc (pp+3, 2*yt-1, pol_eta);		
		
	pol_zeta[0] = z*(1-z)*(1-z);
	for (int k=1;k<=pp-2;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];

	// This one includes constant and linear polynomial
	AutoDiff<3> zz = 2*z-1; 
	IntegratedLegendrePolynomial (pp, zz, pol2_zeta); 
	
	if(usegrad_cell)
	  ii += (pp-1)*pp*(2*pp-1)/6; 
	 	   
	// Type 2a: l.i. combinations of grad-terms   
	// shape = u_i \nabla(v_j) w_k 
	// shape = u_i v_j \nabla(w_k) 
	for(int k=0;k<= pp-2;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
	      { 
		AutoDiff<3> hv1; 
		if(i+j+2>=k)
		  hv1 = pol_xi[i]*pol_zeta[k]; 
		else 
		  hv1 = pol_xi[i]*pol_zeta[i+j+2]*pol2_zeta[k-i-j-2];  
		
		AutoDiff<3> hv = Cross(hv1,pol_eta[j]); 
		for(int l=0;l<3;l++)
		  shape(ii,l) = hv.DValue(l); 
	      }
	
	// Type 2b: shape = v_j w_k \nabla (xt) 
	//          shape = u_i w_k \nabla (yt)
	for(int k = 0;k<= pp-2;k++)
	  for(int j=0;j<=k;j++,ii++) 
	    {
	      AutoDiff<3> vw; 
	      if(j>=k)
		vw = pol_eta[j]*pol_zeta[k];
	      else
		vw = pol_eta[j]*pol_zeta[j]*pol2_zeta[k-j];
	
	      AutoDiff<3> hv = Cross(vw,xt); 
	      for(l=0;l<3;l++)
		shape(ii,l)= hv.DValue(l); 
	    }
	
	for(int  k = 0;k<= pp-2;k++)
	  for (int i=0;i<=k;i++,ii++)
	    {
	      AutoDiff<3> uw; 
	      if(i>=k)
		uw = pol_xi[i]*pol_zeta[k];
	      else
		uw = pol_xi[i]*pol_zeta[i]*pol2_zeta[k-i];
	      
	      AutoDiff<3> hv = Cross(uw,yt); 
	      for(l=0;l<3;l++) shape(ii,l)=hv.DValue(l); 
	    }
	
	// 3rd component spans xi^i eta^j zeta^(k-1), i,j <= k
	// pol_zeta starts linear in zeta 
	// pol_xi and pol_eta quadratic in xi resp. eta 

	pol_zeta[0] = (1-z);
	for (int k=1;k<=pp;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];
	
	for(int k=0;k<= pp-1;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
	      {
		AutoDiff<3> pz; 
		if(i+j+3>=k)
		  pz = pol_zeta[k]*pol_eta[j]*pol_xi[i]; 
		else 
		  pz = pol_zeta[i+j+3]*pol2_zeta[k-i-j-3]*pol_eta[j]*pol_xi[i]; 
		
		AutoDiff<3> hv = Cross(pz,z);
		for(l=0;l<3;l++)
		  shape(ii,l) = hv.DValue(l); 
	      }
      }
#else 
  if (order_inner[0] >= 2)
      {
	int pp = order_inner[0];
	// According H^1 terms: 
	// u_i = L_i+2(2xt-1)
	// v_j = L_j+2(2yt-1) 
	// w_k = z * (1-z)^(k+2)  with 0 <= i,j <= k, 0<= k <= p-2  
	
	T_ORTHOPOL::Calc (pp+3, 2*xt-1, pol_xi);
	T_ORTHOPOL::Calc (pp+3, 2*yt-1, pol_eta);		
		
	pol_zeta[0] = z*(1-z)*(1-z);
	for (int k=1;k<=pp-2;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];

	// This one includes constant and linear polynomial
	AutoDiff<3> zz = 2*z-1; 
	IntegratedLegendrePolynomial (pp, zz, pol2_zeta); 
	
	if(usegrad_cell)
	  ii += (pp-1)*pp*(2*pp-1)/6; 
	 	   
	// Type 2a: l.i. combinations of grad-terms   
	// shape = u_i \nabla(v_j) w_k 
	// shape = u_i v_j \nabla(w_k) 
	for(int k=0;k<= pp-2;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
	      { 
		AutoDiff<3> hv1; 
	
		hv1 = pol_xi[i]*pol_zeta[k]; 
			
		AutoDiff<3> hv = Cross(hv1,pol_eta[j]); 
		for(int l=0;l<3;l++)
		  shape(ii,l) = hv.DValue(l); 
	      }
	
	// Type 2b: shape = v_j w_k \nabla (xt) 
	//          shape = u_i w_k \nabla (yt)
	for(int k = 0;k<= pp-2;k++)
	  for(int j=0;j<=k;j++,ii++) 
	    {
	      AutoDiff<3> vw; 
	  
		vw = pol_eta[j]*pol_zeta[k];
	  	
	      AutoDiff<3> hv = Cross(vw,xt); 
	      for(l=0;l<3;l++)
		shape(ii,l)= hv.DValue(l); 
	    }
	
	for(int  k = 0;k<= pp-2;k++)
	  for (int i=0;i<=k;i++,ii++)
	    {
	      AutoDiff<3> uw; 
	    
		uw = pol_xi[i]*pol_zeta[k];
	    	      
	      AutoDiff<3> hv = Cross(uw,yt); 
	      for(l=0;l<3;l++) shape(ii,l)=hv.DValue(l); 
	    }
	
	// 3rd component spans xi^i eta^j zeta^(k-1), i,j <= k
	// pol_zeta starts linear in zeta 
	// pol_xi and pol_eta quadratic in xi resp. eta 

	pol_zeta[0] = (1-z);
	for (int k=1;k<=pp;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];
	
	for(int k=0;k<= pp-1;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
	      {
		AutoDiff<3> pz; 
	
		  pz = pol_zeta[k]*pol_eta[j]*pol_xi[i]; 
		  
		AutoDiff<3> hv = Cross(pz,z);
		for(l=0;l<3;l++)
		  shape(ii,l) = hv.DValue(l); 
	      }
      }
#endif 
  
  / * (*testout) << " ip " << ip(0) << "\t" << ip(1) << "\t" << ip(2) << endl; 
  (*testout) << " curlshape code  " << endl << shape << endl ; 
  shape -= shape2; 
  (*testout) << " curlshape diff " << endl << shape << endl; 

 
  for(int ii=0;ii<ndof;ii++)
    for(int l=0;l<3;l++)
      if(fabs(shape(ii,l))>=1e-4) 
	{
	*testout << " errorcurl " << endl; 
	*testout << ii << " \t " << l << "\t: " << shape(ii,l) << endl; }
*/ 
  //use numcurl 
    //  shape = shape2; 
   
  }
#endif 
  
  template class  HCurlHighOrderFiniteElement<1>;
  template class  HCurlHighOrderFiniteElement<2>;
  template class  HCurlHighOrderFiniteElement<3>; 
  //template class  HCurlHighOrderTrig<TrigExtensionMonomial>;
  //template class  HCurlHighOrderTrig<TrigExtensionOptimal>; 
  //template class  HCurlHighOrderTrig<TrigExtensionMin>;
  template class  HCurlHighOrderSegm<IntegratedLegendreMonomialExt>; 
  template class  HCurlHighOrderTrig<IntegratedLegendreMonomialExt>;
  template class  HCurlHighOrderTet<IntegratedLegendreMonomialExt>; 
  //template class  HCurlHighOrderTet<TrigExtensionMonomial>;
  //template class  HCurlHighOrderTet<TrigExtensionOptimal>; 
  //template class  HCurlHighOrderTet<TrigExtensionMin>; 
  //template class  HCurlHighOrderPrism<TrigExtensionMonomial>;
  //template class  HCurlHighOrderPrism<TrigExtensionOptimal>; 
  //template class  HCurlHighOrderPrism<TrigExtensionMin>; 
  template class  HCurlHighOrderPrism<IntegratedLegendreMonomialExt>; 
  template class  HCurlHighOrderQuad<IntegratedLegendreMonomialExt>; 
  template class  HCurlHighOrderHex<IntegratedLegendreMonomialExt>; 
  template class  HCurlHighOrderPyr<IntegratedLegendreMonomialExt>; 
}
