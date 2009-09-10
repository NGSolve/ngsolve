/*********************************************************************/
/* File:   hcurlhofe.cpp                                             */
/* Author: Sabine Zaglmayr                                           */
/* Date:   20. Maerz 2003                                            */
/*                                                                   */
/* AutoDiff - revision: J. Schoeberl, March 2009                     */
/*********************************************************************/

#include <fem.hpp>    
#include <hcurlhofe.hpp>

namespace ngfem
{

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

    discontinuous = false;
  }

  template <int D>
  void HCurlHighOrderFiniteElement<D>::
  SetVertexNumbers (FlatArray<int> & avnums)
  {
    for (int i = 0; i < avnums.Size(); i++)
      vnums[i] = avnums[i];
  }
    
  template <int D>
  void HCurlHighOrderFiniteElement<D>::
  SetOrderCell (int oi)
  {
    order_cell = INT<3> (oi,oi,oi); 
  }

  template <int D>
  void HCurlHighOrderFiniteElement<D>::
  SetOrderCell (INT<3> oi)
  {
    order_cell = oi; 
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
  PrintInfo() const
  {
    (*testout) << "order_cell " << order_cell << " order_face ";
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






  /*******************************************/
  /* T_HCurlHOFiniteElement                  */
  /*******************************************/


  template <ELEMENT_TYPE ET>
  T_HCurlHighOrderFiniteElement<ET> :: 
  T_HCurlHighOrderFiniteElement (int aorder)
  {
    for (int i = 0; i < N_EDGE; i++)
      order_edge[i] = aorder;
    for (int i=0; i < N_FACE; i++) 
      order_face[i] = INT<2> (aorder,aorder); 
    if (DIM == 3)
      order_cell = INT<3> (aorder,aorder,aorder);

    for(int i = 0; i < N_EDGE; i++)
      usegrad_edge[i] = 1;
    for(int i=0; i < N_FACE; i++)
      usegrad_face[i] = 1;
    if (DIM == 3)
      usegrad_cell = 1;

    for (int i = 0; i < N_VERTEX; i++)
      vnums[i] = i;
    dimspace = DIM;
    eltype = ET;
  }

  
  template <ELEMENT_TYPE ET>
  void T_HCurlHighOrderFiniteElement<ET> :: 
  CalcShape (const IntegrationPoint & ip, FlatMatrixFixWidth<DIM> shape) const
  {    
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    HCurlShapeAssign<DIM> ds(shape); 
    static_cast<const HCurlHighOrderFE<ET>*> (this) -> T_CalcShape (adp, ds);
  }

  template <ELEMENT_TYPE ET>
  void T_HCurlHighOrderFiniteElement<ET> :: 
  CalcCurlShape (const IntegrationPoint & ip, FlatMatrixFixWidth<DIM_CURL> shape) const
  {  
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    HCurlCurlShapeAssign<DIM> ds(shape); 
    static_cast<const HCurlHighOrderFE<ET>*> (this) -> T_CalcShape (adp, ds);
  }



  template <ELEMENT_TYPE ET>
  void T_HCurlHighOrderFiniteElement<ET> :: 
  CalcMappedShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
                   FlatMatrixFixWidth<DIM> shape) const
  {
    AutoDiff<DIM> adp[DIM];

    for (int i = 0; i < DIM; i++)
      adp[i].Value() = sip.IP()(i);

    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
        adp[i].DValue(j) = sip.GetJacobianInverse()(i,j);

    HCurlShapeAssign<DIM> ds(shape); 
    static_cast<const HCurlHighOrderFE<ET>*> (this) -> T_CalcShape (adp, ds);
  }

  /// compute curl of shape
  template <ELEMENT_TYPE ET>
  void T_HCurlHighOrderFiniteElement<ET> :: 
  CalcMappedCurlShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
                       FlatMatrixFixWidth<DIM_CURL> curlshape) const
  { // not yet tested
    CalcCurlShape (sip.IP(), curlshape);
    if (DIM == 2)
      {
        curlshape /= sip.GetJacobiDet();        
      }
    else
      {
	AutoDiff<DIM> adp[DIM];
	
	for (int i = 0; i < DIM; i++)
	  adp[i].Value() = sip.IP()(i);
	
	for (int i = 0; i < DIM; i++)
	  for (int j = 0; j < DIM; j++)
	    adp[i].DValue(j) = sip.GetJacobianInverse()(i,j);

	HCurlCurlShapeAssign<DIM> ds(curlshape); 
	static_cast<const HCurlHighOrderFE<ET>*> (this) -> T_CalcShape (adp, ds);
	/*
        Mat<DIM> trans = (1.0/sip.GetJacobiDet()) * sip.GetJacobian();
        for (int i = 0; i < ndof; i++)
          {
            Vec<DIM> hs = curlshape.Row(i);
            curlshape.Row(i) = trans * hs;
          }
	*/
      }
  }

  template <ELEMENT_TYPE ET>
  Vec <DIM_CURL_TRAIT<ET_trait<ET>::DIM>::DIM>
  T_HCurlHighOrderFiniteElement<ET> :: 
  EvaluateCurlShape (const IntegrationPoint & ip, 
                     FlatVector<double> x,
                     LocalHeap & lh) const
  {
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    HCurlEvaluateCurl<DIM> ds(x); 
    static_cast<const HCurlHighOrderFE<ET>*> (this) -> T_CalcShape (adp, ds);
    return ds.Sum();
  }

  template <ELEMENT_TYPE ET>
  void T_HCurlHighOrderFiniteElement<ET> :: 
  ComputeNDof()
  {
    ndof = N_EDGE;

    for (int i = 0; i < N_EDGE; i++)
      if(order_edge[i] > 0)
        ndof += usegrad_edge[i]*order_edge[i];
    
    for(int i = 0; i < N_FACE; i++)
      if (FaceType(i) == ET_TRIG)
        {
          if (order_face[i][0] > 1) 
            ndof += ((usegrad_face[i]+1)*order_face[i][0]+2)*(order_face[i][0]-1)/2 ;
        }
      else
        {
          if(order_face[i][0]>=0 && order_face[i][1]>=0)
            ndof +=  (usegrad_face[i]+1)*order_face[i][0]*order_face[i][1] 
              + order_face[i][0] + order_face[i][1]; 
        }

    switch (ET)
      {
      case ET_TET: 
        if(order_cell[0] > 2)
          ndof += ((usegrad_cell + 2) * order_cell[0] + 3) 
            * (order_cell[0]-2) * (order_cell[0]-1) / 6; 
        break;
      case ET_PRISM:
        if(order_cell[2] > 0 && order_cell[0] > 1)
          ndof += ((usegrad_cell+2)*order_cell[2] + 1) * order_cell[0]*(order_cell[0]-1)/2
            + (order_cell[0]-1)*order_cell[2]; 
        break;
      case ET_PYRAMID:
        {
          int pc = order_cell[0]; //SZ: no problem to do anisotropic, but for the moment 
          // is it worth getting crazy :-) 
          if(order_cell[0]>1)
            ndof += usegrad_cell*(pc-1)*pc*(2*pc-1)/6 + pc*(2*pc*pc+3*pc-2)/3; 
          break;
        }
      case ET_HEX:
        if(order_cell[0] >= 0 && order_cell[1]>= 0 && order_cell[2]>=0)
          ndof += (usegrad_cell + 2)* order_cell[0] * order_cell[1] * order_cell[2]
            + order_cell[1]*order_cell[2]  + order_cell[0]*(order_cell[1] + order_cell[2]);  
        break;
      default:
        ;
      }


    order = 0; // max(order_edges,order_face,order_cell);  
    for (int i = 0; i < N_EDGE; i++)
      if (order_edge[i] > order)  order = order_edge[i];

    for(int i=0; i < N_FACE; i++) 
      if (ET_trait<ET>::FaceType(i) == ET_TRIG)
        {
          if (order_face[i][0] > order) 
            order = order_face[i][0]; 
        }
      else
        {
          for(int j=0;j<2;j++)
            if (order_face[i][j] > order) 
              order = order_face[i][j]; 
        }

    if (DIM == 3)
      for (int j = 0; j < 3; j++)
        if (order_cell[j] > order) 
          order = order_cell[j];
    

    // for integration order .. 
    if (ET == ET_PRISM || ET == ET_HEX || ET == ET_PYRAMID || ET == ET_QUAD)
      order++;
    else
      if (order==0) order++;
  }


  template <ELEMENT_TYPE ET>
  void T_HCurlHighOrderFiniteElement<ET> :: 
  GetInternalDofs (Array<int> & idofs) const
  {
    int ni = 0;

    if (discontinuous)
      {
        ni = ndof;
      }
    else
      {
        switch (ET)
          {
          case ET_TRIG:
            if (order_face[0][0] > 1) 
              ni = ((usegrad_face[0]+1)*order_face[0][0]+2)*(order_face[0][0]-1)/2 ;
            break;
          case ET_QUAD:
            if(order_face[0][0]>=0 && order_face[0][1]>=0)
              ni =  (usegrad_face[0]+1)*order_face[0][0]*order_face[0][1] 
                + order_face[0][0] + order_face[0][1]; 
            break;
          case ET_TET: 
            if(order_cell[0] > 2)
              ni = ((usegrad_cell + 2) * order_cell[0] + 3) 
                * (order_cell[0]-2) * (order_cell[0]-1) / 6; 
            break;
          case ET_PRISM:
            if(order_cell[2] > 0 && order_cell[0] > 1)
              ni = ((usegrad_cell+2)*order_cell[2] + 1) * order_cell[0]*(order_cell[0]-1)/2
                + (order_cell[0]-1)*order_cell[2]; 
            break;
          case ET_PYRAMID:
            {
              int pc = order_cell[0];
              if(pc > 1)
                ni = usegrad_cell*(pc-1)*pc*(2*pc-1)/6 + pc*(2*pc*pc+3*pc-2)/3; 
              break;
            }
          case ET_HEX:
            if(order_cell[0] >= 0 && order_cell[1]>= 0 && order_cell[2]>=0)
              ni = (usegrad_cell + 2)* order_cell[0] * order_cell[1] * order_cell[2]
                + order_cell[1]*order_cell[2]  + order_cell[0]*(order_cell[1] + order_cell[2]);  
            break;
          default:
            ;
          }
      }

    idofs.SetSize(ni);
    for (int i = 0; i < ni; i++)
      idofs[i] = ndof-ni+i;
  }



  //------------------------------------------------------------------------
  // HCurlHighOrderSegm
  //------------------------------------------------------------------------
  
  HCurlHighOrderFE<ET_SEGM> :: HCurlHighOrderFE ()
    : HCurlHighOrderFiniteElement<1>(ET_SEGM)
  { 
    // order_cell = INT<3> (aorder,aorder,aorder);
    usegrad_cell = 1; 
    // ComputeNDof(); 
  }

  HCurlHighOrderFE<ET_SEGM> :: HCurlHighOrderFE (int aorder)
    : HCurlHighOrderFiniteElement<1>(ET_SEGM)
  {
    order_cell = INT<3> (aorder,aorder,aorder);
    usegrad_cell = 1; 
    ComputeNDof(); 
  }

  void HCurlHighOrderFE<ET_SEGM> :: ComputeNDof()
  {
    ndof =1; 
    if(usegrad_cell)
      ndof += order_cell[0];    
    order = order_cell[0];
    order++; // integration order 
  }

  
  void HCurlHighOrderFE<ET_SEGM> :: CalcShape (const IntegrationPoint & ip, 
                                               FlatMatrixFixWidth<1> shape) const
  {
    AutoDiff<1> x (ip(0),0); 
    // AutoDiff<1> lami[2] = {1-x,x}; 
    AutoDiff<1> lami[2] = {x,1-x}; 
    
    int es = 0, ee =1;
    if (vnums[es] > vnums[ee]) swap(es,ee);  
    AutoDiff<1> xi = lami[ee] - lami[es]; 
    
    // Nedelec0-shapes
    shape(0,0) = 0.5*(xi.DValue(0)); 

    int ii = 1; 
    if (order_cell[0] >= 1 && usegrad_cell)
      { 
        ArrayMem<AutoDiff<1>, 10> pol_xi(order_cell[0]+2);
        T_ORTHOPOL::Calc (order_cell[0]+1, xi,pol_xi);  
	
        for (int j = 0; j < order_cell[0]; j++)
          shape(ii++,0) = pol_xi[j].DValue(0); 
      }    
  }
  
  //------------------------------------------------------------------------
  // HCurlHighOrderTrig
  //------------------------------------------------------------------------
  
  HCurlHighOrderFE<ET_TRIG> :: HCurlHighOrderFE (int aorder)
    : T_HCurlHighOrderFiniteElement<ET_TRIG> (aorder)
  {
    ComputeNDof();
  }
  
  template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE<ET_TRIG> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1];
    Tx lami[3] = { x, y, 1-x-y };

    ArrayMem<AutoDiff<2>,10> adpol1(order),adpol2(order);	
	
    int ii = 3; 
    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    for (int i = 0; i < 3; i++)
      {
	int es = edges[i][0], ee = edges[i][1];
	if (vnums[es] > vnums[ee])  swap (es, ee);

	//Nedelec low order edge shape function 
        shape[i] = uDv_minus_vDu<2> (lami[es], lami[ee]);

	int p = order_edge[i]; 
	//HO-Edge shapes (Gradient Fields)   
	if(p > 0 && usegrad_edge[i]) 
	  { 
	    AutoDiff<2> xi = lami[ee] - lami[es]; 
	    AutoDiff<2> eta = 1 - lami[ee] - lami[es]; 
	    T_ORTHOPOL::CalcTrigExt(p+1, xi, eta, adpol1); 
	   
	    for(int j = 0; j < p; j++) 
              shape[ii++] = Du<2> (adpol1[j]);
	  }
      }   
     
    //Inner shapes (Face) 
    int p = order_face[0][0];      
    if(p > 1) 
      {
	int fav[3] = { 0, 1, 2 }; 
	//Sort vertices ... v(f0) < v(f1) < v(f2) 
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	  

	AutoDiff<2> xi  = lami[fav[2]]-lami[fav[1]];
	AutoDiff<2> eta = lami[fav[0]]; 

	T_INNERSHAPES::CalcSplitted(p+1, xi, eta, adpol1,adpol2);
	
	// gradients:
	if(usegrad_face[0])
	  for (int j = 0; j < p-1; j++)
	    for (int k = 0; k < p-1-j; k++, ii++)
              shape[ii] = Du<2> (adpol1[j] * adpol2[k]);

	// other combination
	for (int j = 0; j < p-1; j++)
	  for (int k = 0; k < p-1-j; k++, ii++)
            shape[ii] = uDv_minus_vDu<2> (adpol2[k], adpol1[j]);
	
	// rec_pol * Nedelec0 
	for (int j = 0; j < p-1; j++, ii++)
          shape[ii] = wuDv_minus_wvDu<2> (lami[fav[1]], lami[fav[2]], adpol2[j]);
      }
  }


  //------------------------------------------------------------------------
  // HCurlHighOrderQuad
  //------------------------------------------------------------------------
  

  HCurlHighOrderFE<ET_QUAD> :: HCurlHighOrderFE (int aorder)
    : T_HCurlHighOrderFiniteElement<ET_QUAD> (aorder)
  {
    ComputeNDof();
  }
  

  template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE<ET_QUAD> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1];

    AutoDiff<2> lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  

    int ii = 4;
    ArrayMem<AutoDiff<2>, 10> pol_xi(order+2), pol_eta(order+2);

    // edges
    const EDGE * edges = ElementTopology::GetEdges (ET_QUAD);
    for (int i = 0; i < 4; i++)
      {
	int p = order_edge[i]; 
	int es = edges[i][0], ee = edges[i][1];
	if (vnums[es] > vnums[ee]) swap (es, ee);

	AutoDiff<2> xi  = sigma[ee]-sigma[es];
	AutoDiff<2> lam_e = lami[ee]+lami[es];  

	// Nedelec0-shapes
        shape[i] = uDv<2> (0.5 * lam_e, xi); 

	// High Order edges ... Gradient fields 
	if(usegrad_edge[i])
	  {
	    T_ORTHOPOL::Calc (p+1, xi, pol_xi);  
	    for (int j = 0; j < p; j++)
              shape[ii++] = Du<2> (pol_xi[j] * lam_e);
	  }
      }
     
    INT<2> p = order_face[0]; // (order_cell[0],order_cell[1]);
    int fmax = 0; 
    for (int j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax])
	fmax = j;
    
    int f1 = (fmax+3)%4; 
    int f2 = (fmax+1)%4; 
    if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f2 > f1; 

    AutoDiff<2> xi = sigma[fmax]-sigma[f1];  // in [-1,1]
    AutoDiff<2> eta = sigma[fmax]-sigma[f2]; // in [-1,1]
    
    T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
    T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);
    
    //Gradient fields 
    if(usegrad_face[0])
      for (int k = 0; k < p[0]; k++)
	for (int j= 0; j < p[1]; j++)
          shape[ii++] = Du<2> (pol_xi[k]*pol_eta[j]);

    //Rotation of Gradient fields 
    for (int k = 0; k < p[0]; k++)
      for (int j= 0; j < p[1]; j++)
        shape[ii++] = uDv_minus_vDu<2> (pol_eta[j], pol_xi[k]);

    //Missing ones 
    for(int j = 0; j< p[0]; j++)
      shape[ii++] = uDv<2> (0.5*pol_xi[j], eta);

    for(int j = 0; j < p[1]; j++)
      shape[ii++] = uDv<2> (0.5*pol_eta[j], xi); 
  }


  
  //------------------------------------------------------------------------
  //        Tetrahedron
  //------------------------------------------------------------------------
 

  HCurlHighOrderFE<ET_TET> :: HCurlHighOrderFE (int aorder)
    : T_HCurlHighOrderFiniteElement<ET_TET> (aorder)
  {
    ComputeNDof();
  }

  template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE<ET_TET> :: T_CalcShape (Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];
    Tx lami[4] = { x, y, z, 1-x-y-z };

    ArrayMem<AutoDiff<3>,10> adpol1(order+2),adpol2(order+2),adpol3(order+2); 
    int ii = 6; 

    for (int i = 0; i < N_EDGE; i++)
      { 
	int p = order_edge[i]; 

        INT<2> e = GetEdgeSort (i, vnums);	  
	
	//Nedelec low order edge shape function 
        shape[i] = uDv_minus_vDu<3> (lami[e[0]], lami[e[1]]);

	//HO-Edge shape functions (Gradient Fields) 	
	if (p > 0 && usegrad_edge[i]) 
	  {
	    AutoDiff<3> xi = lami[e[1]]-lami[e[0]]; 
	    AutoDiff<3> eta = 1-lami[e[1]]-lami[e[0]]; 
	    T_ORTHOPOL::CalcTrigExt(p+1,xi,eta, adpol1); 

	    for(int j = 0; j < p; j++) 	      
              shape[ii++] = Du<3> (adpol1[j]);
	  }
      }

    // face shape functions
    for(int i = 0; i < N_FACE; i++) 
      if (order_face[i][0] >= 2)
        {
          INT<4> fav = GetFaceSort (i, vnums);
          
          int vop = 6 - fav[0] - fav[1] - fav[2];  	
          int p = order_face[i][0];
          
          AutoDiff<3> xi = lami[fav[2]]-lami[fav[1]];
          AutoDiff<3> eta = lami[fav[0]]; // lo 
          AutoDiff<3> zeta = lami[vop];   // lz 
          
          T_FACESHAPES::CalcSplitted (p+1, xi, eta, zeta, adpol1, adpol2); 
          
          // gradients 
          if (usegrad_face[i])
            for (int j = 0; j <= p-2; j++)
              for (int k = 0; k <= p-2-j; k++, ii++)
                shape[ii] = Du<3> (adpol1[j] * adpol2[k]);
          
          // other combination
          for (int j = 0; j <= p-2; j++)
            for (int k = 0; k <= p-2-j; k++, ii++)
              shape[ii] = uDv_minus_vDu<3> (adpol2[k], adpol1[j]);
          
          // type 3
          for (int j = 0; j <= p-2; j++, ii++)
            shape[ii] = wuDv_minus_wvDu<3> (lami[fav[1]], lami[fav[2]], adpol2[j]);
        }

    
    int p = order_cell[0]; 
    
    T_INNERSHAPES::CalcSplitted(p+1, x-(1-x-y-z), y, z,adpol1, adpol2, adpol3 );
    
    //gradient fields 
    if(usegrad_cell)
      for (int i = 0; i <= p-3; i++)
	for (int j = 0; j <= p-3-i; j++)
	  for (int k = 0; k <= p-3-i-j; k++)
            shape[ii++] = Du<3> (adpol1[i] * adpol2[j] * adpol3[k]);

    // other combinations
    for (int i = 0; i <= p-3; i++)
      for (int j = 0; j <= p-3-i; j++)
	for (int k = 0; k <= p-3-i-j; k++)
          { // not Sabine's original ...
            shape[ii++] = uDv_minus_vDu<3> (adpol1[i], adpol2[j] * adpol3[k]);
            shape[ii++] = uDv_minus_vDu<3> (adpol1[i] * adpol3[k], adpol2[j]);
          }
       
    for (int j= 0; j <= p-3; j++)
      for (int k = 0; k <= p-3-j; k++)
        shape[ii++] = wuDv_minus_wvDu<3> (lami[0], lami[3], adpol2[j] * adpol3[k]);
  }


		        
  //------------------------------------------------------------------------
  //                   Prism
  //------------------------------------------------------------------------

  HCurlHighOrderFE<ET_PRISM> :: HCurlHighOrderFE (int aorder)
    : T_HCurlHighOrderFiniteElement<ET_PRISM> (aorder)
  {
    ComputeNDof();
  }

  template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE<ET_PRISM> :: T_CalcShape (Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];

    AutoDiff<3> lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    AutoDiff<3> muz[6]  = { 1-z, 1-z, 1-z, z, z, z };
       
       
    const EDGE * edges = ElementTopology::GetEdges (ET_PRISM);
    const FACE * faces = ElementTopology::GetFaces (ET_PRISM); 

    ArrayMem<AutoDiff<3>,20> adpolxy1(order+3),adpolxy2(order+3); 
    ArrayMem<AutoDiff<3>,20> adpolz(order+3);   

    int ii = 9;
    
    // horizontal edge shapes
    for (int i = 0; i < 6; i++)
      {
	int p = order_edge[i]; 
        int es = edges[i][0], ee = edges[i][1]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	//Nedelec0
        shape[i] = wuDv_minus_wvDu<3> (lami[es], lami[ee], muz[ee]);
   
	//high order \nabla (P_edge(x,y) * muz)
	if(usegrad_edge[i])
	  {
	    T_ORTHOPOL::CalcTrigExt(p+1, lami[ee]-lami[es],
				    1-lami[es]-lami[ee],adpolxy1);
	    
	    for(int j = 0; j <= p-1; j++)
              shape[ii++] = Du<3> (adpolxy1[j] * muz[ee]);
	  }
      }

    //Vertical Edge Shapes
    for (int i = 6; i < 9; i++)
      {
	int p = order_edge[i]; 
        int es = edges[i][0], ee = edges[i][1];
	if (vnums[es] > vnums[ee]) swap (es, ee);

        shape[i] = wuDv_minus_wvDu<3> (muz[es], muz[ee], lami[ee]);
	
	//high order edges:  \nabla (T_ORTHOPOL^{p+1}(2z-1) * lami(x,y))
	if(usegrad_edge[i])
	  {
	    T_ORTHOPOL::Calc (p+1, muz[ee]-muz[es], adpolz);
	    
	    for (int j = 0; j < p; j++)
              shape[ii++] = Du<3> (adpolz[j] * lami[ee]);
	  }
      }


    // trig face shapes
    for (int i = 0; i < 2; i++)
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

	if(usegrad_face[i])
	  // gradient-fields =>  \nabla( adpolxy1*adpolxy2*muz )
	  for (int j = 0; j <= p-2; j++)
	    for (int k = 0; k <= p-2-j; k++)
              shape[ii++] = Du<3> (adpolxy1[j]*adpolxy2[k] * muz[fav[2]]);
	
	// rotations of grad-fields => grad(uj)*vk*w -  uj*grad(vk)*w 
	for (int j = 0; j <= p-2; j++)
	  for (int k = 0; k <= p-2-j; k++)
            shape[ii++] = wuDv_minus_wvDu<3> (adpolxy2[k], adpolxy1[j], muz[fav[2]]);

	//  Ned0*adpolxy2[j]*muz 
	for (int j = 0; j <= p-2; j++,ii++)
          shape[ii] = wuDv_minus_wvDu<3> (lami[fav[1]], lami[fav[2]], adpolxy2[j]*muz[fav[2]]);
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
	AutoDiff<3> xi = lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
	AutoDiff<3> eta = 1-lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
	AutoDiff<3> zeta = muz[faces[i][fmax]]-muz[faces[i][fz]]; 
	
	int pp = int(max2(p[0],p[1]))+1;
	T_ORTHOPOL::CalcTrigExt(pp,xi,eta,adpolxy1); 
	T_ORTHOPOL::Calc(pp,zeta,adpolz); 

	if(usegrad_face[i])
	  {
	    // Gradientfields nabla(polxy*polz) 
	    if (vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) 
	      for (int k = 0; k <= p[0]-1; k++)
		for (int j = 0; j <= p[1]-1; j++)
                  shape[ii++] = Du<3> (adpolxy1[k] * adpolz[j]);
	    else
	      for (int j = 0; j <= p[0]-1; j++)
		for (int k = 0; k <= p[1]-1; k++)
                  shape[ii++] = Du<3> (adpolxy1[k] * adpolz[j]);
	  }
	  
	// Rotations of GradFields => nabla(polxy)*polz - polxy*nabla(polz)
	if (vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) 
	  for (int k = 0; k <= p[0]-1; k++)
	    for (int j = 0; j <= p[1]-1; j++)
              shape[ii++] = uDv_minus_vDu<3> (adpolz[j], adpolxy1[k]);
	else
	  for (int j = 0; j <= p[0]-1; j++)
	    for (int k = 0; k <= p[1]-1; k++)
              shape[ii++] = uDv_minus_vDu<3> (adpolxy1[k], adpolz[j]);
	
	// Type 3 
	// (ned0_trig)*polz, (ned0_quad)* polxy 

	if(vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) // p = (p_trig,p_z) 
	  {
	    for(int j=0;j<=p[0]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu<3> (muz[faces[i][fz]], muz[faces[i][fmax]], adpolxy1[j]);
	    for(int j=0;j<=p[1]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu<3> (lami[faces[i][ftrig]], lami[faces[i][fmax]], adpolz[j]);
	  }
	else 
	  {
	    for(int j=0;j<=p[0]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu<3> (lami[faces[i][ftrig]], lami[faces[i][fmax]], adpolz[j]);
	    for(int j=0;j<=p[1]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu<3> (muz[faces[i][fz]], muz[faces[i][fmax]], adpolxy1[j]);
	  }
      }
    

    if(order_cell[0] > 1 && order_cell[2] > 0) 
      {
	T_TRIGFACESHAPES::CalcSplitted(order_cell[0]+1,x-y,1-x-y,adpolxy1,adpolxy2);
	T_ORTHOPOL::Calc(order_cell[2]+1,2*z-1,adpolz); 
	
	// gradientfields
	if(usegrad_cell)
	  for(int i=0;i<=order_cell[0]-2;i++)
	    for(int j=0;j<=order_cell[0]-2-i;j++)
	      for(int k=0;k<=order_cell[2]-1;k++)
                shape[ii++] = Du<3> (adpolxy1[i]*adpolxy2[j]*adpolz[k]);

	// Rotations of gradientfields
	for(int i=0;i<=order_cell[0]-2;i++)
	  for(int j=0;j<=order_cell[0]-2-i;j++)
	    for(int k=0;k<=order_cell[2]-1;k++)
	      {
                shape[ii++] = wuDv_minus_wvDu<3> (adpolxy1[i],adpolxy2[j],adpolz[k]);
                shape[ii++] = uDv_minus_vDu<3> (adpolxy1[i],adpolxy2[j]*adpolz[k]);
	      }

	// Type 3 
	// ned0(trig) * polxy2[j]*polz 
	// z.DValue(0) * polxy1[i] * polxy2[j] 
	// double ned_trig[2] = {y.Value(),-x.Value()};  
	for(int j=0;j<=order_cell[0]-2;j++) 
	  for (int k=0;k<=order_cell[2]-1;k++) 
            shape[ii++] = wuDv_minus_wvDu<3> (x,y, adpolxy2[j]*adpolz[k]);

    	for(int i = 0; i <= order_cell[0]-2; i++) 
	  for(int j = 0; j <= order_cell[0]-2-i; j++) 
            shape[ii++] = wuDv_minus_wvDu<3> (z,1-z, adpolxy1[i]*adpolxy2[j]);
      }
  }



  //------------------------------------------------------------------------
  // HCurlHighOrderHex
  //------------------------------------------------------------------------

  HCurlHighOrderFE<ET_HEX> :: HCurlHighOrderFE (int aorder)
    : T_HCurlHighOrderFiniteElement<ET_HEX> (aorder)
  {
    ComputeNDof();
  }
  
  template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE<ET_HEX> :: T_CalcShape (Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];

    AutoDiff<3> lami[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
			 (1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
    AutoDiff<3> sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
			  (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 
     
    int ii = 12; 
    ArrayMem<AutoDiff<3>, 20> pol_xi(order+2),pol_eta(order+2),pol_zeta(order+2);
   
    // edges
    const EDGE * edges = ElementTopology::GetEdges (ET_HEX);
    for (int i = 0; i < 12; i++)
      {
	int es = edges[i][0], ee = edges[i][1];
	int p = order_edge[i]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);

	AutoDiff<3> xi = sigma[ee]-sigma[es]; 
	AutoDiff<3> lam_e = lami[ee]+lami[es]; 
	
	// Nedelec0-shapes
        shape[i] = uDv<3> (0.5*lam_e, xi);

	if(usegrad_edge[i])
	  {
	    // High Order edges ... Gradient fields 
	    T_ORTHOPOL::Calc (p+1, xi, pol_xi);
	    
	    for (int j = 0; j < p; j++)
              shape[ii++] = Du<3> (pol_xi[j] * lam_e);
	  }
      }
    
    //Faces 
    const FACE * faces = ElementTopology::GetFaces (ET_HEX);
    for (int i = 0; i<6; i++)
      {
	INT<2> p = order_face[i];
	
	AutoDiff<3> lam_f = 0;
	for (int j = 0; j < 4; j++)
	  lam_f += lami[faces[i][j]];
	
	int qmax = 0;
	for (int j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] > vnums[faces[i][qmax]])
	    qmax = j;
	
	int q1 = (qmax+3)%4; 
	int q2 = (qmax+1)%4; 

	if(vnums[faces[i][q2]] > vnums[faces[i][q1]])
	  swap(q1,q2);  // fmax > f1 > f2

	int fmax = faces[i][qmax]; 
	int f1 = faces[i][q1]; 
	int f2 = faces[i][q2]; 
	      
	AutoDiff<3> xi = sigma[fmax]-sigma[f1]; 
	AutoDiff<3> eta = sigma[fmax]-sigma[f2]; 
    
	T_ORTHOPOL::Calc(p[0]+1, xi,pol_xi);
	T_ORTHOPOL::Calc(p[1]+1,eta,pol_eta);
	
	//Gradient fields 
	if(usegrad_face[i])
	  for (int k = 0; k < p[0]; k++)
	    for (int j= 0; j < p[1]; j++)
              shape[ii++] = Du<3> (lam_f * pol_xi[k] * pol_eta[j]);
	
	//Rotation of Gradient fields 
	for (int k = 0; k < p[0]; k++)
	  for (int j= 0; j < p[1]; j++)
            shape[ii++] = uDv_minus_vDu<3> (pol_eta[j], lam_f * pol_xi[k]);

	// Missing ones 
	for(int j = 0; j < p[0];j++) 
          shape[ii++] = wuDv_minus_wvDu<3> (0.5, eta, pol_xi[j]*lam_f); 

	for(int j = 0; j < p[1];j++) 
          shape[ii++] = wuDv_minus_wvDu<3> (0.5, xi, pol_eta[j]*lam_f); 
      }
    
    // Element-based shapes
    T_ORTHOPOL::Calc(order_cell[0]+1,2*x-1,pol_xi);
    T_ORTHOPOL::Calc(order_cell[1]+1,2*y-1,pol_eta);
    T_ORTHOPOL::Calc(order_cell[2]+1,2*z-1,pol_zeta); 
    
    //Gradient fields
    if(usegrad_cell)
      for (int i=0; i<order_cell[0]; i++)
	for(int j=0; j<order_cell[1]; j++) 
	  for(int k=0; k<order_cell[2]; k++)
            shape[ii++] = Du<3> (pol_xi[i] * pol_eta[j] * pol_zeta[k]);

    //Rotations of gradient fields
    for (int i=0; i<order_cell[0]; i++)
      for(int j=0; j<order_cell[1]; j++) 
	for(int k=0; k<order_cell[2]; k++)
	  {
            shape[ii++] = uDv_minus_vDu<3> (pol_xi[i] * pol_eta[j], pol_zeta[k]);
            shape[ii++] = uDv_minus_vDu<3> (pol_xi[i], pol_eta[j] * pol_zeta[k]);
	  } 
    
    for(int i = 0; i < order_cell[0]; i++) 
      for(int j = 0; j < order_cell[1]; j++)
        shape[ii++] = wuDv_minus_wvDu<3> (z,1-z,pol_xi[i] * pol_eta[j]);

    for(int i = 0; i < order_cell[0]; i++) 
      for(int k = 0; k < order_cell[2]; k++)
        shape[ii++] = wuDv_minus_wvDu<3> (y,1-y,pol_xi[i] * pol_zeta[k]);

    for(int j = 0; j < order_cell[1]; j++)
      for(int k = 0; k < order_cell[2]; k++)
        shape[ii++] = wuDv_minus_wvDu<3> (x,1-x,pol_eta[j] * pol_zeta[k]);
  }



  //------------------------------------------------------------------------
  //            Pyramid
  //------------------------------------------------------------------------

  HCurlHighOrderFE<ET_PYRAMID> :: HCurlHighOrderFE (int aorder)
    : T_HCurlHighOrderFiniteElement<ET_PYRAMID> (aorder)
  {
    ComputeNDof();
  }

  template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE<ET_PYRAMID> :: T_CalcShape (Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];

    if(z.Value()==1.) z.Value() -=1.e-8; 


    AutoDiff<3> xt = x/(1-z); 
    AutoDiff<3> yt = y/(1-z); 
    AutoDiff<3> sigma[5] = {(1-xt)+(1-yt)+(1-z),xt+(1-yt)+(1-z), xt + yt + (1-z), 
			    (1-xt)+yt+(1-z),z}; 

    AutoDiff<3> lami[5] = {(1-xt)*(1-yt)*(1-z),xt*(1-yt)*(1-z), xt * yt * (1-z), 
			   (1-xt)*yt*(1-z),z}; 

    /*
    AutoDiff<3> sig[5] = {(1-xt)+(1-yt),xt+(1-yt), xt + yt, 
			    (1-xt)+yt,z}; 
    */

    AutoDiff<3> lambda[5] = {(1-xt)*(1-yt),xt*(1-yt), xt * yt, 
			   (1-xt)*yt,z}; 
        
    
       
    const EDGE * edges = ElementTopology::GetEdges (ET_PYRAMID);
    
    ArrayMem<AutoDiff<3>, 20> pol_xi(order+2), pol_eta(order+2), pol_zeta(order+2),pol2_zeta(order+2); 
    
    int ii =8; 
 
    // horizontal edges incl. Nedelec 0
    for (int i = 0; i < 4; i++)
      {
	int es = edges[i][0], ee = edges[i][1]; 
	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	AutoDiff<3> xi  = sigma[ee] - sigma[es];   
	// AutoDiff<3> lam_e = lami[ee] + lami[es];
	AutoDiff<3> lam_t = lambda[ee] + lambda[es]; 
	
        shape[i] = uDv<3> (0.5 * (1-z)*(1-z)*lam_t, xi);

	if(usegrad_edge[i])
	  {
	    T_ORTHOPOL::CalcTrigExt (order_edge[i]+1, xi*(1-z), z, pol_xi);
	    
	    for(int j = 0; j < order_edge[i]; j++)
              shape[ii++] = Du<3> (pol_xi[j]*lam_t);
	  }
      }
    
    // vertical edges incl. Nedelec 0  
    for(int i = 4; i < 8; i++)
      {
	int es = edges[i][0], ee = edges[i][1]; 
	
	if (vnums[es] > vnums[ee]) swap (es, ee);

        shape[i] = uDv_minus_vDu<3> (lami[es], lami[ee]);

	if (usegrad_edge[i])
	  {
	    T_ORTHOPOL::CalcTrigExt (order_edge[i]+1, lami[ee]-lami[es],  
                                     1-lami[es]-lami[ee], pol_xi);
	    for(int j = 0; j < order_edge[i]; j++)
              shape[ii++] = Du<3> (pol_xi[j]);
	  }
      }

    const FACE * faces = ElementTopology::GetFaces (ET_PYRAMID); 

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
	      for(int k=0;k<=p-2-j; k++)
                shape[ii++] = Du<3> (pol_xi[j] * pol_eta[k]);

	  // Type 2:  
	  for(int j=0;j<= p-2; j++)
	    for(int k=0;k<=p-2-j; k++)
              shape[ii++] = uDv_minus_vDu<3> (pol_eta[k], pol_xi[j]);

	  // Type 3: Nedelec-based ones (Ned_0*v_j)
	  for(int j=0;j<=p-2;j++)
            shape[ii++] = wuDv_minus_wvDu<3> (bary[fav[1]], bary[fav[2]], pol_eta[j]);
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
              shape[ii] = Du<3> (pol_xi[k] * pol_eta[j]);

	// Type 2: 
	for (int k = 0; k < px; k++) 
	  for (int j = 0; j < py; j++, ii++) 
            shape[ii] = uDv_minus_vDu<3> (pol_eta[j], pol_xi[k]);

	// Type 3:
	for (int k=0; k< px; k++,ii++)
          shape[ii] = uDv<3> (0.5*pol_xi[k]*fac, eta);

	for (int k=0; k< py; k++,ii++)
          shape[ii] = uDv<3> (0.5*pol_eta[k]*fac, xi);  // shouldn't there be a  *fac  ????
      }


    if (order_cell[0] >= 2)
      {
	int pp = order_cell[0];
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
                  shape[ii] = Du<3> (pol_xi[i]*pol_eta[j]*pol_zeta[k]);
	  }
	
	// Type 2a: l.i. combinations of grad-terms   
	// shape = u_i \nabla(v_j) w_k 
	// shape = u_i v_j \nabla(w_k) 
	for(int k=0;k<= pp-2;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
              shape[ii] = uDv<3> (pol_xi[i]*pol_zeta[k], pol_eta[j]);

	// Type 2b: shape = v_j w_k \nabla (xt) 
	//          shape = u_i w_k \nabla (yt)
	for(int k = 0;k<= pp-2;k++)
	  for(int j=0;j<=k;j++) 
            shape[ii++] = uDv<3> (pol_eta[j]*pol_zeta[k], xt);
	
	for(int  k = 0;k<= pp-2;k++)
	  for (int i=0;i<=k;i++)
            shape[ii++] = uDv<3> (pol_xi[i]*pol_zeta[k], yt);
	
	// 3rd component spans xi^i eta^j zeta^(k-1), i,j <= k
	// pol_zeta starts linear in zeta 
	// pol_xi and pol_eta quadratic in xi resp. eta 
	pol_zeta[0] = (1-z);
	for (int k=1;k<=pp;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];
     
	for(int k=0;k<= pp-1;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
              shape[ii] = uDv<3> (pol_eta[j] * pol_xi[i] * pol_zeta[k], z);
      }
  }
  
  template class  HCurlHighOrderFiniteElement<1>;
  template class  HCurlHighOrderFiniteElement<2>;
  template class  HCurlHighOrderFiniteElement<3>; 
}
