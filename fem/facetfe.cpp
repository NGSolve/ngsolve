#include <fem.hpp>

namespace ngfem {
  
/****************************************************************************
 * FacetVolumeElement 
 ****************************************************************************/

template<int D>
FacetVolumeFiniteElement<D>::FacetVolumeFiniteElement (int adim, ELEMENT_TYPE aeltype) : NodalFiniteElement<D>(adim, aeltype,-1,-1)
{
  for (int i=0; i<8; i++)
    vnums[i] = -1; //HERBERT: warum i
  for (int i=0; i<6; i++)
    facet_order[i] = -1;
  for (int i=0; i<=8; i++)
    first_facet_dof[i] = 0;

  nv = ElementTopology :: GetNVertices(NodalFiniteElement<D>::eltype);
  nf = ElementTopology :: GetNFacets(NodalFiniteElement<D>::eltype);
}
 
template<int D>
void FacetVolumeFiniteElement<D>::SetVertexNumbers (FlatArray<int> & avnums)
{
  for (int i=0; i<avnums.Size(); i++)
    vnums[i] = avnums[i];
/*  cerr << "*** set vertex numbers to ";
  for (int i=0; i<8; i++) cerr << vnums[i] << ", ";
  cerr << endl;*/
}
 
template<int D>
void FacetVolumeFiniteElement<D>::SetOrder(int ao)
{
  NodalFiniteElement<D>::order = ao;
  for (int i=0; i<6; i++)
    facet_order[i] = ao;
  ComputeNDof();
}

template<int D>
void FacetVolumeFiniteElement<D>::SetOrder(FlatArray<int> & ao)
{
  for (int i=0; i<ao.Size(); i++)
    facet_order[i] = ao[i];
  
  NodalFiniteElement<D>::order = facet_order[0];        // integration order (JS, Spet 07)
  for (int i = 1; i < ao.Size(); i++)
    NodalFiniteElement<D>::order = max(NodalFiniteElement<D>::order, ao[i]);

  ComputeNDof();
}

template<int D>
void FacetVolumeFiniteElement<D>::GetFacetDofNrs(int fnr, Array<int>& dnums) const
{
  dnums.SetSize(0);
  for (int i=first_facet_dof[fnr]; i<first_facet_dof[fnr+1]; i++)
    dnums.Append(i);
}


template<int D>
void FacetVolumeFiniteElement<D>::GetVertexNumbers(Array<int> &vn) const
{
  vn.SetSize(nv);
  for (int i=0; i<nv; i++)
    vn[i] = vnums[i];
}

template<int D>
void FacetVolumeFiniteElement<D>::GetFacetOrders(Array<int>& fo) const
{
  fo.SetSize(6);
  for (int i=0; i<6; i++)
    fo[i] = facet_order[i];
}
 

/****************************************************************************
 * FacetFacetElement 
 ****************************************************************************/

/*FacetFacetFiniteElement::FacetFacetFiniteElement (int adim, ELEMENT_TYPE aeltype) : FiniteElement(adim, aeltype, -1, -1)
{
  for (int i=0; i<4; i++)
    vnums[i]=i;
}*/
       
// void FacetFacetFiniteElement::SetVertexNumbers (FlatArray<int> & avnums)
// {
//   for (int i=0; i<avnums.Size(); i++)
//     vnums[i] = avnums[i];
// }





/****************************************************************************
 * FacetFacetSegm 
 ****************************************************************************/

// ----------------------------------------------------------------------------
// void FacetFacetSegm::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
// {
//   double x = ip(0);
// 
//   if (vnums[0] > vnums[1])
//   {
//     x   = 1-x; // keep the right orientation
//   }
//   ArrayMem<double, 20> polx(order+1);
// 
//   LegendrePolynomial(order, 1-2*x, polx);  // basis functions
//   for (int i=0; i<order+1; i++)
//     shape(i) = polx[i];
// }

// void FacetFacetSegm::ComputeNDof() 
// {
//   ndof = order+1;
// }


/****************************************************************************
 * FacetFacetTrig
 ****************************************************************************/

// void FacetFacetTrig::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
// {
//   double x = ip(0), y = ip(1);
// 
// //   cout << "  ++ FacetFacetTrig: order=" << order << endl;
//   // 0 order
//   if ( order==0 )
//   {
//     shape(0) = 1;
//     return;
//   }
//    
//   double lami[3] = {x, y, 1-x-y};
//   shape = 0;
// 
//   // sort according to global orientation: ala h1hofe.cpp
//   int fav[3] = { 0, 1, 2 }; 
//   if(vnums[fav[0]] > vnums[fav[1]]) 
//     swap(fav[0],fav[1]); 
//   if(vnums[fav[1]] > vnums[fav[2]]) 
//     swap(fav[1],fav[2]);
//   if(vnums[fav[0]] > vnums[fav[1]]) 
//     swap(fav[0],fav[1]); 	
//   
//   double xi=lami[fav[0]], eta=lami[fav[1]];
//   
//   // calculate basis ala l2hofe
//   int ii = 0, p=order;
//   ArrayMem<double, 20> polx(p+1), poly(p+1);
//     
//   ScaledLegendrePolynomial (p, 2*xi+eta-1, 1-eta, polx);
//   LegendrePolynomial (p, 2*eta-1, poly);
//     
//   for (int i = 0; i <= p; i++)
//     for (int j = 0; j <= p-i; j++)
//       shape(ii++) = polx[i] * poly[j];
// }

// void FacetFacetTrig::ComputeNDof() 
// {
//   ndof = ((order+2) * (order+1)) / 2;;
// }


/****************************************************************************
 * FacetFacetQuad
 ****************************************************************************/

// void FacetFacetQuad::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
// {
//   double x = ip(0), y = ip(1);
// 
//   // 0 order
//   if ( order==0 )
//   {
//     shape(0) = 1;
//     return;
//   }
//    
// 
//   // sort according to global orientation
//   int smallest=0;
//   for (int i=1; i<4; i++)
//   {
//     if (vnums[i] < vnums[smallest])
//       smallest = i;
//   }
//    
//   switch (smallest)
//   {
//     case 0:
//       if (vnums[1] > vnums[3])
//         swap(x,y);
//       break;
//     case 1:
//       x=1-x;
//       if (vnums[2] < vnums[0])
//         swap(x,y);
//       break;
//     case 2:
//       x=1-x;
//       y=1-y;
//       if (vnums[3] > vnums[1])
//         swap(x,y);
//       break;
//     case 3:
//       y=1-y;
//       if (vnums[0] < vnums[2])
//         swap(x,y);
//       break;
//     default:
//       cerr << "*** error in void FacetFacetQuad::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const" << endl;
//       exit(0);
//   }
// 
//   
//   // calculate basis ala l2hofe
//   int ii = 0, p=order;
//   ArrayMem<double, 20> polx(p+1), poly(p+1);
//     
//   LegendrePolynomial (p, 2*x-1, polx);
//   LegendrePolynomial (p, 2*y-1, poly);
//     
//   for (int i = 0; i <= p; i++)
//     for (int j = 0; j <= p; j++)
//       shape(ii++) = polx[i] * poly[j];
// }

// void FacetFacetQuad::ComputeNDof() 
// {
//   ndof = (order+1) * (order+1);
// }


/****************************************************************************
 * FacetVolumeTrig
 ****************************************************************************/

void FacetVolumeTrig::CalcShape (const IntegrationPoint & ip2d, FlatVector<> shape) const
{  
//   topology: points = 0:()
  IntegrationPoint ip1d;  
  double x=ip2d(0), y=ip2d(1);

  shape = 0.0;
  if (fabs(y)<1e-12) 
  {
    ip1d(0) = x;
//     CalcFacetShape(0, ip1d, shape.Range(first_facet_dof[0], first_facet_dof[1]));
    CalcFacetShape(0, ip1d, shape);
  }
  else if (fabs(x)<1e-12)
  {
    ip1d(0) = 1-y;
//     CalcFacetShape(1, ip1d, shape.Range(first_facet_dof[1], first_facet_dof[2]));
    CalcFacetShape(1, ip1d, shape);
  }
  else if (fabs(1-x-y)<1e-12)
  {
    ip1d(0) = (1-x+y)*0.5;
//     CalcFacetShape(2, ip1d, shape.Range(first_facet_dof[2], first_facet_dof[3]));
    CalcFacetShape(2, ip1d, shape);
  }
} 
 
void FacetVolumeTrig::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
{
  SetFacet(afnr);
  shape=0.0;
  facet.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
}
   
void FacetVolumeTrig::SetFacet(int afnr) const
{
  if (fnr == afnr) return;
  
  FacetVolumeTrig * trig=const_cast<FacetVolumeTrig*>(this);
  trig->fnr = afnr;
  const EDGE * edges = ElementTopology::GetEdges (eltype);
  Array<int> fvnums(2);
  fvnums[0] = vnums[edges[fnr][0]]; 
  fvnums[1] = vnums[edges[fnr][1]]; 
   
  trig->facet.SetVertexNumbers(fvnums);
  trig->facet.SetOrder(facet_order[fnr]);
}  

void FacetVolumeTrig::ComputeNDof() 
{
  ndof = 0;
  for (int i=0; i<nf; i++)
  {
    first_facet_dof[i] = ndof;
    ndof += facet_order[i]+1;
  }
  first_facet_dof[nf] = ndof;
}


/****************************************************************************
 * FacetVolumeQuad
 ****************************************************************************/
// -------------------------------------------------------------------------------
void FacetVolumeQuad::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
{
  // topology: points = 0:(0,0), 1:(1,0), 2:(1,1), 3:(0,1)
  //           edges = 0:(0,1), 1:(2,3), 2:(3,0), 4:(1,2)
  IntegrationPoint ip1d;  
  double x=ip(0), y=ip(1);
  int fnr=-1;

  shape = 0.0;
  if (fabs(y)<1e-12)
  {
    fnr=0;
    ip1d(0) = x;
//     CalcFacetShape(0, ip1d, shape.Range(first_facet_dof[0], first_facet_dof[1]));
    CalcFacetShape(0, ip1d, shape);
  }
  else if (fabs(y)>1-(1e-12))
  {
    fnr=1;
    ip1d(0) = (1-x);
//     CalcFacetShape(1, ip1d, shape.Range(first_facet_dof[1], first_facet_dof[2]));
    CalcFacetShape(1, ip1d, shape);
  }
  else if (fabs(x)<1e-12)
  {
    fnr=2;
    ip1d(0) = 1-y;
//     CalcFacetShape(2, ip1d, shape.Range(first_facet_dof[2], first_facet_dof[3]));
    CalcFacetShape(2, ip1d, shape);
  }
  else if (fabs(x)>1-(1e-12))
  {
    fnr=3;
    ip1d(0) = y;
//     CalcFacetShape(3, ip1d, shape.Range(first_facet_dof[3], first_facet_dof[4]));
    CalcFacetShape(3, ip1d, shape);
  }

}
   
void FacetVolumeQuad::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
{
  SetFacet(afnr);
  shape=0.0;
  facet.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
}
   
void FacetVolumeQuad::SetFacet(int afnr) const
{
  if (fnr == afnr) return;
  
  FacetVolumeQuad * quad=const_cast<FacetVolumeQuad*>(this);
  quad->fnr = afnr;
  const EDGE * edges = ElementTopology::GetEdges (eltype);
  Array<int> fvnums(2);
  fvnums[0] = vnums[edges[fnr][0]]; 
  fvnums[1] = vnums[edges[fnr][1]]; 
   
  quad->facet.SetVertexNumbers(fvnums);
  quad->facet.SetOrder(facet_order[fnr]);
}  


void FacetVolumeQuad::ComputeNDof() 
{
  ndof = 0;
  for (int i=0; i<nf; i++)
  {
    first_facet_dof[i] = ndof;
    ndof += facet_order[i]+1;
  }
  first_facet_dof[nf] = ndof;
}




/****************************************************************************
 * FacetVolumeTet
 ****************************************************************************/

void FacetVolumeTet::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
{
  // topology: points = 0:(1,0,0), 1:(0,1,0), 2:(0,0,1), 3:(0,0,0)
  //           faces = 0:(3,1,2,-1), 1:(3,2,0,-1), 2:(3,0,1,-1), 4:(0,2,1,-1)
  IntegrationPoint ip2d;  
  double x=ip(0), y=ip(1), z=ip(2);
  int fnr=-1;
  shape=0.0;
  
//   cout << "  FacetVolumeTet::CalcShape: ip=" << ip << endl;
  if (fabs(x)<1e-6)
  {
    // (0,0,0)->(1,0); (0,1,0)->(0,1), (0,0,1)->(0,0)
    fnr=0;
    ip2d(0)=1-y-z; ip2d(1)=y;
//     CalcFacetShape(fnr, ip2d,shape.Range(first_facet_dof[fnr], first_facet_dof[fnr+1])); 
    CalcFacetShape(0, ip2d, shape);
  }
  else if (fabs(y)<1e-6)
  {
    // (0,0,0)->(1,0); (0,0,1)->(0,1); (1,0,0)->(0,0)
    fnr=1;
    ip2d(0)=1-x-z; ip2d(1)=z;
//     CalcFacetShape(fnr, ip2d,shape.Range(first_facet_dof[fnr], first_facet_dof[fnr+1])); 
    CalcFacetShape(1, ip2d, shape);
  }
  else if (fabs(z)<1e-6) 
  {
    // (0,0,0)->(1,0); (1,0,0)->(0,1); (0,1,0)->(0,0)
    fnr=2;
    ip2d(0)=1-x-y; ip2d(1)=x;
//     CalcFacetShape(fnr, ip2d,shape.Range(first_facet_dof[fnr], first_facet_dof[fnr+1])); 
    CalcFacetShape(2, ip2d, shape);
  }
  else if (fabs(1-x-y-z)<1e-6)
  {
    // (1,0,0)->(1,0); (0,0,1)->(0,1); (0,1,0)->(0,0)
    fnr=3;
    ip2d(0)=x; ip2d(1)=z;
//     CalcFacetShape(fnr, ip2d,shape.Range(first_facet_dof[fnr], first_facet_dof[fnr+1])); 
    CalcFacetShape(3, ip2d, shape);
  }
  
}
   
void FacetVolumeTet::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
{
  SetFacet(afnr);
  shape=0.0;
  facet.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
}
   
void FacetVolumeTet::SetFacet(int afnr) const
{
//   cerr << "*** Tet::SetFacet " << afnr << " called" << endl; 
  if (fnr == afnr) return;
  
  
  FacetVolumeTet * tet=const_cast<FacetVolumeTet*>(this);
  tet->fnr = afnr;
  Array<int> fvnums(3);  
  const FACE * faces = ElementTopology::GetFaces (eltype);

  fvnums[0] = vnums[faces[fnr][0]]; 
  fvnums[1] = vnums[faces[fnr][1]]; 
  fvnums[2] = vnums[faces[fnr][2]]; 
  
//   cerr << " ** setting vertex numbers << " << fvnums[0] << "," << fvnums[1] << "," << fvnums[2] << endl;
  tet->facet.SetVertexNumbers(fvnums);
  tet->facet.SetOrder(facet_order[fnr]);
}  


void FacetVolumeTet::ComputeNDof() 
{
  ndof = 0;
  for (int i=0; i<nf; i++)
  {
    first_facet_dof[i] = ndof;
    ndof += ( (facet_order[i]+1) * (facet_order[i]+2) ) / 2;
  }
  first_facet_dof[nf] = ndof;
}




/****************************************************************************
 * FacetVolumeHex
 ****************************************************************************/
// -------------------------------------------------------------------------------
void FacetVolumeHex::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
{
  cout << "error: CalcShape not implemented yet for FacetVolumeHex" << endl;
  exit(0);
}
   
void FacetVolumeHex::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
{
  shape=0.0;
  SetFacet(afnr);
  facet.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
}
   
void FacetVolumeHex::SetFacet(int afnr) const
{
  if (fnr == afnr) return;
  
  FacetVolumeHex * hex=const_cast<FacetVolumeHex*>(this);
  hex->fnr = afnr;
  Array<int> fvnums(4);  
  const FACE * faces = ElementTopology::GetFaces (eltype);

  fvnums[0] = vnums[faces[fnr][0]]; 
  fvnums[1] = vnums[faces[fnr][1]]; 
  fvnums[2] = vnums[faces[fnr][2]]; 
  fvnums[3] = vnums[faces[fnr][3]]; 
   
  hex->facet.SetVertexNumbers(fvnums);
  hex->facet.SetOrder(facet_order[fnr]); 
}  

void FacetVolumeHex::ComputeNDof()  
{
  ndof = 0;
  for (int i=0; i<nf; i++)
  {
    first_facet_dof[i] = ndof;
    ndof += (facet_order[i]+1) * (facet_order[i]+1);
  }
  first_facet_dof[nf] = ndof;
}

// const FiniteElement & FacetVolumeHex::GetFacetFE(int fnr, LocalHeap& lh) const
// {
//   FacetFacetQuad* fe = new(lh.Alloc(sizeof(FacetFacetQuad))) FacetFacetQuad();
//    
//   const FACE * faces = ElementTopology::GetFaces (eltype);
//   Array<int> fvnums(4);
//   for (int i=0; i<4; i++)
//     fvnums[i] = vnums[faces[fnr][i]]; 
//   
//   fe->SetVertexNumbers(fvnums);
//   fe->SetOrder(facet_order[fnr]);
//    
//   return *fe;
// }



/****************************************************************************
 * FacetVolumePrism
 ****************************************************************************/
// -------------------------------------------------------------------------------

void FacetVolumePrism::CalcShape (const IntegrationPoint & ip3d, FlatVector<> shape) const
{  
//   topology: points = 0:()
  IntegrationPoint ip2d;  
  double x=ip3d(0), y=ip3d(1), z=ip3d(2);

  shape = 0.0;
  if (fabs(y)<1e-12) // front
  {
    ip2d(0) = x; ip2d(1) = z;
    CalcFacetShape(2, ip2d, shape);
  }
  else if (fabs(x)<1e-12) // left
  {
    ip2d(0) = 1-y; ip2d(1)=z;
    CalcFacetShape(4, ip2d, shape);
  }
  else if (fabs(1-x-y)<1e-12) // right
  {
    ip2d(0) = (1-x+y)*0.5; ip2d(1) = z;
    CalcFacetShape(3, ip2d, shape);
  }
  else if (fabs(z)<1e-12) // botom
  {
    ip2d(0) = 1-x-y; ip2d(1) = y;
    CalcFacetShape(0, ip2d, shape);
  }
  else if (fabs(z)>1-1e-12) // top
  {
    ip2d(0) = x; ip2d(1) = y;
    CalcFacetShape(1, ip2d, shape);
  }
} 
   
void FacetVolumePrism::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
{
  shape=0.0;
  SetFacet(afnr);
  if (afnr < 2) //
    trig.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
  else // quad
    quad.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
}
   
void FacetVolumePrism::SetFacet(int afnr) const
{
  if (qnr == afnr || tnr == afnr) return;
  
  FacetVolumePrism * prism=const_cast<FacetVolumePrism*>(this);
  if (afnr < 2) // triangles
  {
    prism->tnr = afnr;
    Array<int> fvnums(3);  
    const FACE * faces = ElementTopology::GetFaces (eltype);

    fvnums[0] = vnums[faces[afnr][0]]; 
    fvnums[1] = vnums[faces[afnr][1]]; 
    fvnums[2] = vnums[faces[afnr][2]]; 
   
    prism->trig.SetVertexNumbers(fvnums);
    prism->trig.SetOrder(facet_order[tnr]);
  }
  else // quad face
  {
    prism->qnr = afnr;
    Array<int> fvnums(4);  
    const FACE * faces = ElementTopology::GetFaces (eltype);

    fvnums[0] = vnums[faces[afnr][0]]; 
    fvnums[1] = vnums[faces[afnr][1]]; 
    fvnums[2] = vnums[faces[afnr][2]]; 
    fvnums[3] = vnums[faces[afnr][3]]; 
   
    prism->quad.SetVertexNumbers(fvnums);
    prism->quad.SetOrder(facet_order[qnr]);
  }
}  


void FacetVolumePrism::ComputeNDof() 
{
  ndof = 0;
  // triangles
  for (int i=0; i<2; i++)
  {
    first_facet_dof[i] = ndof;
    ndof += ( (facet_order[i]+1) * (facet_order[i]+2) ) / 2;
  }
  //quads
  for (int i=2; i<5; i++)
  {
    first_facet_dof[i] = ndof;
    ndof += (facet_order[i]+1) * (facet_order[i]+1);
  }
  
  first_facet_dof[nf] = ndof;
}

// const FiniteElement & FacetVolumePrism::GetFacetFE(int fnr, LocalHeap& lh) const
// {
//    
//   const FACE * faces = ElementTopology::GetFaces (eltype);
//   if (fnr < 2) // triangle
//   {
//     FacetFacetTrig* fe = new(lh.Alloc(sizeof(FacetFacetTrig))) FacetFacetTrig();
//     Array<int> fvnums(3);
//     for (int i=0; i<3; i++)
//       fvnums[i] = vnums[faces[fnr][i]]; 
//     fe->SetVertexNumbers(fvnums);
//     fe->SetOrder(facet_order[fnr]);
//     return *fe;
//   }
//   else // quad
//   {
//     FacetFacetQuad* fe = new(lh.Alloc(sizeof(FacetFacetQuad))) FacetFacetQuad();
//     Array<int> fvnums(4);
//     for (int i=0; i<4; i++)
//       fvnums[i] = vnums[faces[fnr][i]]; 
//     fe->SetVertexNumbers(fvnums);
//     fe->SetOrder(facet_order[fnr]);
//     return *fe;
//   }
// }


/****************************************************************************
 * FacetVolumePyramid
 ****************************************************************************/
// -------------------------------------------------------------------------------
void FacetVolumePyramid::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
{
  cout << "error: CalcShape not implemented yet for FacetVolumePyramid" << endl;
  exit(0);
}
   
void FacetVolumePyramid::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
{
  shape=0.0;
  SetFacet(afnr);
  if (afnr < 4) // trig
    trig.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
  else
    quad.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
}
   
void FacetVolumePyramid::SetFacet(int afnr) const
{
  if (qnr == afnr || tnr == afnr) return;
  
  FacetVolumePyramid * pyramid=const_cast<FacetVolumePyramid*>(this);
  if (afnr < 4) // triangles
  {
    pyramid->tnr = afnr;
    Array<int> fvnums(3);  const FACE * faces = ElementTopology::GetFaces (eltype);

    fvnums[0] = vnums[faces[tnr][0]]; 
    fvnums[1] = vnums[faces[tnr][1]]; 
    fvnums[2] = vnums[faces[tnr][2]]; 
   
    pyramid->trig.SetVertexNumbers(fvnums);
    pyramid->trig.SetOrder(facet_order[tnr]);
  }
  else // quad face
  {
    pyramid->qnr = afnr;
    Array<int> fvnums(4);  const FACE * faces = ElementTopology::GetFaces (eltype);

    fvnums[0] = vnums[faces[qnr][0]]; 
    fvnums[1] = vnums[faces[qnr][1]]; 
    fvnums[2] = vnums[faces[qnr][2]]; 
    fvnums[3] = vnums[faces[qnr][3]]; 
   
    pyramid->quad.SetVertexNumbers(fvnums);
    pyramid->quad.SetOrder(facet_order[qnr]);
  }
}  


void FacetVolumePyramid::ComputeNDof() 
{
  ndof = 0;
  // triangles
  for (int i=0; i<4; i++)
  {
    first_facet_dof[i] = ndof;
    ndof += ( (facet_order[i]+1) * (facet_order[i]+2) ) / 2;
  }
  //quad - basis
  first_facet_dof[4] = ndof;
  ndof += (facet_order[4]+1) * (facet_order[4]+1);
  
  // final
  first_facet_dof[nf] = ndof;
}

// const FiniteElement & FacetVolumePyramid::GetFacetFE(int fnr, LocalHeap& lh) const
// {
//   const FACE * faces = ElementTopology::GetFaces (eltype);
//   if (fnr < 4) // triangle
//   {
//     FacetFacetTrig* fe = new(lh.Alloc(sizeof(FacetFacetTrig))) FacetFacetTrig();
//     Array<int> fvnums(3);
//     for (int i=0; i<3; i++)
//       fvnums[i] = vnums[faces[fnr][i]]; 
//     fe->SetVertexNumbers(fvnums);
//     fe->SetOrder(facet_order[fnr]);
//     return *fe;
//   }
//   else // facet 5: quad
//   {
//     FacetFacetQuad* fe = new(lh.Alloc(sizeof(FacetFacetQuad))) FacetFacetQuad();
//     Array<int> fvnums(4);
//     for (int i=0; i<4; i++)
//       fvnums[i] = vnums[faces[fnr][i]]; 
//     fe->SetVertexNumbers(fvnums);
//     fe->SetOrder(facet_order[fnr]);
//     return *fe;
//   }
// }



  template class FacetVolumeFiniteElement<1>;
  template class FacetVolumeFiniteElement<2>;
  template class FacetVolumeFiniteElement<3>;



} // namespace


