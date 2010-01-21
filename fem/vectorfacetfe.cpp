/*********************************************************************/
/* File:   vectorfacetfe.cpp                                         */
/* Author: A. Sinwel, (J. Schoeberl)                                 */
/* Date:   2008                                                      */
/*********************************************************************/


#include <fem.hpp>

namespace ngfem {


  using namespace ngfem;

  VectorFacetFacetFiniteElement ::
  VectorFacetFacetFiniteElement (int dim, ELEMENT_TYPE aeltype) :
    FiniteElement ( dim, aeltype, -1, -1 )
  {
    for (int i=0; i<8; i++)
      vnums[i] = -1; 
    order_inner = INT<2> (-1, -1);
  }

  void   VectorFacetFacetFiniteElement ::
  SetVertexNumbers (FlatArray<int> & avnums)
  {
    for ( int i = 0; i < avnums.Size(); i++ )
      vnums[i] = avnums[i];

  }

  void   VectorFacetFacetFiniteElement ::
  SetOrder (int aorder)
  {
    order = aorder;
    order_inner = INT<2> (aorder, aorder);
    ComputeNDof();
  }
  
  void   VectorFacetFacetFiniteElement ::
  SetOrder (INT<2> oi)
  {
    order = max2 (oi[0], oi[1]);
    order_inner = oi;
    ComputeNDof();
  }




  VectorFacetFacetSegm :: VectorFacetFacetSegm (int aorder) :
    VectorFacetFacetFiniteElement (1, ET_SEGM )
  {
    order = aorder;
    order_inner = INT<2> (aorder, aorder);
    ComputeNDof();
  }

  void VectorFacetFacetSegm :: ComputeNDof()
  {
    order = order_inner[0];
    ndof = order + 1;
  }


  /// compute shape
  // FlatMatrix ( ndof * 1 )
  void VectorFacetFacetSegm :: CalcShape (const IntegrationPoint & ip, 
					  FlatMatrix<> shape) const
  {
    AutoDiff<1> x (ip(0),0);
    ArrayMem<double, 10>  polx(order_inner[0]+1);
    // orient
    if ( vnums[0] > vnums[1])
      x = 1-x;
    
    int ii = 0;

    LegendrePolynomial (order_inner[0], 2*x.Value()-1, polx);
    for ( int i = 0; i <= order_inner[0]; i++ )
      shape(ii++,0) = 2 * polx[i] * x.DValue(0);
  }


  VectorFacetFacetTrig :: VectorFacetFacetTrig (int aorder) :
    VectorFacetFacetFiniteElement (2, ET_TRIG)
  {
    order = aorder;
    order_inner = INT<2> ( aorder, aorder);
    ComputeNDof();
  }

  void VectorFacetFacetTrig :: ComputeNDof()
  {
    order = order_inner[0];
    int p = order_inner[0];
    ndof = (p+1)*(p+2);
  }



  /// compute shape
  void VectorFacetFacetTrig :: CalcShape (const IntegrationPoint & ip, 
					  FlatMatrix<> shape) const
  {
    AutoDiff<2> x (ip(0), 0);
    AutoDiff<2> y (ip(1), 1);

    int p = order_inner[0];
    ArrayMem< double, 10> polx(p+1), poly(p+1);
    int ii = 0;

    AutoDiff<2> lami[3] = { x, y, 1-x-y };
    int fav[3] = { 0, 1, 2};
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	


    double xi = lami[fav[0]].Value();
    double eta = lami[fav[1]].Value();

    ScaledLegendrePolynomial (p, 2*xi+eta-1, 1-eta, polx);
    LegendrePolynomial (p, 2*eta-1, poly);

    // *testout << "surface trig, orderinner = " << order_inner[0] << endl;

    for (int i = 0; i <= order_inner[0]; i++)
      for (int j = 0; j <= order_inner[0]-i; j++)
	{
	  double val = polx[i] * poly[j];
	  shape(ii,0) = val * lami[fav[0]].DValue(0);
	  shape(ii,1) = val * lami[fav[0]].DValue(1);
	  ii++;
	  shape(ii,0) = val * lami[fav[1]].DValue(0);
	  shape(ii,1) = val * lami[fav[1]].DValue(1);
	  ii++;
	}
    // *testout << "surface trig " << endl << shape << endl;
  }

  /**
     High order quadrilateral finite element
  */
  VectorFacetFacetQuad :: VectorFacetFacetQuad (int aorder) :
    VectorFacetFacetFiniteElement (2, ET_QUAD)
  {
    order = aorder;
    order_inner = INT<2> (aorder, aorder);
    ComputeNDof();
  }

  void VectorFacetFacetQuad :: ComputeNDof()
  {
    order = max2( order_inner[0], order_inner[2] );
    ndof = 2 * (order_inner[0]+1) * (order_inner[1]+1);
  }

  /// compute shape
  void VectorFacetFacetQuad :: CalcShape (const IntegrationPoint & ip, 
					  FlatMatrix<> shape) const
  {
    AutoDiff<2> x (ip(0), 0);
    AutoDiff<2> y (ip(1), 1);

    // orient: copied from h1
    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    int fmax = 0; 
    for (int j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax]) fmax = j;
    int f1 = (fmax+3)%4; 
    int f2 = (fmax+1)%4; 
    if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f1 > f2; 

    x = sigma[fmax] - sigma[f1]; 
    y = sigma[fmax] - sigma[f2]; 
    
    int ii = 0;
    
    int n = max(order_inner[0],order_inner[1]);
    ArrayMem<double, 20> polx(n+1), poly(n+1);

    LegendrePolynomial (n, x.Value(), polx);
    LegendrePolynomial (n, y.Value(), poly);
    
    for (int i = 0; i <= order_inner[0]; i++)
      for (int j = 0; j <= order_inner[1]; j++)
	{
	  double val = polx[i] * poly[j];
	  shape(ii,0) = val * x.DValue(0);
	  shape(ii,1) = val * y.DValue(1);
	  ii++;
	  shape(ii,0) = val * x.DValue(0);
	  shape(ii,1) = val * y.DValue(1);
	  ii++;
	}
  }


  ///+++++++++++++++++++++++++


  VectorFacetVolumeFiniteElement :: 
  VectorFacetVolumeFiniteElement (int adim, ELEMENT_TYPE aeltype) :
    FiniteElement(adim, aeltype, -1, -1)
  {
    for ( int i = 0; i < 8; i++ )
      vnums[i] = -1;
    for ( int i = 0; i < 6; i++ )
      facet_order[i] = INT<2> (-1, -1);
    for ( int i = 0; i < 7; i++ )
      first_facet_dof[i] = 0;
    switch (eltype)
      {
      case ET_TRIG:
	nv=3; nf=3; break; 
      case ET_QUAD:
	nv=4; nf=4; break;
      case ET_TET:
	nv=4; nf=4; break;
      case ET_PRISM:
	nv=6; nf=5; break;
      case ET_HEX:
	nv=8; nf=6; break;
      default:
	cerr << "*** error in FacetVolumeFiniteElement::FacetVolumeFiniteElement(int, ELEMENT_TYPE):unknown element type" << endl;
	exit(0);
      }
  }

  

  void VectorFacetVolumeFiniteElement :: 
  SetVertexNumbers (FlatArray<int> & avnums)
  {
    for ( int i = 0; i < avnums.Size(); i++ )
      vnums[i] = avnums[i];
  }

  void VectorFacetVolumeFiniteElement :: 
  SetOrder(int ao)
  {
    order = ao;
    for ( int i = 0; i < 6; i++ )
      facet_order[i] = INT<2> (ao, ao);
    ComputeNDof();
  }

  void VectorFacetVolumeFiniteElement :: 
  SetOrder(FlatArray<INT<2> > & ao)
  {
    for ( int i = 0; i < ao.Size(); i++ )
      {
        order = max3 ( order, ao[i][0], ao[i][1] );
	facet_order[i] = ao[i];
      }
    ComputeNDof();
  }

  void VectorFacetVolumeFiniteElement :: 
  SetOrder(FlatArray<int> & ao)
  {
    for ( int i = 0; i < ao.Size(); i++ )
      {
        order = max2 ( order, ao[i] );
	facet_order[i] = INT<2> (ao[i], ao[i]);
      }
    ComputeNDof();
  }

  void VectorFacetVolumeFiniteElement :: 
  GetFacetDofNrs(int afnr, Array<int>& fdnums) const
  {
    fdnums.SetSize(0);
    int first = first_facet_dof[afnr];
    int next = first_facet_dof[afnr+1];
    for ( int i = first; i < next; i++ )
      fdnums.Append(i); 
  }
  
  void VectorFacetVolumeFiniteElement :: 
  GetVertexNumbers(Array<int>& avnums) const
  {
    avnums.SetSize(nv);
    for ( int i = 0; i < nv; i++ )
      avnums[i] = vnums[i];
  }
  
  void VectorFacetVolumeFiniteElement :: GetFacetOrders(Array<INT<2> >& forder) const
  {
    forder.SetSize(6);
    for ( int i = 0; i < 6; i++ )
      forder[i] = facet_order[i];
  }
  

//   // Works only for transformation from 2d to 3d

//   void VectorFacetVolumeFiniteElement :: 
//   TransformFacetToVolumeShape ( int fanr, FlatMatrix<> shape2d, 
// 				FlatMatrix<> shape ) const
//   {
//     int first = first_facet_dof[fanr];
//     int last = first_facet_dof[fanr+1];
//     // the parts of shape/shape2d where the shapes for this face are stored
//     FlatMatrix <> facetshape = FlatMatrix<> (last - first, 3, &shape(first, 0) );
//     FlatMatrix <> facetshape2d = FlatMatrix<> (last - first, 2, &shape2d(0,0) );
//     Mat<2,3> jacobian;

//     const FACE * faces = ElementTopology :: GetFaces ( eltype );
//     const POINT3D * points = ElementTopology :: GetVertices (eltype);
//     const ELEMENT_TYPE facettype = ElementTopology :: GetFacetType (eltype, fanr);
//     int nvert = ElementTopology :: GetNVertices (facettype);

//     int  fav[4] = {faces[fanr][0], faces[fanr][1], faces[fanr][2], faces[fanr][3] };
//     AutoDiff<2> x (0.0, 0);
//     AutoDiff<2> y (0.0, 1);
//     AutoDiff<2> lami[3] = {x, y, 1-x-y};  
//     AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  

//     int jmax = nvert - 1; 
// //     for (int j = 0; j < nvert; j++)
// //       if (vnums[fav[j]] > vnums[fav[jmax]]) jmax = j;
//     int j1 = (jmax+1)%nvert; 
//     int j2 = (jmax+nvert-1)%nvert; 
// //     if(vnums[fav[j1]] > vnums[fav[j2]]) swap(j1,j2);  // fmax > f1 > f2; 
//     if ( nvert == 4 )
//       {
// 	x = sigma[jmax]-sigma[j1]; 
// 	y = sigma[jmax]-sigma[j2]; 
//       }
//     else
//       {
// 	x = lami[j1];
// 	y = lami[j2];
//       }

//     for ( int i = 0; i < 3; i++ )
//       {
// 	jacobian (0,i) = (points[fav[j1]][i] - points[fav[jmax]][i]);
// 	jacobian (1,i) = (points[fav[j2]][i] - points[fav[jmax]][i]);
//       }
//     facetshape = facetshape2d * jacobian;

    
// //     cout << "ok, my facetshape2d is " << endl << facetshape2d << endl
// // 	 << "and my jacobian is " << endl << jacobian << endl 
// //  	 << "the facetshape " << endl << facetshape << endl;
// // 	 << "and the shape " << endl << shape << endl;
    
//   }

  void VectorFacetVolumeTrig :: 
  CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const
  {
    double x=ip(0), y=ip(1);
    int fanr; 

    shape = 0.0;
    if (fabs(y) < 1e-4) 
      {
	fanr = 0;
      }
    else if (fabs(x) < 1e-4)
      {
	fanr = 1;
      }
    else 
      {
	fanr = 2;
      }

    CalcShape ( ip, fanr, shape);
  }

  void VectorFacetVolumeTrig ::
  CalcShape ( const IntegrationPoint & ip, int fanr, FlatMatrix<> shape ) const
  {
    shape = 0.0;
    int first = first_facet_dof[fanr];

    AutoDiff<2> x(ip(0), 0), y(ip(1),1);

    const EDGE * edges = ElementTopology :: GetEdges ( eltype );

    int  fav[2] = {edges[fanr][0], edges[fanr][1] };
    int j1 = 0; 
    int j2 = 1; 
    if(vnums[fav[j1]] > vnums[fav[j2]]) swap(j1,j2); 

    AutoDiff<2> lami[3] = {x, y, 1-x-y};  

    int p = facet_order[fanr][0];
    ArrayMem< double, 10> polx(p+1);
    int ii = first;

    AutoDiff<2> xi = lami[fav[j1]] - lami[fav[j2]];

    LegendrePolynomial (p, xi.Value(), polx);
    for (int i = 0; i <= facet_order[fanr][0]; i++)
      {
	double val = polx[i];
	shape(ii,0) = val * xi.DValue(0);
	shape(ii,1) = val * xi.DValue(1);
	ii++;
      }
  }




//   void VectorFacetVolumeTrig :: CalcFacetShape (int afnr, const IntegrationPoint & ip, 
// 						FlatMatrix<> shape1d) const
//   {
//     SetFacet(afnr);
//     shape1d = 0.0;
//     facet.CalcShape(ip, shape1d);
//     ;
//   }

//   void VectorFacetVolumeTrig :: SetFacet(int afnr) const
//   {
//     if (fnr == afnr) return;
  
//     VectorFacetVolumeTrig * trig=const_cast<VectorFacetVolumeTrig*>(this);
//     trig->fnr = afnr;
//     const EDGE * edges = ElementTopology::GetEdges (eltype);
//     Array<int> fvnums(2);
//     fvnums[0] = vnums[edges[fnr][0]]; 
//     fvnums[1] = vnums[edges[fnr][1]]; 
    
//     trig->facet.SetVertexNumbers(fvnums);
//     trig->facet.SetOrder(facet_order[fnr]);
//   }  
  
    
  void VectorFacetVolumeTrig :: ComputeNDof()
  {
    ndof = 0;
    for ( int i = 0; i < nf; i++ )
      {
	first_facet_dof[i] = ndof;
	ndof += facet_order[i][0] + 1;
      }
    first_facet_dof[nf] = ndof;
  }

  
//   void VectorFacetVolumeTrig :: 
//   TransformFacetToVolumeShape ( int fanr, FlatMatrix<> shape1d, 
// 				FlatMatrix<> shape ) const
//   {
//     int first = first_facet_dof[fanr];
//     int last = first_facet_dof[fanr+1];

//     // the parts in shape/shape1d where the shapes for this face/order are stored
//     FlatMatrix <> facetshape = FlatMatrix<> (last - first, 2, &shape(first, 0) );
//     FlatMatrix <> facetshape1d = FlatMatrix<> (last - first, 1, &shape1d(0, 0) );

//     Mat<1,2> jacobian;
//     const EDGE * edges = ElementTopology :: GetEdges (ET_TRIG);
//     const POINT3D * points = ElementTopology :: GetVertices (ET_TRIG);

//     INT<2> fav (edges[fanr][0], edges[fanr][1]);
// //     if ( vnums[fav[0]] > vnums[fav[1]] )
// //       swap(fav[0], fav[1]);

//     AutoDiff<1> x (0.0, 0);
//     AutoDiff<1> lami[3] = { x, 1-x };

//     jacobian = 0.0;

//     for ( int i = 0; i < 2; i++ )
//       for ( int j = 0; j < 1; j++ )
// 	for ( int k = 0; k < 2; k++ )
// 	  {
// 	    jacobian(j, k) += points[fav[i]][k] * lami[i].DValue(j);
// 	  }

//     facetshape = facetshape1d * jacobian;

//   }
  



  //--------------------------------------------------

  void VectorFacetVolumeQuad :: 
  CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const
  {
    // topology: points = 0:(0,0), 1:(1,0), 2:(1,1), 3:(0,1)
    //           edges = 0:(0,1), 1:(2,3), 2:(3,0), 4:(1,2)
    IntegrationPoint ip1d;  
    double x=ip(0), y=ip(1);
    int fnr=-1;
    
    Matrix<double> shape1d (order+1, 1);

    shape = 0.0;
    if (fabs(y)<1e-12)
      {
	fnr=0;
      }
    else if (fabs(y)>1-(1e-12))
      {
	fnr=1;
      }
    else if (fabs(x)<1e-12)
      {
	fnr=2;
      }
    else if (fabs(x)>1-(1e-12))
      {
	fnr=3;
      }
    CalcShape( ip, fnr, shape);
  }

  void VectorFacetVolumeQuad ::
  CalcShape ( const IntegrationPoint & ip, int fanr, FlatMatrix<> shape ) const
  {

    shape = 0.0;
    int first = first_facet_dof[fanr];
    int last = first_facet_dof[fanr+1];
    // the parts of shape where the shapes for this face are stored
    FlatMatrix <> facetshape = FlatMatrix<> (last - first, 3, &shape(first, 0) );

    AutoDiff<2> x(ip(0), 0), y(ip(1),1);

    const EDGE * faces = ElementTopology :: GetEdges ( eltype );
    // const POINT3D * points = ElementTopology :: GetVertices (eltype);
    // const ELEMENT_TYPE facettype = ElementTopology :: GetFacetType (eltype, fanr);
    // int nvert = ElementTopology :: GetNVertices (facettype);

    int  fav[2] = {faces[fanr][0], faces[fanr][1] };
    int j1 = 0; 
    int j2 = 1; 
    if(vnums[fav[j1]] > vnums[fav[j2]]) swap(j1,j2);  // fmax > f2 > f1; 

    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  

    int p = facet_order[fanr][0];
    ArrayMem< double, 10> polx(p+1);
    int ii = first;

    AutoDiff<2> xi = sigma[fav[j1]] - sigma[fav[j2]];

    LegendrePolynomial (p, 2*xi.Value()-1, polx);
    for (int i = 0; i <= facet_order[fanr][0]; i++)
      {
	double val = polx[i];
	shape(ii,0) = val * xi.DValue(0);
	shape(ii,1) = val * xi.DValue(1);
	ii++;
      }
  }

//   void VectorFacetVolumeQuad :: CalcFacetShape (int afnr, const IntegrationPoint & ip, 
// 						FlatMatrix<> shape1d) const
//   {
//     SetFacet( afnr );
//     shape1d = 0.0;
//     facet.CalcShape (ip, shape1d );
//   }

//   void VectorFacetVolumeQuad :: SetFacet(int afnr) const
//   {
//     if (fnr == afnr) return;
  
//     VectorFacetVolumeQuad * quad=const_cast<VectorFacetVolumeQuad*>(this);
//     quad->fnr = afnr;
//     const EDGE * edges = ElementTopology::GetEdges (eltype);
//     Array<int> fvnums(2);
//     fvnums[0] = vnums[edges[fnr][0]]; 
//     fvnums[1] = vnums[edges[fnr][1]]; 
    
//     quad->facet.SetVertexNumbers(fvnums);
//     quad->facet.SetOrder(facet_order[fnr]);
//   }  
  
    
  void VectorFacetVolumeQuad :: ComputeNDof()
  {
    ndof = 0;
    for ( int i = 0; i < nf; i++ )
      {
	first_facet_dof[i] = ndof;
	ndof += facet_order[i][0] + 1;
      }
    first_facet_dof[nf] = ndof;
  }

  
//   void VectorFacetVolumeQuad :: 
//   TransformFacetToVolumeShape ( int fanr, FlatMatrix<> shape1d, 
// 				FlatMatrix<> shape ) const
//   {
//     int first = first_facet_dof[fanr];
//     int last = first_facet_dof[fanr+1];

//     // the parts of shape/shape1d in which the shape functions for this face are stored
//     FlatMatrix <> facetshape = FlatMatrix<> (last - first, 2, &shape(first, 0) );
//     FlatMatrix <> facetshape1d = FlatMatrix<> (last - first, 1, &shape(0, 0) );

//     Mat<1,2> jacobian;
//     const EDGE * edges = ElementTopology :: GetEdges (ET_QUAD);
//     const POINT3D * points = ElementTopology :: GetVertices (ET_QUAD);

//     INT<2> fav (edges[fanr][0], edges[fanr][1]);
// //     if ( vnums[fav[0]] > vnums[fav[1]] )
// //       swap(fav[0], fav[1]);

//     AutoDiff<1> x (0.0, 0);
//     AutoDiff<1> lami[3] = { x, 1-x };

//     for ( int i = 0; i < 2; i++ )
//       for ( int j = 0; j < 1; j++ )
// 	for ( int k = 0; k < 2; k++ )
// 	  {
// 	    jacobian(j, k) += points[fav[i]][k] * lami[i].DValue(j);
// 	  }

//     facetshape = facetshape1d * jacobian;
//   }


  // #######################################

  void VectorFacetVolumeTet::
  CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const
  {
//     cout << "TET:: CalcShape" << endl;
//     cout << "integration point " << ip << endl;
    // topology: points = 0:(1,0,0), 1:(0,1,0), 2:(0,0,1), 3:(0,0,0)
    //           faces = 0:(3,1,2,-1), 1:(3,2,0,-1), 2:(3,0,1,-1), 4:(0,2,1,-1)
    
    double x=ip(0), y=ip(1), z=ip(2);
    int fanr;

    shape=0.0;


    if (fabs(x)<1e-4)
      {
	// (0,0,0)->(1,0); (0,1,0)->(0,1), (0,0,1)->(0,0)
	fanr=0;
      }
    else if (fabs(y)<1e-4)
      {
	// (0,0,0)->(1,0); (0,0,1)->(0,1); (1,0,0)->(0,0)
	fanr=1;
      }
    else if (fabs(z)<1e-4) 
      {
	// (0,0,0)->(1,0); (1,0,0)->(0,1); (0,1,0)->(0,0)
	fanr=2;
      }
    else if (fabs(1-x-y-z)<1e-1)
      {
	// (1,0,0)->(1,0); (0,0,1)->(0,1); (0,1,0)->(0,0)
	fanr=3;
      }
    else
      return;

    CalcShape ( ip, fanr, shape );

//     CalcFacetShape(fanr, ip2d, shape2d);
//     TransformFacetToVolumeShape (fanr, shape2d, shape);

  }   

  void VectorFacetVolumeTet ::
  CalcShape ( const IntegrationPoint & ip, int fanr, FlatMatrix<> shape ) const
  {
    /*
    IntegrationPoint ip2d;
    double x=ip(0), y=ip(1), z=ip(2);

    int maxmatsize = (order+1)*(order+2);
    Matrix<> shape2d (maxmatsize, 2);

    shape=0.0;

    shape2d = 0.0;

    switch (fanr)
      {
      case 0:
	ip2d(0)=1-y-z; 
	ip2d(1)=y;
	break;
      case 1:
	ip2d(0)=1-x-z; ip2d(1)=z;
	break;
      case 2:
	ip2d(0)=1-x-y; ip2d(1)=x;
	break;      
      case 3:
	ip2d(0)=x; ip2d(1)=z;
	break;
      default: 
	return;
      }

    CalcFacetShape(fanr, ip2d, shape2d);
    TransformFacetToVolumeShape ( fanr, shape2d, shape );
    */
    
    shape = 0.0;
    int first = first_facet_dof[fanr];
    int last = first_facet_dof[fanr+1];
    // the parts of shape where the shapes for this face are stored
    FlatMatrix <> facetshape = FlatMatrix<> (last - first, 3, &shape(first, 0) );

    AutoDiff<3> x(ip(0), 0), y(ip(1),1), z(ip(2),2);

    const FACE * faces = ElementTopology :: GetFaces ( eltype );
    // const POINT3D * points = ElementTopology :: GetVertices (eltype);
    const ELEMENT_TYPE facettype = ElementTopology :: GetFacetType (eltype, fanr);
    int nvert = ElementTopology :: GetNVertices (facettype);

    int  fav[4] = {faces[fanr][0], faces[fanr][1], faces[fanr][2], faces[fanr][3] };
    int jmax = 0; 
    for (int j = 1; j < nvert; j++)
      if (vnums[fav[j]] > vnums[fav[jmax]]) jmax = j;
    int j1 = (jmax+1)%nvert; 
    int j2 = (jmax+nvert-1)%nvert; 
    if(vnums[fav[j1]] > vnums[fav[j2]]) swap(j1,j2);  // fmax > f2 > f1; 

    AutoDiff<3> lami[4] = {x, y, z, 1-x-y-z };

    int p = facet_order[fanr][0];
    ArrayMem< double, 10> polx(p+1), poly(p+1);
    Matrix<> polsy(p+1, p+1);
    int ii = first;

    AutoDiff<3> adxi  = lami[fav[j1]]-lami[fav[jmax]];
    AutoDiff<3> adeta = lami[fav[j2]]-lami[fav[jmax]];

    double xi  = lami[fav[j1]].Value();
    double eta = lami[fav[j2]].Value();

    ScaledLegendrePolynomial (p, 2*xi+eta-1, 1-eta, polx);
    // LegendrePolynomial (p, 2*eta-1, poly);
    DubinerJacobiPolynomials (p, 2*eta-1, 1, 0, polsy);

    /*
    *testout << "facet_order = " << facet_order[fanr][0] << endl;
    *testout << "fanr = " << fanr << endl;
    *testout << "xi, = " << xi << ", eta = " << eta << endl;

    *testout << "polx = " << polx << endl;
    *testout << "poly = " << poly << endl;
    */

    for (int i = 0; i <= facet_order[fanr][0]; i++)
      for (int j = 0; j <= facet_order[fanr][0]-i; j++)
	{
	  // double val = polx[i] * poly[j];
	  double val = polx[i] * polsy(i, j);
	  shape(ii,0) = val * adxi.DValue(0);
	  shape(ii,1) = val * adxi.DValue(1);
	  shape(ii,2) = val * adxi.DValue(2);
	  ii++;
	  shape(ii,0) = val * adeta.DValue(0);
	  shape(ii,1) = val * adeta.DValue(1);
	  shape(ii,2) = val * adeta.DValue(2);
	  ii++;
	}
    
    // *testout << "shape = " << endl << shape << endl;
  }


//   void VectorFacetVolumeTet::
//   CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatMatrix<> shape2d) const
//   {
//     SetFacet( afnr );
//     shape2d = 0.0;
//     facet.CalcShape (ip, shape2d );
//   }
   
//   void VectorFacetVolumeTet::
//   SetFacet(int afnr) const
//   {
//     if (fnr == afnr) return;
    
    
//     VectorFacetVolumeTet * tet=const_cast<VectorFacetVolumeTet*>(this);
//     tet->fnr = afnr;
//     Array<int> fvnums(3);  
//     const FACE * faces = ElementTopology::GetFaces (eltype);
    
//     fvnums[0] = vnums[faces[fnr][0]]; 
//     fvnums[1] = vnums[faces[fnr][1]]; 
//     fvnums[2] = vnums[faces[fnr][2]]; 
  
//     //   cerr << " ** setting vertex numbers << " << fvnums[0] << "," << fvnums[1] << "," << fvnums[2] << endl;
//     tet->facet.SetVertexNumbers(fvnums);
//     tet->facet.SetOrder(facet_order[fnr]);
//   }  


  void VectorFacetVolumeTet :: ComputeNDof() 
  {
    ndof = 0;
    for (int i=0; i<nf; i++)
      {
	first_facet_dof[i] = ndof;
	ndof += ( (facet_order[i][0]+1) * (facet_order[i][0]+2) );
      }
    first_facet_dof[nf] = ndof;
  }


  // --------------------------------------------

  void VectorFacetVolumeHex::
  CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const
  {
    //   topology: points = 0:()
    int fanr = -1;
    double x=ip(0), y=ip(1), z=ip(2);
    
    shape = 0.0;

    if (fabs(z)<1e-6)
      {
	// (0,0) -> (0,0,0), (1,0) -> (0,1,0), (1,1) -> (1,1,0), (0,1) -> (1,0,0)
	fanr=0;
      }
    else if (fabs(z)>1-1e-6)
      {
	// (0,0) -> (0,0,1), (1,0) -> (1,0,0), (1,1) -> (1,1,1), (0,1) -> (0,1,1)
	fanr=1;
      }
    else if (fabs(y)<1e-6) 
      {
	// (0,0) -> (0,0,0), (1,0) -> (1,0,0), (1,1) -> (1,0,1), (0,1) -> (0,0,1)
	fanr=2;
      }
    else if (fabs(x)> 1-1e-6)
      {
	// (0,0) -> (1,0,0), (1,0) -> (1,1,0), (1,1) -> (1,1,1), (0,1) -> (1,0,1)
	fanr=3;
      }
    else if (fabs(y)> 1-1e-6)
      {
	// (0,0) -> (1,1,0), (1,0) -> (0,1,0), (1,1) -> (0,1,1), (0,1) -> (1,1,1)
	fanr=4;
      }
    else if (fabs(x)<1e-6)
      {
	// (0,0) -> (0,1,0), (1,0) -> (0,0,0), (1,1) -> (0,0,1), (0,1) -> (0,1,1)
	fanr=5;
      }

    CalcShape ( ip, fanr, shape );
    ;
  }

  void VectorFacetVolumeHex ::
  CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrix<> shape ) const
  {
    cout << "VectorFacetVolumeHex::CalcShape not implemented!" << endl;
    exit(0);
  }

   
//   void VectorFacetVolumeHex::
//   CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatMatrix<> shape2d) const
//   {
//     shape2d  = 0.0;
//     SetFacet(afnr);
//     facet.CalcShape(ip, shape2d);
//     ;
//   }
   
//   void VectorFacetVolumeHex::SetFacet(int afnr) const
//   {
//     if (fnr == afnr) return;
    
//     VectorFacetVolumeHex * hex=const_cast<VectorFacetVolumeHex*>(this);
//     hex->fnr = afnr;
//     Array<int> fvnums(4);  
//     const FACE * faces = ElementTopology::GetFaces (eltype);
    
//     fvnums[0] = vnums[faces[fnr][0]]; 
//     fvnums[1] = vnums[faces[fnr][1]]; 
//     fvnums[2] = vnums[faces[fnr][2]]; 
//     fvnums[3] = vnums[faces[fnr][3]]; 
    
//     hex->facet.SetVertexNumbers(fvnums);
//     hex->facet.SetOrder(facet_order[fnr]); 
//   }  
  
  void VectorFacetVolumeHex::ComputeNDof()  
  {
    ndof = 0;
    for (int i=0; i<nf; i++)
      {
	first_facet_dof[i] = ndof;
	ndof += 2* (facet_order[i][0]+1) * (facet_order[i][0]+1);
      }
    first_facet_dof[nf] = ndof;
  }


  // -------------------------------------------------------------------------------

  void VectorFacetVolumePrism::
  CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const
  {  
    int fanr = -1;
    double x=ip(0), y=ip(1), z=ip(2);
    
    shape = 0.0;

    if (fabs(y)<1e-12) // front
      {
	fanr = 2;
      }
    else if (fabs(x)<1e-12) // left
      {
	fanr = 4;
      }
    else if (fabs(1-x-y)<1e-12) // right
      {
	fanr = 3;
      }
    else if (fabs(z)<1e-12) // botom
      {
	fanr = 0;
      }
    else if (fabs(z)>1-1e-12) // top
      {
	fanr = 1;
      }

    CalcShape ( ip, fanr, shape );
    ;
  } 
   
  void VectorFacetVolumePrism ::
  CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrix<> shape ) const
  {
    cout << "VectorFacetVolumePrism::CalcShape not implemented" << endl;
    exit(0);
  }

//   void VectorFacetVolumePrism::
//   CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatMatrix<> shape2d) const
//   {
//     shape2d=0.0;
//     SetFacet(afnr);
//     if (afnr < 2) //
//       trig.CalcShape(ip, shape2d);
//     else // quad
//       quad.CalcShape(ip, shape2d);
//   }
   
//   void VectorFacetVolumePrism::SetFacet(int afnr) const
//   {
//     if (qnr == afnr || tnr == afnr) return;
  
//     VectorFacetVolumePrism * prism=const_cast<VectorFacetVolumePrism*>(this);
//     if (afnr < 2) // triangles
//       {
// 	prism->tnr = afnr;
// 	Array<int> fvnums(3);  
// 	const FACE * faces = ElementTopology::GetFaces (eltype);

// 	fvnums[0] = vnums[faces[afnr][0]]; 
// 	fvnums[1] = vnums[faces[afnr][1]]; 
// 	fvnums[2] = vnums[faces[afnr][2]]; 
   
// 	prism->trig.SetVertexNumbers(fvnums);
// 	prism->trig.SetOrder(facet_order[tnr]);
//       }
//     else // quad face
//       {
// 	prism->qnr = afnr;
// 	Array<int> fvnums(4);  
// 	const FACE * faces = ElementTopology::GetFaces (eltype);

// 	fvnums[0] = vnums[faces[afnr][0]]; 
// 	fvnums[1] = vnums[faces[afnr][1]]; 
// 	fvnums[2] = vnums[faces[afnr][2]]; 
// 	fvnums[3] = vnums[faces[afnr][3]]; 
   
// 	prism->quad.SetVertexNumbers(fvnums);
// 	prism->quad.SetOrder(facet_order[qnr]);
//       }
//   }  


  void VectorFacetVolumePrism::ComputeNDof() 
  {
    ndof = 0;
    // triangles
    for (int i=0; i<2; i++)
      {
	first_facet_dof[i] = ndof;
	ndof += ( (facet_order[i][0]+1) * (facet_order[i][0]+2) );
      }
    //quads
    for (int i=2; i<5; i++)
      {
	first_facet_dof[i] = ndof;
	ndof += 2 * (facet_order[i][0]+1) * (facet_order[i][1]+1);
      }
  
    first_facet_dof[nf] = ndof;
  }


  // -------------------------------------------------------------------------------
  void VectorFacetVolumePyramid::
  CalcShape (const IntegrationPoint & ip, FlatMatrix<> shape) const
  {
    ;
    cout << "error in VectorFacetVolumePyramid::CalcShape: not implemented!" << endl;
    exit(0);
  }
   
  void VectorFacetVolumePyramid ::
  CalcShape ( const IntegrationPoint & ip, int facet, FlatMatrix<> shape ) const
  {
    cout << "error in VectorFacetVolumePyramid::CalcShape: not implemented!" << endl;
    exit(0);
  }

//   void VectorFacetVolumePyramid::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatMatrix<> shape) const
//   {
//     ;
//   }
   
//   void VectorFacetVolumePyramid::SetFacet(int afnr) const
//   {
//     if (qnr == afnr || tnr == afnr) return;
  
//     VectorFacetVolumePyramid * pyramid=const_cast<VectorFacetVolumePyramid*>(this);
//     if (afnr < 4) // triangles
//       {
// 	pyramid->tnr = afnr;
// 	Array<int> fvnums(3);  const FACE * faces = ElementTopology::GetFaces (eltype);

// 	fvnums[0] = vnums[faces[tnr][0]]; 
// 	fvnums[1] = vnums[faces[tnr][1]]; 
// 	fvnums[2] = vnums[faces[tnr][2]]; 
   
// 	pyramid->trig.SetVertexNumbers(fvnums);
// 	pyramid->trig.SetOrder(facet_order[tnr]);
//       }
//     else // quad face
//       {
// 	pyramid->qnr = afnr;
// 	Array<int> fvnums(4);  const FACE * faces = ElementTopology::GetFaces (eltype);

// 	fvnums[0] = vnums[faces[qnr][0]]; 
// 	fvnums[1] = vnums[faces[qnr][1]]; 
// 	fvnums[2] = vnums[faces[qnr][2]]; 
// 	fvnums[3] = vnums[faces[qnr][3]]; 
   
// 	pyramid->quad.SetVertexNumbers(fvnums);
// 	pyramid->quad.SetOrder(facet_order[qnr]);
//       }
//   }  


  void VectorFacetVolumePyramid::ComputeNDof() 
  {
    ndof = 0;
    // triangles
    for (int i=0; i<4; i++)
      {
	first_facet_dof[i] = ndof;
	ndof += ( (facet_order[i][0]+1) * (facet_order[i][0]+2) );
      }
    //quad - basis
    first_facet_dof[4] = ndof;
    ndof += 2 * (facet_order[4][0]+1) * (facet_order[4][0]+1);
  
    // final
    first_facet_dof[nf] = ndof;
  }

}



namespace ngfem
{

  using namespace ngfem;


  namespace vectorfacetint {

    class Init
    { 
    public: 
      Init ();

    };
    
    Init::Init()
    {
             
      GetIntegrators().AddBFIntegrator ("massvectorfacet", 3, 1,
					MassEdgeIntegrator<3, VectorFacetVolumeFiniteElement>::Create);
      GetIntegrators().AddBFIntegrator ("massvectorfacet", 2, 1,
					MassEdgeIntegrator<2, VectorFacetVolumeFiniteElement>::Create);
      GetIntegrators().AddBFIntegrator ("robinvectorfacet", 3, 1,
					RobinEdgeIntegrator<3, VectorFacetFacetFiniteElement>::Create);
      GetIntegrators().AddBFIntegrator ("robinvectorfacet", 2, 1,
					RobinEdgeIntegrator<2, VectorFacetFacetFiniteElement>::Create);      
      

      GetIntegrators().AddLFIntegrator ("sourcevectorfacet", 3, 3,
					SourceEdgeIntegrator<3, VectorFacetVolumeFiniteElement>::Create);
      GetIntegrators().AddLFIntegrator ("sourcevectorfacet", 2, 2,
					SourceEdgeIntegrator<2, VectorFacetVolumeFiniteElement>::Create);
      GetIntegrators().AddLFIntegrator ("neumannvectorfacet", 3, 3,
					NeumannEdgeIntegrator<3, VectorFacetFacetFiniteElement>::Create);
      GetIntegrators().AddLFIntegrator ("neumannvectorfacet", 2, 2,
					NeumannEdgeIntegrator<2, VectorFacetFacetFiniteElement>::Create);
    }

    Init init;
  }
}


