/*********************************************************************/
/* File:   vectorfacetfe.cpp                                         */
/* Author: A. Sinwel, (J. Schoeberl)                                 */
/* Date:   2008                                                      */
/*********************************************************************/


#include <fem.hpp>
#include <cassert>

namespace ngfem
{
  /* **************************** Facet Segm ********************************* */


  template<>
  void VectorFacetFacetFE<ET_SEGM>::CalcShape(const IntegrationPoint & ip,
                                                SliceMatrix<> shape) const
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

  template<>
  void VectorFacetFacetFE<ET_SEGM>::ComputeNDof()
  {
    order = order_inner[0];
    ndof = order+1;
  }

  /* **************************** Facet Trig ********************************* */

  template<>
  void VectorFacetFacetFE<ET_TRIG>::CalcShape(const IntegrationPoint & ip,
					  SliceMatrix<> shape) const
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

    AutoDiff<2> adxi  = lami[fav[0]]-lami[fav[2]];
    AutoDiff<2> adeta = lami[fav[1]]-lami[fav[2]];

    LegendrePolynomial::EvalScaled (p, 2*xi+eta-1, 1-eta, polx);

    Matrix<> polsy(p+1, p+1);
    DubinerJacobiPolynomials<1,0> (p, 2*eta-1, polsy);

    for (int i = 0; i <= order_inner[0]; i++)
      for (int j = 0; j <= order_inner[0]-i; j++)
	{
	  double val = polx[i] * polsy(i,j);
	  shape(ii,0) = val * adxi.DValue(0);  // lami[fav[0]].DValue(0);
	  shape(ii,1) = val * adxi.DValue(1);  // lami[fav[0]].DValue(1);
	  ii++;
	  shape(ii,0) = val * adeta.DValue(0);  // lami[fav[1]].DValue(0);
	  shape(ii,1) = val * adeta.DValue(1);  // lami[fav[1]].DValue(1);
	  ii++;
	}
  }

  template<>
  void VectorFacetFacetFE<ET_TRIG>::ComputeNDof()
  {
    order = order_inner[0];
    int p = order_inner[0];
    ndof = (p+1)*(p+2);
  }


  /* **************************** Facet Quad ********************************* */

  template<>
  void VectorFacetFacetFE<ET_QUAD>::CalcShape (const IntegrationPoint & ip,
					  SliceMatrix<> shape) const
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

    AutoDiff<2> xi  = sigma[fmax] - sigma[f1]; 
    AutoDiff<2> eta = sigma[fmax] - sigma[f2]; 
    
    int ii = 0;
    
    int n = max2(order_inner[0],order_inner[1]);
    ArrayMem<double, 20> polx(n+1), poly(n+1);

    LegendrePolynomial (n, xi.Value(), polx);
    LegendrePolynomial (n, eta.Value(), poly);
    
    for (int i = 0; i <= order_inner[0]; i++)
      for (int j = 0; j <= order_inner[1]; j++)
	{
	  double val = polx[i] * poly[j];
	  shape(ii,0) = val * xi.DValue(0);
	  shape(ii,1) = val * xi.DValue(1);
	  ii++;
	  shape(ii,0) = val * eta.DValue(0);
	  shape(ii,1) = val * eta.DValue(1);
	  ii++;
	}
  }

  template<>
  void VectorFacetFacetFE<ET_QUAD>::ComputeNDof()
  {
    order = max2( order_inner[0], order_inner[1] );
    ndof = 2 * (order_inner[0]+1) * (order_inner[1]+1);
  }


  /* **************************** Volume Trig ********************************* */

  template<>
  void VectorFacetVolumeFE<ET_TRIG> ::
  CalcShape ( const IntegrationPoint & ip, int fanr, SliceMatrix<> shape ) const
  {
    for (int i = 0; i < ndof; i++)
      shape(i, 0) = shape(i, 1) = 0;

    int first = first_facet_dof[fanr];

    AutoDiff<2> x(ip(0), 0), y(ip(1),1);

    const EDGE * edges = ElementTopology :: GetEdges (ET_TRIG);

    int fav[2] = { edges[fanr][0], edges[fanr][1] };
    int j1 = 0, j2 = 1;
    if(vnums[fav[j1]] > vnums[fav[j2]]) swap(j1,j2); 

    AutoDiff<2> lami[3] = {x, y, 1-x-y};  

    int p = facet_order[fanr][0];

    AutoDiff<2> xi = lami[fav[j1]] - lami[fav[j2]];
    LegendrePolynomial (p, xi.Value(), 
                        SBLambda([&](int nr, double val)
                                 {
                                   shape(first+nr,0) = val * xi.DValue(0);
                                   shape(first+nr,1) = val * xi.DValue(1);
                                 }));
  }

  template<>
  void VectorFacetVolumeFE<ET_TRIG> ::
  CalcExtraShape ( const IntegrationPoint & ip, int fanr, FlatMatrixFixWidth<2> xshape ) const
  {
    xshape = 0.0;

    AutoDiff<2> x(ip(0), 0), y(ip(1),1);

    const EDGE * edges = ElementTopology :: GetEdges (ET_TRIG);

    int  fav[2] = {edges[fanr][0], edges[fanr][1] };
    int j1 = 0; 
    int j2 = 1; 
    if(vnums[fav[j1]] > vnums[fav[j2]]) swap(j1,j2); 

    AutoDiff<2> lami[3] = {x, y, 1-x-y};  

    int p = facet_order[fanr][0];
    
    ArrayMem< double, 10> polx(p+2);
    // int ii = first;

    AutoDiff<2> xi = lami[fav[j1]] - lami[fav[j2]];

    LegendrePolynomial (p+1, xi.Value(), polx);
    
    double val = polx[p+1];
    xshape(0,0) = val * xi.DValue(0);
    xshape(0,1) = val * xi.DValue(1);
  }

  template<>
  void VectorFacetVolumeFE<ET_TRIG> :: ComputeNDof()
  {
    ndof = 0;
    for ( int i = 0; i < 3; i++ )
      {
	first_facet_dof[i] = ndof;
	ndof += facet_order[i][0] + 1;
      }
    first_facet_dof[3] = ndof;
  }

  /* **************************** Volume Quad ********************************* */

  template<>
  void VectorFacetVolumeFE<ET_QUAD> ::
  CalcShape ( const IntegrationPoint & ip, int fanr, SliceMatrix<> shape ) const
  {
    shape = 0.0;

    AutoDiff<2> x(ip(0), 0), y(ip(1),1);

    const EDGE * faces = ElementTopology :: GetEdges (ET_QUAD);

    int  fav[2] = {faces[fanr][0], faces[fanr][1] };
    int j1 = 0; 
    int j2 = 1; 
    if(vnums[fav[j1]] > vnums[fav[j2]]) swap(j1,j2);  // fmax > f2 > f1; 

    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  

    int p = facet_order[fanr][0];
    ArrayMem< double, 10> polx(p+1);
    int ii = first_facet_dof[fanr];

    AutoDiff<2> xi = sigma[fav[j1]] - sigma[fav[j2]];

    LegendrePolynomial (p, xi.Value(), polx);
    for (int i = 0; i <= facet_order[fanr][0]; i++)
      {
	double val = polx[i];
	shape(ii,0) = val * xi.DValue(0);
	shape(ii,1) = val * xi.DValue(1);
	ii++;
      }
  }


  template<>
  void VectorFacetVolumeFE<ET_QUAD> ::
  CalcExtraShape ( const IntegrationPoint & ip, int fanr, FlatMatrixFixWidth<2> xshape ) const
  {
    xshape = 0.0;

    AutoDiff<2> x(ip(0), 0), y(ip(1),1);

    const EDGE * faces = ElementTopology :: GetEdges (ET_QUAD);

    int  fav[2] = {faces[fanr][0], faces[fanr][1] };
    int j1 = 0; 
    int j2 = 1; 
    if(vnums[fav[j1]] > vnums[fav[j2]]) swap(j1,j2);  // fmax > f2 > f1; 

    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  

    int p = facet_order[fanr][0];
    ArrayMem< double, 10> polx(p+2);

    AutoDiff<2> xi = sigma[fav[j1]] - sigma[fav[j2]];

    LegendrePolynomial (p+1, xi.Value(), polx);
    double val = polx[p+1];
    xshape(0,0) = val * xi.DValue(0);
    xshape(0,1) = val * xi.DValue(1);
  }

  template<>
  void VectorFacetVolumeFE<ET_QUAD> :: ComputeNDof()
  {
    ndof = 0;
    for ( int i = 0; i < 4; i++ )
      {
	first_facet_dof[i] = ndof;
	ndof += facet_order[i][0] + 1;
      }
    first_facet_dof[4] = ndof;
  }

  
  /* **************************** Volume Tet ********************************* */

  template<>
  void VectorFacetVolumeFE<ET_TET> ::
  CalcShape ( const IntegrationPoint & ip, int fanr, SliceMatrix<> shape ) const
  {
    for (int i = 0; i < ndof; i++)
      shape(i,0) = shape(i,1) = shape(i,2) = 0;

    AutoDiff<3> x(ip(0), 0), y(ip(1),1), z(ip(2),2);
    AutoDiff<3> lami[4] = { x, y, z, 1-x-y-z };

    INT<4> fav = ET_T::GetFaceSort (fanr, vnums);

    AutoDiff<3> adxi = lami[fav[0]]-lami[fav[2]];
    AutoDiff<3> adeta = lami[fav[1]]-lami[fav[2]];
    double xi = lami[fav[0]].Value();
    double eta = lami[fav[1]].Value();

    int p = facet_order[fanr][0];
    ArrayMem< double, 10> polx(p+1), poly(p+1);
    Matrix<> polsy(p+1, p+1);
    int ii = first_facet_dof[fanr];

    // ScaledLegendrePolynomial (p, 2*xi+eta-1, 1-eta, polx);
    LegendrePolynomial::EvalScaled (p, 2*xi+eta-1, 1-eta, polx);

    DubinerJacobiPolynomials<1,0> (p, 2*eta-1, polsy);

    for (int i = 0; i <= facet_order[fanr][0]; i++)
      for (int j = 0; j <= facet_order[fanr][0]-i; j++)
	{
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
  }

  template<>
  void VectorFacetVolumeFE<ET_TET> ::
  CalcExtraShape ( const IntegrationPoint & ip, int fanr, FlatMatrixFixWidth<3> xshape ) const
  {
    xshape = 0.0;

    AutoDiff<3> x(ip(0), 0), y(ip(1),1), z(ip(2),2);
    AutoDiff<3> lami[4] = { x, y, z, 1-x-y-z };

    INT<4> fav = ET_trait<ET_TET>::GetFaceSort (fanr, vnums);

    AutoDiff<3> adxi = lami[fav[0]]-lami[fav[2]];
    AutoDiff<3> adeta = lami[fav[1]]-lami[fav[2]];
    double xi = lami[fav[0]].Value();
    double eta = lami[fav[1]].Value();

    int p = facet_order[fanr][0];
    ArrayMem< double, 10> polx(p+2), poly(p+2);
    Matrix<> polsy(p+2, p+2);
    int ii = 0;

    ScaledLegendrePolynomial (p+1, 2*xi+eta-1, 1-eta, polx);
    DubinerJacobiPolynomials<1,0> (p+1, 2*eta-1, polsy);

    for (int i = 0; i <= p+1; i++)
      {
	int j = p+1-i;
	double val = polx[i] * polsy(i, j);
	xshape(ii,0) = val * adxi.DValue(0);
	xshape(ii,1) = val * adxi.DValue(1);
	xshape(ii,2) = val * adxi.DValue(2);
	ii++;
	xshape(ii,0) = val * adeta.DValue(0);
	xshape(ii,1) = val * adeta.DValue(1);
	xshape(ii,2) = val * adeta.DValue(2);
	ii++;
      }
  }


  template<>
  void VectorFacetVolumeFE<ET_TET> :: ComputeNDof()
  {
    ndof = 0;
    for (int i = 0; i < 4; i++)
      {
	first_facet_dof[i] = ndof;
	ndof += (facet_order[i][0]+1) * (facet_order[i][0]+2);
      }
    first_facet_dof[4] = ndof;
  }

  /* **************************** Volume Prism ********************************* */

  template<>
  void VectorFacetVolumeFE<ET_PRISM> ::
  CalcShape ( const IntegrationPoint & ip, int fanr, SliceMatrix<> shape ) const
  {
    AutoDiff<3> x(ip(0), 0), y(ip(1),1), z(ip(2),2);

    AutoDiff<3> lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    AutoDiff<3> muz[6]  = { 1-z, 1-z, 1-z, z, z, z }; 

    AutoDiff<3> sigma[6];
    for (int i = 0; i < 6; i++) sigma[i] = lami[i] + muz[i];

    shape = 0.0;

    // trig face shapes
    if (fanr < 2)
      {
	int p = facet_order[fanr][0];

	INT<4> fav = ET_trait<ET_PRISM>::GetFaceSort (fanr, vnums);

	AutoDiff<3> adxi = lami[fav[0]]-lami[fav[2]];
	AutoDiff<3> adeta = lami[fav[1]]-lami[fav[2]];
	double xi = lami[fav[0]].Value();
	double eta = lami[fav[1]].Value();
	
	ArrayMem< double, 10> polx(p+1), poly(p+1);
	Matrix<> polsy(p+1, p+1);
	
	ScaledLegendrePolynomial (p, 2*xi+eta-1, 1-eta, polx);
	DubinerJacobiPolynomials<1,0> (p, 2*eta-1, polsy);
	
	int ii = first_facet_dof[fanr];
	for (int i = 0; i <= facet_order[fanr][0]; i++)
	  for (int j = 0; j <= facet_order[fanr][0]-i; j++)
	    {
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
      }
    else
      // quad faces
      {
	int p = facet_order[fanr][0];
	 
	INT<4> f = ET_trait<ET_PRISM>::GetFaceSort (fanr, vnums);	  	
	AutoDiff<3> xi  = sigma[f[0]] - sigma[f[1]]; 
	AutoDiff<3> zeta = sigma[f[0]] - sigma[f[3]];

	ArrayMem< double, 10> polx(p+1), poly(p+1);
	LegendrePolynomial (p, xi.Value(), polx);
	LegendrePolynomial (p, zeta.Value(), poly);

	int ii = first_facet_dof[fanr];
	for (int i = 0; i <= p; i++)
	  for (int j = 0; j <= p; j++)
	    {
	      double val = polx[i] * poly[j];
	      shape(ii,0) = val * xi.DValue(0);
	      shape(ii,1) = val * xi.DValue(1);
	      shape(ii,2) = val * xi.DValue(2);
	      ii++;
	      shape(ii,0) = val * zeta.DValue(0);
	      shape(ii,1) = val * zeta.DValue(1);
	      shape(ii,2) = val * zeta.DValue(2);
	      ii++;
	    }	
      }
  }

  template<>
  int VectorFacetVolumeFE<ET_PRISM> ::
  GetNExtraShapes ( int fanr) const
  {
    if (fanr < 2) //trig shape
      return 2*(facet_order[fanr][0]+2);
    else //quad shape
      return 2*(2*facet_order[fanr][0]+3);
  }

  template<>
  void VectorFacetVolumeFE<ET_PRISM> ::
  CalcExtraShape ( const IntegrationPoint & ip, int fanr, FlatMatrixFixWidth<3> xshape ) const
  {
    AutoDiff<3> x(ip(0), 0), y(ip(1),1), z(ip(2),2);

    AutoDiff<3> lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    AutoDiff<3> muz[6]  = { 1-z, 1-z, 1-z, z, z, z }; 

    AutoDiff<3> sigma[6];
    for (int i = 0; i < 6; i++) sigma[i] = lami[i] + muz[i];

    xshape = 0.0;

    // trig face shapes
    if (fanr < 2)
      {
	int p = facet_order[fanr][0];
	INT<4> fav = ET_trait<ET_PRISM>::GetFaceSort (fanr, vnums);

	AutoDiff<3> adxi = lami[fav[0]]-lami[fav[2]];
	AutoDiff<3> adeta = lami[fav[1]]-lami[fav[2]];
	double xi = lami[fav[0]].Value();
	double eta = lami[fav[1]].Value();
	
	ArrayMem< double, 10> polx(p+2), poly(p+2);
	Matrix<> polsy(p+2, p+2);
	
	ScaledLegendrePolynomial (p+1, 2*xi+eta-1, 1-eta, polx);
	DubinerJacobiPolynomials<1,0> (p+1, 2*eta-1, polsy);
	
	int ii = 0;
	for (int i = 0; i <= p+1; i++)
	{
	  int j= p+1-i;
	  double val = polx[i] * polsy(i, j);
	  xshape(ii,0) = val * adxi.DValue(0);
	  xshape(ii,1) = val * adxi.DValue(1);
	  xshape(ii,2) = val * adxi.DValue(2);
	  ii++;
	  xshape(ii,0) = val * adeta.DValue(0);
	  xshape(ii,1) = val * adeta.DValue(1);
	  xshape(ii,2) = val * adeta.DValue(2);
	  ii++;
	}
      }
    else
      // quad faces
      {
	int p = facet_order[fanr][0];
	 
	INT<4> f = ET_trait<ET_PRISM>::GetFaceSort (fanr, vnums);	  	
	AutoDiff<3> xi  = sigma[f[0]] - sigma[f[1]]; 
	AutoDiff<3> zeta = sigma[f[0]] - sigma[f[3]];
	
	ArrayMem< double, 10> polx(p+2), poly(p+2);
	LegendrePolynomial (p+1, xi.Value(), polx);
	LegendrePolynomial (p+1, zeta.Value(), poly);

	int ii = 0;
	for (int i = 0; i <= p+1; i++)
	{
	  for (int j = (i==p+1)?0:p+1; j <= p+1; j++)
	    {
	      double val = polx[i] * poly[j];
	      xshape(ii,0) = val * xi.DValue(0);
	      xshape(ii,1) = val * xi.DValue(1);
	      xshape(ii,2) = val * xi.DValue(2);
	      ii++;
	      xshape(ii,0) = val * zeta.DValue(0);
	      xshape(ii,1) = val * zeta.DValue(1);
	      xshape(ii,2) = val * zeta.DValue(2);
	      ii++;
	    }
	}	
      }
  }

  template<>
  void VectorFacetVolumeFE<ET_PRISM>::ComputeNDof()
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
  
    first_facet_dof[5] = ndof;
  }


  /* **************************** Volume Hex ********************************* */

  template<>
  void VectorFacetVolumeFE<ET_HEX> ::
  CalcShape ( const IntegrationPoint & ip, int facet, SliceMatrix<> shape ) const
  {
    throw Exception("VectorFacetVolumeHex::CalcShape not implemented!");
  }

  template<>
  void VectorFacetVolumeFE<ET_HEX>::ComputeNDof()
  {
    ndof = 0;
    for (int i=0; i<6; i++)
      {
	first_facet_dof[i] = ndof;
	ndof += 2* (facet_order[i][0]+1) * (facet_order[i][0]+1);
      }
    first_facet_dof[6] = ndof;
  }



  /* **************************** Volume Pyramid ********************************* */

  template<>
  void VectorFacetVolumeFE<ET_PYRAMID> ::
  CalcShape ( const IntegrationPoint & ip, int facet, SliceMatrix<> shape ) const
  {
    throw Exception("VectorFacetVolumePyramid::CalcShape: not implemented!");
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

  template<>
  void VectorFacetVolumeFE<ET_PYRAMID>::ComputeNDof()
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
    first_facet_dof[4] = ndof;
  }



  static RegisterBilinearFormIntegrator<RobinEdgeIntegrator<3 /* , VectorFacetFacetFiniteElement<2> */ >  > initrvf3 ("robinvectorfacet", 3, 1);
  static RegisterBilinearFormIntegrator<RobinEdgeIntegrator<2 /* , VectorFacetFacetFiniteElement<1> */ >  > initrvf2 ("robinvectorfacet", 2, 1);
  static RegisterLinearFormIntegrator<NeumannEdgeIntegrator<3 /*, VectorFacetFacetFiniteElement<2> */ >  > initnvf3 ("neumannvectorfacet", 3, 1);
  static RegisterLinearFormIntegrator<NeumannEdgeIntegrator<2 /*, VectorFacetFacetFiniteElement<1> */ >  > initnvf2 ("neumannvectorfacet", 2, 1);
}


