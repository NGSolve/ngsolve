/*********************************************************************/
/* File:   vectorfacetfe.cpp                                         */
/* Author: A. Sinwel, (J. Schoeberl)                                 */
/* Date:   2008                                                      */
/*********************************************************************/


// #include <fem.hpp>
#include "vectorfacetfe.hpp"
#include "hcurl_equations.hpp"
#include <thcurlfe_impl.hpp>
#include "hcurlfe_utils.hpp"
#include <cassert>

namespace ngfem
{

  template <ELEMENT_TYPE ET> template<typename Tx, typename TFA>  
  void VectorFacetVolumeFE<ET>::
  T_CalcShape (Tx hx[DIM], int fnr, TFA & shape) const
  {
    throw ExceptionNOSIMD("VectorFacetVolume::T_CalcShape missing element "+ToString(ET));
  }

  template <> template<typename Tx, typename TFA>  
  void VectorFacetVolumeFE<ET_TRIG>::
  T_CalcShape (Tx hx[DIM], int fanr, TFA & shape) const
  {
    if (fanr == -1) throw Exception("vector-facet element evaluated not at BND");    
    Tx lami[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };  

    int first = first_facet_dof[fanr];
    int p = facet_order[fanr][0];

    INT<2> e = GetVertexOrientedEdge(fanr);
    Tx xi = lami[e[0]] - lami[e[1]];
    
    LegendrePolynomial (p, xi, 
                        SBLambda([&](int nr, auto val)
                                 {
                                   shape[first+nr] = uDv (val, xi);
                                 }));
  }


  template<> template<typename Tx, typename TFA>  
  void VectorFacetVolumeFE<ET_TET> ::
  T_CalcShape (Tx hx[DIM], int fanr, TFA & shape ) const
  {
    if (fanr == -1) throw Exception("vector-facet element evaluated not at BND");
    Tx x = hx[0], y = hx[1], z = hx[2];
    Tx lami[4] = { x, y, z, 1-x-y-z };

    INT<4> fav = ET_T::GetFaceSort (fanr, vnums);

    Tx adxi = lami[fav[0]]-lami[fav[2]];
    Tx adeta = lami[fav[1]]-lami[fav[2]];

    size_t ii = first_facet_dof[fanr];
    DubinerBasis::Eval(facet_order[fanr][0], lami[fav[1]].Value(), lami[fav[0]].Value(),
                       SBLambda([shape,&ii,adxi,adeta] (size_t nr, auto val)
                                {
                                  shape[ii] = uDv(Tx(val), adxi); ii++;
                                  shape[ii] = uDv(Tx(val), adeta); ii++;
                                }));
  }

  
  
  template <ELEMENT_TYPE ET>
  void VectorFacetVolumeFE<ET>::
  CalcMappedShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                   BareSliceMatrix<SIMD<double>> shapes) const
  {
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    shapes.AddSize(DIM*this->GetNDof(), mir.Size()) = 0.0;
    for (size_t i = 0; i < mir.Size(); i++)
      {
        Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
        T_CalcShape (&adp(0), mir[i].IP().FacetNr(),
                     SBLambda ([&] (size_t j, auto s)
                               {
                                 auto shape = s.Value();
                                 for (size_t k = 0; k < DIM; k++)
                                   shapes(j*DIM+k, i) = shape(k);
                               }));
      }
  }

  template <ELEMENT_TYPE ET>
  void VectorFacetVolumeFE<ET>::
  Evaluate (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs,
            BareSliceMatrix<SIMD<double>> values) const
  {
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
        Vec<DIM,SIMD<double>> sum(0.0);
        T_CalcShape (&adp(0), mir[i].IP().FacetNr(),
                     SBLambda ([&] (int j, auto s)
                               {
                                 auto shape = s.Value();
                                 double coef = coefs(j);
                                 Iterate<DIM> ( [&] (auto ii) {
                                     sum(ii.value) += coef * shape(ii.value);
                                   });
                               }));
        for (size_t k = 0; k < DIM; k++)
          values(k,i) = sum(k).Data();
      }
  }
  
  template <ELEMENT_TYPE ET>
  void VectorFacetVolumeFE<ET>::
  AddTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
            BareSliceVector<> coefs) const
  {
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
        Vec<DIM,SIMD<double>> vali = values.Col(i);
        
        T_CalcShape (&adp(0), mir[i].IP().FacetNr(),
                     SBLambda ([vali,coefs] (size_t j, auto s)
                               {
                                 /*
                                 auto shape = s.Value();
                                 SIMD<double> sum = 0.0;
                                 for (int k = 0; k < DIM; k++)
                                   sum += shape(k) * values(k,i);
                                 coefs(j) += HSum(sum);
                                 */
                                 coefs(j) += HSum(InnerProduct(s.Value(), vali));
                               }));
      }
  }

  
  /* **************************** Facet Segm ********************************* */

  template <ELEMENT_TYPE ET> template<typename Tx, typename TFA>  
  void VectorFacetFacetFE<ET>::
  T_CalcShape (TIP<DIM,Tx> tip, TFA & shape) const
  {
    throw ExceptionNOSIMD("VectorFacetFacet::T_CalcShape missing"+ToString(ET));
  }


  
  template<ELEMENT_TYPE ET>
  void VectorFacetFacetFE<ET>::CalcShape(const IntegrationPoint & ip,
                                         BareSliceMatrix<> shape) const
  {
    TIP<DIM,AutoDiff<DIM>> tip = ip;
    T_CalcShape (tip,
                 SBLambda([shape] (size_t i, auto val)
                          {
                            shape.Row(i) = val.Value();
                          }));
  }
  
  template<> template <typename Tx, typename TFA>
  void VectorFacetFacetFE<ET_SEGM>::T_CalcShape(TIP<DIM,Tx> tip,
                                                TFA & shape) const
  {
    // AutoDiff<1> x (ip(0),0);
    /*
    Tx x = tip.x;
    ArrayMem<double, 10>  polx(order_inner[0]+1);
    // orient
    if ( vnums[0] > vnums[1])
      x = 1-x;
    
    int ii = 0;

    LegendrePolynomial (order_inner[0], 2*x.Value()-1, polx);
    for ( int i = 0; i <= order_inner[0]; i++ )
      shape(ii++,0) = 2 * polx[i] * x.DValue(0);
    */
    Tx x = tip.x;
    if ( vnums[0] > vnums[1]) x = 1-x;
    Tx sx = 2*x-1;
    LegendrePolynomial (order_inner[0], sx,
                        SBLambda([&] (size_t i, Tx val)
                                 {
                                   shape[i] = uDv(val, sx);
                                 }));
  }

  template<>
  void VectorFacetFacetFE<ET_SEGM>::ComputeNDof()
  {
    order = order_inner[0];
    ndof = order+1;
  }

  /* **************************** Facet Trig ********************************* */

  template<> template <typename Tx, typename TFA>
  void VectorFacetFacetFE<ET_TRIG>::T_CalcShape(TIP<DIM,Tx> tip, 
                                                TFA &  shape) const
  {
    /*
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
    */
    auto x = tip.x;
    auto y = tip.y;

    int ii = 0;

    Tx lami[3] = { x, y, 1-x-y };
    int fav[3] = { 0, 1, 2};
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	

    AutoDiff<2> adxi  = lami[fav[0]]-lami[fav[2]];
    AutoDiff<2> adeta = lami[fav[1]]-lami[fav[2]];

    DubinerBasis::Eval(order_inner[0], lami[fav[1]].Value(), lami[fav[0]].Value(),
                       SBLambda([&] (size_t nr, Tx val)
                                {
                                  shape[ii] = uDv(Tx(val), adxi); ii++;
                                  shape[ii] = uDv(Tx(val), adeta); ii++;
                                }));
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
                                               BareSliceMatrix<> shape) const
  {
    AutoDiff<2> x (ip(0), 0);
    AutoDiff<2> y (ip(1), 1);

    // orient: copied from h1
    AutoDiff<2> sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    // int fmax = 0; 
    // for (int j = 1; j < 4; j++)
    //   if (vnums[j] > vnums[fmax]) fmax = j;
    // int f1 = (fmax+3)%4; 
    // int f2 = (fmax+1)%4; 
    // if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f1 > f2; 

    // AutoDiff<2> xi  = sigma[fmax] - sigma[f1]; 
    // AutoDiff<2> eta = sigma[fmax] - sigma[f2]; 
    
    INT<4> f = GetFaceSort (0, vnums); 
    AutoDiff<2> xi  = sigma[f[0]] - sigma[f[1]]; 
    AutoDiff<2> eta = sigma[f[0]] - sigma[f[3]]; 

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
  CalcShape ( const IntegrationPoint & ip, int fanr, BareSliceMatrix<> shape ) const
  {
    for (int i = 0; i < ndof; i++)
      shape(i, 0) = shape(i, 1) = 0;

    TIP<DIM,AutoDiff<DIM>> tip = ip;
    T_CalcShape(&tip.x, fanr,
		SBLambda([&](size_t j, auto s) {
		    shape.Row(j) = s.Value();
		  })
		);
  }

  template<>
  template<typename MIP, typename TFA>  
  void VectorFacetVolumeFE<ET_TRIG> :: CalcDualShape2 (const MIP & mip, int fnr, TFA & shape) const
  {
    auto & ip = mip.IP();
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
    T x = ip(0), y = ip(1);
    T lam[3] = { x, y, 1-x-y };
    Vec<2,T> pnts[3] = { { 1, 0 }, { 0, 1 } , { 0, 0 } };

    int first = first_facet_dof[fnr];
    
    if (ip.VB() == BND)
      {
        int p = facet_order[fnr][0];
        
        INT<2> e = GetVertexOrientedEdge(fnr);
        T xi = lam[e[1]]-lam[e[0]];
        Vec<2,T> tauref = pnts[e[1]] - pnts[e[0]];
        auto tau = mip.GetJacobian()*tauref;
	tau /= mip.GetMeasure();

        LegendrePolynomial::Eval
          (p, xi,
           SBLambda([&] (size_t nr, T val)
                    {
                      shape[first+nr] = val*tau;
                    }));
    
      }
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
  CalcShape ( const IntegrationPoint & ip, int fanr, BareSliceMatrix<> shape ) const
  {
    shape.AddSize(ndof, 2) = 0.0;

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

  template<>
  template<typename MIP, typename TFA>  
  void VectorFacetVolumeFE<ET_QUAD> :: CalcDualShape2 (const MIP & mip, int fnr, TFA & shape) const
  {
    throw Exception("calcdualshape2 not implemented for ET_QUAD VectorFacetVolumeFE ");
  }

  
  /* **************************** Volume Tet ********************************* */

  template<>
  void VectorFacetVolumeFE<ET_TET> ::
  CalcShape ( const IntegrationPoint & ip, int fanr, BareSliceMatrix<> shape ) const
  {
    for (int i = 0; i < ndof; i++)
      shape(i,0) = shape(i,1) = shape(i,2) = 0;

    AutoDiff<3> x(ip(0), 0), y(ip(1),1), z(ip(2),2);
    AutoDiff<3> lami[4] = { x, y, z, 1-x-y-z };

    INT<4> fav = ET_T::GetFaceSort (fanr, vnums);

    AutoDiff<3> adxi = lami[fav[0]]-lami[fav[2]];
    AutoDiff<3> adeta = lami[fav[1]]-lami[fav[2]];
    //double xi = lami[fav[0]].Value();
    //double eta = lami[fav[1]].Value();

    //int p = facet_order[fanr][0];
    int ii = first_facet_dof[fanr];

    /*
    cout << "VERY SLOW" << endl;
    ArrayMem< double, 10> polx(p+1), poly(p+1);
    Matrix<> polsy(p+1, p+1);

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
    cout << "old: " << shape.Rows(first_facet_dof[fanr], first_facet_dof[fanr+1]);
    */
    //typedef AutoDiff<3> Tx;
    ii = first_facet_dof[fanr];    
    DubinerBasis::Eval(facet_order[fanr][0], lami[fav[1]].Value(), lami[fav[0]].Value(),
                       SBLambda([&] (size_t nr, double val)
                                {
                                  // shape[ii] = uDv(Tx(val), adxi); ii++;
                                  // shape[ii] = uDv(Tx(val), adeta); ii++;
                                  shape(ii,0) = val * adxi.DValue(0);
                                  shape(ii,1) = val * adxi.DValue(1);
                                  shape(ii,2) = val * adxi.DValue(2);
                                  ii++;
                                  shape(ii,0) = val * adeta.DValue(0);
                                  shape(ii,1) = val * adeta.DValue(1);
                                  shape(ii,2) = val * adeta.DValue(2);
                                  ii++;
                                }));
    // cout << "new: " << shape.Rows(first_facet_dof[fanr], first_facet_dof[fanr+1]);      
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

  template<>
  template<typename MIP, typename TFA>  
  void VectorFacetVolumeFE<ET_TET> :: CalcDualShape2 (const MIP & mip, int fnr, TFA & shape) const
  {
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;        
    auto & ip = mip.IP();
    T x = ip(0), y = ip(1), z = ip(2);

    T lam[4] = { x, y, z, 1-x-y-z };
    Vec<3> pnts[4] = { { 1, 0, 0 }, { 0, 1, 0 } , { 0, 0, 1 }, { 0, 0, 0 } };

    int ii = first_facet_dof[fnr];

    // INT<4> fav = GetFaceSort (fnr, vnums);
    INT<4> fav = GetVertexOrientedFace(fnr);
    Vec<3> adxi = pnts[fav[0]] - pnts[fav[2]];
    Vec<3> adeta = pnts[fav[1]] - pnts[fav[2]];
    T xi = lam[fav[0]];
    T eta = lam[fav[1]];
    
    // Matrix<T> tauhat(3,2);
    Mat<3,2,T> tauhat; 
    tauhat.Col(0) = adxi;//Vec<3,T>(adxi.DValue(0),adxi.DValue(1),adxi.DValue(2));
    tauhat.Col(1) = adeta;//Vec<3,T>(adeta.DValue(0),adeta.DValue(1),adeta.DValue(2));
    // Matrix<T> tau(3,2);
    Mat<3,2,T> tau;
    // (dShat/dS) *F * tauhat
    tau = mip.GetJacobian() * tauhat;
    tau /= mip.GetMeasure();

    DubinerBasis::Eval(order, xi, eta,
		       SBLambda([&] (size_t nr, auto val)
				{
				  shape[ii++] = tau * Vec<2,T>(val, 0);
				  shape[ii++] = tau * Vec<2,T>(0, val);
				}));
  }

  /* **************************** Volume Prism ********************************* */

  template<>
  void VectorFacetVolumeFE<ET_PRISM> ::
  CalcShape ( const IntegrationPoint & ip, int fanr, BareSliceMatrix<> shape ) const
  {
    AutoDiff<3> x(ip(0), 0), y(ip(1),1), z(ip(2),2);

    AutoDiff<3> lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    AutoDiff<3> muz[6]  = { 1-z, 1-z, 1-z, z, z, z }; 

    AutoDiff<3> sigma[6];
    for (int i = 0; i < 6; i++) sigma[i] = lami[i] + muz[i];

    shape.AddSize(ndof, 3) = 0.0;

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

  template<>
  template<typename MIP, typename TFA>  
  void VectorFacetVolumeFE<ET_PRISM> :: CalcDualShape2 (const MIP & mip, int fnr, TFA & shape) const
  {
    throw Exception("calcdualshape2 not implemented for ET_PRISM VectorFacetVolumeFE ");
  }


  /* **************************** Volume Hex ********************************* */

  template<>
  void VectorFacetVolumeFE<ET_HEX> ::
  CalcShape ( const IntegrationPoint & ip, int fanr, BareSliceMatrix<> shape ) const
  {
    AutoDiff<3> x(ip(0), 0), y(ip(1),1), z(ip(2),2);

    // vertex numbering on HEX:
    // { 0, 0, 0 },
    // { 1, 0, 0 },
    // { 1, 1, 0 },
    // { 0, 1, 0 },
    // { 0, 0, 1 },
    // { 1, 0, 1 },
    // { 1, 1, 1 },
    // { 0, 1, 1 }

    AutoDiff<3> mux[8]  = { 1-x,   x,   x, 1-x, 1-x,   x,   x, 1-x };  
    AutoDiff<3> muy[8]  = { 1-y, 1-y,   y,   y, 1-y, 1-y,   y,   y }; 
    AutoDiff<3> muz[8]  = { 1-z, 1-z, 1-z, 1-z,   z,   z,   z,   z };    

    AutoDiff<3> sigma[8];
    for (int i = 0; i < 8; i++) sigma[i] = mux[i] + muy[i] + muz[i];
    
    shape.AddSize(ndof, 3) = 0.0;
    {
      int p = facet_order[fanr][0];
       
      INT<4> f = ET_trait<ET_HEX>::GetFaceSort (fanr, vnums);     
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

  template<>
  template<typename MIP, typename TFA>  
  void VectorFacetVolumeFE<ET_HEX> :: CalcDualShape2 (const MIP & mip, int fnr, TFA & shape) const
  {
    throw Exception("calcdualshape2 not implemented for ET_HEX VectorFacetVolumeFE ");
  }



  /* **************************** Volume Pyramid ********************************* */

  template<>
  void VectorFacetVolumeFE<ET_PYRAMID> ::
  CalcShape ( const IntegrationPoint & ip, int facet, BareSliceMatrix<> shape ) const
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

  template<>
  template<typename MIP, typename TFA>  
  void VectorFacetVolumeFE<ET_PYRAMID> :: CalcDualShape2 (const MIP & mip, int fnr, TFA & shape) const
  {
    throw Exception("calcdualshape2 not implemented for ET_PYRAMID VectorFacetVolumeFE ");
  }

  template <ELEMENT_TYPE ET>
  void VectorFacetVolumeFE<ET> :: CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const
  {
    shape.AddSize(ndof, bmip.DimSpace()) = 0.0;
    Switch<4-DIM>
      (bmip.DimSpace()-DIM,[this,&bmip,shape](auto CODIM)
       {
         auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);
         this->CalcDualShape2 (mip, mip.IP().FacetNr(),
                               SBLambda ([&] (size_t i, auto val)
                                         {
                                           shape.Row(i) = val;
                                         }));
       });
  }

  template <ELEMENT_TYPE ET>
  void VectorFacetVolumeFE<ET> :: CalcDualShape (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> shapes) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,shapes](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         shapes.AddSize(ndof*DIMSPACE, mir.Size()) = 0.0;
         for (size_t i = 0; i < mir.Size(); i++)
           this->CalcDualShape2 (mir[i], mir[i].IP().FacetNr(),
                                 SBLambda ([shapes,i,DIMSPACE] (size_t j, auto val)
                                           {
                                             for (size_t k = 0; k < DIMSPACE; k++)
                                               shapes(j*DIMSPACE+k, i) = val(k);
                                           }));
       });
  }

  template <ELEMENT_TYPE ET>
  void VectorFacetVolumeFE<ET> :: EvaluateDual (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIMSPACE,SIMD<double>> sum (SIMD<double>(0.0));
             this -> CalcDualShape2 (mir[i], mir[i].IP().FacetNr(),
                                     SBLambda([&sum, coefs] (size_t j, auto val)
                                              {
                                                sum += coefs(j) * val;
                                              }));
             for (size_t k = 0; k < DIMSPACE; k++)
               values(k, i) = sum(k);
           }
       });
  }

  template <ELEMENT_TYPE ET>
  void VectorFacetVolumeFE<ET> :: AddDualTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values, BareSliceVector<double> coefs) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIMSPACE,SIMD<double>> value = values.Col(i);
             this -> CalcDualShape2 (mir[i], mir[i].IP().FacetNr(),
                                     SBLambda([value, coefs] (size_t j, auto val)
                                              {
                                                coefs(j) += HSum(InnerProduct(value,val));
                                              }));
           }
       });
  }
  
  template class VectorFacetFacetFE<ET_SEGM>;
  template class VectorFacetFacetFE<ET_TRIG>;
  
  // template class VectorFacetVolumeFE<ET_SEGM>;
  template class VectorFacetVolumeFE<ET_TRIG>;
  template class VectorFacetVolumeFE<ET_QUAD>;
  template class VectorFacetVolumeFE<ET_TET>;
  template class VectorFacetVolumeFE<ET_PRISM>;
  template class VectorFacetVolumeFE<ET_PYRAMID>;
  template class VectorFacetVolumeFE<ET_HEX>;

  static RegisterBilinearFormIntegrator<RobinEdgeIntegrator<3 /* , VectorFacetFacetFiniteElement<2> */ >  > initrvf3 ("robinvectorfacet", 3, 1);
  static RegisterBilinearFormIntegrator<RobinEdgeIntegrator<2 /* , VectorFacetFacetFiniteElement<1> */ >  > initrvf2 ("robinvectorfacet", 2, 1);
  static RegisterLinearFormIntegrator<NeumannEdgeIntegrator<3 /*, VectorFacetFacetFiniteElement<2> */ >  > initnvf3 ("neumannvectorfacet", 3, 1);
  static RegisterLinearFormIntegrator<NeumannEdgeIntegrator<2 /*, VectorFacetFacetFiniteElement<1> */ >  > initnvf2 ("neumannvectorfacet", 2, 1);
}


