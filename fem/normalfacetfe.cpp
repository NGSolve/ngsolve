/*********************************************************************/
/* File:   normalfacetfe.cpp                                         */
/* Author: J. Schoeberl                                              */
/* Date:   2008                                                      */
/*********************************************************************/


#include <fem.hpp>
#include <thdivfe_impl.hpp>
#include "normalfacetfe.hpp"
#include <cassert>

namespace ngfem
{  


  template <int DIM>
  auto GetHDivNormalTIP( const IntegrationPoint & ip);
  
  template<>
  auto GetHDivNormalTIP<1>( const IntegrationPoint & ip)
  {
    TIP<1, AutoDiff<2>> tip(ip.FacetNr(), ip.VB());
    tip.x.Value()= ip(0);
    tip.x.DValue(1) = 1.0;
    tip.x.DValue(0) = 0.0;

    return tip;
  }

  template<>
  auto GetHDivNormalTIP<2>( const IntegrationPoint & ip)
  {
    TIP<2, AutoDiff<3>> tip(ip.FacetNr(), ip.VB());
    tip.x.Value()= ip(0);
    tip.y.Value()= ip(1);

    tip.x.DValue(0) = 0.0;
    tip.x.DValue(1) = 1.0;
    tip.x.DValue(2) = 0.0;

    tip.y.DValue(0) = -1.0;
    tip.y.DValue(1) = 0.0;
    tip.y.DValue(2) = 0.0;

    return tip;
  }

  



  template <ELEMENT_TYPE ET> 
  class NormalFacetVolumeFE_Shape : public NormalFacetVolumeFE<ET>
  {
    using NormalFacetVolumeFE<ET>::order;

    static constexpr int DIM = ngfem::Dim(ET);

  public:
    template<typename Tx, typename TFA>  
      INLINE void T_CalcShape (TIP<DIM,Tx> ip, TFA & shape) const
    { throw Exception ("T_CalcShape not implemented"); }
    
    void CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
    { throw Exception ("dual shape not implemented, H1Ho"); }
  };


  template <> template<typename Tx, typename TFA>  
  void NormalFacetVolumeFE_Shape<ET_TRIG>::
  T_CalcShape (TIP<DIM,Tx> ip, TFA & shape) const
  {
    if (ip.vb != BND) throw Exception("normal-facet element evaluated not at BND");
    Tx lami[3] = { ip.x, ip.y, 1-ip.x-ip.y } ;

    for (int i = 0; i < 3; i++)
      {
        if (ip.facetnr != i)
          {
            for (int j : Range(first_facet_dof[i], first_facet_dof[i+1]))
              shape[j] = Du(Tx(0.0));
            continue;
          }
          
        int first = first_facet_dof[i];
        int p = facet_order[i][0];
        
        INT<2> e = GetVertexOrientedEdge(i);
        Tx xi = lami[e[0]] - lami[e[1]];
        
        LegendrePolynomial (p, xi, 
                            SBLambda([&](int nr, auto val)
                                     {
                                       shape[first+nr] = uDv (val, xi);
                                     }));
      }
  }

  template <> template<typename Tx, typename TFA>  
  void NormalFacetVolumeFE_Shape<ET_QUAD>::
  T_CalcShape (TIP<DIM,Tx> ip, TFA & shape) const
  {
    if (ip.vb != BND) throw Exception("normal-facet element evaluated not at BND");
    
    Tx x = ip.x, y = ip.y;
    //Tx lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  

    for (int i = 0; i < 4; i++)
      {
        if (ip.facetnr != i)
          {
            for (int j : Range(first_facet_dof[i], first_facet_dof[i+1]))
              shape[j] = Du(Tx(0.0));
            continue;
          }
          
        int first = first_facet_dof[i];
        int p = facet_order[i][0];
        
        INT<2> e = GetVertexOrientedEdge(i);
        Tx xi = sigma[e[0]]-sigma[e[1]];
        
        LegendrePolynomial (p, xi, 
                            SBLambda([&](int nr, auto val)
                                     {
                                       shape[first+nr] = uDv (val, xi);
                                     }));
      }
  }

  
  template <> template<typename Tx, typename TFA>  
  void NormalFacetVolumeFE_Shape<ET_TET>::
  T_CalcShape (TIP<DIM,Tx> ip, TFA & shape) const
  {
    if (ip.vb != BND) throw Exception("normal-facet element evaluated not at BND");
    Tx lam[4] = { ip.x, ip.y, ip.z, 1-ip.x-ip.y-ip.z } ;

    for (int i = 0; i < 4; i++)
      {
        if (ip.facetnr != i)
          {
            for (int j : Range(first_facet_dof[i], first_facet_dof[i+1]))
              shape[j] = Du_Cross_Dv (Tx(0.0), Tx(0.0));
            continue;
          }
          
        int first = first_facet_dof[ip.facetnr];
        int p = facet_order[ip.facetnr][0];

        INT<4> f = GetVertexOrientedFace (ip.facetnr);
        auto xi = lam[f[0]]-lam[f[2]];
        auto eta = lam[f[1]]-lam[f[2]];
        DubinerBasis::Eval (p, lam[f[0]], lam[f[1]],
                            SBLambda([&](int nr, auto val)
                                     {
                                       shape[first+nr] = wDu_Cross_Dv (xi, eta, val);
                                     }));
      }
  }

    template <> template<typename Tx, typename TFA>  
  void NormalFacetVolumeFE_Shape<ET_HEX>::
  T_CalcShape (TIP<DIM,Tx> ip, TFA & shape) const
  {
    if (ip.vb != BND) throw Exception("normal-facet element evaluated not at BND");
    
    Tx x = ip.x, y = ip.y, z = ip.z;

    Tx lami[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
                  (1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
    Tx sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
                   (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z};

    const FACE * faces = ElementTopology::GetFaces (ET_HEX);
    for (int i = 0; i < 6; i++)
      {
        if (ip.facetnr != i)
          {
            for (int j : Range(first_facet_dof[i], first_facet_dof[i+1]))
              shape[j] = Du_Cross_Dv (Tx(0.0), Tx(0.0));
            continue;
          }
          
        int first = first_facet_dof[i];
        int p = facet_order[i][0];
        

	Tx lam_f(0);
	for (int j = 0; j < 4; j++)
	  lam_f += lami[faces[i][j]];


        INT<4> f = GetFaceSort (i, vnums);	  
        Tx xi  = sigma[f[0]]-sigma[f[1]];
        Tx eta = sigma[f[0]]-sigma[f[3]];

        shape[first] = wDu_Cross_Dv(eta, xi, -0.25*lam_f);


        ArrayMem<Tx, 20> L_xi(order+2),L_eta(order+2);

        IntegratedLegendrePolynomial::Eval(p+1,xi,L_xi);
        IntegratedLegendrePolynomial::Eval(p+1,eta,L_eta);
        
        /*for(int j = 0; j <= p; j++)
          for(int k = 0; k <= p; k++)
          shape[first + ii++] = L_xi[j]*L_eta[k]*lam_f;

        for (int k = 0; k < p; k++)
          for (int l = 0; l < p; l++, ii++)
            shape[first + ii] = curl_uDvw_minus_Duvw(L_xi[k+2],L_eta[l+2],-lam_f); //divfree
        */

          // could be simpler, but leads to bug in XCode 11.6 (and certainly other):
        for (size_t l = 0; l < p; l++)
          for (size_t k = 0; k < p; k++)
            shape[first + k*p+l+1] = curl_uDvw_minus_Duvw(L_xi[k+2],L_eta[l+2],-lam_f); //divfree
        size_t ii = p*p+1;
        
        for (int k = 0; k < p; k++)
          shape[first + ii++] = Du_Cross_Dv(L_xi[k+2]*lam_f,eta); //divfree

        for (int k = 0; k < p; k++)
          shape[first + ii++] = Du_Cross_Dv(L_eta[k+2]*lam_f,xi); //divfree
      }

  }
    
    



  
  template <ELEMENT_TYPE ET> template<typename Tx, typename TFA>  
  void NormalFacetVolumeFE<ET>::
  T_CalcShape (Tx hx[DIM], int fnr, TFA & shape) const
  {
    throw ExceptionNOSIMD("NormalFacetVolume::T_CalcShape missing element "+ToString(ET));
  }

  template <ELEMENT_TYPE ET> template<typename MIP, typename TFA>  
  void NormalFacetVolumeFE<ET>::
  CalcDualShape2 (MIP mip, int fnr, TFA & shape) const
  {
    throw Exception("NormalFacetVolume::CalcDualShape2 missing element "+ToString(ET));
  }

  template <> template<typename Tx, typename TFA>  
  void NormalFacetVolumeFE<ET_TRIG>::
  T_CalcShape (Tx hx[DIM], int fanr, TFA & shape) const
  {
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

  template <> template<typename MIP, typename TFA>  
  void NormalFacetVolumeFE<ET_TRIG>::
  CalcDualShape2 (MIP mip, int fanr, TFA & shape) const
  {
    auto & ip = mip.IP();
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
    T x = ip(0), y = ip(1);
    T lami[3] = { x, y, 1-x-y };
    Vec<2,T> pnts[3] = { { 1, 0 }, { 0, 1 } , { 0, 0 } };

    int first = first_facet_dof[fanr];
    int p = facet_order[fanr][0];

    if (ip.VB() == BND)
      {
        INT<2> e = GetVertexOrientedEdge(fanr);
        T xi = lami[e[0]] - lami[e[1]];
        Vec<2,T> tauref = pnts[e[0]] - pnts[e[1]];
        Vec<2,T> nvref = Vec<2,T>(tauref[1],-tauref[0]);
        Vec<2,T> nv = L2Norm(tauref) / L2Norm(mip.GetJacobian()*tauref) * Cof(mip.GetJacobian()) * nvref;
        
        LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              shape[first+nr] = val * nv;
                            }));
      }
  }


  template<> template<typename Tx, typename TFA>  
  void NormalFacetVolumeFE<ET_TET> ::
  T_CalcShape (Tx hx[DIM], int fanr, TFA & shape ) const
  {
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

  
#ifdef OLD
  template <ELEMENT_TYPE ET>
  void NormalFacetVolumeFE<ET>::
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
  void NormalFacetVolumeFE<ET>::
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
  void NormalFacetVolumeFE<ET>::
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
#endif

  
  /* **************************** Facet Segm ********************************* */

  template <ELEMENT_TYPE ET> template<typename Tx, typename TFA>  
  void NormalFacetFacetFE<ET>::
  T_CalcShape (TIP<DIM,Tx> tip, TFA & shape) const
  {
    throw ExceptionNOSIMD("NormalFacetFacet::T_CalcShape missing"+ToString(ET));
  }


  
  template<ELEMENT_TYPE ET>
  void NormalFacetFacetFE<ET>::CalcShape(const IntegrationPoint & ip,
                                         FlatVector<> shape) const
  {
    auto tip = GetHDivNormalTIP<DIM>(ip);
    T_CalcShape (tip,
                 SBLambda([shape] (size_t i, auto val)
                          {
                            auto vshape = HDiv2ShapeNew(val);
                            shape(i) = vshape(DIM);
                          }));
  }
  
  template<> template <typename Tx, typename TFA>
  void NormalFacetFacetFE<ET_SEGM>::T_CalcShape(TIP<DIM,Tx> tip,
                                                TFA & shape) const
  {

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
  void NormalFacetFacetFE<ET_SEGM>::ComputeNDof()
  {
    order = order_inner[0];
    ndof = order+1;
  }

  /* **************************** Facet Trig ********************************* */

  template<> template <typename Tx, typename TFA>
  void NormalFacetFacetFE<ET_TRIG>::T_CalcShape(TIP<DIM,Tx> tip, 
                                                TFA &  shape) const
  {
    Tx lam[3] = { tip.x, tip.y, 1-tip.x-tip.y } ;

  
    INT<4> f = GetVertexOrientedFace (0);
    auto xi = lam[f[0]]-lam[f[2]];
    auto eta = lam[f[1]]-lam[f[2]];
    DubinerBasis::Eval (order_inner[0], lam[f[0]], lam[f[1]],
                        SBLambda([&](int nr, auto val)
                                 {
                                   shape[nr] = wDu_Cross_Dv (xi, eta, val);
                                 }));
    
    /*auto x = tip.x;
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
                                  }));*/
  }

  template<>
  void NormalFacetFacetFE<ET_TRIG>::ComputeNDof()
  {
    order = order_inner[0];
    int p = order_inner[0];
    ndof = (p+1)*(p+2)/2;
  }


  /* **************************** Facet Quad ********************************* */

  template<> template <typename Tx, typename TFA>
  void NormalFacetFacetFE<ET_QUAD>::T_CalcShape(TIP<DIM,Tx> tip, 
                                                TFA &  shape) const
  {
    // Tx lam[4] = {(1-tip.x)*(1-tip.y),tip.x*(1-tip.y),tip.x*tip.y,(1-tip.x)*tip.y};  
    Tx sigma[4] = {(1-tip.x)+(1-tip.y),tip.x+(1-tip.y),tip.x+tip.y,(1-tip.x)+tip.y};


    Tx lam_f(0);
    int p = order;
    int ii = 0;
    
    INT<4> f = GetFaceSort (0, vnums);	  
    Tx xi  = sigma[f[0]]-sigma[f[1]];
    Tx eta = sigma[f[0]]-sigma[f[3]];
    lam_f = Cross(xi,eta);
    lam_f.Value() = 1;
    shape[ii++] = wDu_Cross_Dv(eta, xi, -0.25*lam_f);


    ArrayMem<Tx, 20> L_xi(p+2),L_eta(p+2);
    
    IntegratedLegendrePolynomial::Eval(p+1,xi,L_xi);
    IntegratedLegendrePolynomial::Eval(p+1,eta,L_eta);
    
    /*for(int j = 0; j <= p; j++)
      for(int k = 0; k <= p; k++)
      shape[ii++] = L_xi[j]*L_eta[k]*lam_f;*/

    for (int k = 0; k < p; k++)
      for (int l = 0; l < p; l++, ii++)
        shape[ii] = curl_uDvw_minus_Duvw(L_xi[k+2],L_eta[l+2],-lam_f); //divfree
    
    for (int k = 0; k < p; k++)
      shape[ii++] = Du_Cross_Dv(L_xi[k+2]*lam_f,eta); //divfree
    
    for (int k = 0; k < p; k++)
      shape[ii++] = Du_Cross_Dv(L_eta[k+2]*lam_f,xi); //divfree

  }

  template<>
  void NormalFacetFacetFE<ET_QUAD>::ComputeNDof()
  {
    order = max2( order_inner[0], order_inner[1] );
    ndof = (order_inner[0]+1) * (order_inner[1]+1);
  }


  /* **************************** Volume Trig ********************************* */

  template<>
  void NormalFacetVolumeFE<ET_TRIG> ::
  CalcShape ( const IntegrationPoint & ip, int fanr, SliceMatrix<> shape ) const
  {
    throw Exception ("normal facet not readsy 9896");        
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
  void NormalFacetVolumeFE<ET_TRIG> ::
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
  void NormalFacetVolumeFE<ET_TRIG> :: ComputeNDof()
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
  void NormalFacetVolumeFE<ET_QUAD> ::
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
  void NormalFacetVolumeFE<ET_QUAD> ::
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
  void NormalFacetVolumeFE<ET_QUAD> :: ComputeNDof()
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
  void NormalFacetVolumeFE<ET_TET> ::
  CalcShape ( const IntegrationPoint & ip, int fanr, SliceMatrix<> shape ) const
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
  void NormalFacetVolumeFE<ET_TET> ::
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
  void NormalFacetVolumeFE<ET_TET> :: ComputeNDof()
  {
    ndof = 0;
    for (int i = 0; i < 4; i++)
      {
	first_facet_dof[i] = ndof;
	ndof += (facet_order[i][0]+1) * (facet_order[i][0]+2) / 2;
      }
    first_facet_dof[4] = ndof;
  }

  /* **************************** Volume Prism ********************************* */

  template<>
  void NormalFacetVolumeFE<ET_PRISM> ::
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
  int NormalFacetVolumeFE<ET_PRISM> ::
  GetNExtraShapes ( int fanr) const
  {
    if (fanr < 2) //trig shape
      return 2*(facet_order[fanr][0]+2);
    else //quad shape
      return 2*(2*facet_order[fanr][0]+3);
  }

  template<>
  void NormalFacetVolumeFE<ET_PRISM> ::
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
  void NormalFacetVolumeFE<ET_PRISM>::ComputeNDof()
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
  void NormalFacetVolumeFE<ET_HEX> ::
  CalcShape ( const IntegrationPoint & ip, int fanr, SliceMatrix<> shape ) const
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
    
    shape = 0.0;
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
  void NormalFacetVolumeFE<ET_HEX>::ComputeNDof()
  {
    ndof = 0;
    for (int i=0; i<6; i++)
      {
	first_facet_dof[i] = ndof;
	ndof += (facet_order[i][0]+1) * (facet_order[i][0]+1);
      }
    first_facet_dof[6] = ndof;
  }



  /* **************************** Volume Pyramid ********************************* */

  template<>
  void NormalFacetVolumeFE<ET_PYRAMID> ::
  CalcShape ( const IntegrationPoint & ip, int facet, SliceMatrix<> shape ) const
  {
    throw Exception("NormalFacetVolumePyramid::CalcShape: not implemented!");
  }

//   void NormalFacetVolumePyramid::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatMatrix<> shape) const
//   {
//     ;
//   }
   
//   void NormalFacetVolumePyramid::SetFacet(int afnr) const
//   {
//     if (qnr == afnr || tnr == afnr) return;
  
//     NormalFacetVolumePyramid * pyramid=const_cast<NormalFacetVolumePyramid*>(this);
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
  void NormalFacetVolumeFE<ET_PYRAMID>::ComputeNDof()
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

  template<ELEMENT_TYPE ET>
  void NormalFacetVolumeFE<ET> :: CalcDualShape (const BaseMappedIntegrationPoint & bmip, SliceMatrix<> shape) const
    {
      shape = 0.0;
      Switch<4-DIM>
      (bmip.DimSpace()-DIM,[this,&bmip,shape](auto CODIM)
       {
         auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);
         this -> CalcDualShape2 (mip, mip.IP().FacetNr(), SBLambda([shape] (size_t i, auto val) { shape.Row(i) = val; }));
       });
    }

  
  template class NormalFacetFacetFE<ET_SEGM>;
  template class NormalFacetFacetFE<ET_TRIG>;
  template class NormalFacetFacetFE<ET_QUAD>;
  
  // template class NormalFacetVolumeFE<ET_SEGM>;
  template class NormalFacetVolumeFE<ET_TRIG>;
  template class NormalFacetVolumeFE<ET_QUAD>;
  template class NormalFacetVolumeFE<ET_TET>;
  template class NormalFacetVolumeFE<ET_PRISM>;
  template class NormalFacetVolumeFE<ET_PYRAMID>;
  template class NormalFacetVolumeFE<ET_HEX>;

  template class T_HDivFiniteElement<NormalFacetVolumeFE_Shape<ET_TRIG>, ET_TRIG>;
  template class T_HDivFiniteElement<NormalFacetVolumeFE_Shape<ET_QUAD>, ET_QUAD>;
  template class T_HDivFiniteElement<NormalFacetVolumeFE_Shape<ET_TET>, ET_TET>;
  template class T_HDivFiniteElement<NormalFacetVolumeFE_Shape<ET_PRISM>, ET_PRISM>;
  template class T_HDivFiniteElement<NormalFacetVolumeFE_Shape<ET_PYRAMID>, ET_PYRAMID>;
  template class T_HDivFiniteElement<NormalFacetVolumeFE_Shape<ET_HEX>, ET_HEX>;
}


