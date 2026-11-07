#ifndef FILE_HCURLHOFE_IMPL
#define FILE_HCURLHOFE_IMPL

/*********************************************************************/
/* File:   hcurlhofe.hpp                                             */
/* Author: Sabine Zaglmayr, Joachim Schoeber                         */
/* Date:   20. Maerz 2003                                            */
/*                                                                   */
/* AutoCurl - revision: J. Schoeberl, March 2009                     */
/*********************************************************************/

#include "recursive_pol.hpp"
// #include "thdivfe.hpp"
#include "hcurlhofe.hpp"
#include "hcurlfe_utils.hpp"

namespace ngfem
{
  // 1/k, sqrt(k), 1/sqrt(k) for biorth trig coefficients
  struct BiorthTabs
  {
    std::array<double,1024> inv{}, sq{}, isq{};
    BiorthTabs ()
    {
      for (size_t k = 1; k < inv.size(); k++)
        { inv[k] = 1.0/k; sq[k] = sqrt(double(k)); isq[k] = 1.0/sq[k]; }
    }
  };
  inline const BiorthTabs biorth_tabs;

  // biorth trig w~ coefficients U = c0 g - c1 gm1, V = c2 g - c3 gm1, index [i][n], p < 64
  struct BiorthCoefs {
    static constexpr int N = 64;
    std::array<std::array<double,4>,N*N> c{};
    BiorthCoefs () {
      for (int i = 2; i < N; i++)
        for (int n = 0; i+n < N; n++) {
          double t = sqrt(2.*i-1)/sqrt(2.*n+2*i)/(i-1);
          double al = t*(n+2*i)/(n+i+1), be = t*n/(n+i);
          c[i*N+n] = { al*(n+1), be*n, al*i, be*i };
        }
    }
  };
  inline const BiorthCoefs biorth_coefs;

  // dual shape helpers

  // (a,b) with X - P0 = a (P1-P0) + b (P2-P0), X in the plane of the reference face
  template <typename T>
  INLINE void HCurlFaceCoords (const Vec<3,T> & X, Vec<3> P0, Vec<3> P1, Vec<3> P2, T & a, T & b)
  {
    Vec<3> d1 = P1-P0, d2 = P2-P0;
    double g11 = InnerProduct(d1,d1), g12 = InnerProduct(d1,d2), g22 = InnerProduct(d2,d2);
    double idet = 1/(g11*g22-g12*g12);
    T r1(0.0), r2(0.0);
    for (int k = 0; k < 3; k++)
      {
        T rk = X(k)-P0(k);
        r1 += d1(k)*rk;
        r2 += d2(k)*rk;
      }
    a = idet * (g22*r1 - g12*r2);
    b = idet * (g11*r2 - g12*r1);
  }

  // face map 1/(det*measure) J F, F = (t0,t1) reference face tangents
  template <typename T, typename MIP>
  INLINE Mat<3,2,T> HCurlFaceTrafo (const MIP & mip, Vec<3> t0, Vec<3> t1)
  {
    double g11 = InnerProduct(t0,t0), g12 = InnerProduct(t0,t1), g22 = InnerProduct(t1,t1);
    double det = sqrt(g11*g22-g12*g12);
    Mat<3,3,T> jac = mip.GetJacobian();
    T scal = 1/(det*mip.GetMeasure());
    Mat<3,2,T> ftrafo;
    for (int i = 0; i < 3; i++)
      {
        ftrafo(i,0) = scal * (jac(i,0)*t0(0) + jac(i,1)*t0(1) + jac(i,2)*t0(2));
        ftrafo(i,1) = scal * (jac(i,0)*t1(0) + jac(i,1)*t1(1) + jac(i,2)*t1(2));
      }
    return ftrafo;
  }

  // trig face duals, as for the tet; xi, eta = barycentrics of the sorted vertices fav[0], fav[1]
  template <typename T, typename TFA>
  INLINE void HCurlTrigFaceDual (const Mat<3,2,T> & ftrafo, T xi, T eta, int p,
                                 bool usegrad, bool type1, int & ii, TFA & shape)
  {
    if (p < 2) return;
    if (!usegrad)
      {
        // surface rot of P_{p-1}, annihilates the face gradients
        DubinerBasis::Eval(p-1, AutoDiff<2,T>(xi, 0), AutoDiff<2,T>(eta, 1),
                           SBLambda([&] (size_t nr, auto val)
                                    {
                                      if (nr == 0) return;
                                      shape[ii++] = ftrafo*Vec<2,T> (val.DValue(1), -val.DValue(0));
                                    }));
        return;
      }
    DubinerBasis::Eval(p-2, xi, eta,
                       SBLambda([&] (size_t nr, auto val)
                                {
                                  shape[ii++] = ftrafo*Vec<2,T> (val, 0);
                                  if (type1)
                                    shape[ii++] = ftrafo*Vec<2,T>(0, val);
                                  else
                                    shape[ii++] = ftrafo*Vec<2,T> (val*xi, val*eta);
                                }));
    if (!type1)   // P^(1,0)(2xi-1), orthogonal on the face
      JacobiPolynomialFix<1,0>::Eval(p-2, 2*xi-1,
                                     SBLambda([&] (size_t nr, auto val)
                                              {
                                                shape[ii++] = ftrafo*Vec<2,T>(0, val);
                                              }));
  }

  // Legendre L_i(x) with dL_i/dx, x = 2s-1
  template <typename T>
  INLINE void HCurlLegendreAD (int n, T s, FlatArray<AutoDiff<1,T>> l)
  {
    LegendrePolynomial::Eval (n, AutoDiff<1,T>(2*s-1, 0), l);
  }
  // F_i = int_{-1}^x L_i: F_0 = x+1, F_i = (L_{i+1}-L_{i-1})/(2i+1), F_i' = L_i
  template <typename T>
  INLINE T HCurlIntLeg (FlatArray<AutoDiff<1,T>> l, int i)
  {
    if (i == 0) return l[1].Value()+1;
    return (1.0/(2*i+1)) * (l[i+1].Value()-l[i-1].Value());
  }

  /*
    quad duals in coordinates (s,t) in [0,1]^2, x = 2s-1, y = 2t-1
    gradients of H1(p+1-type1): (F_i(x) L_j(y), 0), div onto Q_{p-1-type1}
    curl part of N_[p-1]: rot (L_i(x) L_j(y)), (i,j) != (0,0), div-free
   */
  template <typename T, typename TRAFO, typename TFA>
  INLINE void HCurlQuadDual (const TRAFO & trafo, T s, T t, IVec<2> p,
                             bool usegrad, bool type1, int & ii, TFA & shape)
  {
    int n = max(p[0], p[1]) + 1;
    ArrayMem<AutoDiff<1,T>,20> ls(n+1), lt(n+1);
    HCurlLegendreAD (n, s, ls);
    HCurlLegendreAD (n, t, lt);

    int pg0 = p[0]-type1, pg1 = p[1]-type1;
    if (usegrad && pg0 >= 1 && pg1 >= 1 && !type1)
      // grad of Q_{p+1} bubbles F_{i+1}(x) F_{j+1}(y): blind to div-free fields
      for (int i = 0; i < pg0; i++)
        {
          T bi = HCurlIntLeg<T> (ls, i+1);
          for (int j = 0; j < pg1; j++)
            shape[ii++] = trafo * Vec<2,T>(ls[i+1].Value()*HCurlIntLeg<T> (lt, j+1),
                                           bi*lt[j+1].Value());
        }
    else if (usegrad && pg0 >= 1 && pg1 >= 1)
      for (int i = 0; i < pg0; i++)
        {
          T fi = HCurlIntLeg<T> (ls, i);
          for (int j = 0; j < pg1; j++)
            shape[ii++] = trafo * Vec<2,T>(fi*lt[j].Value(), 0);
        }

    if (p[0] < 1 || p[1] < 1) return;
    for (int i = 0; i < p[0]; i++)
      for (int j = 0; j < p[1]; j++)
        if (i || j)
          shape[ii++] = trafo * Vec<2,T>(ls[i].Value()*lt[j].DValue(0),
                                         -ls[i].DValue(0)*lt[j].Value());
  }



  // declaration of the shapes ...

  

  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET2> class TSHAPES, typename BASE>
  void HCurlHighOrderFE<ET,TSHAPES,BASE> :: ComputeNDof()
  {
    ndof = N_EDGE;

    for (int i = 0; i < N_EDGE; i++)
      if(order_edge[i] > 0)
        ndof += usegrad_edge[i]*order_edge[i];
    
    for (int i = 0; i < N_FACE; i++)
      ndof += hcurl_ndof::Face (FaceType(i), order_face[i], usegrad_face[i], type1);

    if (DIM == 3)
      ndof += hcurl_ndof::Cell (ET, order_cell, usegrad_cell, type1);
    
    TORDER horder = 0; 
    for (int i = 0; i < N_EDGE; i++)
      horder = max2 (horder, order_edge[i]);

    for(int i=0; i < N_FACE; i++) 
      if (ET_trait<ET>::FaceType(i) == ET_TRIG)
        horder = max2 (horder, order_face[i][0]);
      else
        horder = max2 (horder, Max (order_face[i]));

    if (DIM == 3)
      horder = max2 (horder, Max(order_cell));

    // for integration order .. 
    if (ET == ET_PRISM || ET == ET_HEX || ET == ET_PYRAMID || ET == ET_QUAD)
      horder++;
    else
      if (horder==0) horder++;
    order = horder;
  }



  
  template <ELEMENT_TYPE ET> 
  class HCurlHighOrderFE_Shape :  public HCurlHighOrderFE<ET>
  {
    using ET_trait<ET>::DIM;    
    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::N_EDGE;
    using ET_trait<ET>::N_FACE;
    using ET_trait<ET>::N_CELL;
    using ET_trait<ET>::FaceType;
    using Base = HCurlHighOrderFE<ET>;
    using Base::type1;
    using Base::biorth;
    using Base::order;
    using Base::order_edge;
    using Base::order_face;
    using Base::order_cell;
    using Base::vnums;
    using Base::GetEdgeSort;
    using Base::GetFaceSort;
    using Base::usegrad_cell;
    using Base::usegrad_edge;
    using Base::usegrad_face;
    using typename Base::EdgeOrthoPol;
    using typename Base::QuadOrthoPol;
    using typename Base::TrigOrthoPolGrad;
    using typename Base::T_ORTHOPOL;
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (TIP<DIM,Tx> ip, TFA & shape) const;

    template <typename MIP, typename TFA>
    inline void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception(string("CalcDualShape missing for HighOrderHCurl element ")+ElementTopology::GetElementName(ET));
    }
  };



  
  //------------------------------------------------------------------------
  // HCurlHighOrderSegm
  //------------------------------------------------------------------------

  
  template<> template<typename Tx, typename TFA>  
  void HCurlHighOrderFE_Shape<ET_SEGM> :: T_CalcShape (TIP<1,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x;
    Tx lam[2] = { x, 1-x };

    IVec<2> e = GetEdgeSort (0, vnums);	  
    
    //Nedelec low order edge shape function 
    shape[0] = uDv_minus_vDu (lam[e[0]], lam[e[1]]);

    int p = order_edge[0]; //order_cell[0]; 
    //HO-Edge shapes (Gradient Fields)   
    if(p > 0 && usegrad_cell)
      { 
        // LegendrePolynomial::
        size_t ii = 1;
        EdgeOrthoPol::
          EvalScaledMult (p-1, 
                          lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
                          lam[e[0]]*lam[e[1]], 
                          SBLambda ([&](int nr, Tx val)
                                    {
                                      shape[ii++] = Du (val);
                                    }));
      }
  }
  


  
  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_SEGM> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    // tangential moments along the segment, as the edge duals of the trig/quad elements
    auto & ip = mip.IP();
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;
    if (ip.VB() != VOL) return;
    T x = ip(0);
    T lam[2] = { x, 1-x };
    double pnts[2] = { 1, 0 };
    IVec<2> e = GetEdgeSort (0, vnums);
    T xi = lam[e[1]]-lam[e[0]];
    auto tau = mip.GetJacobian() * Vec<1,T>(pnts[e[1]]-pnts[e[0]]);
    tau /= mip.GetMeasure();
    int p = usegrad_cell ? order_edge[0] : 0;
    LegendrePolynomial::Eval
      (p, xi,
       SBLambda([&] (size_t nr, T val)
                {
                  auto vshape = val * tau;
                  shape[nr] = vshape;
                }));
  }

  //------------------------------------------------------------------------
  // HCurlHighOrderTrig
  //------------------------------------------------------------------------
  
  template<> template<typename Tx, typename TFA>  
  INLINE void HCurlHighOrderFE_Shape<ET_TRIG> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x, y = ip.y;
    Tx lam[3] = { x, y, 1-x-y };

	
    int ii = 3; 
    for (int i = 0; i < 3; i++)
      {
        IVec<2> e = GetEdgeSort (i, vnums);	  

	//Nedelec low order edge shape function 
        shape[i] = uDv_minus_vDu (lam[e[0]], lam[e[1]]);

	int p = order_edge[i]; 
	//HO-Edge shapes (Gradient Fields)   
	if(p > 0 && usegrad_edge[i]) 
	  { 
	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalScaledMult (order_edge[i]-1, 
			      lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
			      lam[e[0]]*lam[e[1]], 
                              // adpol1
                              SBLambda ([&](int nr, Tx val)
                                        {
                                          shape[ii++] = Du (val);
                                        }));
	  }
      }   
     
    //Inner shapes (Face) 
    int p = order_face[0][0];      
    if(p > 1) 
      {
	IVec<4> fav = GetFaceSort (0, vnums);

	// Tx xi  = lam[fav[2]]-lam[fav[1]];
	// Tx eta = lam[fav[0]]; 
	int pg = p - 2 - (type1 ? 1 : 0);

	// gradients:	
	if (usegrad_face[0] && pg >=0)
            
	  // DubinerBasis
          TrigOrthoPolGrad::EvalMult
	    (pg, lam[fav[0]], lam[fav[1]], 
	     lam[fav[0]]*lam[fav[1]]*lam[fav[2]], 
	     SBLambda
	     ([&](int nr, Tx val)
	      {
		shape[ii++] = Du (val);
	      }));

        if (biorth)
          {
            // C = fav[0] (collapse), B = fav[1], A = fav[2]
            // f_i ~ lamA lamB P^(1,1)_{i-2}, g_ij ~ lamC IntJac^(2i)_{j-1}
            // w_ij = j g grad f - i f grad g,  i+j <= p+1
            // combine w_n, w_{n-1}: curl ~ Dubiner L_{i-1}(xi) P^(2i-1,0)_n(eta)
            Tx lc = lam[fav[0]], lb = lam[fav[1]], la = lam[fav[2]];
            auto & tb = biorth_tabs;
            JacobiPolynomialFix<1,1>::EvalScaledMult
              (p-2, lb-la, la+lb, la*lb, SBLambda([&](int m, Tx f) LAMBDA_INLINE
              {
                int i = m+2;
                Tx gm1(0.0);
                IntegratedJacobiPolynomialAlpha(2*i).EvalMult
                  (p-2-m, 2*lc-1, lc, SBLambda([&](int n, Tx g) LAMBDA_INLINE
                  {
                    // |curl w~|_ref = 1
                    auto & cf = biorth_coefs.c[i*BiorthCoefs::N+n];
                    shape[ii++] = uDw_minus_wDv (cf[0]*g - cf[1]*gm1, f, cf[2]*g - cf[3]*gm1);
                    gm1 = g;
                  }));
              }));
            // Whitney(A,B) * lamC P^(2,1)_{j-1}, |curl|_ref = 1
            IntegratedJacobiPolynomialAlpha(3).EvalMult
              (p-2, 2*lc-1, lc, SBLambda([&](int j, Tx g) LAMBDA_INLINE
              {
                shape[ii++] = wuDv_minus_wvDu (la, lb, tb.sq[2*j+4]*g);
              }));
          }
	// now ready: compatibility with tets/prisms ! non-gradient inner shapes, same for type 1 and 2
        else if (true) // if (type1)
          {
            DubinerBasis::EvalMult
              (p-2, lam[fav[0]], lam[fav[1]], 
               lam[fav[0]], 
               SBLambda
               ([&](int nr, Tx val)
                {
                  shape[ii++] = wuDv_minus_wvDu (lam[fav[1]], lam[fav[2]], val);
                }));
            
            LegendrePolynomial::EvalMult 
              (p-2, lam[fav[2]]-lam[fav[1]], lam[fav[2]], 
               SBLambda
               ([&] (int j, Tx val)
                {
                  shape[ii++] = wuDv_minus_wvDu (lam[fav[1]], lam[fav[0]], val);
                }));
          }
        else
          {
            Tx xi  = lam[fav[2]]-lam[fav[1]];
            Tx eta = lam[fav[0]];
            
            ArrayMem<Tx,20> adpol1(order),adpol2(order);	
            TrigShapesInnerLegendre::CalcSplitted(p+1, xi, eta, adpol1,adpol2);
	
            // other combination
            for (int j = 0; j < p-1; j++)
              for (int k = 0; k < p-1-j; k++, ii++)
                shape[ii] = uDv_minus_vDu (adpol2[k], adpol1[j]);     
            
            // rec_pol * Nedelec0 
            for (int j = 0; j < p-1; j++, ii++)
              shape[ii] = wuDv_minus_wvDu (lam[fav[1]], lam[fav[2]], adpol2[j]);
          }
      }
  }


  //------------------------------------------------------------------------
  // HCurlHighOrderQuad
  //------------------------------------------------------------------------
  

  template<> template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE_Shape<ET_QUAD> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x, y = ip.y;
    Tx hx[2] = { x, y };
    Tx lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
    Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  

    int ii = 4;
    ArrayMem<Tx, 10> pol_xi(order+2), pol_eta(order+2);

    for (int i = 0; i < 4; i++)
      {
	// int p = order_edge[i]; 
        IVec<2> e = GetEdgeSort (i, vnums);	  
	
	Tx xi  = sigma[e[1]]-sigma[e[0]];
	Tx lam_e = lami[e[0]]+lami[e[1]];  
        Tx bub = 0.25 * lam_e * (1 - xi*xi);

	// Nedelec0-shapes
	shape[i] = uDv (0.5 * lam_e, xi); 

	// High Order edges ... Gradient fields 
	if(usegrad_edge[i])
	  {
            /*
	    LegendrePolynomial::
	      EvalMult (order_edge[i]-1, 
			xi, bub, pol_xi);

	    for (int j = 0; j < p; j++)
              shape[ii++] = Du<2> (pol_xi[j]);
            */
	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalMult (order_edge[i]-1, 
			xi, bub, SBLambda ([&](int i, Tx val)
                                           {
                                             shape[ii++] = Du (val);
                                           }));

	  }
      }



     
    IVec<2> p = order_face[0]; // (order_cell[0],order_cell[1]);
    IVec<2> pg (p[0]-type1, p[1]-type1);

    if (usegrad_face[0] && pg[0] >= 1 && pg[1] >= 1)
      {
        Vec<2,Tx> xi = ET_trait<ET_QUAD>::XiFace(0, hx, vnums);
	Tx bub = 1.0/16 * (1-xi(0)*xi(0))*(1-xi(1)*xi(1));

        QuadOrthoPol::EvalMult1Assign(pg[0]-1, xi(0), bub,
              SBLambda ([&](int i, Tx val) LAMBDA_INLINE
                    {
                      QuadOrthoPol::EvalMult (pg[1]-1, xi(1), val,
                                                    SBLambda([&](int i2, Tx v2)
                                                             {
                                                               shape[ii++] = Du (v2);
                                                             }));
                    }));
      }





    // curl part: the one of N_[p-1], curl onto Q_{p-1}
    if (p[0] < 1 || p[1] < 1) return;
    IVec<2> q (p[0]-1, p[1]-1);

    int fmax = 0;
    for (int j = 1; j < 4; j++)
      if (vnums[j] > vnums[fmax])
	fmax = j;

    int f1 = (fmax+3)%4;
    int f2 = (fmax+1)%4;
    if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f2 > f1;

    Tx xi = sigma[fmax]-sigma[f1];  // in [-1,1]
    Tx eta = sigma[fmax]-sigma[f2]; // in [-1,1]

    T_ORTHOPOL::Calc(q[0]+1, xi,pol_xi);
    T_ORTHOPOL::Calc(q[1]+1,eta,pol_eta);

    //Rotation of Gradient fields
    for (int k = 0; k < q[0]; k++)
      for (int j= 0; j < q[1]; j++)
        shape[ii++] = uDv_minus_vDu (pol_eta[j], pol_xi[k]);

    //Missing ones
    for(int j = 0; j< q[0]; j++)
      shape[ii++] = uDv (0.5*pol_xi[j], eta);

    for(int j = 0; j < q[1]; j++)
      shape[ii++] = uDv (0.5*pol_eta[j], xi);
  }


  
  //------------------------------------------------------------------------
  //        Tetrahedron
  //------------------------------------------------------------------------
 

  template<> template<typename Tx, typename TFA>  
  void HCurlHighOrderFE_Shape<ET_TET> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx x = ip.x, y = ip.y, z = ip.z;    
    Tx lam[4] = { x, y, z, 1-x-y-z };

    ArrayMem<Tx,20> adpol1(order+2),adpol2(order+2),adpol3(order+2); 
    int ii = 6; 

    for (int i = 0; i < N_EDGE; i++)
      { 
	int p = order_edge[i]; 
        IVec<2> e = GetEdgeSort (i, vnums);	  
	
	//Nedelec low order edge shape function 
        shape[i] = uDv_minus_vDu (lam[e[0]], lam[e[1]]);

	//HO-Edge shape functions (Gradient Fields) 	
	if (p > 0 && usegrad_edge[i]) 
	  {     
	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalScaledMult (p-1, 
			      lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
			      lam[e[0]]*lam[e[1]], 
                              SBLambda ([&](int i, Tx val)
                                        {
                                          shape[ii++] = Du (val);
                                        }));
	  }
      }

    // face shape functions
    for(int i = 0; i < N_FACE; i++) 
      if (order_face[i][0] >= 2)
        {
          IVec<4> fav = GetFaceSort (i, vnums);
          
          int vop = 6 - fav[0] - fav[1] - fav[2];  	
          int p = order_face[i][0];
          int pg = p - (type1 ? 1 : 0);
	  
          Tx xi = lam[fav[2]]-lam[fav[1]];
          Tx eta = lam[fav[0]]; // lo 
          Tx zeta = lam[vop];   // lz 
          
          TetShapesFaceLegendre::CalcSplitted (p+1, xi, eta, zeta, adpol1, adpol2); 
          
          // gradients 
          if (usegrad_face[i] && pg >= 2)
            
	    // DubinerBasis
            TrigOrthoPolGrad::EvalScaledMult
	      (pg-2, lam[fav[0]], lam[fav[1]],
	       1-lam[vop],
	       lam[fav[0]]*lam[fav[1]]*lam[fav[2]], 
	       SBLambda
	       ([&](int nr, Tx val)
		{
		  shape[ii++] = Du (val);
		}));
            

	  // non-gradient face shapes
          // if (type1) {
	  if (true) {

	    DubinerBasis::EvalMult
	      (p-2, lam[fav[0]], lam[fav[1]],
	       lam[fav[0]],
	       SBLambda
	       ([&](int nr, Tx val)
		{
		  shape[ii++] =	wuDv_minus_wvDu(lam[fav[1]], lam[fav[2]], val);
		}));
	    
            LegendrePolynomial::EvalMult
	      (p-2, lam[fav[2]]-lam[fav[1]], lam[fav[2]], 
	       SBLambda
	       ([&] (int j, Tx val)
		{
		  shape[ii++] = wuDv_minus_wvDu(lam[fav[1]], lam[fav[0]], val);
		}));
	  }
	  else {

	    // other combination
	    for (int j = 0; j <= p-2; j++)
	      for (int k = 0; k <= p-2-j; k++, ii++)
		shape[ii] = uDv_minus_vDu (adpol2[k], adpol1[j]);
          
	    // type 3
	    for (int j = 0; j <= p-2; j++, ii++)
	      shape[ii] = wuDv_minus_wvDu (lam[fav[1]], lam[fav[2]], adpol2[j]);
	  }
        }
    
    
    int p = order_cell[0];
    int pg = p - (type1 ? 1 : 0);

    // gradient inner shapes
	
    if (usegrad_cell && pg >= 3)

      DubinerBasis3DOrthoBub::EvalMult
	(pg-3, lam[0], lam[1], lam[2], lam[0]*lam[1]*lam[2]*lam[3], 
	 SBLambda
	 ([&](int nr, Tx val)
	  {
	    shape[ii++] = Du(val);
	  }));

    // non-gradient inner shapes, same for type 1 and type 2 Nedelec:

    if (p >= 3)  {

      Tx lam21 = lam[2]*lam[1], lam30 = lam[3]*lam[0];

      DubinerBasis3D::Eval
	(p-3, lam[0], lam[1], lam[2],
	 SBLambda
	 ([&](int nr, Tx val)
	  {
	    shape[ii++] = wuDv_minus_wvDu(lam[3], lam[0], lam21*val);
	    shape[ii++] = wuDv_minus_wvDu(lam[2], lam[1], lam30*val);
	  }));


      DubinerBasis::EvalMult
	(p-3, lam[0], lam[1], lam[0]*lam[1], 
	 SBLambda
	 ([&](int nr, Tx val)
	  {
	    shape[ii++] = wuDv_minus_wvDu(lam[3], lam[2], val);
	  }));
    }
  }


		        
  //------------------------------------------------------------------------
  //                   Prism
  //------------------------------------------------------------------------

  template<> template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE_Shape<ET_PRISM> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;

    // Tx x = hx[0], y = hx[1], z = hx[2];
    Tx x = ip.x, y = ip.y, z = ip.z;
    
    Tx lam[6] = { x, y, 1-x-y, x, y, 1-x-y };
    Tx muz[6]  = { 1-z, 1-z, 1-z, z, z, z };

    Tx sigma[6];
    for (int i = 0; i < 6; i++) sigma[i] = lam[i] + muz[i];

    ArrayMem<Tx,20> adpolxy1(order+3),adpolxy2(order+3); 
    ArrayMem<Tx,20> adpolz(order+3);   

    int ii = 9;
    
    // horizontal edge shapes
    for (int i = 0; i < 6; i++)
      {
	int p = order_edge[i]; 
        IVec<2> e = GetEdgeSort (i, vnums);	  
	
	//Nedelec0
        shape[i] = wuDv_minus_wvDu (lam[e[0]], lam[e[1]], muz[e[1]]);
   
	//high order \nabla (P_edge(x,y) * muz)
	if (p > 0 && usegrad_edge[i])
	  {
	    /*
              T_ORTHOPOL::CalcTrigExt(p+1, lam[e[1]]-lam[e[0]],
              1-lam[e[0]]-lam[e[1]],adpolxy1);
              for(int j = 0; j <= p-1; j++)
              shape[ii++] = Du<3> (adpolxy1[j] * muz[e[1]]);
	    */
            Tx xi = lam[e[1]]-lam[e[0]]; 
            Tx eta = lam[e[0]]+lam[e[1]]; 
            Tx bub = lam[e[0]]*lam[e[1]]*muz[e[1]];

	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalScaledMult (p-1, xi, eta, bub, adpolxy1);
	    for(int j = 0; j <= p-1; j++)
              shape[ii++] = Du (adpolxy1[j]);
	  }
      }

    //Vertical Edge Shapes
    for (int i = 6; i < 9; i++)
      {
	int p = order_edge[i]; 
        IVec<2> e = GetEdgeSort (i, vnums);	  

        shape[i] = wuDv_minus_wvDu (muz[e[0]], muz[e[1]], lam[e[1]]);
	
	//high order edges:  \nabla (T_ORTHOPOL^{p+1}(2z-1) * lam(x,y))
	if(p > 0 && usegrad_edge[i])
	  {
	    // T_ORTHOPOL::Calc (p+1, muz[e[1]]-muz[e[0]], adpolz);
	    // for (int j = 0; j < p; j++)
	    //   shape[ii++] = Du<3> (adpolz[j] * lam[e[1]]);

	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalMult (p-1, 
			muz[e[1]]-muz[e[0]], 
			muz[e[0]]*muz[e[1]]*lam[e[1]], adpolz);
	    
	    for (int j = 0; j < p; j++)
              shape[ii++] = Du (adpolz[j]);
	  }
      }


    const FACE * faces = ElementTopology::GetFaces (ET_PRISM).Data(); 

    // trig face shapes
    for (int i = 0; i < 2; i++)
      {
	int p = order_face[i][0];
	if (p < 2) continue;

	IVec<4> fav = GetFaceSort (i, vnums);

        {
          // gradients
          int pg = p - type1;
          if (usegrad_face[i] && pg >= 2)
            {
              TrigOrthoPolGrad::EvalMult (pg-2, lam[fav[0]], lam[fav[1]],
                                       lam[fav[0]]*lam[fav[1]]*lam[fav[2]]*muz[fav[2]], 
                                       SBLambda ([&](int nr, Tx val)
                                                 {
                                                   shape[ii++] = Du (val);
                                                 }));
            }
        }

        // if (type1) {        
        if (true) {
          
          DubinerBasis::EvalMult
            (p-2, lam[fav[0]], lam[fav[1]],
             lam[fav[0]]*muz[fav[2]],
             SBLambda
             ([&](int nr, Tx val)
             {
               shape[ii++] =  wuDv_minus_wvDu(lam[fav[1]], lam[fav[2]], val);
             }));
	    
          LegendrePolynomial::EvalMult
            (p-2, lam[fav[2]]-lam[fav[1]], lam[fav[2]]*muz[fav[2]], 
             SBLambda
             ([&] (int j, Tx val)
             {
               shape[ii++] = wuDv_minus_wvDu(lam[fav[1]], lam[fav[0]], val);
             }));
        }

        else

          {
            
            Tx xi = lam[fav[2]]-lam[fav[1]];
            Tx eta = lam[fav[0]]; // 1-lam[f2]-lam[f1];
            
            T_TRIGFACESHAPES::CalcSplitted(p+1,xi,eta,adpolxy1,adpolxy2); 
            /*
              if(usegrad_face[i])
              // gradient-fields =>  \nabla( adpolxy1*adpolxy2*muz )
              for (int j = 0; j <= p-2; j++)
              for (int k = 0; k <= p-2-j; k++)
              shape[ii++] = Du<3> (adpolxy1[j]*adpolxy2[k] * muz[fav[2]]);
        */
            
            
            // rotations of grad-fields => grad(uj)*vk*w -  uj*grad(vk)*w 
            for (int j = 0; j <= p-2; j++)
              for (int k = 0; k <= p-2-j; k++)
                shape[ii++] = wuDv_minus_wvDu (adpolxy2[k], adpolxy1[j], muz[fav[2]]);
            
            //  Ned0*adpolxy2[j]*muz 
            for (int j = 0; j <= p-2; j++,ii++)
              shape[ii] = wuDv_minus_wvDu (lam[fav[1]], lam[fav[2]], adpolxy2[j]*muz[fav[2]]);
          }
      }
    

    // quad faces: gradients of H1(p+1-type1), curl part of N_[p-1]
    for (int i = 2; i < 5; i++)
      {
	IVec<2> p = order_face[i];
        IVec<2> pg (p[0]-type1, p[1]-type1);
        IVec<4> f = GetFaceSort (i, vnums);	  

        if (usegrad_face[i] && pg[0] >= 1 && pg[1] >= 1)
          {
            Tx xi  = sigma[f[0]] - sigma[f[1]]; 
            Tx eta = sigma[f[0]] - sigma[f[3]];
        
            Tx scalexi(1.0), scaleeta(1.0);
            if (f[0] / 3 == f[1] / 3)  
              scalexi = lam[f[0]]+lam[f[1]];  // xi is horizontal
            else
              scaleeta = lam[f[0]]+lam[f[3]];
          
            Tx bub = (1.0/16)*(scaleeta*scaleeta-eta*eta)*(scalexi*scalexi-xi*xi);
            QuadOrthoPol::EvalScaled     (pg[0]-1, xi, scalexi, adpolxy1);
            QuadOrthoPol::EvalScaledMult (pg[1]-1, eta, scaleeta, bub, adpolz);

            for (int k = 0; k <= pg[0]-1; k++)
              for (int j = 0; j <= pg[1]-1; j++)
                shape[ii++] = Du (adpolxy1[k] * adpolz[j]);
          }

        if (p[0] < 1 || p[1] < 1) continue;
        IVec<2> q (p[0]-1, p[1]-1);

	int fmax = 0;
	for (int j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;

	int fz = 3-fmax; 
	int ftrig = fmax^1; 
	Tx xi = lam[faces[i][fmax]]-lam[faces[i][ftrig]]; 
	Tx eta = 1-lam[faces[i][fmax]]-lam[faces[i][ftrig]]; 
	Tx zeta = muz[faces[i][fmax]]-muz[faces[i][fz]]; 
	
	int pp = int(max2(q[0],q[1]))+1;
	T_ORTHOPOL::CalcTrigExt(pp,xi,eta,adpolxy1); 
	T_ORTHOPOL::Calc(pp,zeta,adpolz); 

	// Rotations of GradFields => nabla(polxy)*polz - polxy*nabla(polz)
	if (vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) 
	  for (int k = 0; k <= q[0]-1; k++)
	    for (int j = 0; j <= q[1]-1; j++)
              shape[ii++] = uDv_minus_vDu (adpolz[j], adpolxy1[k]);
	else
	  for (int j = 0; j <= q[0]-1; j++)
	    for (int k = 0; k <= q[1]-1; k++)
              shape[ii++] = uDv_minus_vDu (adpolxy1[k], adpolz[j]);
	
	// (ned0_trig)*polz, (ned0_quad)* polxy 
	if(vnums[faces[i][ftrig]] > vnums[faces[i][fz]]) // p = (p_trig,p_z) 
	  {
	    for(int j=0;j<=q[0]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu (muz[faces[i][fz]], muz[faces[i][fmax]], adpolxy1[j]);
	    for(int j=0;j<=q[1]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu (lam[faces[i][ftrig]], lam[faces[i][fmax]], adpolz[j]);
	  }
	else 
	  {
	    for(int j=0;j<=q[0]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu (lam[faces[i][ftrig]], lam[faces[i][fmax]], adpolz[j]);
	    for(int j=0;j<=q[1]-1;j++) 
              shape[ii++] = wuDv_minus_wvDu (muz[faces[i][fz]], muz[faces[i][fmax]], adpolxy1[j]);
	  }
      }

    // cell: gradients of H1(p+1-type1) bubbles, curl part of
    // Ned1_p(trig) x P_pz  +  P_p(trig) x P_{pz-1} e_z
    int p = order_cell[0], pz = order_cell[2];
    int pg = p-type1, pgz = pz-type1;
    if (usegrad_cell && pg >= 2 && pgz >= 1)
      {
        // same basis as H1HighOrderFE
        QuadOrthoPol::EvalMult (pgz-1, 2*z-1, z*(1-z), adpolz);
        TrigOrthoPolGrad::EvalMult (pg-2, x, y, x*y*(1-x-y),
                                SBLambda ([&](int nr, Tx val)
                                          {
                                            for (int k = 0; k < pgz; k++)
                                              shape[ii++] = Du (val * adpolz[k]);
                                          }));
      }

    if (p < 2 || pz < 1) return;
    Tx l2 = 1-x-y;
    // z-bubbles of P_pz
    if (pz >= 2)
      LegendrePolynomial::EvalMult (pz-2, 2*z-1, z*(1-z), adpolz);

    // trig curl bubbles x z-bubbles
    for (int k = 0; k < pz-1; k++)
      {
        Tx w = adpolz[k];
        DubinerBasis::EvalMult (p-2, x, y, x*w,
                                SBLambda ([&](int nr, Tx val)
                                          {
                                            shape[ii++] = wuDv_minus_wvDu (y, l2, val);
                                          }));
        LegendrePolynomial::EvalMult (p-2, l2-y, l2*w,
                                      SBLambda ([&](int j, Tx val)
                                                {
                                                  shape[ii++] = wuDv_minus_wvDu (y, x, val);
                                                }));
      }

    if (p < 3) return;
    // w grad u - u grad w, u trig bubbles of P_p
    DubinerBasis::EvalMult (p-3, x, y, x*y*l2,
                            SBLambda ([&](int nr, Tx u)
                                      {
                                        for (int k = 0; k < pz-1; k++)
                                          shape[ii++] = uDv_minus_vDu (adpolz[k], u);
                                      }));
    // u e_z
    DubinerBasis::EvalMult (p-3, x, y, x*y*l2,
                            SBLambda ([&](int nr, Tx u)
                                      {
                                        shape[ii++] = wuDv_minus_wvDu (z, 1-z, u);
                                      }));
  }



  //------------------------------------------------------------------------
  // HCurlHighOrderHex
  //------------------------------------------------------------------------


  template<> template<typename Tx, typename TFA>  
  void HCurlHighOrderFE_Shape<ET_HEX> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    // Tx x = hx[0], y = hx[1], z = hx[2];
    Tx x = ip.x, y = ip.y, z = ip.z;

    Tx lami[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
                (1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
    Tx sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
                 (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 
    
    int ii = 12; 
    ArrayMem<Tx, 20> pol_xi(order+2),pol_eta(order+2),pol_zeta(order+2);
   
    // edges
    for (int i = 0; i < 12; i++)
      {
	int p = order_edge[i]; 
        IVec<2> e = GetEdgeSort (i, vnums);	  
	
	Tx xi  = sigma[e[1]]-sigma[e[0]];
	Tx lam_e = lami[e[0]]+lami[e[1]];  
        Tx bub = 0.25 * lam_e * (1 - xi*xi);

	// Nedelec0-shapes
	shape[i] = uDv (0.5 * lam_e, xi); 

	// High Order edges ... Gradient fields 
	if(p > 0 && usegrad_edge[i])
	  {
	    //LegendrePolynomial::
            EdgeOrthoPol::
	      EvalMult (p-1, 
			xi, bub, pol_xi);

	    for (int j = 0; j < p; j++)
              shape[ii++] = Du (pol_xi[j]);
	  }
      }
    
    //Faces 
    const FACE * faces = ElementTopology::GetFaces (ET_HEX).Data();
    for (int i = 0; i<6; i++)
      {
	IVec<2> p = order_face[i];
        IVec<2> pg (p[0]-type1, p[1]-type1);

	Tx lam_f(0);
	for (int j = 0; j < 4; j++)
	  lam_f += lami[faces[i][j]];

        if (usegrad_face[i] && pg[0] >= 1 && pg[1] >= 1)
          {
            IVec<4> f = GetFaceSort (i, vnums);	  
            Tx xi  = sigma[f[0]] - sigma[f[1]]; 
            Tx eta = sigma[f[0]] - sigma[f[3]];
        
            Tx bub = lam_f*(1.0/16)*(1.0-eta*eta)*(1.0-xi*xi);
            QuadOrthoPol::Eval     (pg[0]-1, xi, pol_xi);
            QuadOrthoPol::EvalMult (pg[1]-1, eta, bub, pol_eta);
          
            for (int k = 0; k <= pg[0]-1; k++)
              for (int j = 0; j <= pg[1]-1; j++)
                shape[ii++] = Du (pol_xi[k] * pol_eta[j]);
          }

        // curl part of N_[p-1]
        if (p[0] < 1 || p[1] < 1) continue;
        IVec<2> q (p[0]-1, p[1]-1);
	
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
	      
	Tx xi = sigma[fmax]-sigma[f1]; 
	Tx eta = sigma[fmax]-sigma[f2]; 
    
	T_ORTHOPOL::Calc(q[0]+1, xi,pol_xi);
	T_ORTHOPOL::Calc(q[1]+1,eta,pol_eta);

	//Rotation of Gradient fields 
	for (int k = 0; k < q[0]; k++)
	  for (int j= 0; j < q[1]; j++)
            shape[ii++] = uDv_minus_vDu (pol_eta[j], lam_f * pol_xi[k]);

	// Missing ones 
	for(int j = 0; j < q[0];j++) 
          shape[ii++] = wuDv_minus_wvDu (Tx(0.5), eta, pol_xi[j]*lam_f); 

	for(int j = 0; j < q[1];j++) 
          shape[ii++] = wuDv_minus_wvDu (Tx(0.5), xi, pol_eta[j]*lam_f); 
      }


    
    IVec<3> p (order_cell[0], order_cell[1], order_cell[2]);
    IVec<3> pg (p[0]-type1, p[1]-type1, p[2]-type1);
    if (usegrad_cell && pg[0] >= 1 && pg[1] >= 1 && pg[2] >= 1)
      {
        QuadOrthoPol::EvalMult (pg[0]-1, 2*x-1, x*(1-x), pol_xi);
        QuadOrthoPol::EvalMult (pg[1]-1, 2*y-1, y*(1-y), pol_eta);
        QuadOrthoPol::EvalMult (pg[2]-1, 2*z-1, z*(1-z), pol_zeta);
          
        for (int i = 0; i < pg[0]; i++)
          for (int j = 0; j < pg[1]; j++)
            {
              Tx pxy = pol_xi[i] * pol_eta[j];
              for (int k = 0; k < pg[2]; k++)
                shape[ii++] = Du (pxy * pol_zeta[k]);
            }
      }

    // curl part of N_[p-1]
    if (p[0] < 1 || p[1] < 1 || p[2] < 1) return;
    IVec<3> q (p[0]-1, p[1]-1, p[2]-1);
    T_ORTHOPOL::Calc(q[0]+1,2*x-1,pol_xi);
    T_ORTHOPOL::Calc(q[1]+1,2*y-1,pol_eta);
    T_ORTHOPOL::Calc(q[2]+1,2*z-1,pol_zeta); 
    
    //Rotations of gradient fields
    for (int i=0; i<q[0]; i++)
      for(int j=0; j<q[1]; j++) 
	for(int k=0; k<q[2]; k++)
	  {
            shape[ii++] = uDv_minus_vDu (pol_xi[i] * pol_eta[j], pol_zeta[k]);
            shape[ii++] = uDv_minus_vDu (pol_xi[i], pol_eta[j] * pol_zeta[k]);
	  } 
    
    for(int i = 0; i < q[0]; i++) 
      for(int j = 0; j < q[1]; j++)
        shape[ii++] = wuDv_minus_wvDu (z,1-z,pol_xi[i] * pol_eta[j]);

    for(int i = 0; i < q[0]; i++) 
      for(int k = 0; k < q[2]; k++)
        shape[ii++] = wuDv_minus_wvDu (y,1-y,pol_xi[i] * pol_zeta[k]);

    for(int j = 0; j < q[1]; j++)
      for(int k = 0; k < q[2]; k++)
        shape[ii++] = wuDv_minus_wvDu (x,1-x,pol_eta[j] * pol_zeta[k]);
  }



  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_HEX> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;        
    auto & ip = mip.IP();
    T x = ip(0), y = ip(1), z = ip(2);
    T sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
                (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z};
    Vec<3> pnts[8] = 
      { 
	{ 0, 0, 0 },
	{ 1, 0, 0 },
	{ 1, 1, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 1 },
	{ 1, 0, 1 },
	{ 1, 1, 1 },
	{ 0, 1, 1 }
      };
    
    int facetnr = ip.FacetNr();
    int ii = 12;

    if (ip.VB() == BBND)
      { // edge shapes
        for (int i = 0; i < 12; i++)
          {
            int p = order_edge[i] * usegrad_edge[i];
            if (i == facetnr)
              {
                IVec<2> e = GetEdgeSort (i, vnums);
                T xi = sigma[e[1]]-sigma[e[0]];
                Vec<3> tauref = pnts[e[1]] - pnts[e[0]];
                Vec<3,T> tau = mip.GetJacobian()*tauref;
                tau /= mip.GetMeasure();
                LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              Vec<3,T> vshape = val * tau;
                              if (nr==0)
                                shape[i] = vshape;
                              else
                                shape[ii+nr-1] = vshape;
                            }));
              }
            ii += p;
          }
        return;
      }

    for (int i = 0; i < 12; i++)
      ii += order_edge[i] * usegrad_edge[i];

    if (ip.VB() == BND)
      {
        Vec<3,T> X(x, y, z);
        for (int i = 0; i < 6; i++)
          {
            if (i == facetnr)
              {
                IVec<4> f = GetFaceSort (i, vnums);
                T s, t;
                HCurlFaceCoords (X, pnts[f[0]], pnts[f[1]], pnts[f[3]], s, t);
                auto ftrafo = HCurlFaceTrafo<T> (mip, pnts[f[1]]-pnts[f[0]], pnts[f[3]]-pnts[f[0]]);
                HCurlQuadDual (ftrafo, s, t, IVec<2>(order_face[i]), usegrad_face[i], type1, ii, shape);
                return;
              }
            ii += hcurl_ndof::Face (ET_QUAD, order_face[i], usegrad_face[i], type1);
          }
        return;
      }

    for (int i = 0; i < 6; i++)
      ii += hcurl_ndof::Face (ET_QUAD, order_face[i], usegrad_face[i], type1);

    if (ip.VB() == VOL)
      {
        auto trafo = 1/mip.GetMeasure()*mip.GetJacobian();
        IVec<3> p (order_cell[0], order_cell[1], order_cell[2]);
        int n = max(p[0], max(p[1], p[2])) + 1;
        ArrayMem<AutoDiff<1,T>,20> lx(n+1), ly(n+1), lz(n+1);
        HCurlLegendreAD (n, x, lx);
        HCurlLegendreAD (n, y, ly);
        HCurlLegendreAD (n, z, lz);

        // gradients: (F_i L_j L_k, 0, 0), div onto Q_{p-1-type1}
        IVec<3> pg (p[0]-type1, p[1]-type1, p[2]-type1);
        if (usegrad_cell && pg[0] >= 1 && pg[1] >= 1 && pg[2] >= 1 && !type1)
          // grad of Q_{p+1} bubbles: blind to div-free fields
          for (int i = 0; i < pg[0]; i++)
            {
              T bi = HCurlIntLeg<T> (lx, i+1), di = lx[i+1].Value();
              for (int j = 0; j < pg[1]; j++)
                {
                  T bj = HCurlIntLeg<T> (ly, j+1), dj = ly[j+1].Value();
                  for (int k = 0; k < pg[2]; k++)
                    {
                      T bk = HCurlIntLeg<T> (lz, k+1), dk = lz[k+1].Value();
                      shape[ii++] = trafo * Vec<3,T>(di*bj*bk, bi*dj*bk, bi*bj*dk);
                    }
                }
            }
        else if (usegrad_cell && pg[0] >= 1 && pg[1] >= 1 && pg[2] >= 1)
          for (int i = 0; i < pg[0]; i++)
            {
              T fi = HCurlIntLeg<T> (lx, i);
              for (int j = 0; j < pg[1]; j++)
                for (int k = 0; k < pg[2]; k++)
                  shape[ii++] = trafo * Vec<3,T>(fi*ly[j].Value()*lz[k].Value(), 0, 0);
            }

        // curl part of N_[q], q = p-1: div-free subspace of RT_[q-1]
        if (p[0] < 1 || p[1] < 1 || p[2] < 1) return;
        IVec<3> q (p[0]-1, p[1]-1, p[2]-1);
        // M_0 = 1, M_a = F_{a-1}, M_a' = L_{a-1}
        auto M = [&] (FlatArray<AutoDiff<1,T>> l, int a) -> T
          { return a == 0 ? T(1.0) : HCurlIntLeg<T> (l, a-1); };
        auto dM = [&] (FlatArray<AutoDiff<1,T>> l, int a) -> T
          { return a == 0 ? T(0.0) : l[a-1].Value(); };

        // curl (phi e_z), phi = M_a(x) M_b(y) L_c(z), (a,b) != 0
        for (int a = 0; a <= q[0]; a++)
          for (int b = 0; b <= q[1]; b++)
            if (a || b)
              {
                T mxy = M(lx,a)*dM(ly,b), dmxy = dM(lx,a)*M(ly,b);
                for (int c = 0; c < q[2]; c++)
                  shape[ii++] = trafo * Vec<3,T>(mxy*lz[c].Value(), -dmxy*lz[c].Value(), 0);
              }
        // curl (psi e_x), psi = L_i(x) M_j(y) M_l(z), j >= 1
        for (int i = 0; i < q[0]; i++)
          for (int j = 1; j <= q[1]; j++)
            for (int l = 0; l <= q[2]; l++)
              {
                T li = lx[i].Value();
                shape[ii++] = trafo * Vec<3,T>(0, li*M(ly,j)*dM(lz,l), -li*dM(ly,j)*M(lz,l));
              }
      }
  }








  

  //------------------------------------------------------------------------
  //            Pyramid
  //------------------------------------------------------------------------


  template<> template<typename Tx, typename TFA>  
  void  HCurlHighOrderFE_Shape<ET_PYRAMID> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;  

    // Tx x = hx[0], y = hx[1], z = hx[2];
    Tx x = ip.x, y = ip.y, z = ip.z;

    z = ClampPyramidApex (z);   // keeps 1-z > 0 at the apex, exact elsewhere
    Tx xt = x/(1-z);
    Tx yt = y/(1-z);
    Tx sigma[5] = {(1-xt)+(1-yt)+(1-z),xt+(1-yt)+(1-z), xt + yt + (1-z), 
			    (1-xt)+yt+(1-z),z}; 

    Tx lami[5] = {(1-xt)*(1-yt)*(1-z),xt*(1-yt)*(1-z), xt * yt * (1-z), 
			   (1-xt)*yt*(1-z),z}; 

    Tx lambda[5] = {(1-xt)*(1-yt),xt*(1-yt), xt * yt, 
                             (1-xt)*yt,z}; 
        
       
    ArrayMem<Tx, 20> pol_xi(order+2), pol_eta(order+2), pol_zeta(order+2);
    
    int ii =8; 
 
    // horizontal edges incl. Nedelec 0
    for (int i = 0; i < 4; i++)
      {
        int p = order_edge[i];
        IVec<2> e = GetEdgeSort (i, vnums);	  
	
	Tx xi  = sigma[e[1]] - sigma[e[0]];   
	Tx lam_t = lambda[e[1]] + lambda[e[0]]; 
        
        shape[i] = uDv (0.5 * (1-z)*(1-z)*lam_t, xi);

	if(p > 0 && usegrad_edge[i])
	  {
            Tx bub = 0.25*(1-xi*xi)*(1-z)*(1-z)*lam_t;
	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalScaledMult (p-1,
			      xi*(1-z), 1-z, bub,
                              SBLambda ([&](int i, Tx val)
                                        {
                                          shape[ii++] = Du(val);
                                        }));
	  }
      }
    
    // vertical edges incl. Nedelec 0  
    for(int i = 4; i < 8; i++)
      {
        int p = order_edge[i];
        IVec<2> e = GetEdgeSort (i, vnums);	  

        shape[i] = uDv_minus_vDu (lami[e[0]], lami[e[1]]);

	if (p > 0 && usegrad_edge[i])
	  {
            Tx xi = lami[e[1]]-lami[e[0]]; 
            Tx lam_e = lami[e[0]]+lami[e[1]];
            Tx bub = 0.25 * (lam_e*lam_e-xi*xi);

	    // LegendrePolynomial::
            EdgeOrthoPol::
	      EvalScaledMult (p-1,
                              xi, lam_e, bub, 
                              SBLambda ([&](int i, Tx val)
                                        {
                                          shape[ii++] = Du(val);
                                        }));
	  }
      }

    const FACE * faces = ElementTopology::GetFaces (ET_PYRAMID).Data(); 

    // trig face dofs
    for (int i = 0; i < 4; i++)
      if (order_face[i][0] >= 2)
	{
	  int p = order_face[i][0];
	  Tx lam_face = lambda[faces[i][0]] + lambda[faces[i][1]];  
	  Tx bary[3] = 
	    {(sigma[faces[i][0]]-(1-z)-lam_face)*(1-z), 
	     (sigma[faces[i][1]]-(1-z)-lam_face)*(1-z), z}; 
			     
	  int fav[3] = {0, 1, 2};
	  if(vnums[faces[i][fav[0]]] > vnums[faces[i][fav[1]]]) swap(fav[0],fav[1]); 
	  if(vnums[faces[i][fav[1]]] > vnums[faces[i][fav[2]]]) swap(fav[1],fav[2]);
	  if(vnums[faces[i][fav[0]]] > vnums[faces[i][fav[1]]]) swap(fav[0],fav[1]); 	

	  int pg = p - type1;
	  if (usegrad_face[i] && pg >= 2)
            {
              Tx bub = lam_face * bary[fav[0]]*bary[fav[1]]*bary[fav[2]];
              TrigOrthoPolGrad::
                EvalMult (pg-2, bary[fav[0]], bary[fav[1]], bub, 
                          SBLambda ([&](int nr, Tx val)
                                    {
                                      shape[ii++] = Du (val);
                                    }));
            }
          
          /*
	  // phi = pol_xi * pol_eta * lam_face; 
	  // Type 1: Gradient Functions 
	  if(usegrad_face[i])
	    for(int j=0;j<= p-2; j++)
	      for(int k=0;k<=p-2-j; k++)
                shape[ii++] = Du<3> (pol_xi[j] * pol_eta[k]);
          */

          // if (type1) {          
          if (true) {
            
            DubinerBasis::EvalMult
              (p-2, bary[fav[0]], bary[fav[1]],
               bary[fav[0]]*lam_face,
               SBLambda
               ([&](int nr, Tx val)
               {
                 shape[ii++] =  wuDv_minus_wvDu(bary[fav[1]], bary[fav[2]], val);
               }));
	    
            LegendrePolynomial::EvalMult
              (p-2, bary[fav[2]]-bary[fav[1]], bary[fav[2]]*lam_face, 
               SBLambda
               ([&] (int j, Tx val)
               {
                 shape[ii++] = wuDv_minus_wvDu(bary[fav[1]], bary[fav[0]], val);
             }));
          }

          else {
            T_TRIGFACESHAPES::CalcSplitted(p+1, bary[fav[2]]-bary[fav[1]], 
                                           bary[fav[0]],pol_xi,pol_eta);
            
            for(int j=0;j<=p-2;j++) pol_eta[j] *= lam_face;
            
            // Type 2:  
            for(int j=0;j<= p-2; j++)
              for(int k=0;k<=p-2-j; k++)
                shape[ii++] = uDv_minus_vDu (pol_eta[k], pol_xi[j]);
            
            // Type 3: Nedelec-based ones (Ned_0*v_j)
            for(int j=0;j<=p-2;j++)
              shape[ii++] = wuDv_minus_wvDu (bary[fav[1]], bary[fav[2]], pol_eta[j]);
          }
	}


    // quad face 
    if (order_face[4][0] >= 1)
      {
	IVec<4> f = GetFaceSort (4, vnums);	  
	Tx xi  = sigma[f[0]] - sigma[f[1]]; 
	Tx eta = sigma[f[0]] - sigma[f[3]];

        int pg = order_face[4][0] - type1;
        if (usegrad_face[4] && pg >= 1)
          {
            // gradients of H1(pg+1)
            Tx fac(1.0);
            for (int k = 1; k <= pg+1; k++) fac *= (1-z);
            QuadOrthoPol::
              EvalMult (pg-1, xi, fac*0.25*(1-xi*xi), pol_xi);
            QuadOrthoPol::
              EvalMult (pg-1, eta, 0.25*(1-eta*eta), pol_eta);

            for (int k = 0; k <= pg-1; k++) 
              for (int j = 0; j <= pg-1; j++, ii++) 
                shape[ii] = Du (pol_xi[k] * pol_eta[j]);
          }

        // curl part of N_[p-1]
	int px = order_face[4][0]-1;
	int py = order_face[4][0]-1;
	Tx fac(1.0);
	for (int k = 1; k <= px+1; k++) fac *= (1-z);
        
	int fmax = 0;
	for (int l=1; l<4; l++) 
	  if (vnums[l] > vnums[fmax]) fmax = l;  

	int f1 = (fmax+3)%4;
	int f2 = (fmax+1)%4; 
	if(vnums[f1]>vnums[f2]) swap(f1,f2);  // fmax > f2 > f1 

 	xi  = sigma[fmax] - sigma[f2]; 
	eta = sigma[fmax] - sigma[f1];

	T_ORTHOPOL::Calc (px+1, xi, pol_xi);	
	T_ORTHOPOL::Calc (py+1, eta, pol_eta);
	
	for(int k = 0; k < py; k++) pol_eta[k] *= fac;

	// Type 2: 
	for (int k = 0; k < px; k++) 
	  for (int j = 0; j < py; j++) 
            shape[ii++] = uDv_minus_vDu (pol_eta[j], pol_xi[k]);

	// Type 3:
	for (int k = 0; k < px; k++)
          shape[ii++] = uDv (0.5*pol_xi[k]*fac, eta);

	for (int k = 0; k < py; k++)
          shape[ii++] = uDv (0.5*pol_eta[k] /* *fac */, xi); 
      }

    // gradients of the H1(pg+1) cell bubbles, same basis as H1HighOrderFE
    if (int pg = order_cell[0]-type1; usegrad_cell && pg >= 2)
      {
        int ph = pg+1;
        QuadOrthoPol::EvalMult (ph-3, 2*xt-1, xt*(1-xt), pol_xi);
        QuadOrthoPol::EvalMult (ph-3, 2*yt-1, yt*(1-yt), pol_eta);

        Tx bubz = z*(1-z);
        for (int m = 0; m <= ph-3; m++)
          {
            bubz *= 1-z;
            JacobiPolynomial2(m+4, 1).EvalMult (ph-3-m, 2*z-1, bubz, pol_zeta);
            for (int i = 0; i <= m; i++)
              for (int j = 0; j <= m; j++)
                if (i == m || j == m)
                  {
                    Tx pxy = pol_xi[i] * pol_eta[j];
                    for (int n = 0; n <= ph-3-m; n++)
                      shape[ii++] = Du (pxy * pol_zeta[n]);
                  }
          }
      }

    if (order_cell[0] >= 2)
      {
	int pp = order_cell[0];
	// According H^1 terms: 
	// u_i = L_i+2(2xt-1)
	// v_j = L_j+2(2yt-1) 
	// w_k = z * (1-z)^(k+2)  with 0 <= i,j <= k, 0<= k <= p-2  
	
	LegendrePolynomial::EvalMult (pp-1, 2*xt-1, xt*(1-xt), pol_xi);
        LegendrePolynomial::EvalMult (pp-1, 2*yt-1, yt*(1-yt), pol_eta);
	// T_ORTHOPOL::Calc (pp+3, 2*xt-1, pol_xi);
	// T_ORTHOPOL::Calc (pp+3, 2*yt-1, pol_eta);		
		
	pol_zeta[0] = z*(1-z)*(1-z);
	for (int k=1;k<=pp-2;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];

	// Type 2a: l.i. combinations of grad-terms   
	// shape = u_i \nabla(v_j) w_k 
	// shape = u_i v_j \nabla(w_k) 
	for(int k=0;k<= pp-2;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
              shape[ii] = uDv (pol_xi[i]*pol_zeta[k], pol_eta[j]);

	// Type 2b: shape = v_j w_k \nabla (xt) 
	//          shape = u_i w_k \nabla (yt)
	for(int k = 0;k<= pp-2;k++)
	  for(int j=0;j<=k;j++) 
            shape[ii++] = uDv (pol_eta[j]*pol_zeta[k], xt);
	
	for(int  k = 0;k<= pp-2;k++)
	  for (int i=0;i<=k;i++)
            shape[ii++] = uDv (pol_xi[i]*pol_zeta[k], yt);
	
	// 3rd component spans xi^i eta^j zeta^(k-1), i,j <= k
	// pol_zeta starts linear in zeta 
	// pol_xi and pol_eta quadratic in xi resp. eta 
	pol_zeta[0] = (1-z);
	for (int k=1;k<=pp;k++) 
	  pol_zeta[k] = (1-z)*pol_zeta[k-1];
     
	for(int k=0;k<= pp-1;k++)
	  for(int i=0;i<=k;i++)
	    for(int j=0;j<=k;j++,ii++)
              shape[ii] = uDv (pol_eta[j] * pol_xi[i] * pol_zeta[k], z);
      }
  }




  // dual shapes

  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_TRIG> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    // shape = 0;
    auto & ip = mip.IP();
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
    T x = ip(0), y = ip(1);
    T lam[3] = { x, y, 1-x-y };
    Vec<2,T> pnts[3] = { { 1, 0 }, { 0, 1 } , { 0, 0 } };
    int facetnr = ip.FacetNr();

    if (ip.VB() == BND)
      { // facet shapes
        int ii = 3;
        for (int i = 0; i < 3; i++)
          {
            int p = order_edge[i] * usegrad_edge[i];
            if (i == facetnr)
              {
                IVec<2> e = GetEdgeSort (i, vnums);
                T xi = lam[e[1]]-lam[e[0]];
                Vec<2,T> tauref = pnts[e[1]] - pnts[e[0]];
                auto tau = mip.GetJacobian()*tauref;
                tau /= mip.GetMeasure();

                LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              auto vshape = val * tau;
                              if (nr==0)
                                shape[i] = vshape;
                              else
                                shape[ii+nr-1] = vshape;
                            }));
              }
            ii += p;
          }
      }
    if (ip.VB() == VOL)
      { // inner shapes
        int ii = 3;
        for (int i = 0; i < 3; i++)
          ii += order_edge[i] * usegrad_edge[i];


        auto trafo = 1/mip.GetMeasure()*mip.GetJacobian();

        if (biorth)
          {
            // duals B_kl, C_kl of Haubold-Beuchler-Schoeberl (SINUM 2024), (4.12),
            // in collapsed coords xi = (lb-la)/(la+lb), yc = 2lc-1 on (-1,1)^2
            int p = order_face[0][0];
            if (p < 2) return;
            IVec<4> fav = GetFaceSort (0, vnums);
            T lc = lam[fav[0]], lb = lam[fav[1]], la = lam[fav[2]];
            T xs = lb-la, s = la+lb, yc = 2*lc-1;
            Vec<2,T> tex = trafo * Vec<2,T>(0.5 * (pnts[fav[1]]-pnts[fav[2]]));
            Vec<2,T> tey = trafo * Vec<2,T>(0.5 * (pnts[fav[0]]-0.5*(pnts[fav[1]]+pnts[fav[2]])));

            auto & invtab = biorth_tabs.inv;

            // writing grad duals and w duals at their own offsets
            int pg = (usegrad_face[0]) ? p-2-(type1 ? 1 : 0) : -1;
            int ig = ii, iw = ii + (pg+1)*(pg+2)/2;

            // scaled Ls_i and d/dxs on the fly, step i = m+2 uses Ls_{i-2..i}
            typedef AutoDiff<1,T> AD1;
            AD1 Lm1(0.0), Lm2(0.0);
            LegendrePolynomial::EvalScaled1Assign
              (p, AD1(xs,0), AD1(s), SBLambda([&](int i, AD1 Li) LAMBDA_INLINE
              {
                if (i >= 2)
                  {
                    int m = i-2;
                    // B = q1 b, C = q2 c
                    Vec<2,T> b = Lm1.Value() * tex;
                    Vec<2,T> c = (2.0*Li.DValue(0)) * tex - (4.0*s*Lm1.DValue(0)) * tey;
                    Vec<2,T> bg = double(i) * b;

                    // w duals: nrm rot (Ls_{i-1}(xi) Q(yc)), Q = P^(2i-1,0)_n
                    Vec<2,T> dq = (2.0*s*Lm2.DValue(0)) * tex - (4.0*Lm1.DValue(0)) * tey;
                    Vec<2,T> dqp = (4.0*Lm1.Value()) * tex;

                    // P = P^(2i,1)_n = (n+1) IntJac^(2i+1)_n, rest by index shifts
                    T Pm1(0.0), Pam1(0.0);
                    IntegratedJacobiPolynomialAlpha(2*i+1).Eval
                      (p-2-m, yc, SBLambda([&](int n, T R) LAMBDA_INLINE
                      {
                        T P = (n+1)*R;
                        double inv = invtab[2*n+2*i+1];
                        T Pa = inv * ((n+2*i+1)*P - (n+1)*Pm1);   // P^(2i-1,1)_n
                        if (n <= pg-m)
                          {
                            T Pb = inv * ((n+2*i+1)*P + (n+2*i)*Pm1);   // P^(2i,0)_n
                            T a2 = ((n+1)*0.5*invtab[n+2*i]) * Pb;
                            shape[ig++] = Vec<2,T>(Pa*bg(0)+a2*c(0), Pa*bg(1)+a2*c(1));
                          }
                        double nrm = biorth_tabs.sq[2*i-1]*biorth_tabs.sq[2*n+2*i];
                        T Q = (nrm*0.5*invtab[n+i]) * ((n+2*i)*Pa + (n+2*i-1)*Pam1);
                        T dQ = (nrm*0.5*(n+2*i)) * Pm1;   // Q' = (n+2i)/2 P^(2i,1)_{n-1}
                        shape[iw++] = Vec<2,T>(Q*dq(0)+dQ*dqp(0), Q*dq(1)+dQ*dqp(1));
                        Pm1 = P;
                        Pam1 = Pa;
                      }));
                  }
                Lm2 = Lm1;
                Lm1 = Li;
              }));
            ii = iw;

            // (j+1) IntJac^(3)_j = P^(2,1)_j
            IntegratedJacobiPolynomialAlpha(3).Eval
              (p-2, yc, SBLambda([&](int j, T R) LAMBDA_INLINE
              {
                shape[ii++] = ((j+1)*biorth_tabs.isq[2*j+4]*R) * tex;
              }));
            return;
          }

        if (!usegrad_face[0])
          {
            // no gradients: rot of P_{p-1}, annihilates the gradient bubbles
            int p = order_face[0][0];
            if (p >= 2)
              DubinerBasis::Eval(p-1, AutoDiff<2,T>(x, 0), AutoDiff<2,T>(y, 1),
                                 SBLambda([&] (size_t nr, auto val)
                                          {
                                            if (nr == 0) return;
                                            shape[ii++] = trafo * Vec<2,T> (val.DValue(1), -val.DValue(0));
                                          }));
            return;
          }
        DubinerBasis::Eval(order-2, x, y,
                           SBLambda([&] (size_t nr, auto val)
                                    {
                                      shape[ii++] = trafo * Vec<2,T> (val, 0);
                                      if (type1)
                                        shape[ii++] = trafo * Vec<2,T> (0, val);
                                      else
                                        shape[ii++] = trafo * Vec<2,T> (val*x, val*y);
                                    }));
        if (!type1)   // P^(1,0)(2x-1): orthogonal on the trig, Legendre(x) is not
          JacobiPolynomialFix<1,0>::Eval(order-2, 2*x-1,
                                         SBLambda([&] (size_t nr, auto val)
                                                  {
                                                    shape[ii++] = trafo * Vec<2,T>(0,val);
                                                  }));
      }
  }

  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_QUAD> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    auto & ip = mip.IP();
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
    T x = ip(0), y = ip(1);
    Vec<2,T> pnts[4] = { { 0, 0 },
	{ 1, 0 },
	{ 1, 1 },
	{ 0, 1 } };
    T sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    int facetnr = ip.FacetNr();

    if (ip.VB() == BND)
      { // facet shapes
        int ii = 4;
        for (int i = 0; i < 4; i++)
          {
            int p = order_edge[i] * usegrad_edge[i];
            if (i == facetnr)
              {
                IVec<2> e = GetEdgeSort (i, vnums);
                T xi = sigma[e[1]]-sigma[e[0]];
                Vec<2,T> tauref = pnts[e[1]] - pnts[e[0]];
                auto tau = mip.GetJacobian()*tauref;
                tau /= mip.GetMeasure();

                LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              auto vshape = val * tau;
                              if (nr==0)
                                shape[i] = vshape;
                              else
                                shape[ii+nr-1] = vshape;
                            }));
              }
            ii += p;
          }
      }
    if (ip.VB() == VOL)
      { // inner shapes
        int ii = 4;
        for (int i = 0; i < 4; i++)
          ii += order_edge[i] * usegrad_edge[i];

        auto trafo = 1/mip.GetMeasure()*mip.GetJacobian();
        HCurlQuadDual (trafo, x, y, order_face[0], usegrad_face[0], type1, ii, shape);
      }
  }

  

  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_TET> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    // shape = 0;
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;        
    auto & ip = mip.IP();
    T x = ip(0), y = ip(1), z = ip(2);
    T lam[4] = { x, y, z, 1-x-y-z };
    Vec<3> pnts[4] = { { 1, 0, 0 }, { 0, 1, 0 } , { 0, 0, 1 }, { 0, 0, 0 } };
    int facetnr = ip.FacetNr();
    int ii = 6;

    if (ip.VB() == BBND)
      { // edge shapes
        for (int i = 0; i < 6; i++)
          {
            int p = order_edge[i] * usegrad_edge[i];
            if (i == facetnr)
              {
                IVec<2> e = GetEdgeSort (i, vnums);
                T xi = lam[e[1]]-lam[e[0]];
                Vec<3> tauref = pnts[e[1]] - pnts[e[0]];
                Vec<3,T> tau = mip.GetJacobian()*tauref;
                tau /= mip.GetMeasure();
                LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              Vec<3,T> vshape = val * tau;
                              if (nr==0)
                                shape[i] = vshape;
                              else
                                shape[ii+nr-1] = vshape;
                            }));
              }
            ii += p;
          }
      }
    else
      {
        for (int i = 0; i < 6; i++)
          ii += order_edge[i] * usegrad_edge[i];
      }
    if (ip.VB() == BND)
      {
	//AutoDiff<3,T> xa(ip(0), 0), ya(ip(1),1), za(ip(2),2);
	//AutoDiff<3,T> lami[4] = { xa, ya, za, (T)(1.0) };

	for (int f = 0; f < 4; f++)
	  {
	    int p = order_face[f][0];
	     if (f == facetnr)
	       {
		 IVec<4> fav = GetFaceSort (facetnr, vnums);
		 //AutoDiff<3,T> adxi = lami[fav[0]]-lami[fav[2]];
		 //AutoDiff<3,T> adeta = lami[fav[1]]-lami[fav[2]];
		 Vec<3> adxi = pnts[fav[0]] - pnts[fav[2]];
		 Vec<3> adeta = pnts[fav[1]] - pnts[fav[2]];
		 T xi = lam[fav[0]];
		 T eta = lam[fav[1]];
		 
		 Mat<3,2> F;
		 for (int k = 0; k < 3; k++)
		   { F(k,0) = adxi(k); F(k,1) = adeta(k); }
		 Mat<2,2> Ftmp = Trans(F)*F;
		 double det = sqrt(Ftmp(0,0)*Ftmp(1,1)-Ftmp(1,0)*Ftmp(0,1));

		 // face map, once per point: 1/(det*measure) J F
		 Mat<3,3,T> jac = mip.GetJacobian();
		 T scal = 1/(det*mip.GetMeasure());
		 Mat<3,2,T> ftrafo;
		 for (int i = 0; i < 3; i++)
		   for (int j = 0; j < 2; j++)
		     ftrafo(i,j) = scal * (jac(i,0)*F(0,j) + jac(i,1)*F(1,j) + jac(i,2)*F(2,j));

		 if (!usegrad_face[f])
		   {
		     // no gradients: surface rot of P_{p-1}, annihilates the face gradients
		     if (p >= 2)
		       DubinerBasis::Eval(p-1, AutoDiff<2,T>(xi, 0), AutoDiff<2,T>(eta, 1),
					  SBLambda([&] (size_t nr, auto val)
						   {
						     if (nr == 0) return;
						     shape[ii++] = ftrafo*Vec<2,T> (val.DValue(1), -val.DValue(0));
						   }));
		     continue;
		   }
		 
		 DubinerBasis::Eval(order-2, xi, eta,
                                    SBLambda([&] (size_t nr, auto val)
                                             {						 
                                               shape[ii++] = ftrafo*Vec<2,T> (val, 0);
					       if (type1)
						 shape[ii++] = ftrafo*Vec<2,T>(0, val);
					       else
						 shape[ii++] = ftrafo*Vec<2,T> (val*xi, val*eta);
                                             }));
		 if (!type1)   // P^(1,0)(2xi-1), orthogonal on the face
		   JacobiPolynomialFix<1,0>::Eval(order-2, 2*xi-1,
					  SBLambda([&] (size_t nr, auto val)
						   {
						     shape[ii++] = ftrafo*Vec<2,T>(0, val);
						   }));
	       }
	     else
	       ii += usegrad_face[f] ? (p+!type1)*(p-1) : (p+2)*(p-1)/2;
	     
	  }
      }
    else
      {
        for (int i = 0; i < 4; i++)
	  {
	    int p = order_face[i][0];
	    ii += usegrad_face[i] ? (p+!type1)*(p-1) : (p+2)*(p-1)/2;
	  }
      }
    if (ip.VB() == VOL)
      {
	// auto xphys = mip.GetPoint()(0);
        // auto yphys = mip.GetPoint()(1);
	// auto zphys = mip.GetPoint()(2);

	auto trafo = 1/mip.GetMeasure()*mip.GetJacobian();

	if (!usegrad_cell)
	  {
	    // no gradients: div-free P_{p-3}^3, curls of (0,x*a,0), (0,0,x*a), (0,z*b(y,z),0)
	    int p = order_cell[0];
	    AutoDiff<3,T> adx(x, 0), ady(y, 1), adz(z, 2);
	    if (p >= 3)
	      {
		DubinerBasis3D::Eval (p-3, adx, ady, adz,
				      SBLambda([&](size_t j, auto val)
					       {
						 auto f = adx*val;
						 shape[ii++] = trafo*Vec<3,T>(-f.DValue(2), 0, f.DValue(0));
						 shape[ii++] = trafo*Vec<3,T>(f.DValue(1), -f.DValue(0), 0);
					       }));
		DubinerBasis::Eval (p-3, ady, adz,
				    SBLambda([&](size_t j, auto val)
					     {
					       shape[ii++] = trafo*Vec<3,T>(-(adz*val).DValue(2), 0, 0);
					     }));
	      }
	    return;
	  }

	DubinerBasis3D::Eval (order-3, lam[0], lam[1], lam[2],
			      SBLambda([&](size_t j, T val) LAMBDA_INLINE
	 						     {
							       if (type1)
								 shape[ii++] = trafo*Vec<3,T> (0, 0, val);
							       else
								 shape[ii++] = trafo*Vec<3,T>(val*x, val*y, val*z);
	 						       shape[ii++] = trafo*Vec<3,T>(val, 0, 0);
	 						       shape[ii++] = trafo*Vec<3,T>(0, val, 0);
	 						     }));
	
	if (!type1)
	  DubinerBasis::Eval(order-3, x, y,
                           SBLambda([&] (size_t nr, auto val)
                                    {
                                      shape[ii++] = trafo*Vec<3,T> (0, 0, val);
                                    }));
      }
  }

  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_PRISM> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;
    auto & ip = mip.IP();
    T x = ip(0), y = ip(1), z = ip(2);
    T lam[6] = { x, y, 1-x-y, x, y, 1-x-y };
    T muz[6]  = { 1-z, 1-z, 1-z, z, z, z };

    T sigma[6];
    for (int i = 0; i < 6; i++) sigma[i] = lam[i] + muz[i];

    Vec<3> pnts[6] = { { 1, 0, 0 }, { 0, 1, 0 } , { 0, 0, 0 },
                       { 1, 0, 1 }, { 0, 1, 1 }, { 0, 0, 1 }};
    int facetnr = ip.FacetNr();
    int ii = 9;

    if (ip.VB() == BBND)
      {
        for (int i = 0; i < 9; i++)
          {
            int p = order_edge[i] * usegrad_edge[i];
            if (i == facetnr)
              {
                IVec<2> e = GetEdgeSort (i, vnums);
                T xi = sigma[e[1]] - sigma[e[0]];
                Vec<3> tauref = pnts[e[1]] - pnts[e[0]];
                Vec<3,T> tau = mip.GetJacobian()*tauref;
                tau /= mip.GetMeasure();
                LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              Vec<3,T> vshape = val * tau;
                              if (nr==0)
                                shape[i] = vshape;
                              else
                                shape[ii+nr-1] = vshape;
                            }));
              }
            ii += p;
          }
        return;
      }

    for (int i = 0; i < 9; i++)
      ii += order_edge[i] * usegrad_edge[i];

    if (ip.VB() == BND)
      {
        Vec<3,T> X(x, y, z);
        for (int i = 0; i < 5; i++)
          {
            if (i == facetnr)
              {
                IVec<4> f = GetFaceSort (i, vnums);
                if (i < 2)
                  {
                    T xi, eta;
                    HCurlFaceCoords (X, pnts[f[2]], pnts[f[0]], pnts[f[1]], xi, eta);
                    auto ftrafo = HCurlFaceTrafo<T> (mip, pnts[f[0]]-pnts[f[2]], pnts[f[1]]-pnts[f[2]]);
                    HCurlTrigFaceDual (ftrafo, xi, eta, int(order_face[i][0]), usegrad_face[i], type1, ii, shape);
                  }
                else
                  {
                    T s, t;
                    HCurlFaceCoords (X, pnts[f[0]], pnts[f[1]], pnts[f[3]], s, t);
                    auto ftrafo = HCurlFaceTrafo<T> (mip, pnts[f[1]]-pnts[f[0]], pnts[f[3]]-pnts[f[0]]);
                    HCurlQuadDual (ftrafo, s, t, IVec<2>(order_face[i]), usegrad_face[i], type1, ii, shape);
                  }
                return;
              }
            ii += hcurl_ndof::Face (FaceType(i), order_face[i], usegrad_face[i], type1);
          }
        return;
      }

    for (int i = 0; i < 5; i++)
      ii += hcurl_ndof::Face (FaceType(i), order_face[i], usegrad_face[i], type1);

    if (ip.VB() == VOL)
      {
        auto trafo = 1/mip.GetMeasure()*mip.GetJacobian();
        int p = order_cell[0], pz = order_cell[2];
        int pg = p-type1, pgz = pz-type1;
        int n = max(pz, 1) + 1;
        ArrayMem<AutoDiff<1,T>,20> lz(n+1);
        HCurlLegendreAD (n, z, lz);
        // M_{l+1}(z) = int_0^z L_l, d/dz M_{l+1} = L_l(2z-1)
        auto M = [&] (int l) -> T { return (0.5/(2*l+1)) * ((l == 0) ? lz[1].Value()+1 : lz[l+1].Value()-lz[l-1].Value()); };

        // gradients: (0, 0, a(x,y) M_{j+1}(z)), div = a L_j
        if (usegrad_cell && pg >= 2 && pgz >= 1 && !type1)
          {
            // grad of bubbles x y (1-x-y) a(x,y) M_{j+2}(z): blind to div-free fields
            AutoDiff<2,T> adx(x, 0), ady(y, 1);
            AutoDiff<2,T> bub = adx*ady*(1-adx-ady);
            DubinerBasis::Eval (pg-2, adx, ady,
                                SBLambda([&] (size_t nr, auto a)
                                         {
                                           auto ba = bub*a;
                                           for (int j = 0; j < pgz; j++)
                                             {
                                               T mj = M(j+1), dmj = lz[j+1].Value();
                                               shape[ii++] = trafo * Vec<3,T>(ba.DValue(0)*mj, ba.DValue(1)*mj, ba.Value()*dmj);
                                             }
                                         }));
          }
        else if (usegrad_cell && pg >= 2 && pgz >= 1)
          DubinerBasis::Eval (pg-2, x, y,
                              SBLambda([&] (size_t nr, T a)
                                       {
                                         for (int j = 0; j < pgz; j++)
                                           shape[ii++] = trafo * Vec<3,T>(0, 0, a*M(j));
                                       }));

        // curl part: div-free fields
        if (p < 2 || pz < 1) return;
        AutoDiff<2,T> adx(x, 0), ady(y, 1);
        // rot_xy q * L_l(z), q in P_{p-1} without constants
        DubinerBasis::Eval (p-1, adx, ady,
                            SBLambda([&] (size_t nr, auto q)
                                     {
                                       if (nr == 0) return;
                                       for (int l = 0; l < pz-1; l++)
                                         {
                                           T ll = lz[l].Value();
                                           shape[ii++] = trafo * Vec<3,T>(q.DValue(1)*ll, -q.DValue(0)*ll, 0);
                                         }
                                     }));
        if (p < 3) return;
        // (a x L_l, a y L_l, -div_xy(a (x,y)) M_{l+1}), a in P_{p-3}
        DubinerBasis::Eval (p-3, adx, ady,
                            SBLambda([&] (size_t nr, auto a)
                                     {
                                       T av = a.Value();
                                       T dv = 2*av + x*a.DValue(0) + y*a.DValue(1);
                                       for (int l = 0; l < pz-1; l++)
                                         {
                                           T ll = lz[l].Value();
                                           shape[ii++] = trafo * Vec<3,T>(av*x*ll, av*y*ll, -dv*M(l));
                                         }
                                     }));
        // (0, 0, g), g in P_{p-3}
        DubinerBasis::Eval (p-3, x, y,
                            SBLambda([&] (size_t nr, T g)
                                     {
                                       shape[ii++] = trafo * Vec<3,T>(0, 0, g);
                                     }));
      }
  }

  template<> template <typename MIP, typename TFA>
  inline void HCurlHighOrderFE_Shape<ET_PYRAMID> ::
  CalcDualShape2 (const MIP & mip, TFA & shape) const
  {
    typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;
    auto & ip = mip.IP();
    T x = ip(0), y = ip(1), z = ip(2);
    Vec<3,T> X(x, y, z);
    Vec<3> pnts[5] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
    int facetnr = ip.FacetNr();
    int ii = 8;

    if (ip.VB() == BBND)
      {
        for (int i = 0; i < 8; i++)
          {
            int p = order_edge[i] * usegrad_edge[i];
            if (i == facetnr)
              {
                IVec<2> e = GetEdgeSort (i, vnums);
                Vec<3> tauref = pnts[e[1]] - pnts[e[0]];
                T t(0.0);
                for (int k = 0; k < 3; k++)
                  t += (X(k)-pnts[e[0]](k)) * tauref(k);
                T xi = 2*t/InnerProduct(tauref,tauref) - 1;
                Vec<3,T> tau = mip.GetJacobian()*tauref;
                tau /= mip.GetMeasure();
                LegendrePolynomial::Eval
                  (p, xi,
                   SBLambda([&] (size_t nr, T val)
                            {
                              Vec<3,T> vshape = val * tau;
                              if (nr==0)
                                shape[i] = vshape;
                              else
                                shape[ii+nr-1] = vshape;
                            }));
              }
            ii += p;
          }
        return;
      }

    for (int i = 0; i < 8; i++)
      ii += order_edge[i] * usegrad_edge[i];

    if (ip.VB() == BND)
      {
        for (int i = 0; i < 5; i++)
          {
            if (i == facetnr)
              {
                IVec<4> f = GetFaceSort (i, vnums);
                if (i < 4)
                  {
                    T xi, eta;
                    HCurlFaceCoords (X, pnts[f[2]], pnts[f[0]], pnts[f[1]], xi, eta);
                    auto ftrafo = HCurlFaceTrafo<T> (mip, pnts[f[0]]-pnts[f[2]], pnts[f[1]]-pnts[f[2]]);
                    HCurlTrigFaceDual (ftrafo, xi, eta, int(order_face[i][0]), usegrad_face[i], type1, ii, shape);
                  }
                else
                  {
                    T s, t;
                    HCurlFaceCoords (X, pnts[f[0]], pnts[f[1]], pnts[f[3]], s, t);
                    auto ftrafo = HCurlFaceTrafo<T> (mip, pnts[f[1]]-pnts[f[0]], pnts[f[3]]-pnts[f[0]]);
                    HCurlQuadDual (ftrafo, s, t, IVec<2>(order_face[i]), usegrad_face[i], type1, ii, shape);
                  }
                return;
              }
            ii += hcurl_ndof::Face (FaceType(i), order_face[i], usegrad_face[i], type1);
          }
        return;
      }

    if (ip.VB() == VOL)
      {
        // cell: the reference cell shapes themselves, Gram matrix is spd
        int first = this->ndof - hcurl_ndof::Cell (ET_PYRAMID, order_cell, usegrad_cell, type1);
        auto trafo = 1/mip.GetMeasure()*mip.GetJacobian();
        auto cellshapes = SBLambda ([&](size_t i, auto s)
                                    {
                                      if (int(i) < first) return;
                                      Vec<3,T> v = trafo * s.Value();
                                      shape[i] = v;
                                    });
        this->T_CalcShape (GetTIPGrad<3>(ip), cellshapes);
      }
  }
  
  template <ELEMENT_TYPE ET, 
            template <ELEMENT_TYPE ET2> class TSHAPES, 
            typename BASE>
  void HCurlHighOrderFE<ET,TSHAPES,BASE> ::
  CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const
  {
    shape.AddSize(ndof, bmip.DimSpace()) = 0.0;
    Switch<4-DIM>
      (bmip.DimSpace()-DIM,[this,&bmip,shape](auto CODIM)
       {
         auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);
         static_cast<const HCurlHighOrderFE_Shape<ET>*> (this)
           -> CalcDualShape2 (mip, SBLambda([shape] (size_t i, auto val)
                                            {
                                              shape.Row(i) = val;
                                            }));
       });
  }
  

  template <ELEMENT_TYPE ET, 
            template <ELEMENT_TYPE ET2> class TSHAPES, 
            typename BASE>
  void HCurlHighOrderFE<ET,TSHAPES,BASE> ::
  CalcDualShape (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> shapes) const
  {
 
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,shapes](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         shapes.AddSize(ndof*DIMSPACE, mir.Size()) = 0.0;
         for (size_t i = 0; i < mir.Size(); i++)
           static_cast<const HCurlHighOrderFE_Shape<ET>*> (this)
             -> CalcDualShape2 (mir[i], SBLambda([shapes,i,DIMSPACE=DIMSPACE] (size_t j, auto val)
                                                 {
                                                   for (size_t k = 0; k < DIMSPACE; k++)
                                                     shapes(j*DIMSPACE+k, i) = val(k);
                                                 }));
       });
  }

  template <ELEMENT_TYPE ET, 
            template <ELEMENT_TYPE ET2> class TSHAPES, 
            typename BASE>
  void HCurlHighOrderFE<ET,TSHAPES,BASE> ::
  EvaluateDual (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIMSPACE,SIMD<double>> sum (SIMD<double>(0.0));
             static_cast<const HCurlHighOrderFE_Shape<ET>*> (this)
               -> CalcDualShape2 (mir[i], SBLambda([&sum, coefs] (size_t j, auto val)
                                                   {
                                                     sum += coefs(j) * val;
                                                   }));
             for (size_t k = 0; k < DIMSPACE; k++)
               values(k, i) = sum(k);
           }
       });
  }

  template <ELEMENT_TYPE ET, 
            template <ELEMENT_TYPE ET2> class TSHAPES, 
            typename BASE>
  void HCurlHighOrderFE<ET,TSHAPES,BASE> ::
  AddDualTrans (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> values,
                BareSliceVector<double> coefs) const
  {
    Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Vec<DIMSPACE,SIMD<double>> value = values.Col(i);
             // for (size_t k = 0; k < DIMSPACE; k++)
             // value(k) = values(k, i);
             
             static_cast<const HCurlHighOrderFE_Shape<ET>*> (this)
               -> CalcDualShape2 (mir[i], SBLambda([value, coefs] (size_t j, auto val)
                                                   {
                                                     coefs(j) += HSum(InnerProduct(value,val));
                                                   }));
           }
       });
  }
}

#endif
