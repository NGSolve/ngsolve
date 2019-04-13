#ifndef FILE_HCURLDIVFE_IMPL
#define FILE_HCURLDIVFE_IMPL

/**************************************************************************/
/* File:   hcurldivfe_impl.hpp                                            */
/* Author: Philip Lederer                                                 */
/* Date:   2017/2018                                                      */
/**************************************************************************/


namespace ngfem
{
  /* ############### edge basis functions - div-free ############### */
  /* sigma(grad v) = Curl(grad v), where Curl is the 1D to 2D curl operator */
  template <int D, typename T>  class T_Sigma_gradv;
  template <typename T>  class T_Sigma_gradv<2,T>
  {
    AutoDiffDiff<2,T> v;
  public:
    T_Sigma_gradv  (AutoDiffDiff<2,T> av) : v(av){ ; }
    
    Vec<4,T> Shape() {
      return Vec<4,T> (-v.DDValue(0,1), v.DDValue(0,0),
		     -v.DDValue(1,1),v.DDValue(0,1)
		     );
    }

    Vec<2,T> DivShape()
    {      
      return Vec<2,T> (0.0,0.0);     
    }

    Vec<2,T> CurlShape()
    {      
      return Vec<2,T> (0.0,0.0);     
    }
  };

  template <int D, typename T>
  auto Sigma_gradv (AutoDiffDiff<D,T> av) { return T_Sigma_gradv<D,T>(av); }

  /* ############### div-free basis function WITH trace ############### */
  /* For basis functions including the trace */
  
  template <int D, typename T> class T_Sigma_gradu_v;
  template <typename T> class T_Sigma_gradu_v<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_Sigma_gradu_v  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av){ ; }

    Vec<4,T> Shape() {      
      return Vec<4,T> (-u.DDValue(1,0) * v.Value()  -  v.DValue(1)*u.DValue(0),
		     u.DDValue(0,0) * v.Value() + v.DValue(0)*u.DValue(0),
		     -u.DDValue(1,1) * v.Value() - v.DValue(1)*u.DValue(1),
		     u.DDValue(0,1) * v.Value() +  v.DValue(0)*u.DValue(1)
		     );
    }

    Vec<2,T> DivShape()
    {     
      return Vec<2,T> (0.0,0.0);
    }

    Vec<2,T> CurlShape()
    {
      T vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T uxx = u.DDValue(0,0), uxy = u.DDValue(0,1), uyy = u.DDValue(1,1);
      T vx = v.DValue(0), vy = v.DValue(1);

      return Vec<2,T> ( (vy*uxy - vx*uyy) +  (ux * vyy - uy * vxy), (-vy*uxx + vx*uxy) + (-ux*vxy + uy*vxx));     
    }    
  };

  template <int D, typename T>
  auto Sigma_gradu_v (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av ) { return T_Sigma_gradu_v<D,T>(au,av); }

  /* ############### Type 1.1 (for QUAD) - inner basis functions - not div-free ############### */
  /* dev ( grad(u) * Curl(v) ) */
  template <int D, typename T> class T_Gradu_Curlv;
  template <typename T>  class T_Gradu_Curlv<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_Gradu_Curlv  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av){ ; }

    Vec<4,T> Shape() {

      auto trace = (-  v.DValue(1)*u.DValue(0) + v.DValue(0)*u.DValue(1)) / 2.0;
      
      return Vec<4,T> (-  v.DValue(1)*u.DValue(0) - trace,
		      v.DValue(0)*u.DValue(0),
		     - v.DValue(1)*u.DValue(1),
		     v.DValue(0)*u.DValue(1) - trace
		     );
    }

    Vec<2,T> DivShape()
    {
      T vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);
      T uxx = u.DDValue(0,0), uxy = u.DDValue(0,1), uyy = u.DDValue(1,1);      
      T ux = u.DValue(0), uy = u.DValue(1);
      T vx = v.DValue(0), vy = v.DValue(1);
      
      return Vec<2,T> (-vy*uxx + vx*uxy - 0.5*(-vxy*ux - vy*uxx + vxx*uy + vx*uxy),-vy*uxy + vx * uyy - 0.5*(-vyy*ux - vy*uxy + vxy*uy + vx*uyy));
    }

    Vec<2,T> CurlShape()
    {
      T vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T uxx = u.DDValue(0,0), uxy = u.DDValue(0,1), uyy = u.DDValue(1,1);
      T vx = v.DValue(0), vy = v.DValue(1);

      return Vec<2,T> (vyy*ux - vxy*uy,-vxy*ux+vxx*uy);     
    }    
  };

  template <int D, typename T>
  auto Gradu_Curlv (AutoDiffDiff<D,T> au,AutoDiffDiff<D,T> av) { return T_Gradu_Curlv<D, T>(au,av); }

  
  /* ############### Type 2 (QUAD) - inner basis functions - div-free ############### */
  /* u * sigma(grad v) = Curl(grad v), where Curl is the 1D to 2D curl operator */
  template <int D, typename T> class T_u_Sigma_gradv;
  template <typename T> class T_u_Sigma_gradv<2,T>
  {
    AutoDiffDiff<2,T> v,u;
  public:
    T_u_Sigma_gradv  (AutoDiffDiff<2,T> au,AutoDiffDiff<2,T> av) : u(au),v(av){ ; }
    
    Vec<4,T> Shape() {
      
      return Vec<4,T> (-u.Value() * v.DDValue(0,1), u.Value() *  v.DDValue(0,0),
		     -u.Value() * v.DDValue(1,1), u.Value() * v.DDValue(0,1)
		     );
    }

    Vec<2,T> DivShape()
    {      
      return Vec<2,T> (-u.DValue(0) *v.DDValue(0,1) + u.DValue(1) * v.DDValue(0,0) , -u.DValue(0) *v.DDValue(1,1) + u.DValue(1) * v.DDValue(0,1));     
    }

    Vec<2,T> CurlShape()
    {      
      return Vec<2,T> (u.DValue(1) *v.DDValue(0,1) - u.DValue(0) * v.DDValue(1,1) ,-u.DValue(1) *v.DDValue(0,0) + u.DValue(0) * v.DDValue(0,1) );     
    }
  };

  template <int D, typename T>
  auto u_Sigma_gradv (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_u_Sigma_gradv<D,T>(au,av); }


  /* ############### Type 2 - inner basis functions - NOT div-free ############### */
  /* Curl(grad(u)) * v - grad(u) * Curl(v) */
  template <int D, typename T>  class T_Curlgraduv_graducurlv;
  template <typename T>  class T_Curlgraduv_graducurlv<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_Curlgraduv_graducurlv  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av){ ; }

    Vec<4,T> Shape() {

      auto trace = (v.DValue(1)*u.DValue(0) -  v.DValue(0)*u.DValue(1)  )/2.0;
      
      return Vec<4,T> (-u.DDValue(1,0) * v.Value()  +  v.DValue(1)*u.DValue(0) - trace,
		     u.DDValue(0,0) * v.Value() - v.DValue(0)*u.DValue(0),
		     -u.DDValue(1,1) * v.Value() + v.DValue(1)*u.DValue(1),
		     u.DDValue(0,1) * v.Value() -  v.DValue(0)*u.DValue(1) - trace
		     );
    }

    Vec<2,T> DivShape()
    {
      T vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T uxx = u.DDValue(0,0), uxy = u.DDValue(0,1), uyy = u.DDValue(1,1);
      T vx = v.DValue(0), vy = v.DValue(1);
            
      return 0.5 * Vec<2,T> (3 * uxx * vy - 3 * uxy * vx - vxy*ux + vxx*uy, 3* uxy * vy - 3 * uyy * vx - vyy*ux + vxy*uy);

    }

    Vec<2,T> CurlShape()
    {
      T vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T uxx = u.DDValue(0,0), uxy = u.DDValue(0,1), uyy = u.DDValue(1,1);
      T vx = v.DValue(0), vy = v.DValue(1);
      
      return 0.5 *  Vec<2,T> ( (3*vy*uxy - 3*vx*uyy) - (ux * vyy - uy * vxy), (-3*vy*uxx + 3*vx*uxy) - (-ux*vxy + uy*vxx));     
    }
    
  };

  template <int D, typename T>
  auto Curlgraduv_graducurlv (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_Curlgraduv_graducurlv<D,T>(au, av); }

  /* ############### Type 3 - inner basis functions - div-free ############### */
  /*  Curl( [grad(l1) l2 - l1 grad(l2)] * v ) */
  template <int D, typename T> class T_type4;
  template <typename T> class T_type4<2,T>
  {
    AutoDiffDiff<2,T> l1,l2,v;
  public:
    T_type4  (AutoDiffDiff<2,T> lam1, AutoDiffDiff<2,T> lam2, AutoDiffDiff<2,T> av) : l1(lam1), l2(lam2), v(av){ ; }

    Vec<4,T> Shape() {
      T lam1x = l1.DValue(0), lam1y = l1.DValue(1);
      T lam2x = l2.DValue(0), lam2y = l2.DValue(1);
      T vx = v.DValue(0), vy = v.DValue(1);                      

      //auto trace = ( (v.Value() * ( - lam1x * lam2y + lam2x * lam1y) - (lam1x*l2.Value() - lam2x*l1.Value()) * vy)
      //		     + ( v.Value() * ( lam1y * lam2x - lam2y * lam1x) + (lam1y*l2.Value() - lam2y*l1.Value()) * vx))/2.0;
      
      return Vec<4,T> (v.Value() * ( - lam1x * lam2y + lam2x * lam1y) - (lam1x*l2.Value() - lam2x*l1.Value()) * vy,
		      (lam1x*l2.Value() - lam2x*l1.Value()) * vx,
		     -(lam1y*l2.Value() - lam2y*l1.Value()) * vy,
		     v.Value() * ( lam1y * lam2x - lam2y * lam1x) + (lam1y*l2.Value() - lam2y*l1.Value()) * vx
		     ); 
    }

    Vec<2,T> DivShape()
    {     
      return Vec<2,T> (0.0,0.0);     
    }

    Vec<2,T> CurlShape()
    {
      T lam1x = l1.DValue(0), lam1y = l1.DValue(1); 
      T lam2x = l2.DValue(0), lam2y = l2.DValue(1); 
      T vx = v.DValue(0), vy = v.DValue(1), vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);

      return Vec<2,T> ( vyy*(lam1x*l2.Value() - lam2x*l1.Value()) - vxy * (lam1y*l2.Value() - lam2y*l1.Value()) - 3* vy*(-lam1x*lam2y+lam2x*lam1y),
      		      vxx*(lam1y*l2.Value() - lam2y*l1.Value()) - vxy * (lam1x*l2.Value() - lam2x*l1.Value()) + 3* vx*(-lam1x*lam2y+lam2x*lam1y)
      				);
    }    
  }; 

  template <int D, typename T>
  auto type4 (AutoDiffDiff<D,T> al1, AutoDiffDiff<D,T> al2, AutoDiffDiff<D,T> av) { return T_type4<D,T>(al1, al2, av); }
  
  /* GG-bubble 2D from jay */
  
  template <int D, typename T> class T_GGbubble;
  template <typename T> class T_GGbubble<2,T>
  {
    AutoDiffDiff<2,T> S;
    AutoDiffDiff<2,T> b;
  public:
    T_GGbubble  (AutoDiffDiff<2,T> aS, AutoDiffDiff<2,T> ab ) : S(aS), b(ab){ ; }

    Vec<4,T> Shape() {
      auto B = b.Value();
      auto Bx = b.DValue(0);
      auto By = b.DValue(1);
      auto Bxx = b.DDValue(0,0);
      auto Bxy = b.DDValue(0,1);
      auto Byy = b.DDValue(1,1);

      auto trace = 0.5 * (By*S.DValue(0)  - Bx*S.DValue(1));

      return Vec<4,T> (S.DDValue(1,0)*B + By*S.DValue(0) - trace,
		       -S.DDValue(0,0)*B - Bx*S.DValue(0),
		       S.DDValue(1,1)*B + By*S.DValue(1),
		       -S.DDValue(0,1)*B - Bx*S.DValue(1) - trace);
    }

    Vec<2,T> DivShape()
    {
      auto Bx = b.DValue(0);
      auto By = b.DValue(1);
      auto Bxx = b.DDValue(0,0);
      auto Bxy = b.DDValue(0,1);
      auto Byy = b.DDValue(1,1);
      
      return -0.5 * Vec<2,T>(Bxy*S.DValue(0) + By * S.DDValue(0,0) - Bxx*S.DValue(1) - Bx*S.DDValue(0,1),
			   Byy*S.DValue(0) + By * S.DDValue(1,0) - Bxy*S.DValue(1) - Bx*S.DDValue(1,1));
    }

    Vec<2,T> CurlShape()
    {     
      throw Exception("not implemented for curlbubbles");
    }
    
  };

  template <int D, typename T>
  auto GGbubble (AutoDiffDiff<D,T> aS, AutoDiffDiff<D,T> ab ) { return T_GGbubble<D, T>(aS, ab); }
  
  // GG-bubble of Jay
  // computes curl ( curl A B) where A is a skew sym matrix which is L^2 orthogonal on P^k-1
  // and B is the matrix bubble
  
  template <int D, typename T> class T_GGbubble_B1;
  template <typename T> class T_GGbubble_B1<3,T>
  {
    AutoDiffDiff<3,T> S;
    AutoDiffDiff<3,T> b0;
    AutoDiffDiff<3,T> b1;
    AutoDiffDiff<3,T> b2;
    AutoDiffDiff<3,T> b3;
    
  public:
    T_GGbubble_B1  (AutoDiffDiff<3,T> aS, AutoDiffDiff<3,T> ab0, AutoDiffDiff<3,T> ab1, AutoDiffDiff<3,T> ab2, AutoDiffDiff<3,T> ab3) : S(aS), b0(ab0), b1(ab1), b2(ab2), b3(ab3){ ; }

    Vec<9,T> Shape() {
      /////////////////////
      auto b1_x  =  b1.DValue(0);   auto b2_x = b2.DValue(0);
      auto b1_y  =  b1.DValue(1);   auto b2_y = b2.DValue(1); 
      auto b1_z  =  b1.DValue(2);   auto b2_z = b2.DValue(2);
      
      auto b3_x = b3.DValue(0);     auto b0_x = b0.DValue(0);
      auto b3_y = b3.DValue(1);     auto b0_y = b0.DValue(1);
      auto b3_z = b3.DValue(2);     auto b0_z = b0.DValue(2);
      
      auto B0 = b0.Value();  auto B1 = b1.Value();
      auto B2 = b2.Value();  auto B3 = b3.Value();
      /////////////////////
      auto S_x = S.DValue(0); auto S0 = S.Value();
      auto S_y = S.DValue(1);
      auto S_z = S.DValue(2);
    
      auto S_xx = S.DDValue(0,0); auto S_yy = S.DDValue(1,1);
      auto S_xy = S.DDValue(0,1); auto S_yz = S.DDValue(1,2); 
      auto S_xz = S.DDValue(0,2); auto S_zz = S.DDValue(2,2);
      /////////////////////
     
      Vec<9,T> sigmaref;

      auto trace = 1.0/3.0 * (S_x * (b0_z - b0_y) + S_y * (b0_x - b1_z) + S_z * (b1_y - b0_x));

      sigmaref(0) = 0.0;
      sigmaref(1) = 0.0;
      sigmaref(2) = 0.0;
      
      sigmaref(3) = -S_yy*B0 -S_y*b0_y+S_xy*B0 +S_x*b0_y; // dy(a3)
      sigmaref(3)+=  S_yz*B0 +S_y*b0_z-S_xz*B2 -S_x*b2_z; //-dz(a2)
      sigmaref(4) = -S_yz*B1 -S_y*b1_z+S_xz*B0 +S_x*b0_z; // dz(a1)
      sigmaref(4)+=  S_xy*B0 +S_y*b0_x-S_xx*B0 -S_x*b0_x; //-dx(a3)
      sigmaref(5) = -S_xy*B0 -S_y*b0_x+S_xx*B2 +S_x*b2_x; // dx(a2)
      sigmaref(5)+=  S_yy*B1 +S_y*b1_y-S_xy*B0 -S_x*b0_y; //-dy(a1)
      
      sigmaref(6) = -S_yz*B0 -S_z*b0_y+S_xy*B3 +S_x*b3_y; // dy(b3)
      sigmaref(6)+=  S_zz*B0 +S_z*b0_z-S_xz*B0 -S_x*b0_z; //-dz(b2)
      sigmaref(7) = -S_zz*B1 -S_z*b1_z+S_xz*B0 +S_x*b0_z; // dz(b1)
      sigmaref(7)+=  S_xz*B0 +S_z*b0_x-S_xx*B3 -S_x*b3_x; //-dx(b3)
      sigmaref(8) = -S_xz*B0 -S_z*b0_x+S_xx*B0 +S_x*b0_x; // dx(b2)
      sigmaref(8)+=  S_yz*B1 +S_z*b1_y-S_xy*B0 -S_x*b0_y; //-dy(b1)

      sigmaref(0) -= trace;
      sigmaref(4) -= trace;
      sigmaref(8) -= trace;
      
      return sigmaref;
    }

    Vec<3,T> DivShape()
    {
      auto b1_x  =  b1.DValue(0);
      auto b1_y  =  b1.DValue(1); 
      auto b1_z  =  b1.DValue(2);
      
      auto b0_x = b0.DValue(0);
      auto b0_y = b0.DValue(1);
      auto b0_z = b0.DValue(2);

      auto b1_xx  =  b1.DDValue(0,0);
      auto b1_yy  =  b1.DDValue(1,1);         
      auto b1_zz  =  b1.DDValue(2,2);   
      auto b1_yx  =  b1.DDValue(0,1);   
      auto b1_zx  =  b1.DDValue(0,2);
      auto b1_yz  =  b1.DDValue(1,2);

      auto b0_xx  =  b0.DDValue(0,0);
      auto b0_yy  =  b0.DDValue(1,1);         
      auto b0_zz  =  b0.DDValue(2,2);   
      auto b0_yx  =  b0.DDValue(0,1);   
      auto b0_zx  =  b0.DDValue(0,2);
      auto b0_yz  =  b0.DDValue(1,2);
      /////////////////////
      auto S_x = S.DValue(0); 
      auto S_y = S.DValue(1);
      auto S_z = S.DValue(2);
    
      auto S_xx = S.DDValue(0,0); auto S_yy = S.DDValue(1,1);
      auto S_xy = S.DDValue(0,1); auto S_yz = S.DDValue(1,2); 
      auto S_xz = S.DDValue(0,2); auto S_zz = S.DDValue(2,2);
      /////////////////////

      //(S_x * (b0_z - b0_y) + S_y * (b0_x - b1_z) + S_z * (b1_y - b0_x));

      auto div1 = -1.0/3.0 * (S_xx*(b0_z-b0_y)+S_x*(b0_zx-b0_yx) + S_xy*(b0_x-b1_z)+S_y*(b0_xx-b1_zx) + S_xz*(b1_y-b0_x)+S_z*(b1_yx-b0_xx));
      auto div2 = -1.0/3.0 * (S_xy*(b0_z-b0_y)+S_x*(b0_yz-b0_yy) + S_yy*(b0_x-b1_z)+S_y*(b0_yx-b1_yz) + S_yz*(b1_y-b0_x)+S_z*(b1_yy-b0_yx));
      auto div3 = -1.0/3.0 * (S_xz*(b0_z-b0_y)+S_x*(b0_zz-b0_yz) + S_yz*(b0_x-b1_z)+S_y*(b0_zx-b1_zz) + S_zz*(b1_y-b0_x)+S_z*(b1_yz-b0_zx));
      
      //return -1.0/3 * Vec<3,T>(0.0,0.0,0.0);
      return  Vec<3,T> (div1,div2,div3);     
    }

    Vec<3,T> CurlShape()
    {     
      throw Exception("not implemented for GG-bubbles");
    }
    
  };

  template <int D, typename T>
  auto GGbubble_B1 (AutoDiffDiff<D,T> aS, AutoDiffDiff<D,T> ab0, AutoDiffDiff<D,T> ab1, AutoDiffDiff<D,T> ab2,AutoDiffDiff<D,T> ab3) { return T_GGbubble_B1<D, T>(aS, ab0, ab1, ab2, ab3); }

  ////////////////////////////////
  
    template <int D, typename T> class T_GGbubble_B2;
  template <typename T> class T_GGbubble_B2<3,T>
  {
    AutoDiffDiff<3,T> S;
    AutoDiffDiff<3,T> b0;
    AutoDiffDiff<3,T> b1;
    AutoDiffDiff<3,T> b2;
    AutoDiffDiff<3,T> b3;
    
  public:
    T_GGbubble_B2  (AutoDiffDiff<3,T> aS, AutoDiffDiff<3,T> ab0, AutoDiffDiff<3,T> ab1, AutoDiffDiff<3,T> ab2, AutoDiffDiff<3,T> ab3) : S(aS), b0(ab0), b1(ab1), b2(ab2), b3(ab3){ ; }

    Vec<9,T> Shape() {
      /////////////////////
      auto b1_x  =  b1.DValue(0);   auto b2_x = b2.DValue(0);
      auto b1_y  =  b1.DValue(1);   auto b2_y = b2.DValue(1); 
      auto b1_z  =  b1.DValue(2);   auto b2_z = b2.DValue(2);
      
      auto b3_x = b3.DValue(0);     auto b0_x = b0.DValue(0);
      auto b3_y = b3.DValue(1);     auto b0_y = b0.DValue(1);
      auto b3_z = b3.DValue(2);     auto b0_z = b0.DValue(2);
      
      auto B0 = b0.Value();  auto B1 = b1.Value();
      auto B2 = b2.Value();  auto B3 = b3.Value();
      /////////////////////
      auto S_x = S.DValue(0); auto S0 = S.Value();
      auto S_y = S.DValue(1);
      auto S_z = S.DValue(2);
    
      auto S_xx = S.DDValue(0,0); auto S_yy = S.DDValue(1,1);
      auto S_xy = S.DDValue(0,1); auto S_yz = S.DDValue(1,2); 
      auto S_xz = S.DDValue(0,2); auto S_zz = S.DDValue(2,2);
      /////////////////////
      
      Vec<9,T> sigmaref;
      auto trace = 1.0/3.0 * ((b2_z-b0_y)*S_x + (b0_x-b0_z)*S_y + (b0_y-b2_x)*S_z);

      sigmaref(0) = -(-S_yy*B0 -S_y*b0_y+S_xy*B0 +S_x*b0_y); // dy(a3)
      sigmaref(0)+= -( S_yz*B0 +S_y*b0_z-S_xz*B2 -S_x*b2_z); //-dz(a2)
      sigmaref(1) = -(-S_yz*B1 -S_y*b1_z+S_xz*B0 +S_x*b0_z); // dz(a1)
      sigmaref(1)+= -( S_xy*B0 +S_y*b0_x-S_xx*B0 -S_x*b0_x); //-dx(a3)
      sigmaref(2) = -(-S_xy*B0 -S_y*b0_x+S_xx*B2 +S_x*b2_x); // dx(a2)
      sigmaref(2)+= -( S_yy*B1 +S_y*b1_y-S_xy*B0 -S_x*b0_y); //-dy(a1)
      
      sigmaref(3) = 0.0;
      sigmaref(4) = 0.0;
      sigmaref(5) = 0.0;
      
      sigmaref(6) = -S_yz*B0 -S_z*b0_y+S_yy*B3 +S_y*b3_y; // dy(c3)
      sigmaref(6)+=  S_zz*B2 +S_z*b2_z-S_yz*B0 -S_y*b0_z; //-dz(c2)
      sigmaref(7) = -S_zz*B0 -S_z*b0_z+S_yz*B0 +S_y*b0_z; // dz(c1)
      sigmaref(7)+=  S_xz*B0 +S_z*b0_x-S_xy*B3 -S_y*b3_x; //-dx(c3)
      sigmaref(8) = -S_xz*B2 -S_z*b2_x+S_xy*B0 +S_y*b0_x; // dx(c2)
      sigmaref(8)+=  S_yz*B0 +S_z*b0_y-S_yy*B0 -S_y*b0_y; //-dy(c1)

      sigmaref(0) -= trace;
      sigmaref(4) -= trace;
      sigmaref(8) -= trace;
      
      return sigmaref;
    }

    Vec<3,T> DivShape()
    {
      /////////////////////
      auto b1_x  =  b1.DValue(0);   auto b2_x = b2.DValue(0);
      auto b1_y  =  b1.DValue(1);   auto b2_y = b2.DValue(1); 
      auto b1_z  =  b1.DValue(2);   auto b2_z = b2.DValue(2);
      
      auto b3_x = b3.DValue(0);     auto b0_x = b0.DValue(0);
      auto b3_y = b3.DValue(1);     auto b0_y = b0.DValue(1);
      auto b3_z = b3.DValue(2);     auto b0_z = b0.DValue(2);

      ///
      auto b1_xx  =  b1.DDValue(0,0);      auto b2_xx  =  b2.DDValue(0,0);
      auto b1_yy  =  b1.DDValue(1,1);      auto b2_yy  =  b2.DDValue(1,1);         
      auto b1_zz  =  b1.DDValue(2,2);      auto b2_zz  =  b2.DDValue(2,2);   
      auto b1_yx  =  b1.DDValue(0,1);      auto b2_yx  =  b2.DDValue(0,1);   
      auto b1_zx  =  b1.DDValue(0,2);      auto b2_zx  =  b2.DDValue(0,2);
      auto b1_yz  =  b1.DDValue(1,2);      auto b2_yz  =  b2.DDValue(1,2);

      auto b3_xx  =  b3.DDValue(0,0);      auto b0_xx  =  b0.DDValue(0,0);
      auto b3_yy  =  b3.DDValue(1,1);      auto b0_yy  =  b0.DDValue(1,1);         
      auto b3_zz  =  b3.DDValue(2,2);      auto b0_zz  =  b0.DDValue(2,2);   
      auto b3_yx  =  b3.DDValue(0,1);      auto b0_yx  =  b0.DDValue(0,1);   
      auto b3_zx  =  b3.DDValue(0,2);      auto b0_zx  =  b0.DDValue(0,2);
      auto b3_yz  =  b3.DDValue(1,2);      auto b0_yz  =  b0.DDValue(1,2);
      
      /////////////////////
      auto S_x = S.DValue(0); 
      auto S_y = S.DValue(1);
      auto S_z = S.DValue(2);
    
      auto S_xx = S.DDValue(0,0); auto S_yy = S.DDValue(1,1);
      auto S_xy = S.DDValue(0,1); auto S_yz = S.DDValue(1,2); 
      auto S_xz = S.DDValue(0,2); auto S_zz = S.DDValue(2,2);
      /////////////////////

      //return -1.0/3 * Vec<3,T>(0.0,0.0,0.0);
      return -1.0/3.0 * Vec<3,T> (S_xx*(b2_z-b0_y)+S_x*(b2_zx-b0_yx) + S_xy*(b0_x-b0_z)+S_y*(b0_xx-b0_zx) + S_xz*(b0_y-b2_x)+S_z*(b0_yx-b2_xx),
				  S_xy*(b2_z-b0_y)+S_x*(b2_yz-b0_yy) + S_yy*(b0_x-b0_z)+S_y*(b0_yx-b0_yz) + S_yz*(b0_y-b2_x)+S_z*(b0_yy-b2_yx),
				  S_xz*(b2_z-b0_y)+S_x*(b2_zz-b0_yz) + S_yz*(b0_x-b0_z)+S_y*(b0_zx-b0_zz) + S_zz*(b0_y-b2_x)+S_z*(b0_yz-b2_zx));
    }

    Vec<3,T> CurlShape()
    {     
      throw Exception("not implemented for GG-bubbles");
    }
    
  };

  template <int D, typename T>
  auto GGbubble_B2 (AutoDiffDiff<D,T> aS, AutoDiffDiff<D,T> ab0, AutoDiffDiff<D,T> ab1, AutoDiffDiff<D,T> ab2, AutoDiffDiff<D,T> ab3) { return T_GGbubble_B2<D, T>(aS, ab0, ab1, ab2, ab3); }

  //////////////////
  
    template <int D, typename T> class T_GGbubble_B3;
  template <typename T> class T_GGbubble_B3<3,T>
  {
    AutoDiffDiff<3,T> S;
    AutoDiffDiff<3,T> b0;
    AutoDiffDiff<3,T> b1;
    AutoDiffDiff<3,T> b2;
    AutoDiffDiff<3,T> b3;
    
  public:
    T_GGbubble_B3  (AutoDiffDiff<3,T> aS, AutoDiffDiff<3,T> ab0, AutoDiffDiff<3,T> ab1, AutoDiffDiff<3,T> ab2, AutoDiffDiff<3,T> ab3) : S(aS), b0(ab0), b1(ab1), b2(ab2), b3(ab3){ ; }

    Vec<9,T> Shape() {
      /////////////////////
      auto b1_x  =  b1.DValue(0);   auto b2_x = b2.DValue(0);
      auto b1_y  =  b1.DValue(1);   auto b2_y = b2.DValue(1); 
      auto b1_z  =  b1.DValue(2);   auto b2_z = b2.DValue(2);
      
      auto b3_x = b3.DValue(0);     auto b0_x = b0.DValue(0);
      auto b3_y = b3.DValue(1);     auto b0_y = b0.DValue(1);
      auto b3_z = b3.DValue(2);     auto b0_z = b0.DValue(2);
      
      auto B0 = b0.Value();  auto B1 = b1.Value();
      auto B2 = b2.Value();  auto B3 = b3.Value();
      /////////////////////
      auto S_x = S.DValue(0); auto S0 = S.Value();
      auto S_y = S.DValue(1);
      auto S_z = S.DValue(2);
    
      auto S_xx = S.DDValue(0,0); auto S_yy = S.DDValue(1,1);
      auto S_xy = S.DDValue(0,1); auto S_yz = S.DDValue(1,2); 
      auto S_xz = S.DDValue(0,2); auto S_zz = S.DDValue(2,2);
      /////////////////////
      
      Vec<9,T> sigmaref;
      auto trace = 1.0/3.0 *((b0_z-b3_y)*S_x + (b3_x-b0_z)*S_y + (b0_y-b0_x)*S_z);    

      sigmaref(0) = -(-S_yz*B0 -S_z*b0_y+S_xy*B3 +S_x*b3_y); // dy(b3)
      sigmaref(0)+= -( S_zz*B0 +S_z*b0_z-S_xz*B0 -S_x*b0_z); //-dz(b2)
      sigmaref(1) = -(-S_zz*B1 -S_z*b1_z+S_xz*B0 +S_x*b0_z); // dz(b1)
      sigmaref(1)+= -( S_xz*B0 +S_z*b0_x-S_xx*B3 -S_x*b3_x); //-dx(b3)
      sigmaref(2) = -(-S_xz*B0 -S_z*b0_x+S_xx*B0 +S_x*b0_x); // dx(b2)
      sigmaref(2)+= -( S_yz*B1 +S_z*b1_y-S_xy*B0 -S_x*b0_y); //-dy(b1)

      sigmaref(3) = -(-S_yz*B0 -S_z*b0_y+S_yy*B3 +S_y*b3_y); // dy(c3)
      sigmaref(3)+= -( S_zz*B2 +S_z*b2_z-S_yz*B0 -S_y*b0_z); //-dz(c2)
      sigmaref(4) = -(-S_zz*B0 -S_z*b0_z+S_yz*B0 +S_y*b0_z); // dz(c1)
      sigmaref(4)+= -( S_xz*B0 +S_z*b0_x-S_xy*B3 -S_y*b3_x); //-dx(c3)
      sigmaref(5) = -(-S_xz*B2 -S_z*b2_x+S_xy*B0 +S_y*b0_x); // dx(c2)
      sigmaref(5)+= -( S_yz*B0 +S_z*b0_y-S_yy*B0 -S_y*b0_y); //-dy(c1)          
      
      sigmaref(6) = 0.0;
      sigmaref(7) = 0.0;
      sigmaref(8) = 0.0;

      //auto trace = 1.0/3.0 * (sigmaref(0) + sigmaref(4) + sigmaref(8));
      sigmaref(0) -= trace;
      sigmaref(4) -= trace;
      sigmaref(8) -= trace;
      return sigmaref;
    }

    Vec<3,T> DivShape()
    {
      /////////////////////
      auto b1_x  =  b1.DValue(0);   auto b2_x = b2.DValue(0);
      auto b1_y  =  b1.DValue(1);   auto b2_y = b2.DValue(1); 
      auto b1_z  =  b1.DValue(2);   auto b2_z = b2.DValue(2);
      
      auto b3_x = b3.DValue(0);     auto b0_x = b0.DValue(0);
      auto b3_y = b3.DValue(1);     auto b0_y = b0.DValue(1);
      auto b3_z = b3.DValue(2);     auto b0_z = b0.DValue(2);

      ///
      auto b1_xx  =  b1.DDValue(0,0);      auto b2_xx  =  b2.DDValue(0,0);
      auto b1_yy  =  b1.DDValue(1,1);      auto b2_yy  =  b2.DDValue(1,1);         
      auto b1_zz  =  b1.DDValue(2,2);      auto b2_zz  =  b2.DDValue(2,2);   
      auto b1_yx  =  b1.DDValue(0,1);      auto b2_yx  =  b2.DDValue(0,1);   
      auto b1_zx  =  b1.DDValue(0,2);      auto b2_zx  =  b2.DDValue(0,2);
      auto b1_yz  =  b1.DDValue(1,2);      auto b2_yz  =  b2.DDValue(1,2);

      auto b3_xx  =  b3.DDValue(0,0);      auto b0_xx  =  b0.DDValue(0,0);
      auto b3_yy  =  b3.DDValue(1,1);      auto b0_yy  =  b0.DDValue(1,1);         
      auto b3_zz  =  b3.DDValue(2,2);      auto b0_zz  =  b0.DDValue(2,2);   
      auto b3_yx  =  b3.DDValue(0,1);      auto b0_yx  =  b0.DDValue(0,1);   
      auto b3_zx  =  b3.DDValue(0,2);      auto b0_zx  =  b0.DDValue(0,2);
      auto b3_yz  =  b3.DDValue(1,2);      auto b0_yz  =  b0.DDValue(1,2);
      /////////////////////
      auto S_x = S.DValue(0); 
      auto S_y = S.DValue(1);
      auto S_z = S.DValue(2);
    
      auto S_xx = S.DDValue(0,0); auto S_yy = S.DDValue(1,1);
      auto S_xy = S.DDValue(0,1); auto S_yz = S.DDValue(1,2); 
      auto S_xz = S.DDValue(0,2); auto S_zz = S.DDValue(2,2);
      /////////////////////


      //return - 1.0/3 * Vec<3,T>(0.0,0.0,0.0);
      return -1.0/3.0 * Vec<3,T> (S_xx*(b0_z-b3_y)+S_x*(b0_zx-b3_yx) + S_xy*(b3_x-b0_z)+S_y*(b3_xx-b0_zx) + S_xz*(b0_y-b0_x)+S_z*(b0_yx-b0_xx),
				  S_xy*(b0_z-b3_y)+S_x*(b0_yz-b3_yy) + S_yy*(b3_x-b0_z)+S_y*(b3_yx-b0_yz) + S_yz*(b0_y-b0_x)+S_z*(b0_yy-b0_yx),
				  S_xz*(b0_z-b3_y)+S_x*(b0_zz-b3_yz) + S_yz*(b3_x-b0_z)+S_y*(b3_zx-b0_zz) + S_zz*(b0_y-b0_x)+S_z*(b0_yz-b0_zx));
    }

    Vec<3,T> CurlShape()
    {     
      throw Exception("not implemented for GG-bubbles");
    }
    
  };

  template <int D, typename T>
  auto GGbubble_B3 (AutoDiffDiff<D,T> aS, AutoDiffDiff<D,T> ab0, AutoDiffDiff<D,T> ab1, AutoDiffDiff<D,T> ab2, AutoDiffDiff<D,T> ab3) { return T_GGbubble_B3<D, T>(aS, ab0, ab1, ab2, ab3); }
    
   
  /* Face basis functions which are normal-tangential continuous */
  /* calculates [(grad l1) o-times (grad l2 x grad l3)] * legendre */
  /* DivShape assumes that phi_12 =  [(grad l1) o-times (grad l2 x grad l3)] is constant!!! */
  template <typename T>
  class T_Dl1_o_Dl2xDl3_v
  {
    AutoDiff<3,T> l1,l2,l3,v;
  public:
    T_Dl1_o_Dl2xDl3_v  (AutoDiff<3,T> lam1, AutoDiff<3,T> lam2, AutoDiff<3,T> lam3, AutoDiff<3,T> av) : l1(lam1), l2(lam2), l3(lam3), v(av) { ; }
    
    Vec<9,T> Shape() {

      T cross1 = l2.DValue(1)*l3.DValue(2) - l2.DValue(2)*l3.DValue(1);
      T cross2 = -(l2.DValue(0)*l3.DValue(2) - l2.DValue(2)*l3.DValue(0));
      T cross3 = l2.DValue(0)*l3.DValue(1) - l2.DValue(1)*l3.DValue(0);

      Vec<9,T> sigmaref;
      
      for (int i=0; i<3; i++)
	{
	  sigmaref(i*3)= v.Value() * l1.DValue(i) * cross1;
	  sigmaref(i*3+1)= v.Value() * l1.DValue(i) * cross2;
	  sigmaref(i*3+2)= v.Value() * l1.DValue(i) * cross3;
	}

      //T trace_sigma = 1.0/3 * (sigmaref(0) + sigmaref(4) + sigmaref(8));
      T trace_sigma = v.Value()/3 * (l1.DValue(0) * cross1 + l1.DValue(1) * cross2 + l1.DValue(2) * cross3);
      
      sigmaref(0) = sigmaref(0) - trace_sigma;
      sigmaref(4) = sigmaref(4) - trace_sigma;
      sigmaref(8) = sigmaref(8) - trace_sigma;
      
      return sigmaref;  

    }

    Vec<3,T> DivShape()
    {
      T vx = v.DValue(0), vy = v.DValue(1), vz = v.DValue(2);

      T cross1 = l2.DValue(1)*l3.DValue(2) - l2.DValue(2)*l3.DValue(1);
      T cross2 = -(l2.DValue(0)*l3.DValue(2) - l2.DValue(2)*l3.DValue(0));
      T cross3 = l2.DValue(0)*l3.DValue(1) - l2.DValue(1)*l3.DValue(0);

      T trace_sigma = 1.0/3 * (l1.DValue(0) * cross1 + l1.DValue(1) * cross2 + l1.DValue(2) * cross3);

      return Vec<3,T> (vx * l1.DValue(0) * cross1 + vy * l1.DValue(0) * cross2 + vz * l1.DValue(0) * cross3 - vx * trace_sigma  ,
		     vx * l1.DValue(1) * cross1 + vy * l1.DValue(1) * cross2 + vz * l1.DValue(1) * cross3 - vy * trace_sigma,
		     vx * l1.DValue(2) * cross1 + vy * l1.DValue(2) * cross2 + vz * l1.DValue(2) * cross3 - vz * trace_sigma
		     );
    }

    Vec<3,T> CurlShape()
    {     
      throw Exception("not implemented for T_Dl1_o_Dl2xDl3_v");
    }

  };

  
  /* Identity = Inner bubble function (normatl-tangential component is zero) */
  /* calculates I * legendre */
  template <int D, typename T>  class T_Id_v;
  template <typename T>  class T_Id_v<3,T>
  {
     AutoDiff<3,T> v;
  public:
    T_Id_v  (AutoDiff<3,T> av) : v(av) { ; }
    
    Vec<9,T> Shape() {
      T zero = 0.0;
      Vec<9,T> Id_v= zero;

      for (int i=0; i<3; i++)
	Id_v(i*(4))= v.Value();      
      return Id_v;
    }

    Vec<3,T> DivShape()
    {
      Vec<3,T> div_Id_v=0;
      for (int i=0; i<3; i++)
	div_Id_v(i) = v.DValue(i);
      return div_Id_v;
    }

    Vec<3,T> CurlShape()
    {     
      Vec<3, T> curl_Id_v=0;
      curl_Id_v(0) = v.DValue(2) -v.DValue(1);
      curl_Id_v(1) = v.DValue(0) -v.DValue(2);
      curl_Id_v(2) = v.DValue(1) -v.DValue(0); 
          
      return curl_Id_v;
    }

  };

  template <int D, typename T>
  auto Id_v (AutoDiff<D,T> av) { return T_Id_v<D,T>(av); }

}


#endif
  
