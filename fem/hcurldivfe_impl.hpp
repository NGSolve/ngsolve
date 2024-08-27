#ifndef FILE_HCURLDIVFE_IMPL
#define FILE_HCURLDIVFE_IMPL

/**************************************************************************/
/* File:   hcurldivfe_impl.hpp                                            */
/* Author: Philip Lederer                                                 */
/* Date:   2017/2018                                                      */
/**************************************************************************/


#include "shapefunction_utils.hpp"

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
  /* Curl(grad(u) v) = Curl(grad(u)) v + grad(u) o-tiimes Curl(v) */
  /* For basis functions including the trace */
    
  template <int D, typename T> class T_Sigma_gradu_v;
  template <typename T> class T_Sigma_gradu_v<2,T>
  {
    AutoDiffDiff<2,T> u;
    AutoDiffDiff<2,T> v;
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
    AutoDiffDiff<2,T> u;
    AutoDiffDiff<2,T> v;
    double tr;
  public:
    T_Gradu_Curlv  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av, double atr) : u(au), v(av), tr(atr){ ; }

    Vec<4,T> Shape() {

      auto trace = tr * (-  v.DValue(1)*u.DValue(0) + v.DValue(0)*u.DValue(1)) / 2.0;
      
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
      
      return Vec<2,T> (-vy*uxx + vx*uxy - tr * 0.5*(-vxy*ux - vy*uxx + vxx*uy + vx*uxy),-vy*uxy + vx * uyy - tr*0.5*(-vyy*ux - vy*uxy + vxy*uy + vx*uyy));
    }

    Vec<2,T> CurlShape()
    {
      T vxx = v.DDValue(0,0), vxy = v.DDValue(0,1), vyy = v.DDValue(1,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      //T uxx = u.DDValue(0,0), uxy = u.DDValue(0,1), uyy = u.DDValue(1,1);
      //T vx = v.DValue(0), vy = v.DValue(1);

      return Vec<2,T> (vyy*ux - vxy*uy,-vxy*ux+vxx*uy);     
    }    
  };

  template <int D, typename T>
  auto Gradu_Curlv (AutoDiffDiff<D,T> au,AutoDiffDiff<D,T> av, double atr) { return T_Gradu_Curlv<D, T>(au,av,atr); }

  
  /* ############### Type 2 (QUAD) - inner basis functions ############### */
  /* u * sigma(grad v) = u * Curl(grad v), where Curl is the 1D to 2D curl operator */
  template <int D, typename T> class T_u_Sigma_gradv;
  template <typename T> class T_u_Sigma_gradv<2,T>
  {
    AutoDiffDiff<2,T> u;
    AutoDiffDiff<2,T> v;
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
    AutoDiffDiff<2,T> u;
    AutoDiffDiff<2,T> v;
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
      //auto Bxx = b.DDValue(0,0);
      //auto Bxy = b.DDValue(0,1);
      //auto Byy = b.DDValue(1,1);

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
  
    template <typename T>
  INLINE AutoDiffDiff<3,T> Cross (const AutoDiffDiff<3,T> & u,
                                 const AutoDiffDiff<3,T> & v)
  {
    AutoDiffDiff<3,T> ret;
    ret.Value()=0;
    ret.DValue(0) = u.DValue(1)*v.DValue(2)-u.DValue(2)*v.DValue(1);
    ret.DValue(1) = u.DValue(2)*v.DValue(0)-u.DValue(0)*v.DValue(2);
    ret.DValue(2) = u.DValue(0)*v.DValue(1)-u.DValue(1)*v.DValue(0);

    ret.DDValue(0,0) = u.DDValue(0,1)*v.DValue(2)-u.DDValue(0,2)*v.DValue(1) + u.DValue(1)*v.DDValue(0,2)-u.DValue(2)*v.DDValue(0,1);
    ret.DDValue(0,1) = u.DDValue(1,1)*v.DValue(2)-u.DDValue(1,2)*v.DValue(1) + u.DValue(1)*v.DDValue(1,2)-u.DValue(2)*v.DDValue(1,1);
    ret.DDValue(0,2) = u.DDValue(2,1)*v.DValue(2)-u.DDValue(2,2)*v.DValue(1) + u.DValue(1)*v.DDValue(2,2)-u.DValue(2)*v.DDValue(2,1);

    ret.DDValue(1,0) = u.DDValue(0,2)*v.DValue(0)-u.DDValue(0,0)*v.DValue(2) + u.DValue(2)*v.DDValue(0,0)-u.DValue(0)*v.DDValue(0,2);
    ret.DDValue(1,1) = u.DDValue(1,2)*v.DValue(0)-u.DDValue(1,0)*v.DValue(2) + u.DValue(2)*v.DDValue(1,0)-u.DValue(0)*v.DDValue(1,2);
    ret.DDValue(1,2) = u.DDValue(2,2)*v.DValue(0)-u.DDValue(2,0)*v.DValue(2) + u.DValue(2)*v.DDValue(2,0)-u.DValue(0)*v.DDValue(2,2);

    ret.DDValue(2,0) = u.DDValue(0,0)*v.DValue(1)-u.DDValue(0,1)*v.DValue(0) + u.DValue(0)*v.DDValue(0,1)-u.DValue(1)*v.DDValue(0,0);
    ret.DDValue(2,1) = u.DDValue(1,0)*v.DValue(1)-u.DDValue(1,1)*v.DValue(0) + u.DValue(0)*v.DDValue(1,1)-u.DValue(1)*v.DDValue(1,0);
    ret.DDValue(2,2) = u.DDValue(2,0)*v.DValue(1)-u.DDValue(2,1)*v.DValue(0) + u.DValue(0)*v.DDValue(2,1)-u.DValue(1)*v.DDValue(2,0);
    
    return ret;
  }

  template <typename T>
  INLINE Vec<3,AutoDiff<3,T>> Cross (const AutoDiffDiff<3,T> u,
					  const Vec<3,T> v)
  {
    AutoDiff<3,T> hv[3];
    hv[0].Value() = u.DValue(1)*v(2)-u.DValue(2)*v(1);
    hv[1].Value() = u.DValue(2)*v(0)-u.DValue(0)*v(2);
    hv[2].Value() = u.DValue(0)*v(1)-u.DValue(1)*v(0);

    for (int i = 0; i<3; i++)
      {
	hv[0].DValue(i) = u.DDValue(i,1)*v(2)-u.DDValue(i,2)*v(1);
	hv[1].DValue(i) = u.DDValue(i,2)*v(0)-u.DDValue(i,0)*v(2);
	hv[2].DValue(i) = u.DDValue(i,0)*v(1)-u.DDValue(i,1)*v(0);
      }
    
    return Vec<3,AutoDiff<3,T>>(hv[0],hv[1],hv[2]);
  }

  // GG-bubble of Jay
  // computes curl ( curl A B) where A is a skew sym matrix which is L^2 orthogonal on P^k-1
  // and B is the matrix bubble
  
  template <int D, typename T> class T_GGbubble_3D;
  template <typename T> class T_GGbubble_3D<3,T>
  {
    AutoDiffDiff<3,T> q;
    Mat<3,3,T> S;    
    Mat<3,3,T> B;
    AutoDiffDiff<3,T> *curlB;
    
  public:
    T_GGbubble_3D  (AutoDiffDiff<3,T> aq, Mat<3,3,T> aS, Mat<3,3,T> aB, AutoDiffDiff<3,T> acurlB[] ) : q(aq), S(aS), B(aB), curlB(acurlB){ ; }
    Vec<9,T> Shape() {
      Vec<3,AutoDiff<3,T>> grad_q_cross_Si;
      Vec<9,T> sigmaref;
           
      for (int i = 0; i < 3; i++)
	{
	  grad_q_cross_Si = Cross (q, Vec<3,T>(S(i,0),S(i,1),S(i,2)));
	  
	  sigmaref(i*3) = 0;
	  sigmaref(i*3+1) = 0;
	  sigmaref(i*3+2) = 0;
	  for (int j = 0; j < 3; j++)
	    {
	      Vec<3,T> hv1(grad_q_cross_Si(j).DValue(0),grad_q_cross_Si(j).DValue(1),grad_q_cross_Si(j).DValue(2));
	      Vec<3,T> hv2(B(j,0),B(j,1),B(j,2));
	      Vec<3,T> crossres = Cross(hv1,hv2);
	      AutoDiffDiff<3,T> res2 = grad_q_cross_Si(j).Value() * curlB[j];
	      sigmaref(i*3)   += crossres(0) + res2.DValue(0);
	      sigmaref(i*3+1) += crossres(1) + res2.DValue(1);
	      sigmaref(i*3+2) += crossres(2) + res2.DValue(2);
	    }
	}
      
      T sigma_trace = 1/3.0 * (sigmaref(0) + sigmaref(4) + sigmaref(8));	          
      
      sigmaref(0) -= sigma_trace;
      sigmaref(4) -= sigma_trace;
      sigmaref(8) -= sigma_trace;
      
      return sigmaref;      
    }

    Vec<3,T> DivShape()
    {
      T grad_x=0.0;
      T grad_y=0.0;
      T grad_z=0.0;

      Vec<3,AutoDiff<3,T>> grad_q_cross_Si;

      for(int i = 0; i<3; i++)
      	{
	  grad_q_cross_Si = Cross (q, Vec<3,T>(S(i,0),S(i,1),S(i,2)));
	  
	  for(int j = 0; j<3; j++)
	    {
	      grad_x += grad_q_cross_Si(j).DValue(0) * curlB[j].DValue(i) + grad_q_cross_Si(j).Value() * curlB[j].DDValue(i,0);
	      grad_y += grad_q_cross_Si(j).DValue(1) * curlB[j].DValue(i) + grad_q_cross_Si(j).Value() * curlB[j].DDValue(i,1);
	      grad_z += grad_q_cross_Si(j).DValue(2) * curlB[j].DValue(i) + grad_q_cross_Si(j).Value() * curlB[j].DDValue(i,2);
	    }	      
      	}
          
      return -1.0/3.0 * Vec<3,T>(grad_x,grad_y,grad_z);
    }

    Vec<3,T> CurlShape()
    {     
      throw Exception("not implemented for GG-bubbles");
    }
    
  };

  template <int D, typename T>
  auto GGbubble_3D (AutoDiffDiff<D,T> aq, Mat<D,D,T> aS, Mat<D,D,T> aB, AutoDiffDiff<D,T> acurlB[]) { return T_GGbubble_3D<D, T>(aq, aS, aB, acurlB); }

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

    Vec<9,T> CurlShape()
    {
      Vec<9,T> res;
      auto vl1 = Cross (v, l1);
      auto l2l3 = Cross (l2,l3);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          res(3*i+j) = vl1.DValue(i)*l2l3.DValue(j);

      // 1/3 curl (tr A I) = 1/3  hat(grad(v)) * det (grad l1, grad l2, grad l3)  ???
      auto det = -Dot (l1, Cross(l2, l3));
      res(5) += det/3 * v.DValue(0);
      res(7) -= det/3 * v.DValue(0);
      
      res(2) -= det/3 * v.DValue(1);
      res(6) += det/3 * v.DValue(1);
      
      res(1) += det/3 * v.DValue(2);
      res(3) -= det/3 * v.DValue(2);
      return res;
      // throw Exception("not implemented for T_Dl1_o_Dl2xDl3_v");
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
      Vec<3,T> div_Id_v;
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

  //   template <int D, typename T>  class T_Id_v_DD;
  // template <typename T>  class T_Id_v_DD<2,T>
  // {
  //    AutoDiffDiff<2,T> v;
  // public:
  //   T_Id_v_DD  (AutoDiffDiff<2,T> av) : v(av) { ; }
    
  //   Vec<4,T> Shape() {
  //     T zero = 0.0;
  //     Vec<4,T> Id_v= zero;

  //     for (int i=0; i<2; i++)
	//       Id_v(i*(2))= v.Value();      
  //     return Id_v;
  //   }

  //   Vec<2,T> DivShape()
  //   {
  //     Vec<2,T> div_Id_v;
  //     for (int i=0; i<2; i++)
	//         div_Id_v(i) = v.DValue(i);
  //     return div_Id_v;
  //   }

  //   Vec<2,T> CurlShape()
  //   {     
  //     throw Exception("not implemented for T_Id_v and D = 2");
  //   }

  // };
  // template <int D, typename T>
  // auto Id_v_DD (AutoDiffDiff<D,T> av) { return T_Id_v_DD<D,T>(av); }
  

  /* ############### (HEX) - edge basis functions ############### */
  /* calculate legendre * dev((grad l1) o-times (grad l2)) */
  template <typename T>
  class T_dev_Dl1_o_Dl2_v
  {
    AutoDiff<3,T> l1,l2,v;
  public:
    T_dev_Dl1_o_Dl2_v  (AutoDiff<3,T> lam1, AutoDiff<3,T> lam2, AutoDiff<3,T> av) : l1(lam1), l2(lam2), v(av) { ; }
    
    Vec<9,T> Shape() {

      Vec<9,T> sigmaref;

      for (int i=0; i<3; i++)
	{
	  sigmaref(i*3)= v.Value() * l1.DValue(i) * l2.DValue(0);
	  sigmaref(i*3+1)= v.Value() * l1.DValue(i) * l2.DValue(1);
	  sigmaref(i*3+2)= v.Value() * l1.DValue(i) * l2.DValue(2);
	}

      T trace_sigma = v.Value()/3.0 * (l1.DValue(0) * l2.DValue(0) + l1.DValue(1) * l2.DValue(1) + l1.DValue(2) * l2.DValue(2));
      
      sigmaref(0) = sigmaref(0) - trace_sigma;
      sigmaref(4) = sigmaref(4) - trace_sigma;
      sigmaref(8) = sigmaref(8) - trace_sigma;
      
      return sigmaref;  

    }

    Vec<3,T> DivShape()
    {
      T vx = v.DValue(0), vy = v.DValue(1), vz = v.DValue(2);

      T trace_sigma = 1.0/3 * (l1.DValue(0) * l2.DValue(0) + l1.DValue(1) * l2.DValue(1) + l1.DValue(2) * l2.DValue(2));

      return Vec<3,T> (vx * l1.DValue(0) * l2.DValue(0) + vy * l1.DValue(0) * l2.DValue(1) + vz * l1.DValue(0) * l2.DValue(2) - vx * trace_sigma  ,
		       vx * l1.DValue(1) * l2.DValue(0) + vy * l1.DValue(1) * l2.DValue(1) + vz * l1.DValue(1) * l2.DValue(2) - vy * trace_sigma,
		       vx * l1.DValue(2) * l2.DValue(0) + vy * l1.DValue(2) * l2.DValue(1) + vz * l1.DValue(2) * l2.DValue(2) - vz * trace_sigma
		     );
    }

    Vec<3,T> CurlShape()
    {     
      throw Exception("not implemented for T_dev_Dl1_o_Dl2_v");
    }

  };
}


#endif
  
