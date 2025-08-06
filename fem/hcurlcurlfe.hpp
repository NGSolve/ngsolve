#ifndef FILE_HCURLCURLFE
#define FILE_HCURLCURLFE

/*********************************************************************/
/* File:   hcurlcurlfe.hpp                                           */
/* Author: Michael Neunteufel                                        */
/* Date:   June 2018                                                 */
/*********************************************************************/


#include "finiteelement.hpp"
#include "fe_interfaces.hpp"
#include "hcurlfe.hpp"
#include "hcurlfe_utils.hpp"
#include "recursive_pol.hpp"
#include "recursive_pol_trig.hpp"
#include "recursive_pol_tet.hpp"

namespace ngfem
{
  
  template <typename T>
  Mat<3,3,T> TensorCrossProduct(Mat<3,3,T> A, Mat<3,3,T> B)
  {
    // return 0.5 * ( Cof(A+B) - Cof(A-B) ); // more cancelation

    Mat<3,3,T> prod;
    prod.Col(0) = Cross(A.Col(1), B.Col(2)) - Cross(A.Col(2), B.Col(1));
    prod.Col(1) = Cross(A.Col(2), B.Col(0)) - Cross(A.Col(0), B.Col(2));
    prod.Col(2) = Cross(A.Col(0), B.Col(1)) - Cross(A.Col(1), B.Col(0));
    return prod;
  }
  
  template <typename T>
  Mat<3,3,T> TensorCrossProduct(Vec<3,T> v, Mat<3,3,T> A)
  {
    Mat<3,3,T> result;
    for (int j = 0; j < 3; j++)
      result.Col(j) = Cross(v, A.Col(j));
    return result;
  }

  template <typename T>
  Mat<3,3,T> TensorCrossProduct(Mat<3,3,T> A, Vec<3,T> v)
  {
    Mat<3,3,T> result;
    for (int j = 0; j < 3; j++)
      result.Row(j) = Cross(A.Row(j), v);
    return result;
  }
  
  

  template <int DIM>
  class HCurlCurlFiniteElement : public FiniteElement
  {
  public:
    using FiniteElement::FiniteElement;
    using FiniteElement::ndof;
    using FiniteElement::order;
    
    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
                                  BareSliceMatrix<double> shape) const = 0;

    virtual void EvaluateMappedShape (const BaseMappedIntegrationPoint & bmip,
                                      BareSliceVector<double> coefs,
                                      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedCurlShape (const BaseMappedIntegrationPoint & bmip,
      BareSliceMatrix<double> shape) const = 0;

    virtual void CalcMappedIncShape (const BaseMappedIntegrationPoint & bmip,
                                     BareSliceMatrix<double> shape) const = 0;

    virtual void EvaluateMappedIncShape (const BaseMappedIntegrationPoint & bmip,
                                         BareSliceVector<double> coefs,
                                         BareSliceVector<double> inc) const = 0;

    virtual void CalcMappedIncShape (const SIMD_BaseMappedIntegrationRule & bmir,
                                     BareSliceMatrix<SIMD<double>> shape) const = 0;

    virtual void EvaluateIncShape (const SIMD_BaseMappedIntegrationRule & ir,
                                   BareSliceVector<> coefs,
                                   BareSliceMatrix<SIMD<double>> values) const = 0;

    virtual void AddTransIncShape (const SIMD_BaseMappedIntegrationRule & ir,
                                   BareSliceMatrix<SIMD<double>> values,
                                   BareSliceVector<> coefs) const = 0;

    
    virtual void CalcMappedShape (const SIMD<BaseMappedIntegrationPoint> & bmip, 
                                         BareSliceMatrix<SIMD<double>> shapes) const = 0;

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                                         BareSliceMatrix<SIMD<double>> shapes) const = 0;
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                                  BareSliceVector<> coefs,
                                  BareSliceMatrix<SIMD<double>> values) const = 0;

    virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir,
                                  BareSliceMatrix<SIMD<double>> values,
                                  BareSliceVector<> coefs) const = 0;

    virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const = 0;
    virtual void CalcDualShape (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> shape) const = 0;
    virtual void EvaluateDual (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const = 0;
    virtual void AddDualTrans (const SIMD_BaseMappedIntegrationRule& bmir, BareSliceMatrix<SIMD<double>> values, BareSliceVector<double> coefs) const = 0;

  };

  template <int D,typename VEC,typename MAT>
  void VecToSymMat(const VEC & vec, MAT & mat)
  {
    switch(D)
    {
    case 1:
      mat(0) = vec(0);
      break;
    case 2:
      mat(0) = vec(0);
      mat(3) = vec(1);
      mat(1) = mat(2) = vec(2);
      break;
    case 3:
      mat(0) = vec(0);
      mat(4) = vec(1);
      mat(8) = vec(2);
      mat(1) = mat(3) = vec(5);
      mat(2) = mat(6) = vec(4);
      mat(5) = mat(7) = vec(3);
      break;
    }
  }

  template <int H, int W, typename T>
  Mat<H,W,T> DyadProd(Vec<H,T> a, Vec<W,T> b)
  {
    Mat<H,W,T> m;
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        m(i,j) = a(i)*b(j);
    return m;
  }
  
  template <int S, typename T>
  Mat<S,S,T> SymDyadProd(Vec<S,T> a, Vec<S,T> b)
  {
    Mat<S,S,T> m;
    for (int i = 0; i < S; i++)
      for (int j = 0; j < S; j++)
        m(i,j) = a(i)*b(j)+a(j)*b(i);
    return m;
  }


  template <typename T>
  Vec<6, AutoDiff<3,T>> SymDyadProd(AutoDiff<3,T> a, AutoDiff<3,T> b)
  {
    return Vec<6, AutoDiff<3,T>>(2*a.DValue(0)*b.DValue(0),2*a.DValue(1)*b.DValue(1),2*a.DValue(2)*b.DValue(2), a.DValue(1)*b.DValue(2)+a.DValue(2)*b.DValue(1), a.DValue(0)*b.DValue(2)+a.DValue(2)*b.DValue(0),a.DValue(1)*b.DValue(0)+a.DValue(0)*b.DValue(1));
  }

  template <typename T>
  Vec<6, AutoDiff<3,T>> SymDyadProdAD(Vec<3,T> a, Vec<3,T> b)
  {
    return Vec<6, AutoDiff<3,T>>(2*a(0)*b(0),2*a(1)*b(1),2*a(2)*b(2), a(1)*b(2)+a(2)*b(1), a(0)*b(2)+a(2)*b(0),a(1)*b(0)+a(0)*b(1));
  }

  template <typename T>
  Vec<3,AutoDiff<2,T>> SymDyadProd(AutoDiff<2,T> a, AutoDiff<2,T> b)
  {
    return Vec<3,AutoDiff<2,T>>(2*a.DValue(0)*b.DValue(0),2*a.DValue(1)*b.DValue(1),a.DValue(1)*b.DValue(0)+a.DValue(0)*b.DValue(1));
  }

  template <typename T>
  Vec<3,AutoDiff<2,T>> SymDyadProdAD(Vec<2,T> a, Vec<2,T> b)
  {
    return Vec<3,AutoDiff<2,T>>(2*a(0)*b(0),2*a(1)*b(1),a(1)*b(0)+a(0)*b(1));
  }

  template <typename T>
  AutoDiff<1,T> SymDyadProd(AutoDiff<1,T> a, AutoDiff<1,T> b)
  {
    return a.DValue(0)*b.DValue(0);
  }


  //------------------REGGE_SHAPE---------------------
  template <int D, typename T> class T_REGGE_Shape;
  template <typename T> class T_REGGE_Shape<1,T>
  {
    AutoDiff<1,T> u;
  public:
    T_REGGE_Shape  (AutoDiff<1,T> au) : u(au) { ; }
    Vec<1,T> Shape() { return u.Value(); }
    /*0 2
      2 1*/
    Vec<1,T> CurlShape() { return 0.0; }
  };
  
  template <typename T> class T_REGGE_Shape<2,T>
  {
    Vec<3,AutoDiff<2,T>> u;
  public:
    T_REGGE_Shape  (Vec<3,AutoDiff<2,T>> au) : u(au) { ; }
    Vec<3,T> Shape() { return Vec<3,T> (u(0).Value(), u(1).Value(), u(2).Value()); }
    /*0 2
      2 1*/
    Vec<2,T> CurlShape() { return Vec<2,T> (u(2).DValue(0)-u(0).DValue(1), u(1).DValue(0)-u(2).DValue(1)); }
  };
  
  template <typename T> class T_REGGE_Shape<3,T>
  {
    Vec<6,AutoDiff<3,T>> u;
  public:
    T_REGGE_Shape  (Vec<6,AutoDiff<3,T>> au) : u(au) { ; }
    Vec<6,T> Shape() { return Vec<6,T> (u(0).Value(), u(1).Value(), u(2).Value(), u(3).Value(), u(4).Value(), u(5).Value()); }
    /*0 5 4
      5 1 3
      4 3 2*/
    Vec<9,T> CurlShape() { return Vec<9,T> (u(4).DValue(1)-u(5).DValue(2), -u(4).DValue(0)+u(0).DValue(2), u(5).DValue(0)-u(0).DValue(1),
					    u(3).DValue(1)-u(1).DValue(2), -u(3).DValue(0)+u(5).DValue(2), u(1).DValue(0)-u(5).DValue(1),
					    u(2).DValue(1)-u(3).DValue(2), -u(2).DValue(0)+u(4).DValue(2), u(3).DValue(0)-u(4).DValue(1)); }
  };
  //---------------------------------------------------


    // ***************** EpsGrad ****************************** */
  // eps (nabla u)
  
  template <int D, typename T> class T_EpsGrad;
  template <typename T> class T_EpsGrad<2,T>
  {
    AutoDiffDiff<2,T> u;
  public:
    T_EpsGrad  (AutoDiffDiff<2,T> au) : u(au) { ; }
    Vec<3,T> Shape()
    {
      return Vec<3,T> (u.DDValue(0,0), u.DDValue(1,1), u.DDValue(0,1));
    }
    Vec<2,T> CurlShape() { return Vec<2,T> (0.0, 0.0); }
  };
  
  template <int D, typename T>
  auto EpsGrad (AutoDiffDiff<D,T> au) { return T_EpsGrad<D,T>(au); }

  // ***************** wEpsGrad ****************************** */
  // w*eps (nabla u)
  
  template <int D, typename T> class T_wEpsGrad;
  template <typename T> class T_wEpsGrad<2,T>
  {
    AutoDiffDiff<2,T> u;
    AutoDiff<1,T> w;
  public:
    T_wEpsGrad  (AutoDiffDiff<2,T> au, AutoDiff<1,T> aw) : u(au), w(aw) { ; }
    Vec<6,T> Shape()
    {
      return w.Value()*Vec<6,T> (u.DDValue(0,0), u.DDValue(1,1), u.DDValue(2,2), u.DDValue(1,2), u.DDValue(0,2), u.DDValue(0,1));
    }
    Vec<9,T> CurlShape() { return Vec<9,T> (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); }
  };
  
  template <int D, typename T>
  auto wEpsGrad (AutoDiffDiff<D,T> au, AutoDiff<1,T> aw) { return T_wEpsGrad<D,T>(au, aw); }
  
  
  // ***************** Eps_u_Gradv ****************************** */
  // eps (u nabla v)
  
  template <int D, typename T> class T_Eps_u_Gradv;
  template <typename T> class T_Eps_u_Gradv<2,T>
  {
    AutoDiffDiff<2,T> u, v;
  public:
    T_Eps_u_Gradv  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av) { ; }
    Vec<3,T> Shape() { return Vec<3,T> ((u.Value()*v.DDValue(0,0) + u.DValue(0)*v.DValue(0)),
                                        (u.Value()*v.DDValue(1,1) + u.DValue(1)*v.DValue(1)),
                                        u.Value()*v.DDValue(0,1) + 0.5 * (u.DValue(0)*v.DValue(1)+u.DValue(1)*v.DValue(0))); }
    Vec<2,T> CurlShape()
    {
      T uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      T ux = u.DValue(0), uy = u.DValue(1);
      T vxx = v.DDValue(0,0), vyy = v.DDValue(1,1), vxy = v.DDValue(0,1);
      T vx = v.DValue(0), vy = v.DValue(1);
      
      /*return -0.5 * Vec<2,T> (uyy*vx - uxy*vy + uy*vxy - ux*vyy,
      -uxy*vx + uxx*vy - uy*vxx + ux*vxy);*/
      return 0.5 * Vec<2,T>(ux*vxy - uy*vxx - uxy*vx + uxx*vy,
                            ux*vyy + uxy*vy - uyy*vx - uy*vxy);
    }
  };
  
  template <int D, typename T>
  auto Eps_u_Gradv (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_Eps_u_Gradv<D,T>(au, av); }
  
  
  template <int D, typename T> class T_vEpsGradu;
  template <typename T> class T_vEpsGradu<2,T>
  {
    AutoDiffDiff<2,T> u,v;
  public:
    T_vEpsGradu  (AutoDiffDiff<2,T> au, AutoDiffDiff<2,T> av) : u(au), v(av) { ; }
    Vec<3,T> Shape() { return Vec<3,T> (u.DDValue(0,0)*v.Value(),
                                      u.DDValue(1,1)*v.Value(),  (u.DDValue(1,0)*v.Value()));}
    Vec<2,T> CurlShape()
    {
      T uxx = u.DDValue(0,0), uyy = u.DDValue(1,1), uxy = u.DDValue(0,1);
      T vx = v.DValue(0), vy = v.DValue(1);

      //return Vec<2,T> (uyy*vx- uxy*vy, uxx*vy- uxy*vx);
      return Vec<2,T> (uxy*vx - vy*uxx, uyy*vx - uxy*vy);
    }
  };
  
  template <int D, typename T>
  auto vEpsGradu (AutoDiffDiff<D,T> au, AutoDiffDiff<D,T> av) { return T_vEpsGradu<D,T>(au, av); }

  
  template <int D, typename T>
  class ReggeAD
  {
    Mat<D,D,T> value;

  public:
    ReggeAD ()
    {
      value = T(0);
    }
    
    ReggeAD (AutoDiff<D,T> a, AutoDiff<D,T> b)
    {
      Vec<D,T> Da, Db;
      for(int i=0; i<D; i++)
        {
          Da(i) = a.DValue(i);
          Db(i) = b.DValue(i);
        }      
      value = SymDyadProd(Da,Db);
    }

    auto Value() const { return value; }

    Mat<D,D,T> & Value() { return value; }
  };

  template <int D, typename T>
  auto MakeReggeAD(AutoDiff<D,T> a, AutoDiff<D,T> b)
  {
    return ReggeAD<D,T>(a, b);
  }



  template <int D, typename T>
  ReggeAD<D,T> operator* (AutoDiff<D,T> s, ReggeAD<D,T> A)
  {
    ReggeAD<D,T> result;
    result.Value() = s.Value()*A.Value();
    return result;
  }

  
  template <int D, typename T>
  ReggeAD<D,T> operator* (T s, ReggeAD<D,T> A)
  {
    ReggeAD<D,T> result = A;
    result.Value() *= s;
    return result;
  }

  template <int D, typename T>
  ReggeAD<D,T> operator+ (ReggeAD<D,T> A, ReggeAD<D,T> B)
  {
    ReggeAD<D,T> result = A;
    result.Value() += B.Value();
    return result;
  }

  template <int D, typename T>
  ReggeAD<D,T> operator- (ReggeAD<D,T> A, ReggeAD<D,T> B) 
  {
    ReggeAD<D,T> result = A;
    result.Value() -= B.Value();
    return result;
  }
  
  template <int D, typename T> class ReggeADD;

  template <typename T>
  class ReggeADD<3,T>
  {
    Mat<3,3,T> value;
    Mat<3,3,T> curl;
    Mat<3,3,T> inc;

  public:
    ReggeADD ()
    {
      value = T(0);
      curl = T(0);
      inc = T(0);
    }
    
    ReggeADD (AutoDiffDiff<3,T> a, AutoDiffDiff<3,T> b)
    {
      auto Da = Vec<3,T>(a.DValue(0),a.DValue(1),a.DValue(2));
      auto Db = Vec<3,T>(b.DValue(0),b.DValue(1),b.DValue(2));
      value = SymDyadProd(Da,Db);
      // curl(s*v) = nabla s x v + s curl(v) in 3D
      //                     | nabla a_1 x b + a_1 curl(b) |
      // curl(a \otimes b) = | nabla a_2 x b + a_2 curl(b) |
      //                     | nabla a_3 x b + a_3 curl(b) |

      // curl( nabla s ) = 0
      //                                      | nabla d_x a x nabla b |
      // ->  curl( nabla a \otimes nabla b) = | nabla d_y a x nabla b |
      //                                      | nabla d_z a x nabla b |

      Vec<3,T> Ddai [3] = { Vec<3,T>(a.DDValue(0,0), a.DDValue(0,1), a.DDValue(0,2)), Vec<3,T>(a.DDValue(1,0), a.DDValue(1,1), a.DDValue(1,2)), Vec<3,T>(a.DDValue(2,0), a.DDValue(2,1), a.DDValue(2,2)) };
      Vec<3,T> Ddbi [3] = { Vec<3,T>(b.DDValue(0,0), b.DDValue(0,1), b.DDValue(0,2)), Vec<3,T>(b.DDValue(1,0), b.DDValue(1,1), b.DDValue(1,2)), Vec<3,T>(b.DDValue(2,0), b.DDValue(2,1), b.DDValue(2,2)) };

      for (int i = 0; i < 3; i++)
        curl.Row(i) = Cross(Ddai[i], Db) + Cross(Ddbi[i], Da);


      //                                   | nabla d_x a x nabla b |
      //  curl( nabla a \otimes nabla b) = | nabla d_y a x nabla b |
      //                                   | nabla d_z a x nabla b |


      // curl T curl ( nabla a \otimes nabla b):
      //11: d_yd_z(a) d_yd_z(b) - d^2_z(a) d^2_y(b)   - d_y^2(a) d^2_z(b)   + d_yd_z(a)d_yd_z(b)
      //12: d_xd_y(a) d_z^2(b)  - d_xd_z(a) d_yd_z(b) - d_yd_z(a) d_xd_z(b) + d_z^2(a)d_xd_y(b)
      //13: d^2_y(a) d_xd_z(b)  - d_yd_z(a) d_xd_y(b) - d_xd_y(a) d_yd_z(b) + d_xd_z(a)d^2_y(b)
      //22: d_xd_z(a) d_xd_z(b) - d^2_x(a) d^2_z(b)   - d^2_z(a) d^2_x(b)   + d_xd_z(a)d_xd_z(b)
      //23: d_yd_z(a) d^2_x(b)  - d_xd_y(a) d_xd_z(b) - d_xd_z(a) d_xd_y(b) + d^2_x(a)d_yd_z(b)
      //33: d^2_x(a) d^2_y(b)   - d_xd_y(a) d_xd_y(b) - d_xd_y(a) d_xd_y(b) + d^2_y(a)d^2_x(b)

      //  = - hesse(a) x hesse(b) = -eps_imn eps_jlk d_md_l(a) d_nd_k b ??

      /*inc(0,0) = a.DDValue(1,2)*b.DDValue(1,2) - a.DDValue(2,2)*b.DDValue(1,1) - a.DDValue(1,1)*b.DDValue(2,2) + a.DDValue(1,2)*b.DDValue(1,2);
      inc(0,1) = a.DDValue(0,1)*b.DDValue(2,2) - a.DDValue(0,2)*b.DDValue(1,2) - a.DDValue(1,2)*b.DDValue(0,2) + a.DDValue(2,2)*b.DDValue(0,1);
      inc(0,2) = a.DDValue(1,1)*b.DDValue(0,2) - a.DDValue(1,2)*b.DDValue(0,1) - a.DDValue(0,1)*b.DDValue(1,2) + a.DDValue(0,2)*b.DDValue(1,1);
      inc(1,1) = a.DDValue(0,2)*b.DDValue(0,2) - a.DDValue(0,0)*b.DDValue(2,2) - a.DDValue(2,2)*b.DDValue(0,0) + a.DDValue(0,2)*b.DDValue(0,2);
      inc(1,2) = a.DDValue(1,2)*b.DDValue(0,0) - a.DDValue(0,1)*b.DDValue(0,2) - a.DDValue(0,2)*b.DDValue(0,1) + a.DDValue(0,0)*b.DDValue(1,2);
      inc(2,2) = a.DDValue(0,0)*b.DDValue(1,1) - a.DDValue(0,1)*b.DDValue(0,1) - a.DDValue(0,1)*b.DDValue(0,1) + a.DDValue(1,1)*b.DDValue(0,0);
      // curl T curl ( nabla b \otimes nabla a):
      inc(0,0) += b.DDValue(1,2)*a.DDValue(1,2) - b.DDValue(2,2)*a.DDValue(1,1) - b.DDValue(1,1)*a.DDValue(2,2) + b.DDValue(1,2)*a.DDValue(1,2);
      inc(0,1) += b.DDValue(0,1)*a.DDValue(2,2) - b.DDValue(0,2)*a.DDValue(1,2) - b.DDValue(1,2)*a.DDValue(0,2) + b.DDValue(2,2)*a.DDValue(0,1);
      inc(0,2) += b.DDValue(1,1)*a.DDValue(0,2) - b.DDValue(1,2)*a.DDValue(0,1) - b.DDValue(0,1)*a.DDValue(1,2) + b.DDValue(0,2)*a.DDValue(1,1);
      inc(1,1) += b.DDValue(0,2)*a.DDValue(0,2) - b.DDValue(0,0)*a.DDValue(2,2) - b.DDValue(2,2)*a.DDValue(0,0) + b.DDValue(0,2)*a.DDValue(0,2);
      inc(1,2) += b.DDValue(1,2)*a.DDValue(0,0) - b.DDValue(0,1)*a.DDValue(0,2) - b.DDValue(0,2)*a.DDValue(0,1) + b.DDValue(0,0)*a.DDValue(1,2);
      inc(2,2) += b.DDValue(0,0)*a.DDValue(1,1) - b.DDValue(0,1)*a.DDValue(0,1) - b.DDValue(0,1)*a.DDValue(0,1) + b.DDValue(1,1)*a.DDValue(0,0);
      // symmetry
      inc(1,0) = inc(0,1);
      inc(2,0) = inc(0,2);
      inc(2,1) = inc(1,2);*/

      Mat<3,3,T> hesse1, hesse2;
      a.StoreHessian(hesse1.Data());
      b.StoreHessian(hesse2.Data());
      inc = -2*TensorCrossProduct(hesse1,hesse2);
    }

    auto Value() const { return value; }
    auto Curl() const { return curl; }
    auto Inc() const { return inc; }

    Mat<3,3,T> & Value() { return value; }
    Mat<3,3,T> & Curl() { return curl; }
    Mat<3,3,T> & Inc() { return inc; }
  };


  template <typename T>
  class ReggeADD<2,T>
  {
    Mat<2,2,T> value;
    Vec<2,T> curl;
    T inc;

  public:
    ReggeADD ()
    {
      value = T(0);
      curl = T(0);
      inc = T(0);
    }
    
    ReggeADD (AutoDiffDiff<2,T> a, AutoDiffDiff<2,T> b)
    {
      auto Da = Vec<2,T>(a.DValue(0),a.DValue(1));
      auto Db = Vec<2,T>(b.DValue(0),b.DValue(1));
      value = SymDyadProd(Da,Db);
      // curl(s*v) = v* nabla s^perp + s curl(v) in 2D
      //                     | b * nabla a_1^perp + a_1 curl(b) |
      // curl(a \otimes b) = | b * nabla a_2^perp + a_2 curl(b) |

      // curl( nabla s ) = 0
      //                                      | Db * nabla d_x a^perp |
      // ->  curl( nabla a \otimes nabla b) = | Db * nabla d_y a^perp |

      Vec<2,T> Ddai_p [2] = { Vec<2,T>(-a.DDValue(0,1), a.DDValue(0,0)), Vec<2,T>(-a.DDValue(1,1), a.DDValue(1,0)) };
      Vec<2,T> Ddbi_p [2] = { Vec<2,T>(-b.DDValue(0,1), b.DDValue(0,0)), Vec<2,T>(-b.DDValue(1,1), b.DDValue(1,0)) };

      for (int i=0; i<2; i++)
        curl(i) = InnerProduct(Db,Ddai_p[i]) + InnerProduct(Da,Ddbi_p[i]);

      //                                  | Db * nabla d_x a^perp |
      // curl( nabla a \otimes nabla b) = | Db * nabla d_y a^perp |


      // curl T curl ( nabla a \otimes nabla b) = d_x(Db * nabla d_y a^perp) - d_y(Db * nabla d_x a^perp)
      //       =  (d_x Db) * nabla(d_y a)^perp - (d_y Db) * nabla(d_x a)^perp

      inc = InnerProduct(Vec<2,T>(b.DDValue(0,0), b.DDValue(0,1)),Ddai_p[1]) - InnerProduct(Vec<2,T>(b.DDValue(1,0), b.DDValue(1,1)),Ddai_p[0]) + InnerProduct(Vec<2,T>(a.DDValue(0,0), a.DDValue(0,1)),Ddbi_p[1]) - InnerProduct(Vec<2,T>(a.DDValue(1,0), a.DDValue(1,1)),Ddbi_p[0]);
    }

    auto Value() const { return value; }
    auto Curl() const { return curl; }
    auto Inc() const { return inc; }

    Mat<2,2,T> & Value() { return value; }
    Vec<2,T> & Curl() { return curl; }
    T & Inc() { return inc; }
  };
  
  template <int D, typename T>
  auto MakeReggeAD(AutoDiffDiff<D,T> a, AutoDiffDiff<D,T> b)
  {
    return ReggeADD<D,T>(a, b);
  }

  template <typename T>
  ReggeADD<3,T> operator* (AutoDiffDiff<3,T> s, ReggeADD<3,T> A)
  {
    ReggeADD<3,T> result;
    result.Value() = s.Value()*A.Value();

    // s scalar, v vector
    // curl(s*v) = nabla s x v + s curl(v) in 3D
    Vec<3,T> gradient;
    s.StoreGradient(gradient.Data());

    result.Curl() = s.Value()*A.Curl();
    for (int i = 0; i < 3; i++)
      result.Curl().Row(i) += Cross(gradient, Vec<3,T>(A.Value().Row(i)));
    
    // inc(s A) = s inc(A) + 2sym(grad(s) x curl A) + hesse(s) x A, x...Tensor-Cross-Product
    Mat<3,3,T> hesse;
    s.StoreHessian(hesse.Data());
    
    result.Inc() = s.Value()*A.Inc() + TensorCrossProduct(gradient,A.Curl()) + Trans(TensorCrossProduct(gradient,A.Curl())) +  TensorCrossProduct(hesse,A.Value());

    return result;
  }

  
  template <int D, typename T>
  ReggeADD<D,T> operator* (T s, ReggeADD<D,T> A)
  {
    ReggeADD<D,T> result = A;
    result.Value() *= s;
    result.Curl()  *= s;
    result.Inc()   *= s;
    return result;
  }

  template <int D, typename T>
  ReggeADD<D,T> operator+ (ReggeADD<D,T> A, ReggeADD<D,T> B)
  {
    ReggeADD<D,T> result = A;
    result.Value() += B.Value();
    result.Curl()  += B.Curl();
    result.Inc()   += B.Inc();
    return result;
  }

  template <int D, typename T>
  ReggeADD<D,T> operator- (ReggeADD<D,T> A, ReggeADD<D,T> B) 
  {
    ReggeADD<D,T> result = A;
    result.Value() -= B.Value();
    result.Curl()  -= B.Curl();
    result.Inc()   -= B.Inc();
    return result;
  }


  template <typename T>
  ReggeADD<2,T> operator* (AutoDiffDiff<2,T> s, ReggeADD<2,T> A)
  {
    ReggeADD<2,T> result;
    result.Value() = s.Value()*A.Value();

    // s scalar, v vector
    // curl(s*v) = v* nabla s^perp + s curl(v) in 2D
    result.Curl() = A.Value()*Vec<2,T>(-s.DValue(1),s.DValue(0)) + s.Value()*A.Curl();
    
    // inc(sA) = A:( dydy s & -dxdy s \\ -dxdy s & dxdx s) + 2*nabla s^\perp*curl(A) + s*inc(A)
    Mat<2,2,T> hesse;
    hesse(0,0) = s.DDValue(1,1);
    hesse(1,0) = -s.DDValue(1,0);
    hesse(0,1) = -s.DDValue(0,1);
    hesse(1,1) = s.DDValue(0,0);

    result.Inc() = s.Value()*A.Inc() + InnerProduct(hesse,A.Value()) + 2*InnerProduct(Vec<2,T>(-s.DValue(1),s.DValue(0)),A.Curl());  

    return result;
  }



  
  template <ELEMENT_TYPE ET> class HCurlCurlFE;

  
  template <ELEMENT_TYPE ET>
  class T_HCurlCurlFE : public HCurlCurlFiniteElement<ET_trait<ET>::DIM>,
    public VertexOrientedFE<ET>
  {
  protected:
    static constexpr int DIM = ET_trait<ET>::DIM;
    enum { DIM_STRESS = (DIM*(DIM+1))/2 };
    // enum { DIM_DMAT = 7*DIM-12 };
    // enum { DIM_DDMAT = 8*DIM-15 };
    enum { DIM_DMAT = (5*DIM*DIM-11*DIM+6)/2 };
    enum { DIM_DDMAT = (7*DIM*DIM-19*DIM+12)/2 };
    
    using VertexOrientedFE<ET>::vnums;
    using HCurlCurlFiniteElement<ET_trait<ET>::DIM>::ndof;
    using HCurlCurlFiniteElement<ET_trait<ET>::DIM>::order;
    

    int order_edge[ET_trait<ET>::N_EDGE];
    IVec<DIM-1> order_facet[ET_trait<ET>::N_FACET];
    IVec<DIM> order_inner;

    
  public:
    using VertexOrientedFE<ET>::SetVertexNumbers;
    
    T_HCurlCurlFE (int aorder)
    {
      order = aorder;
      for (auto & of : order_facet) of = aorder;
      order_inner = aorder;

    }
    
    virtual ELEMENT_TYPE ElementType() const override { return ET; }
    const HCurlCurlFE<ET> * Cast() const { return static_cast<const HCurlCurlFE<ET>*> (this); } 
    
    INLINE void SetOrderFacet (int nr, IVec<DIM-1,int> order) { order_facet[nr] = order; }
    INLINE void SetOrderEdge (int nr, int order) { order_edge[nr] = order; }
    INLINE void SetOrderInner (IVec<DIM,int> order) { order_inner = order; }

    virtual void ComputeNDof()
    {
      cout << "Error, T_HCurlCurlFE<ET>:: ComputeNDof not available, only for ET == TRIG" << endl;
    }


    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
                                  BareSliceMatrix<double> shapes) const override
    {
      Switch<4-DIM>
        (bmip.DimSpace()-DIM,[this, &bmip, shapes](auto CODIM)
         {
           constexpr auto DIMSPACE = DIM+CODIM.value;
           auto & mip = static_cast<const MappedIntegrationPoint<DIM, DIM+CODIM.value>&> (bmip);
           
           Cast() -> T_CalcShape (GetTIP(mip),SBLambda([shapes,DIMSPACE](int nr,auto val)
                                                   {
                                                     shapes.Row(nr).Range(DIMSPACE*DIMSPACE) = val.Value().AsVector();
                                                   }));
         });
    }


    
    virtual void EvaluateMappedShape (const BaseMappedIntegrationPoint & bmip,
                                      BareSliceVector<double> coefs,
                                      BareSliceMatrix<double> shape) const override
    {
      Switch<4-DIM>
        (bmip.DimSpace()-DIM,[this, &bmip, coefs, shape](auto CODIM)
      {
        auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);

        Mat<DIM+CODIM.value,DIM+CODIM.value> summat(0);
        Cast() -> T_CalcShape (GetTIP(mip), SBLambda ([&summat,coefs] (int nr, auto val)
                                                       {
                                                         summat += coefs(nr) * val.Value();

                                                       }));
        for (int k = 0; k < sqr(DIM+CODIM.value); k++)
          shape(k) = summat(k);
      });

    }

    virtual void CalcMappedIncShape (const BaseMappedIntegrationPoint & bmip,
                                  BareSliceMatrix<double> shapes) const override
    {
      auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM>&> (bmip);
      if constexpr (ET == ET_TET || ET == ET_TRIG || ET == ET_QUAD)
                     Cast() -> T_CalcShape (GetTIPHesse(mip),SBLambda([shapes](int nr,auto val)
                                                                       {
                                                                         if constexpr (DIM==3)
                                                                           shapes.Row(nr).Range(DIM_DDMAT) = val.Inc().AsVector();
                                                                         else
                                                                           shapes.Row(nr).Range(DIM_DDMAT) = val.Inc();
                                                                       }));
      else
        throw Exception("HCurlCurl::CalcMappedIncShape implemented only for TRIG and TET");
      
    }
                            
    virtual void EvaluateMappedIncShape (const BaseMappedIntegrationPoint & bmip,
                                         BareSliceVector<double> coefs,
                                         BareSliceVector<double> inc) const override
    {
      auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM>&> (bmip);

      Mat<DIM*(DIM-1)/2,DIM*(DIM-1)/2> sum = 0.0;
      if constexpr (ET == ET_TET || ET == ET_TRIG || ET == ET_QUAD)
                     Cast() -> T_CalcShape (GetTIPHesse(mip),SBLambda([coefs, &sum](int nr,auto val)
                                                                       {
                                                                         sum += coefs(nr) * Mat<DIM*(DIM-1)/2,DIM*(DIM-1)/2>(val.Inc());
                                                                       }));
      else
        throw Exception("HCurlCurl::EvaluateMappedIncShape implemented only for TRIG and TET");
      
      inc.Range(0,DIM_DDMAT) = sum.AsVector();
    }

    virtual void CalcMappedIncShape (const SIMD_BaseMappedIntegrationRule & bmir,
                                     BareSliceMatrix<SIMD<double>> shapes) const override
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
      for (size_t i = 0; i < mir.Size(); i++)
        {
          if constexpr (ET == ET_TET || ET == ET_TRIG || ET == ET_QUAD)
                     {
                       Cast() -> T_CalcShape (GetTIPHesse(mir[i]),SBLambda([shapes,i](int j,auto val)
                                                                         {
                                                                           if constexpr (DIM==3)
                                                                                          shapes.Rows(j*sqr(DIM), (j+1)*sqr(DIM)).Col(i).Range(0,DIM_DDMAT) = val.Inc().AsVector();
                                                                           else
                                                                             shapes.Rows(j,j+1).Col(i).Range(0,DIM_DDMAT) = val.Inc();
                                                                         }));
                     }
          else
            throw Exception("HCurlCurl::CalcMappedIncShape implemented only for TRIG and TET");

        }
    }

    void EvaluateIncShape (const SIMD_BaseMappedIntegrationRule & bmir,
                           BareSliceVector<> coefs,
                           BareSliceMatrix<SIMD<double>> values) const override
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
      for (size_t i = 0; i < bmir.Size(); i++)
        {
          double *pcoefs = &coefs(0);
          const size_t dist = coefs.Dist();  
          if constexpr (ET == ET_TET && DIM == 3)
                     {
                       Mat<DIM,DIM,SIMD<double>> summat(0);
                       Cast() -> T_CalcShape (GetTIPHesse(mir[i]),
                                               SBLambda ([&summat,&pcoefs,dist] (size_t j, auto val)
                                                         {
                                                           summat += (*pcoefs)*val.Inc();
                                                           pcoefs += dist;
                                                         }));
                       for (size_t k = 0; k < sqr(DIM); k++)
                         values(k,i) = summat(k);
                     }
          if constexpr ((ET == ET_TRIG || ET == ET_QUAD) && DIM == 2)
                     {
                       SIMD<double> summat(0);
                       Cast() -> T_CalcShape (GetTIPHesse(mir[i]),
                                               SBLambda ([&summat,&pcoefs,dist] (size_t j, auto val)
                                                         {
                                                           summat += (*pcoefs)*val.Inc();
                                                           pcoefs += dist;
                                                         }));
                       values(0,i) = summat;
                     }
        }
          
      /*if constexpr (ET == ET_TET || ET == ET_TRIG)
          {
            Switch<1>
              (bmir.DimSpace()-DIM,[values,&bmir,coefs,this](auto CODIM)
              {
                constexpr auto DIMSPACE = DIM+CODIM.value;
                auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
                for (size_t i = 0; i < bmir.Size(); i++)
                  {
                    double *pcoefs = &coefs(0);
                    const size_t dist = coefs.Dist();            
                  

                    if constexpr (DIMSPACE==3)
                       {
                         Mat<DIMSPACE,DIMSPACE,SIMD<double>> summat(0);
                         Cast() -> T_CalcShape (GetTIPHesse(mir[i]),
                                                 SBLambda ([&summat,&pcoefs,dist] (size_t j, auto val)
                                                           {
                                                             summat += (*pcoefs)*val.Inc();
                                                             pcoefs += dist;
                                                           }));
                         for (size_t k = 0; k < sqr(DIMSPACE); k++)
                           values(k,i) = summat(k);
                       }
                    if constexpr (DIMSPACE==2)
                      {
                        SIMD<double> summat(0);
                        Cast() -> T_CalcShape (GetTIPHesse(mir[i]),
                                                SBLambda ([&summat,&pcoefs,dist] (size_t j, auto val)
                                                          {
                                                            summat += (*pcoefs)*val.Inc();
                                                            pcoefs += dist;
                                                          }));
                        values(0,i) = summat;
                        }
                   
                  }
              });
              }*/
    }

    void AddTransIncShape (const SIMD_BaseMappedIntegrationRule & ir,
                                   BareSliceMatrix<SIMD<double>> values,
                                   BareSliceVector<> coefs) const override
    {
      throw ExceptionNOSIMD("HCurlCurl::AddTransIncShape not implemented yet");
    }


    
    virtual void CalcDualShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> shape) const override
    {
      shape.AddSize(ndof, sqr(bmip.DimSpace())) = 0.0;
      Switch<4-DIM>
        (bmip.DimSpace()-DIM,[this, &bmip, shape](auto CODIM)
         {
           auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);

           Cast() -> CalcDualShape2 (mip, SBLambda([shape] (size_t nr, auto val)
                                                   {
                                                     shape.Row(nr) = val.AsVector();
                                                   }));
         });
    }

    virtual void CalcDualShape (const SIMD_BaseMappedIntegrationRule& bmir, BareSliceMatrix<SIMD<double>> shapes) const override
    {
      Switch<4-DIM>
        (bmir.DimSpace()-DIM,[this, &bmir, shapes](auto CODIM)
         {
           constexpr int DIMSPACE = DIM+CODIM.value;
           auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM+CODIM.value>&> (bmir);

           shapes.AddSize(ndof*sqr(DIMSPACE), mir.Size()) = 0.0;
           for (size_t i = 0; i < mir.Size(); i++)
             {
               Cast() -> CalcDualShape2 (mir[i], SBLambda([shapes,i,DIMSPACE] (size_t j, auto val)
                                                          {
                                                            shapes.Rows(j*sqr(DIMSPACE), (j+1)*sqr(DIMSPACE)).Col(i).Range(0,sqr(DIMSPACE)) = val.AsVector();
                                                          }));
             }
         });
    }
    
    virtual void EvaluateDual (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const override
    {
      Switch<4-DIM>
        (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
         {
           constexpr int DIMSPACE = DIM+CODIM.value;
           auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM+CODIM.value>&> (bmir);
           for (size_t i = 0; i < mir.Size(); i++)
             {
               Mat<DIMSPACE,DIMSPACE,SIMD<double>> sum (SIMD<double>(0.0));
               Cast() -> CalcDualShape2 (mir[i], SBLambda([&sum, coefs] (size_t j, auto val)
                                                          {
                                                            sum += coefs(j) * val;
                                                          }));
               for (size_t k = 0; k < sqr(DIMSPACE); k++)
                 values(k, i) = sum(k);
             }
         });
    }

    virtual void AddDualTrans (const SIMD_BaseMappedIntegrationRule& bmir, BareSliceMatrix<SIMD<double>> values, BareSliceVector<double> coefs) const override
    {
      Switch<4-DIM>
        (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
         {
           constexpr int DIMSPACE = DIM+CODIM.value;
           auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM+CODIM.value>&> (bmir);
           for (size_t i = 0; i < mir.Size(); i++)
             {
               Mat<DIMSPACE,DIMSPACE,SIMD<double>> value;
               for (size_t k = 0; k < sqr(DIMSPACE); k++)
                 value(k) = values(k, i);
               
               Cast()-> CalcDualShape2 (mir[i], SBLambda([value, coefs] (size_t j, auto val)
                                                         {
                                                           coefs(j) += HSum(InnerProduct(val,value));
                                                         }));
             }
         });
    }

    virtual void CalcMappedCurlShape (const BaseMappedIntegrationPoint & bmip,
                                     BareSliceMatrix<double> shape) const override
    {
      auto mip = static_cast<const MappedIntegrationPoint<DIM,DIM> &>(bmip);

      if constexpr (ET == ET_TET || ET == ET_TRIG)
                     Cast() -> T_CalcShape (GetTIPHesse(mip),SBLambda([shape](int nr,auto val)
                                                                       {
                                                                         if constexpr (DIM==3)
                                                                           shape.Row(nr).Range(DIM_DMAT) = val.Curl().AsVector();
                                                                         else
                                                                           shape.Row(nr).Range(DIM_DMAT) = val.Curl();
                                                                       }));
      else
        throw Exception("HCurlCurl::CalcMappedCurlShape implemented only for TRIG and TET");
    }


    virtual void CalcMappedShape (const SIMD<BaseMappedIntegrationPoint> & bmip, 
                                         BareSliceMatrix<SIMD<double>> shape) const override
    {
      Switch<4-DIM>
        (bmip.DimSpace()-DIM,[this, &bmip, shape](auto CODIM)
         {
           constexpr auto DIMSPACE = DIM+CODIM.value;
           auto & mip = static_cast<const SIMD<MappedIntegrationPoint<DIM,DIM+CODIM.value>>&> (bmip);
            this->Cast() -> T_CalcShape (GetTIP(mip),
                                          SBLambda ([shape,DIMSPACE] (size_t j, auto val) 
                                                    {
                                                      shape.Rows(j*sqr(DIMSPACE), (j+1)*sqr(DIMSPACE)).Col(0).Range(0,sqr(DIMSPACE)) = val.Value().AsVector();
                                                    }));
         });
    }
      
    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & bmir, 
                                         BareSliceMatrix<SIMD<double>> shapes) const override
    {
      Switch<4-DIM>
        (bmir.DimSpace()-DIM,[this, &bmir, shapes](auto CODIM)
         {
           constexpr auto DIMSPACE = DIM+CODIM.value;
           auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM+CODIM.value>&> (bmir);
           for (size_t i = 0; i < mir.Size(); i++)
             {
               this->Cast() -> T_CalcShape (GetTIP(mir[i]),
                                             SBLambda ([i,shapes,DIMSPACE] (size_t j, auto val) 
                                                       {
                                                         shapes.Rows(j*sqr(DIMSPACE), (j+1)*sqr(DIMSPACE)).Col(i).Range(0,sqr(DIMSPACE)) = val.Value().AsVector();
                                                       }));
             }
         });
    }


    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & bmir,
                                  BareSliceVector<> coefs,
                                  BareSliceMatrix<SIMD<double>> values) const override
    {
      Switch<4-DIM>
        (bmir.DimSpace()-DIM,[values,&bmir, coefs,this](auto CODIM)
       {
         constexpr auto DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < bmir.Size(); i++)
           {
             double *pcoefs = &coefs(0);
             const size_t dist = coefs.Dist();

             Mat<DIMSPACE,DIMSPACE,SIMD<double>> summat(0);
             Cast() -> T_CalcShape (GetTIP(mir[i]),
                                    SBLambda ([&summat,&pcoefs,dist] (size_t j, auto val)
                                              {
                                                summat += (*pcoefs)*val.Value();
                                                pcoefs += dist;
                                              }));
             for (size_t k = 0; k < sqr(DIMSPACE); k++)
               values(k,i) = summat(k);              
           }
        });
    }

    virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & bmir,
                                  BareSliceMatrix<SIMD<double>> values,
                                  BareSliceVector<> coefs) const override
    {
      Switch<4-DIM>
      (bmir.DimSpace()-DIM,[this,&bmir,coefs,values](auto CODIM)
       {
         constexpr int DIMSPACE = DIM+CODIM.value;
         auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
         for (size_t i = 0; i < mir.Size(); i++)
           {
             Mat<DIMSPACE,DIMSPACE,SIMD<double>> vali;
             vali.AsVector()= values.Col(i);

             double *pcoefs = &coefs(0);
             const size_t dist = coefs.Dist();
             Cast()->T_CalcShape (GetTIP(mir[i]),
                                SBLambda ([vali,&pcoefs,dist] (size_t j, auto s)
                                          {
                                            *pcoefs += HSum(InnerProduct(s.Value(), vali));
                                             pcoefs += dist;
                                          }));
           }
       });
      
    }
   
  };



  
#ifdef FILE_HCURLCURLFE_CPP
#define HCURLCURLFE_EXTERN
#else
#define HCURLCURLFE_EXTERN extern
#endif
  
  HCURLCURLFE_EXTERN template class HCurlCurlFiniteElement<2>;
  HCURLCURLFE_EXTERN template class HCurlCurlFiniteElement<3>;

  template <> class HCurlCurlFE<ET_SEGM> : public T_HCurlCurlFE<ET_SEGM> 
  {
    
  public:
    using T_HCurlCurlFE<ET_SEGM> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      ndof += order_inner[0]+1;
      order = max2(order,order_inner[0]);

    }

    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<1,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x;
      Tx lami[2] ={ x, 1-x };
      int ii = 0;

      IVec<2> e = ET_trait<ET_SEGM>::GetEdgeSort (0, vnums);
      Tx ls = lami[e[0]], le = lami[e[1]];

      auto symdyadic = MakeReggeAD(ls, le);

      LegendrePolynomial::Eval(order_inner[0], ls-le, SBLambda([symdyadic, &ii, shape] (size_t nr, auto val)
                            {
                              shape[ii++] = 0.5*val*symdyadic;
                            }));
    }

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      auto & ip = mip.IP();
      typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
      T x = ip(0);
      T lam[2] = { x, 1-x };

      int ii = 0;

      IVec<2> e = ET_trait<ET_SEGM>::GetEdgeSort (0, vnums);
      T xi = lam[e[0]]-lam[e[1]];

      auto tv = mip.GetJacobian()*Vec<1,T>(1);
      auto tt = DyadProd(tv,tv);

      LegendrePolynomial::Eval(order_inner[0], xi, SBLambda([shape,mip,tt,&ii] (size_t nr, T val)
                            {
                              shape[ii++] = 1/mip.GetMeasure()*val*tt;
                            }));

    }
  };

  
  template <> class HCurlCurlFE<ET_TRIG> : public T_HCurlCurlFE<ET_TRIG> 
  {
    
  public:
    using T_HCurlCurlFE<ET_TRIG> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<3; i++)
      {
        ndof += order_facet[i][0]+1;
        order = max2(order, order_facet[i][0]);
      }
      int ninner = 3*order_inner[0]*(order_inner[0]+1)/2 ;
      order = max2(order, order_inner[0]);
      
      ndof += ninner;

    }
    

    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y;
      Tx lami[3] ={ x, y, 1-x-y };
      int ii = 0;

      //    /*int maxorder_facet =
      //      max2(order_facet[0][0],max2(order_facet[1][0],order_facet[2][0]));
      //    ArrayMem<Tx,20> ha(maxorder_facet+1);
      //    ArrayMem<Tx,20> u(order_inner[0]+2), v(order_inner[0]+2);
      
      //    for (int i = 0; i < 3; i++)
      //      {
      //        IVec<2> e = ET_trait<ET_TRIG>::GetEdgeSort(i,vnums);
      //        Tx ls = llami[e[0]], le = llami[e[1]];
      
      //        // edge functions are all curl-free!
      //        IntegratedLegendreMonomialExt::CalcTrigExt(maxorder_facet+2,
      //                                                   le-ls, 1-le-ls, ha);
      
      //        for (int l = 0; l <= order_facet[i][0]; l++)
      //          shape[ii++] = EpsGrad (ha[l]);
      //          }*/
      
      for (int i = 0; i < 3; i++)
        {
          IVec<2> e = ET_trait<ET_TRIG>::GetEdgeSort (i, vnums);
          Tx ls = lami[e[1]], le = lami[e[0]];

          auto symdyadic = MakeReggeAD(ls, le);

          LegendrePolynomial::EvalScaled(order_facet[i][0], ls-le,ls+le, SBLambda([symdyadic, &ii, shape] (size_t nr, auto val)
                            {
                              shape[ii++] = -val*symdyadic;
                            }));
        }


      if (order_inner[0] > 0)
        {
	  IVec<4> f = ET_trait<ET_TRIG>::GetFaceSort(0, vnums); 
	  Tx ls = lami[f[0]], le = lami[f[1]], lt = lami[f[2]];
	  
          auto symdyadic1 = lt*MakeReggeAD(ls, le);
          auto symdyadic2 = ls*MakeReggeAD(lt, le);
          auto symdyadic3 = le*MakeReggeAD(ls, lt);

          
	  DubinerBasis::Eval(order_inner[0]-1, ls,le,
			      SBLambda([symdyadic1,symdyadic2,symdyadic3, &ii, shape] (size_t nr, auto val)
				       {
					 shape[ii++] = 2*val*symdyadic1;
					 shape[ii++] = 2*val*symdyadic2;
                                         shape[ii++] = 2*val*symdyadic3;
				       }));
	}
      
    };

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      auto & ip = mip.IP();
      typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
      T x = ip(0), y = ip(1);
      T lam[3] = { x, y, 1-x-y };
      Vec<2,T> pnts[3] = { { 1, 0 }, { 0, 1 } , { 0, 0 } };
      int facetnr = ip.FacetNr();

      int ii = 0;


      if (ip.VB() == BND)
        { // facet shapes
          for (int i = 0; i < 3; i++)
            {
              int p = order_facet[i][0];
              
              if (i == facetnr)
                {             
                  IVec<2> e = ET_trait<ET_TRIG>::GetEdgeSort (i, vnums);
                  
                  T xi = lam[e[0]]-lam[e[1]];
                  Vec<2,T> tauref = pnts[e[0]] - pnts[e[1]];
                  
                  
                  auto tv = mip.GetJacobian()*tauref;

                  auto tt = DyadProd(tv,tv);
                  LegendrePolynomial::Eval
                    (p, xi,
                     SBLambda([&] (size_t nr, T val)
                              {
                                shape[nr+ii] = 1/mip.GetMeasure()*val*tt;
                              }));
                }
              ii += (p+1);
            }
        }
      else
        {
          for (int i = 0; i < 3; i++)
            ii += order_facet[i][0]+1;
        }
      if (ip.VB() == VOL)
        {
          auto p = order_inner[0]-1;
          if( p >= 0 )
            {
              IVec<4> f =  ET_trait<ET_TRIG>::GetFaceSort(0, vnums);

              DubinerBasis::Eval (p, lam[f[0]], lam[f[1]],
                                   SBLambda([&] (size_t nr, T val)
                                            {
                                              shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<2,2>({{1,0},{0,0}})*Trans(mip.GetJacobian());
                                              shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<2,2>({{0,0},{0,1}})*Trans(mip.GetJacobian());
                                              shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<2,2>({{0,1},{1,0}})*Trans(mip.GetJacobian());
                                            }));
            }
        }
    }

    
  };
  
  template <> class HCurlCurlFE<ET_QUAD> : public T_HCurlCurlFE<ET_QUAD> 
  {
    
  public:
    using T_HCurlCurlFE<ET_QUAD> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i<4; i++)
      {
        ndof += order_facet[i][0]+1;
        order = max2(order, order_facet[i][0]);
      }
      int ninner = order_inner[0]*order_inner[0] + (order_inner[0]+2)*order_inner[0]*2  +1;//+ 2*order_inner[0];
      order = max2(order, order_inner[0]);
      order += 1;
      ndof += ninner;

    }



    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y;
      Tx lx[4] ={ 1-x, x, x, 1-x };
      Tx ly[4] ={ 1-y, 1-y, y, y };
      Tx lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
      Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
      int ii = 0;

      ArrayMem<Tx,20> v(order+2), u(order+2);
      
     
      for (int i = 0; i < 4; i++)
        {
          IVec<2> e = ET_trait<ET_QUAD>::GetEdgeSort (i, vnums);
          Tx xi  = sigma[e[1]]-sigma[e[0]];
          Tx lam_e = lami[e[0]]+lami[e[1]];  
          auto symdyadic = MakeReggeAD(xi, xi);


          //IntLegNoBubble::
            LegendrePolynomial::
            EvalMult (order_facet[i][0], 
                      xi, 0.25*lam_e, SBLambda ([&](int i, auto val)
                                           {
                                             shape[ii++] = val*symdyadic;
                                           }));
        }



      int oi = order_inner[0];

      auto symdyadic = MakeReggeAD(0.5*x,0.5*y); //(0,0.5,  0.5,0) * P(y) * P(x)
      
      Tx eta = ly[2]-ly[1];
      Tx xi = lx[1]-lx[0];
      LegendrePolynomial (oi, eta, v);
      LegendrePolynomial (oi, xi, u);

      for (int i = 0; i <= oi; i++)
        for (int j = 0; j <= oi; j++)
          {
            shape[ii++] = u[i]*v[j]*symdyadic;
          }

      
      auto symdyad = lx[1]*lx[0]*MakeReggeAD(y,y);//x*(1-x)*(0,0,  0,1) * P(y) * P(x)
      for (int i = 0; i < oi; i++)
        for (int j = 0; j <= oi; j++)
          {
            shape[ii++] = u[i]*v[j]*symdyad;
          }

      symdyad = ly[2]*ly[1]*MakeReggeAD(x,x); //y*(1-y)*(1,0,  0,0) * P(x) * P(y)
      
      for (int j = 0; j < oi; j++)
        for (int i = 0; i <= oi; i++)
          {
            shape[ii++] = u[i]*v[j]*symdyad;
          }

      //old version
      //ArrayMem<Tx,20> u(order+2);
      /*for (int i = 0; i < 4; i++)
        {
          IVec<2> e = ET_trait<ET_QUAD>::GetEdgeSort (i, vnums);
          Tx xi = llx[e[1]]+lly[e[1]]-llx[e[0]]-lly[e[0]];
          Tx eta = llx[e[0]]*lly[e[0]]+llx[e[1]]*lly[e[1]];

	  IntegratedLegendreMonomialExt::Calc(order_facet[i][0]+2,xi,u);

          
          for (int l = 0; l <= order_facet[i][0]; l++)
            shape[ii++] = Eps_u_Gradv (eta, u[l]);
        }
      
      IntegratedLegendreMonomialExt::Calc(oi+3,llx[0]-llx[1],u);
      IntegratedLegendreMonomialExt::Calc(oi+3,lly[0]-lly[2],v);

      for(int i = 0; i <= oi-1; i++)
        for(int j = 0; j <= oi-1; j++)
          shape[ii++] = EpsGrad(u[i]*v[j]);

      for(int i = 0; i <= oi+1; i++)
        for(int j = 0; j <= oi-1; j++)
        {
          shape[ii++] = vEpsGradu(u[i],v[j]);
          shape[ii++] = vEpsGradu(v[i],u[j]);
        }
      shape[ii++] = Eps_u_Gradv(lx[0], ly[0]);

      for(int i = 0; i <= oi-1; i++)
      {
        shape[ii++] = Eps_u_Gradv(u[i], ly[0]);
        shape[ii++] = Eps_u_Gradv(v[i], lx[0]);
        }*/

      
      
    };

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      auto & ip = mip.IP();
      typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T; 

      T x = ip(0), y = ip(1);
      T lx[4] = { 1-x, x, x, 1-x };
      T ly[4] = { 1-y, 1-y, y, y };
      T sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};
      
      Vec<2,T> pnts[4] = {  { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
      int facetnr = ip.FacetNr();

      int ii = 0;
      
      ArrayMem<T,20> v(order+2), u(order+2);

      
      if (mip.IP().VB() == BND)
        { // facet shapes
          for (int i = 0; i < 4; i++)
            {
              int p = order_facet[i][0];
              
              if (i == facetnr)
                {             
                  IVec<2> e = ET_trait<ET_QUAD>::GetEdgeSort (i, vnums);
                  
                  //T xi = lam[e[0]]-lam[e[1]];
                  T xi  = sigma[e[1]]-sigma[e[0]];
                  Vec<2,T> tauref = pnts[e[0]] - pnts[e[1]];
                  
                  
                  auto tv = mip.GetJacobian()*tauref;

                  auto tt = DyadProd(tv,tv);
                  LegendrePolynomial::Eval
                    (p, xi,
                     SBLambda([&] (size_t nr, T val)
                              {
                                shape[nr+ii] = 1/mip.GetMeasure()*val*tt;
                                }));
                  /*IVec<2> e = ET_trait<ET_QUAD>::GetEdgeSort (i, vnums);
                  AutoDiff<2,T> xi  = sigma[e[1]]-sigma[e[0]];
                  AutoDiff<2,T> lam_e = lami[e[0]]+lami[e[1]];  
                  Vec<3, AutoDiff<2,T>> symdyadic = SymDyadProd(xi,xi);


                  IntLegNoBubble::
                    EvalMult (order_edge[i], 
                              xi, lam_e, SBLambda ([&](int nr, auto val)
                                                   {
                                                     VecToSymMat<2>(T_REGGE_Shape<2,T>(val*symdyadic).Shape(),tmp);
                                                     shape[nr + ii] = mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
                                                     }));*/
                  /*AutoDiff<2,T> xi  = sigma[e[1]]-sigma[e[0]];
                  AutoDiff<2,T> lam_e = lami[e[0]]+lami[e[1]];  
                  Vec<3, AutoDiff<2,T>> symdyadic = SymDyadProd(xi,xi);


                  IntLegNoBubble::
                    EvalMult (p,xi, lam_e, SBLambda ([&](int nr, auto val)
                                                   {
                                                     VecToSymMat<2>(T_REGGE_Shape<2,T>(val*symdyadic).Shape(),tmp);
                                                     shape[nr + ii] = 1/mip.GetMeasure()*tmp;
                                                     }));*/
                }
              ii += (p+1);
            }
        }
      else
        {
          for (int i = 0; i < 4; i++)
            ii += order_facet[i][0]+1;
        }
      
      if (mip.IP().VB() == VOL)
        {
          auto p = order_inner[0];
         
          T eta = ly[2]-ly[1];
          T xi = lx[1]-lx[0];
          LegendrePolynomial (p, eta, v);
          LegendrePolynomial (p, xi, u);
          
          for (int i = 0; i <= p; i++)
            for (int j = 0; j <= p; j++)
              {
                shape[ii++] = 1/mip.GetMeasure()*u[i]*v[j]*mip.GetJacobian()*Mat<2,2>(Matrix<>({{0,1},{1,0}}))*Trans(mip.GetJacobian());
              }
          
          
          //auto symdyad = lx[1]*lx[0]*SymDyadProd(Vec<2,T>(0,1),Vec<2,T>(0,1));//x*(1-x)*(0,0,  0,1) * P(y) * P(x)
          for (int i = 0; i < p; i++)
            for (int j = 0; j <= p; j++)
              {
                shape[ii++] = 1/mip.GetMeasure()*u[i]*v[j]*mip.GetJacobian()*Mat<2,2>(Matrix<>({{0,0},{0,1}}))*Trans(mip.GetJacobian());
              }
          
          //symdyad = ly[2]*ly[1]*SymDyadProd(Vec<2,T>(1,0),Vec<2,T>(1,0)); //y*(1-y)*(1,0,  0,0) * P(x) * P(y)
          
          for (int j = 0; j < p; j++)
            for (int i = 0; i <= p; i++)
              {
                shape[ii++] = 1/mip.GetMeasure()*u[i]*v[j]*mip.GetJacobian()*Mat<2,2>(Matrix<>({{1,0},{0,0}}))*Trans(mip.GetJacobian());
              }
      
          //IVec<4> f = ET_trait<ET_QUAD>::GetFaceSort(0, vnums);
          
          /*IntegratedLegendreMonomialExt::Calc(p+3,lx[0]-lx[1],u);
          IntegratedLegendreMonomialExt::Calc(p+3,ly[0]-ly[2],v);

          Mat<2,2,T> tmp;
          
          for(int i = 0; i <= p-1; i++)
              for(int j = 0; j <= p-1; j++)
                {
                  VecToSymMat<2>(EpsGrad(u[i]*v[j]).Shape(),tmp);
                  shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
                }
          
          for(int i = 0; i <= p+1; i++)
              for(int j = 0; j <= p-1; j++)
                {
                  VecToSymMat<2>(vEpsGradu(u[i],v[j]).Shape(),tmp);
                  shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
                  VecToSymMat<2>(vEpsGradu(v[i],u[j]).Shape(),tmp);
                  shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
                }

          VecToSymMat<2>(Eps_u_Gradv(lx[0], ly[0]).Shape(),tmp);
          shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
          
          for(int i = 0; i <= p-1; i++)
            {
              VecToSymMat<2>(Eps_u_Gradv(u[i], ly[0]).Shape(),tmp);
              shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
              VecToSymMat<2>(Eps_u_Gradv(v[i], lx[0]).Shape(),tmp);
              shape[ii++] = 1/mip.GetMeasure()*mip.GetJacobian()*tmp*Trans(mip.GetJacobian());
              }*/
        }
    }
  };


  template <> class HCurlCurlFE<ET_PRISM> : public T_HCurlCurlFE<ET_PRISM> 
  {
  public:
    enum { incrorder_xx1 = 0};
    enum { incrorder_zz1 = 0};
    enum { incrorder_xx2 = 0};
    enum { incrorder_zz2 = 0};
    enum { incrorder_xx1_bd = 0};
    enum { incrorder_zz1_bd = 0};
    enum { incrorder_xx2_bd = 0};
    enum { incrorder_zz2_bd = 0};
    using T_HCurlCurlFE<ET_PRISM> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;

      for (int i=0; i < 9; i++)
        {
          ndof += order_edge[i]+1;
          order = max2(order,order_edge[i]);
        }

      for (int i=0; i<2; i++)
        {
          ndof += 3*(order_facet[i][0])*(order_facet[i][0]+1)/2;
          order = max2(order, order_facet[i][0]);
        }

      for (int i=2; i<5; i++)
        {
          ndof += order_facet[i][0]*order_facet[i][0] + (order_facet[i][0]+2)*order_facet[i][0]*2 +1;
          order = max2(order, order_facet[i][0]);
        }
      int p = order_inner[0];
      int ninner =  3*p*(p+1)/2*p + (p-1)*(p)/2*(p+1) + (p+1)*p*(p+1);
      ndof += ninner;

      order = 1+max2(order, p);
    }
    

    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y, z = ip.z;
      Tx lx[6] ={ x, y, 1-x-y, x, y, 1-x-y };
      Tx lz[6] ={ 1-z,1-z,1-z,z,z,z };

      int ii = 0;
      

      const FACE * faces = ElementTopology::GetFaces(ET_PRISM);

      ArrayMem<Tx,20> leg_u(order+2), leg_v(order+3);
      ArrayMem<Tx,20> leg_w(order+2);

      //horizontal edge shapes
      for (int i = 0; i < 6; i++)
        {
          IVec<2> e = ET_trait<ET_PRISM>::GetEdgeSort (i, vnums);
          Tx ls = lx[e[1]], le = lx[e[0]], lm = lz[e[0]];

          auto symdyadic = lm*MakeReggeAD(ls,le);

          LegendrePolynomial::EvalScaled(order_edge[i], ls-le,ls+le, SBLambda([symdyadic, &ii, shape] (size_t nr, auto val)
                            {
                              shape[ii++] = -val*symdyadic;
                            }));
        }



      //vertical edge shapes
      for (int i = 6; i < 9; i++)
        {
          IVec<2> e = ET_trait<ET_PRISM>::GetEdgeSort (i, vnums);
          Tx ls = lx[e[0]], lm1 = lz[e[0]], lm2 = lz[e[1]];
          auto symdyadic = ls*MakeReggeAD(lm1,lm1);
          LegendrePolynomial (order_edge[i],lm1-lm2, leg_v);

          for (int j=0; j <= order_edge[i]; j++)
            shape[ii++] = leg_v[j]*symdyadic;
        }
      


      //horizontal face shaps
      for(int fa = 0; fa < 2; fa++)
        {
          if (order_facet[fa][0] > 0)
            {
              IVec<4> f = ET_trait<ET_PRISM>::GetFaceSort(fa, vnums);
              Tx ls = lx[f[0]], le = lx[f[1]], lt = lx[f[2]], lm = lz[f[0]];
              
              auto symdyadic1 = lm*lt*MakeReggeAD(ls,le);
              auto symdyadic2 = lm*ls*MakeReggeAD(lt,le);
              auto symdyadic3 = lm*le*MakeReggeAD(ls,lt);
              
              DubinerBasis::Eval(order_facet[fa][0]-1, ls,le,
                                 SBLambda([symdyadic1,symdyadic2,symdyadic3, &ii, shape] (size_t nr, auto val)
                                          {
                                            shape[ii++] = val*symdyadic1;
                                            shape[ii++] = val*symdyadic2;
                                            shape[ii++] = val*symdyadic3;
                                          }));
            }
        }


      //vertical face shaps
      for(int fa = 2; fa < 5; fa++)
        {
          int of = order_facet[fa][0];
          
          int fmax = 0;
          for(int j = 1; j < 4; j++)
            if(vnums[faces[fa][j]] > vnums[faces[fa][fmax]]) fmax = j;
          int fz,ftrig;
          fz = 3 - fmax;
          ftrig = fmax^1;          
          fmax = faces[fa][fmax];
          fz = faces[fa][fz];
          ftrig = faces[fa][ftrig];
          
          Tx eta = lz[fz]-lz[fmax];
          Tx xi = lx[ftrig]-lx[fmax];

          LegendrePolynomial (of, eta, leg_v);
          LegendrePolynomial (of, xi, leg_u);

          auto W = uDv_minus_vDu(lx[ftrig],lx[fmax]);
          Tx W_AD;
          W_AD.DValue(0) = W.Value()(0);
          W_AD.DValue(1) = W.Value()(1);
          W_AD.DValue(2) = W.Value()(2);
          auto symdyadic = MakeReggeAD(eta,0.25*W_AD);   //^= (0,1, 1,0) * P(x)*P(y)
          for (int j = 0; j <= of; j++)
            for (int k = 0; k <= of; k++)
              shape[ii++] = leg_v[j]*leg_u[k]*symdyadic;

          
          auto symdyad = 0.25*lx[ftrig]*lx[fmax]*MakeReggeAD(eta,eta);  //^= x*(1-x)*(0,0, 0,1) * P(x) * P(y)
          for (int i = 0; i < of; i++)
            for (int j = 0; j <= of; j++)
                shape[ii++] = leg_u[i]*leg_v[j]*symdyad;

          symdyad = 0.25*lz[fz]*lz[fmax]*MakeReggeAD(lx[ftrig],lx[fmax]);    //^= y*(1-y)*(1,0, 0,0) * P(x)*P(y)
          for (int j = 0; j < of; j++)
            for (int i = 0; i <= of; i++)
                shape[ii++] = leg_u[i]*leg_v[j]*symdyad;
        }
      
      //inner shapes
      int p = order_inner[0];
      if (p > 0)
        {
          
          IVec<4> f = ET_trait<ET_PRISM>::GetFaceSort(0, vnums);

          Tx ls = lx[f[0]], le = lx[f[1]], lt = lx[f[2]], lm = lz[0], ln = lz[3];

          auto symdyadic1 = lm*ln*lt*MakeReggeAD(ls,le);
          auto symdyadic2 = lm*ln*ls*MakeReggeAD(lt,le);
          auto symdyadic3 = lm*ln*le*MakeReggeAD(ls,lt);

          Tx eta = lz[0]-lz[4];
          LegendrePolynomial (p, eta, leg_w);

          // Reg(T) x [0,1]
          DubinerBasis::Eval(p-1, ls,le,
                             SBLambda([symdyadic1,symdyadic2,symdyadic3, &ii, shape,p,leg_w] (size_t nr, auto val)
                                      {
                                        for(int j=0; j < p; j++)
                                          {
                                            shape[ii++] = leg_w[j]*val*symdyadic1;
                                            shape[ii++] = leg_w[j]*val*symdyadic2;
                                            shape[ii++] = leg_w[j]*val*symdyadic3;
                                          }
                                      }));
          

          // H1(T) x [0,1]
          auto symdyadic = ls*le*lt*MakeReggeAD(eta,eta);
          DubinerBasis::Eval(p-2, ls,le,
                             SBLambda([symdyadic, &ii, shape,p,leg_w] (size_t nr, auto val)
                                      {
                                        for(int j=0; j <= p; j++)
                                          {
                                            shape[ii++] = val*leg_w[j]*symdyadic;
                                          }
                                      }));


          // Nedelec_1 x [0,1]
          DubinerBasis::EvalMult(p-2, lx[f[0]], lx[f[1]],lx[f[0]]*lx[f[1]]*lx[f[2]], 
                                 SBLambda([&](int nr, auto val)
                                          {
                                            auto tmp = Du(val);
                                            Tx tmp_AD;
                                            tmp_AD.DValue(0) = tmp.Value()(0);
                                            tmp_AD.DValue(1) = tmp.Value()(1);
                                            tmp_AD.DValue(2) = tmp.Value()(2);
                                            auto symdyadic = MakeReggeAD(tmp_AD,eta);
                                            for(int j=0; j <= p; j++)
                                              shape[ii++] = leg_w[j]*symdyadic;
                                          }));
          
          DubinerBasis::EvalMult(p-1, lx[f[0]], lx[f[1]], lx[f[0]], 
                                 SBLambda([&ii,shape,p,leg_w,eta,f,lx](int nr, auto val)
                                          {
                                            auto tmp = wuDv_minus_wvDu (lx[f[1]], lx[f[2]], val);
                                            Tx tmp_AD;
                                            tmp_AD.DValue(0) = tmp.Value()(0);
                                            tmp_AD.DValue(1) = tmp.Value()(1);
                                            tmp_AD.DValue(2) = tmp.Value()(2);
                                            auto symdyadic = MakeReggeAD(tmp_AD,eta);
                                            for(int j=0; j <= p; j++)
                                              shape[ii++] = leg_w[j]*symdyadic;
                                          }));
          
          LegendrePolynomial::EvalMult(p-1, lx[f[2]]-lx[f[1]], lx[f[2]], 
                                       SBLambda([&ii,shape,p,leg_w,eta,lx,f] (int j, auto val)
                                                {
                                                  auto tmp = wuDv_minus_wvDu (lx[f[1]], lx[f[0]], val);
                                                  Tx tmp_AD;
                                                  tmp_AD.DValue(0) = tmp.Value()(0);
                                                  tmp_AD.DValue(1) = tmp.Value()(1);
                                                  tmp_AD.DValue(2) = tmp.Value()(2);
                                                  auto symdyadic = MakeReggeAD(tmp_AD,eta);
                                                  for(int j=0; j <= p; j++)
                                                    shape[ii++] = leg_w[j]*symdyadic;
                                                }));
          
        }

    }


    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hcurlcurlfe calcdualshape2 not implementend for element type ET_PRISM");
    }

  };


  
  template <> class HCurlCurlFE<ET_TET> : public T_HCurlCurlFE<ET_TET> 
  {
  public:
    using T_HCurlCurlFE<ET_TET> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;

      for (int i=0; i<6; i++)
      {
        ndof += order_edge[i]+1;
        order = max2(order, order_edge[i]);
      }
      
      for (int i=0; i<4; i++)
      {
        ndof += 3*(order_facet[i][0])*(order_facet[i][0]+1)/2;
        order = max2(order, order_facet[i][0]);
      }
     
      int p = order_inner[0];
      int ninner = p > 1 ? 6*(p+1)*(p)*(p-1)/6 : 0;
      ndof += ninner; 

      order = max2(order, p);
    }



    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y, z = ip.z;
      Tx lam[4] = {x, y, z, 1-x-y-z};
      int ii = 0;

      for (int i = 0; i < 6; i++)
        {
          IVec<2> e = ET_trait<ET_TET>::GetEdgeSort (i, vnums);
          Tx ls = lam[e[1]], le = lam[e[0]];

          auto symdyadic = MakeReggeAD(ls, le);
          LegendrePolynomial::EvalScaled(order_edge[i], ls-le,ls+le, SBLambda([symdyadic, &ii, shape] (size_t nr, auto val)
                            {
                              shape[ii++] = -val*symdyadic;
                            }));
        }

      
      for(int fa = 0; fa < 4; fa++)
        {
          if (order_facet[fa][0] > 0)
            {
              IVec<4> f = ET_trait<ET_TET>::GetFaceSort(fa, vnums);
              Tx ls = lam[f[0]], le = lam[f[1]], lt = lam[f[2]];
              
              auto symdyadic1 = lt*MakeReggeAD(ls, le);
              auto symdyadic2 = ls*MakeReggeAD(lt, le);
              auto symdyadic3 = le*MakeReggeAD(ls, lt);
              
              DubinerBasis::Eval(order_facet[fa][0]-1, ls,le,
                                 SBLambda([symdyadic1,symdyadic2,symdyadic3, &ii, shape] (size_t nr, auto val)
                                          {
                                            shape[ii++] = val*symdyadic1;
                                            shape[ii++] = val*symdyadic2;
                                            shape[ii++] = val*symdyadic3;
                                          }));
            }
        }

      if (order_inner[0] > 1)
        {
          Tx li = lam[0], lj = lam[1], lk = lam[2], ll = lam[3];

          auto symdyadic1 = li*lj*MakeReggeAD(lk, ll);
          auto symdyadic2 = lj*lk*MakeReggeAD(ll, li);
          auto symdyadic3 = lk*ll*MakeReggeAD(li, lj);
          auto symdyadic4 = ll*li*MakeReggeAD(lj, lk);
          auto symdyadic5 = li*lk*MakeReggeAD(lj, ll);
          auto symdyadic6 = lj*ll*MakeReggeAD(li, lk);

          
          DubinerBasis3D::Eval (order_inner[0]-2, lam[0], lam[1], lam[2], SBLambda([&ii, shape, symdyadic1, symdyadic2, symdyadic3, symdyadic4, symdyadic5, symdyadic6](size_t j, auto val)
                                                                                   {
                                                                                     shape[ii++] = val*symdyadic1;
                                                                                     shape[ii++] = val*symdyadic2;
                                                                                     shape[ii++] = val*symdyadic3;
                                                                                     shape[ii++] = val*symdyadic4;
                                                                                     shape[ii++] = val*symdyadic5;
                                                                                     shape[ii++] = val*symdyadic6;
                                                                                   }));
        }
      
    }

    

    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      auto & ip = mip.IP();
      typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T;    
      T x = ip(0), y = ip(1), z = ip(2);
      T lam[4] = { x, y, z, 1-x-y-z };
      Vec<3,T> pnts[4] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } , { 0, 0, 0 } };
      int facetnr = ip.FacetNr();

      int ii = 0;

      if (ip.VB() == BBND)
        { // facet shapes
          for (int i = 0; i < 6; i++)
            {
              int p = order_edge[i];
              
              if (i == facetnr)
                {             
                  IVec<2> e = ET_trait<ET_TET>::GetEdgeSort (i, vnums);
                  
                  T xi = lam[e[1]]-lam[e[0]];
                  Vec<3,T> tauref = pnts[e[1]] - pnts[e[0]];
                  Vec<3,T> tau = mip.GetJacobian()*tauref;
                  Mat<3,3,T> tt = DyadProd(tau,tau);
                  LegendrePolynomial::Eval
                    (p, xi,
                     SBLambda([&] (size_t nr, T val)
                              {
                                shape[nr+ii] = 1/mip.GetMeasure()*val*tt;
                              }));
                }
              ii += (p+1);
            }
        }
      else
        {
          for (int i = 0; i < 6; i++)
            ii += order_edge[i]+1;
        }
      if (ip.VB() == BND)
        {
          for (int i = 0; i < 4; i++)
            {
              auto p = order_facet[i][0]-1;
              if( p >= 0 && i == facetnr )
                {
                  IVec<4> fav = ET_trait<ET_TET>:: GetFaceSort(facetnr, vnums);
                  Vec<3,T> adxi = pnts[fav[0]] - pnts[fav[2]];
                  Vec<3,T> adeta = pnts[fav[1]] - pnts[fav[2]];
                  T xi = lam[fav[0]];
                  T eta = lam[fav[1]];
                  
                  Matrix<T> F(3,2);
                  F.Col(0) = adxi;
                  F.Col(1) = adeta;
		 
                  Matrix<T> Ftmp(2,2);
                  Ftmp = Trans(F)*F;
                  auto det = sqrt(Ftmp(0,0)*Ftmp(1,1)-Ftmp(1,0)*Ftmp(0,1));
                                              
                  DubinerBasis::Eval (p, xi, eta,
                                       SBLambda([&] (size_t nr, T val)
                                                {
                                                  shape[ii++] = 1/(det*mip.GetMeasure())*val*Mat<3,3,T>(mip.GetJacobian()*F*Matrix<>({{1,0},{0,0}})*Trans(mip.GetJacobian()*F));
                                                  shape[ii++] = 1/(det*mip.GetMeasure())*val*Mat<3,3,T>(mip.GetJacobian()*F*Matrix<>({{0,0},{0,1}})*Trans(mip.GetJacobian()*F));
                                                  shape[ii++] = 1/(det*mip.GetMeasure())*val*Mat<3,3,T>(mip.GetJacobian()*F*Matrix<>({{0,1},{1,0}})*Trans(mip.GetJacobian()*F));
                                                }));
                }
              else
                ii += 3*(order_facet[i][0])*(order_facet[i][0]+1)/2;
            }
        }
      else
        {
          for (int i = 0; i < 4; i++)
            ii += 3*(order_facet[i][0])*(order_facet[i][0]+1)/2;
        }
      
      if (ip.VB() == VOL && order_inner[0] >= 2)
        {
          DubinerBasis3D::Eval (order_inner[0]-2, lam[0], lam[1], lam[2], SBLambda([&](size_t j, T val)
                               {
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{1,0,0},{0,0,0},{0,0,0}}))*Trans(mip.GetJacobian());
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{0,0,0},{0,1,0},{0,0,0}}))*Trans(mip.GetJacobian());
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{0,0,0},{0,0,0},{0,0,1}}))*Trans(mip.GetJacobian());
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{0,0,0},{0,0,1},{0,1,0}}))*Trans(mip.GetJacobian());
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{0,0,1},{0,0,0},{1,0,0}}))*Trans(mip.GetJacobian());
                                 shape[ii++] = 1/mip.GetMeasure()*val*mip.GetJacobian()*Mat<3,3>(Matrix<>({{0,1,0},{1,0,0},{0,0,0}}))*Trans(mip.GetJacobian());
                               }));
          
        }
    }
  };
  


  template <> class HCurlCurlFE<ET_HEX> : public T_HCurlCurlFE<ET_HEX> 
  {
  public:
    using T_HCurlCurlFE<ET_HEX> :: T_HCurlCurlFE;

    virtual void ComputeNDof()
    {
      order = 0;
      ndof = 0;
      for (int i=0; i < 12; i++)
        {
          ndof += order_edge[i]+1;
          order = max2(order,order_edge[i]);
        }
      for (int i=0; i<6; i++)
      {
        ndof += order_facet[i][0]*order_facet[i][0] + 2*(order_facet[i][0]+2)*order_facet[i][0]+1;
        order = max2(order, order_facet[i][0]);
      }
      int p = order_inner[0];
      ndof += 3*(p*(p+1)*(p+1) + p*p*(p+1) );

      order = 1 + max2(order, p);
    }


    template <typename Tx, typename TFA> 
    void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y, z = ip.z;
      Tx lx[2] ={ 1-x, x};
      Tx ly[2] ={ 1-y, y};
      Tx lz[2] ={ 1-z, z};
      Tx lami[8]={(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
                  (1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
      Tx sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
                   (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z};
      int ii = 0;
      
      const FACE * faces = ElementTopology::GetFaces(ET_HEX);

      ArrayMem<Tx,20> leg_u(order+2), leg_v(order+2), leg_w(order+2);
      
      // edges
      for (int i = 0; i < 12; i++)
        {
          int p = order_edge[i]; 
          IVec<2> e = ET_trait<ET_HEX>::GetEdgeSort (i, vnums);
          Tx xi  = sigma[e[1]]-sigma[e[0]];
          Tx lam_e = lami[e[0]]+lami[e[1]];
          auto symdyadic = MakeReggeAD(xi,xi);
          
          //IntLegNoBubble::
          LegendrePolynomial::
            EvalMult (p, xi, 0.25*lam_e, SBLambda ([&](int i, auto val)
                                              {
                                                shape[ii++] = val*symdyadic;
                                              }));
        }
      
      
      for (int i = 0; i<6; i++)
        {
          int p = order_facet[i][0];
          
          Tx lam_f(0);
          for (int j = 0; j < 4; j++)
            lam_f += lami[faces[i][j]];
          
          IVec<4> f = ET_trait<ET_HEX>::GetFaceSort (i, vnums);	  
          Tx xi  = sigma[f[0]] - sigma[f[1]]; 
          Tx eta = sigma[f[0]] - sigma[f[3]];
          
          LegendrePolynomial (p, eta, leg_u);
          LegendrePolynomial (p, xi, leg_v);
          
          auto symdyadic = 0.25*lam_f*MakeReggeAD(eta,xi);
          for (int j = 0; j <= p; j++)
            for (int k = 0; k <= p; k++)
              shape[ii++] = leg_u[j]*leg_v[k]*symdyadic;

          symdyadic = 0.25*lam_f*(1-eta*eta)*MakeReggeAD(xi,xi);
          for (int j = 0; j < p; j++)
            for (int k = 0; k <= p; k++)
              shape[ii++] = leg_u[j]*leg_v[k]*symdyadic;

          symdyadic = 0.25*lam_f*(1-xi*xi)*MakeReggeAD(eta,eta);
          for (int k = 0; k < p; k++)
            for (int j = 0; j <= p; j++)
              shape[ii++] = leg_u[j]*leg_v[k]*symdyadic;
          
        }

      int p = order_inner[0];
      if (p > 0)
        {
          Tx xi  = sigma[0] - sigma[1];
          Tx eta = sigma[0] - sigma[3];
          Tx nv = sigma[0] - sigma[4];
          
          LegendrePolynomial (p, xi,  leg_u);
          LegendrePolynomial (p, eta, leg_v);
          LegendrePolynomial (p, nv,  leg_w);

          auto symdyadic1 = lz[0]*lz[1]*MakeReggeAD(eta,xi);
          auto symdyadic2 = lx[0]*lx[1]*MakeReggeAD(nv,eta);
          auto symdyadic3 = ly[0]*ly[1]*MakeReggeAD(xi,nv);
          for (int i = 0; i <= p; i++)
            for (int j = 0; j <= p; j++)
              for (int k = 0; k < p; k++)
                {
                  shape[ii++] = leg_u[i]*leg_v[j]*leg_w[k]*symdyadic1;
                  shape[ii++] = leg_v[i]*leg_w[j]*leg_u[k]*symdyadic2;
                  shape[ii++] = leg_w[i]*leg_u[j]*leg_v[k]*symdyadic3;
                }

          symdyadic1 = ly[0]*ly[1]*lz[0]*lz[1]*MakeReggeAD(xi,xi);
          symdyadic2 = lz[0]*lz[1]*lx[0]*lx[1]*MakeReggeAD(eta,eta);
          symdyadic3 = lx[0]*lx[1]*ly[0]*ly[1]*MakeReggeAD(nv,nv);

          for (int i = 0; i <= p; i++)
            for (int j = 0; j < p; j++)
              for (int k = 0; k < p; k++)
                {
                  shape[ii++] = leg_u[i]*leg_v[j]*leg_w[k]*symdyadic1;
                  shape[ii++] = leg_v[i]*leg_w[j]*leg_u[k]*symdyadic2;
                  shape[ii++] = leg_w[i]*leg_u[j]*leg_v[k]*symdyadic3;
                }
        }
    }
    
    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      auto & ip = mip.IP();
      typedef typename std::remove_const<typename std::remove_reference<decltype(mip.IP()(0))>::type>::type T; 

      T x = ip(0), y = ip(1), z = ip(2);
      // T lx[4] = { 1-x, x, x, 1-x };
      // T ly[4] = { 1-y, 1-y, y, y };
      // T lz[4] = { 1-z, 1-z, z, z };
      // T lam[4] = { 1-x-y+x*y, x*(1-y), x*y, y*(1-x) };
      T sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
                   (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z};

      /*Vec<2,AutoDiff<2,T>> adip = ip;
      auto tip = TIP<2,AutoDiffDiff<2,T>>(adip);
      AutoDiffDiff<2,T> xxx = tip.x, yyy = tip.y;
      AutoDiff<2,T> xx(xxx.Value(), &xxx.DValue(0));
      AutoDiff<2,T> yy(yyy.Value(), &yyy.DValue(0));
      AutoDiff<2,T> lami[4] = {(1-xx)*(1-yy),xx*(1-yy),xx*yy,(1-xx)*yy};  
      AutoDiff<2,T> sigma[4] = {(1-xx)+(1-yy),xx+(1-yy),xx+yy,(1-xx)+yy}; */ 
      
      Vec<3,T> pnts[8] = {  { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
      int facetnr = ip.FacetNr();

      int ii = 0;
      
      ArrayMem<T,20> v(order+2), u(order+2), w(order+2);

      if (mip.IP().VB() == BBND)
        { // edge shapes
          for (int i = 0; i < 12; i++)
            {
              int p = order_edge[i];
              
              if (i == facetnr)
                {             
                  IVec<2> e = ET_trait<ET_HEX>::GetEdgeSort (i, vnums);
                  
                  T xi  = sigma[e[1]]-sigma[e[0]];
                  Vec<3,T> tauref = pnts[e[0]] - pnts[e[1]];
                  
                  
                  auto tv = mip.GetJacobian()*tauref;

                  auto tt = DyadProd(tv,tv);
                  LegendrePolynomial::Eval
                    (p, xi,
                     SBLambda([&] (size_t nr, T val)
                              {
                                shape[nr+ii] = 1/mip.GetMeasure()*val*tt;
                                }));
                 
                }
              ii += (p+1);
            }
        }
      else
        {
          for (int i = 0; i < 12; i++)
            ii += order_edge[i]+1;
        }
      if (mip.IP().VB() == BND)
        {
          for (int i = 0; i < 6; i++)
            {
              int p = order_facet[i][0];
              
              if (i == facetnr)
                {
                  IVec<4> f = ET_trait<ET_HEX>::GetFaceSort (i, vnums);	  
                  Vec<3,T> tauref1 = pnts[f[0]] - pnts[f[1]];
                  Vec<3,T> tauref2 = pnts[f[0]] - pnts[f[3]];
                  T xi  = sigma[f[0]] - sigma[f[1]]; 
                  T eta = sigma[f[0]] - sigma[f[3]];
                  //Vec<6, T> symdyadic = SymDyadProd(GetGradient(etaa),GetGradient(xia));
                  auto tv1 = mip.GetJacobian()*tauref1;
                  auto tv2 = mip.GetJacobian()*tauref2;
                  auto symdyadic = SymDyadProd(tv1,tv2);
                  
                  LegendrePolynomial (p, eta, u);
                  LegendrePolynomial (p, xi, v);
                  for (int j = 0; j <= p; j++)
                    for (int k = 0; k <= p; k++)
                      shape[ii + j*(p+1) + k] = u[j]*v[k]*symdyadic;
                      
                  /*                  T eta = ly[2]-ly[1];
          T xi = lx[1]-lx[0];
          LegendrePolynomial (p, eta, v);
          LegendrePolynomial (p, xi, u);
          
          for (int i = 0; i <= p; i++)
            for (int j = 0; j <= p; j++)
              {
                shape[ii++] = 1/mip.GetMeasure()*u[i]*v[j]*mip.GetJacobian()*Mat<2,2>(Matrix<>({{0,1},{1,0}}))*Trans(mip.GetJacobian());
              }
          
          
          //auto symdyad = lx[1]*lx[0]*SymDyadProd(Vec<2,T>(0,1),Vec<2,T>(0,1));//x*(1-x)*(0,0,  0,1) * P(y) * P(x)
          for (int i = 0; i < p; i++)
            for (int j = 0; j <= p; j++)
              {
                shape[ii++] = 1/mip.GetMeasure()*u[i]*v[j]*mip.GetJacobian()*Mat<2,2>(Matrix<>({{0,0},{0,1}}))*Trans(mip.GetJacobian());
              }
          
          //symdyad = ly[2]*ly[1]*SymDyadProd(Vec<2,T>(1,0),Vec<2,T>(1,0)); //y*(1-y)*(1,0,  0,0) * P(x) * P(y)
          
          for (int j = 0; j < p; j++)
            for (int i = 0; i <= p; i++)
              {
                shape[ii++] = 1/mip.GetMeasure()*u[i]*v[j]*mip.GetJacobian()*Mat<2,2>(Matrix<>({{1,0},{0,0}}))*Trans(mip.GetJacobian());
                }*/
                }
              ii += p*p + (p+2)*p*2 + 1;
            }
        }
      else 
        {
          for (int i = 0; i < 6; i++)
            ii += order_facet[i][0]*order_facet[i][0] + (order_facet[i][0]+2)*order_facet[i][0]*2 + 1;
        }

      if (mip.IP().VB() == VOL)
        {
          if (order_inner[0])
            throw Exception ("Hcurlcurlfe calcdualshape2 not implementend for element type ET_HEX for high-order");
        }
    }

  };


  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_SEGM>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_TRIG>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_QUAD>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_TET>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_PRISM>;
  HCURLCURLFE_EXTERN template class T_HCurlCurlFE<ET_HEX>;
}

#endif
