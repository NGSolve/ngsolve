#include <finiteelement.hpp>
#include <analytic_integrals.hpp>


namespace ngsbem
{
  using namespace ngfem;

  /*
    Analytical computation of boundary integrals for the Helmholtz
    equation in three dimensions
    Gumerov+Duraiswami, 21
    
    G = 1 / 4pi / |x-y|
    xhat ... projection of x to plane
    r = |x-y|, rho = |xhat-y|,  r^2 = rho^2 + h^2
    
    using surface div-theorem:
    G = div_s f    with
    f = (r-h)/(rho*rho) * (xhat-y)

    \int_E  f*nu ds = H(l-x,y,z)-H(-x,y,z)
    H =  something ( k0, k1, ...)
    km(x,y,z) = z \int r^m / (x^2+z^2) dx   with r = sqrt(x*x+y*y+z*z)
  */
  
  double sign (double s)
  {
    if (s > 0) return 1;
    if (s < 0) return -1;
    return 0;
  }
  double k0 (double x, double y, double z)
  {
    if (z == 0) return 0;
    return sign(z) * atan2 (x, abs(z));
  }

  double log_r_plus_x (double x, double y, double z)
  {
    double r = sqrt(x*x+y*y+z*z);
    if (x < 0)
      {
        // log(r+x) = log(y*y+z*z) - log(r-x), avoids cancellation.
        double rho2 = y*y+z*z;
        if (rho2 == 0) return 0;
        return log(rho2) - log(r-x);
      }
    return log(r+x);
  }

  double k1 (double x, double y, double z)
  {
    double r = sqrt(x*x+y*y+z*z);
    if (z == 0) return 0;
    return y*sign(z) * atan2 (x*y, abs(z)*r) + z * log_r_plus_x(x, y, z);
  }

  double km1 (double x, double y, double z)
  {
    if (z == 0) return 0;    
    double r = sqrt(x*x+y*y+z*z);
    if (y == 0) return sign(z) * x / (abs(z)*r);
    return sgn(z)/y * atan2 (y*x, abs(z)*r);
  }

  double H (double x, double y, double z)
  {
    return 1/(4*M_PI) * (k1(x,y,z) - y * k0(x,y,z));
  }

  double dHdx (double x, double y, double z)
  {
    double den = x*x+z*z;
    if (den == 0) return 0;
    double r = sqrt(x*x+y*y+z*z);
    return -1/(4*M_PI) * z * (r-y) / den;
  }

  double dHdy (double x, double y, double z)
  {
    return 1/(4*M_PI) * (k0(x,y,z) - y * km1(x,y,z));
  }

  double dHdz (double x, double y, double z)
  {
    double den = x*x+z*z;
    double r = sqrt(x*x+y*y+z*z);
    double tangential = (den == 0) ? 0 : x * (r-y) / den;
    return 1/(4*M_PI) * (tangential - log_r_plus_x(x,y,z));
  }
  
    
  double HDL (double x, double y, double z)
  {
    return 1/(4*M_PI) * (y * km1(x,y,z) - sgn(y) * k0(x,y,z));
  }
  
  
  double LaplaceSL_Triangle (Vec<3> v0, Vec<3> v1, Vec<3> v2, Vec<3> x)
  {
    Vec<3> n = Cross(v1-v0, v2-v0);
    if (InnerProduct(x-v0, n) < 0)
      {
        Swap(v1, v2);
        n *= -1;
      }
    n /= L2Norm(n);

    Vec<3> vi[3] = { v0, v1, v2 };

    auto IntEdge = [=] (Vec<3> ve0, Vec<3> ve1, Vec<3> nu)
    {
      Vec<3> ex = ve1 - ve0; ex /= L2Norm(ex);
      Vec<3> ey = n;
      Vec<3> ez = nu;
      double xp = InnerProduct(x-ve0, ex);
      double yp = InnerProduct(x-ve0, ey);
      double zp = InnerProduct(x-ve0, ez);

      double l = L2Norm(ve0-ve1);
      double anaint = H(l-xp, yp, zp) - H(-xp, yp, zp);
      return anaint;
    };

    double sum = 0;
    for (int i = 0; i < 3; i++)
      {
        auto ve0 = vi[i];
        auto ve1 = vi[(i+1)%3];
        Vec<3> tau = ve1-ve0;
        tau /= L2Norm(tau);
        double intedge = IntEdge (ve0, ve1, Cross(n, tau));
        sum += intedge;
      };
    return sum;
  }


  Vec<3> LaplaceGradSL_Triangle (Vec<3> v0, Vec<3> v1, Vec<3> v2, Vec<3> x)
  {
    Vec<3> n = Cross(v1-v0, v2-v0);
    if (InnerProduct(x-v0, n) < 0)
      {
        Swap(v1, v2);
        n *= -1;
      }
    n /= L2Norm(n);

    Vec<3> vi[3] = { v0, v1, v2 };

    auto IntEdge = [=] (Vec<3> ve0, Vec<3> ve1, Vec<3> nu)
    {
      Vec<3> ex = ve1 - ve0; ex /= L2Norm(ex);
      Vec<3> ey = n;
      Vec<3> ez = nu;
      double xp = InnerProduct(x-ve0, ex);
      double yp = InnerProduct(x-ve0, ey);
      double zp = InnerProduct(x-ve0, ez);

      auto HGrad = [=] (double s)
      {
        return -dHdx(s, yp, zp) * ex
          + dHdy(s, yp, zp) * ey
          + dHdz(s, yp, zp) * ez;
      };

      double l = L2Norm(ve0-ve1);
      return HGrad(l-xp) - HGrad(-xp);
    };

    Vec<3> sum { 0.0, 0.0, 0.0 };
    for (int i = 0; i < 3; i++)
      {
        auto ve0 = vi[i];
        auto ve1 = vi[(i+1)%3];
        Vec<3> tau = ve1-ve0;
        tau /= L2Norm(tau);
        sum += IntEdge (ve0, ve1, Cross(tau, n));
      };
    return sum;
  }


  double LaplaceDL_TriangleOriented (Vec<3> v0, Vec<3> v1, Vec<3> v2, Vec<3> x,
                                     Vec<3> n)
  {
    Vec<3> vi[3] = { v0, v1, v2 };

    // double h = InnerProduct(x-v0, n);
    // Vec<3> xhat = x - h*n;

    auto IntEdge = [=] (Vec<3> ve0, Vec<3> ve1, Vec<3> nu)
    {
      Vec<3> ex = ve1 - ve0; ex /= L2Norm(ex);
      Vec<3> ey = n;
      Vec<3> ez = nu;
      double xp = InnerProduct(x-ve0, ex);
      double yp = InnerProduct(x-ve0, ey);
      double zp = InnerProduct(x-ve0, ez);

      double l = L2Norm(ve0-ve1);
      double anaint = HDL(l-xp, yp, zp) - HDL(-xp, yp, zp);
      return anaint;
    };

    double sum = 0;
    for (int i = 0; i < 3; i++)
      {
        auto ve0 = vi[i];
        auto ve1 = vi[(i+1)%3];
        Vec<3> tau = ve1-ve0;
        tau /= L2Norm(tau);
        double intedge = IntEdge (ve0, ve1, Cross(n, tau));
        sum += intedge;
      };
    return sum;
  }


  double LaplaceDL_Triangle (Vec<3> v0, Vec<3> v1, Vec<3> v2, Vec<3> x,
                             Vec<3> source_normal)
  {
    Vec<3> n = Cross(v1-v0, v2-v0);
    if (InnerProduct(n, source_normal) > 0)
      {
        Swap(v1, v2);
        n *= -1;
      }
    n /= L2Norm(n);
    return LaplaceDL_TriangleOriented(v0, v1, v2, x, n);
  }

}
