#include <fem.hpp>
#include <analytic_integrals.hpp>


namespace ngsbem
{
  using namespace ngfem;

  
  // formulas from Steinbach+Rjasanow, p 246
  double LaplaceSL_TriangleSR (Vec<3> x1, Vec<3> x2, Vec<3> x3, Vec<3> x)
  {
    double ttau = L2Norm(x3-x2);
    Vec<3> r2 = 1/ttau * (x3-x2);

    double stau = L2Norm(Cross(x1-x2,x3-x2)) / L2Norm(x2-x3);
    double tstar = InnerProduct(x1-x2,r2);
    Vec<3> xstar = x2 + tstar * r2;
    Vec<3> r1 = 1/stau * (xstar-x1);

    Vec<3> n = Cross(r1,r2);

    double alpha1 = -tstar/stau;
    double alpha2 = (ttau-tstar) / stau;

    
    double sx = InnerProduct(x-x1,r1);
    double tx = InnerProduct(x-x1,r2);
    double ux = InnerProduct(x-x1,n);

    // cout << "sx = " << sx << ", tx = " << tx << ", ux = " << ux << endl;
    // cout << "al1 = " << alpha1 << ", al2 = " << alpha2 << endl;

    auto F = [=](double s, double alpha)
    {
      // cout << "s = " << s << ", alpha = " << alpha << endl;
      double p = (alpha*tx + sx) / (1+sqr(alpha));
      double q = sqrt ( sqr(ux) + sqr(tx-alpha*sx) / (1+sqr(alpha)));

      /*
      cout << "q = " << q << endl;
      cout << "ux = " << ux << endl;
      cout << "s,p = " << s << ", " << p  << endl;
      cout << "(s-p)*ux = " << (s-p)*ux << endl;
      cout << "alpha s - tx - q = " << (alpha*s-tx-q) << endl;
      cout << "(q-(alpha*sx-tx)/(1+sqr(alpha))) = " << (q-(alpha*sx-tx)/(1+sqr(alpha))) << endl;
      
      cout << "I1 = " << (s-sx) * log (alpha*s-tx+sqrt(sqr(s-sx)+sqr(alpha*s-tx)+sqr(ux))) - s << endl;
      cout << "I2 = " <<  (alpha*sx-tx)/sqrt(1+sqr(alpha)) * (  log(sqrt(1+sqr(alpha))*(s-p)+sqrt((1+sqr(alpha))*sqr(s-p)+sqr(q))) ) << endl;
      cout << "I3 = " << 2*ux*atan ( ( (q-(alpha*sx-tx)/(1+sqr(alpha))) * sqrt((1+sqr(alpha))*sqr(s-p)+q*q) + (alpha*s-tx-q)*q ) /
                                     ( (s-p)*ux )) << endl;

      cout << "atan, y = " << ( (q-(alpha*sx-tx)/(1+sqr(alpha))) * sqrt((1+sqr(alpha))*sqr(s-p)+q*q) + (alpha*s-tx-q)*q )
           << ", x = " <<  ( (s-p)*ux )
           << "atan = " << atan ( ( (q-(alpha*sx-tx)/(1+sqr(alpha))) * sqrt((1+sqr(alpha))*sqr(s-p)+q*q) + (alpha*s-tx-q)*q ) /
                                  ( (s-p)*ux ))
           << ", factor = " << 2*ux
           << endl;
      */
      
      // cout << "I3 = " << 2*ux*atan2 ( ( (q-(alpha*sx-tx)/(1+sqr(alpha))) * sqrt((1+sqr(alpha))*sqr(s-p)+q*q) + (alpha*s-tx-q)*q ) ,
      // ( (s-p)*ux )) << endl;
      
      return (s-sx) * log (alpha*s-tx+sqrt(sqr(s-sx)+sqr(alpha*s-tx)+sqr(ux))) - s
        + (alpha*sx-tx)/sqrt(1+sqr(alpha)) * (  log(sqrt(1+sqr(alpha))*(s-p)+sqrt((1+sqr(alpha))*sqr(s-p)+sqr(q))) )
        + 2*ux*atan ( ( (q-(alpha*sx-tx)/(1+sqr(alpha))) * sqrt((1+sqr(alpha))*sqr(s-p)+q*q) + (alpha*s-tx-q)*q ) /
                       ( (s-p)*ux ));
    };
    
    return 1/(4*M_PI) * (F(stau, alpha2)-F(0,alpha2)-F(stau,alpha1)+F(0,alpha1));
  }

  

  
  
  // from deepseek:
  // Analytical integration for single layer potential on triangle
  // Using the formula from "Fast Multipole Methods for the Helmholtz Equation 
  // in Three Dimensions" by Gumerov and Duraiswami
  double LaplaceSL_TriangleGD (Vec<3> v0, Vec<3> v1, Vec<3> v2, Vec<3> x)
  {
    
    // Vectors from vertices to observation point
    Vec<3> r0 = x - v0;
    Vec<3> r1 = x - v1;
    Vec<3> r2 = x - v2;
        
    double R0 = L2Norm(r0);
    double R1 = L2Norm(r1);
    double R2 = L2Norm(r2);
        
    // Edge vectors
    Vec<3> e0 = v1 - v0;
    Vec<3> e1 = v2 - v1;
    Vec<3> e2 = v0 - v2;
        
    // Edge lengths
    double l0 = L2Norm(e0);
    double l1 = L2Norm(e1);
    double l2 = L2Norm(e2);
        
    // Unit vectors along edges
    Vec<3> u0 = (1.0/l0) * e0;
    Vec<3> u1 = (1.0/l1) * e1;
    Vec<3> u2 = (1.0/l2) * e2;
        
    // Projections of observation point vectors onto edges
    double p0 = InnerProduct(r0,u0);
    double p1 = InnerProduct(r1,u1);
    double p2 = InnerProduct(r2,u2);
    
    // Perpendicular distances
    double h0 = std::sqrt(std::max(0.0, R0*R0 - p0*p0));
    double h1 = std::sqrt(std::max(0.0, R1*R1 - p1*p1));
    double h2 = std::sqrt(std::max(0.0, R2*R2 - p2*p2));
        
    // first guess from deepseek, now it means that is double layer:
    // Angles
    double alpha0 = std::atan2(p0, h0) - std::atan2(p0 - l0, h0);
    double alpha1 = std::atan2(p1, h1) - std::atan2(p1 - l1, h1);
    double alpha2 = std::atan2(p2, h2) - std::atan2(p2 - l2, h2);
    
    // Integration result
    double result = 0.0;
    result += h0 * alpha0;
    result += h1 * alpha1;
    result += h2 * alpha2;
    
    return result / (4*M_PI);


    /*
    // For single layer potential, we need a different formula
    // This is the correct analytical expression for ∫∫ 1/|r| dS
    auto edge_contribution = [](double p, double l, double h, double R) {
        if (h < 1e-12) { // Handle singularity case
            // Special handling when observation point is in the plane
            return p * std::log((p + R) / (p - l + std::sqrt((p - l)*(p - l) + h*h))) 
                   + l * std::log(R + p - l) - R;
        }
        
        double term1 = p * std::log((R + p) / (std::sqrt((p - l)*(p - l) + h*h) + p - l));
        double term2 = l * std::log((R + p - l) / h);
        double term3 = h * (std::atan(p/h) - std::atan((p - l)/h));
        cout << "R = " << R << ", p = " << p << ", R+p-l = " << R+p-l << ", h = " << h << endl;
        cout << "term1 = " << term1 << ", term2 = " << term2 << ", term3 = " << term3 << endl;
        return term1 + term2 + term3;
    };
    */

    /*
    auto edge_contribution = [](double p, double l, double h, double R) {
      double d1 = R + p;
      double d2 = std::sqrt((p - l)*(p - l) + h*h) + p - l;
      double d3 = R + p - l;
      
      // Calculate the sign corrections
      int sign_d1 = (d1 >= 0) ? 1 : -1;
      int sign_d2 = (d2 >= 0) ? 1 : -1;
      int sign_d3 = (d3 >= 0) ? 1 : -1;
      
      double term1 = p * (std::log(std::abs(d1)) - std::log(std::abs(d2)));
      double term2 = l * (std::log(std::abs(d3)) - std::log(h));
      double term3 = h * (std::atan2(p, h) - std::atan2(p - l, h));
      
      // Add imaginary part corrections (which cancel out in the end)
      term1 += p * (std::arg(std::complex<double>(d1, 0)) - std::arg(std::complex<double>(d2, 0)));
      term2 += l * (std::arg(std::complex<double>(d3, 0)) - 0);
      
      return term1 + term2 + term3;
    };
    */
    
    /*
    // Safe edge contribution function
    auto edge_contribution = [](double p, double l, double h, double R) {
        // Handle cases where arguments to log might be negative
        double d1 = R + p;
        double d2 = std::sqrt((p - l)*(p - l) + h*h) + p - l;
        double d3 = R + p - l;
        
        // Ensure positive arguments for logarithms
        double term1 = (std::abs(d1) > 1e-12 && std::abs(d2) > 1e-12) 
                      ? p * std::log(std::abs(d1)/std::abs(d2)) 
                      : 0.0;
        
        double term2 = (std::abs(d3) > 1e-12 && h > 1e-12)
                      ? l * std::log(std::abs(d3)/h)
                      : 0.0;
        
        double term3 = h > 1e-12 
                      ? h * (std::atan2(p, h) - std::atan2(p - l, h))
                      : 0.0;
        
        // Add sign corrections based on actual signs
        if (d1 < 0) term1 += p * M_PI * 1i; // In practice, we need real implementation
        if (d2 < 0) term1 -= p * M_PI * 1i;
        if (d3 < 0) term2 += l * M_PI * 1i;
        
        return term1 + term2 + term3;
    };

    
    double result = 0.0;
    result += edge_contribution(p0, l0, h0, R0);
    result += edge_contribution(p1, l1, h1, R1);
    result += edge_contribution(p2, l2, h2, R2);
    
    return result / (4 * M_PI);
    */
  }




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
    return sign(z) * atan (x / abs(z));
  }

  double k1 (double x, double y, double z)
  {
    double r = sqrt(x*x+y*y+z*z);
    if (z == 0) return 0;
    return y*sign(z) * atan (x*y / (abs(z)*r)) + z * log(r+x);
  }

  double km1 (double x, double y, double z)
  {
    if (z == 0) return 0;    
    double r = sqrt(x*x+y*y+z*z);
    if (y == 0) return sign(z) * x / (abs(z)*r);
    return sgn(z)/y * atan (y*x / (abs(z)*r));
  }

  double H (double x, double y, double z)
  {
    return 1/(4*M_PI) * (k1(x,y,z) - y * k0(x,y,z));
  }
  
    
  double HDL (double x, double y, double z)
  {
    return 1/(4*M_PI) * (y * km1(x,y,z) - sgn(y) * k0(x,y,z));
  }
  
  
  double LaplaceSL_Triangle (Vec<3> v0, Vec<3> v1, Vec<3> v2, Vec<3> x)
  {
    Vec<3> vi[3] = { v0, v1, v2 };
    Vec<3> n = Cross(v1-v0, v2-v0);
    n /= L2Norm(n);

    double h = InnerProduct(x-v0, n);
    Vec<3> xhat = x - h*n;
    h = fabs(h);

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





  

  
  double LaplaceSL_Triangle_exp (Vec<3> v0, Vec<3> v1, Vec<3> v2, Vec<3> x)
  {
    Vec<3> vi[3] = { v0, v1, v2 };
    Vec<3> n = Cross(v1-v0, v2-v0);
    n /= L2Norm(n);

    double h = InnerProduct(x-v0, n);
    Vec<3> xhat = x - h*n;
    h = fabs(h);

    // cout << "x = " << x << ", h = " << h << ", xhat = " << xhat << endl;
    auto f = [=] (Vec<3> y, Vec<3> nu) {
      // double r = L2Norm(x-y);
      // double rho = sqrt(std::max(0., r*r-h*h));
      double rho = L2Norm(xhat-y);
      double r = sqrt(rho*rho+h*h);
      return 1/(4*M_PI) * (r-h)/(rho*rho) * InnerProduct(xhat-y, nu);
    };

    IntegrationRule ir(ET_SEGM, 20);
    auto IntEdge = [=,&ir] (Vec<3> ve0, Vec<3> ve1, Vec<3> nu)
    {
      double fac = L2Norm(ve1-ve0);
      double sum = 0;
      for (auto ip : ir)
        {
          Vec<3> y = ve0 + ip(0) * (ve1-ve0);
          sum += ip.Weight()*fac * f(y, nu);
        }

      // cout << "num edge int = " << sum;

      Vec<3> ex = ve1 - ve0; ex /= L2Norm(ex);
      Vec<3> ey = n;
      Vec<3> ez = nu;
      double xp = InnerProduct(x-ve0, ex);
      double yp = InnerProduct(x-ve0, ey);
      double zp = InnerProduct(x-ve0, ez);
      // cout << "ve0= " << ve0 << ", ve1 = " << ve1 << endl;
      // cout << "x-ve0 = " << x-ve0 << endl;
      // cout << "ex = " << ex << endl;
        
      // cout << "xp = " << xp << endl;
      double l = L2Norm(ve0-ve1);
      // cout << "xp = " << xp << ", yp = " << yp << ", zp = " << zp << endl;
      // cout << "H1 = " << H(l-xp, yp, zp) << ", H2 = " <<  H(-xp, yp, zp) << endl;
      double anaint = H(l-xp, yp, zp) - H(-xp, yp, zp);
      // cout << " =?= " << anaint << " = analytic int" << endl;
      
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






  double LaplaceDL_Triangle (Vec<3> v0, Vec<3> v1, Vec<3> v2, Vec<3> x)
  {
    Vec<3> vi[3] = { v0, v1, v2 };
    Vec<3> n = Cross(v1-v0, v2-v0);
    n /= L2Norm(n);

    double h = InnerProduct(x-v0, n);
    Vec<3> xhat = x - h*n;

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

  



  
  double LaplaceDL_Triangle_exp (Vec<3> v0, Vec<3> v1, Vec<3> v2, Vec<3> x)
  {
    Vec<3> vi[3] = { v0, v1, v2 };
    Vec<3> n = Cross(v1-v0, v2-v0);
    n /= L2Norm(n);

    double h = InnerProduct(x-v0, n);
    Vec<3> xhat = x - h*n;
    // h = fabs(h);

    // cout << "x = " << x << ", h = " << h << ", xhat = " << xhat << endl;
    auto f = [=] (Vec<3> y, Vec<3> nu) {
      // double r = L2Norm(x-y);
      // double rho = sqrt(std::max(0., r*r-h*h));
      double rho = L2Norm(xhat-y);
      double r = sqrt(rho*rho+h*h);
      // return 1/(4*M_PI) * (r-h)/(rho*rho) * InnerProduct(xhat-y, nu);
      return 1/(4*M_PI) * (h/r-sgn(h))/(rho*rho) * InnerProduct(xhat-y, nu);
    };

    IntegrationRule ir(ET_SEGM, 20);
    auto IntEdge = [=,&ir] (Vec<3> ve0, Vec<3> ve1, Vec<3> nu)
    {
      /*
      double fac = L2Norm(ve1-ve0);
      double sum = 0;
      for (auto ip : ir)
        {
          Vec<3> y = ve0 + ip(0) * (ve1-ve0);
          sum += ip.Weight()*fac * f(y, nu);
        }
      */

      // cout << "num edge int = " << sum;

      Vec<3> ex = ve1 - ve0; ex /= L2Norm(ex);
      Vec<3> ey = n;
      Vec<3> ez = nu;
      double xp = InnerProduct(x-ve0, ex);
      double yp = InnerProduct(x-ve0, ey);
      double zp = InnerProduct(x-ve0, ez);
      // cout << "ve0= " << ve0 << ", ve1 = " << ve1 << endl;
      // cout << "x-ve0 = " << x-ve0 << endl;
      // cout << "ex = " << ex << endl;
        
      // cout << "xp = " << xp << endl;
      double l = L2Norm(ve0-ve1);
      // cout << "xp = " << xp << ", yp = " << yp << ", zp = " << zp << endl;
      // cout << "H1 = " << H(l-xp, yp, zp) << ", H2 = " <<  H(-xp, yp, zp) << endl;
      double anaint = HDL(l-xp, yp, zp) - HDL(-xp, yp, zp);
      // cout << " =?= " << anaint << " = analytic int" << endl;
      
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
}
  



