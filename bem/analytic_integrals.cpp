#include <bla.hpp>
#include <analytic_integrals.hpp>


namespace ngsbem
{

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
      double p = (alpha*tx + sx) / (1+sqr(alpha));
      double q = sqrt ( sqr(ux) + sqr(tx-alpha*sx) / (1+sqr(alpha)));

      /*
      cout << "q = " << q << endl;
      cout << "ux = " << ux << endl;
      cout << "s,p = " << s << ", " << p  << endl;
      cout << "(s-p)*ux = " << (s-p)*ux << endl;
      */
      
      /*
      cout << "I1 = " << (s-sx) * log (alpha*s-tx+sqrt(sqr(s-sx)+sqr(alpha*s-tx)+sqr(ux))) - s << endl;
      cout << "I2 = " <<  (alpha*sx-tx)/sqrt(1+sqr(alpha)) * (  log(sqrt(1+sqr(alpha))*(s-p)+sqrt((1+sqr(alpha))*sqr(s-p)+sqr(q))) ) << endl;
      cout << "I3 = " << 2*ux*atan2 ( ( (q-(alpha*sx-tx)/(1+sqr(alpha))) * sqrt((1+sqr(alpha))*sqr(s-p)+q*q) + (alpha*s-tx-q)*q ) ,
                                      ( (s-p)*ux )) << endl;
      */
      // cout << "atan, y = " << ( (q-(alpha*sx-tx)/(1+sqr(alpha))) * sqrt((1+sqr(alpha))*sqr(s-p)+q*q) + (alpha*s-tx-q)*q )
      // << ", x = " <<  ( (s-p)*ux ) << endl;
      
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
  double LaplaceSL_Triangle (Vec<3> v0, Vec<3> v1, Vec<3> v2, Vec<3> x)
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
  }
}
