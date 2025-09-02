#include <bla.hpp>
#include <analytic_integrals.hpp>


namespace ngsbem
{

  // formulas from Steinbach+Rjasanow, p 246
  
  
  double LaplaceSL_Triangle (Vec<3> x1, Vec<3> x2, Vec<3> x3, Vec<3> x)
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

}
