#include <ngstd.hpp>
#include "bspline.hpp"
using namespace ngstd;


namespace ngstd
{
  BSpline :: BSpline (int aorder, 
                      Array<double> at,
                      Array<double> ac)
    : order(aorder), t(at), c(ac) 
  { ; }
  
  
  BSpline BSpline :: Differentiate () const
  {
    if (order <= 1) throw Exception ("cannot differentiate B-spline");
    Array<double> cp(c.Size());
    cp = 0;
    for (int j = 1; j < cp.Size()-order; j++)
      cp[j] = (order-1) * (c[j]-c[j-1]) / (t[j+order-1] - t[j]);
    return BSpline (order-1, Array<double>(t), move(cp));
  }
  
  BSpline BSpline :: Integrate () const
  {
    Array<double> text(t.Size()+2);
    text[0] = t[0];
    text.Range(1, t.Size()+1) = t;
    text[text.Size()-1] = text[text.Size()-2];

    Array<double> ci(c.Size()+2);
    ci = 0;
    for (int j = 0; j < t.Size()-order; j++)
      ci[j+1] = ci[j] + c[j] * (t[j+order] - t[j]) / order;

    cout << "integral, c = " << c << endl << "ci = " << ci << endl;
    return BSpline (order+1, move(text), move(ci));
  }

  
  double BSpline :: Evaluate (double x) const
  {
    for (int m = order-1; m < t.Size()-order+1; m++)
      {
        if ( (t[m] <= x) && (x < t[m+1]))
          {
            Array<double> hc (c);
            for (int p = 1; p < order; p++)
              for (int j = m; j >= m-order+1+p; j--)
                hc[j] = ((x-t[j]) * hc[j] + (t[j+order-p]-x) * hc[j-1])
                  / (t[j+order-p]-t[j]);
            return hc[m];
          }
      }
    return 0;
  }


  ostream & operator<< (ostream & ost, const BSpline & sp)
  {
    ost << "bspline, order = " << sp.order << endl
        << "t = " << sp.t << endl
        << "c = " << sp.c << endl;
    return ost;
  }
  
}


/*
int main ()
{
  auto sp = BSpline (2, 
                     { 0, 0, 1, 2, 3, 4, 5, 6, 6 }, 
                     { 1, 0, 0, 0, 0, 0, 0, 0, 0 } );

  cout << "sp = " << sp << endl;
  auto isp = sp.Integrate();
  cout << "integral(sp)= " << isp << endl;

  ofstream out("spline.out");
  for (double t = -1; t <= 9; t += 0.01)
    {
      out << t << " " << sp.Evaluate(t) << " " << isp.Evaluate(t) << endl;
    }
}
*/
