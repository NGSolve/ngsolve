// #define DEBUG
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
    if (order <= 1) throw Exception ("cannot differentiate B-spline of order <= 1");
    Array<double> cp(c.Size());
    cp = 0;
    if (t[order-1] != t[0])
      cp[0] = (order-1) * (c[0]-0) / (t[0+order-1] - t[0]);
    for (int j = 1; j < t.Size()-order+1; j++)
      if (t[j+order-1] != t[j])
        cp[j] = (order-1) * (c[j]-c[j-1]) / (t[j+order-1] - t[j]);
      else
        cp[j] = 0;
    // throw Exception ("cannot differentiate, B-spline is discontinuous");
    return BSpline (order-1, Array<double>(t), move(cp));
  }
  
  BSpline BSpline :: Integrate () const
  {
    /*
    Array<double> text(t.Size()+2);
    text[0] = t[0];
    text.Range(1, t.Size()+1) = t;
    text[text.Size()-1] = text[text.Size()-2];

    Array<double> ci(c.Size()+2);
    ci = 0;
    for (int j = 0; j < t.Size()-order; j++)
      ci[j+1] = ci[j] + c[j] * (t[j+order] - t[j]) / order;
    // better, but still not correct ...
    int last = t.Size()-order-1;
    for (int j = t.Size()-order; j < ci.Size()-1; j++)
      ci[j+1] = ci[j] + c[last] * (t[last+order] - t[j]) / order;

    cout << "integral, c = " << c << endl << "ci = " << ci << endl;
    return BSpline (order+1, move(text), move(ci));
    */

    Array<double> text(t.Size()+1);
    text.Range(0, t.Size()) = t;
    text[text.Size()-1] = text[text.Size()-2];

    Array<double> ci(c.Size()+1);
    ci = 0;
    double sum = 0;
    for (int j = 0; j < t.Size()-order; j++)
      {
        sum += c[j] * (t[j+order] - t[j]) / order;
        ci[j] = sum;
      }

    for (int j = t.Size()-order; j < t.Size()-1; j++)
      {
        sum += c[t.Size()-order] * (t[t.Size()-1] - t[j]) / order;
        ci[j] = sum;
      }
    ci[t.Size()-1] = ci[t.Size()-2];
    
    // cout << "integral, c = " << c << endl << "ci = " << ci << endl;
    return BSpline (order+1, move(text), move(ci));
  }

  
  double BSpline :: Evaluate (double x) const
  {
    // for (int m = order-1; m < t.Size()-order+1; m++)
    // for (int m = 0; m < t.Size()-order+1; m++)
    for (int m = 0; m < t.Size()-1; m++)
      {
        // cout << "m = " << m << endl;
        if ( (t[m] <= x) && (x < t[m+1]))
          {
            Array<double> hc (c);
            for (int p = 1; p < order; p++)
              for (int j = m; j >= m-order+1+p; j--)
                {
                  // cout << "j = " << j << ", p = " << p << endl;
                  if (j > 0)
                    hc[j] = ((x-t[j]) * hc[j] + (t[j+order-p]-x) * hc[j-1])
                      / (t[j+order-p]-t[j]);
                  else if (j == 0)
                    hc[j] = ((x-t[0]) * hc[j]) / (t[j+order-p]-t[j]);
                }
            return hc[m];
          }
      }
    return 0;
  }

  AutoDiff<1,double> BSpline :: operator() (AutoDiff<1,double> x) const
  {
    // double eps = 1e-5;
    double val = (*this)(x.Value());
    // double valr = (*this)(x.Value()+eps);
    // double vall = (*this)(x.Value()-eps);
    // double dval = (valr-vall) / (2*eps);

    double dval2 = Differentiate()(x.Value());
    // cout << "dval = " << dval << " =?= " << dval2 << endl;

    AutoDiff<1,double> res(val);
    res.DValue(0) = dval2 * x.DValue(0);
    return res;
  }

  AutoDiffDiff<1,double> BSpline :: operator() (AutoDiffDiff<1,double> x) const
  {
    /*
    double eps = 1e-5;
    double val = (*this)(x.Value());
    double valr = (*this)(x.Value()+eps);
    double vall = (*this)(x.Value()-eps);

    double dval = (valr-vall) / (2*eps);
    double ddval = (valr+vall-2*val) / (eps*eps);
    */
    auto diff = Differentiate();
    auto ddiff = diff.Differentiate();
    double val = (*this)(x.Value());
    double dval = diff(x.Value());
    double ddval = ddiff(x.Value());

    AutoDiffDiff<1,double> res(val);
    res.DValue(0) = dval * x.DValue(0);
    res.DDValue(0) = ddval * x.DValue(0)*x.DValue(0) + dval*x.DDValue(0);
    return res;

  }
  
  AutoDiff<1,SIMD<double>> BSpline :: operator() (AutoDiff<1,SIMD<double>> x) const
  {
    // double eps = 1e-5;
    SIMD<double> val = (*this)(x.Value());
    // double valr = (*this)(x.Value()+eps);
    // double vall = (*this)(x.Value()-eps);
    // double dval = (valr-vall) / (2*eps);

    SIMD<double> dval2 = Differentiate()(x.Value());
    // cout << "dval = " << dval << " =?= " << dval2 << endl;

    AutoDiff<1,SIMD<double>> res(val);
    res.DValue(0) = dval2 * x.DValue(0);
    return res;
  }

  AutoDiffDiff<1,SIMD<double>> BSpline :: operator() (AutoDiffDiff<1,SIMD<double>> x) const
  {
    /*
    double eps = 1e-5;
    double val = (*this)(x.Value());
    double valr = (*this)(x.Value()+eps);
    double vall = (*this)(x.Value()-eps);

    double dval = (valr-vall) / (2*eps);
    double ddval = (valr+vall-2*val) / (eps*eps);
    */
    auto diff = Differentiate();
    auto ddiff = diff.Differentiate();
    SIMD<double> val = (*this)(x.Value());
    SIMD<double> dval = diff(x.Value());
    SIMD<double> ddval = ddiff(x.Value());

    AutoDiffDiff<1,SIMD<double>> res(val);
    res.DValue(0) = dval * x.DValue(0);
    res.DDValue(0) = ddval * x.DValue(0)*x.DValue(0) + dval*x.DDValue(0);
    return res;

  }

  SIMD<double> BSpline :: operator() (SIMD<double> x) const
  {
    constexpr int size = SIMD<double>::Size();
    SIMD<double> res;
    for(auto i = 0; i < size; i++) res[i] = Evaluate(x[i]);
    return res;
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
