// #define DEBUG
#include <ngstd.hpp>
#include "bspline.hpp"
using namespace ngstd;

namespace ngstd
{
  BSpline :: BSpline (int aorder, 
                      Array<double> at,
                      Array<double> ac)
    : order(aorder), t(at.Size() + order), c(ac.Size() + order) 
  {
    int j = 0;
    for ( ; j < order ; j++)
      {
        c[j] = 0;
        t[j] = at[0] - order + j;
      }
    c.Range(order,ac.Size()+order) = ac;
    t.Range(order,at.Size()+order) = at;

    if(order >= 2)
      diff = make_shared<BSpline>(Differentiate());
  }
  
  
  BSpline BSpline :: Differentiate () const
  {
    if (order <= 1) throw Exception ("cannot differentiate B-spline of order <= 1");
    //we should create td and cd WITHOUT padding on the leftmost elements
    Array<double> cd(c.Size()-1);
    Array<double> td(t.Size()-2);
    td = t.Range(1,t.Size()-1);
    for (int j = 0; j < c.Size() - 1; ++j)
      if (t[j+order] != t[j+1])
        cd[j] = (order-1) * (c[j+1]-c[j]) / (t[j+order] - t[j+1]);
      else
        cd[j] = 0;
    // throw Exception ("cannot differentiate, B-spline is discontinuous");
    return BSpline (order-1, std::move(td), std::move(cd));
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

    //we should create text and ci WITHOUT padding on the leftmost elements
    const auto origsize = t.Size() - order;
    Array<double> text(origsize+1);
    text.Range(0, origsize) = t.Range(order,t.Size());
    text[text.Size()-1] = text[text.Size()-2];

    Array<double> ci(origsize+1);
    ci = 0;
    double sum = 0;
    for (int j = order; j < t.Size()-order; j++)
      {
        sum += c[j] * (t[j+order] - t[j]) / order;
        ci[j-order] = sum;
      }

    for (int j = t.Size()-order; j < t.Size()-1; j++)
      {
        sum += c[t.Size()-order] * (t[t.Size()-1] - t[j]) / order;
        ci[j-order] = sum;
      }
    ci[origsize] = ci[origsize-1];
    
    // cout << "integral, c = " << c << endl << "ci = " << ci << endl;
    return BSpline (order+1, std::move(text), std::move(ci));
  }

  
  double BSpline :: Evaluate (double x) const
  {
    // static Timer timer_bspline("BSpline::Evaluate");
    // timer_bspline.AddFlops(1);
    // RegionTimer reg (timer_bspline);

    // for (int m = order-1; m < t.Size()-order+1; m++)
    // for (int m = 0; m < t.Size()-order+1; m++)
    for (int m = order; m < t.Size()-1; m++)
      {
        // cout << "m = " << m << endl;
        
        if ( (t[m] <= x) && (x < t[m+1]))
          {
            if(order < 6)
              {
                double val = 0;
                
                Switch<6>(order, [&](auto orderval){
                  if constexpr (orderval.value == 0 ) return;
                  else{ // if constexpr orderval

                    // constexpr auto actualorder = ORDERm1.value_type() + 1;
                    // extern int myvar; myvar = actualorder;
                    // constexpr auto orderval = actualorder;
                    const auto offset = m - orderval.value + 1;
                    double hc[orderval.value];
                    // Vec<orderval, double> hc;

                    {
                      // static Timer timer_bspline_cp("BSpline::Evaluate::Copy");
                      // NgProfiler::AddThreadFlops(timer_bspline_cp, TaskManager::GetThreadId(), 1);
                      // RegionTimer reg_ev (timer_bspline_cp);
                      for(int j = 0; j < orderval.value; j++)
                        {
                          hc[j] = c[j + offset];
                        }
                    }


                    // extern int myvar2; myvar2 = actualorder;
                    {
                      // static Timer timer_bspline_ev("BSpline::Evaluate::Eval");
                      // NgProfiler::AddThreadFlops(timer_bspline_ev, TaskManager::GetThreadId(), 1);
                      // RegionTimer reg_ev (timer_bspline_ev);
                      Iterate<orderval.value-1> ([&] (auto Pm1) {
                        constexpr int p = Pm1.value + 1;
                        for (int jhc = orderval.value - 1; jhc >= p; jhc--)
                          {
                            const auto j = jhc+offset;
                            hc[jhc] = ((x-t[j]) * hc[jhc] +
                                       (t[j+orderval.value-p]-x) * hc[jhc-1])
                              / (t[j+orderval.value-p]-t[j]);
                          }
                      } );
                    }
                    val = hc[orderval.value-1];
                  }
                });
                return val;
              }
            else
              {
                Array<double> hc (c);
                for (int p = 1; p < order; p++)
                  for (int j = m; j >= m-order+1+p; j--)
                    {
                      // cout << "j = " << j << ", p = " << p << endl;
                      //because of padding there is no longer the need to check if j == 0
                      hc[j] = ((x-t[j]) * hc[j] + (t[j+order-p]-x) * hc[j-1])
                        / (t[j+order-p]-t[j]);
                    }
                return hc[m];
              }
          }
      }
    return 0;
  }

  SIMD<double> BSpline :: Evaluate (SIMD<double> x) const
  {
    // static Timer timer_bspline("BSpline::Evaluate SIMD");
    // RegionTimer reg (timer_bspline);
    constexpr int simdSize = SIMD<double>::Size();
    if( order < 6)
      {
        SIMD<int64_t> pos(-1);
        for(int m = order;  m < t.Size() - 1; m++)
          {
            SIMD<double> currPt(t[m]);//leftmost pt of current interval
            SIMD<double> nextPt(t[m+1]);//rightmost pt of current interval
            //shouldUpdate will be set to true
            //if found the correct interval for a given pt
            auto shouldUpdate = (currPt <= x) && (x < nextPt) && (pos < SIMD<int64_t>(0));
            //update pos list
            pos = If (shouldUpdate, SIMD<int64_t>(m), pos);
          }
    
        //pts that were not found should be set to order - 1 so offset = 0
        //and only points where c==0 are taken into account
        pos = If ( pos < SIMD<int64_t>(0), SIMD<int64_t>(order-1),pos);
        //now that the pos vec is correctly set we need to actually
        //calculate the values
        SIMD<double> val(0);//this will store the result
        Switch<6>(order, [&](auto orderval){
          if constexpr (orderval.value == 0 ) return;
          else{ // if constexpr orderval
            const auto offset = pos - SIMD<int64_t>(orderval.value) + 1;
            // SIMD<mask64>offset(offsetd.Data());

            std::array<SIMD<double>,orderval.value> hc;
            // SIMD<double> hc[orderval.value];
            // HTArray<orderval.value, SIMD<double>> hc;
            // Vec<orderval, double> hc;

            {
              // static Timer timer_bspline_cp("BSpline::Evaluate::Copy");
              // NgProfiler::AddThreadFlops(timer_bspline_cp, TaskManager::GetThreadId(), 1);
              // RegionTimer reg_ev (timer_bspline_cp);
              for(int j = 0; j < orderval.value; j++)
                {
                hc[j] = SIMD<double>([&](int i) -> double {
                                         return c[j+offset[i]];});
                }
            }


            // extern int myvar2; myvar2 = actualorder;
            {
              // static Timer timer_bspline_ev("BSpline::Evaluate::Eval");
              // NgProfiler::AddThreadFlops(timer_bspline_ev, TaskManager::GetThreadId(), 1);
              // RegionTimer reg_ev (timer_bspline_ev);

              //this is a workaround attempt of a MSVC non-conformance with the standard
              // auto *dummyhc = &hc;
              Iterate<orderval.value-1> ([&] (auto Pm1) {
                constexpr int p = Pm1.value + 1;
                for (int jhc = orderval.value - 1; jhc >= p; jhc--)
                  {

                    SIMD<double> tj([&](int i) -> double { return t[jhc + offset[i]]; });
                    SIMD<double> tjopd([&](int i) -> double { return t[jhc + offset[i]+orderval.value-p]; });                    
                    /*
                    (*dummyhc)[jhc] = ((x-tj) * (*dummyhc)[jhc] +
                                       (tjopd-x) * (*dummyhc)[jhc-1])
                      / (tjopd-tj);
                    */
                    hc[jhc] = ((x-tj) * hc[jhc] +
                               (tjopd-x) * hc[jhc-1])
                      / (tjopd-tj);
                    
                  }
              } );
            }
            val = hc[orderval.value-1];
          }
        });
        return val;
      }
    else
      {
        /*
        SIMD<double> res;
        for(auto i = 0; i < simdSize; i++) res[i] = Evaluate(x[i]);
        return res;
        */
        double res[simdSize];
        for(auto i = 0; i < simdSize; i++) res[i] = Evaluate(x[i]);
        return SIMD<double>(&res[0]);
      }
    
  }
  
  AutoDiff<1,double> BSpline :: operator() (AutoDiff<1,double> x) const
  {
    if(!diff) throw Exception("BSpline::operator()(AutoDiff) not implemented for order < 2");
    // double eps = 1e-5;
    double val = (*this)(x.Value());
    // double valr = (*this)(x.Value()+eps);
    // double vall = (*this)(x.Value()-eps);
    // double dval = (valr-vall) / (2*eps);

    double dval2 = (*diff)(x.Value());
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
    if(!diff) throw Exception("BSpline::operator()(AutoDiff) not implemented for order < 2");
    // double eps = 1e-5;
    SIMD<double> val = (*this)(x.Value());
    // double valr = (*this)(x.Value()+eps);
    // double vall = (*this)(x.Value()-eps);
    // double dval = (valr-vall) / (2*eps);

    SIMD<double> dval2 = (*diff)(x.Value());
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
    // SIMD<double> res;
    // for(auto i = 0; i < size; i++) res[i] = Evaluate(x[i]);
    // return res;
    double res[size];
    for(auto i = 0; i < size; i++) res[i] = Evaluate(x[i]);
    return SIMD<double>(&res[0]);
  }

  /*
  SIMD<Complex> BSpline :: operator() (SIMD<Complex> x) const
  {
    return { (*this)(x.real()) };
  }
  */
  
  ostream & operator<< (ostream & ost, const BSpline & sp)
  {
    ost << "bspline, order = " << sp.order << endl
        << "t = " << sp.t << endl
        << "c = " << sp.c << endl;
    return ost;
  }

  std::tuple<int,int,bool> BSpline2D :: Search(double x, double y) const
  {
    int ix=0, iy=0;
    if(!extrapolate)
    {
      if(x<px[0] || x>px.Last())
        throw Exception("x out of range: " + ToString(x) + "\t[" + ToString(px[0]) + "," + ToString(px.Last())+"]");

      if(y<py[0] || y>py.Last())
        throw Exception("y out of range: " + ToString(y) + "\t[" + ToString(py[0]) + "," + ToString(py.Last())+"]");
    }

    if(x<px[0]) ix = -1;
    else if(x>px.Last()) ix = px.Size();
    else
      for(auto i : Range(px.Size()-1))
        if(px[i]<=x && x<=px[i+1]) {
          ix = i;
          break;
        }

    if(y<py[0]) iy = -1;
    else if(y>py.Last()) iy = py.Size();
    else
      for(auto i : Range(py.Size()-1))
        if(py[i]<=y && y<=py[i+1]) {
          iy = i;
          break;
        }

    bool outside = min(ix,iy)==-1 || iy==py.Size() || ix==px.Size();
    return {ix, iy, outside};
  }

  double BSpline2D :: Evaluate (double x, double y) const
  {
    auto [ix, iy, outside] = Search(x,y);
    if(outside)
    {
      auto xbnd = min(max(x, px[0]), px.Last());
      auto ybnd = min(max(y, py[0]), py.Last());
      auto dx = AutoDiff<1,double>(xbnd, x-xbnd);
      auto dy = AutoDiff<1,double>(ybnd, y-ybnd);
      auto dval = (*this)(dx, dy);
      return dval.Value()+dval.DValue(0);
    }
    else
    {
      return Interpolate(ix, iy, x, y);
    }
  }

  SIMD<double> BSpline2D :: Evaluate(SIMD<double> x, SIMD<double> y) const
  {
    SIMD<double> res;
    for(auto i : Range(res.Size()))
      reinterpret_cast<double*>(&res)[i] = Evaluate(x[i], y[i]);
    return res;
  }

  AutoDiff<1,double> BSpline2D :: operator() (AutoDiff<1,double> x, AutoDiff<1,double> y) const
  {
    auto [ix, iy, outside] = Search(x.Value(),y.Value());

    if(outside)
    {
      throw Exception("BSpline2D::operator()(AutoDiff<1,double>) out of range");
    }

    return Interpolate(ix, iy, x, y);
  }

  AutoDiffDiff<1,double> BSpline2D :: operator() (AutoDiffDiff<1,double> x, AutoDiffDiff<1,double> y) const
  {
    throw Exception("BSpline2D::AutoDiffDiff eval not implemented");
  }

  AutoDiff<1,SIMD<double>> BSpline2D :: operator() (AutoDiff<1,SIMD<double>> x, AutoDiff<1,SIMD<double>> y) const
  {
    throw ExceptionNOSIMD();
  }

  AutoDiffDiff<1,SIMD<double>> BSpline2D :: operator() (AutoDiffDiff<1,SIMD<double>> x, AutoDiffDiff<1,SIMD<double>> y) const
  {
    throw ExceptionNOSIMD();
  }


  SIMD<double> BSpline2D :: operator() (SIMD<double> x, SIMD<double> y) const
  {
    return Evaluate(x,y);
  }


  ostream & operator<< (ostream & ost, const BSpline2D & sp)
  {
    ost << "bspline2d" << endl
      << "x = " << sp.px << endl
      << "y = " << sp.py << endl
      << "v = " << sp.v << endl;
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

