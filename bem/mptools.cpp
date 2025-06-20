#include "mptools.hpp"
#include "mp_coefficient.hpp"
#include "meshaccess.hpp"


namespace ngsbem
{

  template <typename entry_type>
  entry_type SphericalHarmonics<entry_type> :: Eval (double theta, double phi) const
  {
    static Timer t("mptool sh evaluate"); RegionTimer rg(t);
    
    Matrix legfunc(order+1, order+1);
    NormalizedLegendreFunctions (order, order, cos(theta), legfunc);
    Vector<Complex> exp_imphi(order+1);
    Complex exp_iphi(cos(phi), sin(phi));
    Complex prod = 1.0;
    for (int i = 0; i <= order; i++)
      {
        exp_imphi(i) = prod;
        prod *= -exp_iphi;
      }

    entry_type sum{0.0};
    int ii = 0;
    for (int n = 0; n <= order; n++)
      {
        for (int m = -n; m < 0; m++, ii++)
          sum += conj(exp_imphi(-m)) * legfunc(-m, n) * coefs(ii);
        for (int m = 0; m <= n; m++, ii++)
          sum += exp_imphi(m) * legfunc(m, n) * coefs(ii);
      }

    sum /= sqrt(4*M_PI);
    return sum;
  }


  template <typename entry_type>  
  entry_type SphericalHarmonics<entry_type> :: EvalOrder (int n, double theta, double phi) const
  {
    static Timer t("mptool sh evalorder");

    RegionTimer rg(t);
      
    Matrix legfunc(order+1, order+1);
    NormalizedLegendreFunctions (order, order, cos(theta), legfunc);
    Vector<Complex> exp_imphi(order+1);
    Complex exp_iphi(cos(phi), sin(phi));
    Complex prod = 1.0;
    for (int i = 0; i <= order; i++)
      {
        exp_imphi(i) = prod;
        prod *= -exp_iphi;
      }

    entry_type sum{0.0};
    auto coefsn = CoefsN(n);
    int ii = 0;
    // for (int n = 0; n <= order; n++)
    {
      for (int m = -n; m < 0; m++, ii++)
        sum += conj(exp_imphi(-m)) * legfunc(-m, n) * coefsn(ii);
      for (int m = 0; m <= n; m++, ii++)
        sum += exp_imphi(m) * legfunc(m, n) * coefsn(ii);
    }

    sum /= sqrt(4*M_PI);
    return sum;
  }



  template <typename entry_type>
  void SphericalHarmonics<entry_type> :: EvalOrders (double theta, double phi, FlatVector<entry_type> vals) const
  {
    static Timer ts("mptool sh evalorders small");
    static Timer tl("mptool sh evalorders large");
    // RegionTimer rg(order < 30 ? ts : tl);
    // if (order > 30) tl.Start();
      
    Matrix legfunc(order+1, order+1);
    NormalizedLegendreFunctions (order, order, cos(theta), legfunc);
    Vector<Complex> exp_imphi(order+1);
    Complex exp_iphi(cos(phi), sin(phi));
    Complex prod = 1.0;
    for (int i = 0; i <= order; i++)
      {
        exp_imphi(i) = prod;
        prod *= -exp_iphi;
      }


    for (int n = 0; n <= order; n++)
      {
        auto coefsn = CoefsN(n);
        int ii = 0;
          
        entry_type sum{0.0};
        for (int m = -n; m < 0; m++, ii++)
          sum += conj(exp_imphi(-m)) * legfunc(-m, n) * coefsn(ii);
        for (int m = 0; m <= n; m++, ii++)
          sum += exp_imphi(m) * legfunc(m, n) * coefsn(ii);
        vals(n) = sum;
      }

    vals /= sqrt(4*M_PI);

    // if (order > 30) tl.Stop();      
  }

  
  template <typename entry_type>
  void SphericalHarmonics<entry_type> :: Calc (Vec<3> x, FlatVector<Complex> shapes)
  {
    auto [theta, phi] = Polar(x);

    Matrix legfunc(order+1, order+1);
    NormalizedLegendreFunctions (order, order, cos(theta), legfunc);
    Vector<Complex> exp_imphi(order+1);
    Complex exp_iphi(cos(phi), sin(phi));
    Complex prod = 1.0;
    for (int i = 0; i <= order; i++)
      {
        exp_imphi(i) = prod;
        prod *= -exp_iphi;
      }

    int ii = 0;
    for (int n = 0; n <= order; n++)
      {
        for (int m = -n; m < 0; m++, ii++)
          shapes(ii) = conj(exp_imphi(-m)) * legfunc(-m, n);
        for (int m = 0; m <= n; m++, ii++)
          shapes(ii) = exp_imphi(m) * legfunc(m, n);
      }

    shapes /= sqrt(4*M_PI);
  }

  
  

  // Nail A. Gumerov and Ramani Duraiswami book, formula (2.2.12)
  // add directional derivative divided by kappa to res, both multipoles need same scaling
  template <typename entry_type>  
  void SphericalHarmonics<entry_type> ::
  DirectionalDiffAdd (Vec<3> d, SphericalHarmonics<entry_type> & res, double scale) const
  {
    double fx = d(0);
    double fy = d(1);
    double fz = d(2);
    double invscale = 1./scale;
      
    for (int n = 0; n < order; n++)
      for (int m = -n; m <= n; m++)
        {
          double amn = CalcAmn(m,n);
          double bmn1 = CalcBmn(-m-1, n+1);
          double bmn2 = CalcBmn(m-1,n+1);
            
          res.Coef(n+1, m-1) += invscale * Complex(0.5*fx,0.5*fy)*bmn2 * Coef(n,m);
          res.Coef(n+1, m  ) -= invscale * fz*amn * Coef(n,m);                
          res.Coef(n+1, m+1) += invscale * Complex(0.5*fx,-0.5*fy)*bmn1 * Coef(n,m);
            
          res.Coef(n, m) += scale * Complex(-0.5*fx,0.5*fy)*bmn2 * Coef(n+1,m-1);
          res.Coef(n, m) += scale * fz*amn * Coef(n+1,m);
          res.Coef(n, m) += scale * Complex(-0.5*fx,-0.5*fy)*bmn1 * Coef(n+1,m+1);
        }
  }


  


  template <typename entry_type>
  void SphericalHarmonics<entry_type> :: RotateZ (double alpha)
  {
    // static Timer t("mptool sh RotateZ"); RegionTimer rg(t);
    if (order < 0) return;

    Vector<Complex> exp_imalpha(order+1);
    Complex exp_ialpha(cos(alpha), sin(alpha));
    Complex prod = 1.0;
    for (int i = 0; i <= order; i++)
      {
        exp_imalpha(i) = prod;
        prod *= exp_ialpha;
      }

    int ii = 0;
    for (int n = 0; n <= order; n++)
      {
        for (int m = -n; m < 0; m++, ii++)
          coefs(ii) *= conj(exp_imalpha(-m));
        for (int m = 0; m <= n; m++, ii++)
          coefs(ii) *= exp_imalpha(m);
      }
  }



  template <typename entry_type>    
  void SphericalHarmonics<entry_type> :: RotateY (double alpha)
  {
    LocalHeap lh(8*6*sqr(order) + 8*15*order + 500);
      
    // static Timer t("mptool sh RotateY"); RegionTimer rg(t);
    /*
      static std::map<int, unique_ptr<Timer<>>> timers;
      static std::map<int, unique_ptr<Timer<>>> timerstrafo;      
      if (timers.find(order)==timers.end())
      {
      timers[order] = make_unique<Timer<>> (string("mptools sh RotateY ")+ToString(order));
      timerstrafo[order] = make_unique<Timer<>> (string("mptools sh RotateY trafo ")+ToString(order));
      }
      RegionTimer rg(*timers[order]);
    */
    double s = sin(alpha);
    double c = cos(alpha);

    FlatMatrix<> normalized_leg_func(order+2, order+2, lh);
    NormalizedLegendreFunctions(order+1, order+1, c, normalized_leg_func);

    if (alpha < 0)
      for (int i = 1; i <= order+1; i+=2)
        normalized_leg_func.Row(i) *= -1;
      
    // cout << "leg = " << endl << normalized_leg_func << endl;
    FlatVector<> Dmn(2*order+1, lh);

    for (int n=1; n <= order; n++)
      {
        HeapReset hr(lh);
          
        FlatMatrix<double> trafo(n+1, 2*n+1, lh); 
        /*
          Recursive Computation of Spherical Harmonic Rotation Coefficients of Large Degree
          Nail A. Gumerov and Ramani Duraiswami
          within Excursions in Harmonic Analysis, Volume 3
            
          page 130 
        */

        // Step 2
        // H(0,m)
        trafo.Col(n) = 1.0/sqrt(2*n+1) * normalized_leg_func.Col(n).Range(n+1);
        for (int m = 1; m <= n; m += 2) trafo(m,n) *= -1;
        // Step 3
        // H(1,m)
        FlatVector<double> tmp = 1.0/sqrt(2*n+3) * normalized_leg_func.Col(n+1).Range(n+2) | lh;
        for (int m = 1; m < tmp.Size(); m += 2) tmp(m) *= -1;
        for (int m = 1; m <= n; m++)
          trafo.Col(n+1)(m) = 1/CalcBmn(0,n+1) * (  CalcBmn(-m-1, n+1)*(1-c)/2 * tmp(m+1)
                                                    - CalcBmn(m-1,n+1)*(1+c)/2 * tmp(m-1)
                                                    - CalcAmn(m,n) * s*tmp(m));

        // Step 4
        // diamond - recursion
        for (int mp = -n; mp <= n; mp++)
          Dmn(order+mp) = CalcDmn(mp, n);

        for (int mp = 1; mp < n; mp++)
          {
            double invDmn = 1.0 / Dmn(order+mp);
            for (int m = mp; m < n; m++)
              trafo(m, n+mp+1) = invDmn  * ( Dmn(order+mp-1) *trafo(m  ,n+mp-1)
                                             -Dmn(order+m-1)*trafo(m-1,n+mp)
                                             +Dmn(order+m)  *trafo(m+1,n+mp));
            int m = n;
            trafo(m, n+mp+1) = invDmn * ( Dmn(order+mp-1,n)*trafo(m  ,n+mp-1)
                                          -Dmn(order+m-1,n)*trafo(m-1,n+mp));
          }
      
        // Step 5
        // diamond - recursion, negative      
        for (int mp = 0; mp > -n; mp--)
          {
            double invDmn = 1.0 / Dmn(order+mp-1);              
            for (int m = -mp+1; m < n; m++)
              trafo(m, n+mp-1) = invDmn * (  Dmn(order+mp,n)*trafo(m  ,n+mp+1)
                                             +Dmn(order+m-1,n)*trafo(m-1,n+mp)
                                             -Dmn(order+m  ,n)*trafo(m+1,n+mp));
            int m = n;
            trafo(m, n+mp-1) = invDmn * (  Dmn(order+mp,n)*trafo(m  ,n+mp+1)
                                           +Dmn(order+m-1,n)*trafo(m-1,n+mp));
          }

        // RegionTimer rgtrafo(*timerstrafo[order]);                    
        // Step 6
        // symmetries in m and mp
        for (int m = 0; m <= n; m++)
          {
            auto dst = trafo.Row(m).Range(n+m, n+n+1);
            auto src = trafo.Col(n+m).Range(m, n+1);
            dst = src;
          }
        for (int m = 0; m <= n; m++)
          {
            auto dst = trafo.Row(m).Range(n-n, n-m+1).Reversed();
            auto src = trafo.Col(n-m).Range(m, n+1);
            dst = src;
          }
        /*
          double errortho = L2Norm( Matrix(trafo*Trans(trafo) - Identity(n+1)));
          if (errortho > 1e-10)
          {
          *testout << "n = " << n << " order = " << Order() << ", alpha = " << alpha << ", errortho = " << errortho << endl;
          if (n < 10)
          *testout << trafo*Trans(trafo) << endl;
          }
        */

        FlatVector<entry_type> cn = CoefsN(n);
        FlatVector<entry_type> old = cn | lh;

        cn.Slice(0,1) = Trans(trafo) * old.Range(n, 2*n+1);
        cn.Slice(0,1).Reversed() += Trans(trafo.Rows(1,n+1)) * old.Range(0,n).Reversed();

        for (int m = 1; m <= n; m+=2)
          {
            cn(n+m) *= -1;
            cn(n-m) *= -1;
          }
      }
  }




    
#ifdef VER1
  template <typename TARGET>
  void ShiftZ (double z, MultiPole<TARGET> & target)
  {
    static Timer t("mptool ShiftZ"+ToString(typeid(RADIAL).name())+ToString(typeid(TARGET).name()));
      
    RegionTimer rg(t);
      
    int os = sh.Order();
    int ot = target.SH().Order();

    if (os > 100 && ot > 100 && abs(z)*kappa > 0.3*min(os,ot) && is_same<RADIAL,TARGET>())
      {
        MultiPole<TARGET> tmp {target};
        ShiftZ(z/2, tmp);
        tmp.ShiftZ(z/2, target);
        return;
      }

      
    target.SH().Coefs()=0.0;

    LocalHeap lh( 16*( (os+ot+1)*(os+1) + (os+1 + ot+1) ) + 8*2*(os+ot+1) + 500);

    FlatMatrix<Complex> trafo(os+ot+1, os+1, lh); 
    FlatVector<Complex> hv1(os+1, lh), hv2(ot+1, lh);
    FlatVector<double> amn(os+ot+1, lh);
    FlatVector<double> inv_amn(os+ot+1, lh);
      
    // trafo = Complex(0.0);


    double tscale = target.Scale();
    double inv_tscale = 1.0/tscale;

      
    // (185) from paper 'fast, exact, stable, Gumerov+Duraiswami
    // RADIAL::Eval(os+ot, kappa*abs(z), trafo.Col(0));
    if (typeid(RADIAL) == typeid(TARGET))
      SphericalBessel (os+ot, kappa*abs(z), tscale, trafo.Col(0));
    else
      SphericalHankel1 (os+ot, kappa*abs(z), inv_tscale, trafo.Col(0));

    /*
      if (L2Norm(trafo.Col(0)) > 1e5 || std::isnan(L2Norm(trafo.Col(0))))
      {
      *testout << "large Hankel: " << L2Norm(trafo.Col(0)) << endl;
      *testout << "kappa z = " << kappa*z << ", os = " << os << ", ot = " << ot << endl;
      }
    */
    // if (L2Norm(trafo.Col(0)) > 1e5)
    // throw Exception ("z-shift - coefs large");
      
    if (z < 0)
      for (int l = 1; l < trafo.Height(); l+=2) trafo(l,0) *= -1;
      
    for (int l = 0; l <= os+ot; l++)
      trafo(l,0) *= sqrt(2*l+1);

    if (os > 0)
      {
        for (int l = 1; l < os+ot; l++)   
          trafo(l,1) = -scale/sh.CalcAmn(0,0) * (sh.CalcAmn(0,l)*tscale*trafo(l+1,0)
                                                 -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,0));
        trafo(0,1) = -scale*tscale*trafo(1,0);
      }
      
    for (int n = 1; n < os; n++)
      {
        for (int l = 1; l < os+ot-n; l++)
          trafo(l,n+1) = -scale/sh.CalcAmn(0,n) * (sh.CalcAmn(0,l)*tscale*trafo(l+1,n)
                                                   -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,n)
                                                   -sh.CalcAmn(0,n-1)*scale*trafo(l,n-1));
        trafo(0,n+1) = pow(-scale*tscale,n+1)*trafo(n+1,0);
      }

    cout << "m = " << 0 << endl
         << trafo.Rows(0,ot+1) << endl;

      
    for (int n = 0; n <= os; n++)
      hv1(n) = sh.Coef(n,0);
    hv2 = trafo.Rows(ot+1) * hv1;
    for (int n = 0; n <= ot; n++)
      target.SH().Coef(n,0) = hv2(n);

      
    for (int m = 1; m <= min(os,ot); m++)
      {
        // fill recursive formula (187)
        for (int l = m; l <= os+ot-m; l++)
          trafo(l,m) = scale/sh.CalcBmn(-m, m) * (sh.CalcBmn(-m, l)*inv_tscale*trafo(l-1, m-1)
                                                  -sh.CalcBmn(m-1,l+1)*tscale*trafo(l+1,m-1));


        /*
          cout << "m = " << m << endl;
          cout << "             norm col0 = " << L2Norm(trafo.Col(m).Range(m,os+ot-m+1)) << endl;
        */
          
        for (int l = m-1; l < os+ot-m; l++)
          {
            amn(l) = sh.CalcAmn(m,l);
            inv_amn(l) = scale/amn(l);
          }

        double prod = 1; 
        for (int n = m; n < os; n++)
          {
            for (int l = m+1; l < os+ot-n; l++)
              trafo(l,n+1) = -inv_amn(n) * (amn(l)*tscale*trafo(l+1,n)
                                            -amn(l-1)*inv_tscale*trafo(l-1,n)
                                            -amn(n-1)*scale*trafo(l,n-1));

            prod *= -scale*tscale;
            trafo(m,n+1) = prod*trafo(n+1,m);              
          }

        /*
          cout << "                            norm trafo = "
          << L2Norm(trafo.Rows(m,ot+1).Cols(m,os+1)) 
          << " ortho " << L2Norm( Trans(trafo.Rows(m,ot+1).Cols(m,os+1))*trafo.Rows(m,ot+1).Cols(m,os+1)
          - Identity(os+1-m)) << endl;
        */
        /*
         *testout << "norm trafo = " << L2Norm(trafo.Rows(m, ot+1).Cols(m,os+1)) << endl;
         if ( L2Norm(trafo.Rows(m, ot+1).Cols(m,os+1)) > 1e30)
         {
         *testout << trafo.Rows(m, ot+1).Cols(m,os+1) << endl;
         for (int i = m; i < os+1; i++)
         {
         *testout << "norm col " << i << " = " << L2Norm(trafo.Col(i).Range(m,os+ot-i)) << endl;
         *testout << "col " << i << " = " << trafo.Col(i).Range(m,os+ot-i) << endl;
         }
         throw Exception("large mat");
         }
        */

        cout << "m = " << m << endl
             << trafo.Rows(m,ot+1).Cols(m,os+1) << endl;
          
        for (int n = m; n <= os; n++)
          hv1(n) = sh.Coef(n,m);
        hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
        for (int n = m; n <= ot; n++)
          target.SH().Coef(n,m) = hv2(n);
          
        for (int n = m; n <= os; n++)
          hv1(n) = sh.Coef(n,-m);
        hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
        for (int n = m; n <= ot; n++)
          target.SH().Coef(n,-m) = hv2(n);
      }
  }
#endif


#ifdef VER2
    
  template <typename TARGET>
  void ShiftZ (double z, MultiPole<TARGET> & target)
  {
    static Timer t("mptool ShiftZ"+ToString(typeid(RADIAL).name())+ToString(typeid(TARGET).name()));
    RegionTimer rg(t);
      
    int os = sh.Order();
    int ot = target.SH().Order();
      
    target.SH().Coefs()=0.0;

    LocalHeap lh( 16*( (os+ot+1)*(os+1) + (os+1 + ot+1) ) + 8*2*(os+ot+1) + 500);

    FlatMatrix<Complex> trafo(os+ot+1, os+1, lh); 
    FlatVector<Complex> hv1(os+1, lh), hv2(ot+1, lh);
    FlatVector<double> amn(os+ot+1, lh);
    FlatVector<double> inv_amn(os+ot+1, lh);
      
    // trafo = Complex(0.0);

    double tscale = target.Scale();
    double inv_tscale = 1.0/tscale;

      
    // (185) from paper 'fast, exact, stable, Gumerov+Duraiswami
    // RADIAL::Eval(os+ot, kappa*abs(z), trafo.Col(0));
    if (typeid(RADIAL) == typeid(TARGET))
      SphericalBessel (os+ot, kappa*abs(z), tscale, trafo.Col(0));
    else
      SphericalHankel1 (os+ot, kappa*abs(z), inv_tscale, trafo.Col(0));
      
    if (z < 0)
      for (int l = 1; l < trafo.Height(); l+=2) trafo(l,0) *= -1;
      
    // for (int l = 0; l <= os+ot; l++)
    // trafo(l,0) *= sqrt(2*l+1);

    if (os > 0)
      {
        for (int l = 1; l < os+ot; l++)
          {
            /*
              trafo(l,1) = -scale/sh.CalcAmn(0,0) *
              (sh.CalcAmn(0,l)*tscale*trafo(l+1,0)
              -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,0));
            */

            int m = 0, n = 0;
            double fac = ((2*l+1.0)/(2*n+1.0));
            trafo(l,n+1) = -scale/ ( sqrt((n+1+m)*(n+1-m)) * fac)   *
              (sqrt( (l+1+m)*(l+1-m)) * tscale * trafo(l+1,n)
               -sqrt( (l+m)*(l-m) )  * inv_tscale * trafo(l-1,n)
               -sqrt( (n+m)*(n-m) ) * fac * scale * trafo(l,n-1));
          }
        trafo(0,1) = -scale*tscale*trafo(1,0);
      }
      
    for (int n = 1; n < os; n++)
      {
        for (int l = 1; l < os+ot-n; l++)
          {
            /*
              trafo(l,n+1) = -scale/sh.CalcAmn(0,n) * (sh.CalcAmn(0,l)*tscale*trafo(l+1,n)
              -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,n)
              -sh.CalcAmn(0,n-1)*scale*trafo(l,n-1));
            */
              
            int m = 0;
            double fac = ((2*l+1.0)/(2*n+1.0));
            trafo(l,n+1) = -scale / ( sqrt((n+1+m)*(n+1-m)) * fac)   *
              (sqrt( (l+1+m)*(l+1-m)) * tscale * trafo(l+1,n)
               -sqrt( (l+m)*(l-m) )  * inv_tscale * trafo(l-1,n)
               -sqrt( (n+m)*(n-m) ) * fac * scale * trafo(l,n-1));
          }
        trafo(0,n+1) = pow(-scale*tscale,n+1)*trafo(n+1,0);
      }

      
    Matrix<Complex> scaledtrafo(os+ot+1, os+1);
    for (int l = 0; l <= os+ot; l++)
      for (int n = 0; n <= os; n++)
        scaledtrafo(l,n) = trafo(l,n) * sqrt( (2*l+1)*(2*n+1) );

    // cout << "m = " << 0 << endl
    // << scaledtrafo.Rows(0,ot+1) << endl;
      
    for (int n = 0; n <= os; n++)
      hv1(n) = sh.Coef(n,0) * sqrt(2*n+1);
    hv2 = trafo.Rows(ot+1) * hv1;
    for (int n = 0; n <= ot; n++)
      target.SH().Coef(n,0) = hv2(n) * sqrt(2*n+1);

      
    for (int m = 1; m <= min(os,ot); m++)
      {
        // fill recursive formula (187)
        for (int l = m; l <= os+ot-m; l++)
          {
            trafo(l,m) = scale/sh.CalcBmn(-m, m) *
              (sh.CalcBmn(-m, l)*inv_tscale * trafo(l-1, m-1) * sqrt( (2*l-1)*(2*m-1) )
               -sh.CalcBmn(m-1,l+1)*tscale * trafo(l+1,m-1) * sqrt( (2*l+3)*(2*m-1)) );
            trafo(l,m) /= sqrt( (2*l+1)*(2*m+1) );
          }

        cout << "m = " << m << endl;
        cout << "             norm col0 = " << L2Norm(trafo.Col(m).Range(m,os+ot-m+1)) << endl;
          
        for (int l = m-1; l < os+ot-m; l++)
          {
            amn(l) = sh.CalcAmn(m,l);
            inv_amn(l) = scale/amn(l);
          }

        double prod = 1; 
        for (int n = m; n < os; n++)
          {
            for (int l = m+1; l < os+ot-n; l++)
              {
                /*
                  trafo(l,n+1) = -inv_amn(n) * (amn(l)*tscale*trafo(l+1,n)
                  -amn(l-1)*inv_tscale*trafo(l-1,n)
                  -amn(n-1)*scale*trafo(l,n-1));
                */
                  
                
                double fac = ((2*l+1.0)/(2*n+1.0));
                trafo(l,n+1) = -scale / ( sqrt((n+1+m)*(n+1-m)) * fac)   *
                  (sqrt( (l+1+m)*(l+1-m)) * tscale * trafo(l+1,n)
                   -sqrt( (l+m)*(l-m) )  * inv_tscale * trafo(l-1,n)
                   -sqrt( (n+m)*(n-m) ) * fac * scale * trafo(l,n-1));
              }
            prod *= -scale*tscale;
            trafo(m,n+1) = prod*trafo(n+1,m);              
          }

        double normleft = 0;
        for (int l = m; l <= ot; l++)
          for (int n = m; n <= min(l,os); n++)
            normleft += sqr(abs(trafo(l,n)));
        normleft = sqrt(normleft);
        cout << "                            norm trafo = "
             << L2Norm(trafo.Rows(m,ot+1).Cols(m,os+1)) << ",  normleft = " << normleft << endl;
          
        // << " ortho " << L2Norm( Trans(trafo.Rows(m,ot+1).Cols(m,os+1))*trafo.Rows(m,ot+1).Cols(m,os+1)
        // - Identity(os+1-m)) << endl;
        /*
         *testout << "norm trafo = " << L2Norm(trafo.Rows(m, ot+1).Cols(m,os+1)) << endl;
         if ( L2Norm(trafo.Rows(m, ot+1).Cols(m,os+1)) > 1e30)
         {
         *testout << trafo.Rows(m, ot+1).Cols(m,os+1) << endl;
         for (int i = m; i < os+1; i++)
         {
         *testout << "norm col " << i << " = " << L2Norm(trafo.Col(i).Range(m,os+ot-i)) << endl;
         *testout << "col " << i << " = " << trafo.Col(i).Range(m,os+ot-i) << endl;
         }
         throw Exception("large mat");
         }
        */


        /*
          Matrix<Complex> scaledtrafo(os+ot+1, os+1);
          for (int l = 0; l <= os+ot; l++)
          for (int n = 0; n <= os; n++)
          scaledtrafo(l,n) = trafo(l,n) * sqrt( (2*l+1)*(2*n+1) );
          
          cout << "m = " << m << endl
          << scaledtrafo.Rows(m,ot+1).Cols(m,os+1) << endl;
        */
          
        for (int n = m; n <= os; n++)
          hv1(n) = sh.Coef(n,m) * sqrt(2*n+1);
        hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
        for (int n = m; n <= ot; n++)
          target.SH().Coef(n,m) = hv2(n)*sqrt(2*n+1);
          
        for (int n = m; n <= os; n++)
          hv1(n) = sh.Coef(n,-m) * sqrt(2*n+1);
        hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
        for (int n = m; n <= ot; n++)
          target.SH().Coef(n,-m) = hv2(n) * sqrt(2*n+1);
      }
  }
#endif
  
  /*
    for S->S:  targetscale <= sourcescale <= 1
    for S->R:  targetscale >= 1 >= sourcescale
    for R->R:  targetscale >= sourcescale >= 1

    trafoscaled(k,l) = ss^l / ts^k * trafo(k,l)

    trafo(k,l) = trafo(l,k)  <==>   trafoscale(k,l) = (ss*ts)^(l-k) trafoscale(l,k)
    
    S->S : Col0 ok, Row0 may become very small
   */


#ifdef VER3
  template <typename RADIAL, typename entry_type> template <typename TARGET>
  void MultiPole<RADIAL,entry_type> :: ShiftZ (double z, MultiPole<TARGET,entry_type> & target)
  {
    // static Timer t("mptool ShiftZ"+ToString(typeid(RADIAL).name())+ToString(typeid(TARGET).name()));
    // RegionTimer rg(t);
      
    int os = sh.Order();
    int ot = target.SH().Order();
      
    target.SH().Coefs()=0.0;
    // *testout << "siftz " << typeid(RADIAL).name() << " -> " << typeid(TARGET).name() << endl;
    // *testout << "shiftz, norm input = " << L2Norm(SH().Coefs()) << endl;
    
    LocalHeap lh(( 32*( (os+ot+1)*(os+ot+1) + (os+1 + ot+1) ) + 8*3*(os+ot+1) + 500));

    typedef typename std::conditional<std::is_same<RADIAL,TARGET>::value, double, Complex>::type trafo_type;
    
    FlatMatrix<trafo_type> trafo(os+ot+1, max(os,ot)+1, lh);
    FlatMatrix<trafo_type> oldtrafo(os+ot+1, max(os,ot)+1, lh);     
    FlatVector<entry_type> hv1(os+1, lh), hv2(ot+1, lh);

    trafo = trafo_type(0.0);

    double scale = Scale();
    // double tscale = target.Scale();
    double tscale = target.Scale();
    double inv_tscale = 1.0/tscale;

    /*
    *testout << "SiftZ" +ToString(typeid(RADIAL).name())+ToString(typeid(TARGET).name()) << endl;
    *testout << "z = " << z << endl;
    *testout << "os = " << os << ", ot = " <<  ot << endl;
    *testout << "scale = " << Scale() << ", tscale = " << target.Scale() << endl;
    */

    
    FlatVector<double> amn(os+ot+1, lh);
    FlatVector<double> inv_amn(os+ot+1, lh);

    FlatVector<double> powscale(os+1, lh);
    double prod = 1;
    for (int i = 0; i <= os; i++)
      {
        powscale(i) = prod;
        prod *= -scale*tscale;
      }
    
    // (185) from paper 'fast, exact, stable, Gumerov+Duraiswami
    if constexpr (std::is_same<RADIAL,TARGET>::value)
      SphericalBessel (os+ot, kappa*abs(z), tscale, trafo.Col(0));
    else
      SphericalHankel1 (os+ot, kappa*abs(z), inv_tscale, trafo.Col(0));

    if (z < 0)
      for (int l = 1; l < trafo.Height(); l+=2) trafo(l,0) *= -1;

    for (int l = 0; l <= os+ot; l++)
      trafo(l,0) *= sqrt(2*l+1);


    // *testout << "Norm Col0 = " << L2Norm(trafo.Col(0)) << endl;
    
    if (os > 0)
      {
        for (int l = 1; l < os+ot; l++)   
          trafo(l,1) = -scale/sh.CalcAmn(0,0) * (sh.CalcAmn(0,l)*tscale*trafo(l+1,0)
                                                 -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,0));
        trafo(0,1) = -scale*tscale*trafo(1,0);
      }

    for (int n = 1; n < trafo.Width()-1; n++)
      {
        for (int l = 1; l < os+ot-n; l++)
          trafo(l,n+1) = -scale/sh.CalcAmn(0,n) * (sh.CalcAmn(0,l)*tscale*trafo(l+1,n)
                                                   -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,n)
                                                   -sh.CalcAmn(0,n-1)*scale*trafo(l,n-1));
        trafo(0,n+1) = powscale(n+1)*trafo(n+1,0);
      }

    // use symmetry of matrix (up to scaling)
    // important for small kappa, and R->R
    for (int n = 0; n < trafo.Width(); n++)
      for (int l = n+1; l < trafo.Width(); l++)
        trafo(n,l) = powscale(l-n) * trafo(l,n);

    /*
    if constexpr (std::is_same<RADIAL,MPRegular>::value && std::is_same<TARGET,MPRegular>::value)
      {
        *testout << "RR, old" << endl;
        *testout << "trafo0 = " << trafo.Rows(ot+1).Cols(os+1) << endl;
      }
    */

    
    // *testout << "trafo0" << endl << trafo.Rows(os+ot+1).Cols(os+1) << endl;
    // *testout << "norm trafo col 0 = " << L2Norm(trafo.Rows(ot+1).Col(0)) << endl;
    // *testout << "norm trafo0 = " << L2Norm(trafo.Rows(ot+1).Cols(os+1)) << endl;
    // *testout << "trafo0 = " << trafo.Rows(ot+1).Cols(os+1) << endl;
    
    for (int n = 0; n <= os; n++)
      hv1(n) = sh.Coef(n,0);
    hv2 = trafo.Rows(ot+1).Cols(os+1) * hv1;
    for (int n = 0; n <= ot; n++)
      target.SH().Coef(n,0) = hv2(n);


      
    for (int m = 1; m <= min(os,ot); m++)
      {

        for (int l = m-1; l < os+ot-m; l++)
          {
            amn(l) = sh.CalcAmn(m,l);
            inv_amn(l) = scale/amn(l);
          }

        
        if constexpr (!std::is_same<RADIAL,TARGET>::value)
          {

            trafo.Swap (oldtrafo);
            trafo = trafo_type(0.0);

            // fill recursive formula (187)
            for (int l = m; l <= os+ot-m; l++)
              trafo(l,m) = scale/sh.CalcBmn(-m, m) * (sh.CalcBmn(-m, l)*inv_tscale*oldtrafo(l-1, m-1)  
                                                      -sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,m-1));  
            
            // unstable for large kappa
            for (int n = m; n < os; n++)
              {
                for (int l = n+1; l < os+ot-n; l++)
                  trafo(l,n+1) = -inv_amn(n) * (amn(l)*tscale*trafo(l+1,n)
                                                -amn(l-1)*inv_tscale*trafo(l-1,n)
                                                -amn(n-1)*scale*trafo(l,n-1));
              }


            
            // fix real part by wall-reccurence:
            for (int n = m; n < trafo.Width()-1; n++)
              {
                // int l = 2*order-n-1;
                int l = trafo.Height()-n-2;
                // *testout << "old: " << trafo(l,n+1) << ", new: ";
                trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n) * scale*trafo(l,n-1)
                                                          - sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,n)     
                                                          + sh.CalcBmn(-m,l)  * 1/tscale*oldtrafo(l-1,n) );
                // *testout << trafo(l,n+1) << endl;

                // *testout << "old2: " << trafo(l-1,n) << ", new: ";                
                trafo(l-1,n) = tscale/amn(l-1) * (amn(l)  *   tscale*trafo(l+1,n)
                                                  - amn(n-1)* scale*trafo(l,n-1)
                                                  + amn(n)* 1/scale*trafo(l,n+1));
                // *testout << trafo(l-1,n) << endl;                
              }

            // the same thing 1 row up
            for (int n = m; n < trafo.Width()-2; n++)
              {
                // int l = 2*order-n-2;
                int l = trafo.Height()-n-3;
                  
                trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n)     * scale*trafo(l,n-1)
                                                          - sh.CalcBmn(m-1,l+1) * tscale* oldtrafo(l+1,n)   
                                                          + sh.CalcBmn(-m,l)    * 1/tscale* oldtrafo(l-1,n) ); 
                  
                trafo(l-1,n) = tscale/amn(l-1) * (amn(l)   * tscale*trafo(l+1,n)
                                                  -amn(n-1)* scale*trafo(l,n-1) +
                                                  amn(n)   * 1/scale*trafo(l,n+1)) ;
              }

            for (int l = trafo.Height()-1; l >= m; l--)
              for (int n = m+1; n < min<int>(trafo.Height()-1-l,l); n++)
                {
                  trafo(l-1,n)  = Complex (
                                           (tscale/amn(l-1)* ( amn(l)  * tscale*trafo(l+1,n)
                                                               -amn(n-1)* scale*trafo(l,n-1)
                                                               +amn(n)  * 1/scale*trafo(l,n+1))).real(),
                                           trafo(l-1,n).imag());
                }




            // the imaginary part started from the yn:
            // diagonal down
            for (int n = m; n < trafo.Width()-1; n++)
              {
                int l = n+1;
                // *testout << "old: " << trafo(l,n+1) << ", new: ";
                trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n) * scale*trafo(l,n-1)
                                                          - sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,n)     
                                                          + sh.CalcBmn(-m,l)  * 1/tscale*oldtrafo(l-1,n) );
                // *testout << trafo(l,n+1) << endl;
                // *testout << "old2: " << trafo(l-1,n) << ", new: ";

                trafo(l+1,n) = 
                  1.0 / (amn(l)  *   tscale) * (amn(l-1)/tscale * trafo(l-1,n)
                                                + amn(n-1)* scale*trafo(l,n-1)
                                                - amn(n)* 1/scale*trafo(l,n+1));
                // *testout << trafo(l-1,n) << endl;
              }

            // the next diagonal down
            for (int n = m; n < trafo.Width()-1; n++)
              {
                // int l = 2*order-n-1;
               int l = n+2;
               // *testout << "old: " << trafo(l,n+1) << ", new: ";
               trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n) * scale*trafo(l,n-1)
                                                         - sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,n)     
                                                         + sh.CalcBmn(-m,l)  * 1/tscale*oldtrafo(l-1,n) );
               // *testout << trafo(l,n+1) << endl;
               // *testout << "old2: " << trafo(l-1,n) << ", new: ";

                trafo(l+1,n) = 
                  1.0 / (amn(l)  *   tscale) * (amn(l-1)/tscale * trafo(l-1,n)
                                                + amn(n-1)* scale*trafo(l,n-1)
                                                - amn(n)* 1/scale*trafo(l,n+1));
                // *testout << trafo(l-1,n) << endl;
              }
            
              
            for (int l = m; l < trafo.Height()-1; l++)
              for (int n = m+1; n < min<int>(trafo.Height()-1-l,l-1 ); n++)
                {
                  trafo(l+1,n)  = Complex (trafo(l+1,n).real(),
                                           1/(tscale*amn(l)) * ( amn(l-1)/tscale*trafo(l-1,n)
                                                                 +amn(n-1)* scale*trafo(l,n-1)
                                                                 -amn(n)  * 1/scale*trafo(l,n+1)).imag());
                }

            // *testout << "m = " << m << endl;
            // Matrix<> mat = Real(trafo.Rows(m, ot+1).Cols(m,os+1));
            // *testout << "Norm real = " << L2Norm2(mat) << endl;
            // mat = Imag(trafo.Rows(m, ot+1).Cols(m,os+1)); 
            // *testout << "Norm imag = " << L2Norm2(mat) << endl;
          }

        else

          {
            /*
              for (int n = m; n < trafo.Width()-1; n++)
              {
              for (int l = n+1; l < os+ot-n; l++)
              trafo(l,n+1) = -inv_amn(n) * (amn(l)*tscale*trafo(l+1,n)
              -amn(l-1)*inv_tscale*trafo(l-1,n)
              -amn(n-1)*scale*trafo(l,n-1));
              }
            */

            trafo.Swap (oldtrafo);
            trafo = trafo_type(0.0);
              
            // fill recursive formula (187)
            for (int l = m; l <= os+ot-m; l++)
              trafo(l,m) = scale/sh.CalcBmn(-m, m) * (sh.CalcBmn(-m, l)*inv_tscale*oldtrafo(l-1, m-1)  
                                                      -sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,m-1));  
              
            for (int n = m; n < trafo.Width()-1; n++)
              {
                // int l = 2*order-n-1;
                int l = trafo.Height()-n-2;
                  
                trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n) * scale*trafo(l,n-1)
                                                          - sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,n)     
                                                          + sh.CalcBmn(-m,l)  * 1/tscale*oldtrafo(l-1,n) );
                  
                trafo(l-1,n) = tscale/amn(l-1) * (amn(l)  *   tscale*trafo(l+1,n)
                                                  - amn(n-1)* scale*trafo(l,n-1)
                                                  + amn(n)* 1/scale*trafo(l,n+1));
              }
              
            // the same thing 1 row up
            for (int n = m; n < trafo.Width()-2; n++)
              {
                // int l = 2*order-n-2;
                int l = trafo.Height()-n-3;
                  
                trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n)     * scale*trafo(l,n-1)
                                                          - sh.CalcBmn(m-1,l+1) * tscale* oldtrafo(l+1,n)   
                                                          + sh.CalcBmn(-m,l)    * 1/tscale* oldtrafo(l-1,n) ); 
                  
                trafo(l-1,n) = tscale/amn(l-1) * (amn(l)   * tscale*trafo(l+1,n)
                                                  -amn(n-1)* scale*trafo(l,n-1) +
                                                  amn(n)   * 1/scale*trafo(l,n+1)) ;
              }
              
              
            for (int l = trafo.Height()-1; l >= m; l--)
              for (int n = m+1; n < min<int>(trafo.Height()-1-l,l); n++)
                {
                  trafo(l-1,n) = tscale/amn(l-1)* ( amn(l)  * tscale*trafo(l+1,n)
                                                    -amn(n-1)* scale*trafo(l,n-1)
                                                    +amn(n)  * 1/scale*trafo(l,n+1)) ;
                }
          }
          
          
        // *testout << "m = " << m << endl
        // << trafo << endl;

        for (int n = m; n < os; n++)
          for (int l = n+1; l <= os; l++)
            trafo(n,l) = powscale(l-n) * trafo(l,n);              
        // trafo(n,l) = pow(-scale*tscale, l-n) * trafo(l,n);

        /*
        if (m == 1)
          if constexpr (std::is_same<RADIAL,MPRegular>::value && std::is_same<TARGET,MPRegular>::value)
            {
              *testout << "RR, old" << endl;
              *testout << "trafo1 = " << trafo.Rows(ot+1).Cols(os+ot+1-m) << endl;
            }
        */

        
        // *testout << "norm trafo " << m << " = " 
        // << L2Norm(trafo.Rows(m,ot+1).Cols(m,os+1)) << endl;
          
        for (int n = m; n <= os; n++)
          hv1(n) = sh.Coef(n,m);
        hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
        for (int n = m; n <= ot; n++)
          target.SH().Coef(n,m) = hv2(n);
          
        for (int n = m; n <= os; n++)
          hv1(n) = sh.Coef(n,-m);
        hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
        for (int n = m; n <= ot; n++)
          target.SH().Coef(n,-m) = hv2(n);
      }




    
    if constexpr (std::is_same<RADIAL,MPRegular>::value && std::is_same<TARGET,MPRegular>::value)
      {
        target.SH().Coefs() = 0.0;
        
        int os = sh.Order();
        int ot = target.SH().Order();
        
        target.SH().Coefs()=0.0;
    
        LocalHeap lh(( 32*( (os+ot+1)*(os+ot+1) + (os+1 + ot+1) ) + 8*3*(os+ot+1) + 500));
        
        typedef double trafo_type;
        
        FlatMatrix<trafo_type> trafo(os+ot+1, os+ot+1, lh);
        FlatMatrix<trafo_type> oldtrafo(os+ot+1, os+ot+1, lh);     
        FlatVector<entry_type> hv1(os+1, lh), hv2(ot+1, lh);
        
        trafo = trafo_type(0.0);
        
        double scale = Scale();
        double tscale = target.Scale();
        double inv_tscale = 1.0/tscale;
        
    
        FlatVector<double> amn(os+ot+1, lh);
        FlatVector<double> inv_amn(os+ot+1, lh);
        
        FlatVector<double> powscale(os+1, lh);
        double prod = 1;
        for (int i = 0; i <= os; i++)
          {
            powscale(i) = prod;
            prod /= -scale*tscale;
          }


        
        // (185) from paper 'fast, exact, stable, Gumerov+Duraiswami
        SphericalBessel (os+ot, kappa*abs(z), 1/scale, trafo.Row(0));
        
        if (z > 0)
          for (int l = 1; l < trafo.Width(); l+=2) trafo(0,l) *= -1;
        
        for (int l = 0; l <= os+ot; l++)
          trafo(0,l) *= sqrt(2*l+1);
        

        // *testout << "Norm Row0 = " << L2Norm(trafo.Row(0)) << endl;

        if (ot > 0)
          {
            for (int l = 1; l < os+ot; l++)   
              trafo(1,l) = -1/tscale/sh.CalcAmn(0,0) * (sh.CalcAmn(0,l)/scale*trafo(0,l+1)
                                                        -sh.CalcAmn(0,l-1)*scale*trafo(0,l-1));
            trafo(1,0) = 1.0/(-scale*tscale)*trafo(0,1);
          }
        
        for (int n = 1; n < trafo.Height()-1; n++)
          {
            for (int l = 1; l < os+ot-n; l++)
              trafo(n+1,l) = -1.0/tscale/sh.CalcAmn(0,n) * (sh.CalcAmn(0,l)*1/scale*trafo(n,l+1)
                                                            -sh.CalcAmn(0,l-1)*scale*trafo(n,l-1)
                                                            -sh.CalcAmn(0,n-1)*1/tscale*trafo(n-1,l));
            trafo(n+1,0) = powscale(n+1)*trafo(0,n+1);
          }
        
        // use symmetry of matrix (up to scaling)
        // important for small kappa, and R->R
        for (int n = 0; n < trafo.Height(); n++)
          for (int l = n+1; l < trafo.Height(); l++)
            trafo(l,n) = powscale(l-n) * trafo(n,l);

        
        for (int n = 0; n <= os; n++)
          hv1(n) = sh.Coef(n,0);
        hv2 = trafo.Rows(ot+1).Cols(os+1) * hv1;
        for (int n = 0; n <= ot; n++)
          target.SH().Coef(n,0) = hv2(n);
        // *testout << "trafo0 = " << endl <<  trafo.Rows(ot+1).Cols(os+1) << endl;


        for (int m = 1; m <= min(os,ot); m++)
          {
            
            for (int l = m-1; l < os+ot-m; l++)
              {
                amn(l) = sh.CalcAmn(m,l);
                // inv_amn(l) = scale/amn(l);
              }
            
            trafo.Swap (oldtrafo);
            trafo = trafo_type(0.0);
            
            // fill recursive formula (187)
            for (int l = m; l <= os+ot-m; l++)
              trafo(m,l) = 1/tscale/sh.CalcBmn(-m, m) * (sh.CalcBmn(-m, l)*scale*oldtrafo(m-1,l-1)  
                                                         -sh.CalcBmn(m-1,l+1)*1/scale*oldtrafo(m-1,l+1));  
            

            for (int n = m; n < trafo.Height()-1; n++)
              {
                int l = trafo.Width()-n-2;
                  
                trafo(n+1,l) = 1/tscale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n) * 1/tscale*trafo(n-1,l)
                                                             - sh.CalcBmn(m-1,l+1)*1/scale*oldtrafo(n,l+1)     
                                                             + sh.CalcBmn(-m,l)  * scale*oldtrafo(n,l-1) );
                
                trafo(n,l-1) = 1/scale/amn(l-1) * (amn(l)  *   1/scale*trafo(n,l+1)
                                                   - amn(n-1)* 1/tscale*trafo(n-1,l)
                                                   + amn(n)* tscale*trafo(n+1,l));
              }
            
            // the same thing 1 row up
            for (int n = m; n < trafo.Height()-2; n++)
              {
                // int l = 2*order-n-2;
                int l = trafo.Width()-n-3;
                  
                trafo(n+1,l) = 1/tscale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n)     * 1/tscale*trafo(n-1,l)
                                                             - sh.CalcBmn(m-1,l+1) * 1/scale* oldtrafo(n,l+1)   
                                                             + sh.CalcBmn(-m,l)    * scale* oldtrafo(n,l-1) ); 
                  
                trafo(n,l-1) = 1/scale/amn(l-1) * (amn(l)   * 1/scale*trafo(n,l+1)
                                                   -amn(n-1)* 1/tscale*trafo(n-1,l) +
                                                   amn(n)   * tscale*trafo(n+1,l)) ;
              }
              

            for (int l = trafo.Width()-1; l >= m; l--)
              for (int n = m+1; n < min<int>(trafo.Width()-1-l,l); n++)
                {
                  trafo(n,l-1) = 1/scale/amn(l-1)* ( amn(l)  * 1/scale*trafo(n,l+1)
                                                     -amn(n-1)* 1/tscale*trafo(n-1,l)
                                                     +amn(n)  * tscale*trafo(n+1,l)) ;
                }

            
            
            // *testout << "m = " << m << endl
            // << trafo << endl;

            for (int n = m; n < os; n++)
              for (int l = n+1; l <= os; l++)
                trafo(l,n) = powscale(l-n) * trafo(n,l);

            // trafo(n,l) = pow(-scale*tscale, l-n) * trafo(l,n);

            /*
            if (m == 1)
              if constexpr (std::is_same<RADIAL,MPRegular>::value && std::is_same<TARGET,MPRegular>::value)
                {
                  *testout << "RR, old" << endl;
                  *testout << "trafo1 = " << trafo.Rows(ot+1).Cols(os+ot+1-m) << endl;
                }
            */

            // *testout << "norm trafo " << m << " = " 
            // << L2Norm(trafo.Rows(m,ot+1).Cols(m,os+1)) << endl;
            
            for (int n = m; n <= os; n++)
              hv1(n) = sh.Coef(n,m);
            hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
            for (int n = m; n <= ot; n++)
              target.SH().Coef(n,m) = hv2(n);
            
            for (int n = m; n <= os; n++)
              hv1(n) = sh.Coef(n,-m);
            hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
            for (int n = m; n <= ot; n++)
              target.SH().Coef(n,-m) = hv2(n);
          }
      }
        
    
  }
#endif



#define VER4
#ifdef VER4
  template <typename RADIAL, typename entry_type> template <typename TARGET>
  void MultiPole<RADIAL,entry_type> :: ShiftZ (double z, MultiPole<TARGET,entry_type> & target)
  {
    // static Timer t("mptool ShiftZ"+ToString(typeid(RADIAL).name())+ToString(typeid(TARGET).name()));
    // RegionTimer rg(t);
      
    int os = sh.Order();
    int ot = target.SH().Order();

    double scale = Scale();
    double inv_scale = 1.0/scale;
    double tscale = target.Scale();
    double inv_tscale = 1.0/tscale;

    target.SH().Coefs()=0.0;
    // *testout << "siftz " << typeid(RADIAL).name() << " -> " << typeid(TARGET).name() << endl;
    // *testout << "shiftz, norm input = " << L2Norm(SH().Coefs()) << endl;

    if constexpr (std::is_same<RADIAL,MPSingular>::value && std::is_same<TARGET,MPRegular>::value)
      {
        // sing->reg trafo
        // static int cnt = 0;
        // cnt++;
        // *testout << "sing-reg, cnt = " << cnt << endl;
        // *testout << "norm in = " << L2Norm(SH().Coefs()) << endl;
        LocalHeap lh(( 32*( (os+ot+1)*(os+ot+1) + (os+1 + ot+1) ) + 8*3*(os+ot+1) + 500));
        
        typedef Complex trafo_type;
        
        FlatMatrix<trafo_type> trafo(os+ot+1, max(os,ot)+1, lh);
        FlatMatrix<trafo_type> oldtrafo(os+ot+1, max(os,ot)+1, lh);     
        FlatVector<entry_type> hv1(os+1, lh), hv2(ot+1, lh);
        
        trafo = trafo_type(0.0);
        
        
        /*
         *testout << "SiftZ" +ToString(typeid(RADIAL).name())+ToString(typeid(TARGET).name()) << endl;
         *testout << "z = " << z << endl;
         *testout << "os = " << os << ", ot = " <<  ot << endl;
         *testout << "scale = " << Scale() << ", tscale = " << target.Scale() << endl;
         */
        
        
        FlatVector<double> amn(os+ot+1, lh);
        FlatVector<double> inv_amn(os+ot+1, lh);
        
        FlatVector<double> powscale(os+ot+1, lh);
        double prod = 1;
        for (int i = 0; i <= os; i++)
          {
            powscale(i) = prod;
            prod *= -scale*tscale;
          }
        
        // (185) from paper 'fast, exact, stable, Gumerov+Duraiswami
        // if constexpr (std::is_same<RADIAL,TARGET>::value)
        // SphericalBessel (os+ot, kappa*abs(z), tscale, trafo.Col(0));
        // else
        SphericalHankel1 (os+ot, kappa*abs(z), inv_tscale, trafo.Col(0));
        
        if (z < 0)
          for (int l = 1; l < trafo.Height(); l+=2) trafo(l,0) *= -1;
        
        for (int l = 0; l <= os+ot; l++)
          trafo(l,0) *= sqrt(2*l+1);
        
        
        // *testout << "Norm Col0 = " << L2Norm(trafo.Col(0)) << endl;
        
        if (os > 0)
          {
            for (int l = 1; l < os+ot; l++)   
              trafo(l,1) = -scale/sh.CalcAmn(0,0) * (sh.CalcAmn(0,l)*tscale*trafo(l+1,0)
                                                     -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,0));
            trafo(0,1) = -scale*tscale*trafo(1,0);
          }
        
        for (int n = 1; n < trafo.Width()-1; n++)
          {
            for (int l = 1; l < os+ot-n; l++)
              trafo(l,n+1) = -scale/sh.CalcAmn(0,n) * (sh.CalcAmn(0,l)*tscale*trafo(l+1,n)
                                                       -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,n)
                                                       -sh.CalcAmn(0,n-1)*scale*trafo(l,n-1));
            trafo(0,n+1) = powscale(n+1)*trafo(n+1,0);
          }
        
        // use symmetry of matrix (up to scaling)
        // important for small kappa, and R->R
        for (int n = 0; n < trafo.Width(); n++)
          for (int l = n+1; l < trafo.Width(); l++)
            trafo(n,l) = powscale(l-n) * trafo(l,n);
        
        /*
          if constexpr (std::is_same<RADIAL,MPRegular>::value && std::is_same<TARGET,MPRegular>::value)
          {
          *testout << "RR, old" << endl;
          *testout << "trafo0 = " << trafo.Rows(ot+1).Cols(os+1) << endl;
          }
        */
        
    
        // *testout << "trafo0" << endl << trafo.Rows(os+ot+1).Cols(os+1) << endl;
        // *testout << "norm trafo col 0 = " << L2Norm(trafo.Rows(ot+1).Col(0)) << endl;
        // *testout << "norm trafo0 = " << L2Norm(trafo.Rows(ot+1).Cols(os+1)) << endl;
        // *testout << "norm trafo0,real = " << L2Norm(Real(trafo.Rows(ot+1).Cols(os+1))) << endl;    
        // *testout << "norm trafo0,imag = " << L2Norm(Imag(trafo.Rows(ot+1).Cols(os+1))) << endl;    


        // *testout << "trafo0,real = " << Real(trafo.Rows(ot+1).Cols(os+1)) << endl;    
        // *testout << "trafo0,imag = " << Imag(trafo.Rows(ot+1).Cols(os+1)) << endl;    
// *testout << "trafo0 = " << trafo.Rows(ot+1).Cols(os+1) << endl;
        
        for (int n = 0; n <= os; n++)
          hv1(n) = sh.Coef(n,0);
        hv2 = trafo.Rows(ot+1).Cols(os+1) * hv1;
        for (int n = 0; n <= ot; n++)
          target.SH().Coef(n,0) = hv2(n);
        
        for (int m = 1; m <= min(os,ot); m++)
          {
            for (int l = m-1; l < os+ot-m; l++)
              {
                amn(l) = sh.CalcAmn(m,l);
                inv_amn(l) = scale/amn(l);
              }
            
            
            trafo.Swap (oldtrafo);
            trafo = trafo_type(0.0);  // needs initialization ???
            trafo.Col(m-1) = trafo_type(0.0);  // for wall initial value
            
            // fill recursive formula (187)
            for (int l = m; l <= os+ot-m; l++)
              {
                trafo(l,m) = scale/sh.CalcBmn(-m, m) * (sh.CalcBmn(-m, l)*inv_tscale*oldtrafo(l-1, m-1)  
                                                        -sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,m-1));
              }

            /*
            // unstable for large kappa
            for (int n = m; n < os; n++)
              {
                for (int l = n+1; l < os+ot-n; l++)
                  trafo(l,n+1) = -inv_amn(n) * (amn(l)*tscale*trafo(l+1,n)
                                                -amn(l-1)*inv_tscale*trafo(l-1,n)
                                                -amn(n-1)*scale*trafo(l,n-1));
              }
            */
            
            // fix real part by wall-reccurence:
            for (int n = m; n < trafo.Width()-1; n++)
              {
                int l = trafo.Height()-n-2;
                trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n) * scale*trafo(l,n-1)
                                                          - sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,n)     
                                                          + sh.CalcBmn(-m,l)  * 1/tscale*oldtrafo(l-1,n) );
                
                if (n > m)
                trafo(l-1,n) = tscale/amn(l-1) * (amn(l)  *   tscale*trafo(l+1,n)
                                                  - amn(n-1)* scale*trafo(l,n-1)
                                                  + amn(n)* inv_scale*trafo(l,n+1));
              }

            
            // the same thing 1 row up
            for (int n = m; n < trafo.Width()-2; n++)
              {
                // int l = 2*order-n-2;
                int l = trafo.Height()-n-3;
                
                trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n)     * scale*trafo(l,n-1)
                                                          - sh.CalcBmn(m-1,l+1) * tscale* oldtrafo(l+1,n)   
                                                          + sh.CalcBmn(-m,l)    * 1/tscale* oldtrafo(l-1,n) ); 

                if (n > m)                
                trafo(l-1,n) = tscale/amn(l-1) * (amn(l)   * tscale*trafo(l+1,n)
                                                  -amn(n-1)* scale*trafo(l,n-1) +
                                                  amn(n)   * 1/scale*trafo(l,n+1)) ;
              }
                
            for (int l = trafo.Height()-1; l >= m; l--)
              for (int n = m+1; n < min<int>(trafo.Height()-1-l,l); n++)
                {
                  trafo(l-1,n)  = Complex (
                                           (tscale/amn(l-1)* ( amn(l)  * tscale*trafo(l+1,n)
                                                               -amn(n-1)* scale*trafo(l,n-1)
                                                               +amn(n)  * 1/scale*trafo(l,n+1))).real(),
                                           trafo(l-1,n).imag());
                }
            
            
            // the imaginary part started from the yn:
            // diagonal down
            for (int n = m; n < trafo.Width()-1; n++)
              {
                int l = n+1;
                trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n) * scale*trafo(l,n-1)
                                                          - sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,n)     
                                                          + sh.CalcBmn(-m,l)  * 1/tscale*oldtrafo(l-1,n) );

                if (n > m)
                trafo(l+1,n) = 
                  1.0 / (amn(l)  *   tscale) * (amn(l-1)/tscale * trafo(l-1,n)
                                                + amn(n-1)* scale*trafo(l,n-1)
                                                - amn(n)* 1/scale*trafo(l,n+1));
              }
            
            // the next diagonal down
            for (int n = m; n < trafo.Width()-1; n++)
              {
                int l = n+2;
                trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n) * scale*trafo(l,n-1)
                                                          - sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,n)     
                                                          + sh.CalcBmn(-m,l)  * 1/tscale*oldtrafo(l-1,n) );

                if (n > m)                
                trafo(l+1,n) = 
                  1.0 / (amn(l)  *   tscale) * (amn(l-1)/tscale * trafo(l-1,n)
                                                + amn(n-1)* scale*trafo(l,n-1)
                                                - amn(n)* 1/scale*trafo(l,n+1));
              }
            
            
            for (int l = m; l < trafo.Height()-1; l++)
              for (int n = m+1; n < min<int>(trafo.Height()-1-l,l-1 ); n++)
                {
                  trafo(l+1,n)  = Complex (trafo(l+1,n).real(),
                                           1/(tscale*amn(l)) * ( amn(l-1)/tscale*trafo(l-1,n)
                                                                 +amn(n-1)* scale*trafo(l,n-1)
                                                                 -amn(n)  * 1/scale*trafo(l,n+1)).imag());
                }
            
            for (int n = m; n < os; n++)
              for (int l = n+1; l <= os; l++)
                trafo(n,l) = powscale(l-n) * trafo(l,n);              
              
            for (int n = m; n <= os; n++)
              hv1(n) = sh.Coef(n,m);
            hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
            for (int n = m; n <= ot; n++)
              target.SH().Coef(n,m) = hv2(n);
            
            for (int n = m; n <= os; n++)
              hv1(n) = sh.Coef(n,-m);
            hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
            for (int n = m; n <= ot; n++)
              target.SH().Coef(n,-m) = hv2(n);
          }
      }

    

    

    if constexpr (std::is_same<RADIAL,MPSingular>::value && std::is_same<TARGET,MPSingular>::value)
      {
        // sing->sing transformation
        /*
        static Timer t("mptool ShiftZ sing->sing"); RegionTimer rg(t);
        static Timer tfill("mptool ShiftZ sing->sing fill");
        static Timer tmult("mptool ShiftZ sing->sing mult");
        static Timer ttrans("mptool ShiftZ sing->sing trans");        
        static Timer twall("mptool ShiftZ sing->sing wall");
        static Timer tinit("mptool ShiftZ sing->sing init");
        static Timer tlevel0("mptool ShiftZ sing->sing level0");
        static Timer tcleanup("mptool ShiftZ sing->sing cleanup");                        
        */
        // tlevel0.Start();
        LocalHeap lh(( 32*( (os+ot+1)*(os+ot+1) + (os+1 + ot+1) ) + 8*3*(os+ot+1) + 500));

        typedef double trafo_type;
        
        FlatMatrix<trafo_type> trafo(os+ot+1, max(os,ot)+1, lh);
        FlatMatrix<trafo_type> oldtrafo(os+ot+1, max(os,ot)+1, lh);     
        FlatVector<entry_type> hv1(os+1, lh), hv2(ot+1, lh);
        
        trafo = trafo_type(0.0);
        
        FlatVector<double> amn(os+ot+1, lh);
        FlatVector<double> tscale_inv_amn(os+ot+1, lh);
        
        FlatVector<double> powscale(os+ot+1, lh);
        double prod = 1;
        for (int i = 0; i <= os; i++, prod *= -scale*tscale)
          powscale(i) = prod;

        
        // (185) from paper 'fast, exact, stable, Gumerov+Duraiswami
        SphericalBessel (os+ot, kappa*abs(z), tscale, trafo.Col(0));
        
        if (z < 0)
          for (int l = 1; l < trafo.Height(); l+=2)
            trafo(l,0) *= -1;
        
        for (int l = 0; l <= os+ot; l++)
          trafo(l,0) *= sqrt(2*l+1);
        

        // *testout << "Norm Col0 = " << L2Norm(trafo.Col(0)) << endl;
        
        if (os > 0)
          {
            for (int l = 1; l < os+ot; l++)   
              trafo(l,1) = -scale/sh.CalcAmn(0,0) * (sh.CalcAmn(0,l)*tscale*trafo(l+1,0)
                                                     -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,0));
            trafo(0,1) = -scale*tscale*trafo(1,0);
          }
        
        for (int n = 1; n < trafo.Width()-1; n++)
          {
            for (int l = 1; l < os+ot-n; l++)
              trafo(l,n+1) = -scale/sh.CalcAmn(0,n) * (sh.CalcAmn(0,l)*tscale*trafo(l+1,n)
                                                       -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,n)
                                                       -sh.CalcAmn(0,n-1)*scale*trafo(l,n-1));
            trafo(0,n+1) = powscale(n+1)*trafo(n+1,0);
          }
        
        // use symmetry of matrix (up to scaling)
        // important for small kappa, and R->R
        for (int n = 0; n < trafo.Width(); n++)
          for (int l = n+1; l < trafo.Width(); l++)
            trafo(n,l) = powscale(l-n) * trafo(l,n);
        
        
        for (int n = 0; n <= os; n++)
          hv1(n) = sh.Coef(n,0);
        hv2 = trafo.Rows(ot+1).Cols(os+1) * hv1;
        for (int n = 0; n <= ot; n++)
          target.SH().Coef(n,0) = hv2(n);
        
        // tlevel0.Stop();
      
        for (int m = 1; m <= min(os,ot); m++)
          {
            // tinit.Start();
            for (int l = m-1; l < os+ot-m; l++)
              {
                amn(l) = sh.CalcAmn(m,l);
                tscale_inv_amn(l) = tscale/amn(l);
              }
            
            trafo.Swap (oldtrafo);
            // trafo = trafo_type(0.0);
            // trafo = trafo_type(1.0/0.0);
            // trafo.Col(m-1) = trafo_type(0.0);  // for wall initial value
            
            // fill recursive formula (187)
            for (int l = m; l <= os+ot-m; l++)
              trafo(l,m) = scale/sh.CalcBmn(-m, m) * (sh.CalcBmn(-m, l)*inv_tscale*oldtrafo(l-1, m-1)  
                                                      -sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,m-1));  
            // tinit.Stop();
            // twall.Start();
            {
              int l = os+ot-m-1;
              if (m < trafo.Width()-1)
                trafo(l,m+1) = scale/sh.CalcBmn(-m,m+1)* ( // sh.CalcBmn(m-1,m) * scale*trafo(l,m-1)
                                                           - sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,m)     
                                                           + sh.CalcBmn(-m,l)  * 1/tscale*oldtrafo(l-1,m) );
            }
            for (int n = m+1; n < trafo.Width()-1; n++)
              {
                int l = os+ot-n-1;
                trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n) * scale*trafo(l,n-1)
                                                          - sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,n)     
                                                          + sh.CalcBmn(-m,l)  * 1/tscale*oldtrafo(l-1,n) );
                
                trafo(l-1,n) = tscale_inv_amn(l-1) * (amn(l)  *   tscale*trafo(l+1,n)
                                                      - amn(n-1)* scale*trafo(l,n-1)
                                                      + amn(n)* inv_scale*trafo(l,n+1));
              }
            
            // the same thing 1 row up
            {
              int l = os+ot-m-2;
              if (m < trafo.Width()-2)          
                trafo(l,m+1) = scale/sh.CalcBmn(-m,m+1)* (// sh.CalcBmn(m-1,m)     * scale*trafo(l,m-1)
                                                          - sh.CalcBmn(m-1,l+1) * tscale* oldtrafo(l+1,m)   
                                                          + sh.CalcBmn(-m,l)    * 1/tscale* oldtrafo(l-1,m) );
            }
            for (int n = m+1; n < trafo.Width()-2; n++)
              {
                // int l = 2*order-n-2;
                int l = os+ot-n-2;
                trafo(l,n+1) = scale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n)     * scale*trafo(l,n-1)
                                                          - sh.CalcBmn(m-1,l+1) * tscale* oldtrafo(l+1,n)   
                                                          + sh.CalcBmn(-m,l)    * 1/tscale* oldtrafo(l-1,n) ); 
                
                trafo(l-1,n) = tscale_inv_amn(l-1) * (amn(l)   * tscale*trafo(l+1,n)
                                                      -amn(n-1)* scale*trafo(l,n-1) +
                                                      amn(n)   * 1/scale*trafo(l,n+1)) ;
              }
            
            // twall.Stop();
            // tfill.Start();
            
            for (int l = os+ot-m; l >= m; l--)
              for (int n = m+1; n < min<int>(os+ot-l,l); n++)
                {
                  trafo(l-1,n) = tscale_inv_amn(l-1)* ( amn(l)  * tscale*trafo(l+1,n)
                                                        -amn(n-1)* scale*trafo(l,n-1)
                                                        +amn(n)  * inv_scale*trafo(l,n+1)) ;
                }
            
            // thhfill.Stop();
            // tfill.AddFlops (sqr( (os+ot)/2-m ));

            // ttrans.Start();
            for (int n = m; n < os; n++)
              for (int l = n+1; l <= os; l++)
                trafo(n,l) = powscale(l-n) * trafo(l,n);              
            // ttrans.Stop();
            
            // tmult.Start();
            for (int n = m; n <= os; n++)
              hv1(n) = sh.Coef(n,m);
            hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
            for (int n = m; n <= ot; n++)
              target.SH().Coef(n,m) = hv2(n);
            
            for (int n = m; n <= os; n++)
              hv1(n) = sh.Coef(n,-m);
            hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
            for (int n = m; n <= ot; n++)
              target.SH().Coef(n,-m) = hv2(n);

            // tmult.Stop();
          }
        // RegionTimer cleanup(tcleanup);
      }
    
    
    

    
    if constexpr (std::is_same<RADIAL,MPRegular>::value && std::is_same<TARGET,MPRegular>::value)
      {
        // reg->reg transformation
        
        LocalHeap lh(( 32*( (os+ot+1)*(os+ot+1) + (os+1 + ot+1) ) + 8*3*(os+ot+1) + 500));
        
        typedef double trafo_type;
        
        FlatMatrix<trafo_type> trafo(os+ot+1, os+ot+1, lh);
        FlatMatrix<trafo_type> oldtrafo(os+ot+1, os+ot+1, lh);     
        FlatVector<entry_type> hv1(os+1, lh), hv2(ot+1, lh);
        
        trafo = trafo_type(0.0);
        
    
        FlatVector<double> amn(os+ot+1, lh);
        FlatVector<double> inv_scale_inv_amn(os+ot+1, lh);
        
        FlatVector<double> powscale(os+ot+1, lh);
        double prod = 1;
        for (int i = 0; i <= os; i++)
          {
            powscale(i) = prod;
            prod /= -scale*tscale;
          }

        
        // (185) from paper 'fast, exact, stable, Gumerov+Duraiswami
        SphericalBessel (os+ot, kappa*abs(z), 1/scale, trafo.Row(0));
        
        if (z > 0)
          for (int l = 1; l < trafo.Width(); l+=2) trafo(0,l) *= -1;
        
        for (int l = 0; l <= os+ot; l++)
          trafo(0,l) *= sqrt(2*l+1);
        

        if (ot > 0)
          {
            for (int l = 1; l < os+ot; l++)   
              trafo(1,l) = -1/tscale/sh.CalcAmn(0,0) * (sh.CalcAmn(0,l)/scale*trafo(0,l+1)
                                                        -sh.CalcAmn(0,l-1)*scale*trafo(0,l-1));
            trafo(1,0) = 1.0/(-scale*tscale)*trafo(0,1);
          }
        
        for (int n = 1; n < trafo.Height()-1; n++)
          {
            for (int l = 1; l < os+ot-n; l++)
              trafo(n+1,l) = -1.0/tscale/sh.CalcAmn(0,n) * (sh.CalcAmn(0,l)*1/scale*trafo(n,l+1)
                                                            -sh.CalcAmn(0,l-1)*scale*trafo(n,l-1)
                                                            -sh.CalcAmn(0,n-1)*1/tscale*trafo(n-1,l));
            trafo(n+1,0) = powscale(n+1)*trafo(0,n+1);
          }
        
        // use symmetry of matrix (up to scaling)
        // important for small kappa, and R->R
        for (int n = 0; n < trafo.Height(); n++)
          for (int l = n+1; l < trafo.Height(); l++)
            trafo(l,n) = powscale(l-n) * trafo(n,l);

        
        for (int n = 0; n <= os; n++)
          hv1(n) = sh.Coef(n,0);
        hv2 = trafo.Rows(ot+1).Cols(os+1) * hv1;
        for (int n = 0; n <= ot; n++)
          target.SH().Coef(n,0) = hv2(n);

        for (int m = 1; m <= min(os,ot); m++)
          {
            
            for (int l = m-1; l < os+ot-m; l++)
              {
                amn(l) = sh.CalcAmn(m,l);
                inv_scale_inv_amn(l) = 1.0/scale/amn(l);
              }
            
            trafo.Swap (oldtrafo);
            // trafo = trafo_type(1.0/0.0);
            
            // fill recursive formula (187)
            trafo.Row(m-1) = trafo_type(0.0); // for wall initial value
            for (int l = m; l <= os+ot-m; l++)
              trafo(m,l) = 1/tscale/sh.CalcBmn(-m, m) * (sh.CalcBmn(-m, l)*scale*oldtrafo(m-1,l-1)  
                                                         -sh.CalcBmn(m-1,l+1)*1/scale*oldtrafo(m-1,l+1));  
            

            for (int n = m; n < trafo.Height()-1; n++)
              {
                int l = trafo.Width()-n-2;
                if (l > 0)
                  {
                    trafo(n+1,l) = 1/tscale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n) * 1/tscale*trafo(n-1,l)
                                                                 - sh.CalcBmn(m-1,l+1)*1/scale*oldtrafo(n,l+1)     
                                                                 + sh.CalcBmn(-m,l)  * scale*oldtrafo(n,l-1) );
                    
                    if (n > m)
                      trafo(n,l-1) = inv_scale_inv_amn(l-1) * (amn(l)  *   1/scale*trafo(n,l+1)
                                                               - amn(n-1)* 1/tscale*trafo(n-1,l)
                                                               + amn(n)* tscale*trafo(n+1,l));
                  }
              }
            
            // the same thing 1 row up
            for (int n = m; n < trafo.Height()-2; n++)
              {
                int l = trafo.Width()-n-3;
                if (l > 0)
                  {
                    trafo(n+1,l) = 1/tscale/sh.CalcBmn(-m,n+1)* (sh.CalcBmn(m-1,n)     * 1/tscale*trafo(n-1,l)
                                                                 - sh.CalcBmn(m-1,l+1) * 1/scale* oldtrafo(n,l+1)   
                                                                 + sh.CalcBmn(-m,l)    * scale* oldtrafo(n,l-1) ); 
                    
                    if (n > m)
                      trafo(n,l-1) = inv_scale_inv_amn(l-1) * (amn(l)   * inv_scale*trafo(n,l+1)
                                                               -amn(n-1)* inv_tscale*trafo(n-1,l) +
                                                               amn(n)   * tscale*trafo(n+1,l)) ;
                  }
              }
              

            for (int l = trafo.Width()-1; l >= m; l--)
              for (int n = m+1; n < min<int>(trafo.Width()-1-l,l); n++)
                {
                  trafo(n,l-1) = inv_scale_inv_amn(l-1)* ( amn(l)  * inv_scale*trafo(n,l+1)
                                                           -amn(n-1)* inv_tscale*trafo(n-1,l)
                                                           +amn(n)  * tscale*trafo(n+1,l)) ;
                }


            for (int n = m; n < os; n++)
              for (int l = n+1; l <= os; l++)
                trafo(l,n) = powscale(l-n) * trafo(n,l);

            
            for (int n = m; n <= os; n++)
              hv1(n) = sh.Coef(n,m);
            hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
            for (int n = m; n <= ot; n++)
              target.SH().Coef(n,m) = hv2(n);
            
            for (int n = m; n <= os; n++)
              hv1(n) = sh.Coef(n,-m);
            hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
            for (int n = m; n <= ot; n++)
              target.SH().Coef(n,-m) = hv2(n);
          }
      }
  }
#endif


  
  
  
  template void MultiPole<MPRegular> :: ShiftZ (double z, MultiPole<MPRegular> & target);
  template void MultiPole<MPSingular> :: ShiftZ (double z, MultiPole<MPRegular> & target);  
  template void MultiPole<MPSingular> :: ShiftZ (double z, MultiPole<MPSingular> & target);

  template void MultiPole<MPRegular,Vec<3,Complex>> :: ShiftZ (double z, MultiPole<MPRegular,Vec<3,Complex>> & target);
  template void MultiPole<MPSingular,Vec<3,Complex>> :: ShiftZ (double z, MultiPole<MPRegular,Vec<3,Complex>> & target);  
  template void MultiPole<MPSingular,Vec<3,Complex>> :: ShiftZ (double z, MultiPole<MPSingular,Vec<3,Complex>> & target);

  template void MultiPole<MPRegular,Vec<4,Complex>> :: ShiftZ (double z, MultiPole<MPRegular,Vec<4,Complex>> & target);
  template void MultiPole<MPSingular,Vec<4,Complex>> :: ShiftZ (double z, MultiPole<MPRegular,Vec<4,Complex>> & target);  
  template void MultiPole<MPSingular,Vec<4,Complex>> :: ShiftZ (double z, MultiPole<MPSingular,Vec<4,Complex>> & target);

  template void MultiPole<MPRegular,Vec<6,Complex>> :: ShiftZ (double z, MultiPole<MPRegular,Vec<6,Complex>> & target);
  template void MultiPole<MPSingular,Vec<6,Complex>> :: ShiftZ (double z, MultiPole<MPRegular,Vec<6,Complex>> & target);
  template void MultiPole<MPSingular,Vec<6,Complex>> :: ShiftZ (double z, MultiPole<MPSingular,Vec<6,Complex>> & target);






  // https://fortran-lang.discourse.group/t/looking-for-spherical-bessel-and-hankel-functions-of-first-and-second-kind-and-arbitrary-order/2308/2

  // adapted from fmm3d
  template <typename Tz> 
  void T_besseljs3d (int nterms, Tz z, double scale,
                     FlatVector<Tz> fjs, FlatVector<Tz> fjder)
  {
    /*
      c**********************************************************************
      c
      c PURPOSE:
      c
      c	This subroutine evaluates the first NTERMS spherical Bessel 
      c	functions and if required, their derivatives.
      c	It incorporates a scaling parameter SCALE so that
      c       
      c		fjs_n(z)=j_n(z)/SCALE^n
      c		fjder_n(z)=\frac{\partial fjs_n(z)}{\partial z}
      c
      c INPUT:
      c
      c    nterms (integer): order of expansion of output array fjs 
      c    z     (complex *16): argument of the spherical Bessel functions
      c    scale    (real *8) : scaling factor (discussed above)
      c    ifder  (integer): flag indicating whether to calculate "fjder"
      c		          0	NO
      c		          1	YES
      c OUTPUT:
      c    fjs   (complex *16): array of scaled Bessel functions.
      c    fjder (complex *16): array of derivs of scaled Bessel functions.
      c
      c
    */

  
    // c ... Initializing ...

    // set to asymptotic values if argument is sufficiently small
    if (abs(z) < 1e-200)
      {
        fjs(0) = 1;
        for (int i = 1; i <= nterms; i++)
          fjs(i) = 0.0;

        if (fjder.Size())
          {
            fjder = 0.0;
            fjder(1) = 1.0/(3*scale);
          }
        return;
      }


    //  ... Step 1: recursion up to find ntop, starting from nterms

    Tz zinv=1.0/z;

    Tz fjm1 = 0.0;
    Tz fj0 = 1.0;

    /*
      c
      cc     note max point for upward recurrence is
      c      hard coded to nterms + 1000,
      c      this might cause loss of accuracy for some
      c      arguments in the complex plane for large 
      c      nterms. For example, it is a terrible idea
      c      to use this code for z>>nterms^2
    */
  
    // int lwfjs = nterms + 100000;
    int ntop = nterms+1000;

    for (int i = nterms; ; i++)
      {
        double dcoef = 2*i+1.0;
        Tz fj1 = dcoef*zinv*fj0-fjm1;
        double dd = sqr(abs(fj1));
        if (dd > 1e40)
          {
            ntop=i+1;
            break;
          }
        fjm1 = fj0;
        fj0 = fj1;
        if (i > nterms+100000)
          throw Exception("bessel failed 1");
      }
    
    Array<bool> iscale(ntop+1);
    Vector<Tz> fjtmp(ntop+1);

    /*
      c ... Step 2: Recursion back down to generate the unscaled jfuns:
      c             if magnitude exceeds UPBOUND2, rescale and continue the 
      c	      recursion (saving the order at which rescaling occurred 
      c	      in array iscale.
    */

    iscale = false;

    fjtmp(ntop) = 0.0;
    fjtmp(ntop-1) = 1.0;
    for (int i = ntop-1; i>=1; i--)
      {
        double dcoef = 2*i+1.0;
        fjtmp(i-1) = dcoef*zinv*fjtmp(i)-fjtmp(i+1);
        double dd = sqr(abs(fjtmp(i-1)));
        if (dd > 1e40)
          {
            fjtmp(i) *= 1e-40;
            fjtmp(i-1) *= 1e-40;
            iscale[i] = true;
          }
      }

    /*
      c
      c ...  Step 3: go back up to the top and make sure that all
      c              Bessel functions are scaled by the same factor
      c              (i.e. the net total of times rescaling was invoked
      c              on the way down in the previous loop).
      c              At the same time, add scaling to fjs array.
      c
    */
  
    double scalinv = 1.0/scale;
    double sctot = 1.0;
    for (int i = 1; i <= ntop; i++)
      {
        sctot *= scalinv;
        if (iscale[i-1])
          sctot *= 1e-40;
        fjtmp(i) *= sctot;
      }
  
    //  Determine the normalization parameter:
  
    fj0=sin(z)*zinv;
    Tz fj1=fj0*zinv-cos(z)*zinv;
    if (abs(z) < 1e-2)
      {
        fj1 = 1.0/3 * z - 1.0/30 * z*z*z + 1.0/840 * z*z*z*z*z;
        // |err| < 8/9! * 0.01**7  = 2e-19
      }
    
    double d0=abs(fj0);
    double d1=abs(fj1);
    Tz zscale;
    if (d1 > d0) 
      zscale=fj1/(fjtmp(1)*scale);
    else
      zscale=fj0/fjtmp(0);

    // Scale the jfuns by zscale:
      
    Tz ztmp=zscale;
    for (int i = 0; i <= nterms; i++)
      fjs(i)=fjtmp(i)*ztmp;


    // Finally, calculate the derivatives if desired:

    if (fjder.Size())
      {
        fjder(0) = -fjs(1)*scale;
        for (int i = 1; i <= nterms; i++)
          {
            double dc1=i/(2*i+1.0);
            double dc2=1.0-dc1;
            dc1=dc1*scalinv;
            dc2=dc2*scale;
            fjder(i)=(dc1*fjtmp(i-1)-dc2*fjtmp(i+1))*ztmp;
          }
      }
  }
  

  void besseljs3d (int nterms, double z, double scale, FlatVector<double> fjs, FlatVector<double> fjder)
  {
    T_besseljs3d (nterms, z, scale, fjs, fjder);
  }
  void besseljs3d (int nterms, Complex z, double scale, FlatVector<Complex> fjs, FlatVector<Complex> fjder)
  {
    T_besseljs3d (nterms, z, scale, fjs, fjder);
  }

  /*
  // from A. Barnett 
  // http://www.fresco.org.uk/functions/barnett/index.htm

  spherical bessel functions of first (the j_n) and second (the y_n) kind.

  j0(r) = sin(r)/r
  j1(r) = (sin(r)-r cos(r)) / r**2

  y0(r) = -cos(r)/r
  y1(r) = (-cos(r)-r*sin(r)) / r**2
  */
  void SBESJY (double x, int lmax,
               FlatVector<double> j,
               FlatVector<double> y,
               FlatVector<double> jp,
               FlatVector<double> yp)
  {
    if (x < 1e-8)
      {
        cout << "TODO: special treatment for small x" << endl;
        return;
      }

    double xinv = 1/x;

    if (lmax > 0)
      {
        double twoxi = 2*xinv;
        double sl = lmax*xinv;
        double tk = 2*sl+3*xinv;
        double cf1 = sl;
        double den = 1;
        if (abs(cf1) < 1e-8) cf1 = 1e-8;
        double c = cf1;
        double d = 0;
        for (int l = 1; l < 10000; l++)
          {
            c = tk-1/c;
            d = tk-d;
            if (abs(c) < 1e-8) c = 1e-8;
            if (abs(d) < 1e-8) d = 1e-8;
            d = 1/d;
            double dcf1 = d*c;
            cf1 *= dcf1;
            if (d < 0) den = -den;
            if (abs(dcf1-1) < 1e-10)
              break;
            tk += twoxi;
            // nfp = l;
          }

        j(lmax) = den;
        jp(lmax) = cf1*den;
        for (int l = lmax; l >= 1; l--)
          {
            j(l-1) = (sl+xinv)*j(l) + jp(l);
            sl = sl-xinv;
            jp(l-1) = sl*j(l-1) - j(l);
          }
      }

    double j0 = j(0);
    double jp0 = jp(0);

    // C------ CALCULATE THE L=0 SPHERICAL BESSEL FUNCTIONS DIRECTLY
    j(0) = xinv * sin(x);
    y(0) = -xinv * cos(x);
    jp(0) = -y(0)-xinv*j(0);
    yp(0) = j(0)-xinv*y(0);

    double omega = (abs(j0)>abs(jp0)) ? j(0)/j0 : jp(0) / jp0;    // fix for x \approx 2 pi
    double sl = 0;
    for (int l = 1; l <=lmax; l++)
      {
        j(l) *= omega;
        jp(l) *= omega;
        y(l) = sl*y(l-1) - yp(l-1);
        sl += xinv;
        yp(l) = y(l-1) - (sl+xinv)*y(l);
      }
  }



  /* ************************************** Multipole class ******************************* */

  template <typename RADIAL, typename entry_type>
  entry_type MultiPole<RADIAL,entry_type> :: Eval (Vec<3> x) const
  {
    if (sh.Order() < 0) return entry_type{0.0};

    Vector<Complex> radial(sh.Order()+1);
    Vector<entry_type> shvals(sh.Order()+1);
      
    // RADIAL::Eval(sh.Order(), kappa*L2Norm(x), scale, radial);
    RADIAL::Eval(sh.Order(), kappa, L2Norm(x), rtyp, radial);
    sh.EvalOrders (x, shvals);
      
    entry_type sum{0.0};
    for (int i = 0; i <= sh.Order(); i++)
      sum +=  radial(i) * shvals(i);

    return sum;
  }

  template <typename RADIAL, typename entry_type>
  entry_type MultiPole<RADIAL,entry_type> :: EvalDirectionalDerivative (Vec<3> x, Vec<3> d) const
  {
    if (sh.Order() < 0) return entry_type{0.0};
    MultiPole<RADIAL, entry_type> tmp(Order(), kappa, RTyp());
    this->SH().DirectionalDiffAdd(kappa*d, tmp.SH(), Scale());
    return tmp.Eval(x);
  }


  template <typename RADIAL, typename entry_type>
  void MultiPole<RADIAL,entry_type> :: AddCharge (Vec<3> x, entry_type c)
  {      
    if constexpr (!std::is_same<RADIAL,MPSingular>())
      throw Exception("AddCharge assumes singular MP");
      
    // static Timer t("mptool AddCharge"); RegionTimer rg(t);
    if (L2Norm(x) < 1e-50)
      {
        sh.Coef(0,0) += Complex(0,1)*kappa/sqrt(4*M_PI) * c;
        return;
      }
    // cout << "add charge, kappa rho = " << kappa*L2Norm(x) << ", order = " << sh.Order() << endl;
      
    Vector<Complex> radial(sh.Order()+1);
    Vector<Complex> sh_shapes(sqr (sh.Order()+1));
    
    RADIAL::Eval(sh.Order(), kappa*L2Norm(x), 1.0/Scale(), radial);
    // RADIAL::Eval(sh.Order(), kappa, L2Norm(x), rtyp, radial);
    // MPRegular::Eval(sh.Order(), kappa, L2Norm(x), rtyp, radial);
    
    // cout << "radial = " << Real(radial) << endl;
    sh.Calc(x, sh_shapes);

    for (int i = 0; i <= sh.Order(); i++)
      {
        IntRange r(sqr(i), sqr(i+1));
        // sh.Coefs().Range(r) += c * Complex(0,1)*kappa * radial(i).real()*Conj(sh_shapes.Range(r));
        for (auto j : r)
          sh.Coefs()(j) += Complex(0,1)*kappa * radial(i).real()*Conj(sh_shapes(j)) * c;       
      }
  }

  template <typename RADIAL, typename entry_type>
  void MultiPole<RADIAL, entry_type> :: AddDipole (Vec<3> x, Vec<3> d, entry_type c)
  {
    // static Timer t("mptool AddDipole"); RegionTimer rg(t);      
    /*
      double eps = 1e-4;
      AddCharge(x+eps*d, -c/(2*eps));
      AddCharge(x-eps*d, c/(2*eps));
      return;
    */
      
    if constexpr (!std::is_same<RADIAL,MPSingular>())
      throw Exception("AddDipole assumes singular MP");

    /*
    // book, formula (2.2.20)
    // dipole in origin:
    MultiPole<MPSingular> tmp(1, kappa);
    tmp.SH().Coef(1,1)  += Complex(0,1)*sqr(kappa)*sh.CalcBmn(-1,1)/(2*sqrt(4*M_PI)) * d(0)*c;
    tmp.SH().Coef(1,-1) += Complex(0,1)*sqr(kappa)*sh.CalcBmn(-1,1)/(2*sqrt(4*M_PI)) * d(0)*c;

    tmp.SH().Coef(1,1)  += Complex(1,0)*sqr(kappa)*sh.CalcBmn(-1,1)/(2*sqrt(4*M_PI)) * d(1)*c;
    tmp.SH().Coef(1,-1) -= Complex(1,0)*sqr(kappa)*sh.CalcBmn(-1,1)/(2*sqrt(4*M_PI)) * d(1)*c;
      
    tmp.SH().Coef(1,0) += -Complex(0,1)*kappa*kappa*sh.CalcAmn(0,0)/sqrt(4*M_PI) *d(2)*c;
    tmp.TransformAdd (*this, -x);
    */

    MultiPole<MPSingular, entry_type> tmp(Order(), kappa, RTyp());
    tmp.AddCharge(x, c);
    tmp.SH().DirectionalDiffAdd (kappa*d,  this->SH(), Scale());
  }


  
  template <typename RADIAL, typename entry_type>
  void MultiPole<RADIAL, entry_type> :: AddCurrent (Vec<3> sp, Vec<3> ep, Complex j, int num)
  {
    if constexpr (!std::is_same<RADIAL,MPSingular>() || !std::is_same<entry_type, Vec<3,Complex>>())
      throw Exception("AddCurrent needs a singular vectorial MP");

    Vec<3> tau = ep-sp;
    Vec<3> tau_num = 1.0/num *  tau;
    for (int i = 0; i < num; i++)
      {
        for (int k = 0; k < 3; k++)
          {
            Vec<3> ek{0.0}; ek(k) = 1;
            Vec<3> cp = Cross(tau, ek);
            Vec<3,Complex> source{0.0};
            source(k) = j/double(num);
            if constexpr (std::is_same<entry_type, Vec<3,Complex>>())
              AddDipole (sp+(i+0.5)*tau_num, cp, source);
          }
      }
  }


  
  template <typename entry_type>
  shared_ptr<RegularMLMultiPoleCF<entry_type>>
  SingularMLMultiPoleCF<entry_type> :: CreateRegularExpansion(Vec<3> center, double r) const
  {
    mlmp->CalcMP();
    return make_shared<RegularMLMultiPoleCF<entry_type>> (this->MLMP(), center, r);
  }

  template class SphericalHarmonics<Complex>;
  template class SphericalHarmonics<Vec<3,Complex>>;  
  template class SphericalHarmonics<Vec<4,Complex>>;
  template class SphericalHarmonics<Vec<6,Complex>>;
  
  template class MultiPole<MPSingular>;
  template class MultiPole<MPRegular>;
  template class MultiPole<MPSingular, Vec<3,Complex>>;
  template class MultiPole<MPRegular, Vec<3,Complex>>;    
  template class MultiPole<MPSingular, Vec<4,Complex>>;
  template class MultiPole<MPRegular, Vec<4,Complex>>;    
  template class MultiPole<MPSingular, Vec<6,Complex>>;
  template class MultiPole<MPRegular, Vec<6,Complex>>;


  
  template class SingularMLMultiPoleCF<Complex>;
  template class SingularMLMultiPoleCF<Vec<3,Complex>>;
  
  template<>
  Array<size_t> RegularMLMultiPole<Complex>::nodes_on_level(100);
  template<>  
  Array<size_t> SingularMLMultiPole<Complex>::nodes_on_level(100);
  template<>
  Array<size_t> RegularMLMultiPole<Vec<3,Complex>>::nodes_on_level(100);
  template<>
  Array<size_t> SingularMLMultiPole<Vec<3,Complex>>::nodes_on_level(100);
  template<>
  Array<size_t> RegularMLMultiPole<Vec<4,Complex>>::nodes_on_level(100);
  template<>
  Array<size_t> SingularMLMultiPole<Vec<4,Complex>>::nodes_on_level(100);
  template<>
  Array<size_t> RegularMLMultiPole<Vec<6,Complex>>::nodes_on_level(100);
  template<>
  Array<size_t> SingularMLMultiPole<Vec<6,Complex>>::nodes_on_level(100);


  template class SingularMLMultiPole<Complex>;
  template class SingularMLMultiPole<Vec<3,Complex>>;

  
}
