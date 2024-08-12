#ifndef FILE_MPTOOLS
#define FILE_MPTOOLS

/*
  tools for computing with spherical harmonics and multi-poles
 */


#include <bla.hpp>
#include <coefficient.hpp>


namespace ngfem
{


  
  class SphericalHarmonics
  {
    size_t order;
    Vector<Complex> coefs;

  public:
    SphericalHarmonics (int aorder)
      : order(aorder), coefs(sqr(order+1)) { coefs=0.0; }

    int Order() const { return order; }
    FlatVector<Complex> Coefs() const { return coefs; }

    Complex & Coef(int n, int m)
    {
      return coefs(n*(n+1) + m);  
    }

    auto CoefsN (int n) const
    {
      return FlatVector<Complex> (2*n+1, &coefs(n*(n+1)-n)); 
    }
  
    Complex Eval (Vec<3> x) const
    {
      double phi, theta;
      if (x(0) == 0 && x(1) == 0)
        {
          phi = 0;
          theta = x(2) > 0 ? 0 : M_PI;
        }
      else
        {
          phi = atan2(x(1), x(0));
          theta = acos(x(2)/L2Norm(x));
        }
      return Eval(theta, phi);
    }
  
    Complex Eval (double theta, double phi) const
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

      Complex sum = 0.0;
      int ii = 0;
      for (int n = 0; n <= order; n++)
        {
          for (int m = -n; m < 0; m++, ii++)
            sum += coefs(ii) * conj(exp_imphi(-m)) * legfunc(-m, n);
          for (int m = 0; m <= n; m++, ii++)
            sum += coefs(ii) * exp_imphi(m) * legfunc(m, n);
        }

      sum /= sqrt(4*M_PI);
      return sum;
    }




    Complex EvalOrder (int n, Vec<3> x) const
    {
      double phi, theta;
      if (x(0) == 0 && x(1) == 0)
        {
          phi = 0;
          theta = x(2) > 0 ? 0 : M_PI;
        }
      else
        {
          phi = atan2(x(1), x(0));
          theta = acos(x(2)/L2Norm(x));
        }
      return EvalOrder(n, theta, phi);
    }
  
    Complex EvalOrder (int n, double theta, double phi) const
    {
      static Timer t("mptool sh evalorder"); RegionTimer rg(t);
      
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

      Complex sum = 0.0;
      auto coefsn = CoefsN(n);
      int ii = 0;
      // for (int n = 0; n <= order; n++)
      {
        for (int m = -n; m < 0; m++, ii++)
          sum += coefsn(ii) * conj(exp_imphi(-m)) * legfunc(-m, n);
        for (int m = 0; m <= n; m++, ii++)
          sum += coefsn(ii) * exp_imphi(m) * legfunc(m, n);
      }

      sum /= sqrt(4*M_PI);
      return sum;
    }



    void EvalOrders (Vec<3> x, FlatVector<Complex> vals) const
    {
      double phi, theta;
      if (x(0) == 0 && x(1) == 0)
        {
          phi = 0;
          theta = x(2) > 0 ? 0 : M_PI;
        }
      else
        {
          phi = atan2(x(1), x(0));
          theta = acos(x(2)/L2Norm(x));
        }
      return EvalOrders(theta, phi, vals);
    }
  
    void EvalOrders (double theta, double phi, FlatVector<Complex> vals) const
    {
      static Timer t("mptool sh evalorders"); RegionTimer rg(t);
      
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
          
          Complex sum = 0.0;
          for (int m = -n; m < 0; m++, ii++)
            sum += coefsn(ii) * conj(exp_imphi(-m)) * legfunc(-m, n);
          for (int m = 0; m <= n; m++, ii++)
            sum += coefsn(ii) * exp_imphi(m) * legfunc(m, n);
          vals(n) = sum;
        }

      vals /= sqrt(4*M_PI);
    }




    
    void RotateZ (double alpha)
    {
      static Timer t("mptool sh RotateZ"); RegionTimer rg(t);
      
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




    double CalcAmn (int m, int n)
    {
      if (m < 0) m=-m;
      if (n < m) return 0;
      return sqrt( (n+1.0+m)*(n+1.0-m) / ( (2*n+1)*(2*n+3) ));
    }
  
    double CalcBmn (int m, int n)
    {
      double sgn = (m >= 0) ? 1 : -1;
      if ( (m > n) || (-m > n) ) return 0;
      return sgn * sqrt( (n-m-1.0)*(n-m) / ( (2*n-1.0)*(2*n+1))); 
    }
  
    double CalcDmn (int m, int n)
    {
      double sgn = (m >= 0) ? 1 : -1;
      return sgn/2 * sqrt((n-m)*(n+m+1));
    }
    
    
    void RotateY (double alpha)
    {
      static Timer t("mptool sh RotateY"); RegionTimer rg(t);
      
      double s = sin(alpha);
      double c = cos(alpha);

      Matrix<> normalized_leg_func(order+2, order+2);
      NormalizedLegendreFunctions(order+1, order+1, c, normalized_leg_func);

      if (alpha < 0)
        for (int i = 1; i <= order; i+=2)
          normalized_leg_func.Row(i) *= -1;
      
      // cout << "leg = " << endl << normalized_leg_func << endl;
      Vector<> Dmn(2*order+1), invDmn(2*order+1);
      
      for (int n=1; n <= order; n++)
        {
          Matrix<double> trafo(n+1, 2*n+1);
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
          Vector tmp = 1.0/sqrt(2*n+3) * normalized_leg_func.Col(n+1).Range(n+2);
          for (int m = 1; m < tmp.Size(); m += 2) tmp(m) *= -1;
          for (int m = 1; m <= n; m++)
            trafo.Col(n+1)(m) = 1/CalcBmn(0,n+1) * (  CalcBmn(-m-1, n+1)*(1-c)/2 * tmp(m+1)
                                                      - CalcBmn(m-1,n+1)*(1+c)/2 * tmp(m-1)
                                                      - CalcAmn(m,n) * s*tmp(m));

          // Step 4
          // diamond - recursion
          for (int mp = -n; mp <= n; mp++)
            {
              Dmn(order+mp) = CalcDmn(mp, n);
              invDmn(order+mp) = 1.0/Dmn(order+mp);
            }
          
          for (int mp = 1; mp < n; mp++)
            {
              for (int m = mp; m < n; m++)
                trafo(m, n+mp+1) = invDmn(order+mp,n) * ( Dmn(order+mp-1) *trafo(m  ,n+mp-1)
                                                          -Dmn(order+m-1)*trafo(m-1,n+mp)
                                                          +Dmn(order+m)  *trafo(m+1,n+mp));
              int m = n;
              trafo(m, n+mp+1) = invDmn(order+mp,n) * ( Dmn(order+mp-1,n)*trafo(m  ,n+mp-1)
                                                        -Dmn(order+m-1,n)*trafo(m-1,n+mp));
            }
      
          // Step 5
          // diamond - recursion, negative      
          for (int mp = 0; mp > -n; mp--)
            {
              for (int m = -mp+1; m < n; m++)
                trafo(m, n+mp-1) = invDmn(order+mp-1,n) * (  Dmn(order+mp,n)*trafo(m  ,n+mp+1)
                                                             +Dmn(order+m-1,n)*trafo(m-1,n+mp)
                                                             -Dmn(order+m  ,n)*trafo(m+1,n+mp));
              int m = n;
              trafo(m, n+mp-1) = invDmn(order+mp-1,n) * (  Dmn(order+mp,n)*trafo(m  ,n+mp+1)
                                                           +Dmn(order+m-1,n)*trafo(m-1,n+mp));
            }

      
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




          // cout << "trafo = " << trafo << endl;
          
          // for (int m = 1; m <= n; m+=2)
          // trafo.Row(m) *= -1;
          
          // check orthogonality
          // cout << "trafo: " << endl << trafo << endl;
          // cout << "ortho: " << endl << trafo*Trans(trafo) << endl;
          // auto rdlevel = trafo.Rows(ij+1).Cols(n-ij, n+ij+1);
          /*
            cout << "n = " << ij << ", err ortho = " << L2Norm( Matrix(rdlevel*Trans(rdlevel)) - Identity(ij+1))
            << " row0 = " << L2Norm2(rdlevel.Row(0))-1 
            << " row1 = " << L2Norm2(rdlevel.Row(1))-1 
            << endl;
          */
          FlatVector<Complex> cn = CoefsN(n);
          Vector<Complex> old = cn;
          
          cn = Trans(trafo) * old.Range(n, 2*n+1);
          cn.Slice(0,1).Reversed() += Trans(trafo.Rows(1,n+1)) * old.Range(0,n).Reversed();

          for (int m = 1; m <= n; m+=2)
            {
              cn(n+m) *= -1;
              cn(n-m) *= -1;
            }
        }
    }
    
  };



  // from A. Barnett 
  // http://www.fresco.org.uk/functions/barnett/index.htm

  void SBESJY (double x, int lmax,
               FlatVector<double> j,
               FlatVector<double> y,
               FlatVector<double> jp,
               FlatVector<double> yp)
  {
    if (x < 1e-8)
      {
        cout << "special treatment for small x" << endl;
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
        for (int l = 1; l < 1000; l++)
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
    double den = j(0);

    // C------ CALCULATE THE L=0 SPHERICAL BESSEL FUNCTIONS DIRECTLY
    j(0) = xinv * sin(x);
    y(0) = -xinv * cos(x);
    jp(0) = -y(0)-xinv*j(0);
    yp(0) = j(0)-xinv*y(0);

    double omega = j(0)/den;
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


  
  // Gumerov+Duraiswami book, page 57 (unstable)
  template <typename T>
  void SphericalBessel (int n, double rho, T && values)
  {
    if (rho == 0)
      {
        if (n >= 0) values(0) = 1;
        for (int i = 1; i <= n; i++)
          values(i) = 0;
      }

    Vector j(n+1), y(n+1), jp(n+1), yp(n+1);
    SBESJY (rho, n, j, y, jp, yp);
    values = j;
    /*
    if (n >= 0)
      values(0) = sin(rho)/rho;
    if (n >= 1)
      values(1) = (sin(rho)-rho*cos(rho)) / (rho*rho);

    for (int i = 2; i <= n; i++)
      values(i) = (2*i-1)/rho * values(i-1) - values(i-2);
    */
  }


  template <typename T>
  void SphericalHankel1 (int n, double rho, T && values)
  {
    Complex imag(0,1);
    /*
    if (n >= 0)
      values(0) = exp(imag*rho) / (imag*rho);
    if (n >= 1)
      values(1) = -imag*values(0)*(1.0-1.0/(imag*rho));
    
    for (int i = 2; i <= n; i++)
      values(i) = (2*i-1)/rho * values(i-1) - values(i-2);
    */
    
    if (rho < 1e-8)
      {
        values = Complex(0);
        return;
      }
    Vector j(n+1), y(n+1), jp(n+1), yp(n+1);
    SBESJY (rho, n, j, y, jp, yp);

    values = j + Complex(0,1) * y;
  }




  


  template <typename RADIAL>
  class MultiPole
  {
    SphericalHarmonics sh;
    double kappa;
    // Vec<3> center;
    // double r;
    
  public:
    MultiPole (int aorder, double akappa) // , Vec<3> acenter = Vec<3>{0,0,0}, double ar = 1)
      : sh(aorder), kappa(akappa) { } // , center(acenter), r(ar) { }

    Complex & Coef(int n, int m) { return sh.Coef(n,m); }
    auto & SH() { return sh; }
    Complex Eval (Vec<3> x) const
    {
      // x -= center;
      Vector<Complex> radial(sh.Order()+1);
      Vector<Complex> shvals(sh.Order()+1);
      
      RADIAL::Eval(sh.Order(), kappa*L2Norm(x), radial);
      sh.EvalOrders (x, shvals);
      
      Complex sum = 0;
      for (int i = 0; i <= sh.Order(); i++)
        sum +=  shvals(i) * radial(i);

      return sum;
    }


    template <typename TARGET>
    void ShiftZ (double z, MultiPole<TARGET> & target)
    {
      static Timer t("mptool ShiftZ"); RegionTimer rg(t);

      int os = sh.Order();
      int ot = target.SH().Order();

      target.SH().Coefs()=0.0;
      Matrix<Complex> trafo(os+ot+1, os+1), trafom(os+ot+1, os+1);
      trafo = Complex(0.0);

      // (185) from paper 'fast, exact, stable, Gumerov+Duraiswami
      // RADIAL::Eval(os+ot, kappa*abs(z), trafo.Col(0));
      if (typeid(RADIAL) == typeid(TARGET))
        SphericalBessel (os+ot, kappa*abs(z), trafo.Col(0));
      else
        SphericalHankel1 (os+ot, kappa*abs(z), trafo.Col(0));
      
      if (z < 0)
        for (int l = 1; l < trafo.Height(); l+=2) trafo(l,0) *= -1;
      for (int l = 0; l <= os+ot; l++)
        trafo(l,0) *= sqrt(2*l+1);
      
      for (int l = 1; l < os+ot; l++)   
        trafo(l,1) = -1.0/sh.CalcAmn(0,0) * (sh.CalcAmn(0,l)*trafo(l+1,0)-sh.CalcAmn(0,l-1)*trafo(l-1,0));
      trafo(0,1) = -trafo(1,0);
      
      for (int n = 1; n < os; n++)
        {
          for (int l = 1; l < os+ot-n; l++)
            trafo(l,n+1) = -1.0/sh.CalcAmn(0,n) * (sh.CalcAmn(0,l)*trafo(l+1,n)-sh.CalcAmn(0,l-1)*trafo(l-1,n)-sh.CalcAmn(0,n-1)*trafo(l,n-1));
          trafo(0,n+1) = pow(-1,n+1)*trafo(n+1,0);
        }
      
      Vector<Complex> hv1(os+1), hv2(ot+1);
      for (int n = 0; n <= os; n++)
        hv1(n) = sh.Coef(n,0);
      hv2 = trafo.Rows(ot+1) * hv1;
      for (int n = 0; n <= ot; n++)
        target.SH().Coef(n,0) = hv2(n);
      // cout << "hv2 = " << hv2 << endl;

      for (int m = 1; m <= min(os,ot); m++)
        {
          trafom = Complex(0.0);
          // fill recursive

          // (187)
          for (int l = m; l <= os+ot-m; l++)
            trafom(l,m) = 1/sh.CalcBmn(-m, m) * (sh.CalcBmn(-m, l)*trafo(l-1, m-1)-sh.CalcBmn(m-1,l+1)*trafo(l+1,m-1));
          
          for (int n = m; n < os; n++)
            {
              for (int l = m+1; l < os+ot-n; l++)
                trafom(l,n+1) = -1.0/sh.CalcAmn(m,n) * (sh.CalcAmn(m,l)*trafom(l+1,n)
                                                        -sh.CalcAmn(m,l-1)*trafom(l-1,n)
                                                        -sh.CalcAmn(m,n-1)*trafom(l,n-1));
              trafom(m,n+1) = pow(-1,m+n+1)*trafom(n+1,m);
            }

          // if (m == 1)
          // cout << "trafom = " << trafom.Rows(m,ot+1).Cols(m,os+1) << endl;
          
          for (int n = m; n <= os; n++)
            hv1(n) = sh.Coef(n,m);
          hv2.Range(m,ot+1) = trafom.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
          for (int n = m; n <= ot; n++)
            target.SH().Coef(n,m) = hv2(n);
          // cout << "m = " << m << ", hv1 = " << hv1 << ", hv2 = " << hv2 << endl;
          for (int n = m; n <= os; n++)
            hv1(n) = sh.Coef(n,-m);
          hv2.Range(m,ot+1) = trafom.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
          for (int n = m; n <= ot; n++)
            target.SH().Coef(n,-m) = hv2(n);
          // cout << "m = " << -m << ", hv1 = " << hv1 << ", hv2 = " << hv2 << endl;          
          trafo = trafom; // swap ? 
        }
    }
  };

  
  // hn1 = jn+ i*yn
  class MPSingular
  {
  public:
    template <typename T>
    static void Eval (int order, double r, T && values)
    {
      SphericalHankel1(order, r, values);
    }
  };
  
  // jn
  class MPRegular
  {
  public:    
    template <typename T>
    static void Eval (int order, double r, T && values)
    {
      SphericalBessel (order, r, values);
    }
  };
  

  class SphericalHarmonicsCF : public CoefficientFunction
  {
    SphericalHarmonics sh;
  public:
    SphericalHarmonicsCF (int order)
      : CoefficientFunction(1, true), sh(order) { }
    Complex & Coef(int n, int m) { return sh.Coef(n,m); } 
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("real eval not available"); }

    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const override
    {
      values(0) = sh.Eval(mip.GetPoint());
    }
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const override
    {
      for (int i = 0; i < ir.Size(); i++)
        {
          auto & mip = ir[i];
          values(i,0) = sh.Eval(mip.GetPoint());
        }
    }

    auto & SH() { return sh; }
  };



  template <typename RADIAL>
  class MultiPoleCF : public CoefficientFunction
  {
    MultiPole<RADIAL> mp;
    Vec<3> center;
  public:
    MultiPoleCF (int order, int kappa, Vec<3> acenter)
      : CoefficientFunction(1, true), mp(order, kappa), center(acenter) { }

    Complex & Coef(int n, int m) { return mp.Coef(n,m); } 
    auto & SH() { return mp.SH(); }
    auto & MP() { return mp; }
    Vec<3> Center() const { return center; }
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("real eval not available"); }

    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const override
    {
      values(0) = mp.Eval(mip.GetPoint()-center);
    }

    template <typename TARGET>
    void ShiftZ (double z, MultiPole<TARGET> & target) { mp.ShiftZ(z, target); }

    template <typename TARGET>
    void Transform (MultiPoleCF<TARGET> & target)
    {
      Vec<3> dist = target.Center() - center;
      // cout << "dist = " << dist << endl;
      double len = L2Norm(dist);
      double theta = acos (dist(2) / len);
      double phi = atan2(dist(1), dist(0));
      // cout << "phi = " << phi << ", theta = " << theta << endl;
      MultiPole<RADIAL> tmp{mp};
      tmp.SH().RotateZ(phi);
      tmp.SH().RotateY(theta);
      tmp.ShiftZ(-len, target.MP());
      target.SH().RotateY(-theta); 
      target.SH().RotateZ(-phi);
    }
    
  };
  
}
#endif
