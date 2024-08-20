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
    int order;
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

    auto Polar (Vec<3> x) const
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
      return tuple{theta, phi};
    }
    
    Complex Eval (Vec<3> x) const
    {
      auto [theta, phi] = Polar(x);
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
      auto [theta, phi] = Polar (x);
      return EvalOrder(n, theta, phi);
    }
  
    Complex EvalOrder (int n, double theta, double phi) const
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
      auto [theta, phi] = Polar(x);
      return EvalOrders(theta, phi, vals);
    }
  
    void EvalOrders (double theta, double phi, FlatVector<Complex> vals) const
    {
      static Timer ts("mptool sh evalorders small");
      static Timer tl("mptool sh evalorders large");
      RegionTimer rg(order < 10 ? ts : tl);
      
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


    void Calc (Vec<3> x, FlatVector<Complex> shapes)
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

    
    void RotateZ (double alpha)
    {
      // static Timer t("mptool sh RotateZ"); RegionTimer rg(t);
      
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




    static double CalcAmn (int m, int n)
    {
      if (m < 0) m=-m;
      if (n < m) return 0;
      return sqrt( (n+1.0+m)*(n+1.0-m) / ( (2*n+1)*(2*n+3) ));
    }
  
    static double CalcBmn (int m, int n)
    {
      double sgn = (m >= 0) ? 1 : -1;
      if ( (m > n) || (-m > n) ) return 0;
      return sgn * sqrt( (n-m-1.0)*(n-m) / ( (2*n-1.0)*(2*n+1))); 
    }
  
    static double CalcDmn (int m, int n)
    {
      double sgn = (m >= 0) ? 1 : -1;
      return sgn/2 * sqrt((n-m)*(n+m+1));
    }
    
    
    void RotateY (double alpha)
    {
      LocalHeap lh(8*4*sqr(order) + 10*order + 500);
      static Timer t("mptool sh RotateY"); RegionTimer rg(t);
      static Timer t1("mptool sh RotateY 1");
      static Timer t2("mptool sh RotateY 2");
      static Timer t3("mptool sh RotateY 3");
      static Timer t23("mptool sh RotateY 23");
      static Timer tlp("mptool sh RotateY loop");      
      
      t1.Start();
      double s = sin(alpha);
      double c = cos(alpha);

      FlatMatrix<> normalized_leg_func(order+2, order+2, lh);
      NormalizedLegendreFunctions(order+1, order+1, c, normalized_leg_func);

      if (alpha < 0)
        for (int i = 1; i <= order+1; i+=2)
          normalized_leg_func.Row(i) *= -1;
      
      // cout << "leg = " << endl << normalized_leg_func << endl;
      Vector<> Dmn(2*order+1), invDmn(2*order+1);
      t1.Stop();

      RegionTimer rlp(tlp);
      ArrayMem<double, 1000> trafomem( (order+1)*(2*order+1));
      for (int n=1; n <= order; n++)
        {
          HeapReset hr(lh);
          // RegionTimer r23(t23);
          
          // Matrix<double> trafo(n+1, 2*n+1);
          FlatMatrix<double> trafo(n+1, 2*n+1, trafomem.Data());
          /*
            Recursive Computation of Spherical Harmonic Rotation Coefficients of Large Degree
            Nail A. Gumerov and Ramani Duraiswami
            within Excursions in Harmonic Analysis, Volume 3
            
            page 130 
          */
          // t2.Start();
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


          /*
          double errortho = L2Norm( Matrix(trafo*Trans(trafo) - Identity(n+1)));
          if (errortho > 1e-10)
            {
              *testout << "n = " << n << " order = " << Order() << ", alpha = " << alpha << ", errortho = " << errortho << endl;
              if (n < 10)
                *testout << trafo*Trans(trafo) << endl;
            }
          */

          // t2.Stop();
          // t3.Start();
          
          FlatVector<Complex> cn = CoefsN(n);
          FlatVector<Complex> old = cn | lh;
          
          cn = Trans(trafo) * old.Range(n, 2*n+1);
          cn.Slice(0,1).Reversed() += Trans(trafo.Rows(1,n+1)) * old.Range(0,n).Reversed();

          for (int m = 1; m <= n; m+=2)
            {
              cn(n+m) *= -1;
              cn(n-m) *= -1;
            }
          // t3.Stop();
        }
    }
    
  };


  // https://fortran-lang.discourse.group/t/looking-for-spherical-bessel-and-hankel-functions-of-first-and-second-kind-and-arbitrary-order/2308/2

  // adapted from fmm3d
template <typename Tz> 
void besseljs3d (int nterms, Tz z, double scale,
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


  
  template <typename T>
  void SphericalBessel (int n, double rho, double scale, T && values)
  {
    Vector<double> j(n+1), jp(n+1);
    besseljs3d<double> (n, rho, scale,  j, jp);
    values = j;
  }


  template <typename T>
  void SphericalHankel1 (int n, double rho, double scale, T && values)
  {
    // Complex imag(0,1);
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
    if (scale != 1.0)
      {
        double prod = 1.0;
        for (int i = 0; i <= n; i++)
          {
            values(i) *= prod;
            prod *= scale;
          }
      }
  }




  
  // hn1 = jn+ i*yn
  class MPSingular
  {
  public:
    template <typename T>
    static void Eval (int order, double r, double scale, T && values)
    {
      SphericalHankel1(order, r, scale,  values);
    }
  };
  
  // jn
  class MPRegular
  {
  public:    
    template <typename T>
    static void Eval (int order, double r, double scale, T && values)
    {
      SphericalBessel (order, r, 1.0/scale, values);
    }
  };
  
  


  template <typename RADIAL>
  class MultiPole
  {
    SphericalHarmonics sh;
    double kappa;
    double scale;
  public:
    MultiPole (int aorder, double akappa, double ascale = 1) 
      : sh(aorder), kappa(akappa), scale(ascale) { }

    Complex & Coef(int n, int m) { return sh.Coef(n,m); }
    auto & SH() { return sh; }
    const auto & SH() const { return sh; }
    double Kappa() const { return kappa; }
    double Scale() const { return scale; }    
    int Order() const { return sh.Order(); }
    
    MultiPole<RADIAL> Truncate(int neworder) const
    {
      if (neworder > sh.Order()) neworder=sh.Order();
      MultiPole nmp(neworder, kappa);
      nmp.sh.Coefs() = sh.Coefs().Range(sqr(neworder+1));
      return nmp;
    }

    MultiPole & operator+= (const MultiPole & mp2)
    {
      size_t commonsize = min(SH().Coefs().Size(), mp2.SH().Coefs().Size());
      SH().Coefs().Range(commonsize) += mp2.SH().Coefs().Range(commonsize);
      return *this;
    }
    
    Complex Eval (Vec<3> x) const
    {
      Vector<Complex> radial(sh.Order()+1);
      Vector<Complex> shvals(sh.Order()+1);
      
      RADIAL::Eval(sh.Order(), kappa*L2Norm(x), scale, radial);
      sh.EvalOrders (x, shvals);
      
      Complex sum = 0;
      for (int i = 0; i <= sh.Order(); i++)
        sum +=  shvals(i) * radial(i);

      return sum;
    }

    void AddCharge (Vec<3> x, Complex c)
    {      
      if constexpr (!std::is_same<RADIAL,MPSingular>())
        throw Exception("AddCharge assumes singular MP");
      
      static Timer t("mptool AddCharge"); RegionTimer rg(t);      
      
      if (L2Norm(x) < 1e-10)
        {
          sh.Coef(0,0) += c * Complex(0,1)*kappa/sqrt(4*M_PI);
          return;
        }

      // cout << "add charge, kappa rho = " << kappa*L2Norm(x) << ", order = " << sh.Order() << endl;
      
      Vector<Complex> radial(sh.Order()+1);
      Vector<Complex> sh_shapes(sqr (sh.Order()+1));

      RADIAL::Eval(sh.Order(), kappa*L2Norm(x), 1.0/scale, radial);
      // cout << "radial = " << radial << endl;
      sh.Calc(x, sh_shapes);

      for (int i = 0; i <= sh.Order(); i++)
        {
          IntRange r(sqr(i), sqr(i+1));
          sh.Coefs().Range(r) += c * Complex(0,1)*kappa * radial(i).real()*Conj(sh_shapes.Range(r));
        }
    }

    
    void AddDipole (Vec<3> x, Vec<3> d, Complex c)
    {
      static Timer t("mptool AddDipole"); RegionTimer rg(t);      
      /*
      double eps = 1e-4;
      AddCharge(x+eps*d, -c/(2*eps));
      AddCharge(x-eps*d, c/(2*eps));
      return;
      */
      
      // book, formula (2.2.20)
      if constexpr (!std::is_same<RADIAL,MPSingular>())
        throw Exception("AddCharge assumes singular MP");

      // dipole in origin:
      MultiPole<MPSingular> tmp(1, kappa);
      tmp.SH().Coef(1,1)  += Complex(0,1)*sqr(kappa)*sh.CalcBmn(-1,1)/(2*sqrt(4*M_PI)) * d(0)*c;
      tmp.SH().Coef(1,-1) += Complex(0,1)*sqr(kappa)*sh.CalcBmn(-1,1)/(2*sqrt(4*M_PI)) * d(0)*c;

      tmp.SH().Coef(1,1)  += Complex(1,0)*sqr(kappa)*sh.CalcBmn(-1,1)/(2*sqrt(4*M_PI)) * d(1)*c;
      tmp.SH().Coef(1,-1) -= Complex(1,0)*sqr(kappa)*sh.CalcBmn(-1,1)/(2*sqrt(4*M_PI)) * d(1)*c;
      
      tmp.SH().Coef(1,0) += -Complex(0,1)*kappa*kappa*sh.CalcAmn(0,0)/sqrt(4*M_PI) *d(2)*c;
      tmp.TransformAdd (*this, -x);
    }


    template <typename TARGET>
    void Transform (MultiPole<TARGET> & target, Vec<3> dist) const
    {
      if (target.SH().Order() < 0) return;
      if (SH().Order() < 0)
        {
          target.SH().Coefs() = 0.0;
          return;
        }
      
      static Timer t("mptool Transform "+ToString(typeid(RADIAL).name())+ToString(typeid(TARGET).name()));      
      RegionTimer reg(t);
      double len = L2Norm(dist);
      double theta, phi;

      if (len < 1e-30)
        theta = 0;
      else
        theta = acos (dist(2) / len);

      if (sqr(dist(0))+sqr(dist(1)) < 1e-30)
        phi = 0;
      else
        phi = atan2(dist(1), dist(0));
        
      
      MultiPole<RADIAL> tmp(*this);
      tmp.SH().RotateZ(phi);
      tmp.SH().RotateY(theta);

      tmp.ShiftZ(-len, target);
      
      /*
      // testing ....
      tmp.SH().RotateY(-theta);
      tmp.SH().RotateZ(-phi);
      double err =  L2Norm(tmp.SH().Coefs()-this->SH().Coefs());
      if (err > 1e-5)
        {
          cout << "phi = " << phi << ", theta = " << theta << endl;
          cout << "norm = " << L2Norm(tmp.SH().Coefs()) << "   ";
          cout << "rotation err = " << err << ", order = " << tmp.SH().Order() << endl;
          cout << "vec = " << this->SH().Coefs() << endl;
        }
      */
      
      target.SH().RotateY(-theta);
      target.SH().RotateZ(-phi);
    }
    
    template <typename TARGET>
    void TransformAdd (MultiPole<TARGET> & target, Vec<3> dist) const
    {
      MultiPole<TARGET> tmp{target};
      Transform(tmp, dist);
      target.SH().Coefs() += tmp.SH().Coefs();
    }
    
    template <typename TARGET>
    void ShiftZ (double z, MultiPole<TARGET> & target)
    {
      static Timer t("mptool ShiftZ"+ToString(typeid(RADIAL).name())+ToString(typeid(TARGET).name()));
      RegionTimer rg(t);
      
      int os = sh.Order();
      int ot = target.SH().Order();

      target.SH().Coefs()=0.0;
      Matrix<Complex> trafo(os+ot+1, os+1), trafom(os+ot+1, os+1);
      trafo = Complex(0.0);

      // (185) from paper 'fast, exact, stable, Gumerov+Duraiswami
      // RADIAL::Eval(os+ot, kappa*abs(z), trafo.Col(0));
      if (typeid(RADIAL) == typeid(TARGET))
        SphericalBessel (os+ot, kappa*abs(z), 1.0, trafo.Col(0));
      else
        SphericalHankel1 (os+ot, kappa*abs(z), 1.0, trafo.Col(0));

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
            trafo(l,1) = -1.0/sh.CalcAmn(0,0) * (sh.CalcAmn(0,l)*trafo(l+1,0)-sh.CalcAmn(0,l-1)*trafo(l-1,0));
          trafo(0,1) = -trafo(1,0);
        }
      
      for (int n = 1; n < os; n++)
        {
          for (int l = 1; l < os+ot-n; l++)
            trafo(l,n+1) = -1.0/sh.CalcAmn(0,n) * (sh.CalcAmn(0,l)*trafo(l+1,n)-sh.CalcAmn(0,l-1)*trafo(l-1,n)-sh.CalcAmn(0,n-1)*trafo(l,n-1));
          trafo(0,n+1) = pow(-1,n+1)*trafo(n+1,0);
        }


      double prod = 1;
      for (int i = 0; i <= os; i++, prod *= scale)
        trafo.Col(i) *= prod;
      prod = 1;
      for (int i = 0; i <= os+ot; i++, prod /= target.Scale())
        trafo.Row(i) *= prod;
      
      Vector<Complex> hv1(os+1), hv2(ot+1);
      for (int n = 0; n <= os; n++)
        hv1(n) = sh.Coef(n,0);

      /*
      double prod = 1;
      for (int i = 0; i <= os; i++, prod*=scale)
        hv1(i) *= prod;
      */
      
      hv2 = trafo.Rows(ot+1) * hv1;
      /*
      prod = 1;
      for (int i = 0; i <= ot; i++, prod /= target.Scale())
        hv2(i) *= prod;
      */
      
      for (int n = 0; n <= ot; n++)
        target.SH().Coef(n,0) = hv2(n);
      // cout << "hv2 = " << hv2 << endl;

      for (int m = 1; m <= min(os,ot); m++)
        {
          trafom = Complex(0.0);
          // fill recursive

          // (187)
          for (int l = m; l <= os+ot-m; l++)
            trafom(l,m) = 1/sh.CalcBmn(-m, m) * (sh.CalcBmn(-m, l)*pow(target.Scale(),l-1)*pow(1/scale,m-1)*trafo(l-1, m-1)
                                                 -sh.CalcBmn(m-1,l+1)*pow(target.Scale(),l+1)*pow(1/scale,m-1)*trafo(l+1,m-1));
          
          for (int n = m; n < os; n++)
            {
              for (int l = m+1; l < os+ot-n; l++)
                trafom(l,n+1) = -1.0/sh.CalcAmn(m,n) * (sh.CalcAmn(m,l)*trafom(l+1,n)
                                                        -sh.CalcAmn(m,l-1)*trafom(l-1,n)
                                                        -sh.CalcAmn(m,n-1)*trafom(l,n-1));
              trafom(m,n+1) = pow(-1,m+n+1)*trafom(n+1,m);
            }

          
          double prod = 1;
          for (int i = 0; i <= os; i++, prod *= scale)
            trafom.Col(i) *= prod;
          prod = 1;
          for (int i = 0; i <= os+ot; i++, prod /= target.Scale())
            trafom.Row(i) *= prod;
          
          
          /*
          // if (m == 1)
          cout << "m = " << m << endl;
          cout << "trafom = " << trafom.Rows(m,ot+1).Cols(m,os+1) << endl;
          cout << trafom.Rows(ot+1,os+ot-m+1).Cols(m,os+1) << endl;          
          */
          
          for (int n = m; n <= os; n++)
            hv1(n) = sh.Coef(n,m);
          /*
          prod = 1;
          for (int i = 0; i <= os; i++, prod *= scale)
            hv1(i) *= prod;
          */
          
          hv2.Range(m,ot+1) = trafom.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
          /*
          prod = 1;
          for (int i = 0; i <= ot; i++, prod /= target.Scale())
            hv2(i) *= prod;
          */
          for (int n = m; n <= ot; n++)
            target.SH().Coef(n,m) = hv2(n);
          
          // cout << "m = " << m << ", hv1 = " << hv1 << ", hv2 = " << hv2 << endl;
          for (int n = m; n <= os; n++)
            hv1(n) = sh.Coef(n,-m);
          /*
          prod = 1;
          for (int i = 0; i <= os; i++, prod *= scale)
            hv1(i) *= prod;
          */
          hv2.Range(m,ot+1) = trafom.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
          /*
          prod = 1;
          for (int i = 0; i <= ot; i++, prod /= target.Scale())
            hv2(i) *= prod;
          */
          for (int n = m; n <= ot; n++)
            target.SH().Coef(n,-m) = hv2(n);
          // cout << "m = " << -m << ", hv1 = " << hv1 << ", hv2 = " << hv2 << endl;          
          trafo.Swap (trafom); //  = trafom; // swap ? 
        }
    }
  };

  // ***************** parameters ****************

  static int MPOrder (double rho_kappa)
  {
    return max (10, int(2*rho_kappa));
  }

  class SingularMLMultiPole
  {
    static Array<size_t> nodes_on_level;    
    
    struct Node
    {
      Vec<3> center;
      double r;
      int level;
      std::array<unique_ptr<Node>,8> childs;
      MultiPole<MPSingular> mp;

      Array<tuple<Vec<3>, Complex>> charges;
      Array<tuple<Vec<3>, Vec<3>, Complex>> dipoles;
      
      Node (Vec<3> acenter, double ar, int alevel, int order, double kappa)
        : center(acenter), r(ar), level(alevel), mp(MPOrder(ar*kappa), kappa, min(1.0, r*kappa))
      {
        // cout << "singml, add node, level = " << level << endl;
        if (level < nodes_on_level.Size())
          nodes_on_level[level]++;
      }


      void CreateChilds()
      {
        if (childs[0]) throw Exception("have already childs");
        for (int i = 0; i < 8; i++)
          {
            Vec<3> cc = center;
            cc(0) += (i&1) ? r/2 : -r/2;
            cc(1) += (i&2) ? r/2 : -r/2;
            cc(2) += (i&4) ? r/2 : -r/2;
            childs[i] = make_unique<Node> (cc, r/2, level+1, max(mp.SH().Order()/2, 8), mp.Kappa());
          }
      }
      

      void AddCharge (Vec<3> x, Complex c)
      {
        if (childs[0])
          {
            // directly send to childs:
            int childnum  = 0;
            if (x(0) > center(0)) childnum += 1;
            if (x(1) > center(1)) childnum += 2;
            if (x(2) > center(2)) childnum += 4;
            childs[childnum] -> AddCharge(x, c);
            return;
          }

        charges.Append( tuple{x,c} );

        if (r*mp.Kappa() < 1e-8) return;
        if (charges.Size() < 20 && r*mp.Kappa() < 1)
          return;

        CreateChilds();

        for (auto [x,c] : charges)
          AddCharge (x,c);
        for (auto [x,d,c] : dipoles)
          AddDipole (x,d,c);
        
        charges.SetSize0();
        dipoles.SetSize0();        
      }


      void AddDipole (Vec<3> x, Vec<3> d, Complex c)
      {
        if (childs[0])
          {
            // directly send to childs:

            int childnum  = 0;
            if (x(0) > center(0)) childnum += 1;
            if (x(1) > center(1)) childnum += 2;
            if (x(2) > center(2)) childnum += 4;
            childs[childnum] -> AddDipole(x, d, c);
            return;
          }

        dipoles.Append (tuple{x,d,c});

        if (dipoles.Size() < 20 || r < 1e-8)
          return;
        
        CreateChilds();

        for (auto [x,c] : charges)
          AddCharge (x,c);
        for (auto [x,d,c] : dipoles)
          AddDipole (x,d,c);

        charges.SetSize0();
        dipoles.SetSize0();        
      }

      
      Complex Evaluate(Vec<3> p) const
      {
        Complex sum = 0;
        if (childs[0])
          {
            for (auto & child : childs)
              sum += child->Evaluate(p);
            return sum;
          }

        for (auto [x,c] : charges)
          if (double rho = L2Norm(p-x); rho > 0)
            sum += c*(1/(4*M_PI))*exp(Complex(0,rho*mp.Kappa())) / rho;
        
        for (auto [x,d,c] : dipoles)
          if (double rho = L2Norm(p-x); rho > 0)
            {
              Vec<3> drhodp = 1.0/rho * (p-x);
              Complex dGdrho = c*(1/(4*M_PI))*exp(Complex(0,rho*mp.Kappa())) *
                (Complex(0, mp.Kappa())/rho - 1.0/sqr(rho));
              sum += dGdrho * InnerProduct(drhodp, d);
            }

        return sum;
      }


      void CalcMP()
      {
        mp.SH().Coefs() = 0.0;
        if (childs[0])
          {
            for (auto & child : childs)
              {
                child->CalcMP();
                child->mp.TransformAdd(mp, center-child->center);
              }
          }
        else
          {
            if (charges.Size()+dipoles.Size() == 0)
              {
                mp = MultiPole<MPSingular> (-1, mp.Kappa());
                return;
              }

            for (auto [x,c] : charges)
              mp.AddCharge (x-center,c);
            
            for (auto [x,d,c] : dipoles)
              mp.AddDipole (x-center, d, c);
          }
      }
      
      Complex EvaluateMP(Vec<3> p) const
      {
        if (charges.Size() || dipoles.Size())
          return Evaluate(p);
        
        if (L2Norm(p-center) > 3*r)
          return mp.Eval(p-center);
        
        if (!childs[0]) //  || level==1)
          return Evaluate(p);
          
        Complex sum = 0.0;
        for (auto & child : childs)
          sum += child->EvaluateMP(p);
        return sum;
      }
      
      
      void Print (ostream & ost) const
      {
        ost << "c = " << center << ", r = " << r << endl;
        // for (int i = 0; i < loc_pnts.Size(); i++)
        for (auto [x,c] : charges)
          ost << "xi = " << x << ", ci = " << c << endl;

        for (int i = 0; i < 8; i++)
          if (childs[i]) childs[i] -> Print (ost);
      }

      double Norm () const
      {
        double norm = L2Norm(mp.SH().Coefs());
        if (childs[0])
          for (auto & ch : childs)
            norm += ch->Norm();
        return norm;
      }
    };
    
    Node root;
    bool havemp = false;
    
  public:
    SingularMLMultiPole (Vec<3> center, double r, int order, double kappa)
      : root(center, r, 0, order, kappa)
    {
      nodes_on_level = 0;      
    }

    double Kappa() const { return root.mp.Kappa(); }
    
    void AddCharge(Vec<3> x, Complex c)
    {
      root.AddCharge(x, c);
    }

    void AddDipole(Vec<3> x, Vec<3> d, Complex c)
    {
      root.AddDipole(x, d, c);
    }

    void Print (ostream & ost) const
    {
      root.Print(ost);
    }

    double Norm() const
    {
      return root.Norm();
    }
    
    void CalcMP()
    {
      static Timer t("mptool compute singular MLMP"); RegionTimer rg(t);

      // nodes_on_level = 0;
      
      root.CalcMP();

      int maxlevel = 0;
      for (auto [i,num] : Enumerate(nodes_on_level))
        if (num > 0) maxlevel = i;

      for (int i = 0; i <= maxlevel; i++)
        cout << "sing " <<  i << ": " << nodes_on_level[i] << endl;
      
      havemp = true;
    }

    Complex Evaluate (Vec<3> p) const
    {
      if (havemp)
        return root.EvaluateMP(p);
      else
        return root.Evaluate(p);
    }

    friend class RegularMLMultiPole;
  };


  inline ostream & operator<< (ostream & ost, const SingularMLMultiPole & mlmp)
  {
    mlmp.Print(ost);
    return ost;
  }


  class RegularMLMultiPole
  {
    static Array<size_t> nodes_on_level;
    
    struct Node
    {
      Vec<3> center;
      double r;
      int level;
      std::array<unique_ptr<Node>,8> childs;
      MultiPole<MPRegular> mp;

      Array<const SingularMLMultiPole::Node*> singnodes;

      Node (Vec<3> acenter, double ar, int alevel, int order, double kappa)
        : center(acenter), r(ar), level(alevel),  mp(order, kappa, 1.0/min(1.0, kappa*r))
      {
        if (level < nodes_on_level.Size())
          nodes_on_level[level]++;
      }


      void CreateChilds()
      {
        if (childs[0]) throw Exception("have already childs");
        // create children nodes:
        for (int i = 0; i < 8; i++)
          {
            Vec<3> cc = center;
            cc(0) += (i&1) ? r/2 : -r/2;
            cc(1) += (i&2) ? r/2 : -r/2;
            cc(2) += (i&4) ? r/2 : -r/2;
            childs[i] = make_unique<Node> (cc, r/2, level+1, max(mp.SH().Order()/2, 8), mp.Kappa());
          }
      }
      
      void AddSingularNode (const SingularMLMultiPole::Node & singnode)
      {
        if (singnode.mp.SH().Order() < 0) return;
        if (L2Norm(singnode.mp.SH().Coefs()) == 0) return;
        if (level > 20)
          {
            singnodes.Append(&singnode);            
            return;
          }
        
        // static Timer t("AddSingularNode"); RegionTimer reg(t);
        
        Vec<3> dist = center-singnode.center;

        if (L2Norm(dist)*mp.Kappa() > (mp.Order()+singnode.mp.Order()))
          {
            if (singnode.mp.Order() > 2 * mp.Order() &&
                singnode.childs[0] &&
                singnode.childs[0]->mp.Order() < singnode.mp.Order())
              {
                for (auto & child : singnode.childs)
                  AddSingularNode (*child);
                return;
              }

            static Timer t("mptool transform Helmholtz-criterion"); RegionTimer r(t);
            singnode.mp.TransformAdd(mp, dist);
            return;
          }

        if ( ( (r+singnode.r) * mp.Kappa() < 10) && 
             (L2Norm(dist) > 3*(r + singnode.r)) )
          {
            static Timer t("mptool transform Laplace-criterion"); RegionTimer r(t);       
            
            // int truncorder = max(10, int(L2Norm(dist)*mp.Kappa()));
            int truncorder = 5;
            auto mptrunc = singnode.mp.Truncate(truncorder);
            // mptrunc.TransformAdd(mp, dist);

            MultiPole<MPRegular> tmp{mp.Truncate(truncorder)};
            mptrunc.Transform(tmp, dist);
            mp += tmp;
            return;
          }

        if ( singnode.childs[0]==nullptr
             // || singnode.mp.SH().Order() == singnode.childs[0]->mp.SH().Order()
             )
          {
            singnodes.Append(&singnode);
            return;
          }
        
        if (r > singnode.r)
          {
            if (!childs[0])
              CreateChilds();

            for (auto & ch : childs)
              ch -> AddSingularNode (singnode);
            return;
          }
        else
          {
            for (auto & childsing : singnode.childs)
              AddSingularNode (*childsing);
          }
      }

      void LocalizeExpansion()
      {
        if (mp.Order() > 20 && !childs[0])
          CreateChilds();

        if (childs[0])
          {
            for (auto & ch : childs)
              {
                if (L2Norm(mp.SH().Coefs()) > 0)
                  mp.TransformAdd (ch->mp, ch->center-center);
                ch->LocalizeExpansion();
              }
            mp = MultiPole<MPRegular>(-1, mp.Kappa());
            //mp.SH().Coefs()=0.0;
          }
      }
      
      Complex Evaluate (Vec<3> p) const
      {
        Complex sum = 0.0;
        if (childs[0])
          {
            int childnum = 0;
            if (p(0) > center(0)) childnum += 1;
            if (p(1) > center(1)) childnum += 2;
            if (p(2) > center(2)) childnum += 4;
            sum = childs[childnum]->Evaluate(p);
          }
        else
          sum = mp.Eval(p-center);

        static Timer t("mptool Eval singular nodes");
        t.Start();
        for (auto sn : singnodes)
          sum += sn->EvaluateMP(p);
        t.Stop();
        return sum;
      }

      double Norm() const
      {
        double norm = L2Norm(mp.SH().Coefs());
        if (childs[0])
          for (auto & ch : childs)
            norm += ch->Norm();
        return norm;
      }
    };
    
    Node root;
    shared_ptr<SingularMLMultiPole> singmp;
    
  public:
    RegularMLMultiPole (shared_ptr<SingularMLMultiPole> asingmp, Vec<3> center, double r, int order)
      : root(center, r, 1, order, asingmp->Kappa()), singmp(asingmp)
    {
      if (!singmp->havemp) throw Exception("first call Calc for singular MP");

      nodes_on_level = 0;
      {
        static Timer t("mptool compute regular MLMP"); RegionTimer rg(t);            
        root.AddSingularNode(singmp->root);
        // cout << "norm after S->R conversion: " << root.Norm() << endl;
      }

      
      int maxlevel = 0;
      for (auto [i,num] : Enumerate(nodes_on_level))
        if (num > 0) maxlevel = i;

      for (int i = 0; i <= maxlevel; i++)
        cout << "reg " << i << ": " << nodes_on_level[i] << endl;
      
      {
        static Timer t("mptool expand regular MLMP"); RegionTimer rg(t);                  
        root.LocalizeExpansion();
        // cout << "norm after local expansion: " << root.Norm() << endl;        
      }
    }

    Complex Evaluate (Vec<3> p) const
    {
      if (L2Norm(p-root.center) > root.r) return 0.0;
      return root.Evaluate(p);
    }
  };
  
  Array<size_t> RegularMLMultiPole::nodes_on_level(100);
  Array<size_t> SingularMLMultiPole::nodes_on_level(100);

  // ******************** Coefficient Functions *********************

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
    MultiPoleCF (int order, double kappa, Vec<3> acenter, double scale = 1)
      : CoefficientFunction(1, true), mp(order, kappa, scale), center(acenter) { }

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
      mp.Transform (target.MP(), target.Center()-center);
    }
    
  };

  
  class SingularMLMultiPoleCF : public CoefficientFunction
  {
    shared_ptr<SingularMLMultiPole> mlmp;
  public:
    SingularMLMultiPoleCF (Vec<3> center, double r, int order, double kappa)
      : CoefficientFunction(1, true), mlmp{make_shared<SingularMLMultiPole>(center, r, order, kappa)} { }
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("real eval not available"); }
    
    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const override
    {
      values(0) = mlmp->Evaluate(mip.GetPoint());
    }
    
    shared_ptr<SingularMLMultiPole> MLMP() { return mlmp; }
  };
  

  class RegularMLMultiPoleCF : public CoefficientFunction
  {
    RegularMLMultiPole mlmp;
  public:
    RegularMLMultiPoleCF (shared_ptr<SingularMLMultiPoleCF> asingmp, Vec<3> center, double r, int order)
      : CoefficientFunction(1, true), mlmp(asingmp->MLMP(), center, r, order) { } 
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("real eval not available"); }

    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const override
    {
      values(0) = mlmp.Evaluate(mip.GetPoint());
    }

    RegularMLMultiPole & MLMP() { return mlmp; }
  };

  
}
#endif
