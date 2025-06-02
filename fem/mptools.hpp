#ifndef FILE_MPTOOLS
#define FILE_MPTOOLS

/*
  tools for computing with spherical harmonics and multi-poles
 */


#include <bla.hpp>
#include <coefficient.hpp>
#include <recursive_pol.hpp>


namespace ngcomp
{
  class Region;
}

namespace ngfem
{


  template <typename entry_type = Complex>
  class NGS_DLL_HEADER SphericalHarmonics
  {
    int order;
    Vector<entry_type> coefs;

  public:
    SphericalHarmonics (int aorder)
      : order(aorder), coefs(sqr(order+1)) { coefs=0.0; }

    int Order() const { return order; }
    FlatVector<entry_type> Coefs() const { return coefs; }
    
    entry_type & Coef(int n, int m) { return coefs(n*(n+1) + m); }
    entry_type Coef(int n, int m) const { return coefs(n*(n+1) + m); }    

    auto CoefsN (int n) const
    {
      return coefs.RangeN(n*n, 2*n+1);
    }

    static std::tuple<double,double> Polar (Vec<3> x)
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
      return { theta, phi };
    }
    
    entry_type Eval (Vec<3> x) const
    {
      auto [theta, phi] = Polar(x);
      return Eval(theta, phi);
    }
  
    entry_type Eval (double theta, double phi) const;
    
    entry_type EvalOrder (int n, Vec<3> x) const
    {
      auto [theta, phi] = Polar (x);
      return EvalOrder(n, theta, phi);
    }
  
    entry_type EvalOrder (int n, double theta, double phi) const;

    void EvalOrders (Vec<3> x, FlatVector<entry_type> vals) const
    {
      auto [theta, phi] = Polar(x);
      return EvalOrders(theta, phi, vals);
    }
  
    void EvalOrders (double theta, double phi, FlatVector<entry_type> vals) const;
    
    void Calc (Vec<3> x, FlatVector<Complex> shapes);

    
    void RotateZ (double alpha);
    void RotateY (double alpha);

    
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
    
    // Nail A. Gumerov and Ramani Duraiswami book, formula (2.2.12)
    // add directional derivative divided by kappa to res, both multipoles need same scaling
    void DirectionalDiffAdd (Vec<3> d, SphericalHarmonics<entry_type> & res, double scale = 1) const;

  };


  // https://fortran-lang.discourse.group/t/looking-for-spherical-bessel-and-hankel-functions-of-first-and-second-kind-and-arbitrary-order/2308/2
  NGS_DLL_HEADER  
  void besseljs3d (int nterms, double z, double scale,
                   FlatVector<double> fjs, FlatVector<double> fjder);

  NGS_DLL_HEADER  
  void besseljs3d (int nterms, Complex z, double scale,
                   FlatVector<Complex> fjs, FlatVector<Complex> fjder);

  
  /*
  spherical bessel functions of first (the j_n) and second (the y_n) kind.
  
  j0(r) = sin(r)/r
  j1(r) = (sin(r)-r cos(r)) / r**2

  y0(r) = -cos(r)/r
  y1(r) = (-cos(r)-r*sin(r)) / r**2
  */
  NGS_DLL_HEADER    
  void SBESJY (double x, int lmax,
               FlatVector<double> j,
               FlatVector<double> y,
               FlatVector<double> jp,
               FlatVector<double> yp);

  
  
  template <typename T>
  void SphericalBessel (int n, double rho, double scale, T && values)
  {
    Vector<double> j(n+1), jp(n+1);
    besseljs3d (n, rho, scale,  j, jp);
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
    
    if (rho < 1e-100)
      {
        values = Complex(0);
        return;
      }
    Vector j(n+1), y(n+1), jp(n+1), yp(n+1);
    // SBESJY (rho, n, j, y, jp, yp);

    /*
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
    */

    
    // the bessel-evaluation with scale
    besseljs3d (n, rho, 1/scale,  j, jp);

    // Bessel y directly with the recurrence formula for (y, yp):
    double x = rho;
    double xinv = 1/x;
    y(0) = -xinv * cos(x);
    yp(0) = j(0)-xinv*y(0);

    double sl = 0;
    for (int l = 1; l <= n; l++)
      {
        y(l) = scale * (sl*y(l-1) - yp(l-1));
        sl += xinv;
        yp(l) = scale * y(l-1) - (sl+xinv)*y(l);
      }
    
    for (int i = 0; i <= n; i++)
      values(i) = Complex (j(i), y(i));
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
  
  


  template <typename RADIAL, typename entry_type=Complex>
  class NGS_DLL_HEADER MultiPole
  {
    SphericalHarmonics<entry_type> sh;
    double kappa;
    double scale;
  public:
    MultiPole (int aorder, double akappa, double ascale = 1) 
      : sh(aorder), kappa(akappa), scale(ascale) { }

    entry_type & Coef(int n, int m) { return sh.Coef(n,m); }
    auto & SH() { return sh; }
    const auto & SH() const { return sh; }
    double Kappa() const { return kappa; }
    double Scale() const { return scale; }    
    int Order() const { return sh.Order(); }
    
    MultiPole Truncate(int neworder) const
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
    
    entry_type Eval (Vec<3> x) const;
    entry_type EvalDirectionalDerivative (Vec<3> x, Vec<3> d) const;

    void AddCharge (Vec<3> x, entry_type c);
    void AddDipole (Vec<3> x, Vec<3> d, entry_type c);
    void AddCurrent (Vec<3> ap, Vec<3> ep, Complex j, int num=100);
    
    
    void ChangeScaleTo (double newscale)
    {
      double fac = scale/newscale;
      double prod = 1;
      for (int n = 0; n <= sh.Order(); n++, prod*= fac)
        sh.CoefsN(n) *= prod;
      scale = newscale;
    }

    Vector<double> Spectrum (bool scaled) const
    {
      Vector<double> spec(Order()+1);
      double fac = 1;
      for (int n = 0; n <= Order(); n++)
        {
          spec(n) = fac * L2Norm2(sh.CoefsN(n));
          if (!scaled) fac *= sqr(scale);
        }
      return spec;
    }

    
    template <typename TARGET>
    void Transform (MultiPole<TARGET,entry_type> & target, Vec<3> dist) const
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
        
      
      MultiPole<RADIAL,entry_type> tmp(*this);
      tmp.SH().RotateZ(phi);
      tmp.SH().RotateY(theta);

      tmp.ShiftZ(-len, target);
      
      target.SH().RotateY(-theta);
      target.SH().RotateZ(-phi);
    }
    
    template <typename TARGET>
    void TransformAdd (MultiPole<TARGET,entry_type> & target, Vec<3> dist) const
    {
      if (SH().Order() < 0) return;
      if (target.SH().Order() < 0) return;      
      
      MultiPole<TARGET,entry_type> tmp{target};
      Transform(tmp, dist);
      target.SH().Coefs() += tmp.SH().Coefs();
    }

    template <typename TARGET>
    void ShiftZ (double z, MultiPole<TARGET,entry_type> & target);
    
  };
  
  

  // ***************** parameters ****************

  static constexpr int MPOrder (double rho_kappa)
  {
    return max (20, int(2*rho_kappa));
  }
  static constexpr int maxdirect = 100;


  template <typename entry_type=Complex>
  class SingularMLMultiPole
  {
    static Array<size_t> nodes_on_level;    
    
    struct Node
    {
      Vec<3> center;
      double r;
      int level;
      std::array<unique_ptr<Node>,8> childs;
      MultiPole<MPSingular, entry_type> mp;

      Array<tuple<Vec<3>, entry_type>> charges;
      Array<tuple<Vec<3>, Vec<3>, entry_type>> dipoles;
      Array<tuple<Vec<3>, Vec<3>, Complex,int>> currents;
      int total_sources;
      
      Node (Vec<3> acenter, double ar, int alevel, int order, double kappa)
        : center(acenter), r(ar), level(alevel), mp(MPOrder(ar*kappa), kappa, min(1.0, r*kappa))
          // : center(acenter), r(ar), level(alevel), mp(MPOrder(ar*kappa), kappa, 1.0)
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
      

      void AddCharge (Vec<3> x, entry_type c)
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
        if (charges.Size() < maxdirect && r*mp.Kappa() < 1)
          return;

        CreateChilds();

        for (auto [x,c] : charges)
          AddCharge (x,c);
        for (auto [x,d,c] : dipoles)
          AddDipole (x,d,c);
        for (auto [sp,ep,j,num] : currents)
          AddCurrent (sp,ep,j,num);
        
        charges.SetSize0();
        dipoles.SetSize0();
        currents.SetSize0();
      }


      void AddDipole (Vec<3> x, Vec<3> d, entry_type c)
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

        if (dipoles.Size() < maxdirect || r < 1e-8)
          return;
        
        CreateChilds();

        for (auto [x,c] : charges)
          AddCharge (x,c);
        for (auto [x,d,c] : dipoles)
          AddDipole (x,d,c);
        for (auto [sp,ep,j,num] : currents)
          AddCurrent (sp,ep,j,num);

        charges.SetSize0();
        dipoles.SetSize0();        
        currents.SetSize0();        
      }

      void AddCurrent (Vec<3> sp, Vec<3> ep, Complex j, int num)
      {
        if (childs[0])
          {
            // split line and send to childs
            Array<double> split;
            split.Append(0);
            for (int i = 0; i < 3; i++)
              if (sp(i) < center(i) != ep(i) < center(i))
                split += (center(i)-sp(i)) / (ep(i)-sp(i));  // segment cuts i-th coordinate plane
            split.Append(1);
            BubbleSort(split);

            for (int i = 0; i < split.Size()-1; i++)
              if (split[i+1] > split[i])
                {
                  Vec<3> spi = sp + split[i]*(ep-sp);
                  Vec<3> epi = sp + split[i+1]*(ep-sp);
                  
                  Vec<3> x = 0.5*(spi+epi);
                  
                  int childnum  = 0;
                  if (x(0) > center(0)) childnum += 1;
                  if (x(1) > center(1)) childnum += 2;
                  if (x(2) > center(2)) childnum += 4;
                  childs[childnum] -> AddCurrent(spi, epi, j, num);
                }
            return;
          }
        
        currents.Append (tuple{sp,ep,j,num});

        // if (currents.Size() < maxdirect || r < 1e-8)
        if (currents.Size() < 4 || r < 1e-8)        
          return;
        
        CreateChilds();

        for (auto [x,c] : charges)
          AddCharge (x,c);
        for (auto [x,d,c] : dipoles)
          AddDipole (x,d,c);
        for (auto [sp,ep,j,num] : currents)
          AddCurrent (sp,ep,j,num);

        charges.SetSize0();
        dipoles.SetSize0();
        currents.SetSize0();
      }


      
      
      entry_type Evaluate(Vec<3> p) const
      {
        entry_type sum{0.0};
        if (childs[0])
          {
            for (auto & child : childs)
              sum += child->Evaluate(p);
            return sum;
          }

        for (auto [x,c] : charges)
          if (double rho = L2Norm(p-x); rho > 0)
            sum += (1/(4*M_PI))*exp(Complex(0,rho*mp.Kappa())) / rho * c;
        
        for (auto [x,d,c] : dipoles)
          if (double rho = L2Norm(p-x); rho > 0)
            {
              Vec<3> drhodp = 1.0/rho * (p-x);
              Complex dGdrho = (1/(4*M_PI))*exp(Complex(0,rho*mp.Kappa())) *
                (Complex(0, mp.Kappa())/rho - 1.0/sqr(rho));
              sum += dGdrho * InnerProduct(drhodp, d) * c;
            }

        for (auto [sp,ep,j,num] : currents)
          {
            // should use explizit formula instead ...
            
            Vec<3> tau = ep-sp;
            Vec<3> tau_num = 1.0/num *  tau;
            for (int i = 0; i < num; i++)
              {
                Vec<3> x = sp+(i+0.5)*tau_num;
                
                if (double rho = L2Norm(p-x); rho > 0)
                  {
                    Vec<3> drhodp = 1.0/rho * (p-x);
                    Complex dGdrho = (1/(4*M_PI))*exp(Complex(0,rho*mp.Kappa())) *
                      (Complex(0, mp.Kappa())/rho - 1.0/sqr(rho));

                    if constexpr (std::is_same<entry_type, Vec<3,Complex>>())
                      sum += j*dGdrho * Cross(drhodp, tau_num);
                  }
              }
          }
        
        return sum;
      }

      entry_type EvaluateDeriv(Vec<3> p, Vec<3> d) const
      {
        entry_type sum{0.0};
        if (childs[0])
        {
            for (auto & child : childs)
              sum += child->EvaluateDeriv(p, d);
            return sum;
        }

        if (dipoles.Size())
            throw Exception("EvaluateDeriv not implemented for dipoles in SingularMLMultiPole");

        for (auto [x,c] : charges)
          if (double rho = L2Norm(p-x); rho > 0)
          {
            Vec<3> drhodp = 1.0/rho * (p-x);
            Complex dGdrho = (1/(4*M_PI))*exp(Complex(0,rho*mp.Kappa())) *
            (Complex(0, mp.Kappa())/rho - 1.0/sqr(rho));
            sum += dGdrho * InnerProduct(drhodp, d) * c;
          }
        return sum;
      }

      void CalcTotalSources()
      {
        total_sources = charges.Size() + dipoles.Size();
        for (auto & child : childs)
          if (child)
            {
              child->CalcTotalSources();
              total_sources += child->total_sources;
            }
      }
      
      void CalcMP()
      {
        mp.SH().Coefs() = 0.0;
        if (childs[0])
          {
            if (total_sources < 1000)
              for (auto & child : childs)
                child->CalcMP();
            else
              ParallelFor (8, [&] (int nr)
                           {
                             childs[nr] -> CalcMP();
                           });

            
            for (auto & child : childs)            
              child->mp.TransformAdd(mp, center-child->center);
          }
        else
          {
            if (charges.Size()+dipoles.Size()+currents.Size() == 0)
              {
                mp = MultiPole<MPSingular,entry_type> (-1, mp.Kappa());
                return;
              }

            for (auto [x,c] : charges)
              mp.AddCharge (x-center,c);
            
            for (auto [x,d,c] : dipoles)
              mp.AddDipole (x-center, d, c);

            for (auto [sp,ep,j,num] : currents)
              mp.AddCurrent (sp-center, ep-center, j, num);
          }
      }
      
      entry_type EvaluateMP(Vec<3> p) const
      {
        if (charges.Size() || dipoles.Size())
          return Evaluate(p);
        
        if (L2Norm(p-center) > 3*r)
          return mp.Eval(p-center);
        
        if (!childs[0]) //  || level==1)
          return Evaluate(p);
          
        entry_type sum{0.0};
        for (auto & child : childs)
          sum += child->EvaluateMP(p);
        return sum;
      }

      entry_type EvaluateMPDeriv(Vec<3> p, Vec<3> d) const
      {
        // cout << "EvaluateMPDeriv Singular, p = " << p << ", d = " << d << ", r = " << r << ", center = " << center <<  endl;
        // cout << "Norm: " << L2Norm(p-center) << " > " << 3*r << endl;
        // cout << "charges.Size() = " << charges.Size() << ", dipoles.Size() = " << dipoles.Size() << endl;
        if (charges.Size() || dipoles.Size() || !childs[0])
          return EvaluateDeriv(p, d);

        if (L2Norm(p-center) > 3*r)
          return mp.EvalDirectionalDerivative(p-center, d);

        entry_type sum{0.0};
        for (auto & child : childs)
          sum += child->EvaluateMPDeriv(p, d);
        return sum;
      }

      void Print (ostream & ost, size_t childnr = -1) const
      {
        if (childnr == -1)
          ost << "c = " << center << ", r = " << r << ", level = " << level << endl;
        else
          ost << "c = " << center << ", r = " << r << ", level = " << level << ", childnr = " << childnr << endl;
        // for (int i = 0; i < loc_pnts.Size(); i++)
        for (auto [x,c] : charges)
          ost << "xi = " << x << ", ci = " << c << endl;
        for (auto [x,d,c] : dipoles)
          ost << "xi = " << x << ", di = " << d << ", ci = " << c << endl;

        for (int i = 0; i < 8; i++)
          if (childs[i]) childs[i] -> Print (ost, i);
      }

      double Norm () const
      {
        double norm = L2Norm(mp.SH().Coefs());
        if (childs[0])
          for (auto & ch : childs)
            norm += ch->Norm();
        return norm;
      }
      
      size_t NumCoefficients() const
      {
        size_t num = sqr(mp.SH().Order()+1);
        if (childs[0])
          for (auto & ch : childs)
            num += ch->NumCoefficients();
        return num;
      }
    };
    
    Node root;
    bool havemp = false;
    
  public:
    SingularMLMultiPole (Vec<3> center, double r, int order, double kappa)
      : root(center, r, 0, order, kappa)
    {
      nodes_on_level = 0;
      nodes_on_level[0] = 1;
    }

    double Kappa() const { return root.mp.Kappa(); }
    
    void AddCharge(Vec<3> x, Complex c)
    {
      root.AddCharge(x, c);
    }

    void AddDipole(Vec<3> x, Vec<3> d, entry_type c)
    {
      root.AddDipole(x, d, c);
    }

    void AddCurrent (Vec<3> sp, Vec<3> ep, Complex j, int num)
    {
      if constexpr (!std::is_same<entry_type, Vec<3,Complex>>())
        throw Exception("AddCurrent needs a singular vectorial MP");
      
      root.AddCurrent (sp, ep, j, num);
      /*
         // for testing
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
                root.AddDipole (sp+(i+0.5)*tau_num, cp, source);
            }
        }
      */
    }

    void Print (ostream & ost) const
    {
      root.Print(ost);
    }

    double Norm() const
    {
      return root.Norm();
    }

    size_t NumCoefficients() const
    {
      return root.NumCoefficients();
    }
    
    void CalcMP()
    {
      static Timer t("mptool compute singular MLMP"); RegionTimer rg(t);

      /*
      int maxlevel = 0;
      for (auto [i,num] : Enumerate(nodes_on_level))
        if (num > 0) maxlevel = i;

      for (int i = 0; i <= maxlevel; i++)
        cout << "sing " <<  i << ": " << nodes_on_level[i] << endl;
      */
      root.CalcTotalSources();
      root.CalcMP();
      
      havemp = true;
    }

    entry_type Evaluate (Vec<3> p) const
    {
      if (havemp)
        return root.EvaluateMP(p);
      else
        return root.Evaluate(p);
    }

    template <typename entry_type2>
    friend class RegularMLMultiPole;
  };


  template <typename entry_type>
  inline ostream & operator<< (ostream & ost, const SingularMLMultiPole<entry_type> & mlmp)
  {
    mlmp.Print(ost);
    return ost;
  }


  template <typename elem_type=Complex>
  class NGS_DLL_HEADER RegularMLMultiPole
  {
    static Array<size_t> nodes_on_level;
    
    struct Node
    {
      Vec<3> center;
      double r;
      int level;
      std::array<unique_ptr<Node>,8> childs;
      MultiPole<MPRegular,elem_type> mp;
      Array<Vec<3>> targets;
      int total_targets;

      Array<const typename SingularMLMultiPole<elem_type>::Node*> singnodes;

      Node (Vec<3> acenter, double ar, int alevel, int order, double kappa)
        : center(acenter), r(ar), level(alevel), mp(MPOrder(ar*kappa), kappa, 1.0/min(1.0, 0.25*r*kappa))
          // : center(acenter), r(ar), level(alevel), mp(MPOrder(ar*kappa), kappa, 1.0)
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

      void AddSingularNode (const typename SingularMLMultiPole<elem_type>::Node & singnode, bool allow_refine)
      {
        if (mp.SH().Order() < 0) return;
        if (singnode.mp.SH().Order() < 0) return;
        if (L2Norm(singnode.mp.SH().Coefs()) == 0) return;
        if (level > 20)
          {
            singnodes.Append(&singnode);            
            return;
          }
        
        // static Timer t("AddSingularNode"); RegionTimer reg(t);
        
        Vec<3> dist = center-singnode.center;

        // if (L2Norm(dist)*mp.Kappa() > (mp.Order()+singnode.mp.Order()))
        if (L2Norm(dist) > 2*(r + singnode.r))
          {
            if (singnode.mp.Order() > 2 * mp.Order() &&
                singnode.childs[0] &&
                singnode.childs[0]->mp.Order() < singnode.mp.Order())
              {
                for (auto & child : singnode.childs)
                  AddSingularNode (*child, allow_refine);
                return;
              }

            // static Timer t("mptool transform Helmholtz-criterion"); RegionTimer r(t);
            singnode.mp.TransformAdd(mp, dist);
            return;
          }


        if ( singnode.childs[0]==nullptr )
          {
            singnodes.Append(&singnode);
            return;
          }
        
        if (r > singnode.r)
          {
            if (allow_refine)
              {
                if (!childs[0])
                  CreateChilds();
                
                for (auto & ch : childs)
                  ch -> AddSingularNode (singnode, allow_refine);
              }
            else
              {
                if (total_targets < 1000)
                  {
                    for (auto & ch : childs)
                      if (ch)
                        ch -> AddSingularNode (singnode, allow_refine);
                  }
                else
                  ParallelFor (8, [&] (int nr)
                               {
                                 if (childs[nr])
                                   childs[nr] -> AddSingularNode (singnode, allow_refine);
                               });
                
                if (targets.Size())
                  singnodes.Append(&singnode);
              }
          }
        else
          {
            for (auto & childsing : singnode.childs)
              AddSingularNode (*childsing, allow_refine);
          }
      }

      void LocalizeExpansion(bool allow_refine)
      {
        if (allow_refine)
          if (mp.Order() > 20 && !childs[0])
            CreateChilds();

        if (childs[0])
          {
            for (auto & ch : childs)
              {
                if (L2Norm(mp.SH().Coefs()) > 0)
                  mp.TransformAdd (ch->mp, ch->center-center);
                ch->LocalizeExpansion(allow_refine);
              }
            mp = MultiPole<MPRegular,elem_type>(-1, mp.Kappa());
            //mp.SH().Coefs()=0.0;
          }
      }
      
      elem_type Evaluate (Vec<3> p) const
      {
        // *testout << "eval p = " << p << ", level = " << level << ", center = " << center <<  ", r = " << r << endl;
        elem_type sum{0.0};
        /*
        if (childs[0])
          {
            int childnum = 0;
            if (p(0) > center(0)) childnum += 1;
            if (p(1) > center(1)) childnum += 2;
            if (p(2) > center(2)) childnum += 4;
            sum = childs[childnum]->Evaluate(p);
          }
        */
        int childnum = 0;
        if (p(0) > center(0)) childnum += 1;
        if (p(1) > center(1)) childnum += 2;
        if (p(2) > center(2)) childnum += 4;
        if (childs[childnum])
          sum = childs[childnum]->Evaluate(p);
        else
          sum = mp.Eval(p-center);


        static Timer t("mptool direct evaluate"); RegionTimer r(t);
        for (auto sn : singnodes)
          sum += sn->EvaluateMP(p);

        return sum;
      }

      elem_type EvaluateDirectionalDerivative (Vec<3> p, Vec<3> d) const
      {
        elem_type sum{0.0};
        // cout << "EvaluateDirectionalDerivative RegularMLMP, r = " << r << ", level = " << level << ", center = " << center << endl;
        // cout << "Singnodes: " << singnodes.Size() << ", childs: " << childs[0] << endl;

        int childnum = 0;
        if (p(0) > center(0)) childnum += 1;
        if (p(1) > center(1)) childnum += 2;
        if (p(2) > center(2)) childnum += 4;
        if (childs[childnum])
          sum = childs[childnum]->EvaluateDirectionalDerivative(p, d);
        else
          sum = mp.EvalDirectionalDerivative(p-center, d);

        static Timer t("mptool direct evaluate deriv"); RegionTimer r(t);
        for (auto sn : singnodes)
          sum += sn->EvaluateMPDeriv(p, d);

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

      size_t NumCoefficients() const
      {
        size_t num = sqr(mp.SH().Order()+1);
        if (childs[0])
          for (auto & ch : childs)
            num += ch->NumCoefficients();
        return num;
      }
      
      void AddTarget (Vec<3> x)
      {
        if (childs[0])
          {
            // directly send to childs:
            int childnum  = 0;
            if (x(0) > center(0)) childnum += 1;
            if (x(1) > center(1)) childnum += 2;
            if (x(2) > center(2)) childnum += 4;
            childs[childnum] -> AddTarget( x );
            return;
          }

        targets.Append( x );

        if (r*mp.Kappa() < 1e-8) return;
        if (targets.Size() < maxdirect && r*mp.Kappa() < 1)
          return;

        CreateChilds();

        for (auto t : targets)
          AddTarget (t);
        targets.SetSize0();
      }

      void CalcTotalTargets()
      {
        total_targets = targets.Size();
        for (auto & child : childs)
          if (child)
            {
              child->CalcTotalTargets();
              total_targets += child->total_targets;
            }
      }

      void RemoveEmptyTrees()
      {
        for (auto & child : childs)
          if (child)
            {
              child->RemoveEmptyTrees();
              // if (child->total_targets == 0)
              // child = nullptr;
            }

        if (total_targets == 0)
          mp = MultiPole<MPRegular>(-1, mp.Kappa());
      }
      

      void Print (ostream & ost, size_t childnr = -1) const
      {
        if (childnr == -1)
          ost << "c = " << center << ", r = " << r << ", level = " << level << endl;
        else
          ost << "c = " << center << ", r = " << r << ", level = " << level << ", childnr = " << childnr << endl;
        for (auto x : targets)
          ost << "xi = " << x << endl;

        for (int i = 0; i < 8; i++)
          if (childs[i]) childs[i] -> Print (ost, i);
      }

    };
    
    Node root;
    shared_ptr<SingularMLMultiPole<elem_type>> singmp;
    
  public:
  RegularMLMultiPole (shared_ptr<SingularMLMultiPole<elem_type>> asingmp, Vec<3> center, double r, int order)
      : root(center, r, 0, order, asingmp->Kappa()), singmp(asingmp)
    {
      if (!singmp->havemp) throw Exception("first call Calc for singular MP");

      nodes_on_level = 0;
      nodes_on_level[0] = 1;
      {
        static Timer t("mptool compute regular MLMP"); RegionTimer rg(t);            
        root.AddSingularNode(singmp->root, true);
        // cout << "norm after S->R conversion: " << root.Norm() << endl;
      }


      /*
      int maxlevel = 0;
      for (auto [i,num] : Enumerate(nodes_on_level))
        if (num > 0) maxlevel = i;

      for (int i = 0; i <= maxlevel; i++)
        cout << "reg " << i << ": " << nodes_on_level[i] << endl;
      */
      
      {
        static Timer t("mptool expand regular MLMP"); RegionTimer rg(t);                  
        root.LocalizeExpansion(true);
        // cout << "norm after local expansion: " << root.Norm() << endl;        
      }
    }

    RegularMLMultiPole (Vec<3> center, double r, int order, double kappa)
      : root(center, r, 0, order, kappa)
    {
      nodes_on_level = 0;
      nodes_on_level[0] = 1;
    }
    
    void AddTarget (Vec<3> t)
    {
      root.AddTarget (t);
    }

    void CalcMP(shared_ptr<SingularMLMultiPole<elem_type>> asingmp)
    {
      singmp = asingmp;

      root.CalcTotalTargets();
      root.RemoveEmptyTrees();
        
      root.AddSingularNode(singmp->root, false);

      /*
      int maxlevel = 0;
      for (auto [i,num] : Enumerate(RegularMLMultiPole::nodes_on_level))
        if (num > 0) maxlevel = i;

      for (int i = 0; i <= maxlevel; i++)
        cout << "reg " << i << ": " << RegularMLMultiPole::nodes_on_level[i] << endl;
      */

      root.LocalizeExpansion(false);
    }

    void Print (ostream & ost) const
    {
      root.Print(ost);
    }

    double Norm() const
    {
      return root.Norm();
    }

    size_t NumCoefficients() const
    {
      return root.NumCoefficients();
    }

    elem_type Evaluate (Vec<3> p) const
    {
      // static Timer t("mptool Eval MLMP regular"); RegionTimer r(t);
      if (L2Norm(p-root.center) > root.r) return elem_type{0.0};
      return root.Evaluate(p);
    }

    elem_type EvaluateDirectionalDerivative (Vec<3> p, Vec<3> d) const
    {
        if (L2Norm(p-root.center) > root.r) return elem_type{0.0};
        return root.EvaluateDirectionalDerivative(p, d);
    }

  };

  template <typename elem_type>
  inline ostream & operator<< (ostream & ost, const RegularMLMultiPole<elem_type> & mlmp)
  {
    mlmp.Print(ost);
    // ost << "RegularMLMultiPole" << endl;
    return ost;
  }





  
  // ******************** Coefficient Functions *********************

  
  class SphericalHarmonicsCF : public CoefficientFunction
  {
    SphericalHarmonics<Complex> sh;
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



  template <typename entry_type> class RegularMLMultiPoleCF;
  

  template <typename RADIAL, typename entry_type=Complex>
  class MultiPoleCF : public CoefficientFunction
  {
    MultiPole<RADIAL, entry_type> mp;
    Vec<3> center;
  public:
    MultiPoleCF (int order, double kappa, Vec<3> acenter, double scale = 1)
      : CoefficientFunction(sizeof(entry_type)/sizeof(Complex), true), mp(order, kappa, scale), center(acenter) { }

    entry_type & Coef(int n, int m) { return mp.Coef(n,m); } 
    auto & SH() { return mp.SH(); }
    auto & MP() { return mp; }
    Vec<3> Center() const { return center; }
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("real eval not available"); }

    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const override
    {
      if constexpr (std::is_same<entry_type, Complex>())
        values(0) = mp.Eval(mip.GetPoint()-center);
      else
        values = mp.Eval(mip.GetPoint()-center);
    }

    template <typename TARGET>
    void ShiftZ (double z, MultiPole<TARGET, entry_type> & target) { mp.ShiftZ(z, target); }

    using CoefficientFunction::Transform;        
    template <typename TARGET>
    void Transform (MultiPoleCF<TARGET, entry_type> & target)
    {
      mp.Transform (target.MP(), target.Center()-center);
    }
  };

  template <typename entry_type>
  class SingularMLMultiPoleCF : public CoefficientFunction
  {
    shared_ptr<SingularMLMultiPole<entry_type>> mlmp;
  public:
    SingularMLMultiPoleCF (Vec<3> center, double r, int order, double kappa)
      : CoefficientFunction(sizeof(entry_type)/sizeof(Complex), true), mlmp{make_shared<SingularMLMultiPole<entry_type>>(center, r, order, kappa)} { }
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("real eval not available"); }
    
    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const override
    {
      // values(0) = mlmp->Evaluate(mip.GetPoint());

      if constexpr (std::is_same<entry_type, Complex>())
        values(0) = mlmp->Evaluate(mip.GetPoint());
      else
        values = mlmp->Evaluate(mip.GetPoint());

      
    }
    
    shared_ptr<SingularMLMultiPole<entry_type>> MLMP() const { return mlmp; }
    shared_ptr<RegularMLMultiPoleCF<entry_type>> CreateRegularExpansion(Vec<3> center, double r) const;
  };
  

  template <typename entry_type>
  class RegularMLMultiPoleCF : public CoefficientFunction
  {
    shared_ptr<RegularMLMultiPole<entry_type>> mlmp;
  public:
    RegularMLMultiPoleCF (shared_ptr<SingularMLMultiPoleCF<entry_type>> asingmp, Vec<3> center, double r, int order)
      : CoefficientFunction(sizeof(entry_type)/sizeof(Complex), true), mlmp{make_shared<RegularMLMultiPole<entry_type>>(asingmp->MLMP(), center, r, order)} { } 
    RegularMLMultiPoleCF (shared_ptr<SingularMLMultiPole<entry_type>> asingmp, Vec<3> center, double r, int order)
      : CoefficientFunction(sizeof(entry_type)/sizeof(Complex), true), mlmp{make_shared<RegularMLMultiPole<entry_type>>(asingmp, center, r, order)} { } 
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    { throw Exception("real eval not available"); }

    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const override
    {
      // values(0) = mlmp->Evaluate(mip.GetPoint());

      if constexpr (std::is_same<entry_type, Complex>())
        values(0) = mlmp->Evaluate(mip.GetPoint());
      else
        values = mlmp->Evaluate(mip.GetPoint());
    }

    shared_ptr<RegularMLMultiPole<entry_type>> MLMP() { return mlmp; }
  };

  
}
#endif
