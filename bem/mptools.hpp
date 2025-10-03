#ifndef FILE_MPTOOLS
#define FILE_MPTOOLS

/*
  tools for computing with spherical harmonics and multi-poles
 */


#include <bla.hpp>
#include <coefficient.hpp>
#include <recursive_pol.hpp>


namespace ngsbem
{
  using namespace ngfem;
  
  template<typename T>
  constexpr int VecLength = 1;  // Default: Complex has length 1

  template<int N>
  constexpr int VecLength<Vec<N, Complex>> = N;  // Specialization: Vec<N,Complex> has length N



  constexpr int FMM_SW = 4;

  

  // ************************ SIMD - creation (should end up in simd.hpp) ************* 


  template <int S, typename T, int SW>
  Vec<S,T> HSum (Vec<S,SIMD<T,SW>> v)
  {
    Vec<S,T> res;
    for (int i = 0; i < S; i++)
      res(i) = HSum(v(i));
    // Iterate<S> ([&](auto i) {
    // res.HTData().template Elem<i.value>() = HSum(v.HTData().template Elem<i.value>());
    // });
    return res;
  }
  

  class NGS_DLL_HEADER PrecomputedSqrts
  {
  public:
    Array<double> sqrt_int;
    // Array<double> inv_sqrt_int;
    Array<double> sqrt_n_np1;    // sqrt(n*(n+1))
    Array<double> inv_sqrt_2np1_2np3;  // 1/sqrt( (2n+1)*(2n+3) )
    
    PrecomputedSqrts();
  };
  
  extern NGS_DLL_HEADER PrecomputedSqrts presqrt;
  


  class FMM_Parameters
  {
  public:
    int maxdirect = 100;
    int minorder = 20;    // order = minorder + 2 kappa r 
  };


  
  
  inline std::tuple<double, double, double> SphericalCoordinates(Vec<3> dist){
    double len, theta, phi;
    len = L2Norm(dist);
    if (len < 1e-30)
      theta = 0;
    else
      theta = acos (dist(2) / len);
    if (sqr(dist(0))+sqr(dist(1)) < 1e-30)
      phi = 0;
    else
      phi = atan2(dist(1), dist(0));
    return {len, theta, phi};
  }


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


    void FlipZ ();    
    void RotateZ (double alpha);

    template <typename FUNC>
    void RotateZ (double alpha, FUNC func) const
    {
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
            func(ii, conj(exp_imalpha(-m)));
          for (int m = 0; m <= n; m++, ii++)
            func(ii, exp_imalpha(m));
        };
    };

    template <typename FUNC>
    void RotateZFlip (double alpha, bool flip, FUNC func) const
    {
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

      auto FlipFactor = [] (int n, int m, bool flip)->double
      {
        if (flip)
          return ((n-m)%2) == 1 ? -1 : 1;
        return 1.0;
      };
      
      for (int n = 0; n <= order; n++)
        {
          for (int m = -n; m < 0; m++, ii++)
            func(ii, FlipFactor(n,m,flip)*conj(exp_imalpha(-m)));
          for (int m = 0; m <= n; m++, ii++)
            func(ii, FlipFactor(n,m,flip)*exp_imalpha(m));
        };
    };

    
    
    void RotateY (double alpha, bool parallel = false);

    
    
    static double CalcAmn (int m, int n)
    {
      if (m < 0) m=-m;
      if (n < m) return 0;

      if (2*n+1 < presqrt.sqrt_int.Size())
        return presqrt.sqrt_int[n+1+m]*presqrt.sqrt_int[n+1-m] * presqrt.inv_sqrt_2np1_2np3[n];
      else
        return sqrt( (n+1.0+m)*(n+1.0-m) / ( (2*n+1)*(2*n+3) ));
    }
  
    static double CalcBmn (int m, int n)
    {
      double sgn = (m >= 0) ? 1 : -1;
      if ( (m >= n) || (-m > n) ) return 0;
      if (n <= presqrt.inv_sqrt_2np1_2np3.Size())
        return sgn * presqrt.sqrt_n_np1[n-m-1] * presqrt.inv_sqrt_2np1_2np3[n-1];
      else
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
                   SliceVector<double> fjs, SliceVector<double> fjder = FlatVector<double>(0, nullptr));

  NGS_DLL_HEADER  
  void besseljs3d (int nterms, Complex z, double scale,
                   SliceVector<Complex> fjs, SliceVector<Complex> fjder = FlatVector<Complex>(0, nullptr));

  
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
    besseljs3d (n, rho, scale,  values);
    /*
    Vector<double> j(n+1), jp(n+1);
    besseljs3d (n, rho, scale,  j, jp);
    values = j;
    */
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
  class Singular
  {
  public:
    template <typename T>
    static void Eval (int order, double r, double scale, T && values)
    {
      SphericalHankel1(order, r, scale,  values);
    }

    template <typename T>
    static void Eval (int order, double kappa, double r, double rtyp, T && values)
    {
      double scale = Scale(kappa, rtyp);
      SphericalHankel1(order, r*kappa, scale,  values);
    }

    static double Scale (double kappa, double rtyp)
    {
      // return min(1.0, rtyp*kappa);
      return min(1.0, 0.5*rtyp*kappa);      
    }
  };


  
  // jn
  class Regular
  {
  public:
    template <typename T>
    static void Eval (int order, double r, double scale, T && values)
    {
      SphericalBessel (order, r, 1.0/scale, values);
    }

    template <typename T>
    static void Eval (int order, double kappa, double r, double rtyp, T && values)
    {
      double scale = Scale(kappa, rtyp);
      SphericalBessel (order, r*kappa, 1.0/scale, values);      
    }

    static double Scale (double kappa, double rtyp)
    {
      // return 1.0/ min(1.0, 0.25*rtyp*kappa);
      return 1.0/ min(1.0, 0.5*rtyp*kappa);
    }
    
  };
  
  


  template <typename RADIAL, typename entry_type=Complex>
  class NGS_DLL_HEADER SphericalExpansion
  {
    SphericalHarmonics<entry_type> sh;
    double kappa;
    double rtyp;
  public:

    SphericalExpansion (int aorder, double akappa, double artyp) 
    : sh(aorder), kappa(akappa), rtyp(artyp) { }

  
    entry_type & Coef(int n, int m) { return sh.Coef(n,m); }
    auto & SH() { return sh; }
    const auto & SH() const { return sh; }
    double Kappa() const { return kappa; }
    double Scale() const { return RADIAL::Scale(kappa, rtyp); }
    double RTyp() const { return rtyp; }
    int Order() const { return sh.Order(); }
    
    SphericalExpansion Truncate(int neworder) const
    {
      if (neworder > sh.Order()) neworder=sh.Order();
      SphericalExpansion nmp(neworder, kappa, rtyp);
      nmp.sh.Coefs() = sh.Coefs().Range(sqr(neworder+1));
      return nmp;
    }

    SphericalExpansion & operator+= (const SphericalExpansion & mp2)
    {
      size_t commonsize = min(SH().Coefs().Size(), mp2.SH().Coefs().Size());
      SH().Coefs().Range(commonsize) += mp2.SH().Coefs().Range(commonsize);
      return *this;
    }
    
    entry_type Eval (Vec<3> x) const;
    entry_type EvalDirectionalDerivative (Vec<3> x, Vec<3> d) const;

    void AddCharge (Vec<3> x, entry_type c);
    void AddDipole (Vec<3> x, Vec<3> dir, entry_type c);
    void AddChargeDipole (Vec<3> x, entry_type c, Vec<3> dir, entry_type c2)
    {
      // TODO: add them at once
      AddCharge (x, c);
      AddDipole (x, dir, c2);
    }
    
    void AddPlaneWave (Vec<3> d, entry_type c);    
    void AddCurrent (Vec<3> ap, Vec<3> ep, Complex j, int num=100);
    

    void ChangeRTypTo (double new_rtyp)
    {
      double fac = RADIAL::Scale(kappa, rtyp) / RADIAL::Scale(kappa, new_rtyp);
      double prod = 1;
      for (int n = 0; n <= sh.Order(); n++, prod*= fac)
        sh.CoefsN(n) *= prod;
      rtyp = new_rtyp;
    }
    

    Vector<double> Spectrum (bool scaled) const
    {
      Vector<double> spec(Order()+1);
      double fac = 1;
      for (int n = 0; n <= Order(); n++)
        {
          spec(n) = fac * L2Norm2(sh.CoefsN(n));
          if (!scaled) fac *= sqr(Scale());
        }
      return spec;
    }

    
    template <typename TARGET>
    void Transform (SphericalExpansion<TARGET,entry_type> & target, Vec<3> dist) const
    {
      if (target.SH().Order() < 0) return;
      if (SH().Order() < 0)
        {
          target.SH().Coefs() = 0.0;
          return;
        }
      
      // static Timer t("mptool Transform "+ToString(typeid(RADIAL).name())+ToString(typeid(TARGET).name()));      
      // RegionTimer reg(t);
      
      auto [len, theta, phi] = SphericalCoordinates(dist);
        
      
      // SphericalExpansion<RADIAL,entry_type> tmp{*this};
      SphericalExpansion<RADIAL,entry_type> tmp(Order(), kappa, rtyp);
      tmp.SH().Coefs() = SH().Coefs();
      
      tmp.SH().RotateZ(phi);
      tmp.SH().RotateY(theta);

      tmp.ShiftZ(-len, target);
      
      target.SH().RotateY(-theta);
      target.SH().RotateZ(-phi);
    }
    
    template <typename TARGET>
    void TransformAdd (SphericalExpansion<TARGET,entry_type> & target, Vec<3> dist, bool atomic = false) const
    {
      if (SH().Order() < 0) return;
      if (target.SH().Order() < 0) return;      
      
      SphericalExpansion<TARGET,entry_type> tmp{target};
      Transform(tmp, dist);
      if (!atomic)
        target.SH().Coefs() += tmp.SH().Coefs();
      else
        for (int j = 0; j < target.SH().Coefs().Size(); j++)        
          AtomicAdd(target.SH().Coefs()[j], tmp.SH().Coefs()[j]);        
    }

    template <typename TARGET>
    void ShiftZ (double z, SphericalExpansion<TARGET,entry_type> & target);

    
    template <typename TARGET>
    void In2Out (SphericalExpansion<TARGET,entry_type> & target, double r) const
    {
      Vector<Complex> rad(Order()+1);
      Vector<Complex> radout(target.Order()+1);      
      RADIAL::Eval(Order(), kappa, r, RTyp(), rad);
      TARGET::Eval(target.Order(), kappa, r, target.RTyp(), radout);
      target.SH().Coefs() = 0;
      for (int j = 0; j <= std::min(Order(), target.Order()); j++)
        target.SH().CoefsN(j) = rad(j)/radout(j) * SH().CoefsN(j);
    }
  };
  
  

  // ***************** parameters ****************

  /*
  static constexpr int MPOrder (double rho_kappa)
  {
    // return max (20, int(2*rho_kappa));
    return 20+int(2*rho_kappa);
  }
  static constexpr int maxdirect = 100;
  */


  template <typename SCAL, auto S>
  inline auto VecVector2Matrix (FlatVector<Vec<S,SCAL>> vec)
  {
    return FlatMatrixFixWidth<S,SCAL> (vec.Size(), vec.Data()->Data());
  }
  
  inline auto VecVector2Matrix (FlatVector<Complex> vec)
  {
    return FlatMatrixFixWidth<1,Complex> (vec.Size(), vec.Data());    
  }


  template <typename entry_type=Complex>
  class SingularMLExpansion
  {
    using simd_entry_type = decltype(MakeSimd(declval<std::array<entry_type,FMM_SW>>()));
    static Array<size_t> nodes_on_level;    
    
    struct RecordingSS
    {
      const SphericalExpansion<Singular,entry_type> * mp_source;
      SphericalExpansion<Singular,entry_type> * mp_target;
      Vec<3> dist;
      double len, theta, phi;
      bool flipz;
    public:
      RecordingSS() = default;
      RecordingSS (const SphericalExpansion<Singular,entry_type> * amp_source,
                   SphericalExpansion<Singular,entry_type> * amp_target,
                   Vec<3> adist)
        : mp_source(amp_source), mp_target(amp_target), dist(adist)
      {
        std::tie(len, theta, phi) = SphericalCoordinates(adist);
        // flipz = false;       
        flipz = theta > M_PI/2;
        if (flipz) theta = M_PI-theta;
      }
    };


    static void ProcessBatchSS(FlatArray<RecordingSS*> batch, double len, double theta) {
      constexpr int vec_length = VecLength<entry_type>;
      int batch_size = batch.Size();
      int N = batch_size * vec_length;
      // *testout << "Processing batch of size " << batch.Size() << ", with N = " << N << ", vec_length = " << vec_length << ", Type: " << typeid(entry_type).name() << ", len = " << len << ", theta = " << theta << endl;

      if (N <= 1 || batch_size <= 1) {
        for (auto* rec : batch) {
          rec->mp_source->TransformAdd(*rec->mp_target, rec->dist, true);
        }
      }
      else if (N <= 3) {
        ProcessVectorizedBatchSS<3, vec_length>(batch, len, theta);
      }
      else if (N <= 4) {
        ProcessVectorizedBatchSS<4, vec_length>(batch, len, theta);
      }
      else if (N <= 6) {
        ProcessVectorizedBatchSS<6, vec_length>(batch, len, theta);
      }
      else if (N <= 12) {
        ProcessVectorizedBatchSS<12, vec_length>(batch, len, theta);
      }
      else if (N <= 24) {
        ProcessVectorizedBatchSS<24, vec_length>(batch, len, theta);
      }
      else if (N <= 48) {
        ProcessVectorizedBatchSS<48, vec_length>(batch, len, theta);
      }
      else if (N <= 96) {
        ProcessVectorizedBatchSS<96, vec_length>(batch, len, theta);
      }
      else if (N <= 192) {
        ProcessVectorizedBatchSS<192, vec_length>(batch, len, theta);
      }
      else {
        // Split large batches
        ProcessBatchSS(batch.Range(0, 192 / vec_length), len, theta);
        ProcessBatchSS(batch.Range(192 / vec_length, batch_size), len, theta);
      }
    }

    template<int N, int vec_length>
    static void ProcessVectorizedBatchSS(FlatArray<RecordingSS*> batch, double len, double theta) {

      // *testout << "Processing vectorized S->S batch of size " << batch.Size() << ", with N = " << N << ", vec_length = " << vec_length << ", len = " << len << ", theta = " << theta << endl;
      double kappa = batch[0]->mp_source->Kappa();
      int so = batch[0]->mp_source->Order();
      int to = batch[0]->mp_target->Order();
      SphericalExpansion<Singular, Vec<N,Complex>> vec_source(so, kappa, batch[0]->mp_source->RTyp());
      SphericalExpansion<Singular, Vec<N,Complex>> vec_target(to, kappa, batch[0]->mp_target->RTyp());

      // Copy multipoles into vectorized multipole
      for (int i = 0; i < batch.Size(); i++)
        {
          auto source_i = VecVector2Matrix (batch[i]->mp_source->SH().Coefs());
          auto source_mati = VecVector2Matrix (vec_source.SH().Coefs()).Cols(i*vec_length, (i+1)*vec_length);
          batch[i]->mp_source->SH().RotateZFlip(batch[i]->phi, batch[i]->flipz,
                                            [source_i, source_mati] (size_t ii, Complex factor)
                                            {
                                              source_mati.Row(ii) = factor * source_i.Row(ii);
                                            });
        }
      
      vec_source.SH().RotateY(theta, vec_source.SH().Order() >= 100);
      vec_source.ShiftZ(-len, vec_target);
      vec_target.SH().RotateY(-theta, vec_target.SH().Order() >= 100);

      // Copy vectorized multipole into individual multipoles
      for (int i = 0; i < batch.Size(); i++)
        {
          auto source_mati = VecVector2Matrix (vec_target.SH().Coefs()).Cols(i*vec_length, (i+1)*vec_length);
          auto target_mati = VecVector2Matrix (batch[i]->mp_target->SH().Coefs());
          batch[i]->mp_target->SH().RotateZFlip(-batch[i]->phi, batch[i]->flipz,
                                      [source_mati, target_mati] (size_t ii, Complex factor)
                                      {
                                        AtomicAdd (target_mati.Row(ii), factor * source_mati.Row(ii));
                                      });
      }
    }

    struct Node
    {
      Vec<3> center;
      double r;
      int level;
      std::array<unique_ptr<Node>,8> childs;
      SphericalExpansion<Singular, entry_type> mp;

      Array<tuple<Vec<3>, entry_type>> charges;
      Array<tuple<Vec<3>, Vec<3>, entry_type>> dipoles;
      Array<tuple<Vec<3>, entry_type, Vec<3>, entry_type>> chargedipoles;
      Array<tuple<Vec<3>, Vec<3>, Complex,int>> currents;

      using simd_entry_type = decltype(MakeSimd(declval<std::array<entry_type,FMM_SW>>()));      
      Array<tuple<Vec<3,SIMD<double,FMM_SW>>, simd_entry_type>> simd_charges;
      Array<tuple<Vec<3,SIMD<double,FMM_SW>>, Vec<3,SIMD<double,FMM_SW>>, simd_entry_type>> simd_dipoles;
      Array<tuple<Vec<3,SIMD<double,FMM_SW>>, simd_entry_type,
                  Vec<3,SIMD<double,FMM_SW>>, simd_entry_type>> simd_chargedipoles;
      
      int total_sources;
      const FMM_Parameters & fmm_params;
      std::mutex node_mutex;
      atomic<bool> have_childs{false};
      
      Node (Vec<3> acenter, double ar, int alevel, double akappa, const FMM_Parameters & afmm_params)
      // : center(acenter), r(ar), level(alevel), mp(MPOrder(ar*akappa), akappa, ar), fmm_params(afmm_params)
        : center(acenter), r(ar), level(alevel), mp(afmm_params.minorder+2*ar*akappa, akappa, ar), fmm_params(afmm_params)
      {
        if (level < nodes_on_level.Size())
          nodes_on_level[level]++;
      }

      int GetChildNum (Vec<3> x) const
      {
        int childnum  = 0;
        if (x(0) > center(0)) childnum += 1;
        if (x(1) > center(1)) childnum += 2;
        if (x(2) > center(2)) childnum += 4;
        return childnum;
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
            childs[i] = make_unique<Node> (cc, r/2, level+1, mp.Kappa(), fmm_params);
          }
        have_childs = true;
      }
      

      void SendSourcesToChilds()
      {
        CreateChilds();

        for (auto [x,c] : charges)
          AddCharge (x,c);
        for (auto [x,d,c] : dipoles)
          AddDipole (x,d,c);
        for (auto [x,c,d,c2] : chargedipoles)
          AddChargeDipole (x,c,d,c2);
        for (auto [sp,ep,j,num] : currents)
          AddCurrent (sp,ep,j,num);
        
        charges.DeleteAll();
        dipoles.DeleteAll();
        chargedipoles.DeleteAll();        
        currents.DeleteAll();
      }

      
      void AddCharge (Vec<3> x, entry_type c)
      {
        if (have_childs) // quick check without locking 
          {
            // directly send to childs:
            int childnum = GetChildNum(x);
            childs[childnum] -> AddCharge(x, c);
            return;
          }

        lock_guard<mutex> guard(node_mutex);

        if (have_childs) // test again after locking 
          {
            int childnum  = GetChildNum(x);
            childs[childnum] -> AddCharge(x, c);
            return;
          }

        charges.Append( tuple{x,c} );

        // if (r*mp.Kappa() < 1e-8) return;
        if (level > 20) return;
        if (charges.Size() < fmm_params.maxdirect && r*mp.Kappa() < 5)
          return;
        
        SendSourcesToChilds();
      }


      void AddDipole (Vec<3> x, Vec<3> d, entry_type c)
      {
        if (have_childs)
          {
            // directly send to childs:
            int childnum = GetChildNum(x);
            childs[childnum] -> AddDipole(x, d, c);
            return;
          }

        lock_guard<mutex> guard(node_mutex);

        if (have_childs)
          {
            // directly send to childs:
            int childnum = GetChildNum(x);
            childs[childnum] -> AddDipole(x, d, c);
            return;
          }
        
        dipoles.Append (tuple{x,d,c});
        
        if (level > 20) return;
        if (dipoles.Size() < fmm_params.maxdirect)
          return;

        SendSourcesToChilds();
      }


      void AddChargeDipole (Vec<3> x, entry_type c, Vec<3> dir, entry_type c2)
      {
        if (have_childs)
          {
            // directly send to childs:
            int childnum = GetChildNum(x);
            childs[childnum] -> AddChargeDipole(x, c, dir, c2);
            return;
          }

        lock_guard<mutex> guard(node_mutex);

        if (have_childs)
          {
            // directly send to childs:
            int childnum = GetChildNum(x);
            childs[childnum] -> AddChargeDipole(x, c, dir, c2);
            return;
          }
        
        chargedipoles.Append (tuple{x,c,dir,c2});

        if (chargedipoles.Size() < fmm_params.maxdirect || r < 1e-8)
          return;

        SendSourcesToChilds();

        /*
        AddCharge (x, c);
        AddDipole (x, dir, c2);
        */
      }

      
      // not parallel yet
      void AddCurrent (Vec<3> sp, Vec<3> ep, Complex j, int num)
      {
        if (childs[0])
          {
            // split line and send to childs
            Array<double> split;
            split.Append(0);
            for (int i = 0; i < 3; i++)
              if ((sp(i) < center(i)) != (ep(i) < center(i)))
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

        SendSourcesToChilds();
        /*
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
        */
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

        if (simd_charges.Size())
          {
            // static Timer t("mptool singmp, evaluate, simd charges"); RegionTimer r(t);
            // t.AddFlops (charges.Size());
            
            simd_entry_type vsum{0.0};
            if (mp.Kappa() < 1e-12)
              {
                for (auto [x,c] : simd_charges)
                  {
                    auto rho = L2Norm(p-x);
                    auto kernel = 1/(4*M_PI)/rho;
                    kernel = If(rho > 0.0, kernel, SIMD<double,FMM_SW>(0.0));
                    vsum += kernel * c;

                    /*
                    auto rho2 = L2Norm2(p-x);
                    auto kernel = (1/(4*M_PI)) * rsqrt(rho2);
                    kernel = If(rho2 > 0.0, kernel, SIMD<double,FMM_SW>(0.0));
                    vsum += kernel * c;
                    */
                  }
              }
            else if (mp.Kappa() < 1e-8)
              for (auto [x,c] : simd_charges)
                {
                  auto rho = L2Norm(p-x);
                  auto kernel = (1/(4*M_PI))*SIMD<Complex,FMM_SW> (1,rho*mp.Kappa()) / rho;
                  kernel = If(rho > 0.0, kernel, SIMD<Complex,FMM_SW>(0.0));
                  vsum += kernel * c;
                }
            else
              for (auto [x,c] : simd_charges)
                {
                  auto rho = L2Norm(p-x);
                  auto [si,co] = sincos(rho*mp.Kappa());
                  auto kernel = (1/(4*M_PI))*SIMD<Complex,FMM_SW>(co,si) / rho;
                  kernel = If(rho > 0.0, kernel, SIMD<Complex,FMM_SW>(0.0));
                  vsum += kernel * c;
                }
            
            sum += HSum(vsum);
          }
        else
          {
            if (mp.Kappa() < 1e-8)
              {
                for (auto [x,c] : charges)
                  if (double rho = L2Norm(p-x); rho > 0)
                    sum += (1/(4*M_PI))*Complex(1,rho*mp.Kappa()) / rho * c;
              }
            else
              for (auto [x,c] : charges)
                if (double rho = L2Norm(p-x); rho > 0)
                  sum += (1/(4*M_PI))*exp(Complex(0,rho*mp.Kappa())) / rho * c;
          }

        if (simd_dipoles.Size())
          {
            // static Timer t("mptool singmp, evaluate, simd dipoles"); RegionTimer r(t);
            
            simd_entry_type vsum{0.0};
            for (auto [x,d,c] : simd_dipoles)
              {
                auto rho = L2Norm(p-x);
                auto drhodp = (1.0/rho) * (p-x);
                auto [si,co] = sincos(rho*mp.Kappa());
                auto dGdrho = (1/(4*M_PI))*SIMD<Complex,FMM_SW>(co,si) * 
                  (-1.0/(rho*rho) + SIMD<Complex,FMM_SW>(0, mp.Kappa())/rho);
                auto kernel = dGdrho * InnerProduct(drhodp, d);
                kernel = If(rho > 0.0, kernel, SIMD<Complex,FMM_SW>(0.0));
                vsum += kernel * c;
              }
            sum += HSum(vsum);
          }
        else
          {
            for (auto [x,d,c] : dipoles)
              if (double rho = L2Norm(p-x); rho > 0)
                {
                  Vec<3> drhodp = 1.0/rho * (p-x);
                  Complex dGdrho = (1/(4*M_PI))*exp(Complex(0,rho*mp.Kappa())) *
                    (Complex(0, mp.Kappa())/rho - 1.0/sqr(rho));
                  sum += dGdrho * InnerProduct(drhodp, d) * c;
                }
          }
      
      
      
      if (simd_chargedipoles.Size())
        {
          // static Timer t("mptool singmp, evaluate, simd chargedipoles"); RegionTimer r(t);
          // t.AddFlops (simd_chargedipoles.Size()*FMM_SW);
          
          simd_entry_type vsum{0.0};
          for (auto [x,c,d,c2] : simd_chargedipoles)
            {
              auto rho = L2Norm(p-x);
              auto rhokappa = rho*mp.Kappa();
              auto invrho = If(rho>0.0, 1.0/rho, SIMD<double,FMM_SW>(0.0));
              auto [si,co] = sincos(rhokappa);
              
              auto kernelc = (1/(4*M_PI))*invrho*SIMD<Complex,FMM_SW>(co,si);
              vsum += kernelc * c;   
              
              auto kernel = 
                invrho*invrho * InnerProduct(p-x, d) * 
                kernelc * SIMD<Complex,FMM_SW>(-1.0, rhokappa);
              
              vsum += kernel * c2;
            }
          sum += HSum(vsum);
        }
      else
        {
          // static Timer t("mptool singmp, evaluate, chargedipoles"); RegionTimer r(t);
          // t.AddFlops (chargedipoles.Size());
          
          for (auto [x,c,d,c2] : chargedipoles)
            if (double rho = L2Norm(p-x); rho > 0)
              {
                sum += (1/(4*M_PI))*exp(Complex(0,rho*mp.Kappa())) / rho * c;
                
                Vec<3> drhodp = 1.0/rho * (p-x);
                Complex dGdrho = (1/(4*M_PI))*exp(Complex(0,rho*mp.Kappa())) *
                  (Complex(0, mp.Kappa())/rho - 1.0/sqr(rho));
                
                sum += dGdrho * InnerProduct(drhodp, d) * c2;
              }
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
          {
            static int cnt = 0;
            cnt++;
            if (cnt < 3)
              cout << "we know what we do - evaluateDeriv not implemented for dipoles in SingularMLExpansion" << endl;
            // return sum;
            // throw Exception("EvaluateDeriv not implemented for dipoles in SingularMLExpansion");
          }
        if (chargedipoles.Size())
            throw Exception("EvaluateDeriv not implemented for dipoles in SingularMLExpansion");

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
        total_sources = charges.Size() + dipoles.Size() + chargedipoles.Size();
        for (auto & child : childs)
          if (child)
            {
              child->CalcTotalSources();
              total_sources += child->total_sources;
            }
      }
      
      void CalcMP(Array<RecordingSS> * recording, Array<Node*> * nodes_to_process)
      {
        // mp.SH().Coefs() = 0.0;
        if (childs[0])
          {
            if (total_sources < 1000 || recording)
              for (auto & child : childs)
                child->CalcMP(recording, nodes_to_process);
            else
              ParallelFor (8, [&] (int nr)
                           {
                             childs[nr] -> CalcMP(recording, nodes_to_process);
                           });

            
            for (auto & child : childs){
              if (recording && child->mp.SH().Coefs().Size() > 0)
                *recording += RecordingSS(&child->mp, &mp, center-child->center);
              else
                child->mp.TransformAdd(mp, center-child->center);
            }
          }
        else
          {
            if (charges.Size()+dipoles.Size()+chargedipoles.Size()+currents.Size() == 0)
              {
                mp = SphericalExpansion<Singular,entry_type> (-1, mp.Kappa(), 1.);
                return;
              }

            // make simd charges, comment this block for testing ...
            simd_charges.SetSize( (charges.Size()+FMM_SW-1)/FMM_SW);
            size_t i = 0, ii = 0;
            for ( ; i+FMM_SW <= charges.Size(); i+=FMM_SW, ii++)
              {
                std::array<tuple<Vec<3>,entry_type>, FMM_SW> ca;
                for (int j = 0; j < FMM_SW; j++) ca[j] = charges[i+j];
                simd_charges[ii] = MakeSimd(ca);
              }
            if (i < charges.Size())
              {
                std::array<tuple<Vec<3>,entry_type>, FMM_SW> ca;
                int j = 0;
                for ( ; i+j < charges.Size(); j++) ca[j] = charges[i+j];
                for ( ; j < FMM_SW; j++) ca[j] = tuple( get<0>(ca[0]), entry_type{0.0} );
                simd_charges[ii] = MakeSimd(ca);                
              }

            simd_dipoles.SetSize( (dipoles.Size()+FMM_SW-1)/FMM_SW);
            i = 0, ii = 0;
            for ( ; i+FMM_SW <= dipoles.Size(); i+=FMM_SW, ii++)
              {
                std::array<tuple<Vec<3>,Vec<3>,entry_type>, FMM_SW> di;
                for (int j = 0; j < FMM_SW; j++) di[j] = dipoles[i+j];
                simd_dipoles[ii] = MakeSimd(di);
              }
            if (i < dipoles.Size())
              {
                std::array<tuple<Vec<3>,Vec<3>,entry_type>, FMM_SW> di;
                int j = 0;
                for ( ; i+j < dipoles.Size(); j++) di[j] = dipoles[i+j];
                for ( ; j < FMM_SW; j++) di[j] = tuple( get<0>(di[0]), get<1>(di[0]), entry_type{0.0} );
                simd_dipoles[ii] = MakeSimd(di);
              }


            simd_chargedipoles.SetSize( (chargedipoles.Size()+FMM_SW-1)/FMM_SW);
            i = 0, ii = 0;
            for ( ; i+FMM_SW <= chargedipoles.Size(); i+=FMM_SW, ii++)
              {
                std::array<tuple<Vec<3>,entry_type,Vec<3>,entry_type>, FMM_SW> di;
                for (int j = 0; j < FMM_SW; j++) di[j] = chargedipoles[i+j];
                simd_chargedipoles[ii] = MakeSimd(di);
              }
            if (i < chargedipoles.Size())
              {
                std::array<tuple<Vec<3>,entry_type,Vec<3>,entry_type>, FMM_SW> di;
                int j = 0;
                for ( ; i+j < chargedipoles.Size(); j++) di[j] = chargedipoles[i+j];
                for ( ; j < FMM_SW; j++) di[j] = tuple( get<0>(di[0]), entry_type{0.0}, get<2>(di[0]), entry_type{0.0} );
                simd_chargedipoles[ii] = MakeSimd(di);
              }

            
            if (nodes_to_process)
                *nodes_to_process += this;
            else {
              for (auto [x,c] : charges)
                mp.AddCharge (x-center,c);

              for (auto [x,d,c] : dipoles)
                mp.AddDipole (x-center, d, c);

              for (auto [x,c,d,c2] : chargedipoles)
                mp.AddChargeDipole (x-center, c, d, c2);
              
              for (auto [sp,ep,j,num] : currents)
                mp.AddCurrent (sp-center, ep-center, j, num);
            }
          }
      }
      
      entry_type EvaluateMP(Vec<3> p) const
      {
        if (charges.Size() || dipoles.Size() || chargedipoles.Size())
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
        if (charges.Size() || dipoles.Size() || chargedipoles.Size() || !childs[0])
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
        for (auto [x,c,d,c2] : chargedipoles)
          ost << "xi = " << x << ", c = " << c << ", di = " << d << ", ci = " << c2 << endl;

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

      void TraverseTree (const std::function<void(Node&)> & func)
      {
        func(*this);
        for (auto & child : childs)
          if (child)
            child->TraverseTree(func);
      }
    };
    
    FMM_Parameters fmm_params;
    Node root;    
    bool havemp = false;
    
  public:
    SingularMLExpansion (Vec<3> center, double r, double kappa, FMM_Parameters _params = FMM_Parameters())
      : fmm_params(_params), root(center, r, 0, kappa, fmm_params)
    {
      nodes_on_level = 0;
      nodes_on_level[0] = 1;
    }

    double Kappa() const { return root.mp.Kappa(); }
    
    void AddCharge(Vec<3> x, entry_type c)
    {
      root.AddCharge(x, c);
    }

    void AddDipole(Vec<3> x, Vec<3> d, entry_type c)
    {
      root.AddDipole(x, d, c);
    }

    void AddChargeDipole(Vec<3> x, entry_type c, Vec<3> dir, entry_type c2)
    {
      root.AddChargeDipole(x, c, dir, c2);
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
      static Timer ts2mp("mptool compute singular MLMP - source2mp");
      static Timer tS2S("mptool compute singular MLMP - S->S");
      static Timer trec("mptool comput singular recording");
      static Timer tsort("mptool comput singular sort");      

      /*
      int maxlevel = 0;
      for (auto [i,num] : Enumerate(nodes_on_level))
        if (num > 0) maxlevel = i;

      for (int i = 0; i <= maxlevel; i++)
        cout << "sing " <<  i << ": " << nodes_on_level[i] << endl;
      */
      
      root.CalcTotalSources();

      if constexpr (false)
        // direct evaluation of S->S
        root.CalcMP(nullptr, nullptr);
      else
        {
          
          Array<RecordingSS> recording;
          Array<Node*> nodes_to_process;

          {
            RegionTimer reg(trec);
            root.CalcMP(&recording, &nodes_to_process);
          }
      
          {
            RegionTimer rs2mp(ts2mp);
            ParallelFor(nodes_to_process.Size(), [&](int i)
            {
              auto node = nodes_to_process[i];
              for (auto [x,c]: node->charges)
                node->mp.AddCharge(x-node->center, c);
              for (auto [x,d,c]: node->dipoles)
                node->mp.AddDipole(x-node->center, d, c);
              for (auto [x,c,d,c2]: node->chargedipoles)
                node->mp.AddChargeDipole(x-node->center, c, d, c2);
              for (auto [sp,ep,j,num]: node->currents)
                node->mp.AddCurrent(sp-node->center, ep-node->center, j, num);
            }, TasksPerThread(4));
          }
          
          {
            RegionTimer reg(tsort);
            QuickSort (recording, [] (auto & a, auto & b)
            {
              if (a.len < (1-1e-8) * b.len) return true;
              if (a.len > (1+1e-8) * b.len) return false;
              return a.theta < b.theta;
            });
          }
      
          double current_len = -1e100;
          double current_theta = -1e100;
          Array<RecordingSS*> current_batch;
          Array<Array<RecordingSS*>> batch_group;
          Array<double> group_lengths;
          Array<double> group_thetas;
          for (auto & record : recording)
            {
              bool len_changed = fabs(record.len - current_len) > 1e-8;
              bool theta_changed = fabs(record.theta - current_theta) > 1e-8;
              if ((len_changed || theta_changed) && current_batch.Size() > 0) {
                batch_group.Append(current_batch);
                group_lengths.Append(current_len);
                group_thetas.Append(current_theta);
                current_batch.SetSize(0);
              }
              
              current_len = record.len;
              current_theta = record.theta;
              current_batch.Append(&record);
            }
          
          if (current_batch.Size() > 0) {
            batch_group.Append(current_batch);
            group_lengths.Append(current_len);
            group_thetas.Append(current_theta);
          }

          {
            RegionTimer rS2S(tS2S);
            // ParallelFor(batch_group.Size(), [&](int i) {
            for (int i = 0; i < batch_group.Size(); i++){
              // *testout << "Processing batch " << i << " of size " << batch_group[i].Size() << ", with len = " << group_lengths[i] << ", theta = " << group_thetas[i] << endl;
              int chunk_size = 24;
              if (batch_group[i].Size() < chunk_size)
                ProcessBatchSS(batch_group[i], group_lengths[i], group_thetas[i]);
              else
                ParallelForRange(IntRange(batch_group[i].Size()), [&](IntRange range) {
                  auto sub_batch = batch_group[i].Range(range.First(), range.Next());
                  ProcessBatchSS(sub_batch, group_lengths[i], group_thetas[i]);
                }, TasksPerThread(4));
            }
          }
        }

      // cout << "have singular:" << endl;
      // PrintStatistics (cout);
      havemp = true;
    }

    entry_type Evaluate (Vec<3> p) const
    {
      if (havemp)
        return root.EvaluateMP(p);
      else
        return root.Evaluate(p);
    }


    void PrintStatistics (ostream & ost)
    {
      int levels = 0;
      int cnt = 0;
      root.TraverseTree( [&](Node & node) {
        levels = max(levels, node.level);
        cnt++;
      });
      ost << "levels: " << levels << endl;
      ost << "nodes: " << cnt << endl;

      Array<int> num_on_level(levels+1);
      Array<int> order_on_level(levels+1);
      Array<size_t> coefs_on_level(levels+1);            
      num_on_level = 0;
      order_on_level = 0;
      root.TraverseTree( [&](Node & node) {
        num_on_level[node.level]++;
        order_on_level[node.level] = max(order_on_level[node.level],node.mp.Order());
        coefs_on_level[node.level] += node.mp.SH().Coefs().Size();
      });

      cout << "num on level" << endl;
      for (int i = 0; i < num_on_level.Size(); i++)
        cout << i << ": " << num_on_level[i] << ", order = " << order_on_level[i] << ", coefs " << coefs_on_level[i] << endl;

      size_t totcoefs = 0;
      for (auto n : coefs_on_level)
        totcoefs += n;
      cout << "total mem in coefs: " << sizeof(entry_type)*totcoefs / sqr(1024) << " MB" << endl;
    }


    
    template <typename entry_type2>
    friend class RegularMLExpansion;
  };


  template <typename entry_type>
  inline ostream & operator<< (ostream & ost, const SingularMLExpansion<entry_type> & mlmp)
  {
    mlmp.Print(ost);
    return ost;
  }


  // *********************************** Regular multilevel Expansion
  

  template <typename elem_type=Complex>
  class NGS_DLL_HEADER RegularMLExpansion
  {
    static Array<size_t> nodes_on_level;

    
    struct RecordingRS
    {
      const SphericalExpansion<Singular,elem_type> * mpS;
      SphericalExpansion<Regular,elem_type> * mpR;
      Vec<3> dist;
      double len, theta, phi;
    public:
      RecordingRS() = default;
      RecordingRS (const SphericalExpansion<Singular,elem_type> * ampS,
                   SphericalExpansion<Regular,elem_type> * ampR,
                   Vec<3> adist)
        : mpS(ampS), mpR(ampR), dist(adist)
      {
        std::tie(len, theta, phi) = SphericalCoordinates(dist);
      }
    };

    static void ProcessBatchRS(FlatArray<RecordingRS*> batch, double len, double theta) {
      // static Timer t("ProcessBatchRS"); RegionTimer reg(t, batch.Size());
      constexpr int vec_length = VecLength<elem_type>;
      int batch_size = batch.Size();
      int N = batch_size * vec_length;
      // *testout << "Processing batch of size " << batch.Size() << ", with N = " << N << ", vec_length = " << vec_length << ", Type: " << typeid(elem_type).name() << ", len = " << len << ", theta = " << theta << endl;

      if (N <= 1 || batch_size <= 1) {
        for (auto* rec : batch) {
          rec->mpS->TransformAdd(*rec->mpR, rec->dist);
        }
      }
      else if (N <= 3) {
        ProcessVectorizedBatchRS<3, vec_length>(batch, len, theta);
      }
      else if (N <= 4) {
        ProcessVectorizedBatchRS<4, vec_length>(batch, len, theta);
      }
      else if (N <= 6) {
        ProcessVectorizedBatchRS<6, vec_length>(batch, len, theta);
      }
      else if (N <= 12) {
        ProcessVectorizedBatchRS<12, vec_length>(batch, len, theta);
      }
      else if (N <= 24) {
        ProcessVectorizedBatchRS<24, vec_length>(batch, len, theta);
      }
      else if (N <= 48) {
        ProcessVectorizedBatchRS<48, vec_length>(batch, len, theta);
      }
      else if (N <= 96) {
        ProcessVectorizedBatchRS<96, vec_length>(batch, len, theta);
      }
      else if (N <= 192) {
        ProcessVectorizedBatchRS<192, vec_length>(batch, len, theta);
      }
      else {
        // Split large batches
        /*
        ProcessBatch(batch.Range(0, 192 / vec_length), len, theta);
        ProcessBatch(batch.Range(192 / vec_length, batch_size), len, theta);
        */

        /*
        ParallelFor (2, [&] (int i)
        {
          if (i == 0)
            ProcessBatchRS(batch.Range(0, 192 / vec_length), len, theta);
          else
            ProcessBatchRS(batch.Range(192 / vec_length, batch_size), len, theta);            
        }, 2);
        */

        
        size_t chunksize = 192/vec_length;
        size_t num = (batch.Size()+chunksize-1) / chunksize;
        ParallelFor (num, [&](int i)
        {
          ProcessBatchRS(batch.Range(i*chunksize, min((i+1)*chunksize, batch.Size())), len, theta);          
        }, num);

      }
    }


    template<int N, int vec_length>
    static void ProcessVectorizedBatchRS(FlatArray<RecordingRS*> batch, double len, double theta) {

      // static Timer t("ProcessVectorizedBatch, N = "+ToString(N) + ", vec_len = " + ToString(vec_length));
      // RegionTimer reg(t, batch[0]->mpS->SH().Order());
      // static Timer ttobatch("mptools - copy to batch 2");
      // static Timer tfrombatch("mptools - copy from batch 2");      
      
      // *testout << "Processing vectorized batch of size " << batch.Size() << ", with N = " << N << ", vec_length = " << vec_length << ", len = " << len << ", theta = " << theta << endl;
      SphericalExpansion<Singular, Vec<N,Complex>> vec_source(batch[0]->mpS->Order(), batch[0]->mpS->Kappa(), batch[0]->mpS->RTyp());
      // SphericalExpansion<Singular, elem_type> tmp_source{*batch[0]->mpS};
      SphericalExpansion<Regular, elem_type> tmp_target{*batch[0]->mpR};
      SphericalExpansion<Regular, Vec<N,Complex>> vec_target(batch[0]->mpR->Order(), batch[0]->mpR->Kappa(), batch[0]->mpR->RTyp());

      // Copy multipoles into vectorized multipole
      // ttobatch.Start();
      for (int i = 0; i < batch.Size(); i++)
      {
        auto source_i = VecVector2Matrix (batch[i]->mpS->SH().Coefs());
        auto source_mati = VecVector2Matrix (vec_source.SH().Coefs()).Cols(i*vec_length, (i+1)*vec_length);
        batch[i]->mpS->SH().RotateZ(batch[i]->phi,
            [source_i, source_mati] (size_t ii, Complex factor)
            {
                source_mati.Row(ii) = factor * source_i.Row(ii);
            });
      }

      // ttobatch.Stop();

      vec_source.SH().RotateY(theta);
      vec_source.ShiftZ(-len, vec_target);
      vec_target.SH().RotateY(-theta); 

      // Copy vectorized multipole into individual multipoles
      // tfrombatch.Start();
      for (int i = 0; i < batch.Size(); i++) {
        // auto source_i = VecVector2Matrix (tmp_target.SH().Coefs());
        auto source_mati = VecVector2Matrix (vec_target.SH().Coefs()).Cols(i*vec_length, (i+1)*vec_length);
        auto targeti = VecVector2Matrix(batch[i]->mpR->SH().Coefs());
        
        tmp_target.SH().RotateZ(-batch[i]->phi,
                                [source_mati, targeti] (size_t ii, Complex factor)
                                          {
                                            // source_i.Row(ii) = factor * source_mati.Row(ii);
                                            AtomicAdd (VectorView(targeti.Row(ii)), factor * source_mati.Row(ii));
                                          });
        // for (int j = 0; j < tmp_target.SH().Coefs().Size(); j++)
        // AtomicAdd(batch[i]->mpR->SH().Coefs()[j], tmp_target.SH().Coefs()[j]);
      }
      // tfrombatch.Stop();

    }

    
    struct Node
    {
      Vec<3> center;
      double r;
      int level;
      std::array<unique_ptr<Node>,8> childs;
      SphericalExpansion<Regular,elem_type> mp;
      Array<Vec<3>> targets;
      Array<tuple<Vec<3>,double>> vol_targets;      
      int total_targets;
      std::mutex node_mutex;
      atomic<bool> have_childs{false};

      Array<const typename SingularMLExpansion<elem_type>::Node*> singnodes;
      const FMM_Parameters & params;

      
      Node (Vec<3> acenter, double ar, int alevel, double kappa, const FMM_Parameters & _params)
        : center(acenter), r(ar), level(alevel),
          // mp(MPOrder(ar*kappa), kappa, ar) // 1.0/min(1.0, 0.25*r*kappa))
          mp(-1, kappa, ar), params(_params)
          // : center(acenter), r(ar), level(alevel), mp(MPOrder(ar*kappa), kappa, 1.0)
      {
        if (level < nodes_on_level.Size())
          nodes_on_level[level]++;
      }

      void Allocate()
      {
        // mp = SphericalExpansion<Regular,elem_type>(MPOrder(r*mp.Kappa()), mp.Kappa(), r);
        mp = SphericalExpansion<Regular,elem_type>(params.minorder+2*r*mp.Kappa(), mp.Kappa(), r);        
      }
      
      
      void CreateChilds(bool allocate = false)
      {
        if (childs[0]) throw Exception("have already childs");
        // create children nodes:
        for (int i = 0; i < 8; i++)
          {
            Vec<3> cc = center;
            cc(0) += (i&1) ? r/2 : -r/2;
            cc(1) += (i&2) ? r/2 : -r/2;
            cc(2) += (i&4) ? r/2 : -r/2;
            childs[i] = make_unique<Node> (cc, r/2, level+1, mp.Kappa(), params);
            if (allocate)
              childs[i] -> Allocate();
          }
        have_childs = true;
      }
      
      void AddSingularNode (const typename SingularMLExpansion<elem_type>::Node & singnode, bool allow_refine,
                            Array<RecordingRS> * recording)
      {
        if (mp.SH().Order() < 0) return;
        if (singnode.mp.SH().Order() < 0) return;
        // if (L2Norm(singnode.mp.SH().Coefs()) == 0) return;
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
                  AddSingularNode (*child, allow_refine, recording);
                return;
              }

            // static Timer t("mptool transform Helmholtz-criterion"); RegionTimer r(t);
            if (recording)
              *recording += RecordingRS(&singnode.mp, &mp, dist);
            else
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
                  CreateChilds(true);
                
                for (auto & ch : childs)
                  ch -> AddSingularNode (singnode, allow_refine, recording);
              }
            else
              {
                if (total_targets < 1000 || recording)
                  {
                    for (auto & ch : childs)
                      if (ch)
                        ch -> AddSingularNode (singnode, allow_refine, recording);
                  }
                else
                  ParallelFor (8, [&] (int nr)
                               {
                                 if (childs[nr])
                                   childs[nr] -> AddSingularNode (singnode, allow_refine, recording);
                               });
                
                if (targets.Size()+vol_targets.Size())
                  singnodes.Append(&singnode);
              }
          }
        else
          {
            for (auto & childsing : singnode.childs)
              AddSingularNode (*childsing, allow_refine, recording);
          }
      }

      void LocalizeExpansion(bool allow_refine)
      {
        if (allow_refine)
          if (mp.Order() > 30 && !childs[0])
            CreateChilds(allow_refine);

        if (childs[0])
          {
            if (total_targets < 1000)
              {
                for (int nr = 0; nr < 8; nr++)
                  {
                    if (L2Norm(mp.SH().Coefs()) > 0)
                      mp.TransformAdd (childs[nr]->mp, childs[nr]->center-center);
                    childs[nr]->LocalizeExpansion(allow_refine);
                  }
              }
            else
              ParallelFor(8, [&] (int nr)
              {
                if (L2Norm(mp.SH().Coefs()) > 0)
                  mp.TransformAdd (childs[nr]->mp, childs[nr]->center-center);
                childs[nr]->LocalizeExpansion(allow_refine);
              });
            mp = SphericalExpansion<Regular,elem_type>(-1, mp.Kappa(), 1.);
            //mp.SH().Coefs()=0.0;
          }
      }
      
      elem_type Evaluate (Vec<3> p) const
      {
        elem_type sum{0.0};

        int childnum = 0;
        if (p(0) > center(0)) childnum += 1;
        if (p(1) > center(1)) childnum += 2;
        if (p(2) > center(2)) childnum += 4;
        if (childs[childnum])
          sum = childs[childnum]->Evaluate(p);
        else
          {
            // static Timer t("mptool regmp, evaluate reg"); RegionTimer r(t);          
            sum = mp.Eval(p-center);
          }

        {
          // static Timer t("mptool regmp, evaluate, singnode"); RegionTimer r(t);
          for (auto sn : singnodes)
            sum += sn->EvaluateMP(p);
        }
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

      void TraverseTree (const std::function<void(Node&)> & func)
      {
        func(*this);
        for (auto & child : childs)
          if (child)
            child->TraverseTree(func);
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

      int GetChildNum (Vec<3> x) const
      {
        int childnum  = 0;
        if (x(0) > center(0)) childnum += 1;
        if (x(1) > center(1)) childnum += 2;
        if (x(2) > center(2)) childnum += 4;
        return childnum;
      }

      void AddTarget (Vec<3> x)
      {
        // if (childs[0])
        if (have_childs) // quick check without locking
          {
            // directly send to childs:
            int childnum  = GetChildNum(x);
            childs[childnum] -> AddTarget( x );
            return;
          }

        lock_guard<mutex> guard(node_mutex);

        if (have_childs) // test again after locking
        {
          // directly send to childs:
          int childnum  = GetChildNum(x);
          childs[childnum] -> AddTarget(x);
          return;
        }

        targets.Append( x );

        // if (r*mp.Kappa() < 1e-8) return;
        if (level > 20) return;        
        if (targets.Size() < params.maxdirect && r*mp.Kappa() < 5)
          return;

        CreateChilds();

        for (auto t : targets)
          AddTarget (t);
        for (auto [x,r] : vol_targets)
          AddVolumeTarget (x,r);

        targets.SetSize0();
        vol_targets.SetSize0();
      }


      void AddVolumeTarget (Vec<3> x, double tr)
      {
        if (MaxNorm(x-center) > r+tr) return;

        if (have_childs)
          {
            for (auto & child : childs)
              child->AddVolumeTarget(x, tr);
            return;
          }

        
        lock_guard<mutex> guard(node_mutex);
        
        if (have_childs)
          {
            for (auto & child : childs)
              child->AddVolumeTarget(x, tr);
            return;
          }

        
        vol_targets.Append (tuple(x,tr));

        if (level > 20) return;
        if (vol_targets.Size() < params.maxdirect && (r*mp.Kappa() < 5))
          return;

        CreateChilds();

        for (auto t : targets)
          AddTarget (t);
        for (auto [x,r] : vol_targets)
          AddVolumeTarget (x,r);

        targets.SetSize0();
        vol_targets.SetSize0();
      }

      

      void CalcTotalTargets()
      {
        total_targets = targets.Size() + vol_targets.Size();
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
          mp = SphericalExpansion<Regular,elem_type>(-1, mp.Kappa(),1.);
      }

      void AllocateMemory()
      {
        for (auto & child : childs)
          if (child)
            child->AllocateMemory();

        if (total_targets > 0)
          Allocate();
        // mp = SphericalExpansion<Regular,elem_type>(MPOrder(r*mp.Kappa()), mp.Kappa(), r); // -1, mp.Kappa(),1.);
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

    FMM_Parameters fmm_params;
    Node root;
    shared_ptr<SingularMLExpansion<elem_type>> singmp;
    
  public:
  RegularMLExpansion (shared_ptr<SingularMLExpansion<elem_type>> asingmp, Vec<3> center, double r,
                      const FMM_Parameters & _params)
  : fmm_params(_params), root(center, r, 0, asingmp->Kappa(), fmm_params), singmp(asingmp)
  {
      if (!singmp->havemp) throw Exception("first call Calc for singular MP");
      root.Allocate();
      
      nodes_on_level = 0;
      nodes_on_level[0] = 1;
      {
        static Timer t("mptool compute regular MLMP"); RegionTimer rg(t);
        root.AddSingularNode(singmp->root, true, nullptr);
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

  RegularMLExpansion (Vec<3> center, double r, double kappa, const FMM_Parameters & _params)
  : fmm_params(_params), root(center, r, 0, kappa, fmm_params)
  {
    nodes_on_level = 0;
    nodes_on_level[0] = 1;
  }
  
    void AddTarget (Vec<3> t)
    {
      root.AddTarget (t);
    }

    void AddVolumeTarget (Vec<3> t, double r)
    {
      root.AddVolumeTarget (t, r);
    }

    void CalcMP(shared_ptr<SingularMLExpansion<elem_type>> asingmp, bool onlytargets = true)
    {
      static Timer t("mptool regular MLMP"); RegionTimer rg(t);
      static Timer tremove("removeempty");
      static Timer trec("mptool regular MLMP - recording");
      static Timer tsort("mptool regular MLMP - sort");       
      
      singmp = asingmp;

      
      root.CalcTotalTargets();
      // cout << "before remove empty trees:" << endl;
      // PrintStatistics(cout);

      /*
      tremove.Start();
      if (onlytargets)
        root.RemoveEmptyTrees();
      tremove.Stop();
      */
      
      root.AllocateMemory();

      // cout << "after allocating regular:" << endl;
      // PrintStatistics(cout);

      // cout << "starting S-R converion" << endl;
      // PrintStatistics(cout);


      if constexpr (false)
        {
          root.AddSingularNode(singmp->root, !onlytargets, nullptr);
        }
      else
        {  // use recording
          Array<RecordingRS> recording;
          {
            RegionTimer rrec(trec);
            root.AddSingularNode(singmp->root, !onlytargets, &recording);
          }
          
          // cout << "recorded: " << recording.Size() << endl;
          {
            RegionTimer reg(tsort);
            QuickSort (recording, [] (auto & a, auto & b)
            {
              if (a.len < (1-1e-8) * b.len) return true;
              if (a.len > (1+1e-8) * b.len) return false;
              return a.theta < b.theta;
            });
          }
          
          double current_len = -1e100;
          double current_theta = -1e100;
          Array<RecordingRS*> current_batch;
          Array<Array<RecordingRS*>> batch_group;
          Array<double> group_lengths;
          Array<double> group_thetas;
          for (auto & record : recording)
            {
              bool len_changed = fabs(record.len - current_len) > 1e-8;
              bool theta_changed = fabs(record.theta - current_theta) > 1e-8;
              if ((len_changed || theta_changed) && current_batch.Size() > 0) {
                // ProcessBatch(current_batch, current_len, current_theta);
                batch_group.Append(current_batch);
                group_lengths.Append(current_len);
                group_thetas.Append(current_theta);
                current_batch.SetSize(0);
              }
              
              current_len = record.len;
              current_theta = record.theta;
              current_batch.Append(&record);
            }
          if (current_batch.Size() > 0) {
            // ProcessBatch(current_batch, current_len, current_theta);
            batch_group.Append(current_batch);
            group_lengths.Append(current_len);
            group_thetas.Append(current_theta);
          }
          
          ParallelFor(batch_group.Size(), [&](int i) {
            ProcessBatchRS(batch_group[i], group_lengths[i], group_thetas[i]);
          }, TasksPerThread(4));
        }
          
      
      /*
      int maxlevel = 0;
      for (auto [i,num] : Enumerate(RegularMLExpansion::nodes_on_level))
        if (num > 0) maxlevel = i;

      for (int i = 0; i <= maxlevel; i++)
        cout << "reg " << i << ": " << RegularMLExpansion::nodes_on_level[i] << endl;
      */

      // cout << "starting R-R converion" << endl;
      // PrintStatistics(cout);
      
      static Timer tloc("mptool regular localize expansion"); RegionTimer rloc(tloc);
      root.LocalizeExpansion(!onlytargets);


      // cout << "R-R conversion done" << endl;
      // PrintStatistics(cout);      
    }

    void PrintStatistics (ostream & ost)
    {
      int levels = 0;
      int cnt = 0;
      root.TraverseTree( [&](Node & node) {
        levels = max(levels, node.level);
        cnt++;
      });
      ost << "levels: " << levels << endl;
      ost << "nodes: " << cnt << endl;

      Array<int> num_on_level(levels+1);
      Array<int> order_on_level(levels+1);
      Array<size_t> coefs_on_level(levels+1);            
      num_on_level = 0;
      order_on_level = 0;
      root.TraverseTree( [&](Node & node) {
        num_on_level[node.level]++;
        order_on_level[node.level] = max(order_on_level[node.level],node.mp.Order());
        coefs_on_level[node.level] += node.mp.SH().Coefs().Size();
      });

      cout << "num on level" << endl;
      for (int i = 0; i < num_on_level.Size(); i++)
        cout << i << ": " << num_on_level[i] << ", order = " << order_on_level[i] << ", coefs " << coefs_on_level[i] << endl;

      size_t totcoefs = 0;
      for (auto n : coefs_on_level)
        totcoefs += n;
      cout << "total mem in coefs: " << sizeof(elem_type)*totcoefs / sqr(1024) << " MB" << endl;
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
      // if (L2Norm(p-root.center) > root.r) return elem_type{0.0};
      
      if (MaxNorm(p-root.center) > root.r)
        return singmp->Evaluate(p);
      return root.Evaluate(p);
    }

    elem_type EvaluateDirectionalDerivative (Vec<3> p, Vec<3> d) const
    {
        if (L2Norm(p-root.center) > root.r) return elem_type{0.0};
        return root.EvaluateDirectionalDerivative(p, d);
    }

  };


  template <typename elem_type>
  inline ostream & operator<< (ostream & ost, const RegularMLExpansion<elem_type> & mlmp)
  {
    mlmp.Print(ost);
    // ost << "RegularMLExpansion" << endl;
    return ost;
  }





  
}
#endif
