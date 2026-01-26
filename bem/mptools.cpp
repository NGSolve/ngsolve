#include "mptools.hpp"
#include "mp_coefficient.hpp"
#include "meshaccess.hpp"



namespace ngsbem
{


  PrecomputedSqrts presqrt;

  PrecomputedSqrts :: PrecomputedSqrts()
  {
    size_t N = 10001;
    sqrt_int.SetSize(N);
    // inv_sqrt_int.SetSize(N);
    sqrt_n_np1.SetSize(N);  // // sqrt(n*(n+1))  
    inv_sqrt_2np1_2np3.SetSize(N);  // 1/sqrt( (2n+1)*(2n+3) )
    
    for (size_t n = 0; n < sqrt_int.Size(); n++)
      {
        sqrt_int[n] = sqrt(n);
        // inv_sqrt_int[n] = 1.0/sqrt(n);
        sqrt_n_np1[n] = sqrt(n*(n+1));
        inv_sqrt_2np1_2np3[n] = 1.0/sqrt( (2*n+1)*(2*n+3) );
      }
  }


  
  template <typename entry_type>
  entry_type SphericalHarmonics<entry_type> :: Eval (double theta, double phi) const
  {
    // static Timer t("mptool sh evaluate"); RegionTimer rg(t);

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
    // static Timer ts("mptool sh evalorders small");
    // static Timer tl("mptool sh evalorders large");
    // RegionTimer rg(order < 30 ? ts : tl);
    // if (order > 30) tl.Start();

    ArrayMem<double,1600> mem(sqr(order+1));
    FlatMatrix<double> legfunc(order+1, order+1, mem.Data());
    NormalizedLegendreFunctions (order, order, cos(theta), legfunc);

    VectorMem<40,Complex> exp_imphi(order+1);
    Complex exp_iphi(cos(phi), sin(phi));
    Complex prod = 1.0 / sqrt(4*M_PI);
    for (int i = 0; i <= order; i++)
      {
        exp_imphi(i) = prod;
        prod *= -exp_iphi;
      }

    for (int n = 0; n <= order; n++)
      {
        auto coefsn = CoefsN(n);
        
        entry_type sum1{exp_imphi(0)*legfunc(0,n)*coefsn(n)}, sum2{0.0};
        for (int m = 1; m <= n; m++)
          {
            sum1 += exp_imphi(m) * legfunc(m, n) * coefsn(n+m);
            sum2 += conj(exp_imphi(m)) * legfunc(m, n) * coefsn(n-m);
          }
        
        vals(n) = sum1+sum2;        
      }
  }

  
  template <typename entry_type>
  void SphericalHarmonics<entry_type> :: Calc (Vec<3> x, FlatVector<Complex> shapes)
  {
    auto [theta, phi] = Polar(x);

    ArrayMem<double,1600> mem(sqr(order+1));
    FlatMatrix<double> legfunc(order+1, order+1, mem.Data());
    NormalizedLegendreFunctions (order, order, cos(theta), legfunc);
    
    VectorMem<50,Complex> exp_imphi(order+1);
    Complex exp_iphi(cos(phi), sin(phi));
    
    Complex prod = 1.0/sqrt(4*M_PI);
    for (int i = 0; i <= order; i++)
      {
        exp_imphi(i) = prod;
        prod *= -exp_iphi;
      }

    size_t ii = 0;
    for (int n = 0; n <= order; n++)
      {
        size_t base = ii+n;
        shapes(base) = exp_imphi(0)*legfunc(0,n);
        for (ptrdiff_t m = 1; m <= n; m++)
          {
            shapes(base-m) = conj(exp_imphi(m))*legfunc(m,n);
            shapes(base+m) = exp_imphi(m)*legfunc(m,n);
          }
        ii = base+n+1;
      }
  }

  
  

  // Nail A. Gumerov and Ramani Duraiswami book, formula (2.2.12)
  // add directional derivative divided by kappa to res, both multipoles need same scaling
  template <typename entry_type>  
  void SphericalHarmonics<entry_type> ::
  DirectionalDiffAdd (Vec<3> d, SphericalHarmonics<entry_type> & res, double scale) const
  {
    // static Timer t("mptool Directional Diff Add"); RegionTimer rg(t);
    
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
  void SphericalHarmonics<entry_type> :: FlipZ ()
  {
    throw Exception ("FlipZ not correct!!");
    for (int n = 1; n <= order; n+=2)
      CoefsN(n) *= -1;
  }
  


  template <typename entry_type>
  void SphericalHarmonics<entry_type> :: RotateZ (double alpha)
  {
    // static Timer t("mptool sh RotateZ"); RegionTimer rg(t);
    if (order < 0) return;

    VectorMem<40,Complex> exp_imalpha(order+1);
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
  void SphericalHarmonics<entry_type> :: RotateY (double alpha, bool parallel)
  {
    LocalHeap lh(8*6*sqr(order) + 8*15*order + 2*sizeof(entry_type)*(order+3) + 500, nullptr, parallel);
    // static Timer t("mptool sh RotateY"+ToString(sizeof(entry_type)/sizeof(Complex)));
    // RegionTimer rg(t, order);
    
    
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

    // for (int n=1; n <= order; n++)
    auto transformN = [normalized_leg_func,s,c,this] (int n, LocalHeap & lh)
      {
        HeapReset hr(lh);
          
        FlatVector<> Dmn(2*order+1, lh);
        FlatMatrix<double,RowMajor> trafo(n+1, 2*n+1, lh); 
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
      };


    if (!parallel)
      {
        for (int n = 1; n <= order; n++)
          transformN (n, lh);
      }
    else
      /*
      TaskManager::CreateJob ([this, &lh, transformN] (TaskInfo & ti)
      {
        auto slh = lh.Split();
        for (auto n = ti.task_nr+1; n <= this->order; n += ti.ntasks)
          transformN (n, slh);          
      }, TasksPerThread(4));
      */
      ParallelForRange (IntRange(1,order+1), [&lh, transformN] (IntRange r)
      {
        auto slh = lh.Split();
        for (auto n : r)
          transformN (n, slh);
      }, TasksPerThread(4));
  }




  
  template <typename RADIAL, typename entry_type, typename T_Kappa> template <typename TARGET>
  void SphericalExpansion<RADIAL,entry_type,T_Kappa> :: ShiftZ (double z, SphericalExpansion<TARGET,entry_type,T_Kappa> & target)
  {
    int os = sh.Order();
    int ot = target.SH().Order();

    double scale = Scale();
    double inv_scale = 1.0/scale;
    double tscale = target.Scale();
    double inv_tscale = 1.0/tscale;

    target.SH().Coefs()=0.0;

    LocalHeap lh(( 32*( (os+ot+1)*(os+ot+1) + (os+1 + ot+1) ) + 2*sizeof(entry_type)*(os+ot+3)+ 8*3*(os+ot+1) + 500));

    constexpr bool is_rr = std::is_same<RADIAL,Regular>::value && std::is_same<TARGET,Regular>::value;
    using trafo_type =
        std::conditional_t<
            std::is_same<RADIAL, Singular>::value && std::is_same<TARGET, Regular>::value || std::is_same<T_Kappa, Complex>::value,
            Complex,
            double
        >;

    FlatMatrix<trafo_type> trafo(os+ot+1, max(os,ot)+1, lh);
    FlatMatrix<trafo_type> oldtrafo(os+ot+1, max(os,ot)+1, lh);
    FlatVector<entry_type> hv1(os+1, lh), hv2(ot+1, lh);

    trafo = trafo_type(0.0);

    FlatVector<double> scale_inv_bmn(os+ot+1, lh);
    FlatVector<double> scale_inv_amn(os+ot+1, lh);
    FlatVector<double> powscale(os+ot+1, lh);
    
    // initial values
    if constexpr (std::is_same<RADIAL,Singular>::value && std::is_same<TARGET,Regular>::value)
    {
      SphericalHankel1 (os+ot, kappa*abs(z), inv_tscale, trafo.Col(0));
    }
    if constexpr (is_rr)
    {
      std::swap(scale, inv_tscale);
      std::swap(tscale, inv_scale);
      z = -z;
      SphericalBessel (os+ot, kappa*abs(z), tscale, trafo.Col(0));
    }
    if constexpr (std::is_same<RADIAL,Singular>::value && std::is_same<TARGET,Singular>::value)
    {
      SphericalBessel (os+ot, kappa*abs(z), tscale, trafo.Col(0));
    }

    double prod = 1;
    for (int i = 0; i <= os; i++, prod *= -scale*tscale)
      powscale(i) = prod;
    
    // (185) from paper 'fast, exact, stable, Gumerov+Duraiswami
    if (z < 0)
      for (int l = 1; l < trafo.Height(); l+=2)
        trafo(l,0) *= -1;
    
    for (int l = 0; l <= os+ot; l++)
      trafo(l,0) *= sqrt(2*l+1);

    for (int l = 0; l < os+ot; l++)
      scale_inv_amn(l) = scale/sh.CalcAmn(0,l);

    if (os > 0)
    {
      for (int l = 1; l < os+ot; l++)
        trafo(l,1) = -scale_inv_amn(0) * (sh.CalcAmn(0,l)*tscale*trafo(l+1,0)
              -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,0));
      trafo(0,1) = -scale*tscale*trafo(1,0);
    }

    for (int n = 1; n < trafo.Width()-1; n++)
    {
      for (int l = n; l < os+ot-n; l++)
        trafo(l,n+1) = -scale_inv_amn(n) * (sh.CalcAmn(0,l)*tscale*trafo(l+1,n)
              -sh.CalcAmn(0,l-1)*inv_tscale*trafo(l-1,n)
              -sh.CalcAmn(0,n-1)*scale*trafo(l,n-1));
      trafo(0,n+1) = powscale(n+1)*trafo(n+1,0);
    }
    
    // use symmetry of matrix (up to scaling)
    for (int n = 0; n < trafo.Width(); n++)
      for (int l = n+1; l < trafo.Width(); l++)
        trafo(n,l) = powscale(l-n) * trafo(l,n);

    for (int n = 0; n <= os; n++)
      hv1(n) = sh.Coef(n,0);
    if constexpr (is_rr)
      hv2 = Trans(trafo.Rows(os+1).Cols(ot+1)) * hv1;
    else
      hv2 = trafo.Rows(ot+1).Cols(os+1) * hv1;
    for (int n = 0; n <= ot; n++)
      target.SH().Coef(n,0) = hv2(n);
    
    // recursion over m
    for (int m = 1; m <= min(os,ot); m++)
    {
      for (int l = m-1; l < os+ot-m; l++)
        scale_inv_bmn(l) = scale/sh.CalcBmn(-m,l);
        
      trafo.Swap (oldtrafo);
        
      // fill recursive formula (187)
      for (int l = m; l <= trafo.Height()-1-m; l++)
        trafo(l,m) = scale_inv_bmn(m) * (sh.CalcBmn(-m, l)*inv_tscale*oldtrafo(l-1, m-1)
              -sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,m-1));  
      if (m < trafo.Width()-1)
      {
      for (int l = m; l <= trafo.Height()-1-m; l++)
        trafo(l,m+1) = scale_inv_bmn(m+1)* ( // sh.CalcBmn(m-1,m) * scale*trafo(l,m-1)
              - sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,m)
              + sh.CalcBmn(-m,l)  * inv_tscale*oldtrafo(l-1,m) );
      }
      // fill recursive formula (182)
      for (int n = m+1; n < trafo.Width()-1; n++)
      {
        for (int l = n-1; l <= trafo.Height()-1-m; l++)
          trafo(l,n+1) = scale_inv_bmn(n+1)* (sh.CalcBmn(m-1,n) * scale*trafo(l,n-1)
                       - sh.CalcBmn(m-1,l+1)*tscale*oldtrafo(l+1,n)
                       + sh.CalcBmn(-m,l)  * inv_tscale*oldtrafo(l-1,n) );
      }
        
      for (int n = m; n < os; n++)
        for (int l = n+1; l <= os; l++)
          trafo(n,l) = powscale(l-n) * trafo(l,n);

      for (int n = m; n <= os; n++)
        hv1(n) = sh.Coef(n,m);
      if constexpr (is_rr)
        hv2.Range(m,ot+1) = Trans(trafo.Rows(m,os+1).Cols(m,ot+1)) * hv1.Range(m,os+1);
      else
        hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
      for (int n = m; n <= ot; n++)
        target.SH().Coef(n,m) = hv2(n);
        
      for (int n = m; n <= os; n++)
        hv1(n) = sh.Coef(n,-m);
      if constexpr (is_rr)
        hv2.Range(m,ot+1) = Trans(trafo.Rows(m,os+1).Cols(m,ot+1)) * hv1.Range(m,os+1);
      else
        hv2.Range(m,ot+1) = trafo.Rows(m,ot+1).Cols(m,os+1) * hv1.Range(m,os+1);
      for (int n = m; n <= ot; n++)
        target.SH().Coef(n,-m) = hv2(n);
      }
  }


  
  
  
  template void SphericalExpansion<Regular> :: ShiftZ (double z, SphericalExpansion<Regular> & target);
  template void SphericalExpansion<Singular> :: ShiftZ (double z, SphericalExpansion<Regular> & target);  
  template void SphericalExpansion<Singular> :: ShiftZ (double z, SphericalExpansion<Singular> & target);

  template void SphericalExpansion<Regular,Vec<1,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<1,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<1,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<1,Complex>> & target);  
  template void SphericalExpansion<Singular,Vec<1,Complex>> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<1,Complex>> & target);

  template void SphericalExpansion<Regular,Vec<3,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<3,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<3,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<3,Complex>> & target);  
  template void SphericalExpansion<Singular,Vec<3,Complex>> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<3,Complex>> & target);

  template void SphericalExpansion<Regular,Vec<4,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<4,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<4,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<4,Complex>> & target);  
  template void SphericalExpansion<Singular,Vec<4,Complex>> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<4,Complex>> & target);

  template void SphericalExpansion<Regular,Vec<6,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<6,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<6,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<6,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<6,Complex>> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<6,Complex>> & target);

  template void SphericalExpansion<Regular,Vec<12,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<12,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<12,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<12,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<12,Complex>> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<12,Complex>> & target);

  template void SphericalExpansion<Regular,Vec<24,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<24,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<24,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<24,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<24,Complex>> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<24,Complex>> & target);

  template void SphericalExpansion<Regular,Vec<48,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<48,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<48,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<48,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<48,Complex>> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<48,Complex>> & target);

  template void SphericalExpansion<Regular,Vec<96,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<96,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<96,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<96,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<96,Complex>> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<96,Complex>> & target);

  template void SphericalExpansion<Regular,Vec<192,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<192,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<192,Complex>> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<192,Complex>> & target);
  template void SphericalExpansion<Singular,Vec<192,Complex>> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<192,Complex>> & target);


  template void SphericalExpansion<Regular,Complex,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Complex,Complex> & target);
  template void SphericalExpansion<Singular,Complex,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Complex,Complex> & target);
  template void SphericalExpansion<Singular,Complex,Complex> :: ShiftZ (double z, SphericalExpansion<Singular,Complex,Complex> & target);

  template void SphericalExpansion<Regular,Vec<1,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<1,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<1,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<1,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<1,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<1,Complex>,Complex> & target);

  template void SphericalExpansion<Regular,Vec<3,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<3,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<3,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<3,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<3,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<3,Complex>,Complex> & target);

  template void SphericalExpansion<Regular,Vec<4,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<4,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<4,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<4,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<4,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<4,Complex>,Complex> & target);

  template void SphericalExpansion<Regular,Vec<6,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<6,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<6,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<6,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<6,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<6,Complex>,Complex> & target);

  template void SphericalExpansion<Regular,Vec<12,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<12,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<12,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<12,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<12,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<12,Complex>,Complex> & target);

  template void SphericalExpansion<Regular,Vec<24,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<24,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<24,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<24,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<24,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<24,Complex>,Complex> & target);

  template void SphericalExpansion<Regular,Vec<48,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<48,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<48,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<48,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<48,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<48,Complex>,Complex> & target);

  template void SphericalExpansion<Regular,Vec<96,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<96,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<96,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<96,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<96,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<96,Complex>,Complex> & target);

  template void SphericalExpansion<Regular,Vec<192,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<192,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<192,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Regular,Vec<192,Complex>,Complex> & target);
  template void SphericalExpansion<Singular,Vec<192,Complex>,Complex> :: ShiftZ (double z, SphericalExpansion<Singular,Vec<192,Complex>,Complex> & target);





  // https://fortran-lang.discourse.group/t/looking-for-spherical-bessel-and-hankel-functions-of-first-and-second-kind-and-arbitrary-order/2308/2

  // adapted from fmm3d
  template <typename Tz> 
  void T_besseljs3d (int nterms, Tz z, double scale,
                     SliceVector<Tz> fjs, SliceVector<Tz> fjder)
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
  

  void besseljs3d (int nterms, double z, double scale, SliceVector<double> fjs, SliceVector<double> fjder)
  {
    T_besseljs3d (nterms, z, scale, fjs, fjder);
  }
  void besseljs3d (int nterms, Complex z, double scale, SliceVector<Complex> fjs, SliceVector<Complex> fjder)
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

  template <typename RADIAL, typename entry_type, typename T_Kappa>
  entry_type SphericalExpansion<RADIAL,entry_type,T_Kappa> :: Eval (Vec<3> x) const
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

  template <typename RADIAL, typename entry_type, typename T_Kappa>
  entry_type SphericalExpansion<RADIAL,entry_type,T_Kappa> :: EvalDirectionalDerivative (Vec<3> x, Vec<3> d) const
  {
    if (sh.Order() < 0) return entry_type{0.0};
    SphericalExpansion<RADIAL, entry_type,T_Kappa> tmp(Order(), kappa, RTyp());
    // this->SH().DirectionalDiffAdd(kappa*d, tmp.SH(), Scale());
    this->SH().DirectionalDiffAdd(abs(kappa)*d, tmp.SH(), Scale());
    return tmp.Eval(x);
  }


  template <typename RADIAL, typename entry_type, typename T_Kappa>
  void SphericalExpansion<RADIAL,entry_type,T_Kappa> :: AddCharge (Vec<3> x, entry_type c)
  {
    // static Timer t("mptool AddCharge"); RegionTimer reg(t);
    
    if constexpr (!std::is_same<RADIAL,Singular>())
      throw Exception("AddCharge assumes singular MP");
    
    VectorMem<50,T_Kappa> radial(sh.Order()+1);
    VectorMem<1000,Complex> sh_shapes(sqr (sh.Order()+1));
    
    // SphericalBessel(sh.Order(), kappa*L2Norm(x), Scale(), radial);
    besseljs3d(sh.Order(), kappa*L2Norm(x), Scale(), radial);
    sh.Calc(x, sh_shapes);

    for (int i = 0; i <= sh.Order(); i++)
      {
        // entry_type radc = Complex(0,kappa*radial(i)) * c;
        entry_type radc = kappa*radial(i) * Complex(0,1) * c;
        for (auto j : IntRange(sqr(i), sqr(i+1)))
          sh.Coefs()(j) += Conj(sh_shapes(j)) * radc;
      }
  }

  template <typename RADIAL, typename entry_type, typename T_Kappa>
  void SphericalExpansion<RADIAL, entry_type, T_Kappa> :: AddDipole (Vec<3> x, Vec<3> d, entry_type c)
  {
    // static Timer t("mptool AddDipole"); RegionTimer rg(t);      
    /*
      double eps = 1e-4;
      AddCharge(x+eps*d, -c/(2*eps));
      AddCharge(x-eps*d, c/(2*eps));
      return;
    */
      
    if constexpr (!std::is_same<RADIAL,Singular>())
      throw Exception("AddDipole assumes singular MP");

    SphericalExpansion<Singular, entry_type, T_Kappa> tmp(Order(), kappa, RTyp());
    tmp.AddCharge(x, c);
    // tmp.SH().DirectionalDiffAdd (kappa*d,  this->SH(), Scale());
    tmp.SH().DirectionalDiffAdd (abs(kappa)*d,  this->SH(), Scale());
  }


  template <typename RADIAL, typename entry_type, typename T_Kappa>
  void SphericalExpansion<RADIAL, entry_type, T_Kappa> :: AddPlaneWave (Vec<3> d, entry_type c)
  {
    if constexpr (!std::is_same<RADIAL,Regular>())
      throw Exception("AddPlaneWave assumes regular MP");

    SphericalExpansion<RADIAL, entry_type, T_Kappa> tmp(Order(), kappa, RTyp());
    entry_type fac = kappa / sqrt(M_PI) * c;
    for (int i = 0; i <= Order(); i++)
      {
        tmp.Coef(i,0) += sqrt(2*i+1) * fac;
        fac *= Complex(0,1);
      }

    auto [theta, phi] = SH().Polar(d);
    tmp.SH().RotateY(-theta);
    tmp.SH().RotateZ(-phi);
    this->SH().Coefs() += tmp.SH().Coefs();
  }


  
  
  template <typename RADIAL, typename entry_type, typename T_Kappa>
  void SphericalExpansion<RADIAL, entry_type, T_Kappa> :: AddCurrent (Vec<3> sp, Vec<3> ep, Complex j, int num)
  {
    if constexpr (!std::is_same<RADIAL,Singular>() || !std::is_same<entry_type, Vec<3,Complex>>())
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
  shared_ptr<RegularMLExpansionCF<entry_type>>
  SingularMLExpansionCF<entry_type> :: CreateRegularExpansion(Vec<3> center, double r) const
  {
    mlmp->CalcMP();
    return make_shared<RegularMLExpansionCF<entry_type>> (this->MLExpansion(), center, r);
  }

  template class SphericalHarmonics<Complex>;
  template class SphericalHarmonics<Vec<1,Complex>>;    
  template class SphericalHarmonics<Vec<3,Complex>>;  
  template class SphericalHarmonics<Vec<4,Complex>>;
  template class SphericalHarmonics<Vec<6,Complex>>;
  template class SphericalHarmonics<Vec<12,Complex>>;
  template class SphericalHarmonics<Vec<24,Complex>>;
  template class SphericalHarmonics<Vec<48,Complex>>;
  template class SphericalHarmonics<Vec<96,Complex>>;
  template class SphericalHarmonics<Vec<192,Complex>>;
  
  template class SphericalExpansion<Singular>;
  template class SphericalExpansion<Regular>;
  template class SphericalExpansion<Singular, Vec<1,Complex>>;
  template class SphericalExpansion<Regular, Vec<1,Complex>>;    
  template class SphericalExpansion<Singular, Vec<3,Complex>>;
  template class SphericalExpansion<Regular, Vec<3,Complex>>;    
  template class SphericalExpansion<Singular, Vec<4,Complex>>;
  template class SphericalExpansion<Regular, Vec<4,Complex>>;    
  template class SphericalExpansion<Singular, Vec<6,Complex>>;
  template class SphericalExpansion<Regular, Vec<6,Complex>>;
  template class SphericalExpansion<Singular, Vec<12,Complex>>;
  template class SphericalExpansion<Regular, Vec<12,Complex>>;
  template class SphericalExpansion<Singular, Vec<24,Complex>>;
  template class SphericalExpansion<Regular, Vec<24,Complex>>;
  template class SphericalExpansion<Singular, Vec<48,Complex>>;
  template class SphericalExpansion<Regular, Vec<48,Complex>>;
  template class SphericalExpansion<Singular, Vec<96,Complex>>;
  template class SphericalExpansion<Regular, Vec<96,Complex>>;
  template class SphericalExpansion<Singular, Vec<192,Complex>>;
  template class SphericalExpansion<Regular, Vec<192,Complex>>;

  template class SphericalExpansion<Singular, Complex, Complex>;
  template class SphericalExpansion<Regular, Complex, Complex>;
  template class SphericalExpansion<Singular, Vec<1,Complex>, Complex>;
  template class SphericalExpansion<Regular, Vec<1,Complex>, Complex>;
  template class SphericalExpansion<Singular, Vec<3,Complex>, Complex>;
  template class SphericalExpansion<Regular, Vec<3,Complex>, Complex>;
  template class SphericalExpansion<Singular, Vec<4,Complex>, Complex>;
  template class SphericalExpansion<Regular, Vec<4,Complex>, Complex>;
  template class SphericalExpansion<Singular, Vec<6,Complex>, Complex>;
  template class SphericalExpansion<Regular, Vec<6,Complex>, Complex>;
  template class SphericalExpansion<Singular, Vec<12,Complex>, Complex>;
  template class SphericalExpansion<Regular, Vec<12,Complex>, Complex>;
  template class SphericalExpansion<Singular, Vec<24,Complex>, Complex>;
  template class SphericalExpansion<Regular, Vec<24,Complex>, Complex>;
  template class SphericalExpansion<Singular, Vec<48,Complex>, Complex>;
  template class SphericalExpansion<Regular, Vec<48,Complex>, Complex>;
  template class SphericalExpansion<Singular, Vec<96,Complex>, Complex>;
  template class SphericalExpansion<Regular, Vec<96,Complex>, Complex>;
  template class SphericalExpansion<Singular, Vec<192,Complex>, Complex>;
  template class SphericalExpansion<Regular, Vec<192,Complex>, Complex>;

  
  template class SingularMLExpansionCF<Complex>;
  template class SingularMLExpansionCF<Vec<3,Complex>>;
  
  template<>
  Array<size_t> RegularMLExpansion<Complex>::nodes_on_level(100);
  template<>  
  Array<size_t> SingularMLExpansion<Complex>::nodes_on_level(100);
  template<>
  Array<size_t> RegularMLExpansion<Vec<3,Complex>>::nodes_on_level(100);
  template<>
  Array<size_t> SingularMLExpansion<Vec<3,Complex>>::nodes_on_level(100);
  template<>
  Array<size_t> RegularMLExpansion<Vec<4,Complex>>::nodes_on_level(100);
  template<>
  Array<size_t> SingularMLExpansion<Vec<4,Complex>>::nodes_on_level(100);
  template<>
  Array<size_t> RegularMLExpansion<Vec<6,Complex>>::nodes_on_level(100);
  template<>
  Array<size_t> SingularMLExpansion<Vec<6,Complex>>::nodes_on_level(100);

  template<>
  Array<size_t> RegularMLExpansion<Complex,Complex>::nodes_on_level(100);
  template<>
  Array<size_t> SingularMLExpansion<Complex,Complex>::nodes_on_level(100);
  template<>
  Array<size_t> RegularMLExpansion<Vec<3,Complex>,Complex>::nodes_on_level(100);
  template<>
  Array<size_t> SingularMLExpansion<Vec<3,Complex>,Complex>::nodes_on_level(100);
  template<>
  Array<size_t> RegularMLExpansion<Vec<4,Complex>,Complex>::nodes_on_level(100);
  template<>
  Array<size_t> SingularMLExpansion<Vec<4,Complex>,Complex>::nodes_on_level(100);
  template<>
  Array<size_t> RegularMLExpansion<Vec<6,Complex>,Complex>::nodes_on_level(100);
  template<>
  Array<size_t> SingularMLExpansion<Vec<6,Complex>,Complex>::nodes_on_level(100);


  template class SingularMLExpansion<Complex>;
  template class SingularMLExpansion<Vec<3,Complex>>;

  
}
