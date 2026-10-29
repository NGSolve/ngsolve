#ifndef FILE_BESSELJS3D
#define FILE_BESSELJS3D

/*
  Spherical Bessel functions j_n for the multipole expansions.
  Scalar types are instantiated in mptools.cpp, SIMD types instantiate here.
*/

#include <array>
#include <cstring>
#include <limits>
#include <ngstd.hpp>
#include <bla.hpp>

namespace ngsbem
{
  using namespace ngbla;

  /*
    Scaled spherical Bessel functions of the first kind and their derivatives,
    for 0 <= n <= nterms:

      fjs(n)   = j_n(z) / scale^n
      fjder(n) = d/dz fjs(n)        (only if fjder.Size() > 0)

    z real or complex (double, float, Complex, Complex32), or SIMD<double,N> / SIMD<Complex,N>
    for N arguments at once. Miller's algorithm: an upward recurrence finds the start order
    ntop where the minimal solution is negligible, the downward recurrence of the scaled
    g_n = j_n/scale^n runs from there, and j_0 or j_1 in closed form normalize the result.
    Over- and underflow of g is handled by an extended exponent: g is rescaled by exact
    powers of two and every entry keeps its own count, applied once at the end.
    SIMD: the lanes share the largest ntop, rescaling and exponent counts are per lane.

    Adapted from subroutine besseljs3d(nterms,z,scale,fjs,ifder,fjder) in FMM3D,
    src/Common/besseljs3d.f, https://github.com/flatironinstitute/FMM3D
    (found via https://fortran-lang.discourse.group/t/looking-for-spherical-bessel-and-hankel-functions-of-first-and-second-kind-and-arbitrary-order/2308/2)
  */
  template <typename T>
  NGS_DLL_HEADER
  void besseljs3d (int nterms, T z, double scale,
                   SliceVector<T> fjs, SliceVector<T> fjder = {0, 1, nullptr});

  template <typename T> inline auto BesselNorm2 (T z)
  {
    if constexpr (std::is_same_v<T,double> || std::is_same_v<T,float>) return z*z;
    else return z.real()*z.real()+z.imag()*z.imag();
  }

  template <typename T> inline T BesselLdexp (T v, int e)
  {
    if constexpr (std::is_same_v<T,double> || std::is_same_v<T,float>) return std::ldexp(v, e);
    else return T(std::ldexp(v.real(), e), std::ldexp(v.imag(), e));
  }

  template <typename Tz>
  void besseljs3d_scalar (int nterms, Tz z, double scale,
                          SliceVector<Tz> fjs, SliceVector<Tz> fjder)
  {
    using real_type = typename scal_traits<Tz>::TSCAL_REAL;
    const real_type rscale = scale;
    constexpr bool fp32 = std::is_same_v<real_type,float>;
    const real_type upper2 = fp32 ? 1e20f : 1e40;
    // extended exponent: g is rescaled by exact powers 2^(-+EXPO), the count is applied once at the end
    constexpr int EXPO = fp32 ? 40 : 128;
    const real_type big2 = std::ldexp(real_type(1), 2*EXPO);
    const real_type small2 = std::ldexp(real_type(1), -2*EXPO);
    const real_type down = std::ldexp(real_type(1), -EXPO);
    const real_type up = std::ldexp(real_type(1), EXPO);

    // For FP32, only exact zero may bypass the scaled small-z series.
    if (std::is_same_v<real_type,float> ? z == Tz(0) : abs(z) < 1e-200)
      {
        fjs(0) = 1;
        for (int i = 1; i <= nterms; i++)
          fjs(i) = 0.0;

        if (fjder.Size())
          {
            fjder = 0.0;
            if (nterms >= 1)
              fjder(1) = 1/(3*rscale);
          }
        return;
      }

    if constexpr (std::is_same_v<real_type,float>)
      if (abs(z) < 1e-2f)
        {
          Tz leading = 1;
          Tz dleading = 0;
          Tz zscaled = z/rscale;
          for (int n = 0; n <= nterms; n++)
            {
              Tz term = 1;
              Tz dterm = 0;
              Tz sum = 1;
              Tz dsum = 0;
              for (int k = 0; k < 3; k++)
                {
                  real_type denom = 2*(k+1)*(2*n+2*k+3);
                  Tz factor = -z*z/denom;
                  dterm = dterm*factor-term*real_type(2)*z/denom;
                  term *= factor;
                  sum += term;
                  dsum += dterm;
                }

              fjs(n) = leading*sum;
              if (fjder.Size())
                fjder(n) = dleading*sum+leading*dsum;

              // Advance independently of leading, which may already have underflowed.
              if (n == 0)
                dleading = Tz(1)/(3*rscale);
              else
                dleading *= real_type(n+1)*zscaled/real_type(n*(2*n+3));
              leading *= zscaled/real_type(2*n+3);
            }
          return;
        }


    // Step 1: upward recurrence from nterms until |j| grows by 1e20 (1e10 in float): ntop.
    // For |z| >> nterms this takes about |z| steps; capped at nterms+100000 as in FMM3D.

    Tz zinv=Tz(1.0)/z;

    int ntop;
    {
      Tz fjm1 = 0.0;
      Tz fj0 = 1.0;
      for (int i = nterms; ; i++)
        {
          Tz fj1 = real_type(2*i+1)*zinv*fj0-fjm1;
          if (BesselNorm2(fj1) > upper2)
            {
              ntop=i+1;
              break;
            }
          fjm1 = fj0;
          fj0 = fj1;
          if (i > nterms+100000)
            throw Exception("bessel failed 1");
        }
    }

    /*
      Step 2: downward recurrence for the scaled g_i = j_i/scale^i,
        g_{i-1} = (2i+1) scale/z g_i - scale^2 g_{i+1}.
      Entries 0..nterms go to fjs, nterms+1 to tail. On a rescale the current pair is
      multiplied by 2^(-+EXPO); ex[i] is the rescale count when entry i was stored, so
      the consistent value is stored(i) * 2^(-EXPO (cnt-ex[i])).
    */
    LocalHeap & lh = TLHeap();
    HeapReset hr(lh);
    FlatArray<int> ex(nterms+2, lh);

    Tz sz = rscale*zinv;
    real_type s2 = rscale*rscale;
    Tz tail = 0.0;
    auto store = [&](int k, Tz val, int cnt)
    {
      if (k > nterms+1) return;
      if (k <= nterms) fjs(k) = val; else tail = val;
      ex[k] = cnt;
    };

    Tz gp = 0.0, g = 1.0;   // g_{i+1}, g_i
    int cnt = 0;
    store(ntop, gp, 0);
    store(ntop-1, g, 0);
    for (int i = ntop-1; i >= 1; i--)
      {
        Tz gm = real_type(2*i+1)*sz*g - s2*gp;
        gp = g;
        g = gm;
        auto n2 = BesselNorm2(g);
        if (n2 > big2)
          {
            g *= down;
            gp *= down;
            store(i, gp, ++cnt);
          }
        else if (n2 < small2 && BesselNorm2(gp) < small2)
          {
            // shrinking for scale << |z|; both small, so not a zero crossing
            g *= up;
            gp *= up;
            store(i, gp, --cnt);
          }
        store(i-1, g, cnt);
      }

    auto consistent = [&](int k)
    {
      Tz v = (k <= nterms) ? fjs(k) : tail;
      return cnt != ex[k] ? BesselLdexp(v, -EXPO*(cnt-ex[k])) : v;
    };

    // Step 3: normalize with j_0 or j_1, whichever is larger
  
    Tz fj0=sin(z)*zinv;
    Tz fj1=fj0*zinv-cos(z)*zinv;
    if (abs(z) < 1e-2)
      {
        fj1 = real_type(1.0/3.0) * z - real_type(1.0/30.0) * z*z*z + real_type(1.0/840.0) * z*z*z*z*z;
        // |err| < 8/9! * 0.01**7  = 2e-19
      }
    
    Tz zscale;
    if (abs(fj1) > abs(fj0))
      zscale=fj1/(consistent(1)*rscale);
    else
      zscale=fj0/consistent(0);

    // Scale the jfuns by zscale, one rounding per entry. The shift is constant between
    // rescale points: one power-of-two factor per segment, ldexp only if it is subnormal.
    constexpr int maxe = std::numeric_limits<real_type>::max_exponent-3;
    int dlast = 0;
    real_type fac = 1;
    for (int i = 0; i <= nterms; i++)
      {
        Tz v = fjs(i)*zscale;
        int e = -EXPO*(cnt-ex[i]);
        if (e != dlast)
          {
            dlast = e;
            fac = std::ldexp(real_type(1), e);
          }
        fjs(i) = (e >= -maxe && e <= maxe) ? v*fac : BesselLdexp(v, e);
      }

    // derivatives: j_n' = n/(2n+1) j_{n-1} - (n+1)/(2n+1) j_{n+1}

    if (fjder.Size())
      {
        Tz v = tail*zscale;
        Tz fnext = cnt != ex[nterms+1] ? BesselLdexp(v, -EXPO*(cnt-ex[nterms+1])) : v;
        real_type scalinv = 1/rscale;
        fjder(0) = -(nterms >= 1 ? fjs(1) : fnext)*rscale;
        for (int i = 1; i <= nterms; i++)
          {
            real_type dc1 = real_type(i)/(2*i+1);
            real_type dc2 = 1-dc1;
            Tz fup = i < nterms ? fjs(i+1) : fnext;
            fjder(i) = (dc1*scalinv)*fjs(i-1) - (dc2*rscale)*fup;
          }
      }
  }

  template <typename T> struct BesselSIMDInfo { };
  template <typename S, int N> struct BesselSIMDInfo<SIMD<S,N>>
  {
    using scal = S;
    static constexpr int width = N;
  };

  // SIMD over z, see besseljs3d
  template <typename Tz>
  void besseljs3d_simd (int nterms, Tz z, double scale,
                        SliceVector<Tz> fjs, SliceVector<Tz> fjder)
  {
    constexpr int N = BesselSIMDInfo<Tz>::width;
    constexpr bool cplx = std::is_same_v<typename BesselSIMDInfo<Tz>::scal,Complex>;
    static_assert(cplx || std::is_same_v<typename BesselSIMDInfo<Tz>::scal,double>,
                  "besseljs3d: SIMD<double,N> or SIMD<Complex,N>");
    using TR = SIMD<double,N>;

    constexpr int EXPO = 128;
    constexpr int maxe = std::numeric_limits<double>::max_exponent-3;
    const double rscale = scale;
    const double upper2 = 1e40;
    const TR big2 = std::ldexp(1.0, 2*EXPO), small2 = std::ldexp(1.0, -2*EXPO);
    const TR down = std::ldexp(1.0, -EXPO), up = std::ldexp(1.0, EXPO);

    auto norm2 = [](Tz v) -> TR
    {
      if constexpr (cplx) return v.real()*v.real()+v.imag()*v.imag();
      else return v*v;
    };
    auto sel = [](SIMD<mask64,N> m, Tz a, Tz b) -> Tz
    {
      if constexpr (cplx) return Tz(If(m, a.real(), b.real()), If(m, a.imag(), b.imag()));
      else return If(m, a, b);
    };
    auto any = [](SIMD<mask64,N> m) { return HSum(If(m, TR(1.0), TR(0.0))) != 0; };

    // v * 2^(-EXPO d) lane by lane; exact tells whether all factors are normal
    auto pow2 = [&](TR d, bool & exact) -> TR
    {
      std::array<double,N> f;
      exact = true;
      for (int l = 0; l < N; l++)
        {
          int e = -EXPO*int(d[l]);
          if (e < -maxe || e > maxe) exact = false;
          f[l] = std::ldexp(1.0, e);
        }
      return TR(f.data());
    };
    auto ldexp_lanes = [&](Tz v, TR d) -> Tz
    {
      std::array<double,N> re, im;
      for (int l = 0; l < N; l++)
        {
          int e = -EXPO*int(d[l]);
          if constexpr (cplx)
            {
              re[l] = std::ldexp(v.real()[l], e);
              im[l] = std::ldexp(v.imag()[l], e);
            }
          else
            re[l] = std::ldexp(v[l], e);
        }
      if constexpr (cplx) return Tz(TR(re.data()), TR(im.data()));
      else return TR(re.data());
    };
    auto shift = [&](Tz v, TR d) -> Tz
    {
      bool exact;
      TR fac = pow2(d, exact);
      return exact ? v*fac : ldexp_lanes(v, d);
    };

    // lanes with |z| < 1e-200 run with z = 1 and are overwritten at the end
    TR az2 = norm2(z);
    auto tiny = [&]
    {
      if constexpr (cplx) return (fabs(z.real()) < TR(1e-200)) && (fabs(z.imag()) < TR(1e-200));
      else return fabs(z) < TR(1e-200);
    }();
    Tz zz = sel(tiny, Tz(1.0), z);
    Tz zinv = Tz(1.0)/zz;

    // Step 1: ntop = the largest start order over the lanes
    int ntop;
    {
      Tz fjm1 = 0.0, fj0 = 1.0;
      TR done = 0.0;
      for (int i = nterms; ; i++)
        {
          Tz fj1 = double(2*i+1)*zinv*fj0-fjm1;
          done = If(norm2(fj1) > TR(upper2), TR(1.0), done);
          if (HSum(done) == N)
            {
              ntop = i+1;
              break;
            }
          // finished lanes stay put, no overflow
          auto run = done == TR(0.0);
          fjm1 = sel(run, fj0, fjm1);
          fj0 = sel(run, fj1, fj0);
          if (i > nterms+100000)
            throw Exception("bessel failed 1");
        }
    }

    // Step 2: scaled downward recurrence, per-lane rescaling as in the scalar version
    LocalHeap & lh = TLHeap();
    HeapReset hr(lh);
    FlatArray<TR> ex(nterms+2, lh);

    Tz sz = rscale*zinv;
    double s2 = rscale*rscale;
    Tz tail = 0.0;
    auto store = [&](int k, Tz val, TR cnt)
    {
      if (k > nterms+1) return;
      if (k <= nterms) fjs(k) = val; else tail = val;
      ex[k] = cnt;
    };

    // Rescale checks every K steps only. Per step |g| grows at most by
    // G = (2 ntop+1) scale/|z|_min + scale^2 and shrinks at most by Q = |z|_max/scale,
    // K is chosen such that K steps stay within the double range around 2^(-+EXPO).
    int K;
    {
      double zmin = std::numeric_limits<double>::max(), zmax = 0;
      TR az = sqrt(norm2(zz));
      for (int l = 0; l < N; l++)
        {
          zmin = std::min(zmin, az[l]);
          zmax = std::max(zmax, az[l]);
        }
      double lg = std::max({ std::log2((2*ntop+1)*rscale/zmin + s2), std::log2(zmax/rscale), 1.0 });
      K = std::clamp(int(800/lg), 1, 8);
    }

    Tz gp = 0.0, g = 1.0;
    TR cnt = 0.0;
    store(ntop, gp, cnt);
    store(ntop-1, g, cnt);
    for (int i = ntop-1, k = 1; i >= 1; i--, k++)
      {
        Tz gm = double(2*i+1)*sz*g - s2*gp;
        gp = g;
        g = gm;
        if (k == K || i == 1)
          {
            // K steps can move a lane by more than one factor 2^EXPO: repeat
            k = 0;
            for (bool more = true; more; )
              {
                TR n2 = norm2(g);
                auto mdown = n2 > big2;
                auto mup = (n2 < small2) && (norm2(gp) < small2);
                TR f = If(mdown, down, If(mup, up, TR(1.0)));
                more = any(f != TR(1.0));
                if (!more) break;
                g = g*f;
                gp = gp*f;
                cnt = cnt + If(mdown, TR(1.0), If(mup, TR(-1.0), TR(0.0)));
              }
            if (i <= nterms+1)
              store(i, gp, cnt);
          }
        store(i-1, g, cnt);
      }

    auto consistent = [&](int k)
    {
      return shift((k <= nterms) ? fjs(k) : tail, cnt-ex[k]);
    };

    // Step 3: normalize with j_0 or j_1, whichever is larger
    Tz sn, cs;
    if constexpr (cplx)
      {
        sn = sin(zz);
        cs = cos(zz);
      }
    else
      {
        auto [s,c] = sincos(zz);
        sn = s;
        cs = c;
      }
    Tz fj0 = sn*zinv;
    Tz fj1 = fj0*zinv-cs*zinv;
    Tz ser = (1.0/3.0)*zz - (1.0/30.0)*zz*zz*zz + (1.0/840.0)*zz*zz*zz*zz*zz;
    fj1 = sel(az2 < TR(1e-4), ser, fj1);
    Tz zscale = sel(norm2(fj1) > norm2(fj0), fj1/(consistent(1)*rscale), fj0/consistent(0));

    // one power-of-two factor per segment, recomputed where the counts change
    TR fac = 1.0;
    bool exact = true;
    for (int i = 0; i <= nterms; i++)
      {
        TR d = cnt-ex[i];
        Tz v = fjs(i)*zscale;
        if (i == 0 || std::memcmp(&ex[i], &ex[i-1], sizeof(TR)) != 0)
          fac = pow2(d, exact);
        fjs(i) = exact ? v*fac : ldexp_lanes(v, d);
      }

    if (fjder.Size())
      {
        Tz fnext = shift(tail*zscale, cnt-ex[nterms+1]);
        double scalinv = 1/rscale;
        fjder(0) = -(nterms >= 1 ? fjs(1) : fnext)*rscale;
        for (int i = 1; i <= nterms; i++)
          {
            double dc1 = double(i)/(2*i+1);
            double dc2 = 1-dc1;
            Tz fup = i < nterms ? fjs(i+1) : fnext;
            fjder(i) = (dc1*scalinv)*fjs(i-1) - (dc2*rscale)*fup;
          }
      }

    if (any(tiny))
      for (int i = 0; i <= nterms; i++)
        {
          fjs(i) = sel(tiny, Tz(i == 0 ? 1.0 : 0.0), fjs(i));
          if (fjder.Size())
            fjder(i) = sel(tiny, Tz(i == 1 ? 1/(3*rscale) : 0.0), fjder(i));
        }
  }

  template <typename Tz>
  void besseljs3d (int nterms, Tz z, double scale,
                   SliceVector<Tz> fjs, SliceVector<Tz> fjder)
  {
    if constexpr (requires { BesselSIMDInfo<Tz>::width; })
      besseljs3d_simd (nterms, z, scale, fjs, fjder);
    else
      besseljs3d_scalar (nterms, z, scale, fjs, fjder);
  }

  extern template void besseljs3d (int, double, double, SliceVector<double>, SliceVector<double>);
  extern template void besseljs3d (int, float, double, SliceVector<float>, SliceVector<float>);
  extern template void besseljs3d (int, Complex, double, SliceVector<Complex>, SliceVector<Complex>);
  extern template void besseljs3d (int, Complex32, double, SliceVector<Complex32>, SliceVector<Complex32>);
}

#endif
