#ifndef FILE_BIORTHO_UTILS
#define FILE_BIORTHO_UTILS

// biorthogonal-type trig shapes (Haubold-Beuchler-Schoeberl), shared by HCurl and HDiv

#include "recursive_pol.hpp"
#include "hcurlfe_utils.hpp"

namespace ngfem
{
  // 1/k, sqrt(k), 1/sqrt(k) for biortho trig coefficients
  struct BiorthoTabs
  {
    std::array<double,1024> inv{}, sq{}, isq{};
    BiorthoTabs ()
    {
      for (size_t k = 1; k < inv.size(); k++)
        { inv[k] = 1.0/k; sq[k] = sqrt(double(k)); isq[k] = 1.0/sq[k]; }
    }
  };
  inline const BiorthoTabs biortho_tabs;

  // biortho trig w~ coefficients U = c0 g - c1 gm1, V = c2 g - c3 gm1, index [i][n], p < 64
  struct BiorthoCoefs {
    static constexpr int N = 64;
    std::array<std::array<double,4>,N*N> c{};
    BiorthoCoefs () {
      for (int i = 2; i < N; i++)
        for (int n = 0; i+n < N; n++) {
          double t = sqrt(2.*i-1)/sqrt(2.*n+2*i)/(i-1);
          double al = t*(n+2*i)/(n+i+1), be = t*n/(n+i);
          c[i*N+n] = { al*(n+1), be*n, al*i, be*i };
        }
    }
  };
  inline const BiorthoCoefs biortho_coefs;

  // biortho trig non-gradient shapes, homogeneous in la,lb,lc (t = la+lb+lc), times mz
  // C collapse vertex: g ~ lc IntJac^(2i), f ~ la lb P^(1,1)_{i-2}
  template <typename Tx, typename TZ, typename TFA>
  INLINE void BiorthoTrigShapes (int p, Tx la, Tx lb, Tx lc, TZ mz, int & ii, TFA & shape)
  {
    Tx s = la+lb, t = s+lc;
    JacobiPolynomialFix<1,1>::EvalScaledMult
      (p-2, lb-la, s, la*lb, SBLambda([&](int m, Tx f) LAMBDA_INLINE
      {
        int i = m+2;
        Tx gm1(0.0);
        IntegratedJacobiPolynomialAlpha(2*i).EvalScaledMult
          (p-2-m, 2*lc-t, t, lc, SBLambda([&](int n, Tx g) LAMBDA_INLINE
          {
            // |curl w~|_ref = 1
            auto & cf = biortho_coefs.c[i*BiorthoCoefs::N+n];
            if constexpr (std::is_same_v<TZ,double>)
              shape[ii++] = uDw_minus_wDv (cf[0]*g - cf[1]*gm1, f, cf[2]*g - cf[3]*gm1);
            else
              shape[ii++] = zuDw_minus_wDv (cf[0]*g - cf[1]*gm1, f, cf[2]*g - cf[3]*gm1, mz);
            gm1 = g;
          }));
      }));
    // Whitney(A,B) * lc P^(2,1)_{j-1}, |curl|_ref = 1
    IntegratedJacobiPolynomialAlpha(3).EvalScaledMult
      (p-2, 2*lc-t, t, lc, SBLambda([&](int j, Tx g) LAMBDA_INLINE
      {
        shape[ii++] = wuDv_minus_wvDu (la, lb, (biortho_tabs.sq[2*j+4]*mz)*g);
      }));
  }

  // closed-form duals of the biortho trig shapes (SINUM 2024 (4.12)), for type1 or nograds
  // (default: top grad duals would leave RT_{p-1}); trafo maps 2D coords of pA,pB,pC
  // xs/s = (lb-la)/(la+lb), yc = 2lc-1; D_II diagonal, not normalized
  template <typename T, typename TM, typename TFA>
  INLINE void BiorthoTrigDuals (int p, bool usegrad, bool type1, T la, T lb, T lc, const TM & trafo,
                                Vec<2,T> pA, Vec<2,T> pB, Vec<2,T> pC, int & ii, TFA & shape)
  {
    if (p < 2) return;
    T xs = lb-la, s = la+lb, yc = 2*lc-1;
    constexpr int DIMS = Height<TM>();
    typedef Vec<DIMS,T> TV;
    typedef AutoDiff<1,T> AD1;
    auto lin = [](T a, const TV & u, T b, const TV & v) LAMBDA_INLINE
      { TV r; for (int k = 0; k < DIMS; k++) r(k) = a*u(k)+b*v(k); return r; };
    TV tex = trafo * Vec<2,T>(0.5 * (pB-pA));
    TV tey = trafo * Vec<2,T>(0.5 * (pC-0.5*(pA+pB)));

    // order: grads, w~ (k >= 1), type III (k = 0)
    int pg = usegrad ? p-2-(type1 ? 1 : 0) : -1;
    int ig = ii, iw = ii + (pg+1)*(pg+2)/2, i3 = iw + (p-1)*p/2;
    ArrayMem<T,20> pj(p);

    AD1 Lm1(0.0), Lm2(0.0);   // Ls_{i-1}, Ls_{i-2} with d/dxs
    LegendrePolynomial::EvalScaled1Assign
      (p, AD1(xs,0), AD1(s), SBLambda([&](int i, AD1 Li) LAMBDA_INLINE
      {
        if (i >= 1)
          {
            // rot (Ls_k(xs) P^(2k+1,0)_n(yc)), k = i-1
            int k = i-1;
            TV dq = (2.0*s*Lm2.DValue(0)) * tex - (4.0*Lm1.DValue(0)) * tey;
            TV dqp = (4.0*Lm1.Value()) * tex;
            int & ir = (k == 0) ? i3 : iw;
            JacobiPolynomialAlpha(2*k+1).Eval
              (p-1-k, AD1(yc,0), SBLambda([&](int n, AD1 Q) LAMBDA_INLINE
              {
                if (k > 0 || n > 0)
                  shape[ir++] = lin(Q.Value(), dq, Q.DValue(0), dqp);
              }));
          }
        if (i >= 2 && i-2 <= pg)
          {
            // grads: i B + j/(2(j+2i-1)) C, B ~ Ls_{i-1} P^(2i-1,1)_{j-1}, C ~ P^(2i,0)_{j-1}
            int nmax = pg-(i-2);
            TV bg = (double(i)*Lm1.Value()) * tex;
            TV c = (2.0*Li.DValue(0)) * tex - (4.0*s*Lm1.DValue(0)) * tey;
            JacobiPolynomialAlpha(2*i).Eval (nmax, yc, pj);
            IntegratedJacobiPolynomialAlpha(2*i).Eval
              (nmax, yc, SBLambda([&](int n, T R) LAMBDA_INLINE
              {
                double j = n+1;
                shape[ig++] = lin(j*R, bg, (j*0.5*biortho_tabs.inv[n+2*i])*pj[n], c);
              }));
          }
        Lm2 = Lm1;
        Lm1 = Li;
      }));
    ii = i3;
  }
}

#endif
