/*********************************************************************/
/* File:   p2m_gpu_kernel.h                                          */
/* Purpose: P2M of the singular leaves (charges and dipoles) on the  */
/*          common GPU layer (ngstd/gpukernel.h)                     */
/*********************************************************************/

/*
  SphericalExpansion<Singular>::AddCharge and AddDipole of REAL expansions: per source the basis
  T(n,m) = conj(Y_n^m) j_n(kappa r)/s^n, accumulated with the charge weight i kappa c and, for dipoles, with the
  weight i kappa c2 times the direction components (kind 2) or i kappa times the dipole vector's (kind bit 4:
  per source and component, complex); the dipole moments then get the truncated
  DirectionalDiffAdd as a gather, with the second kappa of AddDipole in its factors kappa/s and kappa*s (so no
  kappa^2 c2 underflows for Laplace, kappa = 1e-16). One group per (leaf, component); every coefficient belongs to
  one work-item, which adds the leaf's sources in their order, so the result is deterministic. P2M_PMAX, P2M_KMAX
  (coefficients per work-item) and P2M_STAGE (sources per stage) are set by GetLeafKernels (mptools_device.cpp), which
  also substitutes REAL and KAPPA_T: REAL, or Complex<REAL> for a complex kappa (then geo holds |kappa| r
  and the kernels multiply by khat = kappa/|kappa|).
*/

#ifdef NGS_GPU_CPU
using std::sin, std::cos, std::sqrt, std::fabs, std::exp;
#endif

#define P2M_W (P2M_PMAX+1)
#define P2M_NC (P2M_W*P2M_W)
#define P2M_TRI (P2M_W*(P2M_PMAX+2)/2)
#define P2M_T(n,m) ((n)*((n)+1)/2+(m))

// the tables, must match P2MTables: NormalizedLegendreFunctions' A(n,m), B(n,m) (W x W), sqrt((2m-1)/(2m)) and
// sqrt(2m+1) (W each), CalcAmn(|m|,n) ((W+1) x (W+1)), CalcBmn(m,n) ((W+1) x (2W+1), m from -W)
#define P2M_TAB_A 0
#define P2M_TAB_B (P2M_W*P2M_W)
#define P2M_TAB_D (2*P2M_W*P2M_W)
#define P2M_TAB_S (P2M_TAB_D+P2M_W)
#define P2M_TAB_AMN (P2M_TAB_S+P2M_W)
#define P2M_TAB_BMN (P2M_TAB_AMN+(P2M_W+1)*(P2M_W+1))

GPU_FUNC REAL P2MAmn (GLOBAL_PTR(const REAL) tab, int m, int n)
{
  return tab[P2M_TAB_AMN + n*(P2M_W+1) + (m < 0 ? -m : m)];
}

GPU_FUNC REAL P2MBmn (GLOBAL_PTR(const REAL) tab, int m, int n)
{
  return tab[P2M_TAB_BMN + n*(2*P2M_W+1) + m + P2M_W];
}

// the kappa dependent parts for a real or complex kappa: i kappa, kappa r, sin, cos, |.|
GPU_FUNC Complex<REAL> IKappa (REAL k) { return Complex<REAL>(0, k); }
GPU_FUNC Complex<REAL> IKappa (Complex<REAL> k) { return Complex<REAL>(-k.im, k.re); }
GPU_FUNC REAL KappaR (REAL, REAL kr) { return kr; }
GPU_FUNC Complex<REAL> KappaR (Complex<REAL> khat, REAL kr) { return khat*kr; }
GPU_FUNC REAL KSin (REAL z) { return sin(z); }
GPU_FUNC REAL KCos (REAL z) { return cos(z); }
// sinh by its series below 1/2 (a fast-math sinh from exp cancels there; to b^13 in double), cosh from exp
GPU_FUNC REAL KSinh (REAL b)
{
  REAL a = fabs(b), b2 = b*b;
  if (a < REAL(0.5))
    {
      if (sizeof(REAL) == 8)
        return b*(REAL(1) + b2*(REAL(1)/6 + b2*(REAL(1)/120 + b2*(REAL(1)/5040 + b2*(REAL(1)/362880
                  + b2*(REAL(1)/39916800 + b2*(REAL(1)/(REAL(39916800)*156))))))));
      return b*(REAL(1) + b2*(REAL(1)/6 + b2*(REAL(1)/120 + b2*(REAL(1)/5040 + b2*(REAL(1)/362880)))));
    }
  REAL e = exp(a), r = REAL(0.5)*(e - REAL(1)/e);
  return b < 0 ? -r : r;
}
GPU_FUNC REAL KCosh (REAL b) { REAL e = exp(fabs(b)); return REAL(0.5)*(e + REAL(1)/e); }
GPU_FUNC Complex<REAL> KSin (Complex<REAL> z) { return Complex<REAL>(sin(z.re)*KCosh(z.im), cos(z.re)*KSinh(z.im)); }
GPU_FUNC Complex<REAL> KCos (Complex<REAL> z) { return Complex<REAL>(cos(z.re)*KCosh(z.im), -sin(z.re)*KSinh(z.im)); }
GPU_FUNC REAL KAbs (REAL z) { return fabs(z); }
GPU_FUNC REAL KAbs (Complex<REAL> z) { return abs(z); }

// f[n] = j_n(z)/s^n, n = 0..p: the ratios f_n/f_{n-1} by backward recursion from p+|z|+16 (40 in double),
// normalized with j_0, or with j_1 near the zeros of j_0 (|z| >= 2; below, j_1 from sin and cos cancels and
// |j_0| > |j_1| anyway)
template <typename Z>
GPU_FUNC void P2MBessel (int p, Z z, REAL s, LOCAL_PTR(Z) f)
{
  if (z == Z(0))
    {
      f[0] = Z(1);
      for (int n = 1; n <= p; n++) f[n] = Z(0);
      return;
    }
  Z s2z = s*s*z, rho = Z(0);
  for (int n = p+(sizeof(REAL) == 8 ? 40 : 16)+int(KAbs(z)); n >= 1; n--)
    {
      Z den = REAL(2*n+1)*s - s2z*rho;
      rho = z / (den != Z(0) ? den : Z(REAL(1e-30)));
      if (n <= p) f[n] = rho;
    }
  Z j0 = KSin(z)/z, j1 = (j0 - KCos(z))/z;
  if (p == 0 || KAbs(z) < 2 || KAbs(j0) >= KAbs(j1))
    {
      f[0] = j0;
      for (int n = 1; n <= p; n++) f[n] = f[n-1]*f[n];
    }
  else
    {
      Z f1 = j1/s;
      f[0] = f1/f[1];
      f[1] = f1;
      for (int n = 2; n <= p; n++) f[n] = f[n-1]*f[n];
    }
}

// sources [lfirst, lfirst+lcount) of leaf GROUP_ID_X into component GROUP_ID_Y of pool slot lslot:
// geo = cos theta, sin theta, cos phi, sin phi, kappa r relative to the leaf center; dir = the dipole direction,
// dipvec[(src*comps+comp)*3+j] the dipole vector
KERNEL(fmm_p2m, GLOBAL_IN(int,lfirst), GLOBAL_IN(int,lcount), GLOBAL_IN(int,lslot),
       GLOBAL_IN(REAL,geo), GLOBAL_IN(REAL,dir), GLOBAL_IN(Complex<REAL>,charge), GLOBAL_IN(Complex<REAL>,dipole),
       GLOBAL_IN(Complex<REAL>,dipvec), GLOBAL_IN(REAL,tab), GLOBAL(Complex<REAL>,pool),
       VALUE(int,p), VALUE(REAL,s), VALUE(KAPPA_T,kappa), VALUE(KAPPA_T,kappa_over_s), VALUE(KAPPA_T,kappa_times_s),
       VALUE(KAPPA_T,khat), VALUE(int,kind), VALUE(int,comps))
{
  int leaf = GROUP_ID_X, comp = GROUP_ID_Y, tid = LOCAL_ID_X, nth = GROUP_SIZE_X;
  int first = lfirst[leaf], count = lcount[leaf];
  int nc = (p+1)*(p+1);
  bool charges = (kind & 1) != 0, dipoles = (kind & 2) != 0;

  const int KW = sizeof(KAPPA_T)/sizeof(REAL);
  const int STAGE = P2M_STAGE*(P2M_TRI + (KW+2)*P2M_W + 8);
  const int MOMENTS = 6*P2M_NC;
  SHARED(REAL, shm, STAGE > MOMENTS ? STAGE : MOMENTS);
  LOCAL_PTR(REAL) leg = shm;                                              // [b][P2M_T(n,m)] P-bar_n^m, with sqrt(2n+1)
  LOCAL_PTR(KAPPA_T) rad = (LOCAL_PTR(KAPPA_T))(shm + P2M_STAGE*P2M_TRI); // [b][n] j_n/s^n
  // [b][m] (-e^{i phi})^m / sqrt(4 pi), then [b] w0, w+, w-, wz
  LOCAL_PTR(Complex<REAL>) ph = (LOCAL_PTR(Complex<REAL>))(shm + P2M_STAGE*(P2M_TRI + KW*P2M_W));
  LOCAL_PTR(Complex<REAL>) wt = ph + P2M_STAGE*P2M_W;

  int jn[P2M_KMAX], jm[P2M_KMAX];
  Complex<REAL> a0[P2M_KMAX], ap[P2M_KMAX], am[P2M_KMAX], az[P2M_KMAX];
  for (int k = 0; k < P2M_KMAX; k++)
    {
      int j = tid + k*nth;
      int n = int(sqrt(REAL(j)));
      while (n*n > j) n--;
      while ((n+1)*(n+1) <= j) n++;
      jn[k] = n; jm[k] = j-n*n-n;
      a0[k] = ap[k] = am[k] = az[k] = Complex<REAL>(0, 0);
    }

  for (int b0 = 0; b0 < count; b0 += P2M_STAGE)
    {
      int nb = count-b0 < P2M_STAGE ? count-b0 : P2M_STAGE;
      // work-item (b, m): column m of source b's Legendre table; m = p: its radial part, m = 0: phases and weights
      for (int t = tid; t < nb*(p+1); t += nth)
        {
          int b = t/(p+1), m = t%(p+1);
          int src = first+b0+b;
          REAL ct = geo[5*src], mst = -geo[5*src+1];
          LOCAL_PTR(REAL) lb = leg + b*P2M_TRI;
          REAL diag = 1;
          for (int k = 1; k <= m; k++) diag = diag*mst*tab[P2M_TAB_D+k];
          REAL y0 = diag;
          lb[P2M_T(m,m)] = y0*tab[P2M_TAB_S+m];
          if (m < p)
            {
              REAL y1 = ct*diag*tab[P2M_TAB_S+m];
              lb[P2M_T(m+1,m)] = y1*tab[P2M_TAB_S+m+1];
              for (int n = m+2; n <= p; n++)
                {
                  REAL y2 = tab[P2M_TAB_A+n*P2M_W+m]*ct*y1 - tab[P2M_TAB_B+n*P2M_W+m]*y0;
                  lb[P2M_T(n,m)] = y2*tab[P2M_TAB_S+n];
                  y0 = y1; y1 = y2;
                }
            }
          if (m == p)
            P2MBessel (p, KappaR(khat, geo[5*src+4]), s, rad + b*P2M_W);
          if (m == 0)
            {
              LOCAL_PTR(Complex<REAL>) pb = ph + b*P2M_W;
              Complex<REAL> em(REAL(0.28209479177387814), 0), me(-geo[5*src+2], -geo[5*src+3]);
              for (int mm = 0; mm <= p; mm++)
                {
                  pb[mm] = em;
                  em = em*me;
                }
              Complex<REAL> w0(0, 0), wp(0, 0), wm(0, 0), wz(0, 0);
              if (charges)
                w0 = IKappa(kappa)*charge[src*comps+comp];
              if (dipoles)
                {
                  Complex<REAL> w = IKappa(kappa)*dipole[src*comps+comp];
                  REAL dx = dir[3*src], dy = dir[3*src+1], dz = dir[3*src+2];
                  wp = Complex<REAL>(REAL(0.5)*dx, REAL(0.5)*dy)*w;
                  wm = Complex<REAL>(REAL(0.5)*dx, REAL(-0.5)*dy)*w;
                  wz = dz*w;
                }
              if (kind & 4)
                {
                  GLOBAL_PTR(const Complex<REAL>) D = dipvec + 3*(src*comps+comp);
                  Complex<REAL> ik = IKappa(kappa), iDy = Complex<REAL>(-D[1].im, D[1].re);
                  wp = wp + ik*(REAL(0.5)*(D[0]+iDy));
                  wm = wm + ik*(REAL(0.5)*(D[0]-iDy));
                  wz = wz + ik*D[2];
                }
              LOCAL_PTR(Complex<REAL>) wb = wt + 4*b;
              wb[0] = w0; wb[1] = wp; wb[2] = wm; wb[3] = wz;
            }
        }
      BARRIER();
      for (int k = 0; k < P2M_KMAX; k++)
        {
          int n = jn[k], m = jm[k];
          if (n > p) continue;
          int ma = m < 0 ? -m : m;
          for (int b = 0; b < nb; b++)
            {
              KAPPA_T lr = leg[b*P2M_TRI+P2M_T(n,ma)] * rad[b*P2M_W+n];
              Complex<REAL> e = ph[b*P2M_W+ma];
              Complex<REAL> T = m >= 0 ? conj(e)*lr : e*lr;
              LOCAL_PTR(Complex<REAL>) wb = wt + 4*b;
              a0[k] = a0[k] + T*wb[0];
              ap[k] = ap[k] + T*wb[1];
              am[k] = am[k] + T*wb[2];
              az[k] = az[k] + T*wb[3];
            }
        }
      BARRIER();
    }

  GLOBAL_PTR(Complex<REAL>) out = pool + (lslot[leaf]*comps+comp)*nc;
  if ((kind & 6) == 0)
    {
      for (int k = 0; k < P2M_KMAX; k++)
        {
          int j = tid + k*nth;
          if (j < nc) out[j] = a0[k];
        }
      return;
    }

  // DirectionalDiffAdd(d, res, s) of kappa times the dipole moments, d = (w+ + w-, (w+ - w-)/i, wz), gathered per coefficient
  LOCAL_PTR(Complex<REAL>) Mp = (LOCAL_PTR(Complex<REAL>))(shm);
  LOCAL_PTR(Complex<REAL>) Mm = Mp + P2M_NC;
  LOCAL_PTR(Complex<REAL>) Mz = Mp + 2*P2M_NC;
  for (int k = 0; k < P2M_KMAX; k++)
    {
      int j = tid + k*nth;
      if (j >= nc) continue;
      Mp[j] = ap[k]; Mm[j] = am[k]; Mz[j] = az[k];
    }
  BARRIER();
  for (int k = 0; k < P2M_KMAX; k++)
    {
      int j = tid + k*nth;
      if (j >= nc) continue;
      int n = jn[k], m = jm[k];
      Complex<REAL> sum = a0[k];
      if (n >= 1)
        {
          int n1 = n-1;
          if (m+1 <= n1) sum = sum + (kappa_over_s*P2MBmn(tab, m, n))*Mp[n1*n1+n1+m+1];
          if (m <= n1 && m >= -n1) sum = sum - (kappa_over_s*P2MAmn(tab, m, n1))*Mz[n1*n1+n1+m];
          if (m-1 >= -n1) sum = sum + (kappa_over_s*P2MBmn(tab, -m, n))*Mm[n1*n1+n1+m-1];
        }
      if (n < p)
        {
          int n2 = n+1;
          sum = sum - (kappa_times_s*P2MBmn(tab, m-1, n2))*Mm[n2*n2+n2+m-1];
          sum = sum + (kappa_times_s*P2MAmn(tab, m, n))*Mz[n2*n2+n2+m];
          sum = sum - (kappa_times_s*P2MBmn(tab, -m-1, n2))*Mp[n2*n2+n2+m+1];
        }
      out[j] = sum;
    }
}
