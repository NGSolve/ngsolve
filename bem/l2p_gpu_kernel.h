/*********************************************************************/
/* File:   l2p_gpu_kernel.h                                          */
/* Purpose: L2P, M2P and P2P at the targets (values, gradients) on   */
/*          the common GPU layer (ngstd/gpukernel.h)                 */
/*********************************************************************/

/*
  SphericalExpansion::Eval of REAL expansions at the targets, sum_n R_n sum_m Y_n^m c_nm (SphericalHarmonics::EvalOrders):
  fmm_l2p with R_n = j_n(kappa r)/s^n for the local expansion of a target's leaf, fmm_m2p with R_n = h_n(kappa r) s^n for
  the far singular nodes of SingularMLExpansion::Node::EvaluateMP; fmm_p2p adds the near sources directly. Compiled after
  p2m_gpu_kernel.h (its tables and P2MBessel); L2P_THREADS bounds the group size (set by GetLeafKernels in
  mptools_device.cpp). Every target value belongs to one work-item per launch: fmm_l2p writes it, the fmm_m2p and fmm_p2p
  launches after it add to it, one after the other. REAL, KAPPA_T, IKappa, KappaR, KSin and KCos as in p2m_gpu_kernel.h.
  mode: bit 0 the value, bit 1 the gradient (EvalDirectionalDerivative, EvaluateMPDeriv, EvaluateDeriv); per target and
  component outs = (mode & 1) + 3 (mode >> 1 & 1) entries, the value first.
*/

// h_0(z) = e^{iz}/(iz) and h_1(z)/h_0(z) = 1/z - i
GPU_FUNC Complex<REAL> HankelH0 (REAL z) { return Complex<REAL>(sin(z)/z, -cos(z)/z); }
GPU_FUNC Complex<REAL> HankelH0 (Complex<REAL> z)
{
  REAL e = exp(-z.im);
  return Complex<REAL>(e*sin(z.re), -e*cos(z.re)) / z;
}
GPU_FUNC Complex<REAL> HankelH1 (REAL z) { return Complex<REAL>(REAL(1)/z, -1); }
GPU_FUNC Complex<REAL> HankelH1 (Complex<REAL> z) { return Complex<REAL>(1)/z - Complex<REAL>(0, 1); }

// h_n(z) s^n, n = 0..p, as SphericalHankel1 (upward recursion)
template <typename Z>
GPU_FUNC void L2PHankel (int p, Z z, REAL s, LOCAL_PTR(Complex<REAL>) h)
{
  Complex<REAL> h0 = HankelH0(z);
  h[0] = h0;
  if (p < 1) return;
  Complex<REAL> h1 = (s*h0)*HankelH1(z);
  h[1] = h1;
  Z sz = s/z;
  REAL s2 = s*s;
  for (int n = 1; n < p; n++)
    {
      Complex<REAL> h2 = (REAL(2*n+1)*sz)*h1 - s2*h0;
      h[n+1] = h2;
      h0 = h1; h1 = h2;
    }
}

// R_n times the Legendre value y and its normalization t
GPU_FUNC Complex<REAL> L2PWeight (REAL r, REAL y, REAL t) { return Complex<REAL>(r*y*t, 0); }
GPU_FUNC Complex<REAL> L2PWeight (Complex<REAL> r, REAL y, REAL t) { return r*(y*t); }

struct L2PGrad3 { Complex<REAL> x, y, z; };

// the derivative expansions of DirectionalDiffAdd (mptools.cpp, Gumerov and Duraiswami 2.2.12) along w+ = (x+iy)/2,
// w- = (x-iy)/2 and z at (n,m), gathered from coef as the P2M gather: kis = kappa/scale on the degree n-1 terms,
// ks = kappa*scale on the degree n+1 terms (the factor kappa of EvalDirectionalDerivative folded in)
template <typename C, typename K>
GPU_FUNC L2PGrad3 L2PDiffCoef (C coef, GLOBAL_PTR(const REAL) tab, int p, int n, int m, K ks, K kis)
{
  L2PGrad3 d = { Complex<REAL>(0, 0), Complex<REAL>(0, 0), Complex<REAL>(0, 0) };
  int ma = m < 0 ? -m : m;
  if (n >= 1)
    {
      int n1 = n-1;
      if (m+1 <= n1) d.x = d.x + (kis*P2MBmn(tab, m, n))*coef[n1*n1+n1+m+1];
      if (m-1 >= -n1) d.y = d.y + (kis*P2MBmn(tab, -m, n))*coef[n1*n1+n1+m-1];
      if (ma <= n1) d.z = d.z - (kis*P2MAmn(tab, m, n1))*coef[n1*n1+n1+m];
    }
  if (n < p)
    {
      int n2 = n+1;
      d.x = d.x - (ks*P2MBmn(tab, -m-1, n2))*coef[n2*n2+n2+m+1];
      d.y = d.y - (ks*P2MBmn(tab, m-1, n2))*coef[n2*n2+n2+m-1];
      d.z = d.z + (ks*P2MAmn(tab, m, n))*coef[n2*n2+n2+m];
    }
  return d;
}

// mode bit 0: sum_n R_n sum_m Y_n^m c_nm, the Legendre columns from cos and -sin theta, me = -e^{i phi}, R real or
// complex; bit 1: its gradient, from the sums P, M, Z of R_n Y_n^m times L2PDiffCoef, grad = ((P+M)/2, i(P-M)/2, Z)
struct L2PSums { Complex<REAL> value; L2PGrad3 grad; };

template <typename RT, typename C, typename K>
GPU_FUNC L2PSums L2PEval (int mode, int p, REAL ct, REAL mst, Complex<REAL> me, LOCAL_PTR(RT) R, C coef,
                          GLOBAL_PTR(const REAL) tab, K ks, K kis)
{
  Complex<REAL> zero(0, 0), em(REAL(0.28209479177387814), 0), value = zero;
  L2PGrad3 sum = { zero, zero, zero };
  REAL diag = 1;
  for (int m = 0; m <= p; m++)
    {
      if (m > 0) { diag = diag*mst*tab[P2M_TAB_D+m]; em = em*me; }
      Complex<REAL> vp = zero, vn = zero;
      L2PGrad3 sp = { zero, zero, zero }, sn = { zero, zero, zero };
      REAL y0 = diag, y1 = ct*diag*tab[P2M_TAB_S+m];
      for (int n = m; n <= p; n++)
        {
          REAL y = y0;
          if (n == m+1) y = y1;
          else if (n > m+1)
            {
              y = tab[P2M_TAB_A+n*P2M_W+m]*ct*y1 - tab[P2M_TAB_B+n*P2M_W+m]*y0;
              y0 = y1; y1 = y;
            }
          Complex<REAL> w = L2PWeight (R[n], y, tab[P2M_TAB_S+n]);
          if (mode & 1)
            {
              vp = vp + w*coef[n*n+n+m];
              if (m > 0) vn = vn + w*coef[n*n+n-m];
            }
          if (mode & 2)
            {
              L2PGrad3 d = L2PDiffCoef (coef, tab, p, n, m, ks, kis);
              sp.x = sp.x + w*d.x; sp.y = sp.y + w*d.y; sp.z = sp.z + w*d.z;
              if (m > 0)
                {
                  d = L2PDiffCoef (coef, tab, p, n, -m, ks, kis);
                  sn.x = sn.x + w*d.x; sn.y = sn.y + w*d.y; sn.z = sn.z + w*d.z;
                }
            }
        }
      if (mode & 1) { value = value + em*vp; if (m > 0) value = value + conj(em)*vn; }
      if (mode & 2)
        {
          sum.x = sum.x + em*sp.x; sum.y = sum.y + em*sp.y; sum.z = sum.z + em*sp.z;
          if (m > 0) { sum.x = sum.x + conj(em)*sn.x; sum.y = sum.y + conj(em)*sn.y; sum.z = sum.z + conj(em)*sn.z; }
        }
    }
  L2PSums res = { value, { REAL(0.5)*(sum.x+sum.y), Complex<REAL>(0, REAL(0.5))*(sum.x-sum.y), sum.z } };
  return res;
}

// the targets [lfirst, lfirst+lcount) of leaf GROUP_ID_X, component GROUP_ID_Y of pool slot lslot:
// out = the local expansion's value and gradient; geo = cos theta, sin theta, cos phi, sin phi, kappa r relative to the
// leaf center
KERNEL(fmm_l2p, GLOBAL_IN(int,lfirst), GLOBAL_IN(int,lcount), GLOBAL_IN(int,lslot), GLOBAL_IN(REAL,geo),
       GLOBAL_IN(REAL,tab), GLOBAL_IN(Complex<REAL>,pool), GLOBAL(Complex<REAL>,out),
       VALUE(int,p), VALUE(REAL,s), VALUE(KAPPA_T,khat), VALUE(KAPPA_T,ks), VALUE(KAPPA_T,kis), VALUE(int,mode),
       VALUE(int,comps))
{
  int leaf = GROUP_ID_X, comp = GROUP_ID_Y, tid = LOCAL_ID_X, nth = GROUP_SIZE_X;
  int first = lfirst[leaf], count = lcount[leaf], nc = (p+1)*(p+1);
  SHARED(REAL, shm, 2*P2M_NC);
  SHARED(REAL, rad, L2P_THREADS*P2M_W*(sizeof(KAPPA_T)/sizeof(REAL)));
  LOCAL_PTR(Complex<REAL>) cf = (LOCAL_PTR(Complex<REAL>))(shm);
  GLOBAL_PTR(const Complex<REAL>) c = pool + (lslot[leaf]*comps+comp)*nc;
  for (int j = tid; j < nc; j += nth) cf[j] = c[j];
  BARRIER();
  LOCAL_PTR(KAPPA_T) R = (LOCAL_PTR(KAPPA_T))(rad) + tid*P2M_W;
  int outs = (mode & 1) + 3*(mode >> 1 & 1);
  for (int t = tid; t < count; t += nth)
    {
      int i = first+t, o = (i*comps+comp)*outs;
      Complex<REAL> me(-geo[5*i+2], -geo[5*i+3]);
      P2MBessel (p, KappaR(khat, geo[5*i+4]), s, R);
      L2PSums r = L2PEval (mode, p, geo[5*i], -geo[5*i+1], me, R, cf, tab, ks, kis);
      if (mode & 1) out[o++] = r.value;
      if (mode & 2) { out[o] = r.grad.x; out[o+1] = r.grad.y; out[o+2] = r.grad.z; }
    }
}

// entry e = GLOBAL_ID_X adds the pairs [efirst, efirst+ecount) to its target etarget, component GROUP_ID_Y:
// pair q is pool slot pslot[q] with pgeo[5q..] = the direction and kappa r of the target relative to the node center
KERNEL(fmm_m2p, GLOBAL_IN(int,etarget), GLOBAL_IN(int,efirst), GLOBAL_IN(int,ecount), GLOBAL_IN(int,pslot),
       GLOBAL_IN(REAL,pgeo), GLOBAL_IN(REAL,tab), GLOBAL_IN(Complex<REAL>,pool), GLOBAL(Complex<REAL>,out),
       VALUE(int,p), VALUE(REAL,s), VALUE(KAPPA_T,khat), VALUE(KAPPA_T,ks), VALUE(KAPPA_T,kis), VALUE(int,mode),
       VALUE(int,comps), VALUE(int,entries))
{
  int e = GLOBAL_ID_X, comp = GROUP_ID_Y;
  int nc = (p+1)*(p+1), outs = (mode & 1) + 3*(mode >> 1 & 1);
  SHARED(REAL, rad, L2P_THREADS*2*P2M_W);
  if (e >= entries) return;
  LOCAL_PTR(Complex<REAL>) R = (LOCAL_PTR(Complex<REAL>))(rad) + LOCAL_ID_X*P2M_W;
  Complex<REAL> zero(0, 0), sum = zero;
  L2PGrad3 gsum = { zero, zero, zero };
  for (int q = efirst[e]; q < efirst[e]+ecount[e]; q++)
    {
      GLOBAL_PTR(const REAL) g = pgeo + 5*q;
      GLOBAL_PTR(const Complex<REAL>) c = pool + (pslot[q]*comps+comp)*nc;
      Complex<REAL> me(-g[2], -g[3]);
      L2PHankel (p, KappaR(khat, g[4]), s, R);
      L2PSums r = L2PEval (mode, p, g[0], -g[1], me, R, c, tab, ks, kis);
      sum = sum + r.value;
      gsum.x = gsum.x + r.grad.x; gsum.y = gsum.y + r.grad.y; gsum.z = gsum.z + r.grad.z;
    }
  int o = (etarget[e]*comps+comp)*outs;
  if (mode & 1)
    { out[o] = out[o] + sum; o++; }
  if (mode & 2)
    {
      out[o] = out[o] + gsum.x; out[o+1] = out[o+1] + gsum.y; out[o+2] = out[o+2] + gsum.z;
    }
}

// e^{i kappa rho}, for a complex kappa with the decay e^{-Im(kappa) rho} (as the host's PhaseFactor)
GPU_FUNC Complex<REAL> P2PPhase (REAL kappa, REAL rho) { REAL kr = kappa*rho; return Complex<REAL>(cos(kr), sin(kr)); }
GPU_FUNC Complex<REAL> P2PPhase (Complex<REAL> kappa, REAL rho)
{
  REAL kr = kappa.re*rho, decay = exp(-kappa.im*rho);
  return Complex<REAL>(cos(kr)*decay, sin(kr)*decay);
}

// |d| as the host's REAL(L2Norm(d)): 0 only for d = 0; very close points without the underflow of the squares
GPU_FUNC REAL P2PDistance (REAL dx, REAL dy, REAL dz)
{
  REAL r2 = dx*dx+dy*dy+dz*dz;
  if (r2 >= REAL(1e-30)) return sqrt(r2);
  REAL m = fabs(dx);
  if (fabs(dy) > m) m = fabs(dy);
  if (fabs(dz) > m) m = fabs(dz);
  if (m == 0) return 0;
  dx /= m; dy /= m; dz /= m;
  return m*sqrt(dx*dx+dy*dy+dz*dz);
}

// (p_hi - x_hi) + (p_lo - x_lo): the hi difference is exact for close points, fast math must not reassociate it with
// the lo parts
GPU_FUNC REAL P2PDiff (REAL p, REAL q, REAL x, REAL xlo)
{
#ifdef __METAL_VERSION__
#pragma clang fp reassociate(off)
#endif
  return (p-x) + (q-xlo);
}

#define P2P_INV4PI REAL(0.07957747154594767)

// the charge kernel of the host's SingularMLExpansion::Node::Evaluate, branch 0: kappa < 1e-12, 1: kappa < 1e-8, 2: else
GPU_FUNC Complex<REAL> P2PCharge (int branch, REAL invr, REAL rho, KAPPA_T kappa)
{
  return branch == 0 ? Complex<REAL>(P2P_INV4PI*invr, 0)
    : branch == 1 ? P2P_INV4PI*(IKappa(kappa)+invr)
    : (P2P_INV4PI*invr)*P2PPhase(kappa, rho);
}

// the near field of the value targets, as the host's SingularMLExpansion::Node::Evaluate (its SIMD branches): entry
// t = GLOBAL_ID_X adds to its target etarget the sources of the leaves [elist, elist+elcount) of llist (first source,
// count), one partial sum per leaf. kind 1 charges, 2 dipoles along sdir, 3 both (ChargeDipoles), bit 4 dipole vectors
// (with charges for bit 1); positions as REAL hi and lo parts; a coincident source adds nothing
KERNEL(fmm_p2p, GLOBAL_IN(int,etarget), GLOBAL_IN(int,elist), GLOBAL_IN(int,elcount), GLOBAL_IN(int,llist),
       GLOBAL_IN(REAL,txyz), GLOBAL_IN(REAL,sxyz), GLOBAL_IN(REAL,sdir),
       GLOBAL_IN(Complex<REAL>,charge), GLOBAL_IN(Complex<REAL>,dipole), GLOBAL_IN(Complex<REAL>,dipvec),
       GLOBAL(Complex<REAL>,out),
       VALUE(KAPPA_T,kappa), VALUE(int,kind), VALUE(int,branch), VALUE(int,comps), VALUE(int,entries), VALUE(int,outs))
{
  int t = GLOBAL_ID_X, comp = GROUP_ID_Y;
  if (t < entries)
    {
      int j = etarget[t], lf = elist[t], lc = elcount[t];
      GLOBAL_PTR(const REAL) pt = txyz + 6*j;
      REAL p0 = pt[0], p1 = pt[1], p2 = pt[2], q0 = pt[3], q1 = pt[4], q2 = pt[5];
      Complex<REAL> acc(0, 0);
      for (int l = lf; l < lf+lc; l++)
        {
          int sf = llist[2*l], sc = llist[2*l+1];
          Complex<REAL> a4[4];
          for (int k = 0; k < 4; k++) a4[k] = Complex<REAL>(0, 0);
          for (int s = sf; s < sf+sc; s++)
            {
              GLOBAL_PTR(const REAL) x = sxyz + 6*s;
              REAL dx = P2PDiff(p0, q0, x[0], x[3]), dy = P2PDiff(p1, q1, x[1], x[4]), dz = P2PDiff(p2, q2, x[2], x[5]);
              REAL rho = P2PDistance (dx, dy, dz);
              if (rho == 0) continue;
              REAL invr = REAL(1)/rho;
              Complex<REAL> term;
              if (kind == 1)
                term = P2PCharge(branch, invr, rho, kappa)*charge[s*comps+comp];
              else if (kind & 4)
                {
                  GLOBAL_PTR(const Complex<REAL>) D = dipvec + 3*(s*comps+comp);
                  Complex<REAL> dd = (dx*D[0] + dy*D[1] + dz*D[2])*invr;
                  term = ((P2P_INV4PI*P2PPhase(kappa, rho))*(IKappa(kappa)*invr + (-invr*invr)))*dd;
                  if (kind & 1) term = term + P2PCharge(branch, invr, rho, kappa)*charge[s*comps+comp];
                }
              else
                {
                  Complex<REAL> e = P2PPhase(kappa, rho);
                  REAL dd = (dx*sdir[3*s]+dy*sdir[3*s+1]+dz*sdir[3*s+2])*invr;
                  if (kind == 2)
                    term = ((P2P_INV4PI*e)*(IKappa(kappa)*invr + (-invr*invr)))*dd*dipole[s*comps+comp];
                  else
                    {
                      Complex<REAL> kc = (P2P_INV4PI*invr)*e;
                      term = kc*charge[s*comps+comp] + ((invr*dd)*kc*(IKappa(kappa)*rho - 1))*dipole[s*comps+comp];
                    }
                }
              a4[(s-sf)&3] = a4[(s-sf)&3] + term;
            }
          acc = acc + ((a4[0]+a4[1]) + (a4[2]+a4[3]));
        }
      int i = (j*comps+comp)*outs;
      out[i] = out[i] + acc;
    }
}

// the gradient near field as the host's SingularMLExpansion::Node::EvaluateDeriv (its SIMD branch): charges only, at
// out + goff of the targets (entries as fmm_p2p)
KERNEL(fmm_p2p_grad, GLOBAL_IN(int,etarget), GLOBAL_IN(int,elist), GLOBAL_IN(int,elcount), GLOBAL_IN(int,llist),
       GLOBAL_IN(REAL,txyz), GLOBAL_IN(REAL,sxyz), GLOBAL_IN(Complex<REAL>,charge), GLOBAL(Complex<REAL>,out),
       VALUE(KAPPA_T,kappa), VALUE(int,comps), VALUE(int,entries), VALUE(int,outs), VALUE(int,goff))
{
  int t = GLOBAL_ID_X, comp = GROUP_ID_Y;
  if (t < entries)
    {
      int j = etarget[t], lf = elist[t], lc = elcount[t];
      GLOBAL_PTR(const REAL) pt = txyz + 6*j;
      REAL p0 = pt[0], p1 = pt[1], p2 = pt[2], q0 = pt[3], q1 = pt[4], q2 = pt[5];
      Complex<REAL> zero(0, 0), gx = zero, gy = zero, gz = zero;
      for (int l = lf; l < lf+lc; l++)
        {
          int sf = llist[2*l], sc = llist[2*l+1];
          Complex<REAL> ax[4], ay[4], az[4];
          for (int k = 0; k < 4; k++) ax[k] = ay[k] = az[k] = zero;
          for (int s = sf; s < sf+sc; s++)
            {
              GLOBAL_PTR(const REAL) x = sxyz + 6*s;
              REAL dx = P2PDiff(p0, q0, x[0], x[3]), dy = P2PDiff(p1, q1, x[1], x[4]), dz = P2PDiff(p2, q2, x[2], x[5]);
              REAL rho = P2PDistance (dx, dy, dz);
              if (rho == 0) continue;
              REAL invr = REAL(1)/rho;
              Complex<REAL> g = ((P2P_INV4PI*P2PPhase(kappa, rho))*(IKappa(kappa)*(invr*invr) + (-(invr*invr*invr))))
                                 *charge[s*comps+comp];
              int k = (s-sf)&3;
              ax[k] = ax[k] + g*dx; ay[k] = ay[k] + g*dy; az[k] = az[k] + g*dz;
            }
          gx = gx + ((ax[0]+ax[1]) + (ax[2]+ax[3]));
          gy = gy + ((ay[0]+ay[1]) + (ay[2]+ay[3]));
          gz = gz + ((az[0]+az[1]) + (az[2]+az[3]));
        }
      int i = (j*comps+comp)*outs + goff;
      out[i] = out[i] + gx; out[i+1] = out[i+1] + gy; out[i+2] = out[i+2] + gz;
    }
}
