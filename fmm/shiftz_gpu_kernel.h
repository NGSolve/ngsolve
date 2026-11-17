/*********************************************************************/
/* File:   shiftz_gpu_kernel.h                                       */
/* Purpose: ShiftZ of a unit of translation groups on the common     */
/*          GPU layer (ngstd/gpukernel.h), layers built in parallel  */
/*          from chained boundaries                                  */
/*********************************************************************/

/*
  Port of SphericalExpansion::ShiftZ (mptools.cpp; Gumerov & Duraiswami, "fast, exact, stable": (185) seed,
  (187) first column, (182) walls, (181) bulk). REAL and TRAFO are substituted by the host: TRAFO = REAL
  (S->S and R->R with real kappa) or Complex<REAL> (S->R or complex kappa).

  Layer m is the H x W matrix trafo(l,n), H = os+ot+1, W = max(os,ot)+1, row-major T[l*W+n]. Its boundary B_m
  is its first column trafo(l,m), l < H, plus the four wall diagonals (Lof(n,d), n), n < W, d = 0..3, with
  Lof = n+d (S->R, upward walls) or os+ot-n-d (S->S, R->R, downward walls); in the boundary buffer at
  base m*(H+4W): col[l] = base+l, band(n,d) = base+H+4n+d. B_m depends only on B_{m-1}, position n of layer m
  only on positions <= n of layer m-1; given B_m the bulk of layer m is independent of every other layer.

  ctab          Amn, 1/Amn, Bmn(-m,n), Bmn(m-1,n), 1/Bmn(-m,n) for m <= min(os,ot), n <= os+ot+1, correctly
                rounded on the host (GetCoef2DTables): five blocks of cblock reals, entry m*cstride+n
  seed          trafo(l,0) of layer 0, l < H, with the z<0 odd-l flip and sqrt(2l+1), host double
  pow_mant/exp  powscale(k) = mant*2^exp, k < W, then the apply's fin(l) and fout(n) (layout: DeviceShiftZ)
  in/out        Complex<REAL>, (os+1)^2 / (ot+1)^2 rows, row-major, row n*(n+1)+m = Coef(n,m)

  One launch serves a unit of G translation groups of one os, ot, variant and layer type, each with its tables
  and workspace at uniform per-group strides; its lanes are columns of in/out (row stride dimx), group after
  group. The apply runs over lane blocks blocks[3i..3i+2] = (group, first lane, lanes).
*/

#ifdef NGS_GPU_CPU
using std::ldexp, std::fma;
#endif

// a*b + c fused per real component: the three-term recurrences are fma chains, so every backend
// evaluates them with the same three roundings (the plain expression has five)
GPU_FUNC REAL ShiftZFma(REAL a, REAL b, REAL c)
{
  return fma(a, b, c);
}
GPU_FUNC Complex<REAL> ShiftZFma(REAL a, Complex<REAL> b, Complex<REAL> c)
{
  return Complex<REAL>(fma(a, b.re, c.re), fma(a, b.im, c.im));
}

// powscale(k)*v with the mantissa/exponent split, per real component
GPU_FUNC REAL ShiftZPow(REAL v, REAL mf, int e)
{
  return e > 0 ? ldexp(v, e-1)*(REAL(2)*mf) : ldexp(v*mf, e);
}
GPU_FUNC Complex<REAL> ShiftZPow(Complex<REAL> v, REAL mf, int e)
{
  return Complex<REAL>(ShiftZPow(v.re, mf, e), ShiftZPow(v.im, mf, e));
}

#define COEF_PARAMS GLOBAL_PTR(const REAL) ctab, int cstride, int cblock
#define COEF_ARGS   ctab, cstride, cblock

GPU_FUNC REAL CoefAmn(COEF_PARAMS, int m, int n)     { return ctab[m*cstride+n]; }
GPU_FUNC REAL CoefInvAmn(COEF_PARAMS, int m, int n)  { return ctab[cblock+m*cstride+n]; }
GPU_FUNC REAL CoefBneg(COEF_PARAMS, int m, int n)    { return ctab[2*cblock+m*cstride+n]; }
GPU_FUNC REAL CoefBpos(COEF_PARAMS, int m, int n)    { return ctab[3*cblock+m*cstride+n]; }
GPU_FUNC REAL CoefInvBneg(COEF_PARAMS, int m, int n) { return ctab[4*cblock+m*cstride+n]; }

GPU_FUNC int ShiftZLof(int variant, int os, int ot, int n, int d)
{
  return variant == 1 ? n+d : os+ot-n-d;
}

// trafo(l,n) of layer m in the boundary buffer; the column n == m is the first column
GPU_FUNC int ShiftZBIndex(int H, int W, int variant, int os, int ot, int m, int l, int n)
{
  return m*(H+4*W) + (n == m ? l : H+4*n+(variant == 1 ? l-n : os+ot-l-n));
}
#define BND(m, l, n) bnd[ShiftZBIndex(H, W, variant, os, ot, m, l, n)]

/*
  One group per translation group g = GROUP_ID_X (grid G). Layer 0 fully into layer0 (column recursion,
  parallel over l; then the symmetric fill), B_0 into the boundary buffer, then the wavefront over steps t:
  layer m is at position n = t-m. Position m is the first column (parallel over l) and the two wall starts;
  position n > m is iteration n of both walls, one thread per layer. Wall entries outside the triangle
  l+n <= os+ot are not computed: the CPU code computes them from stale data and never reads them.
*/
KERNEL(shiftz_boundaries_g, GLOBAL(TRAFO,layer0_all), GLOBAL(TRAFO,bnd_all),
       GLOBAL_IN(TRAFO,seed_all), GLOBAL_IN(REAL,ctab),
       GLOBAL_IN(REAL,pow_mant_all), GLOBAL_IN(int,pow_exp_all), GLOBAL_IN(REAL,scal_all),
       VALUE(int,os), VALUE(int,ot), VALUE(int,variant), VALUE(int,cstride), VALUE(int,cblock))
{
  int tid = LOCAL_ID_X, threads = GROUP_SIZE_X;
  int H = os+ot+1, W = (os > ot ? os : ot)+1, M = os < ot ? os : ot;
  int grp = GROUP_ID_X;
  GLOBAL_PTR(TRAFO) layer0 = layer0_all+grp*(H*W);
  GLOBAL_PTR(TRAFO) bnd = bnd_all+grp*((M+1)*(H+4*W));
  GLOBAL_PTR(const TRAFO) seed = seed_all+grp*H;
  GLOBAL_PTR(const REAL) pow_mant = pow_mant_all+grp*(W+os+ot+2);
  GLOBAL_PTR(const int) pow_exp = pow_exp_all+grp*(W+ot+1);
  GLOBAL_PTR(const REAL) sc = scal_all+grp*9;
  REAL scale = sc[0], tscale = sc[1], inv_scale = sc[2], inv_tscale = sc[3], scale_tscale = sc[4],
       scale_inv_tscale = sc[5], scale2 = sc[6], inv_tscale2 = sc[7], inv_scale_inv_tscale = sc[8];
  bool is_sr = variant == 1;

  // ---- layer 0
  for (int i = tid; i < H*W; i += threads) layer0[i] = TRAFO(REAL(0));
  for (int i = tid; i < (M+1)*(H+4*W); i += threads) bnd[i] = TRAFO(REAL(0));
  DEVICE_BARRIER();
  for (int l = tid; l < H; l += threads) layer0[l*W] = seed[l];
  DEVICE_BARRIER();
  REAL inv_a00 = CoefInvAmn(COEF_ARGS, 0, 0);
  for (int l = 1+tid; l < os+ot; l += threads)
    layer0[l*W+1] = ShiftZFma(-(CoefAmn(COEF_ARGS, 0, l)*inv_a00*scale_tscale), layer0[(l+1)*W],
                              (CoefAmn(COEF_ARGS, 0, l-1)*inv_a00*scale_inv_tscale)*layer0[(l-1)*W]);
  if (tid == 0) layer0[1] = ShiftZPow(layer0[W], pow_mant[1], pow_exp[1]);
  DEVICE_BARRIER();
  for (int n = 1; n < W-1; n++)
  {
    REAL inv_a0n = CoefInvAmn(COEF_ARGS, 0, n);
    REAL a0nm1 = CoefAmn(COEF_ARGS, 0, n-1);
    REAL c3 = a0nm1*inv_a0n*scale2;
    for (int l = n+tid; l < os+ot-n; l += threads)
      layer0[l*W+n+1] = ShiftZFma(-(CoefAmn(COEF_ARGS, 0, l)*inv_a0n*scale_tscale), layer0[(l+1)*W+n],
                          ShiftZFma(CoefAmn(COEF_ARGS, 0, l-1)*inv_a0n*scale_inv_tscale, layer0[(l-1)*W+n],
                            c3*layer0[l*W+n-1]));
    if (tid == 0) layer0[n+1] = ShiftZPow(layer0[(n+1)*W], pow_mant[n+1], pow_exp[n+1]);
    DEVICE_BARRIER();
  }
  // symmetric fill: upper (n,l), n < l < W, from the lower (l,n): disjoint
  for (int i = tid; i < W*W; i += threads)
  {
    int n = i/W, l = i%W;
    if (l > n) layer0[n*W+l] = ShiftZPow(layer0[l*W+n], pow_mant[l-n], pow_exp[l-n]);
  }
  DEVICE_BARRIER();
  // B_0
  for (int l = tid; l < H; l += threads) bnd[l] = layer0[l*W];
  for (int i = tid; i < 4*W; i += threads)
  {
    int n = i/4, d = i%4, l = ShiftZLof(variant, os, ot, n, d);
    if (n >= 1 && l >= 0 && l < H) bnd[H+i] = layer0[l*W+n];
  }
  DEVICE_BARRIER();

  // ---- boundaries B_1..B_M by the wavefront (walls only)
  int tmax = M+W;
  for (int t = 2; t <= tmax; t++)
  {
    // position m of layer m = t/2: the first column (187), parallel over l
    if (t%2 == 0 && t/2 <= M)
    {
      int m = t/2;
      REAL inv_bmm = CoefInvBneg(COEF_ARGS, m, m);
      for (int l = m+tid; l <= H-1-m; l += threads)
        BND(m, l, m) = ShiftZFma(CoefBneg(COEF_ARGS, m, l)*inv_bmm*scale_inv_tscale, BND(m-1, l-1, m-1),
                                 -(CoefBpos(COEF_ARGS, m, l+1)*inv_bmm*scale_tscale)*BND(m-1, l+1, m-1));
    }
    DEVICE_BARRIER();
    // walls: one thread per layer of this step
    for (int m = 1+tid; m <= M; m += threads)
    {
      int n = t-m;
      if (n < m || n > W-2) continue;
      for (int wall = 1; wall <= 2; wall++)
      {
        if (is_sr)
        {
          if (!(n < W-1 && 2*n+wall+1 <= os+ot)) continue;
          int l = n+wall;
          REAL inv_b = CoefInvBneg(COEF_ARGS, m, n+1);
          TRAFO v = ShiftZFma(-(CoefBpos(COEF_ARGS, m, l+1)*inv_b*scale_tscale), BND(m-1, l+1, n),
                              (CoefBneg(COEF_ARGS, m, l)*inv_b*scale_inv_tscale)*BND(m-1, l-1, n));
          if (n == m) { BND(m, l, n+1) = v; continue; }   // wall start (182)
          REAL inv_al = CoefInvAmn(COEF_ARGS, m, l);
          BND(m, l, n+1) = ShiftZFma(CoefBpos(COEF_ARGS, m, n)*inv_b*scale2, BND(m, l, n-1), v);
          BND(m, l+1, n) = ShiftZFma(CoefAmn(COEF_ARGS, m, l-1)*inv_al*inv_tscale2, BND(m, l-1, n),
                             ShiftZFma(CoefAmn(COEF_ARGS, m, n-1)*inv_al*scale_inv_tscale, BND(m, l, n-1),
                               -(CoefAmn(COEF_ARGS, m, n)*inv_al*inv_scale_inv_tscale)*BND(m, l, n+1)));
        }
        else
        {
          if (!(n < W-wall && 2*n+wall+1 <= os+ot)) continue;
          int l = os+ot-n-wall;
          REAL c = scale*CoefInvBneg(COEF_ARGS, m, n+1);
          TRAFO v = ShiftZFma(-(c*CoefBpos(COEF_ARGS, m, l+1)*tscale), BND(m-1, l+1, n),
                              (c*CoefBneg(COEF_ARGS, m, l)*inv_tscale)*BND(m-1, l-1, n));
          if (n == m) { BND(m, l, n+1) = v; continue; }   // wall start (182)
          REAL t_inv_alm1 = tscale*CoefInvAmn(COEF_ARGS, m, l-1);
          BND(m, l, n+1) = ShiftZFma(c*CoefBpos(COEF_ARGS, m, n)*scale, BND(m, l, n-1), v);
          BND(m, l-1, n) = ShiftZFma(t_inv_alm1*CoefAmn(COEF_ARGS, m, l)*tscale, BND(m, l+1, n),
                             ShiftZFma(-(t_inv_alm1*CoefAmn(COEF_ARGS, m, n-1)*scale), BND(m, l, n-1),
                               (t_inv_alm1*CoefAmn(COEF_ARGS, m, n)*inv_scale)*BND(m, l, n+1)));
        }
      }
    }
    DEVICE_BARRIER();
  }
}

/*
  One group per (layer m = m0+GROUP_ID_X >= 1, translation group GROUP_ID_Y), built in slot GROUP_ID_X of its group
  (layer 0 is in layer0): zero, load B_m, bulk (181) row by row (sequential over l, parallel over n) with the ctab
  rows Amn(m,.) and 1/Amn(m,.), symmetric fill with powscale. shiftz_apply_b applies it.
*/
KERNEL(shiftz_layers_g, GLOBAL_IN(TRAFO,bnd_all),
       GLOBAL(TRAFO,layers_all), GLOBAL_IN(REAL,ctab),
       GLOBAL_IN(REAL,pow_mant_all), GLOBAL_IN(int,pow_exp_all), GLOBAL_IN(REAL,scal_all),
       VALUE(int,os), VALUE(int,ot), VALUE(int,variant), VALUE(int,cstride), VALUE(int,cblock),
       VALUE(int,m0), VALUE(int,nslots))
{
  int tid = LOCAL_ID_X, threads = GROUP_SIZE_X;
  int H = os+ot+1, W = (os > ot ? os : ot)+1, M = os < ot ? os : ot;
  int m = m0+GROUP_ID_X;
  bool is_sr = variant == 1, is_rr = variant == 2;
  int grp = GROUP_ID_Y;
  GLOBAL_PTR(const TRAFO) bnd = bnd_all+grp*((M+1)*(H+4*W));
  GLOBAL_PTR(TRAFO) L = layers_all+grp*(nslots*(H*W))+GROUP_ID_X*(H*W);
  GLOBAL_PTR(const REAL) amn = ctab+m*cstride;
  GLOBAL_PTR(const REAL) inv_amn = amn+cblock;
  GLOBAL_PTR(const REAL) pow_mant = pow_mant_all+grp*(W+os+ot+2);
  GLOBAL_PTR(const int) pow_exp = pow_exp_all+grp*(W+ot+1);
  GLOBAL_PTR(const REAL) sc = scal_all+grp*9;
  REAL scale = sc[0], tscale = sc[1], inv_scale = sc[2], scale_inv_tscale = sc[5], inv_tscale2 = sc[7],
       inv_scale_inv_tscale = sc[8];

  for (int i = tid; i < H*W; i += threads) L[i] = TRAFO(REAL(0));
  DEVICE_BARRIER();
  // load B_m
  GLOBAL_PTR(const TRAFO) base = bnd+m*(H+4*W);
  for (int l = tid; l < H; l += threads) L[l*W+m] = base[l];
  for (int i = tid; i < 4*W; i += threads)
  {
    int n = i/4, d = i%4, l = ShiftZLof(variant, os, ot, n, d);
    if (n > m && l >= 0 && l < H) L[l*W+n] = base[H+i];
  }
  DEVICE_BARRIER();
  // bulk (181)
  if (is_sr)
  {
    for (int l = m; l < H-1; l++)
    {
      int nmax = H-1-l < l-1 ? H-1-l : l-1;
      REAL inv_al = inv_amn[l];
      REAL c1 = amn[l-1]*inv_al*inv_tscale2;
      for (int n = m+1+tid; n < nmax; n += threads)
        L[(l+1)*W+n] = ShiftZFma(c1, L[(l-1)*W+n],
                         ShiftZFma(amn[n-1]*inv_al*scale_inv_tscale, L[l*W+n-1],
                           -(amn[n]*inv_al*inv_scale_inv_tscale)*L[l*W+n+1]));
      DEVICE_BARRIER();
    }
  }
  else
  {
    for (int l = os+ot-m; l >= m; l--)
    {
      int nmax = os+ot-l < l ? os+ot-l : l;
      REAL t_inv_alm1 = tscale*inv_amn[l-1];
      REAL c1 = t_inv_alm1*amn[l]*tscale;
      for (int n = m+1+tid; n < nmax; n += threads)
        L[(l-1)*W+n] = ShiftZFma(c1, L[(l+1)*W+n],
                         ShiftZFma(-(t_inv_alm1*amn[n-1]*scale), L[l*W+n-1],
                           (t_inv_alm1*amn[n]*inv_scale)*L[l*W+n+1]));
      DEVICE_BARRIER();
    }
  }
  // symmetric fill (m <= n < l <= symmetry_order)
  int so = is_rr ? ot : os;
  for (int i = tid; i < W*W; i += threads)
  {
    int n = i/W, l = i%W;
    if (n >= m && n < so && l > n && l <= so)
      L[n*W+l] = ShiftZPow(L[l*W+n], pow_mant[l-n], pow_exp[l-n]);
  }
  DEVICE_BARRIER();
}

/*
  Apply +m and -m: one group per (part of a lane block, layer m = m0+GROUP_ID_Y), group x = block*parts + part;
  each (n, lane) item of a block is computed by one work-item. Layer 0 is read from layer0, a layer m >= 1 from
  slot GROUP_ID_Y of its translation group. The groups of layer 0 also zero the rows of out with |m| > os,
  which no layer writes.
*/
KERNEL(shiftz_apply_b, GLOBAL_IN(TRAFO,layer0_all), GLOBAL_IN(TRAFO,layers_all),
       GLOBAL_IN(REAL,pow_mant_all), GLOBAL_IN(int,pow_exp_all), GLOBAL_IN(int,blocks),
       GLOBAL_IN(Complex<REAL>,source), GLOBAL(Complex<REAL>,target),
       VALUE(int,os), VALUE(int,ot), VALUE(int,variant), VALUE(int,m0), VALUE(int,nslots), VALUE(int,dimx),
       VALUE(int,parts))
{
  int tid = LOCAL_ID_X, threads = GROUP_SIZE_X;
  int H = os+ot+1, W = (os > ot ? os : ot)+1;
  int m = m0+GROUP_ID_Y;
  bool is_rr = variant == 2;
  int bi = GROUP_ID_X/parts, part = GROUP_ID_X%parts;
  int grp = blocks[3*bi], lane0 = blocks[3*bi+1], dim = blocks[3*bi+2];
  GLOBAL_PTR(const TRAFO) T = m == 0 ? layer0_all+grp*(H*W) : layers_all+(grp*nslots+GROUP_ID_Y)*(H*W);
  GLOBAL_PTR(const REAL) pow_mant = pow_mant_all+grp*(W+os+ot+2);
  GLOBAL_PTR(const int) pow_exp = pow_exp_all+grp*(W+ot+1);

  // out(n, sm) = fout(n) sum_{l=m}^{os} trafo(n,l) fin(l) in(l, sm), transposed for R->R; lane b is column lane0+b.
  // Eight interleaved partial sums: the S->R sum has up to os+1 terms of comparable size, and one sequential
  // accumulator loses about (os+1)*eps (the CPU float path's matvec accumulates in several lanes too).
#define SHIFTZ_TERM(L) ((is_rr ? T[(L)*W+n] : T[n*W+(L)])*(pow_mant[W+(L)]*source[((L)*((L)+1)+sm)*dimx+lane0+b]))
  int items = (ot-m+1)*dim;
  for (int i = tid+part*threads; i < items; i += threads*parts)
  {
    int n = m+i/dim, b = i%dim;
    for (int s = 0; s < (m == 0 ? 1 : 2); s++)
    {
      int sm = s == 0 ? m : -m;
      Complex<REAL> a0(REAL(0)), a1(REAL(0)), a2(REAL(0)), a3(REAL(0)), a4(REAL(0)), a5(REAL(0)), a6(REAL(0)), a7(REAL(0));
      int l = m;
      for (; l+7 <= os; l += 8)
      {
        a0 += SHIFTZ_TERM(l);   a1 += SHIFTZ_TERM(l+1); a2 += SHIFTZ_TERM(l+2); a3 += SHIFTZ_TERM(l+3);
        a4 += SHIFTZ_TERM(l+4); a5 += SHIFTZ_TERM(l+5); a6 += SHIFTZ_TERM(l+6); a7 += SHIFTZ_TERM(l+7);
      }
      for (; l <= os; l++) a0 += SHIFTZ_TERM(l);
      target[(n*(n+1)+sm)*dimx+lane0+b] = ShiftZPow(((a0+a1)+(a2+a3))+((a4+a5)+(a6+a7)), pow_mant[W+os+1+n], pow_exp[W+n]);
    }
    if (m == 0)
      for (int mm = os+1; mm <= n; mm++)
      {
        target[(n*(n+1)+mm)*dimx+lane0+b] = Complex<REAL>(REAL(0));
        target[(n*(n+1)-mm)*dimx+lane0+b] = Complex<REAL>(REAL(0));
      }
  }
}
