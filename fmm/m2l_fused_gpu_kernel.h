/*********************************************************************/
/* File:   m2l_fused_gpu_kernel.h                                    */
/* Purpose: the S->R translation of a unit (gather, RotateY,         */
/*          ShiftZ, RotateY back, scatter) fused in group memory     */
/*********************************************************************/

/*
  The five steps of transform_gather, rotatey_apply_b, shiftz_apply_b, rotatey_apply_b and transform_scatter for the
  lanes of a lane block: one group per block (group, first lane, lanes <= FUSED_LANES), whose expansions stay in group
  memory, FUSED_ROWS >= (max(os,ot)+1)^2 rows each. The rotation blocks come folded from rot_fold, all ShiftZ layers
  at once from shiftz_boundaries_g and shiftz_layers_g (complex layers, as always for S->R). REAL is substituted by
  the host, FUSED_LANES and FUSED_ROWS are defined by it; compiled after transform_gpu_kernel.h.

  Folded RotateY of degree n: with the trafo block T of rotatey_construct_g, S(m,r) = (T(m,r)+T(m,2n-r))/2 and
  D(m,r) = (T(m,r)-T(m,2n-r))/2, r <= n, and the source coefficients a_m, b_m of m and -m as s_m = a_m+b_m,
  d_m = a_m-b_m (s_0 = d_0 = a_0), the rows r and 2n-r are sgn(r) (sum_m S(m,r) s_m +- sum_m D(m,r) d_m): half the
  multiplications of rotatey_apply_b. The ShiftZ sums keep one accumulator (os+1 <= 28 terms), so the results are
  not bitwise equal to the unfused path. The scatter adds with atomics as transform_scatter (on the CPU device an
  entry receives its lanes in lane order).
*/

#ifdef NGS_GPU_CPU
using std::ldexp;
#endif

// rows (RotateY) and degrees (ShiftZ) per work-item, sharing each load; the register arrays are only indexed by
// compile-time loops (lanes past cnt are computed on zeros and never stored)
#define FUSED_ROTATE_ROWS 2
#define FUSED_SHIFT_DEGREES 4

// degree n's folded block, S then D, each (n+1) x (n+1) [m*(n+1)+r], at this offset in a rotation's fold tables:
// sum_{1 <= k < n} 2(k+1)^2; FoldOffset(order+1) reals per rotation
GPU_FUNC int FoldOffset(int n)
{
  return n*(n+1)*(2*n+1)/3-2;
}

// the folded blocks of the chunk's degrees n = first_degree+GROUP_ID_X of rotation GROUP_ID_Y, from the scratch that
// rotatey_construct_g wrote (same chunk arguments)
KERNEL(rot_fold, GLOBAL_IN(REAL,scratch_all), GLOBAL(REAL,fold_all), VALUE(int,order), VALUE(int,first_degree),
       VALUE(int,work_stride), VALUE(int,nslots))
{
  int n = first_degree+GROUP_ID_X, grp = GROUP_ID_Y, width = 2*n+1, w = n+1;
  GLOBAL_PTR(const REAL) T = scratch_all+(grp*nslots+GROUP_ID_X)*work_stride;
  GLOBAL_PTR(REAL) F = fold_all+grp*FoldOffset(order+1)+FoldOffset(n);
  for (int i = LOCAL_ID_X; i < w*w; i += GROUP_SIZE_X)
  {
    int m = i/w, r = i%w;
    REAL a = T[m*width+r], b = T[m*width+2*n-r];
    F[i] = REAL(0.5)*(a+b);
    F[w*w+i] = REAL(0.5)*(a-b);
  }
}

// the work-items of degree n in FusedRotateY and of order sm in FusedShiftZ
GPU_FUNC int RotateItems(int n) { return (n+FUSED_ROTATE_ROWS)/FUSED_ROTATE_ROWS; }
GPU_FUNC int ShiftItems(int ot, int sm) { return (ot-(sm < 0 ? -sm : sm)+FUSED_SHIFT_DEGREES)/FUSED_SHIFT_DEGREES; }

// out = RotateY(in) for the cnt lanes; in is overwritten by the sums and differences
GPU_FUNC void FusedRotateY(LOCAL_PTR(REAL) in, LOCAL_PTR(REAL) out, GLOBAL_PTR(const REAL) fold, int order, int cnt,
                           int tid, int nth)
{
  int rows = (order+1)*(order+1);
  for (int t = tid; t < cnt*rows; t += nth)
  {
    int l = t/rows, ii = t-l*rows;
    int n = TransformDegree(ii), m = ii-n*n-n;
    if (m <= 0) continue;
    int ia = 2*(l*FUSED_ROWS+ii), ib = 2*(l*FUSED_ROWS+ii-2*m);
    REAL ar = in[ia], ai = in[ia+1], br = in[ib], bi = in[ib+1];
    in[ia] = ar+br; in[ia+1] = ai+bi; in[ib] = ar-br; in[ib+1] = ai-bi;
  }
  BARRIER();
  // items: (degree n, FUSED_ROTATE_ROWS of its rows r <= n)
  int items = 0;
  for (int n = 0; n <= order; n++) items += RotateItems(n);
  for (int it = tid; it < items; it += nth)
  {
    int n = 0, acc = 0;
    while (acc+RotateItems(n) <= it) acc += RotateItems(n++);
    int r0 = (it-acc)*FUSED_ROTATE_ROWS, w = n+1, base = n*n+n;
    if (n == 0)
    {
      for (int l = 0; l < FUSED_LANES; l++)
        if (l < cnt) { out[2*l*FUSED_ROWS] = in[2*l*FUSED_ROWS]; out[2*l*FUSED_ROWS+1] = in[2*l*FUSED_ROWS+1]; }
      continue;
    }
    GLOBAL_PTR(const REAL) F = fold+FoldOffset(n);
    // rows past the degree's end repeat its last row: computed, never stored
    int rr[FUSED_ROTATE_ROWS];
    for (int k = 0; k < FUSED_ROTATE_ROWS; k++) rr[k] = r0+k < w ? r0+k : w-1;
    REAL pr[FUSED_ROTATE_ROWS][FUSED_LANES], pi[FUSED_ROTATE_ROWS][FUSED_LANES];
    REAL qr[FUSED_ROTATE_ROWS][FUSED_LANES], qi[FUSED_ROTATE_ROWS][FUSED_LANES];
    for (int k = 0; k < FUSED_ROTATE_ROWS; k++)
      for (int l = 0; l < FUSED_LANES; l++) { pr[k][l] = 0; pi[k][l] = 0; qr[k][l] = 0; qi[k][l] = 0; }
    for (int m = 0; m <= n; m++)
    {
      REAL sr[FUSED_LANES], si[FUSED_LANES], dr[FUSED_LANES], di[FUSED_LANES];
      for (int l = 0; l < FUSED_LANES; l++)
      {
        bool on = l < cnt;
        sr[l] = on ? in[2*(l*FUSED_ROWS+base+m)] : REAL(0);
        si[l] = on ? in[2*(l*FUSED_ROWS+base+m)+1] : REAL(0);
        dr[l] = m == 0 ? sr[l] : (on ? in[2*(l*FUSED_ROWS+base-m)] : REAL(0));
        di[l] = m == 0 ? si[l] : (on ? in[2*(l*FUSED_ROWS+base-m)+1] : REAL(0));
      }
      for (int k = 0; k < FUSED_ROTATE_ROWS; k++)
      {
        REAL S = F[m*w+rr[k]], D = F[w*w+m*w+rr[k]];
        for (int l = 0; l < FUSED_LANES; l++)
        {
          pr[k][l] += S*sr[l]; pi[k][l] += S*si[l];
          qr[k][l] += D*dr[l]; qi[k][l] += D*di[l];
        }
      }
    }
    for (int k = 0; k < FUSED_ROTATE_ROWS; k++)
    {
      int r = r0+k;
      REAL sgn = ((r-n) % 2) ? REAL(-1) : REAL(1);
      for (int l = 0; l < FUSED_LANES; l++)
      {
        if (r >= w || l >= cnt) continue;
        out[2*(l*FUSED_ROWS+n*n+r)] = sgn*(pr[k][l]+qr[k][l]);
        out[2*(l*FUSED_ROWS+n*n+r)+1] = sgn*(pi[k][l]+qi[k][l]);
        if (r < n)
        {
          out[2*(l*FUSED_ROWS+n*n+2*n-r)] = sgn*(pr[k][l]-qr[k][l]);
          out[2*(l*FUSED_ROWS+n*n+2*n-r)+1] = sgn*(pi[k][l]-qi[k][l]);
        }
      }
    }
  }
}

// powscale(k)*v with the mantissa/exponent split, as ShiftZPow in shiftz_gpu_kernel.h
GPU_FUNC REAL FusedPow(REAL v, REAL mf, int e)
{
  return e > 0 ? ldexp(v, e-1)*(REAL(2)*mf) : ldexp(v*mf, e);
}

// out = ShiftZ(in) for the cnt lanes, as shiftz_apply_b (S->R):
// out(n, sm) = fout(n) sum_{l=|sm|}^{os} T_|sm|(n,l) fin(l) in(l, sm)
GPU_FUNC void FusedShiftZ(LOCAL_PTR(REAL) in, LOCAL_PTR(REAL) out, GLOBAL_PTR(const Complex<REAL>) layer0,
                          GLOBAL_PTR(const Complex<REAL>) layers, GLOBAL_PTR(const REAL) pow_mant,
                          GLOBAL_PTR(const int) pow_exp, int os, int ot, int cnt, int tid, int nth)
{
  int H = os+ot+1, W = (os > ot ? os : ot)+1;
  // items: (sm, FUSED_SHIFT_DEGREES of the target degrees n >= |sm|), sm = -ot..ot
  int items = 0;
  for (int sm = -ot; sm <= ot; sm++) items += ShiftItems(ot, sm);
  for (int it = tid; it < items; it += nth)
  {
    int sm = -ot, acc = 0;
    while (acc+ShiftItems(ot, sm) <= it) acc += ShiftItems(ot, sm++);
    int m = sm < 0 ? -sm : sm;
    int n0 = m+(it-acc)*FUSED_SHIFT_DEGREES;
    if (m > os)
    {
      // no layer writes |m| > os
      for (int k = 0; k < FUSED_SHIFT_DEGREES; k++)
        for (int l = 0; l < FUSED_LANES; l++)
        {
          if (n0+k > ot || l >= cnt) continue;
          int ii = (n0+k)*(n0+k+1)+sm;
          out[2*(l*FUSED_ROWS+ii)] = 0; out[2*(l*FUSED_ROWS+ii)+1] = 0;
        }
      continue;
    }
    GLOBAL_PTR(const Complex<REAL>) T = m == 0 ? layer0 : layers+(m-1)*(H*W);
    int nn[FUSED_SHIFT_DEGREES];
    for (int k = 0; k < FUSED_SHIFT_DEGREES; k++) nn[k] = n0+k <= ot ? n0+k : ot;
    REAL ar[FUSED_SHIFT_DEGREES][FUSED_LANES], ai[FUSED_SHIFT_DEGREES][FUSED_LANES];
    for (int k = 0; k < FUSED_SHIFT_DEGREES; k++)
      for (int l = 0; l < FUSED_LANES; l++) { ar[k][l] = 0; ai[k][l] = 0; }
    for (int j = m; j <= os; j++)
    {
      REAL f = pow_mant[W+j];
      int row = j*(j+1)+sm;
      REAL sr[FUSED_LANES], si[FUSED_LANES];
      for (int l = 0; l < FUSED_LANES; l++)
      {
        sr[l] = l < cnt ? f*in[2*(l*FUSED_ROWS+row)] : REAL(0);
        si[l] = l < cnt ? f*in[2*(l*FUSED_ROWS+row)+1] : REAL(0);
      }
      for (int k = 0; k < FUSED_SHIFT_DEGREES; k++)
      {
        Complex<REAL> t = T[nn[k]*W+j];
        for (int l = 0; l < FUSED_LANES; l++)
        {
          ar[k][l] += t.re*sr[l]-t.im*si[l];
          ai[k][l] += t.re*si[l]+t.im*sr[l];
        }
      }
    }
    for (int k = 0; k < FUSED_SHIFT_DEGREES; k++)
    {
      int n = nn[k], ii = n*(n+1)+sm;
      REAL fm = pow_mant[W+os+1+n];
      int fe = pow_exp[W+n];
      for (int l = 0; l < FUSED_LANES; l++)
      {
        if (n0+k > ot || l >= cnt) continue;
        out[2*(l*FUSED_ROWS+ii)] = FusedPow(ar[k][l], fm, fe);
        out[2*(l*FUSED_ROWS+ii)+1] = FusedPow(ai[k][l], fm, fe);
      }
    }
  }
}

// one group per lane block; the tables of translation group g: fold_all + g*FoldOffset(os+1) (RotateY(theta_g)),
// fold_back_all + g*FoldOffset(ot+1), layer 0 at layer0_all + g*H*W, layer m >= 1 at layers_all + (g*nslots+m-1)*H*W,
// pow_mant/pow_exp as shiftz_apply_b
KERNEL(m2l_fused, GLOBAL_IN(int,blocks), GLOBAL_IN(Complex<REAL>,src_pool), GLOBAL_IN(int,src_idx),
       GLOBAL_IN(Complex<REAL>,phase), GLOBAL_IN(int,flip), GLOBAL_IN(REAL,fold_all), GLOBAL_IN(REAL,fold_back_all),
       GLOBAL_IN(Complex<REAL>,layer0_all), GLOBAL_IN(Complex<REAL>,layers_all),
       GLOBAL_IN(REAL,pow_mant_all), GLOBAL_IN(int,pow_exp_all), GLOBAL_IN(int,tgt_idx),
       GLOBAL_ATOMIC_COMPLEX(REAL,tgt_pool),
       VALUE(int,os), VALUE(int,ot), VALUE(int,phase_stride), VALUE(int,nslots))
{
  int tid = LOCAL_ID_X, nth = GROUP_SIZE_X;
  int grp = blocks[3*GROUP_ID_X], lane0 = blocks[3*GROUP_ID_X+1], cnt = blocks[3*GROUP_ID_X+2];
  int rs = (os+1)*(os+1), rt = (ot+1)*(ot+1), H = os+ot+1, W = (os > ot ? os : ot)+1;
  SHARED(REAL, A, 2*FUSED_LANES*FUSED_ROWS);
  SHARED(REAL, B, 2*FUSED_LANES*FUSED_ROWS);
  for (int t = tid; t < cnt*rs; t += nth)
  {
    int l = t/rs, ii = t-l*rs, b = lane0+l;
    int n = TransformDegree(ii), m = ii-n*n-n;
    Complex<REAL> f = TransformFactor(phase+b*phase_stride, flip[b], n, m, 1);
    Complex<REAL> v = TransformMul(f, src_pool[src_idx[b]*rs+ii]);
    A[2*(l*FUSED_ROWS+ii)] = v.re; A[2*(l*FUSED_ROWS+ii)+1] = v.im;
  }
  BARRIER();
  FusedRotateY(A, B, fold_all+grp*FoldOffset(os+1), os, cnt, tid, nth);
  BARRIER();
  FusedShiftZ(B, A, layer0_all+grp*(H*W), layers_all+grp*nslots*(H*W), pow_mant_all+grp*(W+os+ot+2),
              pow_exp_all+grp*(W+ot+1), os, ot, cnt, tid, nth);
  BARRIER();
  FusedRotateY(A, B, fold_back_all+grp*FoldOffset(ot+1), ot, cnt, tid, nth);
  BARRIER();
  for (int t = tid; t < cnt*rt; t += nth)
  {
    int l = t/rt, ii = t-l*rt, b = lane0+l;
    int n = TransformDegree(ii), m = ii-n*n-n;
    Complex<REAL> f = TransformFactor(phase+b*phase_stride, flip[b], n, m, -1);
    Complex<REAL> v = TransformMul(f, Complex<REAL>(B[2*(l*FUSED_ROWS+ii)], B[2*(l*FUSED_ROWS+ii)+1]));
    int tt = tgt_idx[b]*rt+ii;
    ATOMIC_ADD_COMPLEX(tgt_pool, tt, v);
  }
}
