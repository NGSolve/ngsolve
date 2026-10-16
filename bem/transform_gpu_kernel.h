/*********************************************************************/
/* File:   transform_gpu_kernel.h                                    */
/* Purpose: gather / scatter of a unit of translations (the RotateZ  */
/*          and z-flip steps of the batched Transform) on the common */
/*          GPU layer (ngstd/gpukernel.h)                            */
/*********************************************************************/

/*
  The first and the last step of ProcessDynamicBatchSS/SR/RR, per lane b, ii = n*(n+1)+m, rows = (order+1)^2:
    gather   work(ii,b)                     = f(b,n,m)  * src_pool[src_idx[b]*rows+ii]   RotateZFlip( phi_b, flip_b)
    scatter  tgt_pool[tgt_idx[b]*rows+ii] += f'(b,n,m) * work(ii,b)                      RotateZFlip(-phi_b, flip_b)
  as SphericalHarmonics::RotateZFlip: f = F e^{i m phi_b}, f' = F e^{-i m phi_b}, F = -1 if flip_b and n-m odd
  (FlipFactor), from phase[b*phase_stride+|m|] = e^{i |m| phi_b} (host double, rounded once); work is row-major
  with dim columns. REAL is substituted and TRANSFORM_TILE defined by the host.

  One group per tile of TRANSFORM_TILE rows (x) x TRANSFORM_TILE lanes (y), transposed through group memory (below
  every backend's minimum) so both sides are contiguous; the work-items stride over the tile, any group size works.
  Only multiplications and fmas on the values, nothing fast-math dependent. The scatter adds with atomics, so
  targets may repeat; the additions into a repeated one come in no fixed order (on the CPU device in lane order).
*/

#ifdef NGS_GPU_CPU
using std::fma, std::sqrt;
#endif

// a*b with fused multiply-adds: the same two roundings per component on every backend
GPU_FUNC Complex<REAL> TransformMul(Complex<REAL> a, Complex<REAL> b)
{
  return Complex<REAL>(fma(a.re, b.re, -(a.im*b.im)), fma(a.re, b.im, a.im*b.re));
}

// F * e^{i m phi} (sign = 1) or F * e^{-i m phi} (sign = -1) from the table of e^{i k phi}, k >= 0
GPU_FUNC Complex<REAL> TransformFactor(GLOBAL_PTR(const Complex<REAL>) phase_b, int flip, int n, int m, int sign)
{
  Complex<REAL> f = phase_b[m < 0 ? -m : m];
  if (sign*m < 0) f = Complex<REAL>(f.re, -f.im);
  if (flip != 0 && ((n-m) & 1)) f = Complex<REAL>(-f.re, -f.im);
  return f;
}

// the degree n of row ii = n*(n+1)+m: the float estimate corrected in integers
GPU_FUNC int TransformDegree(int ii)
{
  int n = int(sqrt(float(ii)));
  while (n > 0 && n*n > ii) n--;
  while ((n+1)*(n+1) <= ii) n++;
  return n;
}

KERNEL(transform_gather, GLOBAL_IN(Complex<REAL>,src_pool), GLOBAL_IN(int,src_idx),
       GLOBAL_IN(Complex<REAL>,phase), GLOBAL_IN(int,flip), GLOBAL(Complex<REAL>,work),
       VALUE(int,order), VALUE(int,dim), VALUE(int,phase_stride))
{
  int tid = LOCAL_ID_X, threads = GROUP_SIZE_X;
  int rows = (order+1)*(order+1);
  int r0 = GROUP_ID_X*TRANSFORM_TILE, l0 = GROUP_ID_Y*TRANSFORM_TILE;
  SHARED(REAL, tile_re, TRANSFORM_TILE*(TRANSFORM_TILE+1));
  SHARED(REAL, tile_im, TRANSFORM_TILE*(TRANSFORM_TILE+1));
  // consecutive work-items: consecutive rows of one lane's source expansion
  for (int t = tid; t < TRANSFORM_TILE*TRANSFORM_TILE; t += threads)
  {
    int lane = t/TRANSFORM_TILE, r = t%TRANSFORM_TILE;
    int ii = r0+r, b = l0+lane;
    if (ii >= rows || b >= dim) continue;
    int n = TransformDegree(ii), m = ii-n*n-n;
    Complex<REAL> f = TransformFactor(phase+b*phase_stride, flip[b], n, m, 1);
    Complex<REAL> v = TransformMul(f, src_pool[src_idx[b]*rows+ii]);
    tile_re[lane*(TRANSFORM_TILE+1)+r] = v.re;
    tile_im[lane*(TRANSFORM_TILE+1)+r] = v.im;
  }
  BARRIER();
  // consecutive work-items: consecutive lanes of one row of work
  for (int t = tid; t < TRANSFORM_TILE*TRANSFORM_TILE; t += threads)
  {
    int r = t/TRANSFORM_TILE, lane = t%TRANSFORM_TILE;
    int ii = r0+r, b = l0+lane;
    if (ii >= rows || b >= dim) continue;
    work[ii*dim+b] = Complex<REAL>(tile_re[lane*(TRANSFORM_TILE+1)+r], tile_im[lane*(TRANSFORM_TILE+1)+r]);
  }
}

KERNEL(transform_scatter, GLOBAL_IN(Complex<REAL>,work), GLOBAL_IN(int,tgt_idx),
       GLOBAL_IN(Complex<REAL>,phase), GLOBAL_IN(int,flip), GLOBAL_ATOMIC_COMPLEX(REAL,tgt_pool),
       VALUE(int,order), VALUE(int,dim), VALUE(int,phase_stride))
{
  int tid = LOCAL_ID_X, threads = GROUP_SIZE_X;
  int rows = (order+1)*(order+1);
  int r0 = GROUP_ID_X*TRANSFORM_TILE, l0 = GROUP_ID_Y*TRANSFORM_TILE;
  SHARED(REAL, tile_re, TRANSFORM_TILE*(TRANSFORM_TILE+1));
  SHARED(REAL, tile_im, TRANSFORM_TILE*(TRANSFORM_TILE+1));
  // consecutive work-items: consecutive lanes of one row of work
  for (int t = tid; t < TRANSFORM_TILE*TRANSFORM_TILE; t += threads)
  {
    int r = t/TRANSFORM_TILE, lane = t%TRANSFORM_TILE;
    int ii = r0+r, b = l0+lane;
    if (ii >= rows || b >= dim) continue;
    Complex<REAL> w = work[ii*dim+b];
    tile_re[lane*(TRANSFORM_TILE+1)+r] = w.re;
    tile_im[lane*(TRANSFORM_TILE+1)+r] = w.im;
  }
  BARRIER();
  // consecutive work-items: consecutive rows of one lane's target expansion
  for (int t = tid; t < TRANSFORM_TILE*TRANSFORM_TILE; t += threads)
  {
    int lane = t/TRANSFORM_TILE, r = t%TRANSFORM_TILE;
    int ii = r0+r, b = l0+lane;
    if (ii >= rows || b >= dim) continue;
    int n = TransformDegree(ii), m = ii-n*n-n;
    Complex<REAL> f = TransformFactor(phase+b*phase_stride, flip[b], n, m, -1);
    Complex<REAL> v = TransformMul(f, Complex<REAL>(tile_re[lane*(TRANSFORM_TILE+1)+r], tile_im[lane*(TRANSFORM_TILE+1)+r]));
    int tt = tgt_idx[b]*rows+ii;
    ATOMIC_ADD_COMPLEX(tgt_pool, tt, v);
  }
}

// pool[begin+i] = 0 for i < count: zeroes pool slots inside a recorded program, where no host transfer is allowed
KERNEL(fmm_zero, GLOBAL(Complex<REAL>,pool), VALUE(int,begin), VALUE(int,count))
{
  int stride = NUM_GROUPS_X*GROUP_SIZE_X;
  for (int i = GLOBAL_ID_X; i < count; i += stride)
    pool[begin+i] = Complex<REAL>(REAL(0));
}
