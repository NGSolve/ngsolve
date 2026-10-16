/*********************************************************************/
/* File:   rotatey_gpu_kernel.h                                      */
/* Purpose: RotateY of a unit of spherical harmonics expansions      */
/*          on the common GPU layer (ngstd/gpukernel.h)              */
/*********************************************************************/

/*
  Port of RotateYImpl (mptools.cpp): Gumerov & Duraiswami, "Recursive Computation of Spherical Harmonic
  Rotation Coefficients of Large Degree", Excursions in Harmonic Analysis Vol. 3, p. 130.
  REAL is substituted and ROTATEY_ROWS defined by the host.

  legendre[n*(order+2)+m]  normalized P_n^m(cos theta), 0 <= m <= n <= order+1, host double, odd m negated
                           for theta < 0 (the transpose of the CPU normalized_leg_func(m,n))
  tables                   sqrt(k), 1/sqrt(k) (0 at k=0, never read), 1/sqrt((2n+1)(2n+3)) at offsets 0,
                           table_stride, 2*table_stride (>= 2*order+4), host double cast once: no sqrt and
                           no division here
  scratch                  per degree of a chunk: the (n+1) x (2n+1) trafo block, then the 2n+1 values
                           Dmn(mp,n); work_stride = (order+2)*(2*order+1) REALs per degree
  source/target            Complex<REAL>, (order+1)^2 rows, row-major; row n*(n+1)+m holds Coef(n,m)

  One launch: G rotations of one order, rot_scal[3g..3g+2] = (s, half_minus, half_plus); the apply runs over lane
  blocks[3i..3i+2] = (group, first lane, lanes), lanes are columns of source/target (row stride dimx).
*/

// One group constructs the trafo block of one degree n = first_degree+GROUP_ID_X
// of rotation GROUP_ID_Y; the grid is (degrees of the chunk) x G.
// The host processes degrees in chunks and reuses the scratch buffer.
KERNEL(rotatey_construct_g, GLOBAL_IN(REAL,legendre_all), GLOBAL_IN(REAL,tables),
       GLOBAL_IN(REAL,rot_scal), GLOBAL(REAL,scratch_all), VALUE(int,order),
       VALUE(int,first_degree), VALUE(int,work_stride), VALUE(int,table_stride),
       VALUE(int,nslots))
{
  // the recurrences cancel: fast math must not reassociate them
#ifdef __METAL_VERSION__
#pragma clang fp reassociate(off)
#endif
  int tid = LOCAL_ID_X, threads = GROUP_SIZE_X;
  int leg_width = order+2;
  int grp = GROUP_ID_Y;
  GLOBAL_PTR(const REAL) legendre = legendre_all+grp*(leg_width*leg_width);
  REAL s = rot_scal[3*grp], half_minus = rot_scal[3*grp+1], half_plus = rot_scal[3*grp+2];
  GLOBAL_PTR(const REAL) sqrt_int = tables;
  GLOBAL_PTR(const REAL) inv_sqrt_int = tables+table_stride;
  GLOBAL_PTR(const REAL) inv_sqrt_2np1_2np3 = tables+2*table_stride;
  GLOBAL_PTR(REAL) trafo = scratch_all+(grp*nslots+GROUP_ID_X)*work_stride;
  GLOBAL_PTR(REAL) d = trafo+(order+1)*(2*order+1);

  int n = first_degree+GROUP_ID_X;
  int width = 2*n+1;
  REAL invnorm = inv_sqrt_int[2*n+1];              // 1/sqrt(2n+1)
  REAL invnext = inv_sqrt_int[2*n+3];              // 1/sqrt(2n+3)
  REAL invn = inv_sqrt_2np1_2np3[n];               // 1/sqrt((2n+1)(2n+3))
  // 1/Bmn(0,n+1) = sqrt((2n+1)(2n+3)) / sqrt(n(n+1))
  REAL inv_b0 = sqrt_int[2*n+1]*sqrt_int[2*n+3]*inv_sqrt_int[n]*inv_sqrt_int[n+1];
  // Step 2: H(0,m) from P_n^m, Step 3: H(1,m) from P_{n+1}^m
  for (int m = tid; m <= n; m += threads)
  {
    REAL parity = m%2 ? REAL(-1) : REAL(1);
    trafo[m*width+n] = parity*invnorm*legendre[n*leg_width+m];
    if (m > 0)
    {
      REAL next = -parity*invnext*legendre[(n+1)*leg_width+m+1];
      REAL prev = -parity*invnext*legendre[(n+1)*leg_width+m-1];
      REAL mid = parity*invnext*legendre[(n+1)*leg_width+m];
      REAL b_neg = -sqrt_int[n+m+1]*sqrt_int[n+m+2]*invn;   // Bmn(-m-1,n+1)
      REAL b_pos = sqrt_int[n-m+1]*sqrt_int[n-m+2]*invn;    // Bmn(m-1,n+1)
      REAL a = sqrt_int[n+1+m]*sqrt_int[n+1-m]*invn;        // Amn(m,n)
      trafo[m*width+n+1] = inv_b0*(b_neg*half_minus*next-b_pos*half_plus*prev-a*s*mid);
    }
  }
  // Dmn(m,n) = sgn(m)/2 * sqrt((n-m)(n+m+1))
  for (int j = tid; j < width; j += threads)
  {
    int m = j-n;
    d[j] = (m >= 0 ? REAL(0.5) : REAL(-0.5))*sqrt_int[n-m]*sqrt_int[n+m+1];
  }
  DEVICE_BARRIER();

  // Step 4: diamond recurrence, sequential in mp, parallel in m
  for (int mp = 1; mp < n; mp++)
  {
    REAL invd = REAL(2)*inv_sqrt_int[n-mp]*inv_sqrt_int[n+mp+1];   // 1/Dmn(mp,n)
    for (int m = mp+1+tid; m <= n; m += threads)
    {
      REAL value = d[n+mp-1]*trafo[m*width+n+mp-1]
                   -d[n+m-1]*trafo[(m-1)*width+n+mp];
      if (m < n) value += d[n+m]*trafo[(m+1)*width+n+mp];
      trafo[m*width+n+mp+1] = invd*value;
    }
    DEVICE_BARRIER();
  }
  // Step 5: diamond recurrence, negative mp
  for (int mp = 0; mp > -n; mp--)
  {
    REAL invd = REAL(-2)*inv_sqrt_int[n-mp+1]*inv_sqrt_int[n+mp];   // 1/Dmn(mp-1,n)
    for (int m = -mp+1+tid; m <= n; m += threads)
    {
      REAL value = d[n+mp]*trafo[m*width+n+mp+1]
                   +d[n+m-1]*trafo[(m-1)*width+n+mp];
      if (m < n) value -= d[n+m]*trafo[(m+1)*width+n+mp];
      trafo[m*width+n+mp-1] = invd*value;
    }
    DEVICE_BARRIER();
  }

  // Step 6: the other half by symmetry; all reads are from the half
  // already constructed, so no work-item reads another one's new value.
  for (int i = tid; i < (n+1)*width; i += threads)
  {
    int m = i/width, mp = i%width-n;
    if (mp > m) trafo[i] = trafo[mp*width+n+m];
    if (-mp > m) trafo[i] = trafo[-mp*width+n-m];
  }
}

// Applies the blocks of one chunk: one group per (part of a lane block, degree n = first_degree+GROUP_ID_Y),
// group x = block*parts + part. An item is ROTATEY_ROWS rows of one column: they share each source load and
// keep one sum each, in the m order of a single row, so every output has the bits of a row per work-item.
KERNEL(rotatey_apply_b, GLOBAL_IN(REAL,scratch_all), GLOBAL_IN(int,blocks),
       GLOBAL_IN(Complex<REAL>,source), GLOBAL(Complex<REAL>,target),
       VALUE(int,first_degree), VALUE(int,work_stride), VALUE(int,nslots), VALUE(int,dimx), VALUE(int,parts))
{
  int tid = LOCAL_ID_X, threads = GROUP_SIZE_X;
  int n = first_degree+GROUP_ID_Y;
  int width = 2*n+1;
  int rblocks = (width+ROTATEY_ROWS-1)/ROTATEY_ROWS;
  int bi = GROUP_ID_X/parts, part = GROUP_ID_X%parts;
  int grp = blocks[3*bi], col0 = blocks[3*bi+1], dim = blocks[3*bi+2];
  GLOBAL_PTR(const REAL) trafo = scratch_all+(grp*nslots+GROUP_ID_Y)*work_stride;
  for (int i = tid+part*threads; i < rblocks*dim; i += threads*parts)
  {
    if (n == 1 && i < dim) target[col0+i] = source[col0+i];   // degree 0 is unchanged
    int row0 = (i/dim)*ROTATEY_ROWS, column = col0+i%dim;
    // rows past the block's end repeat its last row: computed, never stored
    int row[ROTATEY_ROWS];
    for (int k = 0; k < ROTATEY_ROWS; k++)
      row[k] = row0+k < width ? row0+k : width-1;
    Complex<REAL> positive[ROTATEY_ROWS], negative[ROTATEY_ROWS];
    for (int k = 0; k < ROTATEY_ROWS; k++)
    {
      positive[k] = Complex<REAL>(REAL(0));
      negative[k] = Complex<REAL>(REAL(0));
    }
    for (int m = 0; m <= n; m++)
    {
      Complex<REAL> s = source[(n*n+n+m)*dimx+column];
      for (int k = 0; k < ROTATEY_ROWS; k++)
        positive[k] += trafo[m*width+row[k]]*s;
    }
    for (int m = 1; m <= n; m++)
    {
      Complex<REAL> s = source[(n*n+n-m)*dimx+column];
      for (int k = 0; k < ROTATEY_ROWS; k++)
        negative[k] += trafo[m*width+2*n-row[k]]*s;
    }
    for (int k = 0; k < ROTATEY_ROWS && row0+k < width; k++)
      target[(n*n+row0+k)*dimx+column] = ((row0+k-n)%2 ? REAL(-1) : REAL(1))*(positive[k]+negative[k]);
  }
}
