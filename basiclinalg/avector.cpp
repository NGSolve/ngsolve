#include <bla.hpp>
using namespace ngbla;


#if defined(__AVX2__)

// namespace ngbla
// {


  // n ... number of SIMDs
  INLINE void MyScal1x4 (int n, 
                         __m256d * a1,
                         __m256d * b1, __m256d * b2, __m256d * b3, __m256d * b4,
                         __m256d & s1)
  {
    __m256d sum11 = _mm256_setzero_pd();
    __m256d sum12 = _mm256_setzero_pd();
    __m256d sum13 = _mm256_setzero_pd();
    __m256d sum14 = _mm256_setzero_pd();

    for (int i = 0; i < n; i++)
      {
        sum11 += a1[i] * b1[i];
        sum12 += a1[i] * b2[i];
        sum13 += a1[i] * b3[i];
        sum14 += a1[i] * b4[i];
      }
  
    s1 = HAdd (sum11, sum12, sum13, sum14);
  }

  // n ... number of SIMDs  
  INLINE void MyScal2x4 (int n, 
                         __m256d * a1, __m256d * a2,
                         __m256d * b1, __m256d * b2, __m256d * b3, __m256d * b4,
                         __m256d & s1, __m256d & s2)
  {
    __m256d sum11 = _mm256_setzero_pd();
    __m256d sum21 = _mm256_setzero_pd();

    __m256d sum12 = _mm256_setzero_pd();
    __m256d sum22 = _mm256_setzero_pd();

    __m256d sum13 = _mm256_setzero_pd();
    __m256d sum23 = _mm256_setzero_pd();

    __m256d sum14 = _mm256_setzero_pd();
    __m256d sum24 = _mm256_setzero_pd();

    for (int i = 0; i < n; i++)
      {
        sum11 += a1[i] * b1[i];
        sum21 += a2[i] * b1[i];
        sum12 += a1[i] * b2[i];
        sum22 += a2[i] * b2[i];
        sum13 += a1[i] * b3[i];
        sum23 += a2[i] * b3[i];
        sum14 += a1[i] * b4[i];
        sum24 += a2[i] * b4[i];
      }
  
    s1 = HAdd (sum11, sum12, sum13, sum14);
    s2 = HAdd (sum21, sum22, sum23, sum24);
  }


  // n ... number of doubles, must be even
  INLINE
  void MyScal4x4 (int n,
                  __m128d * pa1, __m128d * pa2, __m128d * pa3, __m128d * pa4,
                  __m128d * pb1, __m128d * pb2, __m128d * pb3, __m128d * pb4,
                  __m256d & s1, __m256d & s2, __m256d & s3, __m256d & s4)
  {
    __m256d sum11 = _mm256_setzero_pd();
    __m256d sum21 = _mm256_setzero_pd();
    __m256d sum31 = _mm256_setzero_pd();
    __m256d sum41 = _mm256_setzero_pd();
    __m256d sum12 = _mm256_setzero_pd();
    __m256d sum22 = _mm256_setzero_pd();
    __m256d sum32 = _mm256_setzero_pd();
    __m256d sum42 = _mm256_setzero_pd();

    int n2 = n/2;
    for (int i = 0; i < n2; i++)
      {
        __m256d a1 = _mm256_broadcast_pd(pa1+i);
        __m256d a2 = _mm256_broadcast_pd(pa2+i);
        __m256d a3 = _mm256_broadcast_pd(pa3+i);
        __m256d a4 = _mm256_broadcast_pd(pa4+i);

        __m256d b1 = _mm256_broadcast_pd(pb1+i);
        __m256d b2 = _mm256_broadcast_pd(pb2+i);
        __m256d b3 = _mm256_broadcast_pd(pb3+i);
        __m256d b4 = _mm256_broadcast_pd(pb4+i);
        __m256d mb1 = _mm256_blend_pd(b1, b3, 12);
        __m256d mb2 = _mm256_blend_pd(b2, b4, 12);
      
        sum11 += a1 * mb1;
        sum21 += a2 * mb1;
        sum31 += a3 * mb1;
        sum41 += a4 * mb1;
        sum12 += a1 * mb2;
        sum22 += a2 * mb2;
        sum32 += a3 * mb2;
        sum42 += a4 * mb2;
      }
    s1 = _mm256_hadd_pd(sum11, sum12);
    s2 = _mm256_hadd_pd(sum21, sum22);
    s3 = _mm256_hadd_pd(sum31, sum32);
    s4 = _mm256_hadd_pd(sum41, sum42);
  }
        





  // C += A * Trans(B)

  void AddABt (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c)
  {
    int i = 0;
    // clear overhead
    if (a.Width() != 4*a.VWidth())
      {
        int r = 4*a.VWidth()-a.Width();
        __m256i mask = _mm256_cmpgt_epi64(_mm256_set1_epi64x(r),
                                          _mm256_set_epi64x(0,1,2,3));
        /*
          __m256i mask;
          switch (r)
          {
          case 1:
          mask = _mm256_set_epi64x(-1,0,0,0); break;
          case 2:
          mask = _mm256_set_epi64x(-1,-1,0,0); break;
          case 3:
          mask = _mm256_set_epi64x(-1,-1,-1,0); break;
          }
        */
        __m256d zero = _mm256_setzero_pd();
        for (int i = 0; i < a.Height(); i++)
          _mm256_maskstore_pd((double*)&a.Get(i, a.VWidth()-1), mask, zero);
        for (int i = 0; i < b.Height(); i++)
          _mm256_maskstore_pd((double*)&b.Get(i, b.VWidth()-1), mask, zero);
      }
  
    if (a.VWidth() <= 0) return;
  
    for ( ; i < c.Height()-1; i += 2)
      {
        int j = 0;
        for ( ; j < c.Width()-3; j += 4)
          {
            __m256d s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            s1 += _mm256_loadu_pd(&c(i,j));
            s2 += _mm256_loadu_pd(&c(i+1,j));
            _mm256_storeu_pd(&c(i,j), s1);
            _mm256_storeu_pd(&c(i+1,j), s2);
          }
        if (j < c.Width())
          {
            __m256d s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            for (int j2 = 0; j2 < c.Width()-j; j2++)
              {
                c(i,j+j2) += ((double*)(&s1))[j2];
                c(i+1,j+j2) += ((double*)(&s2))[j2];
              }
          }
      }

  
    if (i < c.Height())
      {
        int j = 0;
        for ( ; j < c.Width()-3; j += 4)
          {
            __m256d s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            s1 += _mm256_loadu_pd(&c(i,j));
            _mm256_storeu_pd(&c(i,j), s1);
          }
        if (j < c.Width())
          {
            __m256d s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            for (int j2 = 0; j2 < c.Width()-j; j2++)
              c(i,j+j2) += ((double*)(&s1))[j2];
          }
      }
  }


  void AddABtSymV1 (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c)
  {
    int i = 0;
    // clear overhead
    if (a.Width() != 4*a.VWidth())
      {
        int r = 4*a.VWidth()-a.Width();
        __m256i mask = _mm256_cmpgt_epi64(_mm256_set1_epi64x(r),
                                          _mm256_set_epi64x(0,1,2,3));

        __m256d zero = _mm256_setzero_pd();
        for (int i = 0; i < a.Height(); i++)
          _mm256_maskstore_pd((double*)&a.Get(i, a.VWidth()-1), mask, zero);
        for (int i = 0; i < b.Height(); i++)
          _mm256_maskstore_pd((double*)&b.Get(i, b.VWidth()-1), mask, zero);
      }
  
    if (a.VWidth() <= 0) return;
  
    for ( ; i < c.Height()-1; i += 2)
      {
        int j = 0;
        for ( ; j < i; j += 4)
          {
            __m256d s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            s1 += _mm256_loadu_pd(&c(i,j));
            s2 += _mm256_loadu_pd(&c(i+1,j));
            _mm256_storeu_pd(&c(i,j), s1);
            _mm256_storeu_pd(&c(i+1,j), s2);
          }
        if (j <= i)
          {
            __m256d s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            __m128d s1l = _mm256_extractf128_pd(s1, 0);
            __m128d s2l = _mm256_extractf128_pd(s2, 0);
            s1l += _mm_loadu_pd(&c(i,j));
            s2l += _mm_loadu_pd(&c(i+1,j));
            _mm_storeu_pd(&c(i,j), s1l);
            _mm_storeu_pd(&c(i+1,j), s2l);
          }
      }
  
    if (i < c.Height())
      {
        int j = 0;
        for ( ; j < c.Width()-3; j += 4)
          {
            __m256d s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            s1 += _mm256_loadu_pd(&c(i,j));
            _mm256_storeu_pd(&c(i,j), s1);
          }
        if (j < c.Width())
          {
            __m256d s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            for (int j2 = 0; j2 < c.Width()-j; j2++)
              c(i,j+j2) += ((double*)(&s1))[j2];
          }
      }
  }


  /*
    INLINE void AddABtSym (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c)
    {
    // clear overhead
    if (a.Width() != 4*a.VWidth())
    {
    int r = 4*a.VWidth()-a.Width();
    __m256i mask = _mm256_cmpgt_epi64(_mm256_set1_epi64x(r),
    _mm256_set_epi64x(0,1,2,3));

    __m256d zero = _mm256_setzero_pd();
    for (int i = 0; i < a.Height(); i++)
    _mm256_maskstore_pd((double*)&a.Get(i, a.VWidth()-1), mask, zero);
    for (int i = 0; i < b.Height(); i++)
    _mm256_maskstore_pd((double*)&b.Get(i, b.VWidth()-1), mask, zero);
    }
  
    if (a.VWidth() <= 0) return;
  
    int j = 0;
    for ( ; j < c.Width()-3; j += 4)
    {
    int i = j;
    double * pc = &c(i,j);
    for ( ; i < c.Height()-1; i += 2)
    {
    __m256d s1, s2;
    MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
    &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);

    // s1 += _mm256_loadu_pd(&c(i,j));
    // s2 += _mm256_loadu_pd(&c(i+1,j));
    // _mm256_storeu_pd(&c(i,j), s1);
    // _mm256_storeu_pd(&c(i+1,j), s2);


    s1 += _mm256_loadu_pd(pc);
    _mm256_storeu_pd(pc, s1);
    pc += c.Dist();
    s2 += _mm256_loadu_pd(pc);
    _mm256_storeu_pd(pc, s2);
    pc += c.Dist();          
    }
    if (i < c.Height())
    {
    __m256d s1;
    MyScal1x4 (a.VWidth(), &a.Get(i,0),
    &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1);

    s1 += _mm256_loadu_pd(pc);
    _mm256_storeu_pd(pc, s1);
    }
    }

    for ( ; j < c.Width(); j++)
    for (int i = j; i < c.Height(); i++)
    c(i,j) += InnerProduct(a.Row(i), b.Row(j));
    }
  */

  void AddABtSym (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c)
  {
    // clear overhead
    if (a.Width() != 4*a.VWidth())
      {
        int r = 4*a.VWidth()-a.Width();
        __m256i mask = _mm256_cmpgt_epi64(_mm256_set1_epi64x(r),
                                          _mm256_set_epi64x(0,1,2,3));

        __m256d zero = _mm256_setzero_pd();
        for (int i = 0; i < a.Height(); i++)
          _mm256_maskstore_pd((double*)&a.Get(i, a.VWidth()-1), mask, zero);
        for (int i = 0; i < b.Height(); i++)
          _mm256_maskstore_pd((double*)&b.Get(i, b.VWidth()-1), mask, zero);
      }
  
    if (a.VWidth() <= 0) return;
  
    int j = 0;
    for ( ; j < c.Width()-3; j += 4)
      {
        int i = j;
        double * pc = &c(i,j);
        for ( ; i < c.Height()-3; i += 4)
          {
            __m256d s1, s2, s3, s4;
            MyScal4x4 (4*a.VWidth(),
                       (__m128d*)&a.Get(i,0), (__m128d*)&a.Get(i+1,0), (__m128d*)&a.Get(i+2,0), (__m128d*)&a.Get(i+3,0),
                       (__m128d*)&b.Get(j,0), (__m128d*)&b.Get(j+1,0), (__m128d*)&b.Get(j+2,0), (__m128d*)&b.Get(j+3,0), s1, s2, s3, s4);

            s1 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s1);
            pc += c.Dist();
            s2 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s2);
            pc += c.Dist();          
            s3 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s3);
            pc += c.Dist();
            s4 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s4);
            pc += c.Dist();          
          }

        if (i < c.Height()-1)
          {
            __m256d s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);

            s1 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s1);
            pc += c.Dist();
            s2 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s2);
            pc += c.Dist();
            i += 2;
          }

      if (i < c.Height())
        {
          __m256d s1;
          MyScal1x4 (a.VWidth(), &a.Get(i,0),
                     &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1);

          s1 += _mm256_loadu_pd(pc);
          _mm256_storeu_pd(pc, s1);
        }
      }

    for ( ; j < c.Width(); j++)
      for (int i = j; i < c.Height(); i++)
        c(i,j) += InnerProduct(a.Row(i), b.Row(j));
  }








  // n ... number of doubles
  INLINE
  void MyScal4x4 (int n,
                  double * pa1, double * pa2, double * pa3, double * pa4,
                  Complex * _pb1, Complex * _pb2, Complex * _pb3, Complex * _pb4,
                  __m256d & sum11, __m256d & sum21, __m256d & sum31, __m256d & sum41,
                  __m256d & sum12, __m256d & sum22, __m256d & sum32, __m256d & sum42)
  {
    __m128d * pb1 = reinterpret_cast<__m128d*> (_pb1);
    __m128d * pb2 = reinterpret_cast<__m128d*> (_pb2);
    __m128d * pb3 = reinterpret_cast<__m128d*> (_pb3);
    __m128d * pb4 = reinterpret_cast<__m128d*> (_pb4);


    sum11 = _mm256_setzero_pd();
    sum21 = _mm256_setzero_pd();
    sum31 = _mm256_setzero_pd();
    sum41 = _mm256_setzero_pd();
    sum12 = _mm256_setzero_pd();
    sum22 = _mm256_setzero_pd();
    sum32 = _mm256_setzero_pd();
    sum42 = _mm256_setzero_pd();

    for (int i = 0; i < n; i++)
      {
        __m256d a1 = _mm256_broadcast_sd(pa1+i);
        __m256d a2 = _mm256_broadcast_sd(pa2+i);
        __m256d a3 = _mm256_broadcast_sd(pa3+i);
        __m256d a4 = _mm256_broadcast_sd(pa4+i);

        __m256d b1 = _mm256_broadcast_pd(pb1+i);
        __m256d b2 = _mm256_broadcast_pd(pb2+i);
        __m256d b3 = _mm256_broadcast_pd(pb3+i);
        __m256d b4 = _mm256_broadcast_pd(pb4+i);
        __m256d mb1 = _mm256_blend_pd(b1, b2, 12);
        __m256d mb2 = _mm256_blend_pd(b3, b4, 12);
      
        sum11 += a1 * mb1;
        sum21 += a2 * mb1;
        sum31 += a3 * mb1;
        sum41 += a4 * mb1;
        sum12 += a1 * mb2;
        sum22 += a2 * mb2;
        sum32 += a3 * mb2;
        sum42 += a4 * mb2;
      }
    /*
    s1 = _mm256_hadd_pd(sum11, sum12);
    s2 = _mm256_hadd_pd(sum21, sum22);
    s3 = _mm256_hadd_pd(sum31, sum32);
    s4 = _mm256_hadd_pd(sum41, sum42);
    */
  }

  // n ... number of doubles
  INLINE
  void MyScal2x4 (int n,
                  double * pa1, double * pa2,
                  Complex * _pb1, Complex * _pb2, Complex * _pb3, Complex * _pb4,
                  __m256d & sum11, __m256d & sum21,
                  __m256d & sum12, __m256d & sum22)
  {
    __m128d * pb1 = reinterpret_cast<__m128d*> (_pb1);
    __m128d * pb2 = reinterpret_cast<__m128d*> (_pb2);
    __m128d * pb3 = reinterpret_cast<__m128d*> (_pb3);
    __m128d * pb4 = reinterpret_cast<__m128d*> (_pb4);
    
    sum11 = _mm256_setzero_pd();
    sum21 = _mm256_setzero_pd();
    sum12 = _mm256_setzero_pd();
    sum22 = _mm256_setzero_pd();
    
    for (int i = 0; i < n; i++)
      {
        __m256d a1 = _mm256_broadcast_sd(pa1+i);
        __m256d a2 = _mm256_broadcast_sd(pa2+i);

        __m256d b1 = _mm256_broadcast_pd(pb1+i);
        __m256d b2 = _mm256_broadcast_pd(pb2+i);
        __m256d b3 = _mm256_broadcast_pd(pb3+i);
        __m256d b4 = _mm256_broadcast_pd(pb4+i);
        __m256d mb1 = _mm256_blend_pd(b1, b2, 12);
        __m256d mb2 = _mm256_blend_pd(b3, b4, 12);
      
        sum11 += a1 * mb1;
        sum21 += a2 * mb1;
        sum12 += a1 * mb2;
        sum22 += a2 * mb2;
      }
  }

  // n ... number of doubles
  INLINE
  void MyScal1x4 (int n,
                  double * pa1,
                  Complex * _pb1, Complex * _pb2, Complex * _pb3, Complex * _pb4,
                  __m256d & sum11,
                  __m256d & sum12)
  {
    __m128d * pb1 = reinterpret_cast<__m128d*> (_pb1);
    __m128d * pb2 = reinterpret_cast<__m128d*> (_pb2);
    __m128d * pb3 = reinterpret_cast<__m128d*> (_pb3);
    __m128d * pb4 = reinterpret_cast<__m128d*> (_pb4);
    
    sum11 = _mm256_setzero_pd();
    sum12 = _mm256_setzero_pd();
    
    for (int i = 0; i < n; i++)
      {
        __m256d a1 = _mm256_broadcast_sd(pa1+i);

        __m256d b1 = _mm256_broadcast_pd(pb1+i);
        __m256d b2 = _mm256_broadcast_pd(pb2+i);
        __m256d b3 = _mm256_broadcast_pd(pb3+i);
        __m256d b4 = _mm256_broadcast_pd(pb4+i);
        __m256d mb1 = _mm256_blend_pd(b1, b2, 12);
        __m256d mb2 = _mm256_blend_pd(b3, b4, 12);
      
        sum11 += a1 * mb1;
        sum12 += a1 * mb2;
      }
  }


void AddABtSym (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c)
{
  static Timer t("AddABtSym, double-complex"); RegionTimer reg(t);
  // c += a * Trans(b);
  // return;
  /*
  cout << "a = " << endl << a << endl;
  cout << "b = " << endl << b << endl;
  cout << "c = " << endl << c << endl;
  cout << "c1 = " << endl << c+a*Trans(b) << endl;
  */
  
  if (a.Width() <= 0) return;
  
  int j = 0;
  for ( ; j < c.Width()-3; j += 4)
    {
      int i = j;
      Complex * pc = &c(i,j);
      for ( ; i < c.Height()-3; i += 4)
        {
          __m256d s11, s21, s31, s41;
          __m256d s12, s22, s32, s42;
          MyScal4x4 (a.Width(),
                     &a(i,0), &a(i+1,0), &a(i+2,0), &a(i+3,0),
                     &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0),
                     s11, s21, s31, s41,
                     s12, s22, s32, s42);

          double * pcd = reinterpret_cast<double*> (pc);
          s11 += _mm256_loadu_pd(pcd);
          _mm256_storeu_pd(pcd, s11);
          s12 += _mm256_loadu_pd(pcd+4);
          _mm256_storeu_pd(pcd+4, s12);
          
          pc += c.Dist();
          pcd = reinterpret_cast<double*> (pc);          
          s21 += _mm256_loadu_pd(pcd);
          _mm256_storeu_pd(pcd, s21);
          s22 += _mm256_loadu_pd(pcd+4);
          _mm256_storeu_pd(pcd+4, s22);

          pc += c.Dist();
          pcd = reinterpret_cast<double*> (pc);                    
          s31 += _mm256_loadu_pd(pcd);
          _mm256_storeu_pd(pcd, s31);
          s32 += _mm256_loadu_pd(pcd+4);
          _mm256_storeu_pd(pcd+4, s32);

          pc += c.Dist();
          pcd = reinterpret_cast<double*> (pc);                              
          s41 += _mm256_loadu_pd(pcd);
          _mm256_storeu_pd(pcd, s41);
          s42 += _mm256_loadu_pd(pcd+4);
          _mm256_storeu_pd(pcd+4, s42);
          pc += c.Dist();          
        }

      if (i < c.Height()-1)
        {
          __m256d s11, s21;
          __m256d s12, s22;
          MyScal2x4 (a.Width(), &a(i,0), &a(i+1,0),
                     &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0),
                     s11, s21, s12, s22);

          double * pcd = reinterpret_cast<double*> (pc);          
          s11 += _mm256_loadu_pd(pcd);
          _mm256_storeu_pd(pcd, s11);
          s12 += _mm256_loadu_pd(pcd+4);
          _mm256_storeu_pd(pcd+4, s12);
          
          pc += c.Dist();
          pcd = reinterpret_cast<double*> (pc);          
          s21 += _mm256_loadu_pd(pcd);
          _mm256_storeu_pd(pcd, s21);
          s22 += _mm256_loadu_pd(pcd+4);
          _mm256_storeu_pd(pcd+4, s22);
          
          pc += c.Dist();
          i += 2;
        }

      if (i < c.Height())
        {
          __m256d s11, s12;
          MyScal1x4 (a.Width(), &a(i,0),
                     &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0), s11, s12);
          
          double * pcd = reinterpret_cast<double*> (pc);                    
          s11 += _mm256_loadu_pd(pcd);
          _mm256_storeu_pd(pcd, s11);
          s12 += _mm256_loadu_pd(pcd+4);
          _mm256_storeu_pd(pcd+4, s12);
        }
      }

  for ( ; j < c.Width(); j++)
    for (int i = j; i < c.Height(); i++)
      c(i,j) += InnerProduct(a.Row(i), b.Row(j));

  /*
  cout << "c2 = " << endl << c << endl;
    if (c.Height() > 4)
      exit(1);
  */
}

void AddABtSym (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c)
{
  cout << "AddABTsym complex - complex not implemented" << endl;
}


void AddABt (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c)
{
  // c += a * Trans(b);
  // return;
  
    int i = 0;

    if (a.Width() <= 0) return;
  
    for ( ; i < c.Height()-1; i += 2)
      {
        int j = 0;
        for ( ; j < c.Width()-3; j += 4)
          {
            __m256d s11, s21, s12, s22;
            MyScal2x4 (a.Width(), &a(i,0), &a(i+1,0),
                       &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0),
                       s11, s21, s12, s22);
            s11 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j)));
            s21 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+1,j)));
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j)), s11);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+1,j)), s21);
            s12 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j+2)));
            s22 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+1,j+2)));
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j+2)), s12);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+1,j+2)), s22);
          }
        for ( ; j < c.Width()-1; j += 2)
          {
            __m256d s11, s21, s12, s22;
            MyScal2x4 (a.Width(), &a(i,0), &a(i+1,0),
                       &b(j,0), &b(j+1,0), &b(j+1,0), &b(j+1,0),
                       s11, s21, s12, s22);
            s11 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j)));
            s21 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+1,j)));
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j)), s11);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+1,j)), s21);
          }
        for (  ; j < c.Width(); j++)
          {
            c(i,j) += InnerProduct(a.Row(i), b.Row(j));
            c(i+1,j) += InnerProduct(a.Row(i+1), b.Row(j));
          }
        /*
        if (j < c.Width())
          {
            __m256d s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            for (int j2 = 0; j2 < c.Width()-j; j2++)
              {
                c(i,j+j2) += ((double*)(&s1))[j2];
                c(i+1,j+j2) += ((double*)(&s2))[j2];
              }
          }
        */
      }

  
    if (i < c.Height())
      {
        int j = 0;
        for ( ; j < c.Width()-3; j += 4)
          {
            __m256d s11, s12;
            MyScal1x4 (a.Width(), &a(i,0),
                       &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0), s11, s12);
            s11 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j)));
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j)), s11);
            s12 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j+2)));
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j+2)), s12);
          }
        for ( ; j < c.Width(); j++)
          {
            c(i,j) += InnerProduct(a.Row(i), b.Row(j));            
          }
      }
}

void AddABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c)
{
  cout << "addabt complex-complex not implemented" << endl;
}



// }

#endif
