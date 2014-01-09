#if defined(__SSE3__) || defined(__AVX__)
#include <immintrin.h>
// #include <pmmintrin.h>
#endif

#include <fem.hpp>


/*
  
Computes C += A B^t
with (M) and slice (M2) of A and B is given at compile-time
compute only lower left triangle of C

This is the most expensive operation in element matrix calculation.

*/


#ifdef WIN32
#define __restrict__ __restrict
#endif

namespace ngfem {
  
  typedef std::complex<double> Complex;


  /*
  inline ostream & operator<< (ostream & ost, __m256d d)
  {
    double * p = (double*)&d;
    for (int i = 0; i < 4; i++)
      ost << p[i] << " ";
    return ost;
  }

  double HSum (__m256d d)
  {
    double * p = (double*)&d;
    double sum = 0;
    for (int i = 0; i < 4; i++)
      sum += p[i];
    return sum;
  }
  */


#ifdef __AVX__
  INLINE __m256d HAdd (__m256d v1, __m256d v2, __m256d v3, __m256d v4)
  {
    __m256d hsum1 = _mm256_hadd_pd (v1, v2);
    __m256d hsum2 = _mm256_hadd_pd (v3, v4);

    __m256d hsum = 
      _mm256_add_pd (_mm256_insertf128_pd (hsum1, 
					   _mm256_extractf128_pd (hsum2, 0), 1),
		     _mm256_insertf128_pd (hsum2, 
					   _mm256_extractf128_pd (hsum1, 1), 0));
    return hsum;
  }
#endif


  template <int M> NGS_DLL_HEADER
  void FastMat (int n, int M2, Complex * ba, Complex *  pb, Complex * pc);

  template <int M> NGS_DLL_HEADER
  void FastMat (int n, int M2, double * __restrict__ ba, double *  __restrict__ pb, double * __restrict__ pc);
  
  
  template <int M> 
  void FastMat (int n, 
                int M2,
                double * __restrict__ pa, 
		double * __restrict__ pb, 
		double * __restrict__ pc)
  {
    // static Timer timer (string("Fastmat, M = ")+ToString(M), 2); 
    // RegionTimer reg (timer);  timer.AddFlops (double(M)*n*n/2);

    /* 
    for (int i = 0; i < n; i++)
      for (int j = 0; j <= i; j++)
	{
	  double sum = pc[n*i+j];

	  double * lpa = pa + i * M;
	  double * lpb = pb + j * M;

	  for (int k = 0; k <= M-2; k+=2)
	    sum += lpa[k] * lpb[k] + lpa[k+1] * lpb[k+1];
	  if  (M % 2)
	    sum += lpa[M-1] * lpb[M-1];

	  pc[j+n*i] = sum;
          // pc[i+n*j] = sum;
	}
    */



    /*
    for (int i = 0; i < n-1; i+=2)
      for (int j = 0; j <= i; j+=2)
	{
	  double sum11 = pc[n*i+j];
	  double sum12 = pc[n*i+j+1];
	  double sum21 = pc[n*(i+1)+j];
	  double sum22 = pc[n*(i+1)+j+1];

	  double * lpa1 = pa + i * M;
	  double * lpa2 = pa + (i+1) * M;
	  double * lpb1 = pb + j * M;
	  double * lpb2 = pb + (j+1) * M;
	  
	  for (int k = 0; k < M; k++)
	    {
	      sum11 += lpa1[k] * lpb1[k];
	      sum12 += lpa1[k] * lpb2[k];
	      sum21 += lpa2[k] * lpb1[k];
	      sum22 += lpa2[k] * lpb2[k];
	    }
	  pc[j  +n*i    ] = sum11;
	  pc[j+1+n*i    ] = sum12;
	  pc[j  +n*(i+1)] = sum21;
	  pc[j+1+n*(i+1)] = sum22;
	}

    if (n % 2 == 1)
      {
	int i = n-1;

	for (int j = 0; j <= i; j++)
	  {
	    double sum = pc[n*i+j];
	    
	    double * lpa = pa + i * M;
	    double * lpb = pb + j * M;
	    
	    for (int k = 0; k < M; k++)
	      sum += lpa[k] * lpb[k];
	    
	    pc[j+n*i] = sum;
	  }
      }
    */



#ifndef __SSE3__
    for (int i = 0; i < n-1; i+=2)
      for (int j = 0; j <= i; j+=2)
	{
	  double sum11 = pc[n*i+j];
	  double sum12 = pc[n*i+j+1];
	  double sum21 = pc[n*(i+1)+j];
	  double sum22 = pc[n*(i+1)+j+1];

	  double * lpa1 = pa + i * M2;
	  double * lpa2 = pa + (i+1) * M2;
	  double * lpb1 = pb + j * M2;
	  double * lpb2 = pb + (j+1) * M2;
	  
	  for (int k = 0; k < M-1; k+=2)
	    {
	      sum11 += lpa1[k] * lpb1[k] + lpa1[k+1] * lpb1[k+1];
	      sum12 += lpa1[k] * lpb2[k] + lpa1[k+1] * lpb2[k+1];
	      sum21 += lpa2[k] * lpb1[k] + lpa2[k+1] * lpb1[k+1];
	      sum22 += lpa2[k] * lpb2[k] + lpa2[k+1] * lpb2[k+1];
	    }

	  if (M % 2)
	    {
	      sum11 += lpa1[M-1] * lpb1[M-1];
	      sum12 += lpa1[M-1] * lpb2[M-1];
	      sum21 += lpa2[M-1] * lpb1[M-1];
	      sum22 += lpa2[M-1] * lpb2[M-1];
	    }
	  pc[j  +n*i    ] = sum11;
	  pc[j+1+n*i    ] = sum12;
	  pc[j  +n*(i+1)] = sum21;
	  pc[j+1+n*(i+1)] = sum22;
	}

    if (n % 2 == 1)
      {
	int i = n-1;

	for (int j = 0; j <= i; j++)
	  {
	    double sum = pc[n*i+j];
	    
	    double * lpa = pa + i * M2;
	    double * lpb = pb + j * M2;
	    
	    for (int k = 0; k < M; k++)
	      sum += lpa[k] * lpb[k];
	    
	    pc[j+n*i] = sum;
	  }
      }

#endif



    /*
    for (int i = 0; i < n-1; i+=2)
      for (int j = 0; j <= i; j+=2)
	{
	  // Twins sum11(0,0);
	  __m128d sum11 = _mm_set_pd (0.0, 0.0);
	  Twins sum12(0,0);
	  Twins sum21(0,0);
	  Twins sum22(0,0);

	  double * lpa1 = pa + i * M;
	  double * lpa2 = pa + (i+1) * M;
	  double * lpb1 = pb + j * M;
	  double * lpb2 = pb + (j+1) * M;
	  
	  for (int k = 0; k < M-1; k+=2)
	    {
	      Twins a1(lpa1[k], lpa1[k+1]);
	      Twins a2(lpa2[k], lpa2[k+1]);
	      Twins b1(lpb1[k], lpb1[k+1]);
	      Twins b2(lpb2[k], lpb2[k+1]);

	      // sum11 = _mm_add_pd (sum11, _mm_mul_pd(a1, b1));
	      sum11 += _mm_mul_pd(a1, b1);
	      sum12 += a1*b2;
	      sum21 += a2*b1;
	      sum22 += a2*b2;
	    }

	  Twins sum11t = sum11;
	  pc[j  +n*i    ] += sum11t[0]+sum11t[1];
	  pc[j+1+n*i    ] += sum12[0]+sum12[1];
	  pc[j  +n*(i+1)] += sum21[0]+sum21[1];
	  pc[j+1+n*(i+1)] += sum22[0]+sum22[1];

	  if (M % 2)
	    {
	      pc[j  +n*i    ] += lpa1[M-1] * lpb1[M-1];
	      pc[j+1+n*i    ] += lpa1[M-1] * lpb2[M-1];
	      pc[j  +n*(i+1)] += lpa2[M-1] * lpb1[M-1];
	      pc[j+1+n*(i+1)] += lpa2[M-1] * lpb2[M-1];
	    }
	}


    if (n % 2 == 1)
      {
	int i = n-1;

	for (int j = 0; j <= i; j++)
	  {
	    double sum = pc[n*i+j];
	    
	    double * lpa = pa + i * M;
	    double * lpb = pb + j * M;
	    
	    for (int k = 0; k < M; k++)
	      sum += lpa[k] * lpb[k];
	    
	    pc[j+n*i] = sum;
	  }
      }
    */



#ifdef __SSE3__ 
    /*
    for (int i = 0; i < n-1; i+=2)
      for (int j = 0; j <= i; j+=2)
	{
	  __m128d sum11 = _mm_set1_pd (0.0);
	  __m128d sum12 = _mm_set1_pd (0.0);
	  __m128d sum21 = _mm_set1_pd (0.0);
	  __m128d sum22 = _mm_set1_pd (0.0);

	  double * lpa1 = pa + i * M;
	  double * lpa2 = pa + (i+1) * M;
	  double * lpb1 = pb + j * M;
	  double * lpb2 = pb + (j+1) * M;
	  
	  for (int k = 0; k < M-1; k+=2)
	    {
	      __m128d a1, a2, b1, b2;

	      a1 = _mm_load_pd(lpa1+k);
	      b1 = _mm_load_pd(lpb1+k);

	      if (M % 2) // only 8-byte alignment of matrix rows
		{
		  a2 = _mm_loadu_pd(lpa2+k);
		  b2 = _mm_loadu_pd(lpb2+k);
		}
	      else // 16-byte alignment guaranteed
		{
		  a2 = _mm_load_pd(lpa2+k);
		  b2 = _mm_load_pd(lpb2+k);
		}
	      sum11 += a1*b1; 
	      sum12 += a1*b2;
	      sum21 += a2*b1;
	      sum22 += a2*b2;
	    }

	  __m128d hsum1 = _mm_hadd_pd (sum11, sum12);
	  __m128d hsum2 = _mm_hadd_pd (sum21, sum22);

	  if (M % 2)
	    {
	      __m128d b = _mm_set_pd(lpb2[M-1], lpb1[M-1]);
	      __m128d a1 = _mm_set1_pd(lpa1[M-1]);
	      __m128d a2 = _mm_set1_pd(lpa2[M-1]);
	      hsum1 += a1 * b;
	      hsum2 += a2 * b;
	    }

	  hsum1 += _mm_loadu_pd(&pc[j+n*i]);
	  _mm_storeu_pd (&pc[j+n*i], hsum1);

	  hsum2 += _mm_loadu_pd(&pc[j+n*(i+1)]);
	  _mm_storeu_pd (&pc[j+n*(i+1)], hsum2);
	}


    if (n % 2 == 1)
      {
	int i = n-1;

	for (int j = 0; j <= i; j++)
	  {
	    double sum = pc[n*i+j];
	    
	    double * lpa = pa + i * M;
	    double * lpb = pb + j * M;
	    
	    for (int k = 0; k < M; k++)
	      sum += lpa[k] * lpb[k];
	    
	    pc[j+n*i] = sum;
	  }
      }
    */





#ifdef __AVX__ 



    // 2x4 blocks
    // job 995 calls        1, time 1.3397 sec, MFlops = 7464.08 fastmat - 27
    // job 996 calls        1, time 1.0635 sec, MFlops = 9403.24 fastmat - 30
    // job 997 calls        1, time 1.4389 sec, MFlops = 6949.72 fastmat - 29
    // job 998 calls        1, time 1.0690 sec, MFlops = 9354.53 fastmat - 28

    for (int i = 0; i < n-1; i+=2)
      for (int j = 0; j <= (i & -4); j+=4)
	{
	  __m256d sum11, sum12, sum13, sum14;
	  __m256d sum21, sum22, sum23, sum24;

	  sum11 = sum12 = sum13 = sum14 = _mm256_setzero_pd();
	  sum21 = sum22 = sum23 = sum24 = _mm256_setzero_pd();

	  double * lpa1 = pa + i * M2;
	  double * lpa2 = pa + (i+1) * M2;
	  double * lpb1 = pb + j * M2;
	  double * lpb2 = pb + (j+1) * M2;
	  double * lpb3 = pb + (j+2) * M2;
	  double * lpb4 = pb + (j+3) * M2;
	  
	  
	  __m128i mask = _mm_cmplt_epi32 (_mm_set_epi32(3,2,1,0), _mm_set1_epi32(n-j));
	  __m256i mask64;
	  mask64 = _mm256_insertf128_si256 (mask64, _mm_unpacklo_epi32(mask,mask),0);
	  mask64 = _mm256_insertf128_si256 (mask64, _mm_unpackhi_epi32(mask,mask),1);
 
	  for (int k = 0; k < M-3; k+=4)
	    {
	      __m256d a1, a2, b1, b2, b3, b4;

	      if (M2 % 4)
		{
		  a1 = _mm256_loadu_pd(lpa1+k);
		  a2 = _mm256_loadu_pd(lpa2+k);
		  b1 = _mm256_loadu_pd(lpb1+k);
		  b2 = _mm256_loadu_pd(lpb2+k);
		  b3 = _mm256_loadu_pd(lpb3+k);
		  b4 = _mm256_loadu_pd(lpb4+k);
		}
	      else
		{
		  a1 = _mm256_load_pd(lpa1+k);
		  a2 = _mm256_load_pd(lpa2+k);
		  b1 = _mm256_load_pd(lpb1+k);
		  b2 = _mm256_load_pd(lpb2+k);
		  b3 = _mm256_load_pd(lpb3+k);
		  b4 = _mm256_load_pd(lpb4+k);
		}
	      sum11 += a1*b1; sum12 += a1*b2;
	      sum13 += a1*b3; sum14 += a1*b4;

	      sum21 += a2*b1; sum22 += a2*b2;
	      sum23 += a2*b3; sum24 += a2*b4;
	    }

	  if (M % 4)
	    {
	      int k = M & -4;
	      __m256i mask64;
	      switch (M % 4)
		{
		case 1: mask64 = _mm256_set_epi64x(0,0,0,-1); break;
		case 2: mask64 = _mm256_set_epi64x(0,0,-1,-1); break;
		case 3: mask64 = _mm256_set_epi64x(0,-1,-1,-1); break;
		}

	      __m256d a1 = _mm256_maskload_pd(lpa1+k, mask64);
	      __m256d a2 = _mm256_maskload_pd(lpa2+k, mask64);
	      
	      __m256d b1 = _mm256_maskload_pd(lpb1+k, mask64);
	      __m256d b2 = _mm256_maskload_pd(lpb2+k, mask64);
	      __m256d b3 = _mm256_maskload_pd(lpb3+k, mask64);
	      __m256d b4 = _mm256_maskload_pd(lpb4+k, mask64);

	      sum11 += a1*b1; sum12 += a1*b2;
	      sum13 += a1*b3; sum14 += a1*b4;

	      sum21 += a2*b1; sum22 += a2*b2;
	      sum23 += a2*b3; sum24 += a2*b4;
	    }

	  __m256d hsum1 = HAdd (sum11, sum12, sum13, sum14);
	  __m256d hsum2 = HAdd (sum21, sum22, sum23, sum24);

	  hsum1 += _mm256_maskload_pd(&pc[j+n*i], mask64);
	  hsum2 += _mm256_maskload_pd(&pc[j+n*(i+1)], mask64);

	  _mm256_maskstore_pd (&pc[j+n*i], mask64, hsum1);
	  _mm256_maskstore_pd (&pc[j+n*(i+1)], mask64, hsum2);
	}
    
    if (n % 2)
      {
	int i = n-1;
	int j = 0;
	for ( ; j < n-3; j+=4)
	  {
	    __m256d sum11, sum12, sum13, sum14;

	    sum11 = sum12 = sum13 = sum14 = _mm256_setzero_pd();
	    
	    double * lpa1 = pa + i * M2;
	    double * lpb1 = pb + j * M2;
	    double * lpb2 = pb + (j+1) * M2;
	    double * lpb3 = pb + (j+2) * M2;
	    double * lpb4 = pb + (j+3) * M2;
	    
	    int k = 0;
	    for ( ; k < M-3; k+=4)
	      {
		__m256d a1 = _mm256_loadu_pd(lpa1+k);
		
		__m256d b1 = _mm256_loadu_pd(lpb1+k);
		__m256d b2 = _mm256_loadu_pd(lpb2+k);
		__m256d b3 = _mm256_loadu_pd(lpb3+k);
		__m256d b4 = _mm256_loadu_pd(lpb4+k);
		
		sum11 += a1*b1; sum12 += a1*b2;
		sum13 += a1*b3; sum14 += a1*b4;
	      }

	  __m256d hsum1 = HAdd (sum11, sum12, sum13, sum14);

	  for ( ; k < M; k++)
	    {
	      __m256d b = _mm256_set_pd (lpb4[k],lpb3[k], lpb2[k], lpb1[k]);
	      __m256d a1 = _mm256_set1_pd(lpa1[k]);
	      hsum1 += a1 * b;
	    }

	  hsum1 += _mm256_loadu_pd(&pc[j+n*i]);
	  _mm256_storeu_pd (&pc[j+n*i], hsum1);
	}

	for ( ; j <= i; j++)
	  {
	    double sum = pc[n*i+j];
	    
	    double * lpa = pa + i * M2;
	    double * lpb = pb + j * M2;
	    
	    for (int k = 0; k < M; k++)
	      sum += lpa[k] * lpb[k];
	    
	    pc[j+n*i] = sum;
	  }
      }




    /*

    // 1x4 blocks 
    // job 995 calls        1, time 1.8400 sec, MFlops = 5434.71 fastmat - 27
    // job 996 calls        1, time 1.4012 sec, MFlops = 7136.83 fastmat - 30
    // job 997 calls        1, time 1.9880 sec, MFlops = 5030.22 fastmat - 29
    // job 998 calls        1, time 1.2688 sec, MFlops = 7881.27 fastmat - 28
 
    for (int i = 0; i < n; i++)
      for (int j = 0; j <= i; j+=4)
	{
	  __m256d sum11, sum12, sum13, sum14;

	  sum11 = sum12 = sum13 = sum14 = _mm256_setzero_pd();

	  double * lpa1 = pa + i * M;
	  double * lpb1 = pb + j * M;
	  double * lpb2 = pb + (j+1) * M;
	  double * lpb3 = pb + (j+2) * M;
	  double * lpb4 = pb + (j+3) * M;
	  
	  
	  __m128i mask = _mm_cmplt_epi32 (_mm_set_epi32(3,2,1,0), _mm_set1_epi32(n-j));
	  __m256i mask64;
	  mask64 = _mm256_insertf128_si256 (mask64, _mm_unpacklo_epi32(mask,mask),0);
	  mask64 = _mm256_insertf128_si256 (mask64, _mm_unpackhi_epi32(mask,mask),1);
 
	  for (int k = 0; k < M-3; k+=4)
	    {
	      __m256d a1 = _mm256_loadu_pd(lpa1+k);
	      
	      __m256d b1 = _mm256_loadu_pd(lpb1+k);
	      __m256d b2 = _mm256_loadu_pd(lpb2+k);
	      __m256d b3 = _mm256_loadu_pd(lpb3+k);
	      __m256d b4 = _mm256_loadu_pd(lpb4+k);

	      sum11 += a1*b1; sum12 += a1*b2;
	      sum13 += a1*b3; sum14 += a1*b4;
	    }

	  if (M % 4)
	    {
	      int k = M & -4;
	      __m256i mask64;
	      switch (M % 4)
		{
		case 1: mask64 = _mm256_set_epi64x(0,0,0,-1); break;
		case 2: mask64 = _mm256_set_epi64x(0,0,-1,-1); break;
		case 3: mask64 = _mm256_set_epi64x(0,-1,-1,-1); break;
		}

	      __m256d a1 = _mm256_maskload_pd(lpa1+k, mask64);
	      
	      __m256d b1 = _mm256_maskload_pd(lpb1+k, mask64);
	      __m256d b2 = _mm256_maskload_pd(lpb2+k, mask64);
	      __m256d b3 = _mm256_maskload_pd(lpb3+k, mask64);
	      __m256d b4 = _mm256_maskload_pd(lpb4+k, mask64);

	      sum11 += a1*b1; sum12 += a1*b2;
	      sum13 += a1*b3; sum14 += a1*b4;
	    }

	  __m256d hsum1 = HAdd (sum11, sum12, sum13, sum14);

	  hsum1 += _mm256_maskload_pd(&pc[j+n*i], mask64);
	  _mm256_maskstore_pd (&pc[j+n*i], mask64, hsum1);
	}
    */


    /*
    // 4x4 blocks

    for (int i = 0; i < n-3; i+=4)
      for (int j = 0; j <= i; j+=4)
	{
	  __m256d sum11, sum12, sum13, sum14;
	  __m256d sum21, sum22, sum23, sum24;

	  sum11 = sum12 = sum13 = sum14 = _mm256_setzero_pd();
	  sum21 = sum22 = sum23 = sum24 = _mm256_setzero_pd();

	  double * lpa1 = pa + i * M;
	  double * lpa2 = pa + (i+1) * M;
	  double * lpa3 = pa + (i+2) * M;
	  double * lpa4 = pa + (i+3) * M;
	  double * lpb1 = pb + j * M;
	  double * lpb2 = pb + (j+1) * M;
	  double * lpb3 = pb + (j+2) * M;
	  double * lpb4 = pb + (j+3) * M;
	  
	  
	  __m128i mask = _mm_cmplt_epi32 (_mm_set_epi32(3,2,1,0), _mm_set1_epi32(n-j));
	  __m256i mask64;
	  mask64 = _mm256_insertf128_si256 (mask64, _mm_unpacklo_epi32(mask,mask),0);
	  mask64 = _mm256_insertf128_si256 (mask64, _mm_unpackhi_epi32(mask,mask),1);
 

	  for (int k = 0; k < M-3; k+=4)
	    {
	      __m256d a1 = _mm256_loadu_pd(lpa1+k);
	      __m256d a2 = _mm256_loadu_pd(lpa2+k);
	      __m256d a3 = _mm256_loadu_pd(lpa3+k);
	      __m256d a4 = _mm256_loadu_pd(lpa4+k);
	      
	      __m256d b1 = _mm256_loadu_pd(lpb1+k);
	      __m256d b2 = _mm256_loadu_pd(lpb2+k);
	      __m256d b3 = _mm256_loadu_pd(lpb3+k);
	      __m256d b4 = _mm256_loadu_pd(lpb4+k);

	      sum11 += _mm256_hadd_pd (a1*b1, a2*b1);
	      sum21 += _mm256_hadd_pd (a3*b1, a4*b1);

	      sum12 += _mm256_hadd_pd (a1*b2, a2*b2);
	      sum22 += _mm256_hadd_pd (a3*b2, a4*b2);

	      sum13 += _mm256_hadd_pd (a1*b3, a2*b3);
	      sum23 += _mm256_hadd_pd (a3*b3, a4*b3);

	      sum14 += _mm256_hadd_pd (a1*b4, a2*b4);
	      sum24 += _mm256_hadd_pd (a3*b4, a4*b4);
	    }

	  __m256d hsum1 = HAdd (sum11, sum12, sum13, sum14);
	  __m256d hsum2 = HAdd (sum21, sum22, sum23, sum24);

	  hsum1 += _mm256_maskload_pd(&pc[j+n*i], mask64);
	  hsum2 += _mm256_maskload_pd(&pc[j+n*(i+1)], mask64);

	  _mm256_maskstore_pd (&pc[j+n*i], mask64, hsum1);
	  _mm256_maskstore_pd (&pc[j+n*(i+1)], mask64, hsum2);
	}
*/



#else

    int i = 0;
    for ( ; i < n-3; i+=4)
      for (int j = 0; j <= i+2; j+=2)
	{
	  __m128d sum11, sum12, sum21, sum22;
	  __m128d sum31, sum32, sum41, sum42;

	  sum11 = sum12 = sum21 = sum22 = _mm_setzero_pd();
	  sum31 = sum32 = sum41 = sum42 = _mm_setzero_pd();

	  double * lpa1 = pa + i * M2;
	  double * lpa2 = pa + (i+1) * M2;
	  double * lpa3 = pa + (i+2) * M2;
	  double * lpa4 = pa + (i+3) * M2;
	  double * lpb1 = pb + j * M2;
	  double * lpb2 = pb + (j+1) * M2;
	  
	  for (int k = 0; k < M-1; k+=2)
	    {
	      __m128d a1, a2, a3, a4, b1, b2;
	      
	      a1 = _mm_load_pd(lpa1+k);
	      a3 = _mm_load_pd(lpa3+k);
	      b1 = _mm_load_pd(lpb1+k);

	      if (M2 % 2) // only 8-byte alignment of matrix rows
		{
		  a2 = _mm_loadu_pd(lpa2+k);
		  a4 = _mm_loadu_pd(lpa4+k);
		  b2 = _mm_loadu_pd(lpb2+k);
		}
	      else // 16-byte alignment guaranteed
		{
		  a2 = _mm_load_pd(lpa2+k);
		  a4 = _mm_load_pd(lpa4+k);
		  b2 = _mm_load_pd(lpb2+k);
		}

	      if (k == 0)
		{
		  sum11 = a1*b1; sum12 = a1*b2;
		  sum21 = a2*b1; sum22 = a2*b2;
		  sum31 = a3*b1; sum32 = a3*b2;
		  sum41 = a4*b1; sum42 = a4*b2;
		}
	      else
		{
		  sum11 += a1*b1; sum12 += a1*b2;
		  sum21 += a2*b1; sum22 += a2*b2;
		  sum31 += a3*b1; sum32 += a3*b2;
		  sum41 += a4*b1; sum42 += a4*b2;
		}
	    }

	  __m128d hsum1 = _mm_hadd_pd (sum11, sum12);
	  __m128d hsum2 = _mm_hadd_pd (sum21, sum22);
	  __m128d hsum3 = _mm_hadd_pd (sum31, sum32);
	  __m128d hsum4 = _mm_hadd_pd (sum41, sum42);

	  if (M % 2)
	    {
	      __m128d b = _mm_set_pd(lpb2[M-1], lpb1[M-1]);
	      __m128d a1 = _mm_set1_pd(lpa1[M-1]);
	      __m128d a2 = _mm_set1_pd(lpa2[M-1]);
	      __m128d a3 = _mm_set1_pd(lpa3[M-1]);
	      __m128d a4 = _mm_set1_pd(lpa4[M-1]);

	      if (M == 0)
		{
		  hsum1 = a1 * b;
		  hsum2 = a2 * b;
		  hsum3 = a3 * b;
		  hsum4 = a4 * b;
		}
	      else
		{
		  hsum1 += a1 * b;
		  hsum2 += a2 * b;
		  hsum3 += a3 * b;
		  hsum4 += a4 * b;
		}
	    }

	  hsum1 += _mm_loadu_pd(&pc[j+n*i]);
	  hsum2 += _mm_loadu_pd(&pc[j+n*(i+1)]);
	  hsum3 += _mm_loadu_pd(&pc[j+n*(i+2)]);
	  hsum4 += _mm_loadu_pd(&pc[j+n*(i+3)]);

	  _mm_storeu_pd (&pc[j+n*i], hsum1);
	  _mm_storeu_pd (&pc[j+n*(i+1)], hsum2);
	  _mm_storeu_pd (&pc[j+n*(i+2)], hsum3);
	  _mm_storeu_pd (&pc[j+n*(i+3)], hsum4);
	}


    for ( ; i < n; i++)
      {
	for (int j = 0; j <= i; j++)
	  {
	    double sum = pc[n*i+j];
	    
	    double * lpa = pa + i * M2;
	    double * lpb = pb + j * M2;
	    
	    for (int k = 0; k < M; k++)
	      sum += lpa[k] * lpb[k];
	    
	    pc[j+n*i] = sum;
	  }
      }
#endif


#endif




    // #define xxx
#ifdef xxx

    // war langsamer
    int i;
    for (i = 0; i < n-4; i+=4)
      for (int j = 0; j <= i+3; j++)
	{
	  
          /*
	  for (int k = 0; k <= M-2; k+=2)
            {
              sum1 += lpa1[k] * lpb[k] + lpa1[k+1] * lpb[k+1];
              sum2 += lpa2[k] * lpb[k] + lpa2[k+1] * lpb[k+1];
            }
	  if  (M % 2)
            {
              sum1 += lpa1[M-1] * lpb[M-1];
              sum2 += lpa2[M-1] * lpb[M-1];
            }	  
          */

	  double * lpb = pb + j * M;

	  double sum = pc[n*i+j];
	  double * lpa = pa + i * M;

	  for (int k = 0; k <= M-2; k+=2)
	    sum += lpa[k] * lpb[k] + lpa[k+1] * lpb[k+1];
	  if  (M % 2)
	    sum += lpa[M-1] * lpb[M-1];
	  // for (int k = 0; k < M; k++)
          // sum += lpa[k] * lpb[k];
	  pc[j+n*i] = sum;

	  lpa = pa + (i+1) * M;
	  sum = pc[n*(i+1)+j];
	  for (int k = 0; k <= M-2; k+=2)
	    sum += lpa[k] * lpb[k] + lpa[k+1] * lpb[k+1];
	  if  (M % 2)
	    sum += lpa[M-1] * lpb[M-1];
	  // for (int k = 0; k < M; k++)
          // sum += lpa[k] * lpb[k];
	  pc[j+n*(i+1)] = sum;

	  lpa = pa + (i+2) * M;
	  sum = pc[n*(i+2)+j];
	  for (int k = 0; k <= M-2; k+=2)
	    sum += lpa[k] * lpb[k] + lpa[k+1] * lpb[k+1];
	  if  (M % 2)
	    sum += lpa[M-1] * lpb[M-1];
	  // for (int k = 0; k < M; k++)
          // sum += lpa[k] * lpb[k];
	  pc[j+n*(i+2)] = sum;

	  lpa = pa + (i+3) * M;
	  sum = pc[n*(i+3)+j];
	  for (int k = 0; k <= M-2; k+=2)
	    sum += lpa[k] * lpb[k] + lpa[k+1] * lpb[k+1];
	  if  (M % 2)
	    sum += lpa[M-1] * lpb[M-1];
	  // for (int k = 0; k < M; k++)
          // sum += lpa[k] * lpb[k];
	  pc[j+n*(i+3)] = sum;
	}

    for ( ; i < n; i++)
      for (int j = 0; j <= i; j++)
	{
	  double sum = pc[n*i+j];

	  double * lpa = pa + i * M;
	  double * lpb = pb + j * M;
	  
	  for (int k = 0; k <= M-2; k+=2)
	    sum += lpa[k] * lpb[k] + lpa[k+1] * lpb[k+1];
	  if  (M % 2)
	    sum += lpa[M-1] * lpb[M-1];
	  
	  pc[j+n*i] = sum;
          // pc[i+n*j] = sum;
	}
#endif
  }


  /*
  template <int M>
  void FastMat (int n, 
		double * pa, 
		double * pb, 
		double * pc)
  {
    for (int i = 0; i < n; i++)
      for (int j = 0; j <= i; j++)
	{
	  if (M % 2 == 0 && M >= 4)
	    {
	      __m128d * mpa = (__m128d*) (pa + i * M);
	      __m128d * mpb = (__m128d*) (pb + j * M);
	      
	      __m128d sum1 = _mm_mul_pd ( mpa[0], mpb[0]);
	      __m128d sum2 = _mm_mul_pd ( mpa[1], mpb[1]);
	      
	      for (int k = 2; k < M/2-1; k+=2)
		{
		  sum1  = _mm_add_pd ( sum1, _mm_mul_pd ( mpa[k], mpb[k] ));
		  sum2  = _mm_add_pd ( sum2, _mm_mul_pd ( mpa[k+1], mpb[k+1] ));
		}
	      
	      if (M/2 & 1)
		sum1 = _mm_add_pd ( sum1, _mm_mul_pd ( mpa[M/2-1], mpb[M/2-1] ));
	      
	      sum1 =  _mm_add_pd ( sum1, sum2);
	      sum1  = _mm_add_sd ( sum1, _mm_shuffle_pd ( sum1, sum1, 1 ) );

	      sum1  = _mm_add_sd ( sum1, _mm_load_sd (&pc[n*i+j]));
	      _mm_store_sd ( &pc[n*i+j], sum1);
	      _mm_store_sd ( &pc[n*j+i], sum1);
	    }
	  else
	    {
	      double sum = pc[n*i+j];

	      double * lpa = pa + i * M;
	      double * lpb = pb + j * M;
	  
	      for (int k = 0; k <= M-2; k+=2)
		sum += lpa[k] * lpb[k] + lpa[k+1] * lpb[k+1];
	      if  (M % 2)
		sum += lpa[M-1] * lpb[M-1];
	      
	      pc[j+n*i] = sum;
	      pc[i+n*j] = sum;
	    }
	}
  }
  */








#ifdef __SSE3__

  class SSEComplex
  {
    __m128d data;
  public:
    SSEComplex (__m128d d) { data = d; }
    SSEComplex (double r) { data = _mm_set1_pd (r); }
    SSEComplex (Complex a) { data = _mm_set_pd(a.imag(), a.real()); }

    double real()
    {
      double hd; 
      _mm_store_sd (&hd, data);
      return hd;
    }
    double imag()
    {
      double hd; 
      _mm_storeh_pd (&hd, data);
      return hd;
    }

    operator Complex() const
    { 
      Complex h;
      _mm_storeu_pd ((double*)&h, data);
      return h;
    }
    __m128d Data() const { return data; }
    __m128d & Data() { return data; }
  };


  inline SSEComplex operator+ (SSEComplex a, SSEComplex b)
  {
    return a.Data() + b.Data();
  }

  inline SSEComplex operator+= (SSEComplex & a, SSEComplex b)
  {
    return a.Data() += b.Data();
  }

  inline SSEComplex operator* (SSEComplex a, SSEComplex b)
  {
    SSEComplex ar = a.real();
    SSEComplex ai = a.imag();
    SSEComplex rb = _mm_shuffle_pd (b.Data(), b.Data(), 1);

    return _mm_addsub_pd ( ar.Data()*b.Data(),  
                           ai.Data()*rb.Data() );
  }

  template <int M>
  inline SSEComplex Scal (SSEComplex * pa, SSEComplex * pb)
  {
    SSEComplex sum = 0;
    for (int i = 0; i < M; i++)
      sum += pa[i] * pb[i];
    return sum;
  }
#endif


  template <int M> 
  void FastMat (int n, int M2, Complex * pa, Complex * pb, Complex * pc)
  {
    static Timer timer ("Fastmat, complex", 2);
    RegionTimer reg (timer);
    timer.AddFlops (double(M)*n*n/2);
    
    for (int i = 0; i < n; i++)
      {
	Complex * hpa = pa + i*M2;
	Complex * hpc = pc+n*i;

	for (int j = 0; j < i; j++)
	  {
	    Complex * hpb = pb + j*M2;
	

#ifdef __SSE3___            
            Complex sum = *hpc + Scal<M> ((SSEComplex*)hpa, (SSEComplex*)hpb);
#else
	    Complex sum = *hpc;
	    for (int k = 0; k < M; k++)
	      sum += hpa[k] * hpb[k];
#endif

	    *hpc = sum;
	    pc[i+n*j] = sum;

	    hpb += M;
	    hpc++;
	  }

	Complex * hpb = pb + i*M2;
#ifdef __SSE3___            
        Complex sum = *hpc + Scal<M> ((SSEComplex*)hpa, (SSEComplex*)hpb);
#else
	Complex sum = *hpc;
	for (int k = 0; k < M; k++)
	  sum += hpa[k] * hpb[k];
#endif

	*hpc = sum;

	hpa += M;
      }
  }





  template <int M> 
  void FastMat (int n, int M2, Complex * pa, double * pb, Complex * pc)
  {
    static Timer timer ("Fastmat, complex-double", 2);
    RegionTimer reg (timer);
    timer.AddFlops (double(M)*n*n/2);
    
    for (int i = 0; i < n; i++)
      {
	Complex * hpa = pa + i * M2;
	Complex * hpc = pc+n*i;

	for (int j = 0; j < i; j++)
	  {
	    double * hpb = pb + j * M2;

	    Complex sum = *hpc;
	    for (int k = 0; k < M; k++)
	      sum += hpa[k] * hpb[k];

	    *hpc = sum;
	    pc[i+n*j] = sum;

	    hpb += M;
	    hpc++;
	  }

	double * hpb = pb + i * M2;
	Complex sum = *hpc;
	for (int k = 0; k < M; k++)
	  sum += hpa[k] * hpb[k];

	*hpc = sum;
	hpa += M;
      }
  }





  // template NGS_DLL_HEADER void FastMat<32> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);



  template NGS_DLL_HEADER void FastMat<1> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<2> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<3> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<4> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<5> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<6> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<7> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<8> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<9> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<10> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);



  template NGS_DLL_HEADER void FastMat<25> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<26> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<27> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<28> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<29> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<30> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<40> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);

  /*
  // roundup
  template NGS_DLL_HEADER void FastMat<1,4> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<1,28> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<2,4> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<2,28> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<3,4> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<3,8> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<3,32> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<4,28> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<6,8> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<6,12> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<6,32> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<6,36> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<8,28> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<12,32> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<12,36> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<24,36> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<26,28> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<27,28> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<29,32> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<30,32> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  */



  template NGS_DLL_HEADER void FastMat<12> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<18> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<24> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<32> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<36> (int n, int M2, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template NGS_DLL_HEADER void FastMat<48> (int n, int M2, double * __restrict__  pa, double * __restrict__  pb, double * __restrict__  pc);




  template NGS_DLL_HEADER void FastMat<1> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<2> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<3> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<4> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<5> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<6> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<9> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);


  template NGS_DLL_HEADER void FastMat<25> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<26> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<27> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<28> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<29> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<30> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);


  template NGS_DLL_HEADER void FastMat<12> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<18> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<24> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<32> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<36> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<48> (int n, int M2, Complex * pa, Complex * pb, Complex * pc);








  template NGS_DLL_HEADER void FastMat<1> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<2> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<3> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<4> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<5> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<6> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<8> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<9> (int n, int M2, Complex * pa, double * pb, Complex * pc);


  template NGS_DLL_HEADER void FastMat<25> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<26> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<27> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<28> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<29> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<30> (int n, int M2, Complex * pa, double * pb, Complex * pc);


  template NGS_DLL_HEADER void FastMat<12> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<18> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<24> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<32> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<36> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<48> (int n, int M2, Complex * pa, double * pb, Complex * pc);

  /*
  // roundup
  template NGS_DLL_HEADER void FastMat<1,4> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<1,28> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<2,4> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<2,28> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<3,4> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<3,8> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<3,32> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<4,28> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<6,8> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<6,12> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<6,32> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<6,36> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<8,28> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<12,32> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<12,36> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<24,36> (int n, int M2, Complex * pa, double * pb, Complex * pc);

  template NGS_DLL_HEADER void FastMat<25,28> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<26,28> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<27,28> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<29,32> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  template NGS_DLL_HEADER void FastMat<30,32> (int n, int M2, Complex * pa, double * pb, Complex * pc);
  */



  void FastMatN (int n, int M, double * pa, double * pb, double * pc)
  {
    double * hpa = pa;

    for (int i = 0; i < n; i++)
      {
	double * hpb = pb;
	double * hpc = pc+n*i;

	for (int j = 0; j < i; j++)
	  {
	    double sum = *hpc;
	  
	    for (int k = 0; k < M; k++)
	      sum += hpa[k] * hpb[k];
	  
	    *hpc = sum;
	    pc[i+n*j] = sum;

	    hpb += M;
	    hpc++;
	  }

	double sum = *hpc;
      
	for (int k = 0; k < M; k++)
	  sum += hpa[k] * hpb[k];
      
	*hpc = sum;

	hpa += M;
      }
  }


}
