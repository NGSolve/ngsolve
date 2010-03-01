#define NO_PARALLEL_THREADS

// #include "bla.hpp"
#ifdef SSE
#include <emmintrin.h>
#endif

#include <fem.hpp>

#ifdef WIN32
#define __restrict__ __restrict
#endif

namespace ngfem {
  
  typedef std::complex<double> Complex;



template <int M> NGS_DLL_HEADER
void FastMat (int n, Complex * ba, Complex *  pb, Complex * pc);

template <int M> NGS_DLL_HEADER
void FastMat (int n, double * __restrict__ ba, double *  __restrict__ pb, double * __restrict__ pc);
  

  template <int M>
  void FastMat (int n, 
                double * __restrict__ pa, 
		double * __restrict__ pb, 
		double * __restrict__ pc)
  {
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


  template <int M> 
  void FastMat (int n, Complex * pa, Complex * pb, Complex * pc)
  {
    Complex * hpa = pa;

    for (int i = 0; i < n; i++)
      {
	Complex * hpb = pb;
	Complex * hpc = pc+n*i;

	for (int j = 0; j < i; j++)
	  {
	    Complex sum = *hpc;
	  
	    for (int k = 0; k < M; k++)
	      sum += hpa[k] * hpb[k];
	  
	    *hpc = sum;
	    pc[i+n*j] = sum;

	    hpb += M;
	    hpc++;
	  }

	Complex sum = *hpc;
      
	for (int k = 0; k < M; k++)
	  sum += hpa[k] * hpb[k];
      
	*hpc = sum;

	hpa += M;
      }
  }


  template void FastMat<1> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<2> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<3> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<4> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<5> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<6> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<7> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<8> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<9> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<10> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);



  template void FastMat<25> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<26> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<27> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<28> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<29> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<30> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<40> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);


  template void FastMat<12> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<18> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<24> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<32> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<36> (int n, double * __restrict__ pa, double * __restrict__ pb, double * __restrict__ pc);
  template void FastMat<48> (int n, double * __restrict__  pa, double * __restrict__  pb, double * __restrict__  pc);




  template void FastMat<1> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<2> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<3> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<4> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<5> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<6> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<9> (int n, Complex * pa, Complex * pb, Complex * pc);


  template void FastMat<25> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<26> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<27> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<28> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<29> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<30> (int n, Complex * pa, Complex * pb, Complex * pc);


  template void FastMat<12> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<18> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<24> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<32> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<36> (int n, Complex * pa, Complex * pb, Complex * pc);
  template void FastMat<48> (int n, Complex * pa, Complex * pb, Complex * pc);

  //  template void FastMat<100> (int n, Complex * pa, Complex * pb, Complex * pc);







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
