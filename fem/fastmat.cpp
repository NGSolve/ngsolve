not used anymore


#define FILE_FASTMAT_CPP

#include <fem.hpp>

#include "fastmat.hpp"

namespace ngfem
{
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
