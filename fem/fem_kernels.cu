#include <ngstd.hpp>
#include <bla.hpp>
#include "recursive_pol.hpp"
using namespace ngfem;
using namespace ngs_cuda;

namespace ngfem
{
  __device__ double legendre_coefs[1000][2];
  __device__ double jacobialpha_coefs[100][100][4];
  __device__ int jacobialpha_maxn;
}

__global__ void InitLegendre (int n)
{
  // legendre_coefs = new Vec<2>[n+1];
  legendre_coefs[0][0] = 1;
  legendre_coefs[1][1] = 1; 
  for (int i = 2; i <= 1000; i++)
    {
      legendre_coefs[i][0] = (2.0*i-1)/i;
      legendre_coefs[i][1] = (1.0-i)/i;
    }

  // int alpha = 2*n+2;
  int alpha = 100;
  // jacobialpha_coefs = new Vec<3>[(n+1)*(alpha+1)];

  for (int a = 0; a <= alpha; a++)
    {
      for (int i = 1; i <= n; i++)
	{
	  /*
	  jacobialpha_coefs[a*(n+1)+i][0] = JacobiPolynomialAlpha::CalcA (i, a, 0);
	  jacobialpha_coefs[a*(n+1)+i][1] = JacobiPolynomialAlpha::CalcB (i, a, 0);
	  jacobialpha_coefs[a*(n+1)+i][2] = JacobiPolynomialAlpha::CalcC (i, a, 0);
	  */
	  jacobialpha_coefs[a][i][0] = JacobiPolynomialAlpha::CalcA (i, a, 0);
	  jacobialpha_coefs[a][i][1] = JacobiPolynomialAlpha::CalcB (i, a, 0);
	  jacobialpha_coefs[a][i][2] = JacobiPolynomialAlpha::CalcC (i, a, 0);
	}
      // ALWAYS_INLINE S P1(S x) const { return 0.5 * (2*(al+1)+(al+be+2)*(x-1)); }
      double al = a, be = 0;
      /*
      jacobialpha_coefs[a*(n+1)+1][0] = 0.5 * (al+be+2);
      jacobialpha_coefs[a*(n+1)+1][1] = 0.5 * (2*(al+1)-(al+be+2));
      */
      jacobialpha_coefs[a][1][0] = 0.5 * (al+be+2);
      jacobialpha_coefs[a][1][1] = 0.5 * (2*(al+1)-(al+be+2));
    }
  jacobialpha_maxn = n;
}


class InitFemKernels
{
public:
  InitFemKernels () 
  { 
    cout << "Init device Legendre and Jacobi coefficients ... " << flush;
    InitLegendre<<<1,1>>> (1000); 
    cout << "done" << endl;
  }
};


InitFemKernels ifc;
