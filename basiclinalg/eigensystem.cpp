#include <bla.hpp>

namespace ngbla
{

  void CalcEigenSystem (const FlatMatrix<double> mat1,
			FlatVector<double> lami,
			FlatMatrix<double> eigenvecs)
  {
    int n = mat1.Height();    
    Matrix<double> mat(n, n);
    mat = mat1;

    eigenvecs = 0;
    for (int i = 0; i < n; i++)
      eigenvecs(i,i) = 1;
    
    for (int l = 0; l < 100; l++)
      for (int i = 0; i < n; i++)
	for (int j = 0; j < i; j++)
	  {
	    // find eigensystem of a(i-j,i-j)

	    double a11 = mat(i,i);
	    double a12 = mat(i,j);
	    double a22 = mat(j,j);

	    if (a12*a12 <= 1e-32 * fabs(a11*a22)) continue;

	    /*
	    double lam1, lam2;
	    // y is EV from  a y = lam y
	    // quadratic eq. for eigenvalues lam:
	    // c0 + c1 lam + c2 lam^2 = 0
	  
	    double c0 = a11*a22-a12*a12;
	    double c1 = -a22 - a11;
	    double c2 = 1;
	  
	    lam1 = -c1/(2*c2) + sqrt( sqr(c1/(2*c2)) - c0/c2);
	    lam2 = -c1/(2*c2) - sqrt( sqr(c1/(2*c2)) - c0/c2);
	    cout << "lam1,2 = " << lam1 << ", " << lam2 << endl;

	    lam1 = (a11+a22) / 2 + sqrt ( sqr(a11-a22)/4 + a12*a12);
	    lam2 = (a11+a22) / 2 - sqrt ( sqr(a11-a22)/4 + a12*a12);
	    //	cout << "lam1,2 = " << lam1 << ", " << lam2 << endl;
	    */

	    double p = (a22-a11)/2;
	    double q = a12;
	
	    // compute eigenvectors:
	    double y11, y12, y21, y22, y;

	    y11 = a12;
	    //	y12 = lam1-a11;

	    y12 = 
	      (p >= 0) ?
	      p + sqrt (p*p + q*q) :
	      -q*q / (p - sqrt (p*p+q*q));

	    y = sqrt (y11*y11+y12*y12);
	    y11 /= y;
	    y12 /= y;

	    y21 = a12;
	    //	y22 = lam2-a11;

	    y22 = 
	      (p <= 0) ?
	      p - sqrt (p*p + q*q) :
	      -q*q / (p + sqrt (p*p+q*q));

	    y = sqrt (y21*y21+y22*y22);
	    y21 /= y;
	    y22 /= y;
	
	    /*
	      (*testout) << "evecs = "
	      << "(" << y11 << ", " << y12 << "), "
	      << "(" << y21 << ", " << y22 << ")" << endl;
	      (*testout) << "Y Y = "
	      << (y11*y11+y12*y12) << ", "
	      << (y11*y21+y12*y22) << ", "
	      << (y21*y21+y22*y22) << endl;
	    */
	
	    // V^T A V = V^T G^{-1}  (G^T A G)  G^{-1} V

	    for (int k = 0; k < n; k++)
	      {
		double v1 = mat(k,i);
		double v2 = mat(k,j);
		mat(k,i) = v1 * y11 + v2 * y12;
		mat(k,j) = v1 * y21 + v2 * y22;
	      }

	    for (int k = 0; k < n; k++)
	      {
		double v1 = mat(i,k);
		double v2 = mat(j,k);
		mat(i,k) = v1 * y11 + v2 * y12;
		mat(j,k) = v1 * y21 + v2 * y22;
	      }

	    mat(i,j) = 0;
	    mat(j,i) = 0;

	    for (int k = 0; k < n; k++)
	      {
		double v1 = eigenvecs(i,k);
		double v2 = eigenvecs(j,k);
		eigenvecs(i,k) = v1 * y11 + v2 * y12;
		eigenvecs(j,k) = v1 * y21 + v2 * y22;
	      }
	  }

    for (int i = 0; i < n; i++)
      lami(i) = mat(i,i);
  }
}
