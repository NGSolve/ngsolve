/*********************************************************************/
/* File:   generic_recpol.cpp                                        */
/* Author: Veronika Pillwein and Joachim Schöberl                    */
/* Date:   5. Jan. 2005                                              */
/*********************************************************************/


#ifdef NONE

#include <fem.hpp>
namespace ngfem
{
  using namespace ngfem;

  void RecPol :: MultBubble ()
  {
    for (int i = coefs.Size()-1; i >= 3; i--)
      {
	coefs[i][0] = coefs[i-2][0];
	coefs[i][1] = coefs[i-2][1];
	coefs[i][2] = coefs[i-2][2];
      }
    coefs[2][0] = 0;
    coefs[2][1] = -coefs[0][0];
    coefs[2][2] = coefs[0][1];
    coefs[0][0] = 1;
    coefs[0][1] = 0;
    coefs[0][2] = 0;
    coefs[1][0] = 1;
    coefs[1][1] = 0;
    coefs[1][2] = 0;
  }


  void RecPol :: MultLinear (double a, double b)
  {
    for (int i = coefs.Size()-1; i >= 2; i--)
      {
	coefs[i][0] = coefs[i-1][0];
	coefs[i][1] = coefs[i-1][1];
	coefs[i][2] = coefs[i-1][2];
      }
    coefs[1][0] = a * coefs[0][0];
    coefs[1][1] = b * coefs[0][0];
    coefs[1][2] = 0;

    coefs[0][0] = 1;
    coefs[0][1] = 0;
    coefs[0][2] = 0;
  }


  void GenerateMatrix (const RecPol & pol1, const RecPol & pol2,
		       FlatVector<> pts, FlatVector<> coefs,
		       FlatMatrix<> mat)

  {
    int h = mat.Height();
    int w = mat.Width();

    ArrayMem<double, 20> val1(h);
    ArrayMem<double, 20> val2(w);


    // the standard algorithm
    mat = 0.0;

    if (1 || h < 2 || w < 2)
      {
	for (int i = 0; i < pts.Size(); i++)
	  {
            pol1.Evaluate (h-1, pts(i), val1);
	    pol2.Evaluate (w-1, pts(i), val2);

	    for (int j = 0; j < h; j++)
	      for (int k = 0; k < w; k++)
		mat(j, k) += coefs(i) * val1[j] * val2[k];
	  }
	return;
      }

    // row(0):
    mat = 0.0;
    for (int i = 0; i < pts.Size(); i++)
      {
	pol1.Evaluate (h-1, pts(i), val1);
	pol2.Evaluate (w-1, pts(i), val2);
      
	double fac0 = coefs(i) * val1[0];
	double fac1 = coefs(i) * val1[1];

	for (int k = 0; k < w; k++)
	  {
	    mat(0,k) += fac0 * val2[k];
	    mat(1,k) += fac1 * val2[k];
	  }

	fac0 = coefs(i) * val2[0];
	fac1 = coefs(i) * val2[w-1];

	for (int j = 2; j < h; j++)
	  {
	    mat(j,0) += fac0 * val1[j];
	    mat(j,w-1) += fac1 * val1[j];
	  }
      }

    /*
      for (int j = 2; j < h; j++)
      for (int k = 1; k < w-1; k++)
      mat(j,k) = pol1.A(j) * mat(j-1,k) + pol1.C(j) * mat(j-2,k)
      + pol1.B(j) / pol2.B(k+1) * (mat(j-1,k+1) - pol2.A(k+1)*mat(j-1,k) - pol2.C(k+1)*mat(j-1,k-1));
    */

    for (int j = 2; j < h; j++)
      {
	double * rowj = &mat(j,0);
	double * rowj1 = &mat(j-1,0);
	double * rowj2 = &mat(j-2,0);

	double aj = pol1.A(j);
	double bj = pol1.B(j);
	double cj = pol1.C(j);

	for (int k = 1; k < w-1; k++)
	  rowj[k] = aj * rowj1[k] + cj * rowj2[k]
	    + bj / pol2.B(k+1) * (rowj1[k+1] - pol2.A(k+1)* rowj1[k] - pol2.C(k+1)*rowj1[k-1]);
      }
  }












  void GenerateMatrix (const RecPol & pol1, FlatMatrix<> pol1_vals,
		       const RecPol & pol2, FlatMatrix<> pol2_vals,
		       FlatVector<> coefs,
		       FlatMatrix<> mat)

  {
	      
    int h = mat.Height();
    int w = mat.Width();
    // int wp1 = pol1_vals.Width();
    int wp2 = pol2_vals.Width();

    if (h == 0 || w == 0)
      return;

    // NgProfiler::StartTimer (110);

    ArrayMem<double, 30> prod(coefs.Size());

    /*
    // the standard algorithm
    mat = 0.0;
    for (int j = 0; j < h; j++)
      for (int k = 0; k < w; k++)
	for (int i = 0; i < coefs.Size(); i++)
	  mat(j, k) += coefs(i) * pol1_vals(j,i) * pol2_vals(k,i);
    return;
    */

    // row 0
    for (int i = 0; i < coefs.Size(); i++)
      prod[i] = coefs[i] * pol1_vals(i);
    
    for (int k = 0; k < w; k++)
      {
	double sum = 0;
	int base = k*wp2;
	for (int i = 0; i < coefs.Size(); i++)
	  sum += prod[i] * pol2_vals(base+i);
	mat(0, k) = sum;
      }

    //    NgProfiler::StopTimer (110);

    if (h == 1)
      return;

    //    NgProfiler::StartTimer (111);

    if (w == 1)
      {  
	for (int i = 0; i < coefs.Size(); i++)
	  prod[i] = coefs[i] * pol2_vals(0, i);
	
	for (int k = 0; k < h; k++)
	  {
	    double sum = 0;
	    for (int i = 0; i < coefs.Size(); i++)
	      sum += prod[i] * pol1_vals(k, i); 
	    mat(k, 0) = sum;
	  }
	NgProfiler::StopTimer (110);
	return;
      }


    // col w-1
    for (int i = 0; i < coefs.Size(); i++)
      prod[i] = coefs[i] * pol2_vals(w-1, i);
    for (int j = 1; j < h; j++)
      {
	double sum = 0;
	for (int i = 0; i < coefs.Size(); i++)
	  sum += prod[i] * pol1_vals(j, i);
	mat(j, w-1) = sum;
      }

    //    NgProfiler::StopTimer (111);
    //    NgProfiler::StartTimer (112);

    
    {
      // row 1:
      
      int j = 1;
      double * rowj = &mat(j,0);
      double * rowj1 = &mat(j-1,0);
      
      double aj = pol1.A(j);
      double bj = pol1.B(j);
      // double cj = pol1.C(j);

      // col 0:
      rowj[0] = aj * rowj1[0]
	+ bj / pol2.B(1) * (rowj1[1] - pol2.A(1)* rowj1[0]);

      for (int k = 1; k < w-1; k++)
	rowj[k] = aj * rowj1[k]
	  + bj / pol2.B(k+1) * (rowj1[k+1] - pol2.A(k+1)* rowj1[k] - pol2.C(k+1)*rowj1[k-1]);
    }


    for (int j = 2; j < h; j++)
      {
	double * rowj = &mat(j,0);
	double * rowj1 = &mat(j-1,0);
	double * rowj2 = &mat(j-2,0);

	double aj = pol1.A(j);
	double bj = pol1.B(j);
	double cj = pol1.C(j);

	// col 0
	rowj[0] = aj * rowj1[0] + cj * rowj2[0]
	  + bj / pol2.B(1) * (rowj1[1] - pol2.A(1)* rowj1[0]);

	// cols 1 ... w-2
	for (int k = 1; k < w-1; k++)
	  rowj[k] = aj * rowj1[k] + cj * rowj2[k]
	    + bj / pol2.B(k+1) * (rowj1[k+1] - pol2.A(k+1)* rowj1[k] - pol2.C(k+1)*rowj1[k-1]);
      }

    //    NgProfiler::StopTimer (112);

  }





}


#endif
