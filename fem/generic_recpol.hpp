#ifndef FILE_GENERIC_RECPOL
#define FILE_GENERIC_RECPOL


/*********************************************************************/
/* File:   generic_recpol.hpp                                        */
/* Author: Veronika Pillwein and Joachim Schöberl                    */
/* Date:   5. Jan. 2005                                              */
/*********************************************************************/

namespace ngfem
{


  /**
     Generic recursive polynomial.

     \begin{verbatim}
     p_0 = a_0
     p_1 = (a_1 + b_1 x) p_0
     p_i = (a_i + b_i x) p_{i-1} + c_i p_{i-2}
     \end{verbatim}
  */

  class RecPol
  {
  protected:
    ///
    // Array< double[3] > coefs;
    Array< Vec<3> > coefs;
  
  public:
    /// init coefficient array
    RecPol (int n) : coefs(n+1) 
    { 
      for (int i = 0; i <= n; i++)
	coefs[i][0] = coefs[i][1] = coefs[i][2] = 0.0;
    }

    /// maximal order
    int Order() const { return coefs.Size()-1; }
  
    /// evaluate rec pols up to order n
    template <class S, class T>
    void Evaluate (int n, S x, T & values) const
    {
      // T & values = const_cast<T&> (const_values);
      S v1, v2, v3;

      if (n < 0) return;
      values[0] = v2 = coefs[0][0];

      if (n == 0) return;
      values[1] = v1 = (coefs[1][0] + coefs[1][1] * x) * v2;
      /*
      // hat nichts gebracht
      int i;
      for (i = 2; i < n; i+=2)
        {
          v2 *= coefs[i][2];
          v2 += (coefs[i][0] + coefs[i][1] * x) * v1;
          values[i] = v2;
          
          v1 *= coefs[i+1][2];
          v1 += (coefs[i+1][0] + coefs[i+1][1] * x) * v2;
          values[i+1] = v1;
        }

      if (i <= n)
        {
          v2 *= coefs[i][2];
          v2 += (coefs[i][0] + coefs[i][1] * x) * v1;
          values[i] = v2;
        }
      */
      for (int i = 2; i <= n; i++)
	{
	  v3 = v2; v2 = v1;
	  v1 = (coefs[i][0] + coefs[i][1] * x) * v2 + coefs[i][2] * v3;
	  values[i] = v1;
	}

    }


    /*
      void Evaluate (int n, double x, double * values) const
      {
      double v1, v2, v3;

      if (n >= 0)
      v2 = values[0] = coefs[0][0];
      if (n >= 1)
      v1 = values[1] = (coefs[1][0] + coefs[1][1] * x) * values[0];
    
      for (int i = 2; i <= n; i++)
      {
      v3 = v2; v2 = v1;
      v1 = (coefs[i][0] + coefs[i][1] * x) * v2 + coefs[i][2] * v3;
      values[i] = v1;
      }
      }
    */

    /**
       Multiply polynomial with the bubble $1-x^2$.
       p --> q
       q_0 = 1, q_1 = x
       q_i = (1-x^2) p_{i-2}
    */
    void MultBubble ();

    /**
       Multiply polynomial with the linear factor $a + bx$.
       p --> q
       q_0 = 1
       q_i = (a + b x) p_i-1
    */
    void MultLinear (double a, double b);

  
    void Scale (double s)
    { coefs[0][0] *= s; }


    double A(int i) const { return coefs[i][0]; }
    double B(int i) const { return coefs[i][1]; }
    double C(int i) const { return coefs[i][2]; }
  };


  /// output coefficients 
  inline ostream & operator<< (ostream & ost, const RecPol & pol)
  {
    for (int i = 0; i <= pol.Order(); i++)
      ost << "a(" << i << ") = " << pol.A(i)
	  << ", b(" << i << ") = " << pol.B(i)
	  << ", c(" << i << ") = " << pol.C(i) << endl;
    return ost;
  }


  /**
     Initialize coefficients with x^i
  */
  class RecPolMonomial : public RecPol
  {
  public:
    ///
    RecPolMonomial (int n)
      : RecPol (n)
    {
      coefs[0][0] = 1;
      for (int i = 1; i <= n; i++)
	coefs[i][1] = 1;
    }
  };


  /**
     Initialize coefficients with Jacobi polynomials.
  */
  class RecPolJacobi : public RecPol
  {
  public:
    RecPolJacobi (int n, double alpha, double beta)
      : RecPol (n)
    {
      coefs[0][0] = 1;
      coefs[1][0] = alpha+1 - 0.5 * (alpha + beta + 2);
      coefs[1][1] = 0.5 * (alpha+beta+2);

      for (int i = 2; i <= n; i++)
	{
	  double fac = 1.0 / (2.0*i* (i+alpha+beta)*(2*i-2+alpha+beta));	
	  coefs[i][0] = ( 2*i-1+alpha+beta ) * ( alpha*alpha-beta*beta ) * fac;
	  coefs[i][1] = ( 2*i-2+alpha+beta ) * ( 2*i-1+alpha+beta ) * ( 2*i+alpha+beta ) * fac;
	  coefs[i][2] = -2 * (i-1+alpha) * (i-1+beta) * (2*i+alpha+beta) * fac;
	}
    }
  };


  /**
     Initialize coefficients with Legendre polynomials.
  */
  class RecPolLegendre : public RecPol
  {
  public:
    ///
    RecPolLegendre (int n)
      : RecPol (n)
    {
      coefs[0][0] = 1;
      coefs[1][1] = 1;
      for (int i = 2; i <= n; i++)
	{
	  coefs[i][1] = (2.0*i-1)/i;
	  coefs[i][2] = (1.0-i)/i;
	}
    }
  };










  ///
  extern void GenerateMatrix (const RecPol & pol1, const RecPol & pol2,
			      FlatVector<> pts, FlatVector<> coefs,
			      FlatMatrix<> mat);




  ///
  extern void GenerateMatrix (const RecPol & pol1, FlatMatrix<> pol1_vals,
			      const RecPol & pol2, FlatMatrix<> pol2_vals,
			      FlatVector<> coefs, FlatMatrix<> mat);



}

#endif

