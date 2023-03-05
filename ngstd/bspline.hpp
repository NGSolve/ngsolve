#ifndef FILE_BSPLINE
#define FILE_BSPLINE

/**************************************************************************/
/* File:   bspline.hpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Aug. 2015                                                  */
/**************************************************************************/


namespace ngstd
{

  class BSpline 
  {
    int order;
    Array<double> t;
    Array<double> c;
    
  public:
    BSpline() = default;
    NGS_DLL_HEADER BSpline (int aorder,
             Array<double> at,
             Array<double> ac);
    NGS_DLL_HEADER BSpline (const BSpline &) = default;
    NGS_DLL_HEADER BSpline & operator= (const BSpline &) = default;

    void DoArchive(Archive& ar)
    {
      ar & order & t & c;
    }

    NGS_DLL_HEADER BSpline Differentiate () const;
    NGS_DLL_HEADER BSpline Integrate () const;

    NGS_DLL_HEADER double Evaluate (double x) const;
    NGS_DLL_HEADER SIMD<double> Evaluate(SIMD<double> x) const;
    double operator() (double x) const { return Evaluate(x); }
    // Complex operator() (Complex x) const { return Evaluate(x.real()); }
    /*I had to explicitly write double and SIMD<double>
     here as these functions will be used in GenerateCode and I was
     getting the error "templates must have C++ linkage"*/
    NGS_DLL_HEADER AutoDiff<1,double> operator() (AutoDiff<1,double> x) const;
    NGS_DLL_HEADER AutoDiffDiff<1,double> operator() (AutoDiffDiff<1,double> x) const;
    NGS_DLL_HEADER AutoDiff<1,SIMD<double>> operator() (AutoDiff<1,SIMD<double>> x) const;
    NGS_DLL_HEADER AutoDiffDiff<1,SIMD<double>> operator() (AutoDiffDiff<1,SIMD<double>> x) const;
    
    NGS_DLL_HEADER SIMD<double> operator() (SIMD<double> x) const;
    // NGS_DLL_HEADER SIMD<Complex> operator() (SIMD<Complex> x) const;
    
    friend ostream & operator<< (ostream & ost, const BSpline & sp);
  };

  extern ostream & operator<< (ostream & ost, const BSpline & sp);
}

#endif
