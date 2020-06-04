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
    BSpline (int aorder, 
             Array<double> at,
             Array<double> ac);
    BSpline (const BSpline &) = default;
    BSpline & operator= (const BSpline &) = default;

    void DoArchive(Archive& ar)
    {
      ar & order & t & c;
    }

    BSpline Differentiate () const;
    BSpline Integrate () const;

    double Evaluate (double x) const;
    double operator() (double x) const { return Evaluate(x); }
    AutoDiff<1> operator() (AutoDiff<1> x) const;
    AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const;
    SIMD<double> operator() (SIMD<double> x) const;
    
    friend ostream & operator<< (ostream & ost, const BSpline & sp);
  };

  extern ostream & operator<< (ostream & ost, const BSpline & sp);
}

#endif
