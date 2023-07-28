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

  class BSpline2D
  {
    Array<double> px;
    Array<double> py;
    Array<double> v;
    [[maybe_unused]] int order = 1;
    bool extrapolate = true;

    std::tuple<int,int,bool> Search(double x, double y) const;

    template<typename T>
    T Interpolate(int ix, int iy, T x, T y) const
    {
      double x0 = px[ix];
      double x1 = px[ix+1];
      double y0 = py[iy];
      double y1 = py[iy+1];
      double dx = x1-x0;
      double dy = y1-y0;

      auto fy0 = ((x1-x)*GetV(ix, iy) + (x-x0)*GetV(ix+1, iy))/dx;
      auto fy1 = ((x1-x)*GetV(ix, iy+1) + (x-x0)*GetV(ix+1, iy+1))/dx;
      return ((y1-y)*fy0 + (y-y0)*fy1)/dy;
    }

    double GetV(int ix, int iy) const
    {
      return v[ix*py.Size()+iy];
    }

  public:
    BSpline2D() = default;
    NGS_DLL_HEADER BSpline2D ( Array<double> x, Array<double> y, Array<double> av, int aorder, bool aextrapolate)
      : px(x), py(y), v(av), order(aorder), extrapolate(aextrapolate) { }
    NGS_DLL_HEADER BSpline2D (const BSpline2D &) = default;
    NGS_DLL_HEADER BSpline2D & operator= (const BSpline2D &) = default;

    void DoArchive(Archive& ar)
    {
      ar & px & py & v;
    }

    NGS_DLL_HEADER double Evaluate (double x, double y) const;
    NGS_DLL_HEADER SIMD<double> Evaluate(SIMD<double> x, SIMD<double> y) const;
    double operator() (double x, double y) const { return Evaluate(x, y); }
    NGS_DLL_HEADER AutoDiff<1,double> operator() (AutoDiff<1,double> x, AutoDiff<1,double> y) const;
    NGS_DLL_HEADER AutoDiffDiff<1,double> operator() (AutoDiffDiff<1,double> x, AutoDiffDiff<1,double> y) const;
    NGS_DLL_HEADER AutoDiff<1,SIMD<double>> operator() (AutoDiff<1,SIMD<double>> x, AutoDiff<1,SIMD<double>> y) const;
    NGS_DLL_HEADER AutoDiffDiff<1,SIMD<double>> operator() (AutoDiffDiff<1,SIMD<double>> x, AutoDiffDiff<1,SIMD<double>> y) const;

    NGS_DLL_HEADER SIMD<double> operator() (SIMD<double> x, SIMD<double> y) const;

    friend ostream & operator<< (ostream & ost, const BSpline2D & sp);
  };
  extern ostream & operator<< (ostream & ost, const BSpline2D & sp);

}

#endif
