#ifndef FILE_COEFFICIENT_STDMATH
#define FILE_COEFFICIENT_STDMATH


namespace ngfem
{

  struct GenericSqrt {
    template <typename T> T operator() (T x) const { return sqrt(x); }
    static string Name() { return "sqrt"; }
    void DoArchive(Archive& ar) {}
  };
  
  struct GenericSin {
    template <typename T> T operator() (T x) const { return sin(x); }
    static string Name() { return "sin"; }
    void DoArchive(Archive& ar) {}
  };

  struct GenericCos {
    template <typename T> T operator() (T x) const { return cos(x); }
    static string Name() { return "cos"; }
    void DoArchive(Archive& ar) {}
  };

  struct GenericTan {
    template <typename T> T operator() (T x) const { return tan(x); }
    static string Name() { return "tan"; }
    void DoArchive(Archive& ar) {}
  };

  struct GenericATan {
    template <typename T> T operator() (T x) const { return atan(x); }
    static string Name() { return "atan"; }
    void DoArchive(Archive& ar) {}
  };
  
  struct GenericACos {
    template <typename T> T operator() (T x) const { return acos(x); }
    static string Name() { return "acos"; }
    void DoArchive(Archive& ar) {}
  };
  
  struct GenericASin {
    template <typename T> T operator() (T x) const { return asin(x); }
    static string Name() { return "asin"; }
    void DoArchive(Archive& ar) {}
  };


  

  struct GenericSinh {
    template <typename T> T operator() (T x) const { return sinh(x); }
    static string Name() { return "sinh"; }
    void DoArchive(Archive& ar) {}
  };
  
  struct GenericCosh {
    template <typename T> T operator() (T x) const { return cosh(x); }
    static string Name() { return "cosh"; }
    void DoArchive(Archive& ar) {}
  };

  

  struct GenericExp {
    template <typename T> T operator() (T x) const { return exp(x); }
    static string Name() { return "exp"; }
    void DoArchive(Archive& ar) {}
  };
  
  struct GenericLog {
    template <typename T> T operator() (T x) const { return log(x); }
    static string Name() { return "log"; }
    void DoArchive(Archive& ar) {}
  };


  struct GenericErf {
    template <typename T> T operator() (T x) const { return erf(x); }
    Complex operator() (Complex x) const { throw Exception("no erf for Complex"); }
    SIMD<Complex> operator() (SIMD<Complex> x) const { throw ExceptionNOSIMD("no erf for simd(complex)"); }  
    static string Name() { return "erf"; }
    void DoArchive(Archive& ar) {}
  };



struct GenericFloor {
  template <typename T> T operator() (T x) const { return floor(x); }
  Complex operator() (Complex x) const { throw Exception("no floor for Complex"); }  
  // SIMD<double> operator() (SIMD<double> x) const { throw ExceptionNOSIMD("no floor for simd"); }
  SIMD<Complex> operator() (SIMD<Complex> x) const { throw ExceptionNOSIMD("no floor for simd"); }  
  // AutoDiff<1> operator() (AutoDiff<1> x) const { throw Exception("no floor for AD"); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception("no floor for ADD"); }
  static string Name() { return "floor"; }
  void DoArchive(Archive& ar) {}
};
struct GenericCeil {
  template <typename T> T operator() (T x) const { return ceil(x); }
  Complex operator() (Complex x) const { throw Exception("no ceil for Complex"); }  
  // SIMD<double> operator() (SIMD<double> x) const { throw ExceptionNOSIMD("no ceil for simd"); }
  SIMD<Complex> operator() (SIMD<Complex> x) const { throw ExceptionNOSIMD("no ceil for simd"); }  
  // AutoDiff<1> operator() (AutoDiff<1> x) const { throw Exception("no ceil for AD"); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception("no ceil for ADD"); }
  static string Name() { return "ceil"; }
  void DoArchive(Archive& ar) {}
};


  
  using std::sqrt;
  shared_ptr<CoefficientFunction> sqrt(shared_ptr<CoefficientFunction> x);  

  using std::sin;
  shared_ptr<CoefficientFunction> sin(shared_ptr<CoefficientFunction> x);

  using std::cos;
  shared_ptr<CoefficientFunction> cos(shared_ptr<CoefficientFunction> x);
  
  using std::tan;
  shared_ptr<CoefficientFunction> tan(shared_ptr<CoefficientFunction> x);  

  using std::asin;  
  shared_ptr<CoefficientFunction> asin(shared_ptr<CoefficientFunction> x);

  using std::acos;  
  shared_ptr<CoefficientFunction> acos(shared_ptr<CoefficientFunction> x);
  
  using std::atan;
  shared_ptr<CoefficientFunction> atan(shared_ptr<CoefficientFunction> x);  

  using std::sinh;    
  shared_ptr<CoefficientFunction> sinh(shared_ptr<CoefficientFunction> x);
  
  using std::cosh;
  shared_ptr<CoefficientFunction> cosh(shared_ptr<CoefficientFunction> x);  

  using std::exp;
  shared_ptr<CoefficientFunction> exp(shared_ptr<CoefficientFunction> x);
  
  using std::log;
  shared_ptr<CoefficientFunction> log(shared_ptr<CoefficientFunction> x);  

  using std::erf;
  shared_ptr<CoefficientFunction> erf(shared_ptr<CoefficientFunction> x);  

  using std::floor;
  shared_ptr<CoefficientFunction> floor(shared_ptr<CoefficientFunction> x);

  using std::ceil;        
  shared_ptr<CoefficientFunction> ceil(shared_ptr<CoefficientFunction> x);  
}


#endif
