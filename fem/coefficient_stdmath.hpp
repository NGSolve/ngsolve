#ifndef FILE_COEFFICIENT_STDMATH
#define FILE_COEFFICIENT_STDMATH


namespace ngfem
{
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
  
  
  
  
  shared_ptr<CoefficientFunction> ExpCF(shared_ptr<CoefficientFunction> x);  
  shared_ptr<CoefficientFunction> LogCF(shared_ptr<CoefficientFunction> x);  
  shared_ptr<CoefficientFunction> SinCF(shared_ptr<CoefficientFunction> x);  
  shared_ptr<CoefficientFunction> CosCF(shared_ptr<CoefficientFunction> x);  
  
}


#endif
