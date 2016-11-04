#ifdef NG_PYTHON

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
namespace py = pybind11;
#include <iostream>
#include <sstream>

namespace pybind11 {
template<typename T>
bool CheckCast( py::handle obj ) {
  try{
    obj.cast<T>();
    return true;
  }
  catch (py::cast_error &e) {
    return false;
  }
}


template <typename T>
struct extract
{
  py::handle obj;
  extract( py::handle aobj ) : obj(aobj) {}

  bool check() { return CheckCast<T>(obj); }
  T operator()() { return obj.cast<T>(); }
};
}

struct NGDummyArgument {};

inline void NOOP_Deleter(void *) { ; }

namespace netgen
{

  //////////////////////////////////////////////////////////////////////
  // Lambda to function pointer conversion
  template <typename Function>
  struct function_traits
    : public function_traits<decltype(&Function::operator())> {};

  template <typename ClassType, typename ReturnType, typename... Args>
  struct function_traits<ReturnType(ClassType::*)(Args...) const> {
    typedef ReturnType (*pointer)(Args...);
    typedef ReturnType return_type;
  };

  template <typename Function>
  typename function_traits<Function>::pointer
  FunctionPointer (const Function& lambda) {
    return static_cast<typename function_traits<Function>::pointer>(lambda);
  }


  template <class T>
  inline std::string ToString (const T& t)
  {
    std::stringstream ss;
    ss << t;
    return ss.str();
  }

}

#endif

