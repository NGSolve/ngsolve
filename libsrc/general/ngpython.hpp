#ifdef NG_PYTHON

#include <boost/python.hpp>
namespace bp = boost::python;


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

