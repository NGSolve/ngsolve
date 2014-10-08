#ifdef NG_PYTHON

#include <boost/python.hpp>
namespace bp = boost::python;
#include <iostream>

namespace netgen
{

	class ModuleScope  {
		bp::scope *local_scope;
	public:
		ModuleScope(const std::string name) : local_scope(nullptr) {
			std::string nested_name = name;
			if (bp::scope())
				nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + "." + name);

			bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

			std::cout << "exporting " << nested_name << std::endl;
			bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
			parent.attr(name.c_str()) = module;

			local_scope = new bp::scope(module);
		}

		~ModuleScope() {
			if (local_scope)
				delete (local_scope);
		}
		
	};

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

