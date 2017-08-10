#ifndef PYTHON_NGSTD_HPP___
#define PYTHON_NGSTD_HPP___
#ifdef NGS_PYTHON

#ifdef __clang__
#pragma clang diagnostic push
// #pragma clang diagnostic ignored "-W#pragma-messages"
#pragma clang diagnostic ignored "-Wunused-local-typedefs"
#pragma clang diagnostic ignored "-Wparentheses-equality"
#endif

#include <pybind11/pybind11.h>
#include <pybind11/eval.h>
#include <pybind11/operators.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#include <ngstd.hpp>
#include <thread>
#include <iostream>

namespace py = pybind11;
using namespace pybind11::literals;

using std::string;
using std::cout;
using std::endl;

using namespace ngstd;

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
  catch (py::error_already_set &e) {
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

struct DummyArgument {};

class PythonEnvironment
{
public:

  PythonEnvironment () { ; }

  virtual ~PythonEnvironment() { }
  
  auto operator[] ( const char *s )
  {
    return py::module::import("__main__").attr(s);
  }
  
  virtual void exec(const string s) 
  {
    try{
      PyRun_SimpleString(s.c_str());
    }
    catch (py::error_already_set const &e) {
      cout << "caught python error: " << e.what() << endl;
      PyErr_Print();
    }
  }

  virtual void exec_file(const string fstr) {
    string output;
    ifstream file;
    file.open(fstr.c_str());
    if (file.is_open())
    {
      while (!file.eof())
        {
          std::string line;
          std::getline(file, line);
          output += line.append("\n");
        }
    }
    file.close();
    exec(output);
  }
  
  

  //  py::exec_file(file.c_str(), main_namespace, main_namespace);
  //  return;
  //  try{
  //    py::exec_file(file.c_str(), main_namespace, main_namespace);
  //  }
  //  catch(py::error_already_set const &) {
  //    PyErr_Print();
  //  }
  //  catch (...) {
  //      cout << "caught!" << endl;
  //  }
  //}

  
};



typedef py::gil_scoped_acquire AcquireGIL; 
typedef py::gil_scoped_release ReleaseGIL; 


inline void InitSlice( const py::slice &inds, size_t len, size_t &start, size_t &step, size_t &n ) {
      size_t stop;
      if (!inds.compute(len, &start, &stop, &step, &n))                                          
        throw py::error_already_set();
}



//////////////////////////////////////////////////////////////////////
// Python class name type traits
template <typename T>
struct PyNameTraits {
  static const string & GetName() { static const string name = typeid(T).name(); return name; }
};


template <typename T>
string GetPyName(const char *prefix = 0) {
  string s;
  if(prefix) s = string(prefix);
  s+= PyNameTraits<T>::GetName();
  return s;
}


template<>
struct PyNameTraits<int> {
  static string GetName() { return "I"; }
};

template<>
struct PyNameTraits<float> {
  static string GetName() { return "F"; }
};

template<>
struct PyNameTraits<double> {
  static string GetName() { return "D"; }
};

template<typename T>
struct PyNameTraits<shared_ptr<T>> {
  static string GetName() { return string("sp_")+GetPyName<T>(); }
};



//////////////////////////////////////////////////////////////////////
template <typename T, typename TCLASS = py::class_<T> >
void PyDefToString( py::module &m, TCLASS &c )
{
    c.def("__str__", &ToString<T>);
    c.def("__repr__", &ToString<T>);
}

template <typename T>
class cl_NonElement 
{
public:
  static T Val() { return 0; }
};

template <typename T>
inline T NonElement() { return cl_NonElement<T>::Val(); }


//////////////////////////////////////////////////////////////////////
// read-only bracket operator
template< typename T, typename TELEM, typename TCLASS = py::class_<T>>
void PyDefROBracketOperator( py::module &m, TCLASS &c )
{
    auto Get = [](T& self, int i) { 
      if( i<self.Size() && i>=0 )
        return self[i];
      throw py::index_error();
      return TELEM();
    }; 
    c.def("__getitem__", Get);
    c.def("Get", Get);
}

// read-write bracket operator
template< typename T, typename TELEM, typename TCLASS = py::class_<T>>
void PyDefBracketOperator( py::module &m, TCLASS &c )
{
    PyDefROBracketOperator<T, TELEM>(m, c);
    auto Set = [](T& self, int i, TELEM val) {
      if( i<self.Size() && i>=0 )
        self[i] = val;
      else
        throw py::index_error();
    };
    c.def("__setitem__", Set);
    c.def("Set", Set);
}


//////////////////////////////////////////////////////////////////////
// Export len, bracket operator and iterator protocoll at once
template <typename T, typename TELEM = double, typename TCLASS = py::class_<T> >
void PyDefVector( py::module &m, TCLASS &c )
{
    c.def("__len__",  []( T& v) { return v.Size();}  );
    c.def("__iter__", [] (T &v)
      { return py::make_iterator(v.begin(), v.end()); },
      py::keep_alive<0,1>()
    );
    PyDefBracketOperator<T, TELEM>(m,c);
}


//////////////////////////////////////////////////////////////////
// SymbolTable - template

template <typename T>
class PyRef 
{
  const T & ref;
public:
  PyRef (const T & aref) : ref(aref) { ; }
  const T & Cast () const { return ref; }
};

template <typename T> 
inline ostream & operator<< (ostream & ost, PyRef<T> ref)
{
  ost << (ref.Cast());
  return ost;
}

template<typename T> struct PyTraits { };
template<> struct PyTraits<double> {typedef py::float_ type;};
template<> struct PyTraits<string> {typedef py::str type;};
template<> struct PyTraits<bool> {typedef py::bool_ type;};
template<> struct PyTraits<int> {typedef py::int_ type;};

template<typename T>
Array<T> makeCArray(const py::tuple & obj)
{     
  Array<T> C_vdL(py::len(obj));   
  for (int i = 0; i < py::len(obj); i++)    
    C_vdL[i] = T(typename PyTraits<T>::type(obj[i]));
  return std::move(C_vdL);
}


template<typename T>
Array<T> makeCArray(const py::list & obj)
{     
  Array<T> C_vdL(py::len(obj));   
  for (int i = 0; i < py::len(obj); i++)    
    C_vdL[i] = T(typename PyTraits<T>::type(obj[i]));        
  return std::move(C_vdL);
}

template<typename T>
Array<T> makeCArray(const py::object & obj)
{     
  if (py::isinstance<py::list>(obj))
    return makeCArray<T>(obj.cast<py::list>());
  if (py::isinstance<py::tuple>(obj))
    return makeCArray<T>(obj.cast<py::tuple>());
  throw py::type_error("Cannot convert Python object to C Array");
}

template<typename T>
Array<T> makeCArraySharedPtr(const py::list & obj)
{
  Array<T> C_vdL(py::len(obj));
  for (int i = 0; i < py::len(obj); i++)
    C_vdL[i] = py::extract<T>(obj[i])();
  return std::move(C_vdL);
}
template<typename T>
Array<decltype(std::declval<T>().Get())> makeCArrayUnpackWrapper(const py::list & obj)
{
  Array<decltype(std::declval<T>().Get())> C_vdL(py::len(obj));
  for (int i = 0; i < py::len(obj); i++)
    C_vdL[i] = py::extract<T>(obj[i])().Get();
  return std::move(C_vdL);
}
template<typename T>
struct PyNameTraits<SymbolTable<T>> {
  static string GetName() { return string("SymbolTable_") + GetPyName<T>(); }
};

template<typename T>
struct PyNameTraits<PyRef<T>> {
  static string GetName() { return string("Ref_") + GetPyName<T>(); }
};

template <typename T, typename PY_T = T>
void PyExportSymbolTable (py::module &m)
{
  typedef SymbolTable<T> ST;
  
  string name = GetPyName<ST>();
  py::class_<ST>(m, name.c_str())
    .def("__str__", &ToString<ST>)
    .def("__len__", &ST::Size)
    .def("__contains__", &ST::Used)
    .def("GetName", [](ST & self, int i) { return string(self.GetName(i)); })
    .def("__getitem__", [](ST & self, string name) -> PY_T
                                        {
                                          if (!self.Used(name)) throw py::index_error();
                                          return self[name]; 
                                        })
    .def("__getitem__", [](ST & self, int i) -> PY_T
                                         {
                                           if (i < 0 || i >= self.Size()) throw py::index_error();
                                           return self[i];  
                                         })
    ;
}  


// conversion not possible for shared_ptr<double>, so we have a special treatment:
template <> inline void PyExportSymbolTable<shared_ptr<double>, shared_ptr<double>> (py::module &m)
{
  typedef SymbolTable<shared_ptr<double>> ST;
  
  string name = GetPyName<ST>();
  py::class_<ST>(m, name.c_str())
    .def("__str__", &ToString<ST>)
    .def("__len__", &ST::Size)
    .def("__contains__", &ST::Used)
    .def("GetName", [](ST & self, int i) { return string(self.GetName(i)); })
    .def("__getitem__", [](ST & self, string name)
                                        {
                                          if (!self.Used(name)) throw py::index_error();
                                          return *self[name]; 
                                        })
    .def("__getitem__", [](ST & self, int i)
                                         {
                                           if (i < 0 || i >= self.Size()) throw py::index_error();
                                           return *self[i];  
                                         })
    ;
}

// Parse python kwargs to flags
Flags CreateFlagsFromKwArgs(const py::object& pyclass, const py::kwargs& kwargs, py::list info = py::list());

// replace docu links with plain text for help function
const char* docu_string(const char* str);

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

#endif // NGS_PYTHON
#endif // PYTHON_NGSTD_HPP___
