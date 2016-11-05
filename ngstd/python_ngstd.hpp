#ifndef PYTHON_NGSTD_HPP___
#define PYTHON_NGSTD_HPP___
#ifdef NGS_PYTHON

#pragma clang diagnostic push
// #pragma clang diagnostic ignored "-W#pragma-messages"
#pragma clang diagnostic ignored "-Wunused-local-typedefs"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/complex.h>

#pragma clang diagnostic pop

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
  py::object main_module; 
  py::object main_namespace; 

public:

  PythonEnvironment () { ; }
  PythonEnvironment (py::object _module)
    : main_module(_module), 
      main_namespace(main_module.attr("__dict__")) 
  { ; }

  virtual ~PythonEnvironment() { }
  
  auto operator[] ( const char *s ) -> decltype(main_namespace[s]) 
  {
    return main_namespace[s];
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



extern PythonEnvironment pyenv;

typedef py::gil_scoped_acquire AcquireGIL; 
typedef py::gil_scoped_release ReleaseGIL; 


static void InitSlice( const py::slice &inds, size_t len, size_t &start, size_t &step, size_t &n ) {
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
// Python iterator protocoll

template <typename T, typename TELEM>
class PyIterator
{
  T & v;
  int size;
  int index;
  int startindex;
    
public:
  PyIterator(T & v_, int size_, int startindex_ = 0)
    : v(v_), size(size_), index(startindex_), startindex(startindex_)
    {}

  TELEM Next()
  { 
    if(index<startindex+size) return v[index++];
    else 
      throw py::stop_iteration();
    return NonElement<TELEM>();
  }

  static void Export (py::module &m)
  {
    string name = string("PyIterator")+GetPyName<T>();
    py::class_<PyIterator<T, TELEM> >(m, name.c_str()).def("__next__", &PyIterator<T, TELEM>::Next);
  }
};

//////////////////////////////////////////////////////////////////////
// Export len, bracket operator and iterator protocoll at once
template <typename T, typename TELEM = double, typename TCLASS = py::class_<T> >
void PyDefVector( py::module &m, TCLASS &c )
{
    PyIterator<T, TELEM>::Export(m);

    c.def("__len__",  []( T& v) { return v.Size();}  );
    c.def("__iter__", [](T &v) { return PyIterator<T, TELEM>(v, v.Size(), 0 ); });
    PyDefBracketOperator<T, TELEM>(m,c);
}

template <typename T, typename TELEM = double, typename TCLASS = py::class_<T> >
void PyDefIterable( py::module &m, TCLASS &c )
{
    PyIterator<T, TELEM>::Export();
    c.def(PyDefROBracketOperator<T, TELEM>());
    c.def("__iter__", [](T &v) 
                                       {
                                         return PyIterator<T, TELEM>(v, v.Size(), 0 ); 
                                       });
}

template <typename T>
class PyIterator2
{
  typedef decltype(GetReturnValue(&T::begin)) TITER;
  typedef decltype(GetReturnValue(&TITER::operator*)) TELEM;
  TITER begin, end;
public:
  PyIterator2 (const T & container)
    : begin(container.begin()), end(container.end()) { ; }
  TELEM Next()
  {
    if (! (begin != end))
      throw py::stop_iteration();
      
    auto tmp = begin;
    ++begin;
    return *tmp;
  }
};

// iterable where elements have an increment operator (JS)
template <typename T, typename TCLASS = py::class_<T> >
void PyDefIterable2( py::module &m, TCLASS &c )
{
    string itername = string("PyIterator2_")+GetPyName<T>();
    py::class_<PyIterator2<T>> citer(m, itername.c_str());
    citer.def("__next__", &PyIterator2<T>::Next);
    c.def("__iter__", [](const T & c) { return PyIterator2<T>(c); });
}


// the iterator object copies the container (= range_expr)
// otherwise, it might be destroyed by python too early
template <typename T, typename TCLASS = py::class_<T> >
void PyDefIterable3( py::module &m, TCLASS &c )
{

  typedef decltype(GetReturnValue(&T::begin)) TITER;
  typedef decltype(GetReturnValue(&TITER::operator*)) TELEM;

  class Iterator
  {
    shared_ptr<T> cont2;
    TITER begin, end;
  public:
    Iterator (shared_ptr<T> container)
      : cont2(container), begin(cont2->begin()), end(cont2->end()) 
    { 
      cout << "Iterator from container" << endl;
    }
    
    /*
    Iterator (T && container)
      : cont2(move(container)), begin(cont2->begin()), end(cont2->end())
    { 
      cout << "Iterator from rvalue - container" << endl;
    }
    */

    Iterator (const Iterator & it2)
      : cont2(it2.cont2), begin(cont2->begin()), end(cont2->end()) 
    {
      cout << "copy iterator" << endl; 
    }

    /*
    Iterator (Iterator && it2)
      : cont2(move(it2.cont2)), begin(cont2.begin()), end(cont2.end()) { ; }
    */

    TELEM Next()
    {
      if (! (begin != end))
        throw py::stop_iteration();
        
      auto tmp = begin;
      ++begin;
      return *tmp;
    }
  };

    string itername = string("PyIterator3_")+GetPyName<T>();
    py::class_<Iterator>(m, itername.c_str())
      .def("__next__", &Iterator::Next)
      ;
    c.def("__iter__",
          [](shared_ptr<T> c) 
           { 
             cout << "create python iterator" << endl;
             return Iterator(c); 
           })
      ;
  }

//////////////////////////////////////////////////////////////////////
// Enable numeric expressions for matrix class
inline void PyEnableMatExpr(const char *class_name) {
  string cn(class_name);
  pyenv.exec(cn + ".expr = property(lambda self: MatExpr(self))");
  pyenv.exec(cn + ".data = property(lambda self: None, lambda self, a: Expr(a).AssignTo(self.expr))");
  pyenv.exec(cn + ".__add__ = lambda self,y: self.expr+Expr(y)");
  pyenv.exec(cn + ".__rmul__ = lambda self,y: y*self.expr ");
  pyenv.exec(cn + ".__mul__ = lambda self,y: MatVecExpr(self.expr, Expr(y))");
}

//////////////////////////////////////////////////////////////////////
// Enable numeric expressions for vector class
inline void PyEnableVecExpr(const char *class_name) {

  string cn(class_name);
  pyenv.exec(cn + ".expr = property(lambda self: VecExpr(self))");
  pyenv.exec(cn + ".data = property(lambda self: None, lambda self, a: Expr(a).AssignTo(self.expr))");
  pyenv.exec(cn + ".__add__ = lambda self,y: self.expr+Expr(y)");
  pyenv.exec(cn + ".__rmul__ = lambda self,y: y*self.expr ");
  pyenv.exec(cn + ".__getitem__ = GetSlice");
  pyenv.exec(cn + ".__setitem__ = SetSlice");
}

//////////////////////////////////////////////////////////////////////
// Enable Slicing support
inline void PyEnableSlicing(const char *class_name) {

  string cn(class_name);
  pyenv.exec(cn + ".__getitem__ = GetSlice");
  pyenv.exec(cn + ".__setitem__ = SetSlice");
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

template<typename T> struct PyTraits {};
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
  auto list = py::extract<py::list>(obj);
  if (list.check())
    return makeCArray<T>(list());
  auto tuple = py::extract<py::tuple>(obj);
  if (tuple.check())
    return makeCArray<T>(tuple());
  throw py::type_error("Cannot convert Python object to C Array");
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

template <typename T>
void PyExportSymbolTable (py::module &m)
{
  typedef SymbolTable<T> ST;
  
  string name = GetPyName<ST>();
  py::class_<ST>(m, name.c_str())
    .def("__str__", &ToString<ST>)
    .def("__len__", &ST::Size)
    .def("__contains__", &ST::Used)
    .def("GetName", [](ST & self, int i) { return string(self.GetName(i)); })
    .def("__getitem__", [](ST & self, string name)
                                        {
                                          if (!self.Used(name)) throw py::index_error();
                                          return self[name]; 
                                        })
    .def("__getitem__", [](ST & self, int i)
                                         {
                                           if (i < 0 || i >= self.Size()) throw py::index_error();
                                           return self[i];  
                                         })
    ;
}  


// convertion not possible for shared_ptr<double>, so we have a special treatment:
template <> inline void PyExportSymbolTable<shared_ptr<double>> (py::module &m)
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

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

#endif // NGS_PYTHON
#endif // PYTHON_NGSTD_HPP___
