#ifndef PYTHON_NGSTD_HPP___
#define PYTHON_NGSTD_HPP___
#ifdef NGS_PYTHON


#include <boost/python.hpp>
#include <ngstd.hpp>
#include <thread>
#include <iostream>

namespace bp = boost::python;

using std::string;
using std::cout;
using std::endl;

using namespace ngstd;



class BasePythonEnvironment
{
  // protected:
public:
  BasePythonEnvironment () { ; }

  bp::object main_module; 
  bp::object main_namespace; 

public:
  virtual ~BasePythonEnvironment() { }
  
  auto operator[] ( const char *s ) -> decltype(main_namespace[s]) {
    return main_namespace[s];
  }
  

  virtual void exec(const string s) {
    bp::exec(s.c_str(), main_namespace, main_namespace);
  }

  virtual void exec_file(const char *file) {
    try{
      bp::exec_file(file, main_namespace, main_namespace);
    }
    catch(bp::error_already_set const &) {
      PyErr_Print();
    }
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
};

template <typename Function>
typename function_traits<Function>::pointer
FunctionPointer (const Function& lambda) {
  return static_cast<typename function_traits<Function>::pointer>(lambda);
}


//////////////////////////////////////////////////////////////////////
// Python class name type traits
template <typename T>
struct PyNameTraits {
  static const string & GetName() { static const string name = typeid(T).name(); return name; }
};

template<>
struct PyNameTraits<int> {
  static const string & GetName() { static const string name = "I"; return name; }
};

template<>
struct PyNameTraits<float> {
  static const string & GetName() { static const string name = "F"; return name; }
};

template<>
struct PyNameTraits<double> {
  static const string & GetName() { static const string name = "D"; return name; }
};

template <typename T>
string GetPyName(const char *prefix = 0) {
  string s;
  if(prefix) s = string(prefix);
  s+= PyNameTraits<T>::GetName();
  return s;
}

//////////////////////////////////////////////////////////////////////
template< typename T>
struct PyDefToString : public boost::python::def_visitor<PyDefToString<T> > {
  template <class Tclass>
  void visit(Tclass& c) const {
    c.def("__str__", &PyDefToString<T>::ToString);
    c.def("__repr__", &PyDefToString<T>::ToString);
  }

  static string ToString(T &t ) {
    std::ostringstream s;
    s << t;
    return s.str();
  }
};


//////////////////////////////////////////////////////////////////////
// read-only bracket operator  (Matthias, please polish !)
// enough for iterator
// read-wirt bracket could use inheritance 
template< typename T, typename TELEM = double >
struct PyDefROBracketOperator : public boost::python::def_visitor<PyDefROBracketOperator<T,TELEM> > {
  template <class Tclass>
  void visit(Tclass& c) const {
    c
      .def("__getitem__", &PyDefROBracketOperator<T,TELEM>::Get)
      .def("Get", &PyDefROBracketOperator<T,TELEM>::Get)
      ;
  }

  static TELEM Get(T& self, int i) { 
    if( i<self.Size() && i>=0 )
      return self[i];
    RaiseIndexError();
    // return TELEM();  // problem with deleted default-constructor
  } 

  static void RaiseIndexError() {
    // PythonEnvironment::getInstance().exec("raise IndexError()\n");
    cerr << "python Index error" << endl;
  }

};


template< typename T, typename TELEM = double >
struct PyDefBracketOperator : public boost::python::def_visitor<PyDefBracketOperator<T,TELEM> > {
  template <class Tclass>
  void visit(Tclass& c) const {
    c
      .def("__getitem__", &PyDefBracketOperator<T,TELEM>::Get)
      .def("__setitem__", &PyDefBracketOperator<T,TELEM>::Set)
      .def("Get", &PyDefBracketOperator<T,TELEM>::Get)
      .def("Set", &PyDefBracketOperator<T,TELEM>::Set)
      ;
  }

  static TELEM Get(T& self, int i) { 
    if( i<self.Size() && i>=0 )
      return self[i];
    RaiseIndexError();
    return TELEM();
  } 

  static void Set(T& self, int i, TELEM val) {
    if( i<self.Size() && i>=0 )
      self[i] = val;
    else
      RaiseIndexError();
  }

  static void RaiseIndexError() {
    // PythonEnvironment::getInstance().exec("raise IndexError()\n");
    cerr << "python Index error" << endl;
  }

};

//////////////////////////////////////////////////////////////////////
// Python iterator protocoll
template <typename T, typename TELEM>
class PyIterator {
  T &v;
  int size;
  int index;
  int startindex;
  //     typedef typename std::remove_reference<decltype(v[startindex])>::type TELEM;
    
public:
  PyIterator(T &v_, int size_, int startindex_ = 0) : v(v_), size(size_), index(startindex_), startindex(startindex_) {}

  TELEM Next() { 
    if(index<startindex+size) return v[index++];
    else
      bp::exec("raise StopIteration()\n");
    // cerr << "python Index error" << endl;
    // return TELEM();
  }

  static void Export () {
    bp::class_<PyIterator<T, TELEM> >("PyIterator", bp::no_init).def("__next__", &PyIterator<T, TELEM>::Next);
  }
};




//////////////////////////////////////////////////////////////////////
// Export len, bracket operator and iterator protocoll at once
template< typename T,  typename TELEM = double>
struct PyDefVector : public boost::python::def_visitor<PyDefVector<T,TELEM> > {
  template <class Tclass>
  void visit(Tclass& c) const {
    PyIterator<T, TELEM>::Export();
    c
      .def("__len__", FunctionPointer( []( T& v) { return v.Size();} ) )
      .def(PyDefBracketOperator<T, TELEM>())
      .def("__iter__", FunctionPointer([](T &v) { return PyIterator<T, TELEM>(v, v.Size(), 0 ); }))
      ;
  }
};

// Joachim: Matthias, please polish
template< typename T,  typename TELEM = double>
struct PyDefList : public boost::python::def_visitor<PyDefList<T,TELEM> > {
  template <class Tclass>
  void visit(Tclass& c) const {
    PyIterator<T, TELEM>::Export();
    c
      .def(PyDefROBracketOperator<T, TELEM>())
      .def("__iter__", FunctionPointer([](T &v) { return PyIterator<T, TELEM>(v, v.Size(), 0 ); }))
      ;
  }
};


//////////////////////////////////////////////////////////////////////
// Enable numeric expressions for matrix class
inline void PyEnableMatExpr(BasePythonEnvironment & py_env, const char *class_name) {
  // PythonEnvironment &py_env = PythonEnvironment::getInstance();

  string cn(class_name);
  py_env.exec(cn + ".expr = property(lambda self: MatExpr(self))");
  py_env.exec(cn + ".data = property(lambda self: None, lambda self, a: Expr(a).AssignTo(self.expr))");
  py_env.exec(cn + ".__add__ = lambda self,y: self.expr+Expr(y)");
  py_env.exec(cn + ".__rmul__ = lambda self,y: y*self.expr ");
  py_env.exec(cn + ".__mul__ = lambda self,y: MatVecExpr(self.expr, Expr(y))");
}

//////////////////////////////////////////////////////////////////////
// Enable numeric expressions for vector class
inline void PyEnableVecExpr(BasePythonEnvironment & py_env, const char *class_name) {
  // PythonEnvironment &py_env = PythonEnvironment::getInstance();

  string cn(class_name);
  py_env.exec(cn + ".expr = property(lambda self: VecExpr(self))");
  py_env.exec(cn + ".data = property(lambda self: None, lambda self, a: Expr(a).AssignTo(self.expr))");
  py_env.exec(cn + ".__add__ = lambda self,y: self.expr+Expr(y)");
  py_env.exec(cn + ".__rmul__ = lambda self,y: y*self.expr ");
  py_env.exec(cn + ".__getitem__ = GetSlice");
  py_env.exec(cn + ".__setitem__ = SetSlice");
}

//////////////////////////////////////////////////////////////////////
// Enable Slicing support
inline void PyEnableSlicing(BasePythonEnvironment & py_env, const char *class_name) {
  // PythonEnvironment &py_env = PythonEnvironment::getInstance();

  string cn(class_name);
  py_env.exec(cn + ".__getitem__ = GetSlice");
  py_env.exec(cn + ".__setitem__ = SetSlice");
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



template<typename T>
struct PyNameTraits<SymbolTable<T>> {
  static string GetName()
  { return string("SymbolTable_") + GetPyName<T>(); }
};


template<typename T>
struct PyNameTraits<PyRef<T>> {
  static string GetName()
  { return string("Ref_") + GetPyName<T>(); }
};


template <typename T> struct T_MyRemovePtr 
{ typedef T type; typedef T & ref_type; };

template <typename T> struct T_MyRemovePtr<shared_ptr<T>> 
{ typedef T type; typedef T& ref_type; };

template <class T>
struct cl_remove_pointer
{ static T Val (T v) { return v; } };
template <class T>
struct cl_remove_pointer<T*>
{ static T & Val (T * v) { return *v; }};

template <class T> inline auto MyRemovePtr (T d) -> decltype(cl_remove_pointer<T>::Val(d))
{ return cl_remove_pointer<T>::Val(d); }




template <typename T> void PyExportSymbolTable ()
{
  string name = GetPyName<SymbolTable<T>>();
  bp::class_<SymbolTable<T>>(name.c_str())
    .add_property ("size", &SymbolTable<T>::Size)
    .def(PyDefToString<SymbolTable<T>>())

    .def("__getitem__", FunctionPointer([](SymbolTable<T> & self, bp::str s)
                                        {
                                          string name = bp::extract<string>(s);
                                          if (!self.Used(name))
                                            cerr << "unused" << endl;
                                          return self[name];  // needs boost 1.55 !!
                                        })
         )
    ;



  name = GetPyName<PyRef<SymbolTable<T>>>();
  bp::class_<PyRef<SymbolTable<T>>>(name.c_str(), bp::no_init)
    .add_property ("size", FunctionPointer([](PyRef<SymbolTable<T>> & self)
                                           { return self.Cast().Size(); }))
    .def(PyDefToString<PyRef<SymbolTable<T>>>())
    .def("__getitem__", FunctionPointer([](PyRef<SymbolTable<T>> & self, bp::str s)
                                        -> typename T_MyRemovePtr<T>::ref_type
                                        {
                                          if (!self.Cast().Used(bp::extract<string>(s)))
                                            cerr << "unused" << endl;
                                          return *self.Cast()[bp::extract<string>(s)];
                                        })
         , bp::return_value_policy<bp::reference_existing_object>()
         )
    .def("__len__", FunctionPointer( [](PyRef<SymbolTable<T>> & self) 
                                     { return self.Cast().Size();} ) )
    .def("__getitem__", FunctionPointer([](PyRef<SymbolTable<T>> & self, int nr)
                                        -> typename T_MyRemovePtr<T>::ref_type
                                        {
                                          if (nr < 0 || nr >= self.Cast().Size())
                                            cerr << "unused" << endl;
                                          return *self.Cast()[nr];
                                        })
         , bp::return_value_policy<bp::reference_existing_object>())
    ;
    }



template <typename T> void PyExportSymbolTableStdTypes ()
{
  string name = GetPyName<SymbolTable<T>>();

  bp::class_<SymbolTable<T>>(name.c_str())
    .add_property ("size", &SymbolTable<T>::Size)
    .def(PyDefToString<SymbolTable<T>>())
    .def("__getitem__", FunctionPointer([](SymbolTable<T> & self, bp::str s)
                                        {
                                          if (!self.Used(bp::extract<string>(s)))
                                            cerr << "unused" << endl;
                                          return MyRemovePtr (self[bp::extract<string>(s)]);
                                        }))
    ;
}









#endif // NGS_PYTHON
#endif // PYTHON_NGSTD_HPP___
