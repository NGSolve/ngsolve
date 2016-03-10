#ifndef PYTHON_NGSTD_HPP___
#define PYTHON_NGSTD_HPP___
#ifdef NGS_PYTHON

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#pragma-messages"

#include <boost/python.hpp>

#pragma clang diagnostic pop

#include <ngstd.hpp>
#include <thread>
#include <iostream>


#if BOOST_VERSION >= 106000
  // Boost Python 1.60 does not automatically register shared_ptr<T> with T
  #define REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(type) bp::register_ptr_to_python<type>();
#else
  #define REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(type)
#endif

namespace bp = boost::python;

using std::string;
using std::cout;
using std::endl;

using namespace ngstd;



class PythonEnvironment
{
  bp::object main_module; 
  bp::object main_namespace; 

public:

  PythonEnvironment () { ; }
  PythonEnvironment (bp::object _module)
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
      bp::exec(s.c_str(), main_namespace, main_namespace);
    }
    catch (bp::error_already_set const &) {
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
  
  

  //  bp::exec_file(file.c_str(), main_namespace, main_namespace);
  //  return;
  //  try{
  //    bp::exec_file(file.c_str(), main_namespace, main_namespace);
  //  }
  //  catch(bp::error_already_set const &) {
  //    PyErr_Print();
  //  }
  //  catch (...) {
  //      cout << "caught!" << endl;
  //  }
  //}

  
};



extern PythonEnvironment pyenv;

class AcquireGIL 
{
public:
  inline AcquireGIL(){
    state = PyGILState_Ensure();
  }
  
  inline ~AcquireGIL(){
    PyGILState_Release(state);
  }
private:
  PyGILState_STATE state;
};

class ReleaseGIL
{
public:
  inline ReleaseGIL() {
    m_thread_state = PyEval_SaveThread();
  }
  
  inline ~ReleaseGIL() {
    PyEval_RestoreThread(m_thread_state);
    m_thread_state = NULL;
  }
  
private:
  PyThreadState * m_thread_state;
};












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

template <typename T>
class cl_NonElement 
{
public:
  static T Val() { return 0; }
};

template <typename T>
inline T NonElement() { return cl_NonElement<T>::Val(); }


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

  static TELEM Get(T& self, int i) 
  { 
    if( i<self.Size() && i>=0 )
      return self[i];
    RaiseIndexError();
    return NonElement<TELEM>();
  } 

  static void RaiseIndexError() {
    // PythonEnvironment::getInstance().exec("raise IndexError()\n");
    bp::exec("raise IndexError()\n");
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
  T & v;
  int size;
  int index;
  int startindex;
  //     typedef typename std::remove_reference<decltype(v[startindex])>::type TELEM;
    
public:
  PyIterator(T & v_, int size_, int startindex_ = 0) : v(v_), size(size_), index(startindex_), startindex(startindex_) {}

  TELEM Next() { 
    if(index<startindex+size) return v[index++];
    else 
      bp::exec("raise StopIteration()\n");
    // PyErr_SetNone(PyExc_StopIteration);
    return NonElement<TELEM>();
    // cerr << "python Index error" << endl;
    // return TELEM();
  }

  static void Export () {
    string name = string("PyIterator")+GetPyName<T>();
    bp::class_<PyIterator<T, TELEM> >( name.c_str(), bp::no_init).def("__next__", &PyIterator<T, TELEM>::Next);
    // bp::class_<PyIterator<T, TELEM> >(string("PyIterator")+GetPyName<T>(), bp::no_init).def("__next__", &PyIterator<T, TELEM>::Next);
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
struct PyDefIterable : public boost::python::def_visitor<PyDefIterable<T,TELEM> > {
  template <class Tclass>
  void visit(Tclass& c) const {
    PyIterator<T, TELEM>::Export();
    c
      .def(PyDefROBracketOperator<T, TELEM>())
      .def("__iter__", FunctionPointer([](T &v) 
                                       {
                                         return PyIterator<T, TELEM>(v, v.Size(), 0 ); 
                                       }))
      ;
  }
};



// iterable where elements have an increment operator (JS)
template <typename T>
class PyDefIterable2 : public bp::def_visitor<PyDefIterable2<T>>
{
  // typedef decltype(&T::First) TFUNCTION;
  // typedef typename function_traits<TFUNCTION>::return_type TELEM;

  /*
  typedef decltype(GetReturnValue(&T::First)) TELEM;
  class Iterator
  {
    TELEM first;
    int size, cnt;
  public:
    Iterator (const T & container)
    // : first (container.First()), size(container.Size()), cnt(0) { }
      : first (*container.begin()), size(container.Size()), cnt(0) { }
    TELEM Next()
    {
      if (++cnt > size)
        bp::exec("raise StopIteration()\n");
      return first++;
    }
  };
  */

  /*
  typedef decltype(GetReturnValue(&T::begin)) TITER;
  typedef decltype(GetReturnValue(&TITER::operator*)) TELEM;

  class Iterator
  {
    TITER first;
    int size, cnt;
  public:
    Iterator (const T & container)
    // : first (container.First()), size(container.Size()), cnt(0) { }
      : first (container.begin()), size(container.Size()), cnt(0) { }
    TELEM Next()
    {
      if (++cnt > size)
        bp::exec("raise StopIteration()\n");
      auto tmp = first;
      ++first;
      return *tmp;
    }
  };
  */

  typedef decltype(GetReturnValue(&T::begin)) TITER;
  typedef decltype(GetReturnValue(&TITER::operator*)) TELEM;

  class Iterator
  {
    TITER begin, end;
  public:
    Iterator (const T & container)
      : begin(container.begin()), end(container.end()) { ; }
    TELEM Next()
    {
      if (! (begin != end))
        bp::exec("raise StopIteration()\n");
        
      auto tmp = begin;
      ++begin;
      return *tmp;
    }
  };


public:
  template <class Tclass>
  void visit (Tclass & c) const
  {
    string itername = string("PyIterator2_")+GetPyName<T>();
    bp::class_<Iterator>(itername.c_str(),bp::no_init)
      .def("__next__", &Iterator::Next)
      ;
    c.def("__iter__", FunctionPointer
          ([](const T & c) { return Iterator(c); }))
      ;
  }
};


// the iterator object copies the container (= range_expr)
// otherwise, it might be destroyed by python too early
template <typename T>
class PyDefIterable3 : public bp::def_visitor<PyDefIterable3<T>>
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
        bp::exec("raise StopIteration()\n");
        
      auto tmp = begin;
      ++begin;
      return *tmp;
    }
  };

public:
  template <class Tclass>
  void visit (Tclass & c) const
  {
    string itername = string("PyIterator3_")+GetPyName<T>();
    bp::class_<Iterator>(itername.c_str(),bp::no_init)
      .def("__next__", &Iterator::Next)
      ;
    c.def("__iter__", FunctionPointer
          ([](shared_ptr<T> c) 
           { 
             cout << "create python iterator" << endl;
             return Iterator(c); 
           }))
      ;
  }
};







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


/*
template<typename T>
Array<T> & makeCArray(const bp::object & obj)
{     
  Array<T> * C_vdL = new Array<T>(bp::len(obj));    
  for (int i = 0; i < bp::len(obj); i++)    
    (*C_vdL)[i] = bp::extract<T>(obj[i]);        
  return *C_vdL;
}
*/
template<typename T>
Array<T> makeCArray(const bp::object & obj)
{     
  Array<T> C_vdL(bp::len(obj));   
  for (int i = 0; i < bp::len(obj); i++)    
    C_vdL[i] = bp::extract<T>(obj[i]);        
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

/*
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
*/


template <typename T>
struct PythonDictFromSymbolTable {
  static PyObject* convert(const SymbolTable<T> & st)
    {
      bp::dict res;
      for(int i = 0; i < st.Size(); i++) 
        res[st.GetName(i)] = st[i];
      return bp::incref(res.ptr());
    }
};

template <typename T> void PyExportSymbolTable ()
{
  boost::python::to_python_converter< SymbolTable<T>, PythonDictFromSymbolTable<T> >();
}

// convertion not possible for shared_ptr<double>, so we have a special treatment:
template <> inline void PyExportSymbolTable<shared_ptr<double>> ()
{
  typedef SymbolTable<shared_ptr<double>> ST;
  
  string name = GetPyName<ST>();
  bp::class_<ST>(name.c_str())
    .def("__str__", &ToString<ST>)
    .def("__len__", &ST::Size)
    .def("__contains__", &ST::Used)
    .def("GetName", FunctionPointer([](ST & self, int i) { return string(self.GetName(i)); }))
    .def("__getitem__", FunctionPointer([](ST & self, string name)
                                        {
                                          if (!self.Used(name)) bp::exec("raise KeyError()\n");
                                          return *self[name]; 
                                        }))
    .def("__getitem__", FunctionPointer ([](ST & self, int i)
                                         {
                                           if (i < 0 || i >= self.Size()) bp::exec("raise IndexError()\n");
                                           return *self[i];  
                                         }))
    ;
}  





#endif // NGS_PYTHON
#endif // PYTHON_NGSTD_HPP___
