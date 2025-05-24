#ifndef PYTHON_NGSTD_HPP___
#define PYTHON_NGSTD_HPP___
#ifdef NGS_PYTHON

#ifdef __clang__
#pragma clang diagnostic push
// #pragma clang diagnostic ignored "-W#pragma-messages"
#pragma clang diagnostic ignored "-Wunused-local-typedefs"
#pragma clang diagnostic ignored "-Wparentheses-equality"
#pragma clang diagnostic ignored "-Wunused-value"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/eval.h>
#include <pybind11/operators.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include<pybind11/numpy.h>

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#include <core/python_ngcore.hpp>
#include <ngstd.hpp>
#include <thread>
#include <iostream>

namespace py = pybind11;
using namespace pybind11::literals;

using std::string;
using std::cout;
using std::endl;

namespace ngstd {
  extern bool have_numpy;

  template<typename TClass, typename TFunc, typename... T>
  TClass &  PyDefVectorized(TClass & cls, const char * name, TFunc && f, T && ... args )
  {
      if(have_numpy)
          cls.def(name, py::vectorize(std::forward<TFunc>(f)), std::forward<T>(args)...);
      else
          cls.def(name, std::forward<TFunc>(f), std::forward<T>(args)...);
      return cls;
  }

  template<typename T, typename TOffset>
  auto getNumpyFieldDescriptor(const char *name, TOffset offset) {
      return py::detail::field_descriptor {
                          name, static_cast<py::ssize_t>(offset), sizeof(T),
                          py::format_descriptor<T>::format(),
                          py::detail::npy_format_descriptor<T>::dtype()
      };
  }
}

using namespace ngstd;

namespace pybind11 {
template <typename T>
struct extract
{
  py::handle obj;
  extract( py::handle aobj ) : obj(aobj) {}

  bool check() {
      try
      {
        obj.cast<T>();
      }
      catch(const std::runtime_error &e)
      {
        return false;
      }
      return true;
  }

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
};



typedef py::gil_scoped_acquire AcquireGIL; 
typedef py::gil_scoped_release ReleaseGIL; 


inline void InitSlice( const py::slice &inds, size_t len, size_t &start, size_t &step, size_t &n ) {
      size_t stop;
      if (!inds.compute(len, &start, &stop, &step, &n))                                          
        throw py::error_already_set();
}

template<typename T>
py::array_t<T> MoveToNumpyArray( Array<T> &a )
{
  if(a.Size()) {
      py::capsule free_when_done(&a[0], [](void *f) {
                                 delete [] reinterpret_cast<T *>(f);
                                 });
      a.NothingToDelete();
      return py::array_t<T>(a.Size(), &a[0], free_when_done);
  }
  else
      return py::array_t<T>(0, nullptr);
}


template<typename ... Args>
inline py::object PyRange(Args ... i)
{
  constexpr size_t n = sizeof...(Args);
  static_assert(n>=1 && n<=3, "Wrong number of arguments");
  return py::module::import("__main__").attr("__builtins__").attr("range")(i...);
}

inline py::object PyRange(IntRange i)
{
  return PyRange(i.First(), i.Next());
}

//////////////////////////////////////////////////////////////////////
template <typename T, typename TCLASS = py::class_<T> >
void PyDefToString( py::module &m, TCLASS &c )
{
    c.def("__str__", &ToString<T>);
    c.def("__repr__", &ToString<T>);
}

/*
template <typename T>
class cl_NonElement 
{
public:
  static T Val() { return 0; }
};

template <typename T>
inline T NonElement() { return cl_NonElement<T>::Val(); }
*/

//////////////////////////////////////////////////////////////////////
// read-only bracket operator
template< typename T, typename TELEM, typename TCLASS = py::class_<T>>
void PyDefROBracketOperator( py::module &m, TCLASS &c )
{
    auto Get = [](T& self, int i) {
      if (i < 0) i += self.Size();
      if( i<self.Size() && i>=0 )
        return self[i];
      throw py::index_error();
      return TELEM();
    }; 
    c.def("__getitem__", Get,py::arg("pos"), "Return value at given position");
    c.def("Get", Get, py::arg("pos"), "Return value at given position");
}

// read-write bracket operator
template< typename T, typename TELEM, typename TCLASS = py::class_<T>>
void PyDefBracketOperator( py::module &m, TCLASS &c )
{
    PyDefROBracketOperator<T, TELEM>(m, c);
    auto Set = [](T& self, int i, TELEM val) {
      if (i < 0) i += self.Size();      
      if( i<self.Size() && i>=0 )
        self[i] = val;
      else
        throw py::index_error();
    };
    c.def("__setitem__", Set, py::arg("pos"), py::arg("value"), "Set value at given position");
    c.def("Set", Set, py::arg("pos"), py::arg("value"), "Set value at given position");
}


//////////////////////////////////////////////////////////////////////
// Export len, bracket operator and iterator protocol at once
template <typename T, typename TELEM = double, typename TCLASS = py::class_<T> >
void PyDefVector( py::module &m, TCLASS &c )
{
  c.def("__len__",  []( T& v) { return v.Size();}, "Return length of the array"  );
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

template <typename T>
Table<T> makeCTable (py::list obj)
{
  size_t n = py::len(obj);
  Array<int> entrysize(n);
                       
  for (size_t i = 0; i < n; i++)
    entrysize[i] = py::len(obj[i]);
  
  Table<T> tab(entrysize);
  for (size_t i = 0; i < n; i++)
    {
      const py::object & obji = obj[i];
      tab[i] = makeCArray<T> (obji);
    }
  return tab;
}

template <typename T>
py::tuple MakePyTuple (const BaseArrayObject<T> & ao)
{
  size_t s = ao.Size();
  py::tuple tup(s);
  for (size_t i = 0; i < s; i++)
    tup[i] = ao[i];
  return tup;
}

template <typename T>
py::list MakePyList (const BaseArrayObject<T> & ao)
{
  size_t s = ao.Size();
  py::list l;
  for (size_t i = 0; i < s; i++)
    l.append (ao[i]);
  return l;
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
                                        }, py::arg("name"))
    .def("__getitem__", [](ST & self, int i) -> PY_T
                                         {
                                           if (i < 0 || i >= self.Size()) throw py::index_error();
                                           return self[i];  
                                         }, py::arg("pos"))
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
    .def("GetName", [](ST & self, int i) { return string(self.GetName(i)); }, py::arg("pos"))
    .def("__getitem__", [](ST & self, string name)
                                        {
                                          if (!self.Used(name)) throw py::index_error();
                                          return *self[name]; 
                                        }, py::arg("name"))
    .def("__getitem__", [](ST & self, int i)
                                         {
                                           if (i < 0 || i >= self.Size()) throw py::index_error();
                                           return *self[i];  
                                         }, py::arg("pos"))
    ;
}

// replace docu links with plain text for help function
NGS_DLL_HEADER const char* docu_string(const char* str);

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace ngstd
{
  // MemoryView for pickling without copying
  // doesn't provide memory management, it has to be done from the outside!
  class MemoryView
  {
  private:
    void* ptr;
    size_t size;

  public:
    MemoryView(void* aptr, size_t asize) : ptr(aptr), size(asize) { ; }

    size_t Size() const { return size; }
    void* Ptr() const { return ptr; }
  };
}

#endif // NGS_PYTHON
#endif // PYTHON_NGSTD_HPP___
