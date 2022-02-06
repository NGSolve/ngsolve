#ifndef PYTHON_FEM_HPP___
#define PYTHON_FEM_HPP___
#ifdef NGS_PYTHON

#include <core/register_archive.hpp>
#include <fem.hpp>

namespace ngfem
{
  shared_ptr<CoefficientFunction> MakeCoefficient (py::object val);
  Array<shared_ptr<CoefficientFunction>> MakeCoefficients (py::object py_coef);

  template <typename FUNC>
  void ExportStdMathFunction(py::module &m, string name, string description)
  {
    static RegisterClassForArchive<cl_UnaryOpCF<FUNC>, CoefficientFunction> reguopcf;

    m.def (name.c_str(), [name] (py::object x) -> py::object
              {
                FUNC func;
                py::extract<double> ed(x);
                if (ed.check()) return py::cast(func(ed()));
                if (py::extract<Complex> (x).check())
                  return py::cast(func(py::extract<Complex> (x)()));
                if (py::extract<shared_ptr<CoefficientFunction>>(x).check())
                  {
                    auto coef = py::extract<shared_ptr<CoefficientFunction>>(x)();
                    return py::cast(UnaryOpCF(coef, func, FUNC::Name()));
                  }
                throw py::type_error (string("can't compute math-function, type = ")
                                      + typeid(FUNC).name());
              }, py::arg("x"), description.c_str());
  }

  template <typename FUNC>
  void ExportStdMathFunction_(py::module &m, string name,
                              string description)
  {
    m.def (name.c_str(), [name] (py::object x) -> py::object
              {
                FUNC func;
                py::extract<double> ed(x);
                if (ed.check()) return py::cast(func(ed()));
                if (py::extract<Complex> (x).check())
                  return py::cast(func(py::extract<Complex> (x)()));
                if (py::extract<shared_ptr<CoefficientFunction>>(x).check())
                  {
                    auto coef = py::extract<shared_ptr<CoefficientFunction>>(x)();
                    return py::cast(func(coef));
                  }
                throw py::type_error (string("can't compute math-function, type = ") + name);
              }, py::arg("x"), description.c_str());
  }
  

  template <typename FUNC>
  void ExportStdMathFunction2(py::module &m, string name, string description, string arg0="x", string arg1="y")
  {
    static RegisterClassForArchive<cl_BinaryOpCF<FUNC>, CoefficientFunction> regbinopcf;

    m.def (name.c_str(),
           [name] (py::object x, py::object y) -> py::object
           {
             FUNC func;
             if (py::extract<shared_ptr<CoefficientFunction>>(x).check() || py::extract<shared_ptr<CoefficientFunction>>(y).check())
               {
                 shared_ptr<CoefficientFunction> cx = py::cast<shared_ptr<CoefficientFunction>>(x);
                 shared_ptr<CoefficientFunction> cy = py::cast<shared_ptr<CoefficientFunction>>(y);
                 return py::cast(BinaryOpCF(cx, cy, func, FUNC::Name()));
               }
             py::extract<double> dx(x), dy(y);
             if (dx.check() && dy.check()) return py::cast(func(dx(), dy()));
             py::extract<Complex> cx(x), cy(y);
             if (cx.check() && cy.check()) return py::cast(func(cx(), cy()));
             throw py::type_error (string("can't compute binary math-function")+typeid(FUNC).name());
           }, py::arg(arg0.c_str()), py::arg(arg1.c_str()), description.c_str());
  }


} // namespace ngfem

#endif // NGS_PYTHON
#endif // PYTHON_FEM_HPP___
