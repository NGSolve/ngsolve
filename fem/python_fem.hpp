#ifndef PYTHON_FEM_HPP___
#define PYTHON_FEM_HPP___
#ifdef NGS_PYTHON

#include <core/register_archive.hpp>
// #include <fem.hpp>
#include <coefficient.hpp>
#include <finiteelement.hpp>

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
    static RegisterClassForArchive<cl_UnaryOpCF<FUNC>, CoefficientFunction> reguopcf;

    m.def(name.c_str(), py::vectorize([](double x)
    {
        FUNC func;
        return func(x);
    }), py::arg("x"), description.c_str());
    m.def(name.c_str(), py::vectorize([](Complex x)
    {
        FUNC func;
        return func(x);
    }), py::arg("x"), description.c_str());
    m.def(name.c_str(), [] (shared_ptr<CoefficientFunction> x) -> shared_ptr<CoefficientFunction> {
        FUNC func;
        return func(x);
      }, py::arg("x"), description.c_str());
  }
  

  template <typename FUNC>
  void ExportStdMathFunction2(py::module &m, string name, string description, string arg0="x", string arg1="y")
  {
    static RegisterClassForArchive<cl_BinaryOpCF<FUNC>, CoefficientFunction> regbinopcf;

    m.def(name.c_str(), py::vectorize([](double x, double y)
    {
        FUNC func;
        return func(x,y);
    }), py::arg(arg0.c_str()), py::arg(arg1.c_str()), description.c_str());
    m.def(name.c_str(), py::vectorize([](Complex x, Complex y)
    {
        FUNC func;
        return func(x,y);
    }), py::arg(arg0.c_str()), py::arg(arg1.c_str()), description.c_str());
    m.def(name.c_str(),
          [](shared_ptr<CoefficientFunction> x,
             shared_ptr<CoefficientFunction> y)
          -> shared_ptr<CoefficientFunction>
    {
      FUNC func;
      return BinaryOpCF(x, y, func, FUNC::Name());
    }, py::arg(arg0.c_str()), py::arg(arg1.c_str()), description.c_str());
  }


} // namespace ngfem

#endif // NGS_PYTHON
#endif // PYTHON_FEM_HPP___
