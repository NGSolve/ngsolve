#ifdef NGS_PYTHON
#include "python_ngstd.hpp"

bool rst_docu = false;

PYBIND11_PLUGIN(ngs_docu_flag)
{
  py::module md("ngs_docu_flag","docu_flag_module");
  md.def("Setdocu", []() { rst_docu = true; });
  return md.ptr();
}

#endif //NGS_PYTHON
