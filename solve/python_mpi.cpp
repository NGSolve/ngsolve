not needed anymore

#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include "../ngstd/mpiwrapper.hpp"

using namespace ngstd;

void ExportNgmpi(py::module &m) {
    m.def("SendCommand", [] (string cmd) -> void {
            MyMPI_SendCmd(cmd.c_str());
      }, py::arg("cmd"));

    m.def("Rank", [] () {
        return  MyMPI_GetId(ngs_comm);
    });

    m.def("Barrier", [] () {
        MyMPI_Barrier(ngs_comm);
    }); 
}




#endif
