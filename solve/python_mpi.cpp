#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <mpiwrapper.hpp>

using namespace ngstd;

void ExportNgmpi() {
    std::string nested_name = "ngmpi";
    if( bp::scope() )
      nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".ngmpi");
    
    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << "exporting ngmpi as " << nested_name << endl;
    bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
    parent.attr("ngmpi") = module ;

    bp::scope local_scope(module);

    bp::def("SendCommand", FunctionPointer( [] (string cmd) -> void {
            MyMPI_SendCmd(cmd.c_str());
        }));

    bp::def("Rank", FunctionPointer( [] () {
        return  MyMPI_GetId();
    }));

    bp::def("Barrier", FunctionPointer( [] () {
        MyMPI_Barrier();
    }));
}




#endif
