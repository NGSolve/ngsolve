// a wrapper to load netgen-dll into python

#include <iostream>
#include <../general/ngpython.hpp>

#ifdef WIN32
#define DLL_HEADER __declspec(dllimport)
#else
#define DLL_HEADER
#endif


void DLL_HEADER ExportNetgenMeshing(py::module &m);
void DLL_HEADER ExportMeshVis(py::module &m);
void DLL_HEADER ExportCSG(py::module &m);
void DLL_HEADER ExportCSGVis(py::module &m);
void DLL_HEADER ExportGeom2d(py::module &m);

PYBIND11_PLUGIN(libngpy)
{
    py::module ngpy("libngpy", "pybind netgen module");
    py::module csg = ngpy.def_submodule("_csg", "pybind csg module");
    ExportCSG(csg);
    py::module csgvis = ngpy.def_submodule("csgvis", "pybind csgvis module");
    ExportCSGVis(csgvis);
    py::module meshing = ngpy.def_submodule("_meshing", "pybind meshing module");
    ExportNetgenMeshing(meshing);
    py::module meshvis = ngpy.def_submodule("meshvis", "pybind meshvis module");
    ExportMeshVis(meshvis);
    py::module geom2d = ngpy.def_submodule("_geom2d", "pybind geom2d module");
    ExportGeom2d(geom2d);
    return ngpy.ptr();
}

// Force linking libnglib to libnetgenpy
namespace netgen
{
   void MyBeep (int i);
   void MyDummyToForceLinkingNGLib()
   {
       MyBeep(0);
   }
}
