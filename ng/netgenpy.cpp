#include <iostream>
#include <boost/python.hpp>

void __declspec(dllimport) ExportNetgenMeshing();
void __declspec(dllimport) ExportMeshVis();
void __declspec(dllimport) ExportCSG();
void __declspec(dllimport) ExportCSGVis();
void __declspec(dllimport) ExportGeom2d();

BOOST_PYTHON_MODULE(nglib) 
{
    ExportCSG();
    ExportCSGVis();
    ExportNetgenMeshing();
    ExportMeshVis();
    ExportGeom2d();
}