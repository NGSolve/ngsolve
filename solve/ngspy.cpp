#ifdef NGS_PYTHON

#undef NGS_EXPORTS

#include <solve.hpp>
#include<boost/python.hpp>

void NGS_DLL_HEADER ExportNgstd();
void NGS_DLL_HEADER ExportNgbla();
void NGS_DLL_HEADER ExportNgfem();
void NGS_DLL_HEADER ExportNgla();
void NGS_DLL_HEADER ExportNgcomp();
void NGS_DLL_HEADER ExportNgsolve();


BOOST_PYTHON_MODULE(ngslib)
{

    try
    {
        ExportNgstd();
        ExportNgbla();
        ExportNgfem();
        ExportNgla();
        ExportNgcomp();      
        ExportNgsolve();
    }
    catch (ngstd::Exception & e)
    {
        std::cerr << "\n\nCaught Exception:\n" << e.What() << std::endl;
    }
    catch (...)
    {
        std::cerr << "\n\nCaught Python Exception:\n" << std::endl;
        PyErr_Print();
    }
}

#endif
