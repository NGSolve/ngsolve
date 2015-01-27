#ifdef WIN32

// a wrapper to load netgen-dll into the executable


#include <mydefs.hpp>

DLL_HEADER int NG_main(int argc, char ** argv);

int main(int argc, char ** argv) 
{
    return NG_main(argc, argv);
}


#endif

