
// #define _GLIBCPP_BASIC_FILE_ENCAPSULATION 1


#ifdef HAVE_CONFIG_H
#include <../config.h>
#endif


#include <iostream>
#include <fstream>
// #include <strstream>
#include <sstream>
#include <iomanip>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <climits>
#include <cstring>

#include <new>
#include <exception>
#include <complex>
#include <string>
#include <typeinfo>
#include <memory>

//#ifndef MACOS
//#include <emmintrin.h>
//#endif

/*** threading headers ***/
#ifdef _MSC_VER
# define WIN32_LEAN_AND_MEAN
# include <afxwin.h>
# include <afxmt.h>
# include <windows.h>
# undef WIN32_LEAN_AND_MEAN
# include <winnt.h>
#else
# include <pthread.h>
#endif 





#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/*
extern void* operator new(std::size_t) throw (std::bad_alloc);
extern void* operator new[](std::size_t) throw (std::bad_alloc);
extern void operator delete(void*) throw();
extern void operator delete[](void*) throw();

extern int mem_total_alloc_vector;
*/
