
// #define _GLIBCPP_BASIC_FILE_ENCAPSULATION 1


#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define VERSION "4.9.7"
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


/*** threading headers ***/
#ifdef _MSC_VER
# define WIN32_LEAN_AND_MEAN
# ifdef MSVC_EXPRESS
#  include <pthread.h>
# else
#  include <afxwin.h>
#  include <afxmt.h>
# endif // MSVC_EXPRESS
# include <windows.h>
# undef WIN32_LEAN_AND_MEAN
# include <winnt.h>
#else // Not using MSVC++
# include <pthread.h>
#endif // _MSC_VER





#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _OPENMP
#include <omp.h>
#endif



// gcc compiler 
#ifdef __GNUC__
#define ALWAYS_INLINE __attribute__ ((__always_inline__))
#else
#define ALWAYS_INLINE
#endif

