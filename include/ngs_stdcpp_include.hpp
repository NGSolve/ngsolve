#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef VERSION
#define VERSION "6.0"
#endif

#define NGS_VERSION 600





#ifdef WIN32

// This function or variable may be unsafe. Consider using _ftime64_s instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. See online help for details.
#pragma warning(disable:4996)
#pragma warning(disable:4244)

// multiple inheritance via dominance
#pragma warning(disable:4250)

// bool-int conversion
#pragma warning(disable:4800)

// C++ exception specification ignored except to indicate a function is not __declspec(nothrow)
#pragma warning(disable:4290)

// no suitable definition provided for explicit template instantiation request
#pragma warning(disable:4661)


// needs to have dll-interface to be used by clients of class
#pragma warning(disable:4251)
// why does this apply to inline-only classes ????

// size_t to int conversion:
#pragma warning(disable:4267)

#endif






#ifdef __INTEL_COMPILER
// #pragma warning (disable:175)    // range check 
// #pragma warning (disable:597)    // implicit conversion (2014)

// unknown attribute __leaf__ in /usr/include/x86_64-linux-gnu/sys/sysmacros.h
#pragma warning (disable:1292)  
#endif







#include <iostream>
#include <fstream>
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
#include <initializer_list>
#include <functional>



#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


#ifdef __CUDACC__
#define CUDA
#define HD __host__ __device__
#endif

#ifndef HD
#define HD
#endif


#ifdef __INTEL_COMPILER
#define ALWAYS_INLINE __forceinline
#define INLINE __forceinline inline
#define LAMBDA_INLINE __attribute__ ((__always_inline__))
#else
#ifdef __GNUC__
#define ALWAYS_INLINE __attribute__ ((__always_inline__))
#define INLINE __attribute__ ((__always_inline__)) inline   HD 
#define LAMBDA_INLINE __attribute__ ((__always_inline__))
#define VLA
#else
#define ALWAYS_INLINE
#define INLINE inline
#define LAMBDA_INLINE
#endif
#endif


//#define INLINE __attribute__ ((__always_inline__)) inline
//#define INLINE inline


#ifdef PARALLEL
#include <unistd.h>  // for usleep (only for parallel)
#endif
