#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef VERSION
#define VERSION "5.3-dev"
#endif

#define NGS_VERSION 530


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


#ifndef HD
#define HD
#endif

#ifdef __INTEL_COMPILER
#define ALWAYS_INLINE __forceinline
#define INLINE __forceinline inline
#define LAMBDA_INLINE
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
