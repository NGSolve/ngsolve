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


#ifdef __GNUC__
#if( __GNUC__==8 && __GNUC_MINOR__<=2)
// gcc 8.1/8.2 procudes empty loops in code like
// for (auto vb : {VOL, BND})
#error "This code does not compile with GCC 8.1/8.2, please upgrade your compiler"
#endif
#endif

#if defined(__AVX512F__)
#ifdef __GNUC__
#if( __GNUC__==9 && __GNUC_MINOR__<=2)
// gcc 9.1/9.2 generates wrong code on avx512_skylake platforms:
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=93009
#error "GCC 9.1/9.2 generates wrong code on AVX512 platforms (see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=93009 ). Either build for a different architecture (cmake -DUSE_NATIVE_ARCH=OFF), or use a different compiler (like GCC 8.3 or Clang)"
#endif
#endif
#endif




#ifdef __INTEL_COMPILER
// #pragma warning (disable:175)    // range check 
// #pragma warning (disable:597)    // implicit conversion (2014)

// unknown attribute __leaf__ in /usr/include/x86_64-linux-gnu/sys/sysmacros.h
#pragma warning (disable:1292)  
#endif


#ifdef __clang__
// #pragma clang diagnostic ignored "-Wpotentially-evaluated-expression"
// related to clang Bug 2181:  https://www.nsnam.org/bugzilla/show_bug.cgi?id=2181
#endif

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wvirtual-move-assign"
#pragma GCC diagnostic ignored "-Wattributes"
// this one silences warning: requested alignment 4096 is larger than 256
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
#include <typeindex>
#include <typeinfo>
#include <memory>
#include <initializer_list>
#include <functional>
#include <atomic>
#include <mutex>
#include <list>
#include <array>
#include <optional>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#ifdef __CUDACC__
#define CUDA
#define HD __host__ __device__
#endif

#ifndef HD
#define HD
#endif


#ifdef __INTEL_COMPILER
#ifdef WIN32
#define ALWAYS_INLINE __forceinline
#define INLINE __forceinline inline
#define LAMBDA_INLINE
#else
#define ALWAYS_INLINE __forceinline
#define INLINE __forceinline inline
#define LAMBDA_INLINE __attribute__ ((__always_inline__))
#endif
#else
#ifdef __GNUC__
#define ALWAYS_INLINE __attribute__ ((__always_inline__))
#define INLINE __attribute__ ((__always_inline__)) inline   HD 
#define LAMBDA_INLINE __attribute__ ((__always_inline__))
// #ifndef __clang__
#define VLA
// #endif
#else
#define ALWAYS_INLINE
#define INLINE inline
#define LAMBDA_INLINE
#endif
#endif

// from https://stackoverflow.com/questions/60802864/emulating-gccs-builtin-unreachable-in-visual-studio
#ifdef __GNUC__ // GCC 4.8+, Clang, Intel and other compilers compatible with GCC (-std=c++0x or above)
[[noreturn]] inline __attribute__((always_inline)) void unreachable() {__builtin_unreachable();}
#elif defined(_MSC_VER) // MSVC
[[noreturn]] __forceinline void unreachable() {__assume(false);}
#else // ???
inline void unreachable() {}
#endif


#ifndef __assume
#ifdef __GNUC__
#ifdef __clang__
#define __assume(cond) __builtin_assume(cond)
#else
#define __assume(cond) if (!(cond)) __builtin_unreachable(); else;
#endif
#else
#define __assume(cond)
#endif
#endif

//#define INLINE __attribute__ ((__always_inline__)) inline
//#define INLINE inline


// #ifdef __clang__
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
namespace std
{
  // avoid expensive call to complex mult by using the grammar school implementation
  INLINE std::complex<double> operator* (std::complex<double> a, std::complex<double> b)
  {
    return std::complex<double> (a.real()*b.real()-a.imag()*b.imag(),
                                 a.real()*b.imag()+a.imag()*b.real());
  }
}
#endif

#ifdef PARALLEL
#include <unistd.h>  // for usleep (only for parallel)
#endif
