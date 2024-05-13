#ifndef NGS_STDCPP_INCLUDE_HPP
#define NGS_STDCPP_INCLUDE_HPP

#include <core/ngcore_api.hpp>



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
#include <limits>
#include <cstring>

#include <new>
#include <exception>
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
#include <variant>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// __host__ __device__ for CUDA
#define HD NETGEN_HD

#define ALWAYS_INLINE NETGEN_ALWAYS_INLINE
#define INLINE NETGEN_INLINE
#define LAMBDA_INLINE NETGEN_LAMBDA_INLINE

#ifdef NETGEN_VLA
#define VLA
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


#ifdef NGS_EXPORTS
  #define NGS_DLL_HEADER NGCORE_API_EXPORT
#else
  #define NGS_DLL_HEADER NGCORE_API_IMPORT
#endif

#endif
