#ifndef FILE_MATKERNEL32
#define FILE_MATKERNEL32

#if NETGEN_DEFAULT_SIMD_SIZE==1
#include "matkernel32_1.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==2
#include "matkernel32_2.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==4
#include "matkernel32_4.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==8
#include "matkernel32_8.hpp"
#endif

#endif
