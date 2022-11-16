#if NETGEN_DEFAULT_SIMD_SIZE==1
#include "matkernel_1.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==2
#include "matkernel_2.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==4
#include "matkernel_4.hpp"
#elif NETGEN_DEFAULT_SIMD_SIZE==8
#include "matkernel_8.hpp"
#endif
