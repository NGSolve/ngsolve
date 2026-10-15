#ifndef FILE_NGS_DEFINES
#define FILE_NGS_DEFINES

/**************************************************************************/
/* File:   ngs_defines.hpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   21. Feb. 03                                                    */
/**************************************************************************/




// performs range-checking
// #define DEBUG


// maximal system dimension (3D elasticity = 3, piezo = 4)
// 8 for multiharmonic

#ifndef MAX_SYS_DIM
#define MAX_SYS_DIM 3
#endif

#if MAX_SYS_DIM > 8
#error "instantiation lists support MAX_SYS_DIM <= 8"
#endif


/*
  Instantiation lists over the system dimension:
    #define INST(N) template class JacobiPrecond<Mat<N,N,double>>;
    NGS_FOR_SYS_DIMS(INST)
    #undef INST
  calls INST(N) for N = 1 ... MAX_SYS_DIM, NGS_FOR_SYS_DIMS_FROM_2 for N = 2 ... MAX_SYS_DIM.
*/
#if MAX_SYS_DIM >= 1
#define NGS_IF_SYS_DIM_1(...) __VA_ARGS__
#else
#define NGS_IF_SYS_DIM_1(...)
#endif
#if MAX_SYS_DIM >= 2
#define NGS_IF_SYS_DIM_2(...) __VA_ARGS__
#else
#define NGS_IF_SYS_DIM_2(...)
#endif
#if MAX_SYS_DIM >= 3
#define NGS_IF_SYS_DIM_3(...) __VA_ARGS__
#else
#define NGS_IF_SYS_DIM_3(...)
#endif
#if MAX_SYS_DIM >= 4
#define NGS_IF_SYS_DIM_4(...) __VA_ARGS__
#else
#define NGS_IF_SYS_DIM_4(...)
#endif
#if MAX_SYS_DIM >= 5
#define NGS_IF_SYS_DIM_5(...) __VA_ARGS__
#else
#define NGS_IF_SYS_DIM_5(...)
#endif
#if MAX_SYS_DIM >= 6
#define NGS_IF_SYS_DIM_6(...) __VA_ARGS__
#else
#define NGS_IF_SYS_DIM_6(...)
#endif
#if MAX_SYS_DIM >= 7
#define NGS_IF_SYS_DIM_7(...) __VA_ARGS__
#else
#define NGS_IF_SYS_DIM_7(...)
#endif
#if MAX_SYS_DIM >= 8
#define NGS_IF_SYS_DIM_8(...) __VA_ARGS__
#else
#define NGS_IF_SYS_DIM_8(...)
#endif

#define NGS_FOR_SYS_DIMS_FROM_2(INST)                                   \
  NGS_IF_SYS_DIM_2(INST(2)) NGS_IF_SYS_DIM_3(INST(3))                   \
  NGS_IF_SYS_DIM_4(INST(4)) NGS_IF_SYS_DIM_5(INST(5))                   \
  NGS_IF_SYS_DIM_6(INST(6)) NGS_IF_SYS_DIM_7(INST(7))                   \
  NGS_IF_SYS_DIM_8(INST(8))

#define NGS_FOR_SYS_DIMS(INST)  NGS_IF_SYS_DIM_1(INST(1)) NGS_FOR_SYS_DIMS_FROM_2(INST)




#endif
