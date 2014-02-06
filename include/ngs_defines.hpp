#ifndef FILE_NGS_DEFINES
#define FILE_NGS_DEFINES

/**************************************************************************/
/* File:   ngs_defines.hpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   21. Feb. 03                                                    */
/**************************************************************************/

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
#pragma warning (disable:175)    // range check 
#pragma warning (disable:597)    // implicit conversion (2014)
#endif



// performs range-checking
// #define DEBUG


// maximal system dimension (3D elasticity = 3, piezo = 4)
// 8 for multiharmonic

#ifndef MAX_SYS_DIM
#define MAX_SYS_DIM 3
#endif

#ifndef MAX_CACHEBLOCKS
#define MAX_CACHEBLOCKS 0
#endif



#endif
