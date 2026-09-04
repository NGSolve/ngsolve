#ifndef FILE_TINYBLA
#define FILE_TINYBLA

/*********************************************************************/
/* File:   tinybla.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   3. Aug. 2026                                              */
/*********************************************************************/

#include <string>

// the tinybla kernel source (tinybla.h in the source tree),
// embedded at build time by embed_source.cmake
extern const std::string code_tinybla;

#endif
