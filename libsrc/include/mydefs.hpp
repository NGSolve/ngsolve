#ifndef FILE_MYDEFS
#define FILE_MYDEFS

/**************************************************************************/
/* File:   mydefs.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Mar. 98                                                    */
/**************************************************************************/

/*
  defines for graphics, testmodes, ...
*/


// #define DEBUG

// Philippose - 31/01/2009
// Hack for the Windows Version
// in Linux, "PACKAGE_VERSION" is replaced 
// in the configure/make phases, with the 
// right version number
#ifdef WIN32
#define PACKAGE_VERSION "4.9.6"
#endif


#define noDEMOVERSION
#define noDEVELOP
#define noSTEP
#define noSOLIDGEOM

#define noDEMOAPP
#define noMODELLER

#define noSTAT_STREAM
#define noLOG_STREAM

#endif
