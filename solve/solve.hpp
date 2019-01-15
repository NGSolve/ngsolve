#ifndef FILE_SOLVE
#define FILE_SOLVE

/*********************************************************************/
/* File:   solve.hpp                                                 */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


/* 
   NGS Solvers: BVP, Instationary, ...
*/


#include <ngstd.hpp>
#include <nginterface.h>

#include <comp.hpp>
#include <multigrid.hpp>

struct Tcl_Interp;


/*
#ifdef WIN32
   #define LOCAL_EXPORTS __declspec(dllexport)
#else
   #define LOCAL_EXPORTS 
#endif
*/

/**
   A collection of solvers.
   Numerical procedures NumProc, e.g. for
   the solution of boundary value problems NumProcBVP or
   eigen value problems NumProcEVP, and many many more.
 */
namespace ngsolve
{
  using namespace std;

  using namespace ngstd;
  using namespace ngla;
  using namespace ngfem;
  using namespace ngcomp;
  using namespace ngmg;
}


/*
#include "numproc.hpp"
#include "pde.hpp"
*/


/*
#ifdef WIN32 
// trick from http://social.msdn.microsoft.com/Forums/en/vclanguage/thread/ab642c88-2d2d-4f5d-9fd7-2341442d5a46
// all new/delete allocation is done from ngsolve heap

NGS_DLL_HEADER void * __cdecl my_operator_new_replacement(size_t _count);
NGS_DLL_HEADER void __cdecl my_operator_delete_replacement(void * _ptr);
NGS_DLL_HEADER void * __cdecl my_operator_new_array_replacement(size_t _count);
NGS_DLL_HEADER void __cdecl my_operator_delete_array_replacement(void * _ptr);

#ifndef NGS_EXPORTS
inline void * __cdecl operator new(size_t _count) {
    return my_operator_new_replacement(_count);
}
inline void __cdecl operator delete(void * _ptr) {
    my_operator_delete_replacement(_ptr);
}
inline void * __cdecl operator new[](size_t _count) {
    return my_operator_new_array_replacement(_count);
}
inline void __cdecl operator delete[](void * _ptr) {
    my_operator_delete_array_replacement(_ptr);
}

#endif

#endif
*/




#endif
