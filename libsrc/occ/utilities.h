//  SALOME Utils : general SALOME's definitions and tools
//
//  Copyright (C) 2003  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS 
// 
//  This library is free software; you can redistribute it and/or 
//  modify it under the terms of the GNU Lesser General Public 
//  License as published by the Free Software Foundation; either 
//  version 2.1 of the License. 
// 
//  This library is distributed in the hope that it will be useful, 
//  but WITHOUT ANY WARRANTY; without even the implied warranty of 
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
//  Lesser General Public License for more details. 
// 
//  You should have received a copy of the GNU Lesser General Public 
//  License along with this library; if not, write to the Free Software 
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA 
// 
//  See http://www.opencascade.org/SALOME/ or email : webmaster.salome@opencascade.org 
//
//
//
//  File   : utilities.h
//  Author : Antoine YESSAYAN, Paul RASCLE, EDF
//  Module : SALOME
//  $Header: /cvs/netgen/netgen/libsrc/occ/utilities.h,v 1.3 2008/03/31 14:20:28 wabro Exp $

/* ---  Definition macros file to print informations if _DEBUG_ is defined --- */

#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <iostream>
#include <cstdlib>
// #include "SALOME_Log.hxx"

/* ---  INFOS is always defined (without _DEBUG_): to be used for warnings, with release version --- */

#define INFOS(msg)    {SLog->putMessage(*SLog<<__FILE__<<" ["<<__LINE__<<"] : "<<msg<<endl);}
#define PYSCRIPT(msg) {SLog->putMessage(*SLog<<"---PYSCRIPT--- "<<msg<<endl);}

/* --- To print date and time of compilation of current source --- */

#if defined ( __GNUC__ )
#define COMPILER		"g++" 
#elif defined ( __sun )
#define COMPILER		"CC" 
#elif defined ( __KCC )
#define COMPILER		"KCC" 
#elif defined ( __PGI )
#define COMPILER		"pgCC" 
#elif defined ( __alpha )
#define COMPILER		"cxx" 
#else
#define COMPILER		"undefined" 
#endif

#ifdef INFOS_COMPILATION
#error INFOS_COMPILATION already defined
#endif

#define INFOS_COMPILATION { \
			   SLog->putMessage(\
					   *SLog<<__FILE__<<" ["<< __LINE__<<"] : "\
					   << "COMPILED with " << COMPILER \
					   << ", " << __DATE__ \
					   << " at " << __TIME__ <<endl); }

#ifdef _DEBUG_

/* --- the following MACROS are useful at debug time --- */

#define MYTRACE *SLog << "- Trace " << __FILE__ << " [" << __LINE__ << "] : " 

#define MESSAGE(msg) {SLog->putMessage( MYTRACE <<msg<<endl<<ends); }
#define SCRUTE(var)  {SLog->putMessage( MYTRACE << #var << "=" << var <<endl<<ends); }

#define REPERE *SLog << "   --------------" << endl 
#define BEGIN_OF(msg) {REPERE;MYTRACE<<"Begin of: "     <<msg<<endl;REPERE;} 
#define END_OF(msg)   {REPERE;MYTRACE<<"Normal end of: "<<msg<<endl;REPERE;} 

#define HERE {cout<<flush ;cerr<<"- Trace "<<__FILE__<<" ["<<__LINE__<<"] : "<<flush ;}

#define INTERRUPTION(code) {HERE;cerr<<"INTERRUPTION return code= "<<code<< endl;std::exit(code);}

#ifndef ASSERT
#define ASSERT(condition) \
        if (!(condition)){HERE;cerr<<"CONDITION "<<#condition<<" NOT VERIFIED"<<endl;INTERRUPTION(1);}
#endif /* ASSERT */


#else /* ifdef _DEBUG_*/

#define HERE 
#define SCRUTE(var) {}
#define MESSAGE(msg) {}
#define REPERE
#define BEGIN_OF(msg) {}
#define END_OF(msg) {}

#define INTERRUPTION(code) {}

#ifndef ASSERT
#define ASSERT(condition) {}
#endif /* ASSERT */


#endif /* ifdef _DEBUG_*/

#endif /* ifndef UTILITIES_H */
