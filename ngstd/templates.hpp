nothing here anymore


#ifndef FILE_NGSTD_TEMPLATES
#define FILE_NGSTD_TEMPLATES

/*********************************************************************/
/* File:   templates.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

#include <iostream>
#include <ngs_stdcpp_include.hpp>

namespace ngstd 
{
  using namespace std;

    /*

  
/// min of 2 values
template <class T>
INLINE T min2 (T a, T b)
{
  return (a < b) ? a : b;
}

/// max of 2 values
template <class T>
INLINE T max2 (T a, T b)
{
  return (a > b) ? a : b;
}

/// min of 3 values
template <class T>
INLINE T min3 (T a, T b, T c)
{
  return (a < b) ? (a < c) ? a : c
    : (b < c) ? b : c;
}

/// max of 3 values
template <class T>
INLINE T max3 (T a, T b, T c)
{
  ///
  return (a > b) ? ((a > c) ? a : c)
    : ((b > c) ? b : c);
}
  */


/*
/// sign of value (+1, 0, -1)
template <class T>
INLINE int sgn (T a)
{
  return (a > 0) ? 1 : ( ( a < 0) ? -1 : 0 );
}
*/

  /*
/// square element 
template <class T>
INLINE T sqr (const T a)
{
  return a * a; 
}
  */

   /*
  using ngcore::sqr;
  using ngcore::pow3;
  using ngcore::max2;
  using ngcore::min2;
  using ngcore::max3;
  using ngcore::min3;
  using ngcore::LoadBin;
  using ngcore::SaveBin;
   */
  
    /*

/// element to the third power
template <class T>
INLINE T pow3 (const T a)
{
  return a * a * a; 
}


template <class T>
void SaveBin (ostream & ost, const T & val)
{
  const char * cp = reinterpret_cast<const char*> (&val);
  for (unsigned j = 0; j < sizeof(T); j++)
    ost.put(cp[j]);
}


template <class T>
void LoadBin (istream & ist, T & val)
{
  char * cp = reinterpret_cast<char*> (&val);
  for (unsigned j = 0; j < sizeof(T); j++)
    ist.get(cp[j]);
}
  */



}

#endif
