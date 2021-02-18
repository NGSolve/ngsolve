#ifndef FILE_NGSTD_TEMPLATES
#define FILE_NGSTD_TEMPLATES

/*********************************************************************/
/* File:   templates.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngstd 
{

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


/// sign of value (+1, 0, -1)
template <class T>
INLINE int sgn (T a)
{
  return (a > 0) ? 1 : ( ( a < 0) ? -1 : 0 );
}

/// square element 
template <class T>
INLINE T sqr (const T a)
{
  return a * a; 
}

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





  /*
template <int NUM>
class Cl_Iterate
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC f)
  {
    Cl_Iterate<NUM-1>::Do(f);
    f(IC<NUM>());
  }
};

template <>
class Cl_Iterate<-1>
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC f)  { }
};

template <>
class Cl_Iterate<0>
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC f)  { f(IC<0>()); }
};

template <int NUM, typename FUNC>
INLINE void Iterate (FUNC f)
{
  Cl_Iterate<NUM-1>::Do(f);
}
  */


template <int NUM, typename FUNC>
INLINE void Iterate (FUNC f)
{
  if constexpr (NUM > 1) Iterate<NUM-1> (f);
  if constexpr (NUM >= 1) f(IC<NUM-1>());
}






template <int NUM>
class Cl_Switch
{
public:
  template <typename FUNC>
  static INLINE void Do (size_t nr, FUNC f)
  {
    if (nr == NUM)
      f(IC<NUM>());
    else
      Cl_Switch<NUM-1>::Do(nr, f);
  }
};
  
template <>
class Cl_Switch<-1>
{
public:
  template <typename FUNC>
  static INLINE void Do (size_t /* nr */, FUNC /* f */)  { }
};

template <>
class Cl_Switch<0>
{
public:
  template <typename FUNC>
  static INLINE void Do (size_t /* nr */, FUNC f)
  {
    // if (nr == 0)
    f(IC<0>());
  }
};

template <int NUM, typename FUNC>
INLINE void Switch (size_t nr, FUNC f)
{
  Cl_Switch<NUM-1>::Do(nr, f);
}





}

#endif
