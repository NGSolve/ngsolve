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
inline T min2 (T a, T b)
{
  return (a < b) ? a : b;
}

/// max of 2 values
template <class T>
inline T max2 (T a, T b)
{
  return (a > b) ? a : b;
}

/// min of 3 values
template <class T>
inline T min3 (T a, T b, T c)
{
  return (a < b) ? (a < c) ? a : c
    : (b < c) ? b : c;
}

/// max of 3 values
template <class T>
inline T max3 (T a, T b, T c)
{
  ///
  return (a > b) ? ((a > c) ? a : c)
    : ((b > c) ? b : c);
}


/// swap 2 elements. 
template <class T>
inline void Swap (T & a, T & b)
{
  T temp = a;
  a = b;
  b = temp;
}

/*
template <class T>
inline void swap (T & a, T & b)
{
  T temp = a;
  a = b;
  b = temp;
}
*/

/// sign of value (+1, 0, -1)
template <class T>
inline int sgn (T a)
{
  return (a > 0) ? 1 : ( ( a < 0) ? -1 : 0 );
}

/// square element 
template <class T>
inline T sqr (const T a)
{
  return a * a; 
}

/// element to the third power
template <class T>
inline T pow3 (const T a)
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


/// bubble sort array
template <class T>
inline void BubbleSort (int size, T * data)
{
  T hv;
  for (int i = 0; i < size; i++)
    for (int j = i+1; j < size; j++)
      if (data[i] > data[j])
	{
	  hv = data[i];
	  data[i] = data[j];
	  data[j] = hv;
	}
}
/// bubble sort array
template <class T>
inline void BubbleSort (int size, T * data, int * index)
{
  T hv;
  int hi;
  for (int i = 0; i < size; i++)
    for (int j = i+1; j < size; j++)
      if (data[i] > data[j])
	{
	  hv = data[i];
	  data[i] = data[j];
	  data[j] = hv;

	  hi = index[i];
	  index[i] = index[j];
	  index[j] = hi;
	}
}


/// merge sort array, use help array of same size
template <class T>
void MergeSort (int size, T * data, T * help);
/// merge sort array, use help array of same size
template <class T>
void MergeSort (int size, T * data, T * help, int * index, int * indexhelp);

}

#endif
