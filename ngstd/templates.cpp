/*********************************************************************/
/* File:   templates.hh                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


#include <ngstd.hpp>

namespace ngstd
{
  using namespace ngstd;
  
  template <class T>
  void MergeSort (int size, T * data, T * help)
  {
    if (size < 20)
      {
	BubbleSort (size, data);
	return;
      }
   

    int i, i1, i2;

    int s1 = size/2;
    int s2 = size - s1;

    T * p1 = help;
    T * p2 = help + s1;


    for (i = 0; i < size; i++)
      help[i] = data[i];
  
    MergeSort (s1, p1, data);
    MergeSort (s2, p2, data);


    i1 = 0;
    i2 = 0;
    i = 0;
  
    while (i1 < s1 && i2 < s2)
      {
	if (p1[i1] <= p2[i2])
	  {
	    data[i] = p1[i1];
	    i1++;
	  }
	else
	  {
	    data[i] = p2[i2];
	    i2++;
	  }
	i++;
      }
  
    while (i1 < s1)
      {
	data[i] = p1[i1];
	i++;
	i1++;
      }
    while (i2 < s2)
      {
	data[i] = p2[i2];
	i++;
	i2++;
      }
  }

  template <class T>
  void MergeSort (int size, T * data, T * help, int * index, int * indexhelp)
  {
    if (size < 20)
      {
	BubbleSort (size, data, index);
	return;
      }
   

    int i, i1, i2;

    int s1 = size/2;
    int s2 = size - s1;

    T * p1 = help;
    T * p2 = help + s1;

    int * p1_index = indexhelp;
    int * p2_index = indexhelp + s1;


    for (i = 0; i < size; i++)
      help[i] = data[i];
    for (i = 0; i < size; i++)
      indexhelp[i] = index[i];
  
    MergeSort (s1, p1, data, p1_index, index);
    MergeSort (s2, p2, data, p2_index, index);


    i1 = 0;
    i2 = 0;
    i = 0;
  
    while (i1 < s1 && i2 < s2)
      {
	if (p1[i1] <= p2[i2])
	  {
	    data[i] = p1[i1];
	    index[i] = p1_index[i1];
	    i1++;
	  }
	else
	  {
	    data[i] = p2[i2];
	    index[i] = p2_index[i2];
	    i2++;
	  }
	i++;
      }
  
    while (i1 < s1)
      {
	data[i] = p1[i1];
	index[i] = p1_index[i1];
	i++;
	i1++;
      }
    while (i2 < s2)
      {
	data[i] = p2[i2];
	index[i] = p2_index[i2];
	i++;
	i2++;
      }
  }

  /*
  template <class T>
  void BubbleSort (int size, T * data)
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
  */

  template  void MergeSort<int> (int, int*, int*);
  template  void MergeSort<double> (int, double*, double*, int*, int*);
  // template  void BubbleSort<int> (int, int*);
}
