#ifndef FILE_COMP_UTILS
#define FILE_COMP_UTILS

/*********************************************************************/
/* File:   comp_utils.hpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   18. Jul. 2009                                             */
/*********************************************************************/



namespace ngstd
{

  /*
  class IntRange : public Array<int>::AppendType<IntRange>
  {
    int first, next;
  public: 
    IntRange (int f, int n) : first(f), next(n) {;} 
    int First() const { return first; }
    int Next() const { return next; }
  };

  template <> template<>
  inline int Array<int>::Append (const AppendType<IntRange> & atrange)
  {
    const IntRange & range = static_cast<const IntRange &> (atrange);

    int oldsize = size;
    int s = range.Next() - range.First();
    
    SetSize (oldsize+s);

    for (int i = 0; i < s; i++)
      data[oldsize+i] = range.First()+i;

    return size;
  }
  */

  class IntRange
  {
    int first, next;
  public: 
    IntRange (int f, int n) : first(f), next(n) {;} 
    int First() const { return first; }
    int Next() const { return next; }
    int Size() const { return next-first; }
  };

  inline Array<int> & operator+= (Array<int> & array, const IntRange & range)
  {
    int oldsize = array.Size();
    int s = range.Next() - range.First();
    
    array.SetSize (oldsize+s);

    for (int i = 0; i < s; i++)
      array[oldsize+i] = range.First()+i;

    return array;
  }
  
}


#endif
