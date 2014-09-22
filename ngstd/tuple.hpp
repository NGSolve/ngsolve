#ifndef FILE_NGS_TUPLE
#define FILE_NGS_TUPLE

/**************************************************************************/
/* File:   tuple.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Sep. 14                                                    */
/**************************************************************************/

namespace ngstd
{

template <typename ... ARGS> class Tuple 
{ 
public:
  int Size() const { return 0; }
};

template <typename HEAD, typename ... TAIL>
class Tuple<HEAD, TAIL...> : Tuple<TAIL...>
{
  typedef Tuple<TAIL...> BASE;
  HEAD head;
public:
  Tuple () { ; }
  Tuple (HEAD h, TAIL ... t) : Tuple<TAIL...> (t...), head(h) { ; }

  HEAD Head() const { return head; }
  Tuple<TAIL...> Tail() const { return *this; }

  int Size() const { return BASE::Size()+1; }
};

template <typename ... ARGS>
ostream & operator<< (ostream & ost, Tuple<ARGS...> tup)
{
  return ost;
}

template <typename FIRST, typename ... ARGS>
ostream & operator<< (ostream & ost, Tuple<FIRST, ARGS...> tup)
{
  ost << tup.Head() << ", " << tup.Tail();
  return ost;
}


template <typename ... ARGS>
Tuple<ARGS...> MakeTuple (ARGS ... args)
{
  return Tuple<ARGS...> (args...);
}












}


#endif
