#ifndef FILE_AUTOPTR
#define FILE_AUTOPTR

/**************************************************************************/
/* File:   autoptr.hpp                                                    */
/* Author: STL, Joachim Schoeberl                                         */
/* Date:   29. Dec. 02                                                    */
/**************************************************************************/

namespace ngstd
{

/**
   Pointer to object.
   The object is deleted at the end of the scope of the AutoPtr
 */
template <typename T>
class AutoPtr
{
private:
  /// the pointer
  T * ptr;
public:
  ///
  typedef T* pT;

  /// initialize AutoPtr
  explicit AutoPtr (T * p = 0)  { ptr = p; }
  /// delete object
  ~AutoPtr () { delete ptr; }
  
  /// reference to object
  T & operator*() const { return *ptr; }

  /// reference to object
  T* operator->() const { return ptr; }

  /// reference to the pointer
  T *& Ptr() { return ptr; }

  /// reference to the pointer
  T * Ptr() const { return ptr; }

  /// delete object, and reset pointer
  void Reset(T * p = 0) { if (p != ptr) { delete ptr; ptr = p; } }

  /// is pointer non-zero ?
  operator bool () { return ptr != 0; }

private:
  /// forbid copy constructor
  AutoPtr (AutoPtr &) { ; }
  /// forbid assignement operator
  AutoPtr & operator= (AutoPtr &) { return *this; }
};

}

#endif
