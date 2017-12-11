#ifndef FILE_NGS_AUTOPTR
#define FILE_NGS_AUTOPTR

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
  bool owner;
public:
  ///
  typedef T* pT;

  /// initialize AutoPtr
  explicit AutoPtr (T * p = 0, bool aowner = true)  { ptr = p; owner = aowner;}
  /// delete object
  ~AutoPtr () { if (owner) delete ptr; }
  
  /// reference to object
  T & operator*() const { return *ptr; }

  /// reference to object
  T* operator->() const { return ptr; }

  /// reference to the pointer
  T *& Ptr() { return ptr; }

  /// reference to the pointer
  T * Ptr() const { return ptr; }

  /// delete object, and reset pointer
  void Reset(T * p = 0, bool aowner = true) 
  { if (p != ptr) { if (owner) delete ptr; ptr = p; }; owner = aowner; }

  /// is pointer non-zero ?
  operator bool () { return ptr != 0; }

private:
  /// forbid copy constructor
  AutoPtr (AutoPtr &) { ; }
  /// forbid assignment operator
  AutoPtr & operator= (AutoPtr &) { return *this; }
};

}

#endif
