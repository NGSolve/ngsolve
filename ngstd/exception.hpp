#ifndef FILE_EXCEPTION
#define FILE_EXCEPTION

/**************************************************************************/
/* File:   exception.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   16. Jan. 2002                                                  */
/**************************************************************************/

namespace ngstd
{

#ifdef DEBUG
#define CHECK_RANGE
#endif



/// Base class for all ng exceptions
class NGS_DLL_HEADER Exception 
{
  /// a verbal description of the exception
  string what;
public:
  /// string s describes the exception
  Exception (const string & s);
  ///
  virtual ~Exception ();

  /// append string to description
  void Append (const string & s);
  /// verbal description of exception
  const string & What() const { return what; }
};



/// Out of range exception used for arrays, vectors and matrices
class NGS_DLL_HEADER RangeException : public Exception
{
public:
  /// where it occurs, index, minimal and maximal indices
  RangeException (const string & where, 
		  int ind, int imin, int imax);
};

}

#endif
