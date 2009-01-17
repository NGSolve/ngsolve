/**************************************************************************/
/* File:   exception.cpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   16. Jan. 02                                                    */
/**************************************************************************/

#include <ngstd.hpp>

namespace ngstd
{
  using namespace ngstd;


  Exception :: Exception (const string & s) 
    : what(s)
  { 
    cout << "create ngstd::Exception, what = " << s << endl;
  }

  Exception :: ~Exception () 
  { ; }

  void Exception :: Append (const string & s)
  { 
    what += s; 
  }




  RangeException :: RangeException (const string & where, 
				    int ind, int imin, int imax)
    : Exception ("")
  {
    stringstream str;
    str << where << ": index " << ind << " out of range [" << imin << "," << imax << "]\n";
    Append (str.str());
  }

}
