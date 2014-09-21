/**************************************************************************/
/* File:   exception.cpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   16. Jan. 02                                                    */
/**************************************************************************/

#include <ngstd.hpp>

namespace ngstd
{
  Exception :: Exception ()
  { 
    ;
  }

  Exception :: Exception (const string & s) 
    : what(s)
  { 
    cout << "create ngstd::Exception, what = " << s << endl;
  }

  Exception :: Exception (const char * s) 
    : what(s)
  { 
    cout << "create ngstd::Exception, what = " << s << endl;
  }

  Exception :: ~Exception () 
  { ; }

  Exception & Exception :: Append (const string & s)
  { 
    what += s; 
    return *this;
  }

  Exception & Exception :: Append (const char * s)
  { 
    what += s; 
    return *this;
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
