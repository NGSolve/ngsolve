#ifndef FILE_NGSSTREAM
#define FILE_NGSSTREAM

/**************************************************************************/
/* File:   ngsstream.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   20. Jul. 2011                                                  */
/**************************************************************************/


namespace ngstd
{

 
  // important message
  class IM
  {
    int value;
  public:
    IM (int val) : value(val) { ; }
    int Value () const { return value; }
  };

  class trunc
  {
    double eps;
  public:
    trunc (double aeps) : eps(aeps) { ; }
    double Eps() const { return eps; }
  };
  
  class NGSOStream
  {
    ostream & ost;
    bool active;
    NGS_DLL_HEADER static bool glob_active;
    double trunc;
  public:
    NGSOStream (ostream & aost, bool aactive)
      : ost(aost), active(aactive), trunc(-1) { ; }
    NGSOStream & SetTrunc (double atrunc) { trunc = atrunc; return *this; }
    double GetTrunc () const { return trunc; }
    bool Active () const { return active && glob_active; }
    ostream & GetStream () { return ost; }
    static void SetGlobalActive (bool b) { glob_active = b; }
  };
  
  inline NGSOStream operator<< (ostream & ost, const IM & im)
  {
    return NGSOStream (ost, 
		       (im.Value() <= printmessage_importance));
  }

  /*
    // doesn't work for matrices
  inline NGSOStream operator<< (ostream & ost, trunc tr)
  {
    cout << "set trunc modifier" << endl;
    return NGSOStream (ost, true).SetTrunc (tr.Eps());
  }
  */

  
  template <typename T>
  inline NGSOStream operator<< (NGSOStream ngsost, const T & data)
  {
    if (ngsost.Active())
      ngsost.GetStream() << data;
    return ngsost;
  }

  /*
  inline NGSOStream operator<< (NGSOStream ngsost, const double & data)
  {
    cout << "double out" << endl;
    if (ngsost.Active())
      {
        double hdata = data;
        if (fabs (hdata) < ngsost.GetTrunc()) hdata = 0.0;
        ngsost.GetStream() << hdata;
      }
    return ngsost;
  }
  */
  
  inline NGSOStream operator<< (NGSOStream ngsost, ostream& ( *pf )(ostream&))
  {
    if ( ngsost.Active() )
      ngsost.GetStream() << (*pf);

    return ngsost;
  }
  
  inline NGSOStream operator<< (NGSOStream ngsost, ios& ( *pf )(ios&))
  {
    if ( ngsost.Active() )
      ngsost.GetStream() << (*pf);
    
    return ngsost;
  }

  inline NGSOStream operator<< (NGSOStream ngsost, ios_base& ( *pf )(ios_base&))
  {
    if ( ngsost.Active() )
      ngsost.GetStream() << (*pf);
    
    return ngsost;
  }


}

#endif
