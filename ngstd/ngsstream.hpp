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

  
  class NGSOStream
  {
    ostream & ost;
    bool active;
  public:
    NGSOStream (ostream & aost, bool aactive)
      : ost(aost), active(aactive) { ; }
    bool Active () const { return active; }
    ostream & GetStream () { return ost; }
  };
  
  inline NGSOStream operator<< (ostream & ost, const IM & im)
  {
    return NGSOStream (ost, 
		       (im.Value() <= printmessage_importance)  &&
		       (id == 0) );
  }


  
  template <typename T>
  inline NGSOStream operator<< (NGSOStream ngsost, const T & data)
  {
    if (ngsost.Active())
      ngsost.GetStream() << data;
    return ngsost;
  }
  
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
