#ifndef FILE_XBOOL
#define FILE_XBOOL

/**************************************************************************/
/* File:   xbool.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   14. Nov. 07                                                    */
/**************************************************************************/


namespace ngstd
{
  // an extended bool with values false/maybe/true

  enum TMAYBE { maybe };
  
  class xbool
  {
    uint8_t state;

  public:
    xbool (bool b) : state(b ? 2 : 0) { ; }
    xbool (TMAYBE x) : state(1) { ; }
    xbool () = default;
    xbool (const xbool &) = default;
    
    xbool & operator= (bool b) { state = b ? 2 : 0; return *this; }
    xbool & operator= (TMAYBE x) { state = 1; return *this; }
    
    bool IsTrue () const { return state == 2; }
    bool IsMaybe () const { return state == 1; }
    bool IsFalse () const { return state == 0; }
  };
  
}

#endif

