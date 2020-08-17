#ifndef _STATUSHANDLER_HPP
#define _STATUSHANDLER_HPP

namespace ngstd
{

  /** Access to statusbar. (and more)
   */

  class NGS_DLL_HEADER BaseStatusHandler
  {
  public:
    virtual ~BaseStatusHandler () { ; }
    virtual void PushStatus (const char * str) const;
    virtual void PopStatus () const;
    virtual void SetThreadPercentage (double percent) const;

    virtual void GetStatus (string & str, double & percent) const;

    virtual void SetTerminate(void) const;
    virtual void UnSetTerminate(void) const;
    virtual bool ShouldTerminate(void) const;
  };

}

#endif // _STATUSHANDLER_HPP
