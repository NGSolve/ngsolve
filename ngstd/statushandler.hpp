#ifndef _STATUSHANDLER_HPP
#define _STATUSHANDLER_HPP

namespace ngstd
{
  
  /** Access to statusbar. (and more)
   */

  class NGS_DLL_HEADER BaseStatusHandler
    {
    public:
    static void PushStatus (const std::string& str);
    static void PopStatus ();
    static void SetThreadPercentage (double percent);
    
    static void GetStatus (string & str, double & percent);
    
    static void SetTerminate(void);
    static void UnSetTerminate(void);
    static bool ShouldTerminate(void);
    
    class Region
    {
    public:
      Region(const string& str) { PushStatus(str); }
      ~Region() { PopStatus(); }
    };
  };

}

#endif // _STATUSHANDLER_HPP
