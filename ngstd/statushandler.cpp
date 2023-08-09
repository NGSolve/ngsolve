
#include <nginterface.h>
#include <ngstd.hpp>

namespace ngstd
{
  static mutex m;
  void BaseStatusHandler::PushStatus (const char * str)
  {
    lock_guard lock(m);
    Ng_PushStatus(str);
  }

  void BaseStatusHandler::PopStatus ()
  {
    lock_guard lock(m);
    Ng_PopStatus();
  }

  void BaseStatusHandler::SetThreadPercentage (double percent)
  {
    Ng_SetThreadPercentage(percent);
  }

  void BaseStatusHandler::GetStatus (string & str, double & percent)
  {
    char * s;
    Ng_GetStatus(&s, percent);
    str = s;
  }

  void BaseStatusHandler::SetTerminate(void)
  {
    Ng_SetTerminate();
  }

  void BaseStatusHandler::UnSetTerminate(void)
  {
    Ng_UnSetTerminate();
  }

  bool BaseStatusHandler::ShouldTerminate(void)
  {
    return (Ng_ShouldTerminate() != 0);
  }

} // namespace ngstd

