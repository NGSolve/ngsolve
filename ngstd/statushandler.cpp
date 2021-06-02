
#include <nginterface.h>
#include <ngstd.hpp>

namespace ngstd
{
  static mutex m;
  void BaseStatusHandler::PushStatus (const char * str) const
  {
    lock_guard lock(m);
    Ng_PushStatus(str);
  }

  void BaseStatusHandler::PopStatus () const
  {
    lock_guard lock(m);
    Ng_PopStatus();
  }

  void BaseStatusHandler::SetThreadPercentage (double percent) const
  {
    Ng_SetThreadPercentage(percent);
  }

  void BaseStatusHandler::GetStatus (string & str, double & percent) const
  {
    char * s;
    Ng_GetStatus(&s, percent);
    str = s;
  }

  void BaseStatusHandler::SetTerminate(void) const
  {
    Ng_SetTerminate();
  }

  void BaseStatusHandler::UnSetTerminate(void) const
  {
    Ng_UnSetTerminate();
  }

  bool BaseStatusHandler::ShouldTerminate(void) const
  {
    return (Ng_ShouldTerminate() != 0);
  }
} // namespace ngstd

